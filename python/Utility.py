# -----------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017-2018
# -----------------------------------------------------------------------------

#import pdb; pdb.set_trace()           # enables debugging

import glob, os, sys, errno, math, re, subprocess, traceback
import numpy as np
import scipy
import scipy.misc
import matplotlib
matplotlib.use("Agg")
from Configuration import Configuration


def MakePath(path) -> None:
    """ Takes directory or file name as an input
        and creates the path if it does not exist.
    """
    folder = os.path.dirname(path)
    if folder and (not os.path.exists(folder)):
        try:
            os.makedirs(folder)
            print("Directory created: " + folder)
        except OSError as ex:  # guard against race condition
            if ex.errno != errno.EEXIST:
                sys.exit("Failed to create folder: " + folder)


def MakeFileName(conf, what) -> str:
    """ Function returns the file name for: (1) sensor locations ("sensors"),
        (2) analytic solution ("analytic") or (3) true state field
        ("true_field") given configuration settings. Important, grid resolution
        and the number of time steps (except for sensor locations) are
        encrypted into the file name. This helps to distinguish simulations
        with different settings.
    """
    assert isinstance(conf, Configuration) and isinstance(what, str)
    assert hasattr(conf, 'Nt'), "number of time steps Nt is not specified"

    Nx = round(conf.num_subdomains_x * conf.subdomain_x)
    Ny = round(conf.num_subdomains_y * conf.subdomain_y)
    filename = os.path.join(conf.output_dir,
                            what + "_Nx" + str(Nx) + "_Ny" + str(Ny))
    if what == "sensors":
        filename = filename + ".txt"
    elif what == "analytic":
        filename = filename + "_Nt" + str(conf.Nt) + ".txt"
    elif what == "true_field":
        filename = filename + "_Nt" + str(conf.Nt) + ".bin"
    else:
        sys.exit("unknown entity to make a file name from")

    print("File name: " + filename)
    return filename


def CheckPythonVersion() -> None:
    """ Function checks for the minimum supported Python version.
    """
    if (sys.version_info[0] < 3 or
            (sys.version_info[0] == 3 and sys.version_info[1] < 5)):
        sys.exit("Python of version 3.5+ is expected")


def SwitchToGraphicalBackend() -> object:
    """ Function switches the default "Agg" graphical backend
        to another one capable to show up a window for plotting.
    """
    gui_env = ["TKAgg", "GTKAgg", "Qt4Agg", "WXAgg"]
    ok = False
    for gui in gui_env:
        try:
            print("Trying '" + gui + "' gui ... ", end="", flush=True)
            matplotlib.use(gui, warn=False, force=True)
            from matplotlib import pyplot as plt
            ok = True
            print("found")
            break
        except Exception:
            continue
    if ok:
        print("Graphical backend: '" + str(matplotlib.get_backend()) + "'")
    else:
        print("WARNING: no graphical backend is available")
        return None
    return plt


def ProblemParametersFromFilename(filename, extract_Nt = True,
                                            prefix = None) -> [int,int,int]:
    """ Function extracts problem size (Nx,Ny) from the file name
        and optionally the number of time steps (Nt). The function is useful
        when preparing demonstration from simulation results.
    """
    assert isinstance(filename, str)
    assert (prefix is None) or isinstance(prefix, str)
    # Get the base file name (without path) and check its validity.
    basename = os.path.basename(filename)
    assert (prefix is None) or basename.startswith(prefix), (
                "base file name should start with the prefix: " + prefix)
    # Check file name has the pattern of a file of solution fields.
    assert re.search(r"\w+_Nx\d+_Ny\d+_Nt\d+\.bin", basename), (
                            "base file name does not match the pattern")
    # Extracts Nx, Ny and (optionally) Nt from the file name.
    digits = re.sub("[^0-9]", " ", basename)
    params = [int(s) for s in digits.split() if s.isdigit()]
    Nx, Ny, Nt = int(-1), int(-1), int(-1)
    if extract_Nt:
        assert len(params) == 3, (
                "expecting (Nx,Ny,Nt) in the file name: " + basename)
        Nx = int(params[0])
        Ny = int(params[1])
        Nt = int(params[2])
    else:
        assert len(params) == 2, (
                "expecting (Nx,Ny) in the file name: " + basename)
        Nx = int(params[0])
        Ny = int(params[1])
    return Nx, Ny, Nt


def ReadResultFile(filename) -> [np.ndarray, np.ndarray]:
    """ Function reads the binary file where a number of solutions (full state
        fields) are stacked one after another. The file format has 4-column
        layout (time, abscissa, ordinate, value), all values have float-32
        precision and records are not necessary sorted (in case of C++ output).
    """
    Nx, Ny, Nt = ProblemParametersFromFilename(filename, True, None)
    # Read the solution fields file and sort the records because the parallel
    # C++ application does not guarantee proper ordering. Actually, we do not
    # sort the data as such, rather get the sorting index array.
    with open(filename, "rb") as fid:
        data = np.fromfile(fid, dtype=np.float32)
    data = np.reshape(data, (-1,4))
    Nw = data.shape[0] // (Nx * Ny)             # number of written fields
    assert data.shape[0] == Nx * Ny * Nw, "wrong file size"
    idx = np.lexsort(np.rot90(data[:,0:3]))     # sort by {t,x,y} triples
    # Create the output fields and corresponding timestamps.
    timestamps = np.zeros((Nw,), dtype=int)
    fields = np.zeros((Nw, Nx, Ny), dtype=np.float32)
    # Expected layouts of abscissas and ordinates.
    xpos = np.repeat(np.arange(Nx), Ny).astype(int)
    ypos = np.tile(np.arange(Ny), Nx).astype(int)
    for i in range(Nw):
        # Get indices of a block of records corresponding to a separate field.
        block = idx[i*Nx*Ny : (i+1)*Nx*Ny]
        # Within the block all timestamps are the same.
        timestamps[i] = np.rint(data[block[0], 0])
        assert np.all(np.rint(data[block, 0]) == timestamps[i])
        # Check the layouts of abscissas and ordinates.
        assert np.all(np.rint(data[block, 1]) == xpos)
        assert np.all(np.rint(data[block, 2]) == ypos)
        # Make 2D field from the last column.
        fields[i,:,:] = np.reshape(data[block, 3], (Nx, Ny))
    return timestamps, fields


class PrintProgress:
    """ Class prints a progress-line in terminal window.
    """
    def __init__(self, N):
        """ Constructor sets the range of time/iteration index.
        """
        self.N = round(N)
        assert N > 0

    def Print(self, n) -> None:
        """ Prints progress. The input parameter should fit the range [0..N).
        """
        N = self.N
        assert isinstance(n, int) and 0 <= n and n < N
        if n == 0 or n + 1 == N or (100 * (n - 1)) // N != (100 * n) // N:
            print(".", end="", flush=True)

    def Finalize(self) -> None:
        """ Prints the final new line.
        """
        print("", flush=True)


def Field2Image(field) -> object:
    """ Function converts state field to image. Here we transpose the field
        because Python assumes abscissas in the second dimension. Also, the
        image is flipped vertically so that the origin goes to the left-bottom
        corner.
    """
    assert len(field.shape) == 2
    # Image <- field scaled to [0..1] range.
    vmin = np.amin(field)
    vmax = np.amax(field)
    image = np.abs((np.transpose(field) - vmin) /
                   (vmax - vmin + math.sqrt(np.finfo(float).eps)))
    image = np.flipud(image)            # flip Y
    return image

