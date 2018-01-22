# -----------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017
# -----------------------------------------------------------------------------

import pdb; pdb.set_trace()           # enables debugging

import glob, os, sys, errno, math, re, subprocess, traceback
import numpy as np
import scipy
import scipy.misc
import matplotlib
matplotlib.use("Agg")
from Configuration import Configuration


def RemoveFile(filename) -> None:
    """ Check if a file exists on disk. If exists, delete it.
    """
    if os.path.exists(filename):
        try:
            os.remove(filename)
            print("File removed: " + filename)
        except OSError as ex:
            sys.exit("Error: {} : {}".format(ex.filename, ex.strerror))


def AssurePathExists(path) -> None:
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


def MakeBaseFileName(conf, entity) -> str:
    """ Function returns the file name for sensor locations, analytic solution,
        simulation (solution) field or true field given configuration settings.
        The term "base" means user should supply the proper suffix and
        extension, while the path and the base file name is formed
        by this function.
    """
    assert isinstance(conf, Configuration)
    assert entity == "sensors"  or \
           entity == "analytic" or \
           entity == "field"    or \
           entity == "true-field", \
            "allowed entity strings: " \
            "{'sensors', 'analytic', 'field', 'true-field'}"
    Nx = round(conf.num_subdomains_x * conf.subdomain_x)
    Ny = round(conf.num_subdomains_y * conf.subdomain_y)
    assert hasattr(conf, 'Nt'), "number of time steps Nt is not specified"
    filename = os.path.join(conf.output_dir, entity +
                                "_Nx" + str(Nx) + "_Ny" + str(Ny))
    if entity != "sensors":
        filename = filename + "_Nt" + str(conf.Nt)
    return filename


def MakeFileName(conf, entity, suffix=None) -> str:
    """ Function returns the file name for sensor locations,
        analytic solution, simulation solution or true field
        given configuration settings.
    """
    assert isinstance(conf, Configuration)
    assert entity == "sensors"  or \
           entity == "analytic" or \
           entity == "solution" or \
           entity == "true-field", \
            "allowed entity strings: " \
            "{'sensors', 'analytic', 'simulation', 'true-field'}"
    assert suffix is None or isinstance(suffix, str)
    Nx = round(conf.num_subdomains_x * conf.subdomain_x)
    Ny = round(conf.num_subdomains_y * conf.subdomain_y)
    assert hasattr(conf, 'Nt'), "number of time steps Nt is not specified"
    filename = os.path.join(conf.output_dir, entity +
                                "_Nx" + str(Nx) + "_Ny" + str(Ny))
    if entity != "sensors":
        filename = filename + "_Nt" + str(conf.Nt)

    if suffix is None:
        #assert entity != "true-field", (
        #        "file suffix is always expected for the true field")
        filename = filename + ".txt"
    else:
        if suffix.endswith(".avi"):
            filename = filename + suffix
        elif suffix.endswith(".png") or suffix.endswith(".txt"):
            filename = filename + "_" + suffix
        elif suffix == "*.txt":
            filename = filename + suffix
        else:
            sys.exit("file suffix does not have expected extension")

    print("File name: " + filename)
    return filename


#def ClearOutputDir(conf, ext=None) -> None:
    #""" Function clears the output directory, if it does exist.
        #Optionally, one can remove files of specified type. Note,
        #the extension argument should start from '.', e.g. '.txt'.
    #"""
    #assert isinstance(conf, Configuration)
    #if os.path.exists(conf.output_dir):
        #if ext is None:
            #filelist = [f for f in os.listdir(conf.output_dir)
                        #if f.endswith(".avi") or
                           #f.endswith(".pgm") or
                           #f.endswith(".png") or
                           #f.endswith(".jpg") or
                           #f.endswith(".txt")]
        #else:
            #assert isinstance(ext, str) and ext and ext.startswith(".")
            #filelist = [f for f in os.listdir(conf.output_dir)
                        #if f.endswith(ext)]
        #for f in filelist:
            #os.remove(os.path.join(conf.output_dir, f))


def MakeVideo(conf, filetitle) -> None:
    """ Function creates a video file from a sequence
        of field states written into image files using
        'ffmpeg' utility installed system-wide.
    """
    print("")
    print("")
    print("")
    filename = MakeFileName(conf, filetitle, ".avi")
    if os.path.isfile(filename): os.remove(filename)
    wildcards = os.path.join(conf.output_dir, filetitle + "*.png")
    framerate = 24
    try:
        subprocess.run(["ffmpeg", "-y", "-f", "image2", "-framerate",
                str(framerate), "-pattern_type", "glob", "-i",
                "'" + wildcards + "'", filename], check=True)
    except Exception as error:
        traceback.print_exc()
        print("WARNING: " + str(error.args))
        print("Failed to make video; please, check if 'ffmpeg' was installed")
        print("We can proceed without video ...")

    # Save space by removing image files, which had been encoded in AVI file.
    for f in glob.glob(wildcards):
        os.remove(f)

    print("")
    print("")
    print("")


def CheckPythonVersion() -> None:
    """ Function checks for the minimum supported Python version.
    """
    if (sys.version_info[0] < 3 or
            (sys.version_info[0] == 3 and sys.version_info[1] < 6)):
        sys.exit("Python of version 3.6+ is expected")


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


def WriteFieldAsImage(filename, field) -> object:
    """ Function writes out specified field as an image.
        Here we transpose the field: field(x,y) -> field(y,x)
        because Python assumes abscissas in the second dimension.
        Also, the image is flipped vertically so that the origin
        goes to the left-bottom corner.
    """
    assert len(field.shape) == 2
    # Image <- field scaled to [0..1] range.
    vmin = np.amin(field)
    vmax = np.amax(field)
    image = np.abs((np.transpose(field) - vmin) /
                   (vmax - vmin + math.sqrt(np.finfo(float).eps)))
    # Write image file.
    image = np.flipud(image)            # flip Y
    if filename is not None:
        assert isinstance(filename, str)
        scipy.misc.imsave(filename, image)
    return image


def ReadResultFile(filename):
    """ Function reads the binary file where a number of solutions (full state
        fields) are stacked one after another. The file format has 4-column
        layout (time, abscissa, ordinate, value), all values have float-32
        precision and records are not necessary sorted (in case of C++ output).
    """
    # Check file name has the pattern of a file of solution fields.
    assert re.search(r".*field_Nx\d+_Ny\d+_Nt\d+\.txt", filename), (
            "wrong pattern of the file for solution field")
    # Extracts Nx, Ny and Nt from the file name.
    digits = re.sub("[^0-9]", " ", os.path.basename(filename))
    params = [int(s) for s in digits.split() if s.isdigit()]
    assert len(params) == 3
    Nx = int(params[0])
    Ny = int(params[1])
    Nt = int(params[2])
    # Read the solution fields file and sort the records because the parallel
    # C++ application does not guarantee proper ordering. Actually, we do not
    # sort the data as such, rather get the sorting index array.
    with open(filename, "rb") as fid:
        data = np.fromfile(fid, dtype=np.float32)
    data = np.reshape(data, (-1,4))
    Nw = data.shape[0] // (Nx * Ny)             # number of written records
    assert data.shape[0] == Nx * Ny * Nw, "wrong file size"
    idx = np.lexsort(np.rot90(data[:,0:3]))     # sort by {t,x,y} triples
    # Check the data and form the output fields and corresponding timestamps.
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
    """ Class implements progress printing in terminal window.
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


#if __name__ == "__main__":
#    ReadResultFile("../output/true-field_Nx121_Ny121_Nt840.txt")
