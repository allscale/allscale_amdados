# -----------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017
# -----------------------------------------------------------------------------

import glob, os, sys, errno, math, subprocess, traceback
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
        assert entity != "true-field", (
                "file suffix is always expected for the true field")
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
    assert isinstance(filename, str) and len(field.shape) == 2
    # Image <- field scaled to [0..1] range.
    vmin = np.amin(field)
    vmax = np.amax(field)
    image = np.abs((np.transpose(field) - vmin) /
                   (vmax - vmin + math.sqrt(np.finfo(float).eps)))
    # Write image file.
    image = np.flipud(image)            # flip Y
    scipy.misc.imsave(filename, image)
    return image


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

