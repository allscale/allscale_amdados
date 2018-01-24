# -----------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2018
# -----------------------------------------------------------------------------

#import pdb; pdb.set_trace()           # enables debugging
import numpy as np
from scipy.stats import mstats
import matplotlib
matplotlib.use("Agg")
from matplotlib.ticker import MultipleLocator
import sys, traceback, os, math, argparse
from Utility import *
from Configuration import Configuration

def GetTrueFieldFilename(filename):
    """ Given the name of the solution file, function constructs the name
        of the "true fields" solution file.
    """
    dirname = os.path.dirname(filename)
    basename = os.path.basename(filename)
    if not basename.startswith("field_"):
        sys.exit("\nPlease, use the file named 'field_***.txt' as an input.\n")
    return os.path.join(dirname, "true_" + basename)


def PlotRelativeDifference(field_filename, rel_diff):
    """ Function plots and save in the file the difference between data
        assimilation simulation and analytic one.
    """
    dirname = os.path.dirname(field_filename)
    Nx, Ny, _ = ProblemParametersFromFilename(field_filename, True, "field")
    filename = "rel_diff_Nx" + str(Nx) + "_Ny" + str(Ny) + ".png"

    matplotlib.rcParams["legend.fontsize"] = 10
    fig = plt.figure()
    ax = fig.gca()
    ax.plot(np.arange(1, len(rel_diff) + 1), rel_diff)
    ax.set_xlabel("index of time slice")
    ax.set_ylabel("relative difference")
    plt.title("|analytic - simulation| / |analytic|")
    plt.savefig(os.path.join(dirname, filename))


def PlotSensors(field_filename, params):
    """ Function plots the sensor locations.
    """
    Sx = params["subdomain_x"]
    Sy = params["subdomain_y"]
    # Load file of sensors corresponding to the file of solution fields.
    dirname = os.path.dirname(field_filename)
    Nx, Ny, _ = ProblemParametersFromFilename(field_filename, True, "field")
    sensors_filename = os.path.join(dirname,
                            "sensors_Nx" + str(Nx) + "_Ny" + str(Ny) + ".txt")
    sensors = np.loadtxt(sensors_filename)
    assert Nx % Sx == 0 and Ny % Sy == 0
    # Do plotting with grid spacing at the subdomain boundaries.
    fig = plt.figure()
    ax = fig.gca()
    ax.scatter(sensors[:,0], sensors[:,1], s=max(min(Nx,Ny)*0.01, 1.0))
    ax.grid(which="minor", color="r", linestyle=":")
    ax.xaxis.set_minor_locator(MultipleLocator(Sx))
    ax.yaxis.set_minor_locator(MultipleLocator(Sy))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    plt.title("Sensor locations")
    plt.tight_layout()
    plt.savefig(os.path.join(dirname,
                        "sensors_Nx" + str(Nx) + "_Ny" + str(Ny) + ".png"))


def PlotSchwarzProfile(field_filename, params):
    """ Function plots a picture demonstrating how smoothness of the solution
        improves at the subdomain boundaries as Schwarz iterations progressing.
        If the file of Schwarz iteration profile does not exist (this is only
        available with AMDADOS_DEBUGGING flag defined in C++ application),
        then nothing is done.
    """
    Nschwarz = params["Nschwarz"]
    # Load file of sensors corresponding to the file of solution fields.
    dirname = os.path.dirname(field_filename)
    Nx, Ny, Nt = ProblemParametersFromFilename(field_filename, True, "field")
    profile_filename = os.path.join(dirname,
                    "schwarz_diff_Nx" + str(Nx) + "_Ny" + str(Ny) + ".txt")
    if not os.path.isfile(profile_filename):
        print("WARNING: file " + profile_filename + " does not exist")
        return
    profile = np.loadtxt(profile_filename)
    assert len(profile.shape) == 1 and len(profile) == Nt * Nschwarz
    profile = np.reshape(profile, (Nt, Nschwarz))
    # Compute quantiles among Nt samples of size Nschwarz.
    quantiles = mstats.mquantiles(profile, axis=0)
    # Plotting the picture of the relative difference behaviour.
    fig = plt.figure()
    ax = fig.gca()
    labels = ["25%", "50%", "75%"]
    for i, q in enumerate(quantiles):
        plt.plot(np.arange(1, Nschwarz+1), q, label=labels[i])
    plt.legend()
    ax.set_xlabel("Schwarz iteration")
    ax.set_ylabel("Relative difference")
    plt.title("Decay of relative difference at subdomain boundaries")
    plt.tight_layout()
    plt.grid()
    plt.savefig(os.path.join(dirname,
                    "schwarz_diff_Nx" + str(Nx) + "_Ny" + str(Ny) + ".png"))


if __name__ == "__main__":
    try:
        CheckPythonVersion()
        parser = argparse.ArgumentParser()
        parser.add_argument("--field_file",
                type=str, default=None,
                help="path to solution fields file.")
        opts = parser.parse_args()
        field_file = os.path.expanduser(opts.field_file)
        true_field_file = GetTrueFieldFilename(opts.field_file)
        print("Parameters:")
        print("Solution fields file: " + field_file)
        print("True solution fields file: " + true_field_file)
        print("")

        # Read the data files of the solution fields.
        print("Loading solution fields ...", flush=True)
        true_timestamps, true_fields = ReadResultFile(true_field_file)
        timestamps, fields = ReadResultFile(field_file)
        assert np.all(true_timestamps == timestamps), "mismatch in timestamps"
        assert np.all(true_fields.shape == fields.shape), "mismatch in layouts"
        assert len(timestamps) == fields.shape[0]
        print("done", flush=True)

        # Relative difference between assimilation and true solutions.
        rel_diff = np.zeros((fields.shape[0],))

        # Read some parameters that were stored by ScalabilityTest.py. These
        # parameters must match the configuration used during the simulations.
        params = np.load(
                os.path.join(os.path.dirname(field_file), "params.npy")).item()

        # Plot the separate fields like video.
        plt = SwitchToGraphicalBackend()
        hFig, axarr = plt.subplots(1,2)
        im0, im1 = None, None
        for i in range(fields.shape[0]):
            # Compute the relative difference between two solutions.
            norm_true = np.linalg.norm(true_fields[i,:,:].ravel())
            norm_diff = np.linalg.norm(true_fields[i,:,:].ravel() -
                                            fields[i,:,:].ravel())
            rel_diff[i] = norm_diff / max(norm_true, np.finfo(float).eps)

            # Plot the solution fields.
            true_image = WriteFieldAsImage(None, true_fields[i,:,:])
            image = WriteFieldAsImage(None, fields[i,:,:])
            t = str(timestamps[i])
            if i == 0:
                im0 = axarr[0].imshow(true_image)
            else:
                im0.set_array(true_image)
            axarr[0].set_title("True density, t=" + t,
                                fontsize=10, fontweight="bold")
            if i == 0:
                im1 = axarr[1].imshow(image)
            else:
                im1.set_array(image)
            axarr[1].set_title("Density, t=" + t,
                                fontsize=10, fontweight="bold")
            plt.draw()
            plt.pause(0.20)

        PlotRelativeDifference(field_file, rel_diff)
        PlotSensors(field_file, params)
        PlotSchwarzProfile(field_file, params)

    except AssertionError as error:
        traceback.print_exc()
        print("ERROR: " + str(error.args))
    except ValueError as error:
        traceback.print_exc()
        print("ERROR: " + str(error.args))
    except Exception as error:
        traceback.print_exc()
        print("ERROR: " + str(error.args))
