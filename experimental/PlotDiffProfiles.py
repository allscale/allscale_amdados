# -----------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017
# -----------------------------------------------------------------------------

import pdb; pdb.set_trace()           # enables debugging
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys, traceback, os, glob, getopt, math, argparse, subprocess
import numpy as np
import numpy.linalg
import scipy
import scipy.misc
import scipy.sparse
import scipy.sparse.linalg
from timeit import default_timer as timer
from Configuration import Configuration
from ObservationsGenerator import InitDependentParams, Amdados2D
from Utility import *

# Grid sizes (number of subdomains in each direction) for scalability tests.
# This is a set of too big problems, especially for Python (~ few months):
#Resolutions = [(11,11), (23,25), (43,41), (83,89), (167,163), (367,373)]
# This is a reasonable set of problems (~ 10 days / 80 CPU cores).
Resolutions = [(11,11), (19,17), (23,25), (37,31), (43,41), (83,89)]

# Integration period in seconds.
IntegrationPeriod = 9000

TINY = np.finfo(float).tiny / np.finfo(float).eps**3

def SimplePlot(data, x_label, y_label, plot_title, filename):
    """ Function plots a simple chart y = f(x).
    """
    assert isinstance(data, np.ndarray)
    assert isinstance(x_label, str) and isinstance(y_label, str)
    assert isinstance(plot_title, str) and isinstance(filename, str)
    assert len(data.shape) == 2 and data.shape[1] == 2
    # Note that using plt.subplots below is equivalent to using
    # fig = plt.figure and then ax = fig.add_subplot(111)
    fig, ax = plt.subplots()
    ax.plot(data[:,0], data[:,1])
    ax.set(xlabel=x_label, ylabel=y_label, title=plot_title)
    ax.grid()
    fig.savefig(filename)


def CreateDifferenceProfile(file_list1, file_list2):
    """ Function builds the profile of relative error:
        data assimilation solution vs. ground-truth one.
    """
    assert len(file_list1) == len(file_list2), (
            "it must be the same number of fields outputted by the\n"
            "Python and C++ simulators; please, check your settings")

    # Sorting is important to match files with the same time-stamps
    # and ordered from the oldest to the latest.
    file_list1.sort()
    file_list2.sort()

    Nf = len(file_list1)
    diff_profile = np.zeros((Nf,2))
    for k in range(Nf):
        # Load fields generated at the same time by both simulators.
        f1 = np.loadtxt(file_list1[k])
        f2 = np.loadtxt(file_list2[k])

        # Check layouts.
        assert f1.shape[1] == 3, "3-column field file is expected"
        assert f1.shape == f2.shape, (
            "The Python and C++ simulators outputted fields of\n"
            "different sizes; please, check your settings")

        # Lexicographic sorting guarantees the same grid indexing.
        i1 = np.lexsort(np.rot90(f1[:,0:2]))
        i2 = np.lexsort(np.rot90(f2[:,0:2]))
        assert np.all(np.rint(f1[i1,0:2]) == np.rint(f2[i2,0:2])), (
                "mismatch in grid indexing")

        # Get the columns of values and compute the relative difference.
        f1 = f1[i1,2]
        f2 = f2[i2,2]
        diff_profile[k,0] = float(k) / float(max(Nf - 1, 1))
        diff_profile[k,1] = (np.linalg.norm(f1 - f2) /
                max(max(np.linalg.norm(f1), np.linalg.norm(f2)), TINY))

    # Save the plot of relative error at the current resolution.
    SimplePlot(diff_profile,
               "relative time t/T, t=0..T, T is integration period",
               "relative difference",
               "Data-assimilation solution vs. Ground-truth",
               MakeBaseFileName(conf, "field") + "_rel.error.png")


###############################################################################
# Entry point.
###############################################################################
if __name__ == "__main__":
    try:
        # Read configuration file.
        conf = Configuration("amdados.conf")
        conf.output_dir = "./simulation_results"

        for res_no, res in enumerate(Resolutions):
            # Modify parameters given the current resolution.
            setattr(conf, "num_subdomains_x", int(res[0]))
            setattr(conf, "num_subdomains_y", int(res[1]))
            setattr(conf, "integration_period", int(IntegrationPeriod))
            InitDependentParams(conf)
            conf.PrintParameters()

            # Make the profile of relative differences from the file lists
            # of output fields produced by both simulators.
            true_fields = glob.glob(MakeBaseFileName(conf, "true-field")
                                        + "*.txt")
            fields      = glob.glob(MakeBaseFileName(conf, "field")
                                        + "*.txt")
            CreateDifferenceProfile(true_fields, fields)
    except subprocess.CalledProcessError as error:
        traceback.print_exc()
        if error.output is not None:
            print("ERROR: " + str(error.output))
        else:
            print("CalledProcessError")
    except AssertionError as error:
        traceback.print_exc()
        print("ERROR: " + str(error.args))
    except ValueError as error:
        traceback.print_exc()
        print("ERROR: " + str(error.args))
    except Exception as error:
        traceback.print_exc()
        print("ERROR: " + str(error.args))

