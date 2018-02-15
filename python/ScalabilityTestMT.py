# -----------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017-2018
# -----------------------------------------------------------------------------

""" This script runs several simulations with the same problem size but
    increasing the number of CPU cores. It saves the execution time of each
    simulation in a file that can be used to plot the scalability profile.
    One can modify the problem size and integration period, see the few global
    variables at the beginning of the script.
      Each simulation is twofold. First, we run the Python forward solver that
    generates the ground-truth and observations. The Python code itself uses
    the C++ code running in the special mode for generating sensor locations
    (scenario "sensors"). Second, the C++ data assimilation application is run
    using the observations previously generated by the Python code (scenario
    "simulation").
      The results of all the simulations are accumulated in the output
    directory and can be visualized later on by the script "Visualize.py".
      The configuration file "amdados.conf" is used in all the simulations with
    modification of three parameters: grid sizes (number of subdomains) in both
    dimensions and integration time. The other parameters remain intact. It is
    not recommended to tweak parameters unless you understand what is done.
      If you modified the parameters, please, consider to rerun this script
    as the results in the output directory as not valid any longer.
      The script was designed to fulfil the formal requirements of the
    Allscale project.

"""

#import pdb; pdb.set_trace()           # enables debugging
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys, traceback, os, glob, getopt, math, argparse, subprocess, psutil
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
GridSize = [8,8]

# Integration period in seconds.
IntegrationPeriod = 1000

# Path to the C++ executable.
AMDADOS_EXE = "build/app/amdados"

if __name__ == "__main__":
    try:
        # Read configuration file.
        conf = Configuration("amdados.conf")
        # Create the output directory, if it does not exist.
        if not os.path.isdir(conf.output_dir): os.mkdir(conf.output_dir)
        # Check existence of "amdados" application executable.
        assert os.path.isfile(AMDADOS_EXE), "amdados executable was not found"
        # Save some parameters that will be used for visualization.
        params = {"subdomain_x" : int(round(conf.subdomain_x)),
                  "subdomain_y" : int(round(conf.subdomain_y)),
                  "Nschwarz"    : int(round(conf.schwarz_num_iters))}
        np.save(os.path.join(conf.output_dir, "params.npy"), params)

        # Modify parameters given the current grid size.
        setattr(conf, "num_subdomains_x", int(GridSize[0]))
        setattr(conf, "num_subdomains_y", int(GridSize[1]))
        setattr(conf, "integration_period", int(IntegrationPeriod))
        InitDependentParams(conf)
        conf.PrintParameters()
        config_file = conf.WriteParameterFile("scalability_test.conf")
        subprocess.run("sync", check=True)

        # Python simulator generates the ground-truth and observations.
        Amdados2D(config_file, False)

        # Get the number of "pure" CPU cores (excluding hyper-threading).
        num_cpu = psutil.cpu_count(logical=False)
        assert num_cpu != None and num_cpu > 1
        exe_time_profile = np.zeros((num_cpu,2))

        # For different number of CPU cores but the same problem size ...
        for n in range(1,num_cpu+1):
            # Get the starting time.
            start_time = timer()

            # Run C++ data assimilation application with different number
            # of working threads.
            print("##################################################")
            print("Simulation by 'amdados' application, which will be")
            print("silent if debugging & messaging were disabled ... ")
            print("##################################################")
            proc = subprocess.Popen([AMDADOS_EXE,
                                    "--scenario", "simulation",
                                    "--config", config_file],
                                    env=dict(os.environ, NUM_WORKERS=str(n)))
            proc.wait()
            assert proc.returncode == 0, "amdados returned non-zero status"

            # Get the execution time and corresponding (global) problem size
            # and save the current scalability profile into the file.
            exe_time_profile[n-1,0] = n
            exe_time_profile[n-1,1] = timer() - start_time
            np.savetxt(os.path.join(conf.output_dir, "scalability_mt.txt"),
                       exe_time_profile)

        # Plot and save the scalability profile.
        plt.plot(exe_time_profile[:,0], exe_time_profile[:,1])
        plt.xlabel("number of CPU cores")
        plt.ylabel("time in seconds")
        plt.grid()
        plt.title("Multi-threading scalability")
        plt.savefig(os.path.join(conf.output_dir, "scalability_mt.png"))

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

