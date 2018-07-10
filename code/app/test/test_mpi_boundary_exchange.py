# -----------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017-2018
# -----------------------------------------------------------------------------

""" This script tests boundary exchange mechanism in MPI version of
    Amdados application. Run it from the project root folder:

        python3 code/app/test/test_mpi_boundary_exchange.py
"""
print(__doc__)

#import pdb; pdb.set_trace()           # enables debugging
import sys, traceback, os, subprocess
from timeit import default_timer as timer
sys.path.append("python")
from Configuration import Configuration
from ObservationsGenerator import Amdados2D
from Utility import CheckPythonVersion

# Path to the C++ executables.
AMDADOS_EXE = "build/app/amdados"
MPI_AMDADOS_EXE = "build/mpi_amdados"

# Number of parallel processes (on IBM blade).
NUMPROC = 44

# Temporary (modified) configuration file.
CONFIG = "tmp.conf"

if __name__ == "__main__":
    try:
        CheckPythonVersion()

        # Build and check existence of "amdados" application executables.
        subprocess.run(["bash", "standard.build.sh"], check=True)
        subprocess.run(["bash", "mympi"], check=True)
        subprocess.run("sync", check=True)
        assert os.path.isfile(AMDADOS_EXE), \
                "Allscale amdados executable was not found"
        assert os.path.isfile(MPI_AMDADOS_EXE), \
                "MPI amdados executable was not found"

        # Read configuration file.
        conf = Configuration("amdados.conf")

        # Create the output directory, if it does not exist.
        if not os.path.isdir(conf.output_dir):
            os.mkdir(conf.output_dir)

        # For all the grid and subdomain sizes ...
        for Gx in range(2,10,3):
            for Gy in range(2,10,3):
                for Sx in range(3,17):
                    for Sy in range(3,17):
                        # Modify parameters given the current grid size.
                        setattr(conf, "num_subdomains_x", int(Gx))
                        setattr(conf, "num_subdomains_y", int(Gy))
                        setattr(conf, "subdomain_x", int(Sx))
                        setattr(conf, "subdomain_y", int(Sy))
                        setattr(conf, "integration_period", int(100))
                        config_file = conf.WriteParameterFile(CONFIG)
                        subprocess.run("sync", check=True)

                        # Python simulator generates the ground-truth
                        # and observations.
                        Amdados2D(config_file, False)

                        # Run C++ MPI application in testing mode.
                        print("###########################################")
                        print("Testing with the grid and subdomain sizes: ")
                        print("(Gx,Gy) = (" + str(Gx) + "," + str(Gy) + ")")
                        print("(Sx,Sy) = (" + str(Sx) + "," + str(Sy) + ")")
                        print("###########################################")
                        subprocess.run(["mpirun", "-np", str(NUMPROC),
                                        "-f", "host_file", MPI_AMDADOS_EXE,
                                        "--scenario", "simulation",
                                        "--config", config_file,
                                        "--test", "boundary_exchange"],
                                        check=True)

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

