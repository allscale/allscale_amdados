# -----------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017-2018
# -----------------------------------------------------------------------------

""" This script runs several simulations with increasing problem size
    utilizing all available CPUs. It saves the execution time of each
    simulation in a file that can be used to plot the scalability profile.
    Essential parameters, listed in the first lines, include the set of
    problem sizes and integration period.
      Each simulation is twofold. First, we run the Python forward solver that
    generates the ground-truth and observations ("ObservationsGenerator.py").
    The Python code itself uses the C++ code running in the special mode for
    generating sensor locations (scenario "sensors"). Second, the C++ data
    assimilation application is launched (scenario "simulation") with
    observations generated by "ObservationsGenerator.py".
      The results of all the simulations are accumulated in the output
    directory and can be visualized later on by the script "Visualize.py".
      The configuration file "amdados.conf" is used in all the simulations with
    modification of three parameters: grid sizes (number of subdomains) in both
    dimensions and integration time. Other parameters remain intact. It is not
    recommended to tweak parameters unless their meaning is absolutely clear.
      If you had modified the parameters, please, consider to rerun this script
    because the results in the output directory a not valid any longer.
      The script was designed to fulfil the formal requirements of the
    Allscale project.
"""
print(__doc__)

from timeit import default_timer as timer
import os, cmd
from RandObservationsGenerator import InitDependentParams, Amdados2D
from Utility import *


nthreads = np.arange(0,46,2)
GridSizes = np.zeros([len(nthreads), 2])
nthreads[0] = 1
nthreads[1] = 3
GridSizes[0:13, :] = ([4,2], [6,4], [8,4], [8,6],[8,8], [10,8],[12, 8],[16,7],
                       [16,8], [12,12], [20,8], [16,11], [16,12])
for i in range(13, len(GridSizes)):
        GridSizes[i, :] = [(i)*2, 8]
        problem_size= GridSizes[:,0]*GridSizes[:,1]

execute_time =  np.zeros([len(nthreads), 4])
# Integration period in seconds.
IntegrationPeriod = 25 
IntegrationNsteps = 50
# Path to the C++ executable.
AMDADOS_EXE = "build/mpi_amdados"


if __name__ == "__main__":
    try:
        # Read configuration file.
        conf = Configuration("amdados.conf")
        # Create the output directory, if it does not exist.
        if not os.path.isdir(conf.output_dir):
            os.mkdir(conf.output_dir)
        # Check existence of "amdados" application executable.
        assert os.path.isfile(AMDADOS_EXE), "amdados executable was not found"

        # For all the grid sizes in the list ...
        exe_time_profile = np.zeros((len(GridSizes),2))
        for i in range(0, len(nthreads)):
            grid = GridSizes[i]
            Nproc = nthreads[i]
            # Modify parameters given the current grid size.
            setattr(conf, "num_subdomains_x", int(grid[0]))
            setattr(conf, "num_subdomains_y", int(grid[1]))
            setattr(conf, "integration_period", int(IntegrationPeriod))
            setattr(conf, "integration_nsteps", int(IntegrationNsteps))
            InitDependentParams(conf)
            conf.PrintParameters()
            config_file = conf.WriteParameterFile("scalability_test.conf")
            os.sync()
            # Python simulator generates the ground-truth and observations.
            Amdados2D(config_file, False)

            # Get the starting time.
            start_time = timer()

            # Run C++ data assimilation application.
            print("##################################################")
            print("Simulation by 'amdados' application for series of hpx threads")
            print("Initial simulations for Oceans paper")
            print("##################################################")
            print(AMDADOS_EXE, config_file)
            print("NPROC =", str(Nproc))
            output = subprocess.Popen(["mpirun", "--allow-run-as-root", "-np", str(Nproc),
                "./build/mpi_amdados", "--scenario", "simulation", "--config", config_file], stdout=subprocess.PIPE)
            output.wait()
            assert output.returncode == 0, "amdados returned non-zero status"

            # Strip the execution time from stdout, both total simulation time
            # and throughput (subdomain/s)
            strip_output = str(output.communicate()[0]).split('\\n')
            for line in strip_output:
                if re.search("Simulation took", line):
                    simtime_string = line
                if re.search("Throughput", line):
                    throughput_string = line
            print('simtime read =', simtime_string)
            print('simtime string =', simtime_string.split(' ')[2][:-1])
            print('throughput string =', throughput_string.split(' ')[1])

            simtime_secs = float(simtime_string.split(' ')[2][:-1])
            throughput_secs = float(throughput_string.split(' ')[1])


            # Get the execution time and corresponding (global) problem size
            # and save the current scalability profile into the file.
            problem_size = ( conf.num_subdomains_x * conf.num_subdomains_y  )
            execute_time[i, :] = [problem_size, Nproc, simtime_secs, throughput_secs]
            np.savetxt(os.path.join(conf.output_dir, "scalability_performance_MPI.txt"),
                       execute_time)

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
