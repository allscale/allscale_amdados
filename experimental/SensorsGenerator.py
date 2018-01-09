# -----------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017
# -----------------------------------------------------------------------------

""" Script generates desired number of sensor locations pseudo-randomly
    seeded across the whole domain. Pseudo-random term means we try to avoid
    overlapping or too close locations, i.e. sensors are more or less
    separated in space.
    Here we assume that sensors covers a fraction of nodal points. This
    fraction should not be too small, otherwise data assimilation method would
    not work. On the other hand, it should not be too large, otherwise
    pseuso-random generator can stuck. Ideally,
    fraction < 0.1 * (number of nodal points).
"""

#import pdb; pdb.set_trace()           # enables debugging
import matplotlib
matplotlib.use("Agg")
import traceback, os, sys, random, math, time, argparse
import numpy as np
from timeit import default_timer as timer
from sklearn.neighbors import NearestNeighbors
from Configuration import Configuration
from Utility import *


def RandomSensorLocations(Nx, Ny, Np):
    """ Generates a number of sensor locations pseudo-randomly
        seeded across the whole domain. This is not the best
        possible algorithm. Consider Poisson disk method, etc.
    """
    assert (Nx >= 20) and (Ny >= 20), "domain size(s) is too small"
    assert    Np >= 10,   "number of sensors is too small"
    assert 10*Np < Nx*Ny, "number of sensors is too large"

    mean_dist = max(round(math.sqrt(float(Nx*Ny)/float(Np))), 1)
    print("Mean distance: " + str(mean_dist))
    ini_dist_thr = round(10 * mean_dist)
    print("Initial distance threshold: " + str(ini_dist_thr))

    points = np.zeros((Np,2), dtype=int)    # sensor location points
    new_pt = np.zeros((1,2), dtype=int)     # new point placeholder
    n = 0                                   # number of points generated so far
    random.seed()

    # Generate the first point.
    points[0,0] = random.randint(2, Nx-3)
    points[0,1] = random.randint(2, Ny-3)
    n += 1
    # Create an instance for nearest-neighbour search object.
    neigh = NearestNeighbors(n_neighbors=1, radius=ini_dist_thr)
    neigh.fit(points[0:n, :])

    # Add random points one by one.
    progress = PrintProgress(Np)
    while n < Np:
        attempt = 0
        # TODO: decide how to set "max_attempts".
        max_attempts = max(10, round(math.ceil(math.sqrt(n))))
        max_attempts = 10   # much faster
        total_attempts = 0
        dist_thr = ini_dist_thr
        # Make a number of attempts until the nearest neighbour
        # is not closer than a distance threshold. Reduce the
        # distance threshold as the number of failed attempts
        # increases beyond the limit.
        while True:
            new_pt[0,0] = random.randint(2, Nx-3)
            new_pt[0,1] = random.randint(2, Ny-3)
            dist, idx = neigh.kneighbors(new_pt)
            assert len(dist) == 1 and len(idx) == 1
            if dist[0,0] >= dist_thr:
                points[n,:] = new_pt
                break
            attempt += 1
            if attempt >= max_attempts:
                dist_thr = max(dist_thr // 2, 1)
                attempt = 0
            total_attempts += 1
            if total_attempts >= 10000:
                sys.exit("ERROR: something went really wrong." + os.linesep +
                         "Random sensor locations cannot be generated.")
        n += 1
        neigh.fit(points[0:n, :])
        progress.Print(n-1)

    progress.Finalize()
    assert n == Np
    return points


def WriteOutputFile(filename, points, Nx, Ny):
    """ Function writes the sensor locations to a file.
    """
    with open(filename, "wt") as f:
        f.write("# Line1: domain sizes Nx, Ny; " +
                "Line2: number of points and their dimension" + os.linesep)
        f.write(str(Nx) + " " + str(Ny) + os.linesep)
        f.write(str(points.shape[0]) + " " +
                str(points.shape[1]) + os.linesep)
        for k in range(points.shape[0]):
            f.write("{} {}".format(points[k, 0], points[k, 1]) + os.linesep)
        f.write(os.linesep)
        f.flush()               # make sure the file is really written
    time.sleep(3)               # wait until OS had written the file on disk


def GenerateRandomSensorLocations(conf):
    """ Function generates new random sensor locations and saves
        them into the text file.
    """
    print("###############################")
    print("Generating sensor locations ...")
    print("###############################")
    Nx = round(conf.num_subdomains_x * conf.subdomain_x)
    Ny = round(conf.num_subdomains_y * conf.subdomain_y)
    Np = round(math.ceil(Nx * Ny * conf.sensor_fraction))
    print("Generator of random sensor locations:")
    print("Domain size: " + str(Nx) + " x " + str(Ny))
    print("Number of sensors: " + str(Np))
    print("")

    # Create path for the output, if it does not exist,
    # and remove the old output file, if it exists.
    output_filename = MakeFileName(conf, "sensors")
    AssurePathExists(output_filename)
    RemoveFile(output_filename)

    # Generate random sensor locations.
    start_time = timer()
    points = RandomSensorLocations(Nx, Ny, Np)
    print("")
    print("execution time: " + str(timer() - start_time) + " seconds")
    print("")
    print("")
    WriteOutputFile(output_filename, points, Nx, Ny)
    return points, Nx, Ny


def PlotSensorLocations(plt, points, Nx, Ny):
    """ Function plots sensors locations.
    """
    plt.figure(1)
    plt.plot(points[:, 0], points[:, 1], ".")
    plt.axis([0, Nx, 0, Ny])
    plt.title("Sensor locations")
    plt.tight_layout()
    plt.show()


# Entry point.
if __name__ == "__main__":
    try:
        CheckPythonVersion()

        # Parse command line parameters.
        parser = argparse.ArgumentParser()
        parser.add_argument("--demo", nargs="?", const=True,
                type=bool, default=False,
                help="flag enables visualization of sensor locations.")
        parser.add_argument("--config_file",
                type=str, default="amdados.conf",
                help="path to configuration file.")
        param = parser.parse_args()
        param.config_file = os.path.expanduser(param.config_file)

        conf = Configuration(param.config_file)
        conf.PrintParameters()
        points, Nx, Ny = GenerateRandomSensorLocations(conf)

        if param.demo:
            plt = SwitchToGraphicalBackend()
            if plt is not None:
                PlotSensorLocations(plt, points, Nx, Ny)

    except AssertionError as error:
        traceback.print_exc()
        print("ERROR: " + str(error.args))
    except ValueError as error:
        traceback.print_exc()
        print("ERROR: " + str(error.args))
    except Exception as error:
        traceback.print_exc()
        print("ERROR: " + str(error.args))


