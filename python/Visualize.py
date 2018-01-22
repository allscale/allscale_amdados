# -----------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2018
# -----------------------------------------------------------------------------

#import pdb; pdb.set_trace()           # enables debugging
import matplotlib
matplotlib.use("Agg")
import sys, traceback, os, math, argparse
from Utility import *

if __name__ == "__main__":
    try:
        CheckPythonVersion()
        parser = argparse.ArgumentParser()
        parser.add_argument("--field_file",
                type=str, default=None,
                help="path to solution fields file.")
        param = parser.parse_args()
        param.field_file = os.path.expanduser(param.field_file)
        print("Options:")
        print("Solution fields file: " + param.field_file)
        print("")

        # Read the data file of solution fields.
        timestamps, fields = ReadResultFile(param.field_file)

        # Plot the separate fields like video.
        plt = SwitchToGraphicalBackend()
        hFigure = None
        for i in range(fields.shape[0]):
            image = WriteFieldAsImage(None, fields[i,:,:])
            if hFigure is None:
                hFigure = plt.imshow(image)
            else:
                hFigure.set_data(image)
            plt.title("Density, t=" + str(timestamps[i]),
                        fontsize=10, fontweight="bold")
            plt.pause(0.64)
            plt.draw()

    except AssertionError as error:
        traceback.print_exc()
        print("ERROR: " + str(error.args))
    except ValueError as error:
        traceback.print_exc()
        print("ERROR: " + str(error.args))
    except Exception as error:
        traceback.print_exc()
        print("ERROR: " + str(error.args))
