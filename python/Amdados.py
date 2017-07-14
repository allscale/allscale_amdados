#------------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017
#------------------------------------------------------------------------------

"""
Amdados application.
"""
print(__doc__)

import pdb; pdb.set_trace()

import sys;  sys.path.append("./python")
import os
import getopt
import math
import numpy as np
import scipy
import scipy.misc
import Configuration as Config
import ImplicitEuler
import FlowModel
import Utils
import Observations

import random                   # debugging only

def ParseCmdLine(params):
    """ Parse command-line arguments.
    """
    gen_observations = False
    opts, args = getopt.getopt(params, "hg")
    for opt, arg in opts:
        if opt == "-h":
            print("\nHelp:")
            print("-h print this help")
            print("-g generate observation, otherwise run full-scale application")
            print("\n")
            sys.exit()
        elif opt == "-g":
            gen_observations = True
        else:
            sys.exit("Unknown option: " + opt + ", try -h for help")
    return gen_observations


if __name__ == '__main__':
    try:
        # Parse command-line arguments.
        gen_observations = ParseCmdLine(sys.argv[1:])

        # Read configuration file, create output directories,
        # create and initialize the major application objects.
        conf = Config.ReadConfiguration("amdados.conf")
        Utils.CreateAndCleanOutputDir(conf)
        integrator = ImplicitEuler.ImplicitEuler(conf)
        flow_model = FlowModel.FlowModel(conf)
        field = integrator.Initialize()

        # Either generate (write) observations or read them, but NOT both.
        obvs = Observations.Observations(conf)
        #if not gen_observations:
            #obvs.Read()

        T = float(conf["integration_period"])
        t = float(0)
        dt_base = T / float(conf["integration_nsteps"])
        covar = None

        if     gen_observations: obvs.Write(t, field)
        if not gen_observations: Utils.WriteProgress(conf, t, field)

        # Time-integration loop.
        while True:
            print(".", end='', flush=True)

            dt = dt_base
            dt, field = integrator.Iterate(min(t,T), dt, flow_model, field, covar)

            if gen_observations:
                obvs.Write(t, field)
            else:
                Utils.WriteProgress(conf, t, field)

            if t >= T: break
            t += dt
        print("\n\n\n")
        if not gen_observations: Utils.MakeVideo(conf, "field", field.shape)
    except Exception as error:
        print('ERROR: ' + str(error.args))





        #if True:
            #obvs = Observations.Observations(conf)
            #obvs.Read()
            #for k in range(100):
                #t = random.uniform(obvs.timestamps[0], obvs.timestamps[-1])
                #field = obvs.GetObservation(t)
                #print(str(t) + "    " + str(np.mean(field)))
            #sys.exit()
