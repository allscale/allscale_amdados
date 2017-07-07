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
import math
import numpy as np
import scipy
import scipy.misc
import Configuration as Config
import ImplicitEuler
import FlowModel
import Utils


if __name__ == '__main__':
    try:
        conf = Config.ReadConfiguration("amdados.conf")
        Utils.CreateOutputDir(conf)
        integrator = ImplicitEuler.ImplicitEuler(conf)
        flow_model = FlowModel.FlowModel(conf)
        field = integrator.Initialize()
        T = float(conf["integration_period"])
        t = float(0)
        dt_base = T / float(conf["integration_nsteps"])
        covar = None
        Utils.WriteProgress(conf, field, t)
        while (t < T):
            dt = dt_base
            dt, field = integrator.Iterate(t, dt, flow_model, field, covar)
            Utils.WriteProgress(conf, field, t)
            t += dt
        else:               # finally compute the field exactly at t=T
            dt = dt_base
            dt, field = integrator.Iterate(T, dt, flow_model, field, covar)
            Utils.WriteProgress(conf, field, T, True)
    except Exception as error:
        print('ERROR: ' + str(error.args))

