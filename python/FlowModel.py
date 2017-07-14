#------------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017
#------------------------------------------------------------------------------

import math
import Configuration


class FlowModel:
    """ Class computes flow components on each iteration
        of advection-diffusion problem solver.
        Note, here we assume that flow depends only on time but not on coordinates.
    """
# public:
    def __init__(self, conf):
        """ Constructor.
        """
        self.T = float(conf["integration_period"])
        assert self.T > 0.0
        self.max_vx = float(conf["flow_model_max_vx"])
        self.max_vy = float(conf["flow_model_max_vy"])
        assert self.max_vx >= 0.0 and self.max_vy >= 0.0, "negative max. flow velocities"

# public:
    def Flow(self, t):
        """ Function computes flow components given a time. It also returns
            the modules of maximum flow velocities in each dimension.
        """
        vx = -self.max_vx * math.sin(0.1*float(t)/self.T - math.pi)
        vy = -self.max_vy * math.sin(0.2*float(t)/self.T - math.pi)
        return vx, vy, self.max_vx, self.max_vy

