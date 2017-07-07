#------------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017
#------------------------------------------------------------------------------

import math
import Configuration


class FlowModel:
    """ Class computes flow components on each iteration
        of advection-diffusion problem solver.
    """
# public:
    def __init__(self, conf):
        """ Constructor.
        """
        self.max_vx = float(conf["flow_model_max_vx"])
        self.max_vy = float(conf["flow_model_max_vy"])
        assert self.max_vx >= 0.0 and self.max_vy >= 0.0, "negative max. flow velocities"

# public:
    def Flow(self, t):
        """ Function computes flow components given a time. It also returns
            the modules of maximum flow velocities in each dimension.
        """
        vx = -self.max_vx * math.sin(float(t)/10.0 - math.pi)
        vy = -self.max_vy * math.sin(float(t)/5.0  - math.pi)
        return vx, vy, self.max_vx, self.max_vy

