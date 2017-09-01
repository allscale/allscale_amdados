#------------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017
#------------------------------------------------------------------------------

import numpy as np
import numpy.linalg
import scipy.sparse
import scipy.linalg
import math
import Configuration
import Utils


class KalmanFilter:
    """ Implementation of discrete time Kalman filter.
    """
# public:
    def __init__(self, conf):
        """ Constructor
        """
        self.conf = conf            # reference to the global object for concentration settings
        self.size_x = int(0)        # number of nodal points in x-dimension
        self.size_y = int(0)        # number of nodal points in y-dimension
        self.problem_size = int(0)  # total number of nodal points in the domain
        self.H = None               # observation matrix

# public:
    def Initialize(self):
        """ Function initializes this object, creates a property field
            according to initial condition, and prepares some internal
            variables for the subsequent time-integration process.
            Return: initial state field.
        """
        conf = self.conf
        nx = int(conf["num_domains_x"] * conf["num_elems_x"])   # here 'elems' are nodal points
        ny = int(conf["num_domains_y"] * conf["num_elems_y"])   # here 'elems' are nodal points
        assert nx > 1 and ny > 1

        self.size_x = nx
        self.size_y = ny
        self.problem_size = nx * ny

        # Computes the matrix of observations.
        if True:
            self.H = self.ComputeLowRankH()
        else:
            self.H = scipy.sparse.identity(self.problem_size)

        # Initial state field and covariance matrix.
        field = np.zeros((nx, ny))
        covar = Utils.Symmetrize(self.GetModelCovar())
        return field, covar

# public:
    def GetObservations(self, field):
        """ Given a field, function returns (possibly smaller) vector of observations.
        """
        assert isinstance(field, np.ndarray)
        H = self.H
        field1D = Utils.ColumnVec(field)
        z = H.dot(field1D)
        return z

# public:
    def Iterate(self, observations, field_prior, covar_prior):
        """ Function propagates the Kalman filter one time step ahead,
        """
        assert field_prior.shape == (self.size_x, self.size_y)
        assert isinstance(observations, np.ndarray)
        assert isinstance(field_prior, np.ndarray)
        assert isinstance(covar_prior, np.ndarray)

        # Domain and problem sizes, spacial discretization steps.
        nx = int(self.size_x)
        ny = int(self.size_y)
        problem_size = int(self.problem_size)
        H = self.H
        observations = Utils.ColumnVec(observations)

        # Unroll the field into 1D vector.
        field1D = Utils.ColumnVec(field_prior)

        # y = z - H*x
        y = H.dot(field1D)
        y = observations - y

        # PHt = P_prior * H^t
        # S = H * P_prior * H^t + R
        PHt = H.dot(covar_prior)        # sparse matrix on the left for speed
        PHt = PHt.transpose()
        S = H.dot(PHt)
        S = S + self.ComputeR()
        S = Utils.Symmetrize(S)

        # invS_y = S^{-1} * y
        invS_y = scipy.linalg.solve(S, y)
        invS_y = Utils.ColumnVec(invS_y)

        # invS_HP = S^{-1} * H * P_prior
        HP = PHt.transpose()
        invS_HP = scipy.linalg.solve(S, HP)

        # x = x_prior + K * y = x_prior + P_prior * H^t * S^{-1} * y
        field = PHt.dot(invS_y)
        field = Utils.ColumnVec(field)
        field = field + field1D

        # P = (I - K * H) * P_prior = P_prior - P_prior * H^t * S^{-1} * H * P_prior
        covar = PHt.dot(invS_HP)
        covar = covar_prior - covar

        field = np.reshape(field, (nx, ny))
        covar = Utils.Symmetrize(covar)
        return field, covar

# private:
    def ComputeLowRankH(self):
        """ Function initializes observation matrix H such that observations
            are available at a (small) subset of domain points.
        """
        step = max(int(round(0.5 * float(self.conf["model_covar_radius"]))), int(1));
        if step > 1:
            print("WARNING: the number of observations is less than the number of nodal points")

        # Domain and problem sizes, spacial discretization steps.
        nx = int(self.size_x)
        ny = int(self.size_y)
        problem_size = int(self.problem_size)

        nnz = int(math.ceil(float(nx + step) / float(step)) *
                  math.ceil(float(ny + step) / float(step)))
        rows = np.zeros((nnz,), dtype=int)
        cols = np.zeros((nnz,), dtype=int)
        vals = np.ones((nnz,), dtype=float)

        # Function maps a couple of coordinates to a plain 1D index.
        sub2ind = lambda x, y : int(x * ny + y)

        # Pick up few domain points where observations are available.
        N = int(0)
        for x in range(1,nx-1,step):
            for y in range(1,ny-1,step):
                rows[N] = N
                cols[N] = sub2ind(x,y)
                N += 1
        assert N <= nnz

        # (N x problem_size) matrix of observations.
        H = scipy.sparse.csr_matrix((vals[0:N], (rows[0:N], cols[0:N])), shape=(N, problem_size))
        return H

# private:
    def GetModelCovar(self):
        """ Function computes the model covariance matrix as function of exponential distance.
        """
        # Domain and problem sizes, spacial discretization steps.
        nx = int(self.size_x)
        ny = int(self.size_y)
        problem_size = int(self.problem_size)

        var = float(self.conf["model_covar_var"]);
        radius = max(float(self.conf["model_covar_radius"]), float(1));
        #Rx  = float(self.conf["model_covar_Rx" ]);
        #Ry  = float(self.conf["model_covar_Ry" ]);

        covar = np.zeros((problem_size, problem_size))

        # Function maps a couple of coordinates to a plain 1D index.
        sub2ind = lambda x, y : int(x * ny + y)

        for x1 in range(nx):
            for y1 in range (ny):
                i = sub2ind(x1,y1)
                covar[i,i] = var
                for x2 in range(x1+1,nx):
                    dx = abs(x1 - x2) / radius  # Rx
                    if dx < 4:
                        for y2 in range(y1+1,ny):
                            dy = abs(y1 - y2) / radius      # Ry
                            if dy < 4:
                                j = sub2ind(x2,y2)
                                cov = var * math.exp(-(dx*dx + dy*dy))
                                covar[i,j] = cov
                                covar[j,i] = cov
        return covar

# private:
    def ComputeR(self):
        """ Function computes the observation noise covariance matrix.
        """
        v = np.random.rand(self.H.shape[0]) + 1.0
        R = scipy.sparse.diags(v)
        return R


