#------------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017
#------------------------------------------------------------------------------

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
import Configuration


class ImplicitEuler:
    """ Implicit Euler time integrator for advection-diffusion equation.
    """
# public:
    def __init__(self, conf):
        """ Constructor
        """
        self.conf = conf            # reference to the global object for concentration settings
        self.size_x = int(0)        # number of nodal points in x-dimension
        self.size_y = int(0)        # number of nodal points in y-dimension
        self.problem_size = int(0)  # total number of nodal points in the domain
        self.dx = float(0)          # distance between neighbor nodal points in x-dimension
        self.dy = float(0)          # distance between neighbor nodal points in y-dimension
        self.border = None          # (i,j) indices of the domain's boundary points
        self.row_idx = None         # placeholder for row indexes to build a sparse model matrix
        self.col_idx = None         # placeholder for column indexes to build a sparse model matrix
        self.values = None          # placeholder for values to build a sparse model matrix

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

        self.size_x = nx
        self.size_y = ny
        self.problem_size = nx * ny
        self.dx = float(conf["domain_size_x"]) / float(nx-1)
        self.dy = float(conf["domain_size_y"]) / float(ny-1)
        assert self.dx > 0 and self.dy > 0

        assert float(self.conf["diffusion_coef"]) > 0, "expects positive diffusion coefficient"

        # N x 2 matrix of indices of border points.
        border = np.array( [(x,y) for x in range(nx)
                                  for y in range(ny) if x == 0 or x == nx-1 or
                                                        y == 0 or y == ny-1] )

        # Initial field: all zeros except one point with high concentration.
        field = np.zeros((nx, ny), dtype=float)
        field[round(conf["spot_x"] / self.dx),
              round(conf["spot_y"] / self.dy)] = conf["spot_density"]
        field = self.ApplyBoundaryCondition(border, field)

        # Compute and save infrastructure for assembling sparse model matrix.
        self.border = border
        self.row_idx, self.col_idx, self.values = self.ModelMatrixInfrastruct()
        return field

# public:
    def Iterate(self, t, dt, flow_model, field, covar):
        assert field.shape == (self.size_x, self.size_y)
        field = np.reshape(field, (self.problem_size,))
        assert covar is None
        #assert covar.shape == (self.problem_size, self.problem_size)

        # Domain and problem sizes, spacial discretization steps.
        nx = int(self.size_x)
        ny = int(self.size_y)
        problem_size = int(self.problem_size)
        dx = self.dx
        dy = self.dy
        assert dx > 0 and dy > 0
        D = float(self.conf["diffusion_coef"])

        # By assumption, flow does not depend on coordinates but on time.
        vx, vy, max_vx, max_vy = flow_model.Flow(t)

        # Deduce the optimal time step from the stability criteria.
        tiny = np.finfo(float).tiny / (np.finfo(float).eps)**2
        dt = min(dt, min( min(dx**2, dy**2)/(2.0*D + tiny),
                          1.0/(abs(max_vx)/dx + abs(max_vy)/dy + tiny) ) )
        assert dt > 0

        rho_x = D * dt / dx**2;   v0x = 2.0 * dx / dt;   vx = float(vx) / v0x
        rho_y = D * dt / dy**2;   v0y = 2.0 * dy / dt;   vy = float(vy) / v0y

        # Here we repeat the loops in Initialize(..) except that the sparse
        # matrix entries are computed rather than row/column indices,
        # which are already done at this stage.
        vals = self.values
        N = int(0)
        for x in range(1,nx-1):
            for y in range(1,ny-1):
                vals[N+0] = float(1.0 + 2.0*(rho_x + rho_y))    # mid   (x,y)
                vals[N+1] = float(- vx - rho_x)                 # west  (x-1,y)
                vals[N+2] = float(+ vx - rho_x)                 # east  (x+1,y)
                vals[N+3] = float(- vy - rho_y)                 # south (x,y-1)
                vals[N+4] = float(+ vy - rho_y)                 # north (x,y+1)
                N += 5
        for i in range(self.border.shape[0]):
            vals[N] = 1.0
            N += 1
        assert N == len(self.row_idx)

        # Assemble (inverse) model matrix and propagate the state field
        # forward in time: state_{t+1} = B^{-1} * state_t.
        # Before processing, the state is unrolled into column vector,
        # afterwards it is reshaped back to 2D field.
        field = np.reshape(field, (problem_size,))
        B = csr_matrix((vals, (self.row_idx, self.col_idx)),
                        shape=(problem_size, problem_size))
        field = spsolve(B, field)
        field = np.reshape(field, (nx, ny))
        return dt, field

# private:
    def ApplyBoundaryCondition(self, border, field):
        """
        """
        field[ border[:,0], border[:,1] ] = 0.0
        return field

# private:
    def ModelMatrixInfrastruct(self):
        """ Function creates an infrastructure for quick assembling
            of a sparse model matrix.
        """
        nx = self.size_x
        ny = self.size_y

        # Function maps a couple of coordinates to a plain 1D index.
        sub2ind = lambda x, y : int(x * ny + y)

        nnz = int(4 + 2*(nx-2) + 2*(ny-2) + 5*(nx-2)*(ny-2))
        rows = np.zeros((nnz,), dtype=int)
        cols = np.zeros((nnz,), dtype=int)
        vals = np.zeros((nnz,), dtype=float)

        # Get row and column indices of all internal domain points.
        N = int(0)
        for x in range(1,nx-1):
            for y in range(1,ny-1):
                mid = sub2ind(x,y)
                west = sub2ind(x-1,y)
                east = sub2ind(x+1,y)
                south = sub2ind(x,y-1)
                north = sub2ind(x,y+1)
                rows[N+0] = mid; cols[N+0] = mid
                rows[N+1] = mid; cols[N+1] = west
                rows[N+2] = mid; cols[N+2] = east
                rows[N+3] = mid; cols[N+3] = south
                rows[N+4] = mid; cols[N+4] = north
                N += 5
        assert N < nnz

        # Get row and column indices of all boundary domain points.
        border = self.border
        for i in range(border.shape[0]):
            idx = sub2ind(border[i,0], border[i,1])
            rows[N] = idx
            cols[N] = idx
            N += 1
        assert N == nnz

        return rows, cols, vals

