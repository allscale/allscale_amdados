#------------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017
#------------------------------------------------------------------------------

import numpy as np
import math
import scipy.sparse
import scipy.sparse.linalg
import Configuration
import Utils


class ImplicitEuler:
    """ Process model - implicit Euler time integrator for advection-diffusion equation.
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
        # Get and initialize model parameters.
        conf = self.conf
        nx = int(conf["num_domains_x"] * conf["num_elems_x"])   # here 'elems' are nodal points
        ny = int(conf["num_domains_y"] * conf["num_elems_y"])   # here 'elems' are nodal points
        assert nx > 1 and ny > 1

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

        # Initial field: all zeros except around a point with high concentration.
        field = self.InitialField(conf)
        field = self.ApplyBoundaryCondition(border, field)

        # Compute and save infrastructure for assembling sparse model matrix.
        self.border = border
        self.row_idx, self.col_idx, self.values = self.ModelMatrixInfrastruct()
        return field

# public:
    def Iterate(self, t, dt, flow_model, field, covar):
        """ Function propagates the solution of advection-diffusion equation one time step ahead,
            and updates the state field and covariance matrix according to the formulae:
                field_{t+1} = A_t * field_t
                covar_{t+1} = A_t * covar_t * transpose(A_t) + Q_t
            Also, the time step 'dt' can be adjusted, if it too large, hence the 'next' logical
            timestamp 't+1' corresponds to the physical time 't+dt' after 'dt' adjustment.
        """
        assert field.shape == (self.size_x, self.size_y)

        # Domain and problem sizes, spacial discretization steps.
        nx = int(self.size_x)
        ny = int(self.size_y)
        problem_size = int(self.problem_size)
        dx = self.dx
        dy = self.dy
        D = float(self.conf["diffusion_coef"])

        # By assumption, the flow does not depend on coordinates but on time.
        vx, vy, max_vx, max_vy = flow_model.Flow(t)

        # Deduce the optimal time step from the stability criteria.
        tiny = np.finfo(float).tiny / (np.finfo(float).eps)**2
        dt = min(dt, min( min(dx**2, dy**2)/(2.0*D + tiny),
                          1.0/(abs(max_vx)/dx + abs(max_vy)/dy + tiny) ) )
        assert dt > 0

        # Here we assume the flow velocity is the same everywhere given the time 't'.
        rho_x = D * dt / dx**2;   v0x = 2.0 * dx / dt;   vx = float(vx) / v0x
        rho_y = D * dt / dy**2;   v0y = 2.0 * dy / dt;   vy = float(vy) / v0y

        # Here we repeat the loops in Initialize(..) except that the sparse
        # matrix entries are computed rather than row/column indices.
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
        # Before processing, the state is unrolled into a column vector,
        # afterwards it is reshaped back to 2D field.
        field = Utils.ColumnVec(field)
        B = scipy.sparse.csr_matrix((vals, (self.row_idx, self.col_idx)),
                                    shape=(problem_size, problem_size))

        # XXX
        EnableModelMatrixHack = False
        ApplyModelMatrixHack = False
        if EnableModelMatrixHack:
            neigs = 1
            if ApplyModelMatrixHack:
                B += scipy.sparse.eye(problem_size) * 0.025
            e1 = scipy.sparse.linalg.eigsh(B, k=neigs, which='LM', return_eigenvectors=False)
            e2 = scipy.sparse.linalg.eigsh(B, k=neigs, which='SM', return_eigenvectors=False)
            print('A: largest eigenvalue: ' + str(1.0/e2) + ', smallest eigenvalue: ' + str(1.0/e1))
        # XXX

        field = scipy.sparse.linalg.spsolve(B, field)
        field = np.reshape(field, (nx, ny))
        self.ApplyBoundaryCondition(self.border, field)

        # covar <- A * covar * A^t + Q = B^{-1} * covar * B^{-t} + Q =
        #                                B^{-1} * (B^{-1} * covar^t)^t + Q.
        if covar is not None:
            assert isinstance(covar, np.ndarray)
            assert covar.shape == (self.problem_size, self.problem_size)
            covar = scipy.sparse.linalg.spsolve(B, covar.transpose())
            covar = scipy.sparse.linalg.spsolve(B, covar.transpose())
            assert isinstance(covar, np.ndarray)
            assert covar.shape == (self.problem_size, self.problem_size)
            covar = covar + self.ComputeQ()
        return dt, field, covar

# public:
    def InitialField(self, conf):
        """ Create initial ("true") field with a spike at some point and zeros elsewhere.
            Note, the spike is not very sharp to make the field differentiable.
        """
        nx = self.size_x
        ny = self.size_y
        field = np.zeros((nx, ny), dtype=float)
        sigma = int(1);                             # in logical units (point indices)
        cx = round(float(conf["spot_x"]) / self.dx)
        cy = round(float(conf["spot_y"]) / self.dy)
        a = float(conf["spot_density"]) / (float(sigma**2) * 2.0 * math.pi)
        b = 1.0 / (2.0 * float(sigma)**2)

        for x in range(-4*sigma, 4*sigma+1):
            ix = min(max(x + cx, 0), nx-1)
            for y in range(-4*sigma, 4*sigma+1):
                iy = min(max(y + cy, 0), ny-1)
                field[ix,iy] += a * math.exp(-b * float(x*x + y*y))
        return field

# private:
    def ApplyBoundaryCondition(self, border, field):
        """ Function applies Dirichlet boundary conditions at the boundary.
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

# private:
    def ComputeQ(self):
        """ Function computes the model noise covariance matrix.
        """
        a = float(self.conf["model_noise_Q"])
        v = a * np.random.rand(self.problem_size) + 1.0
        Q = scipy.sparse.diags(v)
        assert Q.shape == (self.problem_size, self.problem_size)
        return Q

