#------------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017
#------------------------------------------------------------------------------

#import pdb; pdb.set_trace()           # enables debugging

import sys                              # ; sys.path.append("./python")
import matplotlib                       # ; matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import traceback
import os
import re
import getopt
import math
import numpy as np
import scipy
import scipy.misc
import scipy.sparse
import scipy.sparse.linalg
from timeit import default_timer as timer

ZERO_BOUNDARY_CONDITIONS = True

def Amdados2D(config_file):
    """ Advection-diffusion PDE solver with data assimilation based on Kalman filter.
    """
    # Initialize parameters, global indices and output directory.
    class Configuration: pass
    conf = Configuration()
    conf = ReadConfiguration(config_file, conf)
    conf = InitDependentParams(conf)
    glo_idx = GlobalIndices(conf)
    for attr, val in vars(conf).items(): print('{} : {}'.format(attr, val))
    CreateAndCleanOutputDir(conf)

    # Run forward simulation and record the "true" solutions.
    true_fields = ComputeTrueFields(conf, glo_idx)
    if conf.generate_observations_only:
        print('\n\n')
        if conf.generate_video: MakeVideo(conf, "true_field")
        return

    # Simulation with data assimilation by the Kalman filter.
    DataAssimilation(conf, glo_idx, true_fields)


###################################################################################################
# Initialization.
###################################################################################################

def ReadConfiguration(filename, conf):
    """ Function reads and parses configuration file of application parameters.
    """
    assert isinstance(filename, str), "expecting file name string"
    assert filename, "empty file name"
    with open(filename) as config_file:
        for line in config_file:
            line_nocomments = re.sub(re.compile("#.*?$"), "", line.strip())
            if line_nocomments:
                name, var = line_nocomments.partition(" ")[::2]
                name = name.strip()
                var = var.strip()
                try:
                    setattr(conf, name, float(var))
                except ValueError:
                    setattr(conf, name, var)        # if not a float, then a string
    return conf


def InitDependentParams(conf):
    """ Function initializes dependent parameters given the primary ones specified by user.
    """
    # Ensure integer values.
    conf.nx = round(conf.num_domains_x * conf.num_elems_x)
    conf.ny = round(conf.num_domains_y * conf.num_elems_y)
    conf.observation_nx = round(conf.observation_nx)
    conf.observation_ny = round(conf.observation_ny)
    conf.integration_nsteps = round(conf.integration_nsteps)

    # Diffusion coefficient must be positive float value.
    D = float(conf.diffusion_coef)
    assert D > 0

    # Deduce space discretization steps.
    conf.problem_size = int(conf.nx * conf.ny)
    dx = float(conf.domain_size_x) / float(conf.nx-1)
    dy = float(conf.domain_size_y) / float(conf.ny-1)
    assert (dx > 0) and (dy > 0)
    conf.dx = dx
    conf.dy = dy

    # Deduce the optimal time step from the stability criteria.
    tiny = np.finfo(float).tiny / (np.finfo(float).eps)**3
    dt_base = float(conf.integration_period) / float(conf.integration_nsteps)
    max_vx = float(conf.flow_model_max_vx)
    max_vy = float(conf.flow_model_max_vy)
    dt = min(dt_base, min( min(dx**2, dy**2)/(2.0*D + tiny),
                           1.0/(abs(max_vx)/dx + abs(max_vy)/dy + tiny) ))
    assert(dt > 0)
    conf.dt = dt
    conf.Nt = round(math.ceil(float(conf.integration_period) / dt))

    # Compute coefficients that will be used in the finite-difference scheme.
    conf.rho_x = float(D * dt / dx**2)
    conf.rho_y = float(D * dt / dy**2)

    conf.v0x = float(2.0 * dx / dt)
    conf.v0y = float(2.0 * dy / dt)
    return conf


def GlobalIndices(conf):
    """ Each nodal point gets a unique global index on the grid.
        Function initializes a 2D array of indices:  index of (x,y) = glo_idx(x,y).
    """
    glo_idx = np.arange(int(conf.problem_size)).reshape((conf.nx, conf.ny))
    return glo_idx


###################################################################################################
# Advection-diffusion PDE stuff.
###################################################################################################

def ApplyBoundaryCondition(conf, field):
    """ Function applies zero Dirichlet boundary condition to the state field.
    """
    assert field.shape[0] == conf.nx and field.shape[1] == conf.ny
    if ZERO_BOUNDARY_CONDITIONS:
        field[0,:] = 0
        field[:,0] = 0
        field[-1,:] = 0
        field[:,-1] = 0
    return field


def InverseModelMatrix(conf, glo_idx, t):
    """ Function initializes inverse matrix of implicit Euler time-integrator:
        B * x_{t+1} = x_{t}, where B = A^{-1} is the matrix returned by this function.
        The matrix must be inverted while iterating forward in time: x_{t+1} = A * x_{t}.
    """
    nx = conf.nx
    ny = conf.ny
    problem_size = conf.problem_size

    if ZERO_BOUNDARY_CONDITIONS:
        #num_nz = 4 + 2*(nx-2) + 2*(ny-2) + 5*(nx-2)*(ny-2)
        num_nz = 5*nx*ny
    else:
        num_nz = 5*nx*ny

    rows = np.zeros((num_nz,), dtype=int)
    cols = np.zeros((num_nz,), dtype=int)
    vals = np.zeros((num_nz,), dtype=float)

    vx, vy = Flow(conf, t)
    rho_x = conf.rho_x
    rho_y = conf.rho_y
    vx = vx / conf.v0x
    vy = vy / conf.v0y

    N = 0
    for x in range(nx):
        for y in range(ny):
            i = glo_idx[x,y]
            if x == 0 or x == nx-1 or y == 0 or y == ny-1:
                if ZERO_BOUNDARY_CONDITIONS:
                    # This is just a Dirichlet boundary condition; insufficient in diffusion case.
                    #rows[N] = i;  cols[N] = i;  vals[N] = 1;  N += 1

                    # Enforce Newmann boundary condition (du/dn = 0).
                    xm = x-1 if x > 0    else 1
                    xp = x+1 if x < nx-1 else nx-2

                    ym = y-1 if y > 0    else 1
                    yp = y+1 if y < ny-1 else ny-2

                    rows[N] = i;  cols[N] = i;     vals[N] = 1 + 2*(rho_x + rho_y);  N += 1
                    rows[N] = i;  cols[N] = glo_idx[xm,y];  vals[N] = - vx - rho_x;  N += 1
                    rows[N] = i;  cols[N] = glo_idx[xp,y];  vals[N] = + vx - rho_x;  N += 1
                    rows[N] = i;  cols[N] = glo_idx[x,ym];  vals[N] = - vy - rho_y;  N += 1
                    rows[N] = i;  cols[N] = glo_idx[x,yp];  vals[N] = + vy - rho_y;  N += 1
                else:
                    rows[N] = i;  cols[N] = i;  vals[N] = 1 + 2*(rho_x + rho_y);  N += 1

                    if x == 0:
                        rows[N] = i;  cols[N] = glo_idx[x,  y];  vals[N] = - 2*vx - rho_x;  N += 1
                        rows[N] = i;  cols[N] = glo_idx[x+1,y];  vals[N] = + 2*vx - rho_x;  N += 1
                    elif x == nx-1:
                        rows[N] = i;  cols[N] = glo_idx[x-1,y];  vals[N] = - 2*vx - rho_x;  N += 1
                        rows[N] = i;  cols[N] = glo_idx[x,  y];  vals[N] = + 2*vx - rho_x;  N += 1
                    else:
                        rows[N] = i;  cols[N] = glo_idx[x-1,y];  vals[N] = - vx - rho_x;    N += 1
                        rows[N] = i;  cols[N] = glo_idx[x+1,y];  vals[N] = + vx - rho_x;    N += 1

                    if y == 0:
                        rows[N] = i;  cols[N] = glo_idx[x,y  ];  vals[N] = - 2*vy - rho_y;  N += 1
                        rows[N] = i;  cols[N] = glo_idx[x,y+1];  vals[N] = + 2*vy - rho_y;  N += 1
                    elif y == ny-1:
                        rows[N] = i;  cols[N] = glo_idx[x,y-1];  vals[N] = - 2*vy - rho_y;  N += 1
                        rows[N] = i;  cols[N] = glo_idx[x,y  ];  vals[N] = + 2*vy - rho_y;  N += 1
                    else:
                        rows[N] = i;  cols[N] = glo_idx[x,y-1];  vals[N] = - vy - rho_y;    N += 1
                        rows[N] = i;  cols[N] = glo_idx[x,y+1];  vals[N] = + vy - rho_y;    N += 1
            else:
                rows[N] = i;  cols[N] = i;               vals[N] = 1 + 2*(rho_x + rho_y);  N += 1
                rows[N] = i;  cols[N] = glo_idx[x-1,y];  vals[N] = - vx - rho_x;           N += 1
                rows[N] = i;  cols[N] = glo_idx[x+1,y];  vals[N] = + vx - rho_x;           N += 1
                rows[N] = i;  cols[N] = glo_idx[x,y-1];  vals[N] = - vy - rho_y;           N += 1
                rows[N] = i;  cols[N] = glo_idx[x,y+1];  vals[N] = + vy - rho_y;           N += 1

    assert N == num_nz
    invA = scipy.sparse.csr_matrix((vals, (rows, cols)), shape=(problem_size, problem_size))
    assert invA.nnz <= num_nz
    return invA


def Flow(conf, t):
    """ Function computes flow components given a time.
    """
    vx = -conf.flow_model_max_vx * math.sin(0.1 * t / conf.integration_period - math.pi)
    vy = -conf.flow_model_max_vy * math.sin(0.2 * t / conf.integration_period - math.pi)
    return vx, vy


def InitialField(conf):
    """ Create initial ("true") field with a spike at some point and zeros elsewhere.
        Note, the spike is not very sharp to make the field differentiable.
    """
    nx = conf.nx
    ny = conf.ny
    field = np.zeros((nx, ny))
    sigma = int(1)                              # in logical units (point indices)
    cx = round(float(conf.spot_x) / conf.dx)
    cy = round(float(conf.spot_y) / conf.dy)
    a = float(conf.spot_density) / (float(sigma)**2 * 2.0 * math.pi)
    b = 1.0 / (2.0 * float(sigma)**2)

    for x in range(-4*sigma, 4*sigma+1):
        ix = min(max(x + cx, 0), nx-1)
        for y in range(-4*sigma, 4*sigma+1):
            iy = min(max(y + cy, 0), ny-1)
            field[ix,iy] += a * math.exp(-b * float(x*x + y*y))

    field = ApplyBoundaryCondition(conf, field)
    return field


def ComputeTrueFields(conf, glo_idx):
    """ Using model matrix A, the function integrates advection-diffusion equation forward in time
        (x_{t+1} = A * x_{t}) and records all the solutions - state fields. These fields are
        considered as the "true" state of the nature and the source of the "true" observations.
    """
    start_time = timer()
    print('\nGenerating observations ...')
    true_fields = np.zeros((conf.nx, conf.ny, conf.Nt))
    field = InitialField(conf)
    writer = None
    if conf.generate_text_output: writer = Writer(conf, conf.analytic_solution)

    hFigure = None
    for k in range(conf.Nt):
        true_fields[:,:,k] = field                          # save the current state field

        # Print the progress.
        print('.', end='', flush=True)
        if (k+1) % 80 == 0: print()
        # Write the field in textual form.
        if conf.generate_text_output: writer.WriteField(field, k, k * conf.dt)
        # Plot the current image.
        vmin = np.amin(field)
        vmax = np.amax(field)
        image = np.abs((np.transpose(field) - vmin) / (vmax - vmin + np.finfo(float).eps))
        if hFigure is None:
            hFigure = plt.imshow(image)
        else:
            hFigure.set_data(image)
        plt.ylim(0, image.shape[0])     # flip Y
        plt.title('True density, t=' + str(k), fontsize=10, fontweight='bold')
        plt.pause(0.05)
        plt.draw()
        # Write image file.
        if conf.generate_video:
            image = np.flipud(image)    # flip Y
            scipy.misc.imsave(conf.output_dir + '/true_field' + format(k,'05d') + '.png', image)

        t = k * conf.dt                                     # physical time
        invA = InverseModelMatrix(conf, glo_idx, t)         # note, we compute A^{-1} not A
        fld1D = np.reshape(field, (conf.problem_size, 1))   # flatten the state field 'x'
        field = scipy.sparse.linalg.spsolve(invA, fld1D)    # x_{t+1} = invA^{-1} * x_t
        field = np.reshape(field, (conf.nx, conf.ny))       # back to 2D field
        field = ApplyBoundaryCondition(conf, field)         # honor the boundary conditions

    if writer is not None: del writer
    plt.close()
    print('\nexecution time: ' + str(timer() - start_time) + ' seconds\n\n')
    return true_fields


###################################################################################################
# Kalman filter stuff.
###################################################################################################

def ComputeQ(conf):
    """ Function computes the model noise covariance matrix.
    """
    v = conf.model_noise_Q * np.random.rand(conf.problem_size) + 1.0
    Q = scipy.sparse.diags(v)
    return Q


def ComputeR(conf):
    """ Function computes the measurement noise covariance matrix.
    """
    v = conf.model_noise_R * np.random.rand(conf.observation_nx * conf.observation_ny) + 1.0
    R = scipy.sparse.diags(v)
    return R


def ComputeH(conf, glo_idx):
    """ Function initializes observation matrix H such that observations
        can be available at a (small) subset of domain points.
    """
    # Domain and problem sizes, spacial discretization steps.
    nx = conf.nx
    ny = conf.ny

    # At least two observations in either dimension, but not more than number of nodal points.
    obs_nx = min(max(round(conf.observation_nx), 2), nx)
    obs_ny = min(max(round(conf.observation_ny), 2), ny)

    num_obs = obs_nx * obs_ny
    if num_obs < conf.problem_size:
        print('WARNING: the number of observations is less than the number of nodal points')
        rows = np.zeros((num_obs,))
        cols = np.zeros((num_obs,))
        vals = np.ones ((num_obs,))

        # Seed the observation points more or less uniformly over the rectangular domain.
        # TODO: avoid observations at the boundary.
        count = 0
        for x in range(obs_nx):
            for y in range(obs_ny):
                rows[count] = count
                cols[count] = glo_idx[round(float(x * nx) / float(obs_nx)),
                                      round(float(y * ny) / float(obs_ny))]
                count += 1

        assert count == num_obs
        H = scipy.sparse.csr_matrix((vals, (rows, cols)), shape=(num_obs, conf.problem_size))
        assert H.nnz == num_obs
    else:
        H = scipy.sparse.eye(conf.problem_size)
    return H


def InitialCovar(conf, glo_idx):
    """ Function computes the initial covariance matrix as function of exponential distance.
    """
    nx = conf.nx
    ny = conf.ny
    problem_size = conf.problem_size

    variance = float(conf.model_ini_var)
    P = np.zeros((problem_size, problem_size))

    # Here we express correlation distances in logical coordinates of nodal points.
    sx = max(nx * float(conf.model_ini_covar_radius) / float(conf.domain_size_x), 1)
    sy = max(ny * float(conf.model_ini_covar_radius) / float(conf.domain_size_y), 1)
    Rx = round(math.ceil(3.0 * sx))
    Ry = round(math.ceil(3.0 * sx))

    for u in range(nx):
        for v in range(ny):
            i = glo_idx[u,v]
            for x in range(u-Rx, u+Rx+1):
                if 0 <= x and x < nx:
                    dx = float(u-x)/sx
                    for y in range(v-Ry, v+Ry+1):
                        if 0 <= y and y < ny:
                            dy = float(v-y)/sy
                            j = glo_idx[x,y]
                            covariance = variance * math.exp(-0.5 * float(dx*dx + dy*dy))
                            P[i,j] = covariance
                            P[j,i] = covariance

    P = (P + np.transpose(P)) * 0.5             # symmetrize
    return P


def DataAssimilation(conf, glo_idx, true_fields):
    """ Function runs the main part of the application - simulation with data assimilation.
    """
    start_time = timer()
    print('\nData assimilation ...')

    # Initially the state field is all zeros because we have no idea what is the "true" one.
    field = np.zeros((conf.nx, conf.ny))
    P = InitialCovar(conf, glo_idx)
    H = ComputeH(conf, glo_idx)
    rel_diff = np.zeros((conf.Nt,))
    writer = None
    if conf.generate_text_output: writer = Writer(conf, conf.simulation_solution)

    # Time integration with data assimilation.
    hFigure = None
    for k in range(conf.Nt):
        # Print the progress.
        print('.', end='', flush=True)
        if (k+1) % 80 == 0: print()
        # Write the field in textual form.
        if conf.generate_text_output: writer.WriteField(field, k, k * conf.dt)
        # Plot the current images (both truth and estimation).
        image = np.hstack((np.transpose(true_fields[:,:,k]),
                           np.zeros((field.shape[1],10)), np.transpose(field)))
        vmin = np.amin(image)
        vmax = np.amax(image)
        image = np.abs((image - vmin) / (vmax - vmin + np.finfo(float).eps))
        if hFigure is None:
            hFigure = plt.imshow(image)
        else:
            hFigure.set_data(image)
        plt.ylim(0, image.shape[0])     # flip Y
        plt.title('True density < vs. > Kalman filtering, t=' + str(k),
                    fontsize=10, fontweight='bold')
        plt.pause(0.05)
        plt.draw()
        # Write image file.
        if conf.generate_video:
            image = np.flipud(image)    # flip Y
            scipy.misc.imsave(conf.output_dir + '/both_fields' + format(k,'05d') + '.png', image)

        Q = ComputeQ(conf)
        R = ComputeR(conf)

        # Compute inverse A (i.e. the inverse model matrix).
        t = k * conf.dt                                     # physical time
        invA = InverseModelMatrix(conf, glo_idx, t)         # note, we compute A^{-1} not A

        # x_prior <- A * x_previous.
        fld1D = np.reshape(field, (conf.problem_size, 1))   # flatten the state field 'x'
        field = scipy.sparse.linalg.spsolve(invA, fld1D)    # x_{t+1} = invA^{-1} * x_t
        field = np.reshape(field, (conf.nx, conf.ny))       # back to 2D field
        field = ApplyBoundaryCondition(conf, field)         # honor the boundary conditions

        # P_prior <- A * P_previous * A^t + Q = invA^{-1} * (invA^{-1} * P_previous)^t + Q.
        P = scipy.sparse.linalg.spsolve(invA, P)
        P = scipy.sparse.linalg.spsolve(invA, P.transpose())
        P += Q
        P = (P + np.transpose(P)) * 0.5             # symmetrize

        # Residual between the "true" observations and the currently simulated ones.
        z = H.dot( np.reshape(true_fields[:,:,k], (conf.problem_size, 1)) )
        y = z - H.dot(fld1D)                        # using flatten 'x_previous'

        # S = H * P_prior * H^t + R
        # N O T E: for some reason P.dot(H.transpose()) does not work, but H.dot(P) does.
        HP = H.dot(P)
        PHt = HP.transpose()
        S = H.dot(PHt)
        S += R
        S = (S + np.transpose(S)) * 0.5             # symmetrize

        # x = x_prior + K * y = x_prior * P * H^t * S^{-1} * y.
        fld1D = np.reshape(field, (conf.problem_size, 1))               # flatten 'x_prior'
        fld1D += PHt.dot( scipy.linalg.solve(S, y) )
        field = np.reshape(fld1D, (conf.nx, conf.ny))                   # back to 2D field

        # P_new = (I - K*H) * P_prior = P_prior - P_prior * H^t * S^{-1} * H * P_prior.
        P -= PHt.dot( scipy.linalg.solve(S, HP) )
        P = (P + np.transpose(P)) * 0.5             # symmetrize

        # Compute relative difference: "true" field vs. data assimulation solution.
        diff = np.linalg.norm((field - true_fields[:,:,k]).ravel())
        norm = np.linalg.norm(true_fields[:,:,k].ravel())
        rel_diff[k] = diff / norm

    if writer is not None: del writer
    plt.close()
    if conf.generate_video:
        MakeVideo(conf, "true_field")
        MakeVideo(conf, "both_fields")
        np.savetxt(conf.output_dir + '/rel_diff.txt', rel_diff)
    print('\nexecution time: ' + str(timer() - start_time) + ' seconds\n\n')


###################################################################################################
# Utilities.
###################################################################################################

class Writer:
    """ Class for writing simulation result in textual form.
    """
    def __init__(self, conf, filename):
        """ Constructor
        """
        assert isinstance(filename, str)
        self.fid = None
        # N O T E: for some strange reason np.savetxt() expects file opened in binary format 'wb'.
        self.fid = open(conf.output_dir + '/' + filename, 'wb')

    def __del__(self):
        """ This method is not a destructor, is just a normal method you can call whenever
            you want to perform any operation, but it is always called before the garbage
            collector destroys the object.
            https://stackoverflow.com/questions/37852560/is-del-really-a-destructor
        """
        if self.fid is not None:
            self.fid.flush()
            self.fid.close()
            self.fid = None

    def WriteField(self, field, discrete_time, physical_time):
        """ Function appends the new field to the output file.
        """
        assert isinstance(discrete_time, int)   and discrete_time >= 0
        assert isinstance(physical_time, float) and physical_time >= 0.0
        nr = field.shape[0]
        nc = field.shape[1]
        rows = np.repeat(np.arange(nr), nc)
        cols = np.tile(np.arange(nc), nr)
        data = np.column_stack((rows, cols, field.flatten()))
        np.savetxt(self.fid, np.reshape(np.array([discrete_time, physical_time]), (1,2)),
                    fmt='%5d %.7f', delimiter='\t', newline='\n');
        np.savetxt(self.fid, data, fmt='%d %d %.10f', delimiter='\t', newline='\n');


def CreateAndCleanOutputDir(conf):
    """ Function creates an output directory inside the current one,
        which is supposed to be the project root folder.
    """
    # Create the output directory if not existent.
    if not os.path.exists(conf.output_dir):
        try:
            os.makedirs(conf.output_dir)
        except OSError as ex:
            if ex.errno != errno.EEXIST:     # guard against race condition
                raise
    # Clear the output directory.
    filelist = [ f for f in os.listdir(conf.output_dir)
        if f.endswith(".avi") or f.endswith(".png") or f.endswith(".jpg") or f.endswith(".txt") ]
    for f in filelist:
        os.remove(conf.output_dir + "/" + f)


def MakeVideo(conf, filetitle):
    """ Function creates a single video file from a sequence of field states
        written into image files.
    """
    print("\n\n\n")
    framerate = 24
    if os.system("ffmpeg -y -f image2 -framerate " + str(framerate) + " -pattern_type glob -i '" +
                    conf.output_dir + "/" + filetitle + "*.png' " +
                    conf.output_dir + "/" + filetitle + ".avi"):
        print("WARNING: unable to write video: ffmpeg failed")
    print("\n\n\n")


if __name__ == '__main__':
    try:
        if len(sys.argv) == 2:
            config_file = sys.argv[1]
            assert isinstance(config_file, str), "script agrument must be a file name string"
        else:
            config_file = 'amdados.conf'
        Amdados2D(config_file)
    except Exception as error:
        traceback.print_exc()
        print('ERROR: ' + str(error.args))


