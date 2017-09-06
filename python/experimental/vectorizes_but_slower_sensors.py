#------------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017
#------------------------------------------------------------------------------

#import pdb; pdb.set_trace()           # enables debugging

#import sys                              # ; sys.path.append("./python")
import matplotlib                       # ; matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import traceback
#import os
#import re
#import getopt
import math
import numpy as np
#import scipy
import scipy.misc
#import scipy.sparse
#import scipy.sparse.linalg
#from timeit import default_timer as timer
import random

# Domain sizes:
Nx = int(127)
Ny = int(97)
# Fraction of points occupied by sensors:
FRACTION = 0.01
# Weight on mutual expulsion of sensors points:
WEIGHT = 2.0
# Number of sensors:
Ns = int(min(max(int(math.ceil(FRACTION * Nx * Ny)), int(1)), int(Nx * Ny)))

def InitialSensorDistribution():
    """ Function generates initial space distribution of sensor points.
    """
    points = np.zeros((2,Ns))
    for k in range(Ns):
        points[0,k] = random.uniform(0.0, 1.0)
        points[1,k] = random.uniform(0.0, 1.0)
    return points


def Evaluate(points):
    """ Function evaluates objective function and its gradient.
    """
    N = points.shape[1]
    EPS = math.sqrt(np.finfo(float).eps)

    # Reciprocal distances to subdomain borders.
    r1 = 1.0 / (       points **2 + EPS)
    r2 = 1.0 / ((1.0 - points)**2 + EPS)

    J = np.sum(r1) + np.sum(r2)
    gradJ = np.zeros((2,N))

    for i in range(N):
        g = np.zeros((2,))
        for j in range(N):
            d = points[:,i] - points[:,j]
            sqdist = float(np.dot(d,d) + EPS)
            J += 1.0 / sqdist
            g -= d / sqdist**2

        gradJ[:,i] = 2.0 * (g - points[:,i] * r1[:,i]**2 + (1.0 - points[:,i]) * r2[:,i]**2)

    return J, gradJ


def DistributeSensors():
    """ Function evenly distributes sensors in the domain.
        It uses gradient descent to minimize objective function, which measures
        the mutual sensor points' expulsion and their expulsion from the borders.
        Experiments showed that quadratic distances (used here) work better
        than absolute values of the distances.
    """
    EPS = math.sqrt(np.finfo(float).eps)
    TINY = np.finfo(float).tiny / np.finfo(float).eps**3
    points = np.random.rand(2,Ns)
    ini_points = np.copy(points)
    J, gradJ = Evaluate(points)
    step = 0.1
    tol = 0.1

    print('Initial J = ' + str(J))
    proceed = True
    while proceed:
        print('.', end='', flush=True)
        points_new = points - step * gradJ
        if any(points_new.ravel() < 0.0) or any(points_new.ravel() > 1.0):
            step *= tol
            continue
        J_new, gradJ_new = Evaluate(points_new)
        if J < J_new:
            step *= tol
            continue
        proceed = (J - J_new > np.finfo(float).eps * J) and (step > TINY)
        np.copyto(points, points_new)
        J = J_new
        np.copyto(gradJ, gradJ_new)
        step *= 2.0
    print('\n  Final J = ' + str(J))
    print('\n\n')

    plt.figure(1)
    plt.subplot(121, autoscale_on=False, aspect=float(Nx)/float(Ny))
    plt.plot(ini_points[0,:] * (Nx-1), ini_points[1,:] * (Ny-1), 'o')
    plt.axis([0, Nx, 0, Ny])
    plt.title('Before')
    plt.subplot(122, autoscale_on=False, aspect=float(Nx)/float(Ny))
    plt.plot(points[0,:] * (Nx-1), points[1,:] * (Ny-1), 'o')
    plt.axis([0, Nx, 0, Ny])
    plt.title('After')
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    try:
        random.seed()
        DistributeSensors()
        #input('Press the "ENTER" key to quit ... ')
    except Exception as error:
        traceback.print_exc()
        print('ERROR: ' + str(error.args))


