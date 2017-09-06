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
import random

Nx = int(127)
Ny = int(97)
FRACTION = 0.01
WEIGHT = 2.0
Ns = int(min(max(int(math.ceil(FRACTION * Nx * Ny)), int(1)), int(Nx * Ny)))

def InitialSensorDistribution():
    """
    """
    x = np.zeros((Ns,))
    y = np.zeros((Ns,))
    for k in range(Ns):
        x[k] = random.uniform(float(1),float(Nx-1))
        y[k] = random.uniform(float(1),float(Ny-1))
    return x, y


def Evaluate(x, y):
    """
    """
    assert x.size == y.size
    N = x.size
    EPS = math.sqrt(np.finfo(float).eps)
    W_div_N = 1 # WEIGHT / N

    J = 0.0
    gradJ = np.zeros((2,N))

    for i in range(N):
        # Reciprocal distances to subdomain borders.
        r_x1 = 1.0 / (       x[i] **2 + EPS)
        r_x2 = 1.0 / ((float(Nx) - x[i])**2 + EPS)
        r_y1 = 1.0 / (       y[i] **2 + EPS)
        r_y2 = 1.0 / ((float(Ny) - y[i])**2 + EPS)

        J += (r_x1**0.5 + r_x2**0.5 +
              r_y1**0.5 + r_y2**0.5)

        gx = 0.0
        gy = 0.0
        for j in range(N):
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            sqdist = float(dx*dx + dy*dy + EPS)
            J  += W_div_N / sqdist
            gx -= dx / sqdist**1.5
            gy -= dy / sqdist**1.5

        gradJ[0,i] = W_div_N * gx - x[i] * r_x1**1.5 + (float(Nx) - x[i]) * r_x2**1.5
        gradJ[1,i] = W_div_N * gy - y[i] * r_y1**1.5 + (float(Ny) - y[i]) * r_y2**1.5

    return J, gradJ


def DistributeSensors():
    """
    """
    EPS = math.sqrt(np.finfo(float).eps)
    x, y = InitialSensorDistribution()
    ini_x = np.copy(x)
    ini_y = np.copy(y)
    J, gradJ = Evaluate(x, y)
    step = 0.1
    tol = 0.1
    print('J = ' + str(J))
    proceed = True
    while proceed:
        print('step = ' + str(step))
        x1 = x - step * gradJ[0,:]
        y1 = y - step * gradJ[1,:]
        if any(x1 < 0) or any(x1 > Nx-1) or any(y1 < 0) or any(y1 > Ny-1):
            step *= tol
            continue
        J1, gradJ1 = Evaluate(x1, y1)
        print('J1 = ' + str(J1))
        if J1 > J:
            step *= tol
            continue
        proceed = (J - J1 > np.finfo(float).eps * J)
        np.copyto(x, x1)
        np.copyto(y, y1)
        J = J1
        np.copyto(gradJ, gradJ1)
        step *= 2.0
        print('.', end='', flush=True)
    print('\n\n')

    #print(x)
    #print(y)

    plt.figure(1)
    plt.subplot(121, autoscale_on=False, aspect=float(Nx)/float(Ny))
    plt.plot(ini_x, ini_y, 'o')
    plt.axis([0, Nx, 0, Ny])
    plt.title('Before')
    plt.subplot(122, autoscale_on=False, aspect=float(Nx)/float(Ny))
    plt.plot(x, y, 'o')
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


