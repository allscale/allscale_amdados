#import pdb; pdb.set_trace()           # enables debugging

import sys
import numpy as np
import matplotlib; matplotlib.use('TkAgg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # Path to kalman filter output file is specified for Linux and Mac OS.
    # For Windows it is not clear what might be the best practice.
    assert len(sys.argv) == 2, "missed path to 'kalman_test.out' file"
    data = np.genfromtxt(sys.argv[1], delimiter=' ')
    assert data.shape[1] == 11, "wrong number of columns in the data matrix"

    t = data[:,0]                   # time-line
    x = data[:,1]                   # true x-coordinate
    y = data[:,2]                   # true y-coordinate
    z = data[:,3]                   # true z-coordinate
    mx = data[:,4]                  # measured noisy x-coordinate
    my = data[:,5]                  # measured noisy y-coordinate
    mz = data[:,6]                  # measured noisy z-coordinate
    x_est = data[:,7]               # Kalman filter estimated x-coordinate
    y_est = data[:,8]               # Kalman filter estimated y-coordinate
    z_est = data[:,9]               # Kalman filter estimated z-coordinate
    cov   = data[:,10]              # mean covariance

    matplotlib.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(x, y, z, label='"true" curve')
    ax.plot(mx, my, mz, label='noisy measurements')
    ax.plot(x_est, y_est, z_est, label='Kalman filter estimation')
    ax.legend()
    plt.show()

