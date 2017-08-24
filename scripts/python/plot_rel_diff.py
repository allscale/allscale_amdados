#------------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017
#------------------------------------------------------------------------------

"""
Script plots the relative difference between the true density field
and data assimilation solution. It expects the file './output/rel_diff.txt'
produced by Amdados2D.py or C++ implementation.
"""
print(__doc__)

import traceback
import matplotlib
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    try:
        data = np.genfromtxt('./output/rel_diff.txt', delimiter=' ')
        assert len(data.shape) == 1 or data.shape[1] == 1, "wrong number of data columns"
        matplotlib.rcParams['legend.fontsize'] = 10
        fig = plt.figure()
        ax = fig.gca()
        ax.plot(data)       # , label='some label text')
        # ax.legend()
        ax.set_xlabel('number of time iterations')
        ax.set_ylabel('relative difference')
        plt.title('||(true field) - (assimilation solution)|| / ||true field||')
        plt.savefig('./output/rel_diff.png')
        plt.show()
    except Exception as error:
        traceback.print_exc()
        print('ERROR: ' + str(error.args))

