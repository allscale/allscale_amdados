#------------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017
#------------------------------------------------------------------------------

"""
Script plots the profile of a relative difference between the true density field
and data assimilation solution as a function of time.
It expects the name of a text file with single value per line.
"""
print(__doc__)

import traceback
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import os, sys

if __name__ == '__main__':
    try:
        assert len(sys.argv) == 2, "script argument must be a file name of a profile"
        filename = sys.argv[1]
        assert isinstance(filename, str), "script agrument must be a file name string"
        data = np.genfromtxt(filename, delimiter=' ')
        assert len(data.shape) == 1 or data.shape[1] == 1, "wrong number of data columns"
        matplotlib.rcParams['legend.fontsize'] = 10
        fig = plt.figure()
        ax = fig.gca()
        ax.plot(data)       # , label='some label text')
        # ax.legend()
        ax.set_xlabel('number of time iterations')
        ax.set_ylabel('relative difference')
        plt.title('||analytic - simulation|| / ||analytic||')
        plt.savefig(os.path.splitext(filename)[0] + '.png')
        plt.show()

        min_i = 1000
        max_i = 1035
        if data.size >= max_i:
            fig = plt.figure()
            ax = fig.gca()
            ax.plot(np.arange(min_i, max_i), data[min_i:max_i])
            # ax.legend()
            ax.set_xlabel('number of time iterations')
            ax.set_ylabel('relative difference')
            plt.title('Subrange of relative differences')
            plt.savefig(os.path.splitext(filename)[0] + '_subrange.png')
            plt.show()
    except Exception as error:
        traceback.print_exc()
        print('ERROR: ' + str(error.args))

