"""
Script converts state fields into pictures given the folder for simulation results.
"""
print(__doc__)

#import pdb; pdb.set_trace()
import numpy as np
import scipy
import scipy.misc
import glob
import sys
import os
from string import digits

NEW_FORMAT = True

if __name__ == '__main__':
    try:
        TINY = np.finfo(np.float64).tiny / np.finfo(np.float64).eps**2
        BORDER = np.uint8(64)
        if len(sys.argv) != 2:
            raise ValueError('please, specify the data folder as a single application argument')

        files = glob.glob(sys.argv[1] + '/*.txt')
        files.sort()
        assert len(files) > 0, ('no state files were found in "' + sys.argv[1] + '"')

        for f in files:
            print('converting file: ' + f + ' ...')

            data = np.loadtxt(f, dtype=np.float64, comments='#')
            assert len(data.shape) == 1, 'wrong data layout'
            dim    = (np.floor(data[0] + 0.5)).astype(int)
            width  = (np.floor(data[1] + 0.5)).astype(int)
            height = (np.floor(data[2] + 0.5)).astype(int)
            assert dim == 2, "two-dimensional data is expected"
            v = np.reshape(data[3:], (width, height))

            maxVal = np.amax(v)
            if maxVal > TINY:
                v = (v * (255.0 / maxVal)).astype(np.uint8)
            else:
                v = np.zeros(v.shape, dtype=np.uint8)
                print('warning: max. value is too small: ' + str(maxVal))
            v[:, 0] = BORDER        # make a visible border
            v[:,-1] = BORDER
            v[ 0,:] = BORDER
            v[-1,:] = BORDER

            # Create a folder based on file title (without digits and extension).
            file_path, file_title = os.path.split(f)                        # get base name
            title = file_title.translate(str.maketrans('', '', digits))     # remove digits
            title = os.path.splitext(os.path.basename(title))[0]            # remove extension
            sub_path = file_path + '/' + title
            if not os.path.isdir(sub_path) or not os.path.exists(sub_path):
                if os.system('mkdir -p ' + sub_path):
                    sys.exit('!mkdir')

            # Save image file in the folder dependent on file title.
            scipy.misc.imsave(sub_path + '/' + file_title + '.png', v) #np.swapaxes(v,0,1))
        print('\n\n\n')
    except ValueError as error:
        print('ERROR: ' + str(error.args))
    except Exception as error:
        print('ERROR: ' + str(error.args))


