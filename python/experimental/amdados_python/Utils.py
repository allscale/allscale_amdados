import os, sys, errno
import math
import numpy as np
import scipy
import scipy.misc
import Configuration


def CreateAndCleanOutputDir(conf):
    """ Function creates an output directory inside the current one,
        which is supposed to be the project root folder.
    """
    output_dir = conf["output_dir"]
    # Create the output directory if not existent.
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        except OSError as ex:
            if ex.errno != errno.EEXIST:     # guard against race condition
                raise
    # Clear the output directory.
    try:
        filelist = [ f for f in os.listdir(output_dir)
                        if f.endswith(".avi") or f.endswith(".png") or f.endswith(".txt") ]
        for f in filelist:
            os.remove(output_dir + "/" + f)
    except Exception as ex:
        sys.exit(str(ex) + "\nfailed to remove files in output directory")


def MakeVideo(conf, filetitle, frame_size):
    """ Function creates a single video file from a sequence of field states
        written into image files.
    """
    print("\n\n\n")
    if conf["video_enable"] == 0: return
    framerate = conf["video_frame_rate"]
    assert isinstance(filetitle, str)
    assert len(frame_size) == 2
    tr_frame_size = (frame_size[1], frame_size[0])      # transposed !!!
    if os.system("ffmpeg -y -f image2 -framerate " + str(framerate) + " -pattern_type glob -i '" +
                    conf["output_dir"] + "/" + filetitle + "*.png' " +
                    "-s " + str(tr_frame_size[0]) + "x" + str(tr_frame_size[1]) + " " +
                    conf["output_dir"] + "/" + filetitle + ".avi"):
        print("WARNING: unable to write video: ffmpeg failed")
    print("\n\n\n")


def WriteProgress(conf, t, field, file_title = "field"):
    """ Function plots or/and save specified property field.
    """
    assert isinstance(t, float) and t >= 0.0
    assert isinstance(field, np.ndarray)
    filename = conf["output_dir"] + "/" + file_title + format(int(round(1000.0*t)), '09d') + ".png"
    tr_frame = np.transpose(field)
    scipy.misc.imsave(filename, tr_frame)


def Symmetrize(x):
    """ Function repairs matrix symmetry that might be lost due to round off errors.
    """
    assert isinstance(x, np.ndarray) and len(x.shape) == 2 and x.shape[0] == x.shape[1]
    x += x.transpose()
    x *= 0.5
    return x


def ColumnVec(x):
    """ Function creates a column vector out of numpy nd-array.
        Note, the shape changing convention is a bit messy with numpy.
    """
    isinstance(x, np.ndarray)
    y = np.reshape(x, (x.size,1))
    return y


