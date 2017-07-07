import os
import errno
import math
import numpy as np
import scipy
import scipy.misc
import Configuration


def CreateOutputDir(conf):
    """ Function creates an output directory inside the current one,
        which is supposed to be the project root folder.
    """
    if not os.path.exists(conf["output_dir"]):
        try:
            os.makedirs(conf["output_dir"])
        except OSError as ex:
            if ex.errno != errno.EEXIST:     # guard against race condition
                raise


def MakeVideo(conf, filetitle, frame_size):
    """ Function creates a single video file from a sequence of field states
        written into image files.
    """
    framerate = 12
    assert isinstance(filetitle, str)
    if os.system("ffmpeg -y -f image2 -framerate " + str(framerate) + " -pattern_type glob -i '" +
                    conf["output_dir"] + "/" + filetitle + "*.png' " +
                    "-s " + str(frame_size[0]) + "x" + str(frame_size[1]) + " " +
                    conf["output_dir"] + "/" + filetitle + ".avi"):
        print("WARNING: unable to write video: ffmpeg failed")


def WriteProgress(conf, field, t, last_frame = False):
    """ Function plots or/and save specified property field.
    """
    assert isinstance(t, float) and t >= 0.0
    assert isinstance(field, np.ndarray)
    # secs = format(math.floor(t), '05d')
    # millisecs = format(round(1000.0 * (t - math.floor(t))), '03d')
    # filename = conf["output_dir"] + "/field_t=" + secs + "_" + millisecs + ".png"
    filename = conf["output_dir"] + "/field" + format(int(round(1000.0*t)), '09d') + ".png"
    tr_frame = np.transpose(field)
    scipy.misc.imsave(filename, tr_frame)
    print(".", end='', flush=True)
    if last_frame:
        print("\n\n\n")
        MakeVideo(conf, "field", tr_frame.shape)
        print("\n\n\n")


