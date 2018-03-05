#------------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017
#------------------------------------------------------------------------------

"""
Script generates an AVI file from a sequence of images.
"""
print(__doc__)

import pdb; pdb.set_trace()
import sys, getopt, os, errno

def MakeVideo(folder, file_title, frame_rate, frame_size = None):
    """ Function creates a single video file from a sequence of field states
        written into image files.
    """
    assert isinstance(folder, str)
    assert isinstance(file_title, str)
    assert isinstance(frame_rate, int)
    if frame_size is not None:
        assert isinstance(frame_size, tuple) and len(frame_size) == 2
        ok = (os.system("ffmpeg -y -f image2 -framerate " + str(frame_rate) +
                        " -pattern_type glob -i '" +
                        folder + "/" + file_title + "*.png' " +
                        "-s " + str(frame_size[0]) + "x" + str(frame_size[1]) + " " +
                        folder + "/" + file_title + ".avi") == 0)
    else:
        ok = (os.system("ffmpeg -y -f image2 -framerate " + str(frame_rate) +
                        " -pattern_type glob -i '" +
                        folder + "/" + file_title + "*.png' " +
                        folder + "/" + file_title + ".avi") == 0)
    if not ok:
        print("WARNING: unable to write video: ffmpeg failed")


def PrintHelp():
    print("\nHelp:")
    print("python3 avi.py -d directory -o file_title -r frame_rate")
    print("All parameters are optional.")
    print("-d source/destination directory, default: output")
    print("-o output file title (i.e. name without extension), default: field")
    print("-r frame rate, default: 12")
    print("\n\n\n");


if __name__ == '__main__':
    try:
        print("\n\n\n")
        folder = "output"
        file_title = "field"
        frame_rate = 12
        frame_size = None
        try:
            opts, args = getopt.getopt(sys.argv[1:], "hd:o:r:")
        except getopt.GetoptError:
            PrintHelp()
            sys.exit(1)
        for opt, arg in opts:
            if opt == "-h":
                PrintHelp()
                sys.exit()
            elif opt == "-d":
                folder = arg
            elif opt == "-o":
                file_title = arg
            elif opt == "-r":
                try:
                    frame_rate = int(arg)
                    assert 1 <= frame_rate and frame_rate <= 100
                except ValueError:
                    sys.exit("frame rate must be a positive integer [1..100]")
            else:
                print("Unknown option: " + opt)
                PrintHelp()
                sys.exit(1)
        print("source/destination directory: " + folder)
        print("file title: " + file_title)
        print("frame rate: " + str(frame_rate))
        print("\n")
        MakeVideo(folder, file_title, frame_rate, frame_size);
        print("\n\n\n")
    except ValueError as error:
        print("ERROR: " + str(error.args))
    except Exception as error:
        print("ERROR: " + str(error.args))


