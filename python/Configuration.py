#------------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017
#------------------------------------------------------------------------------

"""
Module for reading and parsing configuration parameters.
"""
print(__doc__)

import re

# Global instanse of configuration parameters.
__conf_params__ = None

def ReadConfiguration(filename):
    """ Function reads and parses configuration file of application parameters.
    """
    global __conf_params__
    assert isinstance(filename, str), "expecting file name string"
    assert filename, "empty file name"
    __conf_params__ = {}
    with open(filename) as config_file:
        for line in config_file:
            line_nocomments = re.sub(re.compile("#.*?$"), "", line.strip())
            if line_nocomments:
                name, var = line_nocomments.partition(" ")[::2]
                name = name.strip()
                var = var.strip()
                try:
                    __conf_params__[name] = float(var)
                except ValueError:
                    __conf_params__[name] = var         # if not a float, then a string

    print("\n\nConfiguration parameters:\n")
    for key, value in __conf_params__.items():
        print(str(key) + " : " + str(value) + "      [" + str(type(value).__name__) + "]")
    print("\n\n\n")
    return __conf_params__

def GetConfiguration():
    """ Function returns the global instanse of configuration parameters.
    """
    return __conf_params__
