#------------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017
#------------------------------------------------------------------------------

import re

class Configuration:
    """ Class reads and keeps the application parameters.
    """
    def __init__(self, filename):
        """ Constructor reads and parses configuration file of application parameters.
        """
        assert isinstance(filename, str), "expecting file name string"
        assert filename, "empty file name"
        with open(filename) as config_file:
            for line in config_file:
                line_nocomments = re.sub(re.compile("#.*?$"), "", line.strip())
                if line_nocomments:
                    name, var = line_nocomments.partition(" ")[::2]
                    name = name.strip()
                    var = var.strip()
                    try:
                        setattr(self, name, float(var))
                    except ValueError:
                        setattr(self, name, var)    # if not a float, then a string

    def PrintParameters(self):
        for attr, val in vars(self).items(): print('{} : {}'.format(attr, val))

