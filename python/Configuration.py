# -----------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017
# -----------------------------------------------------------------------------

import re, os

class Configuration:
    """ Class reads, writes and keeps the application parameters.
    """

    def __init__(self, filename):
        """ Constructor reads and parses configuration file
            of application parameters.
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
                    try:  # try to convert to float ...
                        setattr(self, name, float(var))
                    except ValueError:  # ... otherwise to string
                        setattr(self, name, str(var))


    def WriteParameterFile(self, filename) -> str:
        """ Function writes modified configuration parameters into a new file.
        """
        assert isinstance(filename, str)
        config_file = os.path.join(self.output_dir, filename)
        with open(config_file, "wt") as f:
            for attr, val in vars(self).items():
                f.write("{} {}\n".format(attr, val))
        return config_file


    def PrintParameters(self):
        """ Prints configuration parameters.
        """
        SEP = "------------------------------------------------------------"
        print(SEP)
        print("Configuration parameters:")
        print(SEP)
        for attr, val in vars(self).items():
            print("{} : {}".format(attr, val))
        print(SEP)
        print("")
