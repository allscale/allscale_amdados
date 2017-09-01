#------------------------------------------------------------------------------
# Author    : Albert Akhriev, albert_akhriev@ie.ibm.com
# Copyright : IBM Research Ireland, 2017
#------------------------------------------------------------------------------

import sys
import numpy as np
import Configuration


class Observations:
    """ Class for writing and reading observations.
    """
# public:
    def __init__(self, conf):
        """ Constructor
        """
        self.mode = None        # reading 'r' or writing 'w'
        self.filename = None    # name of input or output observation file
        self.fields = []        # sequence of (time,field) tuples (writing) or just fields (reading)
        self.timestamps = None; # sequence of timestamps (reading)

        self.filename = conf['output_dir'] + '/' + 'observations.bin'

# public:
    def __del__(self):
        """ This method is not a destructor, is just a normal method you can call whenever
            you want to perform any operation, but it is always called before the garbage
            collector destroys the object.
            https://stackoverflow.com/questions/37852560/is-del-really-a-destructor
        """
        # Write previously accumulated fields into the file.
        if self.mode == 'w' and len(self.fields) > 0:
            fid = None
            try:
                nx = self.fields[0][1].shape[0]
                ny = self.fields[0][1].shape[1]
                fid = open(self.filename, 'wb')
                np.array([np.float32(len(self.fields))]).tofile(fid)    # number of records
                for z in self.fields:
                    t = z[0]
                    field = z[1]
                    assert nx == field.shape[0]         # all the fields have the same size
                    assert ny == field.shape[1]
                    header = np.array([np.float32(5555555), np.float32(t),
                                       np.float32(nx), np.float32(ny)], dtype=np.float32)
                    header.tofile(fid)                              # header
                    field.astype(np.float32).tofile(fid)            # field
                    np.array([np.float32(7777777)]).tofile(fid)     # footer
                fid.flush()
                fid.close()
                print('***** file of observations has been written *****')
            except IOError as ex:
                print(ex)
                if fid is not None: fid.close()
                sys.exit('failed to write to observations file')

# public:
    def Write(self, t, field):
        """ Function accumulates a sequence of fields and their timestamps.
            Upon program termination, the sequence will be written in the file.
        """
        if self.mode is None: self.mode = 'w'
        assert self.mode == 'w', "function Read() was invoked prior to Write(), this is not allowed"
        assert isinstance(t, float) and t >= 0.0
        assert len(field.shape) == 2 and field.dtype == float
        self.fields.append(tuple((np.float32(t), field.astype(np.float32))))

# public:
    def Read(self):
        """ Function reads all the records and save in this object the sequence
            of fields and corresponding sequence of timestamps.
        """
        if self.mode is None: self.mode = 'r'
        assert self.mode == 'r', "function Write() was invoked prior to Read(), this is not allowed"
        self.timestamps = None
        self.fields = None
        nx = 0
        ny = 0
        fid = None
        try:
            fid = open(self.filename, 'rb')
            num_records = np.fromfile(fid, dtype=np.float32, count=1)
            num_records = int(round(num_records[0]))

            for rec in range(num_records):
                header = np.fromfile(fid, dtype=np.float32, count=4)
                if int(round(header[0])) != 5555555:
                    raise IOError("header checkpoint mismatch")

                # Create the output arrays upon reading the first record.
                if rec == 0:
                    nx = int(round(header[2]))
                    ny = int(round(header[3]))
                    self.timestamps = np.zeros((num_records,))
                    self.fields = np.zeros((num_records, nx, ny))
                else:
                    assert nx == int(round(header[2]))
                    assert ny == int(round(header[3]))

                field = np.fromfile(fid, dtype=np.float32, count=nx*ny)
                footer = np.fromfile(fid, dtype=np.float32, count=1)
                if int(round(footer[0])) != 7777777:
                    raise IOError("footer checkpoint mismatch")

                self.timestamps[rec] = float(header[1])
                self.fields[rec,:,:] = field.reshape((nx, ny)).astype(float)

            fid.close()
            assert all(np.sort(self.timestamps, axis=None) == self.timestamps), \
                        "timestamps must be sorted"
        except IOError as ex:
            print(ex)
            if fid is not None: fid.close()
            sys.exit('failed to read from observations file')

# public:
    def GetObservation(self, t):
        """ Function returns a fields of observations at time 't', possibly
            interpolating from two fields adjacent in time.
        """
        assert (self.mode is not None) and self.mode == 'r', "reading mode is expected"
        assert self.timestamps is not None, "observation file must be read beforehand"
        N = len(self.timestamps)
        eps = np.finfo(float).eps
        assert -5*eps <= t and t <= self.timestamps[-1] * (1 + 5*eps), \
                "the timestamp requested is outside observations' range"
        idx = np.searchsorted(self.timestamps, t)
        if idx == 0:
            field = self.fields[0,:,:]
        elif idx == N:
            field = self.fields[-1,:,:]
        else:
            t1 = self.timestamps[idx-1]
            t2 = self.timestamps[idx]
            assert t2 > t1 and t1 <= t and t <= t2
            a = (t2 - t) / float(t2 - t1)
            b = (t - t1) / float(t2 - t1)
            field = a * self.fields[idx-1,:,:] + b * self.fields[idx,:,:]
        assert len(field.shape) == 2
        return field


