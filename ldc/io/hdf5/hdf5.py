""" HDF5 generic I/O 
"""
import h5py
import numpy as np

def str_encode(string):
    """ Encode str to HDF5 file. 
    """
    return string.encode("ascii", "ignore")

def str_decode(string):
    """
    """
    return string.decode()

class HDF5:
    """Provides a suite of I/O routine to load and save any quantities in
    HDF5 file format.
    """

    def __init__(self, filename):
        """ Initialize HDF5 I/O
        """
        self.filename = filename

        
    def save_array(self, arr, name="data", **kwargs):
        """ Write numpy array to hdf5 file
        """
        with h5py.File(self.filename, 'a') as fid:
            if len(arr.shape)==1:
                arr = arr.reshape(len(arr), 1)
                
            dset = fid.create_dataset(name, data=arr, chunks=True,
                                    maxshape=(arr.shape[0], None))
            for k,v in kwargs.items():
                if isinstance(v, str):
                    fid[name].attrs.create(k, str_encode(v))
                else:
                    fid[name].attrs.create(k, v)
                    

    def append_array(self, arr, column_index, name="data"):
        """ Append array to existing file and data set. 
        """
        with h5py.File(self.filename, 'a') as fid:
            fid[name].resize(fid[name].shape[1]+1, axis=1)
            f[name][:,-1] = arr


    def load_array(self, name="data", full_output=True):
        """ Load array from hdf5 file with all its data set attributes
        """
        with h5py.File(self.filename, "r") as fid:
            dset = fid.get(name)
            attr = {}
            for k,v in dset.attrs.items():
                if isinstance(v, str):
                    attr[k] = str_decode(v)
                else:
                    attr[k] = v
            if full_output:
                return np.array(dset), attr
            else:
                return np.array(dset)

    def load_attributes(self, name="data"):
        with h5py.File(self.filename, "r") as fid:
            dset = fid.get(name)
            attr = {}
            for k,v in dset.attrs.items():
                if isinstance(v, str):
                    attr[k] = str_decode(v)
                else:
                    attr[k] = v
            return attr
        
