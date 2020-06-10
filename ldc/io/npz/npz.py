""" Numpy generic I/O 
"""
import numpy as np

def decode(array):
    if array.size==1:
        if array.dtype.kind=='U':
            array = str(array)
        else:
            array = float(array)
    return array

class NPZ:
    """Provides a suite of I/O routine to load and save any quantities in
    Numpy file format.
    """

    def __init__(self, filename):
        """ Initialize NPY I/O
        """
        self.filename = filename
        
    def save_array(self, arr, name="data", **kwargs):
        """ Write numpy array to npy file
        """
        eval("np.savez(self.filename, %s=arr, **kwargs)"%name) 

    def load_array(self, name="data", full_output=True):
        """ Load array from hdf5 file with all its data set attributes
        """
        attr = dict()
        with np.load(self.filename) as data:
            arr = data[name]
            attr = [(k,decode(v)) for k,v in data.items() if k!=name]
        dattr = dict()
        attr = dict(attr)
        return arr, dict(attr)
