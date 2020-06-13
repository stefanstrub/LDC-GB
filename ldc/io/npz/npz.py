"""Provides a suite of I/O routine to load and save any quantities in
Numpy file format.
"""
import numpy as np

def decode(array):
    if array.size==1:
        if array.dtype.kind=='U':
            array = str(array)
        else:
            array = float(array)
    return array

def save_array(filename, arr, name="data", **kwargs):
    """ Write numpy array to npy file
    """
    eval("np.savez(filename, %s=arr, **kwargs)"%name) 

def load_array(filename, name="data", full_output=True):
    """ Load array from hdf5 file with all its data set attributes
    """
    attr = dict()
    with np.load(filename) as data:
        arr = data[name]
        attr = [(k,decode(v)) for k,v in data.items() if k!=name]
    dattr = dict()
    attr = dict(attr)
    return arr, dict(attr)
