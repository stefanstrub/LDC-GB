"""Provides a suite of I/O routine to load and save any quantities in
HDF5 file format.

"""
import os
import h5py
import numpy as np

def str_encode(value):
    """ Encode value to ascii if string 
    """
    if isinstance(value, str) or isinstance(value, bytes):
        return value.encode("ascii", "ignore")
    else:
        return value

def str_decode(value):
    """ Decode value if string 
    """
    if isinstance(value, str) or isinstance(value, bytes):
        return value.decode()
    else:
        return value
    

def save_array(filename, arr, name="data", mode="a", **kwargs):
    """ Write numpy array to hdf5 file

    mode is 'a' (default) or 'w'
    """
    with h5py.File(filename, mode) as fid:
        if len(arr.shape)==1:
            arr = arr.reshape(len(arr), 1)

        dset = fid.create_dataset(name, data=arr, chunks=True,
                                maxshape=(arr.shape[0], None))
        for k,v in kwargs.items():
            fid[name].attrs.create(k, str_encode(v))

def append_array(filename, arr, column_index, name="data"):
    """ Append array to existing file and data set. 
    """
    with h5py.File(filename, 'a') as fid:
        fid[name].resize(fid[name].shape[1]+1, axis=1)
        f[name][:,-1] = arr


def load_array(filename, name="data", full_output=True):
    """ Return array and its attributes from hdf5 file.
    """
    with h5py.File(filename, "r") as fid:
        dset = fid.get(name)
        attr = {}
        for k,v in dset.attrs.items():
            attr[k] = str_decode(v)
        if full_output:
            return np.array(dset).squeeze(), attr
        else:
            return np.array(dset).squeeze()

def load_attributes(filename, name="data"):
    """ Return attributes from data set. 
    """
    with h5py.File(filename, "r") as fid:
        dset = fid.get(name)
        attr = {}
        for k,v in dset.attrs.items():
            attr[k] = str_decode(v)
        return attr

def save_config(filename, cfg, name="config"):
    """ Save configuration file in dedicated group.
    """
    with h5py.File(filename, "a") as fid:
        fid.create_group(name)
        for k,v in cfg.items():
            fid[name].attrs.create(k, str_encode(v))


def load_config(filename, name="config"):
    """ Load configuration from dedicated group. 
    """
    with h5py.File(filename, "r") as fid:
        gp = fid.get(name)
        attr = dict()
        for k,v in gp.attrs.items():
            attr[k] = str_decode(v)
    return attr
