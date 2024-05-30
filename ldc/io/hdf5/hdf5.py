"""Provides a suite of I/O routine to load and save any quantities in
HDF5 file format.

"""
import h5py
import numpy as np
from astropy import units
import ast

def str_encode(value):
    """ Encode value to ascii if string

    >>> str_encode("hello")
    b'hello'
    """
    if isinstance(value, (str)):#, bytes)):
        return value.encode("ascii", "ignore")
    return value

def str_decode(value):
    """ Decode value if string

    >>> str_decode(b'hello')
    'hello'
    """
    if isinstance(value, (bytes)):
        return value.decode()
    return value

def atomic(cfg):
    """Flatten lists and split quantities into data (for values) and
    attributes (for units).

    >>> cfg = {"deltat":[1*units.Unit('s'), 2*units.Unit("s")]}
    >>> data, units = atomic(cfg)
    >>> print(data)
    {'deltat_0': 1.0, 'deltat_1': 2.0}
    >>> print(units)
    {'deltat_0': 's', 'deltat_1': 's'}
    """
    flat = {}
    for k,v in cfg.items():
        if isinstance(v, list) or isinstance(v, tuple):
            for i in range(len(v)):
                flat[k+f"_{i}"] = v[i]
        elif isinstance(v, dict):
            for ki,vi in v.items():
                flat[k+"_"+ki] = vi
        else:
            flat[k] = v
    data = {}
    attr = {}
    for k,v in flat.items():
        if isinstance(v, units.Quantity):
            data[k] = v.value
            attr[k] = str(v.unit)
        else:
            data[k] = v
    return data, attr


def encode_utype(array):
    """ Replace utype column in numpy array by binary format.

    >>> encode_utype(np.rec.fromarrays([["a", "b", "c"], [1, 2, 3]], names=["name", "val"]))
    rec.array([(b'a', 1), (b'b', 2), (b'c', 3)],
              dtype=[('name', 'S1'), ('val', '<i8')])
    """
    sizeof_numpy_unicode_char = np.dtype('U1').itemsize

    if array.dtype.fields:
        new_dtype = [(n, dt[0]) if dt[0].kind != "U"
                     else (n, np.dtype('<S%d'%(dt[0].itemsize//sizeof_numpy_unicode_char)))
                     for n, dt in array.dtype.fields.items()]
        array = array.astype(new_dtype)
    return array

def decode_utype(array):
    """ Replace btype column in numpy array by unicode format.

    >>> decode_utype(np.rec.fromarrays([[b"a", b"b", b"c"], [1, 2, 3]], names=["name", "val"]))
    rec.array([('a', 1), ('b', 2), ('c', 3)],
              dtype=[('name', '<U1'), ('val', '<i8')])
    """
    sizeof_numpy_unicode_char = np.dtype('S1').itemsize

    if array.dtype.fields:
        new_dtype = [(n, dt[0]) if dt[0].kind != "S"
                     else (n, np.dtype('<U%d'%(dt[0].itemsize//sizeof_numpy_unicode_char)))
                     for n, dt in array.dtype.fields.items()]
        array = array.astype(new_dtype)
    return array



def save_array(filename, arr, name="data", mode="a", chunks=False, **kwargs):
    """ Write numpy array to hdf5 file

    mode is 'a' (default) or 'w'

    >>> save_array("test.h5", np.ones((5)), mode="w", author='me')
    >>> print(load_array("test.h5"))
    (array([1., 1., 1., 1., 1.]), {'author': 'me'})
    """
    arr = encode_utype(arr)
    with h5py.File(filename, mode) as fid:
        # TODO PB: check here for Yorsh cat
        if len(arr.shape) == 1:
            arr = arr.reshape(len(arr), 1)
        if chunks:
            fid.create_dataset(name, data=arr, chunks=chunks,
                               maxshape=(arr.shape[0], None))
        else:
            fid.create_dataset(name, data=arr)
        for k, v in kwargs.items():
            fid[name].attrs.create(k, str_encode(v))

def append_array(filename, arr, column_index, name="data"):
    """ Append array to existing file and data set.

    >>> save_array("test.h5", np.ones((5)), mode="w", chunks=True)
    >>> append_array("test.h5", np.ones((5)), column_index=1)
    >>> print(load_array("test.h5"))
    (array([[1., 1.],
           [1., 1.],
           [1., 1.],
           [1., 1.],
           [1., 1.]]), {})
    """
    if column_index == 0:
        save_array(filename, arr, name=name, chunks=True)
    with h5py.File(filename, 'a') as fid:
        fid[name].resize(fid[name].shape[1]+1, axis=1)
        fid[name][:, -1] = arr


def load_array(filename, name="", full_output=True, sl=None):
    """Return array and its attributes from hdf5 file.

    if full_output is True, return array and meta data as dict.
    Otherwise, return array only.

    one can select a subset of data by giving a slice object as sl
    argument.

    """
    with h5py.File(filename, "r") as fid:
        names = [name] if name else list(fid.keys())
        attr = {}
        arrs = []
        for name in names:
            dset = fid.get(name)
            for k, v in dset.attrs.items():
                attr[k] = str_decode(v)
            if sl is not None:
                dset = dset[sl]
            arrs.append(decode_utype(np.array(dset)).squeeze())
        if len(names)>1:
            try: # make a rec array if all arrays share same size
                arr = np.rec.fromarrays(arrs, names=names)
            except ValueError: # make a dict otherwise
                arr = dict(zip(names, arrs))
        else:
            arr = arrs[0]
        if full_output:
            return arr, attr
        return arr

def load_attributes(filename, name="data"):
    """ Return attributes from data set.
    """
    with h5py.File(filename, "r") as fid:
        dset = fid.get(name)
        attr = {}
        for k, v in dset.attrs.items():
            attr[k] = str_decode(v)
        return attr

def save_config(filename, cfg, name="config", mode='a'):
    """ Save configuration file in dedicated group.

    >>> save_config("test.h5", {'author':'me', 'date':'today'})
    >>> print(load_config("test.h5"))
    {'author': 'me', 'date': 'today'}
    """
    data, units = atomic(cfg)
    with h5py.File(filename, mode) as fid:
        fid.create_group(name)
        for k, v in data.items():
            ds = fid[name].create_dataset(k, data=str_encode(v))
            if k in units.keys():
                ds.attrs.create('units', units[k])


def load_config(filename, name="config", deprecated=False):
    """ Load configuration from dedicated group.
    """
    if deprecated:
        return load_config_deprecated(filename, name=name)
    with h5py.File(filename, "r") as fid:
        grp = fid.get(name)
        data = {}
        attr = {}
        for name in grp.keys():
            data[name] = str_decode(grp.get(name)[()])
            if 'units' in grp[name].attrs.keys():
                attr[name] = grp[name].attrs["units"]
    for k,v in attr.items():
        data[k] *= units.Unit(str(v))
    return data

def quantity_decode_deprecated(value):
    """ Decode astropy quantity
    """
    try:
        value = ast.literal_eval(value)
        if isinstance(value, dict) and 'units' in value.keys():
            value = value["value"]*units.Unit(value["units"])
    except:
        pass
    return value

def load_config_deprecated(filename, name="config"):
    """ Load configuration from dedicated group.
    """
    with h5py.File(filename, "r") as fid:
        grp = fid.get(name)
        attr = dict()
        for k, v in grp.attrs.items():
            attr[k] = quantity_decode_deprecated(str_decode(v))
    return attr



def display(filename):
    """ Display file content

    >>> save_array("test.h5", np.ones((5)), mode="w")
    >>> display("test.h5")
    data <HDF5 dataset "data": shape (5, 1), type "<f8">
    """
    def print_attrs(name, obj):
        print(name, obj)
        for key, val in obj.attrs.items():
            print("    %s: %s" % (key, val))
    with h5py.File(filename, "r") as fid:
        fid.visititems(print_attrs)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
