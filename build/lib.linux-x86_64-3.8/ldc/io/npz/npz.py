"""Provides a suite of I/O routine to load and save any quantities in
Numpy file format.
"""
import numpy as np

def is_npy(filename):
    """ npy or npz
    """
    return filename.split(".")[-1] == "npy"

def decode(array):
    """Return float for 1-size array and str for U-type array.

    >>> decode(np.array("hello"))
    'hello'
    >>> decode(np.array([1.0]))
    1.0
    """
    if array.size == 1:
        if array.dtype.kind == 'U':
            array = str(array)
        else:
            array = float(array)
    return array

def save_array(filename, arr, name="data", **kwargs):
    """ Write numpy array to npy file

    >>> save_array("test.npz", np.ones((5)), author='me')
    >>> print(load_array("test.npz"))
    (array([1., 1., 1., 1., 1.]), {'author': 'me'})
    """
    kwargs[name] = arr
    np.savez(filename, **kwargs)

def load_array(filename, name="data", full_output=True):
    """ Load array from hdf5 file with all its data set attributes
    """
    attr = dict()

    if is_npy(filename):
        arr = np.load(filename)
    else:
        with np.load(filename) as data:
            arr = data[name]
            attr = [(k, decode(v)) for k, v in data.items() if k != name]
    if full_output:
        return arr, dict(attr)
    return arr

if __name__ == "__main__":
    import doctest
    doctest.testmod()
