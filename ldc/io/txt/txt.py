"""Provides a suite of I/O routine to load and save any quantities in
ascii file format.
"""
import numpy as np

def save_array(filename, arr, name='data'):
    """ Write numpy array to npy file

    >>> save_array("test.txt.gz", np.ones((5)))
    >>> print(load_array("test.txt.gz")["data"])
    [1. 1. 1. 1. 1.]
    >>> A = np.rec.fromarrays([["s1", "s2"], np.ones((2))], names=["label", "value"])
    >>> save_array("test.txt.gz", A)
    >>> print(load_array("test.txt.gz"))
    [('s1', 1.) ('s2', 1.)]
    """
    header = name
    fmt = '%.18e'
    if isinstance(arr.dtype.names, tuple):
        header = "\t".join(arr.dtype.names)
        fmt = ""
        for t in arr.dtype.descr:
            if 'U' in t[1]:
                fmt += "\t%s"
            else:
                fmt += "\t%.18e"
    np.savetxt(filename, arr, header=header, comments='#', fmt=fmt)

def load_array(filename):
    """ Load array from hdf5 file with all its data set attributes
    """
    arr = np.recfromtxt(filename, encoding='utf-8', names=True)
    return arr

if __name__ == "__main__":
    import doctest
    doctest.testmod()
