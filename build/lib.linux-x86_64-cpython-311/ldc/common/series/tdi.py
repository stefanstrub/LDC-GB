""" Provide TDI facility using xarray dataset
"""

import numpy as np
import xarray as xr
import ldc.io.hdf5 as h5io
from ldc.common.series import TimeSeries, FrequencySeries

# AET matrix
aet_mat = np.array([[-1/np.sqrt(2), 0, 1/np.sqrt(2)],
                    [1/np.sqrt(6), -2/np.sqrt(6), 1/np.sqrt(6)],
                    [1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)]])

#pylint:disable=C0103

def XYZ2AET(X, Y, Z):
    A,E,T = ((Z - X)/np.sqrt(2.0),
            (X - 2.0*Y + Z)/np.sqrt(6.0),
            (X + Y + Z)/np.sqrt(3.0))
    if hasattr(X, "attrs"):
        A.attrs = X.attrs
    return A,E,T


class TDI(xr.Dataset):
    __slots__ = ()
    def XYZ2AET(self):
        """ Convert TDI XYZ into AET

        >>> XYZ = dict([(n,TimeSeries(np.random.randn(50), t0=0, dt=1, name=n)) for n in ["X", "Y", "Z"]])
        >>> tdi = TDI(XYZ)
        >>> tdi.XYZ2AET()
        """
        A,E,T = XYZ2AET(self["X"], self["Y"], self['Z'])
        A.attrs, E.attrs, T.attrs  = [self["X"].attrs]*3
        super().__init__(dict([(k,v) for k,v in zip(["A", "E", "T"], [A, E, T])]))

    def AET2XYZ(self):
        """ Convert TDI AET into XYZ
        """
        [X,Y,Z] = [-1/np.sqrt(2.0)*self["A"] +\
                   1/np.sqrt(6.0)*self["E"] + 1/np.sqrt(3.0)*self["T"],
                   -np.sqrt(2.0/3.0)*self["E"] + 1/np.sqrt(3.0)*self["T"],
                   1/np.sqrt(2.0)*self["A"] + 1/np.sqrt(6.0)*self["E"] +\
                   1/np.sqrt(3.0)*self["T"]]
        X.attrs, Y.attrs, Z.attrs  = [self["A"].attrs]*3
        super().__init__(dict([(k,v) for k,v in zip(["X", "Y", "Z"], [X, Y, Z])]))

    @staticmethod
    def load(filename, name='data', compound=True):
        """ Load tdi from file

        >>> XYZ = dict([(n,TimeSeries(np.random.randn(50), t0=0, dt=1, name=n)) for n in ["X", "Y", "Z"]])
        >>> tdi1 = TDI(XYZ)
        >>> tdi1.save("test.h5", mode='w')
        >>> tdi2 = TDI.load("test.h5")
        >>> tdi1.equals(tdi2)
        True
        """
        if compound:
            arr, attrs = h5io.load_array(filename, name=name)
            if attrs['coord']=='t':
                ds = dict([(k,TimeSeries(arr[k], t0=arr['t'][0], dt=arr['t'][1]-arr['t'][0], units=attrs["units"]))
                      for k in arr.dtype.names if k!='t'])
            elif attrs['coord']=='f':
                df = arr['f'][1] - arr['f'][0]
                ds = dict([(k,FrequencySeries(arr[k], kmin=int(arr['f'][0]//df), df=df, units=attrs["units"]))
                      for k in arr.dtype.names if k!='f'])
        else:
            ds = xr.open_dataset(filename, engine='h5netcdf', group=name)
        return TDI(ds)

    def save(self, filename, mode='a', name='data', compound=True):
        """ Save tdi to file

        >>> XYZ = dict([(n,TimeSeries(np.random.randn(50), t0=0, dt=1, name=n)) for n in ["X", "Y", "Z"]])
        >>> tdi = TDI(XYZ)
        >>> tdi.save("test2.h5", mode='w')
        """
        for k in self.data_vars:
            if self[k].attrs["units"] is None:
                self[k].attrs["units"] = 'dimensionless'

        if not compound:
            self.to_netcdf(filename, format='NETCDF4', engine='h5netcdf', mode=mode, group=name)
        else:
            attr = {}
            for k in self.data_vars:
                attr.update(self[k].attrs)
            attr["coord"] = list(self.coords.keys())[0]
            arr = np.rec.fromarrays([self.coords[k] for k in self.coords.keys()] +\
                                    [self[k] for k in self.keys()],
                                    names=list(self.coords.keys())+list(self.keys()))
            h5io.save_array(filename, arr, name=name, mode=mode, **attr)


def dot(a_mat, b_mat):
    """
    Perform the matrix multiplication of two list of matrices.

    Parameters
    ----------
    a : ndarray
        series of m x n matrices (array of size p x m x n)
    b : ndarray
        series of n x k matrices (array of size p x n x k)

    Returns
    -------
    c : ndarray
        array of size p x m x k containg the dot products of all matrices
        contained in a and b.


    """

    return np.einsum("ijk, ikl -> ijl", a_mat, b_mat)


def dot_vect(a_mat, b_vect):
    """[summary]

    Parameters
    ----------
    a : ndarray
        series of m x n matrices (array of size p x m x n)
    b : ndarray
        series of n x k vectors (array of size p x n)

    Returns
    -------
    c : ndarray
        array of size p x m containg the dot products of all vectors
        contained in a and b.


    """

    return np.einsum("ijk, ik -> ij", a_mat, b_vect)


def transform_covariance(tf, cov):

    return dot(tf, dot(cov, np.swapaxes(tf, 1, 2).conj()))


def delay_operator(fr, delta_t):
    return np.exp(-2j * np.pi  * fr * delta_t)


def compute_tdi_tf(fr, tt, tdi2=True, tdi_var='XYZ'):
    """Compute the TDI transfer function matrix

    Parameters
    ----------
    fr : ndarray
        frequency array [Hz]
    tt : dict
        dictionary of arm light travel times [s]
    gen : str, optional
        TDI generation, by default '2.0'

    Returns
    -------
    tdi_mat_eval
        TDI transfer function matrix for XYZ, array of size nf x 3 x 6
        Assumes that the interferometers are ordered as 1, 2, 3, 1', 2', 3'

    Raises
    ------
    ValueError
        TDI generation can only be 1.5 or 2.0.
    """

    d1 = delay_operator(fr, tt[23])
    d2 = delay_operator(fr, tt[31])
    d3 = delay_operator(fr, tt[12])
    d1p = delay_operator(fr, tt[32])
    d2p = delay_operator(fr, tt[13])
    d3p = delay_operator(fr, tt[21])
    # Build TDI transformation matrix
    zer = np.zeros(fr.size)
    # TDI 1.5
    if not tdi2:
        tdi_mat = np.array(
            [[d2*d2p-1, zer, -(d3*d3p-1) * d2p, -(d3*d3p-1), (d2*d2p-1) * d3, zer],
                [-(d1*d1p-1) * d3p, d3*d3p-1, zer, zer, -(d1*d1p-1), (d3*d3p-1) * d1],
                [zer, - d1p * (d2*d2p-1), d1*d1p-1, (d1*d1p-1) * d2, zer, -(d2*d2p-1)]])
    # TDI 2.0
    else:
        tdi_mat = np.array(
                [[d2p*d2 + d2p*d2*d3*d3p - 1 - d3*d3p*d2p*d2*d2p*d2,
                zer,
                d2p + d2p*d2*d3*d3p*d3*d3p*d2p - d3*d3p*d2p - d3*d3p*d2p*d2*d2p,
                1 + d2p*d2*d3*d3p*d3*d3p - d3*d3p - d3*d3p*d2p*d2,
                d2p*d2*d3 + d2p*d2*d3*d3p*d3 - d3 - d3*d3p*d2p*d2*d2p*d2*d3,
                zer],
                [d3p + d3p*d3*d1*d1p*d1*d1p*d3p - d1*d1p*d3p - d1*d1p*d3p*d3*d3p,
                d3p*d3 + d3p*d3*d1*d1p - 1 - d1*d1p*d3p*d3*d3p*d3,
                zer,
                zer,
                1 + d3p*d3*d1*d1p*d1*d1p - d1*d1p - d1*d1p*d3*d3p,
                d3p*d3*d1 + d3p*d3*d1*d1p*d1 - d1 - d1*d1p*d3p*d3*d3p*d3*d1],
                [zer,
                d1p + d1p*d1*d2*d2p*d2*d2p*d1p - d2*d2p*d1p - d2*d2p*d1p*d1*d1p,
                d1p*d1 + d1p*d1*d2*d2p - 1 - d2*d2p*d1p*d1*d1p*d1,
                d1p*d1*d2 + d1p*d1*d2*d2p*d2 - d2 - d2*d2p*d1p*d1*d1p*d1*d2,
                zer,
                1 + d1p*d1*d2*d2p*d2*d2p - d2*d2p - d2*d2p*d1*d1p]])

    tdi_mat_eval = np.swapaxes(tdi_mat.T, 1, 2)

    if tdi_var == 'AET':
        tdi_mat_eval = dot(aet_mat[np.newaxis, :, :], tdi_mat_eval)

    return tdi_mat_eval


def eta_matrix(fr, tt, noise_type):
    """

    Build matrix that convert single-link measurements to the eta variables,
    to cancel spacraft motion and primed lasers.
    Delay symbols should be ordered as 1, 2, 3, 1p, 2p, 2p


    Parameters
    ----------
    Parameters
    ----------
    fr : ndarray
        frequency array [Hz]
    tt : dict
        dictionary of arm light travel times [s]
    noise_type : string
        Noise source among {"SCI", "TM", "REF"}

    Returns
    -------
    a_eta : ndarray
        transform matrix to eta variables

    """

    d23 = delay_operator(fr, tt[23])
    d31 = delay_operator(fr, tt[31])
    d12 = delay_operator(fr, tt[12])
    d32 = delay_operator(fr, tt[32])
    d13 = delay_operator(fr, tt[13])
    d21 = delay_operator(fr, tt[21])

    zer = np.zeros(fr.size)
    one = np.ones(fr.size)

    if (noise_type == 'TM') | (noise_type == 'REF'):
        # Transformation to xi variables
        a_xi_eps = 0.5 * np.array([[one, zer, zer, zer, d12, zer],
                                   [zer, one, zer, zer, zer, d23],
                                   [zer, zer, one, d31, zer, zer],
                                   [zer, zer, d13, one, zer, zer],
                                   [d21, zer, zer, zer, one, zer],
                                   [zer, d32, zer, zer, zer, one]], dtype=complex)
    if noise_type == 'REF':
        # TM interferometer contribution to eta
        a_eta_tau = 0.5 * np.array([[zer, -d12, zer, zer, d12, zer],
                                    [zer, zer, -d23, zer, zer, d23],
                                    [-d31, zer, zer, d31, zer, zer],
                                    [one, zer, zer, -one, zer, zer],
                                    [zer, one, zer, zer, -one, zer],
                                    [zer, zer, one, zer, zer, -one]], dtype=complex)

    if noise_type == 'SCI':
        mat = np.array([[one, zer, zer, zer, zer, zer],
                        [zer, one, zer, zer, zer, zer],
                        [zer, zer, one, zer, zer, zer],
                        [zer, zer, zer, one, zer, zer],
                        [zer, zer, zer, zer, one, zer],
                        [zer, zer, zer, zer, zer, one]], dtype=complex)
    elif noise_type == 'TM':
        mat = - a_xi_eps
    elif noise_type == 'REF':
        mat = a_eta_tau + a_xi_eps
    # Transform matrices from size 6 x 6 x nf to nf x 6 x 6
    return np.swapaxes(mat.T, 1, 2)


if __name__ == '__main__':
    from timeseries import TimeSeries
    import doctest
    doctest.testmod()
