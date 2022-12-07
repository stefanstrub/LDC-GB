""" Provide TDI facility using xarray dataset
"""

import numpy as np
import xarray as xr
import h5py
import ldc.io.hdf5 as h5io
from ldc.common.series import TimeSeries, FrequencySeries

def XYZ2AET(X, Y, Z):
    return ((Z - X)/np.sqrt(2.0),
            (X - 2.0*Y + Z)/np.sqrt(6.0),
            (X + Y + Z)/np.sqrt(3.0))


class TDI(xr.Dataset):
    __slots__ = ()
    def XYZ2AET(self): 
        """ Convert TDI XYZ into AET

        >>> XYZ = dict([(n,TimeSeries(np.random.randn(50), name=n)) for n in ["X", "Y", "Z"]])
        >>> tdi = TDI(XYZ)
        >>> tdi.XYZ2AET()
        """
        A,E,T = XYZ2AET(self["X"], self["Y"], self['Z'])
        A.attrs, E.attrs, T.attrs  = [self["X"].attrs]*3
        super().__init__(dict([(k,v) for k,v in zip(["A", "E", "T"], [A, E, T])]))

    def AET2XYZ(self):
        """ Convert TDI AET into XYZ
        """
        [X,Y,Z] = [-1/np.sqrt(2.0)*self["A"]+   1/np.sqrt(6.0)*self["E"]+    1/np.sqrt(3.0)*self["T"],
                   -np.sqrt(2.0/3.0)*self["E"]+    1/np.sqrt(3.0)*self["T"],
                   1/np.sqrt(2.0)*self["A"]+   1/np.sqrt(6.0)*self["E"]+    1/np.sqrt(3.0)*self["T"]]
        X.attrs, Y.attrs, Z.attrs  = [self["A"].attrs]*3
        super().__init__(dict([(k,v) for k,v in zip(["X", "Y", "Z"], [X, Y, Z])]))

    @staticmethod
    def load(filename, name='data', compound=True):
        """ Load tdi from file

        >>> XYZ = dict([(n,TimeSeries(np.random.randn(50), name=n)) for n in ["X", "Y", "Z"]])
        >>> tdi1 = TDI(XYZ)
        >>> tdi1.save("test.h5")
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
                ds = dict([(k,FrequencySeries(arr[k], fmin=arr['f'][0], df=arr['f'][1]-arr['f'][0], units=attrs["units"]))
                      for k in arr.dtype.names if k!='f'])
        else:
            ds = xr.open_dataset(filename, engine='h5netcdf', group=name)
        return TDI(ds)

    def save(self, filename, mode='a', name='data', compound=True):
        """ Save tdi to file

        >>> XYZ = dict([(n,TimeSeries(np.random.randn(50), name=n)) for n in ["X", "Y", "Z"]])
        >>> tdi = TDI(XYZ)
        >>> tdi.save("test2.h5")
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
            h5io.save_array(filename, arr, name=name, **attr)


if __name__ == '__main__':
    from timeseries import TimeSeries
    import doctest
    doctest.testmod()
