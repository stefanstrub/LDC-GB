import lisaconstants as constants
import numpy as np

from ldc.waveform.fastGB import get_default_orbits, get_buffersize
from ldc.common.series import TimeSeries, FrequencySeries

try:
    from fastgb import fastgb
except ImportError:
    print("Could not import fastgb. Please pip install fastgb")


class VFastGB(fastgb.FastGB):
    """ Vectorized (and publicly available) fastGB wrapper. 
    """

    def __init__(self, orbits=None, T=6.2914560e7, delta_t=15, t0=0, tinit=0, N=None):
        """ Init orbits and spacecraft positions.
        """
        if orbits is None or isinstance(orbits, dict):
            orbits = get_default_orbits()
        super().__init__(orbits=orbits, T=T, delta_t=delta_t, t0=t0, tinit=tinit, N=N)

    def __reduce__(self):
        return (self.__class__, (None, self.T, self.delta_t, self.t0, self.tinit, self.N))

        
    def get_fd_tdixyz(self, template=None, f0=None, fdot=None, ampl=None,
                      theta=None, beta=None, phi=None, psi=None, incl=None, phi0=None,
                      tdi2=False, xarray=True):
        """ Return TDI X,Y,Z in freq. domain.
        
        f0 in Hz, fdot in Hz/s, ampl in strain,
        theta,phi,psi,incl,phi0 in rad.
        """
        if template is not None:
            params = np.array([template[k] for k in ['Frequency', 'FrequencyDerivative',
                                                     'Amplitude', 'EclipticLatitude',
                                                     'EclipticLongitude', 'Polarization',
                                                     'Inclination', 'InitialPhase']]).T
        else:
            if beta is None:
                beta = 0.5*np.pi-theta
                phi0 = -phi0
            params = np.array([f0, fdot, ampl, beta, phi, psi, incl, phi0]).T
            
        params = params.reshape(-1, 8)
        X, Y, Z, kmin = super().get_fd_tdixyz(params, tdi2=tdi2)
        df = 1/self.T
        if xarray:
            X, Y, Z = (FrequencySeries(X[0,:], df=df, kmin=kmin[0], t0=self.t0, name="X"),
                       FrequencySeries(Y[0,:], df=df, kmin=kmin[0], t0=self.t0, name="Y"),
                       FrequencySeries(Z[0,:], df=df, kmin=kmin[0], t0=self.t0, name="Z"))
            return X, Y, Z
        return np.array([X.squeeze(), Y.squeeze(), Z.squeeze()]), int(kmin)

    def get_td_tdixyz(self, **kwargs):
        """  Return TDI X,Y,Z in time domain.
        """
        fX, fY, fZ = self.get_fd_tdixyz(**kwargs)
        return (fX.ts.ifft(dt=self.delta_t),
                fY.ts.ifft(dt=self.delta_t),
                fZ.ts.ifft(dt=self.delta_t))
