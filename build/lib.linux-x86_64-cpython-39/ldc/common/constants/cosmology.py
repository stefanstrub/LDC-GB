"""
Cosmological constants provided by astropy
"""
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u

#pylint:disable=C0103

#Note that Planck15 is H0 = 67.74, Om0=0.3075
ldc_cosmo = FlatLambdaCDM(H0=67.1, Om0=0.3175)

def DL(zup):
    """ Return distance luminosity and associated units

    >>> print(ldc_cosmo)
    FlatLambdaCDM(H0=67.1 km / (Mpc s), Om0=0.318, Tcmb0=0 K, Neff=3.04, m_nu=None, Ob0=None)
    >>> DL(1)
    (6823.09048017983, Unit("Mpc"))
    """
    quantity = ldc_cosmo.luminosity_distance(zup)
    return quantity.value, quantity.unit

def check_cosmo(dl, z):
    """Check compatibility of distance luminosity and redshift with
    default cosmology.

    >>> check_cosmo(DL(1)[0], 1)
    True
    """
    expected_dl = ldc_cosmo.luminosity_distance(z)* u.Mpc
    dl = dl * u.Mpc
    if np.isscalar(dl.value):
        return np.isclose(expected_dl.value, dl.value)
    else:
        return np.isclose(expected_dl.value, dl.value).sum()==len(dl)
