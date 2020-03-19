import numpy as np
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u

#Note that Planck15 is H0 = 67.74, Om0=0.3075
ldc_cosmo = FlatLambdaCDM(H0=67.1, Om0=0.3175)

def DL(zup):
    quantity = ldc_cosmo.luminosity_distance(zup)
    return quantity.value, quantity.unit

def check_cosmo(dl, z):
    """ dl in Mpc
    """
    expected_dl = ldc_cosmo.luminosity_distance(z)* u.Mpc
    dl = dl * u.Mpc
    if np.isscalar(dl.value):
        return np.isclose(expected_dl.value, dl.value)
    else:
        return np.isclose(expected_dl.value, dl.value).sum()==len(dl)
