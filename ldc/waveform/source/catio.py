""" Load catalog from file in historical formats
"""

import h5py as h5
import numpy as np
import numpy.lib.recfunctions as recf
from ldc.common import constants, tools

YRSID_SI = constants.Nature.SIDEREALYEAR_J2000DAY*24*60*60

def load_gb_catalog(catalog, rename=False):
    """Load GB catalog either in h5 or npy format.
    """
    if catalog.split(".")[-1] in ["h5", "hdf5"]:
        return load_h5_catalog(catalog)
    cat = np.load(catalog)
    units = dict({})
    new_names = dict({})
    for k in cat.dtype.names:
        if '[' in k:
            name, u = k.split("[")
            units[name] = u[:-1]
            new_names[k] = name
        else:
            new_names[k] = k
    if rename:
        cat = recf.rename_fields(cat, new_names)
    return cat, units

def load_h5_catalog(catalog):
    """ Load h5 catalog a la LISAhdf5.
    """
    h5file = h5.File(catalog, mode='r')
    gp = h5file.get("H5LISA/GWSources")
    params = dict()
    sources = list(gp.keys())
    params = list(gp.get(sources[0]).keys())
    D = dict().fromkeys(["Name"]+params)
    for k, v in D.items():
        D[k] = []
    for s in sources:
        src = gp.get(s)
        D["Name"].append(s)
        for ky in list(src.keys()):
            D[str(ky)].append(src.get(ky).value)

    units = dict().fromkeys(["Name"]+params)
    for ky in list(src.keys()):
        p = src.get(ky)
        if 'Units' in p.attrs:
            units[str(ky)] = p.attrs["Units"]
    h5file.close()
    units = to_astropy_units_name(units)
    cat = np.rec.fromarrays(list(D.values()), names=list(D.keys()))
    return cat, units

def to_astropy_units_name(units):
    for k,v in units.items():
        if v=="Radian":
            units[k] = 'rad'
        elif v=='strain':
            units[k] = '1'
    return units

def load_mbhb_catalog(catalog, redshifted_mass=False, non_precessing=False):
    """ Load .dat
    """
    dat = np.genfromtxt(catalog, names=True, dtype=([float]*19 + ['|U3'] + [float]*8))
    if dat.size==1:
        dat = np.array([dat])
    beta = 0.5*np.pi - dat['ecl_colat']
    lam = dat['ecl_long']
    tc = dat['Tc_yrs']*YRSID_SI
    phi0 = dat['phi0']
    mass1 = dat['m1']
    mass2 = dat['m2']
    th1 = dat['thS1']
    th2 = dat['thS2']
    spin1 = dat['a1']
    spin2 = dat['a2']
    psi, incl = tools.aziPolAngleL2PsiIncl(beta, lam, dat['thL'],  dat['phL'])
    psi[psi<0] = psi[psi<0] + 2.0*np.pi
    chi1 = dat['a1']*np.cos(dat['thS1'])
    chi2 = dat['a2']*np.cos(dat['thS2'])

    if redshifted_mass:
        mass1 *= (1. + dat['z'])
        mass2 *= (1. + dat['z'])

    sindex = mass1 < mass2 # switch index
    mass1[sindex], mass2[sindex] = mass2[sindex], mass1[sindex]
    phi0[sindex] += np.pi
    chi1[sindex], chi2[sindex] = chi2[sindex], chi1[sindex]
    th1[sindex], th2[sindex] = th2[sindex], th1[sindex]
    spin1[sindex], spin2[sindex] = spin2[sindex], spin1[sindex]

    if non_precessing:
        cat = np.rec.fromarrays([beta, lam,  chi1, chi2,
                                 mass1, mass2, tc, phi0, psi, incl, dat['DL_Mpc']],
                                names=['EclipticLatitude', 'EclipticLongitude',
                                       'Spin1', 'Spin2','Mass1', 'Mass2',
                                       'CoalescenceTime', 'PhaseAtCoalescence',
                                       'Polarization', 'Inclination','Distance'])
    else:
        cat = np.rec.fromarrays([beta, lam,
                                 th1, th2, spin1, spin2,
                                 mass1, mass2, tc, phi0, dat['thL'], dat['phL'],
                                 dat['z'], dat['DL_Mpc']],
                                names=['EclipticLatitude', 'EclipticLongitude',
                                       'PolarAngleOfSpin1', 'PolarAngleOfSpin2',
                                       'Spin1', 'Spin2', 'Mass1', 'Mass2',
                                       'CoalescenceTime', 'PhaseAtCoalescence',
                                       'InitialPolarAngleL', 'InitialAzimuthalAngleL',
                                       'Redshift', 'Distance'])

    units = {'EclipticLatitude': 'rad', 'EclipticLongitude': 'rad',
             'PolarAngleOfSpin1': 'rad', 'PolarAngleOfSpin2': 'rad',
             'Spin1': '1', 'Spin2': '1', 'Mass1': 'Msun', 'Mass2': 'Msun',
             'CoalescenceTime': 's',
             'PhaseAtCoalescence': 'rad',
             'InitialPolarAngleL': 'rad',
             'InitialAzimuthalAngleL': 'rad',
             'Polarization': 'rad',
             'Inclination': 'rad',
             'Redshift': '1',
             'Distance': 'Mpc' }

    return cat, units
