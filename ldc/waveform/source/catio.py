""" Load catalog from file in historical formats
"""

import h5py as h5
import numpy as np
import pandas as pd
import numpy.lib.recfunctions as recf

import astropy.units as u
from astropy.cosmology import z_at_value

import lisaconstants as constants

from ldc.common import  tools
from ldc.common.constants.cosmology import ldc_cosmo
from ldc.utils.logging import init_logger
import ldc.io.hdf5 as h5io

YRSID_SI = constants.SIDEREALYEAR_J2000DAY*24*60*60

def load_gb_catalog(catalog, logger=None, rename=False):
    """Load GB catalog either in h5 or npy format.
    """
    if logger is None:
        logger = init_logger()
    if catalog.split(".")[-1] in ["h5", "hdf5"]:
        return load_h5_catalog(catalog,logger)
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

def load_h5_catalog(catalog, logger):
    """ Load h5 catalog a la LISAhdf5.
    """
    h5file = h5.File(catalog, mode='r')
    if 'H5LISA' in h5file.keys():
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
                D[str(ky)].append(src.get(ky)[()])
                #D[str(ky)].append(float(np.array(src.get(ky))))

        units = dict().fromkeys(["Name"]+params)
        for ky in list(src.keys()):
            p = src.get(ky)
            if 'Units' in p.attrs:
                units[str(ky)] = p.attrs["Units"]
        h5file.close()
        units = to_astropy_units_name(units)
        cat = np.rec.fromarrays(list(D.values()), names=list(D.keys()))
        return cat, units
    else:
        try:
            # Try to read it as a panda hdf5
            Craw = pd.read_hdf(catalog)
            logger.info(f"Load hdf5 H5panda catalog: {catalog} ...")
            units = dict({})
            for k in Craw.columns:
                if '[' in k:
                    name, u = k.split("[")
                    units[name] = u[:-1]
            cat = np.rec.fromarrays(np.array(Craw).T,names=list(Craw.columns))
            return cat, units
        except:
            return h5io.load_array(catalog)

def to_astropy_units_name(units):
    for k,v in units.items():
        if v=="Radian":
            units[k] = 'rad'
        elif v=='strain':
            units[k] = '1'
    return units

def load_sobbh_catalog(catalog, redshifted_mass=False, non_precessing=False):
    """ Load .dat
    
    .. note :: masses in the catalogue are already redshifted (SB)
    
    """
    dtype = [
        ('Mass1', '<f8'), ('Mass2', '<f8'), 
        ('Chi1', '<f8'), ('Chi2', '<f8'), ('Dist', '<f8'), 
        ('Inclination', '<f8'), ('RA', '<f8'), ('Dec', '<f8'),
        ('Phi', '<f8'), ('Lambda', '<f8'), ('Beta', '<f8'), 
        ('Psi', '<f8'), ('fstart', '<f8'), ('tc', '<f8'),
    ]
    dat = np.genfromtxt(catalog, skip_header=1, dtype=dtype)
    if dat.size==1:
        dat = np.array([dat])
        
    # extract quantities from catalog
    beta = dat['Beta']
    lam = dat['Lambda']
    # tc = dat['tc']*YRSID_SI
    phi0 = dat['Phi']
    mass1 = dat['Mass1']
    mass2 = dat['Mass2']
    
    dist = dat['Dist']
    z = np.array([ 
        z_at_value(ldc_cosmo.luminosity_distance, d * u.Mpc, zmin=0) for d in dist 
    ])
    
    chi1 = dat['Chi1']
    chi2 = dat['Chi2']
    incl = dat['Inclination']
    psi = dat['Psi']
    f0 = dat['fstart']
    
    theta1 = theta2 = 0
    spin1 = chi1 * np.cos(theta1)
    spin2 = chi2 * np.cos(theta2)
    
    # filter with m1 < m2 condition and inverse
    sindex = mass1 < mass2 # switch index
    mass1[sindex], mass2[sindex] = mass2[sindex], mass1[sindex]
    phi0[sindex] += np.pi
    chi1[sindex], chi2[sindex] = chi2[sindex], chi1[sindex]
    spin1[sindex], spin2[sindex] = spin2[sindex], spin1[sindex]
    
    if non_precessing:
        raise NotImplementedError("Non precessing not implemented")
    else:
        params = [
            beta,
            lam,
            mass1,
            mass2,
            spin1,
            spin2,
            incl,
            f0,
            phi0,
            psi,
            z,
            dist
        ]
        names = [
            'EclipticLatitude',
            'EclipticLongitude',
            'Mass1',
            'Mass2',
            'Spin1',
            'Spin2',
            'Inclination',
            'InitialFrequency',
            'InitialPhase',
            'Polarization',
            'Redshift',
            'Distance'
        ]
        cat = np.rec.fromarrays(
            params, names=names
        )
        
    units = {
        'EclipticLatitude': 'rad', 
        'EclipticLongitude': 'rad',
        'Mass1': 'Msun', 
        'Mass2': 'Msun',
        'Spin1': '1', 
        'Spin2': '1', 
        'Inclination': 'rad',
        'InitialFrequency': 'Hz',
        'InitialPhase': 'rad',
        'Polarization': 'rad',
        'Redshift': '1',
        'Distance': 'Mpc'
    }

    return cat, units    
    

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
