""" Select sources to build catalogs
"""
from abc import ABC, abstractmethod
import logging
import re
import sys
import h5py as h5
import numpy as np
import numpy.lib.recfunctions as recf

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import BarycentricMeanEcliptic

from ldc.common import constants, tools
from ldc.waveform.waveform import BBH_IMRPhenomD, GB_fdot

#pylint:disable=W1201
#pylint:disable=C0103

YRSID_SI = constants.Nature.SIDEREALYEAR_J2000DAY*24*60*60
MTsun = constants.Nature.SUN_GM/constants.Nature.VELOCITYOFLIGHT_CONSTANT_VACUUM**3
deg2rad = np.pi/180.

def load_gb_catalog(catalog):
    """Load GB catalog either in h5 or npy format.
    """
    if catalog.split(".")[-1] == "h5":
        return load_h5_catalog(catalog)
    elif catalog.split(".")[-1] in ["dat, data"]:
        return load_txt_catalog(catalog)
    cat = np.load(catalog)
    return cat

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
    h5file.close()
    cat = np.rec.fromarrays(list(D.values()), names=list(D.keys()))
    return cat

def load_mbhb_catalog(catalog):
    """ Load .dat
    """
    dat = np.genfromtxt(catalog, names=True, dtype=([float]*19 + ['|U3'] + [float]*8))

    mass1 = dat['m1']#*(1. + dat['z'])
    mass2 = dat['m2']#*(1. + dat['z'])
    beta = 0.5*np.pi - dat['ecl_colat']
    lam = dat['ecl_long']
    tc = dat['Tc_yrs']*YRSID_SI
    th1 = dat['thS1']
    th2 = dat['thS2']
    phi0 = dat['phi0']
    spin1 = dat['a1']
    spin2 = dat['a2']

    sindex = mass1 < mass2 # switch index
    mass1[sindex], mass2[sindex] = mass2[sindex], mass1[sindex]
    th1[sindex], th2[sindex] = th2[sindex], th1[sindex]
    spin1[sindex], spin2[sindex] = spin2[sindex], spin1[sindex]
    phi0[sindex] += np.pi

    cat = np.rec.fromarrays([beta, lam,
                             dat['thS1'], dat['thS2'], dat['a1'], dat['a2'],
                             mass1, mass2, tc, phi0, dat['thL'], dat['phL'],
                             dat['z'], dat['DL_Mpc']*1e-3],
                            names=['EclipticLatitude', 'EclipticLongitude',
                                   'PolarAngleOfSpin1', 'PolarAngleOfSpin2',
                                   'Spin1', 'Spin2', 'Mass1', 'Mass2',
                                   'CoalescenceTime', 'PhaseAtCoalescence',
                                   'InitialPolarAngleL', 'InitialAzimuthalAngleL',
                                   'Redshift', 'Distance'])
    return cat

def randomize_gaussian(x, randx, xmin, xmax, logger):
    """Add or multiply x by a random quantity, such that it remains
    between xmin and xmax.

    randx can be in ['None', 'gaussian_%fpercent', 'gaussian_%f']
    """
    # analyse randx
    if randx in ['None', '0']:
        return x
    pmul = re.compile(r"gaussian_[-+]?\d*\.\d+|\d+percent")
    padd = re.compile(r"gaussian_[-+]?\d*\.\d+|\d+")
    if pmul.match(randx):
        fmul = float(re.findall(r"\d*\.\d+|\d+", randx)[0])
        fadd = 0
    elif padd.match(randx):
        fadd = float(re.findall(r"\d*\.\d+|\d+", randx)[0])
        fmul = 0
    else:
        logger.error("randomization string unrecognized: %s"%randx)
    logger.info("will use fadd=%f, fmul=%f"%(fadd, fmul))

    i = 0 ; xo = -1e-30
    while not (np.any(xmin < xo) and np.any(xo < xmax)) and i < 1000:
        r = np.random.randn(x.size)
        xo = x * (1 + fmul*r) + fadd*r # 1 of the 2 term is 0
        i += 1
    if i == 1000:
        logger.error("More than 1000 tries in randomization")
    return xo


class SourceMaker(ABC):
    """ Generate source catalogs.

    """
    def __init__(self, source_type, approximant, catalogs=None,
                 logger=None, verbose=True):
        """Initialize source.
        """
        self.source_type = source_type
        self.approximant = approximant
        if catalogs is not None:
            self.catalogs = catalogs
            if isinstance(self.catalogs, str):
                self.catalogs = [self.catalogs]
        if logger is None:
            self.set_logger(verbose=verbose)
        else:
            self.logger = logger
        self.logger.info("Source type is %s"%self.__class__)

    @classmethod
    def type(cls, source_type, approximant, **kwargs):
        """ Return instance corresponding to source type.
        """
        if source_type == "MBHB" and approximant == "IMRPhenomD":
            return MBHBMaker(source_type, approximant, **kwargs)
        elif source_type == "GB":
            return GBMaker(source_type, approximant, **kwargs)
        else:
            raise ValueError("Invalid source_type %s"%source_type)

    @abstractmethod
    def choose_from_catalog(self, nsource, **kwargs):
        """Generate a catalog of nsource sources randomly chosen from
        catalogs.
        """
        pass

    def set_logger(self, verbose=True):
        """ Set a stream logger.
        """
        logger = logging.getLogger()
        if verbose:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.ERROR)
        shandler = logging.StreamHandler(sys.stdout)
        logger.addHandler(shandler)
        self.logger = logger

    def close_logger(self):
        """ Close the logger
        """
        for h in self.logger.handlers:
            self.logger.removeHandler(h)
            h.flush()
            h.close()

class MBHBMaker(SourceMaker, BBH_IMRPhenomD):

    def __init__(self, source_type, approximant, **kwargs):
        """ Initialize object.
        """
        SourceMaker.__init__(self, source_type, approximant, **kwargs)
        BBH_IMRPhenomD.__init__(self, "catalog", source_type, approximant)

    def set_cadence(self, params, base_dt=3, base_Tobs=YRSID_SI):
        """Adjust cadence to source frequency, to ensure a good sampling.

        TODO: this has been deactivated for now.
        """
        Mc = tools.mchirpofm1m2(params["Mass1"], params["Mass2"])
        Ms = (params["Mass1"] + params["Mass2"]) * MTsun
        MfCUT_PhenomD = 0.2-1e-7
        fmax = MfCUT_PhenomD/Ms
        df = 1/(0.5/fmax)
        dt = 1/df
        dt = dt.round(decimals=1)

        dt = np.ones((len(params))) * base_dt
        Tobs = np.ones((len(params))) * base_Tobs
        # dt[dt>base_dt] = base_dt
        # Tobs = base_Tobs/ (base_dt/dt)
        # Tobs = Tobs.round()
        # Tobs[Tobs>base_Tobs] = base_Tobs
        return dt, Tobs

    def draw_random_catalog(self, n=1):
        """Build a completely random catalog.

        TODO: give the possibility to tune the interval for each
        parameter.
        """
        names = list(self.info().keys())
        d = np.rec.fromarrays([np.ones((n))*np.nan]*len(names), names=names)
        d['EclipticLatitude'] = 0.5*np.pi -\
                                np.arccos(np.random.uniform(-1.0, 1.0, size=n))
        d['EclipticLongitude'] = np.random.uniform(0.0, 2.0*np.pi, size=n)
        d['PolarAngleOfSpin1'] = np.arccos(np.random.uniform(-1.0, 1.0, size=n))
        d['PolarAngleOfSpin2'] = np.arccos(np.random.uniform(-1.0, 1.0, size=n))
        d['AzimuthalAngleOfSpin1'] = np.random.uniform(0.0, 2.0*np.pi, size=n)
        d['AzimuthalAngleOfSpin2'] = np.random.uniform(0.0, 2.0*np.pi, size=n)
        d['Spin1'] = np.random.uniform(0.0, 1.0, size=n)
        d['Spin2'] = np.random.uniform(0.0, 1.0, size=n)
        d["Mass1"] = 10.0**np.random.uniform(5.0, 7.0, size=n)
        d["Mass2"] = 10.0**np.random.uniform(5.0, 7.0, size=n)
        d["CoalescenceTime"] = np.random.uniform(0.5, 1.0, size=n)
        d['PhaseAtCoalescence'] = np.random.uniform(0.0, 2.0*np.pi, size=n)
        d['InitialPolarAngleL'] = np.arccos(np.random.uniform(-1.0, 1.0, size=n))
        d['InitialAzimuthalAngleL'] = np.random.uniform(0.0, 2.0*np.pi, size=n)

    def choose_from_catalog(self, nsource, mass_ratio=(1, 10), spin1=(0.01, 0.99),
                            spin2=(0.01, 0.99), coalescence_time=(0.0001, 10),
                            mass_total=(2, 8), **kwargs):
        """Make a random selection of sources.

        If catalogs keyword is provided, sources are randomly chosen
        from them. Otherwise sources are randomly drawn.

        If nsource is negative, number of sources is randomly chosen
        from the lengths of the given ancillary catalogs.
        """

        # load catalogs
        if hasattr(self, 'catalogs'):
            C = [load_mbhb_catalog(cat) for cat in self.catalogs]
            if nsource < 0:
                nsource = np.random.choice([c.size for c in C])
                self.logger.info("number of source is %d"%(nsource))
            C = np.hstack(C)
            self.logger.info("Input catalog has %d sources"%len(C))
            if len(C) < nsource:
                self.logger.error("Number of sources in catalogs (%d) is smaller than requested (%d)"%(len(C), nsource))
                raise ValueError
        else:
            raise AttributeError("Missing catalogs attributes")

        # apply selection criteria
        mratio = C["Mass1"] / C["Mass2"] # mass ratio cut
        selec = (mratio > mass_ratio[0]) & (mratio < mass_ratio[1])
        self.logger.info("mass ratio %d"%(selec.sum()))
        Ms = (C["Mass1"] + C["Mass2"]) * MTsun # total mass cut -> freq. cut
        selec &= (Ms > mass_total[0]) & (Ms < mass_total[1])
        self.logger.info("mass total %d"%(selec.sum()))

        coalescence_time = np.array(coalescence_time)*YRSID_SI
        C["CoalescenceTime"] /= 5. #5 years
        self.logger.info("CoalescenceTime is reset to 1 year before selection")

        for k, c in zip(["Spin1", "Spin2", "CoalescenceTime"], [spin1, spin2, coalescence_time]):
            selec &= (C[k] > c[0]) & (C[k] < c[1])
            self.logger.info("%s %d"%(k, selec.sum()))
        C = C[selec]

        if len(C) < nsource:
            self.logger.error("Number of sources in after selection is smaller than requested")
            raise ValueError

        # random choice of source index
        if 'seed' in kwargs:
            np.random.seed(kwargs["seed"])
        ind = np.random.choice(len(C), nsource, replace=False)
        cadence, obs_duration = self.set_cadence(C[ind])
        C = recf.append_fields(C[ind], ['ObservationDuration', 'Cadence'],
                               [obs_duration, cadence],
                               usemask=False)
        return C


class GBMaker(SourceMaker, GB_fdot):
    """The parameters of an output source is: Name[], Amplitude[], EclipticLatitude[rad],
    EclipticLongitude[rad], Frequency[Hz], FrequencyDerivative[Hz/s], Inclination[rad],
    InitialPhase[rad], Polarization[rad]
    """

    def __init__(self, source_type, approximant, **kwargs):
        """ Initialize object.
        """
        SourceMaker.__init__(self, source_type, approximant, **kwargs)
        GB_fdot.__init__(self, "catalog", source_type, approximant)

    def choose_from_catalog(self, nsource, **kwargs):
        """ Make a random selection of sources.
        """
        C = np.hstack([load_gb_catalog(cat) for cat in self.catalogs])
        self.logger.info("Input catalog has %d sources"%len(C))
        if len(C) < nsource:
            self.logger.error("Number of sources in catalogs is smaller than requested")
        if 'seed' in kwargs:
            np.random.seed(kwargs["seed"])
        ind = np.random.choice(len(C), nsource, replace=False)
        return C[ind]

    def draw_random_catalog(self, nsource, **kwargs):
        """ Make a list of sources randomizing parameters in the catalogs.
        """
        Craw = np.hstack([load_gb_catalog(cat) for cat in self.catalogs])
        self.logger.info("Input catalog has %d sources"%len(Craw))
        if len(Craw) < nsource:
            self.logger.error("Number of sources in catalogs is smaller than requested")
        if 'seed' in kwargs:
            np.random.seed(kwargs["seed"])

        ind = np.random.choice(len(Craw), nsource, replace=False)
        Craw = Craw[ind]

        # Frequency
        if "Period[days]" in Craw.dtype.names:
            P_s = Craw["Period[days]"]*86400. # nsource
        elif "Period[sec]" in Craw.dtype.names:
            P_s = Craw["Period[sec]"]
        else:
            self.logger.error("Period not found in the parameters")
        f = 2./P_s
        f = randomize_gaussian(f, kwargs['random_frequency'], 1.e-6, 1., self.logger)

        # Amplitude
        D_pc = Craw["Distance[kpc]"]*1e3
        DL = D_pc*constants.Nature.PARSEC_METER / \
              constants.Nature.VELOCITYOFLIGHT_CONSTANT_VACUUM
        m1 = randomize_gaussian(Craw["Mass1[MSun]"], kwargs["random_mass"],
                                0.001, 1e3, self.logger)*MTsun
        m2 = randomize_gaussian(Craw["Mass2[MSun]"], kwargs["random_mass"],
                                0.001, 1e3, self.logger)*MTsun
        M = m1 + m2
        eta = m1*m2 / (M*M)
        Mc = M*eta**(3./5.)
        amplitude = 2*(M**(5./3.)*eta/DL)*(np.pi*f)**(2./3.)

        # Frequency derivative
        dfdt = ((96./5.) * Mc**(5./3.) * np.pi**(8./3.) * f**(11./3.))
        if "PeriodDerivative[sec/sec]" in Craw.dtype.names:
            dPdt = Craw["PeriodDerivative[sec/sec]"]
            dfdt = -2.0* dPdt / (P_s*P_s)

        # Sky position
        # sky_gal = [ephem.Galactic(lon*deg2rad, lat*deg2rad, epoch='2000')
        #            for lon,lat in zip(Craw["GalacticLongitude[deg]"],
        #                               Craw["GalacticLatitude[deg]"])]
        # sky_ecl = [ephem.Ecliptic(s) for s in sky_gal]
        # b_ecl = np.array([float(s.lat) for s in sky_ecl]) # in radians
        # l_ecl = np.array([float(s.lon) for s in sky_ecl]) # in radians
        sky_gal = SkyCoord(l=Craw["GalacticLongitude[deg]"]*u.degree,
                           b=Craw["GalacticLatitude[deg]"]*u.degree, frame='galactic')
        sky_ecl = sky_gal.transform_to(BarycentricMeanEcliptic)
        b_ecl = sky_ecl.lat.rad
        l_ecl = sky_ecl.lon.rad

        # Inclination
        if "Inclination[rad]" in Craw.dtype.names:
            inc = Craw["Inclination[rad]"]
        elif kwargs['random_inclination'] != "uniform":
            self.logger.error("Inclination not specified: it should be either in the catalog or set to uniform distribution")
        if kwargs['random_inclination'] == "uniform":
            inc = np.arccos(np.random.uniform(-1., 1., size=nsource))
        else:
            inc = randomize_gaussian(inc, kwargs['random_inclination'],
                                     0.0, np.pi, self.logger)

        # Polarisation and initial phase
        psi = np.random.uniform(0, 2.*np.pi, size=nsource)
        phi0 = np.random.uniform(0, 2.*np.pi, size=nsource)
        name = Craw["ID[]"] if "ID[]" in Craw.dtype.names  else [str(i) for i in np.arange(nsource)]


        C = np.rec.fromarrays([name, amplitude, b_ecl,
                               l_ecl, f, dfdt, inc, phi0, psi],
                              names=['Name', 'Amplitude', 'EclipticLatitude',
                                     'EclipticLongitude', 'Frequency',
                                     'FrequencyDerivative', 'Inclination',
                                     'InitialPhase', 'Polarization'])
        return C
