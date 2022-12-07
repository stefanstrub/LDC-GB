""" Select sources to build catalogs
"""
from abc import ABC, abstractmethod
import logging
import re
import sys
import numpy as np
import numpy.lib.recfunctions as recf

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import BarycentricMeanEcliptic, Galactocentric

import lisaconstants as constants

from ldc.common import tools
from ldc.waveform.waveform import BHB_IMRPhenomD, GB_fdot
from ldc.lisa.noise import simple_snr
from .catio import *

#pylint:disable=W1201
#pylint:disable=C0103

YRSID_SI = constants.SIDEREALYEAR_J2000DAY*24*60*60
MTsun = constants.GM_SUN/constants.SPEED_OF_LIGHT**3
deg2rad = np.pi/180.

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
        if source_type == "SBBH" and approximant == "IMRPhenomD":
            return SOBBHMaker(source_type, approximant, **kwargs)
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


class MBHBMaker(SourceMaker, BHB_IMRPhenomD):

    def __init__(self, source_type, approximant, **kwargs):
        """ Initialize object.
        """
        SourceMaker.__init__(self, source_type, approximant, **kwargs)
        BHB_IMRPhenomD.__init__(self, "catalog", source_type, approximant)

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
                            mass_total=(0, 5000), redshifted_mass=True, non_precessing=False,
                            indices=None,
                            **kwargs):
        """Make a random selection of sources.

        If catalogs keyword is provided, sources are randomly chosen
        from them. Otherwise sources are randomly drawn.

        If nsource is negative, number of sources is randomly chosen
        from the lengths of the given ancillary catalogs.
        """
        # load catalogs
        if hasattr(self, 'catalogs'):
            C = [load_mbhb_catalog(cat, redshifted_mass=redshifted_mass,
                                   non_precessing=non_precessing) for cat in self.catalogs]
            units = C[0][1]
            C = [c[0] for c in C]
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

        coalescence_time = np.array(coalescence_time)*YRSID_SI
        C["CoalescenceTime"] /= 5. #5 years
        self.logger.info("CoalescenceTime is reset to 1 year before selection")
        
        if not indices:
            # apply selection criteria
            mratio = C["Mass1"] / C["Mass2"] # mass ratio cut
            selec = (mratio > mass_ratio[0]) & (mratio < mass_ratio[1])
            self.logger.info("mass ratio %d"%(selec.sum()))
            Ms = (C["Mass1"] + C["Mass2"]) * MTsun # total mass cut -> freq. cut
            selec &= (Ms > mass_total[0]) & (Ms < mass_total[1])
            self.logger.info("mass total %d"%(selec.sum()))

            for k, c in zip(["Spin1", "Spin2", "CoalescenceTime"], [spin1, spin2, coalescence_time]):
                selec &= (C[k] > c[0]) & (C[k] < c[1])
                self.logger.info("%s %d"%(k, selec.sum()))
            C = C[selec]

            if len(C) < nsource:
                self.logger.error("Number of sources in after selection is smaller than requested")
                raise ValueError

            if 'seed' in kwargs:
                np.random.seed(kwargs["seed"])
            ind = np.random.choice(len(C), nsource, replace=False)

        else:
            ind = np.array(indices, dtype=int)
            if np.any(ind>=len(C)):
                self.logger.error("Some indices greater than catalog size")
                raise ValueError

        C = C[ind]
        cadence, obs_duration = self.set_cadence(C)
        C = recf.append_fields(C, ['ObservationDuration', 'Cadence'],
                               [obs_duration, cadence],
                               usemask=False)
        return C, units


class SOBBHMaker(SourceMaker, BHB_IMRPhenomD):
    
    def __init__(self, source_type, approximant, **kwargs):
        """ Initialize object.
        """
        SourceMaker.__init__(self, source_type, approximant, **kwargs)
        BHB_IMRPhenomD.__init__(self, "catalog", source_type, approximant)
        
    def set_cadence(self, params, base_dt=5, base_Tobs=2 * YRSID_SI):
        """Adjust cadence to source frequency, to ensure a good sampling.

        TODO: this has been deactivated for now.
        """
        b = np.ones(params.shape)
        
        dt = base_dt * b
        Tobs = base_Tobs * b
        return dt, Tobs

    def choose_from_catalog(self, nsource, mass_ratio=(1, 10), spin1=(0.01, 0.99),
                            spin2=(0.01, 0.99), coalescence_time=(0.0001, 10),
                            mass_total=(0, 5000), redshifted_mass=True, non_precessing=False,
                            indices=None,
                            **kwargs):
        """Make a random selection of sources.

        If catalogs keyword is provided, sources are randomly chosen
        from them. Otherwise sources are randomly drawn.

        If nsource is negative, number of sources is randomly chosen
        from the lengths of the given ancillary catalogs.
        """
        # load catalogs
        if hasattr(self, 'catalogs'):
            C = [load_sobbh_catalog(cat, redshifted_mass=redshifted_mass,
                                    non_precessing=non_precessing) for cat in self.catalogs]
            units = C[0][1]
            C = [c[0] for c in C]
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

#         coalescence_time = np.array(coalescence_time)*YRSID_SI
#         C["CoalescenceTime"] /= 5. #5 years
#         self.logger.info("CoalescenceTime is reset to 1 year before selection")
 
#         for k, c in zip(["Spin1", "Spin2", "CoalescenceTime"], [spin1, spin2, coalescence_time]):
#             selec &= (C[k] > c[0]) & (C[k] < c[1])
#             self.logger.info("%s %d"%(k, selec.sum()))
        C = C[selec]
 
#         if len(C) < nsource:
#             self.logger.error("Number of sources in after selection is smaller than requested")
#             raise ValueError
 
        # random choice of source index
#         if 'seed' in kwargs:
#             np.random.seed(kwargs["seed"])
#         if not indices:
#             ind = np.random.choice(len(C), nsource, replace=False)
#             #ind = np.arange(len(C))
#         else:
#             ind = np.array(indices, dtype=int)
#             if np.any(ind>=len(C)):
#                 self.logger.error("Some indices greater than catalog size after selection")
#                 raise ValueError
#         C = C[ind]
        
        cadence, obs_duration = self.set_cadence(C)
        C = recf.append_fields(C, ['ObservationDuration', 'Cadence'],
                               [obs_duration, cadence],
                               usemask=False)
        return C, units


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
        C = [load_gb_catalog(cat, self.logger, rename=True) for cat in self.catalogs]
        units = C[0][1]
        C = np.hstack([c[0] for c in C])
        self.logger.info("Input catalog has %d sources"%len(C))
        if len(C) < nsource:
            self.logger.error("Number of sources in catalogs is smaller than requested")
        if 'seed' in kwargs:
            np.random.seed(kwargs["seed"])
        ind = np.random.choice(len(C), nsource, replace=False)
        return C[ind], units

    def get_frequency(self, Craw):
        """ Get frequency parameter from period. 
        """
        P_s = None
        if "Frequency[Hz]" in Craw.dtype.names:
            f = Craw["Frequency[Hz]"]
        elif "Period[days]" in Craw.dtype.names:
            P_s = Craw["Period[days]"]*86400. # nsource
            f = 2./P_s
        elif "Period[sec]" in Craw.dtype.names:
            P_s = Craw["Period[sec]"]
            f = 2./P_s
        else:
            self.logger.error("Frequency and Period not found in the parameters")
        return f, P_s

    def get_sky_position(self, Craw):
        """ Convert in galactic coordinates
        """
        if "GalacticLongitude[deg]" in Craw.dtype.names:
            sky_gal = SkyCoord(l=Craw["GalacticLongitude[deg]"]*u.degree,
                               b=Craw["GalacticLatitude[deg]"]*u.degree, frame='galactic')
        else:
            sky_gal = SkyCoord(x=np.array(Craw['GalacticX[pc]'])*u.pc,\
                     y=np.array(Craw['GalacticY[pc]'])*u.pc,\
                     z=np.array(Craw['GalacticZ[pc]'])*u.pc,\
                     frame=Galactocentric)
        sky_ecl = sky_gal.transform_to(BarycentricMeanEcliptic)
        b_ecl = sky_ecl.lat.rad
        l_ecl = sky_ecl.lon.rad
        return b_ecl, l_ecl, sky_ecl.distance
        
    def get_chirpmass(self, Craw, D_pc, f, **kwargs):
        """ Return chirp mass and total mass
        """
        if "Distance[kpc]" in Craw.dtype.names:
            D_pc = Craw["Distance[kpc]"]*1e3
        DL = D_pc*constants.PARSEC_METER / constants.SPEED_OF_LIGHT
        m1 = randomize_gaussian(Craw["Mass1[MSun]"], kwargs["random_mass"],
                                0.001, 1e3, self.logger)*MTsun
        m2 = randomize_gaussian(Craw["Mass2[MSun]"], kwargs["random_mass"],
                                0.001, 1e3, self.logger)*MTsun
        M = m1 + m2
        eta = m1*m2 / (M*M)
        Mc = M*eta**(3./5.)
        amplitude = 2*(M**(5./3.)*eta/DL)*(np.pi*f)**(2./3.) # Amplitude
        return amplitude, Mc
    
    def draw_random_catalog(self, nsource, **kwargs):
        """ Make a list of sources randomizing parameters in the catalogs.
        """
        Craw = [load_gb_catalog(cat, self.logger) for cat in self.catalogs]
        units = Craw[0][1]
        Craw = np.hstack([c[0] for c in Craw])
        self.logger.info("Input catalog has %d sources"%len(Craw))
        if len(Craw) < nsource:
            self.logger.error(f"Number of sources in catalogs is smaller than requested")

        f, P_s = self.get_frequency(Craw) # Frequency
        b_ecl, l_ecl, D_pc = self.get_sky_position(Craw)  # Sky position
        amplitude, Mc = self.get_chirpmass(Craw, D_pc, f, **kwargs) # mass
        dfdt = ((96./5.) * Mc**(5./3.) * np.pi**(8./3.) * f**(11./3.)) # Frequency derivative

        if "FrequencyDerivative[Hz/sec]" in Craw.dtype.names:
            dfdt = Craw["FrequencyDerivative[Hz/sec]"]
        elif "PeriodDerivative[sec/sec]" in Craw.dtype.names and P_s is not None:
            dPdt = Craw["PeriodDerivative[sec/sec]"]
            dfdt = -2.0* dPdt / (P_s*P_s)
        else:
            self.logger.error("Derivative of Frequency or Period not found in the parameters")
        
        selec = np.ones((len(Craw)), dtype=bool)
        if 'frequency' in kwargs:
            self.logger.info(f"Select on freq. range {kwargs['frequency']}")
            selec &= f>kwargs['frequency'][0]
            selec &= f<kwargs['frequency'][-1]
            self.logger.info(f"keep {selec.sum()} / {len(selec)} sources")

        if 'amplitude' in kwargs:
            self.logger.info(f"Select on amplitude range {kwargs['amplitude']}")
            selec &= amplitude>kwargs['amplitude'][0]
            selec &= amplitude<kwargs['amplitude'][-1]
            self.logger.info(f"keep {selec.sum()} / {len(selec)} sources")

        if 'snr' in kwargs:
            self.logger.info(f"Select on snr range {kwargs['snr']}")
            snr = simple_snr(f, amplitude)
            selec &= snr>kwargs['snr'][0]
            selec &= snr<kwargs['snr'][-1]
            self.logger.info(f"keep {selec.sum()} / {len(selec)} sources")
            
        # Apply selection
        Craw = Craw[selec]
        f, dfdt = f[selec], dfdt[selec]
        b_ecl, l_ecl = b_ecl[selec], l_ecl[selec]
        amplitude = amplitude[selec]

        # Choose source
        if 'seed' in kwargs:
            np.random.seed(kwargs["seed"])
        ind  = np.random.choice(len(Craw), nsource, replace=False)
        Craw = Craw[ind]
        f = randomize_gaussian(f[ind], kwargs['random_frequency'], 1.e-6, 1., self.logger)
        b_ecl, l_ecl = b_ecl[ind], l_ecl[ind]
        amplitude = amplitude[ind]
        dfdt = dfdt[ind]
        
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

        units = {'EclipticLatitude':         'rad',
                 'EclipticLongitude':        'rad',
                 'Amplitude':                '1',
                 'Frequency':                'Hz',
                 'FrequencyDerivative':      'Hz2',
                 'Inclination':              'rad',
                 'Polarization':             'rad',
                 'InitialPhase':             'rad'}

        lbase = [name, amplitude, b_ecl, l_ecl, f, dfdt, inc, phi0, psi]
        nbase = ['Name', 'Amplitude', 'EclipticLatitude',
                  'EclipticLongitude', 'Frequency',
                  'FrequencyDerivative', 'Inclination',
                  'InitialPhase', 'Polarization']
        if 'copy' in kwargs:
            for f in kwargs["copy"]:
                if f in Craw.dtype.names:
                    lbase.append(Craw[f])
                    name, unit = f.split("[")
                    nbase.append(name)
                    units[name] = unit[:-1]
        C = np.rec.fromarrays(lbase, names=nbase)
        return C, units
