import numpy as np
import numpy.lib.recfunctions as rfn
import xarray as xr

from ldc.common.series import FrequencySeries
from ldc.lisa.orbits import AnalyticOrbits

from ldc.common import tools
from ldc.common.series import TDI

import lisaconstants

ASTRONOMICAL_YEAR = lisaconstants.SIDEREALYEAR_J2000DAY*24*60*60
AU_SI = lisaconstants.ASTRONOMICAL_UNIT

class FastBHB :
    def __init__(self, bbh_type, T=6.2914560e7, delta_t=5, approx="IMRPhenomD",
                    orbits=None, modes=None):
        """
        Parameters initilization of the black holes binaries class
            T : duration in second (will be converted in year for lisabeta)
            deltat : time step in second
            approx : approximant for the IMRPhenom waveform computation
            orbits : LISA's orbit type
        """

        #Type of BBH waveform
        if bbh_type.lower() in ["mbhb","sobhb", "sbbh"] :
            self.bbh_type = bbh_type
        else :
            raise NotImplementedError("Only MBHB and SOBHB are available.")

        #Wave duration and time step
        self.Tobs = T
        self.Tobs_yr = self.Tobs/ASTRONOMICAL_YEAR
        self.dt = delta_t

        #Orbits initialization
        if orbits is not None:
            if not isinstance(orbits, AnalyticOrbits):
                raise TypeError('FastBHB approximation requires analytic orbits')
            else:
                init_time = orbits.initial_position*ASTRONOMICAL_YEAR/(2 * np.pi)
                lisa_params = {
                    'OrbitOmega' : (2*np.pi) / ASTRONOMICAL_YEAR,
                    'OrbitPhi0' : orbits.initial_rotation,
                    'Orbitt0' : init_time,
                    'OrbitR' : AU_SI,
                    'OrbitL' : orbits.arm_length}
        else:
            # TODO: to be removed in the end
            import lisabeta.lisa.pyresponse as pyresponse
            lisa_params = pyresponse.LISAconstDict["Proposal"]

        #Waveform parameters used by lisabeta
        self.wvf_pars = {
            "minf": 1e-5,  "maxf": 0.1,
            "timetomerger_max": 1.0,
            "fend": None, "tmin": None, "tmax": self.Tobs_yr,
            "TDI": "TDIAET", "acc": 1e-4,
            "order_fresnel_stencil": 0,
            "approximant": approx,
            "modes": modes,
            "LISAconst": lisa_params,
            "responseapprox": "full",
            "frozenLISA": False,
            "TDIrescaled": False
        }

    @property
    def citation(self):
        return 'arXiv:1806.10734'
        
    def get_waveform(self, template=None, beta=None, lam=None,
                     chi1=None, chi2=None, m1=None, m2=None, Deltat=None,
                     phi=None, psi=None, inc=None, dist=None):
        """ Return the waveform frequency range, amplitude and phase for
        mode (2,2)
        """
        # TODO: There might be a way to remove duplicates with get_fd_tdiaet
        # Check that lisabeta is installed
        # TODO: to be removed in the end
        try:
            import lisabeta.lisa.lisa as lisa
        except ImportError:
            print("Could not import lisabeta.lisa.lisa")
            return None, None, None

        #Rename parameters to match lisabeta notation
        if template is not None:
            params = self.rename_as_lisabeta(template)
        else:
            params = locals()
            params["lambda"] = lam


        #Generate the waveform of the BHB
        if self.bbh_type.lower() == "mbhb" :
            sig = lisa.GenerateLISATDI_SMBH(params, **self.wvf_pars)
        if self.bbh_type.lower() in ["sobhb", "sbbh"]:
            sig = lisa.GenerateLISATDI_SOBH(params, **self.wvf_pars)

        md = (2,2)
        frq = sig[md]['freq']
        amp = sig[md]['amp']
        phase = sig[md]['phase']
        return frq, amp, phase

    def get_tc(self, template):
        """ Return time of coaelescence from initial frequency
        """
        import lisabeta.waveforms.bbh.pyIMRPhenomD as pyIMRPhenomD
        gridfreq0 = np.array([0.1])
        params = self.rename_as_lisabeta(template)
        syst = pyIMRPhenomD.IMRPhenomDh22AmpPhase(gridfreq0, params['m1'], params['m2'], params['chi1'],
                                                  params['chi2'], params['dist'])
        tc = -syst.compute_toff(params['fstart'])
        return tc

    def get_fd_tdiaet(self, template=None, freqs=None, tdi2=False):
        """  Return TDI A,E,T in frequency domain.
        """
        # Check that lisabeta is installed
        # TODO: to be removed in the end
        try:
            import lisabeta.lisa.lisa as lisa
        except ImportError:
            print("Could not import lisabeta.lisa.lisa")
            return None, None, None

        #Get the waveform parameters
        if template is not None:
            params = self.rename_as_lisabeta(template)
        else:
            params = locals()
            params["lambda"] = lam

        #If the simulation relies on TDI2, use AET observables
        if tdi2:
            self.wvf_pars["TDI"] = "TDI2AET"

        #Generate the TDI waveform of the BHB
        if self.bbh_type.lower() == "mbhb" :
            sig = lisa.GenerateLISATDISignal_SMBH(params, **self.wvf_pars)
        if self.bbh_type.lower() in ["sobhb", "sbbh"]:
            sig = lisa.GenerateLISATDISignal_SOBH(params, **self.wvf_pars)

        #Signal info
        df = 1/self.Tobs
        mds = sig['tdi']['modes']

        #Compute frequency range if not given
        if freqs is None :
            fmin = self.wvf_pars['minf']
            fmax = self.wvf_pars['maxf']
            for mod in mds:
                fmin = min(fmin, np.min(sig['tdi'][mod]['freq']))
                fmax = max(fmax, np.max(sig['tdi'][mod]['freq']))
            fmin = np.floor(fmin/df)*df # redefine f0 to be an integer x df
            freqs = np.arange(fmin, fmax, df)
            linearscale = True
        else:
            df = freqs[1]-freqs[0]
            linearscale = np.isclose(np.sum(np.diff(freqs)-df), 0)

        #kmin = int(np.rint(freqs[0]/df))

        #Get TDI signal
        tdifreqseries = lisa.EvaluateTDIFreqseries(sig["tdi"], freqs)

        #Compute the waveform in the frequendy domain
        A = np.zeros(len(freqs), dtype=np.complex128)
        E = np.zeros(len(freqs), dtype=np.complex128)
        T = np.zeros(len(freqs), dtype=np.complex128)
        for md in mds:
            A += np.conj(tdifreqseries[md]['chan1'])
            E += np.conj(tdifreqseries[md]['chan2'])
            T += np.conj(tdifreqseries[md]['chan3'])

        if linearscale:
            A = FrequencySeries(A, fs=freqs, name="A")
            E = FrequencySeries(E, fs=freqs, name="E")
            T = FrequencySeries(T, fs=freqs, name="T")
        else:
            fs = xr.DataArray(freqs, dims=('f'), attrs={'units': '1/s'})
            A = xr.DataArray(A, dims=('f'), coords={'f': fs}, name='A', attrs={'df': df})
            E = xr.DataArray(E, dims=('f'), coords={'f': fs}, name='E', attrs={'df': df})
            T = xr.DataArray(T, dims=('f'), coords={'f': fs}, name='T', attrs={'df': df})
            
        return A,E,T

    def get_fd_tdixyz(self,**kwargs):
        """  Return TDI X,Y,Z in frequency domain.
        """
        #Get the AET combination
        A,E,T = self.get_fd_tdiaet(**kwargs)
        #Create a TDI object
        tdi = TDI(dict(zip(['A', 'E', 'T'], [A, E, T])))
        #Convert to XYZ combination
        tdi.AET2XYZ()
        return tdi["X"],tdi["Y"],tdi["Z"]

    def get_td_tdiaet(self, **kwargs):
        """  Return TDI A,E,T in time domain.
        """
        #Get the waveform in the frequency domain
        fA, fE, fT = self.get_fd_tdiaet(**kwargs)

        #Compute the waveform in the time domain
        tA = fA.ts.ifft(dt=self.dt)
        tE = fE.ts.ifft(dt=self.dt)
        tT = fT.ts.ifft(dt=self.dt)

        return (tA,tE,tT)

    def get_td_tdixyz(self, **kwargs):
        """  Return TDI X,Y,Z in time domain.
        """
        #Get the AET combination
        A,E,T = self.get_td_tdiaet(**kwargs)
        #Create a TDI object
        tdi = TDI(dict(zip(['A', 'E', 'T'], [A, E, T])))
        #Convert to XYZ combination
        tdi.AET2XYZ()
        return (tdi["X"],tdi["Y"],tdi["Z"])

    def rename_as_lisabeta(self, params_i):
        """ Rename fields to match lisabeta parameter names
        """
        #Copy the parameters
        params = params_i.copy()

        #Spin projection
        if 'PolarAngleOfSpin1' in params.keys():
            params["Spin1"] = params['Spin1']*np.cos(params['PolarAngleOfSpin1'])
            params["Spin2"] = params['Spin2']*np.cos(params['PolarAngleOfSpin2'])
            k_param = ['PolarAngleOfSpin1', 'PolarAngleOfSpin2',
                       'ObservationDuration', 'Cadence', 'Redshift']
            if self.bbh_type.lower() == "mbhb" :
                psi, incl = tools.aziPolAngleL2PsiIncl(params["EclipticLatitude"],
                                       params["EclipticLongitude"],
                                       params['InitialPolarAngleL'],
                                       params['InitialAzimuthalAngleL'])
                params['Polarization'] = psi
                params['Inclination'] = incl
                k_param.extend(['InitialPolarAngleL','InitialAzimuthalAngleL'])
            for k in k_param :
                if k in params.keys():
                    params.pop(k)

        #Parameters to match between LDC notation and lisabeta
        dmatch = dict({'Mass1':             "m1",
                       'Mass2':             "m2",
                       'Spin1':             "chi1",
                       'Spin2':             "chi2",
                       'Distance':          'dist',
                       'Inclination':       'inc',
                       'EclipticLongitude': "lambda",
                       'EclipticLatitude':  "beta",
                       "InitialFrequency":  "fstart",
                       'Polarization':      'psi',
                       })
        #Adding additionnal parameters
        if self.bbh_type.lower() == "mbhb" :
            dmatch['CoalescenceTime'] ='Deltat'
            dmatch['PhaseAtCoalescence']='phi'

        if self.bbh_type.lower() in ["sobhb", "sbbh"]:
            dmatch['InitialPhase']='phi'

        #Return only the relevant parameters
        new_params = dict()
        if isinstance(params, dict):
            for k,v in params.items():
                if k in dmatch.keys() :
                    new_params[dmatch[k]] = params[k]
        else:
            new_params = rfn.rename_fields(params, dmatch)
        return new_params
