"""Build TDI from source parameters using various code developed
within LISA.

"""

import numpy as np
from astropy import units as un
import matplotlib.pyplot as plt
import os
from scipy.interpolate import InterpolatedUnivariateSpline
import h5py

from ldc.lisa.projection import ProjectedStrain
from ldc.common.series import TimeSeries, FrequencySeries, TDI
from ldc.lisa.orbits import Orbits
from ldc.common.tools import window
import ldc.waveform.fastGB as fastGB
from ldc.waveform.waveform import LW
from ldc.waveform.lisabeta import FastBHB
from ldc.waveform.waveform import HpHc


import lisagwresponse
import lisaorbits
import lisaconstants

from pytdi.michelson import X1, Y1, Z1
from pytdi.michelson import X2, Y2, Z2

from pytdi import Data

ASTRONOMICAL_YEAR = lisaconstants.ASTRONOMICAL_YEAR


def window(tm, show=False):
    xl = 1000.0
    ind_r = np.argwhere(tm[-1]-tm <= 1000.0)[0][0]
    xr = tm[ind_r]
    kap = 0.005
    winl = 0.5*(1.0 + np.tanh(kap*(tm-xl)))
    winr = 0.5*(1.0 - np.tanh(kap*(tm-xr)))
    if show:
        plt.plot(tm, winl)
        plt.plot(tm, winr)
        plt.grid(True)
        plt.show()
    return (winl*winr)

class TDIFactory():
    """
    """
    
    def __init__(self, param=None, approximant='TD_fdot', source_type="GB", Xonly=True, tdi2=False, dt=5,
                 duration=60*60*24*30*6):
        """ Set source parameters. 
        """
        self.xonly = Xonly
        self.approx = approximant
        self.source_type = source_type
        self.tdi2 = tdi2
        self.dt = dt
        self.t_max = duration
        self.t_min = 0
        
        if param is None:
            param = self.get_default_param()
        if isinstance(param , np.ndarray):
            param = dict(zip(param.dtype.names, np.atleast_1d(param)[0]))
        self.param = param
            
        self.orbits = Orbits.type(dict({'orbit_type':'analytic', 'nominal_arm_length':2.5e9,
                                       "initial_position": 0, "initial_rotation": 0}))
        
    def set_custom_params(self, param):
        self.param = param

    def get_default_param(self):
        if self.approx == "TD_fdot":
            return dict({'Amplitude': 1.07345e-22,
                         'EclipticLatitude': 0.312414*un.rad, #np.pi, #
                         'EclipticLongitude': -2.75291*un.rad, #np.pi/2., #,
                         'Frequency': 0.00135962*un.Hz,
                         'FrequencyDerivative': 8.94581279e-19*un.Unit('Hz2'),
                         'Inclination': 0.523599*un.rad,
                         'InitialPhase': 3.0581565*un.rad,
                         'Polarization': 3.5621656*un.rad})
        elif self.approx == "IMRPhenomD" and self.source_type == "MBHB":
            # return dict({'Mass1': 832628.202,
            #              'Mass2': 700997.2481,
            #              'Spin1': 0.9481998052314212, 
            #              'Spin2': 0.9871324769575264, 
            #              'EclipticLatitude': 0.312414,
            #              'EclipticLongitude': -2.75291,
            #              'Redshift': 2.0,
            #              'Distance': 15974.456786495544,
            #              'Cadence': self.dt,
            #              'ObservationDuration': self.t_max,
            #              'CoalescenceTime': 28086000.0,
            #              'InitialAzimuthalAngleL': 3.9,
            #              'InitialPolarAngleL': 2.3535,
            #              'PhaseAtCoalescence': 3.8,
            #              'PolarAngleOfSpin1': 0.0,
            #              'PolarAngleOfSpin2': 0.0,
            #              })

            return dict({'EclipticLatitude': -0.30300442294174235,
                         'EclipticLongitude': 1.2925183861048521,
                         'PolarAngleOfSpin1': 1.2031361791056812,
                         'PolarAngleOfSpin2': 2.097303543065685,
                         'Spin1': 0.747377,
                         'Spin2': 0.8388,
                         'Mass1': 1323277.47932,
                         'Mass2': 612485.5060299999,
                         'CoalescenceTime': 11526944.921879262,
                         'PhaseAtCoalescence': 1.2201968860015653,
                         'InitialPolarAngleL': 2.6919824500032945,
                         'InitialAzimuthalAngleL': 1.808398497592109,
                         'Redshift': 1.73941,
                         'Distance': 13470.983558972537,
                         'ObservationDuration': self.t_max,#31558149.763545603,
                         'Cadence': self.dt})
        elif self.approx == "IMRPhenomD" and self.source_type == "SBBH":
            ### NOK
            # return dict({"Mass1":                 50.,
            #              "Spin1":                 0.0,
            #              "Mass2":                 40.,
            #              "Spin2":                 0.0,

            #              # "EclipticLatitude":      1.7, # ???
            #              "EclipticLatitude":      0.7,

            #              "EclipticLongitude":     1.0471975511965976,
            #              "Inclination":           1.0471975511965976,
            #              "InitialFrequency":      1.2e-2,
            #              "InitialPhase":          0.7,
            #              "Polarization":          1.2,
            #              "Redshift":              2.0,
            #              "Distance":              15974.456786495544,
            #              'Cadence':               self.dt,                   
            #              'ObservationDuration':   self.t_max,
            #              'PolarAngleOfSpin1': 0.0,
            #              'PolarAngleOfSpin2': 0.0,})
            ### OK
            # return dict({"Mass1":                 46.54269472371296,
            #              "Spin1":                 0.26804245679045086,
            #              "Mass2":                 36.925508361628594,
            #              "Spin2":                 0.12417562853403281,

            #              # "EclipticLatitude":      1.7, # ???
            #              "EclipticLatitude":      0.13566137886653098,

            #              "EclipticLongitude":     0.46209433957127805,
            #              "Inclination":           3.051599397145755,
            #              "InitialFrequency":      0.01284494588573715,
            #              "InitialPhase":          0.0,
            #              "Polarization":          4.447376637333776,
            #              "Redshift":              0.08,
            #              "Distance":              378.5684920200472,
            #              'Cadence':               self.dt,                   
            #              'ObservationDuration':   self.t_max,
            #              'PolarAngleOfSpin1': 0.0,
            #              'PolarAngleOfSpin2': 0.0,})
            ### TESTS
            return dict({"Mass1":                 50.,
                         "Spin1":                 0.,
                         "Mass2":                 40.,
                         "Spin2":                 0.,

                         # "EclipticLatitude":      1.7, # ???
                         "EclipticLatitude":      0.7,

                         "EclipticLongitude":     1.0471975511965976,
                         "Inclination":           1.0471975511965976,
                         "InitialFrequency":      1.2e-2,
                         "InitialPhase":          0.7,
                         "Polarization":          1.2,
                         "Redshift":              2.0,
                         "Distance":              15974.456786495544,
                         'Cadence':               self.dt,                   
                         'ObservationDuration':   self.t_max,
                         'PolarAngleOfSpin1': 0.0,
                         'PolarAngleOfSpin2': 0.0,})
        elif self.approx == "FSEF":
            return dict({
                "MassMBHB": 1e6*un.Msun,
                "MassSOBHB": 10*un.Msun,
                "InitialSemiLatusRect": 12,
                "InitialEccentricity": 0.4,
                'InitialPolarAngleL':  3.9*un.rad,
                'InitialAzimuthalAngleL': 2.3535*un.rad,
                'InitialPhasePhiR': 0.123*un.rad,
                'InitialPhase': 0.456*un.rad,
                "Distance": 400*un.Mpc,
                "Cadence": self.dt,
                "ObservationDuration": self.t_max,
                "eps": 1e-5,
                "EclipticLatitude": 1.7*un.rad,
                "EclipticLongitude": 1.0471975511965976*un.rad,
            })
             
    def ldc(self, param=None, tt_order=0):
        """Compute TDI in time domain from projected strain using LDC toolbox

        """
        GW = HpHc.type("test", self.source_type, self.approx)
        GW.set_param(self.param)
        projector = ProjectedStrain(self.orbits)
        self.yArm = projector.arm_response(self.t_min, self.t_max, self.dt,
                                           GW.split(), tt_order=tt_order)
        self.tvec = np.arange(self.t_min, self.t_max, self.dt)
        if self.xonly:
            tdi = TimeSeries(projector.compute_tdi_x(self.tvec, tdi2=self.tdi2,
                                                     tt_order=tt_order),
                             t0=self.t_min, dt=self.dt)
        else:
            tdi = TDI(
                dict(zip(["X", "Y", "Z"],
                         [TimeSeries(projector.compute_tdi_x(self.tvec, tdi2=self.tdi2,
                                                             tt_order=tt_order),
                                     t0=self.t_min, dt=self.dt),
                          TimeSeries(projector.compute_tdi_y(self.tvec, tdi2=self.tdi2,
                                                             tt_order=tt_order),
                                     t0=self.t_min, dt=self.dt),
                          TimeSeries(projector.compute_tdi_z(self.tvec, tdi2=self.tdi2,
                                                             tt_order=tt_order),
                                     t0=self.t_min, dt=self.dt)])))
        return tdi

    def tofd(self, tdi):
        if isinstance(tdi, TDI):
            pass
        else:
            return tdi.ts.fft(win=window)

    def get_orbits_file(self):
        orbits_fn = 'orbits.h5'
        if os.path.exists(orbits_fn):
            os.remove(orbits_fn)
        tinit = self.orbits.initial_position /(2 * np.pi / ASTRONOMICAL_YEAR)
        m_init = self.orbits.initial_rotation
        L = self.orbits.arm_length
        o = lisaorbits.EqualArmlengthOrbits(L=L, dt=86400, size=10+np.round(self.t_max/86400), tt_order=0,
                                            lambda1=-m_init, m_init1=m_init, tinit=-tinit, t0=-50)
        o.write(orbits_fn)
        return orbits_fn

    def gwresponse_pytdi(self):

        tvec = np.arange(self.t_min, self.t_max, self.dt).astype(np.float64)
        links = ['12', '23', '31', '13', '32', '21']
        orbits_fn = self.get_orbits_file()
        
        if self.approx=="TD_fdot":
            p = self.param
            GB = lisagwresponse.GalacticBinary(
                A=p["Amplitude"], f=p["Frequency"], df=p["FrequencyDerivative"],
                phi0=p["InitialPhase"], iota=p["Inclination"], psi=p["Polarization"],
                orbits=orbits_fn,
                gw_beta=p["EclipticLatitude"], gw_lambda=p["EclipticLongitude"],
                dt=self.dt, t0=self.t_min, size=self.t_max/self.dt)
            #self.yArm = GB.compute_gw_response(links, tvec)

            gwfile = "gw.h5"
            if os.path.exists(gwfile):
                os.remove(gwfile)
            GB.write(path=gwfile) # will recompute response

            data = Data.from_gws(gwfile, orbits_fn)
            if self.tdi2:
                built = X2.build(**data.args)
            else:
                built = X1.build(**data.args)
            data_X = built(data.measurements)
        else:
            raise("Not implemented")
        return TimeSeries(data_X, dt=self.dt, t0=self.t_min)

    def fast(self, option='python'):
        """ Compute fast approximant in Fourier domain. 
        """
        if self.approx == "TD_fdot":
            fast = fastGB.FastGB(delta_t=self.dt, T=self.t_max, orbits=self.orbits) # in seconds
            X,Y,Z = fast.get_fd_tdixyz(template=self.param, oversample=4, tdi2=self.tdi2, radler=(option!='python'))
            tdi = TDI(dict(zip(["X", "Y", "Z"], [X, Y, Z])))
        elif self.approx == "IMRPhenomD":
            fast = FastBHB(self.source_type, T=self.t_max, delta_t=self.dt, approx=self.approx, orbits=self.orbits)
            if self.source_type=='MBHB':
                A,E,T = fast.get_td_tdiaet(template=self.param, tdi2=self.tdi2)
            else:
                A,E,T = fast.get_fd_tdiaet(template=self.param, tdi2=self.tdi2)
            tdi = TDI(dict(zip(["A", "E", "T"], [A, E, T])))
            tdi.AET2XYZ()
        elif self.approximant == "FSEF":
            raise NotImplementedError
        return tdi

    def lw(self):
        """ Compute LW approximant
        """
        lwgb = LW(self.param["Inclination"], self.param["EclipticLatitude"], self.param["EclipticLongitude"],
                  self.param["Polarization"], self.param["InitialPhase"], orbits=self.orbits)
        tvec = np.arange(self.t_min, self.t_max, self.dt).astype(np.float64)

        X, Y, Z = lwgb.get_td_tdixyz(tvec, self.param["Amplitude"], self.param["Frequency"],
                                     self.param["FrequencyDerivative"], source_type=self.source_type)
        return TimeSeries(X, dt=self.dt, t0=self.t_min)

   
if __name__ == '__main__':

    if 0:  # test GB- OK
        facGB = TDIFactory(approximant='TD_fdot', dt=15, duration=60*60*24*365)
        ldcx = facGB.tofd(facGB.ldc())
        gwr_ = facGB.gwresponse_pytdi()
        gwrx = facGB.tofd(gwr_)
        lwx  = facGB.tofd(facGB.lw())
        
        fast = facGB.fast()
        fast2 = facGB.fast(option="C")
        
        plt.figure(figsize=(15,10))
        plt.subplot(211)
        plt.semilogy(ldcx.f, np.abs(ldcx), label='ldc')
        plt.semilogy(fast.f, np.abs(fast.X), label='fast', alpha=0.5)
        plt.semilogy(gwrx.f, np.abs(gwrx), label='gw-response / pytdi', ls='--')
        plt.semilogy(gwrx.f, np.abs(lwx), label='LW', ls='--')

        plt.legend()
        plt.axis([facGB.param["Frequency"]-1e-6, facGB.param["Frequency"]+1e-6, None, None])
        ldc_s = ldcx.sel(f=fast.f, method="nearest")
        ldc_s2 = ldcx.sel(f=fast2.f, method="nearest")
        plt.subplot(212)
        plt.plot(fast.f, np.abs(ldc_s.data-fast.X.data), label="ldc diff with python version")
        plt.plot(fast2.f, np.abs(ldc_s2.data-fast2.X.data), label="ldc diff with C version")
        plt.plot(ldcx.f, np.abs(ldcx.data-gwrx.data), label="ldc diff with gw-response / pytdi")
        plt.plot(ldcx.f, np.abs(ldcx.data-lwx.data), label="ldc diff with LW approx")

        plt.legend()
        plt.axis([facGB.param["Frequency"]-1e-6, facGB.param["Frequency"]+1e-6, None, None])


    if 1: # test MBHB
        facGB = TDIFactory(source_type="MBHB", approximant='IMRPhenomD', dt=5, duration=60*60*24*365)
        ldc = facGB.ldc(tt_order=0)
        #ldc1 = facGB.ldc(tt_order=1)
        fast = facGB.fast()

        plt.figure(figsize=(15,10))
        plt.subplot(211)
        plt.plot(ldc.t, ldc, label='ldc')
        plt.plot(fast.t, fast.X, label='fast')
        plt.legend()
        plt.axis([facGB.param["CoalescenceTime"]-1000,
                  facGB.param["CoalescenceTime"]+500, None, None])
        ldc_s = ldc.sel(t=fast.t, method="nearest")
        #ldc_s1 = ldc1.sel(t=fast.t, method="nearest")
        plt.subplot(212)
        plt.semilogy(fast.t, np.abs(ldc_s.data-fast.X.data), label="ldc00-fast")
        plt.legend()
        plt.grid()
        plt.axis([facGB.param["CoalescenceTime"]-1000,
                  facGB.param["CoalescenceTime"]+500, None, None])

        plt.figure(figsize=(15,10))
        plt.subplot(111)
        plt.semilogy(fast.t, np.abs(ldc_s.data-fast.X.data), label="ldc00-fast")
        #plt.semilogy(fast.t, np.abs(ldc_s1.data-fast.X.data), label="ldc11-fast")
        #plt.semilogy(fast.t, np.abs(ldc_s1.data-ldc_s.data), label="ldc11-ldc00")
        plt.legend()
        plt.grid()
        plt.axis([facGB.param["CoalescenceTime"]-1000,
                  facGB.param["CoalescenceTime"]+500, None, None])
        

    if 0: # test SOBHB
        facGB = TDIFactory(source_type="SBBH", approximant='IMRPhenomD', dt=5, duration=60*60*24*365)
        ldc = facGB.ldc()
        ldc_fd = facGB.tofd(ldc)
        fast = facGB.fast()

        plt.figure()
        plt.subplot(111)
        plt.semilogy(ldc_fd.f, np.abs(ldc_fd), label='ldc')
        plt.plot(fast.f, np.abs(fast.X), label='fast')
        plt.ylabel("|X|")
        plt.legend()

        plt.figure(figsize=(15,10))
        plt.subplot(311)
        plt.plot(ldc_fd.f, ldc_fd.real, label='ldc')
        plt.plot(fast.f, fast.X.real, label='fast')
        plt.ylabel("real X")
        # plt.axis([0.012, 0.01202, None, None])
        plt.subplot(312)
        plt.plot(ldc_fd.f, ldc_fd.imag, label='ldc')
        plt.plot(fast.f, fast.X.imag, label='fast')
        plt.ylabel("imag X")
        # plt.axis([0.012, 0.01202, None, None])
        plt.legend()
        plt.subplot(313)
        plt.plot(ldc_fd.f, np.abs(ldc_fd), label='ldc')
        plt.plot(fast.f, np.abs(fast.X), label='fast')
        plt.ylabel("|X|")
        # plt.axis([0.012, 0.01202, None, None])

        


    if 0: # test EMRI
        facGB = TDIFactory(source_type="EMRI", approximant='FSEF', dt=5, duration=60*60*24*365)
        ldcx = facGB.ldc()
        ldc_fd = facGB.tofd(ldcx)
        #fast = facGB.fast()

        plt.figure(figsize=(15,10))
        plt.subplot(111)
        plt.semilogy(ldc_fd.f, np.abs(ldc_fd), label='ldc')
        plt.legend()

        plt.figure(figsize=(15,10))
        plt.subplot(111)
        plt.plot(ldcx.t, ldcx, label='ldc')
        plt.legend()

    plt.show()
