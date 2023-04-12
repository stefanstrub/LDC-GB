import matplotlib.pyplot as plt
import numpy as np

from ldc.lisa.orbits import Orbits
from ldc.common.series import TDI
import ldc.waveform.fastGB as fastGB
from ldc.lisa.projection import ProjectedStrain
from ldc.waveform.waveform import HpHc
from ldc.common.series import TimeSeries, FrequencySeries
import ldc.io.hdf5 as h5io
from ldc.common.tools import window
from ldc.lisa.noise import get_noise_model
import pyfftw
from TMyFGB import MyFGB
from TMyFGB import MyFGBnb
import os
import MCMC_multichain.MCMC_multichain as mcmc




class TestGB:

    def __init__(self, param, dt=5, tmax=60*60*24*365):
        self.pGB = param

        config = {"initial_position": 0, "initial_rotation": 0, 
                  "nominal_arm_length": 2500000000, "orbit_type": 'analytic'}
        self.lisa_orbits = Orbits.type(config)
        self.dt = dt
        self.t_max = tmax
        self.t_min = 0

    def get_ldcGB(self, freq_domain=True):
        Proj = ProjectedStrain(self.lisa_orbits)
        GW = HpHc.type("debug", "GB", "TD_fdot")
        GW.set_param(self.pGB)
        yArm = Proj.arm_response(self.t_min, self.t_max, self.dt, [GW])
        self.t = np.arange(self.t_min, self.t_max, self.dt)
        tdi_x = TimeSeries(Proj.compute_tdi_x(self.t), dt=self.dt, t0=self.t[0])
        tdi_y = TimeSeries(Proj.compute_tdi_y(self.t), dt=self.dt, t0=self.t[0])
        tdi_z = TimeSeries(Proj.compute_tdi_z(self.t), dt=self.dt, t0=self.t[0])
        if freq_domain:
            tdi = TDI(dict({"X":tdi_x.ts.fft(win=window),
                            "Y":tdi_y.ts.fft(win=window),
                            "Z":tdi_z.ts.fft(win=window)}))
        else:
            tdi = TDI(dict({"X":tdi_x,
                            "Y":tdi_y,
                            "Z":tdi_z}))
        return tdi
            
    def get_fastGB(self, oversample=4, N=4096, option='C'):

        if option=='C':
            GB = fastGB.FastGB(delta_t=self.dt, T=self.t_max, orbits=self.lisa_orbits) # in seconds
            Xs, Ys, Zs = GB.get_fd_tdixyz(template=self.pGB, oversample=oversample, simulator='synthlisa')
            tdi = TDI(dict(zip(["X", "Y", "Z"], [Xs, Ys, Zs])))
            self.df = tdi.f[1]-tdi.f[0]
            self.fastGB = GB
            return tdi
            
        elif option=='python':
            mGB = MyFGB(self.pGB, self.t_max, self.dt, self.lisa_orbits)
            mGB.ConstructSlowPart(N=N)
            fctr = self.t_max/(N*2)
            mfr, mX, mY, mZ = mGB.ComputeXYZ_FD()
            mX = FrequencySeries(fctr*mX, df=mfr[1]-mfr[0], kmin=int(mfr[0]/(mfr[1]-mfr[0])))
            self.mGB = mGB
            return mX

    def set_f(self, freq):
        self.f = freq.values

    def get_sangria_new(self, fn, freq_domain=True, fr=None, remove_last=0):
        tdi, descr = h5io.load_array(fn)
        self.t = tdi['t']
        selec = self.t>=0
        if remove_last>0:
            selec[-remove_last:] = False

        self.t = self.t[selec]
        dt = self.t[1]-self.t[0]
        tdi_sangria_x = TimeSeries(tdi["X"][selec], dt=dt, t0=self.t[0])
        tdi_sangria_y = TimeSeries(tdi["Y"][selec], dt=dt, t0=self.t[0])
        tdi_sangria_z = TimeSeries(tdi["Z"][selec], dt=dt, t0=self.t[0])
        if freq_domain:
            tdi_sangria_fd = TDI(dict({"X":tdi_sangria_x.ts.fft(win=window),
                                       "Y":tdi_sangria_y.ts.fft(win=window),
                                       "Z":tdi_sangria_z.ts.fft(win=window)}))
            #assert self.f[1]-self.f[0]==tdi_sangria_fd.f[1]-tdi_sangria_fd.f[0]
            if hasattr(self, 'f'):
                tdi_sangria_fd = tdi_sangria_fd.sel(f=self.f, method="nearest")
            return tdi_sangria_fd
        return TDI(dict({"X":tdi_sangria_x, "Y":tdi_sangria_y, "Z":tdi_sangria_z}))

        
    def get_sangria(self, fn="LDC2_sangria_gdb-tdi_v1_v3U3MxS.h5", freq_domain=True, fr=None,
                    remove_last=0):
        """ Noiseless Sangria data
        """ 
        tdi, descr = h5io.load_array(fn)
        self.t = tdi['X'][:,0]
        selec = self.t>=0
        if remove_last>0:
            selec[-remove_last:] = False

        self.t = self.t[selec]
        dt = self.t[1]-self.t[0]
        tdi_sangria_x = TimeSeries(tdi["X"][selec,1], dt=dt, t0=self.t[0])
        tdi_sangria_y = TimeSeries(tdi["Y"][selec,1], dt=dt, t0=self.t[0])
        tdi_sangria_z = TimeSeries(tdi["Z"][selec,1], dt=dt, t0=self.t[0])
        if freq_domain:
            tdi_sangria_fd = TDI(dict({"X":tdi_sangria_x.ts.fft(win=window),
                                       "Y":tdi_sangria_y.ts.fft(win=window),
                                       "Z":tdi_sangria_z.ts.fft(win=window)}))
            #assert self.f[1]-self.f[0]==tdi_sangria_fd.f[1]-tdi_sangria_fd.f[0]
            if hasattr(self, 'f'):
                tdi_sangria_fd = tdi_sangria_fd.sel(f=self.f, method="nearest")
            return tdi_sangria_fd
        return TDI(dict({"X":tdi_sangria_x, "Y":tdi_sangria_y, "Z":tdi_sangria_z}))


    def get_lw(self, freq_domain=True):

        from LW_simple import LW
        psi_LW = self.pGB["Polarization"]
        dt_lw = self.t[1] - self.t[0]
        iota = self.pGB['Inclination']
        bet = self.pGB['EclipticLatitude']
        lam = self.pGB['EclipticLongitude']
        phi0 = self.pGB['InitialPhase']
        f0 = self.pGB["Frequency"]
        fdot = self.pGB["FrequencyDerivative"]
        ampl = self.pGB['Amplitude']
        LW_GB = LW(iota, bet, lam, psi_LW, phi0)
        LW_GB.ComputeProjections(self.t, wv="hphc")
        om = 2.0*np.pi*(f0 + fdot*self.t)
        phR, RX, RY, RZ = LW_GB.ComputeResponse(om)
        om = 2.0*np.pi*f0
        omdot = 2.0*np.pi*fdot
        Xlw, Ylw, Zlw = LW_GB.ComputeTDTDI(self.t, ampl, phi0, om, omdot, RX, RY, RZ, src="GB")
        tdi = TDI(dict({"X":TimeSeries(Xlw, dt=dt_lw, t0=self.t[0]),
                        "Y":TimeSeries(Ylw, dt=dt_lw, t0=self.t[0]),
                        "Z":TimeSeries(Zlw, dt=dt_lw, t0=self.t[0])}))
        if freq_domain:
            tdi = TDI(dict({"X":tdi["X"].ts.fft(win=window),
                            "Y":tdi["Y"].ts.fft(win=window),
                            "Z":tdi["Z"].ts.fft(win=window)}))
            if hasattr(self, 'f'):
                tdi = tdi.sel(f=self.f, method="nearest")
        return tdi
    
if __name__ == '__main__':
    
    pGB = dict({"Frequency":0.0109966533684,
                "FrequencyDerivative": 2.95291919174e-14, 
                "EclipticLatitude":-0.150495923831, 
                "EclipticLongitude":4.60685356883, 
                "Amplitude":2.60186425363e-22,
                "Inclination":0.389613740033,
                "Polarization":0.605754423063,
                "InitialPhase":4.54588971721})

    #pGB["Frequency"] = 0.010996652928533865
    #pGB["Inclination"] = np.arccos(0.9528075176942717)
    #pGB["Inclination"] = 0.2961048813900182
    #pGB["Amplitude"] = 10**(-21.57)#391221547315)#-21.56313065014808)
    tgb = TestGB(pGB, dt=5, tmax=60*60*24*365)

    
    tdi = dict()
    #tdi["sangria_td"] = tgb.get_sangria(fn=fn, freq_domain=False) # set time range

    ## Fast GB
    tdi["fastgb_c"] = tgb.get_fastGB(oversample=4, option='C') # XYZ
    tdi["fastgb_py"] = tgb.get_fastGB(N=4096, option='python') # X only
    tgb.set_f(tdi["fastgb_py"].f) # same than tdi["fastgb_py"].f 
    
    ## LW
    #tdi["lw_fd"] = tgb.get_lw(freq_domain=True)
    
    ## sangria
    fn = "LDC2_sangria_gdb-tdi_v1_v3U3MxS.h5"
    tdi["sangria_fd"] = tgb.get_sangria(fn=fn, freq_domain=True)
    tdi["sangria_fd_1"] = tgb.get_sangria(fn=fn, freq_domain=True, remove_last=1)


    fn2 = "/home/maude/data/LDC/sangria/1.7/dgb-tdi-XYZ.h5"
    tdi["sangria_fd_new"] = tgb.get_sangria_new(fn=fn2, freq_domain=True)
    tdi["sangria_fd_new_1"] = tgb.get_sangria_new(fn=fn2, freq_domain=True, remove_last=1)

    #tdi2 = tgb.get_sangria(fn=fn, freq_domain=True, remove_last=1) 
    #tdi["sangria_fd_2"] = tgb.get_sangria(fn=fn, freq_domain=True, remove_last=2)
    #tdi["sangria_fd_3"] = tgb.get_sangria(fn=fn, freq_domain=True, remove_last=3)

    if 0:
        tdi["ldc"] = tgb.get_ldcGB(freq_domain=False)
        N = tdi["ldc"]["X"].size
        tdi1 = TDI(dict({"X":tdi["ldc"].X.ts.fft(win=window, n=N),
                         "Y":tdi["ldc"].Y.ts.fft(win=window, n=N),
                         "Z":tdi["ldc"].Z.ts.fft(win=window, n=N)}))
        tdi2 = TDI(dict({"X":tdi["ldc"].X.ts.fft(win=window, n=N-1),
                         "Y":tdi["ldc"].Y.ts.fft(win=window, n=N-1),
                         "Z":tdi["ldc"].Z.ts.fft(win=window, n=N-1)}))

    if 0:
        tgb2 = TestGB(pGB, dt=15, tmax=10*60*60*24*365)
        tdi3 = tgb2.get_fastGB(oversample=4, option='C') # XYZ
    
        
        plt.figure(figsize=(15,10))
        plt.subplot(111)
        plt.plot(tdi3.f, tdi3.X.values.imag, label='10 years', color='grey', alpha=0.5)
        plt.plot(tdi1.f, tdi1.X.values.imag, label='all', marker='+')#, ls='None')
        plt.plot(tdi2.f, tdi2.X.values.imag, label='minus 1', marker='+')#, ls='None')

    for n, k in [(fn, 'sangria_fd'), (os.path.basename(fn2), 'sangria_fd_new' )]:
        plt.figure()
        plt.title(n)
        #plt.plot(tgb.f, tdi["fastgb_c"].X.values.imag-tdi[k].X.values.imag, color='purple',
        #         label='FastGB - sangria noise free using all samples')
        #plt.plot(tgb.f, tdi["fastgb_c"].X.values.imag-tdi[k+"_1"].X.values.imag, color='orange',
        #         label='FastGB - sangria noise free using all samples minus 1')
        plt.plot(tgb.f, np.abs(tdi["fastgb_py"].values), color='purple', label='FastGB')
        plt.plot(tgb.f, np.abs(tdi[k].X.values), color='orange',
                 label='sangria noise free using all samples', alpha=0.5)
        plt.plot(tgb.f, np.abs(tdi[k+"_1"].X.values), color='k', ls='--',
                 label='sangria noise free using all samples minus 1')
        #plt.plot(tgb.f, np.abs(tdi["fastgb_py"].values), color='k', ls='--', label='pyFastGB')

        #plt.plot(tgb.f, tdi["sangria_fd"].X.values.imag-tdi["lw_fd"].X.imag, color='r', label='sangria - LW')
        #plt.plot(tgb.f, tdi["fastgb_py"].values.imag-tdi["fastgb_c"]["X"].values.imag, alpha=0.5,
        #         color='green', label='myFastGB - FastGB')
        plt.axis([0.010995, 0.01099852, -1.2e-15, 1.2e-15])
        plt.legend()
        plt.ylabel("|X|")
        plt.savefig(n+'.png')
    

stop


if 0:
    Proj = ProjectedStrain(lisa_orbits)
    GW = HpHc.type("debug", "GB", "TD_fdot")
    GW.set_param(pGB)
    yArm = Proj.arm_response(t_min, t_max, dt, [GW])
    tdi_x = TimeSeries(Proj.compute_tdi_x(XYZ.t.values), dt=dt)
    tdi_y = TimeSeries(Proj.compute_tdi_y(XYZ.t.values), dt=dt)
    tdi_z = TimeSeries(Proj.compute_tdi_z(XYZ.t.values), dt=dt)
    tdi_td = TDI(dict({"X":tdi_x, "Y":tdi_y, "Z":tdi_z}))
    #tdi_td.XYZ2AET()
    plt.figure()
    plt.plot(XYZ.t, XYZ["X"], label='fastGB')
    plt.plot(tdi_td.t, tdi_td["X"], alpha=0.5, label='LDC sangria like')
    gw = np.load("XYZ.npy")
    plt.plot(gw["t"], gw["X"], label="GW response + pyTDI", alpha=0.5)
    plt.legend()
    plt.ylabel("X")
    plt.savefig("X.png")
    plt.figure()
    plt.plot(XYZ.t, XYZ["Y"], label='fastGB')
    plt.plot(tdi_td.t, tdi_td["Y"], alpha=0.5, label='LDC sangria like')
    plt.plot(gw["t"], gw["Y"], label="GW response + pyTDI", alpha=0.5)
    plt.legend()
    plt.ylabel("Y")
    plt.savefig("Y.png")
    plt.figure()
    plt.plot(XYZ.t, XYZ["Z"], label='fastGB')
    plt.plot(tdi_td.t, tdi_td["Z"], alpha=0.5, label='LDC sangria like')
    plt.plot(gw["t"], gw["Z"], label="GW response + pyTDI", alpha=0.5)
    plt.legend()
    plt.ylabel("Z")
    plt.savefig("Z.png")
    
    #plt.figure()
    #plt.plot(XYZ2.t, XYZ2["A"])
    #plt.plot(tdi_td.t, tdi_td["A"], alpha=0.5)
    stop




## Timedomain TDI
if 0:
    Proj = ProjectedStrain(lisa_orbits)
    GW = HpHc.type("debug", "GB", "TD_fdot")
    GW.set_param(pGB)
    yArm = Proj.arm_response(t_min, t_max, dt, [GW])
    tdi_x = TimeSeries(Proj.compute_tdi_x(tdi_sangria_x.t.values), dt=dt)
    tdi_y = TimeSeries(Proj.compute_tdi_y(tdi_sangria_x.t.values), dt=dt)
    tdi_z = TimeSeries(Proj.compute_tdi_z(tdi_sangria_x.t.values), dt=dt)
    tdi_fd = TDI(dict({"X":tdi_x.ts.fft(win=window), "Y":tdi_y.ts.fft(win=window), "Z":tdi_z.ts.fft(win=window)}))
    tdi_fd.XYZ2AET()
    subset = tdi_fd.sel(f=AET.f, method="nearest")


## TDI from gw-response / pytdi
if 0:
    XYZ_td = np.load("XYZ.npy")
    tdi_gw = TDI(dict({"X":TimeSeries(XYZ_td["X"], dt=XYZ_td["t"][1]-XYZ_td["t"][0], t0=XYZ_td["t"][0]),
                       "Y":TimeSeries(XYZ_td["Y"], dt=XYZ_td["t"][1]-XYZ_td["t"][0], t0=XYZ_td["t"][0]),
                       "Z":TimeSeries(XYZ_td["Z"], dt=XYZ_td["t"][1]-XYZ_td["t"][0], t0=XYZ_td["t"][0])}))
    #tdi_gw = tdi_gw.sel(t=tdi["X"], method='nearest')
    for k in ["X", "Y", "Z"]:
        tdi_gw[k] = tdi_gw[k].ts.fft(win=window)
    tdi_gw.XYZ2AET()
    subset_gw = tdi_gw.sel(f=AET.f, method="nearest")
    
if 1:
    plt.figure(figsize=(15,10))
    plt.subplot(111)
    #plt.plot(AET.f, AET["X"].imag, label='FastGB')
    #plt.plot(subset_lw.f, subset_lw.values.imag, label='LW')
    
    plt.plot(mX.f, mX.values.imag-subset_sangria.X.values.imag, color='purple', label='myFastGB - sangria noise free')
    plt.plot(tdi_sangria_fd.f, tdi_sangria_fd["X"].imag-tdi_fd_lw.imag, color='r', label='sangria - LW')

    plt.plot(mX.f, mX.values.imag-AET["X"].values.imag, alpha=0.5, color='green', label='myFastGB - FastGB')

    plt.axis([0.010995, 0.01099852, -1.2e-16, 1.2e-16])#-1.1e-15, 1.1e-15])
    plt.legend()
    plt.ylabel("X (imag part)")
    # plt.subplot(212)
    # #plt.plot(AET.f, AET["E"].imag, label='FastGB')
    # #plt.plot(tdi_fd.f, tdi_fd["E"].imag,  label='sangria like')
    # #plt.plot(tdi_sangria_fd.f, tdi_sangria_fd["E"].imag,  label='sangria noise free', ls='--')
    # plt.plot(AET.f, AET["E"].values.imag-subset_sangria.E.values.imag, alpha=0.5, color='grey', label='FastGB - sangria noise free')

    # #plt.plot(tdi_gw.f, tdi_gw["E"].imag,  label='pytdi', ls='--', color='r')
    # #plt.plot(AET.f, AET["E"].values.imag-subset_gw.E.values.imag, ls='--', alpha=0.5, color='orange', label='FastGB - pytdi')
    # #plt.plot(AET.f, subset_sangria.E.values.imag-subset_gw.E.values.imag, alpha=0.5, color='purple', label='sangria noise free - pytdi')
    # #plt.plot(AET.f, AET["E"].values.imag-AET2["E"].values.imag, alpha=0.5, color='r', label='fastGB with - without evol')
    # #plt.plot(AET.f, AET["E"].imag-subset.E.imag, alpha=0.5, color='grey', label='FastGB - sangria like')
    # plt.axis([0.010995, 0.01099852, -4e-17, 4e-17])#-1.1e-15, 1.1e-15])
    # plt.legend()
    # plt.xlabel("Freq [Hz]")
    # plt.ylabel("E (imag part)")
    #plt.savefig("fastGB_vs_window_newf.png")


if 0:
    plt.figure(figsize=(15,10))
    plt.subplot(211)
    plt.plot(AET.f, AET["A"].real, label='FastGB')
    #plt.plot(tdi_fd.f, tdi_fd["A"].real, label='sangria like')
    #plt.plot(tdi_sangria_fd.f, tdi_sangria_fd["A"].real,  label='sangria noise free', ls='--')
    #plt.plot(AET.f, AET["A"].values.real-subset_sangria.A.values.real, alpha=0.5, color='grey', label='FastGB - sangria noise free')

    plt.plot(tdi_gw.f, tdi_gw["A"].real,  label='pytdi', ls='--')
    plt.plot(AET.f, AET["A"].values.real-subset_gw.A.values.real, alpha=0.5, color='grey', label='FastGB - pytdi')

    #plt.plot(AET.f, AET["A"].real-subset.A.real, alpha=0.5, color='grey', label='FastGB - sangria like')
    plt.axis([0.010995, 0.01100, -1.5e-15, 1.5e-15])
    plt.legend()
    plt.ylabel("A (real part)")
    plt.subplot(212)
    plt.plot(AET.f, AET["E"].real, label='FastGB')
    #plt.plot(tdi_fd.f, tdi_fd["E"].real,  label='sangria like')
    #plt.plot(tdi_sangria_fd.f, tdi_sangria_fd["E"].real,  label='sangria noise free', ls='--')
    #plt.plot(AET.f, AET["E"].values.real-subset_sangria.E.values.real, alpha=0.5, color='grey', label='FastGB - sangria noise free')

    plt.plot(tdi_gw.f, tdi_gw["E"].real,  label='pytdi', ls='--')
    plt.plot(AET.f, AET["E"].values.real-subset_gw.E.values.real, alpha=0.5, color='grey', label='FastGB - pytdi')

    
    #plt.plot(AET.f, AET["E"].real-subset.E.real, alpha=0.5, color='grey', label='FastGB - sangria like')
    plt.axis([0.010995, 0.01100, -1.5e-15, 1.5e-15])
    plt.legend()
    plt.xlabel("Freq [Hz]")
    plt.ylabel("E (real part)")
    #plt.savefig("fastGB_vs_window_newf.png")



    
def loglikelihood_dgb(p, full_output=False):
    logA, fr, fdot, sinb, lam, ci, psi, phi0 = p
    inPrior = True

    for i in range(len(p)):
        if (p[i] < prior[i,0] or p[i]>prior[i,1]):
                inPrior=False
    if not inPrior:
        return (-np.inf)
    
    tmpl = dict({"Frequency": fr,
                 "FrequencyDerivative": fdot, 
                 "EclipticLatitude": np.arcsin(sinb), 
                 "EclipticLongitude": lam, 
                 "Amplitude": 10**logA,
                 "Inclination": np.arccos(ci),
                 "Polarization": psi,
                 "InitialPhase": phi0})

    Xs, Ys, Zs = GB.get_fd_tdixyz(template=tmpl, oversample=oversample, simulator='synthlisa')
    tdi = TDI(dict({"X":Xs, 
                    "Y":Ys, 
                    "Z":Zs}))
    tdi.XYZ2AET()
    
    ### overlapping range for the signal
    if (Xs.f[-1] <= fr_min):
        return (-np.inf)
    elif (Xs.f[0] >= fr_max):
        return (-np.inf)

    subset = tdi.sel(f=slice(fr_min, fr_max))#, method="nearest")
    freq = np.array(subset.f)
    At = subset.A.values
    Et = subset.E.values
    
    ### overlapping range for the data
    dsubset = tdi_sangria_fd.sel(f=slice(fr_min, fr_max))
    #dsubset = tdi_fd.sel(f=slice(fr_min, fr_max))
    As = dsubset.A.values
    Es = dsubset.E.values

    At = tdi['A'].interp(f=dsubset.f)
    Et = tdi['E'].interp(f=dsubset.f)

    SA = get_noise_model('SciRDv1', frq=freq).psd(option='A')
    
    sn1 = np.sum(np.real(As*np.conjugate(At) + Es*np.conjugate(Et))/SA ).values
    sn2 = np.sum((np.abs(At)**2 + np.abs(Et)**2)/SA).values
    sn0 = np.sum((np.abs(As)**2 + np.abs(Es)**2)/SA)
    loglik = 4.0*df.values*(sn1 - 0.5*sn2)

    sn1ae = (np.real(As*np.conjugate(At) + Es*np.conjugate(Et))/SA ).values
    sn2ae = ((np.abs(At)**2 + np.abs(Et)**2)/SA).values
    loglikAE = 4.0*df.values*(sn1ae - 0.5*sn2ae)
    
    #loglik = 4.0*df.values*np.abs((sn1 - 0.5*sn2))
    #print(loglik)
    if full_output:
        return loglik, loglikAE, At.f.values
    return loglik

def MakePrior1Src(fr_min, fr_max, f_fctr=1.e-3):
    Amp_bnd = [-24.0, -20.0] ### log10 amplitude
    fr_bnd = np.array([fr_min, fr_max])*f_fctr   ### in Hz
    fdot_bnd = [-1.e-13, 1.e-13]  ### fdot
    sin_bet_bnd = [-1.0, 1.0]
    lam_bnd = [0.0, 2.0*np.pi]
    cos_iota_bnd = [-1.0, 1.0]
    psi_bnd = [-2.0*np.pi, 4.0*np.pi]
    phi0_bnd = [-2.0*np.pi, 4.0*np.pi]
    prior = [Amp_bnd, fr_bnd, fdot_bnd, sin_bet_bnd, lam_bnd, cos_iota_bnd, psi_bnd, phi0_bnd]
    prior = np.array(prior)
    return (prior)

if 0:

    fr_min = 0.010991945713229392
    fr_max = 0.011002061136872575
    
    p = [np.log10(pGB['Amplitude']),  pGB['Frequency'], pGB['FrequencyDerivative'],
         np.sin(pGB['EclipticLatitude']), \
         pGB['EclipticLongitude'], np.cos(pGB['Inclination']), pGB['Polarization'],  pGB['InitialPhase'] ]
    prior = MakePrior1Src(fr_min, fr_max, f_fctr=1.)
    deltas = []
    plt.figure()
    #plt.plot(AET.f, AET["A"].imag, label='FastGB')
    #plt.plot(tdi_sangria_fd.f, tdi_sangria_fd["A"].imag,  label='sangria noise free', ls='--')
    #plt.plot(AET.f, AET["A"].values.imag-subset_sangria.A.values.imag, alpha=0.5, color='grey', label='FastGB - sangria noise free')

    loglik, loglikAE_ref, freq = loglikelihood_dgb(p, full_output=True)
    
    for idelta in [0, 3, 6]:
        new_p = np.cos(pGB["Inclination"])
        delta = new_p*idelta/100. #-5 to 5%
        new_p += delta
        print(new_p)
        p_ = p.copy()
        p_[5] = new_p
        loglik, loglikAE, freq = loglikelihood_dgb(p_, full_output=True)
        plt.plot(freq, loglikAE-loglikAE_ref, label=str(idelta))
        deltas.append(delta)
        p_ = pGB.copy()
        p_["Inclination"] = np.arccos(new_p)
        Xs, Ys, Zs = GB.get_fd_tdixyz(template=p_, oversample=oversample, simulator='synthlisa')
        AET = TDI(dict(zip(["X", "Y", "Z"], [Xs, Ys, Zs])))
        AET.XYZ2AET()
        #plt.plot(AET.f, AET["A"].imag, label=str(idelta))
        #plt.plot(AET.f, AET["A"].values.imag-subset_sangria.A.values.imag, alpha=0.5, label=str(idelta))
        #plt.plot(AET.f, AET["A"].imag-subset.A.imag, alpha=0.5, color='grey', label='FastGB - sangria like')
    #plt.axis([0.010995, 0.01100, -1.5e-15, 1.5e-15])
    plt.axis([0.010995, 0.01100, None, None])

    plt.legend()
    plt.ylabel("A (imag part)")



if 0:

    import corner
     
    fr_min = 0.010991945713229392
    fr_max = 0.011002061136872575

    
    
    p = [np.log10(pGB['Amplitude']),  pGB['Frequency'], pGB['FrequencyDerivative'],
         np.sin(pGB['EclipticLatitude']), \
         pGB['EclipticLongitude'], np.cos(pGB['Inclination']), pGB['Polarization'],  pGB['InitialPhase'] ]
    lbls = [r'$log10(A)$', r'fr', r'$\dot{f}$', r'$\sin(\beta)$', r'$\lambda$', r'$\cos(\iota)$',
            r'$\psi$', r'$phi$']
    prior = MakePrior1Src(fr_min, fr_max, f_fctr=1.)
    
    print ('p=', p)
    print ('prior', prior)
    ll = loglikelihood_dgb(p)
    print ('loglik = ', ll, np.sqrt(2.*ll))

    if 1:
        plt.figure(figsize=(15,15))
        for i_param, precision in [(0, 0.01), (1, 1e-6), (2, 1e-1), (3, 1e0), (4, 1e-1),  (5, 1e0), (6, 1e1)]:
        #for i_param, precision in [(5, 1e0)]:
            deltas = []
            lls = []

            for idelta in range(-10, 10):

                print("-------------")
                print(idelta)

                delta = p[i_param]*idelta*precision/100. #-1 to 1%
                deltas.append(delta)
                p_ = p.copy()
                p_[i_param] += delta
                if i_param==0:
                    print(p_[i_param])
                ll = loglikelihood_dgb(p_)
                lls.append(ll)

            plt.subplot(3,3,i_param+1)
            plt.plot(deltas, lls, marker='+')
            plt.xlabel("delta "+lbls[i_param])
            plt.axvline(x=0)
            plt.savefig("orig_f.png")

    if 0:
        Nchains = 3
        nsrcs = 1
        Nadapt = 1000
        prop_args = {'cov' : None, 'DE_skip': 1000, 'num_modes': 10}

        x1 = [p]*Nchains
        xin = x1

        Ms_dGB = mcmc.MultiChain_MCMC(np.array(prior), Nchains, loglikelihood_dgb, Nadapt, lbls, 1 )
        Ms_dGB.InitialPoints(xin)
        Ms_dGB.InitializeProposals([{'SCAM':70., 'DE':50, 'slice':40}, \
                                    {'SCAM':70., 'DE':40., 'slice':40.},\
                                    {'SCAM':70., 'DE':30., 'slice':40.}], **prop_args)
        #chains = Ms_dGB.runMCMC(500000, verbal=True, printN=10000, verbal_all=True)
        chains = Ms_dGB.runMCMC(50000, verbal=True, printN=1000, verbal_all=True)

        chn_tot = np.array(chains[0].chn)[15000::100,:]
        chn_tot = np.concatenate((chn_tot, np.array(chains[1].chn)[15000::100,:]))
        chn_tot = np.concatenate((chn_tot, np.array(chains[2].chn)[15000::100,:]))
        # print (acor.acor(chn_tot[:, 0])[0])
        # print (np.shape(chn_tot))
        # print (p)
        chn_tot[:, 1] = (chn_tot[:, 1] - p[1])#*Tobs

        fig = corner.corner(chn_tot,  bins=50, hist_kwargs={'density':True, 'lw':3}, labels=lbls, 
                            plot_datapoints=False, fill_contours=False, quantiles=[0.05, 0.5, 0.95],
                            show_titles=True, 
                            color='darkcyan', truths=p, truth_color='k', use_math_test=True,
                            levels=[0.9], title_kwargs={"fontsize": 12})
        plt.show()
