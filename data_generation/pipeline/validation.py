import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from ldc.lisa.orbits import Orbits
from ldc.waveform.waveform import HpHc
from ldc.waveform.waveform import get_fd_tdixyz
from ldc.lisa.projection import ProjectedStrain
from ldc.common.series import TimeSeries
#from lc import run_lisacode
import ldc.io.hdf5 as hdf5io
from ldc.lisa.noise import get_noise_model
from psd import psd as PSD
import scipy

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

def get_cat(key):
    if key == "big-gb":
        cat = np.array([(8837553, 4.688322047828039e-21, -0.7072278874089373,
                         3.5318990119515874,
                         0.01002696164199913, 7.0274293836061735e-15,
                         0.8456506362930373, 0.21987979316696454,
                         2.481152331028798)],
                       dtype=[('Name', '<i8'), ('Amplitude', '<f8'),
                              ('EclipticLatitude', '<f8'),
                              ('EclipticLongitude', '<f8'),
                              ('Frequency', '<f8'), ('FrequencyDerivative', '<f8'),
                              ('Inclination', '<f8'), ('InitialPhase', '<f8'),
                              ('Polarization', '<f8')])
    elif "big-mbhb" in key:
        cat, k = hdf5io.load_array(filename, name="sky/mbhb/cat")
        ik = int(key.split("-")[-1])
        cat = cat[ik:ik+1]
    return cat

def get_GW(cat, key):
    if 'gb' in key:
        GW = HpHc.type("GalBinaries", "GB", "TD_fdot")
    elif 'mbhb' in key:
        GW = HpHc.type("MBHB", "MBHB", "IMRPhenomD")
    pdict = dict(zip(cat.dtype.names, cat[0]))
    GW.set_param(pdict)
    return GW

def get_simple_tdi(cat, key, config, from_file=True):
    filename = 'simple-%s.npy'%key
    if from_file:
        return np.load(filename)
    else:
        GW = get_GW(cat, key)
        orbits = Orbits.type(config)
        P = ProjectedStrain(orbits)
        yArm = P.arm_response(config["t_min"], config["t_max"], config["dt"], [GW],
                              tt_order=config["travel_time_order"])
        trange = np.arange(config["t_min"], config["t_max"], config["dt"])
        X = P.compute_tdi_x(trange)
        np.save(filename, X)
        return X

def get_lisacode(cat, key, config, from_file=True):
    filename = 'lisacode-%s.npy'%key
    if from_file:
        return np.load(filename)
    else:
        GW = get_GW(cat, key)
        X = run_lisacode([GW], config["t_min"], config["t_max"], config["dt"])
        interpolator = spline(X[:,0], X[:,1])
        trange = np.arange(config["t_min"], config["t_max"], config["dt"])
        X = interpolator(trange-251.75)
        np.save(filename, X)
        return X

def get_lisanode(filename, config, name="X", subtract=None):
    X, jk = hdf5io.load_array(filename, name=name)
    if subtract:
        Xs, jk = hdf5io.load_array(subtract, name=name)
        X[:,1] -= Xs[:,1]
    ineg = np.where(X[:,0]>=0)[0][0]
    print(X[ineg:10,0])
    X = X[ineg:-1, 1]
    return X

    
if __name__ == '__main__':

    # load configuration
    dirname = "/home/maude/data/LDC/sangria/1.7"
    filename = os.path.join(dirname, "sangria.h5")
    config = hdf5io.load_config(filename, name="obs/config")

    # replace waveform dt by tdi dt
    tdi_descr = hdf5io.load_config(filename, name="obs/tdi")
    config["dt"] = int(1/(tdi_descr["sampling_frequency"]))
    config["t_min"] = config["t_min"].value
    config["t_max"] = config["t_max"].value
    globals().update(config)
    trange = np.arange(t_min, t_max, dt)
    if 0:
        noise_model = "MRDv1"
        Nmodel = get_noise_model(noise_model, np.logspace(-5, -1, 100))
        Npsd = Nmodel.psd()
        dirname = "/home/maude/data/LDC/sangria/1.6"
        filename = os.path.join(dirname, "sangria.h5")
        #tdi, attr = hdf5io.load_array(filename, name="obs/tdi")
        X = get_lisanode(os.path.join(dirname, "sum-tdi.h5"), config)
        plt.figure(figsize=(8,6))
        f, psdX =  scipy.signal.welch(X, #tdi["X"],
                                      fs=1.0/dt, window='hanning', nperseg=256*256)
        plt.loglog(f, np.sqrt(psdX), label="noise only")
        plt.loglog(Nmodel.freq, np.sqrt(Npsd), label=noise_model, alpha=2)
        plt.legend()
        plt.axis([1e-5, 3e-1, 4e-22, 2e-19])

        
    
    if 1: # mbhb time domain
        key = "big-mbhb-0" #or big-mbhb-9 or big-mbhb-13
        cat = get_cat(key)
        #lisacode = get_lisacode(cat, key, config, from_file=False)
        simple = get_simple_tdi(cat, key, config)#, from_file=False)
        dirname = "/home/maude/data/LDC/sangria/1.7"
        lisanode = get_lisanode(os.path.join(dirname, "mbhb-tdi.h5"), config)
        background = get_lisanode(os.path.join(dirname, "sum-tdi.h5"), config,
                                  subtract=os.path.join(dirname, "mbhb-tdi.h5"))

        plt.figure(figsize=(8,6))
        plt.subplot(111)
        tmin = int((cat["CoalescenceTime"]-600)/dt)
        tmax = int((cat["CoalescenceTime"]+400)/dt)
        plt.plot(trange[tmin:tmax], simple[tmin:tmax], label="simple", color='k')
        #plt.plot(trange[tmin:tmax], lisacode[tmin:tmax], label="lisacode", color='b')
        plt.plot(trange[tmin:tmax], lisanode[tmin:tmax], label="lisanode", color='orange')
        plt.plot(trange[tmin:tmax], background[tmin:tmax], label="background",
                 color='grey', alpha=0.5)
        plt.legend(loc="upper right")

        plt.figure(figsize=(16,6))
        plt.subplot(121)
        #plt.plot(trange[tmin:tmax], lisacode[tmin:tmax]-simple[tmin:tmax],
        #         label="lisacode-simple", color='b')
        plt.plot(trange[tmin:tmax], lisanode[tmin:tmax]-simple[tmin:tmax],
                 label="lisanode-simple", color='orange')
        plt.plot(trange[tmin:tmax], background[tmin:tmax], label="background",
                 color='grey', alpha=0.5)
        plt.legend(loc="upper right")
        
        plt.subplot(122)
        #plt.plot(trange[tmin:tmax], lisacode[tmin:tmax]-simple[tmin:tmax],
        #         label="lisacode-simple", color='b')
        plt.plot(trange[tmin:tmax], lisanode[tmin:tmax]-simple[tmin:tmax],
                 label="lisanode-simple", color='orange')
        plt.plot(trange[tmin:tmax], background[tmin:tmax], label="background",
                 color='grey', alpha=0.5)
        plt.legend(loc="upper right")
        plt.axis([None, None, -1e-21, 1e-21])

    if 0: # gb freq domain
        key = "big-gb"
        cat = get_cat(key)
        lisacode = get_lisacode(cat, key, config)#, from_file=False)
        simple = get_simple_tdi(cat, key, config)#, from_file=False)
        dirname = "/home/maude/data/LDC/sangria/1.7"
        lisanode = get_lisanode(os.path.join(dirname, "dgb-tdi.h5"), config)
        background = get_lisanode(os.path.join(dirname, "sum-tdi.h5"), config,
                                  subtract=os.path.join(dirname, "dgb-tdi.h5"))
        
        simple_Xf = TimeSeries(window(trange)*simple, dt=dt).ts.fft()
        lisacode_Xf = TimeSeries(window(trange)*lisacode, dt=dt).ts.fft()
        lisanode_Xf = TimeSeries(window(trange)*lisanode, dt=dt).ts.fft()
        background_Xf = TimeSeries(window(trange)*background, dt=dt).ts.fft()
        pGB = dict(zip(cat.dtype.names, cat[0]))
        fastGB_Xf, jk, jk = get_fd_tdixyz(dt, t_max, "GB", "TD_fdot", **pGB)
        
        plt.figure(figsize=(8,6))
        plt.subplot(111)
        plt.title("TDI X in freq. domain: abs value")
        plt.plot(simple_Xf.f, np.abs(simple_Xf), label="simple tdi")
        plt.plot(lisacode_Xf.f, np.abs(lisacode_Xf), label="lisacode")
        plt.plot(lisanode_Xf.f, np.abs(lisanode_Xf), label="lisanode")
        plt.plot(background_Xf.f, np.abs(background_Xf), label="background",
                 color='grey', alpha=0.5)
        plt.plot(fastGB_Xf.f, np.abs(fastGB_Xf), label="fastGB", color="k")
        plt.legend(ncol=1, fontsize=8, loc="lower right")
        plt.axis([0.0100255, 0.0100285, None, None])

        plt.figure(figsize=(8,6))
        plt.subplot(111)
        plt.title("TDI X in freq. domain: abs value of the difference")
        plt.plot(simple_Xf.f, np.abs(lisacode_Xf-simple_Xf),
                 label="lisacode-simple")
        plt.plot(simple_Xf.f, np.abs(lisanode_Xf-simple_Xf),
                 label="lisanode-simple")
        plt.plot(simple_Xf.f, np.abs(background_Xf), label="background",
                 color='grey', alpha=0.5)
        plt.plot(simple_Xf.f, np.abs(simple_Xf)*1/100., label="1% simple",
                 color='k', ls='--', alpha=0.5)

        plt.legend(ncol=1, fontsize=8, loc="lower right")
        plt.axis([0.0100255, 0.0100285, 0, 5e-16])
