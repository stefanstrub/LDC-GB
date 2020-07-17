import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from ldc.lisa.orbits import Orbits
from ldc.waveform.waveform import HpHc
from ldc.lisa.projection import ProjectedStrain
from lc import run_lisacode
import ldc.io.hdf5 as hdf5io

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
        cat = np.array([(8837553, 4.68832205e-21, -0.70722789, 3.53189901,
                         0.01002696, 7.02742938e-15, 0.84565064, 0.21987979,
                         2.48115233)],
                       dtype=[('Name', '<i8'), ('Amplitude', '<f8'),
                              ('EclipticLatitude', '<f8'),
                              ('EclipticLongitude', '<f8'),
                              ('Frequency', '<f8'), ('FrequencyDerivative', '<f8'),
                              ('Inclination', '<f8'), ('InitialPhase', '<f8'),
                              ('Polarization', '<f8')])
    elif key == "big-mbhb-1":
        cat = np.array([(-0.30300442, 1.29251839, 1.20313618, 2.09730354, 0.747377,
                         0.8388, 483052., 223583., 11526944.92187926, 1.22019689,
                         2.69198245, 1.8083985, 1.73941, 13.47098356,
                         31558149.7635456, 5.)],
                       dtype=[('EclipticLatitude', '<f8'), ('EclipticLongitude', '<f8'),
                              ('PolarAngleOfSpin1', '<f8'), ('PolarAngleOfSpin2', '<f8'),
                              ('Spin1', '<f8'), ('Spin2', '<f8'), ('Mass1', '<f8'),
                              ('Mass2', '<f8'), ('CoalescenceTime', '<f8'),
                              ('PhaseAtCoalescence', '<f8'),
                              ('InitialPolarAngleL', '<f8'),
                              ('InitialAzimuthalAngleL', '<f8'), ('Redshift', '<f8'),
                              ('Distance', '<f8'), ('ObservationDuration', '<f8'),
                              ('Cadence', '<f8')])
    elif key == "big-mbhb-2":
        cat = np.array([(1.28888303, 2.10716542, 0.62478351, 1.60719115, 0.953675,
                         0.954511, 294021., 260563., 20426222.35914461, 3.43517638,
                         2.06463032, 4.28787516, 3.58897, 32.30154329,
                         31558149.7635456, 5.)],
                       dtype=[('EclipticLatitude', '<f8'), ('EclipticLongitude', '<f8'),
                              ('PolarAngleOfSpin1', '<f8'), ('PolarAngleOfSpin2', '<f8'),
                              ('Spin1', '<f8'), ('Spin2', '<f8'), ('Mass1', '<f8'),
                              ('Mass2', '<f8'), ('CoalescenceTime', '<f8'),
                              ('PhaseAtCoalescence', '<f8'),
                              ('InitialPolarAngleL', '<f8'),
                              ('InitialAzimuthalAngleL', '<f8'), ('Redshift', '<f8'),
                              ('Distance', '<f8'), ('ObservationDuration', '<f8'),
                              ('Cadence', '<f8')])
    elif key == "big-mbhb-3":
        cat = np.array([(-0.56410239, 0.61092685, 0.90897235, 1.18169945, 0.972661,
                         0.972862, 319160., 250435., 4800021.15572853, 4.27592931,
                         2.57753889, 4.09455023, 2.18186, 17.75836794,
                         31558149.7635456, 5.)],
                       dtype=[('EclipticLatitude', '<f8'), ('EclipticLongitude', '<f8'),
                              ('PolarAngleOfSpin1', '<f8'), ('PolarAngleOfSpin2', '<f8'),
                              ('Spin1', '<f8'), ('Spin2', '<f8'), ('Mass1', '<f8'),
                              ('Mass2', '<f8'), ('CoalescenceTime', '<f8'),
                              ('PhaseAtCoalescence', '<f8'),
                              ('InitialPolarAngleL', '<f8'),
                              ('InitialAzimuthalAngleL', '<f8'), ('Redshift', '<f8'),
                              ('Distance', '<f8'), ('ObservationDuration', '<f8'),
                              ('Cadence', '<f8')])

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
        X = interpolator(trange-251.45)
        np.save(filename, X)
        return X

def get_lisanode(filename, config, name="X", subtract=None):
    X, jk = hdf5io.load_array(filename, name=name)
    if subtract:
        Xs, jk = hdf5io.load_array(subtract, name=name)
        X[:,1] -= Xs[:,1]
        
    ineg = np.where(X[:,0]>=0)[0][0]
    X = X[ineg:]
    interpolator = spline(X[:,0], X[:,1])
    trange = np.arange(config["t_min"], config["t_max"], config["dt"])
    X = interpolator(trange)
    return X

    
if __name__ == '__main__':

    # load configuration
    dirname = "/home/maude/data/LDC/sangria/1.3"
    filename = os.path.join(dirname, "sangria.h5")
    config = hdf5io.load_config(filename, name="obs/config")
    config["dt"] = 5
    globals().update(config)

    trange = np.arange(t_min, t_max, dt)

    if 0: # mbhb time domain
        key = "big-mbhb-2"
        cat = get_cat(key)
        lisacode = get_lisacode(cat, key, config)#, from_file=False)
        simple = get_simple_tdi(cat, key, config)#, from_file=False)
        dirname = "/home/maude/data/LDC/sangria/1.3"
        lisanode = get_lisanode(os.path.join(dirname, "mbhb-tdi.h5"), config)
        background = get_lisanode(os.path.join(dirname, "sum-tdi.h5"), config,
                                  subtract=os.path.join(dirname, "mbhb-tdi.h5"))
        
        plt.figure(figsize=(8,6))
        plt.subplot(111)
        tmin = int((cat["CoalescenceTime"]-600)/dt)
        tmax = int((cat["CoalescenceTime"]+400)/dt)
        plt.plot(trange[tmin:tmax], simple[tmin:tmax], label="simple", color='k')
        plt.plot(trange[tmin:tmax], lisacode[tmin:tmax], label="lisacode", color='b')
        plt.plot(trange[tmin:tmax], lisanode[tmin:tmax], label="lisanode", color='orange')
        plt.plot(trange[tmin:tmax], background[tmin:tmax], label="background",
                 color='grey', alpha=0.5)
        plt.legend(loc="upper right")

        plt.figure(figsize=(8,6))
        plt.subplot(111)
        plt.plot(trange[tmin:tmax], lisacode[tmin:tmax]-simple[tmin:tmax],
                 label="lisacode-simple", color='b')
        plt.plot(trange[tmin:tmax], lisanode[tmin:tmax]-simple[tmin:tmax],
                 label="lisanode-simple", color='orange')
        plt.plot(trange[tmin:tmax], simple[tmin:tmax]*10/100., label="10% TDI X",
                 color='k', ls='--', alpha=0.5)
        plt.plot(trange[tmin:tmax], background[tmin:tmax], label="background",
                 color='grey', alpha=0.5)
        plt.legend(loc="upper right")

    if 1: # gb freq domain
        key = "big-gb"
        cat = get_cat(key)
        lisacode = get_lisacode(cat, key, config)#, from_file=False)
        simple = get_simple_tdi(cat, key, config)#, from_file=False)
        dirname = "/home/maude/data/LDC/sangria/1.3"
        lisanode = get_lisanode(os.path.join(dirname, "dgb-tdi.h5"), config)
        background = get_lisanode(os.path.join(dirname, "sum-tdi.h5"), config,
                                  subtract=os.path.join(dirname, "dgb-tdi.h5"))
        
        simple_Xf = np.fft.fft(window(trange)*simple)*dt
        lisacode_Xf = np.fft.fft(window(trange)*lisacode)*dt
        lisanode_Xf = np.fft.fft(window(trange)*lisanode)*dt
        background_Xf = np.fft.fft(window(trange)*background)*dt
        N = len(simple_Xf)
        freq = np.fft.fftfreq(N, d=dt)
        
        plt.figure(figsize=(8,6))
        plt.subplot(111)
        plt.title("TDI X in freq. domain: abs value")
        plt.plot(freq[:N//2], np.abs(simple_Xf[:N//2]), label="simple tdi")
        plt.plot(freq[:N//2], np.abs(lisacode_Xf[:N//2]), label="lisacode")
        plt.plot(freq[:N//2], np.abs(lisanode_Xf[:N//2]), label="lisanode")
        plt.plot(freq[:N//2], np.abs(background_Xf[:N//2]), label="background",
                 color='grey', alpha=0.5)
        plt.legend(ncol=1, fontsize=8, loc="lower right")
        plt.axis([0.0100255, 0.0100285, None, None])

        plt.figure(figsize=(8,6))
        plt.subplot(111)
        plt.title("TDI X in freq. domain: abs value of the difference")
        plt.plot(freq[:N//2], np.abs(lisacode_Xf[:N//2]-simple_Xf[:N//2]),
                 label="lisacode-simple")
        plt.plot(freq[:N//2], np.abs(lisanode_Xf[:N//2]-simple_Xf[:N//2]),
                 label="lisanode-simple")
        plt.plot(freq[:N//2], np.abs(background_Xf[:N//2]), label="background",
                 color='grey', alpha=0.5)
        plt.plot(freq[:N//2], np.abs(simple_Xf[:N//2])*1/100., label="1% simple",
                 color='k', ls='--', alpha=0.5)

        plt.legend(ncol=1, fontsize=8, loc="lower right")
        plt.axis([0.0100255, 0.0100285, None, None])
