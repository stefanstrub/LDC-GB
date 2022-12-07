import numpy as np
import matplotlib.pyplot as plt
from ldc.waveform.waveform import HpHc, NumericHpHc
from ldc.waveform.lisabeta import FastBHB
from scipy.interpolate import InterpolatedUnivariateSpline as spline
import lisabeta.lisa.lisa as lisa
import lisabeta.tools.pytools as pytools
import lisabeta.pyconstants as pyconstants
from ldc.lisa.orbits import Orbits
from ldc.lisa.projection import ProjectedStrain
import lisaconstants as constants
from ldc.common import tools
import ldc.io.hdf5 as h5io
from ldc.common.series import TimeSeries,FrequencySeries
from ldc.common.tools import window

YRSID_SI = constants.SIDEREALYEAR_J2000DAY*24*60*60
MTsun = constants.GM_SUN/constants.SPEED_OF_LIGHT**3

def ifft_positivef_nopad(freqseries):
    n = len(freqseries)
    ncol = freqseries.shape[1] - 1
    deltaf = freqseries[1,0] - freqseries[0,0]
    # Input values for the ifft - has no reason to be real, check we have real and imag
    if not ncol==2:
        raise Exception('Incorrect number of columns in array.')
    # Zero-pad to next power of 2, then 0-pad by factor 2 for the negative frequencies
    # freqseries_pad = zeropad(freqseries, extend=1)
    freqseries_pad = freqseries
    npad = len(freqseries_pad)
    # 1D FD values
    fdvals = freqseries_pad[:,1] + 1j*freqseries_pad[:,2]
    # BEWARE: due to the different convention for the sign of Fourier frequencies, we have to reverse the IFFT input
    # Beware also that the FFT treats effectively the initial time as 0
    # BEWARE: in the reversion of the numpy-convention IFFT input, we have to set aside the 0-frequency term
    fdvals_np = np.zeros_like(fdvals, dtype=complex)
    fdvals_np[0] = fdvals[0]
    fdvals_np[1:] = fdvals[1:][::-1]
    # Inverse FFT
    deltat = 1./(npad*deltaf)
    ifftvals_np = 1./deltat * np.fft.ifft(fdvals_np)
    # Rebuild time series from positive and negative times
    tdvals = np.concatenate((ifftvals_np[npad//2:], ifftvals_np[:npad//2]))
    # Rebuild times
    times = deltat*np.arange(-npad//2, npad//2)
    timeseries = np.array([times, np.real(tdvals), np.imag(tdvals)]).T
    return timeseries

def h22_td_lisabeta(params, tc, T, dt, tstart, tend, minf, maxf,
                    deltaminf=None, deltamaxf=None):
    
    tc_round = dt * int(tc / dt)
    
    params_lisabeta = params.copy()
    params_lisabeta['Deltat'] -= tc_round

    wftdi = lisa.GenerateLISATDI_SMBH(params_lisabeta)
    f, amp, phase = wftdi[(2,2)]['freq'], wftdi[(2,2)]['amp'], wftdi[(2,2)]['phase']
    amp_spline = spline(f, amp)
    phase_spline = spline(f, phase)
    
    # n is chosen to be the first power of 2 so that n*dt >= T
    n = (2**int(np.ceil(np.log2(T/dt)))) # make sure n is a power of 2
    # n = (1 / (dt * df)
    tmax = dt * n
    df = 1./tmax
    freq = df*np.arange(0, n)

    freq_r = freq[(minf <= freq) & (freq <= maxf)]
    
    h22 = amp_spline(freq_r) * np.exp(1j*phase_spline(freq_r))
    plt.figure()
    plt.subplot(111)
    plt.plot(freq_r, np.real(h22))
    if (deltaminf is not None) or (deltamaxf is not None):
        w = pytools.window_planck_vec(freq_r, freq_r[0], freq_r[-1], deltaminf, deltamaxf)
        h22 = w * h22
        plt.plot(freq_r, np.real(h22), alpha=0.5)
        
    istart = int(np.rint(freq_r[0]/df))
    iend = int(np.rint(freq_r[-1]/df))

    if 0:
        h22_fs = FrequencySeries(h22, df=df, kmin=istart)
        h22ts = h22_fs.ts.ifft(dt=dt)
        h22_ts = np.zeros((len(h22ts),2), dtype=float)
        h22_ts[:,0] = h22ts.t.values+tc_round
        h22_ts[:,1] = h22ts.values
        mask = (tstart <= h22_ts[:,0]) & (h22_ts[:,0] <= tend)
        h22_ts = h22_ts[mask]
    else:
        h22_fs = np.zeros((len(freq),3), dtype=float)
        h22_fs[:,0] = freq
        h22_fs[istart:iend+1,1] = np.real(h22)
        h22_fs[istart:iend+1,2] = np.imag(h22)
        
        h22_ts = ifft_positivef_nopad(h22_fs)
    
        h22_ts[:,0] += tc_round
        mask = (tstart <= h22_ts[:,0]) & (h22_ts[:,0] <= tend)
        h22_ts = h22_ts[mask]
    
    return h22_ts

def hphc_td_lisabeta(params, tc, T, dt, tstart, tend, minf, maxf,
                     deltaminf=None, deltamaxf=None):
    inc = params['inc']
    phi = params['phi']
    psi = params['psi']
    
    h22_ts = h22_td_lisabeta(params, tc, T, dt, tstart, tend, minf, maxf, deltaminf=deltaminf, deltamaxf=deltamaxf)
    
    spin_wsh_p = pytools.sYlm(2, 2, inc, phi, s=-2)
    spin_wsh_m = pytools.sYlm(2, -2, inc, phi, s=-2)
    
    t = h22_ts[:,0]
    h22vals = h22_ts[:,1] + 1j*h22_ts[:,2]
    #h22vals = h22_ts[:,1]
    
    h_S = spin_wsh_p * h22vals + spin_wsh_m * np.conj(h22vals) # (-1)^l in the 2nd term, here l=2
    hp_S = np.real(h_S)
    hc_S = -np.imag(h_S)
    
    hp = np.cos(2*psi) * hp_S - np.sin(2*psi) * hc_S
    hc = np.sin(2*psi) * hp_S + np.cos(2*psi) * hc_S
    
    return t, hp, hc


src0 = dict({'Mass1': 832628.202,
             'Mass2': 700997.2481,
             'Spin1': 0.9481998052314212, 
             'Spin2': 0.9871324769575264, 
             'EclipticLatitude': 0.312414,
             'EclipticLongitude': -2.75291,
             'Redshift': 2.0,
             'Distance': 15974.456786495544,
             'Cadence': 5,
             'ObservationDuration': 3.15581498e+07,
             'CoalescenceTime': 28086000.0,
             'InitialAzimuthalAngleL': 3.9,
             'InitialPolarAngleL': 2.3535,
             'PhaseAtCoalescence': 3.8,
             'PolarAngleOfSpin1': 0.0,
             'PolarAngleOfSpin2': 0.0,
             })

src0 = dict(zip(['EclipticLatitude', 'EclipticLongitude',
                 'PolarAngleOfSpin1', 'PolarAngleOfSpin2',
                 'Spin1', 'Spin2', 'Mass1', 'Mass2', 'CoalescenceTime',
                 'PhaseAtCoalescence', 'InitialPolarAngleL',
                 'InitialAzimuthalAngleL', 'Redshift',
                 'Distance', 'ObservationDuration', 'Cadence'] ,
                np.array([-0.56410239, 0.61092685, 0.90897235, 1.18169945,
                          0.972661, 0.972862, 1015522.4376,
                          796849.1091, 4800021.15572853,
                          4.27592931, 2.57753889, 4.09455023,
                          2.18186, 17758.36794127,3.15581498e+07 , 3])))

#src0['InitialPolarAngleL'] = 0.7
#src0['InitialAzimuthalAngleL'] = 2.4

hphc = HpHc.type("test", "MBHB", "IMRPhenomD")
hphc.set_param(src0)

dt = 5 # s
t_min = 0
t_max = src0["CoalescenceTime"] + 1000 #"src0["ObservationDuration"]#CoalescenceTime"] + 1000 
t = np.arange(t_min, t_max, dt)

if 0:

    m1s_si = hphc.m1s*constants.SUN_MASS #  code needs masses in kg
    m2s_si = hphc.m2s*constants.SUN_MASS
    Mc = tools.mchirpofm1m2(hphc.m1s, hphc.m2s)
    f0 = tools.newtonianfoft(Mc, 2.0*hphc.Tobs/YRSID_SI)
    f0 = np.floor(f0/hphc.df)*hphc.df # redefine f0 to be an integer x df
    print(f0)

    FBH = FastBHB("MBHB", approx="IMRPhenomD")
    src0_ = FBH.rename_as_lisabeta(src0)
    tc = src0['CoalescenceTime']
    tc_round = dt * int(tc / dt)
    params_lisabeta = src0_.copy()
    params_lisabeta['Deltat'] -= tc_round
    wftdi = lisa.GenerateLISATDI_SMBH(params_lisabeta)
    f, amp, phase = wftdi[(2,2)]['freq'], wftdi[(2,2)]['amp'], wftdi[(2,2)]['phase']
    amp_spline = spline(f, amp)
    phase_spline = spline(f, phase)
    
    
    frS, phS, ampS = hphc.MBHB_IMRPhenomD_waveform()
    plt.figure()
    plt.subplot(121)
    plt.loglog(frS, ampS, label="ldc")
    plt.loglog(f, amp, label="lb")
    plt.legend()
    plt.title("amp")
    plt.subplot(122)
    plt.loglog(frS, phS, label="ldc")
    plt.loglog(f, phase, label="lb")
    plt.title("phase")
    plt.legend()
    plt.savefig("amp-phase.png")

if 0:

    FBH = FastBHB("MBHB", approx="IMRPhenomD")
    src0_ = FBH.rename_as_lisabeta(src0)
    tc = src0['CoalescenceTime']
    test_t, test_hp, test_hc = hphc_td_lisabeta(src0_, src0['CoalescenceTime'], #hphc.tc,
                                                2*pyconstants.YRSID_SI,
                                                5., t_min, t_max, 5e-5, 2e-2,
                                                deltaminf=1e-5, deltamaxf=4e-3)
    hp, hc = hphc.compute_hphc_td(t, tapering=True) # takes ~40sec
    hp2, hc2 = hphc.compute_hphc_td(t, tapering=False) 

    plt.figure()
    plt.plot(t[::100], hp[::100], label="with tapering")
    plt.plot(t[::100], hp2[::100], alpha=0.5, label="without tapering")
    plt.plot(test_t[::100], test_hp[::100], alpha=0.5, label="hphc_td_lisabeta")
    plt.xlabel("Time [s]")
    plt.ylabel("h+")
    plt.legend()
    # #plt.savefig("hp_.png")

    plt.figure()
    hp_reduced = hp[0:len(test_hp)]
    #hp2_reduced = hp2[0:len(test_hp)]
    plt.plot(test_t[::100], np.abs(test_hp[::100]), alpha=0.5, label="hphc_td_lisabeta")
    plt.semilogy(test_t[::100], np.abs(test_hp[::100]-hp_reduced[::100]),
             alpha=1, label="hphc_td_lisabeta - ldc with tapering")
    #plt.semilogy(test_t[::100], np.abs(test_hp[::100]-hp2_reduced[::100]),
    #         alpha=1, label="hphc_td_lisabeta - ldc without tapering")
    plt.xlabel("Time [s]")
    plt.ylabel("h+ difference")
    plt.legend()
    #plt.savefig("hpdiff.png")
    
if 1:
    config = dict({"nominal_arm_length":2.5e9,
                   "initial_rotation":0,
                   "initial_position":0,
                   "orbit_type":"analytic"})
    orbits = Orbits.type(config)
    P1 = ProjectedStrain(orbits)
    #P2 = ProjectedStrain(orbits)

    hphc = HpHc.type("test", "MBHB", "IMRPhenomD")
    hphc.set_param(src0)
    y1 = P1.arm_response(t_min, t_max, dt, [hphc], tt_order=0, tapering=True)
    #hp,hc = hphc.interp_hphc(t)
    #y2 = P2.arm_response(t_min, t_max, dt, [hphc], tt_order=0, tapering=False)

    FBH = FastBHB("MBHB", approx="IMRPhenomD")
    src0_ = FBH.rename_as_lisabeta(src0)
    tc = src0['CoalescenceTime']
    test_t, test_hp, test_hc = hphc_td_lisabeta(src0_, src0['CoalescenceTime'], #hphc.tc,
                                                2*pyconstants.YRSID_SI,
                                                5., t_min-500, t_max+500, 5e-5, 2e-2,
                                                deltaminf=1e-5, deltamaxf=4e-3)
    #hp_spline = spline(test_t, test_hp)
    #hc_spline = spline(test_t, test_hc)
    #extended_trange = np.arange(t_min-500, t_max+500, dt)
    # istart = np.where(extended_trange==test_t[0])[0][0]
    # istop = np.where(extended_trange==test_t[-1])[0][0]
    
    #hp_ = hp_spline(extended_trange)
    #hc_ = hc_spline(extended_trange)

    # plt.figure()
    # #plt.subplot(211)
    # plt.plot(test_t[100:], test_hp[100:], label="lb")
    # plt.plot(t, hp, label="ldc", alpha=0.5)
    # plt.legend(loc="upper left")
    # plt.figure()
    # #plt.subplot(212)
    # plt.semilogy(t[::100], np.abs(test_hp[100:-100:100]), label="lb")
    # plt.semilogy(t[::100], np.abs(hp[::100]-test_hp[100:-100:100]), label="ldc-lb")
    # plt.legend(loc="upper right")
    # stop
    #plt.plot(extended_trange[::100], hp_[::100])
    
    #hp_[istop:] = hp_[istop];    hp_[0:istart] = hp_[istart]
    #hc_[istop:] = hc_[istop];    hc_[0:istart] = hc_[istart]

    if 1:
        hphc_num = NumericHpHc(test_t, test_hp, test_hc, src0['EclipticLongitude'],
                               src0['EclipticLatitude'])
        P3 = ProjectedStrain(orbits)
        y3 = P3.arm_response(t_min, t_max, dt, [hphc_num], tt_order=0)

    #stop
    # plt.figure()
    # plt.semilogy(t[::100], np.abs(y1[::100,0]), label="with tapering")
    # plt.semilogy(t[::100], np.abs(y2[::100,0]), label="without tapering", alpha=0.5)
    # plt.semilogy(t[::100], np.abs(y3[::100,0]), label="h+,hx lb", alpha=0.5)
    # plt.semilogy(t[::100], np.abs(y3[::100,0]-y1[::100,0]), label="diff", alpha=0.5)
    # plt.legend()
    # plt.xlabel("Time [s]")
    # plt.ylabel("|y1_2|")
    #plt.savefig("y12.png")

    X1 = P1.compute_tdi_x(t, tdi2=False, tt_order=0)
    #X2 = P2.compute_tdi_x(t, tdi2=False, tt_order=0)
    X3 = P3.compute_tdi_x(t, tdi2=False, tt_order=0)

    # plt.figure()
    # plt.semilogy(t[::100], np.abs(X1[::100]), label="with tapering")
    # plt.plot(t[::100], np.abs(X2[::100]), label="without tapering", alpha=0.5)
    # plt.plot(t[::100], np.abs(X3[::100]), label="h+,hx lb", alpha=0.5)
    # plt.legend()
    # plt.xlabel("Time [s]")
    # plt.ylabel("|X|")
    #plt.savefig("tdix.png")

    FBH = FastBHB("MBHB", approx="IMRPhenomD", delta_t=dt, T=t[-1]+dt, orbits=orbits)
    FBH.wvf_pars['minf']=0
    Xlb,Ylb,Zlb = FBH.get_fd_tdixyz(template=src0)

    w = pytools.window_planck_vec(t, 0., src0["CoalescenceTime"]+1200., 1e6, 600.)
    
    X1f = TimeSeries(X1*w,dt=dt, t0=0).ts.fft()
    #X2f = TimeSeries(X2*w,dt=dt).ts.fft()#win=window)
    X3f = TimeSeries(X3*w,dt=dt, t0=0).ts.fft()#win=window)

    X1f_ds = X1f.sel(f=Xlb.f)
    #X2f_ds = X2f.sel(f=Xlb.f)
    X3f_ds = X3f.sel(f=Xlb.f)

    plt.figure()
    np.abs(X1f).plot(label="ldc with tapering")
    #np.abs(X2f).plot(label="ldc without tapering")
    np.abs(X3f).plot(label="lb h+,hx", alpha=0.5)
    np.abs(X1f-X3f).plot(label="with tapering - lb h+,hx")
    #np.abs(X1f-X2f).plot(label="with - without tapering")
    plt.legend(loc='lower right')
    plt.xscale("log")
    plt.yscale("log")
    plt.axis([2e-5, 3e-2, 1e-24, 1e-16])
    plt.grid()
    plt.savefig("tdifdx1.png")

    
    plt.figure()
    np.abs(Xlb).plot(label="lisabeta", color="k")
    np.abs(X1f_ds-Xlb).plot(label="with tapering - lisabeta")
    np.abs(X3f_ds-Xlb).plot(label="h+hx lb - lisabeta", alpha=0.5)
    plt.legend(loc='lower right')
    plt.xscale("log")
    plt.yscale("log")
    plt.axis([2e-5, 3e-2, 1e-24, 1e-16])
    plt.grid()
    plt.savefig("tdifdx2.png")

    
    
if 0:

    i = 1
    for i in range(0,15):
        cat, units = h5io.load_array("/home/maude/data/LDC/sangria/v2/training.h5",
                                     "sky/mbhb/cat")

        src0 = dict(zip(cat.dtype.names, cat[i]))

        # src0 = dict(zip(['EclipticLatitude', 'EclipticLongitude',
        #          'PolarAngleOfSpin1', 'PolarAngleOfSpin2',
        #          'Spin1', 'Spin2', 'Mass1', 'Mass2', 'CoalescenceTime',
        #          'PhaseAtCoalescence', 'InitialPolarAngleL',
        #          'InitialAzimuthalAngleL', 'Redshift',
        #          'Distance', 'ObservationDuration', 'Cadence'] ,
        #                 np.array([-0.56410239, 0.61092685, 0.90897235, 1.18169945,
        #                           0.972661, 0.972862, 1015522.4376,
        #                           796849.1091, 4800021.15572853,
        #                           4.27592931, 2.57753889, 4.09455023,
        #                           2.18186, 17758.36794127,3.15581498e+07 , 3])))
        
        config = dict({"nominal_arm_length":2.5e9,
                       "initial_rotation":0,
                       "initial_position":0,
                       "orbit_type":"analytic"})
        orbits = Orbits.type(config)
        P1 = ProjectedStrain(orbits)
        hphc = HpHc.type("test", "MBHB", "IMRPhenomD")
        hphc.set_param(src0)


        frS, phS, ampS = hphc.MBHB_IMRPhenomD_waveform()
        #stop
        
        y1 = P1.arm_response(t_min, t_max, dt, [hphc], tt_order=0, tapering=True)
        X1 = P1.compute_tdi_x(t, tdi2=False, tt_order=0)

        #hp,hc = hphc.interp_hphc(t)
        #plt.figure()
        #plt.plot(t[::100], hp[::100], label="with tapering")
        #plt.xlabel("Time [s]")
        #plt.ylabel("h+")
        #plt.legend()

        FBH = FastBHB("MBHB", approx="IMRPhenomD", delta_t=dt, T=t_max, orbits=orbits)
        Xlb,Ylb,Zlb = FBH.get_fd_tdixyz(template=src0)

        from ldc.common.series import TimeSeries
        from ldc.common.tools import window
        X1f = TimeSeries(X1,dt=dt).ts.fft(win=window)
        X1f_ds = X1f.interp(f=Xlb.f)

        plt.figure()
        np.abs(X1f_ds-Xlb).plot(label="with tapering - lisabeta")
        np.abs(X1f).plot(label="ldc with tapering")
        np.abs(Xlb).plot(label="lisabeta", color="k")
        plt.axvline(x=frS[0], color="k", ls='--')
        plt.axvline(x=frS[-1], color="k", ls='--')
        plt.legend(loc='lower left')
        plt.xscale("log")
        plt.yscale("log")
        plt.axis([9e-6, 3e-2, 1e-24, 1e-16])
        plt.grid()
        plt.savefig(f"mbhb_{i}.png")

    
stop

FBH = FastBHB("MBHB", approx="IMRPhenomD")
src0_ = FBH.rename_as_lisabeta(src0)
tc = src0['CoalescenceTime']
test_t, test_hp, test_hc = hphc_td_lisabeta(src0_, src0['CoalescenceTime'], #hphc.tc,
                                            2*pyconstants.YRSID_SI,
                                            5., 0., 5e6, 5e-5, 2e-2,
                                            deltaminf=1e-5, deltamaxf=4e-3)
plt.figure()
plt.plot(t[::100], hp[::100])
plt.plot(test_t[::100], test_hp[::100], alpha=0.5)
#plt.figure()
#plt.semilogy(t[::10000], np.abs(hp[::10000]))

