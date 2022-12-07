import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as un
from astropy.cosmology import z_at_value

from lisaconstants import ASTRONOMICAL_YEAR, c

from lisaorbits.orbits import KeplerianOrbits

from few.trajectory.inspiral import EMRIInspiral

from ldc.lisa.orbits import Orbits
from ldc.lisa.projection import ProjectedStrain
from ldc.waveform.waveform import HpHc
from ldc.common.series.timeseries import FrequencySeries, TimeSeries
from ldc.common.series import TDI
from ldc.common.tools import compute_tdi_snr
from ldc.lisa.noise import AnalyticNoise
from ldc.common.constants.cosmology import ldc_cosmo

MODE = "CHECK"
# MODE = "MC"
ORBIT =  'FILE' # 'ANALYTIC'

tobs_yr = 0.5 # yr --- EMRI1
# tobs_yr = 1.5 # yr --- EMRI2
# tobs_yr = 2.1 # yr --- EMRI3
high_ecc = False

DT = 10 # s
TDI2 = False
SNR_WINDOW = [15, 75]
APPROXIMANT = "FSEF"
EPS = 1e-5
L = 2.5e9 # m
CLIGHT = c

tobs_s = tobs_yr * ASTRONOMICAL_YEAR # s

if ORBIT == 'ANALYTIC':
    lisa_orbits = Orbits.type(
        dict({"nominal_arm_length":L*un.m,
            "initial_rotation":0*un.rad,
            "initial_position":0*un.rad,
            "orbit_type":"analytic"})
    )
elif ORBIT == 'FILE':
    dto = 86400
    orbits_fn = "orbits.h5"
    N = int(tobs_s/dto)+2
    t_min_datetime = datetime.datetime(year=1970, month=1, day=1, hour=1)
    lisa_orbits = KeplerianOrbits(L=L, dt=dto, t0=t_min_datetime, size=N, tt_order=2)
    print(lisa_orbits.t0)
    print(lisa_orbits.dt)
    print(lisa_orbits.size)
    if os.path.exists(orbits_fn):
        os.remove(orbits_fn)
    lisa_orbits.write(orbits_fn, mode='x')
    config = {
        'nominal_arm_length': L, 
        'orbit_type':'file', 
        'filename': orbits_fn,
        'lisaorbit_old_version': True,
    }
    lisa_orbits = Orbits.type(config)
    
P = ProjectedStrain(lisa_orbits)
hphc = HpHc.type("demo", "EMRI", APPROXIMANT)


def estimate_snr(hphc, P, T, add_galnoise=False):
    dt = DT
    P.arm_response(0, T, dt, [hphc])
    trange = np.arange(0, T, dt)

    # compute TDI then move TD -> FD
    tdi = TDI( # TD
        dict(zip(["X", "Y", "Z"], [
            TimeSeries(P.compute_tdi_x(trange, tdi2=TDI2), t0=trange[0], dt=dt),
            TimeSeries(P.compute_tdi_y(trange, tdi2=TDI2), t0=trange[0], dt=dt),
            TimeSeries(P.compute_tdi_z(trange, tdi2=TDI2), t0=trange[0], dt=dt)
    ])))
    X = tdi["X"].ts.fft()
    Y = tdi["Y"].ts.fft()
    Z = tdi["Z"].ts.fft()
    XYZ_slow = TDI( # FD
        dict(zip(["X", "Y", "Z"], [X, Y, Z]))
    )

    # compute SNR LDC from TDI XYZ
    # fmin, fmax = X.f[0], X.f[-1]
    if add_galnoise:
        print("Adding galactic confusion noise...")
        # N_ldc = get_noise_model("SciRDv1", X.f, wd=TOBS/ASTRONOMICAL_YEAR)
        noise_model = AnalyticNoise(X.f, model="SciRDv1",
            wd=T/ASTRONOMICAL_YEAR)
    else:
        # N_ldc = get_noise_model("SciRDv1", X.f)
        noise_model = AnalyticNoise(X.f, model="SciRDv1")
    # psd = N_ldc.psd(X.f, option='X', tdi2=TDI2)
    # dict_snr = compute_tdi_snr(XYZ_slow, N_ldc, AET=False, full_output=True)
    dict_snr = compute_tdi_snr(XYZ_slow, noise_model, 
        AET=False, full_output=True)
    snr = np.sqrt(dict_snr["tot2"])

    return snr


def run_mc():
    a = 0.0
    x0 = 1.0 # face-on
    traj = EMRIInspiral(func="SchwarzEccFlux")

    found = False

    while not found:

        ### draw random physical source parameters
        Mz = np.random.uniform(low=1e5, high=3e6) # Msun
        mz = np.random.uniform(low=10, high=200) # Msun
        # cos_incl = np.random.uniform(low=-0.99, high=0.99)
        p0 = np.random.uniform(low=10, high=16.)
        dist = np.random.uniform(low=700, high=5e3) # Mpc
        if high_ecc:
            e0 = np.random.uniform(low=0.5, high=0.7)
        else:
            e0 = np.random.uniform(low=0.01, high=0.7)
        eclLat = np.random.uniform(low=-np.pi, high=np.pi) # beta
        eclLon = np.random.uniform(low=0, high=2*np.pi) # lambda
        theL = np.random.uniform(low=-np.pi/2., high=np.pi/2)
        phiL = np.random.uniform(low=0, high=np.pi)
        # coherent with iota = 0 (face-on)
        # theL = 0
        # eclLat = 3 * np.pi / 2
        # phiL = 0
        # eclLon = np.pi / 2

        phi0 = 0.
        phiR0 = 0.

        # get non redshifted individaul masses
        z = z_at_value(ldc_cosmo.luminosity_distance, dist * un.Mpc).value
        M = Mz / (1+z)
        m = mz / (1+z)

        ### check parameters validity
        if not ((1e5 < M) & (M < 3e6)):
            print(f"Invalid M ({M})")
            continue
        # check mass ratio
        mu = (m*M) / (m+M)
        if not ((mu / M) < 1e-4):
            print(f"Invalid mu/M ({mu/M})")
            continue
        if not ((5 < mu) & (mu < 50)):
            print(f"Invalid mu ({mu})")
            continue
        # check semilatus rectum
        inf_p0 = 10
        sup_p0 = 16 + 2*e0
        if not ((inf_p0 < p0) & (p0 < sup_p0)):
            print(f"Invalid p0 ({p0}) with (inf={inf_p0} | sup={sup_p0})")
            continue

        # integrate trajectories with inf_p0 and sup_p0
        # to get an hint of the time at plunge
        inf_t, _, _, _, _, _, _ = traj(M, mu, a, inf_p0, e0, x0, T=tobs_yr, dt=DT)
        sup_t, _, _, _, _, _, _ = traj(M, mu, a, sup_p0, e0, x0, T=tobs_yr, dt=DT)
        t, _, _, _, _, _, _ = traj(M, mu, a, p0, e0, x0, T=tobs_yr, dt=DT)

        tpl = t[-1]
        inf_tpl, sup_tpl = inf_t[-1], sup_t[-1]
        if (sup_tpl - inf_tpl) < 0.1:
            continue
        print(f"tfplunge: inf={inf_tpl/ASTRONOMICAL_YEAR} yr and sup={sup_tpl/ASTRONOMICAL_YEAR} yr")
        # pack
        d_pars = {
            'MassMBHB': M,
            'MassSOBHB': m,
            'InitialSemiLatusRect': p0,
            'InitialEccentricity':e0,
            'InitialPolarAngleL': theL,
            'InitialAzimuthalAngleL': phiL,
            'InitialPhasePhiR': phiR0,
            'InitialPhase': phi0,
            'Distance': dist,
            'Cadence': DT,
            'ObservationDuration': tobs_s,
            'eps': EPS,
            'EclipticLatitude': eclLat,
            'EclipticLongitude': eclLon,
        }
        
        ### compute SNR
        # # first compute fast SNR average over sky position and polarization
        # snr_avg = estimate_avg_snr(d_pars)
        # if ((SNR_WINDOW[0] < snr_avg) & (snr_avg < SNR_WINDOW[1])):
        #     # compute more precise but slower SNR
        #     print("Computing SNR...")
        #     snr = estimate_snr(d_pars)
        #     if ((SNR_WINDOW[0] < snr) & (snr < SNR_WINDOW[1])):
        #         print("\n##########################")
        #         print(f"Good candidate for params: {d_pars}")
        #         # print(f"Eccentricity @ plunge={epl}")
        #         print(f"tplunge={tpl/ASTRONOMICAL_YEAR} yr")
        #         print(f"SNR={snr}")
        #         print("##########################")
                
        #         found = True
        #     else:
        #         print(f"Invalid SNR ({snr})")
        # else:
        #     print(f"Invalid SNRAVG ({snr_avg})")

        print("Computing SNR...")
        hphc.set_param(d_pars)
        T = d_pars['ObservationDuration']
        # snr_inst = estimate_snr(hphc, P, d_pars, add_galnoise=False)
        snr_inst = estimate_snr(hphc, P, T, add_galnoise=False)
        if ((SNR_WINDOW[0] < snr_inst) & (snr_inst < SNR_WINDOW[1])):
            # compute SNR with instrumental noise + confusion gal noise
            # snr_inst_galnoise = estimate_snr(hphc, P, d_pars, add_galnoise=True)
            snr_inst_galnoise = estimate_snr(hphc, P, T, add_galnoise=True)
            if np.abs(snr_inst - snr_inst_galnoise) < 30:
                print("\n##########################")
                print(f"Good candidate for params: {d_pars}")
                print(f"tplunge={tpl/ASTRONOMICAL_YEAR} yr")
                print(f"SNR (ins. noise only)={snr_inst}")
                print(f"SNR (ins. noise + gal noise)={snr_inst_galnoise}")
                print("##########################")
                
                found = True
        else:
            print(f"Invalid SNR ({snr_inst})")


def check_params(d_pars, plot=False):
    hphc.set_param(d_pars)
    T = d_pars['ObservationDuration']
    # check SNR with instrumetal noise only
    snr_inst = estimate_snr(hphc, P, T, add_galnoise=False)
    # then repeat for instrumatl + confusion noise
    snr_inst_galnoise = estimate_snr(hphc, P, T, add_galnoise=True)
    # display
    print("Checked SNRs:")
    print(f"SNR (ins. noise only)={snr_inst}")
    print(f"SNR (ins. noise + gal noise)={snr_inst_galnoise}")
    # print also info on individual mass for sbbh
    # z = z_at_value(
    #     ldc_cosmo.luminosity_distance, 
    #     d_pars["Distance"] * un.Mpc
    # ).value
    # print(f"indiv. mass SBBH={d_pars['MassSOBHB'] / (1+z)}")
    # plot if any
    if plot:
        # compute missing things
        hphc.set_param(d_pars)
        t = np.arange(0, tobs_s, DT)
        f = np.linspace(1e-5, 1e-1, 10000)
        n_sci = AnalyticNoise(f, model="SciRDv1")
        n_galn = AnalyticNoise(f, model="SciRDv1", wd=tobs_yr)
        n_sci_asd = np.sqrt(n_sci.psd(f, option='X', tdi2=TDI2))
        n_galn_asd = np.sqrt(n_galn.psd(f, option='X', tdi2=TDI2))
        # > hp hc TD
        hp, hc = hphc.compute_hphc_td(t)
        # > hp hc FD
        hp_ts = TimeSeries(hp, t0=t[0], dt=DT)
        hp_fd = hp_ts.ts.fft()
        # > proj strain
        yArm = P.arm_response(0, tobs_s, DT, [hphc])
        # > TDI X
        X = FrequencySeries(P.compute_tdi_x(t, tdi2=TDI2), df=1/tobs_s)
        # > TEST
        f2_hddot = (hp_fd.f[:-2])**2 * np.diff(hp_fd, n=2)
        # plot
        fix, axarr = plt.subplots(2, 2)

        axarr[0, 0].plot(t, hp)
        axarr[0, 0].set_title("TD hplus")
        axarr[0, 0].set_xlabel("Time (s)")
        axarr[0, 0].grid()

        axarr[0, 1].loglog(f, n_sci_asd, label="inst. noise", alpha=0.7)
        axarr[0, 1].loglog(f, n_galn_asd, label=f"inst. + gal. noise (T={tobs_yr} yr)", alpha=0.7)
        axarr[0, 1].loglog(hp_fd.f, np.abs(hp_fd.data), label="EMRI FD", alpha=0.7)
        axarr[0, 1].set_title("FD |hplus_tilde|")
        axarr[0, 1].set_xlabel("Frequency (Hz)")
        axarr[0, 1].legend()
        axarr[0, 1].grid()

        for i in range(yArm.shape[1]):
            axarr[1, 0].semilogy(t, yArm[:, i], label=i+1)
            break
        axarr[1, 0].set_title("TD proj strain")
        axarr[1, 0].set_xlabel("Time (s)")
        axarr[1, 0].legend(loc='upper left')
        axarr[1, 0].grid()

        axarr[1, 1].loglog(f, n_sci_asd, label="inst. noise", alpha=0.7)
        axarr[1, 1].loglog(f, n_galn_asd, label=f"inst. + gal. noise (T={tobs_yr} yr)", alpha=0.7)
        axarr[1, 1].loglog(X.f, np.abs(X.data), label="EMRI TDI")
        axarr[1, 1].loglog(hp_fd.f[:-2], np.abs(f2_hddot), label="TEST")
        axarr[1, 1].set_title("FD |TDI X|")
        axarr[1, 1].set_xlabel("Frequency (Hz)")
        axarr[1, 1].legend()
        axarr[1, 1].grid()

        plt.show()


if __name__ == "__main__":
    if MODE == "MC":
        print("Run Monte-Carlo over EMRI FSEF source parameters ...")
        run_mc()
    
    elif MODE == "CHECK":
        print("Check selected EMRI FSEF source parameters ...")
        # source parameters

        # EMRI 1
        d_pars = {
            'MassMBHB': 634309.0024459653, 
            'MassSOBHB': 22.94506179281826, 
            'InitialSemiLatusRect': 11.606533831613518, 
            'InitialEccentricity': 0.6791984596426677, 
            'InitialPolarAngleL': -0.26662559756010284, 
            'InitialAzimuthalAngleL': 1.3179624458906156, 
            'InitialPhasePhiR': 0.0, 
            'InitialPhase': 0.0, 
            'Distance': 4316.773690290349, 
            'Cadence': 10, 
            'ObservationDuration': 15779074.881772798, 
            'eps': 1e-05, 
            'EclipticLatitude': -3.086939310025947, 
            'EclipticLongitude': 3.4198050647092626,

            # tplunge=0.5000000000000001 yr
            # SNR (ins. noise only)=27.53332301143208
            # SNR (ins. noise + gal noise)=21.86333879246222
            # m_sobhb / (1+z)=13.593252570040262
        }

        
        # EMRI 2

        # EMRI 3
        print(d_pars)
        check_params(d_pars, plot=True)
    else:
        print(f"Unknown mode ({MODE})")
    