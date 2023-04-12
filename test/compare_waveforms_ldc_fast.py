"""Compare hp and hc for SBBH between LDC and Fast waveform generation routines
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from ldc.lisa.orbits import Orbits
from ldc.waveform.lisabeta import FastBHB
from ldc.waveform.waveform import HpHc


def plot_compare_fampphi(
        fr_ldc, fr_fast,
        amp_ldc, amp_fast,
        ph_ldc, ph_fast):
    fig, axs = plt.subplots(2)

    axs[0].plot(fr_ldc, ph_ldc, label='ldc')
    axs[0].plot(fr_fast, ph_fast, label='fast')
    axs[0].set_ylabel("Phi")
    # axs[0].axis([0.012, 0.01202, None, None])
    axs[0].legend()

    axs[1].plot(fr_ldc, amp_ldc, label='ldc')
    axs[1].plot(fr_fast, amp_fast, label='fast')
    axs[1].set_ylabel("Amp")
    # axs[1].axis([0.012, 0.01202, None, None])
    axs[1].legend()

    plt.show()


def plot_diff_ampphi(
        fr_ldc, fr_fast,
        amp_ldc, amp_fast,
        ph_ldc, ph_fast):
    # interpolate
    interp_type = 'cubic'
    f_phi_ = interp1d(fr_ldc, ph_ldc, kind=interp_type)
    f_amp_ = interp1d(fr_ldc, amp_ldc, kind=interp_type)
    fnew = np.linspace(fr_ldc[0], fr_ldc[-1], len(fr_ldc))
    ph_fast_new = f_phi_(fnew)
    amp_fast_new = f_amp_(fnew)

    # proper plot
    fig, axs = plt.subplots(2)

    axs[0].plot(fr_ldc, ph_fast_new-ph_ldc, label='fast-ldc')
    axs[0].set_ylabel("Diff Phi")
    # axs[0].axis([0.012, 0.01202, None, None])
    axs[0].legend()

    axs[1].plot(fr_ldc, amp_fast_new-amp_ldc, label='fast-ldc')
    axs[1].set_ylabel("Diff Amp")
    # axs[1].axis([0.012, 0.01202, None, None])
    axs[1].legend()

    plt.show()


def plot_compare_hphc_fd(
        fr_ldc, amp_ldc, ph_ldc,
        fr_fast, amp_fast, ph_fast,
        mode='FD'):
    # only valid for (2, 2) mode - OK
    h_tilde_ldc = amp_ldc * np.exp(-1j * ph_ldc)
    h_tilde_fast = amp_fast * np.exp(-1j * ph_fast)

    hp_tilde_ldc = np.real(h_tilde_ldc)
    hp_tilde_fast = np.real(h_tilde_fast)
    hc_tilde_ldc = - np.imag(h_tilde_ldc)
    hc_tilde_fast = - np.imag(h_tilde_fast)

    # proper plot
    fig, axs = plt.subplots(2)

    axs[0].plot(fr_ldc, hp_tilde_ldc, label='ldc')
    axs[0].plot(fr_fast, hp_tilde_fast, label='fast')
    axs[0].set_ylabel("hp")
    # axs[0].axis([0.012, 0.01202, None, None])
    axs[0].legend()

    axs[1].plot(fr_ldc, hc_tilde_ldc, label='ldc')
    axs[1].plot(fr_fast, hc_tilde_fast, label='fast')
    axs[1].set_ylabel("hc")
    # axs[1].axis([0.012, 0.01202, None, None])
    axs[1].legend()

    plt.show()


if __name__ == "__main__":

    # general params
    approx = "IMRPhenomD"
    source_type = "SBBH"
    dt = 5  # s
    # t_max = 3600*24*365*0.1 # = 0.1yr
    t_max = 3600*24*365*3  # = 3yr
    t_min = 0

    # source parameters
    params = dict({
         "Mass1":                 50.,
         "Spin1":                 0.0,
         "Mass2":                 40.,
         "Spin2":                 0.0,
         "EclipticLatitude":      1.7,
         "EclipticLongitude":     1.0471975511965976,
         "Inclination":           1.0471975511965976,
         "InitialFrequency":      1.2e-2,
         "InitialPhase":          0.7,
         "Polarization":          1.2,
         "Redshift":              2.0,
         "Distance":              15974.456786495544,
         'Cadence':               dt,
         'ObservationDuration':   t_max-t_min,
         'PolarAngleOfSpin1': 0.0,
         'PolarAngleOfSpin2': 0.0,
    })

    # define orbit.
    orbits = Orbits.type(
        dict({'orbit_type': 'analytic',
              'nominal_arm_length': 2.5e9,
              "initial_position": 0,
              "initial_rotation": 0})
    )

    # Compute hp, hc for LDC
    tvec = np.arange(t_min, t_max, dt)
    SOBHB = HpHc.type("my-sobhb", source_type, approx)
    SOBHB.set_param(params)
    fr_ldc, ph_ldc, amp_ldc = SOBHB.SOBHB_IMRPhenomD_waveform()

    # Compute hp, hc for Fast
    fast = FastBHB(source_type, T=t_max, delta_t=dt,
                   approx=approx, orbits=orbits)
    fr_fast, amp_fast, ph_fast = fast.get_waveform(template=params)

    # Plot f, A, phi for LDC and Fast
    # plot_compare_fampphi(
    #     fr_ldc, fr_fast,
    #     amp_ldc, amp_fast,
    #     ph_ldc, ph_fast
    # )

    # Plot diff f, A, phi for LDC and Fast
    plot_diff_ampphi(
        fr_ldc, fr_fast,
        amp_ldc, amp_fast,
        ph_ldc, ph_fast
    )

    # reconstruct hp_tilde, hc_tilde
    # from f, A, phi for LDC and Fast
    # plot_compare_hphc_fd(
    #     fr_ldc, amp_ldc, ph_ldc,
    #     fr_fast, amp_fast, ph_fast
    # )
