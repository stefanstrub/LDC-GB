import numpy as np
import matplotlib.pyplot as plt

from tdi_factory import TDIFactory

# TODO: factor 2 manually added in tests below
# likely comes from SPA function.

dt = 5 # s
t_max = 60*60*24*365 # s


def test_sbbh_params1(plot=False):
    param = dict({
        "Mass1":                 50.,
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
        'Cadence':               dt,
        'ObservationDuration':   t_max,
        'PolarAngleOfSpin1': 0.0,
        'PolarAngleOfSpin2': 0.0, })
    print(f"src parameters= {param}")

    facGB = TDIFactory(param=param, source_type="SBBH", approximant='IMRPhenomD',
                       dt=dt, duration=t_max)

    ldc = facGB.ldc()
    ldc_fd = facGB.tofd(ldc)
    fast = facGB.fast()

    tdi_x_ldc = ldc_fd
    tdi_x_fast = fast.X

    if plot:
        plt.figure()
        plt.semilogy(ldc_fd.f, np.abs(tdi_x_ldc), label='ldc')
        plt.plot(fast.f, np.abs(tdi_x_fast), label='fast')
        plt.ylabel("|X|")
        plt.legend()
        plt.show()

    # check amplitude
    rdiff_abs = (np.abs(tdi_x_ldc) - np.abs(tdi_x_fast)) / np.abs(tdi_x_fast)

    amp = rdiff_abs.values
    amp[np.isinf(amp)] = np.nan
    print(np.nanmean(amp))
    assert np.nanmean(amp) < 2e-4

    # check phase
    rdiff_ph_r = (tdi_x_ldc.real - tdi_x_fast.real) / tdi_x_fast.real
    rdiff_ph_i = (tdi_x_ldc.imag - tdi_x_fast.imag) / tdi_x_fast.imag

    if plot:
        # plt.plot(fast.f, rdiff_ph_r, label='diff phase real', alpha=0.5)
        # plt.plot(fast.f, rdiff_ph_r, label='diff phase imag', alpha=0.5)
        plt.plot(ldc_fd.f, ldc_fd.real, label='ldc',alpha=0.5)
        plt.plot(fast.f, fast.X.real, label='fast',alpha=0.5)
        plt.ylabel("arg(X)")
        plt.legend()
        plt.show()

    ph_r = rdiff_ph_r.values
    ph_i = rdiff_ph_i.values
    ph_r[np.isinf(ph_r)] = np.nan
    ph_i[np.isinf(ph_i)] = np.nan
    print(np.nanmean(ph_r))
    print(np.nanmean(ph_i))

    assert np.nanmean(ph_r) < 10
    assert np.nanmean(ph_i) < 10


def test_sbbh_params2(plot=False):
    param = dict({
        "Mass1":                 68.2,
        "Spin1":                 0.,
        "Mass2":                 31.4,
        "Spin2":                 0.,

        # "EclipticLatitude":      1.7, # ???
        "EclipticLatitude":      0.1,

        "EclipticLongitude":     1.5,
        "Inclination":           0.9,
        "InitialFrequency":      1.1e-2,
        "InitialPhase":          0.1,
        "Polarization":          0.5,
        "Redshift":              2.0,
        "Distance":              15974.456786495544,
        'Cadence':               dt,
        'ObservationDuration':   t_max,
        'PolarAngleOfSpin1': 0.0,
        'PolarAngleOfSpin2': 0.0, })
    print(f"src parameters= {param}")

    facGB = TDIFactory(param=param, source_type="SBBH", approximant='IMRPhenomD',
                       dt=dt, duration=t_max)

    ldc = facGB.ldc()
    ldc_fd = facGB.tofd(ldc)
    fast = facGB.fast()

    tdi_x_ldc = ldc_fd
    tdi_x_fast = fast.X

    if plot:
        plt.figure()
        plt.semilogy(ldc_fd.f, np.abs(tdi_x_ldc), label='ldc')
        plt.plot(fast.f, np.abs(tdi_x_fast), label='fast')
        plt.ylabel("|X|")
        plt.legend()
        plt.show()

    # check amplitude
    rdiff_abs = (np.abs(tdi_x_ldc) - np.abs(tdi_x_fast)) / np.abs(tdi_x_fast)

    amp = rdiff_abs.values
    amp[np.isinf(amp)] = np.nan
    print(np.nanmean(amp))
    assert np.nanmean(amp) < 2.7e-5

    # check phase
    rdiff_ph_r = (tdi_x_ldc.real - tdi_x_fast.real) / tdi_x_fast.real
    rdiff_ph_i = (tdi_x_ldc.imag - tdi_x_fast.imag) / tdi_x_fast.imag

    if plot:
        plt.plot(fast.f, rdiff_ph_r, label='diff phase real', alpha=0.5)
        plt.plot(fast.f, rdiff_ph_r, label='diff phase imag', alpha=0.5)
        plt.ylabel("arg(X)")
        plt.legend()
        plt.show()

    ph_r = rdiff_ph_r.values
    ph_i = rdiff_ph_i.values
    ph_r[np.isinf(ph_r)] = np.nan
    ph_i[np.isinf(ph_i)] = np.nan
    print(np.nanmean(ph_r))
    print(np.nanmean(ph_i))

    assert np.nanmean(ph_r) < 1
    assert np.nanmean(ph_i) < 1


def test_sbbh_params3(plot=False):
    param = dict({
        "Mass1":                 18.2,
        "Spin1":                 0.,
        "Mass2":                 11.4,
        "Spin2":                 0.,
        "EclipticLatitude":      1.1,
        "EclipticLongitude":     .5,
        "Inclination":           1.9,
        "InitialFrequency":      0.8e-2,
        "InitialPhase":          0.9,
        "Polarization":          0.1,
        "Redshift":              2.0,
        "Distance":              15974.456786495544,
        'Cadence':               dt,
        'ObservationDuration':   t_max,
        'PolarAngleOfSpin1': 0.0,
        'PolarAngleOfSpin2': 0.0, })
    print(f"src parameters= {param}")

    facGB = TDIFactory(param=param, source_type="SBBH", approximant='IMRPhenomD',
                       dt=dt, duration=t_max)

    ldc = facGB.ldc()
    ldc_fd = facGB.tofd(ldc)
    fast = facGB.fast()

    tdi_x_ldc = ldc_fd
    tdi_x_fast = fast.X

    if plot:
        plt.figure()
        plt.semilogy(ldc_fd.f, np.abs(tdi_x_ldc), label='ldc')
        plt.plot(fast.f, np.abs(tdi_x_fast), label='fast')
        plt.ylabel("|X|")
        plt.legend()
        plt.show()

    # check amplitude
    rdiff_abs = (np.abs(tdi_x_ldc) - np.abs(tdi_x_fast)) / np.abs(tdi_x_fast)

    amp = rdiff_abs.values
    amp[np.isinf(amp)] = np.nan
    print(np.nanmean(amp))
    assert np.nanmean(amp) < 2e-2

    # check phase
    rdiff_ph_r = (tdi_x_ldc.real - tdi_x_fast.real) / tdi_x_fast.real
    rdiff_ph_i = (tdi_x_ldc.imag - tdi_x_fast.imag) / tdi_x_fast.imag

    if plot:
        plt.plot(fast.f, rdiff_ph_r, label='diff phase real', alpha=0.5)
        plt.plot(fast.f, rdiff_ph_r, label='diff phase imag', alpha=0.5)
        plt.ylabel("arg(X)")
        plt.legend()
        plt.show()

    ph_r = rdiff_ph_r.values
    ph_i = rdiff_ph_i.values
    ph_r[np.isinf(ph_r)] = np.nan
    ph_i[np.isinf(ph_i)] = np.nan
    print(np.nanmean(ph_r))
    print(np.nanmean(ph_i))

    assert np.nanmean(ph_r) < 1
    assert np.nanmean(ph_i) < 1
