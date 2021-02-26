# Standard useful python module
import numpy as np
import sys, os, re
import time
import copy
# Adding the LDC packages to the python path
sys.path.append("/root/.local/lib/python3.6/site-packages/")
# LDC modules
import tdi
from LISAhdf5 import LISAhdf5,ParsUnits
import LISAConstants as LC

from scipy.interpolate import InterpolatedUnivariateSpline as spline
from scipy import interpolate
import Cosmology
import GenerateFD_SignalTDIs as GenTDIFD
import pyIMRPhenomD
import pyFDresponse as FD_Resp


import corner

# Plotting modules
import matplotlib.pyplot as plt
import matplotlib as mpl

import matplotlib.pyplot as plt
import scipy
import numpy as np
import xarray as xr
from astropy import units as u

import ldc.io.hdf5 as hdfio
from ldc.lisa.noise import get_noise_model
from ldc.lisa.orbits import Orbits
from ldc.lisa.projection import ProjectedStrain
from ldc.common.series import TimeSeries, FrequencySeries, window
import ldc.waveform.fastGB as fastGB
from ldc.common.tools import compute_tdi_snr
from ldc.waveform.waveform import HpHc
import time


### returns array of parameters from hdf5 structure
def GetParams(pGW):
    # print (pGW.get('Mass1')*1.e-6, pGW.get('Mass2')*1.e-6)
    m1 = pGW.get('Mass1') ### Assume masses redshifted
    m2 = pGW.get('Mass2')
    tc = pGW.get('CoalescenceTime')
    chi1 = pGW.get('Spin1')*np.cos(pGW.get('PolarAngleOfSpin1'))
    chi2 = pGW.get('Spin2')*np.cos(pGW.get('PolarAngleOfSpin2'))
    theL = pGW.get('InitialPolarAngleL')
    phiL = pGW.get('InitialAzimuthalAngleL')
    longt = pGW.get('EclipticLongitude')
    lat = pGW.get('EclipticLatitude')
    phi0 = pGW.get('PhaseAtCoalescence')
    z = pGW.get("Redshift")
    DL = Cosmology.DL(z, w=0)[0]
    dist = DL * 1.e6 * LC.pc
    print ("DL = ", DL*1.e-3, "Gpc")
    incl = pGW.get('Inclination')
    # print ("Compare DL:", pGW.getConvert('Distance',LC.convDistance,'mpc'))

    bet, lam, incl, psi = GenTDIFD.GetSkyAndOrientation(p)

    return(m1, m2, tc, chi1, chi2, dist, incl, bet, lam, psi, phi0, DL)

def ComputeMBHBtemplate(p, Tobs, dt, fmin):

    MfCUT_PhenomD = 0.2 - 1e-7 ### for IMRPhenomD 
    Mc, q, tc, chi1, chi2, dist, incl, bet, lam, psi, phi0, DL, m1, m2 = p

    # print ("check 1", Mc, q, tc, chi1, chi2, dist, incl, bet, lam, psi, phi0, DL, m1, m2)

    m1_SI = m1*LC.MsunKG
    m2_SI = m2*LC.MsunKG
    Ms=(m1+m2)*LC.MTsun#*solar mass in sec
    df = 1.0/Tobs
    eta = m1*m2/(m1+m2)**2

    f0 = FD_Resp.funcNewtonianfoft(Mc, Tobs/LC.YRSID_SI)
    if (f0<fmin):
        f0 = fmin
 
    fRef = 0.0 # hardcodded  and defines the waveform in the source frame
    maxf = 0.5/dt
    fmax = min(MfCUT_PhenomD/Ms, maxf)

    acc_sampling = 1.e-5  ## hardcoded tolerance for the interpolation
    freq_PhD = 1/Ms * FD_Resp.WaveformFrequencyGridGeom(eta, Ms*f0, Ms*fmax, acc=acc_sampling)

    phiRef = 0.0 # hardcoded
    wf_PhD_class = pyIMRPhenomD.IMRPhenomDh22AmpPhase(freq_PhD, phiRef, fRef, m1_SI, m2_SI, chi1, chi2, dist)
    wf_PhD = wf_PhD_class.GetWaveform()

    frS = np.array(wf_PhD[0])
    phS = np.array(wf_PhD[2])
    ampS = np.array(wf_PhD[1])

    tfspline = spline(frS, 1/(2.*np.pi)*(phS-phS[0])).derivative()
    tf = tfspline(frS)
    Shift = tf[-1] - tc
    tf = tf - Shift


    index_cuttf = 0
    tfdiff = np.diff(tf)
    while index_cuttf<len(tfdiff)-1 and tfdiff[index_cuttf]>0:
        index_cuttf += 1
    tfr = tf[:index_cuttf+1]
    # print ("cutoff:", index_cuttf, len(tf), tf[index_cuttf], tf[-1])
    frS_r = frS[:index_cuttf+1]
    frspl = spline(tfr, frS_r)
    ind = index_cuttf

    if (tf[0] < 0.0):
        f0 = frspl(0.0)
        freq_PhD = 1/Ms * FD_Resp.WaveformFrequencyGridGeom(eta, Ms*f0, Ms*fmax, acc=acc_sampling)

        # wf_PhD_class = pyIMRPhenomD.IMRPhenomDh22AmpPhase(freq_PhD, phi0, fRef, m1_SI, m2_SI, chi1, chi2, dist)
        wf_PhD_class = pyIMRPhenomD.IMRPhenomDh22AmpPhase(freq_PhD, phiRef, fRef, m1_SI, m2_SI, chi1, chi2, dist)
        wf_PhD = wf_PhD_class.GetWaveform()

        frS = np.array(wf_PhD[0])
        phS = np.array(wf_PhD[2])
        ampS = np.array(wf_PhD[1])

        tfspline = spline(frS, 1/(2.*np.pi)*(phS-phS[0])).derivative()
        tf = tfspline(frS)
        Shift = tf[-1] - tc
        tf = tf - Shift

        # plt.plot(tf, frS)
        # plt.show()

    freq_response = FD_Resp.ResponseFrequencyGrid([frS, ampS, phS])
    tm_response = tfspline(freq_response) - Shift

    order = 0
    epsTfvec = None
    Tfvec = None
    wfTDI =  GenTDIFD.JustLISAFDresponseTDI(freq_response, tm_response, Tfvec, epsTfvec, incl, lam, bet, psi, phi0, t0=0.0, \
                    order_fresnel_stencil=order)

    return(frS, ampS, phS, freq_response, wfTDI)


### computes MBHB template for the given frequency array
### used for slow (conventional) likelihood evaluation
def MBHBtmplFine(par, freq, tdi='XYZ'):
    #        0   1   2   3     4     5     6     7    8    9    10   11  12  13
    # parS = Mc, q, tc, chi1, chi2, dist, incl, bet, lam, psi, phi0, DL, m1, m2
    Mc, q, tc, chi1, chi2, logDL, ci, sb, lam, psi, phi0 = par
    DL = 10.0**logDL
    dist = DL * 1.e6 * LC.pc
    m1=0.5*Mc*((q/((q+1.)**2))**(-3.0/5.0) + ((q/((q+1.)**2))**(-6.0/5.0) - 4.0*(q/((q+1.)**2))**(-1.0/5.0))**(0.5))
    m2=0.5*Mc*((q/((q+1.)**2))**(-3.0/5.0) - ((q/((q+1.)**2))**(-6.0/5.0) - 4.0*(q/((q+1.)**2))**(-1.0/5.0))**(0.5))
    incl = np.arccos(ci)
    bet = np.arcsin(sb)
    par_w = Mc, q, tc, chi1, chi2, dist, incl, bet, lam, psi, phi0, DL, m1, m2

    ### Compute the waveform' ingedients on the coarse grid
    frW, ampW, phW, fr_resp, wfTDI = ComputeMBHBtemplate(par_w, Tobs, del_t, fmin=2.e-5)

    fmin = max(freq[0], frW[0], fr_resp[0])
    fmax = min(freq[-1], frW[-1], fr_resp[-1])
    # print ('fmax= ', fmax, "and", freq[-1], frW[-1], fr_resp[-1])
    ind_in = np.argwhere(freq>fmin)[0][0]
    ind_en = np.argwhere(freq>fmax)[0][0]-1

    ### interpolating on a regular grid and building the waveform 
    amp_spl = spline(frW, ampW)
    ph_spl = spline(frW, phW)
    phRspl = spline(fr_resp, wfTDI['phaseRdelay'])

    phasetimeshift = 2.*np.pi*tc*freq[ind_in:ind_en]

    ph = ph_spl(freq[ind_in:ind_en]) + phRspl(freq[ind_in:ind_en]) + phasetimeshift
    amp = amp_spl(freq[ind_in:ind_en])

    keytrs = ['transferL1', 'transferL2', 'transferL3']

    Nf = len(freq)
    X = np.zeros(Nf, dtype=np.complex128)
    Y = np.zeros(Nf, dtype=np.complex128)
    Z = np.zeros(Nf, dtype=np.complex128)
    fast = amp*np.exp(1.j*ph) # to agree with fft conventions
    for I, ky in enumerate(keytrs):
        transferLRespline = spline(fr_resp, np.real(wfTDI[ky]))
        transferLImspline = spline(fr_resp, np.imag(wfTDI[ky]))
        transferLRe = transferLRespline(freq[ind_in:ind_en])
        transferLIm = transferLImspline(freq[ind_in:ind_en])
        if (ky == 'transferL1'):
                X[ind_in:ind_en] = fast*(transferLRe+1j*transferLIm)
        if (ky == 'transferL2'):
                Y[ind_in:ind_en] = fast*(transferLRe+1j*transferLIm)
        if (ky == 'transferL3'):
                Z[ind_in:ind_en] = fast*(transferLRe+1j*transferLIm)

    # plt.semilogx(freq[ind_in:ind_en], np.real(fast))
    # plt.semilogx(freq[ind_in:ind_en], amp*np.cos(ph), '--')
    # plt.show()
    # sys.exit(0)

    if (tdi=='XYZ'):
        return(np.conjugate(X), np.conjugate(Y), np.conjugate(Z))
    else:
        A = (Z - X)/np.sqrt(2.0)
        E = (X - 2.0*Y + Z)/np.sqrt(6.0)
        T = (X+Y+Z)/np.sqrt(3.0)
        return (np.conjugate(A), np.conjugate(E), np.conjugate(T))


### computes MBHB pseudo-template on the coarse unevenly sampled frequency grid
### used for fast noiseless likelihood evaluation
def MBHBtmplCoarse(par):
    #        0   1   2   3     4     5     6     7    8    9    10   11  12  13
    # parS = Mc, q, tc, chi1, chi2, dist, incl, bet, lam, psi, phi0, DL, m1, m2
    Mc, q, tc, chi1, chi2, logDL, ci, sb, lam, psi, phi0 = par
    DL = 10.0**logDL
    dist = DL * 1.e6 * LC.pc
    m1=0.5*Mc*((q/((q+1.)**2))**(-3.0/5.0) + ((q/((q+1.)**2))**(-6.0/5.0) - 4.0*(q/((q+1.)**2))**(-1.0/5.0))**(0.5))
    m2=0.5*Mc*((q/((q+1.)**2))**(-3.0/5.0) - ((q/((q+1.)**2))**(-6.0/5.0) - 4.0*(q/((q+1.)**2))**(-1.0/5.0))**(0.5))
    incl = np.arccos(ci)
    bet = np.arcsin(sb)
    par_w = Mc, q, tc, chi1, chi2, dist, incl, bet, lam, psi, phi0, DL, m1, m2

    frW, ampW, phW, fr_resp, wfTDI = ComputeMBHBtemplate(par_w, Tobs, del_t, fmin=2.e-5)


    frqs = np.unique(np.concatenate((frW, fr_resp)))
    N = len(frqs)
    amp_spl = spline(frW, ampW)
    ph_spl = spline(frW, phW)
    phRspl = spline(fr_resp, wfTDI['phaseRdelay'])
    amp = amp_spl(frqs)
    phase = ph_spl(frqs)
    phasetimeshift = 2.*np.pi*tc*frqs

    ph = ph_spl(frqs) + phRspl(frqs) + phasetimeshift

    keytrs = ['transferL1', 'transferL2', 'transferL3']

    for I, ky in enumerate(keytrs):
        transferLRespline = spline(fr_resp, np.real(wfTDI[ky]))
        transferLImspline = spline(fr_resp, np.imag(wfTDI[ky]))
        transferLRe = transferLRespline(frqs)
        transferLIm = transferLImspline(frqs)
        if (ky == 'transferL1'):
                X = amp*(transferLRe+1j*transferLIm)
        if (ky == 'transferL2'):
                Y = amp*(transferLRe+1j*transferLIm)
        if (ky == 'transferL3'):
                Z = amp*(transferLRe+1j*transferLIm)

    wvf = {}
    wvf['freq'] = frqs
    wvf['phase'] = ph
    wvf['phPD'] = ph_spl(frqs)
    wvf['phRD'] = phRspl(frqs)
    wvf['ph_shift'] = phasetimeshift
    XYZ = False
    if (XYZ):
        wvf['ampl'] = np.array([X, Y, Z])
    else:
        A = (Z - X)/np.sqrt(2.0)
        E = (X - 2.0*Y + Z)/np.sqrt(6.0)
        T = (X+Y+Z)/np.sqrt(3.0)
        wvf['ampl'] = np.array([A, E, T])

    return (wvf)

def ComputeLogLikSlow(data, par):

    [Mc, q, tc, chi1, chi2, logDL, ci, sb, lam, psi, phi0] = par

    passit = CheckPrior(Mc, q, tc, chi1, chi2, logDL, ci, sb, lam, psi, phi0)

    if (not passit):
        return (-np.inf)
    else:
        freq = data[0]
        Ad = data[1]
        Ed = data[2]
        Et = data[3]

        SA = tdi.noisepsd_AE(freq, model='Proposal', includewd=None)
        ST = tdi.noisepsd_T(freq, model='Proposal', includewd=None)

        print ("check:", np.shape(Ad), np.shape(At), np.shape(Ed), np.shape(Et))

        df = 1.0/Tobs

        print ("chack df:", df, freq[2]-freq[1], freq[1])

        At, Et, Tt = MBHBtmplFine(par, freq, tdi='AET')

        ## we skip T-channel
        SN = np.sum( np.real( Ad * np.conjugate(At) + Ed*np.conjugate(Et) )/SA )
        hh = np.sum( (np.abs(At)**2 + np.abs(Et)**2)/SA )


        return ( 4.0*df*(SN - 0.5*hh) )

def SimpleLogLik(data, template, Sn, df, tdi='XYZ'):

    if (tdi=='XYZ'):
        Xd = data[0]
        Yd = data[1]
        Zd = data[2]

        Xt = template[0]
        Yt = template[1]
        Zt = template[2]

        SNX = np.sum( np.real(Xd*np.conjugate(Xt))/Sn )
        SNY = np.sum( np.real(Yd*np.conjugate(Yt))/Sn )
        SNZ = np.sum( np.real(Zd*np.conjugate(Zt))/Sn )

        # print ('SN = ', 4.0*df*SNX, 4.0*df*SNY, 4.0*df*SNZ)

        XX = np.sum( np.abs(Xt)**2/Sn )
        YY = np.sum( np.abs(Yt)**2/Sn )
        ZZ = np.sum( np.abs(Zt)**2/Sn )

        # print ('hh = ', 4.0*df*XX, 4.0*df*YY, 4.0*df*ZZ)
        llX = 4.0*df*(SNX - 0.5*XX)
        llY = 4.0*df*(SNY - 0.5*YY)
        llZ = 4.0*df*(SNZ - 0.5*ZZ)

        return (llX, llY, llZ)

    else: ### I presume it is A, E
        Ad = data[0]
        Ed = data[1]

        At = template[0]
        Et = template[1]

        SNA = np.sum( np.real(Ad*np.conjugate(At))/Sn )
        SNE = np.sum( np.real(Ed*np.conjugate(Et))/Sn )

        # print ('SN = ', 4.0*df*SNA, 4.0*df*SNE)

        AA = np.sum( np.abs(At)**2/Sn )
        EE = np.sum( np.abs(Et)**2/Sn )

        # print ('hh:', 4.0*df*AA, 4.0*df*EE)

        llA = 4.0*df*(SNA - 0.5*AA)
        llE = 4.0*df*(SNE - 0.5*EE)

        return (llA, llE)


def loglikeStas(sig, par):

    wvfp = MBHBtmplCoarse(par)
    frqs = sig["freq"]
    frqs_com = np.unique(np.concatenate((frqs, wvfp['freq'])))
    ind_rm = []
    for i in range(1, len(frqs_com)):
        if (frqs_com[i] - frqs_com[i-1] < 1.e-10):
            ind_rm.append(i)
    frqs_com = np.delete(frqs_com, ind_rm)

    SA = tdi.noisepsd_AE(frqs_com, model='Proposal', includewd=None)

    splS = spline(sig['freq'], sig["phase"])
    splW = spline(wvfp['freq'], wvfp["phase"])

    SN = []
    hh = []
    for i in range(2):
        splAs_r = spline(sig['freq'], np.real(sig["ampl"][i, :]))
        splAs_i = spline(sig['freq'], np.imag(sig["ampl"][i, :]))

        splAw_r = spline(wvfp['freq'], np.real(wvfp["ampl"][i, :]))
        splAw_i = spline(wvfp['freq'], np.imag(wvfp["ampl"][i, :]))

        SN.append((splAs_r(frqs_com)+1.0j*splAs_i(frqs_com))*(splAw_r(frqs_com)-1.0j*splAw_i(frqs_com))/SA)
        hh.append((splAw_r(frqs_com)**2 + splAw_i(frqs_com)**2)/SA )

    dph = splS(frqs_com) - splW(frqs_com)
    sdph = np.sin(dph)
    cdph = np.cos(dph)
    # splA = spline(frqs_com, np.real(SN[0])*np.cos(dph) - np.imag(SN[0])*np.sin(dph) )
    # splE = spline(frqs_com, np.real(SN[1])*np.cos(dph) - np.imag(SN[1])*np.sin(dph) )
    splA = spline(frqs_com, np.real(SN[0])*cdph - np.imag(SN[0])*sdph )
    splE = spline(frqs_com, np.real(SN[1])*cdph - np.imag(SN[1])*sdph )

    splAA = spline(frqs_com, hh[0])
    splEE = spline(frqs_com, hh[1])

    SNA = 4.0*splA.integral(frqs_com[0], frqs_com[-1])
    SNE = 4.0*splE.integral(frqs_com[0], frqs_com[-1])

    hhA = 4.0*splAA.integral(frqs_com[0], frqs_com[-1])
    hhE = 4.0*splEE.integral(frqs_com[0], frqs_com[-1])

    return(SNA + SNE - 0.5*(hhA+hhE))

DATAPATH = "/home/stefan/LDC/Sangria/data"
sangria_fn = DATAPATH+"/LDC2_sangria_training_v1.h5"
tdi_ts, tdi_descr = hdfio.load_array(sangria_fn, name="obs/tdi")
dt = int(1/(tdi_descr["sampling_frequency"]))

# Build timeseries and frequencyseries object for X,Y,Z
tdi_ts = xr.Dataset(dict([(k,TimeSeries(tdi_ts[k], dt=dt)) for k in ["X", "Y", "Z"]]))
tdi_fs = xr.Dataset(dict([(k,tdi_ts[k].ts.fft(win=window)) for k in ["X", "Y", "Z"]]))



default_units = {'EclipticLatitude':'rad','EclipticLongitude':'rad',
         'PolarAngleOfSpin1':'rad','PolarAngleOfSpin2':'rad',
         'Spin1': '1','Spin2':'1',
         'Mass1':'Msun','Mass2':'Msun',
         'CoalescenceTime': 's','PhaseAtCoalescence':'rad',
         'InitialPolarAngleL':'rad','InitialAzimuthalAngleL':'rad',
         'Cadence': 's','Redshift': '1','Distance': 'Gpc',
         'ObservationDuration':'s'}

mbhb, units = hdfio.load_array(sangria_fn, name="sky/mbhb/cat")
print(units)
if not units:
    units = default_units
config = hdfio.load_config(sangria_fn, name="obs/config")
print(config)
s_index = 0
pMBHB = dict(zip(mbhb.dtype.names, mbhb[s_index]))
units = default_units
#for k,v in pMBHB.items():
#   print(k)
# #  pMBHB[k] *= u.Unit(units[k])

t_max = float(tdi_ts["X"].t[-1]+tdi_ts["X"].attrs["dt"])

pars = GetParams(pMBHB)
[Mc, q, tc, chi1, chi2, logDL, ci, sb, lam, psi, phi0]

sig = MBHBtmplCoarse(pMBHB)