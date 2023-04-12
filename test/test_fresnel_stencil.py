import numpy as np
import numpy.lib.recfunctions as rfn
from scipy.interpolate import InterpolatedUnivariateSpline as spline

from ldc.common.series import FrequencySeries, TDI
import lisaconstants.constants as LC
import pyFDresponse as FD_Resp
# import GenerateFD_SignalTDIs as GenTDIFD
import lisabeta.lisa.lisa as lisa
import lisabeta.lisa.lisa as lll


AU = LC.au
LL = 2.5e9
Reff = 2./np.pi*AU
Omega0 = 2.*np.pi/LC.ASTRONOMICAL_YEAR
armt = LL/LC.c
MTSUN_SI = LC.SUN_SCHWARZSCHILD_RADIUS / (2 * LC.c)


class FrequencyArray(np.ndarray):
    """
    Class to manage array in frequency domain based numpy array class
    """
    def __new__(subtype, data, dtype=None, copy=False, df=None, kmin=None):
        """
        ...
        @param data is ... [required]
        @param dtype is ... [default: None]
        @param copy is ... [default: None]
        @param df is ... [default: None]
        @param kmin is ... [default: None]
        @return
        """
        # make sure we are working with an array, copy the data if requested,
        # then transform the array to our new subclass
        subarr = np.array(data, dtype=dtype, copy=copy)
        subarr = subarr.view(subtype)

        # get df and kmin preferentially from the initialization,
        # then from the data object, otherwise set to None
        subarr.df = df if df is not None else getattr(data, 'df',  None)
        subarr.kmin = int(kmin) if kmin is not None else getattr(data, 'kmin', 0)

        return subarr

    def __array_wrap__(self, out_arr, context=None):
        """
        ...
        @param out_arr is ... [required]
        @param context is ... [default: None]
        @return
        """
        out_arr.df, out_arr.kmin = self.df, self.kmin

        return np.ndarray.__array_wrap__(self, out_arr, context)

    def __getitem__(self, key):
        """
        ...
        @param key is ... [required]
        @return
        """
        return self.view(np.ndarray)[key]

    def __getslice__(self, i, j):
        """
        ...
        @param i is ... [required]
        @param j is ... [required]
        @return
        """
        return self.view(np.ndarray)[i:j]

    # def __array_finalize__(self,obj):
    #    if obj is None: return
    #    self.df   = getattr(obj,'df',  None)
    #    self.kmin = getattr(obj,'kmin',None)

    def __repr__(self):
        """
        ...
        @return
        """
        if self.df is not None:
            return 'Frequency array (f0=%s,df=%s): %s' % (self.kmin * self.df, self.df, self)
        else:
            return str(self)

    def __add__(self, other):
        """
        Combine two FrequencyArrays into a longer one by adding intermediate zeros if necessary
        @param other is ... [required]
        @return
        """
        if isinstance(other, FrequencyArray) and (self.df == other.df):
            beg = min(self.kmin, other.kmin)
            end = max(self.kmin + len(self), other.kmin + len(other))

            ret = np.zeros(end-beg, dtype=np.find_common_type([self.dtype, other.dtype], []))

            ret[(self.kmin - beg):(self.kmin - beg + len(self))] = self
            ret[(other.kmin - beg):(other.kmin - beg + len(other))] += other

            return FrequencyArray(ret, kmin=beg, df=self.df)

        # fall back to simple arrays (may waste memory)
        return np.ndarray.__add__(self, other)

    def __sub__(self, other):
        """
        ...
        same behavior as __add__: TO DO -- consider restricting the result to the extend of the first array
        @param other is ... [required]
        @return
        """
        if isinstance(other, FrequencyArray) and (self.df == other.df):
            beg = min(self.kmin, other.kmin)
            end = max(self.kmin + len(self), other.kmin + len(other))

            ret = np.zeros(end-beg, dtype=np.find_common_type([self.dtype, other.dtype], []))

            ret[(self.kmin - beg):(self.kmin - beg + len(self))] = self
            ret[(other.kmin - beg):(other.kmin - beg + len(other))] -= other

            return FrequencyArray(ret, kmin=beg, df=self.df)

        # fall back to simple arrays (may waste memory)
        return np.ndarray.__sub__(self, other)

    def rsub(self,other):
        """
        Restrict the result to the extent of the first array (useful, e.g., for logL over frequency-limited data)
        @param other is ... [required]
        @return
        """
        if isinstance(other,FrequencyArray) and (self.df == other.df):
            if other.kmin >= self.kmin + len(self) or self.kmin >= other.kmin + len(other):
                return self
            else:
                beg = max(self.kmin,other.kmin)
                end = min(self.kmin + len(self),other.kmin + len(other))

                ret = np.array(self,copy=True,dtype=np.find_common_type([self.dtype,other.dtype],[]))
                ret[(beg - self.kmin):(end - self.kmin)] -= other[(beg - other.kmin):(end - other.kmin)]

                return FrequencyArray(ret,kmin=self.kmin,df=self.df)

        return np.ndarray.__sub__(self,other)

    def __iadd__(self,other):
        """
        The inplace add and sub will work only if the second array is contained in the first one
        also there may be problems with upcasting
        @param other is ... [required]
        @return
        """
        if isinstance(other,FrequencyArray) and (self.df == other.df):
            if (self.kmin <= other.kmin) and (self.kmin + len(self) >= other.kmin + len(other)):
                np.ndarray.__iadd__(self[(other.kmin - self.kmin):(other.kmin - self.kmin + len(other))],other[:])
                return self

        # fall back to simple arrays
        np.ndarray.__iadd__(self,other)
        return self

    def __isub__(self,other):
        """
        ...
        @param other is ... [required]
        @return
        """
        if isinstance(other,FrequencyArray) and (self.df == other.df):
            if (self.kmin <= other.kmin) and (self.kmin + len(self) >= other.kmin + len(other)):
                np.ndarray.__isub__(self[(other.kmin - self.kmin):(other.kmin - self.kmin + len(other))],other[:])
                return self

        # fall back to simple arrays
        np.ndarray.__isub__(self,other)
        return self


    def __mul__(self,other):
        """
        In multiplication, we go for the intersection of arrays (not their union!)
        no intersection return a scalar 0
        @param other is ... [required]
        @return
        """
        if isinstance(other,FrequencyArray) and (self.df == other.df):
            beg = max(self.kmin,other.kmin)
            end = min(self.kmin + len(self),other.kmin + len(other))

            if beg >= end:
                return 0.0
            else:
                ret = np.array(self[(beg - self.kmin):(end - self.kmin)],copy=True,dtype=np.find_common_type([self.dtype,other.dtype],[]))
                ret *= other[(beg - other.kmin):(end - other.kmin)]

                return FrequencyArray(ret,kmin=beg,df=self.df)

        # fall back to simple arrays (may waste memory)
        return np.ndarray.__mul__(self,other)


    def __div__(self,other):
        """
        In division, it's OK if second array is larger, but not if it's smaller (which implies division by zero!)
        @param other is ... [required]
        @return
        """
        if isinstance(other,FrequencyArray) and (self.df == other.df):
            if (other.kmin > self.kmin) or (other.kmin + len(other) < self.kmin + len(self)):
                raise ZeroDivisionError
            else:
                ret = np.array(self,copy=True,dtype=np.find_common_type([self.dtype,other.dtype],[]))
                ret /= other[(self.kmin - other.kmin):(self.kmin - other.kmin + len(self))]

            return FrequencyArray(ret,kmin=self.kmin,df=self.df)

        # fall back to simple arrays
        return np.ndarray.__div__(self,other)


    @property
    def f(self):
        """
        Return the reference frequency array
        """
        return np.linspace(self.kmin * self.df,(self.kmin + len(self) - 1) * self.df,len(self))


    @property
    def fmin(self):
        """
        Return the minimum frequency
        """
        return self.kmin * self.df


    @property
    def fmax(self):
        """
        Return the maximal frequency
        """
        return (self.kmin + len(self)) * self.df


    def ifft(self,dt):
        """
        ...
        @param dt is ...
        """
        n = int(1.0/(dt*self.df))

        ret = np.zeros(int(n/2+1),dtype=self.dtype)
        ret[self.kmin:self.kmin+len(self)] = self[:]
        ret *= n                                        # normalization, ehm, found empirically

        return np.fft.irfft(ret)


    def restrict(self,other):
        """
        Restrict the array to the dimensions of the second, or to dimensions specified as (kmin,len)
        @param Other array
        """
        if isinstance(other,FrequencyArray) and (self.df == other.df):
            kmin, length = other.kmin, len(other)
        elif isinstance(other,(list,tuple)) and len(other) == 2:
            kmin, length = other
        else:
            raise TypeError

        # no need to restrict anything?
        if kmin == self.kmin and length == len(self):
            return other

        ret = FrequencyArray(np.zeros(length,dtype=self.dtype),kmin=kmin,df=self.df)

        beg = max(self.kmin,kmin)
        end = min(self.kmin + len(self),kmin + length)

        ret[(beg - kmin):(end - kmin)] = self[(beg - self.kmin):(end - self.kmin)]

        return ret


    def pad(self,leftpad=1,rightpad=1):
        """
        Pad the array on both sides
        @param leftpad is ...
        @param rightpad is ...
        """
        return self.restrict((self.kmin - int(leftpad)*len(self),int(1+leftpad+rightpad)*len(self)))


def funcestimateforfomorb2nd1(f):
    return (2.*np.pi*f*Reff/LC.c*Omega0**2)


def funcestimateforfomorb2nd2(f):
    return (2.*np.pi*f*Reff/LC.c*Omega0)**2


def funcestimateforfomorb2nd(f):
    return max(funcestimateforfomorb2nd1(f), funcestimateforfomorb2nd2(f))


def funcestimateforfomconst2nd1(f):
    return (4.*np.pi*f*armt*Omega0**2)


def funcestimateforfomconst2nd2(f):
    return (4.*np.pi*f**2*(armt)**2*Omega0**2)


def funcestimateforfomconst2nd3(f):
    return Omega0**2


def funcestimateforfomconst2nd(f):
    return max(funcestimateforfomconst2nd1(f),
               funcestimateforfomconst2nd2(f),
               funcestimateforfomconst2nd3(f))


def funcFOMApprox(key, m1, m2, fstart):
    M = m1 + m2
    q = max(m1/m2, m2/m1)
    nu = q/(1+q)**2
    Msec = M*MTSUN_SI

    def Tffunc(f):
        return 1./8*np.sqrt(5./(3*nu))*Msec*(np.pi*Msec*f)**(-11./6)

    if key == 'Psiorb':
        fom = 1./2*Tffunc(fstart)**2*funcestimateforfomorb2nd(fstart)
    elif key == 'Psiconst':
        fom = 1./2*Tffunc(fstart)**2*funcestimateforfomconst2nd(fstart)
    return fom


def UpSampleTDIs(wfTDI, Tobs, dt):
    # {{{
    w_fr = wfTDI['freq']
    w_ph = wfTDI['phase']
    # get rid of possibly huge constant in the phase before interpolating
    tfspline = spline(w_fr, 1/(2.*np.pi)*(w_ph-w_ph[0])).derivative()
    tfvec = tfspline(w_fr)

    fspl = spline(tfvec, w_fr)
    fend = fspl(Tobs)
    # print ("fend = ", fend)

    length = int(0.5*Tobs/dt) + 1
    df = 1.0/Tobs
    t0 = wfTDI['t0']

    # fmax = 0.5/dt
    freqs = np.arange(length)*df
    fbeg = max(freqs[0], wfTDI['freq'][0])
    fend = min(freqs[-1], fend)
    ibeg = int((fbeg - freqs[0])/df) + 1
    iend = int((fend - freqs[0])/df)
    # print (fbeg, fend)
    fs = freqs[ibeg:iend+1]
    ampspline = spline(wfTDI['freq'], wfTDI['amp'])
    phasespline = spline(wfTDI['freq'], wfTDI['phase'])
    phaseRdelayspline = spline(wfTDI['freq'], wfTDI['phaseRdelay'])
    amp = ampspline(fs)
    phase = phasespline(fs)
    phaseRdelay = phaseRdelayspline(fs)
    phasetimeshift = 2.*np.pi*t0*fs

    keytrs = ['transferL1', 'transferL2', 'transferL3']
    fastPart = amp * np.exp(1j*(phase + phaseRdelay + phasetimeshift))

    keytrs = ['transferL1', 'transferL2', 'transferL3']
    for _, ky in enumerate(keytrs):
        transferLRespline = spline(wfTDI['freq'], np.real(wfTDI[ky]))
        transferLImspline = spline(wfTDI['freq'], np.imag(wfTDI[ky]))
        transferLRe = transferLRespline(fs)
        transferLIm = transferLImspline(fs)
        if (ky == 'transferL1'):
            X = np.conjugate((transferLRe+1j*transferLIm) * fastPart)
        if (ky == 'transferL2'):
            Y = np.conjugate((transferLRe+1j*transferLIm) * fastPart)
        if (ky == 'transferL3'):
            Z = np.conjugate((transferLRe+1j*transferLIm) * fastPart)

    return (fs, X, Y, Z)


def ComputeSOBBH_FDTDI(p, Tobs, dt, buf=None, order=None):
    """
    This function determines the order of frensel stencil and computes SOBBH GW
    signal for a single source
    @param p parameter dictionary as reurned by hdf5 reading function
    @param Tobs observation time (sec)
    @param buf FrequencyArray where the computed signal will be placed
    @param order user defined order of stencil for FD response (forcing order)

    output
    frequency array
    buf FrequencyArray (3xNf) containing X, Y, Z in FD
    wfTDI debugging info dictionary containg Response functions, GW phase,
    amplitude, dopper delay
    """
    m1 = p['m1']
    m2 = p['m2']
    # z = p['z']
    DL = p['dist']

    bet = p['beta']
    lam = p['lambda']
    inc = p['inc']
    psi = p['psi']

    # TODO: HARDCODED
    # chi1 = p['chi1']
    # chi2 = p['chi2']
    a1 = 0
    a2 = 0
    # phi0 = p['phi0']
    phi0 = 0.7

    fstart = p['fstart']

    if (m2 > m1):
        m_tmp = m1
        m1 = m2
        m2 = m_tmp
    # q = m1/m2
    # eta = m1*m2/(m1+m2)**2
    fRef = fstart
    tRef = 0.0
    fny = 0.5/dt
    minf = 1.e-5
    maxf = 1.e-1
    # fend = fny
    df = 1.0/Tobs

    # fom = funcFOMApprox('Psiorb', m1, m2, fstart)
    # [fRef, trajdict, TDItag] = [fstart, FD_Resp.trajdict_MLDC, "TDIXYZ"]
    fRef = fstart
    tobs = 0.0
    t0 = 0.0

    print(f"Order fresnel stencil= {order}")
    assert order is not None

    wfTDI = FD_Resp.LISAGenerateTDI(phi0, fRef, m1, m2, a1, a2, DL, inc, lam, bet, psi, tobs=tobs,
                                    minf=max(minf, fstart), maxf=min(fny, maxf), t0=t0,
                                    settRefAtfRef=True, tRef=tRef, \
                                    trajdict=FD_Resp.trajdict_MLDC, TDItag='TDIXYZ',
                                    nptmin=500, order_fresnel_stencil=order)
    frS, X, Y, Z = UpSampleTDIs(wfTDI, Tobs, dt)

    # frS, X, Y, Z, wfTDI = GenTDIFD.SOBBH_LISAGenerateTDI(phi0, fRef, m1, m2, a1, a2, DL, inc, lam, bet,
    #                                                      psi, Tobs, minf=max(minf, fstart), maxf=fend, t0=t0,
    #                                                      settRefAtfRef=True, tRef=tRef,
    #                                                      trajdict=FD_Resp.trajdict_MLDC,
    #                                                      TDItag='TDIXYZ', nptmin=500,
    #                                                      order_fresnel_stencil=order)

    ### Following block raises Only order_fresnel_stencil=0 implemented for now.
    ### with non zero values of order_fresnel_stencil parameter
    # md = (2, 2)
    # tdisignal = lll.GenerateLISATDISignal_SOBH(p, order_fresnel_stencil=order)
    #
    # wfTDI = tdisignal['wftdi']
    # frq = tdisignal['tdi'][md]['freq']
    # frS = np.arange(0, frq.max(), df)
    # kmin = int(frS[0]/df)
    # tdifreqseries = lisa.EvaluateTDIFreqseries(tdisignal["tdi"], frS)
    #
    # A = FrequencySeries(np.conj(tdifreqseries[md]['chan1']), df=df, kmin=kmin, name="A")
    # E = FrequencySeries(np.conj(tdifreqseries[md]['chan2']), df=df, kmin=kmin, name="E")
    # T = FrequencySeries(np.conj(tdifreqseries[md]['chan3']), df=df, kmin=kmin, name="T")
    # tdi = TDI(dict(zip(["A", "E", "T"], [A, E, T])))
    # tdi.AET2XYZ()
    # X, Y, Z = tdi['X'], tdi['Y'], tdi['Z']

    # return (freqs, buf, wfTDI)
    return (frS, [X, Y, Z], wfTDI)
    # return(frS, [X, Y, Z], wfTDI)


def rename_as_lisabeta(params_i):
    """ Rename fields to match lisabeta parameter names
    """
    # Copy the parameters
    params = params_i.copy()

    # Spin projection
    if 'PolarAngleOfSpin1' in params.keys():
        params["Spin1"] = params['Spin1']*np.cos(params['PolarAngleOfSpin1'])
        params["Spin2"] = params['Spin2']*np.cos(params['PolarAngleOfSpin2'])
        k_param = ['PolarAngleOfSpin1', 'PolarAngleOfSpin2',
                   'ObservationDuration', 'Cadence', 'Redshift']

        for k in k_param:
            if k in params.keys():
                params.pop(k)

    # Parameters to match between LDC notation and lisabeta
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
    # Adding additionnal parameters
    dmatch['InitialPhase'] = 'phi'

    # Return only the relevant parameters
    new_params = dict()
    if isinstance(params, dict):
        for k, v in params.items():
            if k in dmatch.keys():
                new_params[dmatch[k]] = params[k]
    else:
        new_params = rfn.rename_fields(params, dmatch)
    return new_params


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from tdi_factory import TDIFactory

    facGB = TDIFactory(source_type="SBBH", approximant='IMRPhenomD',
                       dt=5, duration=60*60*24*365)
    print(facGB.param)
    ldc = facGB.ldc()
    ldc_fd = facGB.tofd(ldc)
    fast = facGB.fast()

    # instead of computing LDCwaveform use the one with routines above
    params = rename_as_lisabeta(facGB.get_default_param())
    print(params)

    order_fresnel_stencil1 = 1
    freqs1, XYZ1, wfTDI1 = ComputeSOBBH_FDTDI(params, facGB.t_max, facGB.dt,
                                              order=order_fresnel_stencil1)

    order_fresnel_stencil2 = 10
    freqs2, XYZ2, wfTDI2 = ComputeSOBBH_FDTDI(params, facGB.t_max, facGB.dt,
                                              order=order_fresnel_stencil2)

    np.save('ldc_fd_ofs.npy', ldc_fd)
    np.save('fast_ofs.npy', fast.to_array())
    np.save('freqs1_ofs1.npy', freqs1)
    np.save('XYZ1_ofs1.npy', XYZ1)
    np.save('freqs2_ofs10.npy', freqs2)
    np.save('XYZ2_ofs10.npy', XYZ2)

    plt.figure()
    plt.subplot(111)
    # plt.semilogy(ldc_fd.f, np.abs(ldc_fd), label='ldc', alpha=0.5)
    plt.semilogy(fast.f, np.abs(fast.X), label='fast', alpha=0.5)
    plt.plot(freqs1, np.abs(XYZ1[0]), label=f'fast (fresnel) ofs={order_fresnel_stencil1}', alpha=0.5)
    plt.plot(freqs2, np.abs(XYZ2[0]), label=f'fast (fresnel) ofs={order_fresnel_stencil2}', alpha=0.5)
    plt.ylabel("|X|")
    plt.legend()

    plt.figure(figsize=(15, 10))
    plt.subplot(311)
    # plt.plot(ldc_fd.f, ldc_fd.real, label='ldc')
    plt.plot(fast.f, fast.X.real, label='fast')
    plt.plot(freqs1, np.real(XYZ1[0]), label=f'fast (fresnel) ofs={order_fresnel_stencil1}')
    plt.plot(freqs2, np.real(XYZ2[0]), label=f'fast (fresnel) ofs={order_fresnel_stencil2}')
    plt.ylabel("real X")
    # plt.axis([0.012, 0.01202, None, None])
    # plt.axis([0.01208, 0.01210, None, None])
    plt.axis([0.01220, 0.01221, None, None])
    plt.legend()
    plt.subplot(312)
    # plt.plot(ldc_fd.f, ldc_fd.imag, label='ldc')
    plt.plot(fast.f, fast.X.imag, label='fast')
    plt.plot(freqs1, np.imag(XYZ1[0]), label=f'fast (fresnel) ofs={order_fresnel_stencil1}')
    plt.plot(freqs2, np.imag(XYZ2[0]), label=f'fast (fresnel) ofs={order_fresnel_stencil2}')
    plt.ylabel("imag X")
    # plt.axis([0.012, 0.01202, None, None])
    # plt.axis([0.01208, 0.01210, None, None])
    plt.axis([0.01220, 0.01221, None, None])
    plt.legend()
    plt.subplot(313)
    # plt.plot(ldc_fd.f, np.abs(ldc_fd), label='ldc')
    plt.plot(fast.f, np.abs(fast.X), label='fast')
    plt.plot(freqs1, np.abs(XYZ1[0]), label=f'fast (fresnel) ofs={order_fresnel_stencil1}')
    plt.plot(freqs2, np.abs(XYZ2[0]), label=f'fast (fresnel) ofs={order_fresnel_stencil2}')
    plt.ylabel("|X|")
    # plt.axis([0.012, 0.01202, None, None])
    # plt.axis([0.01208, 0.01210, None, None])
    plt.axis([0.01220, 0.01221, None, None])
    plt.legend()

    plt.show()
