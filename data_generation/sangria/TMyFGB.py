import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from scipy import interpolate
import sys
from tqdm import tqdm as tqdm
import corner
import matplotlib.pyplot as plt
import matplotlib as mpl

from ldc.lisa.orbits import Orbits
from ldc.common.series import TDI
import ldc.waveform.fastGB as fastGB
from ldc.lisa.projection import ProjectedStrain
from ldc.waveform.waveform import HpHc
from ldc.common.series import TimeSeries
from ldc.common.tools import window
import ldc.io.hdf5 as h5io
from ldc.common.tools import window
from ldc.lisa.noise import get_noise_model
import lisaconstants as constants

#import MCMC_multichain.MCMC_multichain as mcmc
#import MCMC_multichain.ptMCMC_multichain as ptmcmc

plt.style.use(['seaborn-ticks','seaborn-deep'])
mpl.rcParams.update({'font.size': 16})
plt.rcParams['axes.grid'] = True
plt.rcParams["figure.figsize"] = (11,5)

year = constants.SIDEREALYEAR_J2000DAY*24*60*60
MTsun = constants.GM_SUN/constants.SPEED_OF_LIGHT**3
pc = constants.PARSEC_METER 
clight = constants.SPEED_OF_LIGHT
AU = constants.ASTRONOMICAL_UNIT/clight


import pyFDresponse as FD_Resp

def sinc(x):
    return(np.sinc(x/np.pi))

class miniLW:
    def __init__(self, incl, bet, lam, psi, phi0, arm = 2.5e9):
        self.incl = incl     ### Inclination
        self.bet = bet       ### Ecliptic Latitude
        self.lam = lam       ### Ecliptic Longitude
        self.psi = psi       ### polarization
        self.armt = arm/clight  ### armlength in sec
        self.settRefAtfRef = False

        Y22 = FD_Resp.SpinWeightedSphericalHarmonic(-2, 2, 2, incl, phi0)
        Y2m2 = FD_Resp.SpinWeightedSphericalHarmonic(-2, 2, -2, incl, phi0)
        # self.Yfactorplus = 0.5 * (Y22 + np.conjugate(Y2m2))
        # self.Yfactorcross = -0.5j * (Y22 - np.conjugate(Y2m2))
        self.Yfactorplus = 0.5 * (Y22 + np.conjugate(Y2m2))
        self.Yfactorcross = 0.5j * (Y22 - np.conjugate(Y2m2))

    def ComputeProjections(self, tvec, wv="modes", Ap=None, Ac=None, verbose=False):

        #### We assume the analytic eccentric orbits.
        self.tfvec = tvec
        x, y, z, alpha, a = self.LISAmotion(tvec, verbose=verbose)
        n1 = np.array([x[1,:]-x[2,:], y[1,:]-y[2,:], z[1,:] - z[2,:]])/self.armt
        n2 = np.array([x[2,:]-x[0,:], y[2,:]-y[0,:], z[2,:] - z[0,:]])/self.armt
        n3 = np.array([x[0,:]-x[1,:], y[0,:]-y[1,:], z[0,:] - z[1,:]])/self.armt
        R0 = np.array([a*np.cos(alpha), a*np.sin(alpha), 0])
        self.k = 1.0*np.array([-np.cos(self.bet)*np.cos(self.lam), -np.cos(self.bet)*np.sin(self.lam), -np.sin(self.bet)])
        self.kR0 = np.dot(self.k, R0)
        self.kn1 = np.dot(self.k, n1)
        self.kn2 = np.dot(self.k, n2)
        self.kn3 = np.dot(self.k, n3)

        # print ("LW: k, t, kR0", self.k, tvec[-10:], self.kR0[-10:])


        q1 = np.array([x[0,:]-a*np.cos(alpha), y[0, :]-a*np.sin(alpha), z[0, :]])
        q2 = np.array([x[1,:]-a*np.cos(alpha), y[1, :]-a*np.sin(alpha), z[1, :]])
        q3 = np.array([x[2,:]-a*np.cos(alpha), y[2, :]-a*np.sin(alpha), z[2, :]])

        self.kq1 = np.dot(self.k, q1)
        self.kq2 = np.dot(self.k, q2)
        self.kq3 = np.dot(self.k, q3)

        # print ("LW angles:, lam, bet, psi, incl", self.lam, self.bet, self.psi, self.incl)
        Pl1, Cr1 = self.PolarizationTensor(self.lam, self.bet, self.psi)
        # Pij = 1.0*(Pl1*self.Yfactorplus + Cr1*self.Yfactorcross) ### SB: assumes that h_ij = P_ij A e^{i\Psi} in FD

        Pij = np.conjugate(Pl1*self.Yfactorplus + Cr1*self.Yfactorcross) ### SB because h_ij = P_ij A e^{-i\Psi}
        # print ("polariz", Pl1, Cr1)
        # print ("Pij", Pij)
        if (wv == "hphc"):
            # print ("incl = ", self.incl/np.pi)
            # Ap = (1.0 + (np.cos(self.incl))**2)
            Ap = -(1.0 + (np.cos(self.incl))**2)  ### New conventions
            Ac = -2.0*np.cos(self.incl)
            Pij = Ap*Pl1 - 1.0j*Ac*Cr1 ### SB: here I assume that h_{ij} = Re( Pij * A * exp(i\Phi(t)-phi0) )
            # Pij = Ap*Pl1 + 1.0j*Ac*Cr1 ### SB: here I assume that h_{ij} = Re( Pij * A * exp(-i\Phi(t)) )
            ## This assumes we are working in the time domain and we need to take real part of $p_{ij} A e^{i \phi}$

        N = len(tvec)
        self.A1 = np.zeros(N, dtype='complex128')
        self.A2 = np.zeros(N, dtype='complex128')
        self.A3 = np.zeros(N, dtype='complex128')
        for i in range(N):
            if (wv == 'emri'):
                 Pij = Ap[i]*Pl1 + Ac[i]*Cr1
            self.A1[i] = np.dot(np.dot(n1[:, i], Pij), n1[:, i])
            self.A2[i] = np.dot(np.dot(n2[:, i], Pij), n2[:, i])
            self.A3[i] = np.dot(np.dot(n3[:, i], Pij), n3[:, i])


    def ComputeResponse(self, om, method="full"):
        ## om - single freq, or array of freq: 2*pi*f_gw
        ## method - ["low", "full"], "low" - using the leading order only, "full" - goes beond the leading order

        self.om = om
        self.omL = om*self.armt
        phR = om*self.kR0

        y123 = self.ComputeYslr(send=1, link=2, rec=3, lsign=1.)
        y231 = self.ComputeYslr(send=2, link=3, rec=1, lsign=1.)
        y1m32 = self.ComputeYslr(send=1, link=3, rec=2, lsign=-1.)
        y3m21 = self.ComputeYslr(send=3, link=2, rec=1, lsign=-1.)
        y2m13 = self.ComputeYslr(send=2, link=1, rec=3, lsign=-1.)
        y312 = self.ComputeYslr(send=3, link=1, rec=2, lsign=1.)


        omL = self.omL
        RX = -2.0j*np.sin(omL)*(y231 + (y1m32 - y123)*np.exp(-1.0j*omL) - y3m21)*np.exp(-1.0j*omL)
        RY = -2.0j*np.sin(omL)*(y312 + (y2m13 - y231)*np.exp(-1.0j*omL) - y1m32)*np.exp(-1.0j*omL)
        RZ = -2.0j*np.sin(omL)*(y123 + (y3m21 - y312)*np.exp(-1.0j*omL) - y2m13)*np.exp(-1.0j*omL)

        # print ("LW Y123", y123[-10:])

        ## NOTE the sign of delay phase phR !!!
        return(-phR, RX, RY, RZ)

    def TryResponse(self, om):

        fctr = np.exp(-1.j*om*self.kR0)
        
        y123 = self.ComputeYslr(send=1, link=2, rec=3, lsign=1.)
        y231 = self.ComputeYslr(send=2, link=3, rec=1, lsign=1.)
        y1m32 = self.ComputeYslr(send=1, link=3, rec=2, lsign=-1.)
        y3m21 = self.ComputeYslr(send=3, link=2, rec=1, lsign=-1.)
        y2m13 = self.ComputeYslr(send=2, link=1, rec=3, lsign=-1.)
        y312 = self.ComputeYslr(send=3, link=1, rec=2, lsign=1.)

        omL = self.omL
        Rx = -2.j*(y231 + (y1m32 - y123)*np.exp(-1.0j*omL) - y3m21)*np.exp(-1.0j*om*self.kR0)/self.omL
        self.kR1 = self.kR0 + self.kq1
        self.kR2 = self.kR0 + self.kq2
        self.kR3 = self.kR0 + self.kq3

        # args = 0.5*self.omL*(1.0 - self.kn3)
        # y231 = -self.A3*sinc(args)*np.exp(1.j*(args - self.om*self.kR1))
        
        # args = 0.5*self.omL*(1.0 + self.kn2)
        # y3m21 = -self.A2*sinc(args)*np.exp(1.j*(args - self.om*self.kR1))

        # Rx = -2.j*(y231 - y3m21)*np.exp(-1.0j*om*self.kR0)/self.omL
        

        return (Rx)



    def ComputeTDTDI(self, tm, Amp, phi0, om, omdot, RX, RY, RZ, src="GB"):
        """ So far only for GBs, for other sources I need phase in time domain either numerically or analytic """
        if (src != "GB"):
            raise NotImplementedError

        #### TODO Check if we need to interpolate alreadyprecomputed projections
        xi = tm - self.kR0
        #print ("stas:", np.shape(xi), np.shape(tm), np.shape(self.kR0))

        phase = -phi0 + om*xi + 0.5*omdot*xi**2
        # phase = phi0 + om*xi + 0.5*omdot*xi**2
        #print (np.shape(phase))

        hfast = Amp*np.exp(1.0j*phase)
        #print (np.shape(hfast))

        X = np.real(hfast*RX)
        Y = np.real(hfast*RY)
        Z = np.real(hfast*RZ)
        #print (np.shape(X), np.shape(Y))

        return(X, Y, Z)


    def PolarizationTensor(self, lam, bet, psi):
       sl = np.sin(lam)
       cl = np.cos(lam)
       sb = np.sin(bet)
       cb = np.cos(bet)
       # u = 1.0*np.array([sb*cl, sb*sl, -cb])
       # v = 1.0*np.array([sl, -cl, 0.0])
       v = -1.0*np.array([sb*cl, sb*sl, -cb])
       u = 1.0*np.array([sl, -cl, 0.0])
       Ep = np.zeros((3,3))
       Ec = np.zeros((3,3))
       for i in range(3):
           for j in range(3):
               Ep[i, j] = u[i]*u[j] - v[i]*v[j]   ### Polarization basis in SSB
               #Ep[i, j] = -u[i]*u[j] + v[i]*v[j]
               Ec[i, j] = u[i]*v[j] + u[j]*v[i]

       #psi = 0.5*np.pi - psi
       cp = np.cos(2.0*psi)
       sp = np.sin(2.0*psi)

       # plus = -cp*Ep + sp*Ec
       # cross = -sp*Ep - cp*Ec
       plus = cp*Ep + sp*Ec
       cross = -sp*Ep + cp*Ec   ### polarization basis in the source frame

       return(plus, cross)


    def LISAmotion(self, tm, verbose=False):
        lamb = 0.0
        kappa = 0.0

        N = len(tm)

        a = AU
        e = self.armt/(2.0*np.sqrt(3.0)*a)
        if verbose:
            print ("orbital eccentricity", e)
        nn = np.array([1,2,3])
        # year = YRSID_SI

        Beta = (nn-1)*2.0*np.pi/3.0 + lamb
        alpha = 2.0*np.pi*tm/year + kappa

        x = np.zeros((3, N))
        y = np.zeros((3, N))
        z = np.zeros((3, N))
        for i in range(3):
            x[i, :] = a*np.cos(alpha) + a*e*(np.sin(alpha)*np.cos(alpha)*np.sin(Beta[i]) - (1.0 + (np.sin(alpha))**2)*np.cos(Beta[i]))
            y[i, :] = a*np.sin(alpha) + a*e*(np.sin(alpha)*np.cos(alpha)*np.cos(Beta[i]) - (1.0 + (np.cos(alpha))**2)*np.sin(Beta[i]))
            z[i, :] = -np.sqrt(3.0)*a*e*np.cos(alpha - Beta[i])

        return (x, y, z, alpha, a)

    
    def CompareYslr(self, send, link, rec, lsign=1):

            self.kR1 = self.kR0 + self.kq1
            self.kR2 = self.kR0 + self.kq2
            self.kR3 = self.kR0 + self.kq3

            if (link == 1):
                kn = lsign*self.kn1
                A = self.A1
            elif (link == 2):
                kn = lsign*self.kn2
                A = self.A2
            elif (link == 3):
                kn = lsign*self.kn3
                A = self.A3
            else:
                print ("should never be here")
                raise ValueError
        
            if (rec == 1):
                kR = self.kR1
            elif (rec == 2):
                kR = self.kR2
            elif (rec == 3):
                kR = self.kR3
            else:
                print ("should never be here")
                raise ValueError
            args = 0.5*self.omL*(1.0 - kn)


            # y = -0.5j*self.omL*A*sinc(args)*np.exp(-1.0j*(args + self.om*kR))
            y = -0.25*A*sinc(args)*np.exp(1.0j*(args - self.om*kR))
            # y = args
            print ("A=", A)

            return (y)



    def ComputeYslr(self, send, link, rec, lsign=1):
        ### Does not include fast varying part of the GW signal (only response part)
        if (link == 1):
            kn = lsign*self.kn1
            A = self.A1
        elif (link == 2):
            kn = lsign*self.kn2
            A = self.A2
        elif (link == 3):
            kn = lsign*self.kn3
            A = self.A3
        else:
            print ("should never be here")
            raise ValueError
        if (rec == 1):
            kq = self.kq1
        elif (rec == 2):
            kq = self.kq2
        elif (rec == 3):
            kq = self.kq3
        else:
            print ("should never be here")
            raise ValueError
        args = 0.5*self.omL*(1.0 - kn)

        # print ("A", A)
        # print ("args", args)
        # print ("kq", kq)
        y = -0.5j*self.omL*A*sinc(args)*np.exp(-1.0j*(args + self.om*kq))
        # y = -0.5j*self.omL*A*sinc(args)*np.exp(1.0j*(args - self.om*kq)) ### FIXME
        

        # print ("LW", -0.5j*self.omL, A, sinc(args), np.exp(-1.0j*(args + self.om*kq)) )
        ### FIXME
        #y = 0.5j*omL*A*sinc(args)*np.exp(-1.0j*(args + om*kq))

        #y = y*hf*np.exp(1.j*om*kR0)

        return(y)
    






class MyFGB:
    ''' This is replica of the C-code FastGB. Used for debugging purposes'''
    
    def __init__(self, params, Tobs, dt, orbit):
        '''
        @params params dictionary of parameters
        @Tobs observation time in sec
        @dt cadence
        @orbit orbit object describing LISA's motion
        '''

        self.orbit = orbit
        self.dt = dt
        self.Tobs = Tobs
        self.orbit = orbit

        self.f0 = params["Frequency"]
        self.fdot = params["FrequencyDerivative"]
        self.theta = 0.5*np.pi - params['EclipticLatitude']
        self.phi = params['EclipticLongitude']
        self.ampl = params['Amplitude']
        self.incl = params['Inclination']
        self.psi = params['Polarization']
        self.phi0 = -params['InitialPhase']
        self.q0 = self.f0*self.Tobs

        cosiota= np.cos(self.incl)
        
        cosps = np.cos(2.*self.psi)
        sinps = np.sin(2.*self.psi)

        Aplus = self.ampl*( 1.+ cosiota*cosiota)
        Across = -2.0*self.ampl*cosiota

        self.DP = Aplus*cosps - 1.0j*Across*sinps
        self.DC = -Aplus*sinps - 1.0j*Across*cosps

        self.u = np.zeros(3)
        self.v = np.zeros(3)
        self.k = np.zeros(3)

        sinth = np.sin(self.theta)
        costh = np.cos(self.theta)
        cosph = np.cos(self.phi)
        sinph = np.sin(self.phi)

        self.u[0] = costh*cosph
        self.u[1] = costh*sinph  
        self.u[2] = -sinth
        self.v[0] = sinph        
        self.v[1] = -cosph
        self.k[0] = -sinth*cosph  
        self.k[1] = -sinth*sinph
        self.k[2] = -costh
        
        self.eplus = np.zeros((3,3))
        self.ecross = np.zeros((3,3))

        for ii in range(3):
            for jj in range(3):
                self.eplus[ii, jj] = self.v[ii]*self.v[jj] - self.u[ii]*self.u[jj]
                self.ecross[ii, jj] = self.u[ii]*self.v[jj] + self.v[ii]*self.u[jj]
        
        # print ('DP = ', self.DP, 'DC =', self.DC)
        # print('eplus', self.eplus)
        # print('ecross', self.ecross)
        print ('k = ', self.k)

        
    def ConstructSlowPart(self, N=512):
        
        '''
        @params N number of points in slow part
        '''
        tm = np.linspace(0, self.Tobs, num=N, endpoint=False)
        # tm = np.linspace(0, self.Tobs, num=N)
        self.tm = tm
        
        arm_length =  self.orbit.arm_length
        init_rotation = self.orbit.initial_rotation
        init_position = self.orbit.initial_position

        P1 = self.orbit.compute_position(1, tm)
        P2 = self.orbit.compute_position(2, tm)
        P3 = self.orbit.compute_position(3, tm)
        ### Pi = [3 x Nt] - coordinates vs time

        self.r12 = (P2 - P1)/arm_length ## [3xNt]
        self.r13 = (P3 - P1)/arm_length
        self.r23 = (P3 - P2)/arm_length

        self.r21 = -self.r12
        self.r31 = -self.r13
        self.r32 = -self.r23

        kdotr12 = np.dot(self.k, self.r12) ### should be size Nt
        kdotr13 = np.dot(self.k, self.r13)
        kdotr23 = np.dot(self.k, self.r23)
        
        kdotr21 = -kdotr12
        kdotr32 = -kdotr23
        kdotr31 = -kdotr13

    
        kdotP1 = np.dot(self.k, P1)/clight
        kdotP2 = np.dot(self.k, P2)/clight
        kdotP3 = np.dot(self.k, P3)/clight

        xi =  np.array([tm - kdotP1, tm - kdotP2, tm - kdotP3])

        # print ('tm = ', tm)
        # print ('xi = ', xi)
        # print ('P1', P1/clight, '\n', 'P2', P2/clight, '\n', 'P3', P3/clight)
        # print (f'kdotr12 = {kdotr12}, \n, kdotr13 = {kdotr13}, \n, kdotr23 = {kdotr23}, \n ')
        # print (f'r12 = {self.r12}, \n, r13 = {self.r13}, \n, r23 = {self.r23}, \n ')

        self.Nt = len(tm)

        fi = np.zeros((3, self.Nt))

        for ii in range(3):
            fi[ii, :] = self.f0 + self.fdot*xi[ii, :]

        
        self.fstar = clight/(arm_length*2*np.pi)
	    #Ratio of true frequency to transfer frequency
        fonfs = fi/self.fstar
        

        ### compute transfer f-n

        f0 = self.f0
        self.q = np.rint(self.f0 * self.Tobs) ### index of nearest Fourier bin
        phi0 = self.phi0

        df = 2.*np.pi*(self.q/self.Tobs)

        print ('fonfs = ', fonfs)

        # print ('kdotr12 = ', kdotr12)
        # print ('kdotr23 = ', kdotr23)
        # print ('kdotr31 = ', kdotr31)


        arg12 = 0.5*fonfs[0,:] * (1 + kdotr12)
        arg23 = 0.5*fonfs[1,:] * (1 + kdotr23)
        arg31 = 0.5*fonfs[2,:] * (1 + kdotr31)

        arg13 = 0.5*fonfs[0,:] * (1 + kdotr13)
        arg21 = 0.5*fonfs[1,:] * (1 + kdotr21)
        arg32 = 0.5*fonfs[2,:] * (1 + kdotr32)

        # print ('arg12 = ', arg12)
        # print ('arg23 = ', arg23)
        # print ('arg31 = ', arg31)
        # print ('arg13 = ', arg13)
        # print ('arg21 = ', arg21)
        # print ('arg32 = ', arg32)

        # xi =  np.array([tm - kdotP1, tm - kdotP2, tm - kdotP3])
        # xi = tm - self.kR0
        # om = 2.0*np.pi*(f0 + fdot*tm_LW)
        # phase = -phi0 + om*xi + 0.5*omdot*xi**2

        arg2_1 = 2.0*np.pi*f0*xi[0] + phi0 - df*tm + np.pi*self.fdot*(xi[0]**2)
        arg2_2 = 2.0*np.pi*f0*xi[1] + phi0 - df*tm + np.pi*self.fdot*(xi[1]**2)
        arg2_3 = 2.0*np.pi*f0*xi[2] + phi0 - df*tm + np.pi*self.fdot*(xi[2]**2)


        om = 2.0*np.pi*f0
        self.argS_1 =  phi0 + (om - df)*tm + np.pi*self.fdot*(xi[0]**2)
        self.argS_2 =  phi0 + (om - df)*tm + np.pi*self.fdot*(xi[1]**2)
        self.argS_3 =  phi0 + (om - df)*tm + np.pi*self.fdot*(xi[2]**2)



        # print ("arg2_1", arg2_1)
        # print ("arg2_2", arg2_2)
        # print ("arg2_3", arg2_3)

        Dp12 = np.zeros(self.Nt)
        Dp23 = np.zeros(self.Nt)
        Dp31 = np.zeros(self.Nt)
        Dc12 = np.zeros(self.Nt)
        Dc23 = np.zeros(self.Nt)
        Dc31 = np.zeros(self.Nt)
        for ii in range(self.Nt):
            Dp12[ii] = np.dot(np.dot(self.eplus, self.r12[:,ii]), self.r12[:, ii])  ### sould be dimension Nt
            Dp23[ii] = np.dot(np.dot(self.eplus, self.r23[:,ii]), self.r23[:, ii])  ### sould be dimension Nt
            Dp31[ii] = np.dot(np.dot(self.eplus, self.r31[:,ii]), self.r31[:, ii])  ### sould be dimension Nt

            Dc12[ii] = np.dot(np.dot(self.ecross, self.r12[:,ii]), self.r12[:, ii])  ### sould be dimension Nt
            Dc23[ii] = np.dot(np.dot(self.ecross, self.r23[:,ii]), self.r23[:, ii])  ### sould be dimension Nt
            Dc31[ii] = np.dot(np.dot(self.ecross, self.r31[:,ii]), self.r31[:, ii])  ### sould be dimension Nt

        print ('check', np.shape(Dp12), np.shape(Dc12))
            

        # wfm->kdotr[0][1] += wfm->k[i]*wfm->r12[i];
		# wfm->kdotr[0][2] += wfm->k[i]*wfm->r13[i];
		# wfm->kdotr[1][2] += wfm->k[i]*wfm->r23[i];
        

    
        ### D21=D12, D32=D23, D13=D31

        
        self.y12 = 0.25*np.sin(arg12)/arg12 * np.exp(1.j*(arg12 + arg2_1)) * ( Dp12*self.DP + Dc12*self.DC )
        self.y23 = 0.25*np.sin(arg23)/arg23 * np.exp(1.j*(arg23 + arg2_2)) * ( Dp23*self.DP + Dc23*self.DC )
        self.y31 = 0.25*np.sin(arg31)/arg31 * np.exp(1.j*(arg31 + arg2_3)) * ( Dp31*self.DP + Dc31*self.DC )


        self.y21 = 0.25*np.sin(arg21)/arg21 * np.exp(1.j*(arg21 + arg2_2)) * ( Dp12*self.DP + Dc12*self.DC )
        self.y32 = 0.25*np.sin(arg32)/arg32 * np.exp(1.j*(arg32 + arg2_3)) * ( Dp23*self.DP + Dc23*self.DC )
        self.y13 = 0.25*np.sin(arg13)/arg13 * np.exp(1.j*(arg13 + arg2_1)) * ( Dp31*self.DP + Dc31*self.DC )


        # print ('12',  Dp12*self.DP + Dc12*self.DC )
        # print ('23',  Dp23*self.DP + Dc23*self.DC )
        # print ('31',  Dp31*self.DP + Dc31*self.DC )
        # print ('sinc12', 0.25*np.sin(arg12)/arg12)
        # print ('sinc23', 0.25*np.sin(arg23)/arg23)
        # print ('sinc31', 0.25*np.sin(arg31)/arg31)
        # print ('sinc21', 0.25*np.sin(arg21)/arg21)
        # print ('sinc32', 0.25*np.sin(arg32)/arg32)
        # print ('sinc13', 0.25*np.sin(arg13)/arg13)

        # dphaza_1 = 2.0*np.pi*f0*xi[0] + phi0 - df*tm + np.pi*self.fdot*(xi[0]**2)

        # self.G12 =  0.25*np.sin(arg12)/arg12 * np.exp(1.j*(arg12 -  2.0*np.pi*f0*kdotP1)) * \
                # ( Dp12*self.DP + Dc12*self.DC )/self.ampl
        self.G12 = 0.25*np.sin(arg12)/arg12 * np.exp(1.j*(arg12 -2.0*np.pi*f0*kdotP1))*( Dp12*self.DP + Dc12*self.DC )/self.ampl
        # self.G12 = arg12

        print ('Ampl = ', ( Dp12*self.DP + Dc12*self.DC )/self.ampl)


        # self.G12 = 0.25*np.sin(arg12)/arg12 * np.exp(-1.j*(arg12 + om*kdotP1)) * ( Dp12*self.DP + Dc12*self.DC )
        # self.G23 = 0.25*np.sin(arg23)/arg23 * np.exp(-1.j*(arg23 + om*kdotP2)) * ( Dp23*self.DP + Dc23*self.DC )
        # self.G31 = 0.25*np.sin(arg31)/arg31 * np.exp(-1.j*(arg31 + om*kdotP3)) * ( Dp31*self.DP + Dc31*self.DC )


        # self.G21 = 0.25*np.sin(arg21)/arg21 * np.exp(-1.j*(arg21 + om*kdotP2)) * ( Dp12*self.DP + Dc12*self.DC )
        # self.G32 = 0.25*np.sin(arg32)/arg32 * np.exp(-1.j*(arg32 + om*kdotP3)) * ( Dp23*self.DP + Dc23*self.DC )
        # self.G13 = 0.25*np.sin(arg13)/arg13 * np.exp(-1.j*(arg13 + om*kdotP1)) * ( Dp31*self.DP + Dc31*self.DC )

        A12 = ( Dp12*self.DP + Dc12*self.DC )
        A23 = ( Dp23*self.DP + Dc23*self.DC )
        A31 = ( Dp31*self.DP + Dc31*self.DC )

        self.G12 = 0.25*np.sin(arg12)/arg12 * np.exp(-1.j*(arg12 + om*kdotP1 - self.argS_1)) * A12
        self.G23 = 0.25*np.sin(arg23)/arg23 * np.exp(-1.j*(arg23 + om*kdotP2 - self.argS_2)) * A23
        self.G31 = 0.25*np.sin(arg31)/arg31 * np.exp(-1.j*(arg31 + om*kdotP3 - self.argS_3)) * A31


        self.G21 = 0.25*np.sin(arg21)/arg21 * np.exp(-1.j*(arg21 + om*kdotP2 - self.argS_2)) * A12
        self.G32 = 0.25*np.sin(arg32)/arg32 * np.exp(-1.j*(arg32 + om*kdotP3 - self.argS_3)) * A23
        self.G13 = 0.25*np.sin(arg13)/arg13 * np.exp(-1.j*(arg13 + om*kdotP1 - self.argS_1)) * A31

        # y = -0.5j*self.omL*A*sinc(args)*np.exp(-1.0j*(args + self.om*kq))
        # args = 0.5*self.omL*(1.0 - kn)
        # arg12 = 0.5*fonfs[0,:] * (1 + kdotr12)
        # arg2_1 = 2.0*np.pi*f0*xi[0] + phi0 - df*tm + np.pi*self.fdot*(xi[0]**2)  -> om*k.Ri




        # arg1 = 0.5*wfm->fonfs[i]*(1. + wfm->kdotr[i][j]);

        # arg2 =  PI*2*f0*wfm->xi[i] + phi0 - df*t;

        # sinc = 0.25*sin(arg1)/arg1;

        # tran1r = aevol*(wfm->dplus[i][j]*wfm->DPr + wfm->dcross[i][j]*wfm->DCr);
		# tran1i = aevol*(wfm->dplus[i][j]*wfm->DPi + wfm->dcross[i][j]*wfm->DCi);

        # tran2r = cos(arg1 + arg2);
		# tran2i = sin(arg1 + arg2);

        # wfm->TR[i][j] = sinc*(tran1r*tran2r - tran1i*tran2i);
		# wfm->TI[i][j] = sinc*(tran1r*tran2i + tran1i*tran2r);
    def ComputeXYZslow_FD2(self):

        f = self.f0 + self.fdot*self.tm
        omL = f/self.fstar
        SomL = np.sin(omL)
        
        fctr =np.exp(-1.j*omL)

        fctr2 = 4.0*omL*SomL*fctr/self.ampl


        self.Xsl =  fctr2*(self.G21 - self.G31 + (self.G12 - self.G13)*fctr)        
        self.Ysl =  fctr2*(self.G32 - self.G12 + (self.G23 - self.G21)*fctr)
        self.Zsl =  fctr2*(self.G13 - self.G23 + (self.G31 - self.G32)*fctr)

        Xf_slow = np.fft.fft(self.Xsl)
        Yf_slow = np.fft.fft(self.Ysl)
        Zf_slow = np.fft.fft(self.Zsl)

        # plt.plot(Xf_slow.real)
        # plt.show()

        Xf_slow = self.ampl*Xf_slow
        Yf_slow = self.ampl*Yf_slow
        Zf_slow = self.ampl*Zf_slow

        self.Xtry = 4.0*(self.G21 - self.G31 + (self.G12 - self.G13)*fctr)/self.ampl


        return (Xf_slow, Yf_slow, Zf_slow)


        
    def ComputeXYZslow_TD(self):

        f = self.f0 + self.fdot*self.tm

        print ('f', np.shape(f))
        print ('y', np.shape(self.y12))



        self.Xslow = (self.y12 - self.y13)*np.exp(-3j*f/self.fstar) + (self.y21 - self.y31)*np.exp(-2j*f/self.fstar)+\
             (self.y13 - self.y12)*np.exp(-1j*f/self.fstar) + self.y31 - self.y21


        self.Yslow = (self.y23 - self.y21)*np.exp(-3j*f/self.fstar) + (self.y32 - self.y12)*np.exp(-2j*f/self.fstar)+\
             (self.y21 - self.y23)*np.exp(-1j*f/self.fstar) + self.y12 - self.y32

        
        self.Zslow = (self.y31 - self.y32)*np.exp(-3j*f/self.fstar) + (self.y13 - self.y23)*np.exp(-2j*f/self.fstar)+\
             (self.y32 - self.y31)*np.exp(-1j*f/self.fstar) + self.y23 - self.y13

        

        self.Xslow = 2.0*(f/self.fstar)*self.Xslow
        self.Yslow = 2.0*(f/self.fstar)*self.Yslow
        self.Zslow = 2.0*(f/self.fstar)*self.Zslow

        # self.Xtry = 4.0*(self.G21 - self.G31 + (self.G12 - self.G13)*np.exp(-1.j*f/self.fstar))*np.exp(-1.j*f/self.fstar)/self.ampl
        self.Xtry = 4.0*(self.G21 - self.G31 + (self.G12 - self.G13)*np.exp(-1.j*f/self.fstar))/self.ampl
        # self.Xtry = 4.0*(self.G12 - self.G13 + (self.G21 - self.G31)*np.exp(-1.j*f/self.fstar))/self.ampl
        # self.Xtry = 4.0*(self.G21 - self.G31)*np.exp(-1.j*f/self.fstar)/self.ampl
        # self.Xtry = 4.0*(self.G21 - self.G31)/self.ampl
        # self.Xtry = 4.0*(self.G12 - self.G13)/self.ampl


    def ComputeXYZslow_FD(self):

        y12_fd = np.fft.fft(self.y12)
        y13_fd = np.fft.fft(self.y13)
        y21_fd = np.fft.fft(self.y21)
        y31_fd = np.fft.fft(self.y31)
        y23_fd = np.fft.fft(self.y23)
        y32_fd = np.fft.fft(self.y32)

        M = self.Nt
        f =  (self.q + np.arange(M) - M/2)/self.Tobs

        self.Xslow = (y12_fd - y13_fd)*np.exp(-3j*f/self.fstar) + (y21_fd - y31_fd)*np.exp(-2j*f/self.fstar)+\
             (y13_fd - y12_fd)*np.exp(-1j*f/self.fstar) + y31_fd - y21_fd




    def ComputeXYZ_FD2(self):

        Xf_slow, Yf_slow, Zf_slow = self.ComputeXYZslow_FD2()

        M = len(Xf_slow)

        Xf = np.zeros(M, dtype=np.complex128)
        Yf = np.zeros(M, dtype=np.complex128)
        Zf = np.zeros(M, dtype=np.complex128)

        split = int(M/2)
        Xf[:split] =  Xf_slow[split:]
        Xf[split:] =  Xf_slow[:split]

        freq = (self.q + np.arange(M) - M/2)/self.Tobs

        return (freq, Xf, Yf, Zf)




    def ComputeXYZ_FD(self):

        self.ComputeXYZslow_TD()

        # Delf = (self.q - self.Nt/4)/self.Tobs
        Delf = (self.q - 78)/self.Tobs
        shift = np.exp(-1j*2.*np.pi*Delf*self.tm)
        shift = 1.0

        Xf_slow = np.fft.fft(shift*self.Xslow)
        Yf_slow = np.fft.fft(shift*self.Yslow)
        Zf_slow = np.fft.fft(shift*self.Zslow)

        print (np.shape(Xf_slow))

        phiSL = np.pi/2.0 - self.f0/self.fstar
        fctr = np.exp(1.j*phiSL)
        # fctr = 1.0

        M = len(Xf_slow)

        ### Assuming Nt = even number: 
        ### Xf_slow[0] - zero freq. Zf_slow[-1] - nyquist freq.

        Xf = np.zeros(M, dtype=np.complex128)
        Yf = np.zeros(M, dtype=np.complex128)
        Zf = np.zeros(M, dtype=np.complex128)

        split = int(M/2)
        Xf[:split] = fctr* Xf_slow[split:]
        Xf[split:] = fctr* Xf_slow[:split]

    
        # print ('szs', np.shape(Xf), np.shape(Xf_slow))
        # Xf[0:M-1] = fctr * np.flip(np.conjugate(Xf_slow[1:]))
        # Xf[M-1:] = fctr*Xf_slow

        # Yf[0:M] = fctr * np.flip(np.conjugate(Yf_slow))
        # Yf[1:M+1] = fctr*Yf_slow
        
        # Zf[0:M] = fctr * np.flip(np.conjugate(Zf_slow))
        # Zf[1:M+1] = fctr*Zf_slow

        freq = (self.q + np.arange(M) - M/2)/self.Tobs
        # freq = (self.q + np.arange(M))/self.Tobs

        print ('check fr', self.q,  self.f0, self.q/self.Tobs)
        print (freq[0], freq[int(M/2)], freq[int(M/2)+1])

        return (freq, Xf, Yf, Zf)

        
      
    
def GetParams(p):
    bet = p['EclipticLatitude']
    lam = p['EclipticLongitude']
    Amp = p["Amplitude"]
    f0 = p["Frequency"]
    fdot = p["FrequencyDerivative"]
    iota = p["Inclination"]
    psi = p["Polarization"]
    phi0 = p["InitialPhase"]

    return (bet, lam, Amp, f0, fdot, iota, psi, phi0)
   

def main():

    pGB = dict({"Frequency":0.0109966533684,
                "FrequencyDerivative": 2.95291919174e-14, 
                "EclipticLatitude":-0.150495923831, 
                "EclipticLongitude":4.60685356883, 
                "Amplitude":2.60186425363e-23,
                "Inclination":0.389613740033,
                "Polarization":0.605754423063,
                "InitialPhase":4.54588971721})

    dt = 5.0
    Tobs = 31536000.0
    # Tobs = 31474406.250
    #Nt = 512
    Nt = 2048

    config = {"initial_position": 0, "initial_rotation": 0, 
              "nominal_arm_length": 2500000000, "orbit_type": 'analytic'}
    lisa_orbits = Orbits.type(config)

    mGB = MyFGB(pGB, Tobs, dt, lisa_orbits)


    print('building slow part')
    mGB.ConstructSlowPart()
    print('done')


    # print (mGB.eplus)
    epl = mGB.eplus
    r12 = mGB.r12
    tm = mGB.tm
    Nt =  len(tm)
    # print ('tm = ', tm)
    # np.dot(np.dot(self.eplus, self.r12).T, self.r12)
    # print (Nt, np.shape(r12))
    # plt.plot(tm, r12[0, :])
    # plt.plot(tm, r12[1, :])
    # plt.plot(tm, r12[2, :])
    # plt.show()
    # D1 = np.zeros(Nt)
    # for ii in range(3):
        # for jj in range(3):
            # D1 += epl[ii,jj]*r12[ii, :]*r12[jj,:]


    # D2 = np.zeros(Nt)
    # for ii in range(Nt):
        # D2[ii] = np.dot(np.dot(epl, r12[:,ii]), r12[:, ii])

    # print (np.allclose(D1,D2))

    #print ('y12:', mGB.y12)
    # print ('y13:', mGB.y13)
    # print ('y21:', mGB.y21)

    print (np.shape(mGB.y12))

    ### try to reproduce original FastGB

    #mfr, mX, mY, mZ = mGB.ComputeXYZ_FD()
    mfr, mX, mY, mZ = mGB.ComputeXYZ_FD2()


    ### FastGB....

    GB = fastGB.FastGB(delta_t=dt, T=Tobs, orbits=lisa_orbits) # in seconds
    oversample = 4
    Xs, Ys, Zs = GB.get_fd_tdixyz(template=pGB, oversample=oversample, simulator='synthlisa')

    print ('size of Xs', np.shape(Xs))
    #print (Xs.f.values)
    
    print (pGB["Frequency"], np.rint(pGB["Frequency"]*Tobs))
    # plt.plot(np.real(mX))
    # plt.show()

    ### ------ LW....--------------

   
    bet, lam, Amp, f0, fdot, iota, psi, phi0 = GetParams(pGB)

    psi_LW = psi

    print ('params', bet, lam, Amp, f0, fdot, iota, psi_LW, phi0)

    dt_lw = 15.0
    N_lw = Tobs/dt_lw
    tm_LW = np.linspace(0, Tobs, int(N_lw)+1)
    # tm_LW = np.linspace(0, Tobs, int(N_lw)+1)

    print ('lw N=', len(tm_LW), tm_LW)
    
    # LW_GB = LW.LW(iota, bet, lam, psi_LW, phi0)
    LW_GB = miniLW(iota, bet, lam, psi_LW, phi0)
    LW_GB.ComputeProjections(tm_LW, wv="hphc")

    om = 2.0*np.pi*(f0 + fdot*tm_LW)
    phR, RX, RY, RZ = LW_GB.ComputeResponse(om)

    om = 2.0*np.pi*f0
    omdot = 2.0*np.pi*fdot


    Xlw, Ylw, Zlw = LW_GB.ComputeTDTDI(tm_LW, Amp, phi0, om, omdot, RX, RY, RZ, src="GB")

    # Xlw_f = np.fft.rfft(Xlw)
    # freq = np.linspace(0, 0.5/dt_lw, len(Xlw_f))
    # print ('fr N', len(freq), Tobs*0.5/dt_lw+1)

    tdi_x_lw = TimeSeries(Xlw[:-1], dt=dt_lw, t0=0.0)
    # tdi_x_lw = TimeSeries(Xlw, dt=dt_lw, t0=0.0)

    # print (tdi_x_lw)

    tdi_xfd_lw = tdi_x_lw.ts.fft(win=window)
    
    Xlw_f = tdi_xfd_lw.values
    freq = tdi_xfd_lw.f.values

    print ('lw size', np.shape(Xlw_f))



    # G12 = mGB.G12
    # y12 = LW_GB.CompareYslr(send=2, link=3, rec=1, lsign=1.)

    RX_T = mGB.Xtry 
    RX_LW = LW_GB.TryResponse(om)

    # print (np.shape(G12), np.shape(y12))
    tm_fst = np.linspace(0, Tobs, num=Nt, endpoint=False)
    
    plt.figure(figsize=(14,8))


    print ('for this to work you need to switch off self.argS_i in self.G_ij')
    plt.subplot(211)
    plt.plot(tm_fst, RX_T.real, label='G12')
    plt.plot(tm_LW, RX_LW.real, '--', label='y12')
    # plt.plot(tm_fst, G12.real, label='G12')
    # plt.plot(tm_LW, y12.real, '--', label='y12')
    plt.legend(loc='lower right')
    
    plt.subplot(212)
    plt.plot(tm_fst, RX_T.imag, label='G12')
    plt.plot(tm_LW, RX_LW.imag, '--', label='y12')
    # plt.plot(tm_fst, G12.imag, label='G12')
    # plt.plot(tm_LW, y12.imag, '--', label='y12')
    plt.legend(loc='lower right')

    plt.show()


    # sys.exit(0)
    


    # ----------------------------------



    fctrX = 0.5*Tobs/Nt
    print ('1', dt*np.sqrt(Tobs), '2', fctrX)

    plt.figure(figsize=(15,10))
    plt.subplot(211)
    plt.plot(Xs.f.values, np.real(Xs.values), label='FastGB')
    # plt.plot(mfr, dt*mX.real*np.sqrt(Tobs), label='my')
    plt.plot(mfr, mX.real*fctrX, label='my')
    plt.plot(freq, Xlw_f.real, '--', label='LW')
    plt.axis([0.010995, 0.01100, -1.2e-15, 1.2e-15])
    plt.legend(loc='upper right')
    plt.ylabel("X (real part)")

    plt.subplot(212)
    plt.plot(Xs.f.values, np.imag(Xs.values), label='FastGB')
    # plt.plot(mfr, dt*mX.imag*np.sqrt(Tobs), label='my')
    plt.plot(mfr, fctrX*mX.imag, '--',  label='my')
    plt.axis([0.010995, 0.01100, -1.2e-15, 1.2e-15])
    plt.legend(loc='upper right')
    plt.ylabel("X (imag part)")


    plt.show()
    # sys.exit(0)


    ### ------------- Use Maude's method: ------------

    print ('dt = ', dt, dt_lw)
    print ('check', tdi_x_lw.t.values[:-1],  tdi_x_lw.t.values)
    t_min= 0.0
    # t_max = 31536000.0 +dt_lw
    # t_max = tdi_x_lw.t.values[-1] +dt_lw
    t_max = tm_LW[-1] + dt_lw
    # t_max = 31536000.0
    
    Proj = ProjectedStrain(lisa_orbits)
    GW = HpHc.type("debug", "GB", "TD_fdot")
    GW.set_param(pGB)
    yArm = Proj.arm_response(t_min, t_max, dt_lw, [GW])
    maude_tdi_x = TimeSeries(Proj.compute_tdi_x(tm_LW[:-1]), dt=dt_lw)
    # maude_tdi_x = TimeSeries(Proj.compute_tdi_x(tdi_x_lw.t.values[:-1]), dt=dt_lw)
    # maude_tdi_x = TimeSeries(Proj.compute_tdi_x(tdix_lw.t.values[:2097152]), dt=dt_lw)
    print ('maude time', t_max,  maude_tdi_x.t.values)


    # plt.plot(maude_tdi_x.t, maude_tdi_x.values, label='Maude')
    # plt.plot(tdi_x_lw.t,  tdi_x_lw.values, '--', label='LW')
    # plt.legend(loc='upper right')
    # plt.show()


    tdi_xfd_maude = maude_tdi_x.ts.fft(win=window)
    print ('to compare df:', Tobs,  1.0/Tobs, 1./tdi_x_lw.t.values[-1], 1./(tdi_x_lw.t.values[-1] + dt_lw))
    print ('df1 = ', np.diff(tdi_xfd_maude.f.values))
    print ('df2 = ', np.diff(Xs.f.values))
    
    print (np.shape(tdi_xfd_maude.f), np.shape(tdi_xfd_lw), np.shape(Xlw), np.shape(maude_tdi_x.values), np.shape(mX))
    
    plt.figure(figsize=(15,10))
    plt.subplot(211)
    plt.plot(Xs.f.values, np.real(Xs.values), label='FastGB')
    # plt.plot(mfr, dt*mX.real*np.sqrt(Tobs), label='my')
    plt.plot(mfr, mX.real*fctrX, label='my')
    plt.plot(freq, Xlw_f.real, '--', label='LW')
    plt.plot(tdi_xfd_maude.f, tdi_xfd_maude.values.real, '-.', label='Maude')
    plt.axis([0.010995, 0.01100, -1.2e-15, 1.2e-15])
    plt.legend(loc='upper right')
    plt.ylabel("X (real part)")

    
    plt.subplot(212)
    plt.plot(Xs.f.values, np.imag(Xs.values), label='FastGB')
    # plt.plot(mfr, dt*mX.real*np.sqrt(Tobs), label='my')
    plt.plot(mfr, mX.imag*fctrX, label='my')
    plt.plot(freq, Xlw_f.imag, '--', label='LW')
    plt.plot(tdi_xfd_maude.f, tdi_xfd_maude.values.imag, '-.', label='Maude')
    plt.axis([0.010995, 0.01100, -1.2e-15, 1.2e-15])
    plt.legend(loc='upper right')
    plt.ylabel("X (real part)")




    plt.show()

    sys.exit(0)

    
    dGBh5file = '/Users/stas/Projects/LISA/Sangria/Data/LDC2_sangria_gdb-tdi_v1_orig.h5'
    tdi  = h5io.load_array(dGBh5file)

    
    print ('sangria time', tdi[0]['X'][2:, 0])


    # fctr = 0.5/dt
    fctr = 1.0
    X = fctr*tdi[0]['X'][2:-1, 1]
    Y = fctr*tdi[0]['Y'][2:-1, 1]
    Z = fctr*tdi[0]['Z'][2:-1, 1]

    tdi_sangria_x = TimeSeries(X, dt=dt, t0=0.0)
    tdi_sangria_y = TimeSeries(Y, dt=dt, t0=0.0)
    tdi_sangria_z = TimeSeries(Z, dt=dt, t0=0.0)

    
    fctrXt = np.sqrt(Tobs/dt_lw)

    print ('check', fctrXt)
    Xlw, Ylw, Zlw = LW_GB.ComputeTDTDI(tm_LW, Amp, phi0, om, omdot, RX, RY, RZ, src="GB")

    plt.plot(tdi_sangria_x.t,  tdi_sangria_x.values, label='sangria')
    # plt.plot(tdi_x_lw.t,  fctrXt*tdi_x_lw.values, label='LW')
    plt.plot(tdi_x_lw.t,  Xlw, label='LW')
    # plt.plot(tdi_x_lw.t,  2.*Nt*tdi_x_lw.values/Tobs, label='LW')
    plt.legend()
    plt.show()


    # tdi_sangria_fd = TDI(dict({"X":tdi_sangria_x.ts.fft(),
    tdi_sangria_fd = TDI(dict({"X":tdi_sangria_x.ts.fft(win=window),
                           "Y":tdi_sangria_y.ts.fft(win=window),
                           "Z":tdi_sangria_z.ts.fft(win=window)}))

    subset_sangria = tdi_sangria_fd.sel(f=mfr, method="nearest")

    print (np.diff(mfr)[:3], 1.0/np.diff(mfr)[0])
    print (np.diff(tdi_sangria_fd.f.values)[:3], 1.0/np.diff(tdi_sangria_fd.f.values)[0])

    # sys.exit(0)

    plt.figure(figsize=(14,9))
    plt.subplot(211)
    
    plt.plot(Xs.f.values, np.real(Xs.values), label='FastGB')
    plt.plot(mfr, dt*mX.real*np.sqrt(Tobs), label='my')
    plt.plot(freq, Xlw_f.real,  label='LW')
    plt.plot(tdi_sangria_fd.f, tdi_sangria_fd["X"].real,  label='sangria', ls='--')
    # plt.plot(mfr, dt*mX.real*np.sqrt(Tobs)-subset_sangria.X.values.real, alpha=0.5, color='grey', label='myGB - sangria ')
    # plt.plot(Xs.f.values,  np.real(Xs.values)-subset_sangria.X.values.real, alpha=0.5, color='grey', label='myGB - sangria ')
    plt.axis([0.010995, 0.01100, -1.2e-15, 1.2e-15])
    plt.legend(loc='upper right')
    plt.ylabel("X (real part)")

    plt.subplot(212)
    plt.plot(Xs.f.values, np.imag(Xs.values), label='FastGB')
    plt.plot(mfr, dt*mX.imag*np.sqrt(Tobs), label='my')
    plt.plot(tdi_sangria_fd.f, tdi_sangria_fd["X"].imag,  label='sangria', ls='--')
    plt.plot(mfr, dt*mX.imag*np.sqrt(Tobs)-subset_sangria.X.values.imag, alpha=0.5, color='grey', label='myGB - sangria ')
    plt.axis([0.010995, 0.01100, -1.2e-15, 1.2e-15])
    plt.legend()
    plt.xlabel("Freq [Hz]")
    plt.ylabel("X (imag part)")

    plt.show()


    #plt.plot(mfr, mX.real)
    #plt.plot(mfr, mX.imag, '--')
    #plt.show()


    ### check the links... 

    


    print ('here')


if __name__ == "__main__":

    main()

        









