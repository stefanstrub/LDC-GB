""" Compute the rigid adiabatic TDI. """

import numpy as np
import lisaconstants as constants
from ldc.common import tools

# pylint: disable=C0103
# pylint: disable=R0914

CLIGHT = constants.SPEED_OF_LIGHT
AU = constants.ASTRONOMICAL_UNIT


def polarization_tensor(lam, bet, psi):
    """ Return polarization basis in the source frame.
    """
    sl, cl = np.sin(lam), np.cos(lam)
    sb, cb = np.sin(bet), np.cos(bet)
    v = -1.0*np.array([sb*cl, sb*sl, -cb])
    u = 1.0*np.array([sl, -cl, 0.0])
    u.shape = (3,1)
    v.shape = (3,1)
    Ep = np.dot(u, u.T) - np.dot(v, v.T)
    Ec = np.dot(u, v.T) + np.dot(u, v.T)
    #psi = 0.5*np.pi - psi
    cp, sp = np.cos(2.0*psi), np.sin(2.0*psi)
    plus = cp*Ep + sp*Ec
    cross = -sp*Ep + cp*Ec # polarization basis in the source frame
    return plus, cross


def sinc(x):
    """ sin cardinal of x / pi.
    """
    return np.sinc(x/np.pi)

class LW:
    """ Long wavelength approximation.
    """

    def __init__(self, incl, bet, lam, psi, phi0, orbits=None, ref_at_tref=False):
        """ Initialization of LW container.
        """
        self.incl, self.psi = incl, psi
        self.bet, self.lam = bet, lam
        self.ref_at_tref = ref_at_tref
        self.phi0 = phi0

        if orbits is None:
            orbits = Orbits.type(dict({'orbit_type':'analytic', 'nominal_arm_length':2.5e9,
                                       "initial_position": 0, "initial_rotation": 0}))
        self.orbits = orbits
        self.armt = self.orbits.arm_length/CLIGHT  ### armlength in sec

        Y22 = tools.spinWeightedSphericalHarmonic(-2, 2, 2, incl, phi0)
        Y2m2 = tools.spinWeightedSphericalHarmonic(-2, 2, -2, incl, phi0)
        self.Yfactor_plus = 0.5 * (Y22 + np.conjugate(Y2m2))
        self.Yfactor_cross = 0.5j * (Y22 - np.conjugate(Y2m2))

    def get_td_tdixyz(self, tm, amp, f0, fdot, source_type="GB"):
        """ So far only for GBs, for other sources I need phase in time domain
        either numerically or analytic
        """

        if source_type!="GB":
            raise NotImplementedError

        kR0 = self.compute_projection(tm, wv="hphc", Ap=None, Ac=None)
        ph, RX, RY, RZ = self.compute_response(2*np.pi*(f0 + fdot*tm), kR0)
        om = 2.0*np.pi*f0
        omdot = 2.0*np.pi*fdot
        xi = tm - kR0
        phase = -self.phi0 + om*xi + 0.5*omdot*xi**2
        hfast = amp*np.exp(1.0j*phase)
        X = np.real(hfast*RX)
        Y = np.real(hfast*RY)
        Z = np.real(hfast*RZ)
        return X, Y, Z

    def compute_response(self, om, kR0):
        """ Return combined link response.
        """
        omL = om*self.armt
        phR = om*kR0

        y123 = self.compute_yslr(om, send=1, link=2, rec=3, lsign=1.)
        y231 = self.compute_yslr(om, send=2, link=3, rec=1, lsign=1.)
        y312 = self.compute_yslr(om, send=3, link=1, rec=2, lsign=1.)
        y1m32 = self.compute_yslr(om, send=1, link=3, rec=2, lsign=-1.)
        y3m21 = self.compute_yslr(om, send=3, link=2, rec=1, lsign=-1.)
        y2m13 = self.compute_yslr(om, send=2, link=1, rec=3, lsign=-1.)

        RX = -2.0j*np.sin(omL)*(y231 + (y1m32 - y123)*np.exp(-1.0j*omL) - y3m21)*np.exp(-1.0j*omL)
        RY = -2.0j*np.sin(omL)*(y312 + (y2m13 - y231)*np.exp(-1.0j*omL) - y1m32)*np.exp(-1.0j*omL)
        RZ = -2.0j*np.sin(omL)*(y123 + (y3m21 - y312)*np.exp(-1.0j*omL) - y2m13)*np.exp(-1.0j*omL)

        return -phR, RX, RY, RZ # NOTE the sign of delay phase phR !!!

    def compute_yslr(self, om, send, link, rec, lsign=1):
        """Does not include fast varying part of the GW signal (only response
        part)

        """
        kn = lsign*self.kns[link]
        A = self.As[link]
        kq = self.kqs[rec]
        omL = om*self.armt
        args = 0.5*omL*(1.0 - kn)
        y = -0.5j*omL*A*sinc(args)*np.exp(-1.0j*(args + om*kq))
        ### FIXME
        #y = 0.5j*omL*A*sinc(args)*np.exp(-1.0j*(args + om*kq))
        #y = y*hf*np.exp(1.j*om*kR0)
        return y

    def get_positions(self, tvec):
        """ Return spacecrafts position.
        """
        alpha = self.orbits.alpha(tvec)
        x1,y1,z1 = self.orbits.compute_position(1, tvec)
        x2,y2,z2 = self.orbits.compute_position(2, tvec)
        x3,y3,z3 = self.orbits.compute_position(3, tvec) # Pi = [3 x Nt] - coordinates vs time
        x = np.array([x1, x2, x3])/CLIGHT
        y = np.array([y1, y2, y3])/CLIGHT
        z = np.array([z1, z2, z3])/CLIGHT
        a = AU/CLIGHT  ## AU in sec.
        return x, y, z, a, alpha

    def compute_projection(self, tvec, wv="modes", Ap=None, Ac=None):
        """ Compute projection basis.

        wv can be: 'hphc', 'emri'
        We assume the analytic eccentric orbits.
        """
        self.t = tvec
        x, y, z, a, alpha = self.get_positions(tvec)
        ns = dict()
        ns[1] = np.array([x[1,:]-x[2,:], y[1,:]-y[2,:], z[1,:] - z[2,:]])/self.armt
        ns[2] = np.array([x[2,:]-x[0,:], y[2,:]-y[0,:], z[2,:] - z[0,:]])/self.armt
        ns[3] = np.array([x[0,:]-x[1,:], y[0,:]-y[1,:], z[0,:] - z[1,:]])/self.armt
        R0 = np.array([a*np.cos(alpha), a*np.sin(alpha), 0])
        k = 1.0*np.array([-np.cos(self.bet)*np.cos(self.lam),
                          -np.cos(self.bet)*np.sin(self.lam),
                          -np.sin(self.bet)])
        kR0 = np.dot(k, R0)
        qs = dict()
        self.kqs, self.kns = dict(), dict()
        for i in range(3):
            self.kns[i+1] = np.dot(k, ns[i+1])
            qs[i+1] = np.array([x[i,:]-a*np.cos(alpha), y[i, :]-a*np.sin(alpha), z[i, :]])
            self.kqs[i+1] = np.dot(k, qs[i+1])

        Pl1, Cr1 = polarization_tensor(self.lam, self.bet, self.psi)
        # Pij = 1.0*(Pl1*self.Yfactorplus + Cr1*self.Yfactorcross)
        Pij = np.conjugate(Pl1*self.Yfactor_plus + Cr1*self.Yfactor_cross)# SB because h_ij = P_ij A e^{-i\Psi}

        if wv == "hphc":
            # Ap = (1.0 + (np.cos(self.incl))**2)
            Ap = -(1.0 + (np.cos(self.incl))**2) # new conventions
            Ac = -2.0*np.cos(self.incl)
            Pij = Ap*Pl1 - 1.0j*Ac*Cr1 # SB: assume that h_{ij} = Re( Pij * A * exp(i\Phi(t)-phi0) )
            # Pij = Ap*Pl1 + 1.0j*Ac*Cr1 ### SB: assume that h_{ij} = Re( Pij * A * exp(-i\Phi(t)) )
            ## This assumes we are working in the time domain
            ## and we need to take real part of $p_{ij} A e^{i \phi}$

        N = len(tvec)
        self.As = dict()
        for s in range(1,4):
            self.As[s] = np.zeros(N, dtype='complex128')
            for i in range(N):
                if (wv == 'emri'):
                    Pij = Ap[i]*Pl1 + Ac[i]*Cr1
                self.As[s][i] = np.dot(np.dot(ns[s][:, i], Pij), ns[s][:, i])
        return kR0
