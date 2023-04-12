import lisaconstants
import numpy as np
from numpy import pi, sqrt, cos, sin, exp

#pylint:disable=C0103
#pylint:disable=C0301

MTSUN_SI = lisaconstants.GM_SUN/lisaconstants.SPEED_OF_LIGHT**3
YRSID_SI = lisaconstants.SIDEREALYEAR_J2000DAY*24*60*60

def mchirpofm1m2(m1, m2):
    """ Compute chirp mass from m1 and m2.
    """
    return pow(m1*m2, 3./5) / pow(m1+m2, 1./5)

def newtonianfoft(Mchirp, t):
    """Newtonian estimate of the relation f(deltat) (for the 22 mode
    freq).

    Gives the starting geometric frequency for a given time to merger
    and chirp mass in Hz

    Args:
        Mchirp: chirp mass in solar masses
        t: time in years
    """
    if(t<=0.):
        return 0.
    return 1./np.pi * pow(Mchirp*MTSUN_SI, -5./8) * pow(256.*t*YRSID_SI/5, -3./8)

def newtoniantoffchirp(Mchirp, f):
    """ Newtonian estimate of the relation deltat(f) (for the 22 mode
    freq)

    Gives the time to merger from a starting frequency for a given
    chirp mass in years

    Args:
        Mchirp: input chirp mass in solar masses
        f: frequency in Hz
    """
    return 5./256 * pow(Mchirp*MTSUN_SI, -5./3) * pow(pi*f, -8./3) / YRSID_SI

def deltaMfPowerLaw(eta, Mf, acc):
    """ Used to build frequency grid.
    """
    return 3.8 * np.power(eta * acc, 1./4.) * Mf * np.power(np.pi*Mf, 5./12.)

def deltaMfLog(Mfmin, Mfmax, npt, Mf):
    """ Used to build frequency grid.
    """
    return Mf * (np.power(Mfmax/Mfmin, 1/(npt - 1)) - 1)


def spinWeightedSphericalHarmonic(s, l, m, theta, phi):
    """ Function reproducing XLALSpinWeightedSphericalHarmonic

    Currently only supports s=-2, l=2,3,4,5 modes
    """
    func = "SpinWeightedSphericalHarmonic"
    # Sanity checks
    if ( l < abs(s) ):
        raise ValueError('Error - %s: Invalid mode s=%d, l=%d, m=%d - require |s| <= l\n' % (func, s, l, m))
    if ( l < abs(m) ):
        raise ValueError('Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n' % (func, s, l, m))
    if not ( s == -2 ):
        raise ValueError('Error - %s: Invalid mode s=%d - only s=-2 implemented\n' % (func, s))

    fac = {
        # l=2
        (2,-2) : sqrt( 5.0 / ( 64.0 * pi ) ) * ( 1.0 - cos( theta ))*( 1.0 - cos( theta )),
        (2,-1) : sqrt( 5.0 / ( 16.0 * pi ) ) * sin( theta )*( 1.0 - cos( theta )),
        (2,0) : sqrt( 15.0 / ( 32.0 * pi ) ) * sin( theta )*sin( theta ),
        (2,1) : sqrt( 5.0 / ( 16.0 * pi ) ) * sin( theta )*( 1.0 + cos( theta )),
        (2,2) : sqrt( 5.0 / ( 64.0 * pi ) ) * ( 1.0 + cos( theta ))*( 1.0 + cos( theta )),
        # l=3
        (3,-3) : sqrt(21.0/(2.0*pi))*cos(theta/2.0)*pow(sin(theta/2.0),5.0),
        (3,-2) : sqrt(7.0/(4.0*pi))*(2.0 + 3.0*cos(theta))*pow(sin(theta/2.0),4.0),
        (3,-1) : sqrt(35.0/(2.0*pi))*(sin(theta) + 4.0*sin(2.0*theta) - 3.0*sin(3.0*theta))/32.0,
        (3,0) : (sqrt(105.0/(2.0*pi))*cos(theta)*pow(sin(theta),2.0))/4.0,
        (3,1) : -sqrt(35.0/(2.0*pi))*(sin(theta) - 4.0*sin(2.0*theta) - 3.0*sin(3.0*theta))/32.0,
        (3,2) : sqrt(7.0/(4.0*pi))*(-2.0 + 3.0*cos(theta))*pow(cos(theta/2.0),4.0),
        (3,3) : -sqrt(21.0/(2.0*pi))*pow(cos(theta/2.0),5.0)*sin(theta/2.0),
        # l=4
        (4,-4) : 3.0*sqrt(7.0/pi)*pow(cos(theta/2.0),2.0)*pow(sin(theta/2.0),6.0),
        (4,-3) : 3.0*sqrt(7.0/(2.0*pi))*cos(theta/2.0)*(1.0 + 2.0*cos(theta))*pow(sin(theta/2.0),5.0),
        (4,-2) : (3.0*(9.0 + 14.0*cos(theta) + 7.0*cos(2.0*theta))*pow(sin(theta/2.0),4.0))/(4.0*sqrt(pi)),
        (4,-1) : (3.0*(3.0*sin(theta) + 2.0*sin(2.0*theta) + 7.0*sin(3.0*theta) - 7.0*sin(4.0*theta)))/(32.0*sqrt(2.0*pi)),
        (4,0) : (3.0*sqrt(5.0/(2.0*pi))*(5.0 + 7.0*cos(2.0*theta))*pow(sin(theta),2.0))/16.0,
        (4,1) : (3.0*(3.0*sin(theta) - 2.0*sin(2.0*theta) + 7.0*sin(3.0*theta) + 7.0*sin(4.0*theta)))/(32.0*sqrt(2.0*pi)),
        (4,2) : (3.0*pow(cos(theta/2.0),4.0)*(9.0 - 14.0*cos(theta) + 7.0*cos(2.0*theta)))/(4.0*sqrt(pi)),
        (4,3) : -3.0*sqrt(7.0/(2.0*pi))*pow(cos(theta/2.0),5.0)*(-1.0 + 2.0*cos(theta))*sin(theta/2.0),
        (4,4) : 3.0*sqrt(7.0/pi)*pow(cos(theta/2.0),6.0)*pow(sin(theta/2.0),2.0),
        # l= 5
        (5,-5) : sqrt(330.0/pi)*pow(cos(theta/2.0),3.0)*pow(sin(theta/2.0),7.0),
        (5,-4) : sqrt(33.0/pi)*pow(cos(theta/2.0),2.0)*(2.0 + 5.0*cos(theta))*pow(sin(theta/2.0),6.0),
        (5,-3) : (sqrt(33.0/(2.0*pi))*cos(theta/2.0)*(17.0 + 24.0*cos(theta) + 15.0*cos(2.0*theta))*pow(sin(theta/2.0),5.0))/4.0,
        (5,-2) : (sqrt(11.0/pi)*(32.0 + 57.0*cos(theta) + 36.0*cos(2.0*theta) + 15.0*cos(3.0*theta))*pow(sin(theta/2.0),4.0))/8.0,
        (5,-1) : (sqrt(77.0/pi)*(2.0*sin(theta) + 8.0*sin(2.0*theta) + 3.0*sin(3.0*theta) + 12.0*sin(4.0*theta) - 15.0*sin(5.0*theta)))/256.0,
        (5,0) : (sqrt(1155.0/(2.0*pi))*(5.0*cos(theta) + 3.0*cos(3.0*theta))*pow(sin(theta),2.0))/32.0,
        (5,1) : sqrt(77.0/pi)*(-2.0*sin(theta) + 8.0*sin(2.0*theta) - 3.0*sin(3.0*theta) + 12.0*sin(4.0*theta) + 15.0*sin(5.0*theta))/256.0,
        (5,2) : sqrt(11.0/pi)*pow(cos(theta/2.0),4.0)*(-32.0 + 57.0*cos(theta) - 36.0*cos(2.0*theta) + 15.0*cos(3.0*theta))/8.0,
        (5,3) : -sqrt(33.0/(2.0*pi))*pow(cos(theta/2.0),5.0)*(17.0 - 24.0*cos(theta) + 15.0*cos(2.0*theta))*sin(theta/2.0)/4.0,
        (5,4) : sqrt(33.0/pi)*pow(cos(theta/2.0),6.0)*(-2.0 + 5.0*cos(theta))*pow(sin(theta/2.0),2.0),
        (5,5) : -sqrt(330.0/pi)*pow(cos(theta/2.0),7.0)*pow(sin(theta/2.0),3.0)
        }.get((l,m), None)
    if fac==None:
        raise ValueError('Error - %s: Invalid mode s=%d, l=%d, m=%d - require |m| <= l\n' % (func, s, l, m))

    # Result
    if m==0:
        return fac
    else:
        return fac * exp(1j*m*phi)
