"""Compute the stationary phase approximation method to go from frequency/time
domain to time/frequency domain using an analytical approximation which speeds
up the computation compared to a FFT. It is mainly used for the SOBHB waveform
that are slowly varying.
"""

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

#pylint:disable=C0103

def signal_inverse_spa(signal_fd, sign=1):
    """  Inverse Stationary Phase Approximation for a Fourier-domain signal.

    Allows for chirping or anti-chirping signal with the sign option.
    """
    f, amp_fd, phi_fd = signal_fd.T #Extracting signal in the frequency domain
    spline_phi_fd = InterpolatedUnivariateSpline(f, phi_fd,k=3) #Phase interpolation
    tf      = 1./(2*np.pi) *spline_phi_fd.derivative(n=1)(f) #First order derivative of the interpolation
    Tf2     = sign* 1./(2*np.pi)**2 * spline_phi_fd.derivative(n=2)(f) #Second order derivative
    amp_td = amp_fd / np.sqrt(2*np.pi*Tf2) #Amplitude from frequency to time domain
    phi_td = phi_fd - 2*np.pi*f*tf  + sign*np.pi/4  #Phase from frequency to time domain - Stas version
    # Introduce factor 2 in amplitude
    # to account for negative frequencies
    amp_td *= 2.
    return np.array([tf, amp_td, phi_td]).T


# TODO: commented because not used and "f" variable is not defined.
# def signal_spa(signal_td, sign=1):
#     """ Stationary Phase Approximation for a time-domain signal

#     Allows for chirping or anti-chirping signal with the sign option
#     """
#     t, amp, phi = signal_td.T #Extracting signal in the frequency domain
#     spline_phi_td = InterpolatedUnivariateSpline(t, phi, k=3) #Phase interpolation
#     omega   = spline_phi_td.derivative(n=1)(f) #First derivative
#     omegadot= spline_phi_td.derivative(n=2)(f) #Second derivative
#     f = -omega / (2*np.pi)
#     A_spa = amp * np.sqrt(2*np.pi / (-sign*omegadot)) #Amplitude from time to frequency domain
#     Psi_spa = -phi - 2*np.pi*f*t - sign*np.pi/4 #Phase from time to frequency domain
#     if sign==1: return np.array([f, A_spa, -Psi_spa]).T
#     else:       return np.array([f[::-1], A_spa[::-1], -Psi_spa[::-1]]).T
