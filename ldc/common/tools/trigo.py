import numpy as np

def aziPolAngleL2PsiIncl(bet, lam, theL, phiL):
    """
    Convert polar and azimuthal angles of zS (typically orbital angular momentum L)
    to polarisation and inclination (see doc)
    - bet is the ecliptic latitude of the source in sky [rad]
    - lam is the ecliptic longitude of the source in sky [rad]
    - theL is the polar angle of zS [rad]
    - phiL is the azimuthal angle of zS [rad]
    """
    inc = np.arccos( - np.cos(theL)*np.sin(bet) - np.cos(bet)*np.sin(theL)*np.cos(lam - phiL) )
    down_psi = np.sin(theL)*np.sin(lam - phiL)
    up_psi = -np.sin(bet)*np.sin(theL)*np.cos(lam - phiL) + np.cos(theL)*np.cos(bet)
    psi = np.arctan2(up_psi, down_psi)
    
    return psi, inc
    
