import numpy as np
import lisaconstants
from ldc.common import tools
import ldc.waveform.imrphenomD as imrphenomD

d = dict({'AzimuthalAngleOfSpin1': 4.9, #"radian"),
          'AzimuthalAngleOfSpin2': 5.1, #"radian"),
          'EclipticLatitude': 0.312414, #"radian"),
          'EclipticLongitude': -2.75291,# "radian"),
          'CoalescenceTime': 28086000.0,# 's'), 
          'Distance':  9.14450149011798,# 'Gpc'), 
          'InitialAzimuthalAngleL': 3.9,# 'radian'), 
          'InitialPolarAngleL': 2.3535, #'radian'), 
          'Mass1': 132628.202,# "SolarMass"),
          'Mass2': 30997.2481,# "SolarMass"),
          'PhaseAtCoalescence':  3.8, #'Radian'), 
          'PolarAngleOfSpin1': 0.0 ,#'Radian'),
          'PolarAngleOfSpin2': 0.0 ,#'Radian'),
          'Redshift': 1.2687,# 'dimensionless'),
          'Spin1': 0.9481998052314212, #'MassSquared'),
          'Spin2': 0.9871324769575264, #'MassSquared'),
          'Cadence': 5.,
          'ObservationDuration':95.})#, 's') })


phi0 = 0.           # Orbital phase at fRef (rad)
fRef = 0.           # Reference frequency (Hz)
m1_SI = d["Mass1"]*lisaconstants.SUN_MASS #  code needs masses in kg
m2_SI = d["Mass2"]*lisaconstants.SUN_MASS # Mass of companion 2 (kg)
a1 = np.cos(d["PolarAngleOfSpin1"])*d["Spin1"] # For PhenomD we will use projections
a2 = np.cos(d["PolarAngleOfSpin2"])*d["Spin2"] # Aligned-spin parameter of companion 2
distance = d['Distance']*1e3*1e6*lisaconstants.PARSEC_METER             # Distance of source (m)
pol, incl = tools.aziPolAngleL2PsiIncl(d['EclipticLatitude'],
                                       d["EclipticLongitude"],
                                       d['InitialPolarAngleL'],
                                       d['InitialAzimuthalAngleL'])
freq = np.arange(2e-4, 1e-1, 2e-8) # Frequencies (Hz) on which to evaluate the waveform 

imr_phenomD = imrphenomD.pyIMRPhenomDh22AmpPhase(freq)
frS, ampS, phS = imr_phenomD.get_fd_waveform(phi0, fRef, m1_SI, m2_SI, a1, a2, distance)

df = freq[1]-freq[0]

imr_phenomD = imrphenomD.pyIMRPhenomD()
fr, hpf2, hcf2 = imr_phenomD.get_fd_waveform(phi0, fRef, df, m1_SI, m2_SI, a1, a2,
                                             freq[0], freq[-1], distance, incl)


