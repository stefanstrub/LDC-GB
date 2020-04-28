import numpy as np

def info(self):
    EMRIunits = dict({"AzimuthalAngleOfSpin": 'Radian', 
                      "Distance": 'Gpc', 
                      "EclipticLatitude": "Radian", 
                      "EclipticLongitude": 'Radian', 
                      "InitialAlphaAngle": "Radian", 
                      "InitialAzimuthalOrbitalFrequency": "Hz", 
                      "InitialAzimuthalOrbitalPhase": "Radian", 
                      "InitialEccentricity": 'Unitless', 
                      "InitialTildeGamma": 'Radian', 
                      "LambdaAngle": 'Radian', 
                      "MassOfCompactObject": 'SolarMass', 
                      "MassOfSMBH": 'SolarMass', 
                      "PlungeTime": 'Seconds', 
                      "PolarAngleOfSpin": 'Radian', 
                      "SMBHspin": "MassSquared"})
    return EMRIunits

def check_param(self):
    """ Check parameters and their units
    """
    for k in list(self.info().keys()):
        assert k in self.pnames
    assert self.units["EclipticLatitude"].lower() in ["radian", "rad", "r"]
    assert self.units["EclipticLongitude"].lower() in ["radian", "rad", "r"]
    assert self.units["PolarAngleOfSpin"].lower() in ["radian", "rad", "r"]

def compute_hphc_td(self, t, source_parameters=None, approx_t=False):
    """ Return hp, hx for a time samples in t.

    Source parameters can be updated at the same time. 

    >>> GW = HpHcEMRI("my-emri", "EMRI", "AK")
    >>> hp,hc = GW.compute_hphc_td(np.arange(0,100,10), pEMRI)
    >>> print(hp[0:3], hc[0:3] )
    [-1.66415910e-23 -1.97139923e-23 -2.27346619e-23] [-2.70064032e-23 -2.75318104e-23 -2.79559684e-23]
    """
    if source_parameters != None:
        self.set_param(source_parameters)

    if (self.approximant == 'AK'):
        #Tobs = t[-1]+self.dt
        self.wv.__samples = len(t)
        self.wv.__inittime = int(t[0])
        self.wv.__deltat = int(t[1]-t[0])
        hpSSB, hcSSB = self.wv.waveforms(self.wv.__samples, self.wv.__deltat, self.wv.__inittime)
        tm = t
    else:
        raise NotImplementedError

    ### Interpolate
    if approx_t:
        self.t, self.hp, self.hc = tm, hpSSB, hcSSB
    else:
        self._interp_hphc(tm, hpSSB, hcSSB, t, kind="spline")
    return (self.hp, self.hc)

if __name__ == "__main__":

    import numpy as np
    import doctest
    from hphc import HpHcEMRI

    pEMRI = dict({"AzimuthalAngleOfSpin": 3.3, #'Radian'), 
                  "Distance": 4.066032714225506, #'Gpc'), 
                  "EclipticLatitude": -1.0865079083647613, #"Radian"), 
                  "EclipticLongitude": 5.45915569996,# 'Radian'), 
                  "InitialAlphaAngle": 5.144987211440423,# "Radian"), 
                  "InitialAzimuthalOrbitalFrequency": 0.0006865, #"Hz"), 
                  "InitialAzimuthalOrbitalPhase": 1.774141989722132,# "Radian"), 
                  "InitialEccentricity": 0.3365614512812651,# 'Unitless'), 
                  "InitialTildeGamma": 4.783451533651018,# 'Radian'), 
                  "LambdaAngle": 1.0976000000002033, #'Radian'), 
                  "MassOfCompactObject": 25.89699960277051,# 'SolarMass'), 
                  "MassOfSMBH": 1330605.9163026307,# 'SolarMass'), 
                  "PlungeTime": 41327816.5734, #'s'), 
                  "PolarAngleOfSpin": 2.3754411347141318, #'Radian'), 
                  "SMBHspin": 0.9671,# "MassSquared"),
                  #'Cadence': (5., 's'),
                  })

    doctest.testmod()
