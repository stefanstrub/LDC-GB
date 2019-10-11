"""
Copyright (C) 2006-2011 Gaia Data Processing and Analysis Consortium

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
            THIS IS AN AUTOMATICALLY GENERATED FILE - DO NOT EDIT!

The file has been automatically generated from the contents of the
Gaia Parameter Database at the URL
       https://gaia.esac.esa.int/gpdb/
on 2019-09-30T15:18:18.

Please report any problems arising from the usage of this file to
the Gaia Librarian gaia-helpdesk@cosmos.esa.int
"""


class PdbParam:
    """
    Container class to enclose the contents of the Gaia Parameter
    Database at<br/>
     <code><a href="https://gaia.esac.esa.int/gpdb/">https://gaia.esac.esa.int/gpdb/</a></code><p>
    A hierarchy of nested classes below matches the parameter naming scheme
    detailed in <code><a href="https://dms.cosmos.esa.int/cs/cs/Open/357616">GAIA-JdB-007</a></code><br/>

    author  Gaia SOC, ESA/ESTEC
    version Live
    """

    DBVersion: str = "Live"




    class Nature:

        __A0VSTAR_BMINV: float  = -0.004 # [mag]

        @property
        def A0VSTAR_BMINV(self):
            r"""
            Johnson B minus V (B-V) magnitude of an unreddened A0V star

            #Source: M.S. Bessell, S. Murphy, 2012, 'Spectrophotometric Libraries, Revised Photonic Passbands, and Zero Points for UBVRI, Hipparcos, and Tycho Photometry', PASP, 124, 140<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__A0VSTAR_BMINV


        __A0VSTAR_CALIBRATIONFLUX_LAMBDA: float  = 3.62286e-11 # [W m^-2 nm^-1]

        @property
        def A0VSTAR_CALIBRATIONFLUX_LAMBDA(self):
            r"""
            Flux f_{0\lambda} f an unreddened A0V star with V = 0 mag at the wavelength \lambda_0

            #Basic : false
            #Scalar: true
            #Unit: [W m^-2 nm^-1]
            """

            return self.__A0VSTAR_CALIBRATIONFLUX_LAMBDA


        __A0VSTAR_CALIBRATIONFLUX_NU: float  = 3.65558e-23 # [W m^-2 Hz^-1]

        @property
        def A0VSTAR_CALIBRATIONFLUX_NU(self):
            r"""
            Flux f_{0\nu} of an unreddened A0V star with V = 0 mag at the wavelength \lambda_0 (see also J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1 and R.C. Bohlin, 2014, 'Hubble Space Telescope CALSPEC Flux Standards: Sirius (and Vega)', Astronomical Journal, Volume 147, 127; spectrum named alpha_lyr_mod_002.fits contained in CALSPEC Calibration Database, http://www.stsci.edu/hst/observatory/crds/calspec.html, last modified April 2015)

            #Source: P. Montegriffo, 26 January 2016, 'External calibration for Gaia DR1 integrated photometry', GAIA-C5-TN-OABO-PMN-009, issue 1, revision 0, depends on parameters :Nature:A0VStar_Spectrum, :Nature:Planck_Constant, :Nature:FilterTransmissionCurve_JohnsonCousinsV_002, and :Nature:A0VStar_VMagnitude<br/>
            #Basic : false
            #Scalar: true
            #Unit: [W m^-2 Hz^-1]
            """

            return self.__A0VSTAR_CALIBRATIONFLUX_NU


        __A0VSTAR_CALIBRATIONFLUX_NUMBEROFPHOTONS: float  = 100308492.2552 # [photons s^-1 m^-2 nm^-1]

        @property
        def A0VSTAR_CALIBRATIONFLUX_NUMBEROFPHOTONS(self):
            r"""
            Photon flux N_{0\lambda}(\lambda_0) of an unreddened A0V star with V = 0 mag at the wavelength \lambda_0 (see also J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1). Note that the parameter A0VStar_Spectrum_NumberOfPhotons refers to Pickles' star number 009, and not to the Kurucz Vega spectrum which has been used for flux normalisation/calibration and zero-point definition; the parameter A0VStar_Spectrum_NumberOfPhotons therefore does not precisely have A0VStar_CalibrationFlux photons s^-1 m^-2 nm^-1 at \lambda_0 = A0VStar_CalibrationWavelength

            #Basic : false
            #Scalar: true
            #Unit: [photons s^-1 m^-2 nm^-1]
            """

            return self.__A0VSTAR_CALIBRATIONFLUX_NUMBEROFPHOTONS


        __A0VSTAR_CALIBRATIONFUNCTION: str  = "Nature/A0VStar_CalibrationFunction_002.fits"

        @property
        def A0VSTAR_CALIBRATIONFUNCTION(self):
            r"""
            Calibration function S_V(\lambda). This function can, alternatively, be used to define the zero point of the Johnson V magnitude scale by imposing the requirement that, for any stellar photon flux density N_\lambda (in photons s^-1 m^-2 nm^-1 above the Earth's atmosphere) with V = 0 mag, the integral from 470 to 740 nm (the support interval of the Johnson V band) of N_\lambda times S_V(\lambda) equals N_0 photons s^-1 m^-2. The function S_V(\lambda) and the normalisation constant N_0 depend on the value of Planck's constant (parameter :Nature:Planck_Constant), on the definition of the shape of the Johnson V band (parameter :Nature:FilterTransmissionCurve_JohnsonCousinsV_002), on the monochromatic calibration flux f_{0\lambda} (or f_{0\nu}; parameters :Nature:A0VStar_CalibrationFlux_Lambda and :Nature:A0VStar_CalibrationFlux_Nu) at \lambda_0 (parameter :Nature:A0VStar_CalibrationWavelength), and on the spectrum f_{0\nu}(\lambda) of the general unreddened A0V star (parameter :Nature:A0VStar_Spectrum_Nu_002). First column: wavelength \lambda (in nm; from 470.0 to 740.0). Second column: S_V(\lambda)

            #Source: J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1, Appendix D<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__A0VSTAR_CALIBRATIONFUNCTION


        __A0VSTAR_CALIBRATIONFUNCTION_NORMALISATION: float  = 8630065822.2737 # [photons s^-1 m^-2]

        @property
        def A0VSTAR_CALIBRATIONFUNCTION_NORMALISATION(self):
            r"""
            Calibration (function) normalisation constant N_0. This constant can be used to define the zero point of the Johnson V magnitude scale by imposing the requirement that, for any stellar photon flux density N_\lambda (in photons s^-1 m^-2 nm^-1 above the Earth's atmosphere) with V = 0 mag, the integral from 470 to 740 nm (the support interval of the Johnson V band) of N_\lambda times S_V(\lambda) equals N_0 photons s^-1 m^-2. The function S_V(\lambda) and the normalisation constant N_0 depend on the value of Planck's constant (parameter :Nature:Planck_Constant), on the definition of the shape of the Johnson V band (parameter :Nature:FilterTransmissionCurve_JohnsonCousinsV_002), on the monochromatic calibration flux f_{0\lambda} (or f_{0\nu}; parameters :Nature:A0VStar_CalibrationFlux_Lambda and :Nature:A0VStar_CalibrationFlux_Nu) at \lambda_0 (parameter :Nature:A0VStar_CalibrationWavelength), and on the spectrum f_{0\nu}(\lambda) of the general unreddened A0V star (parameter :Nature:A0VStar_Spectrum_Nu_002)

            #Source: J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1, Appendix D<br/>
            #Basic : true
            #Scalar: true
            #Unit: [photons s^-1 m^-2]
            """

            return self.__A0VSTAR_CALIBRATIONFUNCTION_NORMALISATION


        __A0VSTAR_CALIBRATIONWAVELENGTH: float  = 550.0 # [nm]

        @property
        def A0VSTAR_CALIBRATIONWAVELENGTH(self):
            r"""
            Reference wavelength at which the flux f_{0\lambda} of an unreddened A0V star is calibrated (see also J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)

            #Source: V. Straizys, 1992, 'Multicolor stellar photometry', Pachart Publ. House (Tucson), Table 21<br/>
            #Basic : true
            #Scalar: true
            #Unit: [nm]
            """

            return self.__A0VSTAR_CALIBRATIONWAVELENGTH


        __A0VSTAR_MEANFLUX_NU: float  = 3.58600e-23 # [W m^-2 Hz^-1]

        @property
        def A0VSTAR_MEANFLUX_NU(self):
            r"""
            Weighted mean flux \langle f_{\nu} \rangle over the Johnson V passband (parameter :Nature:FilterTransmissionCurve_JohnsonCousinsV_002) of an unreddened A0V star (parameter :Nature:A0VStar_Spectrum) with V = :Nature:A0VStar_VMagnitude mag

            #Source: P. Montegriffo, 26 January 2016, 'External calibration for Gaia DR1 integrated photometry', GAIA-C5-TN-OABO-PMN-009, issue 1, revision 0, depends on parameters :Nature:A0VStar_Spectrum and :Nature:FilterTransmissionCurve_JohnsonCousinsV_002<br/>
            #Basic : true
            #Scalar: true
            #Unit: [W m^-2 Hz^-1]
            """

            return self.__A0VSTAR_MEANFLUX_NU


        __A0VSTAR_RMINI: float  = -0.001 # [mag]

        @property
        def A0VSTAR_RMINI(self):
            r"""
            Cousins R minus I (R-I) magnitude of an unreddened A0V star

            #Basic : false
            #Scalar: true
            #Unit: [mag]
            """

            return self.__A0VSTAR_RMINI


        __A0VSTAR_SPECTRUM_NU: str  = "Nature/A0VStar_Spectrum_Nu_002.fits"

        @property
        def A0VSTAR_SPECTRUM_NU(self):
            r"""
            Spectrum f_{0\nu}(\lambda) of an unreddened A0V star: high-fidelity, Kurucz-model Vega spectrum (R = 500) with T_eff = 9400 K, log g = 3.90 dex, and [M/H] = -0.5 dex. The Kurucz model has been scaled to fit STIS data (over the interval 554.5-557.0 nm) by a factor 0.994242. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: Eddington flux (in W m^-2 Hz^-1 steradian^-1). Note that the flux at 115.0 nm was obtained using linear interpolation between the available fluxes at 114.9721 and 115.0873 nm (2.521153898563741E-008 and 2.420265114233843E-008, respectively). Note that the flux at 1062.0 nm was obtained using linear interpolation between the available fluxes at 1061.1654 and 1062.2293 nm (2.508019694385564E-005 and 2.504881789424158E-005, respectively)

            #Source: R.C. Bohlin, 2014, 'Hubble Space Telescope CALSPEC Flux Standards: Sirius (and Vega)', Astronomical Journal, Volume 147, 127; spectrum named alpha_lyr_mod_002.fits contained in CALSPEC Calibration Database (http://www.stsci.edu/hst/observatory/crds/calspec.html, last modified April 2015)<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__A0VSTAR_SPECTRUM_NU


        __A0VSTAR_SPECTRUM_NUMBEROFPHOTONS: str  = "Nature/A0VStar_Spectrum_NumberOfPhotons_001.fits"

        @property
        def A0VSTAR_SPECTRUM_NUMBEROFPHOTONS(self):
            r"""
            Photon flux density N_{\lambda}(\lambda) of an unreddened A0V star (Pickles' star number 009) at V = 0 mag. Note that this unreddened A0V star refers to Pickles' star number 009, and not to the Kurucz Vega spectrum which has been used for flux normalisation and zero-point definition. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)

            #Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__A0VSTAR_SPECTRUM_NUMBEROFPHOTONS


        __A0VSTAR_VMAGNITUDE: float  = 0.023 # [mag]

        @property
        def A0VSTAR_VMAGNITUDE(self):
            r"""
            Johnson V magnitude of Vega

            #Source: R.C. Bohlin, 2007, 'HST Stellar Standards with 1\% Accuracy in Absolute Flux', in 'The Future of Photometric, Spectrophotometric and Polarimetric Standardization', ASP Conference Series, Vol. 364, p.315; see also CALSPEC Calibration Database, http://www.stsci.edu/hst/observatory/crds/calspec.html, last modified April 2015<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__A0VSTAR_VMAGNITUDE


        __A0VSTAR_VMING: float  = 0.000 # [mag]

        @property
        def A0VSTAR_VMING(self):
            r"""
            Johnson V minus Gaia G (V-G) magnitude of an unreddened A0V star, applicable to any photometric band G

            #Source: Definition of Gaia G band(s)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__A0VSTAR_VMING


        __A0VSTAR_VMINI: float  = -0.001 # [mag]

        @property
        def A0VSTAR_VMINI(self):
            r"""
            Johnson V minus Cousins I (V-I) magnitude of an unreddened A0V star

            #Source: M.S. Bessell, S. Murphy, 2012, 'Spectrophotometric Libraries, Revised Photonic Passbands, and Zero Points for UBVRI, Hipparcos, and Tycho Photometry', PASP, 124, 140<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__A0VSTAR_VMINI


        __A0VSTAR_VMINR: float  = 0.000 # [mag]

        @property
        def A0VSTAR_VMINR(self):
            r"""
            Johnson V minus Cousins R (V-R) magnitude of an unreddened A0V star

            #Source: M.S. Bessell, S. Murphy, 2012, 'Spectrophotometric Libraries, Revised Photonic Passbands, and Zero Points for UBVRI, Hipparcos, and Tycho Photometry', PASP, 124, 140<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__A0VSTAR_VMINR


        __A3VSTAR_SPECTRUM_NUMBEROFPHOTONS: str  = "Nature/A3VStar_Spectrum_NumberOfPhotons_001.fits"

        @property
        def A3VSTAR_SPECTRUM_NUMBEROFPHOTONS(self):
            r"""
            Photon flux density N_{\lambda}(\lambda) of an unreddened A3V star (Pickles' star number 011) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)

            #Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__A3VSTAR_SPECTRUM_NUMBEROFPHOTONS


        __A5VSTAR_SPECTRUM_NUMBEROFPHOTONS: str  = "Nature/A5VStar_Spectrum_NumberOfPhotons_001.fits"

        @property
        def A5VSTAR_SPECTRUM_NUMBEROFPHOTONS(self):
            r"""
            Photon flux density N_{\lambda}(\lambda) of an unreddened A5V star (Pickles' star number 012) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)

            #Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__A5VSTAR_SPECTRUM_NUMBEROFPHOTONS


        __ABERRATION_CONSTANT_J2000: float  = 20.49122 # [arcsec]

        @property
        def ABERRATION_CONSTANT_J2000(self):
            r"""
            Constant of aberration, nowadays irrelevant as a fundamental constant, at the standard epoch J2000.0. The IAU (1976) System of Astronomical Constants (e.g., T. Lederle, 1980, MitAG, 48, 59, Table 1) lists 20.49552 arcsec

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Equation 3.253-4, page 131<br/>
            #Basic : false
            #Scalar: true
            #Unit: [arcsec]
            """

            return self.__ABERRATION_CONSTANT_J2000


        __ANGSTROM_NANOMETER: float  = 0.1 # [nm]

        @property
        def ANGSTROM_NANOMETER(self):
            r"""
            One Angstrom expressed in units of nm. Note that 'Angstrom' is a non-SI unit which should not be used

            #Basic : true
            #Scalar: true
            #Unit: [nm]
            """

            return self.__ANGSTROM_NANOMETER


        __ARCSECOND_RADIAN: float  = 4.848136811095360e-06 # [rad]

        @property
        def ARCSECOND_RADIAN(self):
            r"""
            One arcsecond in units of radians

            #Basic : false
            #Scalar: true
            #Unit: [rad]
            """

            return self.__ARCSECOND_RADIAN


        __ASTEROID1CERESMASS_SOLARMASS: float  = 4.720e-10

        @property
        def ASTEROID1CERESMASS_SOLARMASS(self):
            r"""
            Ratio of 1 Ceres to solar mass (IAU 2009 CBE value)

            #Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__ASTEROID1CERESMASS_SOLARMASS


        __ASTEROID1CERES_DIAMETER: float  = 933.0 # [km]

        @property
        def ASTEROID1CERES_DIAMETER(self):
            r"""
            Diameter of 1 Ceres

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [km]
            """

            return self.__ASTEROID1CERES_DIAMETER


        __ASTEROID1CERES_LIGHTDEFLECTION_LIMB: float  = 1.0 # [10^-6 arcsec]

        @property
        def ASTEROID1CERES_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of 1 Ceres

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__ASTEROID1CERES_LIGHTDEFLECTION_LIMB


        __ASTEROID1CERES_ORBITALECCENTRICITY_B1950: float  = 0.0780

        @property
        def ASTEROID1CERES_ORBITALECCENTRICITY_B1950(self):
            r"""
            Orbital eccentricity of 1 Ceres (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__ASTEROID1CERES_ORBITALECCENTRICITY_B1950


        __ASTEROID1CERES_ORBITALPERIOD_B1950: float  = 4.607 # [yr]

        @property
        def ASTEROID1CERES_ORBITALPERIOD_B1950(self):
            r"""
            Orbital period of 1 Ceres (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [yr]
            """

            return self.__ASTEROID1CERES_ORBITALPERIOD_B1950


        __ASTEROID1CERES_ORBITALSEMIMAJORAXIS_B1950: float  = 2.769 # [au]

        @property
        def ASTEROID1CERES_ORBITALSEMIMAJORAXIS_B1950(self):
            r"""
            Orbital semi-major axis of 1 Ceres (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [au]
            """

            return self.__ASTEROID1CERES_ORBITALSEMIMAJORAXIS_B1950


        __ASTEROID1CERES_PERIHELIONDISTANCE_B1950: float  = 2.553 # [au]

        @property
        def ASTEROID1CERES_PERIHELIONDISTANCE_B1950(self):
            r"""
            Perihelion distance of 1 Ceres (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic

            #Basic : false
            #Scalar: true
            #Unit: [au]
            """

            return self.__ASTEROID1CERES_PERIHELIONDISTANCE_B1950


        __ASTEROID2PALLASMASS_SOLARMASS: float  = 1.030e-10

        @property
        def ASTEROID2PALLASMASS_SOLARMASS(self):
            r"""
            Ratio of 2 Pallas to solar mass (IAU 2009 CBE value)

            #Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__ASTEROID2PALLASMASS_SOLARMASS


        __ASTEROID2PALLAS_DIAMETER: float  = 525.0 # [km]

        @property
        def ASTEROID2PALLAS_DIAMETER(self):
            r"""
            Diameter of 2 Pallas

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [km]
            """

            return self.__ASTEROID2PALLAS_DIAMETER


        __ASTEROID2PALLAS_LIGHTDEFLECTION_LIMB: float  = 0.0 # [10^-6 arcsec]

        @property
        def ASTEROID2PALLAS_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of 2 Pallas

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__ASTEROID2PALLAS_LIGHTDEFLECTION_LIMB


        __ASTEROID2PALLAS_ORBITALECCENTRICITY_B1950: float  = 0.2347

        @property
        def ASTEROID2PALLAS_ORBITALECCENTRICITY_B1950(self):
            r"""
            Orbital eccentricity of 2 Pallas (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__ASTEROID2PALLAS_ORBITALECCENTRICITY_B1950


        __ASTEROID2PALLAS_ORBITALPERIOD_B1950: float  = 4.611 # [yr]

        @property
        def ASTEROID2PALLAS_ORBITALPERIOD_B1950(self):
            r"""
            Orbital period of 2 Pallas (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [yr]
            """

            return self.__ASTEROID2PALLAS_ORBITALPERIOD_B1950


        __ASTEROID2PALLAS_ORBITALSEMIMAJORAXIS_B1950: float  = 2.770 # [au]

        @property
        def ASTEROID2PALLAS_ORBITALSEMIMAJORAXIS_B1950(self):
            r"""
            Orbital semi-major axis of 2 Pallas (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [au]
            """

            return self.__ASTEROID2PALLAS_ORBITALSEMIMAJORAXIS_B1950


        __ASTEROID2PALLAS_PERIHELIONDISTANCE_B1950: float  = 2.120 # [au]

        @property
        def ASTEROID2PALLAS_PERIHELIONDISTANCE_B1950(self):
            r"""
            Perihelion distance of 2 Pallas (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic

            #Basic : false
            #Scalar: true
            #Unit: [au]
            """

            return self.__ASTEROID2PALLAS_PERIHELIONDISTANCE_B1950


        __ASTEROID3JUNO_DIAMETER: float  = 267.0 # [km]

        @property
        def ASTEROID3JUNO_DIAMETER(self):
            r"""
            Diameter of 3 Juno

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf). Note that a radius of 120 km is found by E.F. Tedesco, 1989, 'Asteroid magnitudes, UBV colors, and IRAS albedos and diameters', in 'Asteroids II', proceedings of the Conference, Tucson, AZ, 8-11 March 1988, eds R.P. Binzel, T. Gehrels, M.S. Matthews, University of Arizona Press, page 1090 (1989aste.conf.1090T)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [km]
            """

            return self.__ASTEROID3JUNO_DIAMETER


        __ASTEROID3JUNO_LIGHTDEFLECTION_LIMB: float  = 0.0 # [10^-6 arcsec]

        @property
        def ASTEROID3JUNO_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of 3 Juno

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__ASTEROID3JUNO_LIGHTDEFLECTION_LIMB


        __ASTEROID3JUNO_MASS: float  = 2.990e+19 # [kg]

        @property
        def ASTEROID3JUNO_MASS(self):
            r"""
            Mass of 3 Juno

            #Source: Value is calculated following J.L. Hilton, 1999, 'US Naval Observatory Ephemerides of the Largest Asteroids', AJ, 117, 1077, who assumes a mean mass density of 3 g cm^-3. A mass of 2.0E19 kg is found on http://nssdc.gsfc.nasa.gov/planetary/factsheet/asteroidfact.html<br/>
            #Basic : false
            #Scalar: true
            #Unit: [kg]
            """

            return self.__ASTEROID3JUNO_MASS


        __ASTEROID3JUNO_ORBITALECCENTRICITY_B1950: float  = 0.0258

        @property
        def ASTEROID3JUNO_ORBITALECCENTRICITY_B1950(self):
            r"""
            Orbital eccentricity of 3 Juno (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__ASTEROID3JUNO_ORBITALECCENTRICITY_B1950


        __ASTEROID3JUNO_ORBITALPERIOD_B1950: float  = 4.359 # [yr]

        @property
        def ASTEROID3JUNO_ORBITALPERIOD_B1950(self):
            r"""
            Orbital period of 3 Juno (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [yr]
            """

            return self.__ASTEROID3JUNO_ORBITALPERIOD_B1950


        __ASTEROID3JUNO_ORBITALSEMIMAJORAXIS_B1950: float  = 2.668 # [au]

        @property
        def ASTEROID3JUNO_ORBITALSEMIMAJORAXIS_B1950(self):
            r"""
            Orbital semi-major axis of 3 Juno (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [au]
            """

            return self.__ASTEROID3JUNO_ORBITALSEMIMAJORAXIS_B1950


        __ASTEROID3JUNO_PERIHELIONDISTANCE_B1950: float  = 2.599 # [au]

        @property
        def ASTEROID3JUNO_PERIHELIONDISTANCE_B1950(self):
            r"""
            Perihelion distance of 3 Juno (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic

            #Basic : false
            #Scalar: true
            #Unit: [au]
            """

            return self.__ASTEROID3JUNO_PERIHELIONDISTANCE_B1950


        __ASTEROID4VESTAMASS_SOLARMASS: float  = 1.350e-10

        @property
        def ASTEROID4VESTAMASS_SOLARMASS(self):
            r"""
            Ratio of 4 Vesta to solar mass (IAU 2009 CBE value)

            #Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__ASTEROID4VESTAMASS_SOLARMASS


        __ASTEROID4VESTA_DIAMETER: float  = 510.0 # [km]

        @property
        def ASTEROID4VESTA_DIAMETER(self):
            r"""
            Diameter of 4 Vesta

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [km]
            """

            return self.__ASTEROID4VESTA_DIAMETER


        __ASTEROID4VESTA_LIGHTDEFLECTION_LIMB: float  = 1.0 # [10^-6 arcsec]

        @property
        def ASTEROID4VESTA_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of 4 Vesta

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__ASTEROID4VESTA_LIGHTDEFLECTION_LIMB


        __ASTEROID4VESTA_ORBITALECCENTRICITY_B1950: float  = 0.0906

        @property
        def ASTEROID4VESTA_ORBITALECCENTRICITY_B1950(self):
            r"""
            Orbital eccentricity of 4 Vesta (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__ASTEROID4VESTA_ORBITALECCENTRICITY_B1950


        __ASTEROID4VESTA_ORBITALPERIOD_B1950: float  = 3.629 # [yr]

        @property
        def ASTEROID4VESTA_ORBITALPERIOD_B1950(self):
            r"""
            Orbital period of 4 Vesta (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [yr]
            """

            return self.__ASTEROID4VESTA_ORBITALPERIOD_B1950


        __ASTEROID4VESTA_ORBITALSEMIMAJORAXIS_B1950: float  = 2.361 # [au]

        @property
        def ASTEROID4VESTA_ORBITALSEMIMAJORAXIS_B1950(self):
            r"""
            Orbital semi-major axis of 4 Vesta (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [au]
            """

            return self.__ASTEROID4VESTA_ORBITALSEMIMAJORAXIS_B1950


        __ASTEROID4VESTA_PERIHELIONDISTANCE_B1950: float  = 2.147 # [au]

        @property
        def ASTEROID4VESTA_PERIHELIONDISTANCE_B1950(self):
            r"""
            Perihelion distance of 4 Vesta (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic

            #Basic : false
            #Scalar: true
            #Unit: [au]
            """

            return self.__ASTEROID4VESTA_PERIHELIONDISTANCE_B1950


        __ASTEROIDBELTMASS_SOLARMASS: float  = 1.400e-09

        @property
        def ASTEROIDBELTMASS_SOLARMASS(self):
            r"""
            Total mass of the main asteroid belt, in units of solar masses

            #Source: E. Pitjeva, 2003, 'The Dynamic Estimation of the Mass of the Main Asteroid Belt', in 'Physical Properties and Morphology of Small Solar-System Bodies', XXV-th General Assembly of the IAU, Joint Discussion 19, 23 July 2003, Sidney, Australia (2003IAUJD..19E..22P)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__ASTEROIDBELTMASS_SOLARMASS


        __ASTEROID_VELOCITYDISTRIBUTIONSTANDARDDEVIATION_AC: float  = 13.0 # [mas s^-1]

        @property
        def ASTEROID_VELOCITYDISTRIBUTIONSTANDARDDEVIATION_AC(self):
            r"""
            The velocity distribution of main-belt asteroids (MBOs) is very roughly Gaussian with zero mean and a standard deviation of 13.0 mas s^-1 across-scan (for a solar-aspect angle of 45 degrees)

            #Source: F. Mignard, 2002, 'Observations of solar-system objects with Gaia. I. Detection of NEOs', A&A, 393, 727, Section 4.4 (2002A&A...393..727M). See also E. Hoeg, F. Arenou, P. Hjorth, U.G. Joergensen, F. Mignard, S. Wolff, 28 February 2003, 'Faint objects and NEOs with Gaia', GAIA-CUO-118, issue 1, revision 0 and F. Mignard, 22 June 2001, 'Observation of main-belt asteroids with Gaia', GAIA-FM-009, issue 1, revision 0. Current value, for a solar-aspect angle of 45 degrees, from F. Mignard, priv. comm., 10 August 2005<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mas s^-1]
            """

            return self.__ASTEROID_VELOCITYDISTRIBUTIONSTANDARDDEVIATION_AC


        __ASTEROID_VELOCITYDISTRIBUTIONSTANDARDDEVIATION_AL: float  = 7.0 # [mas s^-1]

        @property
        def ASTEROID_VELOCITYDISTRIBUTIONSTANDARDDEVIATION_AL(self):
            r"""
            The velocity distribution of main-belt asteroids (MBOs) is very roughly Gaussian with zero mean and a standard deviation of 7.0 mas s^-1 along-scan (for a solar-aspect angle of 45 degrees)

            #Source: F. Mignard, 2002, 'Observations of solar-system objects with Gaia. I. Detection of NEOs', A&A, 393, 727, Section 4.4 (2002A&A...393..727M). See also E. Hoeg, F. Arenou, P. Hjorth, U.G. Joergensen, F. Mignard, S. Wolff, 28 February 2003, 'Faint objects and NEOs with Gaia', GAIA-CUO-118, issue 1, revision 0 and F. Mignard, 22 June 2001, 'Observation of main-belt asteroids with Gaia', GAIA-FM-009, issue 1, revision 0. Current value, for a solar-aspect angle of 45 degrees, from F. Mignard, priv. comm., 10 August 2005<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mas s^-1]
            """

            return self.__ASTEROID_VELOCITYDISTRIBUTIONSTANDARDDEVIATION_AL


        __ASTRONOMICALUNIT_METER: float  = 149597870700.0 # [m]

        @property
        def ASTRONOMICALUNIT_METER(self):
            r"""
            Astronomical unit (au) length. The au is a conventional unit of length and is a defining constant. The numerical value is in agreement with the value adopted in IAU 2009 Resolution B2. The definition applies to all time scales such as TCB, TDB, TCG, TT, etc.

            #Source: IAU, August 2012, 'Re-definition of the astronomical unit of length', IAU 2012 Resolution B2 adopted at the XXVIII-th General Assembly of the IAU<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__ASTRONOMICALUNIT_METER


        __ATOMICMASS_CONSTANT: float  = 1.6605390404e-27 # [kg]

        @property
        def ATOMICMASS_CONSTANT(self):
            r"""
            Atomic mass constant (also known as atomic mass unit [amu]; 1 amu is defined as 1/12-th of the mass of a 12-C atom). Note: best-measured value equals 1.660539040E-27 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))

            #Basic : false
            #Scalar: true
            #Unit: [kg]
            """

            return self.__ATOMICMASS_CONSTANT


        __AVOGADRO_CONSTANT: float  = 6.0221408570e+23 # [mol^-1]

        @property
        def AVOGADRO_CONSTANT(self):
            r"""
            Avogadro's constant

            #Source: P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mol^-1]
            """

            return self.__AVOGADRO_CONSTANT


        __B0ISTAR_SPECTRUM_NUMBEROFPHOTONS: str  = "Nature/B0IStar_Spectrum_NumberOfPhotons_001.fits"

        @property
        def B0ISTAR_SPECTRUM_NUMBEROFPHOTONS(self):
            r"""
            Photon flux density N_{\lambda}(\lambda) of an unreddened B0I star (Pickles' star number 114) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)

            #Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__B0ISTAR_SPECTRUM_NUMBEROFPHOTONS


        __B1VSTAR_SPECTRUM_NUMBEROFPHOTONS: str  = "Nature/B1VStar_Spectrum_NumberOfPhotons_001.fits"

        @property
        def B1VSTAR_SPECTRUM_NUMBEROFPHOTONS(self):
            r"""
            Photon flux density N_{\lambda}(\lambda) of an unreddened B1V star (Pickles' star number 004) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)

            #Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__B1VSTAR_SPECTRUM_NUMBEROFPHOTONS


        __B1VSTAR_SPECTRUM_NUMBEROFPHOTONSHIGHRESOLUTION: str  = "Nature/B1VStar_Spectrum_NumberOfPhotonsHighResolution_001.fits"

        @property
        def B1VSTAR_SPECTRUM_NUMBEROFPHOTONSHIGHRESOLUTION(self):
            r"""
            High-resolution photon-flux density N_{\lambda}(\lambda) of an unreddened B1V star at V = 15 mag. The data refer to a high-resolution Kurucz-model spectrum with the following properties: effective temperature T_eff = 25500 K, logarithm of surface gravity log g = 4.0, metallicity [Fe/H] = 0.0, alpha-elements [\alpha/Fe] = 0.0, rotational velocity v sini = 50 km s^-1, micro-turbulence velocity = 2.0 km s^-1, length of convective bubble divided by pressure scale height = 0.50, no convective overshooting, macro-turbulence velocity = 2.0 km s^-1, and resolving power R = \lambda / \delta \lambda = 250,000. First column: wavelength \lambda (in nm; from 830.1673264 to 889.8217922). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1). The 34698 lines have an average wavelength step of 0.00172 nm; the spectrum extent is thus 59.7 nm

            #Source: ESA, 20 June 2005, 'Photon-flux distributions for reference stars', GAIA-EST-TN-00539, issue 1, revision 0, based on D. Katz, priv. comm., 11 May 2005<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__B1VSTAR_SPECTRUM_NUMBEROFPHOTONSHIGHRESOLUTION


        __BOHRRADIUS_CONSTANT: float  = 5.29177210564e-11 # [m]

        @property
        def BOHRRADIUS_CONSTANT(self):
            r"""
            Bohr radius. Note: best-measured value equals 0.52917721067E-10 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__BOHRRADIUS_CONSTANT


        __BOLTZMANN_CONSTANT: float  = 1.380648510e-23 # [J K^-1]

        @property
        def BOLTZMANN_CONSTANT(self):
            r"""
            Boltzmann's constant. Note: best-measured value equals 1.38064852E-23 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))

            #Basic : false
            #Scalar: true
            #Unit: [J K^-1]
            """

            return self.__BOLTZMANN_CONSTANT


        __CAIITRIPLET_EQUIVALENTWIDTH_1SUN: float  = 0.146 # [nm]

        @property
        def CAIITRIPLET_EQUIVALENTWIDTH_1SUN(self):
            r"""
            Solar value of the equivalent width of the first line of the CaII-triplet (3p^{6}3d ^{2}D_{3/2} - 3p^{6}4p ^{2}P_{3/2} transition). A relative intensity of 130 is mentioned in Y. Ralchenko, A.E. Kramida, J. Reader, and the NIST ASD Team, 4 October 2010, 'NIST Atomic Spectra Database', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.nist.gov/pml/data/asd.cfm (Web Version 4.0)

            #Source: C.E. Moore, M.G.J. Minnaert, J. Houtgast, 1966, 'The solar spectrum 2935 AA to 8770 AA', US National Bureau of Standards, Monograph 61<br/>
            #Basic : true
            #Scalar: true
            #Unit: [nm]
            """

            return self.__CAIITRIPLET_EQUIVALENTWIDTH_1SUN


        __CAIITRIPLET_EQUIVALENTWIDTH_2SUN: float  = 0.367 # [nm]

        @property
        def CAIITRIPLET_EQUIVALENTWIDTH_2SUN(self):
            r"""
            Solar value of the equivalent width of the second line of the CaII-triplet (3p^{6}3d ^{2}D_{5/2} - 3p^{6}4p ^{2}P_{3/2} transition). A relative intensity of 170 is mentioned in Y. Ralchenko, A.E. Kramida, J. Reader, and the NIST ASD Team, 4 October 2010, 'NIST Atomic Spectra Database', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.nist.gov/pml/data/asd.cfm (Web Version 4.0)

            #Source: C.E. Moore, M.G.J. Minnaert, J. Houtgast, 1966, 'The solar spectrum 2935 AA to 8770 AA', US National Bureau of Standards, Monograph 61<br/>
            #Basic : true
            #Scalar: true
            #Unit: [nm]
            """

            return self.__CAIITRIPLET_EQUIVALENTWIDTH_2SUN


        __CAIITRIPLET_EQUIVALENTWIDTH_3SUN: float  = 0.260 # [nm]

        @property
        def CAIITRIPLET_EQUIVALENTWIDTH_3SUN(self):
            r"""
            Solar value of the equivalent width of the third line of the CaII-triplet (3p^{6}3d ^{2}D_{3/2} - 3p^{6}4p ^{2}P_{1/2} transition). A relative intensity of 160 is mentioned in Y. Ralchenko, A.E. Kramida, J. Reader, and the NIST ASD Team, 4 October 2010, 'NIST Atomic Spectra Database', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.nist.gov/pml/data/asd.cfm (Web Version 4.0)

            #Source: C.E. Moore, M.G.J. Minnaert, J. Houtgast, 1966, 'The solar spectrum 2935 AA to 8770 AA', US National Bureau of Standards, Monograph 61<br/>
            #Basic : true
            #Scalar: true
            #Unit: [nm]
            """

            return self.__CAIITRIPLET_EQUIVALENTWIDTH_3SUN


        __CAIITRIPLET_OSCILLATORSTRENGTH_1: float  = -1.318

        @property
        def CAIITRIPLET_OSCILLATORSTRENGTH_1(self):
            r"""
            Oscillator strength of the first line of the CaII-triplet (3p^{6}3d ^{2}D_{3/2} - 3p^{6}4p ^{2}P_{3/2} transition). The estimated accuracy is better than 25%

            #Source: Y. Ralchenko, A.E. Kramida, J. Reader, and the NIST ASD Team, 4 October 2010, 'NIST Atomic Spectra Database', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.nist.gov/pml/data/asd.cfm (Web Version 4.0)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__CAIITRIPLET_OSCILLATORSTRENGTH_1


        __CAIITRIPLET_OSCILLATORSTRENGTH_2: float  = -0.36

        @property
        def CAIITRIPLET_OSCILLATORSTRENGTH_2(self):
            r"""
            Oscillator strength of the second line of the CaII-triplet (3p^{6}3d ^{2}D_{5/2} - 3p^{6}4p ^{2}P_{3/2} transition). The estimated accuracy is better than 25%

            #Source: Y. Ralchenko, A.E. Kramida, J. Reader, and the NIST ASD Team, 4 October 2010, 'NIST Atomic Spectra Database', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.nist.gov/pml/data/asd.cfm (Web Version 4.0)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__CAIITRIPLET_OSCILLATORSTRENGTH_2


        __CAIITRIPLET_OSCILLATORSTRENGTH_3: float  = -0.622

        @property
        def CAIITRIPLET_OSCILLATORSTRENGTH_3(self):
            r"""
            Oscillator strength of the third line of the CaII-triplet (3p^{6}3d ^{2}D_{3/2} - 3p^{6}4p ^{2}P_{1/2} transition). The estimated accuracy is better than 25%

            #Source: Y. Ralchenko, A.E. Kramida, J. Reader, and the NIST ASD Team, 4 October 2010, 'NIST Atomic Spectra Database', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.nist.gov/pml/data/asd.cfm (Web Version 4.0)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__CAIITRIPLET_OSCILLATORSTRENGTH_3


        __CAIITRIPLET_WAVELENGTH_1: float  = 850.035 # [nm]

        @property
        def CAIITRIPLET_WAVELENGTH_1(self):
            r"""
            Rest-wavelength in vacuum of the first line of the CaII-triplet (3p^{6}3d ^{2}D_{3/2} - 3p^{6}4p ^{2}P_{3/2} transition), as calculated from the difference between the energy of the upper and lower level of the transition

            #Source: Y. Ralchenko, A.E. Kramida, J. Reader, and the NIST ASD Team, 4 October 2010, 'NIST Atomic Spectra Database', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.nist.gov/pml/data/asd.cfm (Web Version 4.0)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [nm]
            """

            return self.__CAIITRIPLET_WAVELENGTH_1


        __CAIITRIPLET_WAVELENGTH_2: float  = 854.444 # [nm]

        @property
        def CAIITRIPLET_WAVELENGTH_2(self):
            r"""
            Rest-wavelength in vacuum of the second line of the CaII-triplet (3p^{6}3d ^{2}D_{5/2} - 3p^{6}4p ^{2}P_{3/2} transition), as calculated from the difference between the energy of the upper and lower level of the transition

            #Source: Y. Ralchenko, A.E. Kramida, J. Reader, and the NIST ASD Team, 4 October 2010, 'NIST Atomic Spectra Database', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.nist.gov/pml/data/asd.cfm (Web Version 4.0)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [nm]
            """

            return self.__CAIITRIPLET_WAVELENGTH_2


        __CAIITRIPLET_WAVELENGTH_3: float  = 866.452 # [nm]

        @property
        def CAIITRIPLET_WAVELENGTH_3(self):
            r"""
            Rest-wavelength in vacuum of the third line of the CaII-triplet (3p^{6}3d ^{2}D_{3/2} - 3p^{6}4p ^{2}P_{1/2} transition), as calculated from the difference between the energy of the upper and lower level of the transition

            #Source: Y. Ralchenko, A.E. Kramida, J. Reader, and the NIST ASD Team, 4 October 2010, 'NIST Atomic Spectra Database', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.nist.gov/pml/data/asd.cfm (Web Version 4.0)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [nm]
            """

            return self.__CAIITRIPLET_WAVELENGTH_3


        __CHARON_ENCOMPASSINGSPHERERADIUS: float  = 6.050e+05 # [m]

        @property
        def CHARON_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around Charon which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__CHARON_ENCOMPASSINGSPHERERADIUS


        __CHARON_GM: float  = 1.0320e+11 # [m^3 s^-2]

        @property
        def CHARON_GM(self):
            r"""
            GM of Charon

            #Source: R.A. Jacobson, 2007, 'Constants used in the PLU017 ephemeris', priv. comm.; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__CHARON_GM


        __CHARON_GEOMETRICALBEDO: float  = 0.372

        @property
        def CHARON_GEOMETRICALBEDO(self):
            r"""
            Geometric albedo of Charon (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)

            #Source: K. Reinsch, V. Burwitz, and M.C. Festou, 1994, 'Albedo Maps of Pluto and Improved Physical Parameters of the Pluto-Charon System', Icarus, 108, 209-218; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__CHARON_GEOMETRICALBEDO


        __CHARON_LIGHTDEFLECTION_LIMB: float  = 1.0 # [10^-6 arcsec]

        @property
        def CHARON_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of Charon

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__CHARON_LIGHTDEFLECTION_LIMB


        __CHARON_MASS: float  = 1.4705e+21 # [kg]

        @property
        def CHARON_MASS(self):
            r"""
            Mass of Charon (do not use for high-precision (orbit) calculations)

            #Basic : false
            #Scalar: true
            #Unit: [kg]
            """

            return self.__CHARON_MASS


        __CHARON_MASSDENSITY_MEAN: float  = 1.585 # [g cm^-3]

        @property
        def CHARON_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of Charon

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__CHARON_MASSDENSITY_MEAN


        __CHARON_RADIUS: float  = 6.050e+05 # [m]

        @property
        def CHARON_RADIUS(self):
            r"""
            Radius of Charon

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__CHARON_RADIUS


        __COSMICRAY_CATALOGUE_AFCCD: str  = "Nature/CosmicRay_Catalogue_AFCCD_001.fits" # [e^-]

        @property
        def COSMICRAY_CATALOGUE_AFCCD(self):
            r"""
            Catalogue containing CCD images, in units of photo-electrons, of typical (galactic) cosmic-ray events for an AF CCD (used in BAM, WFS, SM, and AF; this CCD is also used in BP albeit with a different anti-reflection coating). Cosmic rays will be present as constant background all through the mission. Note that cosmic rays and solar protons are distinct particles, collectively refered to as prompt-particle events (PPEs). The catalogue contains 12389 events. The structure of the FITS file is as follows: the first FITS-file extension contains a list of events containing event number, number of pixels across-scan in the image, and number of pixels along-scan in the image. The following extensions contain the individual images ('pixel matrices'), in units of photo-electron counts, one image per extension

            #Source: A. Short (ESA), priv. comm., 12 May 2006<br/>
            #Basic : true
            #Scalar: false
            #Unit: [e^-]
            """

            return self.__COSMICRAY_CATALOGUE_AFCCD


        __COSMICRAY_CATALOGUE_REDENHANCEDCCD: str  = "Nature/CosmicRay_Catalogue_RedEnhancedCCD_001.fits" # [e^-]

        @property
        def COSMICRAY_CATALOGUE_REDENHANCEDCCD(self):
            r"""
            Catalogue containing CCD images, in units of photo-electrons, of typical (galactic) cosmic-ray events for a red-enhanced CCD (used in RP and RVS). Cosmic rays will be present as constant background all through the mission. Note that cosmic rays and solar protons are distinct particles, collectively refered to as prompt-particle events (PPEs). The catalogue contains 4718 events. The structure of the FITS file is as follows: the first FITS-file extension contains a list of events containing event number, number of pixels across-scan in the image, and number of pixels along-scan in the image. The following extensions contain the individual images ('pixel matrices'), in units of photo-electron counts, one image per extension

            #Source: A. Short (ESA), priv. comm., 5 May 2004<br/>
            #Basic : true
            #Scalar: false
            #Unit: [e^-]
            """

            return self.__COSMICRAY_CATALOGUE_REDENHANCEDCCD


        __COSMICRAY_FLUX_L2: float  = 5.0 # [particles cm^-2 s^-1]

        @property
        def COSMICRAY_FLUX_L2(self):
            r"""
            Typical expected (galactic) cosmic-ray flux at L2, in units of particles cm^-2 s^-1. This flux will be present as constant background all through the mission. Note that cosmic rays and solar protons are distinct particles, collectively refered to as prompt-particle events (PPEs). An isotropic prompt-particle event flux N, in units of events cm^-2 s^-1, generates 2 N A / 4 events s^-1 CCD^-1, where A denotes the active-pixel area of the CCD in units of cm^2 (including any reduction as a result of TDI-gating), the factor 2 results from considering 'inflow' through both the illuminated and the non-illuminated faces of the CCD, and the factor 4 results from the 'flat geometry' of the CCD (see J.H.J. de Bruijne, A. Short, 7 September 2005, 'prompt-particle events: from fluxes to count rates', GAIA-JdB-026, issue 1, revision 0)

            #Source: A. Short (ESA), priv. comm., 20 December 2005<br/>
            #Basic : true
            #Scalar: true
            #Unit: [particles cm^-2 s^-1]
            """

            return self.__COSMICRAY_FLUX_L2


        class DE405:

            __ASTEROID1CERESMASS_SOLARMASS: float  = 4.70e-10

            @property
            def ASTEROID1CERESMASS_SOLARMASS(self):
                r"""
                Ratio of 1 Ceres to solar mass (DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__ASTEROID1CERESMASS_SOLARMASS


            __ASTEROID2PALLASMASS_SOLARMASS: float  = 1.00e-10

            @property
            def ASTEROID2PALLASMASS_SOLARMASS(self):
                r"""
                Ratio of 2 Pallas to solar mass (DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__ASTEROID2PALLASMASS_SOLARMASS


            __ASTEROID4VESTAMASS_SOLARMASS: float  = 1.30e-10

            @property
            def ASTEROID4VESTAMASS_SOLARMASS(self):
                r"""
                Ratio of 4 Vesta to solar mass (DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__ASTEROID4VESTAMASS_SOLARMASS


            __ASTEROID_MASSDENSITY_MEANCLASSC: float  = 1.8 # [g cm^-3]

            @property
            def ASTEROID_MASSDENSITY_MEANCLASSC(self):
                r"""
                Mean mass density of C-class asteroids (DE405 value). In JPL's DE ephemerides, masses for the 300 most massive asteroids (except Ceres, Pallas, and Vesta) are derived using the relation GM = 6.27E-22 Radius^3 \rho, where Radius is the known asteroid radius in km and \rho is the mean mass density in g cm^-3 (GM in this relation is in 'solar-system units', i.e., Sun_GM = Gauss_Constant * Gauss_Constant [au^3/2 day^-1 M_Sun^-1/2]; note that the reference erroneously lists a prefactor 6.27E20; typo confirmed by E.M. Standish, priv. comm., 18 September 2003). The mean density \rho is assumed to be constant within each of the three taxonomic asteroid classes C, S, and M

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048, Section IV<br/>
                #Basic : true
                #Scalar: true
                #Unit: [g cm^-3]
                """

                return self.__ASTEROID_MASSDENSITY_MEANCLASSC


            __ASTEROID_MASSDENSITY_MEANCLASSM: float  = 5.0 # [g cm^-3]

            @property
            def ASTEROID_MASSDENSITY_MEANCLASSM(self):
                r"""
                Mean mass density of M-class asteroids (DE405 value). In JPL's DE ephemerides, masses for the 300 most massive asteroids (except Ceres, Pallas, and Vesta) are derived using the relation GM = 6.27E-22 Radius^3 \rho, where Radius is the known asteroid radius in km and \rho is the mean mass density in g cm^-3 (GM in this relation is in 'solar-system units', i.e., Sun_GM = Gauss_Constant * Gauss_Constant [au^3/2 day^-1 M_Sun^-1/2]; note that the reference erroneously lists a prefactor 6.27E20; typo confirmed by E.M. Standish, priv. comm., 18 September 2003). The mean density \rho is assumed to be constant within each of the three taxonomic asteroid classes C, S, and M

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048, Section IV<br/>
                #Basic : true
                #Scalar: true
                #Unit: [g cm^-3]
                """

                return self.__ASTEROID_MASSDENSITY_MEANCLASSM


            __ASTEROID_MASSDENSITY_MEANCLASSS: float  = 2.4 # [g cm^-3]

            @property
            def ASTEROID_MASSDENSITY_MEANCLASSS(self):
                r"""
                Mean mass density of S-class asteroids (DE405 value). In JPL's DE ephemerides, masses for the 300 most massive asteroids (except Ceres, Pallas, and Vesta) are derived using the relation GM = 6.27E-22 Radius^3 \rho, where Radius is the known asteroid radius in km and \rho is the mean mass density in g cm^-3 (GM in this relation is in 'solar-system units', i.e., Sun_GM = Gauss_Constant * Gauss_Constant [au^3/2 day^-1 M_Sun^-1/2]; note that the reference erroneously lists a prefactor 6.27E20; typo confirmed by E.M. Standish, priv. comm., 18 September 2003). The mean density \rho is assumed to be constant within each of the three taxonomic asteroid classes C, S, and M

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048, Section IV<br/>
                #Basic : true
                #Scalar: true
                #Unit: [g cm^-3]
                """

                return self.__ASTEROID_MASSDENSITY_MEANCLASSS


            __ASTRONOMICALUNIT_METER: float  = 1.4959787301053391e+11 # [m]

            @property
            def ASTRONOMICALUNIT_METER(self):
                r"""
                Astronomical unit (au) length (TCB-compatible value; DE405 value; see S.A. Klioner, 2008, 'Relativistic scaling of astronomical quantities and the system of astronomical units', A&A, 478, 951-958)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m]
                """

                return self.__ASTRONOMICALUNIT_METER


            __ASTRONOMICALUNIT_SECOND: float  = 4.9900479154326796e+02 # [s]

            @property
            def ASTRONOMICALUNIT_SECOND(self):
                r"""
                Astronomical unit (au) light time (TCB-compatible value; DE405 value; see S.A. Klioner, 2008, 'Relativistic scaling of astronomical quantities and the system of astronomical units', A&A, 478, 951-958)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : false
                #Scalar: true
                #Unit: [s]
                """

                return self.__ASTRONOMICALUNIT_SECOND


            __ASTRONOMICALUNIT_TDBMETER: float  = 1.4959787301053392e+11 # [m (TDB)]

            @property
            def ASTRONOMICALUNIT_TDBMETER(self):
                r"""
                Astronomical unit (au) length (TDB-compatible value; DE405 value). Do not use this parameter but use the TCB-compatible value from parameter :Nature:DE405:AstronomicalUnit_Meter instead

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m (TDB)]
                """

                return self.__ASTRONOMICALUNIT_TDBMETER


            __ASTRONOMICALUNIT_TDBSECOND: float  = 4.9900478380610000e+02 # [s (TDB)]

            @property
            def ASTRONOMICALUNIT_TDBSECOND(self):
                r"""
                Astronomical unit (au) light time (TDB-compatible value; DE405 value). Do not use this parameter but use the TCB-compatible value from parameter :Nature:DE405:AstronomicalUnit_Second instead

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: [s (TDB)]
                """

                return self.__ASTRONOMICALUNIT_TDBSECOND


            __EARTHTOMOON_MASSRATIO: float  = 81.30056

            @property
            def EARTHTOMOON_MASSRATIO(self):
                r"""
                Ratio of Earth to Moon mass (DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__EARTHTOMOON_MASSRATIO


            __EARTH_EQUATORIALRADIUS: float  = 6378137.0 # [m]

            @property
            def EARTH_EQUATORIALRADIUS(self):
                r"""
                Equatorial radius of the Earth (DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: [m]
                """

                return self.__EARTH_EQUATORIALRADIUS


            __EARTH_GM: float  = 3.986004576184e+14 # [m^3 s^-2]

            @property
            def EARTH_GM(self):
                r"""
                Geocentric gravitational constant (TCB-compatible value; DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__EARTH_GM


            __EARTH_GM_TDB: float  = 3.986004514380e+14 # [m^3 s^-2 (TDB)]

            @property
            def EARTH_GM_TDB(self):
                r"""
                Geocentric gravitational constant (TDB-compatible value; DE405 value). Do not use this parameter but use the TCB-compatible value from parameter :Nature:DE405:Earth_GM instead

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2 (TDB)]
                """

                return self.__EARTH_GM_TDB


            __EARTH_JSUB2DOT: float  = 0.00e-01 # [cy^-1]

            @property
            def EARTH_JSUB2DOT(self):
                r"""
                Secular (long-term) variation of the dynamical form-factor J_2 of the Earth (also known as oblateness and as Stokes' second-degree zonal harmonic of the geopotential) due to the post-glacial rebound of the mantle (DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: [cy^-1]
                """

                return self.__EARTH_JSUB2DOT


            __EARTH_LOVENUMBER_20: float  = 0.34

            @property
            def EARTH_LOVENUMBER_20(self):
                r"""
                Love number k_20 of harmonic (2,0) of the Earth's harmonic potential expansion (rigid-Earth tide / slow zonal tides; DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__EARTH_LOVENUMBER_20


            __EARTH_LOVENUMBER_21: float  = 0.30

            @property
            def EARTH_LOVENUMBER_21(self):
                r"""
                Love number k_21 of harmonic (2,1) of the Earth's harmonic potential expansion (tidal deformation / diurnal tides; DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__EARTH_LOVENUMBER_21


            __EARTH_LOVENUMBER_22: float  = 0.30

            @property
            def EARTH_LOVENUMBER_22(self):
                r"""
                Love number k_22 of harmonic (2,1) of the Earth's harmonic potential expansion (rotational deformation / semi-diurnal tides; DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__EARTH_LOVENUMBER_22


            __EARTH_POTENTIALEXPANSION_C = [ -0.001082626,  0.0,  0.0,  0.000002533,  0.0,  0.0,  0.0,  0.000001616,  0.0,  0.0,  0.0,  0.0 ]

            @property
            def EARTH_POTENTIALEXPANSION_C(self):
                r"""
                Harmonic potential coefficients of the Earth (DE405 values). The vector elements denote C_nm, with (n,m) = (2,0), (2,1), (2,2), (3,0), (3,1), (3,2), (3,3), (4,0), (4,1), (4,2), (4,3), and (4,4). A zonal harmonic J_n is a spherical harmonic of the form P_n(cos\theta), i.e., one which reduces to a Legendre polynomial of degree n. A tesseral harmonic C_nm/S_nm is a spherical harmonic of the form cos/sin(m\phi) P_n^m(cos\theta), where P_n^m is a Legendre function of degree n and order m. Special notations include -C_20 = J_2 = Stokes' second degree zonal harmonic (oblateness), -C_30 = J_3 = Stokes' third degree zonal harmonic, and -C_40 = J_4 = Stokes' fourth degree zonal harmonic

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__EARTH_POTENTIALEXPANSION_C


            __EARTH_POTENTIALEXPANSION_DEGREE: int  = 4

            @property
            def EARTH_POTENTIALEXPANSION_DEGREE(self):
                r"""
                Degree of harmonic expansion of the Earth's potential (DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__EARTH_POTENTIALEXPANSION_DEGREE


            __EARTH_POTENTIALEXPANSION_S = [ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 ]

            @property
            def EARTH_POTENTIALEXPANSION_S(self):
                r"""
                Harmonic potential coefficients of the Earth (DE405 values). The vector elements denote S_nm, with (n,m) = (2,0), (2,1), (2,2), (3,0), (3,1), (3,2), (3,3), (4,0), (4,1), (4,2), (4,3), and (4,4). A zonal harmonic J_n is a spherical harmonic of the form P_n(cos\theta), i.e., one which reduces to a Legendre polynomial of degree n. A tesseral harmonic C_nm/S_nm is a spherical harmonic of the form cos/sin(m\phi) P_n^m(cos\theta), where P_n^m is a Legendre function of degree n and order m

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__EARTH_POTENTIALEXPANSION_S


            __EARTH_TIMEDELAY_20: float  = 0.0 # [day]

            @property
            def EARTH_TIMEDELAY_20(self):
                r"""
                Time delay \tau_20 used to compute tidal effects for harmonic (2,0) of the Earth's harmonic potential expansion (rigid-Earth tide / slow zonal tides; DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: [day]
                """

                return self.__EARTH_TIMEDELAY_20


            __EARTH_TIMEDELAY_21: float  = 0.01290895939 # [day]

            @property
            def EARTH_TIMEDELAY_21(self):
                r"""
                Time delay \tau_21 used to compute tidal effects for harmonic (2,1) of the Earth's harmonic potential expansion (tidal deformation / diurnal tides; DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: [day]
                """

                return self.__EARTH_TIMEDELAY_21


            __EARTH_TIMEDELAY_22: float  = 0.00694178558 # [day]

            @property
            def EARTH_TIMEDELAY_22(self):
                r"""
                Time delay \tau_22 used to compute tidal effects for harmonic (2,1) of the Earth's harmonic potential expansion (rotational deformation / semi-diurnal tides; DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: [day]
                """

                return self.__EARTH_TIMEDELAY_22


            __MOONTOEARTH_MASSRATIO: float  = 0.0123000383

            @property
            def MOONTOEARTH_MASSRATIO(self):
                r"""
                Ratio of Moon to Earth mass (DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : false
                #Scalar: true
                #Unit: []
                """

                return self.__MOONTOEARTH_MASSRATIO


            __MOON_EQUATORIALRADIUS: float  = 1738000.0 # [m]

            @property
            def MOON_EQUATORIALRADIUS(self):
                r"""
                Equatorial radius of the Moon (DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: [m]
                """

                return self.__MOON_EQUATORIALRADIUS


            __MOON_LOVENUMBER_2: float  = 0.0299221167

            @property
            def MOON_LOVENUMBER_2(self):
                r"""
                Love number k_2 of the Moon's harmonic potential expansion, assumed to be the same for all harmonic coefficients of order 2 (DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__MOON_LOVENUMBER_2


            __MOON_POTENTIALEXPANSION_C = [ -99.99,  -99.99,  -99.99,  -0.000008785470,  0.000030803810,  0.000004879807,  0.000001770176,  1.45383E-7,  -0.000007177801,  -0.000001439518,  -8.5479E-8,  -1.54904E-7 ]

            @property
            def MOON_POTENTIALEXPANSION_C(self):
                r"""
                Harmonic potential coefficients of the Moon (DE405 values). The vector elements denote C_nm, with (n,m) = (2,0), (2,1), (2,2), (3,0), (3,1), (3,2), (3,3), (4,0), (4,1), (4,2), (4,3), and (4,4). A zonal harmonic J_n is a spherical harmonic of the form P_n(cos\theta), i.e., one which reduces to a Legendre polynomial of degree n. A tesseral harmonic C_nm/S_nm is a spherical harmonic of the form cos/sin(m\phi) P_n^m(cos\theta), where P_n^m is a Legendre function of degree n and order m. Special notations include -C_20 = J_2 = Stokes' second degree zonal harmonic (oblateness), -C_30 = J_3 = Stokes' third degree zonal harmonic, and -C_40 = J_4 = Stokes' fourth degree zonal harmonic. The second-degree lunar gravity field is time variable and the time-variable harmonic coefficients are computed in the DE405 ephemeris from the time-variable moment-of-inertia tensor. The numerical values of the C_20, C_21, and C_22 coefficients reported here (-99.99) are hence spurious

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__MOON_POTENTIALEXPANSION_C


            __MOON_POTENTIALEXPANSION_DEGREE: int  = 4

            @property
            def MOON_POTENTIALEXPANSION_DEGREE(self):
                r"""
                Degree of harmonic expansion of the Moon's potential (DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__MOON_POTENTIALEXPANSION_DEGREE


            __MOON_POTENTIALEXPANSION_S = [ 0.0,  -99.99,  -99.99,  0.0,  0.000004259329,  0.000001695516,  -2.70970E-7,  0.0,  0.000002947434,  -0.000002884372,  -7.88967E-7,  5.6404E-8 ]

            @property
            def MOON_POTENTIALEXPANSION_S(self):
                r"""
                Harmonic potential coefficients of the Moon (DE405 values). The vector elements denote S_nm, with (n,m) = (2,0), (2,1), (2,2), (3,0), (3,1), (3,2), (3,3), (4,0), (4,1), (4,2), (4,3), and (4,4). A zonal harmonic J_n is a spherical harmonic of the form P_n(cos\theta), i.e., one which reduces to a Legendre polynomial of degree n. A tesseral harmonic C_nm/S_nm is a spherical harmonic of the form cos/sin(m\phi) P_n^m(cos\theta), where P_n^m is a Legendre function of degree n and order m. The second-degree lunar gravity field is time variable and the time-variable harmonic coefficients are computed in the DE405 ephemeris from the time-variable moment-of-inertia tensor. The numerical values of the S_21 and S_22 coefficients reported here (-99.99) are hence spurious

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__MOON_POTENTIALEXPANSION_S


            __MOON_TIMEDELAY_2: float  = 0.1667165558 # [day]

            @property
            def MOON_TIMEDELAY_2(self):
                r"""
                Time delay \tau_2 used to compute tidal effects for the Moon's solid-body tide, assumed to be the same for all harmonic coefficients of order 2 (DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: [day]
                """

                return self.__MOON_TIMEDELAY_2


            __SUNTOEARTHSYSTEM_MASSRATIO: float  = 328900.561400

            @property
            def SUNTOEARTHSYSTEM_MASSRATIO(self):
                r"""
                Ratio of Sun to Earth-system mass (DE405 value). The planetary mass includes the contribution of its satellite, the Moon

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__SUNTOEARTHSYSTEM_MASSRATIO


            __SUNTOEARTH_MASSRATIO: float  = 332946.050895

            @property
            def SUNTOEARTH_MASSRATIO(self):
                r"""
                Ratio of Sun to Earth mass (DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : false
                #Scalar: true
                #Unit: []
                """

                return self.__SUNTOEARTH_MASSRATIO


            __SUNTOJUPITERSYSTEM_MASSRATIO: float  = 1047.3486

            @property
            def SUNTOJUPITERSYSTEM_MASSRATIO(self):
                r"""
                Ratio of Sun to Jupiter-system mass (DE405 value). The planetary mass includes the contribution of its satellites

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048. Numerical value: J.K. Campbell, S.P. Synnott, 1985, 'Gravity field of the Jovian system from Pioneer and Voyager tracking data', AJ, 90, 364; this reference lists GM = 126712767 km^3 s^-2<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__SUNTOJUPITERSYSTEM_MASSRATIO


            __SUNTOMARSSYSTEM_MASSRATIO: float  = 3098708.0

            @property
            def SUNTOMARSSYSTEM_MASSRATIO(self):
                r"""
                Ratio of Sun to Mars-system mass (DE405 value). The planetary mass includes the contribution of its satellites

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048. Numerical value: G.W. Null, 1969, 'A solution for the mass and dynamical oblateness of Mars using Mariner-IV Doppler data', Bull. Am. Astron. Soc., 1, 356<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__SUNTOMARSSYSTEM_MASSRATIO


            __SUNTOMERCURYSYSTEM_MASSRATIO: float  = 6023600.0

            @property
            def SUNTOMERCURYSYSTEM_MASSRATIO(self):
                r"""
                Ratio of Sun to Mercury(-system) mass (DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048. Numerical value: J.D. Anderson, et al., 1987, 'The mass, gravity field, and ephemeris of Mercury', Icarus, 71, 337; this reference lists GM = 22032.09 km^3 s^-2<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__SUNTOMERCURYSYSTEM_MASSRATIO


            __SUNTONEPTUNESYSTEM_MASSRATIO: float  = 19412.24

            @property
            def SUNTONEPTUNESYSTEM_MASSRATIO(self):
                r"""
                Ratio of Sun to Neptune-system mass (DE405 value). The planetary mass includes the contribution of its satellites

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048. Numerical value: R.A. Jacobson, et al., 1991, 'The orbits of Triton and Nereid from spacecraft and Earth-based observations', A&A, 247, 565; this reference lists GM = 6836535 km^3 s^-2<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__SUNTONEPTUNESYSTEM_MASSRATIO


            __SUNTOPLUTOSYSTEM_MASSRATIO: float  = 135200000.0

            @property
            def SUNTOPLUTOSYSTEM_MASSRATIO(self):
                r"""
                Ratio of Sun to Pluto-system mass (DE405 value). The 'planetary' mass includes the contribution of its satellite, Charon

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048. Numerical value: D.J. Tholen, M.W. Buie, 1997, 'The Orbit of Charon', Icarus, 125, 245, although these authors list 1.3522E8 rather than 1.3521E8; this reference lists GM = 981.5 km^3 s^-2<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__SUNTOPLUTOSYSTEM_MASSRATIO


            __SUNTOSATURNSYSTEM_MASSRATIO: float  = 3497.898

            @property
            def SUNTOSATURNSYSTEM_MASSRATIO(self):
                r"""
                Ratio of Sun to Saturn-system mass (DE405 value). The planetary mass includes the contribution of its satellites

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048. Numerical value: J.K. Campbell, J.D. Anderson, 1989, 'Gravity field of the Saturnian system from Pioneer and Voyager tracking data', AJ, 97, 1485; this reference lists GM = 37940630 km^3 s^-2<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__SUNTOSATURNSYSTEM_MASSRATIO


            __SUNTOURANUSSYSTEM_MASSRATIO: float  = 22902.98

            @property
            def SUNTOURANUSSYSTEM_MASSRATIO(self):
                r"""
                Ratio of Sun to Uranus-system mass (DE405 value). The planetary mass includes the contribution of its satellites

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048. Numerical value: R.A. Jacobson, et al., 1992, AJ, 'The masses of Uranus and its major satellites from Voyager tracking data and Earth-based Uranian satellite data', 103, 2068; this reference lists GM = 5794548.6 km^3 s^-2<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__SUNTOURANUSSYSTEM_MASSRATIO


            __SUNTOVENUSSYSTEM_MASSRATIO: float  = 408523.71

            @property
            def SUNTOVENUSSYSTEM_MASSRATIO(self):
                r"""
                Ratio of Sun to Venus(-system) mass (DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048. Numerical value: W.L. Sjogren, et al., 1990, 'Venus - A total mass estimate', Geophysical Research Letters, 17, 1485; this reference lists GM = 324858.60 km^3 s^-2<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__SUNTOVENUSSYSTEM_MASSRATIO


            __SUN_GM: float  = 1.327124482489e+20 # [m^3 s^-2]

            @property
            def SUN_GM(self):
                r"""
                Heliocentric gravitational constant (TCB-compatible value; DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__SUN_GM


            __SUN_GM_TDB: float  = 1.327124461912e+20 # [m^3 s^-2 (TDB)]

            @property
            def SUN_GM_TDB(self):
                r"""
                Heliocentric gravitational constant (TDB-compatible value; DE405 value). Do not use this parameter but use the TCB-compatible value from parameter :Nature:DE405:Sun_GM instead

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2 (TDB)]
                """

                return self.__SUN_GM_TDB


            __SUN_JSUB2: float  = 2.0e-07

            @property
            def SUN_JSUB2(self):
                r"""
                Dynamical form-factor of the Sun (Stokes' second-degree zonal harmonic of the solar potential; DE405 value)

                #Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__SUN_JSUB2


        class DE410:

            __ASTEROID1CERESMASS_SOLARMASS: float  = 4.690e-10

            @property
            def ASTEROID1CERESMASS_SOLARMASS(self):
                r"""
                Ratio of 1 Ceres to solar mass (DE410 value)

                #Source: E.M. Standish, 24 April 2003, 'JPL Planetary Ephemeris DE410', JPL IOM 312.N-03-009, Table III<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__ASTEROID1CERESMASS_SOLARMASS


            __ASTEROID2PALLASMASS_SOLARMASS: float  = 1.050e-10

            @property
            def ASTEROID2PALLASMASS_SOLARMASS(self):
                r"""
                Ratio of 2 Pallas to solar mass (DE410 value)

                #Source: E.M. Standish, 24 April 2003, 'JPL Planetary Ephemeris DE410', JPL IOM 312.N-03-009, Table III<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__ASTEROID2PALLASMASS_SOLARMASS


            __ASTEROID4VESTAMASS_SOLARMASS: float  = 1.360e-10

            @property
            def ASTEROID4VESTAMASS_SOLARMASS(self):
                r"""
                Ratio of 4 Vesta to solar mass (DE410 value)

                #Source: E.M. Standish, 24 April 2003, 'JPL Planetary Ephemeris DE410', JPL IOM 312.N-03-009, Table III<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__ASTEROID4VESTAMASS_SOLARMASS


            __ASTEROIDRINGMASS_SOLARMASS: float  = 1.032e-10

            @property
            def ASTEROIDRINGMASS_SOLARMASS(self):
                r"""
                Ratio of the Krasinsky asteroid ring to solar mass (originally expressed in terms of M_Ceres; DE410 value). Following G.A. Krasinsky, E.V. Pitjeva, M.V. Vasilyev, E.I. Yagudina, 1 February 2002, 'Hidden Mass in the Asteroid Belt', Icarus, 158, 98-105, the gravitational effect of all but the 300 heaviest asteroids can be modeled as an acceleration caused by a solid ring of this mass in the ecliptic plane (see also E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048, Section IV)

                #Source: E.M. Standish, 24 April 2003, 'JPL Planetary Ephemeris DE410', JPL IOM 312.N-03-009, Table III<br/>
                #Basic : false
                #Scalar: true
                #Unit: []
                """

                return self.__ASTEROIDRINGMASS_SOLARMASS


            __ASTEROIDRING_ORBITALSEMIMAJORAXIS: float  = 2.8 # [au]

            @property
            def ASTEROIDRING_ORBITALSEMIMAJORAXIS(self):
                r"""
                Barycentric distance (orbital semi-major axis) of the Krasinsky asteroid ring (DE410 value). Following G.A. Krasinsky, E.V. Pitjeva, M.V. Vasilyev, E.I. Yagudina, 1 February 2002, 'Hidden Mass in the Asteroid Belt', Icarus, 158, 98-105, the gravitational effect of all but the 300 heaviest asteroids can be modeled as an acceleration caused by a solid ring with this barycentric distance/radius in the ecliptic plane (see also E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048, Section IV)

                #Source: E.M. Standish, 24 April 2003, 'JPL Planetary Ephemeris DE410', JPL IOM 312.N-03-009, Table III<br/>
                #Basic : true
                #Scalar: true
                #Unit: [au]
                """

                return self.__ASTEROIDRING_ORBITALSEMIMAJORAXIS


            __ASTEROID_MASSDENSITY_MEANCLASSC: float  = 1.55 # [g cm^-3]

            @property
            def ASTEROID_MASSDENSITY_MEANCLASSC(self):
                r"""
                Mean mass density of C-class asteroids (DE410 value). In JPL's DE ephemerides, masses for the 300 most massive asteroids (except Ceres, Pallas, and Vesta) are derived using the relation GM = 6.27E-22 Radius^3 \rho, where Radius is the known asteroid radius in km and \rho is the mean mass density in g cm^-3 (GM in this relation is in 'solar-system units', i.e., Sun_GM = Gauss_Constant * Gauss_Constant [au^3/2 day^-1 M_Sun^-1/2]; note that the reference erroneously lists a prefactor 6.27E20; typo confirmed by E.M. Standish, priv. comm., 18 September 2003). The mean density \rho is assumed to be constant within each of the three taxonomic asteroid classes C, S, and M (see E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048, Section IV)

                #Source: E.M. Standish, 24 April 2003, 'JPL Planetary Ephemeris DE410', JPL IOM 312.N-03-009, Table III<br/>
                #Basic : true
                #Scalar: true
                #Unit: [g cm^-3]
                """

                return self.__ASTEROID_MASSDENSITY_MEANCLASSC


            __ASTEROID_MASSDENSITY_MEANCLASSM: float  = 4.5 # [g cm^-3]

            @property
            def ASTEROID_MASSDENSITY_MEANCLASSM(self):
                r"""
                Mean mass density of M-class asteroids (DE410 value). In JPL's DE ephemerides, masses for the 300 most massive asteroids (except Ceres, Pallas, and Vesta) are derived using the relation GM = 6.27E-22 Radius^3 \rho, where Radius is the known asteroid radius in km and \rho is the mean mass density in g cm^-3 (GM in this relation is in 'solar-system units', i.e., Sun_GM = Gauss_Constant * Gauss_Constant [au^3/2 day^-1 M_Sun^-1/2]; note that the reference erroneously lists a prefactor 6.27E20; typo confirmed by E.M. Standish, priv. comm., 18 September 2003). The mean density \rho is assumed to be constant within each of the three taxonomic asteroid classes C, S, and M (see E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048, Section IV)

                #Source: E.M. Standish, 24 April 2003, 'JPL Planetary Ephemeris DE410', JPL IOM 312.N-03-009, Table III<br/>
                #Basic : true
                #Scalar: true
                #Unit: [g cm^-3]
                """

                return self.__ASTEROID_MASSDENSITY_MEANCLASSM


            __ASTEROID_MASSDENSITY_MEANCLASSS: float  = 2.13 # [g cm^-3]

            @property
            def ASTEROID_MASSDENSITY_MEANCLASSS(self):
                r"""
                Mean mass density of S-class asteroids (DE410 value). In JPL's DE ephemerides, masses for the 300 most massive asteroids (except Ceres, Pallas, and Vesta) are derived using the relation GM = 6.27E-22 Radius^3 \rho, where Radius is the known asteroid radius in km and \rho is the mean mass density in g cm^-3 (GM in this relation is in 'solar-system units', i.e., Sun_GM = Gauss_Constant * Gauss_Constant [au^3/2 day^-1 M_Sun^-1/2]; note that the reference erroneously lists a prefactor 6.27E20; typo confirmed by E.M. Standish, priv. comm., 18 September 2003). The mean density \rho is assumed to be constant within each of the three taxonomic asteroid classes C, S, and M (see E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048, Section IV)

                #Source: E.M. Standish, 24 April 2003, 'JPL Planetary Ephemeris DE410', JPL IOM 312.N-03-009, Table III<br/>
                #Basic : true
                #Scalar: true
                #Unit: [g cm^-3]
                """

                return self.__ASTEROID_MASSDENSITY_MEANCLASSS


            __SUN_JSUB2: float  = 2.90e-07

            @property
            def SUN_JSUB2(self):
                r"""
                Dynamical form-factor of the Sun (Stokes' second-degree zonal harmonic of the solar potential; DE410 value)

                #Source: E.M. Standish, 24 April 2003, 'JPL Planetary Ephemeris DE410', JPL IOM 312.N-03-009, Table III<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__SUN_JSUB2


        __DAY_SECOND: float  = 86400.0 # [s]

        @property
        def DAY_SECOND(self):
            r"""
            Number of seconds per day

            #Source: IAU definition<br/>
            #Basic : false
            #Scalar: true
            #Unit: [s]
            """

            return self.__DAY_SECOND


        __DEGREE_RADIAN: float  = 1.745329251994330e-02 # [rad]

        @property
        def DEGREE_RADIAN(self):
            r"""
            One degree in units of radians

            #Basic : false
            #Scalar: true
            #Unit: [rad]
            """

            return self.__DEGREE_RADIAN


        __EMBC_ORBITALECCENTRICITY_J2000: float  = 0.01671123

        @property
        def EMBC_ORBITALECCENTRICITY_J2000(self):
            r"""
            Mean orbital eccentricity of the Earth-Moon barycentre (EMBC) orbit, at the standard epoch J2000.0. The mean orbital eccentricity is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. In this fit, each orbital element is allowed to vary linearly with time (the resulting evolution of the orbital eccentricity is -0.00004392 radians per century). The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0. This solution fits the DE405 orbit of the Earth-Moon barycentre to about 22 arcsec. DE405 is based upon the International Celestial Reference Frame (ICRF)

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__EMBC_ORBITALECCENTRICITY_J2000


        __EMBC_ORBITALINCLINATION_J2000: float  = -0.00001531 # [deg]

        @property
        def EMBC_ORBITALINCLINATION_J2000(self):
            r"""
            Mean orbital inclination of the Earth-Moon barycentre (EMBC) orbit, at the standard epoch J2000.0. The mean orbital inclination is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. In this fit, each orbital element is allowed to vary linearly with time (the resulting evolution of the orbital inclination is -0.01294668 degrees per century). The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0. This solution fits the DE405 orbit of the Earth-Moon barycentre to about 22 arcsec. DE405 is based upon the International Celestial Reference Frame (ICRF)

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__EMBC_ORBITALINCLINATION_J2000


        __EMBC_ORBITALSEMIMAJORAXIS_J2000: float  = 1.00000261 # [au]

        @property
        def EMBC_ORBITALSEMIMAJORAXIS_J2000(self):
            r"""
            Mean orbital semi-major axis of the Earth-Moon barycentre (EMBC) orbit, at the standard epoch J2000.0. The mean orbital semi-major axis is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. In this fit, each orbital element is allowed to vary linearly with time (the resulting evolution of the orbital semi-major axis is 0.00000562 au per century). The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0. This solution fits the DE405 orbit of the Earth-Moon barycentre to about 22 arcsec. DE405 is based upon the International Celestial Reference Frame (ICRF)

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: [au]
            """

            return self.__EMBC_ORBITALSEMIMAJORAXIS_J2000


        __EARTHELLIPSOID_GM: float  = 3.9860044180e+14 # [m^3 s^-2]

        @property
        def EARTHELLIPSOID_GM(self):
            r"""
            Geocentric gravitational constant (TCB-compatible value), including the Earth's atmosphere but excluding the mass of the Moon

            #Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__EARTHELLIPSOID_GM


        __EARTHELLIPSOID_INVERSEFLATTENING_ZEROFREQUENCYTIDE: float  = 298.25642

        @property
        def EARTHELLIPSOID_INVERSEFLATTENING_ZEROFREQUENCYTIDE(self):
            r"""
            Inverse of the geometrical flattening factor f of the Earth (f = (a-b)/a; zero-frequency-tide value)

            #Source: Numerical value is zero-frequency-tide value from E. Groten, 2000, 'Report of Special Commission 3 of IAG', in Proceedings of IAU Colloquium 180, 'Towards models and constants for sub-microarcsecond astrometry', eds K.J. Johnston, D.D. McCarthy, B.J. Luzum, G.H. Kaplan, page 337<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__EARTHELLIPSOID_INVERSEFLATTENING_ZEROFREQUENCYTIDE


        __EARTHELLIPSOID_JSUB2_ZEROFREQUENCYTIDE: float  = 1.08263590e-03

        @property
        def EARTHELLIPSOID_JSUB2_ZEROFREQUENCYTIDE(self):
            r"""
            Dynamical form-factor of the Earth, i.e., second-degree zonal harmonic of the geopotential including the indirect tidal distortion on J_2, i.e., in the zero-frequency-tide system JGM-3. The (long-term) rate of change of this parameter equals -3.001E-9 cy^-1

            #Source: Numerical value is zero-frequency-tide value from E. Groten, 2000, 'Report of Special Commission 3 of IAG', in Proceedings of IAU Colloquium 180, 'Towards models and constants for sub-microarcsecond astrometry', eds K.J. Johnston, D.D. McCarthy, B.J. Luzum, G.H. Kaplan, page 337. See also the IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__EARTHELLIPSOID_JSUB2_ZEROFREQUENCYTIDE


        __EARTHELLIPSOID_RSUB0: float  = 6.36367256e+06 # [m]

        @property
        def EARTHELLIPSOID_RSUB0(self):
            r"""
            Geopotential scale factor

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__EARTHELLIPSOID_RSUB0


        __EARTHELLIPSOID_SEMIMAJORAXIS_ZEROFREQUENCYTIDE: float  = 6378136.6 # [m]

        @property
        def EARTHELLIPSOID_SEMIMAJORAXIS_ZEROFREQUENCYTIDE(self):
            r"""
            Semi-major axis of the Earth reference ellipsoid (zero-frequency-tide value)

            #Source: Numerical value is zero-frequency-tide value from E. Groten, 2000, 'Report of Special Commission 3 of IAG', in Proceedings of IAU Colloquium 180, 'Towards models and constants for sub-microarcsecond astrometry', eds K.J. Johnston, D.D. McCarthy, B.J. Luzum, G.H. Kaplan, page 337. See also the IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__EARTHELLIPSOID_SEMIMAJORAXIS_ZEROFREQUENCYTIDE


        __EARTHELLIPSOID_SPINRATE_NOMINAL: float  = 7.2921150e-05 # [rad s^-1]

        @property
        def EARTHELLIPSOID_SPINRATE_NOMINAL(self):
            r"""
            Nominal mean angular velocity of the Earth

            #Source: Numerical value is from E. Groten, 2000, 'Report of Special Commission 3 of IAG', in Proceedings of IAU Colloquium 180, 'Towards models and constants for sub-microarcsecond astrometry', eds K.J. Johnston, D.D. McCarthy, B.J. Luzum, G.H. Kaplan, page 337. See also the IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [rad s^-1]
            """

            return self.__EARTHELLIPSOID_SPINRATE_NOMINAL


        __EARTHELLIPSOID_WSUB0: float  = 6.263685600e+07 # [m^2 s^-2]

        @property
        def EARTHELLIPSOID_WSUB0(self):
            r"""
            Potential of the geoid

            #Source: Numerical value is from E. Groten, 2000, 'Report of Special Commission 3 of IAG', in Proceedings of IAU Colloquium 180, 'Towards models and constants for sub-microarcsecond astrometry', eds K.J. Johnston, D.D. McCarthy, B.J. Luzum, G.H. Kaplan, page 337. See also the IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^2 s^-2]
            """

            return self.__EARTHELLIPSOID_WSUB0


        __EARTHSYSTEM_ASTROMETRICSIGNATURE_10PARSEC: float  = 0.304 # [10^-6 arcsec]

        @property
        def EARTHSYSTEM_ASTROMETRICSIGNATURE_10PARSEC(self):
            r"""
            Astrometric signature of the Sun induced by the Earth system for an observer located at a distance of 10 pc from the Sun

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.7, Equation 1.22, page 11<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__EARTHSYSTEM_ASTROMETRICSIGNATURE_10PARSEC


        __EARTHSYSTEM_ORBITALPERIOD: float  = 1.0000174 # [yr]

        @property
        def EARTHSYSTEM_ORBITALPERIOD(self):
            r"""
            Sidereal orbital period

            #Source: Values derived from the mean longitude rates in Table 5.8.1 in P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, page 316; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [yr]
            """

            return self.__EARTHSYSTEM_ORBITALPERIOD


        __EARTHSYSTEM_RADIALVELOCITYSIGNATURE: float  = 0.091 # [m s^-1]

        @property
        def EARTHSYSTEM_RADIALVELOCITYSIGNATURE(self):
            r"""
            Radial-velocity amplitude of the Sun induced by the Earth system for 'an edge-on observer' (i.e., an observer in the orbital plane of the Earth system)

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.4, Equation 1.18, page 9<br/>
            #Basic : false
            #Scalar: true
            #Unit: [m s^-1]
            """

            return self.__EARTHSYSTEM_RADIALVELOCITYSIGNATURE


        __EARTHTOMOON_MASSRATIO: float  = 81.30057

        @property
        def EARTHTOMOON_MASSRATIO(self):
            r"""
            Ratio of Earth to Moon mass

            #Basic : false
            #Scalar: true
            #Unit: []
            """

            return self.__EARTHTOMOON_MASSRATIO


        __EARTH_ECCENTRICITY: float  = 8.181930088e-02

        @property
        def EARTH_ECCENTRICITY(self):
            r"""
            Eccentricity e of the Earth (its shape, not its orbit)

            #Basic : false
            #Scalar: true
            #Unit: []
            """

            return self.__EARTH_ECCENTRICITY


        __EARTH_ENCOMPASSINGSPHERERADIUS: float  = 6378137.0 # [m]

        @property
        def EARTH_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around the Earth which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__EARTH_ENCOMPASSINGSPHERERADIUS


        __EARTH_EQUATORIALRADIUS: float  = 6.37813660e+06 # [m]

        @property
        def EARTH_EQUATORIALRADIUS(self):
            r"""
            Equatorial radius of the Earth

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__EARTH_EQUATORIALRADIUS


        __EARTH_EQUATORIALRADIUS_NOMINAL: float  = 6.37810e+06 # [m]

        @property
        def EARTH_EQUATORIALRADIUS_NOMINAL(self):
            r"""
            Nominal equatorial radius of the Earth (zero-frequency-tide value), in SI units. This nominal value shall be understood as conversion factor only

            #Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__EARTH_EQUATORIALRADIUS_NOMINAL


        __EARTH_FLATTENING: float  = 3.352819698e-03

        @property
        def EARTH_FLATTENING(self):
            r"""
            Geometrical flattening factor f of the Earth (f = (a-b)/a); this quantity is also refered to as ellipticity, but is not identical to eccentricity

            #Basic : false
            #Scalar: true
            #Unit: []
            """

            return self.__EARTH_FLATTENING


        __EARTH_FLUXREDUCTION_MAXIMUM: float  = 0.008 # [%]

        @property
        def EARTH_FLUXREDUCTION_MAXIMUM(self):
            r"""
            Maximum reduction of the solar flux for an observer external to the solar system during a transit of Earth

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.2, Equation 2.4, page 14<br/>
            #Basic : false
            #Scalar: true
            #Unit: [%]
            """

            return self.__EARTH_FLUXREDUCTION_MAXIMUM


        __EARTH_GM: float  = 3.9860044180e+14 # [m^3 s^-2]

        @property
        def EARTH_GM(self):
            r"""
            Geocentric gravitational constant (TCB-compatible value), including the Earth's atmosphere but excluding the mass of the Moon

            #Basic : false
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__EARTH_GM


        __EARTH_GM_NOMINAL: float  = 3.9860040e+14 # [m^3 s^-2]

        @property
        def EARTH_GM_NOMINAL(self):
            r"""
            Nominal GM of the Earth, in SI units. This nominal value shall be understood as conversion factor only

            #Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__EARTH_GM_NOMINAL


        __EARTH_GEOMETRICALBEDO: float  = 0.367

        @property
        def EARTH_GEOMETRICALBEDO(self):
            r"""
            Geometric albedo of the Earth (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__EARTH_GEOMETRICALBEDO


        __EARTH_MASS: float  = 5.97237e+24 # [kg]

        @property
        def EARTH_MASS(self):
            r"""
            Earth mass, including its atmosphere but excluding the mass of the Moon

            #Basic : false
            #Scalar: true
            #Unit: [kg]
            """

            return self.__EARTH_MASS


        __EARTH_MASSDENSITY_MEAN: float  = 5.514 # [g cm^-3]

        @property
        def EARTH_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of the Earth

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__EARTH_MASSDENSITY_MEAN


        __EARTH_NORTHROTATIONALPOLE_DECLINATION: float  = 90.00 # [deg]

        @property
        def EARTH_NORTHROTATIONALPOLE_DECLINATION(self):
            r"""
            IAU-recommended value for the declination \delta_0 of the north pole of rotation of Earth. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value shall be used for comparative purposes only

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__EARTH_NORTHROTATIONALPOLE_DECLINATION


        __EARTH_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE: float  = -0.00001525 # [deg day^-1]

        @property
        def EARTH_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the declination \delta_0 of the north pole of rotation of Earth. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value shall be used for comparative purposes only

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : false
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__EARTH_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE


        __EARTH_NORTHROTATIONALPOLE_RIGHTASCENSION: float  = 0.00 # [deg]

        @property
        def EARTH_NORTHROTATIONALPOLE_RIGHTASCENSION(self):
            r"""
            IAU-recommended value for the right ascension \alpha_0 of the north pole of rotation of Earth. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value shall be used for comparative purposes only

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__EARTH_NORTHROTATIONALPOLE_RIGHTASCENSION


        __EARTH_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE: float  = -0.00001755 # [deg day^-1]

        @property
        def EARTH_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the right ascension \alpha_0 of the north pole of rotation of Earth. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value shall be used for comparative purposes only

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : false
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__EARTH_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE


        __EARTH_ORBITALSPEED_MEAN: float  = 2.97827e+04 # [m s^-1]

        @property
        def EARTH_ORBITALSPEED_MEAN(self):
            r"""
            Mean orbital speed of the Earth (mean velocity over an unperturbed elliptic orbit). The equation is accurate to 4-th order in EMBC_OrbitalEccentricity_J2000

            #Basic : false
            #Scalar: true
            #Unit: [m s^-1]
            """

            return self.__EARTH_ORBITALSPEED_MEAN


        __EARTH_POLARRADIUS: float  = 6.35675186e+06 # [m]

        @property
        def EARTH_POLARRADIUS(self):
            r"""
            Mean polar radius of the Earth

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__EARTH_POLARRADIUS


        __EARTH_POLARRADIUS_NOMINAL: float  = 6.35680e+06 # [m]

        @property
        def EARTH_POLARRADIUS_NOMINAL(self):
            r"""
            Nominal polar radius of the Earth (zero-frequency-tide value), in SI units. This nominal value shall be understood as conversion factor only

            #Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__EARTH_POLARRADIUS_NOMINAL


        __EARTH_PRIMEMERIDIAN_EPHEMERISPOSITION: float  = 190.147 # [deg]

        @property
        def EARTH_PRIMEMERIDIAN_EPHEMERISPOSITION(self):
            r"""
            IAU-recommended value for the ephemeris position of the prime meridian of Earth. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde. The numerical value shall be used for comparative purposes only

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__EARTH_PRIMEMERIDIAN_EPHEMERISPOSITION


        __EARTH_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE: float  = 360.9856235 # [deg day^-1]

        @property
        def EARTH_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch of 1.5 January 2000 = JD2451545.0 TDB) of the ephemeris position of the prime meridian of Earth. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde. The numerical value shall be used for comparative purposes only

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__EARTH_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE


        __EARTH_RADIUS_MEAN: float  = 6.371008e+06 # [m]

        @property
        def EARTH_RADIUS_MEAN(self):
            r"""
            Mean radius of the Earth

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__EARTH_RADIUS_MEAN


        __EARTH_SURFACEGRAVITY: float  = 9.798 # [m s^-2]

        @property
        def EARTH_SURFACEGRAVITY(self):
            r"""
            Surface gravity of the Earth

            #Basic : false
            #Scalar: true
            #Unit: [m s^-2]
            """

            return self.__EARTH_SURFACEGRAVITY


        __EARTH_SURFACEGRAVITY_MEAN: float  = 9.784 # [m s^-2]

        @property
        def EARTH_SURFACEGRAVITY_MEAN(self):
            r"""
            Mean surface gravity of the Earth. The value for the International Standard Atmopshere is 9.80665 m s^-2

            #Source: F. Budnik (ESA), 8 March 2013, 'Gaia FDS-SOC Orbit ICD', GAIA-ESC-ICD-0012, issue 2, revision 0, Annex A. Reference documents: F. Kleijer (Netherlands Geodetic Commission, Delft), 1 April 2004, 'Troposphere Modeling and Filtering for Precise GPS Leveling', Publications on Geodesy 56, ISBN 90 6132 284 7 (http://www.ncg.knaw.nl/Publicaties/Geodesy/pdf/56Kleijer.pdf and http://repository.tudelft.nl/view/ir/uuid%3Aea1f0cf0-4e48-421b-b7ae-4ae3e36d1880/), J. Saastamoinen, 1 January 1972, 'Atmospheric correction for the troposphere and stratosphere in radio ranging of satellites' in 'The use of artificial satellites for geodesy', editors S.W. Henrikson et al., Geophysical Monograph Series, 15, 247-251, and B.R. Bean and E.J. Dutton, 1 March 1966, 'Radio Meteorology', National Bureau of Standards Monograph, 92, 1-44<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m s^-2]
            """

            return self.__EARTH_SURFACEGRAVITY_MEAN


        __EARTH_TRANSITPROBABILITY: float  = 0.469 # [%]

        @property
        def EARTH_TRANSITPROBABILITY(self):
            r"""
            Geometric transit probability (Earth transiting the Sun) for an observer external to the solar system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.1, Equation 2.2, page 14<br/>
            #Basic : false
            #Scalar: true
            #Unit: [%]
            """

            return self.__EARTH_TRANSITPROBABILITY


        __EARTH_TRANSITTIME_MAXIMUM: float  = 0.55 # [day]

        @property
        def EARTH_TRANSITTIME_MAXIMUM(self):
            r"""
            Maximum transit time of Earth (transiting the Sun) for an observer external to the solar system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.3, Equation 2.5, page 15<br/>
            #Basic : false
            #Scalar: true
            #Unit: [day]
            """

            return self.__EARTH_TRANSITTIME_MAXIMUM


        __EARTH_VONEZEROMAGNITUDE: float  = -3.86 # [mag]

        @property
        def EARTH_VONEZEROMAGNITUDE(self):
            r"""
            V(1,0) magnitude of the Earth (i.e., the visual magnitude of the planet reduced to a distance of 1 au from both the Sun and Earth and phase angle zero). This parameter is also refered to as absolute magnitude in planetary sciences

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__EARTH_VONEZEROMAGNITUDE


        __EARTH_VOLUMETRICRADIUS: float  = 6.371000e+06 # [m]

        @property
        def EARTH_VOLUMETRICRADIUS(self):
            r"""
            Mean volumetric radius of the Earth

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__EARTH_VOLUMETRICRADIUS


        __ELECTRIC_CONSTANT: float  = 8.854187817620390e-12 # [F m^-1]

        @property
        def ELECTRIC_CONSTANT(self):
            r"""
            Electric constant (defining constant)

            #Basic : false
            #Scalar: true
            #Unit: [F m^-1]
            """

            return self.__ELECTRIC_CONSTANT


        __ELECTRON_CLASSICALRADIUS: float  = 2.81794032201e-15 # [m]

        @property
        def ELECTRON_CLASSICALRADIUS(self):
            r"""
            Classical electron radius. Note: best-measured value equals 2.8179403227E-15 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__ELECTRON_CLASSICALRADIUS


        __ELECTRON_CROSSSECTION_THOMSONSCATTERING: float  = 6.65245871237e-29 # [m^2]

        @property
        def ELECTRON_CROSSSECTION_THOMSONSCATTERING(self):
            r"""
            Thomson free-electron-scattering absorption coefficient (cross section per electron). Note: best-measured value equals 0.66524587158E-28 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))

            #Source: E.g., R. Kippenhahn, A. Weigert, 1991, 'Stellar structure and evolution' (corrected 2-nd printing), Springer Verlag, Berlin, Section 17, Equation 17.1, page 137<br/>
            #Basic : false
            #Scalar: true
            #Unit: [m^2]
            """

            return self.__ELECTRON_CROSSSECTION_THOMSONSCATTERING


        __ELECTRON_MASS: float  = 9.109383560e-31 # [kg]

        @property
        def ELECTRON_MASS(self):
            r"""
            Electron mass

            #Source: P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [kg]
            """

            return self.__ELECTRON_MASS


        __ELEMENTARYCHARGE_CONSTANT: float  = 1.60217662080e-19 # [C]

        @property
        def ELEMENTARYCHARGE_CONSTANT(self):
            r"""
            Elementary charge

            #Source: P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [C]
            """

            return self.__ELEMENTARYCHARGE_CONSTANT


        __ERG_JOULE: float  = 1.0e-07 # [J]

        @property
        def ERG_JOULE(self):
            r"""
            One erg expressed in units of J. Note that 'erg' is a non-SI unit which should not be used

            #Basic : true
            #Scalar: true
            #Unit: [J]
            """

            return self.__ERG_JOULE


        __EXTINCTION_CURVE: str  = "Nature/Extinction_Curve_002.fits"

        @property
        def EXTINCTION_CURVE(self):
            r"""
            Interstellar extinction curve (see also J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1). First column: wavelength \lambda (in nm; \lambda = 200.0 to 1100.0). Second column: normalised interstellar extinction A(\lambda) / A(\lambda_{ref}) with \lambda_{ref} = 1000 / 1.82 = 549.45 nm. Note that it has become customary within the Gaia community to denote A(\lambda_{ref}) for \lambda_{ref} = 1000 / 1.82 = 549.45 nm as A(550 nm)

            #Source: J.A. Cardelli, G.C. Clayton, J.S. Mathis, 1989, 'The relationship between infrared, optical, and ultraviolet extinction', Astrophysical Journal (ApJ), 345, 245<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__EXTINCTION_CURVE


        __EXTINCTION_GALACTICPLANE_TYPICAL: float  = 1.0 # [mag kpc^-1]

        @property
        def EXTINCTION_GALACTICPLANE_TYPICAL(self):
            r"""
            Typical extinction in the Johnson V band (A_V) per kpc in the Galactic plane; values ranging from 0.5 to 1.5 mag kpc^-1 are considered 'normal' (e.g., H. Jonch-Sorensen, 1994, 'CCD uvby-beta photometry of faint stars. 2: Reddening in six fields in the Galaxy', A&A, 292, 92)

            #Source: Typical value ('common lore')<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag kpc^-1]
            """

            return self.__EXTINCTION_GALACTICPLANE_TYPICAL


        __EXTINCTION_TOTALTOSELECTIVERATIO: float  = 3.10

        @property
        def EXTINCTION_TOTALTOSELECTIVERATIO(self):
            r"""
            Ratio of total to selective absorption (typical value). One has: R_V = A_V / E(B-V), where A_V is the total extinction in the Johnson V band and E(B-V) is the colour excess (see also J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)

            #Source: J.A. Cardelli, G.C. Clayton, J.S. Mathis, 1989, 'The relationship between infrared, optical, and ultraviolet extinction', Astrophysical Journal (ApJ), 345, 245<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__EXTINCTION_TOTALTOSELECTIVERATIO


        __F2VSTAR_SPECTRUM_NUMBEROFPHOTONS: str  = "Nature/F2VStar_Spectrum_NumberOfPhotons_001.fits"

        @property
        def F2VSTAR_SPECTRUM_NUMBEROFPHOTONS(self):
            r"""
            Photon flux density N_{\lambda}(\lambda) of an unreddened F2V star (Pickles' star number 015) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)

            #Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__F2VSTAR_SPECTRUM_NUMBEROFPHOTONS


        __F6VSTAR_SPECTRUM_NUMBEROFPHOTONS: str  = "Nature/F6VStar_Spectrum_NumberOfPhotons_001.fits"

        @property
        def F6VSTAR_SPECTRUM_NUMBEROFPHOTONS(self):
            r"""
            Photon flux density N_{\lambda}(\lambda) of an unreddened F6V star (Pickles' star number 018) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)

            #Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__F6VSTAR_SPECTRUM_NUMBEROFPHOTONS


        __F8VSTAR_SPECTRUM_NUMBEROFPHOTONS: str  = "Nature/F8VStar_Spectrum_NumberOfPhotons_001.fits"

        @property
        def F8VSTAR_SPECTRUM_NUMBEROFPHOTONS(self):
            r"""
            Photon flux density N_{\lambda}(\lambda) of an unreddened F8V star (Pickles' star number 020) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)

            #Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__F8VSTAR_SPECTRUM_NUMBEROFPHOTONS


        __FILTERTRANSMISSIONCURVE_JOHNSONCOUSINSB: str  = "Nature/FilterTransmissionCurve_JohnsonCousinsB_002.fits"

        @property
        def FILTERTRANSMISSIONCURVE_JOHNSONCOUSINSB(self):
            r"""
            Johnson B band filter profile (normalised standard passband in the Johnson-Cousins UBVRI photometric system). First column: wavelength \lambda (in nm; from 360.0 to 560.0). Second column: normalised transmittance

            #Source: M.S. Bessell, S. Murphy, 2012, 'Spectrophotometric Libraries, Revised Photonic Passbands, and Zero Points for UBVRI, Hipparcos, and Tycho Photometry', PASP, 124, 140<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__FILTERTRANSMISSIONCURVE_JOHNSONCOUSINSB


        __FILTERTRANSMISSIONCURVE_JOHNSONCOUSINSI: str  = "Nature/FilterTransmissionCurve_JohnsonCousinsI_002.fits"

        @property
        def FILTERTRANSMISSIONCURVE_JOHNSONCOUSINSI(self):
            r"""
            Cousins I band filter profile (normalised standard passband in the Johnson-Cousins UBVRI photometric system). First column: wavelength \lambda (in nm; from 700.0 to 920.0). Second column: normalised transmittance

            #Source: M.S. Bessell, S. Murphy, 2012, 'Spectrophotometric Libraries, Revised Photonic Passbands, and Zero Points for UBVRI, Hipparcos, and Tycho Photometry', PASP, 124, 140<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__FILTERTRANSMISSIONCURVE_JOHNSONCOUSINSI


        __FILTERTRANSMISSIONCURVE_JOHNSONCOUSINSR: str  = "Nature/FilterTransmissionCurve_JohnsonCousinsR_002.fits"

        @property
        def FILTERTRANSMISSIONCURVE_JOHNSONCOUSINSR(self):
            r"""
            Cousins R band filter profile (normalised standard passband in the Johnson-Cousins UBVRI photometric system). First column: wavelength \lambda (in nm; from 550.0 to 910.0). Second column: normalised transmittance

            #Source: M.S. Bessell, S. Murphy, 2012, 'Spectrophotometric Libraries, Revised Photonic Passbands, and Zero Points for UBVRI, Hipparcos, and Tycho Photometry', PASP, 124, 140<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__FILTERTRANSMISSIONCURVE_JOHNSONCOUSINSR


        __FILTERTRANSMISSIONCURVE_JOHNSONCOUSINSV: str  = "Nature/FilterTransmissionCurve_JohnsonCousinsV_002.fits"

        @property
        def FILTERTRANSMISSIONCURVE_JOHNSONCOUSINSV(self):
            r"""
            Johnson V band filter profile (normalised standard passband in the Johnson-Cousins UBVRI photometric system). First column: wavelength \lambda (in nm; from 470.0 to 740.0). Second column: normalised transmittance

            #Source: M.S. Bessell, S. Murphy, 2012, 'Spectrophotometric Libraries, Revised Photonic Passbands, and Zero Points for UBVRI, Hipparcos, and Tycho Photometry', PASP, 124, 140<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__FILTERTRANSMISSIONCURVE_JOHNSONCOUSINSV


        __FINESTRUCTURE_CONSTANT: float  = 7.29735256621e-03

        @property
        def FINESTRUCTURE_CONSTANT(self):
            r"""
            Fine structure constant. Note: best-measured value equals 7.2973525664E-3 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))

            #Basic : false
            #Scalar: true
            #Unit: []
            """

            return self.__FINESTRUCTURE_CONSTANT


        __FORESHORTENING_CONSTANT: float  = 1.0227121650e-09 # [mas^-1 km^-1 yr^-1 s]

        @property
        def FORESHORTENING_CONSTANT(self):
            r"""
            Foreshortening constant A_z^-1 (see ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, page 25, Table 1.2.2)

            #Basic : false
            #Scalar: true
            #Unit: [mas^-1 km^-1 yr^-1 s]
            """

            return self.__FORESHORTENING_CONSTANT


        __FOURPISTERADIAN_DEGREESQUARE: float  = 41252.9612494193 # [deg^2]

        @property
        def FOURPISTERADIAN_DEGREESQUARE(self):
            r"""
            Surface area of unit sphere (4 Pi steradians) in units of square degrees

            #Basic : false
            #Scalar: true
            #Unit: [deg^2]
            """

            return self.__FOURPISTERADIAN_DEGREESQUARE


        __FUSEDSILICA_LINEARTHERMALCOEFFICIENTOFEXPANSION_293K: float  = 0.5 # [ppm K^-1]

        @property
        def FUSEDSILICA_LINEARTHERMALCOEFFICIENTOFEXPANSION_293K(self):
            r"""
            Linear thermal expansion coefficient of synthetic fused Silica at 293 K. The quoted value is specified to be valid over the temperature range T = 293 - 373 K; a value for T = 170 K is not available

            #Source: Schott Lithotec AG, 5 August 2010, 'Optical glass data sheets (Lithosil-Q)', http://www.schott.com/advanced_optics/english/download/- schott_optical_glass_august_2010_en.pdf; see also http://www.schott.com/advanced_optics/english/download/- schott_fused_silica_jan_2010_en_brochure.pdf<br/>
            #Basic : true
            #Scalar: true
            #Unit: [ppm K^-1]
            """

            return self.__FUSEDSILICA_LINEARTHERMALCOEFFICIENTOFEXPANSION_293K


        __FUSEDSILICA_MASSDENSITY: float  = 2.2 # [g cm^-3]

        @property
        def FUSEDSILICA_MASSDENSITY(self):
            r"""
            Density of synthetic fused Silica

            #Source: Schott Lithotec AG, 5 August 2010, 'Optical glass data sheets (Lithosil-Q)', http://www.schott.com/advanced_optics/english/download/- schott_optical_glass_august_2010_en.pdf; see also http://www.schott.com/advanced_optics/english/download/- schott_fused_silica_jan_2010_en_brochure.pdf<br/>
            #Basic : true
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__FUSEDSILICA_MASSDENSITY


        __FUSEDSILICA_SELLMEIERCOEFFICIENT_B1: float  = 1.10053898145

        @property
        def FUSEDSILICA_SELLMEIERCOEFFICIENT_B1(self):
            r"""
            Sellmeier coefficient B_1 of synthetic fused Silica, which is dimensionless, at T = 120 K. The Sellmeier equation is an empirical relation between the refractive index n and wavelength \lambda for transparent media in the form: n^2(\lambda) = 1 + \frac{B_1 \lambda^2}{\lambda^2 - C_1} + \frac{B_2 \lambda^2}{\lambda^2 - C_2} + \frac{B_3 \lambda^2}{\lambda^2 - C_3}, where B_1, B_2, B_3 and C_1, C_2, C_3 are experimentally-determined Sellmeier coefficients. These coefficients are defined for \lambda measured in 10^-6 m. The wavelength \lambda is the vacuum wavelength and not that in the material itself, which is \lambda / n(\lambda)

            #Source: EADS-Astrium, 4 July 2008, 'Code V Models Description', GAIA.ASF.TCN.PLM.00287, issue 1, revision 0. For temperature-dependent Sellmeier coefficients, see G. Ghosh, M. Endo, T. Iwasaki, 1 August 1994, 'Temperature-dependent Sellmeier coefficients and chromatic dispersions for some optical fiber glasses', Journal of Lightwave Technology, Volume 12, Number 8, pages 1338-1342 (1994JLwT...12.1338G) and http://www.schott.com/advanced_optics/english/download/- schott_tie-29_refractive_index_v3_jan_2007_en.pdf<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__FUSEDSILICA_SELLMEIERCOEFFICIENT_B1


        __FUSEDSILICA_SELLMEIERCOEFFICIENT_B2: float  = 0.00144043087

        @property
        def FUSEDSILICA_SELLMEIERCOEFFICIENT_B2(self):
            r"""
            Sellmeier coefficient B_2 of synthetic fused Silica, which is dimensionless, at T = 120 K. The Sellmeier equation is an empirical relation between the refractive index n and wavelength \lambda for transparent media in the form: n^2(\lambda) = 1 + \frac{B_1 \lambda^2}{\lambda^2 - C_1} + \frac{B_2 \lambda^2}{\lambda^2 - C_2} + \frac{B_3 \lambda^2}{\lambda^2 - C_3}, where B_1, B_2, B_3 and C_1, C_2, C_3 are experimentally-determined Sellmeier coefficients. These coefficients are defined for \lambda measured in 10^-6 m. The wavelength \lambda is the vacuum wavelength and not that in the material itself, which is \lambda / n(\lambda)

            #Source: EADS-Astrium, 4 July 2008, 'Code V Models Description', GAIA.ASF.TCN.PLM.00287, issue 1, revision 0. For temperature-dependent Sellmeier coefficients, see G. Ghosh, M. Endo, T. Iwasaki, 1 August 1994, 'Temperature-dependent Sellmeier coefficients and chromatic dispersions for some optical fiber glasses', Journal of Lightwave Technology, Volume 12, Number 8, pages 1338-1342 (1994JLwT...12.1338G) and http://www.schott.com/advanced_optics/english/download/- schott_tie-29_refractive_index_v3_jan_2007_en.pdf<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__FUSEDSILICA_SELLMEIERCOEFFICIENT_B2


        __FUSEDSILICA_SELLMEIERCOEFFICIENT_B3: float  = 0.77782846144

        @property
        def FUSEDSILICA_SELLMEIERCOEFFICIENT_B3(self):
            r"""
            Sellmeier coefficient B_3 of synthetic fused Silica, which is dimensionless, at T = 120 K. The Sellmeier equation is an empirical relation between the refractive index n and wavelength \lambda for transparent media in the form: n^2(\lambda) = 1 + \frac{B_1 \lambda^2}{\lambda^2 - C_1} + \frac{B_2 \lambda^2}{\lambda^2 - C_2} + \frac{B_3 \lambda^2}{\lambda^2 - C_3}, where B_1, B_2, B_3 and C_1, C_2, C_3 are experimentally-determined Sellmeier coefficients. These coefficients are defined for \lambda measured in 10^-6 m. The wavelength \lambda is the vacuum wavelength and not that in the material itself, which is \lambda / n(\lambda)

            #Source: EADS-Astrium, 4 July 2008, 'Code V Models Description', GAIA.ASF.TCN.PLM.00287, issue 1, revision 0. For temperature-dependent Sellmeier coefficients, see G. Ghosh, M. Endo, T. Iwasaki, 1 August 1994, 'Temperature-dependent Sellmeier coefficients and chromatic dispersions for some optical fiber glasses', Journal of Lightwave Technology, Volume 12, Number 8, pages 1338-1342 (1994JLwT...12.1338G) and http://www.schott.com/advanced_optics/english/download/- schott_tie-29_refractive_index_v3_jan_2007_en.pdf<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__FUSEDSILICA_SELLMEIERCOEFFICIENT_B3


        __FUSEDSILICA_SELLMEIERCOEFFICIENT_C1: float  = 0.00787874390 # [(10^-6 m)^2]

        @property
        def FUSEDSILICA_SELLMEIERCOEFFICIENT_C1(self):
            r"""
            Sellmeier coefficient C_1 of synthetic fused Silica, in units of (10^-6 m)^2, at T = 120 K. It is emphasised that this C-value is already squared, thus complying with the denominator units of the Sellmeier equation, in (10^-6 m)^2. The Sellmeier equation is an empirical relation between the refractive index n and wavelength \lambda for transparent media in the form: n^2(\lambda) = 1 + \frac{B_1 \lambda^2}{\lambda^2 - C_1} + \frac{B_2 \lambda^2}{\lambda^2 - C_2} + \frac{B_3 \lambda^2}{\lambda^2 - C_3}, where B_1, B_2, B_3 and C_1, C_2, C_3 are experimentally-determined Sellmeier coefficients. These coefficients are defined for \lambda measured in 10^-6 m. The wavelength \lambda is the vacuum wavelength and not that in the material itself, which is \lambda / n(\lambda)

            #Source: EADS-Astrium, 4 July 2008, 'Code V Models Description', GAIA.ASF.TCN.PLM.00287, issue 1, revision 0. For temperature-dependent Sellmeier coefficients, see G. Ghosh, M. Endo, T. Iwasaki, 1 August 1994, 'Temperature-dependent Sellmeier coefficients and chromatic dispersions for some optical fiber glasses', Journal of Lightwave Technology, Volume 12, Number 8, pages 1338-1342 (1994JLwT...12.1338G) and http://www.schott.com/advanced_optics/english/download/- schott_tie-29_refractive_index_v3_jan_2007_en.pdf<br/>
            #Basic : false
            #Scalar: true
            #Unit: [(10^-6 m)^2]
            """

            return self.__FUSEDSILICA_SELLMEIERCOEFFICIENT_C1


        __FUSEDSILICA_SELLMEIERCOEFFICIENT_C2: float  = 0.07320427965 # [(10^-6 m)^2]

        @property
        def FUSEDSILICA_SELLMEIERCOEFFICIENT_C2(self):
            r"""
            Sellmeier coefficient C_2 of synthetic fused Silica, in units of (10^-6 m)^2, at T = 120 K. It is emphasised that this C-value is already squared, thus complying with the denominator units of the Sellmeier equation, in (10^-6 m)^2. The Sellmeier equation is an empirical relation between the refractive index n and wavelength \lambda for transparent media in the form: n^2(\lambda) = 1 + \frac{B_1 \lambda^2}{\lambda^2 - C_1} + \frac{B_2 \lambda^2}{\lambda^2 - C_2} + \frac{B_3 \lambda^2}{\lambda^2 - C_3}, where B_1, B_2, B_3 and C_1, C_2, C_3 are experimentally-determined Sellmeier coefficients. These coefficients are defined for \lambda measured in 10^-6 m. The wavelength \lambda is the vacuum wavelength and not that in the material itself, which is \lambda / n(\lambda)

            #Source: EADS-Astrium, 4 July 2008, 'Code V Models Description', GAIA.ASF.TCN.PLM.00287, issue 1, revision 0. For temperature-dependent Sellmeier coefficients, see G. Ghosh, M. Endo, T. Iwasaki, 1 August 1994, 'Temperature-dependent Sellmeier coefficients and chromatic dispersions for some optical fiber glasses', Journal of Lightwave Technology, Volume 12, Number 8, pages 1338-1342 (1994JLwT...12.1338G) and http://www.schott.com/advanced_optics/english/download/- schott_tie-29_refractive_index_v3_jan_2007_en.pdf<br/>
            #Basic : false
            #Scalar: true
            #Unit: [(10^-6 m)^2]
            """

            return self.__FUSEDSILICA_SELLMEIERCOEFFICIENT_C2


        __FUSEDSILICA_SELLMEIERCOEFFICIENT_C3: float  = 85.64043680329 # [(10^-6 m)^2]

        @property
        def FUSEDSILICA_SELLMEIERCOEFFICIENT_C3(self):
            r"""
            Sellmeier coefficient C_3 of synthetic fused Silica, in units of (10^-6 m)^2, at T = 120 K. It is emphasised that this C-value is already squared, thus complying with the denominator units of the Sellmeier equation, in (10^-6 m)^2. The Sellmeier equation is an empirical relation between the refractive index n and wavelength \lambda for transparent media in the form: n^2(\lambda) = 1 + \frac{B_1 \lambda^2}{\lambda^2 - C_1} + \frac{B_2 \lambda^2}{\lambda^2 - C_2} + \frac{B_3 \lambda^2}{\lambda^2 - C_3}, where B_1, B_2, B_3 and C_1, C_2, C_3 are experimentally-determined Sellmeier coefficients. These coefficients are defined for \lambda measured in 10^-6 m. The wavelength \lambda is the vacuum wavelength and not that in the material itself, which is \lambda / n(\lambda)

            #Source: EADS-Astrium, 4 July 2008, 'Code V Models Description', GAIA.ASF.TCN.PLM.00287, issue 1, revision 0. For temperature-dependent Sellmeier coefficients, see G. Ghosh, M. Endo, T. Iwasaki, 1 August 1994, 'Temperature-dependent Sellmeier coefficients and chromatic dispersions for some optical fiber glasses', Journal of Lightwave Technology, Volume 12, Number 8, pages 1338-1342 (1994JLwT...12.1338G) and http://www.schott.com/advanced_optics/english/download/- schott_tie-29_refractive_index_v3_jan_2007_en.pdf<br/>
            #Basic : false
            #Scalar: true
            #Unit: [(10^-6 m)^2]
            """

            return self.__FUSEDSILICA_SELLMEIERCOEFFICIENT_C3


        __FUSEDSILICA_TRANSMISSIVITY_10MM: str  = "Nature/FusedSilica_Transmissivity_10mm_001.fits"

        @property
        def FUSEDSILICA_TRANSMISSIVITY_10MM(self):
            r"""
            Typical transmission of synthetic fused Silica, including Fresnel reflection losses for an uncoated surface, for a 10-mm path length. First column: wavelength \lambda (in nm; from 200.0 to 1250.0). Second column: typical transmission. Explanatory, supplementary information: this parameter provides the bulk transmission of 10 mm of synthetic fused Silica, including Fresnel reflection losses. Fresnel diffraction losses, however, are only applicable in the absence of an anti-relfection coating. Since all Gaia prisms (BP/RP and RVS) do have anti-relfection coatings, the bulk transmission (i) first has to be corrected for (the absence of) Fresnel diffraction losses, and (ii) subsequently has to be scaled for the proper path length in the fused Silica. Ad (i): the Fresnel reflectivity per surface equals R = (n-1)^2 (n+1)^-2, where n is the index of refraction of fused Silica (which is a function of wavelength and temperature). The Fresnel loss per surface is hence 1 - R. The index of refraction n can be calculated from the Sellmeier coefficients (see parameters :Nature:FusedSilica_SellmeierCoefficient_*). The corrected bulk transmission of 10 mm of fused Silica exceeds 99.9% above 250 nm. Ad (ii): let us denote the transmission curve of 10 mm synthetic fused Silica, corrected for (the absence of) Fresnel diffraction losses, by \eta_10. The transmission curve for d mm path length (see parameters :Satellite:BP:Prism_Thickness, :Satellite:RP:Prism_Thickness, and :Satellite:RVS:Prism_Thickness) can then be calculated - following the Bouguer-Lambert law - as \eta_d = \eta_10^(d/10). As secondary effect, prism wedge angles could be included, effectively increasing the path lengths d

            #Source: Schott Lithotec AG, 5 August 2010, 'Optical glass data sheets (Lithosil-Q)', http://www.schott.com/advanced_optics/english/download/- schott_optical_glass_august_2010_en.pdf; see also http://www.schott.com/advanced_optics/english/download/- schott_fused_silica_jan_2010_en_brochure.pdf<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__FUSEDSILICA_TRANSMISSIVITY_10MM


        __G2VSTAR_SPECTRUM_NUMBEROFPHOTONS: str  = "Nature/G2VStar_Spectrum_NumberOfPhotons_001.fits"

        @property
        def G2VSTAR_SPECTRUM_NUMBEROFPHOTONS(self):
            r"""
            Photon flux density N_{\lambda}(\lambda) of an unreddened G2V star (Pickles' star number 026) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)

            #Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__G2VSTAR_SPECTRUM_NUMBEROFPHOTONS


        __G2VSTAR_SPECTRUM_NUMBEROFPHOTONSHIGHRESOLUTION: str  = "Nature/G2VStar_Spectrum_NumberOfPhotonsHighResolution_001.fits"

        @property
        def G2VSTAR_SPECTRUM_NUMBEROFPHOTONSHIGHRESOLUTION(self):
            r"""
            High-resolution photon-flux density N_{\lambda}(\lambda) of an unreddened G2V star at V = 15 mag. The data refer to a high-resolution Kurucz-model spectrum with the following properties: effective temperature T_eff = 5800 K, logarithm of surface gravity log g = 4.5, metallicity [Fe/H] = 0.0, alpha-elements [\alpha/Fe] = 0.0, rotational velocity v sini = 5 km s^-1, micro-turbulence velocity = 2.0 km s^-1, length of convective bubble divided by pressure scale height = 0.50, no convective overshooting, macro-turbulence velocity = 2.0 km s^-1, and resolving power R = \lambda / \delta \lambda = 250,000. First column: wavelength \lambda (in nm; from 830.1673264 to 889.8217922). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1). The 34698 lines have an average wavelength step of 0.00172 nm; the spectrum extent is thus 59.7 nm

            #Source: ESA, 20 June 2005, 'Photon-flux distributions for reference stars', GAIA-EST-TN-00539, issue 1, revision 0, based on D. Katz, priv. comm., 11 May 2005<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__G2VSTAR_SPECTRUM_NUMBEROFPHOTONSHIGHRESOLUTION


        __G8IIISTAR_SPECTRUM_NUMBEROFPHOTONS: str  = "Nature/G8IIIStar_Spectrum_NumberOfPhotons_001.fits"

        @property
        def G8IIISTAR_SPECTRUM_NUMBEROFPHOTONS(self):
            r"""
            Photon flux density N_{\lambda}(\lambda) of an unreddened G8III star (Pickles' star number 076) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)

            #Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__G8IIISTAR_SPECTRUM_NUMBEROFPHOTONS


        __HUBBLE_CONSTANT: float  = 69.32 # [km s^-1 Mpc^-1]

        @property
        def HUBBLE_CONSTANT(self):
            r"""
            Hubble constant (uncertainty is 0.80 km s^-1 Mpc^-1)

            #Source: C.L. Bennett, et al., 1 October 2013, 'Nine-Year Wilkinson Microwave Anisotropy Probe (WMAP) Observations: Final Maps and Results', Astrophysical Journal Supplement, Volume 208, 20<br/>
            #Basic : true
            #Scalar: true
            #Unit: [km s^-1 Mpc^-1]
            """

            return self.__HUBBLE_CONSTANT


        class IAU2000A:

            __PRECESSIONNUTATIONMODEL_ARGUMENTMULTIPLIERS_LUNISOLAR: str  = "Nature/IAU2000A/PrecessionNutationModel_ArgumentMultipliers_LuniSolar_001.fits"

            @property
            def PRECESSIONNUTATIONMODEL_ARGUMENTMULTIPLIERS_LUNISOLAR(self):
                r"""
                The IAU 2000A Precession-Nutation Model, developed by Mathews et al. (2002; MHB), is given by a series for nutation in longitude (\Delta\psi) and obliquity (\Delta\epsilon) - referred to the mean ecliptic of date, with time measured in Julian centuries of TDB from epoch J2000.0 - plus the contribution of the corrections to the IAU 1976 precession rates, plus the frame bias in longitude and obliquity. The 'total nutation' includes all components, with the exception of the free core nutation (FCN). The nutation series - providing the direction of the celestial pole in the GCRS with an accuracy of 0.2 mas - includes N_k = 678 luni-solar terms and N_k = 687 planetary terms, which are expressed as 'in-phase' components (A_k, A^\prime_k, B_k, and B^\prime_k) and 'out-of-phase' components (A^\prime\prime_k, A^\prime\prime\prime_k, B^\prime\prime_k, and B^\prime\prime\prime_k) with their time variations: \Delta\psi = Sum_{k=1}^{N_k} (A_k + A^\prime_k * t) * SIN(ARGUMENT) + (A^\prime\prime_k + A^\prime\prime\prime_k * t) * COS(ARGUMENT) and \Delta\epsilon = Sum_{k=1}^{N_k} (B_k + B^\prime_k * t) * COS(ARGUMENT) + (B^\prime\prime_k + B^\prime\prime\prime_k * t) * SIN(ARGUMENT). Each of the N_k = 678 luni-solar terms in the nutation series is characterised by a set of five integers N_j which determines the ARGUMENT for the term as a linear combination of the five Fundamental Arguments F_j, namely the Delaunay variables (l = mean anomaly of the Moon, l^\prime = mean anomaly of the Sun, F = L - \Omega [with l the mean longitude of the Moon], D = mean elongation of the Moon from the Sun, \Omega = mean longitude of the ascending node the Moon): ARGUMENT = Sum_{j=1}^{5} N_j * F_j, where the values (N_1, ..., N_5) of the multipliers characterise the term. The F_j are functions of time, and the angular frequency of the nutation described by the term is given by \omega = d(ARGUMENT) / dt. The N_k = 687 planetary nutation terms differ from the luni-solar terms described above only in that ARGUMENT = SUM_{j=1}^{14} N^\prime_j * F^\prime_j, where F^\prime_6 through F^\prime_13 are the mean longitudes of the planets Mercury through Neptune (including F^\prime_8 for the Earth), and F^\prime_14 is the general precession in longitude p_a. Over time scales involved in nutation studies, the frequency \omega is effectively time-independent, and one may write, for the k-th term in the nutation series, ARGUMENT = \omega_k + \alpha_k. This parameter provides the j = 1, ..., 5 argument multipliers N_j for the N_k = 678 luni-solar terms. These multipliers are dimensionless

                #Source: Mathews, P.M., Herring, T.A., Buffett, B.A., 2002, 'Modeling of nutation-precession: new nutation series for non-rigid Earth, and insights into the Earth's interior', Journal of Geophysical Research (Solid Earth), Volume 107, Issue B4, 2068, 2002JGRB..107.2068M. Reference document: G. Petit, B. Luzum, 21 October 2010,  'IERS Conventions (2010)', IERS Technical Note 36, Chapter 5 (http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html)<br/>
                #Basic : true
                #Scalar: false
                #Unit: []
                """

                return self.__PRECESSIONNUTATIONMODEL_ARGUMENTMULTIPLIERS_LUNISOLAR


            __PRECESSIONNUTATIONMODEL_ARGUMENTMULTIPLIERS_PLANETARY: str  = "Nature/IAU2000A/PrecessionNutationModel_ArgumentMultipliers_Planetary_001.fits"

            @property
            def PRECESSIONNUTATIONMODEL_ARGUMENTMULTIPLIERS_PLANETARY(self):
                r"""
                The IAU 2000A Precession-Nutation Model, developed by Mathews et al. (2002; MHB), is given by a series for nutation in longitude (\Delta\psi) and obliquity (\Delta\epsilon) - referred to the mean ecliptic of date, with time measured in Julian centuries of TDB from epoch J2000.0 - plus the contribution of the corrections to the IAU 1976 precession rates, plus the frame bias in longitude and obliquity. The 'total nutation' includes all components, with the exception of the free core nutation (FCN). The nutation series - providing the direction of the celestial pole in the GCRS with an accuracy of 0.2 mas - includes N_k = 678 luni-solar terms and N_k = 687 planetary terms, which are expressed as 'in-phase' components (A_k, A^\prime_k, B_k, and B^\prime_k) and 'out-of-phase' components (A^\prime\prime_k, A^\prime\prime\prime_k, B^\prime\prime_k, and B^\prime\prime\prime_k) with their time variations: \Delta\psi = Sum_{k=1}^{N_k} (A_k + A^\prime_k * t) * SIN(ARGUMENT) + (A^\prime\prime_k + A^\prime\prime\prime_k * t) * COS(ARGUMENT) and \Delta\epsilon = Sum_{k=1}^{N_k} (B_k + B^\prime_k * t) * COS(ARGUMENT) + (B^\prime\prime_k + B^\prime\prime\prime_k * t) * SIN(ARGUMENT). Each of the N_k = 678 luni-solar terms in the nutation series is characterised by a set of five integers N_j which determines the ARGUMENT for the term as a linear combination of the five Fundamental Arguments F_j, namely the Delaunay variables (l = mean anomaly of the Moon, l^\prime = mean anomaly of the Sun, F = L - \Omega [with l the mean longitude of the Moon], D = mean elongation of the Moon from the Sun, \Omega = mean longitude of the ascending node the Moon): ARGUMENT = Sum_{j=1}^{5} N_j * F_j, where the values (N_1, ..., N_5) of the multipliers characterise the term. The F_j are functions of time, and the angular frequency of the nutation described by the term is given by \omega = d(ARGUMENT) / dt. The N_k = 687 planetary nutation terms differ from the luni-solar terms described above only in that ARGUMENT = SUM_{j=1}^{14} N^\prime_j * F^\prime_j, where F^\prime_6 through F^\prime_13 are the mean longitudes of the planets Mercury through Neptune (including F^\prime_8 for the Earth), and F^\prime_14 is the general precession in longitude p_a. Over time scales involved in nutation studies, the frequency \omega is effectively time-independent, and one may write, for the k-th term in the nutation series, ARGUMENT = \omega_k + \alpha_k. This parameter provides the j = 1, ..., 14 argument multipliers N^\prime_j for the N_k = 687 planetary terms. These multipliers are dimensionless

                #Source: Mathews, P.M., Herring, T.A., Buffett, B.A., 2002, 'Modeling of nutation-precession: new nutation series for non-rigid Earth, and insights into the Earth's interior', Journal of Geophysical Research (Solid Earth), Volume 107, Issue B4, 2068, 2002JGRB..107.2068M. Reference document: G. Petit, B. Luzum, 21 October 2010,  'IERS Conventions (2010)', IERS Technical Note 36, Chapter 5 (http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html)<br/>
                #Basic : true
                #Scalar: false
                #Unit: []
                """

                return self.__PRECESSIONNUTATIONMODEL_ARGUMENTMULTIPLIERS_PLANETARY


            __PRECESSIONNUTATIONMODEL_COEFFICIENTS_LUNISOLAR: str  = "Nature/IAU2000A/PrecessionNutationModel_Coefficients_LuniSolar_001.fits"

            @property
            def PRECESSIONNUTATIONMODEL_COEFFICIENTS_LUNISOLAR(self):
                r"""
                The IAU 2000A Precession-Nutation Model, developed by Mathews et al. (2002; MHB), is given by a series for nutation in longitude (\Delta\psi) and obliquity (\Delta\epsilon) - referred to the mean ecliptic of date, with time measured in Julian centuries of TDB from epoch J2000.0 - plus the contribution of the corrections to the IAU 1976 precession rates, plus the frame bias in longitude and obliquity. The 'total nutation' includes all components, with the exception of the free core nutation (FCN). The nutation series - providing the direction of the celestial pole in the GCRS with an accuracy of 0.2 mas - includes N_k = 678 luni-solar terms and N_k = 687 planetary terms, which are expressed as 'in-phase' components (A_k, A^\prime_k, B_k, and B^\prime_k) and 'out-of-phase' components (A^\prime\prime_k, A^\prime\prime\prime_k, B^\prime\prime_k, and B^\prime\prime\prime_k) with their time variations: \Delta\psi = Sum_{k=1}^{N_k} (A_k + A^\prime_k * t) * SIN(ARGUMENT) + (A^\prime\prime_k + A^\prime\prime\prime_k * t) * COS(ARGUMENT) and \Delta\epsilon = Sum_{k=1}^{N_k} (B_k + B^\prime_k * t) * COS(ARGUMENT) + (B^\prime\prime_k + B^\prime\prime\prime_k * t) * SIN(ARGUMENT). Each of the N_k = 678 luni-solar terms in the nutation series is characterised by a set of five integers N_j which determines the ARGUMENT for the term as a linear combination of the five Fundamental Arguments F_j, namely the Delaunay variables (l = mean anomaly of the Moon, l^\prime = mean anomaly of the Sun, F = L - \Omega [with l the mean longitude of the Moon], D = mean elongation of the Moon from the Sun, \Omega = mean longitude of the ascending node the Moon): ARGUMENT = Sum_{j=1}^{5} N_j * F_j, where the values (N_1, ..., N_5) of the multipliers characterise the term. The F_j are functions of time, and the angular frequency of the nutation described by the term is given by \omega = d(ARGUMENT) / dt. The N_k = 687 planetary nutation terms differ from the luni-solar terms described above only in that ARGUMENT = SUM_{j=1}^{14} N^\prime_j * F^\prime_j, where F^\prime_6 through F^\prime_13 are the mean longitudes of the planets Mercury through Neptune (including F^\prime_8 for the Earth), and F^\prime_14 is the general precession in longitude p_a. Over time scales involved in nutation studies, the frequency \omega is effectively time-independent, and one may write, for the k-th term in the nutation series, ARGUMENT = \omega_k + \alpha_k. This parameter provides the 4 'in-phase' and 4 'out-of-phase' components (nutation coefficients for longitude and obliquity) for the N_k = 678 luni-solar terms in the order A_k, A^\prime_k, B_k, B^\prime_k, A^\prime\prime_k, A^\prime\prime\prime_k, B^\prime\prime_k, B^\prime\prime\prime_k. Units are 10^-3 arcsec and 10^-3 arcsec per century for the coefficients and their time variations, respectively

                #Source: Mathews, P.M., Herring, T.A., Buffett, B.A., 2002, 'Modeling of nutation-precession: new nutation series for non-rigid Earth, and insights into the Earth's interior', Journal of Geophysical Research (Solid Earth), Volume 107, Issue B4, 2068, 2002JGRB..107.2068M. Reference document: G. Petit, B. Luzum, 21 October 2010,  'IERS Conventions (2010)', IERS Technical Note 36, Chapter 5 (http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html)<br/>
                #Basic : true
                #Scalar: false
                #Unit: []
                """

                return self.__PRECESSIONNUTATIONMODEL_COEFFICIENTS_LUNISOLAR


            __PRECESSIONNUTATIONMODEL_COEFFICIENTS_PLANETARY: str  = "Nature/IAU2000A/PrecessionNutationModel_Coefficients_Planetary_001.fits"

            @property
            def PRECESSIONNUTATIONMODEL_COEFFICIENTS_PLANETARY(self):
                r"""
                The IAU 2000A Precession-Nutation Model, developed by Mathews et al. (2002; MHB), is given by a series for nutation in longitude (\Delta\psi) and obliquity (\Delta\epsilon) - referred to the mean ecliptic of date, with time measured in Julian centuries of TDB from epoch J2000.0 - plus the contribution of the corrections to the IAU 1976 precession rates, plus the frame bias in longitude and obliquity. The 'total nutation' includes all components, with the exception of the free core nutation (FCN). The nutation series - providing the direction of the celestial pole in the GCRS with an accuracy of 0.2 mas - includes N_k = 678 luni-solar terms and N_k = 687 planetary terms, which are expressed as 'in-phase' components (A_k, A^\prime_k, B_k, and B^\prime_k) and 'out-of-phase' components (A^\prime\prime_k, A^\prime\prime\prime_k, B^\prime\prime_k, and B^\prime\prime\prime_k) with their time variations: \Delta\psi = Sum_{k=1}^{N_k} (A_k + A^\prime_k * t) * SIN(ARGUMENT) + (A^\prime\prime_k + A^\prime\prime\prime_k * t) * COS(ARGUMENT) and \Delta\epsilon = Sum_{k=1}^{N_k} (B_k + B^\prime_k * t) * COS(ARGUMENT) + (B^\prime\prime_k + B^\prime\prime\prime_k * t) * SIN(ARGUMENT). Each of the N_k = 678 luni-solar terms in the nutation series is characterised by a set of five integers N_j which determines the ARGUMENT for the term as a linear combination of the five Fundamental Arguments F_j, namely the Delaunay variables (l = mean anomaly of the Moon, l^\prime = mean anomaly of the Sun, F = L - \Omega [with l the mean longitude of the Moon], D = mean elongation of the Moon from the Sun, \Omega = mean longitude of the ascending node the Moon): ARGUMENT = Sum_{j=1}^{5} N_j * F_j, where the values (N_1, ..., N_5) of the multipliers characterise the term. The F_j are functions of time, and the angular frequency of the nutation described by the term is given by \omega = d(ARGUMENT) / dt. The N_k = 687 planetary nutation terms differ from the luni-solar terms described above only in that ARGUMENT = SUM_{j=1}^{14} N^\prime_j * F^\prime_j, where F^\prime_6 through F^\prime_13 are the mean longitudes of the planets Mercury through Neptune (including F^\prime_8 for the Earth), and F^\prime_14 is the general precession in longitude p_a. Over time scales involved in nutation studies, the frequency \omega is effectively time-independent, and one may write, for the k-th term in the nutation series, ARGUMENT = \omega_k + \alpha_k. This parameter provides the 2 'in-phase' and 2 'out-of-phase' components (nutation coefficients for longitude and obliquity) for the N_k = 687 planetary terms in the order A_k, A^\prime\prime_k, B_k, B^\prime\prime_k. Units are 10^-3 arcsec

                #Source: Mathews, P.M., Herring, T.A., Buffett, B.A., 2002, 'Modeling of nutation-precession: new nutation series for non-rigid Earth, and insights into the Earth's interior', Journal of Geophysical Research (Solid Earth), Volume 107, Issue B4, 2068, 2002JGRB..107.2068M. Reference document: G. Petit, B. Luzum, 21 October 2010,  'IERS Conventions (2010)', IERS Technical Note 36, Chapter 5 (http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html)<br/>
                #Basic : true
                #Scalar: false
                #Unit: []
                """

                return self.__PRECESSIONNUTATIONMODEL_COEFFICIENTS_PLANETARY


        __ICRS_LONGITUDEOFASCENDINGNODE: float  = 32.93192 # [deg]

        @property
        def ICRS_LONGITUDEOFASCENDINGNODE(self):
            r"""
            Galactic longitude of the ascending node of the Galactic plane on the equator of the ICRS (numerical value should be regarded as exact)

            #Source: ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, page 91<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__ICRS_LONGITUDEOFASCENDINGNODE


        __ICRS_NORTHGALACTICPOLE_DECLINATION: float  = 27.12825 # [deg]

        @property
        def ICRS_NORTHGALACTICPOLE_DECLINATION(self):
            r"""
            Declination of the north Galactic pole in the ICRS system (numerical value should be regarded as exact)

            #Source: ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, page 91<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__ICRS_NORTHGALACTICPOLE_DECLINATION


        __ICRS_NORTHGALACTICPOLE_RIGHTASCENSION: float  = 192.85948 # [deg]

        @property
        def ICRS_NORTHGALACTICPOLE_RIGHTASCENSION(self):
            r"""
            Right ascension of the north Galactic pole in the ICRS system (numerical value should be regarded as exact)

            #Source: ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, page 91<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__ICRS_NORTHGALACTICPOLE_RIGHTASCENSION


        class INPOP10e:

            __ASTEROID1013TOMBECKA_GM: float  = 7.9531899440686237e+06 # [m^3 s^-2]

            @property
            def ASTEROID1013TOMBECKA_GM(self):
                r"""
                GM of asteroid 1013 Tombecka (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID1013TOMBECKA_GM


            __ASTEROID1021FLAMMARIO_GM: float  = 7.3151076468356557e+07 # [m^3 s^-2]

            @property
            def ASTEROID1021FLAMMARIO_GM(self):
                r"""
                GM of asteroid 1021 Flammario (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID1021FLAMMARIO_GM


            __ASTEROID1036GANYMED_GM: float  = 7.7604237885203171e+06 # [m^3 s^-2]

            @property
            def ASTEROID1036GANYMED_GM(self):
                r"""
                GM of asteroid 1036 Ganymed (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID1036GANYMED_GM


            __ASTEROID105ARTEMIS_GM: float  = 4.0425884340051366e+08 # [m^3 s^-2]

            @property
            def ASTEROID105ARTEMIS_GM(self):
                r"""
                GM of asteroid 105 Artemis (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID105ARTEMIS_GM


            __ASTEROID106DIONE_GM: float  = 5.1363068033089922e+08 # [m^3 s^-2]

            @property
            def ASTEROID106DIONE_GM(self):
                r"""
                GM of asteroid 106 Dione (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID106DIONE_GM


            __ASTEROID107CAMILLA_GM: float  = 4.5300167118439350e+08 # [m^3 s^-2]

            @property
            def ASTEROID107CAMILLA_GM(self):
                r"""
                GM of asteroid 107 Camilla (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID107CAMILLA_GM


            __ASTEROID109FELICITAS_GM: float  = 2.1334036970797304e+07 # [m^3 s^-2]

            @property
            def ASTEROID109FELICITAS_GM(self):
                r"""
                GM of asteroid 109 Felicitas (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID109FELICITAS_GM


            __ASTEROID10HYGIEA_GM: float  = 5.8390041443806954e+09 # [m^3 s^-2]

            @property
            def ASTEROID10HYGIEA_GM(self):
                r"""
                GM of asteroid 10 Hygiea (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID10HYGIEA_GM


            __ASTEROID111ATE_GM: float  = 5.9580019641313180e+08 # [m^3 s^-2]

            @property
            def ASTEROID111ATE_GM(self):
                r"""
                GM of asteroid 111 Ate (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID111ATE_GM


            __ASTEROID112IPHIGENIA_GM: float  = 9.1960926871076204e+07 # [m^3 s^-2]

            @property
            def ASTEROID112IPHIGENIA_GM(self):
                r"""
                GM of asteroid 112 Iphigenia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID112IPHIGENIA_GM


            __ASTEROID117LOMIA_GM: float  = 8.0437920943187555e+08 # [m^3 s^-2]

            @property
            def ASTEROID117LOMIA_GM(self):
                r"""
                GM of asteroid 117 Lomia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID117LOMIA_GM


            __ASTEROID11PARTHENOPE_GM: float  = 5.0048256226049669e+08 # [m^3 s^-2]

            @property
            def ASTEROID11PARTHENOPE_GM(self):
                r"""
                GM of asteroid 11 Parthenope (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID11PARTHENOPE_GM


            __ASTEROID120LACHESIS_GM: float  = 2.0019174836054619e+08 # [m^3 s^-2]

            @property
            def ASTEROID120LACHESIS_GM(self):
                r"""
                GM of asteroid 120 Lachesis (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID120LACHESIS_GM


            __ASTEROID121HERMIONE_GM: float  = 2.0870018202750795e+09 # [m^3 s^-2]

            @property
            def ASTEROID121HERMIONE_GM(self):
                r"""
                GM of asteroid 121 Hermione (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID121HERMIONE_GM


            __ASTEROID126VELLEDA_GM: float  = 2.2017513093567601e+07 # [m^3 s^-2]

            @property
            def ASTEROID126VELLEDA_GM(self):
                r"""
                GM of asteroid 126 Velleda (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID126VELLEDA_GM


            __ASTEROID128NEMESIS_GM: float  = 9.1024212903081425e+08 # [m^3 s^-2]

            @property
            def ASTEROID128NEMESIS_GM(self):
                r"""
                GM of asteroid 128 Nemesis (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID128NEMESIS_GM


            __ASTEROID12VICTORIA_GM: float  = 6.9482386231968752e+07 # [m^3 s^-2]

            @property
            def ASTEROID12VICTORIA_GM(self):
                r"""
                GM of asteroid 12 Victoria (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID12VICTORIA_GM


            __ASTEROID132AETHRA_GM: float  = 1.9253478609675224e+07 # [m^3 s^-2]

            @property
            def ASTEROID132AETHRA_GM(self):
                r"""
                GM of asteroid 132 Aethra (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID132AETHRA_GM


            __ASTEROID134SOPHROSYNE_GM: float  = 1.3455129358325903e+08 # [m^3 s^-2]

            @property
            def ASTEROID134SOPHROSYNE_GM(self):
                r"""
                GM of asteroid 134 Sophrosyne (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID134SOPHROSYNE_GM


            __ASTEROID135HERTHA_GM: float  = 1.2167072819223394e+08 # [m^3 s^-2]

            @property
            def ASTEROID135HERTHA_GM(self):
                r"""
                GM of asteroid 135 Hertha (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID135HERTHA_GM


            __ASTEROID139JUEWA_GM: float  = 2.8269837868259157e+08 # [m^3 s^-2]

            @property
            def ASTEROID139JUEWA_GM(self):
                r"""
                GM of asteroid 139 Juewa (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID139JUEWA_GM


            __ASTEROID13EGERIA_GM: float  = 6.2548618644966828e+08 # [m^3 s^-2]

            @property
            def ASTEROID13EGERIA_GM(self):
                r"""
                GM of asteroid 13 Egeria (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID13EGERIA_GM


            __ASTEROID141LUMEN_GM: float  = 5.5025545915405380e+08 # [m^3 s^-2]

            @property
            def ASTEROID141LUMEN_GM(self):
                r"""
                GM of asteroid 141 Lumen (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID141LUMEN_GM


            __ASTEROID156XANTHIPPE_GM: float  = 4.3257295425004363e+08 # [m^3 s^-2]

            @property
            def ASTEROID156XANTHIPPE_GM(self):
                r"""
                GM of asteroid 156 Xanthippe (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID156XANTHIPPE_GM


            __ASTEROID15EUNOMIA_GM: float  = 2.1020267614940911e+09 # [m^3 s^-2]

            @property
            def ASTEROID15EUNOMIA_GM(self):
                r"""
                GM of asteroid 15 Eunomia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID15EUNOMIA_GM


            __ASTEROID164EVA_GM: float  = 1.9354802800126533e+05 # [m^3 s^-2]

            @property
            def ASTEROID164EVA_GM(self):
                r"""
                GM of asteroid 164 Eva (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID164EVA_GM


            __ASTEROID168SIBYLLA_GM: float  = 7.9887495652663551e+08 # [m^3 s^-2]

            @property
            def ASTEROID168SIBYLLA_GM(self):
                r"""
                GM of asteroid 168 Sibylla (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID168SIBYLLA_GM


            __ASTEROID1694KAISER_GM: float  = 5.1974702708939677e+06 # [m^3 s^-2]

            @property
            def ASTEROID1694KAISER_GM(self):
                r"""
                GM of asteroid 1694 Kaiser (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID1694KAISER_GM


            __ASTEROID16PSYCHE_GM: float  = 1.6739220218687147e+09 # [m^3 s^-2]

            @property
            def ASTEROID16PSYCHE_GM(self):
                r"""
                GM of asteroid 16 Psyche (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID16PSYCHE_GM


            __ASTEROID171OPHELIA_GM: float  = 3.8845612002754150e+08 # [m^3 s^-2]

            @property
            def ASTEROID171OPHELIA_GM(self):
                r"""
                GM of asteroid 171 Ophelia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID171OPHELIA_GM


            __ASTEROID172BAUCIS_GM: float  = 5.9473568787444189e+07 # [m^3 s^-2]

            @property
            def ASTEROID172BAUCIS_GM(self):
                r"""
                GM of asteroid 172 Baucis (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID172BAUCIS_GM


            __ASTEROID173INO_GM: float  = 8.9487138876601197e+08 # [m^3 s^-2]

            @property
            def ASTEROID173INO_GM(self):
                r"""
                GM of asteroid 173 Ino (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID173INO_GM


            __ASTEROID179KLYTAEMNESTRA_GM: float  = 1.1462526623299412e+08 # [m^3 s^-2]

            @property
            def ASTEROID179KLYTAEMNESTRA_GM(self):
                r"""
                GM of asteroid 179 Klytaemnestra (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID179KLYTAEMNESTRA_GM


            __ASTEROID17THETIS_GM: float  = 4.8917509952675321e+08 # [m^3 s^-2]

            @property
            def ASTEROID17THETIS_GM(self):
                r"""
                GM of asteroid 17 Thetis (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID17THETIS_GM


            __ASTEROID185EUNIKE_GM: float  = 5.8481152014919869e+08 # [m^3 s^-2]

            @property
            def ASTEROID185EUNIKE_GM(self):
                r"""
                GM of asteroid 185 Eunike (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID185EUNIKE_GM


            __ASTEROID187LAMBERTA_GM: float  = 1.1232253601955330e+07 # [m^3 s^-2]

            @property
            def ASTEROID187LAMBERTA_GM(self):
                r"""
                GM of asteroid 187 Lamberta (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID187LAMBERTA_GM


            __ASTEROID194PROKNE_GM: float  = 7.4333048018458559e+08 # [m^3 s^-2]

            @property
            def ASTEROID194PROKNE_GM(self):
                r"""
                GM of asteroid 194 Prokne (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID194PROKNE_GM


            __ASTEROID19FORTUNA_GM: float  = 6.4925120044129564e+08 # [m^3 s^-2]

            @property
            def ASTEROID19FORTUNA_GM(self):
                r"""
                GM of asteroid 19 Fortuna (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID19FORTUNA_GM


            __ASTEROID1CERES_GM: float  = 6.2012183942528447e+10 # [m^3 s^-2]

            @property
            def ASTEROID1CERES_GM(self):
                r"""
                GM of asteroid 1 Ceres (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID1CERES_GM


            __ASTEROID200DYNAMENE_GM: float  = 7.6156064781391007e+07 # [m^3 s^-2]

            @property
            def ASTEROID200DYNAMENE_GM(self):
                r"""
                GM of asteroid 200 Dynamene (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID200DYNAMENE_GM


            __ASTEROID204KALLISTO_GM: float  = 2.8036593061651360e+07 # [m^3 s^-2]

            @property
            def ASTEROID204KALLISTO_GM(self):
                r"""
                GM of asteroid 204 Kallisto (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID204KALLISTO_GM


            __ASTEROID209DIDO_GM: float  = 1.0005158700123591e+09 # [m^3 s^-2]

            @property
            def ASTEROID209DIDO_GM(self):
                r"""
                GM of asteroid 209 Dido (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID209DIDO_GM


            __ASTEROID20MASSALIA_GM: float  = 3.8450306541675409e+08 # [m^3 s^-2]

            @property
            def ASTEROID20MASSALIA_GM(self):
                r"""
                GM of asteroid 20 Massalia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID20MASSALIA_GM


            __ASTEROID210ISABELLA_GM: float  = 1.5915077609293021e+08 # [m^3 s^-2]

            @property
            def ASTEROID210ISABELLA_GM(self):
                r"""
                GM of asteroid 210 Isabella (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID210ISABELLA_GM


            __ASTEROID211ISOLDA_GM: float  = 7.1779392920251512e+08 # [m^3 s^-2]

            @property
            def ASTEROID211ISOLDA_GM(self):
                r"""
                GM of asteroid 211 Isolda (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID211ISOLDA_GM


            __ASTEROID212MEDEA_GM: float  = 3.1495425151492577e+08 # [m^3 s^-2]

            @property
            def ASTEROID212MEDEA_GM(self):
                r"""
                GM of asteroid 212 Medea (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID212MEDEA_GM


            __ASTEROID217EUDORA_GM: float  = 7.1074477882268534e+07 # [m^3 s^-2]

            @property
            def ASTEROID217EUDORA_GM(self):
                r"""
                GM of asteroid 217 Eudora (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID217EUDORA_GM


            __ASTEROID21LUTETIA_GM: float  = 1.1503721169099261e+08 # [m^3 s^-2]

            @property
            def ASTEROID21LUTETIA_GM(self):
                r"""
                GM of asteroid 21 Lutetia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID21LUTETIA_GM


            __ASTEROID22KALLIOPE_GM: float  = 1.1112751634270690e+09 # [m^3 s^-2]

            @property
            def ASTEROID22KALLIOPE_GM(self):
                r"""
                GM of asteroid 22 Kalliope (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID22KALLIOPE_GM


            __ASTEROID234BARBARA_GM: float  = 2.0492014129900661e+07 # [m^3 s^-2]

            @property
            def ASTEROID234BARBARA_GM(self):
                r"""
                GM of asteroid 234 Barbara (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID234BARBARA_GM


            __ASTEROID23THALIA_GM: float  = 2.0499213757716355e+08 # [m^3 s^-2]

            @property
            def ASTEROID23THALIA_GM(self):
                r"""
                GM of asteroid 23 Thalia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID23THALIA_GM


            __ASTEROID247EUKRATE_GM: float  = 5.9420762173955001e+08 # [m^3 s^-2]

            @property
            def ASTEROID247EUKRATE_GM(self):
                r"""
                GM of asteroid 247 Eukrate (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID247EUKRATE_GM


            __ASTEROID250BETTINA_GM: float  = 1.2408181646515161e+08 # [m^3 s^-2]

            @property
            def ASTEROID250BETTINA_GM(self):
                r"""
                GM of asteroid 250 Bettina (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID250BETTINA_GM


            __ASTEROID253MATHILDE_GM: float  = 8.4021633178130707e+07 # [m^3 s^-2]

            @property
            def ASTEROID253MATHILDE_GM(self):
                r"""
                GM of asteroid 253 Mathilde (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID253MATHILDE_GM


            __ASTEROID25PHOCAEA_GM: float  = 1.0366196678124746e+08 # [m^3 s^-2]

            @property
            def ASTEROID25PHOCAEA_GM(self):
                r"""
                GM of asteroid 25 Phocaea (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID25PHOCAEA_GM


            __ASTEROID266ALINE_GM: float  = 3.1756085199183367e+08 # [m^3 s^-2]

            @property
            def ASTEROID266ALINE_GM(self):
                r"""
                GM of asteroid 266 Aline (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID266ALINE_GM


            __ASTEROID26PROSERPINA_GM: float  = 2.0834242686152029e+08 # [m^3 s^-2]

            @property
            def ASTEROID26PROSERPINA_GM(self):
                r"""
                GM of asteroid 26 Proserpina (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID26PROSERPINA_GM


            __ASTEROID29AMPHITRITE_GM: float  = 9.5914046612273862e+08 # [m^3 s^-2]

            @property
            def ASTEROID29AMPHITRITE_GM(self):
                r"""
                GM of asteroid 29 Amphitrite (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID29AMPHITRITE_GM


            __ASTEROID2PALLAS_GM: float  = 1.3623519829068894e+10 # [m^3 s^-2]

            @property
            def ASTEROID2PALLAS_GM(self):
                r"""
                GM of asteroid 2 Pallas (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID2PALLAS_GM


            __ASTEROID304OLGA_GM: float  = 7.6417754608396792e+07 # [m^3 s^-2]

            @property
            def ASTEROID304OLGA_GM(self):
                r"""
                GM of asteroid 304 Olga (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID304OLGA_GM


            __ASTEROID308POLYXO_GM: float  = 5.1288139004370254e+08 # [m^3 s^-2]

            @property
            def ASTEROID308POLYXO_GM(self):
                r"""
                GM of asteroid 308 Polyxo (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID308POLYXO_GM


            __ASTEROID313CHALDAEA_GM: float  = 1.3568618932908778e+08 # [m^3 s^-2]

            @property
            def ASTEROID313CHALDAEA_GM(self):
                r"""
                GM of asteroid 313 Chaldaea (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID313CHALDAEA_GM


            __ASTEROID31EUPHROSYNE_GM: float  = 1.7562517866270849e+09 # [m^3 s^-2]

            @property
            def ASTEROID31EUPHROSYNE_GM(self):
                r"""
                GM of asteroid 31 Euphrosyne (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID31EUPHROSYNE_GM


            __ASTEROID322PHAEO_GM: float  = 8.6933740062976113e+07 # [m^3 s^-2]

            @property
            def ASTEROID322PHAEO_GM(self):
                r"""
                GM of asteroid 322 Phaeo (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID322PHAEO_GM


            __ASTEROID324BAMBERGA_GM: float  = 6.3292657018921930e+08 # [m^3 s^-2]

            @property
            def ASTEROID324BAMBERGA_GM(self):
                r"""
                GM of asteroid 324 Bamberga (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID324BAMBERGA_GM


            __ASTEROID335ROBERTA_GM: float  = 1.7274312784818528e+08 # [m^3 s^-2]

            @property
            def ASTEROID335ROBERTA_GM(self):
                r"""
                GM of asteroid 335 Roberta (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID335ROBERTA_GM


            __ASTEROID336LACADIERA_GM: float  = 8.1456992030180356e+07 # [m^3 s^-2]

            @property
            def ASTEROID336LACADIERA_GM(self):
                r"""
                GM of asteroid 336 Lacadiera (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID336LACADIERA_GM


            __ASTEROID33POLYHYMNIA_GM: float  = 2.8960743674231378e+08 # [m^3 s^-2]

            @property
            def ASTEROID33POLYHYMNIA_GM(self):
                r"""
                GM of asteroid 33 Polyhymnia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID33POLYHYMNIA_GM


            __ASTEROID346HERMENTARIA_GM: float  = 2.9556035531795935e+08 # [m^3 s^-2]

            @property
            def ASTEROID346HERMENTARIA_GM(self):
                r"""
                GM of asteroid 346 Hermentaria (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID346HERMENTARIA_GM


            __ASTEROID34CIRCE_GM: float  = 1.9273254973801046e+08 # [m^3 s^-2]

            @property
            def ASTEROID34CIRCE_GM(self):
                r"""
                GM of asteroid 34 Circe (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID34CIRCE_GM


            __ASTEROID350ORNAMENTA_GM: float  = 4.0527274829271427e+08 # [m^3 s^-2]

            @property
            def ASTEROID350ORNAMENTA_GM(self):
                r"""
                GM of asteroid 350 Ornamenta (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID350ORNAMENTA_GM


            __ASTEROID356LIGURIA_GM: float  = 5.5379028002660318e+08 # [m^3 s^-2]

            @property
            def ASTEROID356LIGURIA_GM(self):
                r"""
                GM of asteroid 356 Liguria (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID356LIGURIA_GM


            __ASTEROID36ATALANTE_GM: float  = 2.8796815326611882e+08 # [m^3 s^-2]

            @property
            def ASTEROID36ATALANTE_GM(self):
                r"""
                GM of asteroid 36 Atalante (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID36ATALANTE_GM


            __ASTEROID37FIDES_GM: float  = 3.1097050608961499e+08 # [m^3 s^-2]

            @property
            def ASTEROID37FIDES_GM(self):
                r"""
                GM of asteroid 37 Fides (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID37FIDES_GM


            __ASTEROID381MYRRHA_GM: float  = 4.2872473963951988e+08 # [m^3 s^-2]

            @property
            def ASTEROID381MYRRHA_GM(self):
                r"""
                GM of asteroid 381 Myrrha (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID381MYRRHA_GM


            __ASTEROID387AQUITANIA_GM: float  = 2.4837629612901240e+08 # [m^3 s^-2]

            @property
            def ASTEROID387AQUITANIA_GM(self):
                r"""
                GM of asteroid 387 Aquitania (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID387AQUITANIA_GM


            __ASTEROID388CHARYBDIS_GM: float  = 3.2080410583920504e+08 # [m^3 s^-2]

            @property
            def ASTEROID388CHARYBDIS_GM(self):
                r"""
                GM of asteroid 388 Charybdis (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID388CHARYBDIS_GM


            __ASTEROID3JUNO_GM: float  = 1.5650376033611044e+09 # [m^3 s^-2]

            @property
            def ASTEROID3JUNO_GM(self):
                r"""
                GM of asteroid 3 Juno (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID3JUNO_GM


            __ASTEROID404ARSINOE_GM: float  = 8.3335002699994517e+07 # [m^3 s^-2]

            @property
            def ASTEROID404ARSINOE_GM(self):
                r"""
                GM of asteroid 404 Arsinoe (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID404ARSINOE_GM


            __ASTEROID405THIA_GM: float  = 1.8291269394614959e+08 # [m^3 s^-2]

            @property
            def ASTEROID405THIA_GM(self):
                r"""
                GM of asteroid 405 Thia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID405THIA_GM


            __ASTEROID409ASPASIA_GM: float  = 1.4485207257722234e+05 # [m^3 s^-2]

            @property
            def ASTEROID409ASPASIA_GM(self):
                r"""
                GM of asteroid 409 Aspasia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID409ASPASIA_GM


            __ASTEROID410CHLORIS_GM: float  = 4.6130313447746411e+08 # [m^3 s^-2]

            @property
            def ASTEROID410CHLORIS_GM(self):
                r"""
                GM of asteroid 410 Chloris (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID410CHLORIS_GM


            __ASTEROID419AURELIA_GM: float  = 2.1890666809971182e+08 # [m^3 s^-2]

            @property
            def ASTEROID419AURELIA_GM(self):
                r"""
                GM of asteroid 419 Aurelia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID419AURELIA_GM


            __ASTEROID420BERTHOLDA_GM: float  = 6.8901048902653521e+08 # [m^3 s^-2]

            @property
            def ASTEROID420BERTHOLDA_GM(self):
                r"""
                GM of asteroid 420 Bertholda (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID420BERTHOLDA_GM


            __ASTEROID442EICHSFELDIA_GM: float  = 7.2696133676778026e+07 # [m^3 s^-2]

            @property
            def ASTEROID442EICHSFELDIA_GM(self):
                r"""
                GM of asteroid 442 Eichsfeldia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID442EICHSFELDIA_GM


            __ASTEROID444GYPTIS_GM: float  = 7.0727216281453015e+08 # [m^3 s^-2]

            @property
            def ASTEROID444GYPTIS_GM(self):
                r"""
                GM of asteroid 444 Gyptis (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID444GYPTIS_GM


            __ASTEROID445EDNA_GM: float  = 1.6203293575514477e+08 # [m^3 s^-2]

            @property
            def ASTEROID445EDNA_GM(self):
                r"""
                GM of asteroid 445 Edna (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID445EDNA_GM


            __ASTEROID44NYSA_GM: float  = 8.6199506289667382e+07 # [m^3 s^-2]

            @property
            def ASTEROID44NYSA_GM(self):
                r"""
                GM of asteroid 44 Nysa (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID44NYSA_GM


            __ASTEROID451PATIENTIA_GM: float  = 1.9885686160983943e+09 # [m^3 s^-2]

            @property
            def ASTEROID451PATIENTIA_GM(self):
                r"""
                GM of asteroid 451 Patientia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID451PATIENTIA_GM


            __ASTEROID455BRUCHSALIA_GM: float  = 1.4702122822900930e+08 # [m^3 s^-2]

            @property
            def ASTEROID455BRUCHSALIA_GM(self):
                r"""
                GM of asteroid 455 Bruchsalia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID455BRUCHSALIA_GM


            __ASTEROID469ARGENTINA_GM: float  = 5.2207383300069546e+05 # [m^3 s^-2]

            @property
            def ASTEROID469ARGENTINA_GM(self):
                r"""
                GM of asteroid 469 Argentina (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID469ARGENTINA_GM


            __ASTEROID46HESTIA_GM: float  = 4.6782985130525508e+08 # [m^3 s^-2]

            @property
            def ASTEROID46HESTIA_GM(self):
                r"""
                GM of asteroid 46 Hestia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID46HESTIA_GM


            __ASTEROID481EMITA_GM: float  = 1.1562903987899217e+08 # [m^3 s^-2]

            @property
            def ASTEROID481EMITA_GM(self):
                r"""
                GM of asteroid 481 Emita (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID481EMITA_GM


            __ASTEROID488KREUSA_GM: float  = 6.8444731733902973e+08 # [m^3 s^-2]

            @property
            def ASTEROID488KREUSA_GM(self):
                r"""
                GM of asteroid 488 Kreusa (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID488KREUSA_GM


            __ASTEROID4VESTA_GM: float  = 1.7288981939943646e+10 # [m^3 s^-2]

            @property
            def ASTEROID4VESTA_GM(self):
                r"""
                GM of asteroid 4 Vesta (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID4VESTA_GM


            __ASTEROID503EVELYN_GM: float  = 1.3326002607937277e+08 # [m^3 s^-2]

            @property
            def ASTEROID503EVELYN_GM(self):
                r"""
                GM of asteroid 503 Evelyn (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID503EVELYN_GM


            __ASTEROID505CAVA_GM: float  = 2.6643963557510264e+08 # [m^3 s^-2]

            @property
            def ASTEROID505CAVA_GM(self):
                r"""
                GM of asteroid 505 Cava (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID505CAVA_GM


            __ASTEROID511DAVIDA_GM: float  = 2.4220172669698846e+09 # [m^3 s^-2]

            @property
            def ASTEROID511DAVIDA_GM(self):
                r"""
                GM of asteroid 511 Davida (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID511DAVIDA_GM


            __ASTEROID516AMHERSTIA_GM: float  = 4.6419532383144037e+07 # [m^3 s^-2]

            @property
            def ASTEROID516AMHERSTIA_GM(self):
                r"""
                GM of asteroid 516 Amherstia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID516AMHERSTIA_GM


            __ASTEROID51NEMAUSA_GM: float  = 1.1850701657405471e+06 # [m^3 s^-2]

            @property
            def ASTEROID51NEMAUSA_GM(self):
                r"""
                GM of asteroid 51 Nemausa (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID51NEMAUSA_GM


            __ASTEROID52EUROPA_GM: float  = 1.4176132635840474e+09 # [m^3 s^-2]

            @property
            def ASTEROID52EUROPA_GM(self):
                r"""
                GM of asteroid 52 Europa (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID52EUROPA_GM


            __ASTEROID532HERCULINA_GM: float  = 1.5330914367975053e+09 # [m^3 s^-2]

            @property
            def ASTEROID532HERCULINA_GM(self):
                r"""
                GM of asteroid 532 Herculina (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID532HERCULINA_GM


            __ASTEROID53KALYPSO_GM: float  = 2.6631974354533394e+08 # [m^3 s^-2]

            @property
            def ASTEROID53KALYPSO_GM(self):
                r"""
                GM of asteroid 53 Kalypso (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID53KALYPSO_GM


            __ASTEROID54ALEXANDRA_GM: float  = 1.1137607808140530e+09 # [m^3 s^-2]

            @property
            def ASTEROID54ALEXANDRA_GM(self):
                r"""
                GM of asteroid 54 Alexandra (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID54ALEXANDRA_GM


            __ASTEROID55PANDORA_GM: float  = 7.2565502740078594e+07 # [m^3 s^-2]

            @property
            def ASTEROID55PANDORA_GM(self):
                r"""
                GM of asteroid 55 Pandora (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID55PANDORA_GM


            __ASTEROID568CHERUSKIA_GM: float  = 1.6092033080277349e+08 # [m^3 s^-2]

            @property
            def ASTEROID568CHERUSKIA_GM(self):
                r"""
                GM of asteroid 568 Cheruskia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID568CHERUSKIA_GM


            __ASTEROID569MISA_GM: float  = 9.4974536474941922e+07 # [m^3 s^-2]

            @property
            def ASTEROID569MISA_GM(self):
                r"""
                GM of asteroid 569 Misa (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID569MISA_GM


            __ASTEROID56MELETE_GM: float  = 3.5510132170246012e+08 # [m^3 s^-2]

            @property
            def ASTEROID56MELETE_GM(self):
                r"""
                GM of asteroid 56 Melete (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID56MELETE_GM


            __ASTEROID583KLOTILDE_GM: float  = 1.3306434325297216e+08 # [m^3 s^-2]

            @property
            def ASTEROID583KLOTILDE_GM(self):
                r"""
                GM of asteroid 583 Klotilde (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID583KLOTILDE_GM


            __ASTEROID584SEMIRAMIS_GM: float  = 3.8506494714933353e+07 # [m^3 s^-2]

            @property
            def ASTEROID584SEMIRAMIS_GM(self):
                r"""
                GM of asteroid 584 Semiramis (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID584SEMIRAMIS_GM


            __ASTEROID591IRMGARD_GM: float  = 1.5125998853964005e+06 # [m^3 s^-2]

            @property
            def ASTEROID591IRMGARD_GM(self):
                r"""
                GM of asteroid 591 Irmgard (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID591IRMGARD_GM


            __ASTEROID593TITANIA_GM: float  = 1.0449214411384316e+08 # [m^3 s^-2]

            @property
            def ASTEROID593TITANIA_GM(self):
                r"""
                GM of asteroid 593 Titania (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID593TITANIA_GM


            __ASTEROID595POLYXENA_GM: float  = 3.1721169220784370e+08 # [m^3 s^-2]

            @property
            def ASTEROID595POLYXENA_GM(self):
                r"""
                GM of asteroid 595 Polyxena (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID595POLYXENA_GM


            __ASTEROID599LUISA_GM: float  = 6.6785968807632834e+07 # [m^3 s^-2]

            @property
            def ASTEROID599LUISA_GM(self):
                r"""
                GM of asteroid 599 Luisa (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID599LUISA_GM


            __ASTEROID602MARIANNA_GM: float  = 2.5659694460527015e+08 # [m^3 s^-2]

            @property
            def ASTEROID602MARIANNA_GM(self):
                r"""
                GM of asteroid 602 Marianna (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID602MARIANNA_GM


            __ASTEROID618ELFRIEDE_GM: float  = 5.1050504107617317e+05 # [m^3 s^-2]

            @property
            def ASTEROID618ELFRIEDE_GM(self):
                r"""
                GM of asteroid 618 Elfriede (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID618ELFRIEDE_GM


            __ASTEROID61DANAE_GM: float  = 4.1297923333635475e+07 # [m^3 s^-2]

            @property
            def ASTEROID61DANAE_GM(self):
                r"""
                GM of asteroid 61 Danae (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID61DANAE_GM


            __ASTEROID626NOTBURGA_GM: float  = 2.5001067131733260e+08 # [m^3 s^-2]

            @property
            def ASTEROID626NOTBURGA_GM(self):
                r"""
                GM of asteroid 626 Notburga (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID626NOTBURGA_GM


            __ASTEROID62ERATO_GM: float  = 2.1232338569167166e+08 # [m^3 s^-2]

            @property
            def ASTEROID62ERATO_GM(self):
                r"""
                GM of asteroid 62 Erato (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID62ERATO_GM


            __ASTEROID63AUSONIA_GM: float  = 4.5857887479153414e+05 # [m^3 s^-2]

            @property
            def ASTEROID63AUSONIA_GM(self):
                r"""
                GM of asteroid 63 Ausonia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID63AUSONIA_GM


            __ASTEROID65CYBELE_GM: float  = 5.5875034960511513e+08 # [m^3 s^-2]

            @property
            def ASTEROID65CYBELE_GM(self):
                r"""
                GM of asteroid 65 Cybele (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID65CYBELE_GM


            __ASTEROID667DENISE_GM: float  = 7.8762291169638808e+07 # [m^3 s^-2]

            @property
            def ASTEROID667DENISE_GM(self):
                r"""
                GM of asteroid 667 Denise (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID667DENISE_GM


            __ASTEROID675LUDMILLA_GM: float  = 7.0108090353291896e+08 # [m^3 s^-2]

            @property
            def ASTEROID675LUDMILLA_GM(self):
                r"""
                GM of asteroid 675 Ludmilla (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID675LUDMILLA_GM


            __ASTEROID67ASIA_GM: float  = 4.7960213071519115e+07 # [m^3 s^-2]

            @property
            def ASTEROID67ASIA_GM(self):
                r"""
                GM of asteroid 67 Asia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID67ASIA_GM


            __ASTEROID690WRATISLAVIA_GM: float  = 5.9712951363080574e+08 # [m^3 s^-2]

            @property
            def ASTEROID690WRATISLAVIA_GM(self):
                r"""
                GM of asteroid 690 Wratislavia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID690WRATISLAVIA_GM


            __ASTEROID694EKARD_GM: float  = 2.8164391089152126e+03 # [m^3 s^-2]

            @property
            def ASTEROID694EKARD_GM(self):
                r"""
                GM of asteroid 694 Ekard (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID694EKARD_GM


            __ASTEROID6HEBE_GM: float  = 9.4013859632630443e+08 # [m^3 s^-2]

            @property
            def ASTEROID6HEBE_GM(self):
                r"""
                GM of asteroid 6 Hebe (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID6HEBE_GM


            __ASTEROID704INTERAMNIA_GM: float  = 2.5503367929164454e+09 # [m^3 s^-2]

            @property
            def ASTEROID704INTERAMNIA_GM(self):
                r"""
                GM of asteroid 704 Interamnia (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID704INTERAMNIA_GM


            __ASTEROID70PANOPAEA_GM: float  = 2.7708718966786989e+08 # [m^3 s^-2]

            @property
            def ASTEROID70PANOPAEA_GM(self):
                r"""
                GM of asteroid 70 Panopaea (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID70PANOPAEA_GM


            __ASTEROID718ERIDA_GM: float  = 4.0557482554503977e+07 # [m^3 s^-2]

            @property
            def ASTEROID718ERIDA_GM(self):
                r"""
                GM of asteroid 718 Erida (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID718ERIDA_GM


            __ASTEROID735MARGHANNA_GM: float  = 1.0038522798255246e+08 # [m^3 s^-2]

            @property
            def ASTEROID735MARGHANNA_GM(self):
                r"""
                GM of asteroid 735 Marghanna (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID735MARGHANNA_GM


            __ASTEROID739MANDEVILLE_GM: float  = 3.0413247983619986e+08 # [m^3 s^-2]

            @property
            def ASTEROID739MANDEVILLE_GM(self):
                r"""
                GM of asteroid 739 Mandeville (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID739MANDEVILLE_GM


            __ASTEROID747WINCHESTER_GM: float  = 9.6009858605926026e+07 # [m^3 s^-2]

            @property
            def ASTEROID747WINCHESTER_GM(self):
                r"""
                GM of asteroid 747 Winchester (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID747WINCHESTER_GM


            __ASTEROID751FAINA_GM: float  = 3.2994347243673589e+08 # [m^3 s^-2]

            @property
            def ASTEROID751FAINA_GM(self):
                r"""
                GM of asteroid 751 Faina (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID751FAINA_GM


            __ASTEROID752SULAMITIS_GM: float  = 6.0508535365628257e+07 # [m^3 s^-2]

            @property
            def ASTEROID752SULAMITIS_GM(self):
                r"""
                GM of asteroid 752 Sulamitis (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID752SULAMITIS_GM


            __ASTEROID75EURYDIKE_GM: float  = 4.2715759272180030e+07 # [m^3 s^-2]

            @property
            def ASTEROID75EURYDIKE_GM(self):
                r"""
                GM of asteroid 75 Eurydike (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID75EURYDIKE_GM


            __ASTEROID78DIANA_GM: float  = 3.4003238496349289e+08 # [m^3 s^-2]

            @property
            def ASTEROID78DIANA_GM(self):
                r"""
                GM of asteroid 78 Diana (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID78DIANA_GM


            __ASTEROID790PRETORIA_GM: float  = 8.6438330145248586e+08 # [m^3 s^-2]

            @property
            def ASTEROID790PRETORIA_GM(self):
                r"""
                GM of asteroid 790 Pretoria (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID790PRETORIA_GM


            __ASTEROID791ANI_GM: float  = 2.7128482223705070e+08 # [m^3 s^-2]

            @property
            def ASTEROID791ANI_GM(self):
                r"""
                GM of asteroid 791 Ani (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID791ANI_GM


            __ASTEROID7IRIS_GM: float  = 8.3629208510440051e+08 # [m^3 s^-2]

            @property
            def ASTEROID7IRIS_GM(self):
                r"""
                GM of asteroid 7 Iris (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID7IRIS_GM


            __ASTEROID804HISPANIA_GM: float  = 1.3669667630778693e+08 # [m^3 s^-2]

            @property
            def ASTEROID804HISPANIA_GM(self):
                r"""
                GM of asteroid 804 Hispania (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID804HISPANIA_GM


            __ASTEROID80SAPPHO_GM: float  = 1.1784221437378471e+08 # [m^3 s^-2]

            @property
            def ASTEROID80SAPPHO_GM(self):
                r"""
                GM of asteroid 80 Sappho (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID80SAPPHO_GM


            __ASTEROID814TAURIS_GM: float  = 3.2159462197573089e+08 # [m^3 s^-2]

            @property
            def ASTEROID814TAURIS_GM(self):
                r"""
                GM of asteroid 814 Tauris (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID814TAURIS_GM


            __ASTEROID85IO_GM: float  = 2.9069348822079933e+08 # [m^3 s^-2]

            @property
            def ASTEROID85IO_GM(self):
                r"""
                GM of asteroid 85 Io (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID85IO_GM


            __ASTEROID88THISBE_GM: float  = 9.4067567843322719e+08 # [m^3 s^-2]

            @property
            def ASTEROID88THISBE_GM(self):
                r"""
                GM of asteroid 88 Thisbe (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID88THISBE_GM


            __ASTEROID8FLORA_GM: float  = 4.4557515312994348e+08 # [m^3 s^-2]

            @property
            def ASTEROID8FLORA_GM(self):
                r"""
                GM of asteroid 8 Flora (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID8FLORA_GM


            __ASTEROID914PALISANA_GM: float  = 5.7331462009756002e+07 # [m^3 s^-2]

            @property
            def ASTEROID914PALISANA_GM(self):
                r"""
                GM of asteroid 914 Palisana (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID914PALISANA_GM


            __ASTEROID93MINERVA_GM: float  = 5.0039749795125115e+08 # [m^3 s^-2]

            @property
            def ASTEROID93MINERVA_GM(self):
                r"""
                GM of asteroid 93 Minerva (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID93MINERVA_GM


            __ASTEROID949HEL_GM: float  = 8.0964451809677710e+07 # [m^3 s^-2]

            @property
            def ASTEROID949HEL_GM(self):
                r"""
                GM of asteroid 949 Hel (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID949HEL_GM


            __ASTEROID94AURORA_GM: float  = 1.9949958990308900e+09 # [m^3 s^-2]

            @property
            def ASTEROID94AURORA_GM(self):
                r"""
                GM of asteroid 94 Aurora (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID94AURORA_GM


            __ASTEROID97KLOTHO_GM: float  = 6.1873651052720885e+06 # [m^3 s^-2]

            @property
            def ASTEROID97KLOTHO_GM(self):
                r"""
                GM of asteroid 97 Klotho (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID97KLOTHO_GM


            __ASTEROID9METIS_GM: float  = 5.5771967924306053e+08 # [m^3 s^-2]

            @property
            def ASTEROID9METIS_GM(self):
                r"""
                GM of asteroid 9 Metis (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROID9METIS_GM


            __ASTEROIDRING_GM: float  = 8.9697172580075571e+09 # [m^3 s^-2]

            @property
            def ASTEROIDRING_GM(self):
                r"""
                GM of the asteroid ring (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__ASTEROIDRING_GM


            __ASTEROIDRING_ORBITALSEMIMAJORAXIS: float  = 3.1477080248116020e+00 # [au]

            @property
            def ASTEROIDRING_ORBITALSEMIMAJORAXIS(self):
                r"""
                Barycentric distance (orbital semi-major axis) of the asteroid ring (INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : true
                #Scalar: true
                #Unit: [au]
                """

                return self.__ASTEROIDRING_ORBITALSEMIMAJORAXIS


            __ASTRONOMICALUNIT_METER: float  = 1.4959787070000000e+11 # [m]

            @property
            def ASTRONOMICALUNIT_METER(self):
                r"""
                Astronomical unit (au) length in m

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : true
                #Scalar: true
                #Unit: [m]
                """

                return self.__ASTRONOMICALUNIT_METER


            __EARTHSYSTEM_GM: float  = 4.0350325101102718e+14 # [m^3 s^-2]

            @property
            def EARTHSYSTEM_GM(self):
                r"""
                GM of the Earth-system (TCB-compatible value). The gravitational constant includes the contribution of its satellite, the Moon

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__EARTHSYSTEM_GM


            __EARTHTOMOON_MASSRATIO: float  = 8.1300570000000000e+01

            @property
            def EARTHTOMOON_MASSRATIO(self):
                r"""
                Ratio of Earth to Moon mass (INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__EARTHTOMOON_MASSRATIO


            __EARTH_EQUATORIALRADIUS: float  = 6.3781366988942710e+06 # [m]

            @property
            def EARTH_EQUATORIALRADIUS(self):
                r"""
                Mean equatorial radius of the Earth (INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : true
                #Scalar: true
                #Unit: [m]
                """

                return self.__EARTH_EQUATORIALRADIUS


            __EARTH_GM: float  = 3.9860045081157502e+14 # [m^3 s^-2]

            @property
            def EARTH_GM(self):
                r"""
                Geocentric gravitational constant (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__EARTH_GM


            __JUPITERSYSTEM_GM: float  = 1.2671276453465735e+17 # [m^3 s^-2]

            @property
            def JUPITERSYSTEM_GM(self):
                r"""
                GM of the Jupiter-system (TCB-compatible value). The gravitational constant includes the contribution of its satellites

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__JUPITERSYSTEM_GM


            __JUPITER_GM: float  = 1.2671276453465734e+17 # [m^3 s^-2]

            @property
            def JUPITER_GM(self):
                r"""
                GM of Jupiter (TCB-compatible value; this is a low-accuracy parameter used in the relativistic model, ignoring the contribution of planetary satellites)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__JUPITER_GM


            __MARSSYSTEM_GM: float  = 4.2828375886337909e+13 # [m^3 s^-2]

            @property
            def MARSSYSTEM_GM(self):
                r"""
                GM of the Mars-system (TCB-compatible value). The gravitational constant includes the contribution of its satellites

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__MARSSYSTEM_GM


            __MARS_GM: float  = 4.2828375886337906e+13 # [m^3 s^-2]

            @property
            def MARS_GM(self):
                r"""
                GM of Mars (TCB-compatible value; this is a low-accuracy parameter used in the relativistic model, ignoring the contribution of planetary satellites)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__MARS_GM


            __MERCURYSYSTEM_GM: float  = 2.2032080834196276e+13 # [m^3 s^-2]

            @property
            def MERCURYSYSTEM_GM(self):
                r"""
                GM of the Mercury(-system) (TCB-compatible value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__MERCURYSYSTEM_GM


            __MERCURY_GM: float  = 2.2032080834196277e+13 # [m^3 s^-2]

            @property
            def MERCURY_GM(self):
                r"""
                GM of Mercury (TCB-compatible value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__MERCURY_GM


            __MOONTOEARTH_MASSRATIO: float  = 0.0123000368

            @property
            def MOONTOEARTH_MASSRATIO(self):
                r"""
                Ratio of Moon to Earth mass (INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: []
                """

                return self.__MOONTOEARTH_MASSRATIO


            __MOON_EQUATORIALRADIUS: float  = 1.7380000269480340e+06 # [m]

            @property
            def MOON_EQUATORIALRADIUS(self):
                r"""
                Mean equatorial radius of the Moon (INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : true
                #Scalar: true
                #Unit: [m]
                """

                return self.__MOON_EQUATORIALRADIUS


            __MOON_GM: float  = 4.9028001994521693e+12 # [m^3 s^-2]

            @property
            def MOON_GM(self):
                r"""
                Selenocentric gravitational constant (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__MOON_GM


            __NEPTUNESYSTEM_GM: float  = 6.8365271283644811e+15 # [m^3 s^-2]

            @property
            def NEPTUNESYSTEM_GM(self):
                r"""
                GM of the Neptune-system (TCB-compatible value). The gravitational constant includes the contribution of its satellites

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__NEPTUNESYSTEM_GM


            __NEPTUNE_GM: float  = 6.8365271283644810e+15 # [m^3 s^-2]

            @property
            def NEPTUNE_GM(self):
                r"""
                GM of Neptune (TCB-compatible value; this is a low-accuracy parameter used in the relativistic model, ignoring the contribution of planetary satellites)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__NEPTUNE_GM


            __PPN_BETA: float  = 1.0

            @property
            def PPN_BETA(self):
                r"""
                General relativistic standard PPN parameter \beta (INPOP10e value, fixed to 1)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__PPN_BETA


            __PPN_GAMMA: float  = 1.0

            @property
            def PPN_GAMMA(self):
                r"""
                General relativistic standard PPN parameter \gamma (INPOP10e value, fixed to 1)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__PPN_GAMMA


            __PLUTOSYSTEM_GM: float  = 9.7178245029026624e+11 # [m^3 s^-2]

            @property
            def PLUTOSYSTEM_GM(self):
                r"""
                GM of the Pluto-system (TCB-compatible value). The gravitational constant includes the contribution of its satellite, Charon

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__PLUTOSYSTEM_GM


            __PLUTO_GM: float  = 9.7178245029026624e+11 # [m^3 s^-2]

            @property
            def PLUTO_GM(self):
                r"""
                GM of Pluto (TCB-compatible value; this is a low-accuracy parameter used in the relativistic model, ignoring the contribution of Charon)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__PLUTO_GM


            __SATURNSYSTEM_GM: float  = 3.7940585442640140e+16 # [m^3 s^-2]

            @property
            def SATURNSYSTEM_GM(self):
                r"""
                GM of the Saturn-system (TCB-compatible value). The gravitational constant includes the contribution of its satellites

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__SATURNSYSTEM_GM


            __SATURN_GM: float  = 3.7940585442640144e+16 # [m^3 s^-2]

            @property
            def SATURN_GM(self):
                r"""
                GM of Saturn (TCB-compatible value; this is a low-accuracy parameter used in the relativistic model, ignoring the contribution of planetary satellites)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__SATURN_GM


            __SUNTOEARTHSYSTEM_MASSRATIO: float  = 328900.552289

            @property
            def SUNTOEARTHSYSTEM_MASSRATIO(self):
                r"""
                Ratio of Sun to Earth-system mass (INPOP10e value). The planetary mass includes the contribution of its satellite, the Moon

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: []
                """

                return self.__SUNTOEARTHSYSTEM_MASSRATIO


            __SUNTOJUPITERSYSTEM_MASSRATIO: float  = 1047.348644

            @property
            def SUNTOJUPITERSYSTEM_MASSRATIO(self):
                r"""
                Ratio of Sun to Jupiter-system mass (INPOP10e value). The planetary mass includes the contribution of its satellites

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: []
                """

                return self.__SUNTOJUPITERSYSTEM_MASSRATIO


            __SUNTOMARSSYSTEM_MASSRATIO: float  = 3098704.0

            @property
            def SUNTOMARSSYSTEM_MASSRATIO(self):
                r"""
                Ratio of Sun to Mars-system mass (INPOP10e value). The planetary mass includes the contribution of its satellites

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: []
                """

                return self.__SUNTOMARSSYSTEM_MASSRATIO


            __SUNTOMERCURYSYSTEM_MASSRATIO: float  = 6023600.0

            @property
            def SUNTOMERCURYSYSTEM_MASSRATIO(self):
                r"""
                Ratio of Sun to Mercury(-system) mass (INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: []
                """

                return self.__SUNTOMERCURYSYSTEM_MASSRATIO


            __SUNTONEPTUNESYSTEM_MASSRATIO: float  = 19412.26

            @property
            def SUNTONEPTUNESYSTEM_MASSRATIO(self):
                r"""
                Ratio of Sun to Neptune-system mass (INPOP10e value). The planetary mass includes the contribution of its satellites

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: []
                """

                return self.__SUNTONEPTUNESYSTEM_MASSRATIO


            __SUNTOPLUTOSYSTEM_MASSRATIO: float  = 136566000.0

            @property
            def SUNTOPLUTOSYSTEM_MASSRATIO(self):
                r"""
                Ratio of Sun to Pluto-system mass (INPOP10e value). The 'planetary' mass includes the contribution of its satellite, Charon

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: []
                """

                return self.__SUNTOPLUTOSYSTEM_MASSRATIO


            __SUNTOSATURNSYSTEM_MASSRATIO: float  = 3497.902

            @property
            def SUNTOSATURNSYSTEM_MASSRATIO(self):
                r"""
                Ratio of Sun to Saturn-system mass (INPOP10e value). The planetary mass includes the contribution of its satellites

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: []
                """

                return self.__SUNTOSATURNSYSTEM_MASSRATIO


            __SUNTOURANUSSYSTEM_MASSRATIO: float  = 22902.98

            @property
            def SUNTOURANUSSYSTEM_MASSRATIO(self):
                r"""
                Ratio of Sun to Uranus-system mass (INPOP10e value). The planetary mass includes the contribution of its satellites

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: []
                """

                return self.__SUNTOURANUSSYSTEM_MASSRATIO


            __SUNTOVENUSSYSTEM_MASSRATIO: float  = 408523.72

            @property
            def SUNTOVENUSSYSTEM_MASSRATIO(self):
                r"""
                Ratio of Sun to Venus(-system) mass (INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: []
                """

                return self.__SUNTOVENUSSYSTEM_MASSRATIO


            __SUN_EQUATORIALRADIUS: float  = 6.9600001079161780e+08 # [m]

            @property
            def SUN_EQUATORIALRADIUS(self):
                r"""
                Mean equatorial radius of the Sun (INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : true
                #Scalar: true
                #Unit: [m]
                """

                return self.__SUN_EQUATORIALRADIUS


            __SUN_GM: float  = 1.3271244210789467e+20 # [m^3 s^-2]

            @property
            def SUN_GM(self):
                r"""
                Heliocentric gravitational constant (TCB-compatible value; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__SUN_GM


            __SUN_JSUB2: float  = 1.8000000000000000e-07

            @property
            def SUN_JSUB2(self):
                r"""
                Dynamical form-factor of the Sun (Stokes' second-degree zonal harmonic of the solar potential; INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : true
                #Scalar: true
                #Unit: []
                """

                return self.__SUN_JSUB2


            __SUN_NORTHROTATIONALPOLE_DECLINATION: float  = 6.3870000000000000e+01 # [deg]

            @property
            def SUN_NORTHROTATIONALPOLE_DECLINATION(self):
                r"""
                Declination \delta_0 of the north pole of rotation of the Sun (INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : true
                #Scalar: true
                #Unit: [deg]
                """

                return self.__SUN_NORTHROTATIONALPOLE_DECLINATION


            __SUN_NORTHROTATIONALPOLE_RIGHTASCENSION: float  = 2.8613000000000000e+02 # [deg]

            @property
            def SUN_NORTHROTATIONALPOLE_RIGHTASCENSION(self):
                r"""
                Right ascension \alpha_0 of the north pole of rotation of the Sun (INPOP10e value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : true
                #Scalar: true
                #Unit: [deg]
                """

                return self.__SUN_NORTHROTATIONALPOLE_RIGHTASCENSION


            __URANUSSYSTEM_GM: float  = 5.7945490985393442e+15 # [m^3 s^-2]

            @property
            def URANUSSYSTEM_GM(self):
                r"""
                GM of the Uranus-system (TCB-compatible value). The gravitational constant includes the contribution of its satellites

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__URANUSSYSTEM_GM


            __URANUS_GM: float  = 5.7945490985393440e+15 # [m^3 s^-2]

            @property
            def URANUS_GM(self):
                r"""
                GM of Uranus (TCB-compatible value; this is a low-accuracy parameter used in the relativistic model, ignoring the contribution of planetary satellites)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__URANUS_GM


            __VENUSSYSTEM_GM: float  = 3.2485859679756975e+14 # [m^3 s^-2]

            @property
            def VENUSSYSTEM_GM(self):
                r"""
                GM of the Venus(-system) (TCB-compatible value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__VENUSSYSTEM_GM


            __VENUS_GM: float  = 3.2485859679756975e+14 # [m^3 s^-2]

            @property
            def VENUS_GM(self):
                r"""
                GM of Venus (TCB-compatible value)

                #Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0<br/>
                #Basic : false
                #Scalar: true
                #Unit: [m^3 s^-2]
                """

                return self.__VENUS_GM


        __INVERSEFINESTRUCTURE_CONSTANT: float  = 137.035999142

        @property
        def INVERSEFINESTRUCTURE_CONSTANT(self):
            r"""
            Inverse of the fine structure constant. Note: best-measured value equals 137.035999139 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))

            #Basic : false
            #Scalar: true
            #Unit: []
            """

            return self.__INVERSEFINESTRUCTURE_CONSTANT


        __JULIANCENTURY_JULIANYEAR: float  = 36525.0 # [day]

        @property
        def JULIANCENTURY_JULIANYEAR(self):
            r"""
            Number of days per Julian century

            #Source: IAU definition<br/>
            #Basic : true
            #Scalar: true
            #Unit: [day]
            """

            return self.__JULIANCENTURY_JULIANYEAR


        __JULIANDATE_J2000: float  = 2451545.0 # [JD]

        @property
        def JULIANDATE_J2000(self):
            r"""
            Julian date of the standard epoch J2000 = J2000.0, i.e., calendar date 2000 January 1, 12:00:00 h = 2000-01-01T12:00:00 TT

            #Source: Definition (e.g., ESA, 1997, 'The Hipparcos and Tycho Catalogues', Volume 1, page 27)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [JD]
            """

            return self.__JULIANDATE_J2000


        __JULIANDATE_J2010: float  = 2455197.5 # [JD]

        @property
        def JULIANDATE_J2010(self):
            r"""
            Julian date of the Gaia-specific reference epoch J2010 = J2010.0 = JD2455197.5 = 2010-01-01T00:00:00

            #Source: U. Bastian, 5 July 2007, 'Reference systems, conventions, and notations for Gaia', GAIA-CA-SP-ARI-BAS-003, issue 6, revision 1, Section 3.5<br/>
            #Basic : false
            #Scalar: true
            #Unit: [JD]
            """

            return self.__JULIANDATE_J2010


        __JULIANYEAR_DAY: float  = 365.25 # [day]

        @property
        def JULIANYEAR_DAY(self):
            r"""
            Number of days per Julian year

            #Source: IAU definition<br/>
            #Basic : false
            #Scalar: true
            #Unit: [day]
            """

            return self.__JULIANYEAR_DAY


        __JUPITERSYSTEM_ASTROMETRICSIGNATURE_10PARSEC: float  = 497.0 # [10^-6 arcsec]

        @property
        def JUPITERSYSTEM_ASTROMETRICSIGNATURE_10PARSEC(self):
            r"""
            Astrometric signature of the Sun induced by the Jupiter system for an observer located at a distance of 10 pc from the Sun

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.7, Equation 1.22, page 11<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__JUPITERSYSTEM_ASTROMETRICSIGNATURE_10PARSEC


        __JUPITERSYSTEM_MASS: float  = 1.89858e+27 # [kg]

        @property
        def JUPITERSYSTEM_MASS(self):
            r"""
            Jupiter-system mass (IAU 2009 CBE value). The planetary mass includes the contribution of its satellites

            #Basic : false
            #Scalar: true
            #Unit: [kg]
            """

            return self.__JUPITERSYSTEM_MASS


        __JUPITERSYSTEM_ORBITALECCENTRICITY_J2000: float  = 0.04838624

        @property
        def JUPITERSYSTEM_ORBITALECCENTRICITY_J2000(self):
            r"""
            Mean orbital eccentricity of Jupiter, at the standard epoch J2000.0. The mean orbital eccentricity is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__JUPITERSYSTEM_ORBITALECCENTRICITY_J2000


        __JUPITERSYSTEM_ORBITALINCLINATION_J2000: float  = 1.30439695 # [deg]

        @property
        def JUPITERSYSTEM_ORBITALINCLINATION_J2000(self):
            r"""
            Mean orbital inclination of Jupiter, at the standard epoch J2000.0. The mean orbital inclination is is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__JUPITERSYSTEM_ORBITALINCLINATION_J2000


        __JUPITERSYSTEM_ORBITALPERIOD: float  = 11.862615 # [yr]

        @property
        def JUPITERSYSTEM_ORBITALPERIOD(self):
            r"""
            Sidereal orbital period

            #Source: Values derived from the mean longitude rates in Table 5.8.1 in P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, page 316; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [yr]
            """

            return self.__JUPITERSYSTEM_ORBITALPERIOD


        __JUPITERSYSTEM_ORBITALSEMIMAJORAXIS_J2000: float  = 5.20288700 # [au]

        @property
        def JUPITERSYSTEM_ORBITALSEMIMAJORAXIS_J2000(self):
            r"""
            Mean orbital semi-major axis of Jupiter, at the standard epoch J2000.0. The mean orbital semi-major axis is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: [au]
            """

            return self.__JUPITERSYSTEM_ORBITALSEMIMAJORAXIS_J2000


        __JUPITERSYSTEM_RADIALVELOCITYSIGNATURE: float  = 12.5 # [m s^-1]

        @property
        def JUPITERSYSTEM_RADIALVELOCITYSIGNATURE(self):
            r"""
            Radial-velocity amplitude of the Sun induced by the Jupiter system for 'an edge-on observer' (i.e., an observer in the orbital plane of the Jupiter system)

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.4, Equation 1.18, page 9<br/>
            #Basic : false
            #Scalar: true
            #Unit: [m s^-1]
            """

            return self.__JUPITERSYSTEM_RADIALVELOCITYSIGNATURE


        __JUPITER_ENCOMPASSINGSPHERERADIUS: float  = 7.14917e+07 # [m]

        @property
        def JUPITER_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around Jupiter which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__JUPITER_ENCOMPASSINGSPHERERADIUS


        __JUPITER_EQUATORIALRADIUS: float  = 7.14917e+07 # [m]

        @property
        def JUPITER_EQUATORIALRADIUS(self):
            r"""
            Equatorial radius of Jupiter

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__JUPITER_EQUATORIALRADIUS


        __JUPITER_EQUATORIALRADIUS_NOMINAL: float  = 7.14920e+07 # [m]

        @property
        def JUPITER_EQUATORIALRADIUS_NOMINAL(self):
            r"""
            Nominal equatorial radius of Jupiter (one-bar value), in SI units. This nominal value shall be understood as conversion factor only

            #Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__JUPITER_EQUATORIALRADIUS_NOMINAL


        __JUPITER_FLATTENING: float  = 6.487440e-02

        @property
        def JUPITER_FLATTENING(self):
            r"""
            Geometrical flattening factor f of Jupiter (f = (a-b)/a)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__JUPITER_FLATTENING


        __JUPITER_FLUXREDUCTION_MAXIMUM: float  = 1.010 # [%]

        @property
        def JUPITER_FLUXREDUCTION_MAXIMUM(self):
            r"""
            Maximum reduction of the solar flux for an observer external to the solar system during a transit of Jupiter

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.2, Equation 2.4, page 14<br/>
            #Basic : false
            #Scalar: true
            #Unit: [%]
            """

            return self.__JUPITER_FLUXREDUCTION_MAXIMUM


        __JUPITER_GM_NOMINAL: float  = 1.26686530e+17 # [m^3 s^-2]

        @property
        def JUPITER_GM_NOMINAL(self):
            r"""
            Nominal GM of Jupiter, in SI units. This nominal value shall be understood as conversion factor only

            #Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__JUPITER_GM_NOMINAL


        __JUPITER_GEOMETRICALBEDO: float  = 0.52

        @property
        def JUPITER_GEOMETRICALBEDO(self):
            r"""
            Geometric albedo of Jupiter (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__JUPITER_GEOMETRICALBEDO


        __JUPITER_JSUB2: float  = 0.014697

        @property
        def JUPITER_JSUB2(self):
            r"""
            Dynamical form-factor of Jupiter (oblateness or Stokes' second-degree zonal harmonic of the potential)

            #Source: P.R. Weissman, L.-A. McFadden, T.V. Johnson (eds.), 1999, 'Encyclopedia of the Solar System (first edition)', Academic Press, page 342<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__JUPITER_JSUB2


        __JUPITER_LIGHTDEFLECTION_LIMB: float  = 16635.0 # [10^-6 arcsec]

        @property
        def JUPITER_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of Jupiter

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__JUPITER_LIGHTDEFLECTION_LIMB


        __JUPITER_MASS: float  = 1.898130e+27 # [kg]

        @property
        def JUPITER_MASS(self):
            r"""
            Mass of Jupiter (do not use for high-precision (orbit) calculations)

            #Source: R.A. Jacobson, 2005, 'Jovian Satellite ephemeris - JUP230', priv. comm.; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [kg]
            """

            return self.__JUPITER_MASS


        __JUPITER_MASSDENSITY_MEAN: float  = 1.3262 # [g cm^-3]

        @property
        def JUPITER_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of Jupiter

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__JUPITER_MASSDENSITY_MEAN


        __JUPITER_NORTHROTATIONALPOLE_DECLINATION: float  = 64.495303 # [deg]

        @property
        def JUPITER_NORTHROTATIONALPOLE_DECLINATION(self):
            r"""
            IAU-recommended value for the declination \delta_0 of the north pole of rotation of Jupiter. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__JUPITER_NORTHROTATIONALPOLE_DECLINATION


        __JUPITER_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE: float  = 0.000000066064 # [deg day^-1]

        @property
        def JUPITER_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the declination \delta_0 of the north pole of rotation of Jupiter. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : false
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__JUPITER_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE


        __JUPITER_NORTHROTATIONALPOLE_RIGHTASCENSION: float  = 268.056595 # [deg]

        @property
        def JUPITER_NORTHROTATIONALPOLE_RIGHTASCENSION(self):
            r"""
            IAU-recommended value for the right ascension \alpha_0 of the north pole of rotation of Jupiter. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__JUPITER_NORTHROTATIONALPOLE_RIGHTASCENSION


        __JUPITER_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE: float  = -0.000000177933 # [deg day^-1]

        @property
        def JUPITER_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the right ascension \alpha_0 of the north pole of rotation of Jupiter. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : false
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__JUPITER_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE


        __JUPITER_POLARRADIUS: float  = 6.68537e+07 # [m]

        @property
        def JUPITER_POLARRADIUS(self):
            r"""
            Polar radius of Jupiter

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__JUPITER_POLARRADIUS


        __JUPITER_POLARRADIUS_NOMINAL: float  = 6.68540e+07 # [m]

        @property
        def JUPITER_POLARRADIUS_NOMINAL(self):
            r"""
            Nominal polar radius of Jupiter (one-bar value), in SI units. This nominal value shall be understood as conversion factor only

            #Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__JUPITER_POLARRADIUS_NOMINAL


        __JUPITER_PRIMEMERIDIAN_EPHEMERISPOSITION: float  = 284.95 # [deg]

        @property
        def JUPITER_PRIMEMERIDIAN_EPHEMERISPOSITION(self):
            r"""
            IAU-recommended value for the ephemeris position of the prime meridian of Jupiter. The prime meridian refers to the rotation of the magnetic field System III. System I (W_{I} = 67.1 deg + 877.900 deg day^-1) refers to the mean atmospheric equatorial rotation. System II (W_{II} = 43.3 deg + 870.270 deg day^-1) refers to the mean atmospheric rotation north of the south component of the north equatorial belt, and south of the north component of the south equatorial belt. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__JUPITER_PRIMEMERIDIAN_EPHEMERISPOSITION


        __JUPITER_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE: float  = 870.5360000 # [deg day^-1]

        @property
        def JUPITER_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch of 1.5 January 2000 = JD2451545.0 TDB) of the ephemeris position of the prime meridian of Jupiter. The prime meridian refers to the rotation of the magnetic field System III. System I (W_{I} = 67.1 deg + 877.900 deg day^-1) refers to the mean atmospheric equatorial rotation. System II (W_{II} = 43.3 deg + 870.270 deg day^-1) refers to the mean atmospheric rotation north of the south component of the north equatorial belt, and south of the north component of the south equatorial belt. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__JUPITER_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE


        __JUPITER_TRANSITPROBABILITY: float  = 0.098 # [%]

        @property
        def JUPITER_TRANSITPROBABILITY(self):
            r"""
            Geometric transit probability (Jupiter transiting the Sun) for an observer external to the solar system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.1, Equation 2.2, page 14<br/>
            #Basic : false
            #Scalar: true
            #Unit: [%]
            """

            return self.__JUPITER_TRANSITPROBABILITY


        __JUPITER_TRANSITTIME_MAXIMUM: float  = 1.36 # [day]

        @property
        def JUPITER_TRANSITTIME_MAXIMUM(self):
            r"""
            Maximum transit time of Jupiter (transiting the Sun) for an observer external to the solar system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.3, Equation 2.5, page 15<br/>
            #Basic : false
            #Scalar: true
            #Unit: [day]
            """

            return self.__JUPITER_TRANSITTIME_MAXIMUM


        __JUPITER_VONEZEROMAGNITUDE: float  = -9.40 # [mag]

        @property
        def JUPITER_VONEZEROMAGNITUDE(self):
            r"""
            V(1,0) magnitude of Jupiter (i.e., the visual magnitude of the planet reduced to a distance of 1 au from both the Sun and Earth and phase angle zero). This parameter is also refered to as absolute magnitude in planetary sciences

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__JUPITER_VONEZEROMAGNITUDE


        __JUPITER_VOLUMETRICRADIUS: float  = 6.99110e+07 # [m]

        @property
        def JUPITER_VOLUMETRICRADIUS(self):
            r"""
            Mean volumetric radius of Jupiter

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__JUPITER_VOLUMETRICRADIUS


        __K1IIIMPSTAR_SPECTRUM_NUMBEROFPHOTONS: str  = "Nature/K1IIIMPStar_Spectrum_NumberOfPhotons_001.fits"

        @property
        def K1IIIMPSTAR_SPECTRUM_NUMBEROFPHOTONS(self):
            r"""
            Photon flux density N_{\lambda}(\lambda) of an unreddened K1III metal-poor (MP) star (Pickles' star number 082) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)

            #Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__K1IIIMPSTAR_SPECTRUM_NUMBEROFPHOTONS


        __K1IIIMPSTAR_SPECTRUM_NUMBEROFPHOTONSHIGHRESOLUTION: str  = "Nature/K1IIIMPStar_Spectrum_NumberOfPhotonsHighResolution_001.fits"

        @property
        def K1IIIMPSTAR_SPECTRUM_NUMBEROFPHOTONSHIGHRESOLUTION(self):
            r"""
            High-resolution photon-flux density N_{\lambda}(\lambda) of an unreddened K1III metal-poor (MP) star at V = 15 mag. The data refer to a high-resolution Kurucz-model spectrum with the following properties: effective temperature T_eff = 4500 K, logarithm of surface gravity log g = 2.0, metallicity [Fe/H] = -1.5, alpha-elements [\alpha/Fe] = +0.4, rotational velocity v sini = 5 km s^-1, micro-turbulence velocity = 2.0 km s^-1, length of convective bubble divided by pressure scale height = 0.50, no convective overshooting, macro-turbulence velocity = 2.0 km s^-1, and resolving power R = \lambda / \delta \lambda = 250,000. First column: wavelength \lambda (in nm; from 830.1673264 to 889.8217922). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1). The 34698 lines have an average wavelength step of 0.00172 nm; the spectrum extent is thus 59.7 nm

            #Source: ESA, 20 June 2005, 'Photon-flux distributions for reference stars', GAIA-EST-TN-00539, issue 1, revision 0, based on D. Katz, priv. comm., 11 May 2005<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__K1IIIMPSTAR_SPECTRUM_NUMBEROFPHOTONSHIGHRESOLUTION


        __K3IIISTAR_SPECTRUM_NUMBEROFPHOTONS: str  = "Nature/K3IIIStar_Spectrum_NumberOfPhotons_001.fits"

        @property
        def K3IIISTAR_SPECTRUM_NUMBEROFPHOTONS(self):
            r"""
            Photon flux density N_{\lambda}(\lambda) of an unreddened K3III star (Pickles' star number 087) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)

            #Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__K3IIISTAR_SPECTRUM_NUMBEROFPHOTONS


        __K3VSTAR_SPECTRUM_NUMBEROFPHOTONS: str  = "Nature/K3VStar_Spectrum_NumberOfPhotons_001.fits"

        @property
        def K3VSTAR_SPECTRUM_NUMBEROFPHOTONS(self):
            r"""
            Photon flux density N_{\lambda}(\lambda) of an unreddened K3V star (Pickles' star number 034) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)

            #Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__K3VSTAR_SPECTRUM_NUMBEROFPHOTONS


        __L2_ALPHA: float  = 0.0100447147

        @property
        def L2_ALPHA(self):
            r"""
            Central auxiliary parameter in GAIA-FM-011, issue 1, revision 0

            #Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Equation 22<br/>
            #Basic : false
            #Scalar: true
            #Unit: []
            """

            return self.__L2_ALPHA


        __L2_BETA: float  = 3.18722929

        @property
        def L2_BETA(self):
            r"""
            Axis ratio of the 'horizontal' elliptic motion around L2 (Equations 70-71 in GAIA-FM-011, issue 1, revision 0)

            #Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Equation 66; see also the addition to GAIA-FM-011, issue 1, revision 0, dated 25 September 2003<br/>
            #Basic : false
            #Scalar: true
            #Unit: []
            """

            return self.__L2_BETA


        __L2_CAPITALOMEGA: float  = 3.940522185

        @property
        def L2_CAPITALOMEGA(self):
            r"""
            Auxiliary variable in GAIA-FM-011, issue 1, revision 0

            #Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Equations 34 and 38; see also the addition to GAIA-FM-011, issue 1, revision 0, dated 25 September 2003<br/>
            #Basic : false
            #Scalar: true
            #Unit: []
            """

            return self.__L2_CAPITALOMEGA


        __L2_MU: float  = 3.040423402e-06

        @property
        def L2_MU(self):
            r"""
            Reduced mass \mu of the Sun and Earth-Moon system as used in GAIA-FM-011, issue 1, revision 0. Note that the INPOP10e value 328900.552289 gives \mu = 3.04042347E-6

            #Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Section 2<br/>
            #Basic : false
            #Scalar: true
            #Unit: []
            """

            return self.__L2_MU


        __L2_OMEGA: float  = 1.985074856 # [{sidereal year}^-1]

        @property
        def L2_OMEGA(self):
            r"""
            Frequency of the vertical oscillation around L2

            #Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Equation 40; see also the addition to GAIA-FM-011, issue 1, revision 0, dated 25 September 2003<br/>
            #Basic : false
            #Scalar: true
            #Unit: [{sidereal year}^-1]
            """

            return self.__L2_OMEGA


        __L2_OMEGA_PERIOD: float  = 184.00 # [day]

        @property
        def L2_OMEGA_PERIOD(self):
            r"""
            Synodic period of the vertical oscillation around L2 in units of days

            #Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Equation 40<br/>
            #Basic : false
            #Scalar: true
            #Unit: [day]
            """

            return self.__L2_OMEGA_PERIOD


        __L2_ORBITALECCENTRICITY_J2000: float  = 0.01671123

        @property
        def L2_ORBITALECCENTRICITY_J2000(self):
            r"""
            Mean orbital eccentricity of the L2 orbit of the Sun and Earth-Moon system, at the standard epoch J2000.0. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0. DE405 is based upon the International Celestial Reference Frame (ICRF)

            #Source: See, e.g., F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Section 3.2, pages 6-7<br/>
            #Basic : false
            #Scalar: true
            #Unit: []
            """

            return self.__L2_ORBITALECCENTRICITY_J2000


        __L2_ORBITALSEMIMAJORAXIS_J2000: float  = 1.01008088 # [au]

        @property
        def L2_ORBITALSEMIMAJORAXIS_J2000(self):
            r"""
            Mean orbital semi-major axis of the L2 orbit of the Sun and Earth-Moon system, at the standard epoch J2000.0. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0. DE405 is based upon the International Celestial Reference Frame (ICRF)

            #Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Section 3<br/>
            #Basic : false
            #Scalar: true
            #Unit: [au]
            """

            return self.__L2_ORBITALSEMIMAJORAXIS_J2000


        __L2_P: float  = 27.0

        @property
        def L2_P(self):
            r"""
            The quantity p_n in GAIA-FM-011, issue 1, revision 0; this choice guarantees an eclipse-free orbit around L2 for more than 6 years

            #Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Tables 1 and 2 and Sections 6.2 and 6.4<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__L2_P


        __L2_PENUMBRARADIUS: float  = 13923.0 # [km]

        @property
        def L2_PENUMBRARADIUS(self):
            r"""
            Maximum radius of the penumbra of the Earth during a solar eclipse as seen from a point 1.0E5 km 'behind' L2. GAIA-FM-011, issue 1, revision 0 uses a rounded value of 14000 km

            #Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Section 4 (symbol \sigma_2) and Equation 86 (symbol s)<br/>
            #Basic : false
            #Scalar: true
            #Unit: [km]
            """

            return self.__L2_PENUMBRARADIUS


        __L2_Q: float  = 26.0

        @property
        def L2_Q(self):
            r"""
            The quantity q_n in GAIA-FM-011, issue 1, revision 0; this choice guarantees an eclipse-free orbit around L2 for more than 6 years

            #Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Tables 1 and 2 and Sections 6.2 and 6.4<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__L2_Q


        __L2_QTIMESEPSILON: float  = 5.78e-02

        @property
        def L2_QTIMESEPSILON(self):
            r"""
            The 'residual quantity' p - q a = q \epsilon, with a = \sigma / \omega, in GAIA-FM-011, issue 1, revision 0

            #Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Table 1 and Section 6.4<br/>
            #Basic : false
            #Scalar: true
            #Unit: []
            """

            return self.__L2_QTIMESEPSILON


        __L2_RHO: float  = 0.01007824044

        @property
        def L2_RHO(self):
            r"""
            Normalised separation between L2 and the Earth-Moon barycentre

            #Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Equation 23; see also the addition to GAIA-FM-011, issue 1, revision 0, dated 25 September 2003<br/>
            #Basic : false
            #Scalar: true
            #Unit: []
            """

            return self.__L2_RHO


        __L2_SIGMA: float  = 2.057014191 # [{sidereal year}^-1]

        @property
        def L2_SIGMA(self):
            r"""
            Frequency of the oscillation around L2 in the horizontal plane

            #Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Equations 49 and 52; see also the addition to GAIA-FM-011, issue 1, revision 0, dated 25 September 2003<br/>
            #Basic : false
            #Scalar: true
            #Unit: [{sidereal year}^-1]
            """

            return self.__L2_SIGMA


        __L2_SIGMA_PERIOD: float  = 177.57 # [day]

        @property
        def L2_SIGMA_PERIOD(self):
            r"""
            Synodic period of the 'horizontal' elliptic motion around L2 (Equations 70-71 in GAIA-FM-011, issue 1, revision 0) in units of days

            #Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Section 5.4<br/>
            #Basic : false
            #Scalar: true
            #Unit: [day]
            """

            return self.__L2_SIGMA_PERIOD


        __LSR_ABERRATION_MAXIMUM: float  = 4.2 # [10^-6 arcsec yr^-1]

        @property
        def LSR_ABERRATION_MAXIMUM(self):
            r"""
            Maximum (change in) aberration brought about by the acceleration of the LSR relative to the Galactic centre (resulting, if not corrected for, in spurious (apparent) proper motions for extra-Galactic sources in some regions of the sky)

            #Source: J. Kovalevsky, 2003, 'Aberration in proper motions', A&A, 404, 743. See also ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 1.8.10, page 110<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec yr^-1]
            """

            return self.__LSR_ABERRATION_MAXIMUM


        __LSR_ANGULARVELOCITY: float  = 27.19 # [km s^-1 kpc^-1]

        @property
        def LSR_ANGULARVELOCITY(self):
            r"""
            Angular velocity of circular rotation at the solar Galactocentric radius. Note: Transactions of the IAU, Volume XIXB, 1985, page 254: \Omega_0 = 220/8.5 = 25.88 km s^-1 kpc^-1

            #Basic : false
            #Scalar: true
            #Unit: [km s^-1 kpc^-1]
            """

            return self.__LSR_ANGULARVELOCITY


        __LSR_CIRCULARVELOCITY: float  = 217.520 # [km s^-1]

        @property
        def LSR_CIRCULARVELOCITY(self):
            r"""
            Velocity of circular rotation at the solar Galactocentric radius (local circular speed). Note: Transactions of the IAU, Volume XIXB, 1985, page 254: V_0 = 220 km s^-1

            #Basic : false
            #Scalar: true
            #Unit: [km s^-1]
            """

            return self.__LSR_CIRCULARVELOCITY


        __LSR_GALACTICROTATIONPERIOD: float  = 225.95 # [Myr]

        @property
        def LSR_GALACTICROTATIONPERIOD(self):
            r"""
            Period of rotation around the Galactic centre at the solar Galactocentric radius for a circular orbit

            #Basic : false
            #Scalar: true
            #Unit: [Myr]
            """

            return self.__LSR_GALACTICROTATIONPERIOD


        __LSR_GALACTOCENTRICRADIUS: float  = 8.0 # [kpc]

        @property
        def LSR_GALACTOCENTRICRADIUS(self):
            r"""
            Distance from the Local Standard of Rest (LSR, or: Sun) to the Galactic centre (solar Galactocentric radius)

            #Source: Z.M. Malkin, 28 February 2012, 'The current best estimate of the Galactocentric distance of the Sun based on comparison of different statistical techniques', http://adsabs.harvard.edu/abs/2012arXiv1202.6128M; see also Z.M. Malkin, 2013, Astronomicheskii Zhurnal, Volume 90, Number 2, pages 152-157 and Z.M. Malkin, 1 February 2013, 'Analysis of determinations of the distance between the sun and the galactic center', Astronomy Reports, Volume 57, Issue 2, pages 128-133<br/>
            #Basic : true
            #Scalar: true
            #Unit: [kpc]
            """

            return self.__LSR_GALACTOCENTRICRADIUS


        __LSUBB_CONSTANT: float  = 1.5505197680e-08

        @property
        def LSUBB_CONSTANT(self):
            r"""
            Average value of 1-d(TT)/d(TCB), defined as a defining constant, based on the '2006 best estimate' of LSubC_Constant + LSubG_Constant - LSubC_Constant * LSubG_Constant

            #Source: IAU, August 2006, 'Re-definition of Barycentric Dynamical Time, TDB', IAU 2006 Resolution 3 adopted at the XXVI-th General Assembly of the IAU. See also the IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__LSUBB_CONSTANT


        __LSUBC_CONSTANT: float  = 1.480826867410e-08

        @property
        def LSUBC_CONSTANT(self):
            r"""
            Average value of 1-d(TCG)/d(TCB)

            #Source: A. Irwin, T. Fukushima, 1999, A&A, 348, 642. See also the IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__LSUBC_CONSTANT


        __LSUBG_CONSTANT: float  = 6.9692901340e-10

        @property
        def LSUBG_CONSTANT(self):
            r"""
            Value of 1-d(TT)/d(TCG) (defining constant, based on best-estimate EarthEllipsoid_WSub0 = 6.26368556E7 m^2 s^-2)

            #Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html). See also IAU, August 2000, 'Re-definition of Terrestrial Time, TT', IAU 2000 Resolution B1.9 adopted at the XXIV-th General Assembly of the IAU<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__LSUBG_CONSTANT


        __LIGHTYEAR_ASTRONOMICALUNIT: float  = 63241.077084 # [au]

        @property
        def LIGHTYEAR_ASTRONOMICALUNIT(self):
            r"""
            Light year expressed in au

            #Basic : false
            #Scalar: true
            #Unit: [au]
            """

            return self.__LIGHTYEAR_ASTRONOMICALUNIT


        __LIGHTYEAR_METER: float  = 9.4607304726e+15 # [m]

        @property
        def LIGHTYEAR_METER(self):
            r"""
            Light year expressed in meters

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__LIGHTYEAR_METER


        __LIGHTYEAR_PARSEC: float  = 0.306601 # [pc]

        @property
        def LIGHTYEAR_PARSEC(self):
            r"""
            Light year expressed in parsecs

            #Basic : false
            #Scalar: true
            #Unit: [pc]
            """

            return self.__LIGHTYEAR_PARSEC


        __M0IIISTAR_SPECTRUM_NUMBEROFPHOTONS: str  = "Nature/M0IIIStar_Spectrum_NumberOfPhotons_001.fits"

        @property
        def M0IIISTAR_SPECTRUM_NUMBEROFPHOTONS(self):
            r"""
            Photon flux density N_{\lambda}(\lambda) of an unreddened M0III star (Pickles' star number 095) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)

            #Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__M0IIISTAR_SPECTRUM_NUMBEROFPHOTONS


        __M0VSTAR_SPECTRUM_NUMBEROFPHOTONS: str  = "Nature/M0VStar_Spectrum_NumberOfPhotons_001.fits"

        @property
        def M0VSTAR_SPECTRUM_NUMBEROFPHOTONS(self):
            r"""
            Photon flux density N_{\lambda}(\lambda) of an unreddened M0V star (Pickles' star number 038) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)

            #Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__M0VSTAR_SPECTRUM_NUMBEROFPHOTONS


        __M6VSTAR_SPECTRUM_NUMBEROFPHOTONS: str  = "Nature/M6VStar_Spectrum_NumberOfPhotons_001.fits"

        @property
        def M6VSTAR_SPECTRUM_NUMBEROFPHOTONS(self):
            r"""
            Photon flux density N_{\lambda}(\lambda) of an unreddened M6V star (Pickles' star number 045) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)

            #Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__M6VSTAR_SPECTRUM_NUMBEROFPHOTONS


        __M7IIISTAR_SPECTRUM_NUMBEROFPHOTONS: str  = "Nature/M7IIIStar_Spectrum_NumberOfPhotons_001.fits"

        @property
        def M7IIISTAR_SPECTRUM_NUMBEROFPHOTONS(self):
            r"""
            Photon flux density N_{\lambda}(\lambda) of an unreddened M7III star (Pickles' star number 102) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)

            #Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__M7IIISTAR_SPECTRUM_NUMBEROFPHOTONS


        __MAGNETIC_CONSTANT: float  = 1.256637061435917e-06 # [N A^-2]

        @property
        def MAGNETIC_CONSTANT(self):
            r"""
            Magnetic constant (defining constant)

            #Basic : false
            #Scalar: true
            #Unit: [N A^-2]
            """

            return self.__MAGNETIC_CONSTANT


        __MARSSYSTEM_ASTROMETRICSIGNATURE_10PARSEC: float  = 0.049 # [10^-6 arcsec]

        @property
        def MARSSYSTEM_ASTROMETRICSIGNATURE_10PARSEC(self):
            r"""
            Astrometric signature of the Sun induced by the Mars system for an observer located at a distance of 10 pc from the Sun

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.7, Equation 1.22, page 11<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__MARSSYSTEM_ASTROMETRICSIGNATURE_10PARSEC


        __MARSSYSTEM_MASS: float  = 6.41712e+23 # [kg]

        @property
        def MARSSYSTEM_MASS(self):
            r"""
            Mars-system mass (IAU 2009 CBE value). The planetary mass includes the contribution of its satellites

            #Basic : false
            #Scalar: true
            #Unit: [kg]
            """

            return self.__MARSSYSTEM_MASS


        __MARSSYSTEM_ORBITALECCENTRICITY_J2000: float  = 0.09339410

        @property
        def MARSSYSTEM_ORBITALECCENTRICITY_J2000(self):
            r"""
            Mean orbital eccentricity of Mars, at the standard epoch J2000.0. The mean orbital eccentricity is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__MARSSYSTEM_ORBITALECCENTRICITY_J2000


        __MARSSYSTEM_ORBITALINCLINATION_J2000: float  = 1.84969142 # [deg]

        @property
        def MARSSYSTEM_ORBITALINCLINATION_J2000(self):
            r"""
            Mean orbital inclination of Mars, at the standard epoch J2000.0. The mean orbital inclination is is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__MARSSYSTEM_ORBITALINCLINATION_J2000


        __MARSSYSTEM_ORBITALPERIOD: float  = 1.8808476 # [yr]

        @property
        def MARSSYSTEM_ORBITALPERIOD(self):
            r"""
            Sidereal orbital period

            #Source: Values derived from the mean longitude rates in Table 5.8.1 in P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, page 316; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [yr]
            """

            return self.__MARSSYSTEM_ORBITALPERIOD


        __MARSSYSTEM_ORBITALSEMIMAJORAXIS_J2000: float  = 1.52371034 # [au]

        @property
        def MARSSYSTEM_ORBITALSEMIMAJORAXIS_J2000(self):
            r"""
            Mean orbital semi-major axis of Mars, at the standard epoch J2000.0. The mean orbital semi-major axis is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: [au]
            """

            return self.__MARSSYSTEM_ORBITALSEMIMAJORAXIS_J2000


        __MARSSYSTEM_RADIALVELOCITYSIGNATURE: float  = 0.008 # [m s^-1]

        @property
        def MARSSYSTEM_RADIALVELOCITYSIGNATURE(self):
            r"""
            Radial-velocity amplitude of the Sun induced by the Mars system for 'an edge-on observer' (i.e., an observer in the orbital plane of the Mars system)

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.4, Equation 1.18, page 9<br/>
            #Basic : false
            #Scalar: true
            #Unit: [m s^-1]
            """

            return self.__MARSSYSTEM_RADIALVELOCITYSIGNATURE


        __MARS_ENCOMPASSINGSPHERERADIUS: float  = 3.396190e+06 # [m]

        @property
        def MARS_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around Mars which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MARS_ENCOMPASSINGSPHERERADIUS


        __MARS_EQUATORIALRADIUS: float  = 3.396190e+06 # [m]

        @property
        def MARS_EQUATORIALRADIUS(self):
            r"""
            Equatorial radius of Mars

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MARS_EQUATORIALRADIUS


        __MARS_FLATTENING: float  = 5.89790e-03

        @property
        def MARS_FLATTENING(self):
            r"""
            Geometrical flattening factor f of Mars (f = (a-b)/a). Mars has a significant dichotomy in shape between the northern and southern hemispheres

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__MARS_FLATTENING


        __MARS_FLUXREDUCTION_MAXIMUM: float  = 0.002 # [%]

        @property
        def MARS_FLUXREDUCTION_MAXIMUM(self):
            r"""
            Maximum reduction of the solar flux for an observer external to the solar system during a transit of Mars

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.2, Equation 2.4, page 14<br/>
            #Basic : false
            #Scalar: true
            #Unit: [%]
            """

            return self.__MARS_FLUXREDUCTION_MAXIMUM


        __MARS_GEOMETRICALBEDO: float  = 0.150

        @property
        def MARS_GEOMETRICALBEDO(self):
            r"""
            Geometric albedo of Mars (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__MARS_GEOMETRICALBEDO


        __MARS_JSUB2: float  = 1.9640e-03

        @property
        def MARS_JSUB2(self):
            r"""
            Dynamical form-factor of Mars (oblateness or Stokes' second-degree zonal harmonic of the potential)

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__MARS_JSUB2


        __MARS_LIGHTDEFLECTION_LIMB: float  = 116.0 # [10^-6 arcsec]

        @property
        def MARS_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of Mars

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__MARS_LIGHTDEFLECTION_LIMB


        __MARS_MASS: float  = 6.416930e+23 # [kg]

        @property
        def MARS_MASS(self):
            r"""
            Mass of Mars (do not use for high-precision (orbit) calculations)

            #Source: R.A. Jacobson, 2008, 'Ephemerides of the Martian Satellites - MAR080', JPL IOM 343.R-08-006; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [kg]
            """

            return self.__MARS_MASS


        __MARS_MASSDENSITY_MEAN: float  = 3.9340 # [g cm^-3]

        @property
        def MARS_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of Mars

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__MARS_MASSDENSITY_MEAN


        __MARS_NORTHROTATIONALPOLE_DECLINATION: float  = 52.88650 # [deg]

        @property
        def MARS_NORTHROTATIONALPOLE_DECLINATION(self):
            r"""
            IAU-recommended value for the declination \delta_0 of the north pole of rotation of Mars. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__MARS_NORTHROTATIONALPOLE_DECLINATION


        __MARS_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE: float  = -0.0000016674 # [deg day^-1]

        @property
        def MARS_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the declination \delta_0 of the north pole of rotation of Mars. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : false
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__MARS_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE


        __MARS_NORTHROTATIONALPOLE_RIGHTASCENSION: float  = 317.68143 # [deg]

        @property
        def MARS_NORTHROTATIONALPOLE_RIGHTASCENSION(self):
            r"""
            IAU-recommended value for the right ascension \alpha_0 of the north pole of rotation of Mars. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__MARS_NORTHROTATIONALPOLE_RIGHTASCENSION


        __MARS_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE: float  = -0.0000029049 # [deg day^-1]

        @property
        def MARS_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the right ascension \alpha_0 of the north pole of rotation of Mars. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : false
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__MARS_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE


        __MARS_POLARRADIUS: float  = 3.376160e+06 # [m]

        @property
        def MARS_POLARRADIUS(self):
            r"""
            Polar radius of Mars. Mars has a significant dichotomy in shape between the northern and southern hemispheres: the average polar radius is listed as 3.37620E6 m in B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MARS_POLARRADIUS


        __MARS_PRIMEMERIDIAN_EPHEMERISPOSITION: float  = 176.630 # [deg]

        @property
        def MARS_PRIMEMERIDIAN_EPHEMERISPOSITION(self):
            r"""
            IAU-recommended value for the ephemeris position of the prime meridian of Mars. The 0-deg meridian is defined by the crater Airy-0. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__MARS_PRIMEMERIDIAN_EPHEMERISPOSITION


        __MARS_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE: float  = 350.89198226 # [deg day^-1]

        @property
        def MARS_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch of 1.5 January 2000 = JD2451545.0 TDB) of the ephemeris position of the prime meridian of Mars. The 0-deg meridian is defined by the crater Airy-0. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__MARS_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE


        __MARS_TRANSITPROBABILITY: float  = 0.307 # [%]

        @property
        def MARS_TRANSITPROBABILITY(self):
            r"""
            Geometric transit probability (Mars transiting the Sun) for an observer external to the solar system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.1, Equation 2.2, page 14<br/>
            #Basic : false
            #Scalar: true
            #Unit: [%]
            """

            return self.__MARS_TRANSITPROBABILITY


        __MARS_TRANSITTIME_MAXIMUM: float  = 0.67 # [day]

        @property
        def MARS_TRANSITTIME_MAXIMUM(self):
            r"""
            Maximum transit time of Mars (transiting the Sun) for an observer external to the solar system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.3, Equation 2.5, page 15<br/>
            #Basic : false
            #Scalar: true
            #Unit: [day]
            """

            return self.__MARS_TRANSITTIME_MAXIMUM


        __MARS_VONEZEROMAGNITUDE: float  = -1.52 # [mag]

        @property
        def MARS_VONEZEROMAGNITUDE(self):
            r"""
            V(1,0) magnitude of Mars (i.e., the visual magnitude of the planet reduced to a distance of 1 au from both the Sun and Earth and phase angle zero). This parameter is also refered to as absolute magnitude in planetary sciences

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__MARS_VONEZEROMAGNITUDE


        __MARS_VOLUMETRICRADIUS: float  = 3.389500e+06 # [m]

        @property
        def MARS_VOLUMETRICRADIUS(self):
            r"""
            Mean volumetric radius of Mars

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__MARS_VOLUMETRICRADIUS


        __MERCURYSYSTEM_ASTROMETRICSIGNATURE_10PARSEC: float  = 0.006 # [10^-6 arcsec]

        @property
        def MERCURYSYSTEM_ASTROMETRICSIGNATURE_10PARSEC(self):
            r"""
            Astrometric signature of the Sun induced by Mercury for an observer located at a distance of 10 pc from the Sun

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.7, Equation 1.22, page 11<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__MERCURYSYSTEM_ASTROMETRICSIGNATURE_10PARSEC


        __MERCURYSYSTEM_MASS: float  = 3.3011e+23 # [kg]

        @property
        def MERCURYSYSTEM_MASS(self):
            r"""
            Mercury(-system) mass (IAU 2009 CBE value)

            #Basic : false
            #Scalar: true
            #Unit: [kg]
            """

            return self.__MERCURYSYSTEM_MASS


        __MERCURYSYSTEM_ORBITALECCENTRICITY_J2000: float  = 0.20563593

        @property
        def MERCURYSYSTEM_ORBITALECCENTRICITY_J2000(self):
            r"""
            Mean orbital eccentricity of Mercury, at the standard epoch J2000.0. The mean orbital eccentricity is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__MERCURYSYSTEM_ORBITALECCENTRICITY_J2000


        __MERCURYSYSTEM_ORBITALINCLINATION_J2000: float  = 7.00497902 # [deg]

        @property
        def MERCURYSYSTEM_ORBITALINCLINATION_J2000(self):
            r"""
            Mean orbital inclination of Mercury, at the standard epoch J2000.0. The mean orbital inclination is is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__MERCURYSYSTEM_ORBITALINCLINATION_J2000


        __MERCURYSYSTEM_ORBITALPERIOD: float  = 0.2408467 # [yr]

        @property
        def MERCURYSYSTEM_ORBITALPERIOD(self):
            r"""
            Sidereal orbital period

            #Source: Values derived from the mean longitude rates in Table 5.8.1 in P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, page 316; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [yr]
            """

            return self.__MERCURYSYSTEM_ORBITALPERIOD


        __MERCURYSYSTEM_ORBITALSEMIMAJORAXIS_J2000: float  = 0.38709927 # [au]

        @property
        def MERCURYSYSTEM_ORBITALSEMIMAJORAXIS_J2000(self):
            r"""
            Mean orbital semi-major axis of Mercury, at the standard epoch J2000.0. The mean orbital semi-major axis is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: [au]
            """

            return self.__MERCURYSYSTEM_ORBITALSEMIMAJORAXIS_J2000


        __MERCURYSYSTEM_RADIALVELOCITYSIGNATURE: float  = 0.008 # [m s^-1]

        @property
        def MERCURYSYSTEM_RADIALVELOCITYSIGNATURE(self):
            r"""
            Radial-velocity amplitude of the Sun induced by Mercury for 'an edge-on observer' (i.e., an observer in the orbital plane of Mercury)

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.4, Equation 1.18, page 9<br/>
            #Basic : false
            #Scalar: true
            #Unit: [m s^-1]
            """

            return self.__MERCURYSYSTEM_RADIALVELOCITYSIGNATURE


        __MERCURY_ENCOMPASSINGSPHERERADIUS: float  = 2.4397e+06 # [m]

        @property
        def MERCURY_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around Mercury which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MERCURY_ENCOMPASSINGSPHERERADIUS


        __MERCURY_EQUATORIALRADIUS: float  = 2.43970e+06 # [m]

        @property
        def MERCURY_EQUATORIALRADIUS(self):
            r"""
            Equatorial radius of Mercury

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MERCURY_EQUATORIALRADIUS


        __MERCURY_FLATTENING: float  = 0.0

        @property
        def MERCURY_FLATTENING(self):
            r"""
            Geometrical flattening factor f of Mercury (f = (a-b)/a)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__MERCURY_FLATTENING


        __MERCURY_FLUXREDUCTION_MAXIMUM: float  = 0.001 # [%]

        @property
        def MERCURY_FLUXREDUCTION_MAXIMUM(self):
            r"""
            Maximum reduction of the solar flux for an observer external to the solar system during a transit of Mercury

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.2, Equation 2.4, page 14<br/>
            #Basic : false
            #Scalar: true
            #Unit: [%]
            """

            return self.__MERCURY_FLUXREDUCTION_MAXIMUM


        __MERCURY_GEOMETRICALBEDO: float  = 0.106

        @property
        def MERCURY_GEOMETRICALBEDO(self):
            r"""
            Geometric albedo of Mercury (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__MERCURY_GEOMETRICALBEDO


        __MERCURY_JSUB2: float  = 6.00e-05

        @property
        def MERCURY_JSUB2(self):
            r"""
            Dynamical form-factor of Mercury (oblateness or Stokes' second-degree zonal harmonic of the potential)

            #Source: J.D. Anderson, G. Colombo, P.B. Esposito, E.L. Lau, G.B. Trager, 1 September 1987, 'The mass, gravity field, and ephemeris of Mercury', Icarus, 71, 337-349<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__MERCURY_JSUB2


        __MERCURY_LIGHTDEFLECTION_LIMB: float  = 83.0 # [10^-6 arcsec]

        @property
        def MERCURY_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of Mercury

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__MERCURY_LIGHTDEFLECTION_LIMB


        __MERCURY_MASS: float  = 3.301040e+23 # [kg]

        @property
        def MERCURY_MASS(self):
            r"""
            Mass of Mercury (do not use for high-precision (orbit) calculations)

            #Source: J.D. Anderson, et al., 1987, 'The mass, gravity field, and ephemeris of Mercury', Icarus, 71, 337-349; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [kg]
            """

            return self.__MERCURY_MASS


        __MERCURY_MASSDENSITY_MEAN: float  = 5.427 # [g cm^-3]

        @property
        def MERCURY_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of Mercury

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__MERCURY_MASSDENSITY_MEAN


        __MERCURY_NORTHROTATIONALPOLE_DECLINATION: float  = 61.4143 # [deg]

        @property
        def MERCURY_NORTHROTATIONALPOLE_DECLINATION(self):
            r"""
            IAU-recommended value for the declination \delta_0 of the north pole of rotation of Mercury. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__MERCURY_NORTHROTATIONALPOLE_DECLINATION


        __MERCURY_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE: float  = -0.0000001342 # [deg day^-1]

        @property
        def MERCURY_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the declination \delta_0 of the north pole of rotation of Mercury. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : false
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__MERCURY_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE


        __MERCURY_NORTHROTATIONALPOLE_RIGHTASCENSION: float  = 281.0097 # [deg]

        @property
        def MERCURY_NORTHROTATIONALPOLE_RIGHTASCENSION(self):
            r"""
            IAU-recommended value for the right ascension \alpha_0 of the north pole of rotation of Mercury. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__MERCURY_NORTHROTATIONALPOLE_RIGHTASCENSION


        __MERCURY_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE: float  = -0.0000008980 # [deg day^-1]

        @property
        def MERCURY_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the right ascension \alpha_0 of the north pole of rotation of Mercury. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : false
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__MERCURY_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE


        __MERCURY_POLARRADIUS: float  = 2.43970e+06 # [m]

        @property
        def MERCURY_POLARRADIUS(self):
            r"""
            Polar radius of Mercury

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MERCURY_POLARRADIUS


        __MERCURY_PRIMEMERIDIAN_EPHEMERISPOSITION: float  = 329.5469 # [deg]

        @property
        def MERCURY_PRIMEMERIDIAN_EPHEMERISPOSITION(self):
            r"""
            IAU-recommended value for the ephemeris position of the prime meridian of Mercury. The 20-deg meridian is defined by the crater Hun Kal. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__MERCURY_PRIMEMERIDIAN_EPHEMERISPOSITION


        __MERCURY_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE: float  = 6.1385025 # [deg day^-1]

        @property
        def MERCURY_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch of 1.5 January 2000 = JD2451545.0 TDB) of the ephemeris position of the prime meridian of Mercury. The 20-deg meridian is defined by the crater Hun Kal. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__MERCURY_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE


        __MERCURY_TRANSITPROBABILITY: float  = 1.206 # [%]

        @property
        def MERCURY_TRANSITPROBABILITY(self):
            r"""
            Geometric transit probability (Mercury transiting the Sun) for an observer external to the solar system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.1, Equation 2.2, page 14<br/>
            #Basic : false
            #Scalar: true
            #Unit: [%]
            """

            return self.__MERCURY_TRANSITPROBABILITY


        __MERCURY_TRANSITTIME_MAXIMUM: float  = 0.34 # [day]

        @property
        def MERCURY_TRANSITTIME_MAXIMUM(self):
            r"""
            Maximum transit time of Mercury (transiting the Sun) for an observer external to the solar system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.3, Equation 2.5, page 15<br/>
            #Basic : false
            #Scalar: true
            #Unit: [day]
            """

            return self.__MERCURY_TRANSITTIME_MAXIMUM


        __MERCURY_VONEZEROMAGNITUDE: float  = -0.60 # [mag]

        @property
        def MERCURY_VONEZEROMAGNITUDE(self):
            r"""
            V(1,0) magnitude of Mercury (i.e., the visual magnitude of the planet reduced to a distance of 1 au from both the Sun and Earth and phase angle zero). This parameter is also refered to as absolute magnitude in planetary sciences

            #Source: J.L. Hilton, 2005, 'Improving the Visual Magnitudes of the Planets in The Astronomical Almanac. I. Mercury and Venus', AJ, 129, 2902-2906; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__MERCURY_VONEZEROMAGNITUDE


        __MERCURY_VOLUMETRICRADIUS: float  = 2.43970e+06 # [m]

        @property
        def MERCURY_VOLUMETRICRADIUS(self):
            r"""
            Mean volumetric radius of Mercury

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__MERCURY_VOLUMETRICRADIUS


        __MICROARCSECOND_RADIAN: float  = 4.848136811095360e-12 # [rad]

        @property
        def MICROARCSECOND_RADIAN(self):
            r"""
            One micro-arcsecond in units of radians

            #Basic : false
            #Scalar: true
            #Unit: [rad]
            """

            return self.__MICROARCSECOND_RADIAN


        __MICROMETEOROID_FLUX_L2LARGEPARTICLES: float  = 2.60e-18 # [particles m^-2 s^-1 hemisphere^-1]

        @property
        def MICROMETEOROID_FLUX_L2LARGEPARTICLES(self):
            r"""
            The typical expected micro-meteoroid flux at L2, in units of particles cm^-2 s^-1 hemisphere^-1, resulting from particles with masses greater than m (in units of kg and for m > 2.8E-11 kg) equals 2.6E-18 * m^(-7/6)

            #Source: L. Lindegren, 13 July 2000, 'Effects of micro-meteoroids on Gaia attitude', GAIA-LL-031, issue 1, revision 0. See also K. Yamakoshi, 1994, 'Extraterrestrial dust; Laboratory studies of interplanetary dust', Astrophysics and Space Science Library, 181, Kluwer Academic Publishers, Dordrecht (1994edls.book.....Y)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [particles m^-2 s^-1 hemisphere^-1]
            """

            return self.__MICROMETEOROID_FLUX_L2LARGEPARTICLES


        __MICROMETEOROID_FLUX_L2SMALLPARTICLES: float  = 2.80e-11 # [particles m^-2 s^-1 hemisphere^-1]

        @property
        def MICROMETEOROID_FLUX_L2SMALLPARTICLES(self):
            r"""
            The typical expected micro-meteoroid flux at L2, in units of particles cm^-2 s^-1 hemisphere^-1, resulting from particles with masses greater than m (in units of kg and for m < 2.8E-11 kg) equals 2.8E-11 * m^(-1/2)

            #Source: L. Lindegren, 13 July 2000, 'Effects of micro-meteoroids on Gaia attitude', GAIA-LL-031, issue 1, revision 0. See also K. Yamakoshi, 1994, 'Extraterrestrial dust; Laboratory studies of interplanetary dust', Astrophysics and Space Science Library, 181, Kluwer Academic Publishers, Dordrecht (1994edls.book.....Y)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [particles m^-2 s^-1 hemisphere^-1]
            """

            return self.__MICROMETEOROID_FLUX_L2SMALLPARTICLES


        __MILLIARCSECOND_RADIAN: float  = 4.848136811095360e-09 # [rad]

        @property
        def MILLIARCSECOND_RADIAN(self):
            r"""
            One milli-arcsecond in units of radians

            #Basic : false
            #Scalar: true
            #Unit: [rad]
            """

            return self.__MILLIARCSECOND_RADIAN


        __MOLARGAS_CONSTANT: float  = 8.3144598 # [J mol^-1 K^-1]

        @property
        def MOLARGAS_CONSTANT(self):
            r"""
            Molar gas constant

            #Source: P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [J mol^-1 K^-1]
            """

            return self.__MOLARGAS_CONSTANT


        __MOONJ1_ENCOMPASSINGSPHERERADIUS: float  = 1.82149e+06 # [m]

        @property
        def MOONJ1_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around J1 (Io) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONJ1_ENCOMPASSINGSPHERERADIUS


        __MOONJ1_GM: float  = 5.9599160e+12 # [m^3 s^-2]

        @property
        def MOONJ1_GM(self):
            r"""
            GM of J1 (Io)

            #Source: R.A. Jacobson, 2003, 'Constants used in the JUP230 ephemeris', priv. comm.; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__MOONJ1_GM


        __MOONJ1_LIGHTDEFLECTION_LIMB: float  = 30.0 # [10^-6 arcsec]

        @property
        def MOONJ1_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of J1 (Io)

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__MOONJ1_LIGHTDEFLECTION_LIMB


        __MOONJ1_MASSDENSITY_MEAN: float  = 3.528 # [g cm^-3]

        @property
        def MOONJ1_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of J1 (Io)

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__MOONJ1_MASSDENSITY_MEAN


        __MOONJ1_RADIUS: float  = 1.821490e+06 # [m]

        @property
        def MOONJ1_RADIUS(self):
            r"""
            Radius of J1 (Io)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONJ1_RADIUS


        __MOONJ2_ENCOMPASSINGSPHERERADIUS: float  = 1.56080e+06 # [m]

        @property
        def MOONJ2_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around J2 (Europa) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONJ2_ENCOMPASSINGSPHERERADIUS


        __MOONJ2_GM: float  = 3.2027390e+12 # [m^3 s^-2]

        @property
        def MOONJ2_GM(self):
            r"""
            GM of J2 (Europa)

            #Source: R.A. Jacobson, 2003, 'Constants used in the JUP230 ephemeris', priv. comm.; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__MOONJ2_GM


        __MOONJ2_LIGHTDEFLECTION_LIMB: float  = 19.0 # [10^-6 arcsec]

        @property
        def MOONJ2_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of J2 (Europa)

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__MOONJ2_LIGHTDEFLECTION_LIMB


        __MOONJ2_MASSDENSITY_MEAN: float  = 3.013 # [g cm^-3]

        @property
        def MOONJ2_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of J2 (Europa)

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__MOONJ2_MASSDENSITY_MEAN


        __MOONJ2_RADIUS: float  = 1.56080e+06 # [m]

        @property
        def MOONJ2_RADIUS(self):
            r"""
            Radius of J2 (Europa)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONJ2_RADIUS


        __MOONJ3_ENCOMPASSINGSPHERERADIUS: float  = 2.63120e+06 # [m]

        @property
        def MOONJ3_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around J3 (Ganymede) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONJ3_ENCOMPASSINGSPHERERADIUS


        __MOONJ3_GM: float  = 9.8878340e+12 # [m^3 s^-2]

        @property
        def MOONJ3_GM(self):
            r"""
            GM of J3 (Ganymede)

            #Source: R.A. Jacobson, 2003, 'Constants used in the JUP230 ephemeris', priv. comm.; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__MOONJ3_GM


        __MOONJ3_LIGHTDEFLECTION_LIMB: float  = 34.0 # [10^-6 arcsec]

        @property
        def MOONJ3_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of J3 (Ganymede)

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__MOONJ3_LIGHTDEFLECTION_LIMB


        __MOONJ3_MASSDENSITY_MEAN: float  = 1.942 # [g cm^-3]

        @property
        def MOONJ3_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of J3 (Ganymede)

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__MOONJ3_MASSDENSITY_MEAN


        __MOONJ3_RADIUS: float  = 2.63120e+06 # [m]

        @property
        def MOONJ3_RADIUS(self):
            r"""
            Radius of J3 (Ganymede)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONJ3_RADIUS


        __MOONJ4_ENCOMPASSINGSPHERERADIUS: float  = 2.41030e+06 # [m]

        @property
        def MOONJ4_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around J4 (Callisto) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONJ4_ENCOMPASSINGSPHERERADIUS


        __MOONJ4_GM: float  = 7.1792890e+12 # [m^3 s^-2]

        @property
        def MOONJ4_GM(self):
            r"""
            GM of J4 (Callisto)

            #Source: R.A. Jacobson, 2003, 'Constants used in the JUP230 ephemeris', priv. comm.; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__MOONJ4_GM


        __MOONJ4_LIGHTDEFLECTION_LIMB: float  = 27.0 # [10^-6 arcsec]

        @property
        def MOONJ4_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of J4 (Callisto)

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__MOONJ4_LIGHTDEFLECTION_LIMB


        __MOONJ4_MASSDENSITY_MEAN: float  = 1.834 # [g cm^-3]

        @property
        def MOONJ4_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of J4 (Callisto)

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__MOONJ4_MASSDENSITY_MEAN


        __MOONJ4_RADIUS: float  = 2.41030e+06 # [m]

        @property
        def MOONJ4_RADIUS(self):
            r"""
            Radius of J4 (Callisto)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONJ4_RADIUS


        __MOONN1_ENCOMPASSINGSPHERERADIUS: float  = 1.35260e+06 # [m]

        @property
        def MOONN1_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around N1 (Triton) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONN1_ENCOMPASSINGSPHERERADIUS


        __MOONN1_GM: float  = 1.42760e+12 # [m^3 s^-2]

        @property
        def MOONN1_GM(self):
            r"""
            GM of N1 (Triton)

            #Source: R.A. Jacobson, 2009, 'The Orbits of the Neptunian Satellites and the Orientation of the Pole of Neptune', AJ, 137, 4322; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__MOONN1_GM


        __MOONN1_LIGHTDEFLECTION_LIMB: float  = 10.0 # [10^-6 arcsec]

        @property
        def MOONN1_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of N1 (Triton)

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__MOONN1_LIGHTDEFLECTION_LIMB


        __MOONN1_MASSDENSITY_MEAN: float  = 2.064 # [g cm^-3]

        @property
        def MOONN1_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of N1 (Triton)

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__MOONN1_MASSDENSITY_MEAN


        __MOONN1_RADIUS: float  = 1.35260e+06 # [m]

        @property
        def MOONN1_RADIUS(self):
            r"""
            Radius of N1 (Triton)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONN1_RADIUS


        __MOONS3_ENCOMPASSINGSPHERERADIUS: float  = 5.3100e+05 # [m]

        @property
        def MOONS3_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around S3 (Tethys) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONS3_ENCOMPASSINGSPHERERADIUS


        __MOONS3_GM: float  = 4.120670e+10 # [m^3 s^-2]

        @property
        def MOONS3_GM(self):
            r"""
            GM of S3 (Tethys)

            #Source: R.A. Jacobson, et al., 2008, 'The Gravity Gield of the Saturnian System and the Orbits of the Major Saturnian satellites', Presented at the Saturn After Cassini-Huygens Symposium held at Imperial College London, UK, 28 July - 1 August 2008; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__MOONS3_GM


        __MOONS3_LIGHTDEFLECTION_LIMB: float  = 1.0 # [10^-6 arcsec]

        @property
        def MOONS3_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of S3 (Tethys)

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__MOONS3_LIGHTDEFLECTION_LIMB


        __MOONS3_MASSDENSITY_MEAN: float  = 0.984 # [g cm^-3]

        @property
        def MOONS3_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of S3 (Tethys)

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__MOONS3_MASSDENSITY_MEAN


        __MOONS3_RADIUS: float  = 5.3100e+05 # [m]

        @property
        def MOONS3_RADIUS(self):
            r"""
            Radius of S3 (Tethys)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONS3_RADIUS


        __MOONS4_ENCOMPASSINGSPHERERADIUS: float  = 5.614e+05 # [m]

        @property
        def MOONS4_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around S4 (Dione) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONS4_ENCOMPASSINGSPHERERADIUS


        __MOONS4_GM: float  = 7.311460e+10 # [m^3 s^-2]

        @property
        def MOONS4_GM(self):
            r"""
            GM of S4 (Dione)

            #Source: R.A. Jacobson, et al., 2008, 'The Gravity Gield of the Saturnian System and the Orbits of the Major Saturnian satellites', Presented at the Saturn After Cassini-Huygens Symposium held at Imperial College London, UK, 28 July - 1 August 2008; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__MOONS4_GM


        __MOONS4_LIGHTDEFLECTION_LIMB: float  = 1.0 # [10^-6 arcsec]

        @property
        def MOONS4_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of S4 (Dione)

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__MOONS4_LIGHTDEFLECTION_LIMB


        __MOONS4_MASSDENSITY_MEAN: float  = 1.478 # [g cm^-3]

        @property
        def MOONS4_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of S4 (Dione)

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__MOONS4_MASSDENSITY_MEAN


        __MOONS4_RADIUS: float  = 5.6140e+05 # [m]

        @property
        def MOONS4_RADIUS(self):
            r"""
            Radius of S4 (Dione)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONS4_RADIUS


        __MOONS5_ENCOMPASSINGSPHERERADIUS: float  = 7.635e+05 # [m]

        @property
        def MOONS5_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around S5 (Rhea) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONS5_ENCOMPASSINGSPHERERADIUS


        __MOONS5_GM: float  = 1.5394260e+11 # [m^3 s^-2]

        @property
        def MOONS5_GM(self):
            r"""
            GM of S5 (Rhea)

            #Source: R.A. Jacobson, et al., 2008, 'The Gravity Gield of the Saturnian System and the Orbits of the Major Saturnian satellites', Presented at the Saturn After Cassini-Huygens Symposium held at Imperial College London, UK, 28 July - 1 August 2008; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__MOONS5_GM


        __MOONS5_LIGHTDEFLECTION_LIMB: float  = 2.0 # [10^-6 arcsec]

        @property
        def MOONS5_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of S5 (Rhea)

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__MOONS5_LIGHTDEFLECTION_LIMB


        __MOONS5_MASSDENSITY_MEAN: float  = 1.237 # [g cm^-3]

        @property
        def MOONS5_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of S5 (Rhea)

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__MOONS5_MASSDENSITY_MEAN


        __MOONS5_RADIUS: float  = 7.6350e+05 # [m]

        @property
        def MOONS5_RADIUS(self):
            r"""
            Radius of S5 (Rhea)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONS5_RADIUS


        __MOONS6_ENCOMPASSINGSPHERERADIUS: float  = 2.5747e+06 # [m]

        @property
        def MOONS6_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around S6 (Titan) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONS6_ENCOMPASSINGSPHERERADIUS


        __MOONS6_GM: float  = 8.97813820e+12 # [m^3 s^-2]

        @property
        def MOONS6_GM(self):
            r"""
            GM of S6 (Titan)

            #Source: R.A. Jacobson, et al., 2008, 'The Gravity Gield of the Saturnian System and the Orbits of the Major Saturnian satellites', Presented at the Saturn After Cassini-Huygens Symposium held at Imperial College London, UK, 28 July - 1 August 2008; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__MOONS6_GM


        __MOONS6_LIGHTDEFLECTION_LIMB: float  = 32.0 # [10^-6 arcsec]

        @property
        def MOONS6_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of S6 (Titan)

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__MOONS6_LIGHTDEFLECTION_LIMB


        __MOONS6_MASSDENSITY_MEAN: float  = 1.882 # [g cm^-3]

        @property
        def MOONS6_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of S6 (Titan)

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__MOONS6_MASSDENSITY_MEAN


        __MOONS6_RADIUS: float  = 2.574730e+06 # [m]

        @property
        def MOONS6_RADIUS(self):
            r"""
            Radius of S6 (Titan)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONS6_RADIUS


        __MOONS8_ENCOMPASSINGSPHERERADIUS: float  = 7.343e+05 # [m]

        @property
        def MOONS8_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around S8 (Iapetus) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONS8_ENCOMPASSINGSPHERERADIUS


        __MOONS8_GM: float  = 1.2050380e+11 # [m^3 s^-2]

        @property
        def MOONS8_GM(self):
            r"""
            GM of S8 (Iapetus)

            #Source: R.A. Jacobson, et al., 2008, 'The Gravity Gield of the Saturnian System and the Orbits of the Major Saturnian satellites', Presented at the Saturn After Cassini-Huygens Symposium held at Imperial College London, UK, 28 July - 1 August 2008; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__MOONS8_GM


        __MOONS8_LIGHTDEFLECTION_LIMB: float  = 2.0 # [10^-6 arcsec]

        @property
        def MOONS8_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of S8 (Iapetus)

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__MOONS8_LIGHTDEFLECTION_LIMB


        __MOONS8_MASSDENSITY_MEAN: float  = 1.089 # [g cm^-3]

        @property
        def MOONS8_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of S8 (Iapetus)

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__MOONS8_MASSDENSITY_MEAN


        __MOONS8_RADIUS: float  = 7.3430e+05 # [m]

        @property
        def MOONS8_RADIUS(self):
            r"""
            Radius of S8 (Iapetus)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONS8_RADIUS


        __MOONTOEARTH_MASSRATIO: float  = 0.0123000371

        @property
        def MOONTOEARTH_MASSRATIO(self):
            r"""
            Ratio of Moon to Earth mass

            #Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__MOONTOEARTH_MASSRATIO


        __MOONU1_ENCOMPASSINGSPHERERADIUS: float  = 5.7890e+05 # [m]

        @property
        def MOONU1_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around U1 (Ariel) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONU1_ENCOMPASSINGSPHERERADIUS


        __MOONU1_GM: float  = 8.640e+10 # [m^3 s^-2]

        @property
        def MOONU1_GM(self):
            r"""
            GM of U1 (Ariel)

            #Source: R.A. Jacobson, 2007, 'The Gravity Field of the Uranian System and the Orbits of the Uranian Satellites and Rings', BAAS, 39, 453; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__MOONU1_GM


        __MOONU1_LIGHTDEFLECTION_LIMB: float  = 1.0 # [10^-6 arcsec]

        @property
        def MOONU1_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of U1 (Ariel)

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__MOONU1_LIGHTDEFLECTION_LIMB


        __MOONU1_MASSDENSITY_MEAN: float  = 1.593 # [g cm^-3]

        @property
        def MOONU1_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of U1 (Ariel)

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__MOONU1_MASSDENSITY_MEAN


        __MOONU1_RADIUS: float  = 5.7890e+05 # [m]

        @property
        def MOONU1_RADIUS(self):
            r"""
            Radius of U1 (Ariel)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONU1_RADIUS


        __MOONU2_ENCOMPASSINGSPHERERADIUS: float  = 5.8470e+05 # [m]

        @property
        def MOONU2_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around U2 (Umbriel) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONU2_ENCOMPASSINGSPHERERADIUS


        __MOONU2_GM: float  = 8.150e+10 # [m^3 s^-2]

        @property
        def MOONU2_GM(self):
            r"""
            GM of U2 (Umbriel)

            #Source: R.A. Jacobson, 2007, 'The Gravity Field of the Uranian System and the Orbits of the Uranian Satellites and Rings', BAAS, 39, 453; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__MOONU2_GM


        __MOONU2_LIGHTDEFLECTION_LIMB: float  = 1.0 # [10^-6 arcsec]

        @property
        def MOONU2_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of U2 (Umbriel)

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__MOONU2_LIGHTDEFLECTION_LIMB


        __MOONU2_MASSDENSITY_MEAN: float  = 1.458 # [g cm^-3]

        @property
        def MOONU2_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of U2 (Umbriel)

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__MOONU2_MASSDENSITY_MEAN


        __MOONU2_RADIUS: float  = 5.8470e+05 # [m]

        @property
        def MOONU2_RADIUS(self):
            r"""
            Radius of U2 (Umbriel)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONU2_RADIUS


        __MOONU3_ENCOMPASSINGSPHERERADIUS: float  = 7.8890e+05 # [m]

        @property
        def MOONU3_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around U3 (Titania) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONU3_ENCOMPASSINGSPHERERADIUS


        __MOONU3_GM: float  = 2.2820e+11 # [m^3 s^-2]

        @property
        def MOONU3_GM(self):
            r"""
            GM of U3 (Titania)

            #Source: R.A. Jacobson, 2007, 'The Gravity Field of the Uranian System and the Orbits of the Uranian Satellites and Rings', BAAS, 39, 453; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__MOONU3_GM


        __MOONU3_LIGHTDEFLECTION_LIMB: float  = 3.0 # [10^-6 arcsec]

        @property
        def MOONU3_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of U3 (Titania)

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__MOONU3_LIGHTDEFLECTION_LIMB


        __MOONU3_MASSDENSITY_MEAN: float  = 1.663 # [g cm^-3]

        @property
        def MOONU3_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of U3 (Titania)

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__MOONU3_MASSDENSITY_MEAN


        __MOONU3_RADIUS: float  = 7.8890e+05 # [m]

        @property
        def MOONU3_RADIUS(self):
            r"""
            Radius of U3 (Titania)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONU3_RADIUS


        __MOONU4_ENCOMPASSINGSPHERERADIUS: float  = 7.6140e+05 # [m]

        @property
        def MOONU4_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around U4 (Oberon) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONU4_ENCOMPASSINGSPHERERADIUS


        __MOONU4_GM: float  = 1.9240e+11 # [m^3 s^-2]

        @property
        def MOONU4_GM(self):
            r"""
            GM of U4 (Oberon)

            #Source: R.A. Jacobson, 2007, 'The Gravity Field of the Uranian System and the Orbits of the Uranian Satellites and Rings', BAAS, 39, 453; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__MOONU4_GM


        __MOONU4_LIGHTDEFLECTION_LIMB: float  = 2.0 # [10^-6 arcsec]

        @property
        def MOONU4_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of U4 (Oberon)

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__MOONU4_LIGHTDEFLECTION_LIMB


        __MOONU4_MASSDENSITY_MEAN: float  = 1.559 # [g cm^-3]

        @property
        def MOONU4_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of U4 (Oberon)

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__MOONU4_MASSDENSITY_MEAN


        __MOONU4_RADIUS: float  = 7.6140e+05 # [m]

        @property
        def MOONU4_RADIUS(self):
            r"""
            Radius of U4 (Oberon)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOONU4_RADIUS


        __MOON_DIURNALPARALLAX: float  = 3422.595 # [arcsec]

        @property
        def MOON_DIURNALPARALLAX(self):
            r"""
            Lunar diurnal parallax. Formally, this parameter is defined as the ratio of a fictitious mean equatorial radius of the Earth to the perturbed mean distance of the Moon; the ratio F_2 of the perturbed mean distance to the Moon (the perturbation being due to the Sun) to the two-body mean distance of the Moon (with the Sun not present and constant mean motion) equals 0.999093141975298 (see T.D. Moyer, 15 May 1971, 'Mathematical formulation of the Double-Precision Orbit Determination Programme (DPODP)', NASA JPL Technical Report 32-1527, pages 25-26)

            #Basic : false
            #Scalar: true
            #Unit: [arcsec]
            """

            return self.__MOON_DIURNALPARALLAX


        __MOON_ENCOMPASSINGSPHERERADIUS: float  = 1738000.0 # [m]

        @property
        def MOON_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around the Moon which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOON_ENCOMPASSINGSPHERERADIUS


        __MOON_GM: float  = 4.9028002e+12 # [m^3 s^-2]

        @property
        def MOON_GM(self):
            r"""
            Selenocentric gravitational constant (TCB-compatible value)

            #Basic : false
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__MOON_GM


        __MOON_GEOMETRICALBEDO: float  = 0.12

        @property
        def MOON_GEOMETRICALBEDO(self):
            r"""
            Geometric albedo of the Moon (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__MOON_GEOMETRICALBEDO


        __MOON_MASS: float  = 7.3460e+22 # [kg]

        @property
        def MOON_MASS(self):
            r"""
            Lunar mass (do not use for high-precision (orbit) calculations)

            #Basic : false
            #Scalar: true
            #Unit: [kg]
            """

            return self.__MOON_MASS


        __MOON_MASSDENSITY_MEAN: float  = 3.344 # [g cm^-3]

        @property
        def MOON_MASSDENSITY_MEAN(self):
            r"""
            Mean Moon mass density

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__MOON_MASSDENSITY_MEAN


        __MOON_OPPOSITIONVMAGNITUDE: float  = -12.74 # [mag]

        @property
        def MOON_OPPOSITIONVMAGNITUDE(self):
            r"""
            Johnson V band mean opposition magnitude of the Moon

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__MOON_OPPOSITIONVMAGNITUDE


        __MOON_ORBITALECCENTRICITY_J2000: float  = 0.0554

        @property
        def MOON_ORBITALECCENTRICITY_J2000(self):
            r"""
            Eccentricity of Lunar orbit (mean ecliptic orbital elements, at the standard epoch J2000.0, based on JPL's Planetary and Lunar Ephemerides DE405/LE405)

            #Source: E.M. Standish, 2001, 'Approximate Mean Ecliptic Elements of the Lunar Orbit', JPL IOM 312.F-01-004 (http://ssd.jpl.nasa.gov/?sat_elem)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__MOON_ORBITALECCENTRICITY_J2000


        __MOON_ORBITALINCLINATION_J2000: float  = 5.16 # [deg]

        @property
        def MOON_ORBITALINCLINATION_J2000(self):
            r"""
            Inclination of Lunar orbit with respect to the ecliptic (mean ecliptic orbital elements, at the standard epoch J2000.0, based on JPL's Planetary and Lunar Ephemerides DE405/LE405)

            #Source: E.M. Standish, 2001, 'Approximate Mean Ecliptic Elements of the Lunar Orbit', JPL IOM 312.F-01-004 (http://ssd.jpl.nasa.gov/?sat_elem)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__MOON_ORBITALINCLINATION_J2000


        __MOON_ORBITALSEMIMAJORAXIS_J2000: float  = 3.844000e+08 # [m]

        @property
        def MOON_ORBITALSEMIMAJORAXIS_J2000(self):
            r"""
            Semi-major axis of Lunar orbit (mean ecliptic orbital elements, at the standard epoch J2000.0, based on JPL's Planetary and Lunar Ephemerides DE405/LE405)

            #Source: E.M. Standish, 2001, 'Approximate Mean Ecliptic Elements of the Lunar Orbit', JPL IOM 312.F-01-004 (http://ssd.jpl.nasa.gov/?sat_elem)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOON_ORBITALSEMIMAJORAXIS_J2000


        __MOON_PRECESSIONPERIOD_J2000ARGUMENTOFPERIAPSIS: float  = 5.997 # [yr]

        @property
        def MOON_PRECESSIONPERIOD_J2000ARGUMENTOFPERIAPSIS(self):
            r"""
            Precession period of the argument of periapsis of Lunar orbit, i.e., apsidal period (mean ecliptic orbital elements, at the standard epoch J2000.0, based on JPL's Planetary and Lunar Ephemerides DE405/LE405)

            #Source: E.M. Standish, 2001, 'Approximate Mean Ecliptic Elements of the Lunar Orbit', JPL IOM 312.F-01-004 (http://ssd.jpl.nasa.gov/?sat_elem)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [yr]
            """

            return self.__MOON_PRECESSIONPERIOD_J2000ARGUMENTOFPERIAPSIS


        __MOON_PRECESSIONPERIOD_J2000LONGITUDEOFASCENDINGNODE: float  = 18.600 # [yr]

        @property
        def MOON_PRECESSIONPERIOD_J2000LONGITUDEOFASCENDINGNODE(self):
            r"""
            Precession period of the longitude of the ascending node of Lunar orbit, i.e., nodal period (mean ecliptic orbital elements, at the standard epoch J2000.0, based on JPL's Planetary and Lunar Ephemerides DE405/LE405)

            #Source: E.M. Standish, 2001, 'Approximate Mean Ecliptic Elements of the Lunar Orbit', JPL IOM 312.F-01-004 (http://ssd.jpl.nasa.gov/?sat_elem)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [yr]
            """

            return self.__MOON_PRECESSIONPERIOD_J2000LONGITUDEOFASCENDINGNODE


        __MOON_SIDEREALPERIOD_J2000: float  = 27.322 # [day]

        @property
        def MOON_SIDEREALPERIOD_J2000(self):
            r"""
            Sidereal period of Lunar orbit (mean ecliptic orbital elements, at the standard epoch J2000.0, based on JPL's Planetary and Lunar Ephemerides DE405/LE405)

            #Source: E.M. Standish, 2001, 'Approximate Mean Ecliptic Elements of the Lunar Orbit', JPL IOM 312.F-01-004 (http://ssd.jpl.nasa.gov/?sat_elem)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [day]
            """

            return self.__MOON_SIDEREALPERIOD_J2000


        __MOON_VOLUMETRICRADIUS: float  = 1.73740e+06 # [m]

        @property
        def MOON_VOLUMETRICRADIUS(self):
            r"""
            Mean volumetric radius of the Moon

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__MOON_VOLUMETRICRADIUS


        __NAPIER_CONSTANT: float  = 2.71828182845904523536028747135

        @property
        def NAPIER_CONSTANT(self):
            r"""
            Napier's constant (also known as Neper's constant), i.e., base of the natural logarithm. Although the symbol 'e' refers to Euler, the Napier constant should not be confused with the Euler(-Mascheroni) constant \gamma = 0.5772156649... Note that double-precision, floating-point numbers in any programming language which follows the IEEE standard (true for C, C++, Java, and most, if not all, others) have only 16 significant digits (64 bits); the representation here, using 30 significant digits, is thus amply sufficient

            #Source: Well-known mathematical constant; numerical value can be extracted, e.g., from Mathematica 4.0 for Solaris (Wolfram Research, Inc.) using 'N[Exp[1],30]'<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__NAPIER_CONSTANT


        __NEAREARTHASTEROID1036GANYMED_DIAMETER: float  = 38.5 # [km]

        @property
        def NEAREARTHASTEROID1036GANYMED_DIAMETER(self):
            r"""
            Diameter of near-Earth asteroid 1036 Ganymed

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 16 'Prominent Aten-, Apollo-, and Amor-class near-Earth asteroids with well-determined orbits and V(1,0) < 18 mag' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [km]
            """

            return self.__NEAREARTHASTEROID1036GANYMED_DIAMETER


        __NEAREARTHASTEROID1580BETULIA_DIAMETER: float  = 7.4 # [km]

        @property
        def NEAREARTHASTEROID1580BETULIA_DIAMETER(self):
            r"""
            Diameter of near-Earth asteroid 1580 Betulia

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 16 'Prominent Aten-, Apollo-, and Amor-class near-Earth asteroids with well-determined orbits and V(1,0) < 18 mag' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [km]
            """

            return self.__NEAREARTHASTEROID1580BETULIA_DIAMETER


        __NEAREARTHASTEROID1627IVAR_DIAMETER: float  = 8.1 # [km]

        @property
        def NEAREARTHASTEROID1627IVAR_DIAMETER(self):
            r"""
            Diameter of near-Earth asteroid 1627 Ivar

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 16 'Prominent Aten-, Apollo-, and Amor-class near-Earth asteroids with well-determined orbits and V(1,0) < 18 mag' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [km]
            """

            return self.__NEAREARTHASTEROID1627IVAR_DIAMETER


        __NEAREARTHASTEROID1685TORO_DIAMETER: float  = 5.2 # [km]

        @property
        def NEAREARTHASTEROID1685TORO_DIAMETER(self):
            r"""
            Diameter of near-Earth asteroid 1685 Toro

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 16 'Prominent Aten-, Apollo-, and Amor-class near-Earth asteroids with well-determined orbits and V(1,0) < 18 mag' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [km]
            """

            return self.__NEAREARTHASTEROID1685TORO_DIAMETER


        __NEAREARTHASTEROID1866SISYPHUS_DIAMETER: float  = 8.2 # [km]

        @property
        def NEAREARTHASTEROID1866SISYPHUS_DIAMETER(self):
            r"""
            Diameter of near-Earth asteroid 1866 Sisyphus

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 16 'Prominent Aten-, Apollo-, and Amor-class near-Earth asteroids with well-determined orbits and V(1,0) < 18 mag' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [km]
            """

            return self.__NEAREARTHASTEROID1866SISYPHUS_DIAMETER


        __NEAREARTHASTEROID3200PHAETON_DIAMETER: float  = 6.9 # [km]

        @property
        def NEAREARTHASTEROID3200PHAETON_DIAMETER(self):
            r"""
            Diameter of near-Earth asteroid 3200 Phaeton

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 16 'Prominent Aten-, Apollo-, and Amor-class near-Earth asteroids with well-determined orbits and V(1,0) < 18 mag' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [km]
            """

            return self.__NEAREARTHASTEROID3200PHAETON_DIAMETER


        __NEAREARTHASTEROID3552DONQUIXOTE_DIAMETER: float  = 18.7 # [km]

        @property
        def NEAREARTHASTEROID3552DONQUIXOTE_DIAMETER(self):
            r"""
            Diameter of near-Earth asteroid 3552 Don Quixote

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 16 'Prominent Aten-, Apollo-, and Amor-class near-Earth asteroids with well-determined orbits and V(1,0) < 18 mag' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [km]
            """

            return self.__NEAREARTHASTEROID3552DONQUIXOTE_DIAMETER


        __NEAREARTHASTEROID433EROS_DIAMETER: float  = 22.0 # [km]

        @property
        def NEAREARTHASTEROID433EROS_DIAMETER(self):
            r"""
            Diameter of near-Earth asteroid 433 Eros

            #Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 16 'Prominent Aten-, Apollo-, and Amor-class near-Earth asteroids with well-determined orbits and V(1,0) < 18 mag' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [km]
            """

            return self.__NEAREARTHASTEROID433EROS_DIAMETER


        __NEAREARTHOBJECT_VELOCITYDISTRIBUTIONSTANDARDDEVIATION_AC: float  = 30.0 # [mas s^-1]

        @property
        def NEAREARTHOBJECT_VELOCITYDISTRIBUTIONSTANDARDDEVIATION_AC(self):
            r"""
            The velocity distribution of near-Earth objects (NEOs) is approximately Gaussian with zero mean and a standard deviation of 30.0 mas s^-1 across-scan (for a solar-aspect angle of 45 degrees)

            #Source: F. Mignard, 2002, 'Observations of solar-system objects with Gaia. I. Detection of NEOs', A&A, 393, 727, Section 4.4 (2002A&A...393..727M). See also E. Hoeg, F. Arenou, P. Hjorth, U.G. Joergensen, F. Mignard, S. Wolff, 28 February 2003, 'Faint objects and NEOs with Gaia', GAIA-CUO-118, issue 1, revision 0. Current value, for a solar-aspect angle of 45 degrees, from F. Mignard, priv. comm., 10 August 2005<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mas s^-1]
            """

            return self.__NEAREARTHOBJECT_VELOCITYDISTRIBUTIONSTANDARDDEVIATION_AC


        __NEAREARTHOBJECT_VELOCITYDISTRIBUTIONSTANDARDDEVIATION_AL: float  = 22.5 # [mas s^-1]

        @property
        def NEAREARTHOBJECT_VELOCITYDISTRIBUTIONSTANDARDDEVIATION_AL(self):
            r"""
            The velocity distribution of near-Earth objects (NEOs) is approximately Gaussian with zero mean and a standard deviation of 22.5 mas s^-1 along-scan (for a solar-aspect angle of 45 degrees)

            #Source: F. Mignard, 2002, 'Observations of solar-system objects with Gaia. I. Detection of NEOs', A&A, 393, 727, Section 4.4 (2002A&A...393..727M). See also E. Hoeg, F. Arenou, P. Hjorth, U.G. Joergensen, F. Mignard, S. Wolff, 28 February 2003, 'Faint objects and NEOs with Gaia', GAIA-CUO-118, issue 1, revision 0. Current value, for a solar-aspect angle of 45 degrees, from F. Mignard, priv. comm., 10 August 2005<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mas s^-1]
            """

            return self.__NEAREARTHOBJECT_VELOCITYDISTRIBUTIONSTANDARDDEVIATION_AL


        __NEPTUNESYSTEM_ASTROMETRICSIGNATURE_10PARSEC: float  = 155.0 # [10^-6 arcsec]

        @property
        def NEPTUNESYSTEM_ASTROMETRICSIGNATURE_10PARSEC(self):
            r"""
            Astrometric signature of the Sun induced by the Neptune system for an observer located at a distance of 10 pc from the Sun

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.7, Equation 1.22, page 11<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__NEPTUNESYSTEM_ASTROMETRICSIGNATURE_10PARSEC


        __NEPTUNESYSTEM_MASS: float  = 1.02434e+26 # [kg]

        @property
        def NEPTUNESYSTEM_MASS(self):
            r"""
            Neptune-system mass (IAU 2009 CBE value). The planetary mass includes the contribution of its satellites

            #Basic : false
            #Scalar: true
            #Unit: [kg]
            """

            return self.__NEPTUNESYSTEM_MASS


        __NEPTUNESYSTEM_ORBITALECCENTRICITY_J2000: float  = 0.00859048

        @property
        def NEPTUNESYSTEM_ORBITALECCENTRICITY_J2000(self):
            r"""
            Mean orbital eccentricity of Neptune, at the standard epoch J2000.0. The mean orbital eccentricity is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__NEPTUNESYSTEM_ORBITALECCENTRICITY_J2000


        __NEPTUNESYSTEM_ORBITALINCLINATION_J2000: float  = 1.77004347 # [deg]

        @property
        def NEPTUNESYSTEM_ORBITALINCLINATION_J2000(self):
            r"""
            Mean orbital inclination of Neptune, at the standard epoch J2000.0. The mean orbital inclination is is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__NEPTUNESYSTEM_ORBITALINCLINATION_J2000


        __NEPTUNESYSTEM_ORBITALPERIOD: float  = 164.79132 # [yr]

        @property
        def NEPTUNESYSTEM_ORBITALPERIOD(self):
            r"""
            Sidereal orbital period

            #Source: Values derived from the mean longitude rates in Table 5.8.1 in P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, page 316; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [yr]
            """

            return self.__NEPTUNESYSTEM_ORBITALPERIOD


        __NEPTUNESYSTEM_ORBITALSEMIMAJORAXIS_J2000: float  = 30.06992276 # [au]

        @property
        def NEPTUNESYSTEM_ORBITALSEMIMAJORAXIS_J2000(self):
            r"""
            Mean orbital semi-major axis of Neptune, at the standard epoch J2000.0. The mean orbital semi-major axis is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: [au]
            """

            return self.__NEPTUNESYSTEM_ORBITALSEMIMAJORAXIS_J2000


        __NEPTUNESYSTEM_RADIALVELOCITYSIGNATURE: float  = 0.3 # [m s^-1]

        @property
        def NEPTUNESYSTEM_RADIALVELOCITYSIGNATURE(self):
            r"""
            Radial-velocity amplitude of the Sun induced by the Neptune system for 'an edge-on observer' (i.e., an observer in the orbital plane of the Neptune system)

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.4, Equation 1.18, page 9<br/>
            #Basic : false
            #Scalar: true
            #Unit: [m s^-1]
            """

            return self.__NEPTUNESYSTEM_RADIALVELOCITYSIGNATURE


        __NEPTUNE_ENCOMPASSINGSPHERERADIUS: float  = 2.47640e+07 # [m]

        @property
        def NEPTUNE_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around Neptune which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__NEPTUNE_ENCOMPASSINGSPHERERADIUS


        __NEPTUNE_EQUATORIALRADIUS: float  = 2.47640e+07 # [m]

        @property
        def NEPTUNE_EQUATORIALRADIUS(self):
            r"""
            Equatorial radius of Neptune

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__NEPTUNE_EQUATORIALRADIUS


        __NEPTUNE_FLATTENING: float  = 1.710e-02

        @property
        def NEPTUNE_FLATTENING(self):
            r"""
            Geometrical flattening factor f of Neptune (f = (a-b)/a)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__NEPTUNE_FLATTENING


        __NEPTUNE_FLUXREDUCTION_MAXIMUM: float  = 0.125 # [%]

        @property
        def NEPTUNE_FLUXREDUCTION_MAXIMUM(self):
            r"""
            Maximum reduction of the solar flux for an observer external to the solar system during a transit of Neptune

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.2, Equation 2.4, page 14<br/>
            #Basic : false
            #Scalar: true
            #Unit: [%]
            """

            return self.__NEPTUNE_FLUXREDUCTION_MAXIMUM


        __NEPTUNE_GEOMETRICALBEDO: float  = 0.41

        @property
        def NEPTUNE_GEOMETRICALBEDO(self):
            r"""
            Geometric albedo of Neptune (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__NEPTUNE_GEOMETRICALBEDO


        __NEPTUNE_JSUB2: float  = 0.003538

        @property
        def NEPTUNE_JSUB2(self):
            r"""
            Dynamical form-factor of Neptune (oblateness or Stokes' second-degree zonal harmonic of the potential)

            #Source: P.R. Weissman, L.-A. McFadden, T.V. Johnson (eds.), 1999, 'Encyclopedia of the Solar System (first edition)', Academic Press, page 342<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__NEPTUNE_JSUB2


        __NEPTUNE_LIGHTDEFLECTION_LIMB: float  = 2548.0 # [10^-6 arcsec]

        @property
        def NEPTUNE_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of Neptune

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__NEPTUNE_LIGHTDEFLECTION_LIMB


        __NEPTUNE_MASS: float  = 1.024100e+26 # [kg]

        @property
        def NEPTUNE_MASS(self):
            r"""
            Mass of Neptune (do not use for high-precision (orbit) calculations)

            #Source: R.A. Jacobson, 2008, 'The orbits of the Neptunian satellites and the orientation of the pole of Neptune', BAAS, 40, 296; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [kg]
            """

            return self.__NEPTUNE_MASS


        __NEPTUNE_MASSDENSITY_MEAN: float  = 1.638 # [g cm^-3]

        @property
        def NEPTUNE_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of Neptune

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__NEPTUNE_MASSDENSITY_MEAN


        __NEPTUNE_NORTHROTATIONALPOLE_DECLINATION: float  = 42.950 # [deg]

        @property
        def NEPTUNE_NORTHROTATIONALPOLE_DECLINATION(self):
            r"""
            IAU-recommended value for the declination \delta_0 of the north pole of rotation of Neptune. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value is based on an approximate formula, accurate to first order in 'd' (see below); the true equation is \delta_0 = 43.46 - 0.51 * cos(357.85 + 52.316 * d / 36525), where d is the number of Julian days calculated from the standard epoch. The numerical accuracy of this equation is 0.1 deg

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : false
            #Scalar: true
            #Unit: [deg]
            """

            return self.__NEPTUNE_NORTHROTATIONALPOLE_DECLINATION


        __NEPTUNE_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE: float  = -0.0000004783 # [deg day^-1]

        @property
        def NEPTUNE_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the declination \delta_0 of the north pole of rotation of Neptune. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value is based on an approximate formula, accurate to first order in 'd' (see below); the true equation is \delta_0 = 43.46 - 0.51 * cos(357.85 + 52.316 * d / 36525), where d is the number of Julian days calculated from the standard epoch. The numerical accuracy of this equation is 0.1 deg

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : false
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__NEPTUNE_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE


        __NEPTUNE_NORTHROTATIONALPOLE_RIGHTASCENSION: float  = 299.334 # [deg]

        @property
        def NEPTUNE_NORTHROTATIONALPOLE_RIGHTASCENSION(self):
            r"""
            IAU-recommended value for the right ascension \alpha_0 of the north pole of rotation of Neptune. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value is based on an approximate formula, accurate to first order in 'd' (see below); the true equation is \alpha_0 = 299.36 + 0.70 * sin(357.85 + 52.316 * d / 36525), where d is the number of Julian days calculated from the standard epoch. The numerical accuracy of this equation is 0.1 deg

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : false
            #Scalar: true
            #Unit: [deg]
            """

            return self.__NEPTUNE_NORTHROTATIONALPOLE_RIGHTASCENSION


        __NEPTUNE_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE: float  = 0.0000174869 # [deg day^-1]

        @property
        def NEPTUNE_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the right ascension \alpha_0 of the north pole of rotation of Neptune. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value is based on an approximate formula, accurate to first order in 'd' (see below); the true equation is \alpha_0 = 299.36 + 0.70 * sin(357.85 + 52.316 * d / 36525), where d is the number of Julian days calculated from the standard epoch. The numerical accuracy of this equation is 0.1 deg

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : false
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__NEPTUNE_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE


        __NEPTUNE_POLARRADIUS: float  = 2.43405e+07 # [m]

        @property
        def NEPTUNE_POLARRADIUS(self):
            r"""
            Polar radius of Neptune

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__NEPTUNE_POLARRADIUS


        __NEPTUNE_PRIMEMERIDIAN_EPHEMERISPOSITION: float  = 253.198 # [deg]

        @property
        def NEPTUNE_PRIMEMERIDIAN_EPHEMERISPOSITION(self):
            r"""
            IAU-recommended value for the ephemeris position of the prime meridian of Neptune. The prime meridian refers to the rotation of the magnetic field System III. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde. The numerical value is based on an approximate formula, accurate to first order in 'd' (see below); the true equation is W = 253.18 + 536.3128492 * d - 0.48 * sin(357.85 + 52.316 * d / 36525), where d is the number of Julian days calculated from the standard epoch. The numerical accuracy of this equation is 0.1 deg

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : false
            #Scalar: true
            #Unit: [deg]
            """

            return self.__NEPTUNE_PRIMEMERIDIAN_EPHEMERISPOSITION


        __NEPTUNE_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE: float  = 536.3128372090 # [deg day^-1]

        @property
        def NEPTUNE_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch of 1.5 January 2000 = JD2451545.0 TDB) of the ephemeris position of the prime meridian of Neptune. The prime meridian refers to the rotation of the magnetic field System III. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde. The numerical value is based on an approximate formula, accurate to first order in 'd' (see below); the true equation is W = 253.18 + 536.3128492 * d - 0.48 * sin(357.85 + 52.316 * d / 36525), where d is the number of Julian days calculated from the standard epoch. The numerical accuracy of this equation is 0.1 deg

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : false
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__NEPTUNE_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE


        __NEPTUNE_TRANSITPROBABILITY: float  = 0.016 # [%]

        @property
        def NEPTUNE_TRANSITPROBABILITY(self):
            r"""
            Geometric transit probability (Neptune transiting the Sun) for an observer external to the solar system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.1, Equation 2.2, page 14<br/>
            #Basic : false
            #Scalar: true
            #Unit: [%]
            """

            return self.__NEPTUNE_TRANSITPROBABILITY


        __NEPTUNE_TRANSITTIME_MAXIMUM: float  = 3.07 # [day]

        @property
        def NEPTUNE_TRANSITTIME_MAXIMUM(self):
            r"""
            Maximum transit time of Neptune (transiting the Sun) for an observer external to the solar system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.3, Equation 2.5, page 15<br/>
            #Basic : false
            #Scalar: true
            #Unit: [day]
            """

            return self.__NEPTUNE_TRANSITTIME_MAXIMUM


        __NEPTUNE_VONEZEROMAGNITUDE: float  = -6.87 # [mag]

        @property
        def NEPTUNE_VONEZEROMAGNITUDE(self):
            r"""
            V(1,0) magnitude of Neptune (i.e., the visual magnitude of the planet reduced to a distance of 1 au from both the Sun and Earth and phase angle zero). This parameter is also refered to as absolute magnitude in planetary sciences

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__NEPTUNE_VONEZEROMAGNITUDE


        __NEPTUNE_VOLUMETRICRADIUS: float  = 2.46220e+07 # [m]

        @property
        def NEPTUNE_VOLUMETRICRADIUS(self):
            r"""
            Mean volumetric radius of Neptune

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__NEPTUNE_VOLUMETRICRADIUS


        __NEWTON_CONSTANT: float  = 6.674080e-11 # [m^3 kg^-1 s^-2]

        @property
        def NEWTON_CONSTANT(self):
            r"""
            Newton's universal constant of gravitation

            #Source: P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0). See also the IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^3 kg^-1 s^-2]
            """

            return self.__NEWTON_CONSTANT


        __NOMINALSUN_MEANLONGITUDERATE_J2000: float  = 0.98560903 # [deg day^-1]

        @property
        def NOMINALSUN_MEANLONGITUDERATE_J2000(self):
            r"""
            Mean (geometric) longitude rate of the nominal Sun for use in simulations of the NSL (mean ecliptic orbital elements, at the standard epoch J2000.0). Note that a value of 1295977422.83429 / (1.0E3 * 365.25 * 3600.0) = 0.98560911 degrees day^-1 is given in Section 5.8.3 of J.L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touze, G. Francou, J. Laskar, 1994, 'Numerical expressions for precession formulae and mean elements for the Moon and the planets', A&A, 282, 663 (1994A&A...282..663S)

            #Source: F. Mignard, priv. comm., 14 July 2004, based on the Fortran-90 subroutine SCANNING (version 4.1, April 2004); see also GAIA-FM-010 and GAIA-FM-017<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__NOMINALSUN_MEANLONGITUDERATE_J2000


        __NOMINALSUN_MEANLONGITUDE_J2000: float  = 280.4665 # [deg]

        @property
        def NOMINALSUN_MEANLONGITUDE_J2000(self):
            r"""
            Mean (geometric) longitude of the nominal Sun for use in simulations of the NSL (mean ecliptic orbital elements, at the standard epoch J2000.0); subtract aberration (about 20 arcsec; see parameter :Nature:Aberration_Constant_J2000) to obtain the apparent longitude of the nominal Sun. Note that a value of 280.46645683 degrees is given in Section 5.8.3 of J.L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touze, G. Francou, J. Laskar, 1994, 'Numerical expressions for precession formulae and mean elements for the Moon and the planets', A&A, 282, 663 (1994A&A...282..663S)

            #Source: F. Mignard, priv. comm., 14 July 2004, based on the Fortran-90 subroutine SCANNING (version 4.1, April 2004); see also GAIA-FM-010 and GAIA-FM-017<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__NOMINALSUN_MEANLONGITUDE_J2000


        __NOMINALSUN_ORBITALECCENTRICITY_J2000: float  = 0.01671

        @property
        def NOMINALSUN_ORBITALECCENTRICITY_J2000(self):
            r"""
            Mean orbital eccentricity of the nominal Sun for use in simulations of the NSL (mean ecliptic orbital elements, at the standard epoch J2000.0). See also the parameter :Nature:EMBC_OrbitalEccentricity_J2000. Note that a value of 0.0167086342 is given in Section 5.8.3 of J.L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touze, G. Francou, J. Laskar, 1994, 'Numerical expressions for precession formulae and mean elements for the Moon and the planets', A&A, 282, 663 (1994A&A...282..663S)

            #Source: F. Mignard, priv. comm., 14 July 2004, based on the Fortran-90 subroutine SCANNING (version 4.1, April 2004); see also GAIA-FM-010 and GAIA-FM-017<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__NOMINALSUN_ORBITALECCENTRICITY_J2000


        __NOMINALSUN_ORBITALMEANANOMALYRATE_J2000: float  = 0.98560020 # [deg day^-1]

        @property
        def NOMINALSUN_ORBITALMEANANOMALYRATE_J2000(self):
            r"""
            Orbital mean anomaly rate of the nominal Sun for use in simulations of the NSL (mean ecliptic orbital elements, at the standard epoch J2000.0). Note that a value of (1295977422.83429 - 11612.35290) / (1.0E3 * 365.25 * 3600.0) = 0.98560028 degrees day^-1 is given in Section 5.8.3 of J.L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touze, G. Francou, J. Laskar, 1994, 'Numerical expressions for precession formulae and mean elements for the Moon and the planets', A&A, 282, 663 (1994A&A...282..663S)

            #Source: F. Mignard, priv. comm., 14 July 2004, based on the Fortran-90 subroutine SCANNING (version 4.1, April 2004); see also GAIA-FM-010 and GAIA-FM-017<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__NOMINALSUN_ORBITALMEANANOMALYRATE_J2000


        __NOMINALSUN_ORBITALMEANANOMALY_J2000: float  = 357.529 # [deg]

        @property
        def NOMINALSUN_ORBITALMEANANOMALY_J2000(self):
            r"""
            Orbital mean anomaly of the nominal Sun for use in simulations of the NSL (mean ecliptic orbital elements, at the standard epoch J2000.0). Note that a value of 100.46645683 - 102.93734808 + 360 = 357.52910875 degrees is given in Section 5.8.3 of J.L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touze, G. Francou, J. Laskar, 1994, 'Numerical expressions for precession formulae and mean elements for the Moon and the planets', A&A, 282, 663 (1994A&A...282..663S)

            #Source: F. Mignard, priv. comm., 14 July 2004, based on the Fortran-90 subroutine SCANNING (version 4.1, April 2004); see also GAIA-FM-010 and GAIA-FM-017<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__NOMINALSUN_ORBITALMEANANOMALY_J2000


        __NUTATION_CONSTANT_J2000: float  = 9.2025 # [arcsec]

        @property
        def NUTATION_CONSTANT_J2000(self):
            r"""
            Constant of nutation, at the standard epoch J2000.0, nowadays irrelevant as a fundamental constant

            #Source: P.K. Seidelmann, May 1982, '1980 IAU theory of nutation. The final report of the IAU Working Group on Nutation', Celestial Mechanics, 27, pages 79-106 (1982CeMec..27...79S)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [arcsec]
            """

            return self.__NUTATION_CONSTANT_J2000


        __OBLIQUITYOFECLIPTIC_J2000: float  = 84381.41100 # [arcsec]

        @property
        def OBLIQUITYOFECLIPTIC_J2000(self):
            r"""
            Obliquity of the ecliptic with respect to the ICRS reference plane, at the standard epoch J2000.0. Note that the ICRS origin is shifted in the equatorial plane from \Gamma by \phi = 0.05542 arcsec, positive from \Gamma to the ICRS origin. Note that the value of the obliquity of the ecliptic in the inertial sense, i.e., with respect to the CIP equator (Celestial Intermediate Pole, formerly the Mean Celestial Ephemeris Pole or MCEP) equals 84381.406 arcsec (from the IAU (2009) System of Astronomical Constants: IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)

            #Source: J. Chapront, M. Chapront-Touze, G. Francou, 2002, 'A new determination of lunar orbital parameters, precession constant, and tidal acceleration from LLR measurements', A&A, 387, 700<br/>
            #Basic : true
            #Scalar: true
            #Unit: [arcsec]
            """

            return self.__OBLIQUITYOFECLIPTIC_J2000


        __OORT_CONSTANT_A: float  = 14.82 # [km s^-1 kpc^-1]

        @property
        def OORT_CONSTANT_A(self):
            r"""
            Oort constant A

            #Source: M. Feast, P. Whitelock, 1997, MNRAS, 291, 683<br/>
            #Basic : true
            #Scalar: true
            #Unit: [km s^-1 kpc^-1]
            """

            return self.__OORT_CONSTANT_A


        __OORT_CONSTANT_B: float  = -12.37 # [km s^-1 kpc^-1]

        @property
        def OORT_CONSTANT_B(self):
            r"""
            Oort constant B

            #Source: M. Feast, P. Whitelock, 1997, MNRAS, 291, 683<br/>
            #Basic : true
            #Scalar: true
            #Unit: [km s^-1 kpc^-1]
            """

            return self.__OORT_CONSTANT_B


        __PPN_ALPHA_1: float  = 0.0

        @property
        def PPN_ALPHA_1(self):
            r"""
            General relativistic standard PPN parameter \alpha_1, quantifying prefered-frame effects

            #Source: E.g., C.M. Will, 2006, 'The Confrontation between General Relativity and Experiment', Living Reviews in Relativity, Section 3.2 'The parametrized post-Newtonian formalism', http://www.livingreviews.org/lrr-2006-3<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__PPN_ALPHA_1


        __PPN_ALPHA_2: float  = 0.0

        @property
        def PPN_ALPHA_2(self):
            r"""
            General relativistic standard PPN parameter \alpha_2, quantifying prefered-frame effects

            #Source: E.g., C.M. Will, 2006, 'The Confrontation between General Relativity and Experiment', Living Reviews in Relativity, Section 3.2 'The parametrized post-Newtonian formalism', http://www.livingreviews.org/lrr-2006-3<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__PPN_ALPHA_2


        __PPN_ALPHA_3: float  = 0.0

        @property
        def PPN_ALPHA_3(self):
            r"""
            General relativistic standard PPN parameter \alpha_3, quantifying prefered-frame effects

            #Source: E.g., C.M. Will, 2006, 'The Confrontation between General Relativity and Experiment', Living Reviews in Relativity, Section 3.2 'The parametrized post-Newtonian formalism', http://www.livingreviews.org/lrr-2006-3<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__PPN_ALPHA_3


        __PPN_BETA: float  = 1.0

        @property
        def PPN_BETA(self):
            r"""
            General relativistic standard PPN parameter \beta, quantifying non-linearity in the gravitational superposition law

            #Source: E.g., C.M. Will, 2006, 'The Confrontation between General Relativity and Experiment', Living Reviews in Relativity, Section 3.2 'The parametrized post-Newtonian formalism', http://www.livingreviews.org/lrr-2006-3<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__PPN_BETA


        __PPN_GAMMA: float  = 1.0

        @property
        def PPN_GAMMA(self):
            r"""
            General relativistic standard PPN parameter \gamma, quantifying space curvature per unit rest mass

            #Source: E.g., C.M. Will, 2006, 'The Confrontation between General Relativity and Experiment', Living Reviews in Relativity, Section 3.2 'The parametrized post-Newtonian formalism', http://www.livingreviews.org/lrr-2006-3<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__PPN_GAMMA


        __PPN_XI: float  = 0.0

        @property
        def PPN_XI(self):
            r"""
            General relativistic standard PPN parameter \xi, quantifying prefered-location effects

            #Source: E.g., C.M. Will, 2006, 'The Confrontation between General Relativity and Experiment', Living Reviews in Relativity, Section 3.2 'The parametrized post-Newtonian formalism', http://www.livingreviews.org/lrr-2006-3<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__PPN_XI


        __PPN_ZETA_1: float  = 0.0

        @property
        def PPN_ZETA_1(self):
            r"""
            General relativistic standard PPN parameter \zeta_1, quantifying violation of conservation of total momentum

            #Source: E.g., C.M. Will, 2006, 'The Confrontation between General Relativity and Experiment', Living Reviews in Relativity, Section 3.2 'The parametrized post-Newtonian formalism', http://www.livingreviews.org/lrr-2006-3<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__PPN_ZETA_1


        __PPN_ZETA_2: float  = 0.0

        @property
        def PPN_ZETA_2(self):
            r"""
            General relativistic standard PPN parameter \zeta_2, quantifying violation of conservation of total momentum

            #Source: E.g., C.M. Will, 2006, 'The Confrontation between General Relativity and Experiment', Living Reviews in Relativity, Section 3.2 'The parametrized post-Newtonian formalism', http://www.livingreviews.org/lrr-2006-3<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__PPN_ZETA_2


        __PPN_ZETA_3: float  = 0.0

        @property
        def PPN_ZETA_3(self):
            r"""
            General relativistic standard PPN parameter \zeta_3, quantifying violation of conservation of total momentum

            #Source: E.g., C.M. Will, 2006, 'The Confrontation between General Relativity and Experiment', Living Reviews in Relativity, Section 3.2 'The parametrized post-Newtonian formalism', http://www.livingreviews.org/lrr-2006-3<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__PPN_ZETA_3


        __PPN_ZETA_4: float  = 0.0

        @property
        def PPN_ZETA_4(self):
            r"""
            General relativistic standard PPN parameter \zeta_4, quantifying violation of conservation of total momentum

            #Source: E.g., C.M. Will, 2006, 'The Confrontation between General Relativity and Experiment', Living Reviews in Relativity, Section 3.2 'The parametrized post-Newtonian formalism', http://www.livingreviews.org/lrr-2006-3<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__PPN_ZETA_4


        __PARSEC_ASTRONOMICALUNIT: float  = 206264.806247096 # [au]

        @property
        def PARSEC_ASTRONOMICALUNIT(self):
            r"""
            Parsec expressed in au

            #Basic : false
            #Scalar: true
            #Unit: [au]
            """

            return self.__PARSEC_ASTRONOMICALUNIT


        __PARSEC_METER: float  = 3.0856775814913674e+16 # [m]

        @property
        def PARSEC_METER(self):
            r"""
            Parsec expressed in m

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__PARSEC_METER


        __PI_CONSTANT: float  = 3.14159265358979323846264338328

        @property
        def PI_CONSTANT(self):
            r"""
            The constant Pi (also known as Archimedes' constant). Note that double-precision, floating-point numbers in any programming language which follows the IEEE standard (true for C, C++, Java, and most, if not all, others) have only 16 significant digits (64 bits); the representation here, using 30 significant digits, is thus amply sufficient

            #Source: Well-known mathematical constant; numerical value can be extracted, e.g., from Mathematica 4.0 for Solaris (Wolfram Research, Inc.) using 'N[Pi,30]'<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__PI_CONSTANT


        __PLANCK_CONSTANT: float  = 6.6260700400e-34 # [J s]

        @property
        def PLANCK_CONSTANT(self):
            r"""
            Planck's constant

            #Source: P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [J s]
            """

            return self.__PLANCK_CONSTANT


        __PLUTOSYSTEM_ASTROMETRICSIGNATURE_10PARSEC: float  = 0.029 # [10^-6 arcsec]

        @property
        def PLUTOSYSTEM_ASTROMETRICSIGNATURE_10PARSEC(self):
            r"""
            Astrometric signature of the Sun induced by the Pluto system for an observer located at a distance of 10 pc from the Sun

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.7, Equation 1.22, page 11<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__PLUTOSYSTEM_ASTROMETRICSIGNATURE_10PARSEC


        __PLUTOSYSTEM_MASS: float  = 1.4561e+22 # [kg]

        @property
        def PLUTOSYSTEM_MASS(self):
            r"""
            Pluto-system mass (IAU 2009 CBE value). The 'planetary' mass includes the contribution of its satellite, Charon

            #Basic : false
            #Scalar: true
            #Unit: [kg]
            """

            return self.__PLUTOSYSTEM_MASS


        __PLUTOSYSTEM_ORBITALECCENTRICITY_J2000: float  = 0.24882730

        @property
        def PLUTOSYSTEM_ORBITALECCENTRICITY_J2000(self):
            r"""
            Mean orbital eccentricity of Pluto, at the standard epoch J2000.0. The mean orbital eccentricity is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__PLUTOSYSTEM_ORBITALECCENTRICITY_J2000


        __PLUTOSYSTEM_ORBITALINCLINATION_J2000: float  = 17.14001206 # [deg]

        @property
        def PLUTOSYSTEM_ORBITALINCLINATION_J2000(self):
            r"""
            Mean orbital inclination of Pluto, at the standard epoch J2000.0. The mean orbital inclination is is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__PLUTOSYSTEM_ORBITALINCLINATION_J2000


        __PLUTOSYSTEM_ORBITALPERIOD: float  = 247.92065 # [yr]

        @property
        def PLUTOSYSTEM_ORBITALPERIOD(self):
            r"""
            Sidereal orbital period

            #Source: Values derived from the mean longitude rates in Table 5.8.1 in P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, page 316; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [yr]
            """

            return self.__PLUTOSYSTEM_ORBITALPERIOD


        __PLUTOSYSTEM_ORBITALSEMIMAJORAXIS_J2000: float  = 39.48211675 # [au]

        @property
        def PLUTOSYSTEM_ORBITALSEMIMAJORAXIS_J2000(self):
            r"""
            Mean orbital semi-major axis of Pluto, at the standard epoch J2000.0. The mean orbital semi-major axis is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: [au]
            """

            return self.__PLUTOSYSTEM_ORBITALSEMIMAJORAXIS_J2000


        __PLUTOSYSTEM_RADIALVELOCITYSIGNATURE: float  = 3.58e-05 # [m s^-1]

        @property
        def PLUTOSYSTEM_RADIALVELOCITYSIGNATURE(self):
            r"""
            Radial-velocity amplitude of the Sun induced by the Pluto system for 'an edge-on observer' (i.e., an observer in the orbital plane of the Pluto system)

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.4, Equation 1.18, page 9<br/>
            #Basic : false
            #Scalar: true
            #Unit: [m s^-1]
            """

            return self.__PLUTOSYSTEM_RADIALVELOCITYSIGNATURE


        __PLUTO_ENCOMPASSINGSPHERERADIUS: float  = 1.1950e+06 # [m]

        @property
        def PLUTO_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around Pluto which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__PLUTO_ENCOMPASSINGSPHERERADIUS


        __PLUTO_EQUATORIALRADIUS: float  = 1.1950e+06 # [m]

        @property
        def PLUTO_EQUATORIALRADIUS(self):
            r"""
            Equatorial radius of Pluto

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__PLUTO_EQUATORIALRADIUS


        __PLUTO_FLATTENING: float  = 0.0

        @property
        def PLUTO_FLATTENING(self):
            r"""
            Geometrical flattening factor f of Pluto (f = (a-b)/a)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__PLUTO_FLATTENING


        __PLUTO_FLUXREDUCTION_MAXIMUM: float  = 0.0003 # [%]

        @property
        def PLUTO_FLUXREDUCTION_MAXIMUM(self):
            r"""
            Maximum reduction of the solar flux for an observer external to the solar system during a transit of Pluto

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.2, Equation 2.4, page 14<br/>
            #Basic : false
            #Scalar: true
            #Unit: [%]
            """

            return self.__PLUTO_FLUXREDUCTION_MAXIMUM


        __PLUTO_GEOMETRICALBEDO: float  = 0.3

        @property
        def PLUTO_GEOMETRICALBEDO(self):
            r"""
            Geometric albedo of Pluto (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__PLUTO_GEOMETRICALBEDO


        __PLUTO_JSUB2: float  = 0.0

        @property
        def PLUTO_JSUB2(self):
            r"""
            Dynamical form-factor of Pluto (oblateness or Stokes' second-degree zonal harmonic of the potential); value is unknown

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__PLUTO_JSUB2


        __PLUTO_LIGHTDEFLECTION_LIMB: float  = 7.0 # [10^-6 arcsec]

        @property
        def PLUTO_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of Pluto

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__PLUTO_LIGHTDEFLECTION_LIMB


        __PLUTO_MASS: float  = 1.3090e+22 # [kg]

        @property
        def PLUTO_MASS(self):
            r"""
            Mass of Pluto (do not use for high-precision (orbit) calculations)

            #Source: R.A. Jacobson, 2007, 'The orbits of the satellites of Pluto - Ephemeris PLU017', priv. comm.; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [kg]
            """

            return self.__PLUTO_MASS


        __PLUTO_MASSDENSITY_MEAN: float  = 1.83 # [g cm^-3]

        @property
        def PLUTO_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of Pluto (rough estimate)

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__PLUTO_MASSDENSITY_MEAN


        __PLUTO_POLARRADIUS: float  = 1.1950e+06 # [m]

        @property
        def PLUTO_POLARRADIUS(self):
            r"""
            Polar radius of Pluto

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__PLUTO_POLARRADIUS


        __PLUTO_TRANSITPROBABILITY: float  = 0.012 # [%]

        @property
        def PLUTO_TRANSITPROBABILITY(self):
            r"""
            Geometric transit probability (Pluto transiting the Sun) for an observer external to the solar system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.1, Equation 2.2, page 14<br/>
            #Basic : false
            #Scalar: true
            #Unit: [%]
            """

            return self.__PLUTO_TRANSITPROBABILITY


        __PLUTO_TRANSITTIME_MAXIMUM: float  = 3.40 # [day]

        @property
        def PLUTO_TRANSITTIME_MAXIMUM(self):
            r"""
            Maximum transit time of Pluto (transiting the Sun) for an observer external to the solar system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.3, Equation 2.5, page 15<br/>
            #Basic : false
            #Scalar: true
            #Unit: [day]
            """

            return self.__PLUTO_TRANSITTIME_MAXIMUM


        __PLUTO_VONEZEROMAGNITUDE: float  = -1.0 # [mag]

        @property
        def PLUTO_VONEZEROMAGNITUDE(self):
            r"""
            V(1,0) magnitude of Pluto (i.e., the visual magnitude of the 'planet' reduced to a distance of 1 au from both the Sun and Earth and phase angle zero). This parameter is also refered to as absolute magnitude in planetary sciences

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__PLUTO_VONEZEROMAGNITUDE


        __PLUTO_VOLUMETRICRADIUS: float  = 1.1950e+06 # [m]

        @property
        def PLUTO_VOLUMETRICRADIUS(self):
            r"""
            Mean volumetric radius of Pluto

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__PLUTO_VOLUMETRICRADIUS


        __PRECESSIONLONGITUDE_CONSTANT_J2000: float  = 5028.796195 # [arcsec cy^-1]

        @property
        def PRECESSIONLONGITUDE_CONSTANT_J2000(self):
            r"""
            Speed of general precession in ecliptic longitude, in arcsec per Julian century, at the standard epoch J2000.0, nowadays irrelevant as a fundamental constant

            #Source: N. Capitaine, P.T. Wallace, J. Chapront, 2003, 'Expressions for IAU 2000 precession quantities', A&A, 412, 567-586, 'P03 solution'<br/>
            #Basic : true
            #Scalar: true
            #Unit: [arcsec cy^-1]
            """

            return self.__PRECESSIONLONGITUDE_CONSTANT_J2000


        __PRECESSION_CONSTANT_J2000M: float  = 3.075887 # [s yr^-1]

        @property
        def PRECESSION_CONSTANT_J2000M(self):
            r"""
            Precession constant m = p cos(\epsilon_0), in s per Julian year, at the standard epoch J2000.0, nowadays irrelevant as a fundamental constant. The precession in right ascension \alpha equals m + n sin(\alpha) tan(\delta); the precession in declination \delta equals n cos(\alpha)

            #Basic : false
            #Scalar: true
            #Unit: [s yr^-1]
            """

            return self.__PRECESSION_CONSTANT_J2000M


        __PRECESSION_CONSTANT_J2000N: float  = 20.003394 # [arcsec yr^-1]

        @property
        def PRECESSION_CONSTANT_J2000N(self):
            r"""
            Precession constant n = p sin(\epsilon_0), in arcsec per Julian year, at the standard epoch J2000.0, nowadays irrelevant as a fundamental constant. The precession in right ascension \alpha equals m + n sin(\alpha) tan(\delta); the precession in declination \delta equals n cos(\alpha)

            #Basic : false
            #Scalar: true
            #Unit: [arcsec yr^-1]
            """

            return self.__PRECESSION_CONSTANT_J2000N


        __PROPERMOTION_CONSTANT: float  = 4.740470464 # [km yr s^-1]

        @property
        def PROPERMOTION_CONSTANT(self):
            r"""
            Proper-motion conversion constant A_v = 4.74... km yr s^-1 (see ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, page 25, Table 1.2.2)

            #Basic : false
            #Scalar: true
            #Unit: [km yr s^-1]
            """

            return self.__PROPERMOTION_CONSTANT


        __QSO_SPECTRUM_LAMBDA: str  = "Nature/QSO_Spectrum_Lambda_001.fits"

        @property
        def QSO_SPECTRUM_LAMBDA(self):
            r"""
            Composite (unweighted-mean) zero-redshift quasar (QSO) spectrum, based on 949 LBQS QSOs, 191 2MASS AGN, 37 Hamburg-ESO QSOs, and 18 PG QSOs. First column: wavelength \lambda (in nm; from 96.10 to 932.85). Second column: unnormalised flux density f_\lambda (in W m^-2 nm^-1)

            #Source: P.J. Francis, et al., 1991, 'A high signal-to-noise ratio composite quasar spectrum', Astrophysical Journal (ApJ), 373, 465<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__QSO_SPECTRUM_LAMBDA


        __RADIAN_ARCSECOND: float  = 2.062648062470964e+05 # [arcsec]

        @property
        def RADIAN_ARCSECOND(self):
            r"""
            One radian in units of arcseconds

            #Basic : false
            #Scalar: true
            #Unit: [arcsec]
            """

            return self.__RADIAN_ARCSECOND


        __RADIAN_DEGREE: float  = 5.729577951308232e+01 # [deg]

        @property
        def RADIAN_DEGREE(self):
            r"""
            One radian in units of degrees

            #Basic : false
            #Scalar: true
            #Unit: [deg]
            """

            return self.__RADIAN_DEGREE


        __RADIAN_MICROARCSECOND: float  = 2.062648062470964e+11 # [micro-arcsec]

        @property
        def RADIAN_MICROARCSECOND(self):
            r"""
            One radian in units of micro-arcseconds

            #Basic : false
            #Scalar: true
            #Unit: [micro-arcsec]
            """

            return self.__RADIAN_MICROARCSECOND


        __RADIAN_MILLIARCSECOND: float  = 2.062648062470964e+08 # [milli-arcsec]

        @property
        def RADIAN_MILLIARCSECOND(self):
            r"""
            One radian in units of milli-arcseconds

            #Basic : false
            #Scalar: true
            #Unit: [milli-arcsec]
            """

            return self.__RADIAN_MILLIARCSECOND


        __RADIATION_CONSTANT: float  = 7.5657229e-16 # [J m^-3 K^-4]

        @property
        def RADIATION_CONSTANT(self):
            r"""
            Radiation constant, also known as radiation-density constant, linking the energy density u (= 4 Pi I / c) of black-body radiation and temperature T via u = a T^4

            #Source: E.g., H. Karttunen, et al., 1987, 'Fundamental Astronomy', Springer Verlag, Berlin, Section 11.2, page 247 or R. Kippenhahn, A. Weigert, 1991, 'Stellar structure and evolution' (corrected 2-nd printing), Springer Verlag, Berlin, Section 3.1, page 16, and Section 5.1.2, page 28<br/>
            #Basic : false
            #Scalar: true
            #Unit: [J m^-3 K^-4]
            """

            return self.__RADIATION_CONSTANT


        __RADIATION_CONSTANT_FIRST: float  = 3.7417717901e-16 # [W m^2]

        @property
        def RADIATION_CONSTANT_FIRST(self):
            r"""
            First radiation constant. Note: best-measured value equals 3.741771790E-16 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))

            #Basic : false
            #Scalar: true
            #Unit: [W m^2]
            """

            return self.__RADIATION_CONSTANT_FIRST


        __RADIATION_CONSTANT_SECOND: float  = 1.438777364e-02 # [m K]

        @property
        def RADIATION_CONSTANT_SECOND(self):
            r"""
            Second radiation constant. Note: best-measured value equals 1.43877736E-2 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))

            #Basic : false
            #Scalar: true
            #Unit: [m K]
            """

            return self.__RADIATION_CONSTANT_SECOND


        __REFERENCEEPOCH_TCB = "JD2443144.5003725 TCB"

        @property
        def REFERENCEEPOCH_TCB(self):
            r"""
            The origin of TCB is defined in terms of TAI: the reading of TCB on 1 January 1977, 00:00:00 TAI = JD2443144.5 TAI must be 1 January 1977, 00:00:32.184 TCB = JD2443144.5003725 TCB. This origin has been arbitrarily set so that TCB coincides with TT at the geocentre on 1 January 1977, 00:00:00 TAI

            #Source: IAU, July 1991, 'Recommendations from the Working Group on Reference Systems', IAU 1991 Resolution A4, Recommendation 4, adopted at the XXI-st General Assembly of the IAU. See, for example, M. Soffel, et al., 1 December 2003, 'The IAU 2000 resolutions for astrometry, celestial mechanics, and metrology in the relativistic framework: explanatory supplement', AJ, 126, 2687-2706<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__REFERENCEEPOCH_TCB


        __REFERENCEEPOCH_TCG = "JD2443144.5003725 TCG"

        @property
        def REFERENCEEPOCH_TCG(self):
            r"""
            The origin of TCG is defined in terms of TAI: the reading of TCG on 1 January 1977, 00:00:00 TAI = JD2443144.5 TAI must be 1 January 1977, 00:00:32.184 TCG = JD2443144.5003725 TCG. This origin has been arbitrarily set so that TCG coincides with TT at the geocentre on 1 January 1977, 00:00:00 TAI

            #Source: IAU, July 1991, 'Recommendations from the Working Group on Reference Systems', IAU 1991 Resolution A4, Recommendation 4, adopted at the XXI-st General Assembly of the IAU. See, for example, M. Soffel, et al., 1 December 2003, 'The IAU 2000 resolutions for astrometry, celestial mechanics, and metrology in the relativistic framework: explanatory supplement', AJ, 126, 2687-2706<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__REFERENCEEPOCH_TCG


        __REFERENCEEPOCH_TDBSUBZERO: float  = -6.550e-05 # [s]

        @property
        def REFERENCEEPOCH_TDBSUBZERO(self):
            r"""
            IAU 2006 Resolution B3, entitled 'Re-definition of Barycentric Dynamical Time, TDB', recommends that 'TDB be defined as the following linear transformation of TCB: TDB = TCB - L_B x ( JD_TCB - T_0 ) x 86400 + TDB_0, where T_0 = 2443144.5003725 (parameter :Nature:ReferenceEpoch_TCB), and L_B = 1.550519768E-8 (parameter :Nature:LSubB_Constant) and TDB_0 = -6.55E-5 s are defining constants'. The number 86400 is equal to parameter :Nature:Day_Second

            #Source: IAU, August 2006, 'Re-definition of Barycentric Dynamical Time, TDB', IAU 2006 Resolution 3 adopted at the XXVI-th General Assembly of the IAU. See also the IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [s]
            """

            return self.__REFERENCEEPOCH_TDBSUBZERO


        __RYDBERG_CONSTANT: float  = 10973731.570551 # [m^-1]

        @property
        def RYDBERG_CONSTANT(self):
            r"""
            Rydberg constant. Note: best-measured value equals 10973731.568508 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))

            #Basic : false
            #Scalar: true
            #Unit: [m^-1]
            """

            return self.__RYDBERG_CONSTANT


        __SSBC_VELOCITYCMB = [ -26.29,  -244.96,  275.93 ] # [km s^-1]

        @property
        def SSBC_VELOCITYCMB(self):
            r"""
            Velocity vector of the Solar-System BaryCentre (SSBC) with respect to the Cosmic Microwave Background (CMB), in units of km s^-1. The vector elements refer to Galactic (U,V,W) coordinates

            #Source: G. Hinshaw, et al., 11 February 2009, 'Five-Year Wilkinson Microwave Anisotropy Probe (WMAP) Observations: Data Processing, Sky Maps, and Basic Results', Astrophysical Journal Supplement (ApJS), Volume 180, pages 225-245<br/>
            #Basic : true
            #Scalar: true
            #Unit: [km s^-1]
            """

            return self.__SSBC_VELOCITYCMB


        __SATURNSYSTEM_ASTROMETRICSIGNATURE_10PARSEC: float  = 273.0 # [10^-6 arcsec]

        @property
        def SATURNSYSTEM_ASTROMETRICSIGNATURE_10PARSEC(self):
            r"""
            Astrometric signature of the Sun induced by the Saturn system for an observer located at a distance of 10 pc from the Sun

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.7, Equation 1.22, page 11<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__SATURNSYSTEM_ASTROMETRICSIGNATURE_10PARSEC


        __SATURNSYSTEM_MASS: float  = 5.68477e+26 # [kg]

        @property
        def SATURNSYSTEM_MASS(self):
            r"""
            Saturn-system mass (IAU 2009 CBE value). The planetary mass includes the contribution of its satellites

            #Basic : false
            #Scalar: true
            #Unit: [kg]
            """

            return self.__SATURNSYSTEM_MASS


        __SATURNSYSTEM_ORBITALECCENTRICITY_J2000: float  = 0.05386179

        @property
        def SATURNSYSTEM_ORBITALECCENTRICITY_J2000(self):
            r"""
            Mean orbital eccentricity of Saturn, at the standard epoch J2000.0. The mean orbital eccentricity is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__SATURNSYSTEM_ORBITALECCENTRICITY_J2000


        __SATURNSYSTEM_ORBITALINCLINATION_J2000: float  = 2.48599187 # [deg]

        @property
        def SATURNSYSTEM_ORBITALINCLINATION_J2000(self):
            r"""
            Mean orbital inclination of Saturn, at the standard epoch J2000.0. The mean orbital inclination is is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__SATURNSYSTEM_ORBITALINCLINATION_J2000


        __SATURNSYSTEM_ORBITALPERIOD: float  = 29.447498 # [yr]

        @property
        def SATURNSYSTEM_ORBITALPERIOD(self):
            r"""
            Sidereal orbital period

            #Source: Values derived from the mean longitude rates in Table 5.8.1 in P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, page 316; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [yr]
            """

            return self.__SATURNSYSTEM_ORBITALPERIOD


        __SATURNSYSTEM_ORBITALSEMIMAJORAXIS_J2000: float  = 9.53667594 # [au]

        @property
        def SATURNSYSTEM_ORBITALSEMIMAJORAXIS_J2000(self):
            r"""
            Mean orbital semi-major axis of Saturn, at the standard epoch J2000.0. The mean orbital semi-major axis is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: [au]
            """

            return self.__SATURNSYSTEM_ORBITALSEMIMAJORAXIS_J2000


        __SATURNSYSTEM_RADIALVELOCITYSIGNATURE: float  = 2.8 # [m s^-1]

        @property
        def SATURNSYSTEM_RADIALVELOCITYSIGNATURE(self):
            r"""
            Radial-velocity amplitude of the Sun induced by the Saturn system for 'an edge-on observer' (i.e., an observer in the orbital plane of the Saturn system)

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.4, Equation 1.18, page 9<br/>
            #Basic : false
            #Scalar: true
            #Unit: [m s^-1]
            """

            return self.__SATURNSYSTEM_RADIALVELOCITYSIGNATURE


        __SATURN_ENCOMPASSINGSPHERERADIUS: float  = 6.02680e+07 # [m]

        @property
        def SATURN_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around Saturn which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__SATURN_ENCOMPASSINGSPHERERADIUS


        __SATURN_EQUATORIALRADIUS: float  = 6.02680e+07 # [m]

        @property
        def SATURN_EQUATORIALRADIUS(self):
            r"""
            Equatorial radius of Saturn

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__SATURN_EQUATORIALRADIUS


        __SATURN_FLATTENING: float  = 9.796240e-02

        @property
        def SATURN_FLATTENING(self):
            r"""
            Geometrical flattening factor f of Saturn (f = (a-b)/a)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__SATURN_FLATTENING


        __SATURN_FLUXREDUCTION_MAXIMUM: float  = 0.701 # [%]

        @property
        def SATURN_FLUXREDUCTION_MAXIMUM(self):
            r"""
            Maximum reduction of the solar flux for an observer external to the solar system during a transit of Saturn

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.2, Equation 2.4, page 14<br/>
            #Basic : false
            #Scalar: true
            #Unit: [%]
            """

            return self.__SATURN_FLUXREDUCTION_MAXIMUM


        __SATURN_GEOMETRICALBEDO: float  = 0.47

        @property
        def SATURN_GEOMETRICALBEDO(self):
            r"""
            Geometric albedo of Saturn (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__SATURN_GEOMETRICALBEDO


        __SATURN_JSUB2: float  = 0.016331

        @property
        def SATURN_JSUB2(self):
            r"""
            Dynamical form-factor of Saturn (oblateness or Stokes' second-degree zonal harmonic of the potential)

            #Source: P.R. Weissman, L.-A. McFadden, T.V. Johnson (eds.), 1999, 'Encyclopedia of the Solar System (first edition)', Academic Press, page 342<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__SATURN_JSUB2


        __SATURN_LIGHTDEFLECTION_LIMB: float  = 5980.0 # [10^-6 arcsec]

        @property
        def SATURN_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of Saturn

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__SATURN_LIGHTDEFLECTION_LIMB


        __SATURN_MASS: float  = 5.683190e+26 # [kg]

        @property
        def SATURN_MASS(self):
            r"""
            Mass of Saturn (do not use for high-precision (orbit) calculations)

            #Source: R.A. Jacobson, et al., 2006, 'The gravity field of the Saturnian system from satellite observations and spacecraft tracking data', AJ, 132, 2520-2526; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [kg]
            """

            return self.__SATURN_MASS


        __SATURN_MASSDENSITY_MEAN: float  = 0.6871 # [g cm^-3]

        @property
        def SATURN_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of Saturn

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__SATURN_MASSDENSITY_MEAN


        __SATURN_NORTHROTATIONALPOLE_DECLINATION: float  = 83.537 # [deg]

        @property
        def SATURN_NORTHROTATIONALPOLE_DECLINATION(self):
            r"""
            IAU-recommended value for the declination \delta_0 of the north pole of rotation of Saturn. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__SATURN_NORTHROTATIONALPOLE_DECLINATION


        __SATURN_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE: float  = -0.00000011 # [deg day^-1]

        @property
        def SATURN_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the declination \delta_0 of the north pole of rotation of Saturn. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : false
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__SATURN_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE


        __SATURN_NORTHROTATIONALPOLE_RIGHTASCENSION: float  = 40.589 # [deg]

        @property
        def SATURN_NORTHROTATIONALPOLE_RIGHTASCENSION(self):
            r"""
            IAU-recommended value for the right ascension \alpha_0 of the north pole of rotation of Saturn. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__SATURN_NORTHROTATIONALPOLE_RIGHTASCENSION


        __SATURN_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE: float  = -0.00000099 # [deg day^-1]

        @property
        def SATURN_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the right ascension \alpha_0 of the north pole of rotation of Saturn. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : false
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__SATURN_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE


        __SATURN_POLARRADIUS: float  = 5.43640e+07 # [m]

        @property
        def SATURN_POLARRADIUS(self):
            r"""
            Polar radius of Saturn

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__SATURN_POLARRADIUS


        __SATURN_PRIMEMERIDIAN_EPHEMERISPOSITION: float  = 38.90 # [deg]

        @property
        def SATURN_PRIMEMERIDIAN_EPHEMERISPOSITION(self):
            r"""
            IAU-recommended value for the ephemeris position of the prime meridian of Saturn. The prime meridian refers to the rotation of the magnetic field System III. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__SATURN_PRIMEMERIDIAN_EPHEMERISPOSITION


        __SATURN_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE: float  = 810.7939024 # [deg day^-1]

        @property
        def SATURN_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch of 1.5 January 2000 = JD2451545.0 TDB) of the ephemeris position of the prime meridian of Saturn. The prime meridian refers to the rotation of the magnetic field System III. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__SATURN_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE


        __SATURN_TRANSITPROBABILITY: float  = 0.053 # [%]

        @property
        def SATURN_TRANSITPROBABILITY(self):
            r"""
            Geometric transit probability (Saturn transiting the Sun) for an observer external to the solar system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.1, Equation 2.2, page 14<br/>
            #Basic : false
            #Scalar: true
            #Unit: [%]
            """

            return self.__SATURN_TRANSITPROBABILITY


        __SATURN_TRANSITTIME_MAXIMUM: float  = 1.81 # [day]

        @property
        def SATURN_TRANSITTIME_MAXIMUM(self):
            r"""
            Maximum transit time of Saturn (transiting the Sun) for an observer external to the solar system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.3, Equation 2.5, page 15<br/>
            #Basic : false
            #Scalar: true
            #Unit: [day]
            """

            return self.__SATURN_TRANSITTIME_MAXIMUM


        __SATURN_VONEZEROMAGNITUDE: float  = -8.88 # [mag]

        @property
        def SATURN_VONEZEROMAGNITUDE(self):
            r"""
            V(1,0) magnitude of Saturn (i.e., the visual magnitude of the planet reduced to a distance of 1 au from both the Sun and Earth and phase angle zero). This parameter is also refered to as absolute magnitude in planetary sciences. The magnitude refers to the planetary disk only

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__SATURN_VONEZEROMAGNITUDE


        __SATURN_VOLUMETRICRADIUS: float  = 5.82320e+07 # [m]

        @property
        def SATURN_VOLUMETRICRADIUS(self):
            r"""
            Mean volumetric radius of Saturn

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__SATURN_VOLUMETRICRADIUS


        __SHAPIRO_CONSTANT: float  = 9.8510e-06 # [s]

        @property
        def SHAPIRO_CONSTANT(self):
            r"""
            Shapiro's time-delay constant for the Sun

            #Source: See, e.g., L. Lindegren, D. Dravins, 2003, 'The fundamental definition of <radial velocity>', A&A, 401, 1185, Section 4.3 or P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Equation 5.3211-1, page 295<br/>
            #Basic : false
            #Scalar: true
            #Unit: [s]
            """

            return self.__SHAPIRO_CONSTANT


        __SIC_BULKMODULUS: float  = 206.0 # [GPa]

        @property
        def SIC_BULKMODULUS(self):
            r"""
            The bulk modulus of SiC (Boostec) under isotropic stress; also known as compression modulus

            #Basic : false
            #Scalar: true
            #Unit: [GPa]
            """

            return self.__SIC_BULKMODULUS


        __SIC_COMPRESSIBILITY: float  = 4.857e-12 # [m^2 N^-1]

        @property
        def SIC_COMPRESSIBILITY(self):
            r"""
            The compressibility of SiC (Boostec)

            #Basic : false
            #Scalar: true
            #Unit: [m^2 N^-1]
            """

            return self.__SIC_COMPRESSIBILITY


        __SIC_CRYOGENICLINEARSCALEFACTOR_293K: float  = 1.00021971

        @property
        def SIC_CRYOGENICLINEARSCALEFACTOR_293K(self):
            r"""
            SiC cryogenic linear scale factor \xi between 293 K and 120 K. The assumed linear thermal expansion coefficient of SiC at 293 K is 1.27 ppm K^-1 (see also parameter SiC_LinearThermalCoefficientOfExpansion_293K)

            #Source: A. Mora, 21 June 2011, 'Conversion between image and object space coordinates', GAIA-CH-TN-ESAC-AMF-008, issue 2, revision 1, Section 3<br/>
            #Basic : false
            #Scalar: true
            #Unit: []
            """

            return self.__SIC_CRYOGENICLINEARSCALEFACTOR_293K


        __SIC_LAMECONSTANT_FIRST: float  = 85.0 # [GPa]

        @property
        def SIC_LAMECONSTANT_FIRST(self):
            r"""
            First Lame constant of SiC (Boostec)

            #Basic : false
            #Scalar: true
            #Unit: [GPa]
            """

            return self.__SIC_LAMECONSTANT_FIRST


        __SIC_LAMECONSTANT_SECOND: float  = 181.0 # [GPa]

        @property
        def SIC_LAMECONSTANT_SECOND(self):
            r"""
            Second Lame constant of SiC (Boostec), also known as shear modulus or rigidity

            #Basic : false
            #Scalar: true
            #Unit: [GPa]
            """

            return self.__SIC_LAMECONSTANT_SECOND


        __SIC_LINEARTHERMALCOEFFICIENTOFEXPANSION_100K: float  = 1.0 # [ppm K^-1]

        @property
        def SIC_LINEARTHERMALCOEFFICIENTOFEXPANSION_100K(self):
            r"""
            Average linear thermal expansion coefficient of isotropic homogeneous SiC (Boostec) at 100 K

            #Basic : false
            #Scalar: true
            #Unit: [ppm K^-1]
            """

            return self.__SIC_LINEARTHERMALCOEFFICIENTOFEXPANSION_100K


        __SIC_LINEARTHERMALCOEFFICIENTOFEXPANSION_170K: float  = 1.1 # [ppm K^-1]

        @property
        def SIC_LINEARTHERMALCOEFFICIENTOFEXPANSION_170K(self):
            r"""
            Average linear thermal expansion coefficient of isotropic homogeneous SiC (Boostec) at 170 K

            #Basic : false
            #Scalar: true
            #Unit: [ppm K^-1]
            """

            return self.__SIC_LINEARTHERMALCOEFFICIENTOFEXPANSION_170K


        __SIC_LINEARTHERMALCOEFFICIENTOFEXPANSION_293K: float  = 0.7 # [ppm K^-1]

        @property
        def SIC_LINEARTHERMALCOEFFICIENTOFEXPANSION_293K(self):
            r"""
            Average linear thermal expansion coefficient of isotropic homogeneous SiC (Boostec) at 293 K

            #Basic : false
            #Scalar: true
            #Unit: [ppm K^-1]
            """

            return self.__SIC_LINEARTHERMALCOEFFICIENTOFEXPANSION_293K


        __SIC_MASSDENSITY: float  = 3.16 # [g cm^-3]

        @property
        def SIC_MASSDENSITY(self):
            r"""
            Density of SiC (Boostec)

            #Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9<br/>
            #Basic : true
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__SIC_MASSDENSITY


        __SIC_PWAVESPEED: float  = 11.90 # [km s^-1]

        @property
        def SIC_PWAVESPEED(self):
            r"""
            The P-wave speed in SiC (Boostec); a P-wave (pressure wave) is a longitudinal wave in an elastic medium in which the restoring force is provided by the medium's bulk modulus

            #Basic : false
            #Scalar: true
            #Unit: [km s^-1]
            """

            return self.__SIC_PWAVESPEED


        __SIC_POISSONRATIO: float  = 0.16

        @property
        def SIC_POISSONRATIO(self):
            r"""
            Poisson ratio of SiC (Boostec)

            #Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__SIC_POISSONRATIO


        __SIC_SWAVESPEED: float  = 7.57 # [km s^-1]

        @property
        def SIC_SWAVESPEED(self):
            r"""
            The S-wave speed in SiC (Boostec); an S-wave is a wave in an elastic medium in which the restoring force is provided by shear

            #Basic : false
            #Scalar: true
            #Unit: [km s^-1]
            """

            return self.__SIC_SWAVESPEED


        __SIC_SPECIFICHEATATCONSTANTPRESSURE_100K: float  = 100.0 # [J K^-1 kg^-1]

        @property
        def SIC_SPECIFICHEATATCONSTANTPRESSURE_100K(self):
            r"""
            Specific heat of SiC (Boostec) at 100 K

            #Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9<br/>
            #Basic : true
            #Scalar: true
            #Unit: [J K^-1 kg^-1]
            """

            return self.__SIC_SPECIFICHEATATCONSTANTPRESSURE_100K


        __SIC_SPECIFICHEATATCONSTANTPRESSURE_170K: float  = 320.0 # [J K^-1 kg^-1]

        @property
        def SIC_SPECIFICHEATATCONSTANTPRESSURE_170K(self):
            r"""
            Specific heat of SiC (Boostec) at 170 K

            #Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9<br/>
            #Basic : true
            #Scalar: true
            #Unit: [J K^-1 kg^-1]
            """

            return self.__SIC_SPECIFICHEATATCONSTANTPRESSURE_170K


        __SIC_SPECIFICHEATATCONSTANTPRESSURE_293K: float  = 680.0 # [J K^-1 kg^-1]

        @property
        def SIC_SPECIFICHEATATCONSTANTPRESSURE_293K(self):
            r"""
            Specific heat of SiC (Boostec) at 293 K

            #Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9<br/>
            #Basic : true
            #Scalar: true
            #Unit: [J K^-1 kg^-1]
            """

            return self.__SIC_SPECIFICHEATATCONSTANTPRESSURE_293K


        __SIC_SPECIFICSTIFNESS: float  = 133.0 # [10^6 N m kg^-1]

        @property
        def SIC_SPECIFICSTIFNESS(self):
            r"""
            Specific stiffness of SiC (Boostec); also known as specific modulus

            #Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^6 N m kg^-1]
            """

            return self.__SIC_SPECIFICSTIFNESS


        __SIC_STRENGTH_AVERAGE: float  = 390.0 # [10^6 Pa]

        @property
        def SIC_STRENGTH_AVERAGE(self):
            r"""
            Average strength of SiC (Boostec) in the four-point bending test

            #Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9<br/>
            #Basic : true
            #Scalar: true
            #Unit: [10^6 Pa]
            """

            return self.__SIC_STRENGTH_AVERAGE


        __SIC_THERMALCONDUCTIVITY_100K: float  = 180.0 # [W m^-1 K^-1]

        @property
        def SIC_THERMALCONDUCTIVITY_100K(self):
            r"""
            Thermal conductivity of isotropic homogeneous SiC (Boostec) at 100 K

            #Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9<br/>
            #Basic : true
            #Scalar: true
            #Unit: [W m^-1 K^-1]
            """

            return self.__SIC_THERMALCONDUCTIVITY_100K


        __SIC_THERMALCONDUCTIVITY_170K: float  = 220.0 # [W m^-1 K^-1]

        @property
        def SIC_THERMALCONDUCTIVITY_170K(self):
            r"""
            Thermal conductivity of isotropic homogeneous SiC (Boostec) at 170 K

            #Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9<br/>
            #Basic : true
            #Scalar: true
            #Unit: [W m^-1 K^-1]
            """

            return self.__SIC_THERMALCONDUCTIVITY_170K


        __SIC_THERMALCONDUCTIVITY_293K: float  = 190.0 # [W m^-1 K^-1]

        @property
        def SIC_THERMALCONDUCTIVITY_293K(self):
            r"""
            Thermal conductivity of isotropic homogeneous SiC (Boostec) at 293 K

            #Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9<br/>
            #Basic : true
            #Scalar: true
            #Unit: [W m^-1 K^-1]
            """

            return self.__SIC_THERMALCONDUCTIVITY_293K


        __SIC_THERMALDIFFUSIVITY_100K: float  = 570.0 # [m^2 s^-1]

        @property
        def SIC_THERMALDIFFUSIVITY_100K(self):
            r"""
            Thermal diffusivity of SiC (Boostec) at 100 K

            #Basic : false
            #Scalar: true
            #Unit: [m^2 s^-1]
            """

            return self.__SIC_THERMALDIFFUSIVITY_100K


        __SIC_THERMALDIFFUSIVITY_170K: float  = 218.0 # [m^2 s^-1]

        @property
        def SIC_THERMALDIFFUSIVITY_170K(self):
            r"""
            Thermal diffusivity of SiC (Boostec) at 170 K

            #Basic : false
            #Scalar: true
            #Unit: [m^2 s^-1]
            """

            return self.__SIC_THERMALDIFFUSIVITY_170K


        __SIC_THERMALDIFFUSIVITY_293K: float  = 88.0 # [m^2 s^-1]

        @property
        def SIC_THERMALDIFFUSIVITY_293K(self):
            r"""
            Thermal diffusivity of SiC (Boostec) at 293 K

            #Basic : false
            #Scalar: true
            #Unit: [m^2 s^-1]
            """

            return self.__SIC_THERMALDIFFUSIVITY_293K


        __SIC_VOLUMETRICTHERMALCOEFFICIENTOFEXPANSION_100K: float  = 3.1 # [ppm K^-1]

        @property
        def SIC_VOLUMETRICTHERMALCOEFFICIENTOFEXPANSION_100K(self):
            r"""
            Average volumetric thermal expansion coefficient of isotropic homogeneous SiC (Boostec) at 100 K

            #Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9<br/>
            #Basic : true
            #Scalar: true
            #Unit: [ppm K^-1]
            """

            return self.__SIC_VOLUMETRICTHERMALCOEFFICIENTOFEXPANSION_100K


        __SIC_VOLUMETRICTHERMALCOEFFICIENTOFEXPANSION_170K: float  = 3.4 # [ppm K^-1]

        @property
        def SIC_VOLUMETRICTHERMALCOEFFICIENTOFEXPANSION_170K(self):
            r"""
            Average volumetric thermal expansion coefficient of isotropic homogeneous SiC (Boostec) at 170 K

            #Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9<br/>
            #Basic : true
            #Scalar: true
            #Unit: [ppm K^-1]
            """

            return self.__SIC_VOLUMETRICTHERMALCOEFFICIENTOFEXPANSION_170K


        __SIC_VOLUMETRICTHERMALCOEFFICIENTOFEXPANSION_293K: float  = 2.2 # [ppm K^-1]

        @property
        def SIC_VOLUMETRICTHERMALCOEFFICIENTOFEXPANSION_293K(self):
            r"""
            Average volumetric thermal expansion coefficient of isotropic homogeneous SiC (Boostec) at 293 K

            #Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9<br/>
            #Basic : true
            #Scalar: true
            #Unit: [ppm K^-1]
            """

            return self.__SIC_VOLUMETRICTHERMALCOEFFICIENTOFEXPANSION_293K


        __SIC_WEIBULLMODULUS_FIREDSLASHSINTERED: float  = 12.0

        @property
        def SIC_WEIBULLMODULUS_FIREDSLASHSINTERED(self):
            r"""
            The Weibull modulus of fired/sintered SiC (Boostec). The Weibull modulus is the slope of the linear plot of log(log(P^-1)) versus log(\sigma), where P is the probability of survival under a given stress \sigma. The numerical value is a lower limit for the Weibull modulus

            #Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__SIC_WEIBULLMODULUS_FIREDSLASHSINTERED


        __SIC_WEIBULLMODULUS_MACHINED: float  = 10.0

        @property
        def SIC_WEIBULLMODULUS_MACHINED(self):
            r"""
            The Weibull modulus of machined SiC (Boostec). The Weibull modulus is the slope of the linear plot of log(log(P^-1)) versus log(\sigma), where P is the probability of survival under a given stress \sigma. The numerical value is a lower limit for the Weibull modulus

            #Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__SIC_WEIBULLMODULUS_MACHINED


        __SIC_YOUNGMODULUS: float  = 420.0 # [GPa]

        @property
        def SIC_YOUNGMODULUS(self):
            r"""
            Young's modulus of SiC (Boostec); a measure for the resistance of a material to elastic (recoverable) deformation under load. Also known as elastic modulus and tension modulus

            #Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9<br/>
            #Basic : true
            #Scalar: true
            #Unit: [GPa]
            """

            return self.__SIC_YOUNGMODULUS


        __SI_BANDGAP_0K: float  = 1.166 # [eV]

        @property
        def SI_BANDGAP_0K(self):
            r"""
            Silicon bandgap, in eV, at 0 K

            #Source: B. van Zeghbroeck, 1 January 2007, 'Principles of Semiconductor Devices', http://ece-www.colorado.edu/~bart/book/<br/>
            #Basic : true
            #Scalar: true
            #Unit: [eV]
            """

            return self.__SI_BANDGAP_0K


        __SI_BANDGAP_160K: float  = 1.151 # [eV]

        @property
        def SI_BANDGAP_160K(self):
            r"""
            Silicon bandgap, in eV, at 160 K. The energy bandgap of Silicon (and semi-conductors in general) decreases with increasing temperature. This is explained as follows: an increased temperature enhances the amplitude of atomic vibrations due to the increased thermal energy; this causes the interatomic spacing to increase; this decreases the potential seen by the electrons in the material; this, finally, reduces the size of the energy bandgap. The temperature dependence of the energy bandgap has been experimentally determined to be E_g(T) = E_g(0) - \alpha * T^2 * (T + \beta)^-1, where T is the temperature in K,  E_g(0) is the bandgap energy in eV at zero absolute temperature (parameter Si_Bandgap_0K), and \alpha and \beta are material constants (parameters Si_Constant_Alpha and Si_Constant_Beta for Silicon)

            #Basic : false
            #Scalar: true
            #Unit: [eV]
            """

            return self.__SI_BANDGAP_160K


        __SI_CONSTANT_ALPHA: float  = 4.730e-04 # [eV K^-1]

        @property
        def SI_CONSTANT_ALPHA(self):
            r"""
            Silicon constant \alpha, in eV K^-1. The energy bandgap of Silicon (and semi-conductors in general) decreases with increasing temperature. This is explained as follows: an increased temperature enhances the amplitude of atomic vibrations due to the increased thermal energy; this causes the interatomic spacing to increase; this decreases the potential seen by the electrons in the material; this, finally, reduces the size of the energy bandgap. The temperature dependence of the energy bandgap has been experimentally determined to be E_g(T) = E_g(0) - \alpha * T^2 * (T + \beta)^-1, where T is the temperature in K,  E_g(0) is the bandgap energy in eV at zero absolute temperature (parameter Si_Bandgap_0K), and \alpha and \beta are material constants (parameters Si_Constant_Alpha and Si_Constant_Beta for Silicon)

            #Source: B. van Zeghbroeck, 1 January 2007, 'Principles of Semiconductor Devices', http://ece-www.colorado.edu/~bart/book/<br/>
            #Basic : true
            #Scalar: true
            #Unit: [eV K^-1]
            """

            return self.__SI_CONSTANT_ALPHA


        __SI_CONSTANT_BETA: float  = 636.0 # [K]

        @property
        def SI_CONSTANT_BETA(self):
            r"""
            Silicon constant \beta, in K. The energy bandgap of Silicon (and semi-conductors in general) decreases with increasing temperature. This is explained as follows: an increased temperature enhances the amplitude of atomic vibrations due to the increased thermal energy; this causes the interatomic spacing to increase; this decreases the potential seen by the electrons in the material; this, finally, reduces the size of the energy bandgap. The temperature dependence of the energy bandgap has been experimentally determined to be E_g(T) = E_g(0) - \alpha * T^2 * (T + \beta)^-1, where T is the temperature in K,  E_g(0) is the bandgap energy in eV at zero absolute temperature (parameter Si_Bandgap_0K), and \alpha and \beta are material constants (parameters Si_Constant_Alpha and Si_Constant_Beta for Silicon)

            #Source: B. van Zeghbroeck, 1 January 2007, 'Principles of Semiconductor Devices', http://ece-www.colorado.edu/~bart/book/<br/>
            #Basic : true
            #Scalar: true
            #Unit: [K]
            """

            return self.__SI_CONSTANT_BETA


        __SI_CUTOFFWAVELENGTH_160K: float  = 1077.39 # [nm]

        @property
        def SI_CUTOFFWAVELENGTH_160K(self):
            r"""
            Silicon cut-off wavelength, in nm, at 160 K

            #Basic : false
            #Scalar: true
            #Unit: [nm]
            """

            return self.__SI_CUTOFFWAVELENGTH_160K


        __SI_DIFFUSIONCOEFFICIENT: float  = 0.0039 # [m^2 s^-1]

        @property
        def SI_DIFFUSIONCOEFFICIENT(self):
            r"""
            Silicon diffusion coefficient

            #Source: J.R. Janesick, 2001, 'Scientific CCDs', SPIE, Bellingham, Washington, Example 4.17, page 348<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^2 s^-1]
            """

            return self.__SI_DIFFUSIONCOEFFICIENT


        __SI_DIFFUSIONLENGTH: float  = 0.0006 # [m]

        @property
        def SI_DIFFUSIONLENGTH(self):
            r"""
            Diffusion length in the epitaxial Silicon

            #Source: J.R. Janesick, 2001, 'Scientific CCDs', SPIE, Bellingham, Washington, Example 4.17, page 348<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__SI_DIFFUSIONLENGTH


        __SI_OPTICALABSORPTIONCOEFFICIENT_160K: str  = "Nature/Si_OpticalAbsorptionCoefficient_160K_001.fits"

        @property
        def SI_OPTICALABSORPTIONCOEFFICIENT_160K(self):
            r"""
            Silicon optical absorption coefficient as a function of (photon) wavelength. First column: wavelength \lambda (in nm; from 200.0 to 1100.0). Second column: Silicon optical absorption coefficient \alpha (in [10^-6 m]^-1) at T = 160 K. Third column: Silicon photon absorption depth L_A = \alpha^-1 (in 10^-6 m) at T = 160 K

            #Source: K. Rajkanan, R. Singh, J. Shewchun, 1979, 'Absorption coefficient of Silicon for solar cell calculations', Solid-State Electronics, 22, 793-795; analytical phenomenological model presented in this reference has an accuracy of about 20%; formal validity ranges are 20-500 K in temperature T and 1.1-4.0 eV, i.e., 310-1127 nm, in photon energy (and wavelength)<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__SI_OPTICALABSORPTIONCOEFFICIENT_160K


        __SI_PHOTONABSORPTIONDEPTH_160K: str  = "Nature/Si_PhotonAbsorptionDepth_160K_001.fits"

        @property
        def SI_PHOTONABSORPTIONDEPTH_160K(self):
            r"""
            Silicon photon absorption depth as a function of (photon) wavelength. First column: wavelength \lambda (in nm; from 200.0 to 1100.0). Second column: Silicon optical absorption coefficient \alpha (in [10^-6 m]^-1) at T = 160 K. Third column: Silicon photon absorption depth L_A = \alpha^-1 (in 10^-6 m) at T = 160 K

            #Source: K. Rajkanan, R. Singh, J. Shewchun, 1979, 'Absorption coefficient of Silicon for solar cell calculations', Solid-State Electronics, 22, 793-795; analytical phenomenological model presented in this reference has an accuracy of about 20%; formal validity ranges are 20-500 K in temperature T and 1.1-4.0 eV, i.e., 310-1127 nm, in photon energy (and wavelength)<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__SI_PHOTONABSORPTIONDEPTH_160K


        __SIDEREALYEAR_J2000DAY: float  = 365.256363004 # [day]

        @property
        def SIDEREALYEAR_J2000DAY(self):
            r"""
            Number of days per sidereal year

            #Source: J.L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touze, G. Francou, J. Laskar, 1994, 'Numerical expressions for precession formulae and mean elements for the Moon and the planets', A&A, 282, 663 (1994A&A...282..663S)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [day]
            """

            return self.__SIDEREALYEAR_J2000DAY


        __SKY_NUMBEROFGALAXIES_G = [ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.009,  0.038,  0.104,  0.241,  0.509,  0.983,  1.822,  3.326,  6.048,  11.075,  21.145,  41.602,  82.794,  166.012 ] # [10^6 objects]

        @property
        def SKY_NUMBEROFGALAXIES_G(self):
            r"""
            Cumulative number of galaxies, integrated over the full sky, as function of G magnitude for the following limits: up to G = 4.5 mag, up to G = 5.0 mag, ..., up to G = 20.5 mag, and up to G = 21.0 mag (numerical values from parameter Sky_ObjectDensity_003, extracted by B. Holl on 4 October 2012 using GT SkyDensity.class). See parameter Sky_NumberOfStars_G for star counts and Sky_NumberOfObjects_G for object counts

            #Source: The data used to produce the standard Gaia Galaxy model have been provided by the teams of the Observatoire de Besancon (A. Robin, C. Reyle, et al.) and the Observatorio Astronomico di Torino (R. Drimmel, et al.); the harmonisation has been provided by the SWG (X. Luri). See U. Lammers, 22 April 2005, 'Gaia Standard Galaxy Model Access Software (GSGMAS); User Guide', GAIA-UL-010, issue 1, revision 0. See also R. Drimmel, et al., 19 July 2005, 'Recommendations on the use of estimated star counts for Gaia studies', SWG-RD-004<br/>
            #Basic : true
            #Scalar: true
            #Unit: [10^6 objects]
            """

            return self.__SKY_NUMBEROFGALAXIES_G


        __SKY_NUMBEROFGALAXIES_GRVS = [ 0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.010,  0.040,  0.111,  0.255,  0.537,  1.032,  1.912,  3.490,  6.349,  11.607 ] # [10^6 objects]

        @property
        def SKY_NUMBEROFGALAXIES_GRVS(self):
            r"""
            Cumulative number of galaxies, integrated over the full sky, as function of G_RVS (= C1M861 = RVF) magnitude for the following limits: up to G_RVS = 4.5 mag, up to G_RVS = 5.0 mag, ..., up to G_RVS = 17.5 mag, and up to G_RVS = 18.0 mag (numerical values from parameter Sky_ObjectDensity_003, extracted by B. Holl on 4 October 2012 using GT SkyDensity.class). See parameter Sky_NumberOfStars_GRVS for star counts and Sky_NumberOfObjects_GRVS for object counts

            #Source: The data used to produce the standard Gaia Galaxy model have been provided by the teams of the Observatoire de Besancon (A. Robin, C. Reyle, et al.) and the Observatorio Astronomico di Torino (R. Drimmel, et al.); the harmonisation has been provided by the SWG (X. Luri). See U. Lammers, 22 April 2005, 'Gaia Standard Galaxy Model Access Software (GSGMAS); User Guide', GAIA-UL-010, issue 1, revision 0. See also R. Drimmel, et al., 19 July 2005, 'Recommendations on the use of estimated star counts for Gaia studies', SWG-RD-004<br/>
            #Basic : true
            #Scalar: true
            #Unit: [10^6 objects]
            """

            return self.__SKY_NUMBEROFGALAXIES_GRVS


        __SKY_NUMBEROFOBJECTS_G = [ 0.003,  0.007,  0.014,  0.026,  0.044,  0.070,  0.107,  0.160,  0.241,  0.356,  0.522,  0.768,  1.141,  1.715,  2.640,  4.053,  6.267,  9.708,  14.749,  22.060,  32.556,  47.488,  68.545,  97.837,  137.944,  192.408,  265.629,  363.261,  494.546,  672.011,  911.605,  1233.940,  1667.145,  2247.560 ] # [10^6 objects]

        @property
        def SKY_NUMBEROFOBJECTS_G(self):
            r"""
            Cumulative number of objects (stars + galaxies; see parameters Sky_NumberOfStars_G and Sky_NumberOfGalaxies_G), integrated over the full sky, as function of G magnitude for the following limits: up to G = 4.5 mag, up to G = 5.0 mag, ..., up to G = 20.5 mag, and up to G = 21.0 mag (numerical values from parameter Sky_ObjectDensity_003, extracted by B. Holl on 4 October 2012 using GT SkyDensity.class)

            #Source: The data used to produce the standard Gaia Galaxy model have been provided by the teams of the Observatoire de Besancon (A. Robin, C. Reyle, et al.) and the Observatorio Astronomico di Torino (R. Drimmel, et al.); the harmonisation has been provided by the SWG (X. Luri). See U. Lammers, 22 April 2005, 'Gaia Standard Galaxy Model Access Software (GSGMAS); User Guide', GAIA-UL-010, issue 1, revision 0. See also R. Drimmel, et al., 19 July 2005, 'Recommendations on the use of estimated star counts for Gaia studies', SWG-RD-004<br/>
            #Basic : true
            #Scalar: true
            #Unit: [10^6 objects]
            """

            return self.__SKY_NUMBEROFOBJECTS_G


        __SKY_NUMBEROFOBJECTS_GRVS = [ 0.018,  0.042,  0.072,  0.111,  0.160,  0.227,  0.320,  0.460,  0.675,  0.995,  1.471,  2.196,  3.293,  4.902,  7.280,  10.876,  16.184,  23.915,  35.157,  51.387,  74.697,  107.743,  154.027,  217.589,  303.945,  420.851,  580.228,  800.547 ] # [10^6 objects]

        @property
        def SKY_NUMBEROFOBJECTS_GRVS(self):
            r"""
            Cumulative number of objects (stars + galaxies; see parameters Sky_NumberOfStars_GRVS and Sky_NumberOfGalaxies_GRVS), integrated over the full sky, as function of G_RVS (= C1M861 = RVF) magnitude for the following limits: up to G_RVS = 4.5 mag, up to G_RVS = 5.0 mag, ..., up to G_RVS = 17.5 mag, and up to G_RVS = 18.0 mag (numerical values from parameter Sky_ObjectDensity_003, extracted by B. Holl on 4 October 2012 using GT SkyDensity.class)

            #Source: The data used to produce the standard Gaia Galaxy model have been provided by the teams of the Observatoire de Besancon (A. Robin, C. Reyle, et al.) and the Observatorio Astronomico di Torino (R. Drimmel, et al.); the harmonisation has been provided by the SWG (X. Luri). See U. Lammers, 22 April 2005, 'Gaia Standard Galaxy Model Access Software (GSGMAS); User Guide', GAIA-UL-010, issue 1, revision 0. See also R. Drimmel, et al., 19 July 2005, 'Recommendations on the use of estimated star counts for Gaia studies', SWG-RD-004<br/>
            #Basic : true
            #Scalar: true
            #Unit: [10^6 objects]
            """

            return self.__SKY_NUMBEROFOBJECTS_GRVS


        __SKY_NUMBEROFSTARS_DISKGIANT: float  = 92.0 # [10^6 stars]

        @property
        def SKY_NUMBEROFSTARS_DISKGIANT(self):
            r"""
            Predicted number of observed stars (i.e., G <= 20.00 mag) in the disk (giants)

            #Source: ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 6.4.4 and Table 6.3, pages 239-240 (Galaxy model from J. Torra, et al., 1999, 'Predicting Gaia Observations from a Star-Count Model', Baltic Astronomy, 8, 171 and extinction law from J. Hakkila, et al., 1997, 'A Computerised Model of Large-Scale Visual Interstellar Extinction', AJ, 114, 2043)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [10^6 stars]
            """

            return self.__SKY_NUMBEROFSTARS_DISKGIANT


        __SKY_NUMBEROFSTARS_DISKMSPLUSWD: float  = 779.0 # [10^6 stars]

        @property
        def SKY_NUMBEROFSTARS_DISKMSPLUSWD(self):
            r"""
            Predicted number of observed stars (i.e., G <= 20.00 mag) in the disk (main sequence stars plus white dwarfs)

            #Source: ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 6.4.4 and Table 6.3, pages 239-240 (Galaxy model from J. Torra, et al., 1999, 'Predicting Gaia Observations from a Star-Count Model', Baltic Astronomy, 8, 171 and extinction law from J. Hakkila, et al., 1997, 'A Computerised Model of Large-Scale Visual Interstellar Extinction', AJ, 114, 2043)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [10^6 stars]
            """

            return self.__SKY_NUMBEROFSTARS_DISKMSPLUSWD


        __SKY_NUMBEROFSTARS_G = [ 0.003,  0.007,  0.014,  0.026,  0.044,  0.070,  0.107,  0.160,  0.241,  0.356,  0.522,  0.768,  1.141,  1.715,  2.640,  4.053,  6.267,  9.708,  14.749,  22.060,  32.547,  47.450,  68.441,  97.595,  137.435,  191.424,  263.806,  359.935,  488.497,  660.937,  890.460,  1192.338,  1584.350,  2081.548 ] # [10^6 objects]

        @property
        def SKY_NUMBEROFSTARS_G(self):
            r"""
            Cumulative number of stars, integrated over the full sky, as function of G magnitude for the following limits: up to G = 4.5 mag, up to G = 5.0 mag, ..., up to G = 20.5 mag, and up to G = 21.0 mag (numerical values from parameter Sky_ObjectDensity_003, extracted by B. Holl on 4 October 2012 using GT SkyDensity.class). See parameter Sky_NumberOfGalaxies_G for galaxy counts and Sky_NumberOfObjects_G for object counts

            #Source: The data used to produce the standard Gaia Galaxy model have been provided by the teams of the Observatoire de Besancon (A. Robin, C. Reyle, et al.) and the Observatorio Astronomico di Torino (R. Drimmel, et al.); the harmonisation has been provided by the SWG (X. Luri). See U. Lammers, 22 April 2005, 'Gaia Standard Galaxy Model Access Software (GSGMAS); User Guide', GAIA-UL-010, issue 1, revision 0. See also R. Drimmel, et al., 19 July 2005, 'Recommendations on the use of estimated star counts for Gaia studies', SWG-RD-004<br/>
            #Basic : true
            #Scalar: true
            #Unit: [10^6 objects]
            """

            return self.__SKY_NUMBEROFSTARS_G


        __SKY_NUMBEROFSTARS_GRVS = [ 0.018,  0.042,  0.072,  0.111,  0.160,  0.227,  0.320,  0.460,  0.675,  0.995,  1.471,  2.196,  3.293,  4.902,  7.280,  10.876,  16.184,  23.915,  35.147,  51.347,  74.587,  107.488,  153.490,  216.557,  302.033,  417.361,  573.879,  788.939 ] # [10^6 objects]

        @property
        def SKY_NUMBEROFSTARS_GRVS(self):
            r"""
            Cumulative number of stars, integrated over the full sky, as function of G_RVS (= C1M861 = RVF) magnitude for the following limits: up to G_RVS = 4.5 mag, up to G_RVS = 5.0 mag, ..., up to G_RVS = 17.5 mag, and up to G_RVS = 18.0 mag (numerical values from parameter Sky_ObjectDensity_003, extracted by B. Holl on 4 October 2012 using GT SkyDensity.class). See parameter Sky_NumberOfGalaxies_GRVS for galaxy counts and Sky_NumberOfObjects_GRVS for object counts

            #Source: The data used to produce the standard Gaia Galaxy model have been provided by the teams of the Observatoire de Besancon (A. Robin, C. Reyle, et al.) and the Observatorio Astronomico di Torino (R. Drimmel, et al.); the harmonisation has been provided by the SWG (X. Luri). See U. Lammers, 22 April 2005, 'Gaia Standard Galaxy Model Access Software (GSGMAS); User Guide', GAIA-UL-010, issue 1, revision 0. See also R. Drimmel, et al., 19 July 2005, 'Recommendations on the use of estimated star counts for Gaia studies', SWG-RD-004<br/>
            #Basic : true
            #Scalar: true
            #Unit: [10^6 objects]
            """

            return self.__SKY_NUMBEROFSTARS_GRVS


        __SKY_NUMBEROFSTARS_SPHEROID: float  = 67.0 # [10^6 stars]

        @property
        def SKY_NUMBEROFSTARS_SPHEROID(self):
            r"""
            Predicted number of observed stars (i.e., G <= 20.00 mag) in the spheroid (including the bulge)

            #Source: ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 6.4.4 and Table 6.3, pages 239-240 (Galaxy model from J. Torra, et al., 1999, 'Predicting Gaia Observations from a Star-Count Model', Baltic Astronomy, 8, 171 and extinction law from J. Hakkila, et al., 1997, 'A Computerised Model of Large-Scale Visual Interstellar Extinction', AJ, 114, 2043)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [10^6 stars]
            """

            return self.__SKY_NUMBEROFSTARS_SPHEROID


        __SKY_NUMBEROFSTARS_THICKDISK: float  = 97.0 # [10^6 stars]

        @property
        def SKY_NUMBEROFSTARS_THICKDISK(self):
            r"""
            Predicted number of observed stars (i.e., G <= 20.00 mag) in the thick disk

            #Source: ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 6.4.4 and Table 6.3, pages 239-240 (Galaxy model from J. Torra, et al., 1999, 'Predicting Gaia Observations from a Star-Count Model', Baltic Astronomy, 8, 171 and extinction law from J. Hakkila, et al., 1997, 'A Computerised Model of Large-Scale Visual Interstellar Extinction', AJ, 114, 2043)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [10^6 stars]
            """

            return self.__SKY_NUMBEROFSTARS_THICKDISK


        __SKY_NUMBEROFSTARS_TOTAL: float  = 1.0350e+09 # [stars]

        @property
        def SKY_NUMBEROFSTARS_TOTAL(self):
            r"""
            Predicted total number of observed stars (i.e., G <= 20.00 mag)

            #Source: See also ESA, 21 May 2013, 'Gaia mission requirements document (MRD)', GAIA-EST-RD-00553, issue 3, revision 1, Requirement SCI-130<br/>
            #Basic : false
            #Scalar: true
            #Unit: [stars]
            """

            return self.__SKY_NUMBEROFSTARS_TOTAL


        __SKY_OBJECTDENSITY: str  = "Nature/Sky_ObjectDensity_003.fits"

        @property
        def SKY_OBJECTDENSITY(self):
            r"""
            Standard Gaia Galaxy model, providing predictions of object densities in the three Gaia-2-relevant photometric bands G, GS, and G_RVS (= C1M861 = RVF) as a function of limiting magnitude for any point on the celestial sphere. Object densities are defined on a Hierarchical Triangular Mesh of level 6 comprised of almost-equal-area cells of approximately 1 square degree in size. The model is a synthesis of count predictions from two distinct sources: (i) the Bescancon Galaxy model (http://model.obs-besancon.fr/), in conjunction with the extinction law of R. Drimmel, for stars; and (ii) the GSC-II catalogue (http://gsss.stsci.edu/Catalogs/GSC/GSC2/GSC2.htm) for stars and galaxies. The combined result can be summarised as follows: (a) Bright stars: 4 <= G <= 12.5: GSC-II (Tycho-2 catalogue included); 4 <= GS <= 12.5: GSC-II (Tycho-2 catalogue included); 4 <= G_RVS <= 12.5: GSC-II (Tycho-2 catalogue included); (b) Faint stars: 12.5 <= G <= 21: Bescancon; 12.5 <= GS <= 21: Bescancon; 12.5 <= G_RVS <= 18: Bescancon; (c) Stars in special regions (LMC/SMC): counts for full range of magnitudes, for circle with radius 10 deg (LMC) and 7 deg (SMC), taken from the GSC-II catalogue; (d) galaxies: G <= 21: GSC-II; GS <= 21: GSC-II; G_RVS <= 18: GSC-II. The following Bintables are defined: (1) STAR-G-MAGGRID: magnitude grid for column DENSITY in following extension STAR-G-DENSITY (the magnitude grids STAR-G-MAGGRID and STAR-C1M861-MAGGRID are defined through parameters Sky_ObjectDensity_MagMin, Sky_ObjectDensity_MagMax, and Sky_ObjectDensity_MagBinWidth); (2) STAR-G-DENSITY: table with star magnitude counts in column DENSITY per HTM cell in the G band; each row corresponds to one cell. Table columns are: HTMID: unique HTM cell designator; AREA: area of cell [square deg]; ALPHA: right ascension of cell centre [h]; DELTA: declination of cell centre [deg]; LGAL: Galactic longitude of cell centre [deg]; BGAL: Galactic latitude of cell centre [deg]; DENSITY: non-cumulative number of objects per square degree in cell in 0.5 magnitude-bin interval, i.e., DENSITY[i] is the number of stars per square degree in the magnitude interval [STAR-G-MAGGRID[i], STAR-G-MAGGRID[i]+0.5]; (3) GAL-G-MAGGRID: same as STAR-G-MAGGRID, but for galaxies (the magnitude grids STAR-G-MAGGRID and STAR-C1M861-MAGGRID are defined through parameters Sky_ObjectDensity_MagMin, Sky_ObjectDensity_MagMax, and Sky_ObjectDensity_MagBinWidth); (4) GAL-G-DENSITY: same as STAR-G-DENSITY, but for galaxies; (5) STAR-GS-MAGGRID: same as STAR-G-MAGGRID, but in GS band; (6) STAR-GS-DENSITY: same as STAR-G-DENSITY, but in GS band; (7) GAL-GS-MAGGRID: same as STAR-GS-MAGGRID, but for galaxies; (8) GAL-GS-DENSITY: same as STAR-GS-DENSITY, but for galaxies; (9) STAR-C1M861-MAGGRID: same as STAR-G-MAGGRID, but in G_RVS band; (10) STAR-C1M861-DENSITY: same as STAR-G-DENSITY, but in G_RVS band; (11) GAL-C1M861-MAGGRID: same as STAR-C1M861-MAGGRID, but for galaxies; and (12) GAL-C1M861-DENSITY: same as STAR-C1M861-DENSITY, but for galaxies

            #Source: The data used to produce the standard Gaia Galaxy model have been provided by the teams of the Observatoire de Besancon (A. Robin, C. Reyle, et al.) and the Observatorio Astronomico di Torino (R. Drimmel, et al.); the harmonisation has been provided by the SWG (X. Luri). See U. Lammers, 22 April 2005, 'Gaia Standard Galaxy Model Access Software (GSGMAS); User Guide', GAIA-UL-010, issue 1, revision 0. See also R. Drimmel, et al., 19 July 2005, 'Recommendations on the use of estimated star counts for Gaia studies', SWG-RD-004<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__SKY_OBJECTDENSITY


        __SKY_OBJECTDENSITY_MAGBINWIDTH: float  = 0.5 # [mag]

        @property
        def SKY_OBJECTDENSITY_MAGBINWIDTH(self):
            r"""
            Magnitude step (bin size) of the magnitude grid (STAR-G-MAGGRID and STAR-C1M861-MAGGRID) underlying the object counts in parameter Sky_ObjectDensity

            #Source: Extracted from parameter Sky_ObjectDensity_003<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__SKY_OBJECTDENSITY_MAGBINWIDTH


        __SKY_OBJECTDENSITY_MAGMAX: float  = 21.0 # [mag]

        @property
        def SKY_OBJECTDENSITY_MAGMAX(self):
            r"""
            Maximum magnitude of the magnitude grid (STAR-G-MAGGRID and STAR-C1M861-MAGGRID) underlying the object counts in parameter Sky_ObjectDensity

            #Source: Extracted from parameter Sky_ObjectDensity_003<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__SKY_OBJECTDENSITY_MAGMAX


        __SKY_OBJECTDENSITY_MAGMIN: float  = 4.0 # [mag]

        @property
        def SKY_OBJECTDENSITY_MAGMIN(self):
            r"""
            Minimum magnitude of the magnitude grid (STAR-G-MAGGRID and STAR-C1M861-MAGGRID) underlying the object counts in parameter Sky_ObjectDensity

            #Source: Extracted from parameter Sky_ObjectDensity_003<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__SKY_OBJECTDENSITY_MAGMIN


        __SKY_STARDENSITY_AVERAGE: float  = 25089.0 # [stars deg^-2]

        @property
        def SKY_STARDENSITY_AVERAGE(self):
            r"""
            Sky-averaged density of stars with G <= 20.00 mag. Note that a value of 25.1E3 stars deg^-2 is given in ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 6.4.5 and Table 6.6, pages 240-242

            #Source: See also ESA, 21 May 2013, 'Gaia mission requirements document (MRD)', GAIA-EST-RD-00553, issue 3, revision 1, Section 4.1<br/>
            #Basic : false
            #Scalar: true
            #Unit: [stars deg^-2]
            """

            return self.__SKY_STARDENSITY_AVERAGE


        __SKY_STARDENSITY_DESIGN: float  = 600000.0 # [stars deg^-2]

        @property
        def SKY_STARDENSITY_DESIGN(self):
            r"""
            Density of stars with G <= 20.00 mag below which Gaia is designed to operate nominally: nominal observations and all astrometric requirements shall be achieved in the two superimposed fields of view computed using the design object density in one instrument field of view plus the typical object density in the other instrument field of view (Requirement SCI-150). Requirement SCI-160 guarantees that a strategy to observe high-density sky regions (e.g., Baade's window), with stellar densities up to the worst-case star density, will be implemented (if higher densities than the worst-case density are encountered, the brightest objects up to the worst-case density will be observed). One option is to have several successive transits of the same sky region at a reduced precession rate using a modified scanning law (MSL; see L. Lindegren, 7 April 2005, 'Multi-pass scanning across Baade`s Window', GAIA-LL-058, issue 1, revision 0 and L. Lindegren, 22 August 2010, 'Reference Scanning Law for Gaia', GAIA-C3-TN-LU-LL-085-01)

            #Source: ESA, 21 May 2013, 'Gaia mission requirements document (MRD)', GAIA-EST-RD-00553, issue 3, revision 1, Section 4.1, Requirements SCI-150 and SCI-160<br/>
            #Basic : true
            #Scalar: true
            #Unit: [stars deg^-2]
            """

            return self.__SKY_STARDENSITY_DESIGN


        __SKY_STARDENSITY_LATITUDE = [ 101.6,  79.8,  31.2,  11.4,  3.8 ] # [1E3 stars deg^-2]

        @property
        def SKY_STARDENSITY_LATITUDE(self):
            r"""
            Density of stars with G <= 20.00 mag in the Galactic-latitude ranges 0-5 degrees, 5-10 degrees, 10-20 degrees, 20-30 degrees, and 30-90 degrees

            #Source: ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 6.4.5 and Table 6.6, pages 240-242<br/>
            #Basic : true
            #Scalar: true
            #Unit: [1E3 stars deg^-2]
            """

            return self.__SKY_STARDENSITY_LATITUDE


        __SKY_STARDENSITY_RVSREDUCTIONFACTOR: float  = 6.0

        @property
        def SKY_STARDENSITY_RVSREDUCTIONFACTOR(self):
            r"""
            The brighter star-selection limiting magnitude of RVS compared to SM/AF/BP/RP (G_RVS (= C1M861 = RVF) = 17.00 mag versus G = 20.00 mag, respectively) is assumed to correspond to a factor 6 reduction in star density and the number of stars

            #Source: D. Katz, M. Cropper, J.-M. Desert, et al., 3 November 2006, 'CU6 - Spectroscopic processing - preliminary functional and data-flow analysis', GAIA-C6-TN-OPM-DK-001, issue 3, revision 0, Section 4<br/>
            #Basic : false
            #Scalar: true
            #Unit: []
            """

            return self.__SKY_STARDENSITY_RVSREDUCTIONFACTOR


        __SKY_STARDENSITY_TYPICAL: float  = 150000.0 # [stars deg^-2]

        @property
        def SKY_STARDENSITY_TYPICAL(self):
            r"""
            Typical density of stars with G <= 20.00 mag. This value is a typical value encountered in the Galactic plane and, as such, is more representative than the sky-averaged value

            #Source: J.H.J. de Bruijne, 31 October 2003, 'PDHE load assumptions: properties of the sky', GAIA-JdB-009. See also ESA, 21 May 2013, 'Gaia mission requirements document (MRD)', GAIA-EST-RD-00553, issue 3, revision 1, Section 4.1 and ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 3.6.6, page 176<br/>
            #Basic : true
            #Scalar: true
            #Unit: [stars deg^-2]
            """

            return self.__SKY_STARDENSITY_TYPICAL


        __SKY_STARDENSITY_WORSTCASE: float  = 3.0e+06 # [stars deg^-2]

        @property
        def SKY_STARDENSITY_WORSTCASE(self):
            r"""
            Density of stars with G <= 20.00 mag (worst-case, localised value, e.g., in Baade's window)

            #Source: ESA, 21 May 2013, 'Gaia mission requirements document (MRD)', GAIA-EST-RD-00553, issue 3, revision 1, Section 4.1. See also ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 6.5.2, page 244<br/>
            #Basic : true
            #Scalar: true
            #Unit: [stars deg^-2]
            """

            return self.__SKY_STARDENSITY_WORSTCASE


        __SKY_SURFACEBRIGHTNESS_AVERAGE: float  = 22.5 # [mag arcsec^-2]

        @property
        def SKY_SURFACEBRIGHTNESS_AVERAGE(self):
            r"""
            Johnson V band sky-background surface brightness seen from space; average value

            #Source: ESA, 21 May 2013, 'Gaia mission requirements document (MRD)', GAIA-EST-RD-00553, issue 3, revision 1, Requirement SCI-090. See also ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 6.4.3, page 239<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag arcsec^-2]
            """

            return self.__SKY_SURFACEBRIGHTNESS_AVERAGE


        __SKY_SURFACEBRIGHTNESS_WORSTCASE: float  = 21.0 # [mag arcsec^-2]

        @property
        def SKY_SURFACEBRIGHTNESS_WORSTCASE(self):
            r"""
            Johnson V band sky-background surface brightness seen from space; worst-case value (in the ecliptic)

            #Source: ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 6.4.3, page 239<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag arcsec^-2]
            """

            return self.__SKY_SURFACEBRIGHTNESS_WORSTCASE


        __SOLARPROTON_CATALOGUE_AFCCD: str  = "Nature/SolarProton_Catalogue_AFCCD_001.fits" # [e^-]

        @property
        def SOLARPROTON_CATALOGUE_AFCCD(self):
            r"""
            Catalogue containing CCD images, in units of photo-electrons, of typical solar-proton events for an AF CCD (used in BAM, WFS, SM, and AF; this CCD is also used in BP albeit with a different anti-reflection coating). Normally, outside periods of solar activity (solar flares), the solar-proton flux will be negligibly small. During solar activity (solar flares), solar-proton fluxes vary strongly with time, reaching peak levels of 1E6 particles cm^-2 s^-1 or higher. Note that cosmic rays and solar protons are distinct particles, collectively refered to as prompt-particle events (PPEs). The catalogue contains 12954 events. The structure of the FITS file is as follows: the first FITS-file extension contains a list of events containing event number, number of pixels across-scan in the image, and number of pixels along-scan in the image. The following extensions contain the individual images ('pixel matrices'), in units of photo-electron counts, one image per extension

            #Source: A. Short (ESA), priv. comm., 12 May 2006<br/>
            #Basic : true
            #Scalar: false
            #Unit: [e^-]
            """

            return self.__SOLARPROTON_CATALOGUE_AFCCD


        __SOLARPROTON_CATALOGUE_REDENHANCEDCCD: str  = "Nature/SolarProton_Catalogue_RedEnhancedCCD_001.fits" # [e^-]

        @property
        def SOLARPROTON_CATALOGUE_REDENHANCEDCCD(self):
            r"""
            Catalogue containing CCD images, in units of photo-electrons, of typical solar-proton events for a red-enhanced CCD (used in RP and RVS). Normally, outside periods of solar activity (solar flares), the solar-proton flux will be negligibly small. During solar activity (solar flares), solar-proton fluxes vary strongly with time, reaching peak levels of 1E6 particles cm^-2 s^-1 or higher. Note that cosmic rays and solar protons are distinct particles, collectively refered to as prompt-particle events (PPEs). The catalogue contains 4950 events. The structure of the FITS file is as follows: the first FITS-file extension contains a list of events containing event number, number of pixels across-scan in the image, and number of pixels along-scan in the image. The following extensions contain the individual images ('pixel matrices'), in units of photo-electron counts, one image per extension

            #Source: A. Short (ESA), priv. comm., 1 September 2008<br/>
            #Basic : true
            #Scalar: false
            #Unit: [e^-]
            """

            return self.__SOLARPROTON_CATALOGUE_REDENHANCEDCCD


        __SOLARPROTON_FLUX_L2: float  = 1300.0 # [particles cm^-2 s^-1]

        @property
        def SOLARPROTON_FLUX_L2(self):
            r"""
            Maximum-sustainable solar-proton flux at L2, in units of particles cm^-2 s^-1. Normally, outside periods of solar activity (solar flares), the solar-proton flux will be negligibly small. During solar activity (solar flares), solar-proton fluxes vary strongly with time, reaching peak levels of 1E6 particles cm^-2 s^-1 or higher. Note that cosmic rays and solar protons are distinct particles, collectively refered to as prompt-particle events (PPEs). A PPE rate of 1300 particles cm^-2 s^-1 is assumed to correspond to the operational limit below which Gaia functions nominally and above which the spacecraft will be in AOCS TranSition Mode (TSM). An isotropic prompt-particle event flux N, in units of events cm^-2 s^-1, generates 2 N A / 4 events s^-1 CCD^-1, where A denotes the active-pixel area of the CCD in units of cm^2 (including any reduction as a result of TDI-gating), the factor 2 results from considering 'inflow' through both the illuminated and the non-illuminated faces of the CCD, and the factor 4 results from the 'flat geometry' of the CCD (see J.H.J. de Bruijne, A. Short, 7 September 2005, 'prompt-particle events: from fluxes to count rates', GAIA-JdB-026, issue 1, revision 0)

            #Source: EADS-Astrium, 3 March 2011, 'Dead-time budget', GAIA.ASF.TCN.SAT.00133, issue 5, revision 1<br/>
            #Basic : true
            #Scalar: true
            #Unit: [particles cm^-2 s^-1]
            """

            return self.__SOLARPROTON_FLUX_L2


        __SPACE_TEMPERATURE_L2: float  = 3.0 # [K]

        @property
        def SPACE_TEMPERATURE_L2(self):
            r"""
            Space sink temperature at L2

            #Source: European Cooperation for Space Standards (ECSS), 15 November 2008, 'Space environment standard', ECSS-E-ST-10-04C, Section 6.2.1 (http://www.ecss.nl/)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [K]
            """

            return self.__SPACE_TEMPERATURE_L2


        __SPECIFICGAS_CONSTANT_DRYAIR: float  = 287.060 # [m^2 s^-2 K^-1]

        @property
        def SPECIFICGAS_CONSTANT_DRYAIR(self):
            r"""
            Specific gas constant for dry air, defined as the molar gas constant (:Nature:MolarGas_Constant) divided by the molar mass (which is 0.0289644 kg mol^-1 for the International Standard Atmopshere)

            #Source: F. Kleijer (Netherlands Geodetic Commission, Delft), 1 April 2004, 'Troposphere Modeling and Filtering for Precise GPS Leveling', Publications on Geodesy 56, ISBN 90 6132 284 7 (http://www.ncg.knaw.nl/Publicaties/Geodesy/pdf/56Kleijer.pdf and http://repository.tudelft.nl/view/ir/uuid%3Aea1f0cf0-4e48-421b-b7ae-4ae3e36d1880/)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^2 s^-2 K^-1]
            """

            return self.__SPECIFICGAS_CONSTANT_DRYAIR


        __SPECIFICGAS_CONSTANT_WATERVAPOUR: float  = 461.525 # [m^2 s^-2 K^-1]

        @property
        def SPECIFICGAS_CONSTANT_WATERVAPOUR(self):
            r"""
            Specific gas constant for water vapour, defined as the molar gas constant (:Nature:MolarGas_Constant) divided by the molar mass (which is 0.0180152 kg mol^-1)

            #Source: F. Kleijer (Netherlands Geodetic Commission, Delft), 1 April 2004, 'Troposphere Modeling and Filtering for Precise GPS Leveling', Publications on Geodesy 56, ISBN 90 6132 284 7 (http://www.ncg.knaw.nl/Publicaties/Geodesy/pdf/56Kleijer.pdf and http://repository.tudelft.nl/view/ir/uuid%3Aea1f0cf0-4e48-421b-b7ae-4ae3e36d1880/)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^2 s^-2 K^-1]
            """

            return self.__SPECIFICGAS_CONSTANT_WATERVAPOUR


        __STEFANBOLTZMANN_CONSTANT: float  = 5.6703666e-08 # [W m^-2 K^-4]

        @property
        def STEFANBOLTZMANN_CONSTANT(self):
            r"""
            Stefan-Boltzmann constant. Note: best-measured value equals 5.670367E-8 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))

            #Basic : false
            #Scalar: true
            #Unit: [W m^-2 K^-4]
            """

            return self.__STEFANBOLTZMANN_CONSTANT


        __STERADIAN_DEGREESQUARE: float  = 3282.80635001174 # [deg^2]

        @property
        def STERADIAN_DEGREESQUARE(self):
            r"""
            One steradian in units of square degrees

            #Basic : false
            #Scalar: true
            #Unit: [deg^2]
            """

            return self.__STERADIAN_DEGREESQUARE


        __SUNTOEARTHSYSTEM_MASSRATIO: float  = 328900.559616

        @property
        def SUNTOEARTHSYSTEM_MASSRATIO(self):
            r"""
            Ratio of Sun to Earth-system mass. The planetary mass includes the contribution of its satellite, the Moon

            #Basic : false
            #Scalar: true
            #Unit: []
            """

            return self.__SUNTOEARTHSYSTEM_MASSRATIO


        __SUNTOEARTH_MASSRATIO: float  = 332946.048701

        @property
        def SUNTOEARTH_MASSRATIO(self):
            r"""
            Ratio of Sun to Earth mass. The Earth mass includes the Earth's atmosphere but excludes the mass of the Moon

            #Basic : false
            #Scalar: true
            #Unit: []
            """

            return self.__SUNTOEARTH_MASSRATIO


        __SUNTOERISSYSTEM_MASSRATIO: float  = 1.1910e+08

        @property
        def SUNTOERISSYSTEM_MASSRATIO(self):
            r"""
            Ratio of Sun to Eris-system mass (IAU 2009 CBE value). The 'planetary' mass includes the contribution of its satellite, Dysnomia

            #Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__SUNTOERISSYSTEM_MASSRATIO


        __SUNTOJUPITERSYSTEM_MASSRATIO: float  = 1.0473486440e+03

        @property
        def SUNTOJUPITERSYSTEM_MASSRATIO(self):
            r"""
            Ratio of Sun to Jupiter-system mass (IAU 2009 CBE value). The planetary mass includes the contribution of its satellites

            #Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__SUNTOJUPITERSYSTEM_MASSRATIO


        __SUNTOMARSSYSTEM_MASSRATIO: float  = 3.098703590e+06

        @property
        def SUNTOMARSSYSTEM_MASSRATIO(self):
            r"""
            Ratio of Sun to Mars-system mass (IAU 2009 CBE value). The planetary mass includes the contribution of its satellites

            #Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__SUNTOMARSSYSTEM_MASSRATIO


        __SUNTOMERCURYSYSTEM_MASSRATIO: float  = 6.02360e+06

        @property
        def SUNTOMERCURYSYSTEM_MASSRATIO(self):
            r"""
            Ratio of Sun to Mercury(-system) mass (IAU 2009 CBE value)

            #Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__SUNTOMERCURYSYSTEM_MASSRATIO


        __SUNTONEPTUNESYSTEM_MASSRATIO: float  = 1.9412260e+04

        @property
        def SUNTONEPTUNESYSTEM_MASSRATIO(self):
            r"""
            Ratio of Sun to Neptune-system mass (IAU 2009 CBE value). The planetary mass includes the contribution of its satellites

            #Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__SUNTONEPTUNESYSTEM_MASSRATIO


        __SUNTOPLUTOSYSTEM_MASSRATIO: float  = 1.365660e+08

        @property
        def SUNTOPLUTOSYSTEM_MASSRATIO(self):
            r"""
            Ratio of Sun to Pluto-system mass (IAU 2009 CBE value). The 'planetary' mass includes the contribution of its satellite, Charon

            #Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__SUNTOPLUTOSYSTEM_MASSRATIO


        __SUNTOSATURNSYSTEM_MASSRATIO: float  = 3.49790180e+03

        @property
        def SUNTOSATURNSYSTEM_MASSRATIO(self):
            r"""
            Ratio of Sun to Saturn-system mass (IAU 2009 CBE value). The planetary mass includes the contribution of its satellites

            #Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__SUNTOSATURNSYSTEM_MASSRATIO


        __SUNTOURANUSSYSTEM_MASSRATIO: float  = 2.2902980e+04

        @property
        def SUNTOURANUSSYSTEM_MASSRATIO(self):
            r"""
            Ratio of Sun to Uranus-system mass (IAU 2009 CBE value). The planetary mass includes the contribution of its satellites

            #Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__SUNTOURANUSSYSTEM_MASSRATIO


        __SUNTOVENUSSYSTEM_MASSRATIO: float  = 4.085237190e+05

        @property
        def SUNTOVENUSSYSTEM_MASSRATIO(self):
            r"""
            Ratio of Sun to Venus(-system) mass (IAU 2009 CBE value)

            #Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__SUNTOVENUSSYSTEM_MASSRATIO


        __SUN_ABSOLUTEVMAGNITUDE: float  = 4.812 # [mag]

        @property
        def SUN_ABSOLUTEVMAGNITUDE(self):
            r"""
            Johnson V absolute magnitude of the Sun

            #Basic : false
            #Scalar: true
            #Unit: [mag]
            """

            return self.__SUN_ABSOLUTEVMAGNITUDE


        __SUN_ALPHA: float  = 1.68

        @property
        def SUN_ALPHA(self):
            r"""
            Mixing length parameter \alpha of the Sun

            #Source: L. Girardi, A. Bressan, G. Bertelli, C. Chiosi, 2000, 'Evolutionary tracks and isochrones for low- and intermediate-mass stars; from M = 0.15 to 7 M_Sun and from Z = 0.0004 to 0.03', A&AS, 141, 371<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__SUN_ALPHA


        __SUN_BMINV: float  = 0.678 # [mag]

        @property
        def SUN_BMINV(self):
            r"""
            Johnson B-V colour of the Sun

            #Source: Derived from :Nature:Planck_Constant, :Nature:VelocityOfLight_Constant_Vacuum, :Nature:A0VStar_CalibrationWavelength, :Nature:A0VStar_Spectrum_Nu_002, :Nature:A0VStar_CalibrationFunction_002, :Nature:Sun_Spectrum_Nu_001, :Nature:A0VStar_VMinI, :Nature:A0VStar_VMinR, :Nature:A0VStar_BMinV, :Nature:A0VStar_RMinI, :Nature:FilterTransmissionCurve_JohnsonCousinsB_002, :Nature:FilterTransmissionCurve_JohnsonCousinsV_002, :Nature:FilterTransmissionCurve_JohnsonCousinsR_002, and :Nature:FilterTransmissionCurve_JohnsonCousinsI_002. See also I. Ramirez, et al., 2012, 'The UBV(RI)C Colors of the Sun', Astrophysical Journal (ApJ), 752, 5, J. Holmberg, C. Flynn, L. Portinari, 2006, 'The colours of the Sun', MNRAS, 367, 449, and B.J. Taylor, 1998, 'The colours of the Sun', proceedings of IAU Symposium 189 on 'Fundamental Stellar Properties: The Interaction between Observation and Theory', eds T.R. Bedding, A.J. Booth, J. Davis, Kluwer, Dordrecht, ISBN 0-7923-4651-3, page 83 (1998IAUS..189...83T)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__SUN_BMINV


        __SUN_BOLOMETRICCORRECTIONVMAGNITUDE: float  = -0.072 # [mag]

        @property
        def SUN_BOLOMETRICCORRECTIONVMAGNITUDE(self):
            r"""
            Johnson V bolometric correction of the Sun

            #Basic : false
            #Scalar: true
            #Unit: [mag]
            """

            return self.__SUN_BOLOMETRICCORRECTIONVMAGNITUDE


        __SUN_BOLOMETRICMAGNITUDE: float  = 4.740 # [mag]

        @property
        def SUN_BOLOMETRICMAGNITUDE(self):
            r"""
            Absolute bolometric magnitude of the Sun

            #Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B2 on Recommended Zero Points for the Absolute and Apparent Bolometric Magnitude Scales', arXiv:1510.06262 (see also https://www.iau.org/administration/resolutions/general_assemblies/)<br/>
            #Basic : false
            #Scalar: true
            #Unit: [mag]
            """

            return self.__SUN_BOLOMETRICMAGNITUDE


        __SUN_DIURNALPARALLAX: float  = 8.794143 # [arcsec]

        @property
        def SUN_DIURNALPARALLAX(self):
            r"""
            Solar diurnal parallax

            #Basic : false
            #Scalar: true
            #Unit: [arcsec]
            """

            return self.__SUN_DIURNALPARALLAX


        __SUN_EFFECTIVETEMPERATURE: float  = 5772.0 # [K]

        @property
        def SUN_EFFECTIVETEMPERATURE(self):
            r"""
            Effective (black-body) temperature of the Sun

            #Basic : false
            #Scalar: true
            #Unit: [K]
            """

            return self.__SUN_EFFECTIVETEMPERATURE


        __SUN_EFFECTIVETEMPERATURE_NOMINAL: float  = 5772.0 # [K]

        @property
        def SUN_EFFECTIVETEMPERATURE_NOMINAL(self):
            r"""
            Nominal effective (black-body) temperature of the Sun, in SI units. This nominal value shall be understood as conversion factor only

            #Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [K]
            """

            return self.__SUN_EFFECTIVETEMPERATURE_NOMINAL


        __SUN_ENCOMPASSINGSPHERERADIUS: float  = 6.960e+08 # [m]

        @property
        def SUN_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around the Sun which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__SUN_ENCOMPASSINGSPHERERADIUS


        __SUN_ENERGYFLUX_ASTRONOMICALUNIT: float  = 1361.0 # [W m^-2]

        @property
        def SUN_ENERGYFLUX_ASTRONOMICALUNIT(self):
            r"""
            Energy flux of the Sun at a distance of 1 au (also refered to as solar constant or total solar irradiance). This is the cycle-23-averaged, measured value. Due to orbital modulation, this value changes by \pm 3.4% during a year; this value changes by \pm 0.1% during a solar cycle; during solar maximum, sunspots can induce variations over one solar rotation at the level of a few tenths of a percent

            #Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [W m^-2]
            """

            return self.__SUN_ENERGYFLUX_ASTRONOMICALUNIT


        __SUN_ENERGYFLUX_L2: float  = 1334.0 # [W m^-2]

        @property
        def SUN_ENERGYFLUX_L2(self):
            r"""
            Energy flux of the Sun at L2. Due to orbital modulation, this value changes by \pm 3.4% during a year; this value changes by \pm 0.1% during a solar cycle; during solar maximum, sunspots can induce variations over one solar rotation at the level of a few tenths of a percent

            #Basic : false
            #Scalar: true
            #Unit: [W m^-2]
            """

            return self.__SUN_ENERGYFLUX_L2


        __SUN_ENERGYFLUX_NOMINAL: float  = 1361.0 # [W m^-2]

        @property
        def SUN_ENERGYFLUX_NOMINAL(self):
            r"""
            Nominal energy flux of the Sun at a distance of 1 au (also refered to as solar constant or total solar irradiance), in SI units. This nominal value shall be understood as conversion factor only

            #Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [W m^-2]
            """

            return self.__SUN_ENERGYFLUX_NOMINAL


        __SUN_EQUATORIALRADIUS: float  = 6.956580e+08 # [m]

        @property
        def SUN_EQUATORIALRADIUS(self):
            r"""
            Solar (equatorial) radius (photosphere)

            #Source: M. Haberreiter, W. Schmutz, A.G. Kosovichev, 2008, 'Solving the Discrepancy between the Seismic and Photospheric Solar Radius', Astrophysical Journal (ApJ), 675, L53<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__SUN_EQUATORIALRADIUS


        __SUN_EQUATORIALRADIUS_APPARENT: float  = 959.17 # [arcsec]

        @property
        def SUN_EQUATORIALRADIUS_APPARENT(self):
            r"""
            Solar apparent (equatorial) radius at unit distance

            #Basic : false
            #Scalar: true
            #Unit: [arcsec]
            """

            return self.__SUN_EQUATORIALRADIUS_APPARENT


        __SUN_EQUATORIALRADIUS_NOMINAL: float  = 6.9570e+08 # [m]

        @property
        def SUN_EQUATORIALRADIUS_NOMINAL(self):
            r"""
            Nominal solar (equatorial) radius (photosphere), in SI units. This nominal value shall be understood as conversion factor only

            #Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__SUN_EQUATORIALRADIUS_NOMINAL


        __SUN_FREEPHOTONPATH_MEAN: float  = 0.021 # [m]

        @property
        def SUN_FREEPHOTONPATH_MEAN(self):
            r"""
            Coarse estimate of the solar (mean) value of the mean free photon path (assuming complete ionisation)

            #Source: E.g., R. Kippenhahn, A. Weigert, 1991, 'Stellar structure and evolution' (corrected 2-nd printing), Springer Verlag, Berlin, Section 5, Equation 5.1, page 27<br/>
            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__SUN_FREEPHOTONPATH_MEAN


        __SUN_GM: float  = 1.3271244210789466e+20 # [m^3 s^-2]

        @property
        def SUN_GM(self):
            r"""
            Heliocentric gravitational constant (TCB-compatible value). Note that IAU 2012 Resolution B2 adopted at the XXVIII-th General Assembly of the IAU recommends that this parameter be determined observationally in SI units

            #Basic : false
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__SUN_GM


        __SUN_GM_NOMINAL: float  = 1.32712440e+20 # [m^3 s^-2]

        @property
        def SUN_GM_NOMINAL(self):
            r"""
            Nominal heliocentric gravitational constant, in SI units. This nominal value shall be understood as conversion factor only

            #Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m^3 s^-2]
            """

            return self.__SUN_GM_NOMINAL


        __SUN_LIGHTDEFLECTION_LIMB: float  = 1751293.0 # [10^-6 arcsec]

        @property
        def SUN_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle, for an observer at 1 au from the Sun, of a Solar-limb-grazing light ray due to the spherically symmetric part of the gravitational field of the Sun

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__SUN_LIGHTDEFLECTION_LIMB


        __SUN_LIGHTDEFLECTION_RIGHTANGLES: float  = 4072.0 # [10^-6 arcsec]

        @property
        def SUN_LIGHTDEFLECTION_RIGHTANGLES(self):
            r"""
            Post-Newtonian deflection angle, for an observer at 1 au from the Sun, of a light ray arriving at right angles to the solar direction due to the spherically symmetric part of the gravitational field of the Sun

            #Source: ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 3, page 331; cf. S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__SUN_LIGHTDEFLECTION_RIGHTANGLES


        __SUN_LUMINOSITY: float  = 3.8275e+26 # [W]

        @property
        def SUN_LUMINOSITY(self):
            r"""
            Luminosity of the Sun. This value changes by \pm 0.1% during a solar cycle; during solar maximum, sunspots can induce variations over one solar rotation at the level of a few tenths of a percent

            #Basic : false
            #Scalar: true
            #Unit: [W]
            """

            return self.__SUN_LUMINOSITY


        __SUN_LUMINOSITY_NOMINAL: float  = 3.8280e+26 # [W]

        @property
        def SUN_LUMINOSITY_NOMINAL(self):
            r"""
            Nominal luminosity of the Sun, in SI units. This nominal value shall be understood as conversion factor only

            #Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [W]
            """

            return self.__SUN_LUMINOSITY_NOMINAL


        __SUN_MASS: float  = 1.98848e+30 # [kg]

        @property
        def SUN_MASS(self):
            r"""
            Solar mass

            #Basic : false
            #Scalar: true
            #Unit: [kg]
            """

            return self.__SUN_MASS


        __SUN_MASSDENSITY_MEAN: float  = 1.410 # [g cm^-3]

        @property
        def SUN_MASSDENSITY_MEAN(self):
            r"""
            Mean solar mass density

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__SUN_MASSDENSITY_MEAN


        __SUN_MEANMOLECULARWEIGHT: float  = 0.6092

        @property
        def SUN_MEANMOLECULARWEIGHT(self):
            r"""
            Solar value of the mean molecular weight (assuming complete ionisation)

            #Source: E.g., H. Karttunen, et al., 1987, 'Fundamental Astronomy', Springer Verlag, Berlin, Section 11.2, Equation 11.8, page 245<br/>
            #Basic : false
            #Scalar: true
            #Unit: []
            """

            return self.__SUN_MEANMOLECULARWEIGHT


        __SUN_MEANMOLECULARWEIGHT_PERFREEELECTRON: float  = 1.1651

        @property
        def SUN_MEANMOLECULARWEIGHT_PERFREEELECTRON(self):
            r"""
            Solar value of the mean molecular weight per free electron (assuming complete ionisation)

            #Source: E.g., H. Karttunen, et al., 1987, 'Fundamental Astronomy', Springer Verlag, Berlin, Section 11.2, page 246<br/>
            #Basic : false
            #Scalar: true
            #Unit: []
            """

            return self.__SUN_MEANMOLECULARWEIGHT_PERFREEELECTRON


        __SUN_NORTHROTATIONALPOLE_DECLINATION: float  = 63.87 # [deg]

        @property
        def SUN_NORTHROTATIONALPOLE_DECLINATION(self):
            r"""
            IAU-recommended value for the declination \delta_0 of the north pole of rotation of the Sun. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value shall be used for comparative purposes only

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__SUN_NORTHROTATIONALPOLE_DECLINATION


        __SUN_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE: float  = 0.0 # [deg day^-1]

        @property
        def SUN_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the declination \delta_0 of the north pole of rotation of the Sun. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value shall be used for comparative purposes only

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__SUN_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE


        __SUN_NORTHROTATIONALPOLE_RIGHTASCENSION: float  = 286.13 # [deg]

        @property
        def SUN_NORTHROTATIONALPOLE_RIGHTASCENSION(self):
            r"""
            IAU-recommended value for the right ascension \alpha_0 of the north pole of rotation of the Sun. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value shall be used for comparative purposes only

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__SUN_NORTHROTATIONALPOLE_RIGHTASCENSION


        __SUN_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE: float  = 0.0 # [deg day^-1]

        @property
        def SUN_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the right ascension \alpha_0 of the north pole of rotation of the Sun. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value shall be used for comparative purposes only

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__SUN_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE


        __SUN_ORBITALSEMIMAJORAXIS_EARTHSYSTEM: float  = 455.0 # [km]

        @property
        def SUN_ORBITALSEMIMAJORAXIS_EARTHSYSTEM(self):
            r"""
            Orbital semi-major axis of the Sun's orbit around the (hypothetical) common barycentre of the Earth-system-Sun system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.1, Equation 1.6, page 6<br/>
            #Basic : false
            #Scalar: true
            #Unit: [km]
            """

            return self.__SUN_ORBITALSEMIMAJORAXIS_EARTHSYSTEM


        __SUN_ORBITALSEMIMAJORAXIS_JUPITERSYSTEM: float  = 743154.0 # [km]

        @property
        def SUN_ORBITALSEMIMAJORAXIS_JUPITERSYSTEM(self):
            r"""
            Orbital semi-major axis of the Sun's orbit around the (hypothetical) common barycentre of the Jupiter-system-Sun system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.1, Equation 1.6, page 6<br/>
            #Basic : false
            #Scalar: true
            #Unit: [km]
            """

            return self.__SUN_ORBITALSEMIMAJORAXIS_JUPITERSYSTEM


        __SUN_ORBITALSEMIMAJORAXIS_MARSSYSTEM: float  = 74.0 # [km]

        @property
        def SUN_ORBITALSEMIMAJORAXIS_MARSSYSTEM(self):
            r"""
            Orbital semi-major axis of the Sun's orbit around the (hypothetical) common barycentre of the Mars-system-Sun system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.1, Equation 1.6, page 6<br/>
            #Basic : false
            #Scalar: true
            #Unit: [km]
            """

            return self.__SUN_ORBITALSEMIMAJORAXIS_MARSSYSTEM


        __SUN_ORBITALSEMIMAJORAXIS_MERCURYSYSTEM: float  = 10.0 # [km]

        @property
        def SUN_ORBITALSEMIMAJORAXIS_MERCURYSYSTEM(self):
            r"""
            Orbital semi-major axis of the Sun's orbit around the (hypothetical) common barycentre of the Mercury-Sun system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.1, Equation 1.6, page 6<br/>
            #Basic : false
            #Scalar: true
            #Unit: [km]
            """

            return self.__SUN_ORBITALSEMIMAJORAXIS_MERCURYSYSTEM


        __SUN_ORBITALSEMIMAJORAXIS_NEPTUNESYSTEM: float  = 231730.0 # [km]

        @property
        def SUN_ORBITALSEMIMAJORAXIS_NEPTUNESYSTEM(self):
            r"""
            Orbital semi-major axis of the Sun's orbit around the (hypothetical) common barycentre of the Neptune-system-Sun system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.1, Equation 1.6, page 6<br/>
            #Basic : false
            #Scalar: true
            #Unit: [km]
            """

            return self.__SUN_ORBITALSEMIMAJORAXIS_NEPTUNESYSTEM


        __SUN_ORBITALSEMIMAJORAXIS_PLUTOSYSTEM: float  = 43.0 # [km]

        @property
        def SUN_ORBITALSEMIMAJORAXIS_PLUTOSYSTEM(self):
            r"""
            Orbital semi-major axis of the Sun's orbit around the (hypothetical) common barycentre of the Pluto-system-Sun system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.1, Equation 1.6, page 6<br/>
            #Basic : false
            #Scalar: true
            #Unit: [km]
            """

            return self.__SUN_ORBITALSEMIMAJORAXIS_PLUTOSYSTEM


        __SUN_ORBITALSEMIMAJORAXIS_SATURNSYSTEM: float  = 407863.0 # [km]

        @property
        def SUN_ORBITALSEMIMAJORAXIS_SATURNSYSTEM(self):
            r"""
            Orbital semi-major axis of the Sun's orbit around the (hypothetical) common barycentre of the Saturn-system-Sun system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.1, Equation 1.6, page 6<br/>
            #Basic : false
            #Scalar: true
            #Unit: [km]
            """

            return self.__SUN_ORBITALSEMIMAJORAXIS_SATURNSYSTEM


        __SUN_ORBITALSEMIMAJORAXIS_URANUSSYSTEM: float  = 125340.0 # [km]

        @property
        def SUN_ORBITALSEMIMAJORAXIS_URANUSSYSTEM(self):
            r"""
            Orbital semi-major axis of the Sun's orbit around the (hypothetical) common barycentre of the Uranus-system-Sun system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.1, Equation 1.6, page 6<br/>
            #Basic : false
            #Scalar: true
            #Unit: [km]
            """

            return self.__SUN_ORBITALSEMIMAJORAXIS_URANUSSYSTEM


        __SUN_ORBITALSEMIMAJORAXIS_VENUSSYSTEM: float  = 265.0 # [km]

        @property
        def SUN_ORBITALSEMIMAJORAXIS_VENUSSYSTEM(self):
            r"""
            Orbital semi-major axis of the Sun's orbit around the (hypothetical) common barycentre of the Venus-Sun system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.1, Equation 1.6, page 6<br/>
            #Basic : false
            #Scalar: true
            #Unit: [km]
            """

            return self.__SUN_ORBITALSEMIMAJORAXIS_VENUSSYSTEM


        __SUN_PRIMEMERIDIAN_EPHEMERISPOSITION: float  = 84.176 # [deg]

        @property
        def SUN_PRIMEMERIDIAN_EPHEMERISPOSITION(self):
            r"""
            IAU-recommended value for the ephemeris position of the prime meridian of the Sun. The location of the prime meridian is specified by the angle that is measured along the Sun's equator in an easterly direction with respect to the Sun's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the Sun's equator on the standard equator to the point B where the prime meridian crosses the Sun's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the Sun's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the Sun, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the Sun. If W increases with time, the body has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde. The numerical value shall be used for comparative purposes only

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__SUN_PRIMEMERIDIAN_EPHEMERISPOSITION


        __SUN_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE: float  = 14.1844000 # [deg day^-1]

        @property
        def SUN_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch of 1.5 January 2000 = JD2451545.0 TDB) of the ephemeris position of the prime meridian of the Sun. The location of the prime meridian is specified by the angle that is measured along the Sun's equator in an easterly direction with respect to the Sun's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the Sun's equator on the standard equator to the point B where the prime meridian crosses the Sun's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the Sun's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the Sun, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the Sun. If W increases with time, the body has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde. The numerical value shall be used for comparative purposes only

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__SUN_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE


        __SUN_RMINI: float  = 0.344 # [mag]

        @property
        def SUN_RMINI(self):
            r"""
            Cousins R-I colour of the Sun

            #Source: Derived from :Nature:Planck_Constant, :Nature:VelocityOfLight_Constant_Vacuum, :Nature:A0VStar_CalibrationWavelength, :Nature:A0VStar_Spectrum_Nu_002, :Nature:A0VStar_CalibrationFunction_002, :Nature:Sun_Spectrum_Nu_001, :Nature:A0VStar_VMinI, :Nature:A0VStar_VMinR, :Nature:A0VStar_BMinV, :Nature:A0VStar_RMinI, :Nature:FilterTransmissionCurve_JohnsonCousinsB_002, :Nature:FilterTransmissionCurve_JohnsonCousinsV_002, :Nature:FilterTransmissionCurve_JohnsonCousinsR_002, and :Nature:FilterTransmissionCurve_JohnsonCousinsI_002. See also I. Ramirez, et al., 2012, 'The UBV(RI)C Colors of the Sun', Astrophysical Journal (ApJ), 752, 5, J. Holmberg, C. Flynn, L. Portinari, 2006, 'The colours of the Sun', MNRAS, 367, 449, and B.J. Taylor, 1998, 'The colours of the Sun', proceedings of IAU Symposium 189 on 'Fundamental Stellar Properties: The Interaction between Observation and Theory', eds T.R. Bedding, A.J. Booth, J. Davis, Kluwer, Dordrecht, ISBN 0-7923-4651-3, page 83 (1998IAUS..189...83T)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__SUN_RMINI


        __SUN_RADIATIONPRESSURE_ASTRONOMICALUNIT: float  = 4.540e-06 # [Pa]

        @property
        def SUN_RADIATIONPRESSURE_ASTRONOMICALUNIT(self):
            r"""
            Mean value of the solar-radiation pressure at a distance of 1 au. Due to orbital modulation, this value changes by \pm 3.4% during a year; this value changes by \pm 0.1% during a solar cycle; during solar maximum, sunspots can induce variations over one solar rotation at the level of a few tenths of a percent

            #Basic : false
            #Scalar: true
            #Unit: [Pa]
            """

            return self.__SUN_RADIATIONPRESSURE_ASTRONOMICALUNIT


        __SUN_RADIATIONPRESSURE_L2: float  = 4.450e-06 # [Pa]

        @property
        def SUN_RADIATIONPRESSURE_L2(self):
            r"""
            Mean value of the solar-radiation pressure at L2. Due to orbital modulation, this value changes by \pm 3.4% during a year; this value changes by \pm 0.1% during a solar cycle; during solar maximum, sunspots can induce variations over one solar rotation at the level of a few tenths of a percent

            #Basic : false
            #Scalar: true
            #Unit: [Pa]
            """

            return self.__SUN_RADIATIONPRESSURE_L2


        __SUN_ROSSELANDMEANOPACITY_THOMSONSCATTERING: float  = 0.0344 # [m^2 kg^-1]

        @property
        def SUN_ROSSELANDMEANOPACITY_THOMSONSCATTERING(self):
            r"""
            Solar (mean) value of the Rosseland mean opacity for Thomson free-electron-scattering

            #Source: E.g., R. Kippenhahn, A. Weigert, 1991, 'Stellar structure and evolution' (corrected 2-nd printing), Springer Verlag, Berlin, Section 17, Equation 17.2, page 137<br/>
            #Basic : false
            #Scalar: true
            #Unit: [m^2 kg^-1]
            """

            return self.__SUN_ROSSELANDMEANOPACITY_THOMSONSCATTERING


        __SUN_SPECTRUM_NU: str  = "Nature/Sun_Spectrum_Nu_001.fits"

        @property
        def SUN_SPECTRUM_NU(self):
            r"""
            Spectrum f_{\odot\nu}(\lambda) of the Sun: Kurucz/ATLAS9 spectrum (CDROM 19). First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: Eddington flux (in W m^-2 Hz^-1 steradian^-1). Note that the flux at 115.0 nm was obtained using linear interpolation between the available fluxes at 114.5 and 115.5 nm (0.6745E-15 and 0.8131E-15, respectively). Note that the flux at 1062.0 nm was obtained using linear interpolation between the available fluxes at 1057.5 and 1062.5 nm (0.8892E-05 and 0.8861E-05, respectively)

            #Source: C. Jordi, priv. comm., 17 February 2003; see also http://gaia.am.ub.es/PWG/sun.mod. Note that data beyond the current wavelength limits are available<br/>
            #Basic : true
            #Scalar: false
            #Unit: []
            """

            return self.__SUN_SPECTRUM_NU


        __SUN_SURFACEGRAVITY: float  = 274.2 # [m s^-2]

        @property
        def SUN_SURFACEGRAVITY(self):
            r"""
            Surface gravity of the Sun

            #Basic : false
            #Scalar: true
            #Unit: [m s^-2]
            """

            return self.__SUN_SURFACEGRAVITY


        __SUN_VMAGNITUDE: float  = -26.760 # [mag]

        @property
        def SUN_VMAGNITUDE(self):
            r"""
            Johnson V magnitude of the Sun (observed)

            #Source: M.S. Bessell, F. Castelli, B. Plez, 1998, 'Model atmospheres, broad-band colors, bolometric corrections, and temperature calibrations for O-M stars', A&A, 333, 231, Appendices C and D (erratum 1998, A&A, 337, 321)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__SUN_VMAGNITUDE


        __SUN_VMINI: float  = 0.711 # [mag]

        @property
        def SUN_VMINI(self):
            r"""
            Johnson-Cousins V-I colour of the Sun

            #Source: Derived from :Nature:Planck_Constant, :Nature:VelocityOfLight_Constant_Vacuum, :Nature:A0VStar_CalibrationWavelength, :Nature:A0VStar_Spectrum_Nu_002, :Nature:A0VStar_CalibrationFunction_002, :Nature:Sun_Spectrum_Nu_001, :Nature:A0VStar_VMinI, :Nature:A0VStar_VMinR, :Nature:A0VStar_BMinV, :Nature:A0VStar_RMinI, :Nature:FilterTransmissionCurve_JohnsonCousinsB_002, :Nature:FilterTransmissionCurve_JohnsonCousinsV_002, :Nature:FilterTransmissionCurve_JohnsonCousinsR_002, and :Nature:FilterTransmissionCurve_JohnsonCousinsI_002. See also I. Ramirez, et al., 2012, 'The UBV(RI)C Colors of the Sun', Astrophysical Journal (ApJ), 752, 5, J. Holmberg, C. Flynn, L. Portinari, 2006, 'The colours of the Sun', MNRAS, 367, 449, and B.J. Taylor, 1998, 'The colours of the Sun', proceedings of IAU Symposium 189 on 'Fundamental Stellar Properties: The Interaction between Observation and Theory', eds T.R. Bedding, A.J. Booth, J. Davis, Kluwer, Dordrecht, ISBN 0-7923-4651-3, page 83 (1998IAUS..189...83T)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__SUN_VMINI


        __SUN_VMINR: float  = 0.367 # [mag]

        @property
        def SUN_VMINR(self):
            r"""
            Johnson-Cousins V-R colour of the Sun

            #Source: Derived from :Nature:Planck_Constant, :Nature:VelocityOfLight_Constant_Vacuum, :Nature:A0VStar_CalibrationWavelength, :Nature:A0VStar_Spectrum_Nu_002, :Nature:A0VStar_CalibrationFunction_002, :Nature:Sun_Spectrum_Nu_001, :Nature:A0VStar_VMinI, :Nature:A0VStar_VMinR, :Nature:A0VStar_BMinV, :Nature:A0VStar_RMinI, :Nature:FilterTransmissionCurve_JohnsonCousinsB_002, :Nature:FilterTransmissionCurve_JohnsonCousinsV_002, :Nature:FilterTransmissionCurve_JohnsonCousinsR_002, and :Nature:FilterTransmissionCurve_JohnsonCousinsI_002. See also I. Ramirez, et al., 2012, 'The UBV(RI)C Colors of the Sun', Astrophysical Journal (ApJ), 752, 5, J. Holmberg, C. Flynn, L. Portinari, 2006, 'The colours of the Sun', MNRAS, 367, 449, and B.J. Taylor, 1998, 'The colours of the Sun', proceedings of IAU Symposium 189 on 'Fundamental Stellar Properties: The Interaction between Observation and Theory', eds T.R. Bedding, A.J. Booth, J. Davis, Kluwer, Dordrecht, ISBN 0-7923-4651-3, page 83 (1998IAUS..189...83T)<br/>
            #Basic : false
            #Scalar: true
            #Unit: [mag]
            """

            return self.__SUN_VMINR


        __SUN_VELOCITYLSR: float  = 18.04 # [km s^-1]

        @property
        def SUN_VELOCITYLSR(self):
            r"""
            Velocity of the Sun with respect to the local standard of rest (LSR)

            #Basic : false
            #Scalar: true
            #Unit: [km s^-1]
            """

            return self.__SUN_VELOCITYLSR


        __SUN_VELOCITYLSR_U: float  = 11.10 # [km s^-1]

        @property
        def SUN_VELOCITYLSR_U(self):
            r"""
            Velocity of the Sun with respect to the local standard of rest (LSR); U-component, i.e., radially inwards (towards the Galactic centre). The random error is +0.69 and -0.75 km s^-1; the systematic error is 1.0 km s^-1

            #Source: R. Schoenrich, J.J. Binney, W. Dehnen, 1 April 2010, 'Local kinematics and the local standard of rest', MNRAS, 403, 1829-1833 (2010MNRAS.403.1829S)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [km s^-1]
            """

            return self.__SUN_VELOCITYLSR_U


        __SUN_VELOCITYLSR_V: float  = 12.24 # [km s^-1]

        @property
        def SUN_VELOCITYLSR_V(self):
            r"""
            Velocity of the Sun with respect to the local standard of rest (LSR); V-component, i.e., in the direction of Galactic rotation. The random error is +0.47 and -0.47 km s^-1; the systematic error is 2.0 km s^-1

            #Source: R. Schoenrich, J.J. Binney, W. Dehnen, 1 April 2010, 'Local kinematics and the local standard of rest', MNRAS, 403, 1829-1833 (2010MNRAS.403.1829S)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [km s^-1]
            """

            return self.__SUN_VELOCITYLSR_V


        __SUN_VELOCITYLSR_W: float  = 7.25 # [km s^-1]

        @property
        def SUN_VELOCITYLSR_W(self):
            r"""
            Velocity of the Sun with respect to the local standard of rest (LSR); W-component, i.e., vertically upwards (towards the north Galactic pole). The random error is +0.37 and -0.36 km s^-1; the systematic error is 0.5 km s^-1

            #Source: R. Schoenrich, J.J. Binney, W. Dehnen, 1 April 2010, 'Local kinematics and the local standard of rest', MNRAS, 403, 1829-1833 (2010MNRAS.403.1829S)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [km s^-1]
            """

            return self.__SUN_VELOCITYLSR_W


        __SUN_X: float  = 0.7166

        @property
        def SUN_X(self):
            r"""
            Hydrogen abundance of the Sun by mass

            #Basic : false
            #Scalar: true
            #Unit: []
            """

            return self.__SUN_X


        __SUN_Y: float  = 0.2659

        @property
        def SUN_Y(self):
            r"""
            Helium abundance of the Sun by mass

            #Source: N. Grevesse, A. Noels, 1993, in 'Association Vaudoise des chercheurs en physique: la formation des elements chimiques', eds B. Hauck, S. Plantani, D. Raboud and N. Grevesse, A. Noels, 1993, in 'Origin and evolution of the elements', eds N. Prantzos, E. Vangioni-Flam, M. Casse, Cambridge University Press, Cambridge, page 15<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__SUN_Y


        __SUN_Z: float  = 0.0175

        @property
        def SUN_Z(self):
            r"""
            Metal abundance of the Sun by mass

            #Source: N. Grevesse, A. Noels, 1993, in 'Association Vaudoise des chercheurs en physique: la formation des elements chimiques', eds B. Hauck, S. Plantani, D. Raboud and N. Grevesse, A. Noels, 1993, in 'Origin and evolution of the elements', eds N. Prantzos, E. Vangioni-Flam, M. Casse, Cambridge University Press, Cambridge, page 15<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__SUN_Z


        __TAIMINUTC_CONSTANT = [ [2441317.5,  2441499.5,  10.0],  [2441499.5,  2441683.5,  11.0],  [2441683.5,  2442048.5,  12.0],  [2442048.5,  2442413.5,  13.0],  [2442413.5,  2442778.5,  14.0],  [2442778.5,  2443144.5,  15.0],  [2443144.5,  2443509.5,  16.0],  [2443509.5,  2443874.5,  17.0],  [2443874.5,  2444239.5,  18.0],  [2444239.5,  2444786.5,  19.0],  [2444786.5,  2445151.5,  20.0],  [2445151.5,  2445516.5,  21.0],  [2445516.5,  2446247.5,  22.0],  [2446247.5,  2447161.5,  23.0],  [2447161.5,  2447892.5,  24.0],  [2447892.5,  2448257.5,  25.0],  [2448257.5,  2448804.5,  26.0],  [2448804.5,  2449169.5,  27.0],  [2449169.5,  2449534.5,  28.0],  [2449534.5,  2450083.5,  29.0],  [2450083.5,  2450630.5,  30.0],  [2450630.5,  2451179.5,  31.0],  [2451179.5,  2453736.5,  32.0],  [2453736.5,  2454832.5,  33.0],  [2454832.5,  2456109.5,  34.0],  [2456109.5,  2457204.5,  35.0],  [2457204.5,  2457754.5,  36.0],  [2457754.5,  2500000.5,  37.0] ] # [s]

        @property
        def TAIMINUTC_CONSTANT(self):
            r"""
            The difference of TAI and UTC as function of time, since 1972, in units of seconds. UTC differs from TAI by an integral number of seconds ('leap seconds'), in such a way that UT1-UTC stays smaller than 0.9 s in absolute value. The decision to introduce a leap second in UTC to meet this condition is the responsability of the IERS. Announcements are made through IERS Bulletin C (https://hpiers.obspm.fr/eoppc/bul/bulc/). Positive leap seconds make the difference TAI-UTC more positive. The vector elements of this parameter define triplets {STARTTIME, ENDTIME, LEAPSECOND} where STARTTIME denotes the start time of the validity interval (in JD UTC), ENDTIME denotes the end time of the validity interval (in JD UTC), and LEAPSECOND denotes the difference TAI - UTC in units of seconds. Note that the ENDTIME of the final triplet (JD2500000.5 UTC) indicates that the end time of the current validity interval is undefined

            #Source: IERS Bulletin C (https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat; see also http://maia.usno.navy.mil/ser7/tai-utc.dat)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [s]
            """

            return self.__TAIMINUTC_CONSTANT


        __TTMINTAI_CONSTANT_NOMINAL: float  = 32.184 # [s]

        @property
        def TTMINTAI_CONSTANT_NOMINAL(self):
            r"""
            Nominal value of the (constant) offset between TAI and TT: TT(TAI) = TAI + 32.184 s. This offset was chosen to give continuity with the previously used, but now obsolete, Ephemeris Time (see P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7). Deviations of this value, which are attributable to physical defects of atomic time standards, are probably between the limits \pm 10 \mus

            #Source: ESA, 1997, 'The Hipparcos and Tycho Catalogues', Volume 1, page 23<br/>
            #Basic : true
            #Scalar: true
            #Unit: [s]
            """

            return self.__TTMINTAI_CONSTANT_NOMINAL


        __TEMPERATURE_CONSTANT: float  = 5039.8 # [K]

        @property
        def TEMPERATURE_CONSTANT(self):
            r"""
            Temperature constant encountered in physics of stellar atmospheres

            #Source: E.g., E. Bohm-Vitense, 1997, 'Introduction to stellar astrophysics; Volume 2: stellar atmospheres', Cambridge University Press, page 74<br/>
            #Basic : false
            #Scalar: true
            #Unit: [K]
            """

            return self.__TEMPERATURE_CONSTANT


        __TRANSFORMATIONMATRIX_ECLIPTICTOICRS = [ 0.9999999999999639,  2.465125329E-7,  -1.068762105E-7,  -2.686837421E-7,  0.9174821334228226,  -0.3977769913529863,  0E-16,  0.3977769913530006,  0.9174821334228557 ]

        @property
        def TRANSFORMATIONMATRIX_ECLIPTICTOICRS(self):
            r"""
            Transformation matrix which transforms the unit-direction vector r_ecl, expressed in ecliptic coordinates, into the unit-direction vector r_equ in equatorial coordinates (ICRS): r_equ = A_K times r_ecl (see also ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, Section 1.5.3, inverse of Equation 1.5.12). Note that the ICRS origin is shifted in the equatorial plane from \Gamma by \phi = 0.05542 arcsec, positive from \Gamma to the ICRS origin (see J. Chapront, M. Chapront-Touze, G. Francou, 2002, 'A new determination of lunar orbital parameters, precession constant, and tidal acceleration from LLR measurements', A&A, 387, 700). The ICRS has an unambiguous definition with an origin in the ICRF equator defined by the realisation of the ICRF. The ecliptic system is less well-defined, potentially depending on additional conventions in dynamical theories. The transformation quantified here corresponds to the inertial mean ecliptic with obliquity (see parameter :Nature:ObliquityOfEcliptic_J2000) and \Gamma defined by reference to the ICRS equator (other possibilities include the mean equator for J2000 or one of the JPL ephemerides equators). Both the obliquity and the position of \Gamma on the ICRS equator with respect to the ICRS origin have been obtained from LLR measurements. The transformation quantified here has no time dependence (there is no secular variation of the obliquity and no precession): it simply defines the relative situation of the various planes at J2000.0

            #Source: ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, Section 1.5.3, Equation 1.5.7 (generalised matrix A_K)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__TRANSFORMATIONMATRIX_ECLIPTICTOICRS


        __TRANSFORMATIONMATRIX_GALACTICTOICRS = [ -0.0548755604162154,  0.4941094278755837,  -0.8676661490190047,  -0.8734370902348850,  -0.4448296299600112,  -0.1980763734312015,  -0.4838350155487132,  0.7469822444972189,  0.4559837761750669 ]

        @property
        def TRANSFORMATIONMATRIX_GALACTICTOICRS(self):
            r"""
            Transformation matrix which transforms the unit-direction vector r_gal, expressed in Galactic coordinates, into the unit-direction vector r_equ in equatorial coordinates (ICRS): r_equ = A_G times r_gal (see ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, Section 1.5.3, inverse of Equation 1.5.13)

            #Source: ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, Section 1.5.3, Equation 1.5.11 (matrix A_G)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__TRANSFORMATIONMATRIX_GALACTICTOICRS


        __TRANSFORMATIONMATRIX_ICRSTOECLIPTIC = [ 0.9999999999999639,  -2.686837421E-7,  0E-16,  2.465125329E-7,  0.9174821334228226,  0.3977769913530006,  -1.068762105E-7,  -0.3977769913529863,  0.9174821334228557 ]

        @property
        def TRANSFORMATIONMATRIX_ICRSTOECLIPTIC(self):
            r"""
            Transformation matrix which transforms the unit-direction vector r_equ, expressed in equatorial coordinates (ICRS), into the unit-direction vector r_ecl in ecliptic coordinates: r_ecl = A_K^T times r_equ (see also ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, Section 1.5.3, Equation 1.5.12; A_K^T denotes the transpose of matrix A_K). Note that the ICRS origin is shifted in the equatorial plane from \Gamma by \phi = 0.05542 arcsec, positive from \Gamma to the ICRS origin (see J. Chapront, M. Chapront-Touze, G. Francou, 2002, 'A new determination of lunar orbital parameters, precession constant, and tidal acceleration from LLR measurements', A&A, 387, 700). The ICRS has an unambiguous definition with an origin in the ICRF equator defined by the realisation of the ICRF. The ecliptic system is less well-defined, potentially depending on additional conventions in dynamical theories. The transformation quantified here corresponds to the inertial mean ecliptic with obliquity (see parameter :Nature:ObliquityOfEcliptic_J2000) and \Gamma defined by reference to the ICRS equator (other possibilities include the mean equator for J2000 or one of the JPL ephemerides equators). Both the obliquity and the position of \Gamma on the ICRS equator with respect to the ICRS origin have been obtained from LLR measurements. The transformation quantified here has no time dependence (there is no secular variation of the obliquity and no precession): it simply defines the relative situation of the various planes at J2000.0

            #Source: ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, Section 1.5.3, transpose of Equation 1.5.7 (transpose of generalised matrix A_K)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__TRANSFORMATIONMATRIX_ICRSTOECLIPTIC


        __TRANSFORMATIONMATRIX_ICRSTOGALACTIC = [ -0.0548755604162154,  -0.8734370902348850,  -0.4838350155487132,  0.4941094278755837,  -0.4448296299600112,  0.7469822444972189,  -0.8676661490190047,  -0.1980763734312015,  0.4559837761750669 ]

        @property
        def TRANSFORMATIONMATRIX_ICRSTOGALACTIC(self):
            r"""
            Transformation matrix which transforms the unit-direction vector r_equ, expressed in equatorial coordinates (ICRS), into the unit-direction vector r_gal in Galactic coordinates: r_gal = A_G^T times r_equ (see ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, Section 1.5.3, Equation 1.5.13; A_G^T denotes the transpose of matrix A_G)

            #Source: ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, Section 1.5.3, transpose of Equation 1.5.11 (transpose of matrix A_G)<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__TRANSFORMATIONMATRIX_ICRSTOGALACTIC


        __TROPICALYEAR_J2000DAY: float  = 365.242190402 # [day]

        @property
        def TROPICALYEAR_J2000DAY(self):
            r"""
            Number of days per tropical year

            #Source: J.L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touze, G. Francou, J. Laskar, 1994, 'Numerical expressions for precession formulae and mean elements for the Moon and the planets', A&A, 282, 663 (1994A&A...282..663S)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [day]
            """

            return self.__TROPICALYEAR_J2000DAY


        __URANUSSYSTEM_ASTROMETRICSIGNATURE_10PARSEC: float  = 84.0 # [10^-6 arcsec]

        @property
        def URANUSSYSTEM_ASTROMETRICSIGNATURE_10PARSEC(self):
            r"""
            Astrometric signature of the Sun induced by the Uranus system for an observer located at a distance of 10 pc from the Sun

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.7, Equation 1.22, page 11<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__URANUSSYSTEM_ASTROMETRICSIGNATURE_10PARSEC


        __URANUSSYSTEM_MASS: float  = 8.68217e+25 # [kg]

        @property
        def URANUSSYSTEM_MASS(self):
            r"""
            Uranus-system mass (IAU 2009 CBE value). The planetary mass includes the contribution of its satellites

            #Basic : false
            #Scalar: true
            #Unit: [kg]
            """

            return self.__URANUSSYSTEM_MASS


        __URANUSSYSTEM_ORBITALECCENTRICITY_J2000: float  = 0.04725744

        @property
        def URANUSSYSTEM_ORBITALECCENTRICITY_J2000(self):
            r"""
            Mean orbital eccentricity of Uranus, at the standard epoch J2000.0. The mean orbital eccentricity is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__URANUSSYSTEM_ORBITALECCENTRICITY_J2000


        __URANUSSYSTEM_ORBITALINCLINATION_J2000: float  = 0.77263783 # [deg]

        @property
        def URANUSSYSTEM_ORBITALINCLINATION_J2000(self):
            r"""
            Mean orbital inclination of Uranus, at the standard epoch J2000.0. The mean orbital inclination is is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__URANUSSYSTEM_ORBITALINCLINATION_J2000


        __URANUSSYSTEM_ORBITALPERIOD: float  = 84.016846 # [yr]

        @property
        def URANUSSYSTEM_ORBITALPERIOD(self):
            r"""
            Sidereal orbital period

            #Source: Values derived from the mean longitude rates in Table 5.8.1 in P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, page 316; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [yr]
            """

            return self.__URANUSSYSTEM_ORBITALPERIOD


        __URANUSSYSTEM_ORBITALSEMIMAJORAXIS_J2000: float  = 19.18916464 # [au]

        @property
        def URANUSSYSTEM_ORBITALSEMIMAJORAXIS_J2000(self):
            r"""
            Mean orbital semi-major axis of Uranus, at the standard epoch J2000.0. The mean orbital semi-major axis is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: [au]
            """

            return self.__URANUSSYSTEM_ORBITALSEMIMAJORAXIS_J2000


        __URANUSSYSTEM_RADIALVELOCITYSIGNATURE: float  = 0.3 # [m s^-1]

        @property
        def URANUSSYSTEM_RADIALVELOCITYSIGNATURE(self):
            r"""
            Radial-velocity amplitude of the Sun induced by the Uranus system for 'an edge-on observer' (i.e., an observer in the orbital plane of the Uranus system)

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.4, Equation 1.18, page 9<br/>
            #Basic : false
            #Scalar: true
            #Unit: [m s^-1]
            """

            return self.__URANUSSYSTEM_RADIALVELOCITYSIGNATURE


        __URANUS_ENCOMPASSINGSPHERERADIUS: float  = 2.55588e+07 # [m]

        @property
        def URANUS_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around Uranus which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__URANUS_ENCOMPASSINGSPHERERADIUS


        __URANUS_EQUATORIALRADIUS: float  = 2.55588e+07 # [m]

        @property
        def URANUS_EQUATORIALRADIUS(self):
            r"""
            Equatorial radius of Uranus

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__URANUS_EQUATORIALRADIUS


        __URANUS_FLATTENING: float  = 2.292730e-02

        @property
        def URANUS_FLATTENING(self):
            r"""
            Geometrical flattening factor f of Uranus (f = (a-b)/a)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__URANUS_FLATTENING


        __URANUS_FLUXREDUCTION_MAXIMUM: float  = 0.133 # [%]

        @property
        def URANUS_FLUXREDUCTION_MAXIMUM(self):
            r"""
            Maximum reduction of the solar flux for an observer external to the solar system during a transit of Uranus

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.2, Equation 2.4, page 14<br/>
            #Basic : false
            #Scalar: true
            #Unit: [%]
            """

            return self.__URANUS_FLUXREDUCTION_MAXIMUM


        __URANUS_GEOMETRICALBEDO: float  = 0.51

        @property
        def URANUS_GEOMETRICALBEDO(self):
            r"""
            Geometric albedo of Uranus (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__URANUS_GEOMETRICALBEDO


        __URANUS_JSUB2: float  = 0.003516

        @property
        def URANUS_JSUB2(self):
            r"""
            Dynamical form-factor of Uranus (oblateness or Stokes' second-degree zonal harmonic of the potential)

            #Source: P.R. Weissman, L.-A. McFadden, T.V. Johnson (eds.), 1999, 'Encyclopedia of the Solar System (first edition)', Academic Press, page 342<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__URANUS_JSUB2


        __URANUS_LIGHTDEFLECTION_LIMB: float  = 2097.0 # [10^-6 arcsec]

        @property
        def URANUS_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of Uranus

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__URANUS_LIGHTDEFLECTION_LIMB


        __URANUS_MASS: float  = 8.681030e+25 # [kg]

        @property
        def URANUS_MASS(self):
            r"""
            Mass of Uranus (do not use for high-precision (orbit) calculations)

            #Source: R.A. Jacobson, 2007, 'The gravity field of the Uranian system and the orbits of the Uranian satellites and rings', BAAS, 39, 453; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [kg]
            """

            return self.__URANUS_MASS


        __URANUS_MASSDENSITY_MEAN: float  = 1.270 # [g cm^-3]

        @property
        def URANUS_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of Uranus

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__URANUS_MASSDENSITY_MEAN


        __URANUS_NORTHROTATIONALPOLE_DECLINATION: float  = -15.175 # [deg]

        @property
        def URANUS_NORTHROTATIONALPOLE_DECLINATION(self):
            r"""
            IAU-recommended value for the declination \delta_0 of the north pole of rotation of Uranus. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__URANUS_NORTHROTATIONALPOLE_DECLINATION


        __URANUS_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE: float  = 0.0 # [deg day^-1]

        @property
        def URANUS_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the declination \delta_0 of the north pole of rotation of Uranus. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__URANUS_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE


        __URANUS_NORTHROTATIONALPOLE_RIGHTASCENSION: float  = 257.311 # [deg]

        @property
        def URANUS_NORTHROTATIONALPOLE_RIGHTASCENSION(self):
            r"""
            IAU-recommended value for the right ascension \alpha_0 of the north pole of rotation of Uranus. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__URANUS_NORTHROTATIONALPOLE_RIGHTASCENSION


        __URANUS_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE: float  = 0.0 # [deg day^-1]

        @property
        def URANUS_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the right ascension \alpha_0 of the north pole of rotation of Uranus. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__URANUS_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE


        __URANUS_POLARRADIUS: float  = 2.49728e+07 # [m]

        @property
        def URANUS_POLARRADIUS(self):
            r"""
            Polar radius of Uranus

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__URANUS_POLARRADIUS


        __URANUS_PRIMEMERIDIAN_EPHEMERISPOSITION: float  = 203.81 # [deg]

        @property
        def URANUS_PRIMEMERIDIAN_EPHEMERISPOSITION(self):
            r"""
            IAU-recommended value for the ephemeris position of the prime meridian of Uranus. The prime meridian refers to the rotation of the magnetic field System III. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__URANUS_PRIMEMERIDIAN_EPHEMERISPOSITION


        __URANUS_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE: float  = -501.1600928 # [deg day^-1]

        @property
        def URANUS_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch of 1.5 January 2000 = JD2451545.0 TDB) of the ephemeris position of the prime meridian of Uranus. The prime meridian refers to the rotation of the magnetic field System III. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__URANUS_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE


        __URANUS_TRANSITPROBABILITY: float  = 0.025 # [%]

        @property
        def URANUS_TRANSITPROBABILITY(self):
            r"""
            Geometric transit probability (Uranus transiting the Sun) for an observer external to the solar system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.1, Equation 2.2, page 14<br/>
            #Basic : false
            #Scalar: true
            #Unit: [%]
            """

            return self.__URANUS_TRANSITPROBABILITY


        __URANUS_TRANSITTIME_MAXIMUM: float  = 2.45 # [day]

        @property
        def URANUS_TRANSITTIME_MAXIMUM(self):
            r"""
            Maximum transit time of Uranus (transiting the Sun) for an observer external to the solar system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.3, Equation 2.5, page 15<br/>
            #Basic : false
            #Scalar: true
            #Unit: [day]
            """

            return self.__URANUS_TRANSITTIME_MAXIMUM


        __URANUS_VONEZEROMAGNITUDE: float  = -7.19 # [mag]

        @property
        def URANUS_VONEZEROMAGNITUDE(self):
            r"""
            V(1,0) magnitude of Uranus (i.e., the visual magnitude of the planet reduced to a distance of 1 au from both the Sun and Earth and phase angle zero). This parameter is also refered to as absolute magnitude in planetary sciences

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__URANUS_VONEZEROMAGNITUDE


        __URANUS_VOLUMETRICRADIUS: float  = 2.53620e+07 # [m]

        @property
        def URANUS_VOLUMETRICRADIUS(self):
            r"""
            Mean volumetric radius of Uranus

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__URANUS_VOLUMETRICRADIUS


        __VELOCITYOFLIGHT_CONSTANT_VACUUM: float  = 299792458.0 # [m s^-1]

        @property
        def VELOCITYOFLIGHT_CONSTANT_VACUUM(self):
            r"""
            Velocity of light in vacuum (defining constant)

            #Source: P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0). See also the IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m s^-1]
            """

            return self.__VELOCITYOFLIGHT_CONSTANT_VACUUM


        __VENUSSYSTEM_ASTROMETRICSIGNATURE_10PARSEC: float  = 0.177 # [10^-6 arcsec]

        @property
        def VENUSSYSTEM_ASTROMETRICSIGNATURE_10PARSEC(self):
            r"""
            Astrometric signature of the Sun induced by Venus for an observer located at a distance of 10 pc from the Sun

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.7, Equation 1.22, page 11<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__VENUSSYSTEM_ASTROMETRICSIGNATURE_10PARSEC


        __VENUSSYSTEM_MASS: float  = 4.86747e+24 # [kg]

        @property
        def VENUSSYSTEM_MASS(self):
            r"""
            Venus(-system) mass (IAU 2009 CBE value)

            #Basic : false
            #Scalar: true
            #Unit: [kg]
            """

            return self.__VENUSSYSTEM_MASS


        __VENUSSYSTEM_ORBITALECCENTRICITY_J2000: float  = 0.00677672

        @property
        def VENUSSYSTEM_ORBITALECCENTRICITY_J2000(self):
            r"""
            Mean orbital eccentricity of Venus, at the standard epoch J2000.0. The mean orbital eccentricity is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__VENUSSYSTEM_ORBITALECCENTRICITY_J2000


        __VENUSSYSTEM_ORBITALINCLINATION_J2000: float  = 3.39467605 # [deg]

        @property
        def VENUSSYSTEM_ORBITALINCLINATION_J2000(self):
            r"""
            Mean orbital inclination of Venus, at the standard epoch J2000.0. The mean orbital inclination is is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__VENUSSYSTEM_ORBITALINCLINATION_J2000


        __VENUSSYSTEM_ORBITALPERIOD: float  = 0.61519726 # [yr]

        @property
        def VENUSSYSTEM_ORBITALPERIOD(self):
            r"""
            Sidereal orbital period

            #Source: Values derived from the mean longitude rates in Table 5.8.1 in P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, page 316; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [yr]
            """

            return self.__VENUSSYSTEM_ORBITALPERIOD


        __VENUSSYSTEM_ORBITALSEMIMAJORAXIS_J2000: float  = 0.72333566 # [au]

        @property
        def VENUSSYSTEM_ORBITALSEMIMAJORAXIS_J2000(self):
            r"""
            Mean orbital semi-major axis of Venus, at the standard epoch J2000.0. The mean orbital semi-major axis is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document

            #Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem<br/>
            #Basic : true
            #Scalar: true
            #Unit: [au]
            """

            return self.__VENUSSYSTEM_ORBITALSEMIMAJORAXIS_J2000


        __VENUSSYSTEM_RADIALVELOCITYSIGNATURE: float  = 0.086 # [m s^-1]

        @property
        def VENUSSYSTEM_RADIALVELOCITYSIGNATURE(self):
            r"""
            Radial-velocity amplitude of the Sun induced by Venus for 'an edge-on observer' (i.e., an observer in the orbital plane of the Venus)

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.4, Equation 1.18, page 9<br/>
            #Basic : false
            #Scalar: true
            #Unit: [m s^-1]
            """

            return self.__VENUSSYSTEM_RADIALVELOCITYSIGNATURE


        __VENUS_ENCOMPASSINGSPHERERADIUS: float  = 6.051800e+06 # [m]

        @property
        def VENUS_ENCOMPASSINGSPHERERADIUS(self):
            r"""
            Radius of the smallest hypothetical sphere around Venus which would encompass the body (this is a low-accuracy parameter used in the relativistic model)

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__VENUS_ENCOMPASSINGSPHERERADIUS


        __VENUS_EQUATORIALRADIUS: float  = 6.05180e+06 # [m]

        @property
        def VENUS_EQUATORIALRADIUS(self):
            r"""
            Equatorial radius of Venus

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__VENUS_EQUATORIALRADIUS


        __VENUS_FLATTENING: float  = 0.0

        @property
        def VENUS_FLATTENING(self):
            r"""
            Geometrical flattening factor f of Venus (f = (a-b)/a)

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__VENUS_FLATTENING


        __VENUS_FLUXREDUCTION_MAXIMUM: float  = 0.008 # [%]

        @property
        def VENUS_FLUXREDUCTION_MAXIMUM(self):
            r"""
            Maximum reduction of the solar flux for an observer external to the solar system during a transit of Venus

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.2, Equation 2.4, page 14<br/>
            #Basic : false
            #Scalar: true
            #Unit: [%]
            """

            return self.__VENUS_FLUXREDUCTION_MAXIMUM


        __VENUS_GEOMETRICALBEDO: float  = 0.65

        @property
        def VENUS_GEOMETRICALBEDO(self):
            r"""
            Geometric albedo of Venus (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__VENUS_GEOMETRICALBEDO


        __VENUS_JSUB2: float  = 2.70e-05

        @property
        def VENUS_JSUB2(self):
            r"""
            Dynamical form-factor of Venus (oblateness or Stokes' second-degree zonal harmonic of the potential)

            #Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8<br/>
            #Basic : true
            #Scalar: true
            #Unit: []
            """

            return self.__VENUS_JSUB2


        __VENUS_LIGHTDEFLECTION_LIMB: float  = 493.0 # [10^-6 arcsec]

        @property
        def VENUS_LIGHTDEFLECTION_LIMB(self):
            r"""
            Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of Venus

            #Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0<br/>
            #Basic : false
            #Scalar: true
            #Unit: [10^-6 arcsec]
            """

            return self.__VENUS_LIGHTDEFLECTION_LIMB


        __VENUS_MASS: float  = 4.867320e+24 # [kg]

        @property
        def VENUS_MASS(self):
            r"""
            Mass of Venus (do not use for high-precision (orbit) calculations)

            #Source: A.S. Konopliv, et al., 1999, 'Venus gravity: 180-th degree and order model', Icarus, 139, 3-18; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [kg]
            """

            return self.__VENUS_MASS


        __VENUS_MASSDENSITY_MEAN: float  = 5.243 # [g cm^-3]

        @property
        def VENUS_MASSDENSITY_MEAN(self):
            r"""
            Mean mass density of Venus

            #Basic : false
            #Scalar: true
            #Unit: [g cm^-3]
            """

            return self.__VENUS_MASSDENSITY_MEAN


        __VENUS_NORTHROTATIONALPOLE_DECLINATION: float  = 67.16 # [deg]

        @property
        def VENUS_NORTHROTATIONALPOLE_DECLINATION(self):
            r"""
            IAU-recommended value for the declination \delta_0 of the north pole of rotation of Venus. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__VENUS_NORTHROTATIONALPOLE_DECLINATION


        __VENUS_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE: float  = 0.0 # [deg day^-1]

        @property
        def VENUS_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the declination \delta_0 of the north pole of rotation of Venus. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__VENUS_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE


        __VENUS_NORTHROTATIONALPOLE_RIGHTASCENSION: float  = 272.76 # [deg]

        @property
        def VENUS_NORTHROTATIONALPOLE_RIGHTASCENSION(self):
            r"""
            IAU-recommended value for the right ascension \alpha_0 of the north pole of rotation of Venus. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__VENUS_NORTHROTATIONALPOLE_RIGHTASCENSION


        __VENUS_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE: float  = 0.0 # [deg day^-1]

        @property
        def VENUS_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the right ascension \alpha_0 of the north pole of rotation of Venus. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__VENUS_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE


        __VENUS_POLARRADIUS: float  = 6.05180e+06 # [m]

        @property
        def VENUS_POLARRADIUS(self):
            r"""
            Polar radius of Venus

            #Basic : false
            #Scalar: true
            #Unit: [m]
            """

            return self.__VENUS_POLARRADIUS


        __VENUS_PRIMEMERIDIAN_EPHEMERISPOSITION: float  = 160.20 # [deg]

        @property
        def VENUS_PRIMEMERIDIAN_EPHEMERISPOSITION(self):
            r"""
            IAU-recommended value for the ephemeris position of the prime meridian of Venus. The 0-deg meridian is defined by the central peak in the crater Ariadne. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg]
            """

            return self.__VENUS_PRIMEMERIDIAN_EPHEMERISPOSITION


        __VENUS_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE: float  = -1.4813688 # [deg day^-1]

        @property
        def VENUS_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE(self):
            r"""
            IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch of 1.5 January 2000 = JD2451545.0 TDB) of the ephemeris position of the prime meridian of Venus. The 0-deg meridian is defined by the central peak in the crater Ariadne. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg day^-1]
            """

            return self.__VENUS_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE


        __VENUS_TRANSITPROBABILITY: float  = 0.648 # [%]

        @property
        def VENUS_TRANSITPROBABILITY(self):
            r"""
            Geometric transit probability (Venus transiting the Sun) for an observer external to the solar system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.1, Equation 2.2, page 14<br/>
            #Basic : false
            #Scalar: true
            #Unit: [%]
            """

            return self.__VENUS_TRANSITPROBABILITY


        __VENUS_TRANSITTIME_MAXIMUM: float  = 0.46 # [day]

        @property
        def VENUS_TRANSITTIME_MAXIMUM(self):
            r"""
            Maximum transit time of Venus (transiting the Sun) for an observer external to the solar system

            #Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.3, Equation 2.5, page 15<br/>
            #Basic : false
            #Scalar: true
            #Unit: [day]
            """

            return self.__VENUS_TRANSITTIME_MAXIMUM


        __VENUS_VONEZEROMAGNITUDE: float  = -4.47 # [mag]

        @property
        def VENUS_VONEZEROMAGNITUDE(self):
            r"""
            V(1,0) magnitude of Venus (i.e., the visual magnitude of the planet reduced to a distance of 1 au from both the Sun and Earth and phase angle zero). This parameter is also refered to as absolute magnitude in planetary sciences

            #Source: J.L. Hilton, 2005, 'Improving the Visual Magnitudes of the Planets in The Astronomical Almanac. I. Mercury and Venus', AJ, 129, 2902-2906; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [mag]
            """

            return self.__VENUS_VONEZEROMAGNITUDE


        __VENUS_VOLUMETRICRADIUS: float  = 6.05180e+06 # [m]

        @property
        def VENUS_VOLUMETRICRADIUS(self):
            r"""
            Mean volumetric radius of Venus

            #Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m]
            """

            return self.__VENUS_VOLUMETRICRADIUS


        __WIEN_CONSTANT: float  = 2.89777290e-03 # [m K]

        @property
        def WIEN_CONSTANT(self):
            r"""
            Wien's displacement-law constant (for \lambda_max)

            #Source: P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0)<br/>
            #Basic : true
            #Scalar: true
            #Unit: [m K]
            """

            return self.__WIEN_CONSTANT


        __ZEROCELSIUS_KELVIN: float  = 273.15 # [K]

        @property
        def ZEROCELSIUS_KELVIN(self):
            r"""
            Zero degrees Celsius (ice point) expressed in degrees Kelvin

            #Basic : false
            #Scalar: true
            #Unit: [K]
            """

            return self.__ZEROCELSIUS_KELVIN


        __ZEROKELVIN_CELSIUS: float  = -273.15 # [deg C]

        @property
        def ZEROKELVIN_CELSIUS(self):
            r"""
            Zero degrees Kelvin expressed in degrees Celsius. The triple point of water is the only realizable defining fixed point common to the Kelvin Thermodynamic Temperature Scale (KTTS) and the International Temperature Scale of 1990 (ITS-90); the assigned value of the triple point of water on these scales is 273.16 K (0.01 C)

            #Source: The International Temperature Scale of 1990 (ITS-90); http://www.its-90.com/its-90.html<br/>
            #Basic : true
            #Scalar: true
            #Unit: [deg C]
            """

            return self.__ZEROKELVIN_CELSIUS
