//-----------------------------------------------------------------------------
// Copyright (C) 2006-2011 Gaia Data Processing and Analysis Consortium
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//-----------------------------------------------------------------------------
//
//             THIS IS AN AUTOMATICALLY GENERATED FILE - DO NOT EDIT!
//
//    The file has been automatically generated from the contents of the
//    Gaia Parameter Database at the URL
//        https://gaia.esac.esa.int/gpdb/
//    on 2019-09-30T15:15:21.
//
//    Please report any problems arising from the usage of this file to
//    the Gaia Librarian gaia-helpdesk@cosmos.esa.int
//
#ifndef PDBPARAM_HPP
#define PDBPARAM_HPP

//
// Namespace to enclose the contents of the Gaia Parameter Database at
// https://gaia.esac.esa.int/gpdb/
// A hierarchy of nested classes below matches the parameter naming scheme
// detailed in <code><a href="https://dms.cosmos.esa.int/cs/cs/Open/357616">GAIA-JdB-007</a></code><br/>
//
//
// Author: Gaia SOC, ESA/ESTEC
// Version: Live
//

#include <string>

namespace PdbParam {

static const char *const DBVersion = "Live";

class Nature {
public:
    // Johnson B minus V (B-V) magnitude of an unreddened A0V star
    // Source: M.S. Bessell, S. Murphy, 2012, 'Spectrophotometric Libraries, Revised Photonic Passbands, and Zero Points for UBVRI, Hipparcos, and Tycho Photometry', PASP, 124, 140
    // Basic : true
    // Scalar: true
    static const double A0VSTAR_BMINV = -0.004; // [mag]

    // Flux f_{0\lambda} f an unreddened A0V star with V = 0 mag at the wavelength \lambda_0
    // Basic : false
    // Scalar: true
    static const double A0VSTAR_CALIBRATIONFLUX_LAMBDA = 3.62286e-11; // [W m^-2 nm^-1]

    // Flux f_{0\nu} of an unreddened A0V star with V = 0 mag at the wavelength \lambda_0 (see also J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1 and R.C. Bohlin, 2014, 'Hubble Space Telescope CALSPEC Flux Standards: Sirius (and Vega)', Astronomical Journal, Volume 147, 127; spectrum named alpha_lyr_mod_002.fits contained in CALSPEC Calibration Database, http://www.stsci.edu/hst/observatory/crds/calspec.html, last modified April 2015)
    // Source: P. Montegriffo, 26 January 2016, 'External calibration for Gaia DR1 integrated photometry', GAIA-C5-TN-OABO-PMN-009, issue 1, revision 0, depends on parameters :Nature:A0VStar_Spectrum, :Nature:Planck_Constant, :Nature:FilterTransmissionCurve_JohnsonCousinsV_002, and :Nature:A0VStar_VMagnitude
    // Basic : false
    // Scalar: true
    static const double A0VSTAR_CALIBRATIONFLUX_NU = 3.65558e-23; // [W m^-2 Hz^-1]

    // Photon flux N_{0\lambda}(\lambda_0) of an unreddened A0V star with V = 0 mag at the wavelength \lambda_0 (see also J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1). Note that the parameter A0VStar_Spectrum_NumberOfPhotons refers to Pickles' star number 009, and not to the Kurucz Vega spectrum which has been used for flux normalisation/calibration and zero-point definition; the parameter A0VStar_Spectrum_NumberOfPhotons therefore does not precisely have A0VStar_CalibrationFlux photons s^-1 m^-2 nm^-1 at \lambda_0 = A0VStar_CalibrationWavelength
    // Basic : false
    // Scalar: true
    static const double A0VSTAR_CALIBRATIONFLUX_NUMBEROFPHOTONS = 100308492.2552; // [photons s^-1 m^-2 nm^-1]

    // Calibration function S_V(\lambda). This function can, alternatively, be used to define the zero point of the Johnson V magnitude scale by imposing the requirement that, for any stellar photon flux density N_\lambda (in photons s^-1 m^-2 nm^-1 above the Earth's atmosphere) with V = 0 mag, the integral from 470 to 740 nm (the support interval of the Johnson V band) of N_\lambda times S_V(\lambda) equals N_0 photons s^-1 m^-2. The function S_V(\lambda) and the normalisation constant N_0 depend on the value of Planck's constant (parameter :Nature:Planck_Constant), on the definition of the shape of the Johnson V band (parameter :Nature:FilterTransmissionCurve_JohnsonCousinsV_002), on the monochromatic calibration flux f_{0\lambda} (or f_{0\nu}; parameters :Nature:A0VStar_CalibrationFlux_Lambda and :Nature:A0VStar_CalibrationFlux_Nu) at \lambda_0 (parameter :Nature:A0VStar_CalibrationWavelength), and on the spectrum f_{0\nu}(\lambda) of the general unreddened A0V star (parameter :Nature:A0VStar_Spectrum_Nu_002). First column: wavelength \lambda (in nm; from 470.0 to 740.0). Second column: S_V(\lambda)
    // Source: J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1, Appendix D
    // Basic : true
    // Scalar: false
    static const char  *const A0VSTAR_CALIBRATIONFUNCTION() { return "Nature/A0VStar_CalibrationFunction_002.fits"; }

    // Calibration (function) normalisation constant N_0. This constant can be used to define the zero point of the Johnson V magnitude scale by imposing the requirement that, for any stellar photon flux density N_\lambda (in photons s^-1 m^-2 nm^-1 above the Earth's atmosphere) with V = 0 mag, the integral from 470 to 740 nm (the support interval of the Johnson V band) of N_\lambda times S_V(\lambda) equals N_0 photons s^-1 m^-2. The function S_V(\lambda) and the normalisation constant N_0 depend on the value of Planck's constant (parameter :Nature:Planck_Constant), on the definition of the shape of the Johnson V band (parameter :Nature:FilterTransmissionCurve_JohnsonCousinsV_002), on the monochromatic calibration flux f_{0\lambda} (or f_{0\nu}; parameters :Nature:A0VStar_CalibrationFlux_Lambda and :Nature:A0VStar_CalibrationFlux_Nu) at \lambda_0 (parameter :Nature:A0VStar_CalibrationWavelength), and on the spectrum f_{0\nu}(\lambda) of the general unreddened A0V star (parameter :Nature:A0VStar_Spectrum_Nu_002)
    // Source: J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1, Appendix D
    // Basic : true
    // Scalar: true
    static const double A0VSTAR_CALIBRATIONFUNCTION_NORMALISATION = 8630065822.2737; // [photons s^-1 m^-2]

    // Reference wavelength at which the flux f_{0\lambda} of an unreddened A0V star is calibrated (see also J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)
    // Source: V. Straizys, 1992, 'Multicolor stellar photometry', Pachart Publ. House (Tucson), Table 21
    // Basic : true
    // Scalar: true
    static const double A0VSTAR_CALIBRATIONWAVELENGTH = 550.0; // [nm]

    // Weighted mean flux \langle f_{\nu} \rangle over the Johnson V passband (parameter :Nature:FilterTransmissionCurve_JohnsonCousinsV_002) of an unreddened A0V star (parameter :Nature:A0VStar_Spectrum) with V = :Nature:A0VStar_VMagnitude mag
    // Source: P. Montegriffo, 26 January 2016, 'External calibration for Gaia DR1 integrated photometry', GAIA-C5-TN-OABO-PMN-009, issue 1, revision 0, depends on parameters :Nature:A0VStar_Spectrum and :Nature:FilterTransmissionCurve_JohnsonCousinsV_002
    // Basic : true
    // Scalar: true
    static const double A0VSTAR_MEANFLUX_NU = 3.58600e-23; // [W m^-2 Hz^-1]

    // Cousins R minus I (R-I) magnitude of an unreddened A0V star
    // Basic : false
    // Scalar: true
    static const double A0VSTAR_RMINI = -0.001; // [mag]

    // Spectrum f_{0\nu}(\lambda) of an unreddened A0V star: high-fidelity, Kurucz-model Vega spectrum (R = 500) with T_eff = 9400 K, log g = 3.90 dex, and [M/H] = -0.5 dex. The Kurucz model has been scaled to fit STIS data (over the interval 554.5-557.0 nm) by a factor 0.994242. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: Eddington flux (in W m^-2 Hz^-1 steradian^-1). Note that the flux at 115.0 nm was obtained using linear interpolation between the available fluxes at 114.9721 and 115.0873 nm (2.521153898563741E-008 and 2.420265114233843E-008, respectively). Note that the flux at 1062.0 nm was obtained using linear interpolation between the available fluxes at 1061.1654 and 1062.2293 nm (2.508019694385564E-005 and 2.504881789424158E-005, respectively)
    // Source: R.C. Bohlin, 2014, 'Hubble Space Telescope CALSPEC Flux Standards: Sirius (and Vega)', Astronomical Journal, Volume 147, 127; spectrum named alpha_lyr_mod_002.fits contained in CALSPEC Calibration Database (http://www.stsci.edu/hst/observatory/crds/calspec.html, last modified April 2015)
    // Basic : true
    // Scalar: false
    static const char  *const A0VSTAR_SPECTRUM_NU() { return "Nature/A0VStar_Spectrum_Nu_002.fits"; }

    // Photon flux density N_{\lambda}(\lambda) of an unreddened A0V star (Pickles' star number 009) at V = 0 mag. Note that this unreddened A0V star refers to Pickles' star number 009, and not to the Kurucz Vega spectrum which has been used for flux normalisation and zero-point definition. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)
    // Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863
    // Basic : true
    // Scalar: false
    static const char  *const A0VSTAR_SPECTRUM_NUMBEROFPHOTONS() { return "Nature/A0VStar_Spectrum_NumberOfPhotons_001.fits"; }

    // Johnson V magnitude of Vega
    // Source: R.C. Bohlin, 2007, 'HST Stellar Standards with 1\% Accuracy in Absolute Flux', in 'The Future of Photometric, Spectrophotometric and Polarimetric Standardization', ASP Conference Series, Vol. 364, p.315; see also CALSPEC Calibration Database, http://www.stsci.edu/hst/observatory/crds/calspec.html, last modified April 2015
    // Basic : true
    // Scalar: true
    static const double A0VSTAR_VMAGNITUDE = 0.023; // [mag]

    // Johnson V minus Gaia G (V-G) magnitude of an unreddened A0V star, applicable to any photometric band G
    // Source: Definition of Gaia G band(s)
    // Basic : true
    // Scalar: true
    static const double A0VSTAR_VMING = 0.000; // [mag]

    // Johnson V minus Cousins I (V-I) magnitude of an unreddened A0V star
    // Source: M.S. Bessell, S. Murphy, 2012, 'Spectrophotometric Libraries, Revised Photonic Passbands, and Zero Points for UBVRI, Hipparcos, and Tycho Photometry', PASP, 124, 140
    // Basic : true
    // Scalar: true
    static const double A0VSTAR_VMINI = -0.001; // [mag]

    // Johnson V minus Cousins R (V-R) magnitude of an unreddened A0V star
    // Source: M.S. Bessell, S. Murphy, 2012, 'Spectrophotometric Libraries, Revised Photonic Passbands, and Zero Points for UBVRI, Hipparcos, and Tycho Photometry', PASP, 124, 140
    // Basic : true
    // Scalar: true
    static const double A0VSTAR_VMINR = 0.000; // [mag]

    // Photon flux density N_{\lambda}(\lambda) of an unreddened A3V star (Pickles' star number 011) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)
    // Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863
    // Basic : true
    // Scalar: false
    static const char  *const A3VSTAR_SPECTRUM_NUMBEROFPHOTONS() { return "Nature/A3VStar_Spectrum_NumberOfPhotons_001.fits"; }

    // Photon flux density N_{\lambda}(\lambda) of an unreddened A5V star (Pickles' star number 012) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)
    // Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863
    // Basic : true
    // Scalar: false
    static const char  *const A5VSTAR_SPECTRUM_NUMBEROFPHOTONS() { return "Nature/A5VStar_Spectrum_NumberOfPhotons_001.fits"; }

    // Constant of aberration, nowadays irrelevant as a fundamental constant, at the standard epoch J2000.0. The IAU (1976) System of Astronomical Constants (e.g., T. Lederle, 1980, MitAG, 48, 59, Table 1) lists 20.49552 arcsec
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Equation 3.253-4, page 131
    // Basic : false
    // Scalar: true
    static const double ABERRATION_CONSTANT_J2000 = 20.49122; // [arcsec]

    // One Angstrom expressed in units of nm. Note that 'Angstrom' is a non-SI unit which should not be used
    // Basic : true
    // Scalar: true
    static const double ANGSTROM_NANOMETER = 0.1; // [nm]

    // One arcsecond in units of radians
    // Basic : false
    // Scalar: true
    static const double ARCSECOND_RADIAN = 4.848136811095360e-06; // [rad]

    // Ratio of 1 Ceres to solar mass (IAU 2009 CBE value)
    // Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double ASTEROID1CERESMASS_SOLARMASS = 4.720e-10;

    // Diameter of 1 Ceres
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double ASTEROID1CERES_DIAMETER = 933.; // [km]

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of 1 Ceres
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double ASTEROID1CERES_LIGHTDEFLECTION_LIMB = 1.; // [10^-6 arcsec]

    // Orbital eccentricity of 1 Ceres (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double ASTEROID1CERES_ORBITALECCENTRICITY_B1950 = 0.0780;

    // Orbital period of 1 Ceres (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double ASTEROID1CERES_ORBITALPERIOD_B1950 = 4.607; // [yr]

    // Orbital semi-major axis of 1 Ceres (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double ASTEROID1CERES_ORBITALSEMIMAJORAXIS_B1950 = 2.769; // [au]

    // Perihelion distance of 1 Ceres (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic
    // Basic : false
    // Scalar: true
    static const double ASTEROID1CERES_PERIHELIONDISTANCE_B1950 = 2.553; // [au]

    // Ratio of 2 Pallas to solar mass (IAU 2009 CBE value)
    // Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double ASTEROID2PALLASMASS_SOLARMASS = 1.030e-10;

    // Diameter of 2 Pallas
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double ASTEROID2PALLAS_DIAMETER = 525.; // [km]

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of 2 Pallas
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double ASTEROID2PALLAS_LIGHTDEFLECTION_LIMB = 0.; // [10^-6 arcsec]

    // Orbital eccentricity of 2 Pallas (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double ASTEROID2PALLAS_ORBITALECCENTRICITY_B1950 = 0.2347;

    // Orbital period of 2 Pallas (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double ASTEROID2PALLAS_ORBITALPERIOD_B1950 = 4.611; // [yr]

    // Orbital semi-major axis of 2 Pallas (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double ASTEROID2PALLAS_ORBITALSEMIMAJORAXIS_B1950 = 2.770; // [au]

    // Perihelion distance of 2 Pallas (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic
    // Basic : false
    // Scalar: true
    static const double ASTEROID2PALLAS_PERIHELIONDISTANCE_B1950 = 2.120; // [au]

    // Diameter of 3 Juno
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf). Note that a radius of 120 km is found by E.F. Tedesco, 1989, 'Asteroid magnitudes, UBV colors, and IRAS albedos and diameters', in 'Asteroids II', proceedings of the Conference, Tucson, AZ, 8-11 March 1988, eds R.P. Binzel, T. Gehrels, M.S. Matthews, University of Arizona Press, page 1090 (1989aste.conf.1090T)
    // Basic : true
    // Scalar: true
    static const double ASTEROID3JUNO_DIAMETER = 267.; // [km]

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of 3 Juno
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double ASTEROID3JUNO_LIGHTDEFLECTION_LIMB = 0.; // [10^-6 arcsec]

    // Mass of 3 Juno
    // Source: Value is calculated following J.L. Hilton, 1999, 'US Naval Observatory Ephemerides of the Largest Asteroids', AJ, 117, 1077, who assumes a mean mass density of 3 g cm^-3. A mass of 2.0E19 kg is found on http://nssdc.gsfc.nasa.gov/planetary/factsheet/asteroidfact.html
    // Basic : false
    // Scalar: true
    static const double ASTEROID3JUNO_MASS = 2.990e+19; // [kg]

    // Orbital eccentricity of 3 Juno (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double ASTEROID3JUNO_ORBITALECCENTRICITY_B1950 = 0.0258;

    // Orbital period of 3 Juno (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double ASTEROID3JUNO_ORBITALPERIOD_B1950 = 4.359; // [yr]

    // Orbital semi-major axis of 3 Juno (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double ASTEROID3JUNO_ORBITALSEMIMAJORAXIS_B1950 = 2.668; // [au]

    // Perihelion distance of 3 Juno (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic
    // Basic : false
    // Scalar: true
    static const double ASTEROID3JUNO_PERIHELIONDISTANCE_B1950 = 2.599; // [au]

    // Ratio of 4 Vesta to solar mass (IAU 2009 CBE value)
    // Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double ASTEROID4VESTAMASS_SOLARMASS = 1.350e-10;

    // Diameter of 4 Vesta
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double ASTEROID4VESTA_DIAMETER = 510.; // [km]

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of 4 Vesta
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double ASTEROID4VESTA_LIGHTDEFLECTION_LIMB = 1.; // [10^-6 arcsec]

    // Orbital eccentricity of 4 Vesta (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double ASTEROID4VESTA_ORBITALECCENTRICITY_B1950 = 0.0906;

    // Orbital period of 4 Vesta (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double ASTEROID4VESTA_ORBITALPERIOD_B1950 = 3.629; // [yr]

    // Orbital semi-major axis of 4 Vesta (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 15 'Prominent minor planets or asteroids' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double ASTEROID4VESTA_ORBITALSEMIMAJORAXIS_B1950 = 2.361; // [au]

    // Perihelion distance of 4 Vesta (epoch 1 October 1989). The orbital elements are referenced to the B1950 equinox and ecliptic
    // Basic : false
    // Scalar: true
    static const double ASTEROID4VESTA_PERIHELIONDISTANCE_B1950 = 2.147; // [au]

    // Total mass of the main asteroid belt, in units of solar masses
    // Source: E. Pitjeva, 2003, 'The Dynamic Estimation of the Mass of the Main Asteroid Belt', in 'Physical Properties and Morphology of Small Solar-System Bodies', XXV-th General Assembly of the IAU, Joint Discussion 19, 23 July 2003, Sidney, Australia (2003IAUJD..19E..22P)
    // Basic : true
    // Scalar: true
    static const double ASTEROIDBELTMASS_SOLARMASS = 1.400e-09;

    // The velocity distribution of main-belt asteroids (MBOs) is very roughly Gaussian with zero mean and a standard deviation of 13.0 mas s^-1 across-scan (for a solar-aspect angle of 45 degrees)
    // Source: F. Mignard, 2002, 'Observations of solar-system objects with Gaia. I. Detection of NEOs', A&A, 393, 727, Section 4.4 (2002A&A...393..727M). See also E. Hoeg, F. Arenou, P. Hjorth, U.G. Joergensen, F. Mignard, S. Wolff, 28 February 2003, 'Faint objects and NEOs with Gaia', GAIA-CUO-118, issue 1, revision 0 and F. Mignard, 22 June 2001, 'Observation of main-belt asteroids with Gaia', GAIA-FM-009, issue 1, revision 0. Current value, for a solar-aspect angle of 45 degrees, from F. Mignard, priv. comm., 10 August 2005
    // Basic : true
    // Scalar: true
    static const double ASTEROID_VELOCITYDISTRIBUTIONSTANDARDDEVIATION_AC = 13.0; // [mas s^-1]

    // The velocity distribution of main-belt asteroids (MBOs) is very roughly Gaussian with zero mean and a standard deviation of 7.0 mas s^-1 along-scan (for a solar-aspect angle of 45 degrees)
    // Source: F. Mignard, 2002, 'Observations of solar-system objects with Gaia. I. Detection of NEOs', A&A, 393, 727, Section 4.4 (2002A&A...393..727M). See also E. Hoeg, F. Arenou, P. Hjorth, U.G. Joergensen, F. Mignard, S. Wolff, 28 February 2003, 'Faint objects and NEOs with Gaia', GAIA-CUO-118, issue 1, revision 0 and F. Mignard, 22 June 2001, 'Observation of main-belt asteroids with Gaia', GAIA-FM-009, issue 1, revision 0. Current value, for a solar-aspect angle of 45 degrees, from F. Mignard, priv. comm., 10 August 2005
    // Basic : true
    // Scalar: true
    static const double ASTEROID_VELOCITYDISTRIBUTIONSTANDARDDEVIATION_AL = 7.0; // [mas s^-1]

    // Astronomical unit (au) length. The au is a conventional unit of length and is a defining constant. The numerical value is in agreement with the value adopted in IAU 2009 Resolution B2. The definition applies to all time scales such as TCB, TDB, TCG, TT, etc.
    // Source: IAU, August 2012, 'Re-definition of the astronomical unit of length', IAU 2012 Resolution B2 adopted at the XXVIII-th General Assembly of the IAU
    // Basic : true
    // Scalar: true
    static const double ASTRONOMICALUNIT_METER = 149597870700.; // [m]

    // Atomic mass constant (also known as atomic mass unit [amu]; 1 amu is defined as 1/12-th of the mass of a 12-C atom). Note: best-measured value equals 1.660539040E-27 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))
    // Basic : false
    // Scalar: true
    static const double ATOMICMASS_CONSTANT = 1.6605390404e-27; // [kg]

    // Avogadro's constant
    // Source: P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0)
    // Basic : true
    // Scalar: true
    static const double AVOGADRO_CONSTANT = 6.0221408570e+23; // [mol^-1]

    // Photon flux density N_{\lambda}(\lambda) of an unreddened B0I star (Pickles' star number 114) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)
    // Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863
    // Basic : true
    // Scalar: false
    static const char  *const B0ISTAR_SPECTRUM_NUMBEROFPHOTONS() { return "Nature/B0IStar_Spectrum_NumberOfPhotons_001.fits"; }

    // Photon flux density N_{\lambda}(\lambda) of an unreddened B1V star (Pickles' star number 004) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)
    // Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863
    // Basic : true
    // Scalar: false
    static const char  *const B1VSTAR_SPECTRUM_NUMBEROFPHOTONS() { return "Nature/B1VStar_Spectrum_NumberOfPhotons_001.fits"; }

    // High-resolution photon-flux density N_{\lambda}(\lambda) of an unreddened B1V star at V = 15 mag. The data refer to a high-resolution Kurucz-model spectrum with the following properties: effective temperature T_eff = 25500 K, logarithm of surface gravity log g = 4.0, metallicity [Fe/H] = 0.0, alpha-elements [\alpha/Fe] = 0.0, rotational velocity v sini = 50 km s^-1, micro-turbulence velocity = 2.0 km s^-1, length of convective bubble divided by pressure scale height = 0.50, no convective overshooting, macro-turbulence velocity = 2.0 km s^-1, and resolving power R = \lambda / \delta \lambda = 250,000. First column: wavelength \lambda (in nm; from 830.1673264 to 889.8217922). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1). The 34698 lines have an average wavelength step of 0.00172 nm; the spectrum extent is thus 59.7 nm
    // Source: ESA, 20 June 2005, 'Photon-flux distributions for reference stars', GAIA-EST-TN-00539, issue 1, revision 0, based on D. Katz, priv. comm., 11 May 2005
    // Basic : true
    // Scalar: false
    static const char  *const B1VSTAR_SPECTRUM_NUMBEROFPHOTONSHIGHRESOLUTION() { return "Nature/B1VStar_Spectrum_NumberOfPhotonsHighResolution_001.fits"; }

    // Bohr radius. Note: best-measured value equals 0.52917721067E-10 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))
    // Basic : false
    // Scalar: true
    static const double BOHRRADIUS_CONSTANT = 5.29177210564e-11; // [m]

    // Boltzmann's constant. Note: best-measured value equals 1.38064852E-23 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))
    // Basic : false
    // Scalar: true
    static const double BOLTZMANN_CONSTANT = 1.380648510e-23; // [J K^-1]

    // Solar value of the equivalent width of the first line of the CaII-triplet (3p^{6}3d ^{2}D_{3/2} - 3p^{6}4p ^{2}P_{3/2} transition). A relative intensity of 130 is mentioned in Y. Ralchenko, A.E. Kramida, J. Reader, and the NIST ASD Team, 4 October 2010, 'NIST Atomic Spectra Database', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.nist.gov/pml/data/asd.cfm (Web Version 4.0)
    // Source: C.E. Moore, M.G.J. Minnaert, J. Houtgast, 1966, 'The solar spectrum 2935 AA to 8770 AA', US National Bureau of Standards, Monograph 61
    // Basic : true
    // Scalar: true
    static const double CAIITRIPLET_EQUIVALENTWIDTH_1SUN = 0.146; // [nm]

    // Solar value of the equivalent width of the second line of the CaII-triplet (3p^{6}3d ^{2}D_{5/2} - 3p^{6}4p ^{2}P_{3/2} transition). A relative intensity of 170 is mentioned in Y. Ralchenko, A.E. Kramida, J. Reader, and the NIST ASD Team, 4 October 2010, 'NIST Atomic Spectra Database', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.nist.gov/pml/data/asd.cfm (Web Version 4.0)
    // Source: C.E. Moore, M.G.J. Minnaert, J. Houtgast, 1966, 'The solar spectrum 2935 AA to 8770 AA', US National Bureau of Standards, Monograph 61
    // Basic : true
    // Scalar: true
    static const double CAIITRIPLET_EQUIVALENTWIDTH_2SUN = 0.367; // [nm]

    // Solar value of the equivalent width of the third line of the CaII-triplet (3p^{6}3d ^{2}D_{3/2} - 3p^{6}4p ^{2}P_{1/2} transition). A relative intensity of 160 is mentioned in Y. Ralchenko, A.E. Kramida, J. Reader, and the NIST ASD Team, 4 October 2010, 'NIST Atomic Spectra Database', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.nist.gov/pml/data/asd.cfm (Web Version 4.0)
    // Source: C.E. Moore, M.G.J. Minnaert, J. Houtgast, 1966, 'The solar spectrum 2935 AA to 8770 AA', US National Bureau of Standards, Monograph 61
    // Basic : true
    // Scalar: true
    static const double CAIITRIPLET_EQUIVALENTWIDTH_3SUN = 0.260; // [nm]

    // Oscillator strength of the first line of the CaII-triplet (3p^{6}3d ^{2}D_{3/2} - 3p^{6}4p ^{2}P_{3/2} transition). The estimated accuracy is better than 25%
    // Source: Y. Ralchenko, A.E. Kramida, J. Reader, and the NIST ASD Team, 4 October 2010, 'NIST Atomic Spectra Database', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.nist.gov/pml/data/asd.cfm (Web Version 4.0)
    // Basic : true
    // Scalar: true
    static const double CAIITRIPLET_OSCILLATORSTRENGTH_1 = -1.318;

    // Oscillator strength of the second line of the CaII-triplet (3p^{6}3d ^{2}D_{5/2} - 3p^{6}4p ^{2}P_{3/2} transition). The estimated accuracy is better than 25%
    // Source: Y. Ralchenko, A.E. Kramida, J. Reader, and the NIST ASD Team, 4 October 2010, 'NIST Atomic Spectra Database', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.nist.gov/pml/data/asd.cfm (Web Version 4.0)
    // Basic : true
    // Scalar: true
    static const double CAIITRIPLET_OSCILLATORSTRENGTH_2 = -0.36;

    // Oscillator strength of the third line of the CaII-triplet (3p^{6}3d ^{2}D_{3/2} - 3p^{6}4p ^{2}P_{1/2} transition). The estimated accuracy is better than 25%
    // Source: Y. Ralchenko, A.E. Kramida, J. Reader, and the NIST ASD Team, 4 October 2010, 'NIST Atomic Spectra Database', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.nist.gov/pml/data/asd.cfm (Web Version 4.0)
    // Basic : true
    // Scalar: true
    static const double CAIITRIPLET_OSCILLATORSTRENGTH_3 = -0.622;

    // Rest-wavelength in vacuum of the first line of the CaII-triplet (3p^{6}3d ^{2}D_{3/2} - 3p^{6}4p ^{2}P_{3/2} transition), as calculated from the difference between the energy of the upper and lower level of the transition
    // Source: Y. Ralchenko, A.E. Kramida, J. Reader, and the NIST ASD Team, 4 October 2010, 'NIST Atomic Spectra Database', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.nist.gov/pml/data/asd.cfm (Web Version 4.0)
    // Basic : true
    // Scalar: true
    static const double CAIITRIPLET_WAVELENGTH_1 = 850.035; // [nm]

    // Rest-wavelength in vacuum of the second line of the CaII-triplet (3p^{6}3d ^{2}D_{5/2} - 3p^{6}4p ^{2}P_{3/2} transition), as calculated from the difference between the energy of the upper and lower level of the transition
    // Source: Y. Ralchenko, A.E. Kramida, J. Reader, and the NIST ASD Team, 4 October 2010, 'NIST Atomic Spectra Database', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.nist.gov/pml/data/asd.cfm (Web Version 4.0)
    // Basic : true
    // Scalar: true
    static const double CAIITRIPLET_WAVELENGTH_2 = 854.444; // [nm]

    // Rest-wavelength in vacuum of the third line of the CaII-triplet (3p^{6}3d ^{2}D_{3/2} - 3p^{6}4p ^{2}P_{1/2} transition), as calculated from the difference between the energy of the upper and lower level of the transition
    // Source: Y. Ralchenko, A.E. Kramida, J. Reader, and the NIST ASD Team, 4 October 2010, 'NIST Atomic Spectra Database', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.nist.gov/pml/data/asd.cfm (Web Version 4.0)
    // Basic : true
    // Scalar: true
    static const double CAIITRIPLET_WAVELENGTH_3 = 866.452; // [nm]

    // Radius of the smallest hypothetical sphere around Charon which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double CHARON_ENCOMPASSINGSPHERERADIUS = 6.050e+05; // [m]

    // GM of Charon
    // Source: R.A. Jacobson, 2007, 'Constants used in the PLU017 ephemeris', priv. comm.; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double CHARON_GM = 1.0320e+11; // [m^3 s^-2]

    // Geometric albedo of Charon (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)
    // Source: K. Reinsch, V. Burwitz, and M.C. Festou, 1994, 'Albedo Maps of Pluto and Improved Physical Parameters of the Pluto-Charon System', Icarus, 108, 209-218; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double CHARON_GEOMETRICALBEDO = 0.372;

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of Charon
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double CHARON_LIGHTDEFLECTION_LIMB = 1.; // [10^-6 arcsec]

    // Mass of Charon (do not use for high-precision (orbit) calculations)
    // Basic : false
    // Scalar: true
    static const double CHARON_MASS = 1.4705e+21; // [kg]

    // Mean mass density of Charon
    // Basic : false
    // Scalar: true
    static const double CHARON_MASSDENSITY_MEAN = 1.585; // [g cm^-3]

    // Radius of Charon
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double CHARON_RADIUS = 6.050e+05; // [m]

    // Catalogue containing CCD images, in units of photo-electrons, of typical (galactic) cosmic-ray events for an AF CCD (used in BAM, WFS, SM, and AF; this CCD is also used in BP albeit with a different anti-reflection coating). Cosmic rays will be present as constant background all through the mission. Note that cosmic rays and solar protons are distinct particles, collectively refered to as prompt-particle events (PPEs). The catalogue contains 12389 events. The structure of the FITS file is as follows: the first FITS-file extension contains a list of events containing event number, number of pixels across-scan in the image, and number of pixels along-scan in the image. The following extensions contain the individual images ('pixel matrices'), in units of photo-electron counts, one image per extension
    // Source: A. Short (ESA), priv. comm., 12 May 2006
    // Basic : true
    // Scalar: false
    static const char  *const COSMICRAY_CATALOGUE_AFCCD() { return "Nature/CosmicRay_Catalogue_AFCCD_001.fits"; } // [e^-]

    // Catalogue containing CCD images, in units of photo-electrons, of typical (galactic) cosmic-ray events for a red-enhanced CCD (used in RP and RVS). Cosmic rays will be present as constant background all through the mission. Note that cosmic rays and solar protons are distinct particles, collectively refered to as prompt-particle events (PPEs). The catalogue contains 4718 events. The structure of the FITS file is as follows: the first FITS-file extension contains a list of events containing event number, number of pixels across-scan in the image, and number of pixels along-scan in the image. The following extensions contain the individual images ('pixel matrices'), in units of photo-electron counts, one image per extension
    // Source: A. Short (ESA), priv. comm., 5 May 2004
    // Basic : true
    // Scalar: false
    static const char  *const COSMICRAY_CATALOGUE_REDENHANCEDCCD() { return "Nature/CosmicRay_Catalogue_RedEnhancedCCD_001.fits"; } // [e^-]

    // Typical expected (galactic) cosmic-ray flux at L2, in units of particles cm^-2 s^-1. This flux will be present as constant background all through the mission. Note that cosmic rays and solar protons are distinct particles, collectively refered to as prompt-particle events (PPEs). An isotropic prompt-particle event flux N, in units of events cm^-2 s^-1, generates 2 N A / 4 events s^-1 CCD^-1, where A denotes the active-pixel area of the CCD in units of cm^2 (including any reduction as a result of TDI-gating), the factor 2 results from considering 'inflow' through both the illuminated and the non-illuminated faces of the CCD, and the factor 4 results from the 'flat geometry' of the CCD (see J.H.J. de Bruijne, A. Short, 7 September 2005, 'prompt-particle events: from fluxes to count rates', GAIA-JdB-026, issue 1, revision 0)
    // Source: A. Short (ESA), priv. comm., 20 December 2005
    // Basic : true
    // Scalar: true
    static const double COSMICRAY_FLUX_L2 = 5.; // [particles cm^-2 s^-1]

    class DE405 {
    public:
        // Ratio of 1 Ceres to solar mass (DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double ASTEROID1CERESMASS_SOLARMASS = 4.70e-10;

        // Ratio of 2 Pallas to solar mass (DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double ASTEROID2PALLASMASS_SOLARMASS = 1.00e-10;

        // Ratio of 4 Vesta to solar mass (DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double ASTEROID4VESTAMASS_SOLARMASS = 1.30e-10;

        // Mean mass density of C-class asteroids (DE405 value). In JPL's DE ephemerides, masses for the 300 most massive asteroids (except Ceres, Pallas, and Vesta) are derived using the relation GM = 6.27E-22 Radius^3 \rho, where Radius is the known asteroid radius in km and \rho is the mean mass density in g cm^-3 (GM in this relation is in 'solar-system units', i.e., Sun_GM = Gauss_Constant * Gauss_Constant [au^3/2 day^-1 M_Sun^-1/2]; note that the reference erroneously lists a prefactor 6.27E20; typo confirmed by E.M. Standish, priv. comm., 18 September 2003). The mean density \rho is assumed to be constant within each of the three taxonomic asteroid classes C, S, and M
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048, Section IV
        // Basic : true
        // Scalar: true
        static const double ASTEROID_MASSDENSITY_MEANCLASSC = 1.8; // [g cm^-3]

        // Mean mass density of M-class asteroids (DE405 value). In JPL's DE ephemerides, masses for the 300 most massive asteroids (except Ceres, Pallas, and Vesta) are derived using the relation GM = 6.27E-22 Radius^3 \rho, where Radius is the known asteroid radius in km and \rho is the mean mass density in g cm^-3 (GM in this relation is in 'solar-system units', i.e., Sun_GM = Gauss_Constant * Gauss_Constant [au^3/2 day^-1 M_Sun^-1/2]; note that the reference erroneously lists a prefactor 6.27E20; typo confirmed by E.M. Standish, priv. comm., 18 September 2003). The mean density \rho is assumed to be constant within each of the three taxonomic asteroid classes C, S, and M
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048, Section IV
        // Basic : true
        // Scalar: true
        static const double ASTEROID_MASSDENSITY_MEANCLASSM = 5.0; // [g cm^-3]

        // Mean mass density of S-class asteroids (DE405 value). In JPL's DE ephemerides, masses for the 300 most massive asteroids (except Ceres, Pallas, and Vesta) are derived using the relation GM = 6.27E-22 Radius^3 \rho, where Radius is the known asteroid radius in km and \rho is the mean mass density in g cm^-3 (GM in this relation is in 'solar-system units', i.e., Sun_GM = Gauss_Constant * Gauss_Constant [au^3/2 day^-1 M_Sun^-1/2]; note that the reference erroneously lists a prefactor 6.27E20; typo confirmed by E.M. Standish, priv. comm., 18 September 2003). The mean density \rho is assumed to be constant within each of the three taxonomic asteroid classes C, S, and M
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048, Section IV
        // Basic : true
        // Scalar: true
        static const double ASTEROID_MASSDENSITY_MEANCLASSS = 2.4; // [g cm^-3]

        // Astronomical unit (au) length (TCB-compatible value; DE405 value; see S.A. Klioner, 2008, 'Relativistic scaling of astronomical quantities and the system of astronomical units', A&A, 478, 951-958)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : false
        // Scalar: true
        static const double ASTRONOMICALUNIT_METER = 1.4959787301053391e+11; // [m]

        // Astronomical unit (au) light time (TCB-compatible value; DE405 value; see S.A. Klioner, 2008, 'Relativistic scaling of astronomical quantities and the system of astronomical units', A&A, 478, 951-958)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : false
        // Scalar: true
        static const double ASTRONOMICALUNIT_SECOND = 4.9900479154326796e+02; // [s]

        // Astronomical unit (au) length (TDB-compatible value; DE405 value). Do not use this parameter but use the TCB-compatible value from parameter :Nature:DE405:AstronomicalUnit_Meter instead
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : false
        // Scalar: true
        static const double ASTRONOMICALUNIT_TDBMETER = 1.4959787301053392e+11; // [m (TDB)]

        // Astronomical unit (au) light time (TDB-compatible value; DE405 value). Do not use this parameter but use the TCB-compatible value from parameter :Nature:DE405:AstronomicalUnit_Second instead
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double ASTRONOMICALUNIT_TDBSECOND = 4.9900478380610000e+02; // [s (TDB)]

        // Ratio of Earth to Moon mass (DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double EARTHTOMOON_MASSRATIO = 81.30056;

        // Equatorial radius of the Earth (DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double EARTH_EQUATORIALRADIUS = 6378137.; // [m]

        // Geocentric gravitational constant (TCB-compatible value; DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : false
        // Scalar: true
        static const double EARTH_GM = 3.986004576184e+14; // [m^3 s^-2]

        // Geocentric gravitational constant (TDB-compatible value; DE405 value). Do not use this parameter but use the TCB-compatible value from parameter :Nature:DE405:Earth_GM instead
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : false
        // Scalar: true
        static const double EARTH_GM_TDB = 3.986004514380e+14; // [m^3 s^-2 (TDB)]

        // Secular (long-term) variation of the dynamical form-factor J_2 of the Earth (also known as oblateness and as Stokes' second-degree zonal harmonic of the geopotential) due to the post-glacial rebound of the mantle (DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double EARTH_JSUB2DOT = 0.00e-01; // [cy^-1]

        // Love number k_20 of harmonic (2,0) of the Earth's harmonic potential expansion (rigid-Earth tide / slow zonal tides; DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double EARTH_LOVENUMBER_20 = 0.34;

        // Love number k_21 of harmonic (2,1) of the Earth's harmonic potential expansion (tidal deformation / diurnal tides; DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double EARTH_LOVENUMBER_21 = 0.30;

        // Love number k_22 of harmonic (2,1) of the Earth's harmonic potential expansion (rotational deformation / semi-diurnal tides; DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double EARTH_LOVENUMBER_22 = 0.30;

        // Harmonic potential coefficients of the Earth (DE405 values). The vector elements denote C_nm, with (n,m) = (2,0), (2,1), (2,2), (3,0), (3,1), (3,2), (3,3), (4,0), (4,1), (4,2), (4,3), and (4,4). A zonal harmonic J_n is a spherical harmonic of the form P_n(cos\theta), i.e., one which reduces to a Legendre polynomial of degree n. A tesseral harmonic C_nm/S_nm is a spherical harmonic of the form cos/sin(m\phi) P_n^m(cos\theta), where P_n^m is a Legendre function of degree n and order m. Special notations include -C_20 = J_2 = Stokes' second degree zonal harmonic (oblateness), -C_30 = J_3 = Stokes' third degree zonal harmonic, and -C_40 = J_4 = Stokes' fourth degree zonal harmonic
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double  *const EARTH_POTENTIALEXPANSION_C ()  { static double _v[12] = { -0.001082626,  0.,  0.,  0.000002533,  0.,  0.,  0.,  0.000001616,  0.,  0.,  0.,  0. }; return _v; }

        // Degree of harmonic expansion of the Earth's potential (DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const unsigned EARTH_POTENTIALEXPANSION_DEGREE = 4;

        // Harmonic potential coefficients of the Earth (DE405 values). The vector elements denote S_nm, with (n,m) = (2,0), (2,1), (2,2), (3,0), (3,1), (3,2), (3,3), (4,0), (4,1), (4,2), (4,3), and (4,4). A zonal harmonic J_n is a spherical harmonic of the form P_n(cos\theta), i.e., one which reduces to a Legendre polynomial of degree n. A tesseral harmonic C_nm/S_nm is a spherical harmonic of the form cos/sin(m\phi) P_n^m(cos\theta), where P_n^m is a Legendre function of degree n and order m
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double  *const EARTH_POTENTIALEXPANSION_S ()  { static double _v[12] = { 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0. }; return _v; }

        // Time delay \tau_20 used to compute tidal effects for harmonic (2,0) of the Earth's harmonic potential expansion (rigid-Earth tide / slow zonal tides; DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double EARTH_TIMEDELAY_20 = 0.; // [day]

        // Time delay \tau_21 used to compute tidal effects for harmonic (2,1) of the Earth's harmonic potential expansion (tidal deformation / diurnal tides; DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double EARTH_TIMEDELAY_21 = 0.01290895939; // [day]

        // Time delay \tau_22 used to compute tidal effects for harmonic (2,1) of the Earth's harmonic potential expansion (rotational deformation / semi-diurnal tides; DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double EARTH_TIMEDELAY_22 = 0.00694178558; // [day]

        // Ratio of Moon to Earth mass (DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : false
        // Scalar: true
        static const double MOONTOEARTH_MASSRATIO = 0.0123000383;

        // Equatorial radius of the Moon (DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double MOON_EQUATORIALRADIUS = 1738000.; // [m]

        // Love number k_2 of the Moon's harmonic potential expansion, assumed to be the same for all harmonic coefficients of order 2 (DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double MOON_LOVENUMBER_2 = 0.0299221167;

        // Harmonic potential coefficients of the Moon (DE405 values). The vector elements denote C_nm, with (n,m) = (2,0), (2,1), (2,2), (3,0), (3,1), (3,2), (3,3), (4,0), (4,1), (4,2), (4,3), and (4,4). A zonal harmonic J_n is a spherical harmonic of the form P_n(cos\theta), i.e., one which reduces to a Legendre polynomial of degree n. A tesseral harmonic C_nm/S_nm is a spherical harmonic of the form cos/sin(m\phi) P_n^m(cos\theta), where P_n^m is a Legendre function of degree n and order m. Special notations include -C_20 = J_2 = Stokes' second degree zonal harmonic (oblateness), -C_30 = J_3 = Stokes' third degree zonal harmonic, and -C_40 = J_4 = Stokes' fourth degree zonal harmonic. The second-degree lunar gravity field is time variable and the time-variable harmonic coefficients are computed in the DE405 ephemeris from the time-variable moment-of-inertia tensor. The numerical values of the C_20, C_21, and C_22 coefficients reported here (-99.99) are hence spurious
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double  *const MOON_POTENTIALEXPANSION_C ()  { static double _v[12] = { -99.99,  -99.99,  -99.99,  -0.000008785470,  0.000030803810,  0.000004879807,  0.000001770176,  1.45383E-7,  -0.000007177801,  -0.000001439518,  -8.5479E-8,  -1.54904E-7 }; return _v; }

        // Degree of harmonic expansion of the Moon's potential (DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const unsigned MOON_POTENTIALEXPANSION_DEGREE = 4;

        // Harmonic potential coefficients of the Moon (DE405 values). The vector elements denote S_nm, with (n,m) = (2,0), (2,1), (2,2), (3,0), (3,1), (3,2), (3,3), (4,0), (4,1), (4,2), (4,3), and (4,4). A zonal harmonic J_n is a spherical harmonic of the form P_n(cos\theta), i.e., one which reduces to a Legendre polynomial of degree n. A tesseral harmonic C_nm/S_nm is a spherical harmonic of the form cos/sin(m\phi) P_n^m(cos\theta), where P_n^m is a Legendre function of degree n and order m. The second-degree lunar gravity field is time variable and the time-variable harmonic coefficients are computed in the DE405 ephemeris from the time-variable moment-of-inertia tensor. The numerical values of the S_21 and S_22 coefficients reported here (-99.99) are hence spurious
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double  *const MOON_POTENTIALEXPANSION_S ()  { static double _v[12] = { 0.,  -99.99,  -99.99,  0.,  0.000004259329,  0.000001695516,  -2.70970E-7,  0.,  0.000002947434,  -0.000002884372,  -7.88967E-7,  5.6404E-8 }; return _v; }

        // Time delay \tau_2 used to compute tidal effects for the Moon's solid-body tide, assumed to be the same for all harmonic coefficients of order 2 (DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double MOON_TIMEDELAY_2 = 0.1667165558; // [day]

        // Ratio of Sun to Earth-system mass (DE405 value). The planetary mass includes the contribution of its satellite, the Moon
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double SUNTOEARTHSYSTEM_MASSRATIO = 328900.561400;

        // Ratio of Sun to Earth mass (DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : false
        // Scalar: true
        static const double SUNTOEARTH_MASSRATIO = 332946.050895;

        // Ratio of Sun to Jupiter-system mass (DE405 value). The planetary mass includes the contribution of its satellites
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048. Numerical value: J.K. Campbell, S.P. Synnott, 1985, 'Gravity field of the Jovian system from Pioneer and Voyager tracking data', AJ, 90, 364; this reference lists GM = 126712767 km^3 s^-2
        // Basic : true
        // Scalar: true
        static const double SUNTOJUPITERSYSTEM_MASSRATIO = 1047.3486;

        // Ratio of Sun to Mars-system mass (DE405 value). The planetary mass includes the contribution of its satellites
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048. Numerical value: G.W. Null, 1969, 'A solution for the mass and dynamical oblateness of Mars using Mariner-IV Doppler data', Bull. Am. Astron. Soc., 1, 356
        // Basic : true
        // Scalar: true
        static const double SUNTOMARSSYSTEM_MASSRATIO = 3098708.;

        // Ratio of Sun to Mercury(-system) mass (DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048. Numerical value: J.D. Anderson, et al., 1987, 'The mass, gravity field, and ephemeris of Mercury', Icarus, 71, 337; this reference lists GM = 22032.09 km^3 s^-2
        // Basic : true
        // Scalar: true
        static const double SUNTOMERCURYSYSTEM_MASSRATIO = 6023600.;

        // Ratio of Sun to Neptune-system mass (DE405 value). The planetary mass includes the contribution of its satellites
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048. Numerical value: R.A. Jacobson, et al., 1991, 'The orbits of Triton and Nereid from spacecraft and Earth-based observations', A&A, 247, 565; this reference lists GM = 6836535 km^3 s^-2
        // Basic : true
        // Scalar: true
        static const double SUNTONEPTUNESYSTEM_MASSRATIO = 19412.24;

        // Ratio of Sun to Pluto-system mass (DE405 value). The 'planetary' mass includes the contribution of its satellite, Charon
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048. Numerical value: D.J. Tholen, M.W. Buie, 1997, 'The Orbit of Charon', Icarus, 125, 245, although these authors list 1.3522E8 rather than 1.3521E8; this reference lists GM = 981.5 km^3 s^-2
        // Basic : true
        // Scalar: true
        static const double SUNTOPLUTOSYSTEM_MASSRATIO = 135200000.;

        // Ratio of Sun to Saturn-system mass (DE405 value). The planetary mass includes the contribution of its satellites
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048. Numerical value: J.K. Campbell, J.D. Anderson, 1989, 'Gravity field of the Saturnian system from Pioneer and Voyager tracking data', AJ, 97, 1485; this reference lists GM = 37940630 km^3 s^-2
        // Basic : true
        // Scalar: true
        static const double SUNTOSATURNSYSTEM_MASSRATIO = 3497.898;

        // Ratio of Sun to Uranus-system mass (DE405 value). The planetary mass includes the contribution of its satellites
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048. Numerical value: R.A. Jacobson, et al., 1992, AJ, 'The masses of Uranus and its major satellites from Voyager tracking data and Earth-based Uranian satellite data', 103, 2068; this reference lists GM = 5794548.6 km^3 s^-2
        // Basic : true
        // Scalar: true
        static const double SUNTOURANUSSYSTEM_MASSRATIO = 22902.98;

        // Ratio of Sun to Venus(-system) mass (DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048. Numerical value: W.L. Sjogren, et al., 1990, 'Venus - A total mass estimate', Geophysical Research Letters, 17, 1485; this reference lists GM = 324858.60 km^3 s^-2
        // Basic : true
        // Scalar: true
        static const double SUNTOVENUSSYSTEM_MASSRATIO = 408523.71;

        // Heliocentric gravitational constant (TCB-compatible value; DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : false
        // Scalar: true
        static const double SUN_GM = 1.327124482489e+20; // [m^3 s^-2]

        // Heliocentric gravitational constant (TDB-compatible value; DE405 value). Do not use this parameter but use the TCB-compatible value from parameter :Nature:DE405:Sun_GM instead
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : false
        // Scalar: true
        static const double SUN_GM_TDB = 1.327124461912e+20; // [m^3 s^-2 (TDB)]

        // Dynamical form-factor of the Sun (Stokes' second-degree zonal harmonic of the solar potential; DE405 value)
        // Source: E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048
        // Basic : true
        // Scalar: true
        static const double SUN_JSUB2 = 2.0e-07;

    };

    class DE410 {
    public:
        // Ratio of 1 Ceres to solar mass (DE410 value)
        // Source: E.M. Standish, 24 April 2003, 'JPL Planetary Ephemeris DE410', JPL IOM 312.N-03-009, Table III
        // Basic : true
        // Scalar: true
        static const double ASTEROID1CERESMASS_SOLARMASS = 4.690e-10;

        // Ratio of 2 Pallas to solar mass (DE410 value)
        // Source: E.M. Standish, 24 April 2003, 'JPL Planetary Ephemeris DE410', JPL IOM 312.N-03-009, Table III
        // Basic : true
        // Scalar: true
        static const double ASTEROID2PALLASMASS_SOLARMASS = 1.050e-10;

        // Ratio of 4 Vesta to solar mass (DE410 value)
        // Source: E.M. Standish, 24 April 2003, 'JPL Planetary Ephemeris DE410', JPL IOM 312.N-03-009, Table III
        // Basic : true
        // Scalar: true
        static const double ASTEROID4VESTAMASS_SOLARMASS = 1.360e-10;

        // Ratio of the Krasinsky asteroid ring to solar mass (originally expressed in terms of M_Ceres; DE410 value). Following G.A. Krasinsky, E.V. Pitjeva, M.V. Vasilyev, E.I. Yagudina, 1 February 2002, 'Hidden Mass in the Asteroid Belt', Icarus, 158, 98-105, the gravitational effect of all but the 300 heaviest asteroids can be modeled as an acceleration caused by a solid ring of this mass in the ecliptic plane (see also E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048, Section IV)
        // Source: E.M. Standish, 24 April 2003, 'JPL Planetary Ephemeris DE410', JPL IOM 312.N-03-009, Table III
        // Basic : false
        // Scalar: true
        static const double ASTEROIDRINGMASS_SOLARMASS = 1.032e-10;

        // Barycentric distance (orbital semi-major axis) of the Krasinsky asteroid ring (DE410 value). Following G.A. Krasinsky, E.V. Pitjeva, M.V. Vasilyev, E.I. Yagudina, 1 February 2002, 'Hidden Mass in the Asteroid Belt', Icarus, 158, 98-105, the gravitational effect of all but the 300 heaviest asteroids can be modeled as an acceleration caused by a solid ring with this barycentric distance/radius in the ecliptic plane (see also E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048, Section IV)
        // Source: E.M. Standish, 24 April 2003, 'JPL Planetary Ephemeris DE410', JPL IOM 312.N-03-009, Table III
        // Basic : true
        // Scalar: true
        static const double ASTEROIDRING_ORBITALSEMIMAJORAXIS = 2.8; // [au]

        // Mean mass density of C-class asteroids (DE410 value). In JPL's DE ephemerides, masses for the 300 most massive asteroids (except Ceres, Pallas, and Vesta) are derived using the relation GM = 6.27E-22 Radius^3 \rho, where Radius is the known asteroid radius in km and \rho is the mean mass density in g cm^-3 (GM in this relation is in 'solar-system units', i.e., Sun_GM = Gauss_Constant * Gauss_Constant [au^3/2 day^-1 M_Sun^-1/2]; note that the reference erroneously lists a prefactor 6.27E20; typo confirmed by E.M. Standish, priv. comm., 18 September 2003). The mean density \rho is assumed to be constant within each of the three taxonomic asteroid classes C, S, and M (see E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048, Section IV)
        // Source: E.M. Standish, 24 April 2003, 'JPL Planetary Ephemeris DE410', JPL IOM 312.N-03-009, Table III
        // Basic : true
        // Scalar: true
        static const double ASTEROID_MASSDENSITY_MEANCLASSC = 1.55; // [g cm^-3]

        // Mean mass density of M-class asteroids (DE410 value). In JPL's DE ephemerides, masses for the 300 most massive asteroids (except Ceres, Pallas, and Vesta) are derived using the relation GM = 6.27E-22 Radius^3 \rho, where Radius is the known asteroid radius in km and \rho is the mean mass density in g cm^-3 (GM in this relation is in 'solar-system units', i.e., Sun_GM = Gauss_Constant * Gauss_Constant [au^3/2 day^-1 M_Sun^-1/2]; note that the reference erroneously lists a prefactor 6.27E20; typo confirmed by E.M. Standish, priv. comm., 18 September 2003). The mean density \rho is assumed to be constant within each of the three taxonomic asteroid classes C, S, and M (see E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048, Section IV)
        // Source: E.M. Standish, 24 April 2003, 'JPL Planetary Ephemeris DE410', JPL IOM 312.N-03-009, Table III
        // Basic : true
        // Scalar: true
        static const double ASTEROID_MASSDENSITY_MEANCLASSM = 4.5; // [g cm^-3]

        // Mean mass density of S-class asteroids (DE410 value). In JPL's DE ephemerides, masses for the 300 most massive asteroids (except Ceres, Pallas, and Vesta) are derived using the relation GM = 6.27E-22 Radius^3 \rho, where Radius is the known asteroid radius in km and \rho is the mean mass density in g cm^-3 (GM in this relation is in 'solar-system units', i.e., Sun_GM = Gauss_Constant * Gauss_Constant [au^3/2 day^-1 M_Sun^-1/2]; note that the reference erroneously lists a prefactor 6.27E20; typo confirmed by E.M. Standish, priv. comm., 18 September 2003). The mean density \rho is assumed to be constant within each of the three taxonomic asteroid classes C, S, and M (see E.M. Standish, 26 August 1998, 'JPL Planetary and Lunar Ephemerides, DE405/LE405', JPL IOM 312.F-98-048, Section IV)
        // Source: E.M. Standish, 24 April 2003, 'JPL Planetary Ephemeris DE410', JPL IOM 312.N-03-009, Table III
        // Basic : true
        // Scalar: true
        static const double ASTEROID_MASSDENSITY_MEANCLASSS = 2.13; // [g cm^-3]

        // Dynamical form-factor of the Sun (Stokes' second-degree zonal harmonic of the solar potential; DE410 value)
        // Source: E.M. Standish, 24 April 2003, 'JPL Planetary Ephemeris DE410', JPL IOM 312.N-03-009, Table III
        // Basic : true
        // Scalar: true
        static const double SUN_JSUB2 = 2.90e-07;

    };

    // Number of seconds per day
    // Source: IAU definition
    // Basic : false
    // Scalar: true
    static const double DAY_SECOND = 86400.; // [s]

    // One degree in units of radians
    // Basic : false
    // Scalar: true
    static const double DEGREE_RADIAN = 1.745329251994330e-02; // [rad]

    // Mean orbital eccentricity of the Earth-Moon barycentre (EMBC) orbit, at the standard epoch J2000.0. The mean orbital eccentricity is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. In this fit, each orbital element is allowed to vary linearly with time (the resulting evolution of the orbital eccentricity is -0.00004392 radians per century). The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0. This solution fits the DE405 orbit of the Earth-Moon barycentre to about 22 arcsec. DE405 is based upon the International Celestial Reference Frame (ICRF)
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double EMBC_ORBITALECCENTRICITY_J2000 = 0.01671123;

    // Mean orbital inclination of the Earth-Moon barycentre (EMBC) orbit, at the standard epoch J2000.0. The mean orbital inclination is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. In this fit, each orbital element is allowed to vary linearly with time (the resulting evolution of the orbital inclination is -0.01294668 degrees per century). The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0. This solution fits the DE405 orbit of the Earth-Moon barycentre to about 22 arcsec. DE405 is based upon the International Celestial Reference Frame (ICRF)
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double EMBC_ORBITALINCLINATION_J2000 = -0.00001531; // [deg]

    // Mean orbital semi-major axis of the Earth-Moon barycentre (EMBC) orbit, at the standard epoch J2000.0. The mean orbital semi-major axis is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. In this fit, each orbital element is allowed to vary linearly with time (the resulting evolution of the orbital semi-major axis is 0.00000562 au per century). The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0. This solution fits the DE405 orbit of the Earth-Moon barycentre to about 22 arcsec. DE405 is based upon the International Celestial Reference Frame (ICRF)
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double EMBC_ORBITALSEMIMAJORAXIS_J2000 = 1.00000261; // [au]

    // Geocentric gravitational constant (TCB-compatible value), including the Earth's atmosphere but excluding the mass of the Moon
    // Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double EARTHELLIPSOID_GM = 3.9860044180e+14; // [m^3 s^-2]

    // Inverse of the geometrical flattening factor f of the Earth (f = (a-b)/a; zero-frequency-tide value)
    // Source: Numerical value is zero-frequency-tide value from E. Groten, 2000, 'Report of Special Commission 3 of IAG', in Proceedings of IAU Colloquium 180, 'Towards models and constants for sub-microarcsecond astrometry', eds K.J. Johnston, D.D. McCarthy, B.J. Luzum, G.H. Kaplan, page 337
    // Basic : true
    // Scalar: true
    static const double EARTHELLIPSOID_INVERSEFLATTENING_ZEROFREQUENCYTIDE = 298.25642;

    // Dynamical form-factor of the Earth, i.e., second-degree zonal harmonic of the geopotential including the indirect tidal distortion on J_2, i.e., in the zero-frequency-tide system JGM-3. The (long-term) rate of change of this parameter equals -3.001E-9 cy^-1
    // Source: Numerical value is zero-frequency-tide value from E. Groten, 2000, 'Report of Special Commission 3 of IAG', in Proceedings of IAU Colloquium 180, 'Towards models and constants for sub-microarcsecond astrometry', eds K.J. Johnston, D.D. McCarthy, B.J. Luzum, G.H. Kaplan, page 337. See also the IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double EARTHELLIPSOID_JSUB2_ZEROFREQUENCYTIDE = 1.08263590e-03;

    // Geopotential scale factor
    // Basic : false
    // Scalar: true
    static const double EARTHELLIPSOID_RSUB0 = 6.36367256e+06; // [m]

    // Semi-major axis of the Earth reference ellipsoid (zero-frequency-tide value)
    // Source: Numerical value is zero-frequency-tide value from E. Groten, 2000, 'Report of Special Commission 3 of IAG', in Proceedings of IAU Colloquium 180, 'Towards models and constants for sub-microarcsecond astrometry', eds K.J. Johnston, D.D. McCarthy, B.J. Luzum, G.H. Kaplan, page 337. See also the IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double EARTHELLIPSOID_SEMIMAJORAXIS_ZEROFREQUENCYTIDE = 6378136.6; // [m]

    // Nominal mean angular velocity of the Earth
    // Source: Numerical value is from E. Groten, 2000, 'Report of Special Commission 3 of IAG', in Proceedings of IAU Colloquium 180, 'Towards models and constants for sub-microarcsecond astrometry', eds K.J. Johnston, D.D. McCarthy, B.J. Luzum, G.H. Kaplan, page 337. See also the IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double EARTHELLIPSOID_SPINRATE_NOMINAL = 7.2921150e-05; // [rad s^-1]

    // Potential of the geoid
    // Source: Numerical value is from E. Groten, 2000, 'Report of Special Commission 3 of IAG', in Proceedings of IAU Colloquium 180, 'Towards models and constants for sub-microarcsecond astrometry', eds K.J. Johnston, D.D. McCarthy, B.J. Luzum, G.H. Kaplan, page 337. See also the IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double EARTHELLIPSOID_WSUB0 = 6.263685600e+07; // [m^2 s^-2]

    // Astrometric signature of the Sun induced by the Earth system for an observer located at a distance of 10 pc from the Sun
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.7, Equation 1.22, page 11
    // Basic : false
    // Scalar: true
    static const double EARTHSYSTEM_ASTROMETRICSIGNATURE_10PARSEC = 0.304; // [10^-6 arcsec]

    // Sidereal orbital period
    // Source: Values derived from the mean longitude rates in Table 5.8.1 in P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, page 316; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double EARTHSYSTEM_ORBITALPERIOD = 1.0000174; // [yr]

    // Radial-velocity amplitude of the Sun induced by the Earth system for 'an edge-on observer' (i.e., an observer in the orbital plane of the Earth system)
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.4, Equation 1.18, page 9
    // Basic : false
    // Scalar: true
    static const double EARTHSYSTEM_RADIALVELOCITYSIGNATURE = 0.091; // [m s^-1]

    // Ratio of Earth to Moon mass
    // Basic : false
    // Scalar: true
    static const double EARTHTOMOON_MASSRATIO = 81.30057;

    // Eccentricity e of the Earth (its shape, not its orbit)
    // Basic : false
    // Scalar: true
    static const double EARTH_ECCENTRICITY = 8.181930088e-02;

    // Radius of the smallest hypothetical sphere around the Earth which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double EARTH_ENCOMPASSINGSPHERERADIUS = 6378137.; // [m]

    // Equatorial radius of the Earth
    // Basic : false
    // Scalar: true
    static const double EARTH_EQUATORIALRADIUS = 6.37813660e+06; // [m]

    // Nominal equatorial radius of the Earth (zero-frequency-tide value), in SI units. This nominal value shall be understood as conversion factor only
    // Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)
    // Basic : true
    // Scalar: true
    static const double EARTH_EQUATORIALRADIUS_NOMINAL = 6.37810e+06; // [m]

    // Geometrical flattening factor f of the Earth (f = (a-b)/a); this quantity is also refered to as ellipticity, but is not identical to eccentricity
    // Basic : false
    // Scalar: true
    static const double EARTH_FLATTENING = 3.352819698e-03;

    // Maximum reduction of the solar flux for an observer external to the solar system during a transit of Earth
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.2, Equation 2.4, page 14
    // Basic : false
    // Scalar: true
    static const double EARTH_FLUXREDUCTION_MAXIMUM = 0.008; // [%]

    // Geocentric gravitational constant (TCB-compatible value), including the Earth's atmosphere but excluding the mass of the Moon
    // Basic : false
    // Scalar: true
    static const double EARTH_GM = 3.9860044180e+14; // [m^3 s^-2]

    // Nominal GM of the Earth, in SI units. This nominal value shall be understood as conversion factor only
    // Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)
    // Basic : true
    // Scalar: true
    static const double EARTH_GM_NOMINAL = 3.9860040e+14; // [m^3 s^-2]

    // Geometric albedo of the Earth (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double EARTH_GEOMETRICALBEDO = 0.367;

    // Earth mass, including its atmosphere but excluding the mass of the Moon
    // Basic : false
    // Scalar: true
    static const double EARTH_MASS = 5.97237e+24; // [kg]

    // Mean mass density of the Earth
    // Basic : false
    // Scalar: true
    static const double EARTH_MASSDENSITY_MEAN = 5.514; // [g cm^-3]

    // IAU-recommended value for the declination \delta_0 of the north pole of rotation of Earth. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value shall be used for comparative purposes only
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double EARTH_NORTHROTATIONALPOLE_DECLINATION = 90.00; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the declination \delta_0 of the north pole of rotation of Earth. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value shall be used for comparative purposes only
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : false
    // Scalar: true
    static const double EARTH_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE = -0.00001525; // [deg day^-1]

    // IAU-recommended value for the right ascension \alpha_0 of the north pole of rotation of Earth. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value shall be used for comparative purposes only
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double EARTH_NORTHROTATIONALPOLE_RIGHTASCENSION = 0.00; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the right ascension \alpha_0 of the north pole of rotation of Earth. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value shall be used for comparative purposes only
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : false
    // Scalar: true
    static const double EARTH_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE = -0.00001755; // [deg day^-1]

    // Mean orbital speed of the Earth (mean velocity over an unperturbed elliptic orbit). The equation is accurate to 4-th order in EMBC_OrbitalEccentricity_J2000
    // Basic : false
    // Scalar: true
    static const double EARTH_ORBITALSPEED_MEAN = 2.97827e+04; // [m s^-1]

    // Mean polar radius of the Earth
    // Basic : false
    // Scalar: true
    static const double EARTH_POLARRADIUS = 6.35675186e+06; // [m]

    // Nominal polar radius of the Earth (zero-frequency-tide value), in SI units. This nominal value shall be understood as conversion factor only
    // Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)
    // Basic : true
    // Scalar: true
    static const double EARTH_POLARRADIUS_NOMINAL = 6.35680e+06; // [m]

    // IAU-recommended value for the ephemeris position of the prime meridian of Earth. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde. The numerical value shall be used for comparative purposes only
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double EARTH_PRIMEMERIDIAN_EPHEMERISPOSITION = 190.147; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch of 1.5 January 2000 = JD2451545.0 TDB) of the ephemeris position of the prime meridian of Earth. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde. The numerical value shall be used for comparative purposes only
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double EARTH_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE = 360.9856235; // [deg day^-1]

    // Mean radius of the Earth
    // Basic : false
    // Scalar: true
    static const double EARTH_RADIUS_MEAN = 6.371008e+06; // [m]

    // Surface gravity of the Earth
    // Basic : false
    // Scalar: true
    static const double EARTH_SURFACEGRAVITY = 9.798; // [m s^-2]

    // Mean surface gravity of the Earth. The value for the International Standard Atmopshere is 9.80665 m s^-2
    // Source: F. Budnik (ESA), 8 March 2013, 'Gaia FDS-SOC Orbit ICD', GAIA-ESC-ICD-0012, issue 2, revision 0, Annex A. Reference documents: F. Kleijer (Netherlands Geodetic Commission, Delft), 1 April 2004, 'Troposphere Modeling and Filtering for Precise GPS Leveling', Publications on Geodesy 56, ISBN 90 6132 284 7 (http://www.ncg.knaw.nl/Publicaties/Geodesy/pdf/56Kleijer.pdf and http://repository.tudelft.nl/view/ir/uuid%3Aea1f0cf0-4e48-421b-b7ae-4ae3e36d1880/), J. Saastamoinen, 1 January 1972, 'Atmospheric correction for the troposphere and stratosphere in radio ranging of satellites' in 'The use of artificial satellites for geodesy', editors S.W. Henrikson et al., Geophysical Monograph Series, 15, 247-251, and B.R. Bean and E.J. Dutton, 1 March 1966, 'Radio Meteorology', National Bureau of Standards Monograph, 92, 1-44
    // Basic : true
    // Scalar: true
    static const double EARTH_SURFACEGRAVITY_MEAN = 9.784; // [m s^-2]

    // Geometric transit probability (Earth transiting the Sun) for an observer external to the solar system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.1, Equation 2.2, page 14
    // Basic : false
    // Scalar: true
    static const double EARTH_TRANSITPROBABILITY = 0.469; // [%]

    // Maximum transit time of Earth (transiting the Sun) for an observer external to the solar system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.3, Equation 2.5, page 15
    // Basic : false
    // Scalar: true
    static const double EARTH_TRANSITTIME_MAXIMUM = 0.55; // [day]

    // V(1,0) magnitude of the Earth (i.e., the visual magnitude of the planet reduced to a distance of 1 au from both the Sun and Earth and phase angle zero). This parameter is also refered to as absolute magnitude in planetary sciences
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double EARTH_VONEZEROMAGNITUDE = -3.86; // [mag]

    // Mean volumetric radius of the Earth
    // Basic : false
    // Scalar: true
    static const double EARTH_VOLUMETRICRADIUS = 6.371000e+06; // [m]

    // Electric constant (defining constant)
    // Basic : false
    // Scalar: true
    static const double ELECTRIC_CONSTANT = 8.854187817620390e-12; // [F m^-1]

    // Classical electron radius. Note: best-measured value equals 2.8179403227E-15 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))
    // Basic : false
    // Scalar: true
    static const double ELECTRON_CLASSICALRADIUS = 2.81794032201e-15; // [m]

    // Thomson free-electron-scattering absorption coefficient (cross section per electron). Note: best-measured value equals 0.66524587158E-28 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))
    // Source: E.g., R. Kippenhahn, A. Weigert, 1991, 'Stellar structure and evolution' (corrected 2-nd printing), Springer Verlag, Berlin, Section 17, Equation 17.1, page 137
    // Basic : false
    // Scalar: true
    static const double ELECTRON_CROSSSECTION_THOMSONSCATTERING = 6.65245871237e-29; // [m^2]

    // Electron mass
    // Source: P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0)
    // Basic : true
    // Scalar: true
    static const double ELECTRON_MASS = 9.109383560e-31; // [kg]

    // Elementary charge
    // Source: P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0)
    // Basic : true
    // Scalar: true
    static const double ELEMENTARYCHARGE_CONSTANT = 1.60217662080e-19; // [C]

    // One erg expressed in units of J. Note that 'erg' is a non-SI unit which should not be used
    // Basic : true
    // Scalar: true
    static const double ERG_JOULE = 1.0e-07; // [J]

    // Interstellar extinction curve (see also J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1). First column: wavelength \lambda (in nm; \lambda = 200.0 to 1100.0). Second column: normalised interstellar extinction A(\lambda) / A(\lambda_{ref}) with \lambda_{ref} = 1000 / 1.82 = 549.45 nm. Note that it has become customary within the Gaia community to denote A(\lambda_{ref}) for \lambda_{ref} = 1000 / 1.82 = 549.45 nm as A(550 nm)
    // Source: J.A. Cardelli, G.C. Clayton, J.S. Mathis, 1989, 'The relationship between infrared, optical, and ultraviolet extinction', Astrophysical Journal (ApJ), 345, 245
    // Basic : true
    // Scalar: false
    static const char  *const EXTINCTION_CURVE() { return "Nature/Extinction_Curve_002.fits"; }

    // Typical extinction in the Johnson V band (A_V) per kpc in the Galactic plane; values ranging from 0.5 to 1.5 mag kpc^-1 are considered 'normal' (e.g., H. Jonch-Sorensen, 1994, 'CCD uvby-beta photometry of faint stars. 2: Reddening in six fields in the Galaxy', A&A, 292, 92)
    // Source: Typical value ('common lore')
    // Basic : true
    // Scalar: true
    static const double EXTINCTION_GALACTICPLANE_TYPICAL = 1.0; // [mag kpc^-1]

    // Ratio of total to selective absorption (typical value). One has: R_V = A_V / E(B-V), where A_V is the total extinction in the Johnson V band and E(B-V) is the colour excess (see also J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)
    // Source: J.A. Cardelli, G.C. Clayton, J.S. Mathis, 1989, 'The relationship between infrared, optical, and ultraviolet extinction', Astrophysical Journal (ApJ), 345, 245
    // Basic : true
    // Scalar: true
    static const double EXTINCTION_TOTALTOSELECTIVERATIO = 3.10;

    // Photon flux density N_{\lambda}(\lambda) of an unreddened F2V star (Pickles' star number 015) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)
    // Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863
    // Basic : true
    // Scalar: false
    static const char  *const F2VSTAR_SPECTRUM_NUMBEROFPHOTONS() { return "Nature/F2VStar_Spectrum_NumberOfPhotons_001.fits"; }

    // Photon flux density N_{\lambda}(\lambda) of an unreddened F6V star (Pickles' star number 018) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)
    // Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863
    // Basic : true
    // Scalar: false
    static const char  *const F6VSTAR_SPECTRUM_NUMBEROFPHOTONS() { return "Nature/F6VStar_Spectrum_NumberOfPhotons_001.fits"; }

    // Photon flux density N_{\lambda}(\lambda) of an unreddened F8V star (Pickles' star number 020) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)
    // Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863
    // Basic : true
    // Scalar: false
    static const char  *const F8VSTAR_SPECTRUM_NUMBEROFPHOTONS() { return "Nature/F8VStar_Spectrum_NumberOfPhotons_001.fits"; }

    // Johnson B band filter profile (normalised standard passband in the Johnson-Cousins UBVRI photometric system). First column: wavelength \lambda (in nm; from 360.0 to 560.0). Second column: normalised transmittance
    // Source: M.S. Bessell, S. Murphy, 2012, 'Spectrophotometric Libraries, Revised Photonic Passbands, and Zero Points for UBVRI, Hipparcos, and Tycho Photometry', PASP, 124, 140
    // Basic : true
    // Scalar: false
    static const char  *const FILTERTRANSMISSIONCURVE_JOHNSONCOUSINSB() { return "Nature/FilterTransmissionCurve_JohnsonCousinsB_002.fits"; }

    // Cousins I band filter profile (normalised standard passband in the Johnson-Cousins UBVRI photometric system). First column: wavelength \lambda (in nm; from 700.0 to 920.0). Second column: normalised transmittance
    // Source: M.S. Bessell, S. Murphy, 2012, 'Spectrophotometric Libraries, Revised Photonic Passbands, and Zero Points for UBVRI, Hipparcos, and Tycho Photometry', PASP, 124, 140
    // Basic : true
    // Scalar: false
    static const char  *const FILTERTRANSMISSIONCURVE_JOHNSONCOUSINSI() { return "Nature/FilterTransmissionCurve_JohnsonCousinsI_002.fits"; }

    // Cousins R band filter profile (normalised standard passband in the Johnson-Cousins UBVRI photometric system). First column: wavelength \lambda (in nm; from 550.0 to 910.0). Second column: normalised transmittance
    // Source: M.S. Bessell, S. Murphy, 2012, 'Spectrophotometric Libraries, Revised Photonic Passbands, and Zero Points for UBVRI, Hipparcos, and Tycho Photometry', PASP, 124, 140
    // Basic : true
    // Scalar: false
    static const char  *const FILTERTRANSMISSIONCURVE_JOHNSONCOUSINSR() { return "Nature/FilterTransmissionCurve_JohnsonCousinsR_002.fits"; }

    // Johnson V band filter profile (normalised standard passband in the Johnson-Cousins UBVRI photometric system). First column: wavelength \lambda (in nm; from 470.0 to 740.0). Second column: normalised transmittance
    // Source: M.S. Bessell, S. Murphy, 2012, 'Spectrophotometric Libraries, Revised Photonic Passbands, and Zero Points for UBVRI, Hipparcos, and Tycho Photometry', PASP, 124, 140
    // Basic : true
    // Scalar: false
    static const char  *const FILTERTRANSMISSIONCURVE_JOHNSONCOUSINSV() { return "Nature/FilterTransmissionCurve_JohnsonCousinsV_002.fits"; }

    // Fine structure constant. Note: best-measured value equals 7.2973525664E-3 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))
    // Basic : false
    // Scalar: true
    static const double FINESTRUCTURE_CONSTANT = 7.29735256621e-03;

    // Foreshortening constant A_z^-1 (see ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, page 25, Table 1.2.2)
    // Basic : false
    // Scalar: true
    static const double FORESHORTENING_CONSTANT = 1.0227121650e-09; // [mas^-1 km^-1 yr^-1 s]

    // Surface area of unit sphere (4 Pi steradians) in units of square degrees
    // Basic : false
    // Scalar: true
    static const double FOURPISTERADIAN_DEGREESQUARE = 41252.9612494193; // [deg^2]

    // Linear thermal expansion coefficient of synthetic fused Silica at 293 K. The quoted value is specified to be valid over the temperature range T = 293 - 373 K; a value for T = 170 K is not available
    // Source: Schott Lithotec AG, 5 August 2010, 'Optical glass data sheets (Lithosil-Q)', http://www.schott.com/advanced_optics/english/download/- schott_optical_glass_august_2010_en.pdf; see also http://www.schott.com/advanced_optics/english/download/- schott_fused_silica_jan_2010_en_brochure.pdf
    // Basic : true
    // Scalar: true
    static const double FUSEDSILICA_LINEARTHERMALCOEFFICIENTOFEXPANSION_293K = 0.5; // [ppm K^-1]

    // Density of synthetic fused Silica
    // Source: Schott Lithotec AG, 5 August 2010, 'Optical glass data sheets (Lithosil-Q)', http://www.schott.com/advanced_optics/english/download/- schott_optical_glass_august_2010_en.pdf; see also http://www.schott.com/advanced_optics/english/download/- schott_fused_silica_jan_2010_en_brochure.pdf
    // Basic : true
    // Scalar: true
    static const double FUSEDSILICA_MASSDENSITY = 2.2; // [g cm^-3]

    // Sellmeier coefficient B_1 of synthetic fused Silica, which is dimensionless, at T = 120 K. The Sellmeier equation is an empirical relation between the refractive index n and wavelength \lambda for transparent media in the form: n^2(\lambda) = 1 + \frac{B_1 \lambda^2}{\lambda^2 - C_1} + \frac{B_2 \lambda^2}{\lambda^2 - C_2} + \frac{B_3 \lambda^2}{\lambda^2 - C_3}, where B_1, B_2, B_3 and C_1, C_2, C_3 are experimentally-determined Sellmeier coefficients. These coefficients are defined for \lambda measured in 10^-6 m. The wavelength \lambda is the vacuum wavelength and not that in the material itself, which is \lambda / n(\lambda)
    // Source: EADS-Astrium, 4 July 2008, 'Code V Models Description', GAIA.ASF.TCN.PLM.00287, issue 1, revision 0. For temperature-dependent Sellmeier coefficients, see G. Ghosh, M. Endo, T. Iwasaki, 1 August 1994, 'Temperature-dependent Sellmeier coefficients and chromatic dispersions for some optical fiber glasses', Journal of Lightwave Technology, Volume 12, Number 8, pages 1338-1342 (1994JLwT...12.1338G) and http://www.schott.com/advanced_optics/english/download/- schott_tie-29_refractive_index_v3_jan_2007_en.pdf
    // Basic : true
    // Scalar: true
    static const double FUSEDSILICA_SELLMEIERCOEFFICIENT_B1 = 1.10053898145;

    // Sellmeier coefficient B_2 of synthetic fused Silica, which is dimensionless, at T = 120 K. The Sellmeier equation is an empirical relation between the refractive index n and wavelength \lambda for transparent media in the form: n^2(\lambda) = 1 + \frac{B_1 \lambda^2}{\lambda^2 - C_1} + \frac{B_2 \lambda^2}{\lambda^2 - C_2} + \frac{B_3 \lambda^2}{\lambda^2 - C_3}, where B_1, B_2, B_3 and C_1, C_2, C_3 are experimentally-determined Sellmeier coefficients. These coefficients are defined for \lambda measured in 10^-6 m. The wavelength \lambda is the vacuum wavelength and not that in the material itself, which is \lambda / n(\lambda)
    // Source: EADS-Astrium, 4 July 2008, 'Code V Models Description', GAIA.ASF.TCN.PLM.00287, issue 1, revision 0. For temperature-dependent Sellmeier coefficients, see G. Ghosh, M. Endo, T. Iwasaki, 1 August 1994, 'Temperature-dependent Sellmeier coefficients and chromatic dispersions for some optical fiber glasses', Journal of Lightwave Technology, Volume 12, Number 8, pages 1338-1342 (1994JLwT...12.1338G) and http://www.schott.com/advanced_optics/english/download/- schott_tie-29_refractive_index_v3_jan_2007_en.pdf
    // Basic : true
    // Scalar: true
    static const double FUSEDSILICA_SELLMEIERCOEFFICIENT_B2 = 0.00144043087;

    // Sellmeier coefficient B_3 of synthetic fused Silica, which is dimensionless, at T = 120 K. The Sellmeier equation is an empirical relation between the refractive index n and wavelength \lambda for transparent media in the form: n^2(\lambda) = 1 + \frac{B_1 \lambda^2}{\lambda^2 - C_1} + \frac{B_2 \lambda^2}{\lambda^2 - C_2} + \frac{B_3 \lambda^2}{\lambda^2 - C_3}, where B_1, B_2, B_3 and C_1, C_2, C_3 are experimentally-determined Sellmeier coefficients. These coefficients are defined for \lambda measured in 10^-6 m. The wavelength \lambda is the vacuum wavelength and not that in the material itself, which is \lambda / n(\lambda)
    // Source: EADS-Astrium, 4 July 2008, 'Code V Models Description', GAIA.ASF.TCN.PLM.00287, issue 1, revision 0. For temperature-dependent Sellmeier coefficients, see G. Ghosh, M. Endo, T. Iwasaki, 1 August 1994, 'Temperature-dependent Sellmeier coefficients and chromatic dispersions for some optical fiber glasses', Journal of Lightwave Technology, Volume 12, Number 8, pages 1338-1342 (1994JLwT...12.1338G) and http://www.schott.com/advanced_optics/english/download/- schott_tie-29_refractive_index_v3_jan_2007_en.pdf
    // Basic : true
    // Scalar: true
    static const double FUSEDSILICA_SELLMEIERCOEFFICIENT_B3 = 0.77782846144;

    // Sellmeier coefficient C_1 of synthetic fused Silica, in units of (10^-6 m)^2, at T = 120 K. It is emphasised that this C-value is already squared, thus complying with the denominator units of the Sellmeier equation, in (10^-6 m)^2. The Sellmeier equation is an empirical relation between the refractive index n and wavelength \lambda for transparent media in the form: n^2(\lambda) = 1 + \frac{B_1 \lambda^2}{\lambda^2 - C_1} + \frac{B_2 \lambda^2}{\lambda^2 - C_2} + \frac{B_3 \lambda^2}{\lambda^2 - C_3}, where B_1, B_2, B_3 and C_1, C_2, C_3 are experimentally-determined Sellmeier coefficients. These coefficients are defined for \lambda measured in 10^-6 m. The wavelength \lambda is the vacuum wavelength and not that in the material itself, which is \lambda / n(\lambda)
    // Source: EADS-Astrium, 4 July 2008, 'Code V Models Description', GAIA.ASF.TCN.PLM.00287, issue 1, revision 0. For temperature-dependent Sellmeier coefficients, see G. Ghosh, M. Endo, T. Iwasaki, 1 August 1994, 'Temperature-dependent Sellmeier coefficients and chromatic dispersions for some optical fiber glasses', Journal of Lightwave Technology, Volume 12, Number 8, pages 1338-1342 (1994JLwT...12.1338G) and http://www.schott.com/advanced_optics/english/download/- schott_tie-29_refractive_index_v3_jan_2007_en.pdf
    // Basic : false
    // Scalar: true
    static const double FUSEDSILICA_SELLMEIERCOEFFICIENT_C1 = 0.00787874390; // [(10^-6 m)^2]

    // Sellmeier coefficient C_2 of synthetic fused Silica, in units of (10^-6 m)^2, at T = 120 K. It is emphasised that this C-value is already squared, thus complying with the denominator units of the Sellmeier equation, in (10^-6 m)^2. The Sellmeier equation is an empirical relation between the refractive index n and wavelength \lambda for transparent media in the form: n^2(\lambda) = 1 + \frac{B_1 \lambda^2}{\lambda^2 - C_1} + \frac{B_2 \lambda^2}{\lambda^2 - C_2} + \frac{B_3 \lambda^2}{\lambda^2 - C_3}, where B_1, B_2, B_3 and C_1, C_2, C_3 are experimentally-determined Sellmeier coefficients. These coefficients are defined for \lambda measured in 10^-6 m. The wavelength \lambda is the vacuum wavelength and not that in the material itself, which is \lambda / n(\lambda)
    // Source: EADS-Astrium, 4 July 2008, 'Code V Models Description', GAIA.ASF.TCN.PLM.00287, issue 1, revision 0. For temperature-dependent Sellmeier coefficients, see G. Ghosh, M. Endo, T. Iwasaki, 1 August 1994, 'Temperature-dependent Sellmeier coefficients and chromatic dispersions for some optical fiber glasses', Journal of Lightwave Technology, Volume 12, Number 8, pages 1338-1342 (1994JLwT...12.1338G) and http://www.schott.com/advanced_optics/english/download/- schott_tie-29_refractive_index_v3_jan_2007_en.pdf
    // Basic : false
    // Scalar: true
    static const double FUSEDSILICA_SELLMEIERCOEFFICIENT_C2 = 0.07320427965; // [(10^-6 m)^2]

    // Sellmeier coefficient C_3 of synthetic fused Silica, in units of (10^-6 m)^2, at T = 120 K. It is emphasised that this C-value is already squared, thus complying with the denominator units of the Sellmeier equation, in (10^-6 m)^2. The Sellmeier equation is an empirical relation between the refractive index n and wavelength \lambda for transparent media in the form: n^2(\lambda) = 1 + \frac{B_1 \lambda^2}{\lambda^2 - C_1} + \frac{B_2 \lambda^2}{\lambda^2 - C_2} + \frac{B_3 \lambda^2}{\lambda^2 - C_3}, where B_1, B_2, B_3 and C_1, C_2, C_3 are experimentally-determined Sellmeier coefficients. These coefficients are defined for \lambda measured in 10^-6 m. The wavelength \lambda is the vacuum wavelength and not that in the material itself, which is \lambda / n(\lambda)
    // Source: EADS-Astrium, 4 July 2008, 'Code V Models Description', GAIA.ASF.TCN.PLM.00287, issue 1, revision 0. For temperature-dependent Sellmeier coefficients, see G. Ghosh, M. Endo, T. Iwasaki, 1 August 1994, 'Temperature-dependent Sellmeier coefficients and chromatic dispersions for some optical fiber glasses', Journal of Lightwave Technology, Volume 12, Number 8, pages 1338-1342 (1994JLwT...12.1338G) and http://www.schott.com/advanced_optics/english/download/- schott_tie-29_refractive_index_v3_jan_2007_en.pdf
    // Basic : false
    // Scalar: true
    static const double FUSEDSILICA_SELLMEIERCOEFFICIENT_C3 = 85.64043680329; // [(10^-6 m)^2]

    // Typical transmission of synthetic fused Silica, including Fresnel reflection losses for an uncoated surface, for a 10-mm path length. First column: wavelength \lambda (in nm; from 200.0 to 1250.0). Second column: typical transmission. Explanatory, supplementary information: this parameter provides the bulk transmission of 10 mm of synthetic fused Silica, including Fresnel reflection losses. Fresnel diffraction losses, however, are only applicable in the absence of an anti-relfection coating. Since all Gaia prisms (BP/RP and RVS) do have anti-relfection coatings, the bulk transmission (i) first has to be corrected for (the absence of) Fresnel diffraction losses, and (ii) subsequently has to be scaled for the proper path length in the fused Silica. Ad (i): the Fresnel reflectivity per surface equals R = (n-1)^2 (n+1)^-2, where n is the index of refraction of fused Silica (which is a function of wavelength and temperature). The Fresnel loss per surface is hence 1 - R. The index of refraction n can be calculated from the Sellmeier coefficients (see parameters :Nature:FusedSilica_SellmeierCoefficient_*). The corrected bulk transmission of 10 mm of fused Silica exceeds 99.9% above 250 nm. Ad (ii): let us denote the transmission curve of 10 mm synthetic fused Silica, corrected for (the absence of) Fresnel diffraction losses, by \eta_10. The transmission curve for d mm path length (see parameters :Satellite:BP:Prism_Thickness, :Satellite:RP:Prism_Thickness, and :Satellite:RVS:Prism_Thickness) can then be calculated - following the Bouguer-Lambert law - as \eta_d = \eta_10^(d/10). As secondary effect, prism wedge angles could be included, effectively increasing the path lengths d
    // Source: Schott Lithotec AG, 5 August 2010, 'Optical glass data sheets (Lithosil-Q)', http://www.schott.com/advanced_optics/english/download/- schott_optical_glass_august_2010_en.pdf; see also http://www.schott.com/advanced_optics/english/download/- schott_fused_silica_jan_2010_en_brochure.pdf
    // Basic : true
    // Scalar: false
    static const char  *const FUSEDSILICA_TRANSMISSIVITY_10MM() { return "Nature/FusedSilica_Transmissivity_10mm_001.fits"; }

    // Photon flux density N_{\lambda}(\lambda) of an unreddened G2V star (Pickles' star number 026) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)
    // Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863
    // Basic : true
    // Scalar: false
    static const char  *const G2VSTAR_SPECTRUM_NUMBEROFPHOTONS() { return "Nature/G2VStar_Spectrum_NumberOfPhotons_001.fits"; }

    // High-resolution photon-flux density N_{\lambda}(\lambda) of an unreddened G2V star at V = 15 mag. The data refer to a high-resolution Kurucz-model spectrum with the following properties: effective temperature T_eff = 5800 K, logarithm of surface gravity log g = 4.5, metallicity [Fe/H] = 0.0, alpha-elements [\alpha/Fe] = 0.0, rotational velocity v sini = 5 km s^-1, micro-turbulence velocity = 2.0 km s^-1, length of convective bubble divided by pressure scale height = 0.50, no convective overshooting, macro-turbulence velocity = 2.0 km s^-1, and resolving power R = \lambda / \delta \lambda = 250,000. First column: wavelength \lambda (in nm; from 830.1673264 to 889.8217922). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1). The 34698 lines have an average wavelength step of 0.00172 nm; the spectrum extent is thus 59.7 nm
    // Source: ESA, 20 June 2005, 'Photon-flux distributions for reference stars', GAIA-EST-TN-00539, issue 1, revision 0, based on D. Katz, priv. comm., 11 May 2005
    // Basic : true
    // Scalar: false
    static const char  *const G2VSTAR_SPECTRUM_NUMBEROFPHOTONSHIGHRESOLUTION() { return "Nature/G2VStar_Spectrum_NumberOfPhotonsHighResolution_001.fits"; }

    // Photon flux density N_{\lambda}(\lambda) of an unreddened G8III star (Pickles' star number 076) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)
    // Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863
    // Basic : true
    // Scalar: false
    static const char  *const G8IIISTAR_SPECTRUM_NUMBEROFPHOTONS() { return "Nature/G8IIIStar_Spectrum_NumberOfPhotons_001.fits"; }

    // Hubble constant (uncertainty is 0.80 km s^-1 Mpc^-1)
    // Source: C.L. Bennett, et al., 1 October 2013, 'Nine-Year Wilkinson Microwave Anisotropy Probe (WMAP) Observations: Final Maps and Results', Astrophysical Journal Supplement, Volume 208, 20
    // Basic : true
    // Scalar: true
    static const double HUBBLE_CONSTANT = 69.32; // [km s^-1 Mpc^-1]

    class IAU2000A {
    public:
        // The IAU 2000A Precession-Nutation Model, developed by Mathews et al. (2002; MHB), is given by a series for nutation in longitude (\Delta\psi) and obliquity (\Delta\epsilon) - referred to the mean ecliptic of date, with time measured in Julian centuries of TDB from epoch J2000.0 - plus the contribution of the corrections to the IAU 1976 precession rates, plus the frame bias in longitude and obliquity. The 'total nutation' includes all components, with the exception of the free core nutation (FCN). The nutation series - providing the direction of the celestial pole in the GCRS with an accuracy of 0.2 mas - includes N_k = 678 luni-solar terms and N_k = 687 planetary terms, which are expressed as 'in-phase' components (A_k, A^\prime_k, B_k, and B^\prime_k) and 'out-of-phase' components (A^\prime\prime_k, A^\prime\prime\prime_k, B^\prime\prime_k, and B^\prime\prime\prime_k) with their time variations: \Delta\psi = Sum_{k=1}^{N_k} (A_k + A^\prime_k * t) * SIN(ARGUMENT) + (A^\prime\prime_k + A^\prime\prime\prime_k * t) * COS(ARGUMENT) and \Delta\epsilon = Sum_{k=1}^{N_k} (B_k + B^\prime_k * t) * COS(ARGUMENT) + (B^\prime\prime_k + B^\prime\prime\prime_k * t) * SIN(ARGUMENT). Each of the N_k = 678 luni-solar terms in the nutation series is characterised by a set of five integers N_j which determines the ARGUMENT for the term as a linear combination of the five Fundamental Arguments F_j, namely the Delaunay variables (l = mean anomaly of the Moon, l^\prime = mean anomaly of the Sun, F = L - \Omega [with l the mean longitude of the Moon], D = mean elongation of the Moon from the Sun, \Omega = mean longitude of the ascending node the Moon): ARGUMENT = Sum_{j=1}^{5} N_j * F_j, where the values (N_1, ..., N_5) of the multipliers characterise the term. The F_j are functions of time, and the angular frequency of the nutation described by the term is given by \omega = d(ARGUMENT) / dt. The N_k = 687 planetary nutation terms differ from the luni-solar terms described above only in that ARGUMENT = SUM_{j=1}^{14} N^\prime_j * F^\prime_j, where F^\prime_6 through F^\prime_13 are the mean longitudes of the planets Mercury through Neptune (including F^\prime_8 for the Earth), and F^\prime_14 is the general precession in longitude p_a. Over time scales involved in nutation studies, the frequency \omega is effectively time-independent, and one may write, for the k-th term in the nutation series, ARGUMENT = \omega_k + \alpha_k. This parameter provides the j = 1, ..., 5 argument multipliers N_j for the N_k = 678 luni-solar terms. These multipliers are dimensionless
        // Source: Mathews, P.M., Herring, T.A., Buffett, B.A., 2002, 'Modeling of nutation-precession: new nutation series for non-rigid Earth, and insights into the Earth's interior', Journal of Geophysical Research (Solid Earth), Volume 107, Issue B4, 2068, 2002JGRB..107.2068M. Reference document: G. Petit, B. Luzum, 21 October 2010,  'IERS Conventions (2010)', IERS Technical Note 36, Chapter 5 (http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html)
        // Basic : true
        // Scalar: false
        static const char  *const PRECESSIONNUTATIONMODEL_ARGUMENTMULTIPLIERS_LUNISOLAR() { return "Nature/IAU2000A/PrecessionNutationModel_ArgumentMultipliers_LuniSolar_001.fits"; }

        // The IAU 2000A Precession-Nutation Model, developed by Mathews et al. (2002; MHB), is given by a series for nutation in longitude (\Delta\psi) and obliquity (\Delta\epsilon) - referred to the mean ecliptic of date, with time measured in Julian centuries of TDB from epoch J2000.0 - plus the contribution of the corrections to the IAU 1976 precession rates, plus the frame bias in longitude and obliquity. The 'total nutation' includes all components, with the exception of the free core nutation (FCN). The nutation series - providing the direction of the celestial pole in the GCRS with an accuracy of 0.2 mas - includes N_k = 678 luni-solar terms and N_k = 687 planetary terms, which are expressed as 'in-phase' components (A_k, A^\prime_k, B_k, and B^\prime_k) and 'out-of-phase' components (A^\prime\prime_k, A^\prime\prime\prime_k, B^\prime\prime_k, and B^\prime\prime\prime_k) with their time variations: \Delta\psi = Sum_{k=1}^{N_k} (A_k + A^\prime_k * t) * SIN(ARGUMENT) + (A^\prime\prime_k + A^\prime\prime\prime_k * t) * COS(ARGUMENT) and \Delta\epsilon = Sum_{k=1}^{N_k} (B_k + B^\prime_k * t) * COS(ARGUMENT) + (B^\prime\prime_k + B^\prime\prime\prime_k * t) * SIN(ARGUMENT). Each of the N_k = 678 luni-solar terms in the nutation series is characterised by a set of five integers N_j which determines the ARGUMENT for the term as a linear combination of the five Fundamental Arguments F_j, namely the Delaunay variables (l = mean anomaly of the Moon, l^\prime = mean anomaly of the Sun, F = L - \Omega [with l the mean longitude of the Moon], D = mean elongation of the Moon from the Sun, \Omega = mean longitude of the ascending node the Moon): ARGUMENT = Sum_{j=1}^{5} N_j * F_j, where the values (N_1, ..., N_5) of the multipliers characterise the term. The F_j are functions of time, and the angular frequency of the nutation described by the term is given by \omega = d(ARGUMENT) / dt. The N_k = 687 planetary nutation terms differ from the luni-solar terms described above only in that ARGUMENT = SUM_{j=1}^{14} N^\prime_j * F^\prime_j, where F^\prime_6 through F^\prime_13 are the mean longitudes of the planets Mercury through Neptune (including F^\prime_8 for the Earth), and F^\prime_14 is the general precession in longitude p_a. Over time scales involved in nutation studies, the frequency \omega is effectively time-independent, and one may write, for the k-th term in the nutation series, ARGUMENT = \omega_k + \alpha_k. This parameter provides the j = 1, ..., 14 argument multipliers N^\prime_j for the N_k = 687 planetary terms. These multipliers are dimensionless
        // Source: Mathews, P.M., Herring, T.A., Buffett, B.A., 2002, 'Modeling of nutation-precession: new nutation series for non-rigid Earth, and insights into the Earth's interior', Journal of Geophysical Research (Solid Earth), Volume 107, Issue B4, 2068, 2002JGRB..107.2068M. Reference document: G. Petit, B. Luzum, 21 October 2010,  'IERS Conventions (2010)', IERS Technical Note 36, Chapter 5 (http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html)
        // Basic : true
        // Scalar: false
        static const char  *const PRECESSIONNUTATIONMODEL_ARGUMENTMULTIPLIERS_PLANETARY() { return "Nature/IAU2000A/PrecessionNutationModel_ArgumentMultipliers_Planetary_001.fits"; }

        // The IAU 2000A Precession-Nutation Model, developed by Mathews et al. (2002; MHB), is given by a series for nutation in longitude (\Delta\psi) and obliquity (\Delta\epsilon) - referred to the mean ecliptic of date, with time measured in Julian centuries of TDB from epoch J2000.0 - plus the contribution of the corrections to the IAU 1976 precession rates, plus the frame bias in longitude and obliquity. The 'total nutation' includes all components, with the exception of the free core nutation (FCN). The nutation series - providing the direction of the celestial pole in the GCRS with an accuracy of 0.2 mas - includes N_k = 678 luni-solar terms and N_k = 687 planetary terms, which are expressed as 'in-phase' components (A_k, A^\prime_k, B_k, and B^\prime_k) and 'out-of-phase' components (A^\prime\prime_k, A^\prime\prime\prime_k, B^\prime\prime_k, and B^\prime\prime\prime_k) with their time variations: \Delta\psi = Sum_{k=1}^{N_k} (A_k + A^\prime_k * t) * SIN(ARGUMENT) + (A^\prime\prime_k + A^\prime\prime\prime_k * t) * COS(ARGUMENT) and \Delta\epsilon = Sum_{k=1}^{N_k} (B_k + B^\prime_k * t) * COS(ARGUMENT) + (B^\prime\prime_k + B^\prime\prime\prime_k * t) * SIN(ARGUMENT). Each of the N_k = 678 luni-solar terms in the nutation series is characterised by a set of five integers N_j which determines the ARGUMENT for the term as a linear combination of the five Fundamental Arguments F_j, namely the Delaunay variables (l = mean anomaly of the Moon, l^\prime = mean anomaly of the Sun, F = L - \Omega [with l the mean longitude of the Moon], D = mean elongation of the Moon from the Sun, \Omega = mean longitude of the ascending node the Moon): ARGUMENT = Sum_{j=1}^{5} N_j * F_j, where the values (N_1, ..., N_5) of the multipliers characterise the term. The F_j are functions of time, and the angular frequency of the nutation described by the term is given by \omega = d(ARGUMENT) / dt. The N_k = 687 planetary nutation terms differ from the luni-solar terms described above only in that ARGUMENT = SUM_{j=1}^{14} N^\prime_j * F^\prime_j, where F^\prime_6 through F^\prime_13 are the mean longitudes of the planets Mercury through Neptune (including F^\prime_8 for the Earth), and F^\prime_14 is the general precession in longitude p_a. Over time scales involved in nutation studies, the frequency \omega is effectively time-independent, and one may write, for the k-th term in the nutation series, ARGUMENT = \omega_k + \alpha_k. This parameter provides the 4 'in-phase' and 4 'out-of-phase' components (nutation coefficients for longitude and obliquity) for the N_k = 678 luni-solar terms in the order A_k, A^\prime_k, B_k, B^\prime_k, A^\prime\prime_k, A^\prime\prime\prime_k, B^\prime\prime_k, B^\prime\prime\prime_k. Units are 10^-3 arcsec and 10^-3 arcsec per century for the coefficients and their time variations, respectively
        // Source: Mathews, P.M., Herring, T.A., Buffett, B.A., 2002, 'Modeling of nutation-precession: new nutation series for non-rigid Earth, and insights into the Earth's interior', Journal of Geophysical Research (Solid Earth), Volume 107, Issue B4, 2068, 2002JGRB..107.2068M. Reference document: G. Petit, B. Luzum, 21 October 2010,  'IERS Conventions (2010)', IERS Technical Note 36, Chapter 5 (http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html)
        // Basic : true
        // Scalar: false
        static const char  *const PRECESSIONNUTATIONMODEL_COEFFICIENTS_LUNISOLAR() { return "Nature/IAU2000A/PrecessionNutationModel_Coefficients_LuniSolar_001.fits"; }

        // The IAU 2000A Precession-Nutation Model, developed by Mathews et al. (2002; MHB), is given by a series for nutation in longitude (\Delta\psi) and obliquity (\Delta\epsilon) - referred to the mean ecliptic of date, with time measured in Julian centuries of TDB from epoch J2000.0 - plus the contribution of the corrections to the IAU 1976 precession rates, plus the frame bias in longitude and obliquity. The 'total nutation' includes all components, with the exception of the free core nutation (FCN). The nutation series - providing the direction of the celestial pole in the GCRS with an accuracy of 0.2 mas - includes N_k = 678 luni-solar terms and N_k = 687 planetary terms, which are expressed as 'in-phase' components (A_k, A^\prime_k, B_k, and B^\prime_k) and 'out-of-phase' components (A^\prime\prime_k, A^\prime\prime\prime_k, B^\prime\prime_k, and B^\prime\prime\prime_k) with their time variations: \Delta\psi = Sum_{k=1}^{N_k} (A_k + A^\prime_k * t) * SIN(ARGUMENT) + (A^\prime\prime_k + A^\prime\prime\prime_k * t) * COS(ARGUMENT) and \Delta\epsilon = Sum_{k=1}^{N_k} (B_k + B^\prime_k * t) * COS(ARGUMENT) + (B^\prime\prime_k + B^\prime\prime\prime_k * t) * SIN(ARGUMENT). Each of the N_k = 678 luni-solar terms in the nutation series is characterised by a set of five integers N_j which determines the ARGUMENT for the term as a linear combination of the five Fundamental Arguments F_j, namely the Delaunay variables (l = mean anomaly of the Moon, l^\prime = mean anomaly of the Sun, F = L - \Omega [with l the mean longitude of the Moon], D = mean elongation of the Moon from the Sun, \Omega = mean longitude of the ascending node the Moon): ARGUMENT = Sum_{j=1}^{5} N_j * F_j, where the values (N_1, ..., N_5) of the multipliers characterise the term. The F_j are functions of time, and the angular frequency of the nutation described by the term is given by \omega = d(ARGUMENT) / dt. The N_k = 687 planetary nutation terms differ from the luni-solar terms described above only in that ARGUMENT = SUM_{j=1}^{14} N^\prime_j * F^\prime_j, where F^\prime_6 through F^\prime_13 are the mean longitudes of the planets Mercury through Neptune (including F^\prime_8 for the Earth), and F^\prime_14 is the general precession in longitude p_a. Over time scales involved in nutation studies, the frequency \omega is effectively time-independent, and one may write, for the k-th term in the nutation series, ARGUMENT = \omega_k + \alpha_k. This parameter provides the 2 'in-phase' and 2 'out-of-phase' components (nutation coefficients for longitude and obliquity) for the N_k = 687 planetary terms in the order A_k, A^\prime\prime_k, B_k, B^\prime\prime_k. Units are 10^-3 arcsec
        // Source: Mathews, P.M., Herring, T.A., Buffett, B.A., 2002, 'Modeling of nutation-precession: new nutation series for non-rigid Earth, and insights into the Earth's interior', Journal of Geophysical Research (Solid Earth), Volume 107, Issue B4, 2068, 2002JGRB..107.2068M. Reference document: G. Petit, B. Luzum, 21 October 2010,  'IERS Conventions (2010)', IERS Technical Note 36, Chapter 5 (http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html)
        // Basic : true
        // Scalar: false
        static const char  *const PRECESSIONNUTATIONMODEL_COEFFICIENTS_PLANETARY() { return "Nature/IAU2000A/PrecessionNutationModel_Coefficients_Planetary_001.fits"; }

    };

    // Galactic longitude of the ascending node of the Galactic plane on the equator of the ICRS (numerical value should be regarded as exact)
    // Source: ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, page 91
    // Basic : true
    // Scalar: true
    static const double ICRS_LONGITUDEOFASCENDINGNODE = 32.93192; // [deg]

    // Declination of the north Galactic pole in the ICRS system (numerical value should be regarded as exact)
    // Source: ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, page 91
    // Basic : true
    // Scalar: true
    static const double ICRS_NORTHGALACTICPOLE_DECLINATION = 27.12825; // [deg]

    // Right ascension of the north Galactic pole in the ICRS system (numerical value should be regarded as exact)
    // Source: ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, page 91
    // Basic : true
    // Scalar: true
    static const double ICRS_NORTHGALACTICPOLE_RIGHTASCENSION = 192.85948; // [deg]

    class INPOP10e {
    public:
        // GM of asteroid 1013 Tombecka (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID1013TOMBECKA_GM = 7.9531899440686237e+06; // [m^3 s^-2]

        // GM of asteroid 1021 Flammario (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID1021FLAMMARIO_GM = 7.3151076468356557e+07; // [m^3 s^-2]

        // GM of asteroid 1036 Ganymed (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID1036GANYMED_GM = 7.7604237885203171e+06; // [m^3 s^-2]

        // GM of asteroid 105 Artemis (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID105ARTEMIS_GM = 4.0425884340051366e+08; // [m^3 s^-2]

        // GM of asteroid 106 Dione (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID106DIONE_GM = 5.1363068033089922e+08; // [m^3 s^-2]

        // GM of asteroid 107 Camilla (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID107CAMILLA_GM = 4.5300167118439350e+08; // [m^3 s^-2]

        // GM of asteroid 109 Felicitas (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID109FELICITAS_GM = 2.1334036970797304e+07; // [m^3 s^-2]

        // GM of asteroid 10 Hygiea (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID10HYGIEA_GM = 5.8390041443806954e+09; // [m^3 s^-2]

        // GM of asteroid 111 Ate (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID111ATE_GM = 5.9580019641313180e+08; // [m^3 s^-2]

        // GM of asteroid 112 Iphigenia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID112IPHIGENIA_GM = 9.1960926871076204e+07; // [m^3 s^-2]

        // GM of asteroid 117 Lomia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID117LOMIA_GM = 8.0437920943187555e+08; // [m^3 s^-2]

        // GM of asteroid 11 Parthenope (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID11PARTHENOPE_GM = 5.0048256226049669e+08; // [m^3 s^-2]

        // GM of asteroid 120 Lachesis (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID120LACHESIS_GM = 2.0019174836054619e+08; // [m^3 s^-2]

        // GM of asteroid 121 Hermione (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID121HERMIONE_GM = 2.0870018202750795e+09; // [m^3 s^-2]

        // GM of asteroid 126 Velleda (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID126VELLEDA_GM = 2.2017513093567601e+07; // [m^3 s^-2]

        // GM of asteroid 128 Nemesis (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID128NEMESIS_GM = 9.1024212903081425e+08; // [m^3 s^-2]

        // GM of asteroid 12 Victoria (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID12VICTORIA_GM = 6.9482386231968752e+07; // [m^3 s^-2]

        // GM of asteroid 132 Aethra (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID132AETHRA_GM = 1.9253478609675224e+07; // [m^3 s^-2]

        // GM of asteroid 134 Sophrosyne (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID134SOPHROSYNE_GM = 1.3455129358325903e+08; // [m^3 s^-2]

        // GM of asteroid 135 Hertha (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID135HERTHA_GM = 1.2167072819223394e+08; // [m^3 s^-2]

        // GM of asteroid 139 Juewa (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID139JUEWA_GM = 2.8269837868259157e+08; // [m^3 s^-2]

        // GM of asteroid 13 Egeria (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID13EGERIA_GM = 6.2548618644966828e+08; // [m^3 s^-2]

        // GM of asteroid 141 Lumen (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID141LUMEN_GM = 5.5025545915405380e+08; // [m^3 s^-2]

        // GM of asteroid 156 Xanthippe (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID156XANTHIPPE_GM = 4.3257295425004363e+08; // [m^3 s^-2]

        // GM of asteroid 15 Eunomia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID15EUNOMIA_GM = 2.1020267614940911e+09; // [m^3 s^-2]

        // GM of asteroid 164 Eva (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID164EVA_GM = 1.9354802800126533e+05; // [m^3 s^-2]

        // GM of asteroid 168 Sibylla (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID168SIBYLLA_GM = 7.9887495652663551e+08; // [m^3 s^-2]

        // GM of asteroid 1694 Kaiser (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID1694KAISER_GM = 5.1974702708939677e+06; // [m^3 s^-2]

        // GM of asteroid 16 Psyche (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID16PSYCHE_GM = 1.6739220218687147e+09; // [m^3 s^-2]

        // GM of asteroid 171 Ophelia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID171OPHELIA_GM = 3.8845612002754150e+08; // [m^3 s^-2]

        // GM of asteroid 172 Baucis (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID172BAUCIS_GM = 5.9473568787444189e+07; // [m^3 s^-2]

        // GM of asteroid 173 Ino (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID173INO_GM = 8.9487138876601197e+08; // [m^3 s^-2]

        // GM of asteroid 179 Klytaemnestra (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID179KLYTAEMNESTRA_GM = 1.1462526623299412e+08; // [m^3 s^-2]

        // GM of asteroid 17 Thetis (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID17THETIS_GM = 4.8917509952675321e+08; // [m^3 s^-2]

        // GM of asteroid 185 Eunike (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID185EUNIKE_GM = 5.8481152014919869e+08; // [m^3 s^-2]

        // GM of asteroid 187 Lamberta (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID187LAMBERTA_GM = 1.1232253601955330e+07; // [m^3 s^-2]

        // GM of asteroid 194 Prokne (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID194PROKNE_GM = 7.4333048018458559e+08; // [m^3 s^-2]

        // GM of asteroid 19 Fortuna (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID19FORTUNA_GM = 6.4925120044129564e+08; // [m^3 s^-2]

        // GM of asteroid 1 Ceres (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID1CERES_GM = 6.2012183942528447e+10; // [m^3 s^-2]

        // GM of asteroid 200 Dynamene (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID200DYNAMENE_GM = 7.6156064781391007e+07; // [m^3 s^-2]

        // GM of asteroid 204 Kallisto (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID204KALLISTO_GM = 2.8036593061651360e+07; // [m^3 s^-2]

        // GM of asteroid 209 Dido (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID209DIDO_GM = 1.0005158700123591e+09; // [m^3 s^-2]

        // GM of asteroid 20 Massalia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID20MASSALIA_GM = 3.8450306541675409e+08; // [m^3 s^-2]

        // GM of asteroid 210 Isabella (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID210ISABELLA_GM = 1.5915077609293021e+08; // [m^3 s^-2]

        // GM of asteroid 211 Isolda (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID211ISOLDA_GM = 7.1779392920251512e+08; // [m^3 s^-2]

        // GM of asteroid 212 Medea (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID212MEDEA_GM = 3.1495425151492577e+08; // [m^3 s^-2]

        // GM of asteroid 217 Eudora (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID217EUDORA_GM = 7.1074477882268534e+07; // [m^3 s^-2]

        // GM of asteroid 21 Lutetia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID21LUTETIA_GM = 1.1503721169099261e+08; // [m^3 s^-2]

        // GM of asteroid 22 Kalliope (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID22KALLIOPE_GM = 1.1112751634270690e+09; // [m^3 s^-2]

        // GM of asteroid 234 Barbara (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID234BARBARA_GM = 2.0492014129900661e+07; // [m^3 s^-2]

        // GM of asteroid 23 Thalia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID23THALIA_GM = 2.0499213757716355e+08; // [m^3 s^-2]

        // GM of asteroid 247 Eukrate (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID247EUKRATE_GM = 5.9420762173955001e+08; // [m^3 s^-2]

        // GM of asteroid 250 Bettina (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID250BETTINA_GM = 1.2408181646515161e+08; // [m^3 s^-2]

        // GM of asteroid 253 Mathilde (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID253MATHILDE_GM = 8.4021633178130707e+07; // [m^3 s^-2]

        // GM of asteroid 25 Phocaea (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID25PHOCAEA_GM = 1.0366196678124746e+08; // [m^3 s^-2]

        // GM of asteroid 266 Aline (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID266ALINE_GM = 3.1756085199183367e+08; // [m^3 s^-2]

        // GM of asteroid 26 Proserpina (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID26PROSERPINA_GM = 2.0834242686152029e+08; // [m^3 s^-2]

        // GM of asteroid 29 Amphitrite (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID29AMPHITRITE_GM = 9.5914046612273862e+08; // [m^3 s^-2]

        // GM of asteroid 2 Pallas (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID2PALLAS_GM = 1.3623519829068894e+10; // [m^3 s^-2]

        // GM of asteroid 304 Olga (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID304OLGA_GM = 7.6417754608396792e+07; // [m^3 s^-2]

        // GM of asteroid 308 Polyxo (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID308POLYXO_GM = 5.1288139004370254e+08; // [m^3 s^-2]

        // GM of asteroid 313 Chaldaea (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID313CHALDAEA_GM = 1.3568618932908778e+08; // [m^3 s^-2]

        // GM of asteroid 31 Euphrosyne (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID31EUPHROSYNE_GM = 1.7562517866270849e+09; // [m^3 s^-2]

        // GM of asteroid 322 Phaeo (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID322PHAEO_GM = 8.6933740062976113e+07; // [m^3 s^-2]

        // GM of asteroid 324 Bamberga (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID324BAMBERGA_GM = 6.3292657018921930e+08; // [m^3 s^-2]

        // GM of asteroid 335 Roberta (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID335ROBERTA_GM = 1.7274312784818528e+08; // [m^3 s^-2]

        // GM of asteroid 336 Lacadiera (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID336LACADIERA_GM = 8.1456992030180356e+07; // [m^3 s^-2]

        // GM of asteroid 33 Polyhymnia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID33POLYHYMNIA_GM = 2.8960743674231378e+08; // [m^3 s^-2]

        // GM of asteroid 346 Hermentaria (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID346HERMENTARIA_GM = 2.9556035531795935e+08; // [m^3 s^-2]

        // GM of asteroid 34 Circe (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID34CIRCE_GM = 1.9273254973801046e+08; // [m^3 s^-2]

        // GM of asteroid 350 Ornamenta (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID350ORNAMENTA_GM = 4.0527274829271427e+08; // [m^3 s^-2]

        // GM of asteroid 356 Liguria (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID356LIGURIA_GM = 5.5379028002660318e+08; // [m^3 s^-2]

        // GM of asteroid 36 Atalante (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID36ATALANTE_GM = 2.8796815326611882e+08; // [m^3 s^-2]

        // GM of asteroid 37 Fides (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID37FIDES_GM = 3.1097050608961499e+08; // [m^3 s^-2]

        // GM of asteroid 381 Myrrha (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID381MYRRHA_GM = 4.2872473963951988e+08; // [m^3 s^-2]

        // GM of asteroid 387 Aquitania (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID387AQUITANIA_GM = 2.4837629612901240e+08; // [m^3 s^-2]

        // GM of asteroid 388 Charybdis (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID388CHARYBDIS_GM = 3.2080410583920504e+08; // [m^3 s^-2]

        // GM of asteroid 3 Juno (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID3JUNO_GM = 1.5650376033611044e+09; // [m^3 s^-2]

        // GM of asteroid 404 Arsinoe (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID404ARSINOE_GM = 8.3335002699994517e+07; // [m^3 s^-2]

        // GM of asteroid 405 Thia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID405THIA_GM = 1.8291269394614959e+08; // [m^3 s^-2]

        // GM of asteroid 409 Aspasia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID409ASPASIA_GM = 1.4485207257722234e+05; // [m^3 s^-2]

        // GM of asteroid 410 Chloris (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID410CHLORIS_GM = 4.6130313447746411e+08; // [m^3 s^-2]

        // GM of asteroid 419 Aurelia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID419AURELIA_GM = 2.1890666809971182e+08; // [m^3 s^-2]

        // GM of asteroid 420 Bertholda (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID420BERTHOLDA_GM = 6.8901048902653521e+08; // [m^3 s^-2]

        // GM of asteroid 442 Eichsfeldia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID442EICHSFELDIA_GM = 7.2696133676778026e+07; // [m^3 s^-2]

        // GM of asteroid 444 Gyptis (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID444GYPTIS_GM = 7.0727216281453015e+08; // [m^3 s^-2]

        // GM of asteroid 445 Edna (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID445EDNA_GM = 1.6203293575514477e+08; // [m^3 s^-2]

        // GM of asteroid 44 Nysa (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID44NYSA_GM = 8.6199506289667382e+07; // [m^3 s^-2]

        // GM of asteroid 451 Patientia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID451PATIENTIA_GM = 1.9885686160983943e+09; // [m^3 s^-2]

        // GM of asteroid 455 Bruchsalia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID455BRUCHSALIA_GM = 1.4702122822900930e+08; // [m^3 s^-2]

        // GM of asteroid 469 Argentina (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID469ARGENTINA_GM = 5.2207383300069546e+05; // [m^3 s^-2]

        // GM of asteroid 46 Hestia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID46HESTIA_GM = 4.6782985130525508e+08; // [m^3 s^-2]

        // GM of asteroid 481 Emita (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID481EMITA_GM = 1.1562903987899217e+08; // [m^3 s^-2]

        // GM of asteroid 488 Kreusa (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID488KREUSA_GM = 6.8444731733902973e+08; // [m^3 s^-2]

        // GM of asteroid 4 Vesta (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID4VESTA_GM = 1.7288981939943646e+10; // [m^3 s^-2]

        // GM of asteroid 503 Evelyn (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID503EVELYN_GM = 1.3326002607937277e+08; // [m^3 s^-2]

        // GM of asteroid 505 Cava (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID505CAVA_GM = 2.6643963557510264e+08; // [m^3 s^-2]

        // GM of asteroid 511 Davida (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID511DAVIDA_GM = 2.4220172669698846e+09; // [m^3 s^-2]

        // GM of asteroid 516 Amherstia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID516AMHERSTIA_GM = 4.6419532383144037e+07; // [m^3 s^-2]

        // GM of asteroid 51 Nemausa (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID51NEMAUSA_GM = 1.1850701657405471e+06; // [m^3 s^-2]

        // GM of asteroid 52 Europa (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID52EUROPA_GM = 1.4176132635840474e+09; // [m^3 s^-2]

        // GM of asteroid 532 Herculina (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID532HERCULINA_GM = 1.5330914367975053e+09; // [m^3 s^-2]

        // GM of asteroid 53 Kalypso (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID53KALYPSO_GM = 2.6631974354533394e+08; // [m^3 s^-2]

        // GM of asteroid 54 Alexandra (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID54ALEXANDRA_GM = 1.1137607808140530e+09; // [m^3 s^-2]

        // GM of asteroid 55 Pandora (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID55PANDORA_GM = 7.2565502740078594e+07; // [m^3 s^-2]

        // GM of asteroid 568 Cheruskia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID568CHERUSKIA_GM = 1.6092033080277349e+08; // [m^3 s^-2]

        // GM of asteroid 569 Misa (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID569MISA_GM = 9.4974536474941922e+07; // [m^3 s^-2]

        // GM of asteroid 56 Melete (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID56MELETE_GM = 3.5510132170246012e+08; // [m^3 s^-2]

        // GM of asteroid 583 Klotilde (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID583KLOTILDE_GM = 1.3306434325297216e+08; // [m^3 s^-2]

        // GM of asteroid 584 Semiramis (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID584SEMIRAMIS_GM = 3.8506494714933353e+07; // [m^3 s^-2]

        // GM of asteroid 591 Irmgard (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID591IRMGARD_GM = 1.5125998853964005e+06; // [m^3 s^-2]

        // GM of asteroid 593 Titania (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID593TITANIA_GM = 1.0449214411384316e+08; // [m^3 s^-2]

        // GM of asteroid 595 Polyxena (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID595POLYXENA_GM = 3.1721169220784370e+08; // [m^3 s^-2]

        // GM of asteroid 599 Luisa (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID599LUISA_GM = 6.6785968807632834e+07; // [m^3 s^-2]

        // GM of asteroid 602 Marianna (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID602MARIANNA_GM = 2.5659694460527015e+08; // [m^3 s^-2]

        // GM of asteroid 618 Elfriede (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID618ELFRIEDE_GM = 5.1050504107617317e+05; // [m^3 s^-2]

        // GM of asteroid 61 Danae (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID61DANAE_GM = 4.1297923333635475e+07; // [m^3 s^-2]

        // GM of asteroid 626 Notburga (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID626NOTBURGA_GM = 2.5001067131733260e+08; // [m^3 s^-2]

        // GM of asteroid 62 Erato (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID62ERATO_GM = 2.1232338569167166e+08; // [m^3 s^-2]

        // GM of asteroid 63 Ausonia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID63AUSONIA_GM = 4.5857887479153414e+05; // [m^3 s^-2]

        // GM of asteroid 65 Cybele (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID65CYBELE_GM = 5.5875034960511513e+08; // [m^3 s^-2]

        // GM of asteroid 667 Denise (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID667DENISE_GM = 7.8762291169638808e+07; // [m^3 s^-2]

        // GM of asteroid 675 Ludmilla (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID675LUDMILLA_GM = 7.0108090353291896e+08; // [m^3 s^-2]

        // GM of asteroid 67 Asia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID67ASIA_GM = 4.7960213071519115e+07; // [m^3 s^-2]

        // GM of asteroid 690 Wratislavia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID690WRATISLAVIA_GM = 5.9712951363080574e+08; // [m^3 s^-2]

        // GM of asteroid 694 Ekard (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID694EKARD_GM = 2.8164391089152126e+03; // [m^3 s^-2]

        // GM of asteroid 6 Hebe (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID6HEBE_GM = 9.4013859632630443e+08; // [m^3 s^-2]

        // GM of asteroid 704 Interamnia (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID704INTERAMNIA_GM = 2.5503367929164454e+09; // [m^3 s^-2]

        // GM of asteroid 70 Panopaea (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID70PANOPAEA_GM = 2.7708718966786989e+08; // [m^3 s^-2]

        // GM of asteroid 718 Erida (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID718ERIDA_GM = 4.0557482554503977e+07; // [m^3 s^-2]

        // GM of asteroid 735 Marghanna (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID735MARGHANNA_GM = 1.0038522798255246e+08; // [m^3 s^-2]

        // GM of asteroid 739 Mandeville (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID739MANDEVILLE_GM = 3.0413247983619986e+08; // [m^3 s^-2]

        // GM of asteroid 747 Winchester (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID747WINCHESTER_GM = 9.6009858605926026e+07; // [m^3 s^-2]

        // GM of asteroid 751 Faina (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID751FAINA_GM = 3.2994347243673589e+08; // [m^3 s^-2]

        // GM of asteroid 752 Sulamitis (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID752SULAMITIS_GM = 6.0508535365628257e+07; // [m^3 s^-2]

        // GM of asteroid 75 Eurydike (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID75EURYDIKE_GM = 4.2715759272180030e+07; // [m^3 s^-2]

        // GM of asteroid 78 Diana (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID78DIANA_GM = 3.4003238496349289e+08; // [m^3 s^-2]

        // GM of asteroid 790 Pretoria (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID790PRETORIA_GM = 8.6438330145248586e+08; // [m^3 s^-2]

        // GM of asteroid 791 Ani (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID791ANI_GM = 2.7128482223705070e+08; // [m^3 s^-2]

        // GM of asteroid 7 Iris (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID7IRIS_GM = 8.3629208510440051e+08; // [m^3 s^-2]

        // GM of asteroid 804 Hispania (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID804HISPANIA_GM = 1.3669667630778693e+08; // [m^3 s^-2]

        // GM of asteroid 80 Sappho (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID80SAPPHO_GM = 1.1784221437378471e+08; // [m^3 s^-2]

        // GM of asteroid 814 Tauris (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID814TAURIS_GM = 3.2159462197573089e+08; // [m^3 s^-2]

        // GM of asteroid 85 Io (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID85IO_GM = 2.9069348822079933e+08; // [m^3 s^-2]

        // GM of asteroid 88 Thisbe (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID88THISBE_GM = 9.4067567843322719e+08; // [m^3 s^-2]

        // GM of asteroid 8 Flora (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID8FLORA_GM = 4.4557515312994348e+08; // [m^3 s^-2]

        // GM of asteroid 914 Palisana (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID914PALISANA_GM = 5.7331462009756002e+07; // [m^3 s^-2]

        // GM of asteroid 93 Minerva (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID93MINERVA_GM = 5.0039749795125115e+08; // [m^3 s^-2]

        // GM of asteroid 949 Hel (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID949HEL_GM = 8.0964451809677710e+07; // [m^3 s^-2]

        // GM of asteroid 94 Aurora (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID94AURORA_GM = 1.9949958990308900e+09; // [m^3 s^-2]

        // GM of asteroid 97 Klotho (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID97KLOTHO_GM = 6.1873651052720885e+06; // [m^3 s^-2]

        // GM of asteroid 9 Metis (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROID9METIS_GM = 5.5771967924306053e+08; // [m^3 s^-2]

        // GM of the asteroid ring (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double ASTEROIDRING_GM = 8.9697172580075571e+09; // [m^3 s^-2]

        // Barycentric distance (orbital semi-major axis) of the asteroid ring (INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : true
        // Scalar: true
        static const double ASTEROIDRING_ORBITALSEMIMAJORAXIS = 3.1477080248116020e+00; // [au]

        // Astronomical unit (au) length in m
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : true
        // Scalar: true
        static const double ASTRONOMICALUNIT_METER = 1.4959787070000000e+11; // [m]

        // GM of the Earth-system (TCB-compatible value). The gravitational constant includes the contribution of its satellite, the Moon
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double EARTHSYSTEM_GM = 4.0350325101102718e+14; // [m^3 s^-2]

        // Ratio of Earth to Moon mass (INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : true
        // Scalar: true
        static const double EARTHTOMOON_MASSRATIO = 8.1300570000000000e+01;

        // Mean equatorial radius of the Earth (INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : true
        // Scalar: true
        static const double EARTH_EQUATORIALRADIUS = 6.3781366988942710e+06; // [m]

        // Geocentric gravitational constant (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double EARTH_GM = 3.9860045081157502e+14; // [m^3 s^-2]

        // GM of the Jupiter-system (TCB-compatible value). The gravitational constant includes the contribution of its satellites
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double JUPITERSYSTEM_GM = 1.2671276453465735e+17; // [m^3 s^-2]

        // GM of Jupiter (TCB-compatible value; this is a low-accuracy parameter used in the relativistic model, ignoring the contribution of planetary satellites)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double JUPITER_GM = 1.2671276453465734e+17; // [m^3 s^-2]

        // GM of the Mars-system (TCB-compatible value). The gravitational constant includes the contribution of its satellites
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double MARSSYSTEM_GM = 4.2828375886337909e+13; // [m^3 s^-2]

        // GM of Mars (TCB-compatible value; this is a low-accuracy parameter used in the relativistic model, ignoring the contribution of planetary satellites)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double MARS_GM = 4.2828375886337906e+13; // [m^3 s^-2]

        // GM of the Mercury(-system) (TCB-compatible value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double MERCURYSYSTEM_GM = 2.2032080834196276e+13; // [m^3 s^-2]

        // GM of Mercury (TCB-compatible value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double MERCURY_GM = 2.2032080834196277e+13; // [m^3 s^-2]

        // Ratio of Moon to Earth mass (INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double MOONTOEARTH_MASSRATIO = 0.0123000368;

        // Mean equatorial radius of the Moon (INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : true
        // Scalar: true
        static const double MOON_EQUATORIALRADIUS = 1.7380000269480340e+06; // [m]

        // Selenocentric gravitational constant (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double MOON_GM = 4.9028001994521693e+12; // [m^3 s^-2]

        // GM of the Neptune-system (TCB-compatible value). The gravitational constant includes the contribution of its satellites
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double NEPTUNESYSTEM_GM = 6.8365271283644811e+15; // [m^3 s^-2]

        // GM of Neptune (TCB-compatible value; this is a low-accuracy parameter used in the relativistic model, ignoring the contribution of planetary satellites)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double NEPTUNE_GM = 6.8365271283644810e+15; // [m^3 s^-2]

        // General relativistic standard PPN parameter \beta (INPOP10e value, fixed to 1)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : true
        // Scalar: true
        static const double PPN_BETA = 1.;

        // General relativistic standard PPN parameter \gamma (INPOP10e value, fixed to 1)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : true
        // Scalar: true
        static const double PPN_GAMMA = 1.;

        // GM of the Pluto-system (TCB-compatible value). The gravitational constant includes the contribution of its satellite, Charon
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double PLUTOSYSTEM_GM = 9.7178245029026624e+11; // [m^3 s^-2]

        // GM of Pluto (TCB-compatible value; this is a low-accuracy parameter used in the relativistic model, ignoring the contribution of Charon)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double PLUTO_GM = 9.7178245029026624e+11; // [m^3 s^-2]

        // GM of the Saturn-system (TCB-compatible value). The gravitational constant includes the contribution of its satellites
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double SATURNSYSTEM_GM = 3.7940585442640140e+16; // [m^3 s^-2]

        // GM of Saturn (TCB-compatible value; this is a low-accuracy parameter used in the relativistic model, ignoring the contribution of planetary satellites)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double SATURN_GM = 3.7940585442640144e+16; // [m^3 s^-2]

        // Ratio of Sun to Earth-system mass (INPOP10e value). The planetary mass includes the contribution of its satellite, the Moon
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double SUNTOEARTHSYSTEM_MASSRATIO = 328900.552289;

        // Ratio of Sun to Jupiter-system mass (INPOP10e value). The planetary mass includes the contribution of its satellites
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double SUNTOJUPITERSYSTEM_MASSRATIO = 1047.348644;

        // Ratio of Sun to Mars-system mass (INPOP10e value). The planetary mass includes the contribution of its satellites
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double SUNTOMARSSYSTEM_MASSRATIO = 3098704.;

        // Ratio of Sun to Mercury(-system) mass (INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double SUNTOMERCURYSYSTEM_MASSRATIO = 6023600.;

        // Ratio of Sun to Neptune-system mass (INPOP10e value). The planetary mass includes the contribution of its satellites
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double SUNTONEPTUNESYSTEM_MASSRATIO = 19412.26;

        // Ratio of Sun to Pluto-system mass (INPOP10e value). The 'planetary' mass includes the contribution of its satellite, Charon
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double SUNTOPLUTOSYSTEM_MASSRATIO = 136566000.;

        // Ratio of Sun to Saturn-system mass (INPOP10e value). The planetary mass includes the contribution of its satellites
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double SUNTOSATURNSYSTEM_MASSRATIO = 3497.902;

        // Ratio of Sun to Uranus-system mass (INPOP10e value). The planetary mass includes the contribution of its satellites
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double SUNTOURANUSSYSTEM_MASSRATIO = 22902.98;

        // Ratio of Sun to Venus(-system) mass (INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double SUNTOVENUSSYSTEM_MASSRATIO = 408523.72;

        // Mean equatorial radius of the Sun (INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : true
        // Scalar: true
        static const double SUN_EQUATORIALRADIUS = 6.9600001079161780e+08; // [m]

        // Heliocentric gravitational constant (TCB-compatible value; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double SUN_GM = 1.3271244210789467e+20; // [m^3 s^-2]

        // Dynamical form-factor of the Sun (Stokes' second-degree zonal harmonic of the solar potential; INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : true
        // Scalar: true
        static const double SUN_JSUB2 = 1.8000000000000000e-07;

        // Declination \delta_0 of the north pole of rotation of the Sun (INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : true
        // Scalar: true
        static const double SUN_NORTHROTATIONALPOLE_DECLINATION = 6.3870000000000000e+01; // [deg]

        // Right ascension \alpha_0 of the north pole of rotation of the Sun (INPOP10e value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : true
        // Scalar: true
        static const double SUN_NORTHROTATIONALPOLE_RIGHTASCENSION = 2.8613000000000000e+02; // [deg]

        // GM of the Uranus-system (TCB-compatible value). The gravitational constant includes the contribution of its satellites
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double URANUSSYSTEM_GM = 5.7945490985393442e+15; // [m^3 s^-2]

        // GM of Uranus (TCB-compatible value; this is a low-accuracy parameter used in the relativistic model, ignoring the contribution of planetary satellites)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double URANUS_GM = 5.7945490985393440e+15; // [m^3 s^-2]

        // GM of the Venus(-system) (TCB-compatible value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double VENUSSYSTEM_GM = 3.2485859679756975e+14; // [m^3 s^-2]

        // GM of Venus (TCB-compatible value)
        // Source: A. Fienga, H. Manche, J. Laskar, M. Gastineau, A. Verma, 7 March 2013, 'DPAC INPOP final release: INPOP10e', GAIA-CA-TN-IMC-AF-002, issue 1, revision 0
        // Basic : false
        // Scalar: true
        static const double VENUS_GM = 3.2485859679756975e+14; // [m^3 s^-2]

    };

    // Inverse of the fine structure constant. Note: best-measured value equals 137.035999139 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))
    // Basic : false
    // Scalar: true
    static const double INVERSEFINESTRUCTURE_CONSTANT = 137.035999142;

    // Number of days per Julian century
    // Source: IAU definition
    // Basic : true
    // Scalar: true
    static const double JULIANCENTURY_JULIANYEAR = 36525.; // [day]

    // Julian date of the standard epoch J2000 = J2000.0, i.e., calendar date 2000 January 1, 12:00:00 h = 2000-01-01T12:00:00 TT
    // Source: Definition (e.g., ESA, 1997, 'The Hipparcos and Tycho Catalogues', Volume 1, page 27)
    // Basic : true
    // Scalar: true
    static const double JULIANDATE_J2000 = 2451545.0; // [JD]

    // Julian date of the Gaia-specific reference epoch J2010 = J2010.0 = JD2455197.5 = 2010-01-01T00:00:00
    // Source: U. Bastian, 5 July 2007, 'Reference systems, conventions, and notations for Gaia', GAIA-CA-SP-ARI-BAS-003, issue 6, revision 1, Section 3.5
    // Basic : false
    // Scalar: true
    static const double JULIANDATE_J2010 = 2455197.5; // [JD]

    // Number of days per Julian year
    // Source: IAU definition
    // Basic : false
    // Scalar: true
    static const double JULIANYEAR_DAY = 365.25; // [day]

    // Astrometric signature of the Sun induced by the Jupiter system for an observer located at a distance of 10 pc from the Sun
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.7, Equation 1.22, page 11
    // Basic : false
    // Scalar: true
    static const double JUPITERSYSTEM_ASTROMETRICSIGNATURE_10PARSEC = 497.; // [10^-6 arcsec]

    // Jupiter-system mass (IAU 2009 CBE value). The planetary mass includes the contribution of its satellites
    // Basic : false
    // Scalar: true
    static const double JUPITERSYSTEM_MASS = 1.89858e+27; // [kg]

    // Mean orbital eccentricity of Jupiter, at the standard epoch J2000.0. The mean orbital eccentricity is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double JUPITERSYSTEM_ORBITALECCENTRICITY_J2000 = 0.04838624;

    // Mean orbital inclination of Jupiter, at the standard epoch J2000.0. The mean orbital inclination is is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double JUPITERSYSTEM_ORBITALINCLINATION_J2000 = 1.30439695; // [deg]

    // Sidereal orbital period
    // Source: Values derived from the mean longitude rates in Table 5.8.1 in P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, page 316; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double JUPITERSYSTEM_ORBITALPERIOD = 11.862615; // [yr]

    // Mean orbital semi-major axis of Jupiter, at the standard epoch J2000.0. The mean orbital semi-major axis is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double JUPITERSYSTEM_ORBITALSEMIMAJORAXIS_J2000 = 5.20288700; // [au]

    // Radial-velocity amplitude of the Sun induced by the Jupiter system for 'an edge-on observer' (i.e., an observer in the orbital plane of the Jupiter system)
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.4, Equation 1.18, page 9
    // Basic : false
    // Scalar: true
    static const double JUPITERSYSTEM_RADIALVELOCITYSIGNATURE = 12.5; // [m s^-1]

    // Radius of the smallest hypothetical sphere around Jupiter which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double JUPITER_ENCOMPASSINGSPHERERADIUS = 7.14917e+07; // [m]

    // Equatorial radius of Jupiter
    // Basic : false
    // Scalar: true
    static const double JUPITER_EQUATORIALRADIUS = 7.14917e+07; // [m]

    // Nominal equatorial radius of Jupiter (one-bar value), in SI units. This nominal value shall be understood as conversion factor only
    // Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)
    // Basic : true
    // Scalar: true
    static const double JUPITER_EQUATORIALRADIUS_NOMINAL = 7.14920e+07; // [m]

    // Geometrical flattening factor f of Jupiter (f = (a-b)/a)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double JUPITER_FLATTENING = 6.487440e-02;

    // Maximum reduction of the solar flux for an observer external to the solar system during a transit of Jupiter
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.2, Equation 2.4, page 14
    // Basic : false
    // Scalar: true
    static const double JUPITER_FLUXREDUCTION_MAXIMUM = 1.010; // [%]

    // Nominal GM of Jupiter, in SI units. This nominal value shall be understood as conversion factor only
    // Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)
    // Basic : true
    // Scalar: true
    static const double JUPITER_GM_NOMINAL = 1.26686530e+17; // [m^3 s^-2]

    // Geometric albedo of Jupiter (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double JUPITER_GEOMETRICALBEDO = 0.52;

    // Dynamical form-factor of Jupiter (oblateness or Stokes' second-degree zonal harmonic of the potential)
    // Source: P.R. Weissman, L.-A. McFadden, T.V. Johnson (eds.), 1999, 'Encyclopedia of the Solar System (first edition)', Academic Press, page 342
    // Basic : true
    // Scalar: true
    static const double JUPITER_JSUB2 = 0.014697;

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of Jupiter
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double JUPITER_LIGHTDEFLECTION_LIMB = 16635.; // [10^-6 arcsec]

    // Mass of Jupiter (do not use for high-precision (orbit) calculations)
    // Source: R.A. Jacobson, 2005, 'Jovian Satellite ephemeris - JUP230', priv. comm.; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double JUPITER_MASS = 1.898130e+27; // [kg]

    // Mean mass density of Jupiter
    // Basic : false
    // Scalar: true
    static const double JUPITER_MASSDENSITY_MEAN = 1.3262; // [g cm^-3]

    // IAU-recommended value for the declination \delta_0 of the north pole of rotation of Jupiter. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double JUPITER_NORTHROTATIONALPOLE_DECLINATION = 64.495303; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the declination \delta_0 of the north pole of rotation of Jupiter. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : false
    // Scalar: true
    static const double JUPITER_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE = 0.000000066064; // [deg day^-1]

    // IAU-recommended value for the right ascension \alpha_0 of the north pole of rotation of Jupiter. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double JUPITER_NORTHROTATIONALPOLE_RIGHTASCENSION = 268.056595; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the right ascension \alpha_0 of the north pole of rotation of Jupiter. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : false
    // Scalar: true
    static const double JUPITER_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE = -0.000000177933; // [deg day^-1]

    // Polar radius of Jupiter
    // Basic : false
    // Scalar: true
    static const double JUPITER_POLARRADIUS = 6.68537e+07; // [m]

    // Nominal polar radius of Jupiter (one-bar value), in SI units. This nominal value shall be understood as conversion factor only
    // Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)
    // Basic : true
    // Scalar: true
    static const double JUPITER_POLARRADIUS_NOMINAL = 6.68540e+07; // [m]

    // IAU-recommended value for the ephemeris position of the prime meridian of Jupiter. The prime meridian refers to the rotation of the magnetic field System III. System I (W_{I} = 67.1 deg + 877.900 deg day^-1) refers to the mean atmospheric equatorial rotation. System II (W_{II} = 43.3 deg + 870.270 deg day^-1) refers to the mean atmospheric rotation north of the south component of the north equatorial belt, and south of the north component of the south equatorial belt. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double JUPITER_PRIMEMERIDIAN_EPHEMERISPOSITION = 284.95; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch of 1.5 January 2000 = JD2451545.0 TDB) of the ephemeris position of the prime meridian of Jupiter. The prime meridian refers to the rotation of the magnetic field System III. System I (W_{I} = 67.1 deg + 877.900 deg day^-1) refers to the mean atmospheric equatorial rotation. System II (W_{II} = 43.3 deg + 870.270 deg day^-1) refers to the mean atmospheric rotation north of the south component of the north equatorial belt, and south of the north component of the south equatorial belt. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double JUPITER_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE = 870.5360000; // [deg day^-1]

    // Geometric transit probability (Jupiter transiting the Sun) for an observer external to the solar system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.1, Equation 2.2, page 14
    // Basic : false
    // Scalar: true
    static const double JUPITER_TRANSITPROBABILITY = 0.098; // [%]

    // Maximum transit time of Jupiter (transiting the Sun) for an observer external to the solar system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.3, Equation 2.5, page 15
    // Basic : false
    // Scalar: true
    static const double JUPITER_TRANSITTIME_MAXIMUM = 1.36; // [day]

    // V(1,0) magnitude of Jupiter (i.e., the visual magnitude of the planet reduced to a distance of 1 au from both the Sun and Earth and phase angle zero). This parameter is also refered to as absolute magnitude in planetary sciences
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double JUPITER_VONEZEROMAGNITUDE = -9.40; // [mag]

    // Mean volumetric radius of Jupiter
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double JUPITER_VOLUMETRICRADIUS = 6.99110e+07; // [m]

    // Photon flux density N_{\lambda}(\lambda) of an unreddened K1III metal-poor (MP) star (Pickles' star number 082) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)
    // Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863
    // Basic : true
    // Scalar: false
    static const char  *const K1IIIMPSTAR_SPECTRUM_NUMBEROFPHOTONS() { return "Nature/K1IIIMPStar_Spectrum_NumberOfPhotons_001.fits"; }

    // High-resolution photon-flux density N_{\lambda}(\lambda) of an unreddened K1III metal-poor (MP) star at V = 15 mag. The data refer to a high-resolution Kurucz-model spectrum with the following properties: effective temperature T_eff = 4500 K, logarithm of surface gravity log g = 2.0, metallicity [Fe/H] = -1.5, alpha-elements [\alpha/Fe] = +0.4, rotational velocity v sini = 5 km s^-1, micro-turbulence velocity = 2.0 km s^-1, length of convective bubble divided by pressure scale height = 0.50, no convective overshooting, macro-turbulence velocity = 2.0 km s^-1, and resolving power R = \lambda / \delta \lambda = 250,000. First column: wavelength \lambda (in nm; from 830.1673264 to 889.8217922). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1). The 34698 lines have an average wavelength step of 0.00172 nm; the spectrum extent is thus 59.7 nm
    // Source: ESA, 20 June 2005, 'Photon-flux distributions for reference stars', GAIA-EST-TN-00539, issue 1, revision 0, based on D. Katz, priv. comm., 11 May 2005
    // Basic : true
    // Scalar: false
    static const char  *const K1IIIMPSTAR_SPECTRUM_NUMBEROFPHOTONSHIGHRESOLUTION() { return "Nature/K1IIIMPStar_Spectrum_NumberOfPhotonsHighResolution_001.fits"; }

    // Photon flux density N_{\lambda}(\lambda) of an unreddened K3III star (Pickles' star number 087) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)
    // Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863
    // Basic : true
    // Scalar: false
    static const char  *const K3IIISTAR_SPECTRUM_NUMBEROFPHOTONS() { return "Nature/K3IIIStar_Spectrum_NumberOfPhotons_001.fits"; }

    // Photon flux density N_{\lambda}(\lambda) of an unreddened K3V star (Pickles' star number 034) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)
    // Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863
    // Basic : true
    // Scalar: false
    static const char  *const K3VSTAR_SPECTRUM_NUMBEROFPHOTONS() { return "Nature/K3VStar_Spectrum_NumberOfPhotons_001.fits"; }

    // Central auxiliary parameter in GAIA-FM-011, issue 1, revision 0
    // Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Equation 22
    // Basic : false
    // Scalar: true
    static const double L2_ALPHA = 0.0100447147;

    // Axis ratio of the 'horizontal' elliptic motion around L2 (Equations 70-71 in GAIA-FM-011, issue 1, revision 0)
    // Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Equation 66; see also the addition to GAIA-FM-011, issue 1, revision 0, dated 25 September 2003
    // Basic : false
    // Scalar: true
    static const double L2_BETA = 3.18722929;

    // Auxiliary variable in GAIA-FM-011, issue 1, revision 0
    // Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Equations 34 and 38; see also the addition to GAIA-FM-011, issue 1, revision 0, dated 25 September 2003
    // Basic : false
    // Scalar: true
    static const double L2_CAPITALOMEGA = 3.940522185;

    // Reduced mass \mu of the Sun and Earth-Moon system as used in GAIA-FM-011, issue 1, revision 0. Note that the INPOP10e value 328900.552289 gives \mu = 3.04042347E-6
    // Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Section 2
    // Basic : false
    // Scalar: true
    static const double L2_MU = 3.040423402e-06;

    // Frequency of the vertical oscillation around L2
    // Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Equation 40; see also the addition to GAIA-FM-011, issue 1, revision 0, dated 25 September 2003
    // Basic : false
    // Scalar: true
    static const double L2_OMEGA = 1.985074856; // [{sidereal year}^-1]

    // Synodic period of the vertical oscillation around L2 in units of days
    // Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Equation 40
    // Basic : false
    // Scalar: true
    static const double L2_OMEGA_PERIOD = 184.00; // [day]

    // Mean orbital eccentricity of the L2 orbit of the Sun and Earth-Moon system, at the standard epoch J2000.0. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0. DE405 is based upon the International Celestial Reference Frame (ICRF)
    // Source: See, e.g., F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Section 3.2, pages 6-7
    // Basic : false
    // Scalar: true
    static const double L2_ORBITALECCENTRICITY_J2000 = 0.01671123;

    // Mean orbital semi-major axis of the L2 orbit of the Sun and Earth-Moon system, at the standard epoch J2000.0. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0. DE405 is based upon the International Celestial Reference Frame (ICRF)
    // Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Section 3
    // Basic : false
    // Scalar: true
    static const double L2_ORBITALSEMIMAJORAXIS_J2000 = 1.01008088; // [au]

    // The quantity p_n in GAIA-FM-011, issue 1, revision 0; this choice guarantees an eclipse-free orbit around L2 for more than 6 years
    // Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Tables 1 and 2 and Sections 6.2 and 6.4
    // Basic : true
    // Scalar: true
    static const double L2_P = 27.;

    // Maximum radius of the penumbra of the Earth during a solar eclipse as seen from a point 1.0E5 km 'behind' L2. GAIA-FM-011, issue 1, revision 0 uses a rounded value of 14000 km
    // Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Section 4 (symbol \sigma_2) and Equation 86 (symbol s)
    // Basic : false
    // Scalar: true
    static const double L2_PENUMBRARADIUS = 13923.; // [km]

    // The quantity q_n in GAIA-FM-011, issue 1, revision 0; this choice guarantees an eclipse-free orbit around L2 for more than 6 years
    // Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Tables 1 and 2 and Sections 6.2 and 6.4
    // Basic : true
    // Scalar: true
    static const double L2_Q = 26.;

    // The 'residual quantity' p - q a = q \epsilon, with a = \sigma / \omega, in GAIA-FM-011, issue 1, revision 0
    // Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Table 1 and Section 6.4
    // Basic : false
    // Scalar: true
    static const double L2_QTIMESEPSILON = 5.78e-02;

    // Normalised separation between L2 and the Earth-Moon barycentre
    // Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Equation 23; see also the addition to GAIA-FM-011, issue 1, revision 0, dated 25 September 2003
    // Basic : false
    // Scalar: true
    static const double L2_RHO = 0.01007824044;

    // Frequency of the oscillation around L2 in the horizontal plane
    // Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Equations 49 and 52; see also the addition to GAIA-FM-011, issue 1, revision 0, dated 25 September 2003
    // Basic : false
    // Scalar: true
    static const double L2_SIGMA = 2.057014191; // [{sidereal year}^-1]

    // Synodic period of the 'horizontal' elliptic motion around L2 (Equations 70-71 in GAIA-FM-011, issue 1, revision 0) in units of days
    // Source: F. Mignard, 15 March 2002, 'Considerations on the orbit of Gaia for simulations', GAIA-FM-011, issue 1, revision 0, Section 5.4
    // Basic : false
    // Scalar: true
    static const double L2_SIGMA_PERIOD = 177.57; // [day]

    // Maximum (change in) aberration brought about by the acceleration of the LSR relative to the Galactic centre (resulting, if not corrected for, in spurious (apparent) proper motions for extra-Galactic sources in some regions of the sky)
    // Source: J. Kovalevsky, 2003, 'Aberration in proper motions', A&A, 404, 743. See also ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 1.8.10, page 110
    // Basic : false
    // Scalar: true
    static const double LSR_ABERRATION_MAXIMUM = 4.2; // [10^-6 arcsec yr^-1]

    // Angular velocity of circular rotation at the solar Galactocentric radius. Note: Transactions of the IAU, Volume XIXB, 1985, page 254: \Omega_0 = 220/8.5 = 25.88 km s^-1 kpc^-1
    // Basic : false
    // Scalar: true
    static const double LSR_ANGULARVELOCITY = 27.19; // [km s^-1 kpc^-1]

    // Velocity of circular rotation at the solar Galactocentric radius (local circular speed). Note: Transactions of the IAU, Volume XIXB, 1985, page 254: V_0 = 220 km s^-1
    // Basic : false
    // Scalar: true
    static const double LSR_CIRCULARVELOCITY = 217.520; // [km s^-1]

    // Period of rotation around the Galactic centre at the solar Galactocentric radius for a circular orbit
    // Basic : false
    // Scalar: true
    static const double LSR_GALACTICROTATIONPERIOD = 225.95; // [Myr]

    // Distance from the Local Standard of Rest (LSR, or: Sun) to the Galactic centre (solar Galactocentric radius)
    // Source: Z.M. Malkin, 28 February 2012, 'The current best estimate of the Galactocentric distance of the Sun based on comparison of different statistical techniques', http://adsabs.harvard.edu/abs/2012arXiv1202.6128M; see also Z.M. Malkin, 2013, Astronomicheskii Zhurnal, Volume 90, Number 2, pages 152-157 and Z.M. Malkin, 1 February 2013, 'Analysis of determinations of the distance between the sun and the galactic center', Astronomy Reports, Volume 57, Issue 2, pages 128-133
    // Basic : true
    // Scalar: true
    static const double LSR_GALACTOCENTRICRADIUS = 8.0; // [kpc]

    // Average value of 1-d(TT)/d(TCB), defined as a defining constant, based on the '2006 best estimate' of LSubC_Constant + LSubG_Constant - LSubC_Constant * LSubG_Constant
    // Source: IAU, August 2006, 'Re-definition of Barycentric Dynamical Time, TDB', IAU 2006 Resolution 3 adopted at the XXVI-th General Assembly of the IAU. See also the IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double LSUBB_CONSTANT = 1.5505197680e-08;

    // Average value of 1-d(TCG)/d(TCB)
    // Source: A. Irwin, T. Fukushima, 1999, A&A, 348, 642. See also the IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double LSUBC_CONSTANT = 1.480826867410e-08;

    // Value of 1-d(TT)/d(TCG) (defining constant, based on best-estimate EarthEllipsoid_WSub0 = 6.26368556E7 m^2 s^-2)
    // Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html). See also IAU, August 2000, 'Re-definition of Terrestrial Time, TT', IAU 2000 Resolution B1.9 adopted at the XXIV-th General Assembly of the IAU
    // Basic : true
    // Scalar: true
    static const double LSUBG_CONSTANT = 6.9692901340e-10;

    // Light year expressed in au
    // Basic : false
    // Scalar: true
    static const double LIGHTYEAR_ASTRONOMICALUNIT = 63241.077084; // [au]

    // Light year expressed in meters
    // Basic : false
    // Scalar: true
    static const double LIGHTYEAR_METER = 9.4607304726e+15; // [m]

    // Light year expressed in parsecs
    // Basic : false
    // Scalar: true
    static const double LIGHTYEAR_PARSEC = 0.306601; // [pc]

    // Photon flux density N_{\lambda}(\lambda) of an unreddened M0III star (Pickles' star number 095) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)
    // Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863
    // Basic : true
    // Scalar: false
    static const char  *const M0IIISTAR_SPECTRUM_NUMBEROFPHOTONS() { return "Nature/M0IIIStar_Spectrum_NumberOfPhotons_001.fits"; }

    // Photon flux density N_{\lambda}(\lambda) of an unreddened M0V star (Pickles' star number 038) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)
    // Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863
    // Basic : true
    // Scalar: false
    static const char  *const M0VSTAR_SPECTRUM_NUMBEROFPHOTONS() { return "Nature/M0VStar_Spectrum_NumberOfPhotons_001.fits"; }

    // Photon flux density N_{\lambda}(\lambda) of an unreddened M6V star (Pickles' star number 045) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)
    // Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863
    // Basic : true
    // Scalar: false
    static const char  *const M6VSTAR_SPECTRUM_NUMBEROFPHOTONS() { return "Nature/M6VStar_Spectrum_NumberOfPhotons_001.fits"; }

    // Photon flux density N_{\lambda}(\lambda) of an unreddened M7III star (Pickles' star number 102) at V = 0 mag. First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: photon flux density above the Earth's atmosphere (in photons m^-2 s^-1 nm^-1; the photometric zero point is described in J.H.J. de Bruijne, 26 May 2003, 'Stellar fluxes: transformations and calibrations', GAIA-JdB-005, issue 1, revision 1)
    // Source: A.J. Pickles, 1998, 'A stellar spectral flux library: 1150-25000 AA', PASP, 110, 863
    // Basic : true
    // Scalar: false
    static const char  *const M7IIISTAR_SPECTRUM_NUMBEROFPHOTONS() { return "Nature/M7IIIStar_Spectrum_NumberOfPhotons_001.fits"; }

    // Magnetic constant (defining constant)
    // Basic : false
    // Scalar: true
    static const double MAGNETIC_CONSTANT = 1.256637061435917e-06; // [N A^-2]

    // Astrometric signature of the Sun induced by the Mars system for an observer located at a distance of 10 pc from the Sun
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.7, Equation 1.22, page 11
    // Basic : false
    // Scalar: true
    static const double MARSSYSTEM_ASTROMETRICSIGNATURE_10PARSEC = 0.049; // [10^-6 arcsec]

    // Mars-system mass (IAU 2009 CBE value). The planetary mass includes the contribution of its satellites
    // Basic : false
    // Scalar: true
    static const double MARSSYSTEM_MASS = 6.41712e+23; // [kg]

    // Mean orbital eccentricity of Mars, at the standard epoch J2000.0. The mean orbital eccentricity is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double MARSSYSTEM_ORBITALECCENTRICITY_J2000 = 0.09339410;

    // Mean orbital inclination of Mars, at the standard epoch J2000.0. The mean orbital inclination is is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double MARSSYSTEM_ORBITALINCLINATION_J2000 = 1.84969142; // [deg]

    // Sidereal orbital period
    // Source: Values derived from the mean longitude rates in Table 5.8.1 in P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, page 316; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double MARSSYSTEM_ORBITALPERIOD = 1.8808476; // [yr]

    // Mean orbital semi-major axis of Mars, at the standard epoch J2000.0. The mean orbital semi-major axis is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double MARSSYSTEM_ORBITALSEMIMAJORAXIS_J2000 = 1.52371034; // [au]

    // Radial-velocity amplitude of the Sun induced by the Mars system for 'an edge-on observer' (i.e., an observer in the orbital plane of the Mars system)
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.4, Equation 1.18, page 9
    // Basic : false
    // Scalar: true
    static const double MARSSYSTEM_RADIALVELOCITYSIGNATURE = 0.008; // [m s^-1]

    // Radius of the smallest hypothetical sphere around Mars which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double MARS_ENCOMPASSINGSPHERERADIUS = 3.396190e+06; // [m]

    // Equatorial radius of Mars
    // Basic : false
    // Scalar: true
    static const double MARS_EQUATORIALRADIUS = 3.396190e+06; // [m]

    // Geometrical flattening factor f of Mars (f = (a-b)/a). Mars has a significant dichotomy in shape between the northern and southern hemispheres
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double MARS_FLATTENING = 5.89790e-03;

    // Maximum reduction of the solar flux for an observer external to the solar system during a transit of Mars
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.2, Equation 2.4, page 14
    // Basic : false
    // Scalar: true
    static const double MARS_FLUXREDUCTION_MAXIMUM = 0.002; // [%]

    // Geometric albedo of Mars (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double MARS_GEOMETRICALBEDO = 0.150;

    // Dynamical form-factor of Mars (oblateness or Stokes' second-degree zonal harmonic of the potential)
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8
    // Basic : true
    // Scalar: true
    static const double MARS_JSUB2 = 1.9640e-03;

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of Mars
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double MARS_LIGHTDEFLECTION_LIMB = 116.; // [10^-6 arcsec]

    // Mass of Mars (do not use for high-precision (orbit) calculations)
    // Source: R.A. Jacobson, 2008, 'Ephemerides of the Martian Satellites - MAR080', JPL IOM 343.R-08-006; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double MARS_MASS = 6.416930e+23; // [kg]

    // Mean mass density of Mars
    // Basic : false
    // Scalar: true
    static const double MARS_MASSDENSITY_MEAN = 3.9340; // [g cm^-3]

    // IAU-recommended value for the declination \delta_0 of the north pole of rotation of Mars. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double MARS_NORTHROTATIONALPOLE_DECLINATION = 52.88650; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the declination \delta_0 of the north pole of rotation of Mars. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : false
    // Scalar: true
    static const double MARS_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE = -0.0000016674; // [deg day^-1]

    // IAU-recommended value for the right ascension \alpha_0 of the north pole of rotation of Mars. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double MARS_NORTHROTATIONALPOLE_RIGHTASCENSION = 317.68143; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the right ascension \alpha_0 of the north pole of rotation of Mars. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : false
    // Scalar: true
    static const double MARS_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE = -0.0000029049; // [deg day^-1]

    // Polar radius of Mars. Mars has a significant dichotomy in shape between the northern and southern hemispheres: the average polar radius is listed as 3.37620E6 m in B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : false
    // Scalar: true
    static const double MARS_POLARRADIUS = 3.376160e+06; // [m]

    // IAU-recommended value for the ephemeris position of the prime meridian of Mars. The 0-deg meridian is defined by the crater Airy-0. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double MARS_PRIMEMERIDIAN_EPHEMERISPOSITION = 176.630; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch of 1.5 January 2000 = JD2451545.0 TDB) of the ephemeris position of the prime meridian of Mars. The 0-deg meridian is defined by the crater Airy-0. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double MARS_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE = 350.89198226; // [deg day^-1]

    // Geometric transit probability (Mars transiting the Sun) for an observer external to the solar system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.1, Equation 2.2, page 14
    // Basic : false
    // Scalar: true
    static const double MARS_TRANSITPROBABILITY = 0.307; // [%]

    // Maximum transit time of Mars (transiting the Sun) for an observer external to the solar system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.3, Equation 2.5, page 15
    // Basic : false
    // Scalar: true
    static const double MARS_TRANSITTIME_MAXIMUM = 0.67; // [day]

    // V(1,0) magnitude of Mars (i.e., the visual magnitude of the planet reduced to a distance of 1 au from both the Sun and Earth and phase angle zero). This parameter is also refered to as absolute magnitude in planetary sciences
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double MARS_VONEZEROMAGNITUDE = -1.52; // [mag]

    // Mean volumetric radius of Mars
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double MARS_VOLUMETRICRADIUS = 3.389500e+06; // [m]

    // Astrometric signature of the Sun induced by Mercury for an observer located at a distance of 10 pc from the Sun
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.7, Equation 1.22, page 11
    // Basic : false
    // Scalar: true
    static const double MERCURYSYSTEM_ASTROMETRICSIGNATURE_10PARSEC = 0.006; // [10^-6 arcsec]

    // Mercury(-system) mass (IAU 2009 CBE value)
    // Basic : false
    // Scalar: true
    static const double MERCURYSYSTEM_MASS = 3.3011e+23; // [kg]

    // Mean orbital eccentricity of Mercury, at the standard epoch J2000.0. The mean orbital eccentricity is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double MERCURYSYSTEM_ORBITALECCENTRICITY_J2000 = 0.20563593;

    // Mean orbital inclination of Mercury, at the standard epoch J2000.0. The mean orbital inclination is is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double MERCURYSYSTEM_ORBITALINCLINATION_J2000 = 7.00497902; // [deg]

    // Sidereal orbital period
    // Source: Values derived from the mean longitude rates in Table 5.8.1 in P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, page 316; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double MERCURYSYSTEM_ORBITALPERIOD = 0.2408467; // [yr]

    // Mean orbital semi-major axis of Mercury, at the standard epoch J2000.0. The mean orbital semi-major axis is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double MERCURYSYSTEM_ORBITALSEMIMAJORAXIS_J2000 = 0.38709927; // [au]

    // Radial-velocity amplitude of the Sun induced by Mercury for 'an edge-on observer' (i.e., an observer in the orbital plane of Mercury)
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.4, Equation 1.18, page 9
    // Basic : false
    // Scalar: true
    static const double MERCURYSYSTEM_RADIALVELOCITYSIGNATURE = 0.008; // [m s^-1]

    // Radius of the smallest hypothetical sphere around Mercury which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double MERCURY_ENCOMPASSINGSPHERERADIUS = 2.4397e+06; // [m]

    // Equatorial radius of Mercury
    // Basic : false
    // Scalar: true
    static const double MERCURY_EQUATORIALRADIUS = 2.43970e+06; // [m]

    // Geometrical flattening factor f of Mercury (f = (a-b)/a)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double MERCURY_FLATTENING = 0.;

    // Maximum reduction of the solar flux for an observer external to the solar system during a transit of Mercury
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.2, Equation 2.4, page 14
    // Basic : false
    // Scalar: true
    static const double MERCURY_FLUXREDUCTION_MAXIMUM = 0.001; // [%]

    // Geometric albedo of Mercury (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double MERCURY_GEOMETRICALBEDO = 0.106;

    // Dynamical form-factor of Mercury (oblateness or Stokes' second-degree zonal harmonic of the potential)
    // Source: J.D. Anderson, G. Colombo, P.B. Esposito, E.L. Lau, G.B. Trager, 1 September 1987, 'The mass, gravity field, and ephemeris of Mercury', Icarus, 71, 337-349
    // Basic : true
    // Scalar: true
    static const double MERCURY_JSUB2 = 6.00e-05;

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of Mercury
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double MERCURY_LIGHTDEFLECTION_LIMB = 83.; // [10^-6 arcsec]

    // Mass of Mercury (do not use for high-precision (orbit) calculations)
    // Source: J.D. Anderson, et al., 1987, 'The mass, gravity field, and ephemeris of Mercury', Icarus, 71, 337-349; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double MERCURY_MASS = 3.301040e+23; // [kg]

    // Mean mass density of Mercury
    // Basic : false
    // Scalar: true
    static const double MERCURY_MASSDENSITY_MEAN = 5.427; // [g cm^-3]

    // IAU-recommended value for the declination \delta_0 of the north pole of rotation of Mercury. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double MERCURY_NORTHROTATIONALPOLE_DECLINATION = 61.4143; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the declination \delta_0 of the north pole of rotation of Mercury. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : false
    // Scalar: true
    static const double MERCURY_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE = -0.0000001342; // [deg day^-1]

    // IAU-recommended value for the right ascension \alpha_0 of the north pole of rotation of Mercury. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double MERCURY_NORTHROTATIONALPOLE_RIGHTASCENSION = 281.0097; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the right ascension \alpha_0 of the north pole of rotation of Mercury. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : false
    // Scalar: true
    static const double MERCURY_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE = -0.0000008980; // [deg day^-1]

    // Polar radius of Mercury
    // Basic : false
    // Scalar: true
    static const double MERCURY_POLARRADIUS = 2.43970e+06; // [m]

    // IAU-recommended value for the ephemeris position of the prime meridian of Mercury. The 20-deg meridian is defined by the crater Hun Kal. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double MERCURY_PRIMEMERIDIAN_EPHEMERISPOSITION = 329.5469; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch of 1.5 January 2000 = JD2451545.0 TDB) of the ephemeris position of the prime meridian of Mercury. The 20-deg meridian is defined by the crater Hun Kal. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double MERCURY_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE = 6.1385025; // [deg day^-1]

    // Geometric transit probability (Mercury transiting the Sun) for an observer external to the solar system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.1, Equation 2.2, page 14
    // Basic : false
    // Scalar: true
    static const double MERCURY_TRANSITPROBABILITY = 1.206; // [%]

    // Maximum transit time of Mercury (transiting the Sun) for an observer external to the solar system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.3, Equation 2.5, page 15
    // Basic : false
    // Scalar: true
    static const double MERCURY_TRANSITTIME_MAXIMUM = 0.34; // [day]

    // V(1,0) magnitude of Mercury (i.e., the visual magnitude of the planet reduced to a distance of 1 au from both the Sun and Earth and phase angle zero). This parameter is also refered to as absolute magnitude in planetary sciences
    // Source: J.L. Hilton, 2005, 'Improving the Visual Magnitudes of the Planets in The Astronomical Almanac. I. Mercury and Venus', AJ, 129, 2902-2906; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double MERCURY_VONEZEROMAGNITUDE = -0.60; // [mag]

    // Mean volumetric radius of Mercury
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double MERCURY_VOLUMETRICRADIUS = 2.43970e+06; // [m]

    // One micro-arcsecond in units of radians
    // Basic : false
    // Scalar: true
    static const double MICROARCSECOND_RADIAN = 4.848136811095360e-12; // [rad]

    // The typical expected micro-meteoroid flux at L2, in units of particles cm^-2 s^-1 hemisphere^-1, resulting from particles with masses greater than m (in units of kg and for m > 2.8E-11 kg) equals 2.6E-18 * m^(-7/6)
    // Source: L. Lindegren, 13 July 2000, 'Effects of micro-meteoroids on Gaia attitude', GAIA-LL-031, issue 1, revision 0. See also K. Yamakoshi, 1994, 'Extraterrestrial dust; Laboratory studies of interplanetary dust', Astrophysics and Space Science Library, 181, Kluwer Academic Publishers, Dordrecht (1994edls.book.....Y)
    // Basic : true
    // Scalar: true
    static const double MICROMETEOROID_FLUX_L2LARGEPARTICLES = 2.60e-18; // [particles m^-2 s^-1 hemisphere^-1]

    // The typical expected micro-meteoroid flux at L2, in units of particles cm^-2 s^-1 hemisphere^-1, resulting from particles with masses greater than m (in units of kg and for m < 2.8E-11 kg) equals 2.8E-11 * m^(-1/2)
    // Source: L. Lindegren, 13 July 2000, 'Effects of micro-meteoroids on Gaia attitude', GAIA-LL-031, issue 1, revision 0. See also K. Yamakoshi, 1994, 'Extraterrestrial dust; Laboratory studies of interplanetary dust', Astrophysics and Space Science Library, 181, Kluwer Academic Publishers, Dordrecht (1994edls.book.....Y)
    // Basic : true
    // Scalar: true
    static const double MICROMETEOROID_FLUX_L2SMALLPARTICLES = 2.80e-11; // [particles m^-2 s^-1 hemisphere^-1]

    // One milli-arcsecond in units of radians
    // Basic : false
    // Scalar: true
    static const double MILLIARCSECOND_RADIAN = 4.848136811095360e-09; // [rad]

    // Molar gas constant
    // Source: P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0)
    // Basic : true
    // Scalar: true
    static const double MOLARGAS_CONSTANT = 8.3144598; // [J mol^-1 K^-1]

    // Radius of the smallest hypothetical sphere around J1 (Io) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double MOONJ1_ENCOMPASSINGSPHERERADIUS = 1.82149e+06; // [m]

    // GM of J1 (Io)
    // Source: R.A. Jacobson, 2003, 'Constants used in the JUP230 ephemeris', priv. comm.; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONJ1_GM = 5.9599160e+12; // [m^3 s^-2]

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of J1 (Io)
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double MOONJ1_LIGHTDEFLECTION_LIMB = 30.; // [10^-6 arcsec]

    // Mean mass density of J1 (Io)
    // Basic : false
    // Scalar: true
    static const double MOONJ1_MASSDENSITY_MEAN = 3.528; // [g cm^-3]

    // Radius of J1 (Io)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONJ1_RADIUS = 1.821490e+06; // [m]

    // Radius of the smallest hypothetical sphere around J2 (Europa) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double MOONJ2_ENCOMPASSINGSPHERERADIUS = 1.56080e+06; // [m]

    // GM of J2 (Europa)
    // Source: R.A. Jacobson, 2003, 'Constants used in the JUP230 ephemeris', priv. comm.; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONJ2_GM = 3.2027390e+12; // [m^3 s^-2]

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of J2 (Europa)
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double MOONJ2_LIGHTDEFLECTION_LIMB = 19.; // [10^-6 arcsec]

    // Mean mass density of J2 (Europa)
    // Basic : false
    // Scalar: true
    static const double MOONJ2_MASSDENSITY_MEAN = 3.013; // [g cm^-3]

    // Radius of J2 (Europa)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONJ2_RADIUS = 1.56080e+06; // [m]

    // Radius of the smallest hypothetical sphere around J3 (Ganymede) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double MOONJ3_ENCOMPASSINGSPHERERADIUS = 2.63120e+06; // [m]

    // GM of J3 (Ganymede)
    // Source: R.A. Jacobson, 2003, 'Constants used in the JUP230 ephemeris', priv. comm.; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONJ3_GM = 9.8878340e+12; // [m^3 s^-2]

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of J3 (Ganymede)
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double MOONJ3_LIGHTDEFLECTION_LIMB = 34.; // [10^-6 arcsec]

    // Mean mass density of J3 (Ganymede)
    // Basic : false
    // Scalar: true
    static const double MOONJ3_MASSDENSITY_MEAN = 1.942; // [g cm^-3]

    // Radius of J3 (Ganymede)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONJ3_RADIUS = 2.63120e+06; // [m]

    // Radius of the smallest hypothetical sphere around J4 (Callisto) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double MOONJ4_ENCOMPASSINGSPHERERADIUS = 2.41030e+06; // [m]

    // GM of J4 (Callisto)
    // Source: R.A. Jacobson, 2003, 'Constants used in the JUP230 ephemeris', priv. comm.; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONJ4_GM = 7.1792890e+12; // [m^3 s^-2]

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of J4 (Callisto)
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double MOONJ4_LIGHTDEFLECTION_LIMB = 27.; // [10^-6 arcsec]

    // Mean mass density of J4 (Callisto)
    // Basic : false
    // Scalar: true
    static const double MOONJ4_MASSDENSITY_MEAN = 1.834; // [g cm^-3]

    // Radius of J4 (Callisto)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONJ4_RADIUS = 2.41030e+06; // [m]

    // Radius of the smallest hypothetical sphere around N1 (Triton) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double MOONN1_ENCOMPASSINGSPHERERADIUS = 1.35260e+06; // [m]

    // GM of N1 (Triton)
    // Source: R.A. Jacobson, 2009, 'The Orbits of the Neptunian Satellites and the Orientation of the Pole of Neptune', AJ, 137, 4322; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONN1_GM = 1.42760e+12; // [m^3 s^-2]

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of N1 (Triton)
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double MOONN1_LIGHTDEFLECTION_LIMB = 10.; // [10^-6 arcsec]

    // Mean mass density of N1 (Triton)
    // Basic : false
    // Scalar: true
    static const double MOONN1_MASSDENSITY_MEAN = 2.064; // [g cm^-3]

    // Radius of N1 (Triton)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONN1_RADIUS = 1.35260e+06; // [m]

    // Radius of the smallest hypothetical sphere around S3 (Tethys) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double MOONS3_ENCOMPASSINGSPHERERADIUS = 5.3100e+05; // [m]

    // GM of S3 (Tethys)
    // Source: R.A. Jacobson, et al., 2008, 'The Gravity Gield of the Saturnian System and the Orbits of the Major Saturnian satellites', Presented at the Saturn After Cassini-Huygens Symposium held at Imperial College London, UK, 28 July - 1 August 2008; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONS3_GM = 4.120670e+10; // [m^3 s^-2]

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of S3 (Tethys)
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double MOONS3_LIGHTDEFLECTION_LIMB = 1.; // [10^-6 arcsec]

    // Mean mass density of S3 (Tethys)
    // Basic : false
    // Scalar: true
    static const double MOONS3_MASSDENSITY_MEAN = 0.984; // [g cm^-3]

    // Radius of S3 (Tethys)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONS3_RADIUS = 5.3100e+05; // [m]

    // Radius of the smallest hypothetical sphere around S4 (Dione) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double MOONS4_ENCOMPASSINGSPHERERADIUS = 5.614e+05; // [m]

    // GM of S4 (Dione)
    // Source: R.A. Jacobson, et al., 2008, 'The Gravity Gield of the Saturnian System and the Orbits of the Major Saturnian satellites', Presented at the Saturn After Cassini-Huygens Symposium held at Imperial College London, UK, 28 July - 1 August 2008; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONS4_GM = 7.311460e+10; // [m^3 s^-2]

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of S4 (Dione)
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double MOONS4_LIGHTDEFLECTION_LIMB = 1.; // [10^-6 arcsec]

    // Mean mass density of S4 (Dione)
    // Basic : false
    // Scalar: true
    static const double MOONS4_MASSDENSITY_MEAN = 1.478; // [g cm^-3]

    // Radius of S4 (Dione)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONS4_RADIUS = 5.6140e+05; // [m]

    // Radius of the smallest hypothetical sphere around S5 (Rhea) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double MOONS5_ENCOMPASSINGSPHERERADIUS = 7.635e+05; // [m]

    // GM of S5 (Rhea)
    // Source: R.A. Jacobson, et al., 2008, 'The Gravity Gield of the Saturnian System and the Orbits of the Major Saturnian satellites', Presented at the Saturn After Cassini-Huygens Symposium held at Imperial College London, UK, 28 July - 1 August 2008; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONS5_GM = 1.5394260e+11; // [m^3 s^-2]

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of S5 (Rhea)
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double MOONS5_LIGHTDEFLECTION_LIMB = 2.; // [10^-6 arcsec]

    // Mean mass density of S5 (Rhea)
    // Basic : false
    // Scalar: true
    static const double MOONS5_MASSDENSITY_MEAN = 1.237; // [g cm^-3]

    // Radius of S5 (Rhea)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONS5_RADIUS = 7.6350e+05; // [m]

    // Radius of the smallest hypothetical sphere around S6 (Titan) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double MOONS6_ENCOMPASSINGSPHERERADIUS = 2.5747e+06; // [m]

    // GM of S6 (Titan)
    // Source: R.A. Jacobson, et al., 2008, 'The Gravity Gield of the Saturnian System and the Orbits of the Major Saturnian satellites', Presented at the Saturn After Cassini-Huygens Symposium held at Imperial College London, UK, 28 July - 1 August 2008; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONS6_GM = 8.97813820e+12; // [m^3 s^-2]

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of S6 (Titan)
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double MOONS6_LIGHTDEFLECTION_LIMB = 32.; // [10^-6 arcsec]

    // Mean mass density of S6 (Titan)
    // Basic : false
    // Scalar: true
    static const double MOONS6_MASSDENSITY_MEAN = 1.882; // [g cm^-3]

    // Radius of S6 (Titan)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONS6_RADIUS = 2.574730e+06; // [m]

    // Radius of the smallest hypothetical sphere around S8 (Iapetus) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double MOONS8_ENCOMPASSINGSPHERERADIUS = 7.343e+05; // [m]

    // GM of S8 (Iapetus)
    // Source: R.A. Jacobson, et al., 2008, 'The Gravity Gield of the Saturnian System and the Orbits of the Major Saturnian satellites', Presented at the Saturn After Cassini-Huygens Symposium held at Imperial College London, UK, 28 July - 1 August 2008; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONS8_GM = 1.2050380e+11; // [m^3 s^-2]

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of S8 (Iapetus)
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double MOONS8_LIGHTDEFLECTION_LIMB = 2.; // [10^-6 arcsec]

    // Mean mass density of S8 (Iapetus)
    // Basic : false
    // Scalar: true
    static const double MOONS8_MASSDENSITY_MEAN = 1.089; // [g cm^-3]

    // Radius of S8 (Iapetus)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONS8_RADIUS = 7.3430e+05; // [m]

    // Ratio of Moon to Earth mass
    // Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double MOONTOEARTH_MASSRATIO = 0.0123000371;

    // Radius of the smallest hypothetical sphere around U1 (Ariel) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double MOONU1_ENCOMPASSINGSPHERERADIUS = 5.7890e+05; // [m]

    // GM of U1 (Ariel)
    // Source: R.A. Jacobson, 2007, 'The Gravity Field of the Uranian System and the Orbits of the Uranian Satellites and Rings', BAAS, 39, 453; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONU1_GM = 8.640e+10; // [m^3 s^-2]

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of U1 (Ariel)
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double MOONU1_LIGHTDEFLECTION_LIMB = 1.; // [10^-6 arcsec]

    // Mean mass density of U1 (Ariel)
    // Basic : false
    // Scalar: true
    static const double MOONU1_MASSDENSITY_MEAN = 1.593; // [g cm^-3]

    // Radius of U1 (Ariel)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONU1_RADIUS = 5.7890e+05; // [m]

    // Radius of the smallest hypothetical sphere around U2 (Umbriel) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double MOONU2_ENCOMPASSINGSPHERERADIUS = 5.8470e+05; // [m]

    // GM of U2 (Umbriel)
    // Source: R.A. Jacobson, 2007, 'The Gravity Field of the Uranian System and the Orbits of the Uranian Satellites and Rings', BAAS, 39, 453; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONU2_GM = 8.150e+10; // [m^3 s^-2]

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of U2 (Umbriel)
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double MOONU2_LIGHTDEFLECTION_LIMB = 1.; // [10^-6 arcsec]

    // Mean mass density of U2 (Umbriel)
    // Basic : false
    // Scalar: true
    static const double MOONU2_MASSDENSITY_MEAN = 1.458; // [g cm^-3]

    // Radius of U2 (Umbriel)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONU2_RADIUS = 5.8470e+05; // [m]

    // Radius of the smallest hypothetical sphere around U3 (Titania) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double MOONU3_ENCOMPASSINGSPHERERADIUS = 7.8890e+05; // [m]

    // GM of U3 (Titania)
    // Source: R.A. Jacobson, 2007, 'The Gravity Field of the Uranian System and the Orbits of the Uranian Satellites and Rings', BAAS, 39, 453; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONU3_GM = 2.2820e+11; // [m^3 s^-2]

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of U3 (Titania)
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double MOONU3_LIGHTDEFLECTION_LIMB = 3.; // [10^-6 arcsec]

    // Mean mass density of U3 (Titania)
    // Basic : false
    // Scalar: true
    static const double MOONU3_MASSDENSITY_MEAN = 1.663; // [g cm^-3]

    // Radius of U3 (Titania)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONU3_RADIUS = 7.8890e+05; // [m]

    // Radius of the smallest hypothetical sphere around U4 (Oberon) which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double MOONU4_ENCOMPASSINGSPHERERADIUS = 7.6140e+05; // [m]

    // GM of U4 (Oberon)
    // Source: R.A. Jacobson, 2007, 'The Gravity Field of the Uranian System and the Orbits of the Uranian Satellites and Rings', BAAS, 39, 453; see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONU4_GM = 1.9240e+11; // [m^3 s^-2]

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of U4 (Oberon)
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double MOONU4_LIGHTDEFLECTION_LIMB = 2.; // [10^-6 arcsec]

    // Mean mass density of U4 (Oberon)
    // Basic : false
    // Scalar: true
    static const double MOONU4_MASSDENSITY_MEAN = 1.559; // [g cm^-3]

    // Radius of U4 (Oberon)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOONU4_RADIUS = 7.6140e+05; // [m]

    // Lunar diurnal parallax. Formally, this parameter is defined as the ratio of a fictitious mean equatorial radius of the Earth to the perturbed mean distance of the Moon; the ratio F_2 of the perturbed mean distance to the Moon (the perturbation being due to the Sun) to the two-body mean distance of the Moon (with the Sun not present and constant mean motion) equals 0.999093141975298 (see T.D. Moyer, 15 May 1971, 'Mathematical formulation of the Double-Precision Orbit Determination Programme (DPODP)', NASA JPL Technical Report 32-1527, pages 25-26)
    // Basic : false
    // Scalar: true
    static const double MOON_DIURNALPARALLAX = 3422.595; // [arcsec]

    // Radius of the smallest hypothetical sphere around the Moon which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double MOON_ENCOMPASSINGSPHERERADIUS = 1738000.; // [m]

    // Selenocentric gravitational constant (TCB-compatible value)
    // Basic : false
    // Scalar: true
    static const double MOON_GM = 4.9028002e+12; // [m^3 s^-2]

    // Geometric albedo of the Moon (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8
    // Basic : true
    // Scalar: true
    static const double MOON_GEOMETRICALBEDO = 0.12;

    // Lunar mass (do not use for high-precision (orbit) calculations)
    // Basic : false
    // Scalar: true
    static const double MOON_MASS = 7.3460e+22; // [kg]

    // Mean Moon mass density
    // Basic : false
    // Scalar: true
    static const double MOON_MASSDENSITY_MEAN = 3.344; // [g cm^-3]

    // Johnson V band mean opposition magnitude of the Moon
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8
    // Basic : true
    // Scalar: true
    static const double MOON_OPPOSITIONVMAGNITUDE = -12.74; // [mag]

    // Eccentricity of Lunar orbit (mean ecliptic orbital elements, at the standard epoch J2000.0, based on JPL's Planetary and Lunar Ephemerides DE405/LE405)
    // Source: E.M. Standish, 2001, 'Approximate Mean Ecliptic Elements of the Lunar Orbit', JPL IOM 312.F-01-004 (http://ssd.jpl.nasa.gov/?sat_elem)
    // Basic : true
    // Scalar: true
    static const double MOON_ORBITALECCENTRICITY_J2000 = 0.0554;

    // Inclination of Lunar orbit with respect to the ecliptic (mean ecliptic orbital elements, at the standard epoch J2000.0, based on JPL's Planetary and Lunar Ephemerides DE405/LE405)
    // Source: E.M. Standish, 2001, 'Approximate Mean Ecliptic Elements of the Lunar Orbit', JPL IOM 312.F-01-004 (http://ssd.jpl.nasa.gov/?sat_elem)
    // Basic : true
    // Scalar: true
    static const double MOON_ORBITALINCLINATION_J2000 = 5.16; // [deg]

    // Semi-major axis of Lunar orbit (mean ecliptic orbital elements, at the standard epoch J2000.0, based on JPL's Planetary and Lunar Ephemerides DE405/LE405)
    // Source: E.M. Standish, 2001, 'Approximate Mean Ecliptic Elements of the Lunar Orbit', JPL IOM 312.F-01-004 (http://ssd.jpl.nasa.gov/?sat_elem)
    // Basic : true
    // Scalar: true
    static const double MOON_ORBITALSEMIMAJORAXIS_J2000 = 3.844000e+08; // [m]

    // Precession period of the argument of periapsis of Lunar orbit, i.e., apsidal period (mean ecliptic orbital elements, at the standard epoch J2000.0, based on JPL's Planetary and Lunar Ephemerides DE405/LE405)
    // Source: E.M. Standish, 2001, 'Approximate Mean Ecliptic Elements of the Lunar Orbit', JPL IOM 312.F-01-004 (http://ssd.jpl.nasa.gov/?sat_elem)
    // Basic : true
    // Scalar: true
    static const double MOON_PRECESSIONPERIOD_J2000ARGUMENTOFPERIAPSIS = 5.997; // [yr]

    // Precession period of the longitude of the ascending node of Lunar orbit, i.e., nodal period (mean ecliptic orbital elements, at the standard epoch J2000.0, based on JPL's Planetary and Lunar Ephemerides DE405/LE405)
    // Source: E.M. Standish, 2001, 'Approximate Mean Ecliptic Elements of the Lunar Orbit', JPL IOM 312.F-01-004 (http://ssd.jpl.nasa.gov/?sat_elem)
    // Basic : true
    // Scalar: true
    static const double MOON_PRECESSIONPERIOD_J2000LONGITUDEOFASCENDINGNODE = 18.600; // [yr]

    // Sidereal period of Lunar orbit (mean ecliptic orbital elements, at the standard epoch J2000.0, based on JPL's Planetary and Lunar Ephemerides DE405/LE405)
    // Source: E.M. Standish, 2001, 'Approximate Mean Ecliptic Elements of the Lunar Orbit', JPL IOM 312.F-01-004 (http://ssd.jpl.nasa.gov/?sat_elem)
    // Basic : true
    // Scalar: true
    static const double MOON_SIDEREALPERIOD_J2000 = 27.322; // [day]

    // Mean volumetric radius of the Moon
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 3 April 2009, 'Planetary Satellite Physical Parameters', http://ssd.jpl.nasa.gov/?sat_phys_par
    // Basic : true
    // Scalar: true
    static const double MOON_VOLUMETRICRADIUS = 1.73740e+06; // [m]

    // Napier's constant (also known as Neper's constant), i.e., base of the natural logarithm. Although the symbol 'e' refers to Euler, the Napier constant should not be confused with the Euler(-Mascheroni) constant \gamma = 0.5772156649... Note that double-precision, floating-point numbers in any programming language which follows the IEEE standard (true for C, C++, Java, and most, if not all, others) have only 16 significant digits (64 bits); the representation here, using 30 significant digits, is thus amply sufficient
    // Source: Well-known mathematical constant; numerical value can be extracted, e.g., from Mathematica 4.0 for Solaris (Wolfram Research, Inc.) using 'N[Exp[1],30]'
    // Basic : true
    // Scalar: true
    static const double NAPIER_CONSTANT = 2.71828182845904523536028747135;

    // Diameter of near-Earth asteroid 1036 Ganymed
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 16 'Prominent Aten-, Apollo-, and Amor-class near-Earth asteroids with well-determined orbits and V(1,0) < 18 mag' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double NEAREARTHASTEROID1036GANYMED_DIAMETER = 38.5; // [km]

    // Diameter of near-Earth asteroid 1580 Betulia
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 16 'Prominent Aten-, Apollo-, and Amor-class near-Earth asteroids with well-determined orbits and V(1,0) < 18 mag' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double NEAREARTHASTEROID1580BETULIA_DIAMETER = 7.4; // [km]

    // Diameter of near-Earth asteroid 1627 Ivar
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 16 'Prominent Aten-, Apollo-, and Amor-class near-Earth asteroids with well-determined orbits and V(1,0) < 18 mag' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double NEAREARTHASTEROID1627IVAR_DIAMETER = 8.1; // [km]

    // Diameter of near-Earth asteroid 1685 Toro
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 16 'Prominent Aten-, Apollo-, and Amor-class near-Earth asteroids with well-determined orbits and V(1,0) < 18 mag' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double NEAREARTHASTEROID1685TORO_DIAMETER = 5.2; // [km]

    // Diameter of near-Earth asteroid 1866 Sisyphus
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 16 'Prominent Aten-, Apollo-, and Amor-class near-Earth asteroids with well-determined orbits and V(1,0) < 18 mag' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double NEAREARTHASTEROID1866SISYPHUS_DIAMETER = 8.2; // [km]

    // Diameter of near-Earth asteroid 3200 Phaeton
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 16 'Prominent Aten-, Apollo-, and Amor-class near-Earth asteroids with well-determined orbits and V(1,0) < 18 mag' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double NEAREARTHASTEROID3200PHAETON_DIAMETER = 6.9; // [km]

    // Diameter of near-Earth asteroid 3552 Don Quixote
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 16 'Prominent Aten-, Apollo-, and Amor-class near-Earth asteroids with well-determined orbits and V(1,0) < 18 mag' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double NEAREARTHASTEROID3552DONQUIXOTE_DIAMETER = 18.7; // [km]

    // Diameter of near-Earth asteroid 433 Eros
    // Source: C.F. Yoder, 1995, 'Astrometric and Geodetic Properties of Earth and the solar System', in 'Global Earth Physics; A Handbook of Physical Constants', AGU Reference Shelf 1, American Geophysical Union, Table 16 'Prominent Aten-, Apollo-, and Amor-class near-Earth asteroids with well-determined orbits and V(1,0) < 18 mag' (http://www.agu.org/books/rf/v001/RF001p0001/RF001p0001.pdf)
    // Basic : true
    // Scalar: true
    static const double NEAREARTHASTEROID433EROS_DIAMETER = 22.; // [km]

    // The velocity distribution of near-Earth objects (NEOs) is approximately Gaussian with zero mean and a standard deviation of 30.0 mas s^-1 across-scan (for a solar-aspect angle of 45 degrees)
    // Source: F. Mignard, 2002, 'Observations of solar-system objects with Gaia. I. Detection of NEOs', A&A, 393, 727, Section 4.4 (2002A&A...393..727M). See also E. Hoeg, F. Arenou, P. Hjorth, U.G. Joergensen, F. Mignard, S. Wolff, 28 February 2003, 'Faint objects and NEOs with Gaia', GAIA-CUO-118, issue 1, revision 0. Current value, for a solar-aspect angle of 45 degrees, from F. Mignard, priv. comm., 10 August 2005
    // Basic : true
    // Scalar: true
    static const double NEAREARTHOBJECT_VELOCITYDISTRIBUTIONSTANDARDDEVIATION_AC = 30.0; // [mas s^-1]

    // The velocity distribution of near-Earth objects (NEOs) is approximately Gaussian with zero mean and a standard deviation of 22.5 mas s^-1 along-scan (for a solar-aspect angle of 45 degrees)
    // Source: F. Mignard, 2002, 'Observations of solar-system objects with Gaia. I. Detection of NEOs', A&A, 393, 727, Section 4.4 (2002A&A...393..727M). See also E. Hoeg, F. Arenou, P. Hjorth, U.G. Joergensen, F. Mignard, S. Wolff, 28 February 2003, 'Faint objects and NEOs with Gaia', GAIA-CUO-118, issue 1, revision 0. Current value, for a solar-aspect angle of 45 degrees, from F. Mignard, priv. comm., 10 August 2005
    // Basic : true
    // Scalar: true
    static const double NEAREARTHOBJECT_VELOCITYDISTRIBUTIONSTANDARDDEVIATION_AL = 22.5; // [mas s^-1]

    // Astrometric signature of the Sun induced by the Neptune system for an observer located at a distance of 10 pc from the Sun
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.7, Equation 1.22, page 11
    // Basic : false
    // Scalar: true
    static const double NEPTUNESYSTEM_ASTROMETRICSIGNATURE_10PARSEC = 155.; // [10^-6 arcsec]

    // Neptune-system mass (IAU 2009 CBE value). The planetary mass includes the contribution of its satellites
    // Basic : false
    // Scalar: true
    static const double NEPTUNESYSTEM_MASS = 1.02434e+26; // [kg]

    // Mean orbital eccentricity of Neptune, at the standard epoch J2000.0. The mean orbital eccentricity is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double NEPTUNESYSTEM_ORBITALECCENTRICITY_J2000 = 0.00859048;

    // Mean orbital inclination of Neptune, at the standard epoch J2000.0. The mean orbital inclination is is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double NEPTUNESYSTEM_ORBITALINCLINATION_J2000 = 1.77004347; // [deg]

    // Sidereal orbital period
    // Source: Values derived from the mean longitude rates in Table 5.8.1 in P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, page 316; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double NEPTUNESYSTEM_ORBITALPERIOD = 164.79132; // [yr]

    // Mean orbital semi-major axis of Neptune, at the standard epoch J2000.0. The mean orbital semi-major axis is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double NEPTUNESYSTEM_ORBITALSEMIMAJORAXIS_J2000 = 30.06992276; // [au]

    // Radial-velocity amplitude of the Sun induced by the Neptune system for 'an edge-on observer' (i.e., an observer in the orbital plane of the Neptune system)
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.4, Equation 1.18, page 9
    // Basic : false
    // Scalar: true
    static const double NEPTUNESYSTEM_RADIALVELOCITYSIGNATURE = 0.3; // [m s^-1]

    // Radius of the smallest hypothetical sphere around Neptune which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double NEPTUNE_ENCOMPASSINGSPHERERADIUS = 2.47640e+07; // [m]

    // Equatorial radius of Neptune
    // Basic : false
    // Scalar: true
    static const double NEPTUNE_EQUATORIALRADIUS = 2.47640e+07; // [m]

    // Geometrical flattening factor f of Neptune (f = (a-b)/a)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double NEPTUNE_FLATTENING = 1.710e-02;

    // Maximum reduction of the solar flux for an observer external to the solar system during a transit of Neptune
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.2, Equation 2.4, page 14
    // Basic : false
    // Scalar: true
    static const double NEPTUNE_FLUXREDUCTION_MAXIMUM = 0.125; // [%]

    // Geometric albedo of Neptune (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double NEPTUNE_GEOMETRICALBEDO = 0.41;

    // Dynamical form-factor of Neptune (oblateness or Stokes' second-degree zonal harmonic of the potential)
    // Source: P.R. Weissman, L.-A. McFadden, T.V. Johnson (eds.), 1999, 'Encyclopedia of the Solar System (first edition)', Academic Press, page 342
    // Basic : true
    // Scalar: true
    static const double NEPTUNE_JSUB2 = 0.003538;

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of Neptune
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double NEPTUNE_LIGHTDEFLECTION_LIMB = 2548.; // [10^-6 arcsec]

    // Mass of Neptune (do not use for high-precision (orbit) calculations)
    // Source: R.A. Jacobson, 2008, 'The orbits of the Neptunian satellites and the orientation of the pole of Neptune', BAAS, 40, 296; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double NEPTUNE_MASS = 1.024100e+26; // [kg]

    // Mean mass density of Neptune
    // Basic : false
    // Scalar: true
    static const double NEPTUNE_MASSDENSITY_MEAN = 1.638; // [g cm^-3]

    // IAU-recommended value for the declination \delta_0 of the north pole of rotation of Neptune. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value is based on an approximate formula, accurate to first order in 'd' (see below); the true equation is \delta_0 = 43.46 - 0.51 * cos(357.85 + 52.316 * d / 36525), where d is the number of Julian days calculated from the standard epoch. The numerical accuracy of this equation is 0.1 deg
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : false
    // Scalar: true
    static const double NEPTUNE_NORTHROTATIONALPOLE_DECLINATION = 42.950; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the declination \delta_0 of the north pole of rotation of Neptune. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value is based on an approximate formula, accurate to first order in 'd' (see below); the true equation is \delta_0 = 43.46 - 0.51 * cos(357.85 + 52.316 * d / 36525), where d is the number of Julian days calculated from the standard epoch. The numerical accuracy of this equation is 0.1 deg
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : false
    // Scalar: true
    static const double NEPTUNE_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE = -0.0000004783; // [deg day^-1]

    // IAU-recommended value for the right ascension \alpha_0 of the north pole of rotation of Neptune. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value is based on an approximate formula, accurate to first order in 'd' (see below); the true equation is \alpha_0 = 299.36 + 0.70 * sin(357.85 + 52.316 * d / 36525), where d is the number of Julian days calculated from the standard epoch. The numerical accuracy of this equation is 0.1 deg
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : false
    // Scalar: true
    static const double NEPTUNE_NORTHROTATIONALPOLE_RIGHTASCENSION = 299.334; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the right ascension \alpha_0 of the north pole of rotation of Neptune. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value is based on an approximate formula, accurate to first order in 'd' (see below); the true equation is \alpha_0 = 299.36 + 0.70 * sin(357.85 + 52.316 * d / 36525), where d is the number of Julian days calculated from the standard epoch. The numerical accuracy of this equation is 0.1 deg
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : false
    // Scalar: true
    static const double NEPTUNE_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE = 0.0000174869; // [deg day^-1]

    // Polar radius of Neptune
    // Basic : false
    // Scalar: true
    static const double NEPTUNE_POLARRADIUS = 2.43405e+07; // [m]

    // IAU-recommended value for the ephemeris position of the prime meridian of Neptune. The prime meridian refers to the rotation of the magnetic field System III. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde. The numerical value is based on an approximate formula, accurate to first order in 'd' (see below); the true equation is W = 253.18 + 536.3128492 * d - 0.48 * sin(357.85 + 52.316 * d / 36525), where d is the number of Julian days calculated from the standard epoch. The numerical accuracy of this equation is 0.1 deg
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : false
    // Scalar: true
    static const double NEPTUNE_PRIMEMERIDIAN_EPHEMERISPOSITION = 253.198; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch of 1.5 January 2000 = JD2451545.0 TDB) of the ephemeris position of the prime meridian of Neptune. The prime meridian refers to the rotation of the magnetic field System III. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde. The numerical value is based on an approximate formula, accurate to first order in 'd' (see below); the true equation is W = 253.18 + 536.3128492 * d - 0.48 * sin(357.85 + 52.316 * d / 36525), where d is the number of Julian days calculated from the standard epoch. The numerical accuracy of this equation is 0.1 deg
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : false
    // Scalar: true
    static const double NEPTUNE_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE = 536.3128372090; // [deg day^-1]

    // Geometric transit probability (Neptune transiting the Sun) for an observer external to the solar system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.1, Equation 2.2, page 14
    // Basic : false
    // Scalar: true
    static const double NEPTUNE_TRANSITPROBABILITY = 0.016; // [%]

    // Maximum transit time of Neptune (transiting the Sun) for an observer external to the solar system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.3, Equation 2.5, page 15
    // Basic : false
    // Scalar: true
    static const double NEPTUNE_TRANSITTIME_MAXIMUM = 3.07; // [day]

    // V(1,0) magnitude of Neptune (i.e., the visual magnitude of the planet reduced to a distance of 1 au from both the Sun and Earth and phase angle zero). This parameter is also refered to as absolute magnitude in planetary sciences
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double NEPTUNE_VONEZEROMAGNITUDE = -6.87; // [mag]

    // Mean volumetric radius of Neptune
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double NEPTUNE_VOLUMETRICRADIUS = 2.46220e+07; // [m]

    // Newton's universal constant of gravitation
    // Source: P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0). See also the IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double NEWTON_CONSTANT = 6.674080e-11; // [m^3 kg^-1 s^-2]

    // Mean (geometric) longitude rate of the nominal Sun for use in simulations of the NSL (mean ecliptic orbital elements, at the standard epoch J2000.0). Note that a value of 1295977422.83429 / (1.0E3 * 365.25 * 3600.0) = 0.98560911 degrees day^-1 is given in Section 5.8.3 of J.L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touze, G. Francou, J. Laskar, 1994, 'Numerical expressions for precession formulae and mean elements for the Moon and the planets', A&A, 282, 663 (1994A&A...282..663S)
    // Source: F. Mignard, priv. comm., 14 July 2004, based on the Fortran-90 subroutine SCANNING (version 4.1, April 2004); see also GAIA-FM-010 and GAIA-FM-017
    // Basic : true
    // Scalar: true
    static const double NOMINALSUN_MEANLONGITUDERATE_J2000 = 0.98560903; // [deg day^-1]

    // Mean (geometric) longitude of the nominal Sun for use in simulations of the NSL (mean ecliptic orbital elements, at the standard epoch J2000.0); subtract aberration (about 20 arcsec; see parameter :Nature:Aberration_Constant_J2000) to obtain the apparent longitude of the nominal Sun. Note that a value of 280.46645683 degrees is given in Section 5.8.3 of J.L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touze, G. Francou, J. Laskar, 1994, 'Numerical expressions for precession formulae and mean elements for the Moon and the planets', A&A, 282, 663 (1994A&A...282..663S)
    // Source: F. Mignard, priv. comm., 14 July 2004, based on the Fortran-90 subroutine SCANNING (version 4.1, April 2004); see also GAIA-FM-010 and GAIA-FM-017
    // Basic : true
    // Scalar: true
    static const double NOMINALSUN_MEANLONGITUDE_J2000 = 280.4665; // [deg]

    // Mean orbital eccentricity of the nominal Sun for use in simulations of the NSL (mean ecliptic orbital elements, at the standard epoch J2000.0). See also the parameter :Nature:EMBC_OrbitalEccentricity_J2000. Note that a value of 0.0167086342 is given in Section 5.8.3 of J.L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touze, G. Francou, J. Laskar, 1994, 'Numerical expressions for precession formulae and mean elements for the Moon and the planets', A&A, 282, 663 (1994A&A...282..663S)
    // Source: F. Mignard, priv. comm., 14 July 2004, based on the Fortran-90 subroutine SCANNING (version 4.1, April 2004); see also GAIA-FM-010 and GAIA-FM-017
    // Basic : true
    // Scalar: true
    static const double NOMINALSUN_ORBITALECCENTRICITY_J2000 = 0.01671;

    // Orbital mean anomaly rate of the nominal Sun for use in simulations of the NSL (mean ecliptic orbital elements, at the standard epoch J2000.0). Note that a value of (1295977422.83429 - 11612.35290) / (1.0E3 * 365.25 * 3600.0) = 0.98560028 degrees day^-1 is given in Section 5.8.3 of J.L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touze, G. Francou, J. Laskar, 1994, 'Numerical expressions for precession formulae and mean elements for the Moon and the planets', A&A, 282, 663 (1994A&A...282..663S)
    // Source: F. Mignard, priv. comm., 14 July 2004, based on the Fortran-90 subroutine SCANNING (version 4.1, April 2004); see also GAIA-FM-010 and GAIA-FM-017
    // Basic : true
    // Scalar: true
    static const double NOMINALSUN_ORBITALMEANANOMALYRATE_J2000 = 0.98560020; // [deg day^-1]

    // Orbital mean anomaly of the nominal Sun for use in simulations of the NSL (mean ecliptic orbital elements, at the standard epoch J2000.0). Note that a value of 100.46645683 - 102.93734808 + 360 = 357.52910875 degrees is given in Section 5.8.3 of J.L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touze, G. Francou, J. Laskar, 1994, 'Numerical expressions for precession formulae and mean elements for the Moon and the planets', A&A, 282, 663 (1994A&A...282..663S)
    // Source: F. Mignard, priv. comm., 14 July 2004, based on the Fortran-90 subroutine SCANNING (version 4.1, April 2004); see also GAIA-FM-010 and GAIA-FM-017
    // Basic : true
    // Scalar: true
    static const double NOMINALSUN_ORBITALMEANANOMALY_J2000 = 357.529; // [deg]

    // Constant of nutation, at the standard epoch J2000.0, nowadays irrelevant as a fundamental constant
    // Source: P.K. Seidelmann, May 1982, '1980 IAU theory of nutation. The final report of the IAU Working Group on Nutation', Celestial Mechanics, 27, pages 79-106 (1982CeMec..27...79S)
    // Basic : true
    // Scalar: true
    static const double NUTATION_CONSTANT_J2000 = 9.2025; // [arcsec]

    // Obliquity of the ecliptic with respect to the ICRS reference plane, at the standard epoch J2000.0. Note that the ICRS origin is shifted in the equatorial plane from \Gamma by \phi = 0.05542 arcsec, positive from \Gamma to the ICRS origin. Note that the value of the obliquity of the ecliptic in the inertial sense, i.e., with respect to the CIP equator (Celestial Intermediate Pole, formerly the Mean Celestial Ephemeris Pole or MCEP) equals 84381.406 arcsec (from the IAU (2009) System of Astronomical Constants: IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Source: J. Chapront, M. Chapront-Touze, G. Francou, 2002, 'A new determination of lunar orbital parameters, precession constant, and tidal acceleration from LLR measurements', A&A, 387, 700
    // Basic : true
    // Scalar: true
    static const double OBLIQUITYOFECLIPTIC_J2000 = 84381.41100; // [arcsec]

    // Oort constant A
    // Source: M. Feast, P. Whitelock, 1997, MNRAS, 291, 683
    // Basic : true
    // Scalar: true
    static const double OORT_CONSTANT_A = 14.82; // [km s^-1 kpc^-1]

    // Oort constant B
    // Source: M. Feast, P. Whitelock, 1997, MNRAS, 291, 683
    // Basic : true
    // Scalar: true
    static const double OORT_CONSTANT_B = -12.37; // [km s^-1 kpc^-1]

    // General relativistic standard PPN parameter \alpha_1, quantifying prefered-frame effects
    // Source: E.g., C.M. Will, 2006, 'The Confrontation between General Relativity and Experiment', Living Reviews in Relativity, Section 3.2 'The parametrized post-Newtonian formalism', http://www.livingreviews.org/lrr-2006-3
    // Basic : true
    // Scalar: true
    static const double PPN_ALPHA_1 = 0.;

    // General relativistic standard PPN parameter \alpha_2, quantifying prefered-frame effects
    // Source: E.g., C.M. Will, 2006, 'The Confrontation between General Relativity and Experiment', Living Reviews in Relativity, Section 3.2 'The parametrized post-Newtonian formalism', http://www.livingreviews.org/lrr-2006-3
    // Basic : true
    // Scalar: true
    static const double PPN_ALPHA_2 = 0.;

    // General relativistic standard PPN parameter \alpha_3, quantifying prefered-frame effects
    // Source: E.g., C.M. Will, 2006, 'The Confrontation between General Relativity and Experiment', Living Reviews in Relativity, Section 3.2 'The parametrized post-Newtonian formalism', http://www.livingreviews.org/lrr-2006-3
    // Basic : true
    // Scalar: true
    static const double PPN_ALPHA_3 = 0.;

    // General relativistic standard PPN parameter \beta, quantifying non-linearity in the gravitational superposition law
    // Source: E.g., C.M. Will, 2006, 'The Confrontation between General Relativity and Experiment', Living Reviews in Relativity, Section 3.2 'The parametrized post-Newtonian formalism', http://www.livingreviews.org/lrr-2006-3
    // Basic : true
    // Scalar: true
    static const double PPN_BETA = 1.;

    // General relativistic standard PPN parameter \gamma, quantifying space curvature per unit rest mass
    // Source: E.g., C.M. Will, 2006, 'The Confrontation between General Relativity and Experiment', Living Reviews in Relativity, Section 3.2 'The parametrized post-Newtonian formalism', http://www.livingreviews.org/lrr-2006-3
    // Basic : true
    // Scalar: true
    static const double PPN_GAMMA = 1.;

    // General relativistic standard PPN parameter \xi, quantifying prefered-location effects
    // Source: E.g., C.M. Will, 2006, 'The Confrontation between General Relativity and Experiment', Living Reviews in Relativity, Section 3.2 'The parametrized post-Newtonian formalism', http://www.livingreviews.org/lrr-2006-3
    // Basic : true
    // Scalar: true
    static const double PPN_XI = 0.;

    // General relativistic standard PPN parameter \zeta_1, quantifying violation of conservation of total momentum
    // Source: E.g., C.M. Will, 2006, 'The Confrontation between General Relativity and Experiment', Living Reviews in Relativity, Section 3.2 'The parametrized post-Newtonian formalism', http://www.livingreviews.org/lrr-2006-3
    // Basic : true
    // Scalar: true
    static const double PPN_ZETA_1 = 0.;

    // General relativistic standard PPN parameter \zeta_2, quantifying violation of conservation of total momentum
    // Source: E.g., C.M. Will, 2006, 'The Confrontation between General Relativity and Experiment', Living Reviews in Relativity, Section 3.2 'The parametrized post-Newtonian formalism', http://www.livingreviews.org/lrr-2006-3
    // Basic : true
    // Scalar: true
    static const double PPN_ZETA_2 = 0.;

    // General relativistic standard PPN parameter \zeta_3, quantifying violation of conservation of total momentum
    // Source: E.g., C.M. Will, 2006, 'The Confrontation between General Relativity and Experiment', Living Reviews in Relativity, Section 3.2 'The parametrized post-Newtonian formalism', http://www.livingreviews.org/lrr-2006-3
    // Basic : true
    // Scalar: true
    static const double PPN_ZETA_3 = 0.;

    // General relativistic standard PPN parameter \zeta_4, quantifying violation of conservation of total momentum
    // Source: E.g., C.M. Will, 2006, 'The Confrontation between General Relativity and Experiment', Living Reviews in Relativity, Section 3.2 'The parametrized post-Newtonian formalism', http://www.livingreviews.org/lrr-2006-3
    // Basic : true
    // Scalar: true
    static const double PPN_ZETA_4 = 0.;

    // Parsec expressed in au
    // Basic : false
    // Scalar: true
    static const double PARSEC_ASTRONOMICALUNIT = 206264.806247096; // [au]

    // Parsec expressed in m
    // Basic : false
    // Scalar: true
    static const double PARSEC_METER = 3.0856775814913674e+16; // [m]

    // The constant Pi (also known as Archimedes' constant). Note that double-precision, floating-point numbers in any programming language which follows the IEEE standard (true for C, C++, Java, and most, if not all, others) have only 16 significant digits (64 bits); the representation here, using 30 significant digits, is thus amply sufficient
    // Source: Well-known mathematical constant; numerical value can be extracted, e.g., from Mathematica 4.0 for Solaris (Wolfram Research, Inc.) using 'N[Pi,30]'
    // Basic : true
    // Scalar: true
    static const double PI_CONSTANT = 3.14159265358979323846264338328;

    // Planck's constant
    // Source: P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0)
    // Basic : true
    // Scalar: true
    static const double PLANCK_CONSTANT = 6.6260700400e-34; // [J s]

    // Astrometric signature of the Sun induced by the Pluto system for an observer located at a distance of 10 pc from the Sun
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.7, Equation 1.22, page 11
    // Basic : false
    // Scalar: true
    static const double PLUTOSYSTEM_ASTROMETRICSIGNATURE_10PARSEC = 0.029; // [10^-6 arcsec]

    // Pluto-system mass (IAU 2009 CBE value). The 'planetary' mass includes the contribution of its satellite, Charon
    // Basic : false
    // Scalar: true
    static const double PLUTOSYSTEM_MASS = 1.4561e+22; // [kg]

    // Mean orbital eccentricity of Pluto, at the standard epoch J2000.0. The mean orbital eccentricity is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double PLUTOSYSTEM_ORBITALECCENTRICITY_J2000 = 0.24882730;

    // Mean orbital inclination of Pluto, at the standard epoch J2000.0. The mean orbital inclination is is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double PLUTOSYSTEM_ORBITALINCLINATION_J2000 = 17.14001206; // [deg]

    // Sidereal orbital period
    // Source: Values derived from the mean longitude rates in Table 5.8.1 in P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, page 316; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double PLUTOSYSTEM_ORBITALPERIOD = 247.92065; // [yr]

    // Mean orbital semi-major axis of Pluto, at the standard epoch J2000.0. The mean orbital semi-major axis is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double PLUTOSYSTEM_ORBITALSEMIMAJORAXIS_J2000 = 39.48211675; // [au]

    // Radial-velocity amplitude of the Sun induced by the Pluto system for 'an edge-on observer' (i.e., an observer in the orbital plane of the Pluto system)
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.4, Equation 1.18, page 9
    // Basic : false
    // Scalar: true
    static const double PLUTOSYSTEM_RADIALVELOCITYSIGNATURE = 3.58e-05; // [m s^-1]

    // Radius of the smallest hypothetical sphere around Pluto which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double PLUTO_ENCOMPASSINGSPHERERADIUS = 1.1950e+06; // [m]

    // Equatorial radius of Pluto
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double PLUTO_EQUATORIALRADIUS = 1.1950e+06; // [m]

    // Geometrical flattening factor f of Pluto (f = (a-b)/a)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double PLUTO_FLATTENING = 0.;

    // Maximum reduction of the solar flux for an observer external to the solar system during a transit of Pluto
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.2, Equation 2.4, page 14
    // Basic : false
    // Scalar: true
    static const double PLUTO_FLUXREDUCTION_MAXIMUM = 0.0003; // [%]

    // Geometric albedo of Pluto (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double PLUTO_GEOMETRICALBEDO = 0.3;

    // Dynamical form-factor of Pluto (oblateness or Stokes' second-degree zonal harmonic of the potential); value is unknown
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8
    // Basic : true
    // Scalar: true
    static const double PLUTO_JSUB2 = 0.;

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of Pluto
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double PLUTO_LIGHTDEFLECTION_LIMB = 7.; // [10^-6 arcsec]

    // Mass of Pluto (do not use for high-precision (orbit) calculations)
    // Source: R.A. Jacobson, 2007, 'The orbits of the satellites of Pluto - Ephemeris PLU017', priv. comm.; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double PLUTO_MASS = 1.3090e+22; // [kg]

    // Mean mass density of Pluto (rough estimate)
    // Basic : false
    // Scalar: true
    static const double PLUTO_MASSDENSITY_MEAN = 1.83; // [g cm^-3]

    // Polar radius of Pluto
    // Basic : false
    // Scalar: true
    static const double PLUTO_POLARRADIUS = 1.1950e+06; // [m]

    // Geometric transit probability (Pluto transiting the Sun) for an observer external to the solar system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.1, Equation 2.2, page 14
    // Basic : false
    // Scalar: true
    static const double PLUTO_TRANSITPROBABILITY = 0.012; // [%]

    // Maximum transit time of Pluto (transiting the Sun) for an observer external to the solar system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.3, Equation 2.5, page 15
    // Basic : false
    // Scalar: true
    static const double PLUTO_TRANSITTIME_MAXIMUM = 3.40; // [day]

    // V(1,0) magnitude of Pluto (i.e., the visual magnitude of the 'planet' reduced to a distance of 1 au from both the Sun and Earth and phase angle zero). This parameter is also refered to as absolute magnitude in planetary sciences
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double PLUTO_VONEZEROMAGNITUDE = -1.0; // [mag]

    // Mean volumetric radius of Pluto
    // Basic : false
    // Scalar: true
    static const double PLUTO_VOLUMETRICRADIUS = 1.1950e+06; // [m]

    // Speed of general precession in ecliptic longitude, in arcsec per Julian century, at the standard epoch J2000.0, nowadays irrelevant as a fundamental constant
    // Source: N. Capitaine, P.T. Wallace, J. Chapront, 2003, 'Expressions for IAU 2000 precession quantities', A&A, 412, 567-586, 'P03 solution'
    // Basic : true
    // Scalar: true
    static const double PRECESSIONLONGITUDE_CONSTANT_J2000 = 5028.796195; // [arcsec cy^-1]

    // Precession constant m = p cos(\epsilon_0), in s per Julian year, at the standard epoch J2000.0, nowadays irrelevant as a fundamental constant. The precession in right ascension \alpha equals m + n sin(\alpha) tan(\delta); the precession in declination \delta equals n cos(\alpha)
    // Basic : false
    // Scalar: true
    static const double PRECESSION_CONSTANT_J2000M = 3.075887; // [s yr^-1]

    // Precession constant n = p sin(\epsilon_0), in arcsec per Julian year, at the standard epoch J2000.0, nowadays irrelevant as a fundamental constant. The precession in right ascension \alpha equals m + n sin(\alpha) tan(\delta); the precession in declination \delta equals n cos(\alpha)
    // Basic : false
    // Scalar: true
    static const double PRECESSION_CONSTANT_J2000N = 20.003394; // [arcsec yr^-1]

    // Proper-motion conversion constant A_v = 4.74... km yr s^-1 (see ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, page 25, Table 1.2.2)
    // Basic : false
    // Scalar: true
    static const double PROPERMOTION_CONSTANT = 4.740470464; // [km yr s^-1]

    // Composite (unweighted-mean) zero-redshift quasar (QSO) spectrum, based on 949 LBQS QSOs, 191 2MASS AGN, 37 Hamburg-ESO QSOs, and 18 PG QSOs. First column: wavelength \lambda (in nm; from 96.10 to 932.85). Second column: unnormalised flux density f_\lambda (in W m^-2 nm^-1)
    // Source: P.J. Francis, et al., 1991, 'A high signal-to-noise ratio composite quasar spectrum', Astrophysical Journal (ApJ), 373, 465
    // Basic : true
    // Scalar: false
    static const char  *const QSO_SPECTRUM_LAMBDA() { return "Nature/QSO_Spectrum_Lambda_001.fits"; }

    // One radian in units of arcseconds
    // Basic : false
    // Scalar: true
    static const double RADIAN_ARCSECOND = 2.062648062470964e+05; // [arcsec]

    // One radian in units of degrees
    // Basic : false
    // Scalar: true
    static const double RADIAN_DEGREE = 5.729577951308232e+01; // [deg]

    // One radian in units of micro-arcseconds
    // Basic : false
    // Scalar: true
    static const double RADIAN_MICROARCSECOND = 2.062648062470964e+11; // [micro-arcsec]

    // One radian in units of milli-arcseconds
    // Basic : false
    // Scalar: true
    static const double RADIAN_MILLIARCSECOND = 2.062648062470964e+08; // [milli-arcsec]

    // Radiation constant, also known as radiation-density constant, linking the energy density u (= 4 Pi I / c) of black-body radiation and temperature T via u = a T^4
    // Source: E.g., H. Karttunen, et al., 1987, 'Fundamental Astronomy', Springer Verlag, Berlin, Section 11.2, page 247 or R. Kippenhahn, A. Weigert, 1991, 'Stellar structure and evolution' (corrected 2-nd printing), Springer Verlag, Berlin, Section 3.1, page 16, and Section 5.1.2, page 28
    // Basic : false
    // Scalar: true
    static const double RADIATION_CONSTANT = 7.5657229e-16; // [J m^-3 K^-4]

    // First radiation constant. Note: best-measured value equals 3.741771790E-16 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))
    // Basic : false
    // Scalar: true
    static const double RADIATION_CONSTANT_FIRST = 3.7417717901e-16; // [W m^2]

    // Second radiation constant. Note: best-measured value equals 1.43877736E-2 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))
    // Basic : false
    // Scalar: true
    static const double RADIATION_CONSTANT_SECOND = 1.438777364e-02; // [m K]

    // The origin of TCB is defined in terms of TAI: the reading of TCB on 1 January 1977, 00:00:00 TAI = JD2443144.5 TAI must be 1 January 1977, 00:00:32.184 TCB = JD2443144.5003725 TCB. This origin has been arbitrarily set so that TCB coincides with TT at the geocentre on 1 January 1977, 00:00:00 TAI
    // Source: IAU, July 1991, 'Recommendations from the Working Group on Reference Systems', IAU 1991 Resolution A4, Recommendation 4, adopted at the XXI-st General Assembly of the IAU. See, for example, M. Soffel, et al., 1 December 2003, 'The IAU 2000 resolutions for astrometry, celestial mechanics, and metrology in the relativistic framework: explanatory supplement', AJ, 126, 2687-2706
    // Basic : true
    // Scalar: true
    static const char  *const REFERENCEEPOCH_TCB() { return "JD2443144.5003725 TCB"; }

    // The origin of TCG is defined in terms of TAI: the reading of TCG on 1 January 1977, 00:00:00 TAI = JD2443144.5 TAI must be 1 January 1977, 00:00:32.184 TCG = JD2443144.5003725 TCG. This origin has been arbitrarily set so that TCG coincides with TT at the geocentre on 1 January 1977, 00:00:00 TAI
    // Source: IAU, July 1991, 'Recommendations from the Working Group on Reference Systems', IAU 1991 Resolution A4, Recommendation 4, adopted at the XXI-st General Assembly of the IAU. See, for example, M. Soffel, et al., 1 December 2003, 'The IAU 2000 resolutions for astrometry, celestial mechanics, and metrology in the relativistic framework: explanatory supplement', AJ, 126, 2687-2706
    // Basic : true
    // Scalar: true
    static const char  *const REFERENCEEPOCH_TCG() { return "JD2443144.5003725 TCG"; }

    // IAU 2006 Resolution B3, entitled 'Re-definition of Barycentric Dynamical Time, TDB', recommends that 'TDB be defined as the following linear transformation of TCB: TDB = TCB - L_B x ( JD_TCB - T_0 ) x 86400 + TDB_0, where T_0 = 2443144.5003725 (parameter :Nature:ReferenceEpoch_TCB), and L_B = 1.550519768E-8 (parameter :Nature:LSubB_Constant) and TDB_0 = -6.55E-5 s are defining constants'. The number 86400 is equal to parameter :Nature:Day_Second
    // Source: IAU, August 2006, 'Re-definition of Barycentric Dynamical Time, TDB', IAU 2006 Resolution 3 adopted at the XXVI-th General Assembly of the IAU. See also the IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double REFERENCEEPOCH_TDBSUBZERO = -6.550e-05; // [s]

    // Rydberg constant. Note: best-measured value equals 10973731.568508 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))
    // Basic : false
    // Scalar: true
    static const double RYDBERG_CONSTANT = 10973731.570551; // [m^-1]

    // Velocity vector of the Solar-System BaryCentre (SSBC) with respect to the Cosmic Microwave Background (CMB), in units of km s^-1. The vector elements refer to Galactic (U,V,W) coordinates
    // Source: G. Hinshaw, et al., 11 February 2009, 'Five-Year Wilkinson Microwave Anisotropy Probe (WMAP) Observations: Data Processing, Sky Maps, and Basic Results', Astrophysical Journal Supplement (ApJS), Volume 180, pages 225-245
    // Basic : true
    // Scalar: true
    static const double  *const SSBC_VELOCITYCMB ()  { static double _v[3] = { -26.29,  -244.96,  275.93 }; return _v; } // [km s^-1]

    // Astrometric signature of the Sun induced by the Saturn system for an observer located at a distance of 10 pc from the Sun
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.7, Equation 1.22, page 11
    // Basic : false
    // Scalar: true
    static const double SATURNSYSTEM_ASTROMETRICSIGNATURE_10PARSEC = 273.; // [10^-6 arcsec]

    // Saturn-system mass (IAU 2009 CBE value). The planetary mass includes the contribution of its satellites
    // Basic : false
    // Scalar: true
    static const double SATURNSYSTEM_MASS = 5.68477e+26; // [kg]

    // Mean orbital eccentricity of Saturn, at the standard epoch J2000.0. The mean orbital eccentricity is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double SATURNSYSTEM_ORBITALECCENTRICITY_J2000 = 0.05386179;

    // Mean orbital inclination of Saturn, at the standard epoch J2000.0. The mean orbital inclination is is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double SATURNSYSTEM_ORBITALINCLINATION_J2000 = 2.48599187; // [deg]

    // Sidereal orbital period
    // Source: Values derived from the mean longitude rates in Table 5.8.1 in P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, page 316; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double SATURNSYSTEM_ORBITALPERIOD = 29.447498; // [yr]

    // Mean orbital semi-major axis of Saturn, at the standard epoch J2000.0. The mean orbital semi-major axis is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double SATURNSYSTEM_ORBITALSEMIMAJORAXIS_J2000 = 9.53667594; // [au]

    // Radial-velocity amplitude of the Sun induced by the Saturn system for 'an edge-on observer' (i.e., an observer in the orbital plane of the Saturn system)
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.4, Equation 1.18, page 9
    // Basic : false
    // Scalar: true
    static const double SATURNSYSTEM_RADIALVELOCITYSIGNATURE = 2.8; // [m s^-1]

    // Radius of the smallest hypothetical sphere around Saturn which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double SATURN_ENCOMPASSINGSPHERERADIUS = 6.02680e+07; // [m]

    // Equatorial radius of Saturn
    // Basic : false
    // Scalar: true
    static const double SATURN_EQUATORIALRADIUS = 6.02680e+07; // [m]

    // Geometrical flattening factor f of Saturn (f = (a-b)/a)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double SATURN_FLATTENING = 9.796240e-02;

    // Maximum reduction of the solar flux for an observer external to the solar system during a transit of Saturn
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.2, Equation 2.4, page 14
    // Basic : false
    // Scalar: true
    static const double SATURN_FLUXREDUCTION_MAXIMUM = 0.701; // [%]

    // Geometric albedo of Saturn (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double SATURN_GEOMETRICALBEDO = 0.47;

    // Dynamical form-factor of Saturn (oblateness or Stokes' second-degree zonal harmonic of the potential)
    // Source: P.R. Weissman, L.-A. McFadden, T.V. Johnson (eds.), 1999, 'Encyclopedia of the Solar System (first edition)', Academic Press, page 342
    // Basic : true
    // Scalar: true
    static const double SATURN_JSUB2 = 0.016331;

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of Saturn
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double SATURN_LIGHTDEFLECTION_LIMB = 5980.; // [10^-6 arcsec]

    // Mass of Saturn (do not use for high-precision (orbit) calculations)
    // Source: R.A. Jacobson, et al., 2006, 'The gravity field of the Saturnian system from satellite observations and spacecraft tracking data', AJ, 132, 2520-2526; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double SATURN_MASS = 5.683190e+26; // [kg]

    // Mean mass density of Saturn
    // Basic : false
    // Scalar: true
    static const double SATURN_MASSDENSITY_MEAN = 0.6871; // [g cm^-3]

    // IAU-recommended value for the declination \delta_0 of the north pole of rotation of Saturn. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double SATURN_NORTHROTATIONALPOLE_DECLINATION = 83.537; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the declination \delta_0 of the north pole of rotation of Saturn. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : false
    // Scalar: true
    static const double SATURN_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE = -0.00000011; // [deg day^-1]

    // IAU-recommended value for the right ascension \alpha_0 of the north pole of rotation of Saturn. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double SATURN_NORTHROTATIONALPOLE_RIGHTASCENSION = 40.589; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the right ascension \alpha_0 of the north pole of rotation of Saturn. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : false
    // Scalar: true
    static const double SATURN_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE = -0.00000099; // [deg day^-1]

    // Polar radius of Saturn
    // Basic : false
    // Scalar: true
    static const double SATURN_POLARRADIUS = 5.43640e+07; // [m]

    // IAU-recommended value for the ephemeris position of the prime meridian of Saturn. The prime meridian refers to the rotation of the magnetic field System III. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double SATURN_PRIMEMERIDIAN_EPHEMERISPOSITION = 38.90; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch of 1.5 January 2000 = JD2451545.0 TDB) of the ephemeris position of the prime meridian of Saturn. The prime meridian refers to the rotation of the magnetic field System III. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double SATURN_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE = 810.7939024; // [deg day^-1]

    // Geometric transit probability (Saturn transiting the Sun) for an observer external to the solar system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.1, Equation 2.2, page 14
    // Basic : false
    // Scalar: true
    static const double SATURN_TRANSITPROBABILITY = 0.053; // [%]

    // Maximum transit time of Saturn (transiting the Sun) for an observer external to the solar system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.3, Equation 2.5, page 15
    // Basic : false
    // Scalar: true
    static const double SATURN_TRANSITTIME_MAXIMUM = 1.81; // [day]

    // V(1,0) magnitude of Saturn (i.e., the visual magnitude of the planet reduced to a distance of 1 au from both the Sun and Earth and phase angle zero). This parameter is also refered to as absolute magnitude in planetary sciences. The magnitude refers to the planetary disk only
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double SATURN_VONEZEROMAGNITUDE = -8.88; // [mag]

    // Mean volumetric radius of Saturn
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double SATURN_VOLUMETRICRADIUS = 5.82320e+07; // [m]

    // Shapiro's time-delay constant for the Sun
    // Source: See, e.g., L. Lindegren, D. Dravins, 2003, 'The fundamental definition of <radial velocity>', A&A, 401, 1185, Section 4.3 or P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Equation 5.3211-1, page 295
    // Basic : false
    // Scalar: true
    static const double SHAPIRO_CONSTANT = 9.8510e-06; // [s]

    // The bulk modulus of SiC (Boostec) under isotropic stress; also known as compression modulus
    // Basic : false
    // Scalar: true
    static const double SIC_BULKMODULUS = 206.; // [GPa]

    // The compressibility of SiC (Boostec)
    // Basic : false
    // Scalar: true
    static const double SIC_COMPRESSIBILITY = 4.857e-12; // [m^2 N^-1]

    // SiC cryogenic linear scale factor \xi between 293 K and 120 K. The assumed linear thermal expansion coefficient of SiC at 293 K is 1.27 ppm K^-1 (see also parameter SiC_LinearThermalCoefficientOfExpansion_293K)
    // Source: A. Mora, 21 June 2011, 'Conversion between image and object space coordinates', GAIA-CH-TN-ESAC-AMF-008, issue 2, revision 1, Section 3
    // Basic : false
    // Scalar: true
    static const double SIC_CRYOGENICLINEARSCALEFACTOR_293K = 1.00021971;

    // First Lame constant of SiC (Boostec)
    // Basic : false
    // Scalar: true
    static const double SIC_LAMECONSTANT_FIRST = 85.; // [GPa]

    // Second Lame constant of SiC (Boostec), also known as shear modulus or rigidity
    // Basic : false
    // Scalar: true
    static const double SIC_LAMECONSTANT_SECOND = 181.; // [GPa]

    // Average linear thermal expansion coefficient of isotropic homogeneous SiC (Boostec) at 100 K
    // Basic : false
    // Scalar: true
    static const double SIC_LINEARTHERMALCOEFFICIENTOFEXPANSION_100K = 1.0; // [ppm K^-1]

    // Average linear thermal expansion coefficient of isotropic homogeneous SiC (Boostec) at 170 K
    // Basic : false
    // Scalar: true
    static const double SIC_LINEARTHERMALCOEFFICIENTOFEXPANSION_170K = 1.1; // [ppm K^-1]

    // Average linear thermal expansion coefficient of isotropic homogeneous SiC (Boostec) at 293 K
    // Basic : false
    // Scalar: true
    static const double SIC_LINEARTHERMALCOEFFICIENTOFEXPANSION_293K = 0.7; // [ppm K^-1]

    // Density of SiC (Boostec)
    // Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9
    // Basic : true
    // Scalar: true
    static const double SIC_MASSDENSITY = 3.16; // [g cm^-3]

    // The P-wave speed in SiC (Boostec); a P-wave (pressure wave) is a longitudinal wave in an elastic medium in which the restoring force is provided by the medium's bulk modulus
    // Basic : false
    // Scalar: true
    static const double SIC_PWAVESPEED = 11.90; // [km s^-1]

    // Poisson ratio of SiC (Boostec)
    // Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9
    // Basic : true
    // Scalar: true
    static const double SIC_POISSONRATIO = 0.16;

    // The S-wave speed in SiC (Boostec); an S-wave is a wave in an elastic medium in which the restoring force is provided by shear
    // Basic : false
    // Scalar: true
    static const double SIC_SWAVESPEED = 7.57; // [km s^-1]

    // Specific heat of SiC (Boostec) at 100 K
    // Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9
    // Basic : true
    // Scalar: true
    static const double SIC_SPECIFICHEATATCONSTANTPRESSURE_100K = 100.; // [J K^-1 kg^-1]

    // Specific heat of SiC (Boostec) at 170 K
    // Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9
    // Basic : true
    // Scalar: true
    static const double SIC_SPECIFICHEATATCONSTANTPRESSURE_170K = 320.; // [J K^-1 kg^-1]

    // Specific heat of SiC (Boostec) at 293 K
    // Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9
    // Basic : true
    // Scalar: true
    static const double SIC_SPECIFICHEATATCONSTANTPRESSURE_293K = 680.; // [J K^-1 kg^-1]

    // Specific stiffness of SiC (Boostec); also known as specific modulus
    // Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9
    // Basic : false
    // Scalar: true
    static const double SIC_SPECIFICSTIFNESS = 133.; // [10^6 N m kg^-1]

    // Average strength of SiC (Boostec) in the four-point bending test
    // Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9
    // Basic : true
    // Scalar: true
    static const double SIC_STRENGTH_AVERAGE = 390.; // [10^6 Pa]

    // Thermal conductivity of isotropic homogeneous SiC (Boostec) at 100 K
    // Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9
    // Basic : true
    // Scalar: true
    static const double SIC_THERMALCONDUCTIVITY_100K = 180.; // [W m^-1 K^-1]

    // Thermal conductivity of isotropic homogeneous SiC (Boostec) at 170 K
    // Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9
    // Basic : true
    // Scalar: true
    static const double SIC_THERMALCONDUCTIVITY_170K = 220.; // [W m^-1 K^-1]

    // Thermal conductivity of isotropic homogeneous SiC (Boostec) at 293 K
    // Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9
    // Basic : true
    // Scalar: true
    static const double SIC_THERMALCONDUCTIVITY_293K = 190.; // [W m^-1 K^-1]

    // Thermal diffusivity of SiC (Boostec) at 100 K
    // Basic : false
    // Scalar: true
    static const double SIC_THERMALDIFFUSIVITY_100K = 570.; // [m^2 s^-1]

    // Thermal diffusivity of SiC (Boostec) at 170 K
    // Basic : false
    // Scalar: true
    static const double SIC_THERMALDIFFUSIVITY_170K = 218.; // [m^2 s^-1]

    // Thermal diffusivity of SiC (Boostec) at 293 K
    // Basic : false
    // Scalar: true
    static const double SIC_THERMALDIFFUSIVITY_293K = 88.; // [m^2 s^-1]

    // Average volumetric thermal expansion coefficient of isotropic homogeneous SiC (Boostec) at 100 K
    // Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9
    // Basic : true
    // Scalar: true
    static const double SIC_VOLUMETRICTHERMALCOEFFICIENTOFEXPANSION_100K = 3.1; // [ppm K^-1]

    // Average volumetric thermal expansion coefficient of isotropic homogeneous SiC (Boostec) at 170 K
    // Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9
    // Basic : true
    // Scalar: true
    static const double SIC_VOLUMETRICTHERMALCOEFFICIENTOFEXPANSION_170K = 3.4; // [ppm K^-1]

    // Average volumetric thermal expansion coefficient of isotropic homogeneous SiC (Boostec) at 293 K
    // Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9
    // Basic : true
    // Scalar: true
    static const double SIC_VOLUMETRICTHERMALCOEFFICIENTOFEXPANSION_293K = 2.2; // [ppm K^-1]

    // The Weibull modulus of fired/sintered SiC (Boostec). The Weibull modulus is the slope of the linear plot of log(log(P^-1)) versus log(\sigma), where P is the probability of survival under a given stress \sigma. The numerical value is a lower limit for the Weibull modulus
    // Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9
    // Basic : true
    // Scalar: true
    static const double SIC_WEIBULLMODULUS_FIREDSLASHSINTERED = 12.;

    // The Weibull modulus of machined SiC (Boostec). The Weibull modulus is the slope of the linear plot of log(log(P^-1)) versus log(\sigma), where P is the probability of survival under a given stress \sigma. The numerical value is a lower limit for the Weibull modulus
    // Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9
    // Basic : true
    // Scalar: true
    static const double SIC_WEIBULLMODULUS_MACHINED = 10.;

    // Young's modulus of SiC (Boostec); a measure for the resistance of a material to elastic (recoverable) deformation under load. Also known as elastic modulus and tension modulus
    // Source: F. Safa (EADS-Astrium), 23 April 2002, 'SiC versus C-SiC trade-off', pages 8-9
    // Basic : true
    // Scalar: true
    static const double SIC_YOUNGMODULUS = 420.; // [GPa]

    // Silicon bandgap, in eV, at 0 K
    // Source: B. van Zeghbroeck, 1 January 2007, 'Principles of Semiconductor Devices', http://ece-www.colorado.edu/~bart/book/
    // Basic : true
    // Scalar: true
    static const double SI_BANDGAP_0K = 1.166; // [eV]

    // Silicon bandgap, in eV, at 160 K. The energy bandgap of Silicon (and semi-conductors in general) decreases with increasing temperature. This is explained as follows: an increased temperature enhances the amplitude of atomic vibrations due to the increased thermal energy; this causes the interatomic spacing to increase; this decreases the potential seen by the electrons in the material; this, finally, reduces the size of the energy bandgap. The temperature dependence of the energy bandgap has been experimentally determined to be E_g(T) = E_g(0) - \alpha * T^2 * (T + \beta)^-1, where T is the temperature in K,  E_g(0) is the bandgap energy in eV at zero absolute temperature (parameter Si_Bandgap_0K), and \alpha and \beta are material constants (parameters Si_Constant_Alpha and Si_Constant_Beta for Silicon)
    // Basic : false
    // Scalar: true
    static const double SI_BANDGAP_160K = 1.151; // [eV]

    // Silicon constant \alpha, in eV K^-1. The energy bandgap of Silicon (and semi-conductors in general) decreases with increasing temperature. This is explained as follows: an increased temperature enhances the amplitude of atomic vibrations due to the increased thermal energy; this causes the interatomic spacing to increase; this decreases the potential seen by the electrons in the material; this, finally, reduces the size of the energy bandgap. The temperature dependence of the energy bandgap has been experimentally determined to be E_g(T) = E_g(0) - \alpha * T^2 * (T + \beta)^-1, where T is the temperature in K,  E_g(0) is the bandgap energy in eV at zero absolute temperature (parameter Si_Bandgap_0K), and \alpha and \beta are material constants (parameters Si_Constant_Alpha and Si_Constant_Beta for Silicon)
    // Source: B. van Zeghbroeck, 1 January 2007, 'Principles of Semiconductor Devices', http://ece-www.colorado.edu/~bart/book/
    // Basic : true
    // Scalar: true
    static const double SI_CONSTANT_ALPHA = 4.730e-04; // [eV K^-1]

    // Silicon constant \beta, in K. The energy bandgap of Silicon (and semi-conductors in general) decreases with increasing temperature. This is explained as follows: an increased temperature enhances the amplitude of atomic vibrations due to the increased thermal energy; this causes the interatomic spacing to increase; this decreases the potential seen by the electrons in the material; this, finally, reduces the size of the energy bandgap. The temperature dependence of the energy bandgap has been experimentally determined to be E_g(T) = E_g(0) - \alpha * T^2 * (T + \beta)^-1, where T is the temperature in K,  E_g(0) is the bandgap energy in eV at zero absolute temperature (parameter Si_Bandgap_0K), and \alpha and \beta are material constants (parameters Si_Constant_Alpha and Si_Constant_Beta for Silicon)
    // Source: B. van Zeghbroeck, 1 January 2007, 'Principles of Semiconductor Devices', http://ece-www.colorado.edu/~bart/book/
    // Basic : true
    // Scalar: true
    static const double SI_CONSTANT_BETA = 636.; // [K]

    // Silicon cut-off wavelength, in nm, at 160 K
    // Basic : false
    // Scalar: true
    static const double SI_CUTOFFWAVELENGTH_160K = 1077.39; // [nm]

    // Silicon diffusion coefficient
    // Source: J.R. Janesick, 2001, 'Scientific CCDs', SPIE, Bellingham, Washington, Example 4.17, page 348
    // Basic : true
    // Scalar: true
    static const double SI_DIFFUSIONCOEFFICIENT = 0.0039; // [m^2 s^-1]

    // Diffusion length in the epitaxial Silicon
    // Source: J.R. Janesick, 2001, 'Scientific CCDs', SPIE, Bellingham, Washington, Example 4.17, page 348
    // Basic : true
    // Scalar: true
    static const double SI_DIFFUSIONLENGTH = 0.0006; // [m]

    // Silicon optical absorption coefficient as a function of (photon) wavelength. First column: wavelength \lambda (in nm; from 200.0 to 1100.0). Second column: Silicon optical absorption coefficient \alpha (in [10^-6 m]^-1) at T = 160 K. Third column: Silicon photon absorption depth L_A = \alpha^-1 (in 10^-6 m) at T = 160 K
    // Source: K. Rajkanan, R. Singh, J. Shewchun, 1979, 'Absorption coefficient of Silicon for solar cell calculations', Solid-State Electronics, 22, 793-795; analytical phenomenological model presented in this reference has an accuracy of about 20%; formal validity ranges are 20-500 K in temperature T and 1.1-4.0 eV, i.e., 310-1127 nm, in photon energy (and wavelength)
    // Basic : true
    // Scalar: false
    static const char  *const SI_OPTICALABSORPTIONCOEFFICIENT_160K() { return "Nature/Si_OpticalAbsorptionCoefficient_160K_001.fits"; }

    // Silicon photon absorption depth as a function of (photon) wavelength. First column: wavelength \lambda (in nm; from 200.0 to 1100.0). Second column: Silicon optical absorption coefficient \alpha (in [10^-6 m]^-1) at T = 160 K. Third column: Silicon photon absorption depth L_A = \alpha^-1 (in 10^-6 m) at T = 160 K
    // Source: K. Rajkanan, R. Singh, J. Shewchun, 1979, 'Absorption coefficient of Silicon for solar cell calculations', Solid-State Electronics, 22, 793-795; analytical phenomenological model presented in this reference has an accuracy of about 20%; formal validity ranges are 20-500 K in temperature T and 1.1-4.0 eV, i.e., 310-1127 nm, in photon energy (and wavelength)
    // Basic : true
    // Scalar: false
    static const char  *const SI_PHOTONABSORPTIONDEPTH_160K() { return "Nature/Si_PhotonAbsorptionDepth_160K_001.fits"; }

    // Number of days per sidereal year
    // Source: J.L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touze, G. Francou, J. Laskar, 1994, 'Numerical expressions for precession formulae and mean elements for the Moon and the planets', A&A, 282, 663 (1994A&A...282..663S)
    // Basic : true
    // Scalar: true
    static const double SIDEREALYEAR_J2000DAY = 365.256363004; // [day]

    // Cumulative number of galaxies, integrated over the full sky, as function of G magnitude for the following limits: up to G = 4.5 mag, up to G = 5.0 mag, ..., up to G = 20.5 mag, and up to G = 21.0 mag (numerical values from parameter Sky_ObjectDensity_003, extracted by B. Holl on 4 October 2012 using GT SkyDensity.class). See parameter Sky_NumberOfStars_G for star counts and Sky_NumberOfObjects_G for object counts
    // Source: The data used to produce the standard Gaia Galaxy model have been provided by the teams of the Observatoire de Besancon (A. Robin, C. Reyle, et al.) and the Observatorio Astronomico di Torino (R. Drimmel, et al.); the harmonisation has been provided by the SWG (X. Luri). See U. Lammers, 22 April 2005, 'Gaia Standard Galaxy Model Access Software (GSGMAS); User Guide', GAIA-UL-010, issue 1, revision 0. See also R. Drimmel, et al., 19 July 2005, 'Recommendations on the use of estimated star counts for Gaia studies', SWG-RD-004
    // Basic : true
    // Scalar: true
    static const double  *const SKY_NUMBEROFGALAXIES_G ()  { static double _v[34] = { 0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.009,  0.038,  0.104,  0.241,  0.509,  0.983,  1.822,  3.326,  6.048,  11.075,  21.145,  41.602,  82.794,  166.012 }; return _v; } // [10^6 objects]

    // Cumulative number of galaxies, integrated over the full sky, as function of G_RVS (= C1M861 = RVF) magnitude for the following limits: up to G_RVS = 4.5 mag, up to G_RVS = 5.0 mag, ..., up to G_RVS = 17.5 mag, and up to G_RVS = 18.0 mag (numerical values from parameter Sky_ObjectDensity_003, extracted by B. Holl on 4 October 2012 using GT SkyDensity.class). See parameter Sky_NumberOfStars_GRVS for star counts and Sky_NumberOfObjects_GRVS for object counts
    // Source: The data used to produce the standard Gaia Galaxy model have been provided by the teams of the Observatoire de Besancon (A. Robin, C. Reyle, et al.) and the Observatorio Astronomico di Torino (R. Drimmel, et al.); the harmonisation has been provided by the SWG (X. Luri). See U. Lammers, 22 April 2005, 'Gaia Standard Galaxy Model Access Software (GSGMAS); User Guide', GAIA-UL-010, issue 1, revision 0. See also R. Drimmel, et al., 19 July 2005, 'Recommendations on the use of estimated star counts for Gaia studies', SWG-RD-004
    // Basic : true
    // Scalar: true
    static const double  *const SKY_NUMBEROFGALAXIES_GRVS ()  { static double _v[28] = { 0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.000,  0.010,  0.040,  0.111,  0.255,  0.537,  1.032,  1.912,  3.490,  6.349,  11.607 }; return _v; } // [10^6 objects]

    // Cumulative number of objects (stars + galaxies; see parameters Sky_NumberOfStars_G and Sky_NumberOfGalaxies_G), integrated over the full sky, as function of G magnitude for the following limits: up to G = 4.5 mag, up to G = 5.0 mag, ..., up to G = 20.5 mag, and up to G = 21.0 mag (numerical values from parameter Sky_ObjectDensity_003, extracted by B. Holl on 4 October 2012 using GT SkyDensity.class)
    // Source: The data used to produce the standard Gaia Galaxy model have been provided by the teams of the Observatoire de Besancon (A. Robin, C. Reyle, et al.) and the Observatorio Astronomico di Torino (R. Drimmel, et al.); the harmonisation has been provided by the SWG (X. Luri). See U. Lammers, 22 April 2005, 'Gaia Standard Galaxy Model Access Software (GSGMAS); User Guide', GAIA-UL-010, issue 1, revision 0. See also R. Drimmel, et al., 19 July 2005, 'Recommendations on the use of estimated star counts for Gaia studies', SWG-RD-004
    // Basic : true
    // Scalar: true
    static const double  *const SKY_NUMBEROFOBJECTS_G ()  { static double _v[34] = { 0.003,  0.007,  0.014,  0.026,  0.044,  0.070,  0.107,  0.160,  0.241,  0.356,  0.522,  0.768,  1.141,  1.715,  2.640,  4.053,  6.267,  9.708,  14.749,  22.060,  32.556,  47.488,  68.545,  97.837,  137.944,  192.408,  265.629,  363.261,  494.546,  672.011,  911.605,  1233.940,  1667.145,  2247.560 }; return _v; } // [10^6 objects]

    // Cumulative number of objects (stars + galaxies; see parameters Sky_NumberOfStars_GRVS and Sky_NumberOfGalaxies_GRVS), integrated over the full sky, as function of G_RVS (= C1M861 = RVF) magnitude for the following limits: up to G_RVS = 4.5 mag, up to G_RVS = 5.0 mag, ..., up to G_RVS = 17.5 mag, and up to G_RVS = 18.0 mag (numerical values from parameter Sky_ObjectDensity_003, extracted by B. Holl on 4 October 2012 using GT SkyDensity.class)
    // Source: The data used to produce the standard Gaia Galaxy model have been provided by the teams of the Observatoire de Besancon (A. Robin, C. Reyle, et al.) and the Observatorio Astronomico di Torino (R. Drimmel, et al.); the harmonisation has been provided by the SWG (X. Luri). See U. Lammers, 22 April 2005, 'Gaia Standard Galaxy Model Access Software (GSGMAS); User Guide', GAIA-UL-010, issue 1, revision 0. See also R. Drimmel, et al., 19 July 2005, 'Recommendations on the use of estimated star counts for Gaia studies', SWG-RD-004
    // Basic : true
    // Scalar: true
    static const double  *const SKY_NUMBEROFOBJECTS_GRVS ()  { static double _v[28] = { 0.018,  0.042,  0.072,  0.111,  0.160,  0.227,  0.320,  0.460,  0.675,  0.995,  1.471,  2.196,  3.293,  4.902,  7.280,  10.876,  16.184,  23.915,  35.157,  51.387,  74.697,  107.743,  154.027,  217.589,  303.945,  420.851,  580.228,  800.547 }; return _v; } // [10^6 objects]

    // Predicted number of observed stars (i.e., G <= 20.00 mag) in the disk (giants)
    // Source: ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 6.4.4 and Table 6.3, pages 239-240 (Galaxy model from J. Torra, et al., 1999, 'Predicting Gaia Observations from a Star-Count Model', Baltic Astronomy, 8, 171 and extinction law from J. Hakkila, et al., 1997, 'A Computerised Model of Large-Scale Visual Interstellar Extinction', AJ, 114, 2043)
    // Basic : true
    // Scalar: true
    static const double SKY_NUMBEROFSTARS_DISKGIANT = 92.; // [10^6 stars]

    // Predicted number of observed stars (i.e., G <= 20.00 mag) in the disk (main sequence stars plus white dwarfs)
    // Source: ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 6.4.4 and Table 6.3, pages 239-240 (Galaxy model from J. Torra, et al., 1999, 'Predicting Gaia Observations from a Star-Count Model', Baltic Astronomy, 8, 171 and extinction law from J. Hakkila, et al., 1997, 'A Computerised Model of Large-Scale Visual Interstellar Extinction', AJ, 114, 2043)
    // Basic : true
    // Scalar: true
    static const double SKY_NUMBEROFSTARS_DISKMSPLUSWD = 779.; // [10^6 stars]

    // Cumulative number of stars, integrated over the full sky, as function of G magnitude for the following limits: up to G = 4.5 mag, up to G = 5.0 mag, ..., up to G = 20.5 mag, and up to G = 21.0 mag (numerical values from parameter Sky_ObjectDensity_003, extracted by B. Holl on 4 October 2012 using GT SkyDensity.class). See parameter Sky_NumberOfGalaxies_G for galaxy counts and Sky_NumberOfObjects_G for object counts
    // Source: The data used to produce the standard Gaia Galaxy model have been provided by the teams of the Observatoire de Besancon (A. Robin, C. Reyle, et al.) and the Observatorio Astronomico di Torino (R. Drimmel, et al.); the harmonisation has been provided by the SWG (X. Luri). See U. Lammers, 22 April 2005, 'Gaia Standard Galaxy Model Access Software (GSGMAS); User Guide', GAIA-UL-010, issue 1, revision 0. See also R. Drimmel, et al., 19 July 2005, 'Recommendations on the use of estimated star counts for Gaia studies', SWG-RD-004
    // Basic : true
    // Scalar: true
    static const double  *const SKY_NUMBEROFSTARS_G ()  { static double _v[34] = { 0.003,  0.007,  0.014,  0.026,  0.044,  0.070,  0.107,  0.160,  0.241,  0.356,  0.522,  0.768,  1.141,  1.715,  2.640,  4.053,  6.267,  9.708,  14.749,  22.060,  32.547,  47.450,  68.441,  97.595,  137.435,  191.424,  263.806,  359.935,  488.497,  660.937,  890.460,  1192.338,  1584.350,  2081.548 }; return _v; } // [10^6 objects]

    // Cumulative number of stars, integrated over the full sky, as function of G_RVS (= C1M861 = RVF) magnitude for the following limits: up to G_RVS = 4.5 mag, up to G_RVS = 5.0 mag, ..., up to G_RVS = 17.5 mag, and up to G_RVS = 18.0 mag (numerical values from parameter Sky_ObjectDensity_003, extracted by B. Holl on 4 October 2012 using GT SkyDensity.class). See parameter Sky_NumberOfGalaxies_GRVS for galaxy counts and Sky_NumberOfObjects_GRVS for object counts
    // Source: The data used to produce the standard Gaia Galaxy model have been provided by the teams of the Observatoire de Besancon (A. Robin, C. Reyle, et al.) and the Observatorio Astronomico di Torino (R. Drimmel, et al.); the harmonisation has been provided by the SWG (X. Luri). See U. Lammers, 22 April 2005, 'Gaia Standard Galaxy Model Access Software (GSGMAS); User Guide', GAIA-UL-010, issue 1, revision 0. See also R. Drimmel, et al., 19 July 2005, 'Recommendations on the use of estimated star counts for Gaia studies', SWG-RD-004
    // Basic : true
    // Scalar: true
    static const double  *const SKY_NUMBEROFSTARS_GRVS ()  { static double _v[28] = { 0.018,  0.042,  0.072,  0.111,  0.160,  0.227,  0.320,  0.460,  0.675,  0.995,  1.471,  2.196,  3.293,  4.902,  7.280,  10.876,  16.184,  23.915,  35.147,  51.347,  74.587,  107.488,  153.490,  216.557,  302.033,  417.361,  573.879,  788.939 }; return _v; } // [10^6 objects]

    // Predicted number of observed stars (i.e., G <= 20.00 mag) in the spheroid (including the bulge)
    // Source: ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 6.4.4 and Table 6.3, pages 239-240 (Galaxy model from J. Torra, et al., 1999, 'Predicting Gaia Observations from a Star-Count Model', Baltic Astronomy, 8, 171 and extinction law from J. Hakkila, et al., 1997, 'A Computerised Model of Large-Scale Visual Interstellar Extinction', AJ, 114, 2043)
    // Basic : true
    // Scalar: true
    static const double SKY_NUMBEROFSTARS_SPHEROID = 67.; // [10^6 stars]

    // Predicted number of observed stars (i.e., G <= 20.00 mag) in the thick disk
    // Source: ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 6.4.4 and Table 6.3, pages 239-240 (Galaxy model from J. Torra, et al., 1999, 'Predicting Gaia Observations from a Star-Count Model', Baltic Astronomy, 8, 171 and extinction law from J. Hakkila, et al., 1997, 'A Computerised Model of Large-Scale Visual Interstellar Extinction', AJ, 114, 2043)
    // Basic : true
    // Scalar: true
    static const double SKY_NUMBEROFSTARS_THICKDISK = 97.; // [10^6 stars]

    // Predicted total number of observed stars (i.e., G <= 20.00 mag)
    // Source: See also ESA, 21 May 2013, 'Gaia mission requirements document (MRD)', GAIA-EST-RD-00553, issue 3, revision 1, Requirement SCI-130
    // Basic : false
    // Scalar: true
    static const double SKY_NUMBEROFSTARS_TOTAL = 1.0350e+09; // [stars]

    // Standard Gaia Galaxy model, providing predictions of object densities in the three Gaia-2-relevant photometric bands G, GS, and G_RVS (= C1M861 = RVF) as a function of limiting magnitude for any point on the celestial sphere. Object densities are defined on a Hierarchical Triangular Mesh of level 6 comprised of almost-equal-area cells of approximately 1 square degree in size. The model is a synthesis of count predictions from two distinct sources: (i) the Bescancon Galaxy model (http://model.obs-besancon.fr/), in conjunction with the extinction law of R. Drimmel, for stars; and (ii) the GSC-II catalogue (http://gsss.stsci.edu/Catalogs/GSC/GSC2/GSC2.htm) for stars and galaxies. The combined result can be summarised as follows: (a) Bright stars: 4 <= G <= 12.5: GSC-II (Tycho-2 catalogue included); 4 <= GS <= 12.5: GSC-II (Tycho-2 catalogue included); 4 <= G_RVS <= 12.5: GSC-II (Tycho-2 catalogue included); (b) Faint stars: 12.5 <= G <= 21: Bescancon; 12.5 <= GS <= 21: Bescancon; 12.5 <= G_RVS <= 18: Bescancon; (c) Stars in special regions (LMC/SMC): counts for full range of magnitudes, for circle with radius 10 deg (LMC) and 7 deg (SMC), taken from the GSC-II catalogue; (d) galaxies: G <= 21: GSC-II; GS <= 21: GSC-II; G_RVS <= 18: GSC-II. The following Bintables are defined: (1) STAR-G-MAGGRID: magnitude grid for column DENSITY in following extension STAR-G-DENSITY (the magnitude grids STAR-G-MAGGRID and STAR-C1M861-MAGGRID are defined through parameters Sky_ObjectDensity_MagMin, Sky_ObjectDensity_MagMax, and Sky_ObjectDensity_MagBinWidth); (2) STAR-G-DENSITY: table with star magnitude counts in column DENSITY per HTM cell in the G band; each row corresponds to one cell. Table columns are: HTMID: unique HTM cell designator; AREA: area of cell [square deg]; ALPHA: right ascension of cell centre [h]; DELTA: declination of cell centre [deg]; LGAL: Galactic longitude of cell centre [deg]; BGAL: Galactic latitude of cell centre [deg]; DENSITY: non-cumulative number of objects per square degree in cell in 0.5 magnitude-bin interval, i.e., DENSITY[i] is the number of stars per square degree in the magnitude interval [STAR-G-MAGGRID[i], STAR-G-MAGGRID[i]+0.5]; (3) GAL-G-MAGGRID: same as STAR-G-MAGGRID, but for galaxies (the magnitude grids STAR-G-MAGGRID and STAR-C1M861-MAGGRID are defined through parameters Sky_ObjectDensity_MagMin, Sky_ObjectDensity_MagMax, and Sky_ObjectDensity_MagBinWidth); (4) GAL-G-DENSITY: same as STAR-G-DENSITY, but for galaxies; (5) STAR-GS-MAGGRID: same as STAR-G-MAGGRID, but in GS band; (6) STAR-GS-DENSITY: same as STAR-G-DENSITY, but in GS band; (7) GAL-GS-MAGGRID: same as STAR-GS-MAGGRID, but for galaxies; (8) GAL-GS-DENSITY: same as STAR-GS-DENSITY, but for galaxies; (9) STAR-C1M861-MAGGRID: same as STAR-G-MAGGRID, but in G_RVS band; (10) STAR-C1M861-DENSITY: same as STAR-G-DENSITY, but in G_RVS band; (11) GAL-C1M861-MAGGRID: same as STAR-C1M861-MAGGRID, but for galaxies; and (12) GAL-C1M861-DENSITY: same as STAR-C1M861-DENSITY, but for galaxies
    // Source: The data used to produce the standard Gaia Galaxy model have been provided by the teams of the Observatoire de Besancon (A. Robin, C. Reyle, et al.) and the Observatorio Astronomico di Torino (R. Drimmel, et al.); the harmonisation has been provided by the SWG (X. Luri). See U. Lammers, 22 April 2005, 'Gaia Standard Galaxy Model Access Software (GSGMAS); User Guide', GAIA-UL-010, issue 1, revision 0. See also R. Drimmel, et al., 19 July 2005, 'Recommendations on the use of estimated star counts for Gaia studies', SWG-RD-004
    // Basic : true
    // Scalar: false
    static const char  *const SKY_OBJECTDENSITY() { return "Nature/Sky_ObjectDensity_003.fits"; }

    // Magnitude step (bin size) of the magnitude grid (STAR-G-MAGGRID and STAR-C1M861-MAGGRID) underlying the object counts in parameter Sky_ObjectDensity
    // Source: Extracted from parameter Sky_ObjectDensity_003
    // Basic : true
    // Scalar: true
    static const double SKY_OBJECTDENSITY_MAGBINWIDTH = 0.5; // [mag]

    // Maximum magnitude of the magnitude grid (STAR-G-MAGGRID and STAR-C1M861-MAGGRID) underlying the object counts in parameter Sky_ObjectDensity
    // Source: Extracted from parameter Sky_ObjectDensity_003
    // Basic : true
    // Scalar: true
    static const double SKY_OBJECTDENSITY_MAGMAX = 21.0; // [mag]

    // Minimum magnitude of the magnitude grid (STAR-G-MAGGRID and STAR-C1M861-MAGGRID) underlying the object counts in parameter Sky_ObjectDensity
    // Source: Extracted from parameter Sky_ObjectDensity_003
    // Basic : true
    // Scalar: true
    static const double SKY_OBJECTDENSITY_MAGMIN = 4.0; // [mag]

    // Sky-averaged density of stars with G <= 20.00 mag. Note that a value of 25.1E3 stars deg^-2 is given in ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 6.4.5 and Table 6.6, pages 240-242
    // Source: See also ESA, 21 May 2013, 'Gaia mission requirements document (MRD)', GAIA-EST-RD-00553, issue 3, revision 1, Section 4.1
    // Basic : false
    // Scalar: true
    static const double SKY_STARDENSITY_AVERAGE = 25089.; // [stars deg^-2]

    // Density of stars with G <= 20.00 mag below which Gaia is designed to operate nominally: nominal observations and all astrometric requirements shall be achieved in the two superimposed fields of view computed using the design object density in one instrument field of view plus the typical object density in the other instrument field of view (Requirement SCI-150). Requirement SCI-160 guarantees that a strategy to observe high-density sky regions (e.g., Baade's window), with stellar densities up to the worst-case star density, will be implemented (if higher densities than the worst-case density are encountered, the brightest objects up to the worst-case density will be observed). One option is to have several successive transits of the same sky region at a reduced precession rate using a modified scanning law (MSL; see L. Lindegren, 7 April 2005, 'Multi-pass scanning across Baade`s Window', GAIA-LL-058, issue 1, revision 0 and L. Lindegren, 22 August 2010, 'Reference Scanning Law for Gaia', GAIA-C3-TN-LU-LL-085-01)
    // Source: ESA, 21 May 2013, 'Gaia mission requirements document (MRD)', GAIA-EST-RD-00553, issue 3, revision 1, Section 4.1, Requirements SCI-150 and SCI-160
    // Basic : true
    // Scalar: true
    static const double SKY_STARDENSITY_DESIGN = 600000.; // [stars deg^-2]

    // Density of stars with G <= 20.00 mag in the Galactic-latitude ranges 0-5 degrees, 5-10 degrees, 10-20 degrees, 20-30 degrees, and 30-90 degrees
    // Source: ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 6.4.5 and Table 6.6, pages 240-242
    // Basic : true
    // Scalar: true
    static const double  *const SKY_STARDENSITY_LATITUDE ()  { static double _v[5] = { 101.6,  79.8,  31.2,  11.4,  3.8 }; return _v; } // [1E3 stars deg^-2]

    // The brighter star-selection limiting magnitude of RVS compared to SM/AF/BP/RP (G_RVS (= C1M861 = RVF) = 17.00 mag versus G = 20.00 mag, respectively) is assumed to correspond to a factor 6 reduction in star density and the number of stars
    // Source: D. Katz, M. Cropper, J.-M. Desert, et al., 3 November 2006, 'CU6 - Spectroscopic processing - preliminary functional and data-flow analysis', GAIA-C6-TN-OPM-DK-001, issue 3, revision 0, Section 4
    // Basic : false
    // Scalar: true
    static const double SKY_STARDENSITY_RVSREDUCTIONFACTOR = 6.0;

    // Typical density of stars with G <= 20.00 mag. This value is a typical value encountered in the Galactic plane and, as such, is more representative than the sky-averaged value
    // Source: J.H.J. de Bruijne, 31 October 2003, 'PDHE load assumptions: properties of the sky', GAIA-JdB-009. See also ESA, 21 May 2013, 'Gaia mission requirements document (MRD)', GAIA-EST-RD-00553, issue 3, revision 1, Section 4.1 and ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 3.6.6, page 176
    // Basic : true
    // Scalar: true
    static const double SKY_STARDENSITY_TYPICAL = 150000.; // [stars deg^-2]

    // Density of stars with G <= 20.00 mag (worst-case, localised value, e.g., in Baade's window)
    // Source: ESA, 21 May 2013, 'Gaia mission requirements document (MRD)', GAIA-EST-RD-00553, issue 3, revision 1, Section 4.1. See also ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 6.5.2, page 244
    // Basic : true
    // Scalar: true
    static const double SKY_STARDENSITY_WORSTCASE = 3.0e+06; // [stars deg^-2]

    // Johnson V band sky-background surface brightness seen from space; average value
    // Source: ESA, 21 May 2013, 'Gaia mission requirements document (MRD)', GAIA-EST-RD-00553, issue 3, revision 1, Requirement SCI-090. See also ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 6.4.3, page 239
    // Basic : true
    // Scalar: true
    static const double SKY_SURFACEBRIGHTNESS_AVERAGE = 22.5; // [mag arcsec^-2]

    // Johnson V band sky-background surface brightness seen from space; worst-case value (in the ecliptic)
    // Source: ESA, July 2000, 'Gaia; Composition, Formation and Evolution of the Galaxy; Concept and Technology Study Report', ESA-SCI(2000)4, Section 6.4.3, page 239
    // Basic : true
    // Scalar: true
    static const double SKY_SURFACEBRIGHTNESS_WORSTCASE = 21.0; // [mag arcsec^-2]

    // Catalogue containing CCD images, in units of photo-electrons, of typical solar-proton events for an AF CCD (used in BAM, WFS, SM, and AF; this CCD is also used in BP albeit with a different anti-reflection coating). Normally, outside periods of solar activity (solar flares), the solar-proton flux will be negligibly small. During solar activity (solar flares), solar-proton fluxes vary strongly with time, reaching peak levels of 1E6 particles cm^-2 s^-1 or higher. Note that cosmic rays and solar protons are distinct particles, collectively refered to as prompt-particle events (PPEs). The catalogue contains 12954 events. The structure of the FITS file is as follows: the first FITS-file extension contains a list of events containing event number, number of pixels across-scan in the image, and number of pixels along-scan in the image. The following extensions contain the individual images ('pixel matrices'), in units of photo-electron counts, one image per extension
    // Source: A. Short (ESA), priv. comm., 12 May 2006
    // Basic : true
    // Scalar: false
    static const char  *const SOLARPROTON_CATALOGUE_AFCCD() { return "Nature/SolarProton_Catalogue_AFCCD_001.fits"; } // [e^-]

    // Catalogue containing CCD images, in units of photo-electrons, of typical solar-proton events for a red-enhanced CCD (used in RP and RVS). Normally, outside periods of solar activity (solar flares), the solar-proton flux will be negligibly small. During solar activity (solar flares), solar-proton fluxes vary strongly with time, reaching peak levels of 1E6 particles cm^-2 s^-1 or higher. Note that cosmic rays and solar protons are distinct particles, collectively refered to as prompt-particle events (PPEs). The catalogue contains 4950 events. The structure of the FITS file is as follows: the first FITS-file extension contains a list of events containing event number, number of pixels across-scan in the image, and number of pixels along-scan in the image. The following extensions contain the individual images ('pixel matrices'), in units of photo-electron counts, one image per extension
    // Source: A. Short (ESA), priv. comm., 1 September 2008
    // Basic : true
    // Scalar: false
    static const char  *const SOLARPROTON_CATALOGUE_REDENHANCEDCCD() { return "Nature/SolarProton_Catalogue_RedEnhancedCCD_001.fits"; } // [e^-]

    // Maximum-sustainable solar-proton flux at L2, in units of particles cm^-2 s^-1. Normally, outside periods of solar activity (solar flares), the solar-proton flux will be negligibly small. During solar activity (solar flares), solar-proton fluxes vary strongly with time, reaching peak levels of 1E6 particles cm^-2 s^-1 or higher. Note that cosmic rays and solar protons are distinct particles, collectively refered to as prompt-particle events (PPEs). A PPE rate of 1300 particles cm^-2 s^-1 is assumed to correspond to the operational limit below which Gaia functions nominally and above which the spacecraft will be in AOCS TranSition Mode (TSM). An isotropic prompt-particle event flux N, in units of events cm^-2 s^-1, generates 2 N A / 4 events s^-1 CCD^-1, where A denotes the active-pixel area of the CCD in units of cm^2 (including any reduction as a result of TDI-gating), the factor 2 results from considering 'inflow' through both the illuminated and the non-illuminated faces of the CCD, and the factor 4 results from the 'flat geometry' of the CCD (see J.H.J. de Bruijne, A. Short, 7 September 2005, 'prompt-particle events: from fluxes to count rates', GAIA-JdB-026, issue 1, revision 0)
    // Source: EADS-Astrium, 3 March 2011, 'Dead-time budget', GAIA.ASF.TCN.SAT.00133, issue 5, revision 1
    // Basic : true
    // Scalar: true
    static const double SOLARPROTON_FLUX_L2 = 1300.; // [particles cm^-2 s^-1]

    // Space sink temperature at L2
    // Source: European Cooperation for Space Standards (ECSS), 15 November 2008, 'Space environment standard', ECSS-E-ST-10-04C, Section 6.2.1 (http://www.ecss.nl/)
    // Basic : true
    // Scalar: true
    static const double SPACE_TEMPERATURE_L2 = 3.; // [K]

    // Specific gas constant for dry air, defined as the molar gas constant (:Nature:MolarGas_Constant) divided by the molar mass (which is 0.0289644 kg mol^-1 for the International Standard Atmopshere)
    // Source: F. Kleijer (Netherlands Geodetic Commission, Delft), 1 April 2004, 'Troposphere Modeling and Filtering for Precise GPS Leveling', Publications on Geodesy 56, ISBN 90 6132 284 7 (http://www.ncg.knaw.nl/Publicaties/Geodesy/pdf/56Kleijer.pdf and http://repository.tudelft.nl/view/ir/uuid%3Aea1f0cf0-4e48-421b-b7ae-4ae3e36d1880/)
    // Basic : true
    // Scalar: true
    static const double SPECIFICGAS_CONSTANT_DRYAIR = 287.060; // [m^2 s^-2 K^-1]

    // Specific gas constant for water vapour, defined as the molar gas constant (:Nature:MolarGas_Constant) divided by the molar mass (which is 0.0180152 kg mol^-1)
    // Source: F. Kleijer (Netherlands Geodetic Commission, Delft), 1 April 2004, 'Troposphere Modeling and Filtering for Precise GPS Leveling', Publications on Geodesy 56, ISBN 90 6132 284 7 (http://www.ncg.knaw.nl/Publicaties/Geodesy/pdf/56Kleijer.pdf and http://repository.tudelft.nl/view/ir/uuid%3Aea1f0cf0-4e48-421b-b7ae-4ae3e36d1880/)
    // Basic : true
    // Scalar: true
    static const double SPECIFICGAS_CONSTANT_WATERVAPOUR = 461.525; // [m^2 s^-2 K^-1]

    // Stefan-Boltzmann constant. Note: best-measured value equals 5.670367E-8 (P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0))
    // Basic : false
    // Scalar: true
    static const double STEFANBOLTZMANN_CONSTANT = 5.6703666e-08; // [W m^-2 K^-4]

    // One steradian in units of square degrees
    // Basic : false
    // Scalar: true
    static const double STERADIAN_DEGREESQUARE = 3282.80635001174; // [deg^2]

    // Ratio of Sun to Earth-system mass. The planetary mass includes the contribution of its satellite, the Moon
    // Basic : false
    // Scalar: true
    static const double SUNTOEARTHSYSTEM_MASSRATIO = 328900.559616;

    // Ratio of Sun to Earth mass. The Earth mass includes the Earth's atmosphere but excludes the mass of the Moon
    // Basic : false
    // Scalar: true
    static const double SUNTOEARTH_MASSRATIO = 332946.048701;

    // Ratio of Sun to Eris-system mass (IAU 2009 CBE value). The 'planetary' mass includes the contribution of its satellite, Dysnomia
    // Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double SUNTOERISSYSTEM_MASSRATIO = 1.1910e+08;

    // Ratio of Sun to Jupiter-system mass (IAU 2009 CBE value). The planetary mass includes the contribution of its satellites
    // Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double SUNTOJUPITERSYSTEM_MASSRATIO = 1.0473486440e+03;

    // Ratio of Sun to Mars-system mass (IAU 2009 CBE value). The planetary mass includes the contribution of its satellites
    // Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double SUNTOMARSSYSTEM_MASSRATIO = 3.098703590e+06;

    // Ratio of Sun to Mercury(-system) mass (IAU 2009 CBE value)
    // Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double SUNTOMERCURYSYSTEM_MASSRATIO = 6.02360e+06;

    // Ratio of Sun to Neptune-system mass (IAU 2009 CBE value). The planetary mass includes the contribution of its satellites
    // Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double SUNTONEPTUNESYSTEM_MASSRATIO = 1.9412260e+04;

    // Ratio of Sun to Pluto-system mass (IAU 2009 CBE value). The 'planetary' mass includes the contribution of its satellite, Charon
    // Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double SUNTOPLUTOSYSTEM_MASSRATIO = 1.365660e+08;

    // Ratio of Sun to Saturn-system mass (IAU 2009 CBE value). The planetary mass includes the contribution of its satellites
    // Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double SUNTOSATURNSYSTEM_MASSRATIO = 3.49790180e+03;

    // Ratio of Sun to Uranus-system mass (IAU 2009 CBE value). The planetary mass includes the contribution of its satellites
    // Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double SUNTOURANUSSYSTEM_MASSRATIO = 2.2902980e+04;

    // Ratio of Sun to Venus(-system) mass (IAU 2009 CBE value)
    // Source: The IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double SUNTOVENUSSYSTEM_MASSRATIO = 4.085237190e+05;

    // Johnson V absolute magnitude of the Sun
    // Basic : false
    // Scalar: true
    static const double SUN_ABSOLUTEVMAGNITUDE = 4.812; // [mag]

    // Mixing length parameter \alpha of the Sun
    // Source: L. Girardi, A. Bressan, G. Bertelli, C. Chiosi, 2000, 'Evolutionary tracks and isochrones for low- and intermediate-mass stars; from M = 0.15 to 7 M_Sun and from Z = 0.0004 to 0.03', A&AS, 141, 371
    // Basic : true
    // Scalar: true
    static const double SUN_ALPHA = 1.68;

    // Johnson B-V colour of the Sun
    // Source: Derived from :Nature:Planck_Constant, :Nature:VelocityOfLight_Constant_Vacuum, :Nature:A0VStar_CalibrationWavelength, :Nature:A0VStar_Spectrum_Nu_002, :Nature:A0VStar_CalibrationFunction_002, :Nature:Sun_Spectrum_Nu_001, :Nature:A0VStar_VMinI, :Nature:A0VStar_VMinR, :Nature:A0VStar_BMinV, :Nature:A0VStar_RMinI, :Nature:FilterTransmissionCurve_JohnsonCousinsB_002, :Nature:FilterTransmissionCurve_JohnsonCousinsV_002, :Nature:FilterTransmissionCurve_JohnsonCousinsR_002, and :Nature:FilterTransmissionCurve_JohnsonCousinsI_002. See also I. Ramirez, et al., 2012, 'The UBV(RI)C Colors of the Sun', Astrophysical Journal (ApJ), 752, 5, J. Holmberg, C. Flynn, L. Portinari, 2006, 'The colours of the Sun', MNRAS, 367, 449, and B.J. Taylor, 1998, 'The colours of the Sun', proceedings of IAU Symposium 189 on 'Fundamental Stellar Properties: The Interaction between Observation and Theory', eds T.R. Bedding, A.J. Booth, J. Davis, Kluwer, Dordrecht, ISBN 0-7923-4651-3, page 83 (1998IAUS..189...83T)
    // Basic : true
    // Scalar: true
    static const double SUN_BMINV = 0.678; // [mag]

    // Johnson V bolometric correction of the Sun
    // Basic : false
    // Scalar: true
    static const double SUN_BOLOMETRICCORRECTIONVMAGNITUDE = -0.072; // [mag]

    // Absolute bolometric magnitude of the Sun
    // Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B2 on Recommended Zero Points for the Absolute and Apparent Bolometric Magnitude Scales', arXiv:1510.06262 (see also https://www.iau.org/administration/resolutions/general_assemblies/)
    // Basic : false
    // Scalar: true
    static const double SUN_BOLOMETRICMAGNITUDE = 4.740; // [mag]

    // Solar diurnal parallax
    // Basic : false
    // Scalar: true
    static const double SUN_DIURNALPARALLAX = 8.794143; // [arcsec]

    // Effective (black-body) temperature of the Sun
    // Basic : false
    // Scalar: true
    static const double SUN_EFFECTIVETEMPERATURE = 5772.; // [K]

    // Nominal effective (black-body) temperature of the Sun, in SI units. This nominal value shall be understood as conversion factor only
    // Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)
    // Basic : true
    // Scalar: true
    static const double SUN_EFFECTIVETEMPERATURE_NOMINAL = 5772.; // [K]

    // Radius of the smallest hypothetical sphere around the Sun which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double SUN_ENCOMPASSINGSPHERERADIUS = 6.960e+08; // [m]

    // Energy flux of the Sun at a distance of 1 au (also refered to as solar constant or total solar irradiance). This is the cycle-23-averaged, measured value. Due to orbital modulation, this value changes by \pm 3.4% during a year; this value changes by \pm 0.1% during a solar cycle; during solar maximum, sunspots can induce variations over one solar rotation at the level of a few tenths of a percent
    // Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)
    // Basic : true
    // Scalar: true
    static const double SUN_ENERGYFLUX_ASTRONOMICALUNIT = 1361.; // [W m^-2]

    // Energy flux of the Sun at L2. Due to orbital modulation, this value changes by \pm 3.4% during a year; this value changes by \pm 0.1% during a solar cycle; during solar maximum, sunspots can induce variations over one solar rotation at the level of a few tenths of a percent
    // Basic : false
    // Scalar: true
    static const double SUN_ENERGYFLUX_L2 = 1334.; // [W m^-2]

    // Nominal energy flux of the Sun at a distance of 1 au (also refered to as solar constant or total solar irradiance), in SI units. This nominal value shall be understood as conversion factor only
    // Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)
    // Basic : true
    // Scalar: true
    static const double SUN_ENERGYFLUX_NOMINAL = 1361.; // [W m^-2]

    // Solar (equatorial) radius (photosphere)
    // Source: M. Haberreiter, W. Schmutz, A.G. Kosovichev, 2008, 'Solving the Discrepancy between the Seismic and Photospheric Solar Radius', Astrophysical Journal (ApJ), 675, L53
    // Basic : true
    // Scalar: true
    static const double SUN_EQUATORIALRADIUS = 6.956580e+08; // [m]

    // Solar apparent (equatorial) radius at unit distance
    // Basic : false
    // Scalar: true
    static const double SUN_EQUATORIALRADIUS_APPARENT = 959.17; // [arcsec]

    // Nominal solar (equatorial) radius (photosphere), in SI units. This nominal value shall be understood as conversion factor only
    // Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)
    // Basic : true
    // Scalar: true
    static const double SUN_EQUATORIALRADIUS_NOMINAL = 6.9570e+08; // [m]

    // Coarse estimate of the solar (mean) value of the mean free photon path (assuming complete ionisation)
    // Source: E.g., R. Kippenhahn, A. Weigert, 1991, 'Stellar structure and evolution' (corrected 2-nd printing), Springer Verlag, Berlin, Section 5, Equation 5.1, page 27
    // Basic : false
    // Scalar: true
    static const double SUN_FREEPHOTONPATH_MEAN = 0.021; // [m]

    // Heliocentric gravitational constant (TCB-compatible value). Note that IAU 2012 Resolution B2 adopted at the XXVIII-th General Assembly of the IAU recommends that this parameter be determined observationally in SI units
    // Basic : false
    // Scalar: true
    static const double SUN_GM = 1.3271244210789466e+20; // [m^3 s^-2]

    // Nominal heliocentric gravitational constant, in SI units. This nominal value shall be understood as conversion factor only
    // Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)
    // Basic : true
    // Scalar: true
    static const double SUN_GM_NOMINAL = 1.32712440e+20; // [m^3 s^-2]

    // Post-Newtonian deflection angle, for an observer at 1 au from the Sun, of a Solar-limb-grazing light ray due to the spherically symmetric part of the gravitational field of the Sun
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64
    // Basic : false
    // Scalar: true
    static const double SUN_LIGHTDEFLECTION_LIMB = 1751293.; // [10^-6 arcsec]

    // Post-Newtonian deflection angle, for an observer at 1 au from the Sun, of a light ray arriving at right angles to the solar direction due to the spherically symmetric part of the gravitational field of the Sun
    // Source: ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 3, page 331; cf. S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64
    // Basic : false
    // Scalar: true
    static const double SUN_LIGHTDEFLECTION_RIGHTANGLES = 4072.; // [10^-6 arcsec]

    // Luminosity of the Sun. This value changes by \pm 0.1% during a solar cycle; during solar maximum, sunspots can induce variations over one solar rotation at the level of a few tenths of a percent
    // Basic : false
    // Scalar: true
    static const double SUN_LUMINOSITY = 3.8275e+26; // [W]

    // Nominal luminosity of the Sun, in SI units. This nominal value shall be understood as conversion factor only
    // Source: E. Mamajek, et al., 2015, 'IAU 2015 Resolution B3 on Recommended Nominal Conversion Constants for Selected Solar and Planetary Properties', arXiv:1510.07674 (see also https://www.iau.org/administration/resolutions/general_assemblies/)
    // Basic : true
    // Scalar: true
    static const double SUN_LUMINOSITY_NOMINAL = 3.8280e+26; // [W]

    // Solar mass
    // Basic : false
    // Scalar: true
    static const double SUN_MASS = 1.98848e+30; // [kg]

    // Mean solar mass density
    // Basic : false
    // Scalar: true
    static const double SUN_MASSDENSITY_MEAN = 1.410; // [g cm^-3]

    // Solar value of the mean molecular weight (assuming complete ionisation)
    // Source: E.g., H. Karttunen, et al., 1987, 'Fundamental Astronomy', Springer Verlag, Berlin, Section 11.2, Equation 11.8, page 245
    // Basic : false
    // Scalar: true
    static const double SUN_MEANMOLECULARWEIGHT = 0.6092;

    // Solar value of the mean molecular weight per free electron (assuming complete ionisation)
    // Source: E.g., H. Karttunen, et al., 1987, 'Fundamental Astronomy', Springer Verlag, Berlin, Section 11.2, page 246
    // Basic : false
    // Scalar: true
    static const double SUN_MEANMOLECULARWEIGHT_PERFREEELECTRON = 1.1651;

    // IAU-recommended value for the declination \delta_0 of the north pole of rotation of the Sun. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value shall be used for comparative purposes only
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double SUN_NORTHROTATIONALPOLE_DECLINATION = 63.87; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the declination \delta_0 of the north pole of rotation of the Sun. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value shall be used for comparative purposes only
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double SUN_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE = 0.; // [deg day^-1]

    // IAU-recommended value for the right ascension \alpha_0 of the north pole of rotation of the Sun. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value shall be used for comparative purposes only
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double SUN_NORTHROTATIONALPOLE_RIGHTASCENSION = 286.13; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the right ascension \alpha_0 of the north pole of rotation of the Sun. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees. The numerical value shall be used for comparative purposes only
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double SUN_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE = 0.; // [deg day^-1]

    // Orbital semi-major axis of the Sun's orbit around the (hypothetical) common barycentre of the Earth-system-Sun system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.1, Equation 1.6, page 6
    // Basic : false
    // Scalar: true
    static const double SUN_ORBITALSEMIMAJORAXIS_EARTHSYSTEM = 455.; // [km]

    // Orbital semi-major axis of the Sun's orbit around the (hypothetical) common barycentre of the Jupiter-system-Sun system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.1, Equation 1.6, page 6
    // Basic : false
    // Scalar: true
    static const double SUN_ORBITALSEMIMAJORAXIS_JUPITERSYSTEM = 743154.; // [km]

    // Orbital semi-major axis of the Sun's orbit around the (hypothetical) common barycentre of the Mars-system-Sun system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.1, Equation 1.6, page 6
    // Basic : false
    // Scalar: true
    static const double SUN_ORBITALSEMIMAJORAXIS_MARSSYSTEM = 74.; // [km]

    // Orbital semi-major axis of the Sun's orbit around the (hypothetical) common barycentre of the Mercury-Sun system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.1, Equation 1.6, page 6
    // Basic : false
    // Scalar: true
    static const double SUN_ORBITALSEMIMAJORAXIS_MERCURYSYSTEM = 10.; // [km]

    // Orbital semi-major axis of the Sun's orbit around the (hypothetical) common barycentre of the Neptune-system-Sun system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.1, Equation 1.6, page 6
    // Basic : false
    // Scalar: true
    static const double SUN_ORBITALSEMIMAJORAXIS_NEPTUNESYSTEM = 231730.; // [km]

    // Orbital semi-major axis of the Sun's orbit around the (hypothetical) common barycentre of the Pluto-system-Sun system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.1, Equation 1.6, page 6
    // Basic : false
    // Scalar: true
    static const double SUN_ORBITALSEMIMAJORAXIS_PLUTOSYSTEM = 43.; // [km]

    // Orbital semi-major axis of the Sun's orbit around the (hypothetical) common barycentre of the Saturn-system-Sun system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.1, Equation 1.6, page 6
    // Basic : false
    // Scalar: true
    static const double SUN_ORBITALSEMIMAJORAXIS_SATURNSYSTEM = 407863.; // [km]

    // Orbital semi-major axis of the Sun's orbit around the (hypothetical) common barycentre of the Uranus-system-Sun system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.1, Equation 1.6, page 6
    // Basic : false
    // Scalar: true
    static const double SUN_ORBITALSEMIMAJORAXIS_URANUSSYSTEM = 125340.; // [km]

    // Orbital semi-major axis of the Sun's orbit around the (hypothetical) common barycentre of the Venus-Sun system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.1, Equation 1.6, page 6
    // Basic : false
    // Scalar: true
    static const double SUN_ORBITALSEMIMAJORAXIS_VENUSSYSTEM = 265.; // [km]

    // IAU-recommended value for the ephemeris position of the prime meridian of the Sun. The location of the prime meridian is specified by the angle that is measured along the Sun's equator in an easterly direction with respect to the Sun's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the Sun's equator on the standard equator to the point B where the prime meridian crosses the Sun's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the Sun's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the Sun, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the Sun. If W increases with time, the body has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde. The numerical value shall be used for comparative purposes only
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double SUN_PRIMEMERIDIAN_EPHEMERISPOSITION = 84.176; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch of 1.5 January 2000 = JD2451545.0 TDB) of the ephemeris position of the prime meridian of the Sun. The location of the prime meridian is specified by the angle that is measured along the Sun's equator in an easterly direction with respect to the Sun's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the Sun's equator on the standard equator to the point B where the prime meridian crosses the Sun's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the Sun's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the Sun, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the Sun. If W increases with time, the body has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde. The numerical value shall be used for comparative purposes only
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double SUN_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE = 14.1844000; // [deg day^-1]

    // Cousins R-I colour of the Sun
    // Source: Derived from :Nature:Planck_Constant, :Nature:VelocityOfLight_Constant_Vacuum, :Nature:A0VStar_CalibrationWavelength, :Nature:A0VStar_Spectrum_Nu_002, :Nature:A0VStar_CalibrationFunction_002, :Nature:Sun_Spectrum_Nu_001, :Nature:A0VStar_VMinI, :Nature:A0VStar_VMinR, :Nature:A0VStar_BMinV, :Nature:A0VStar_RMinI, :Nature:FilterTransmissionCurve_JohnsonCousinsB_002, :Nature:FilterTransmissionCurve_JohnsonCousinsV_002, :Nature:FilterTransmissionCurve_JohnsonCousinsR_002, and :Nature:FilterTransmissionCurve_JohnsonCousinsI_002. See also I. Ramirez, et al., 2012, 'The UBV(RI)C Colors of the Sun', Astrophysical Journal (ApJ), 752, 5, J. Holmberg, C. Flynn, L. Portinari, 2006, 'The colours of the Sun', MNRAS, 367, 449, and B.J. Taylor, 1998, 'The colours of the Sun', proceedings of IAU Symposium 189 on 'Fundamental Stellar Properties: The Interaction between Observation and Theory', eds T.R. Bedding, A.J. Booth, J. Davis, Kluwer, Dordrecht, ISBN 0-7923-4651-3, page 83 (1998IAUS..189...83T)
    // Basic : true
    // Scalar: true
    static const double SUN_RMINI = 0.344; // [mag]

    // Mean value of the solar-radiation pressure at a distance of 1 au. Due to orbital modulation, this value changes by \pm 3.4% during a year; this value changes by \pm 0.1% during a solar cycle; during solar maximum, sunspots can induce variations over one solar rotation at the level of a few tenths of a percent
    // Basic : false
    // Scalar: true
    static const double SUN_RADIATIONPRESSURE_ASTRONOMICALUNIT = 4.540e-06; // [Pa]

    // Mean value of the solar-radiation pressure at L2. Due to orbital modulation, this value changes by \pm 3.4% during a year; this value changes by \pm 0.1% during a solar cycle; during solar maximum, sunspots can induce variations over one solar rotation at the level of a few tenths of a percent
    // Basic : false
    // Scalar: true
    static const double SUN_RADIATIONPRESSURE_L2 = 4.450e-06; // [Pa]

    // Solar (mean) value of the Rosseland mean opacity for Thomson free-electron-scattering
    // Source: E.g., R. Kippenhahn, A. Weigert, 1991, 'Stellar structure and evolution' (corrected 2-nd printing), Springer Verlag, Berlin, Section 17, Equation 17.2, page 137
    // Basic : false
    // Scalar: true
    static const double SUN_ROSSELANDMEANOPACITY_THOMSONSCATTERING = 0.0344; // [m^2 kg^-1]

    // Spectrum f_{\odot\nu}(\lambda) of the Sun: Kurucz/ATLAS9 spectrum (CDROM 19). First column: wavelength \lambda (in nm; from 115.0 to 1062.0). Second column: Eddington flux (in W m^-2 Hz^-1 steradian^-1). Note that the flux at 115.0 nm was obtained using linear interpolation between the available fluxes at 114.5 and 115.5 nm (0.6745E-15 and 0.8131E-15, respectively). Note that the flux at 1062.0 nm was obtained using linear interpolation between the available fluxes at 1057.5 and 1062.5 nm (0.8892E-05 and 0.8861E-05, respectively)
    // Source: C. Jordi, priv. comm., 17 February 2003; see also http://gaia.am.ub.es/PWG/sun.mod. Note that data beyond the current wavelength limits are available
    // Basic : true
    // Scalar: false
    static const char  *const SUN_SPECTRUM_NU() { return "Nature/Sun_Spectrum_Nu_001.fits"; }

    // Surface gravity of the Sun
    // Basic : false
    // Scalar: true
    static const double SUN_SURFACEGRAVITY = 274.2; // [m s^-2]

    // Johnson V magnitude of the Sun (observed)
    // Source: M.S. Bessell, F. Castelli, B. Plez, 1998, 'Model atmospheres, broad-band colors, bolometric corrections, and temperature calibrations for O-M stars', A&A, 333, 231, Appendices C and D (erratum 1998, A&A, 337, 321)
    // Basic : true
    // Scalar: true
    static const double SUN_VMAGNITUDE = -26.760; // [mag]

    // Johnson-Cousins V-I colour of the Sun
    // Source: Derived from :Nature:Planck_Constant, :Nature:VelocityOfLight_Constant_Vacuum, :Nature:A0VStar_CalibrationWavelength, :Nature:A0VStar_Spectrum_Nu_002, :Nature:A0VStar_CalibrationFunction_002, :Nature:Sun_Spectrum_Nu_001, :Nature:A0VStar_VMinI, :Nature:A0VStar_VMinR, :Nature:A0VStar_BMinV, :Nature:A0VStar_RMinI, :Nature:FilterTransmissionCurve_JohnsonCousinsB_002, :Nature:FilterTransmissionCurve_JohnsonCousinsV_002, :Nature:FilterTransmissionCurve_JohnsonCousinsR_002, and :Nature:FilterTransmissionCurve_JohnsonCousinsI_002. See also I. Ramirez, et al., 2012, 'The UBV(RI)C Colors of the Sun', Astrophysical Journal (ApJ), 752, 5, J. Holmberg, C. Flynn, L. Portinari, 2006, 'The colours of the Sun', MNRAS, 367, 449, and B.J. Taylor, 1998, 'The colours of the Sun', proceedings of IAU Symposium 189 on 'Fundamental Stellar Properties: The Interaction between Observation and Theory', eds T.R. Bedding, A.J. Booth, J. Davis, Kluwer, Dordrecht, ISBN 0-7923-4651-3, page 83 (1998IAUS..189...83T)
    // Basic : true
    // Scalar: true
    static const double SUN_VMINI = 0.711; // [mag]

    // Johnson-Cousins V-R colour of the Sun
    // Source: Derived from :Nature:Planck_Constant, :Nature:VelocityOfLight_Constant_Vacuum, :Nature:A0VStar_CalibrationWavelength, :Nature:A0VStar_Spectrum_Nu_002, :Nature:A0VStar_CalibrationFunction_002, :Nature:Sun_Spectrum_Nu_001, :Nature:A0VStar_VMinI, :Nature:A0VStar_VMinR, :Nature:A0VStar_BMinV, :Nature:A0VStar_RMinI, :Nature:FilterTransmissionCurve_JohnsonCousinsB_002, :Nature:FilterTransmissionCurve_JohnsonCousinsV_002, :Nature:FilterTransmissionCurve_JohnsonCousinsR_002, and :Nature:FilterTransmissionCurve_JohnsonCousinsI_002. See also I. Ramirez, et al., 2012, 'The UBV(RI)C Colors of the Sun', Astrophysical Journal (ApJ), 752, 5, J. Holmberg, C. Flynn, L. Portinari, 2006, 'The colours of the Sun', MNRAS, 367, 449, and B.J. Taylor, 1998, 'The colours of the Sun', proceedings of IAU Symposium 189 on 'Fundamental Stellar Properties: The Interaction between Observation and Theory', eds T.R. Bedding, A.J. Booth, J. Davis, Kluwer, Dordrecht, ISBN 0-7923-4651-3, page 83 (1998IAUS..189...83T)
    // Basic : false
    // Scalar: true
    static const double SUN_VMINR = 0.367; // [mag]

    // Velocity of the Sun with respect to the local standard of rest (LSR)
    // Basic : false
    // Scalar: true
    static const double SUN_VELOCITYLSR = 18.04; // [km s^-1]

    // Velocity of the Sun with respect to the local standard of rest (LSR); U-component, i.e., radially inwards (towards the Galactic centre). The random error is +0.69 and -0.75 km s^-1; the systematic error is 1.0 km s^-1
    // Source: R. Schoenrich, J.J. Binney, W. Dehnen, 1 April 2010, 'Local kinematics and the local standard of rest', MNRAS, 403, 1829-1833 (2010MNRAS.403.1829S)
    // Basic : true
    // Scalar: true
    static const double SUN_VELOCITYLSR_U = 11.10; // [km s^-1]

    // Velocity of the Sun with respect to the local standard of rest (LSR); V-component, i.e., in the direction of Galactic rotation. The random error is +0.47 and -0.47 km s^-1; the systematic error is 2.0 km s^-1
    // Source: R. Schoenrich, J.J. Binney, W. Dehnen, 1 April 2010, 'Local kinematics and the local standard of rest', MNRAS, 403, 1829-1833 (2010MNRAS.403.1829S)
    // Basic : true
    // Scalar: true
    static const double SUN_VELOCITYLSR_V = 12.24; // [km s^-1]

    // Velocity of the Sun with respect to the local standard of rest (LSR); W-component, i.e., vertically upwards (towards the north Galactic pole). The random error is +0.37 and -0.36 km s^-1; the systematic error is 0.5 km s^-1
    // Source: R. Schoenrich, J.J. Binney, W. Dehnen, 1 April 2010, 'Local kinematics and the local standard of rest', MNRAS, 403, 1829-1833 (2010MNRAS.403.1829S)
    // Basic : true
    // Scalar: true
    static const double SUN_VELOCITYLSR_W = 7.25; // [km s^-1]

    // Hydrogen abundance of the Sun by mass
    // Basic : false
    // Scalar: true
    static const double SUN_X = 0.7166;

    // Helium abundance of the Sun by mass
    // Source: N. Grevesse, A. Noels, 1993, in 'Association Vaudoise des chercheurs en physique: la formation des elements chimiques', eds B. Hauck, S. Plantani, D. Raboud and N. Grevesse, A. Noels, 1993, in 'Origin and evolution of the elements', eds N. Prantzos, E. Vangioni-Flam, M. Casse, Cambridge University Press, Cambridge, page 15
    // Basic : true
    // Scalar: true
    static const double SUN_Y = 0.2659;

    // Metal abundance of the Sun by mass
    // Source: N. Grevesse, A. Noels, 1993, in 'Association Vaudoise des chercheurs en physique: la formation des elements chimiques', eds B. Hauck, S. Plantani, D. Raboud and N. Grevesse, A. Noels, 1993, in 'Origin and evolution of the elements', eds N. Prantzos, E. Vangioni-Flam, M. Casse, Cambridge University Press, Cambridge, page 15
    // Basic : true
    // Scalar: true
    static const double SUN_Z = 0.0175;

    // The difference of TAI and UTC as function of time, since 1972, in units of seconds. UTC differs from TAI by an integral number of seconds ('leap seconds'), in such a way that UT1-UTC stays smaller than 0.9 s in absolute value. The decision to introduce a leap second in UTC to meet this condition is the responsability of the IERS. Announcements are made through IERS Bulletin C (https://hpiers.obspm.fr/eoppc/bul/bulc/). Positive leap seconds make the difference TAI-UTC more positive. The vector elements of this parameter define triplets {STARTTIME, ENDTIME, LEAPSECOND} where STARTTIME denotes the start time of the validity interval (in JD UTC), ENDTIME denotes the end time of the validity interval (in JD UTC), and LEAPSECOND denotes the difference TAI - UTC in units of seconds. Note that the ENDTIME of the final triplet (JD2500000.5 UTC) indicates that the end time of the current validity interval is undefined
    // Source: IERS Bulletin C (https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat; see also http://maia.usno.navy.mil/ser7/tai-utc.dat)
    // Basic : true
    // Scalar: true
    static const double  (  *const TAIMINUTC_CONSTANT ()  )[3] { static double _v[28][3] = { {2441317.5,  2441499.5,  10.0},  {2441499.5,  2441683.5,  11.0},  {2441683.5,  2442048.5,  12.0},  {2442048.5,  2442413.5,  13.0},  {2442413.5,  2442778.5,  14.0},  {2442778.5,  2443144.5,  15.0},  {2443144.5,  2443509.5,  16.0},  {2443509.5,  2443874.5,  17.0},  {2443874.5,  2444239.5,  18.0},  {2444239.5,  2444786.5,  19.0},  {2444786.5,  2445151.5,  20.0},  {2445151.5,  2445516.5,  21.0},  {2445516.5,  2446247.5,  22.0},  {2446247.5,  2447161.5,  23.0},  {2447161.5,  2447892.5,  24.0},  {2447892.5,  2448257.5,  25.0},  {2448257.5,  2448804.5,  26.0},  {2448804.5,  2449169.5,  27.0},  {2449169.5,  2449534.5,  28.0},  {2449534.5,  2450083.5,  29.0},  {2450083.5,  2450630.5,  30.0},  {2450630.5,  2451179.5,  31.0},  {2451179.5,  2453736.5,  32.0},  {2453736.5,  2454832.5,  33.0},  {2454832.5,  2456109.5,  34.0},  {2456109.5,  2457204.5,  35.0},  {2457204.5,  2457754.5,  36.0},  {2457754.5,  2500000.5,  37.0} }; return _v; } // [s]

    // Nominal value of the (constant) offset between TAI and TT: TT(TAI) = TAI + 32.184 s. This offset was chosen to give continuity with the previously used, but now obsolete, Ephemeris Time (see P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7). Deviations of this value, which are attributable to physical defects of atomic time standards, are probably between the limits \pm 10 \mus
    // Source: ESA, 1997, 'The Hipparcos and Tycho Catalogues', Volume 1, page 23
    // Basic : true
    // Scalar: true
    static const double TTMINTAI_CONSTANT_NOMINAL = 32.184; // [s]

    // Temperature constant encountered in physics of stellar atmospheres
    // Source: E.g., E. Bohm-Vitense, 1997, 'Introduction to stellar astrophysics; Volume 2: stellar atmospheres', Cambridge University Press, page 74
    // Basic : false
    // Scalar: true
    static const double TEMPERATURE_CONSTANT = 5039.8; // [K]

    // Transformation matrix which transforms the unit-direction vector r_ecl, expressed in ecliptic coordinates, into the unit-direction vector r_equ in equatorial coordinates (ICRS): r_equ = A_K times r_ecl (see also ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, Section 1.5.3, inverse of Equation 1.5.12). Note that the ICRS origin is shifted in the equatorial plane from \Gamma by \phi = 0.05542 arcsec, positive from \Gamma to the ICRS origin (see J. Chapront, M. Chapront-Touze, G. Francou, 2002, 'A new determination of lunar orbital parameters, precession constant, and tidal acceleration from LLR measurements', A&A, 387, 700). The ICRS has an unambiguous definition with an origin in the ICRF equator defined by the realisation of the ICRF. The ecliptic system is less well-defined, potentially depending on additional conventions in dynamical theories. The transformation quantified here corresponds to the inertial mean ecliptic with obliquity (see parameter :Nature:ObliquityOfEcliptic_J2000) and \Gamma defined by reference to the ICRS equator (other possibilities include the mean equator for J2000 or one of the JPL ephemerides equators). Both the obliquity and the position of \Gamma on the ICRS equator with respect to the ICRS origin have been obtained from LLR measurements. The transformation quantified here has no time dependence (there is no secular variation of the obliquity and no precession): it simply defines the relative situation of the various planes at J2000.0
    // Source: ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, Section 1.5.3, Equation 1.5.7 (generalised matrix A_K)
    // Basic : true
    // Scalar: true
    static const double  *const TRANSFORMATIONMATRIX_ECLIPTICTOICRS ()  { static double _v[9] = { 0.9999999999999639,  2.465125329E-7,  -1.068762105E-7,  -2.686837421E-7,  0.9174821334228226,  -0.3977769913529863,  0E-16,  0.3977769913530006,  0.9174821334228557 }; return _v; }

    // Transformation matrix which transforms the unit-direction vector r_gal, expressed in Galactic coordinates, into the unit-direction vector r_equ in equatorial coordinates (ICRS): r_equ = A_G times r_gal (see ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, Section 1.5.3, inverse of Equation 1.5.13)
    // Source: ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, Section 1.5.3, Equation 1.5.11 (matrix A_G)
    // Basic : true
    // Scalar: true
    static const double  *const TRANSFORMATIONMATRIX_GALACTICTOICRS ()  { static double _v[9] = { -0.0548755604162154,  0.4941094278755837,  -0.8676661490190047,  -0.8734370902348850,  -0.4448296299600112,  -0.1980763734312015,  -0.4838350155487132,  0.7469822444972189,  0.4559837761750669 }; return _v; }

    // Transformation matrix which transforms the unit-direction vector r_equ, expressed in equatorial coordinates (ICRS), into the unit-direction vector r_ecl in ecliptic coordinates: r_ecl = A_K^T times r_equ (see also ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, Section 1.5.3, Equation 1.5.12; A_K^T denotes the transpose of matrix A_K). Note that the ICRS origin is shifted in the equatorial plane from \Gamma by \phi = 0.05542 arcsec, positive from \Gamma to the ICRS origin (see J. Chapront, M. Chapront-Touze, G. Francou, 2002, 'A new determination of lunar orbital parameters, precession constant, and tidal acceleration from LLR measurements', A&A, 387, 700). The ICRS has an unambiguous definition with an origin in the ICRF equator defined by the realisation of the ICRF. The ecliptic system is less well-defined, potentially depending on additional conventions in dynamical theories. The transformation quantified here corresponds to the inertial mean ecliptic with obliquity (see parameter :Nature:ObliquityOfEcliptic_J2000) and \Gamma defined by reference to the ICRS equator (other possibilities include the mean equator for J2000 or one of the JPL ephemerides equators). Both the obliquity and the position of \Gamma on the ICRS equator with respect to the ICRS origin have been obtained from LLR measurements. The transformation quantified here has no time dependence (there is no secular variation of the obliquity and no precession): it simply defines the relative situation of the various planes at J2000.0
    // Source: ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, Section 1.5.3, transpose of Equation 1.5.7 (transpose of generalised matrix A_K)
    // Basic : true
    // Scalar: true
    static const double  *const TRANSFORMATIONMATRIX_ICRSTOECLIPTIC ()  { static double _v[9] = { 0.9999999999999639,  -2.686837421E-7,  0E-16,  2.465125329E-7,  0.9174821334228226,  0.3977769913530006,  -1.068762105E-7,  -0.3977769913529863,  0.9174821334228557 }; return _v; }

    // Transformation matrix which transforms the unit-direction vector r_equ, expressed in equatorial coordinates (ICRS), into the unit-direction vector r_gal in Galactic coordinates: r_gal = A_G^T times r_equ (see ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, Section 1.5.3, Equation 1.5.13; A_G^T denotes the transpose of matrix A_G)
    // Source: ESA, 1997, 'The Hipparcos and Tycho Catalogues', ESA SP-1200, Volume 1, Section 1.5.3, transpose of Equation 1.5.11 (transpose of matrix A_G)
    // Basic : true
    // Scalar: true
    static const double  *const TRANSFORMATIONMATRIX_ICRSTOGALACTIC ()  { static double _v[9] = { -0.0548755604162154,  -0.8734370902348850,  -0.4838350155487132,  0.4941094278755837,  -0.4448296299600112,  0.7469822444972189,  -0.8676661490190047,  -0.1980763734312015,  0.4559837761750669 }; return _v; }

    // Number of days per tropical year
    // Source: J.L. Simon, P. Bretagnon, J. Chapront, M. Chapront-Touze, G. Francou, J. Laskar, 1994, 'Numerical expressions for precession formulae and mean elements for the Moon and the planets', A&A, 282, 663 (1994A&A...282..663S)
    // Basic : true
    // Scalar: true
    static const double TROPICALYEAR_J2000DAY = 365.242190402; // [day]

    // Astrometric signature of the Sun induced by the Uranus system for an observer located at a distance of 10 pc from the Sun
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.7, Equation 1.22, page 11
    // Basic : false
    // Scalar: true
    static const double URANUSSYSTEM_ASTROMETRICSIGNATURE_10PARSEC = 84.; // [10^-6 arcsec]

    // Uranus-system mass (IAU 2009 CBE value). The planetary mass includes the contribution of its satellites
    // Basic : false
    // Scalar: true
    static const double URANUSSYSTEM_MASS = 8.68217e+25; // [kg]

    // Mean orbital eccentricity of Uranus, at the standard epoch J2000.0. The mean orbital eccentricity is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double URANUSSYSTEM_ORBITALECCENTRICITY_J2000 = 0.04725744;

    // Mean orbital inclination of Uranus, at the standard epoch J2000.0. The mean orbital inclination is is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double URANUSSYSTEM_ORBITALINCLINATION_J2000 = 0.77263783; // [deg]

    // Sidereal orbital period
    // Source: Values derived from the mean longitude rates in Table 5.8.1 in P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, page 316; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double URANUSSYSTEM_ORBITALPERIOD = 84.016846; // [yr]

    // Mean orbital semi-major axis of Uranus, at the standard epoch J2000.0. The mean orbital semi-major axis is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double URANUSSYSTEM_ORBITALSEMIMAJORAXIS_J2000 = 19.18916464; // [au]

    // Radial-velocity amplitude of the Sun induced by the Uranus system for 'an edge-on observer' (i.e., an observer in the orbital plane of the Uranus system)
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.4, Equation 1.18, page 9
    // Basic : false
    // Scalar: true
    static const double URANUSSYSTEM_RADIALVELOCITYSIGNATURE = 0.3; // [m s^-1]

    // Radius of the smallest hypothetical sphere around Uranus which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double URANUS_ENCOMPASSINGSPHERERADIUS = 2.55588e+07; // [m]

    // Equatorial radius of Uranus
    // Basic : false
    // Scalar: true
    static const double URANUS_EQUATORIALRADIUS = 2.55588e+07; // [m]

    // Geometrical flattening factor f of Uranus (f = (a-b)/a)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double URANUS_FLATTENING = 2.292730e-02;

    // Maximum reduction of the solar flux for an observer external to the solar system during a transit of Uranus
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.2, Equation 2.4, page 14
    // Basic : false
    // Scalar: true
    static const double URANUS_FLUXREDUCTION_MAXIMUM = 0.133; // [%]

    // Geometric albedo of Uranus (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double URANUS_GEOMETRICALBEDO = 0.51;

    // Dynamical form-factor of Uranus (oblateness or Stokes' second-degree zonal harmonic of the potential)
    // Source: P.R. Weissman, L.-A. McFadden, T.V. Johnson (eds.), 1999, 'Encyclopedia of the Solar System (first edition)', Academic Press, page 342
    // Basic : true
    // Scalar: true
    static const double URANUS_JSUB2 = 0.003516;

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of Uranus
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double URANUS_LIGHTDEFLECTION_LIMB = 2097.; // [10^-6 arcsec]

    // Mass of Uranus (do not use for high-precision (orbit) calculations)
    // Source: R.A. Jacobson, 2007, 'The gravity field of the Uranian system and the orbits of the Uranian satellites and rings', BAAS, 39, 453; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double URANUS_MASS = 8.681030e+25; // [kg]

    // Mean mass density of Uranus
    // Basic : false
    // Scalar: true
    static const double URANUS_MASSDENSITY_MEAN = 1.270; // [g cm^-3]

    // IAU-recommended value for the declination \delta_0 of the north pole of rotation of Uranus. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double URANUS_NORTHROTATIONALPOLE_DECLINATION = -15.175; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the declination \delta_0 of the north pole of rotation of Uranus. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double URANUS_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE = 0.; // [deg day^-1]

    // IAU-recommended value for the right ascension \alpha_0 of the north pole of rotation of Uranus. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double URANUS_NORTHROTATIONALPOLE_RIGHTASCENSION = 257.311; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the right ascension \alpha_0 of the north pole of rotation of Uranus. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double URANUS_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE = 0.; // [deg day^-1]

    // Polar radius of Uranus
    // Basic : false
    // Scalar: true
    static const double URANUS_POLARRADIUS = 2.49728e+07; // [m]

    // IAU-recommended value for the ephemeris position of the prime meridian of Uranus. The prime meridian refers to the rotation of the magnetic field System III. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double URANUS_PRIMEMERIDIAN_EPHEMERISPOSITION = 203.81; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch of 1.5 January 2000 = JD2451545.0 TDB) of the ephemeris position of the prime meridian of Uranus. The prime meridian refers to the rotation of the magnetic field System III. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double URANUS_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE = -501.1600928; // [deg day^-1]

    // Geometric transit probability (Uranus transiting the Sun) for an observer external to the solar system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.1, Equation 2.2, page 14
    // Basic : false
    // Scalar: true
    static const double URANUS_TRANSITPROBABILITY = 0.025; // [%]

    // Maximum transit time of Uranus (transiting the Sun) for an observer external to the solar system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.3, Equation 2.5, page 15
    // Basic : false
    // Scalar: true
    static const double URANUS_TRANSITTIME_MAXIMUM = 2.45; // [day]

    // V(1,0) magnitude of Uranus (i.e., the visual magnitude of the planet reduced to a distance of 1 au from both the Sun and Earth and phase angle zero). This parameter is also refered to as absolute magnitude in planetary sciences
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double URANUS_VONEZEROMAGNITUDE = -7.19; // [mag]

    // Mean volumetric radius of Uranus
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double URANUS_VOLUMETRICRADIUS = 2.53620e+07; // [m]

    // Velocity of light in vacuum (defining constant)
    // Source: P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0). See also the IAU (2009) System of Astronomical Constants (IAU, August 2009, 'IAU 2009 Astronomical Constants', IAU 2009 Resolution B2 adopted at the XXVII-th General Assembly of the IAU. See also IAU, 10 August 2009, 'IAU WG on NSFA Current Best Estimates', http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
    // Basic : true
    // Scalar: true
    static const double VELOCITYOFLIGHT_CONSTANT_VACUUM = 299792458.; // [m s^-1]

    // Astrometric signature of the Sun induced by Venus for an observer located at a distance of 10 pc from the Sun
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.7, Equation 1.22, page 11
    // Basic : false
    // Scalar: true
    static const double VENUSSYSTEM_ASTROMETRICSIGNATURE_10PARSEC = 0.177; // [10^-6 arcsec]

    // Venus(-system) mass (IAU 2009 CBE value)
    // Basic : false
    // Scalar: true
    static const double VENUSSYSTEM_MASS = 4.86747e+24; // [kg]

    // Mean orbital eccentricity of Venus, at the standard epoch J2000.0. The mean orbital eccentricity is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double VENUSSYSTEM_ORBITALECCENTRICITY_J2000 = 0.00677672;

    // Mean orbital inclination of Venus, at the standard epoch J2000.0. The mean orbital inclination is is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double VENUSSYSTEM_ORBITALINCLINATION_J2000 = 3.39467605; // [deg]

    // Sidereal orbital period
    // Source: Values derived from the mean longitude rates in Table 5.8.1 in P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, page 316; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double VENUSSYSTEM_ORBITALPERIOD = 0.61519726; // [yr]

    // Mean orbital semi-major axis of Venus, at the standard epoch J2000.0. The mean orbital semi-major axis is associated with a mean orbit solution from a 250-year least-squares fit of JPL's DE405 ephemeris to a Keplerian orbit covering the interval of years 1800 AD - 2050 AD. The orbital element is referenced to the mean ecliptic and equinox of J2000 at the standard epoch J2000.0; linear rates and fit errors can be found in the reference document
    // Source: E.M. Standish, 16 February 2006, 'Keplerian Elements for Approximate Positions of the Major Planets', http://ssd.jpl.nasa.gov/?planets\#elem
    // Basic : true
    // Scalar: true
    static const double VENUSSYSTEM_ORBITALSEMIMAJORAXIS_J2000 = 0.72333566; // [au]

    // Radial-velocity amplitude of the Sun induced by Venus for 'an edge-on observer' (i.e., an observer in the orbital plane of the Venus)
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 1.4, Equation 1.18, page 9
    // Basic : false
    // Scalar: true
    static const double VENUSSYSTEM_RADIALVELOCITYSIGNATURE = 0.086; // [m s^-1]

    // Radius of the smallest hypothetical sphere around Venus which would encompass the body (this is a low-accuracy parameter used in the relativistic model)
    // Basic : false
    // Scalar: true
    static const double VENUS_ENCOMPASSINGSPHERERADIUS = 6.051800e+06; // [m]

    // Equatorial radius of Venus
    // Basic : false
    // Scalar: true
    static const double VENUS_EQUATORIALRADIUS = 6.05180e+06; // [m]

    // Geometrical flattening factor f of Venus (f = (a-b)/a)
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double VENUS_FLATTENING = 0.;

    // Maximum reduction of the solar flux for an observer external to the solar system during a transit of Venus
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.2, Equation 2.4, page 14
    // Basic : false
    // Scalar: true
    static const double VENUS_FLUXREDUCTION_MAXIMUM = 0.008; // [%]

    // Geometric albedo of Venus (i.e., the ratio of the body's brightness at zero phase angle to the brightness of a perfectly diffusing disk with the same position and apparent size as the body)
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8 on page 706; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double VENUS_GEOMETRICALBEDO = 0.65;

    // Dynamical form-factor of Venus (oblateness or Stokes' second-degree zonal harmonic of the potential)
    // Source: P.K. Seidelmann, 1992, 'Explanatory Supplement to the Astronomical Almanac', University Science Books, Mill Valley, Ca., ISBN 0-935702-68-7, Table 15.8
    // Basic : true
    // Scalar: true
    static const double VENUS_JSUB2 = 2.70e-05;

    // Post-Newtonian deflection angle of a limb-grazing light ray due to the spherically symmetric part of the gravitational field of Venus
    // Source: S.A. Klioner, 2003, 'A Practical Relativistic Model for Microarcsecond Astrometry in Space', AJ, 125, 1580, Equation 64 with cot(x) = cos(x) / sin(x) approximated by x^-1; see also J.H.J. de Bruijne, 19 February 2002, 'Gravitational light deflection', GAIA-JdB-001, issue 1, revision 0
    // Basic : false
    // Scalar: true
    static const double VENUS_LIGHTDEFLECTION_LIMB = 493.; // [10^-6 arcsec]

    // Mass of Venus (do not use for high-precision (orbit) calculations)
    // Source: A.S. Konopliv, et al., 1999, 'Venus gravity: 180-th degree and order model', Icarus, 139, 3-18; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double VENUS_MASS = 4.867320e+24; // [kg]

    // Mean mass density of Venus
    // Basic : false
    // Scalar: true
    static const double VENUS_MASSDENSITY_MEAN = 5.243; // [g cm^-3]

    // IAU-recommended value for the declination \delta_0 of the north pole of rotation of Venus. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double VENUS_NORTHROTATIONALPOLE_DECLINATION = 67.16; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the declination \delta_0 of the north pole of rotation of Venus. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double VENUS_NORTHROTATIONALPOLE_DECLINATIONRATEOFCHANGE = 0.; // [deg day^-1]

    // IAU-recommended value for the right ascension \alpha_0 of the north pole of rotation of Venus. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double VENUS_NORTHROTATIONALPOLE_RIGHTASCENSION = 272.76; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch) of the right ascension \alpha_0 of the north pole of rotation of Venus. The pair (\alpha_0, \delta_0) denotes standard equatorial coordinates with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB). The north pole is that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double VENUS_NORTHROTATIONALPOLE_RIGHTASCENSIONRATEOFCHANGE = 0.; // [deg day^-1]

    // Polar radius of Venus
    // Basic : false
    // Scalar: true
    static const double VENUS_POLARRADIUS = 6.05180e+06; // [m]

    // IAU-recommended value for the ephemeris position of the prime meridian of Venus. The 0-deg meridian is defined by the central peak in the crater Ariadne. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double VENUS_PRIMEMERIDIAN_EPHEMERISPOSITION = 160.20; // [deg]

    // IAU-recommended value for the rate of change (in degrees per Julian day, calculated from the standard epoch of 1.5 January 2000 = JD2451545.0 TDB) of the ephemeris position of the prime meridian of Venus. The 0-deg meridian is defined by the central peak in the crater Ariadne. The location of the prime meridian is specified by the angle that is measured along the planet's equator in an easterly direction with respect to the planet's north pole from the node Q (located at right ascension 90 deg + \alpha_0, where \alpha_0 denotes the right ascension of the north pole of rotation) of the planet's equator on the standard equator to the point B where the prime meridian crosses the planet's equator. The right ascension of the point Q is 90 deg + \alpha_0 and the inclination of the planet's equator to the standard equator is 90 deg - \delta_0, where \delta_0 denotes the declination of the north pole of rotation. (The pair (\alpha_0, \delta_0) denotes the standard equatorial coordinates, with equinox J2000 at epoch J2000 (the standard epoch is 1.5 January 2000 = JD2451545.0 TDB), of the north pole of rotation, which itself is defined as that pole of rotation that lies on the north side of the invariable plane of the solar system; the approximate coordinates of the north pole of the invariable plane are \alpha_0 = 273.85 deg and \delta_0 = 66.99 degrees.) Because the prime meridian is assumed to rotate uniformly with the planet, W accordingly varies linearly with time. In addition, \alpha_0, \delta_0, and W may vary with time due to a precession of the axis of rotation of the planet. If W increases with time, the planet has a direct (or prograde) rotation; if W decreases with time, the rotation is said to be retrograde
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A)
    // Basic : true
    // Scalar: true
    static const double VENUS_PRIMEMERIDIAN_EPHEMERISPOSITIONRATEOFCHANGE = -1.4813688; // [deg day^-1]

    // Geometric transit probability (Venus transiting the Sun) for an observer external to the solar system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.1, Equation 2.2, page 14
    // Basic : false
    // Scalar: true
    static const double VENUS_TRANSITPROBABILITY = 0.648; // [%]

    // Maximum transit time of Venus (transiting the Sun) for an observer external to the solar system
    // Source: E.g., A. Johansen, 26 March 2002, 'Detection of planetary transits with the Gaia satellite', GAIA-CUO-106, issue 1, revision 3, Section 2.3, Equation 2.5, page 15
    // Basic : false
    // Scalar: true
    static const double VENUS_TRANSITTIME_MAXIMUM = 0.46; // [day]

    // V(1,0) magnitude of Venus (i.e., the visual magnitude of the planet reduced to a distance of 1 au from both the Sun and Earth and phase angle zero). This parameter is also refered to as absolute magnitude in planetary sciences
    // Source: J.L. Hilton, 2005, 'Improving the Visual Magnitudes of the Planets in The Astronomical Almanac. I. Mercury and Venus', AJ, 129, 2902-2906; see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double VENUS_VONEZEROMAGNITUDE = -4.47; // [mag]

    // Mean volumetric radius of Venus
    // Source: B.A. Archinal, et al., 1 February 2011, 'Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009', Celestial Mechanics and Dynamical Astronomy, 109, 101-135 (2011CeMDA.109..101A); see also D.K. Yeomans (NASA/JPL), 5 November 2008, 'Planets and Pluto: Physical Characteristics', http://ssd.jpl.nasa.gov/?planet_phys_par
    // Basic : true
    // Scalar: true
    static const double VENUS_VOLUMETRICRADIUS = 6.05180e+06; // [m]

    // Wien's displacement-law constant (for \lambda_max)
    // Source: P.J. Mohr, B.N. Taylor, D.B. Newell, 9 July 2015, 'The 2014 CODATA Recommended Values of the Fundamental Physical Constants', National Institute of Standards and Technology, Gaithersburg, MD 20899-8401; http://www.codata.org/ and http://physics.nist.gov/constants (Web Version 7.0)
    // Basic : true
    // Scalar: true
    static const double WIEN_CONSTANT = 2.89777290e-03; // [m K]

    // Zero degrees Celsius (ice point) expressed in degrees Kelvin
    // Basic : false
    // Scalar: true
    static const double ZEROCELSIUS_KELVIN = 273.15; // [K]

    // Zero degrees Kelvin expressed in degrees Celsius. The triple point of water is the only realizable defining fixed point common to the Kelvin Thermodynamic Temperature Scale (KTTS) and the International Temperature Scale of 1990 (ITS-90); the assigned value of the triple point of water on these scales is 273.16 K (0.01 C)
    // Source: The International Temperature Scale of 1990 (ITS-90); http://www.its-90.com/its-90.html
    // Basic : true
    // Scalar: true
    static const double ZEROKELVIN_CELSIUS = -273.15; // [deg C]

};

}
#endif
