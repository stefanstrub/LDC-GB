//
// LISA Constants.
//
// This header provides values sanctioned by the LISA Consortium for physical constants and mission parameters.
//
// LISA Constants is intended to be consistently used by other pieces of software related to the simulation of
// the instrument, of gravitational wave signals, and others.
//
// Authors:
//    Jean-Baptiste Bayle <j2b.bayle@gmail.com>
//    Aurelien Hees <aurelien.hees@obspm.fr>
//    Maude Lejeune <lejeune@apc.in2p3.fr>
//
#ifndef LISACONSTANTS_H
#define LISACONSTANTS_H
#include <stdbool.h>

/**
 Speed of light in a vacuum
 Unit: m/s.

**/
static const double LISA_SPEED_OF_LIGHT = 299792458.0;

/**
 Speed of light in a vacuum
 Unit: m/s.

**/
static const double LISA_c = 299792458.0;

/**
 Number of days per sidereal year
 Unit: day.

**/
static const double LISA_SIDEREALYEAR_J2000DAY = 365.256363004;

/**
 Number of days per tropical year
 Unit: day.

**/
static const double LISA_TROPICALYEAR_J2000DAY = 365.242190402;

/**
 Astronomical year
 Unit: s.

**/
static const double LISA_ASTRONOMICAL_YEAR = 31558149.763545595;

/**
 Astronomical unit
 Unit: m.

**/
static const double LISA_ASTRONOMICAL_UNIT = 149597870700.0;

/**
 Astronomical unit
 Unit: m.

**/
static const double LISA_au = 149597870700.0;

/**
 Sun gravitational parameter
 Unit: m^3/s^2.

**/
static const double LISA_GM_SUN = 1.327124400419394e+20;

/**
 Sun Schwarzschild radius
 Unit: m.

**/
static const double LISA_SUN_SCHWARZSCHILD_RADIUS = 2953.2500770335273;

/**
 Parsec expressed in meters
 Unit: m.

**/
static const double LISA_PARSEC_METER = 3.085677581491367e+16;

/**
 Newton's universal constant of gravitation
 Unit: m^3/kg/s^2.

**/
static const double LISA_NEWTON_CONSTANT = 6.67408e-11;

/**
 Solar mass
 Unit: kg.

**/
static const double LISA_SUN_MASS = 1.98848e+30;

#endif