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
#pragma once
#include <string>

namespace LisaConstants {

class Constants {
public:

/**
 Speed of light in a vacuum
 Unit: m/s.

**/
static constexpr double SPEED_OF_LIGHT = 299792458.0;

/**
 Speed of light in a vacuum
 Unit: m/s.

**/
static constexpr double c = 299792458.0;

/**
 Number of days per sidereal year
 Unit: day.

**/
static constexpr double SIDEREALYEAR_J2000DAY = 365.256363004;

/**
 Number of days per tropical year
 Unit: day.

**/
static constexpr double TROPICALYEAR_J2000DAY = 365.242190402;

/**
 Astronomical year
 Unit: s.

**/
static constexpr double ASTRONOMICAL_YEAR = 31558149.763545595;

/**
 Astronomical unit
 Unit: m.

**/
static constexpr double ASTRONOMICAL_UNIT = 149597870700.0;

/**
 Astronomical unit
 Unit: m.

**/
static constexpr double au = 149597870700.0;

/**
 Sun gravitational parameter
 Unit: m^3/s^2.

**/
static constexpr double GM_SUN = 1.327124400419394e+20;

/**
 Sun Schwarzschild radius
 Unit: m.

**/
static constexpr double SUN_SCHWARZSCHILD_RADIUS = 2953.2500770335273;

/**
 Parsec expressed in meters
 Unit: m.

**/
static constexpr double PARSEC_METER = 3.085677581491367e+16;

/**
 Newton's universal constant of gravitation
 Unit: m^3/kg/s^2.

**/
static constexpr double NEWTON_CONSTANT = 6.67408e-11;

/**
 Solar mass
 Unit: kg.

**/
static constexpr double SUN_MASS = 1.98848e+30;

};
}