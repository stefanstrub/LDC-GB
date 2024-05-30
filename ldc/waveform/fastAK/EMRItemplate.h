/*
Copyright (C) 2005  S. Babak

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/


/******************  CVS info ************************
#define CVSTAG "$Name:  $"
#define CVSHEADER "$Header: /afs/aeiw/cvsroot/waves/people/stas/LISAWP/base/include/EMRItemplate.hh,v 1.2 2008/04/18 11:48:08 stba Exp $"
*/

#ifndef EMRITEMPLATEHH
#define EMRITEMPLATEHH

#include <vector>
#include <math.h>
#include <complex>
#include "lisaconstants.hpp" // LISA constants

using namespace LisaConstants;

#define LISAWP_C_SI Constants::SPEED_OF_LIGHT /* Speed of light in vacuo, m s^-1 */
#define GMSUN  Constants::GM_SUN
#define LISAWP_MTSUN_SI GMSUN/(LISAWP_C_SI*LISAWP_C_SI*LISAWP_C_SI) /* Geometrized solar mass, s */
#define LISAWP_PC_SI Constants::PARSEC_METER /* Parsec, m */

#define LISAWP_PI         3.1415926535897932384626433832795029L  /* pi */
#define LISAWP_TWOPI      6.2831853071795864769252867665590058L  /* 2*pi */
#define LISAWP_PI_2       1.5707963267948966192313216916397514L  /* pi/2 */
#define LISAWP_PI_4       0.7853981633974483096156608458198757L  /* pi/4 */
#define LISAWP_1_PI       0.3183098861837906715377675267450287L  /* 1/pi */
#define LISAWP_2_PI       0.6366197723675813430755350534900574L  /* 2/pi */
#define LISAWP_AU_SI Constants::ASTRONOMICAL_UNIT /* Astronomical unit, m */
//#define LISA_arm    2.5e9 /* LISA armlength in m */
// #define year  31457280.0
#define LISAWP_YRSID_SI (double)(Constants::SIDEREALYEAR_J2000DAY*24*60*60)  /* Sidereal year, s */


  /** Template class
  * @author S. Babak, 2006
  */

class EMRItemplate{

  public:

  EMRItemplate();

  /** Costructor
  * @param Mtotal total mass in solar masses
  * @param eta reduced mass ratio
  * @param lambda in radian
  * @param beta in radian
  * @param Tc time to coalescence in seconds
  */

  EMRItemplate(double Mass, double mu, double spin, double lambda, double quadmom);

//  void SetInitialConditions(double nu, double e, double Phi, double gamma, double alpha);

  void SetPosition(double thetaS, double phiS, double thetaK, double phiK, double D);

  EMRItemplate(const EMRItemplate& p);

  EMRItemplate operator=(const EMRItemplate& p);

//  EMRItemplate operator=(const EMRItemplate&);

  double M;
  double Mt;
  double m;
  double mt;
  double a;
  double lam;

  double e0;
  double nu0;
  double Phi0;
  double gamma0;
  double alpha0;
//  double tStart; // start of observation
  double t0;  // instance of time when initial conditions are defined
  double fgam0;
  double falph0;
  double qm;

  double tPl;
  double e_pl;
  double nu_pl;
  double Phi_pl;
  double alpha_pl;
  double gamma_pl;

  double thS;   // co-latittude !!!
  double phS;
  double thK;
  double phK;
  double stS, stK, ctS, ctK; // cos and sin of thetaS and thetaK
  double cpS, spS, cpK, spK; // cos and sin of phiS, and phiK

  double dist; // distance in seconds
  double Ampl; // dimensionless amplitude: mu/D
  bool SkySet;

  double SNR;
  double LogL;


};


#endif