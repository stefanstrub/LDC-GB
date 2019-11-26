#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include "PdbParam.hpp"

using namespace std;
using namespace PdbParam;

#define BARYCENTER_ANGULAR_VELOCITY 2*Nature::PI_CONSTEXPRANT/Nature::SIDEREALYEAR_J2000DAY //1.9909865927683785e-07 //2pi/yrsid_si 
#define ASTRONOMICAL_YEAR Nature::SIDEREALYEAR_J2000DAY*24*60*60 //31558149.7635456
#define CLIGHT Nature::VELOCITYOFLIGHT_CONSTEXPRANT_VACUUM
#define HALF_SCHWARZSCHILD_RADIUS 1.47664E3
#define POST_NEWTONIAN_CONSTANT 1.0

double norm(array<double, 3> v);
double dot_product(array<double, 3> v1, array<double, 3> v2);
double travel_time(array<double, 3> r_i, array<double, 3> r_j,array<double, 3> v_j, int order);

