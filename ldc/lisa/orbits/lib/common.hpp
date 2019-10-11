#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>

using namespace std;


#define AU_IN_M 149597870700.0
#define BARYCENTER_ANGULAR_VELOCITY 1.9909865927683785e-07 //2pi/yrsid_si 
#define ASTRONOMICAL_YEAR 31558149.7635456
#define CLIGHT 299792458.0
#define HALF_SCHWARZSCHILD_RADIUS 1.47664E3
#define POST_NEWTONIAN_CONSTANT 1.0

double norm(array<double, 3> v);
double dot_product(array<double, 3> v1, array<double, 3> v2);
double travel_time(array<double, 3> r_i, array<double, 3> r_j,array<double, 3> v_j, int order);

