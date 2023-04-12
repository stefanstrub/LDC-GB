#include "orbits.hpp"
#include "c_wrapper.h"

extern "C" {
  AnalyticOrbits* newAnalyticOrbits(double arm_length, double init_position,
				    double init_rotation) {
    return new AnalyticOrbits(arm_length, init_position, init_rotation);
  }
  
  double AnalyticOrbits_get_armlength(AnalyticOrbits* v) {
    return v->spacecraft_separation;
  }
  double AnalyticOrbits_get_position_x(AnalyticOrbits* v, int sci, double t) {
    double x = 0;
    v->position_x(sci, &t, &x, 1);
    return x;
  }
  double AnalyticOrbits_get_position_y(AnalyticOrbits* v, int sci, double t) {
    double y = 0;
    v->position_y(sci, &t, &y, 1);
    return y;
  }
  double AnalyticOrbits_get_position_z(AnalyticOrbits* v, int sci, double t) {
    double z = 0;
    v->position_z(sci, &t, &z, 1);
    return z;
  }
  
  void deleteAnalyticOrbits(AnalyticOrbits* v) {
    delete v;
  }
}
