#ifndef ORBITS_H
#define ORBITS_H

#include <array>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include "common.hpp"


class AnalyticOrbits{

public:

  double spacecraft_separation = 2.5E9;

  int number_of_spacecrafts = 3;
  int number_of_arms = 6;
  double inter_spacecraft_phase = 2 * M_PI / 3;

  array<double, 3> rotation;
  array<double, 3> c_rotation;
  array<double, 3> s_rotation;

  double eccentricity;
  double init_time;

  AnalyticOrbits();
  AnalyticOrbits(double arm_length, double init_position, double init_rotation);
  double alpha(double time);
  void position_x(int sci, double* time, double* x, int nt) ;
  void position_y(int sci, double* time, double* y, int nt) ;
  void position_z(int sci, double* time, double* z, int nt) ;
  array<double, 3> position(int sci, double time);

  void velocity_x(int sci, double* time, double* vx, int nt) ;
  void velocity_y(int sci, double* time, double* vy, int nt) ;
  void velocity_z(int sci, double* time, double* vz, int nt) ;
  array<double, 3> velocity(int sci, double time);

  void get_travel_time(int sci, int scj, double* rec_time, double* tt, int nt, int order);
  void get_travel_time(int sci, int scj, double* rec_time,
		       double* x_emitter, double* y_emitter, double* z_emitter,
		       double* x_receiver, double* y_receiver, double* z_receiver,
		       double* tt, int nt, int order);
};

#endif
