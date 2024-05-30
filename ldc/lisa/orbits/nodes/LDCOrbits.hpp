#pragma once

#include <array>
#include <cmath>
#include "Node.hpp"
#include "orbits.hpp"


using namespace std;

class LDCOrbits: public Node, public AnalyticOrbits {

public:

  // MARK: Parameters cf AnalyticOrbits
  double arm_length = 2.5e9;
  double init_position = 0;
  double init_rotation = 0;
  
  // MARK: Outputs
  
  /** Positions along X of spacecraft, in m */
  array<double, 3> x;
  
  /** Positions along Y of spacecraft, in m */
  array<double, 3> y;
  
  /** Positions along Z of spacecraft, in m */
  array<double, 3> z;
  
  /** Velocities along X of spacecraft, in m */
  array<double, 3> vx;
  
  /** Velocities along Y of spacecraft, in m */
  array<double, 3> vy;
  
  /** Velocities along Z of spacecraft, in m */
  array<double, 3> vz;

  void prepare() {
    initialize_orbit_parameters();
  }

  void fire(double time) {
    // Update position and velocity for each spacecraft
    for (int spacecraft_index = 0; spacecraft_index < 3; spacecraft_index += 1) {
      compute_spacecraft_state(spacecraft_index, time);
    }
  }
  
protected:


  
  /**
     Compute internal orbital parameters.
     
     These are computed once before similation starts from the node's params.
     They are stored as part of the internal (private) state of the node.
  */
  void initialize_orbit_parameters() {
    //AnalyticOrbits(arm_length, init_position, init_rotation);
    // TODO avoid code duplication here.
    for (int i = 0; i < number_of_spacecrafts; i += 1) {
      rotation[i] = i * inter_spacecraft_phase + init_rotation;
      c_rotation[i] = cos(rotation[i]);
      s_rotation[i] = sin(rotation[i]);
    }
    eccentricity = arm_length/(2*sqrt(3)*AU_IN_M);
    init_time = init_position * ASTRONOMICAL_YEAR / (2 * M_PI);
  }
  
  void compute_spacecraft_state(int spacecraft_index, double time) {

    array<double, 3> p;
    p = position(spacecraft_index+1, time);
    // Store spacecraft position
    x[spacecraft_index] = p[0];
    y[spacecraft_index] = p[1];
    z[spacecraft_index] = p[2];
    
    // Store spacecraft velocity
    array<double, 3> v;
    v = velocity(spacecraft_index+1, time);
    vx[spacecraft_index] = v[0];
    vy[spacecraft_index] = v[1];
    vz[spacecraft_index] = v[2];
  }
  
};
