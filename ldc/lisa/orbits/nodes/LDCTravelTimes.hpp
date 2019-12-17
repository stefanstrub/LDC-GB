#pragma once

#include <array>
#include <cmath>
#include "Node.hpp"
#include "orbits.hpp"

using namespace std;
using namespace boost;

class LDCTravelTimes: public Node {

public:

  // MARK: Parameters
  int order = 0;
  
  // MARK: Inputs
  
  /** Positions along X of spacecraft, in m */
  array<circular_buffer<double>, 3> x;
  
  /** Positions along Y of spacecraft, in m */
  array<circular_buffer<double>, 3> y;
  
  /** Positions along Z of spacecraft, in m */
  array<circular_buffer<double>, 3> z;
  
  /** Velocities along X of spacecraft, in m */
  array<circular_buffer<double>, 3> vx;
  
  /** Velocities along Y of spacecraft, in m */
  array<circular_buffer<double>, 3> vy;
  
  /** Velocities along Z of spacecraft, in m */
  array<circular_buffer<double>, 3> vz;
  
  // MARK: Outputs
  
  /** Matrix of travel times, in s */
  array<array<double, 3>, 3> travel_times;
  
  
  void prepare() {
    for (int spacecraft_index = 0; spacecraft_index < 3; spacecraft_index += 1) {
      x[spacecraft_index] = circular_buffer<double>(1, 0.0);
      y[spacecraft_index] = circular_buffer<double>(1, 0.0);
      z[spacecraft_index] = circular_buffer<double>(1, 0.0);
      vx[spacecraft_index] = circular_buffer<double>(1, 0.0);
      vy[spacecraft_index] = circular_buffer<double>(1, 0.0);
      vz[spacecraft_index] = circular_buffer<double>(1, 0.0);
    }
  }
  
  void fire(double time) {
    
    for (int spacecraft1 = 0; spacecraft1 < 3; spacecraft1 += 1) {
      auto& from_spacecraft = travel_times[spacecraft1];
      for (int spacecraft2 = 0; spacecraft2 < 3; spacecraft2 += 1) {
	from_spacecraft[spacecraft2] = travel_time_between(time, spacecraft1, spacecraft2);
      }
    }
  }
  
protected:

    /**
     Compute travel time between two spacecraft.

     \param spacecraft1 Emitter spacecraft.
     \param spacecraft2 Receiver spacecraft.
     \return Light travel time between both spacecraft.
     */
    double travel_time_between(double time, int spacecraft1, int spacecraft2) {

        // No-op if same spacecraft
        if (spacecraft1 == spacecraft2) {
            return 0.0;
        }

        array<double, 3> emitter_position = { x[spacecraft1][0],
					      y[spacecraft1][0],
					      z[spacecraft1][0] };
        array<double, 3> receiver_position = { x[spacecraft2][0],
					       y[spacecraft2][0],
					       z[spacecraft2][0] };
	array<double, 3> receiver_velocity = {
	  vx[spacecraft2][0], vy[spacecraft2][0], vz[spacecraft2][0]};

	double tt = travel_time(emitter_position, receiver_position, receiver_velocity, order);
        return tt;
    }


};
