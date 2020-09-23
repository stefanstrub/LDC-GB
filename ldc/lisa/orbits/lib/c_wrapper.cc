#include "orbits.hpp"
#include "c_wrapper.h"

extern "C" {
        AnalyticOrbits* newAnalyticOrbits(double arm_length, double init_position, double init_rotation) {
	  return new AnalyticOrbits(arm_length, init_position, init_rotation);
        }

  double MyClass_int_get(AnalyticOrbits* v, double t) {
                return v->alpha(t);
        }

        void deleteAnalyticOrbits(AnalyticOrbits* v) {
                delete v;
        }
}
