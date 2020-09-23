#ifndef __C_WRAPPER_H
#define __C_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct AnalyticOrbits AnalyticOrbits;

AnalyticOrbits* newAnalyticOrbits(double arm_length, double init_position, double init_rotation);

  double MyClass_int_get(AnalyticOrbits* v, double t);

void deleteAnalyticOrbits(AnalyticOrbits* v);

#ifdef __cplusplus
}
#endif

#endif
