#ifndef __C_WRAPPER_H
#define __C_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct AnalyticOrbits AnalyticOrbits;

AnalyticOrbits* newAnalyticOrbits(double arm_length, double init_position, double init_rotation);

  double AnalyticOrbits_get_armlength(AnalyticOrbits* v);
  double AnalyticOrbits_get_position_x(AnalyticOrbits* v, int spacecraft_index, double t);
  double AnalyticOrbits_get_position_y(AnalyticOrbits* v, int spacecraft_index, double t);
  double AnalyticOrbits_get_position_z(AnalyticOrbits* v, int spacecraft_index, double t);
  
  
void deleteAnalyticOrbits(AnalyticOrbits* v);

#ifdef __cplusplus
}
#endif

#endif
