#include "common.hpp"

double norm(array<double, 3> p){
  double res = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
  return res;
}

double dot_product(array<double, 3> v1, array<double, 3> v2){
  double res = 0;
  for (int i=0; i < 3; i += 1)
    res += v1[i]*v2[i];
      return res;
}

double travel_time(array<double, 3> r_i, array<double, 3> r_j, array<double, 3> v_j, int order) {
  array<double, 3> r_ij;
  for (int i=0; i < 3; i+= 1)
    r_ij[i] = r_j[i] - r_i[i];

  double n = norm(r_ij);
  double tt = n / (double)CLIGHT;
  
  if (order>0){
    for (int i=0; i < 3; i+= 1)
      v_j[i] /= CLIGHT;
      
    array<double, 3> n_ij;
    for (int i=0; i < 3; i+= 1)
      n_ij[i] = r_ij[i] / n;

    if (order>1){
      double c_21 = 0.5 * tt * (dot_product(v_j, v_j) + pow(dot_product(n_ij, v_j), 2.0));
      double c_22 = -0.5 * tt*tt * CLIGHT * HALF_SCHWARZSCHILD_RADIUS * \
	dot_product(n_ij, r_j) / pow(dot_product(r_j, r_j), 1.5);
      
      double numerator = sqrt(pow(CLIGHT*tt + dot_product(n_ij, r_i), 2.0) + \
			      dot_product(r_i, r_i) - pow(dot_product(r_i,n_ij), 2.0))+\
	CLIGHT * tt + dot_product(n_ij, r_i);
      double denominator = dot_product(n_ij, r_i) + sqrt(dot_product(r_i, r_i));
      
      double c_23 = (HALF_SCHWARZSCHILD_RADIUS / CLIGHT) *		\
	(1.0 + POST_NEWTONIAN_CONSTANT) * log(numerator / denominator);
      tt += (c_21 + c_22 + c_23) + tt*dot_product(n_ij, v_j);
    }
    else{
      tt *= (1 + dot_product(n_ij, v_j));
    }
  }
  return tt;
}
