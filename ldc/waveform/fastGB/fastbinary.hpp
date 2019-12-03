#ifndef FASTBINARY_H

#include "PdbParam.hpp"

using namespace PdbParam;

// constants
const double pi = Nature::PI_CONSTEXPRANT;
const double sq3 = 1.73205080757;
const double clight = Nature::VELOCITYOFLIGHT_CONSTEXPRANT_VACUUM;  // m/s
const double AU = Nature::ASTRONOMICALUNIT_METER;   //  m
const double year = Nature::SIDEREALYEAR_J2000DAY*24*60*60; // siderial year in s
const double fm = 1.0/year;

const double kappa = 0.0;           // initial azimuthal position of the guiding center
const double lambda = 0.0;          // initial orientation of the LISA constellation

extern "C" {
    #include <fftw3.h>
}


class FastBinary {
  private:
    long N, M;
    double T, dt;

    double *u,*v,*k;    // Gravitational Wave basis vectors
    double *kdotx, **kdotr;  // Dot products
    // Distance, gravitational wave frequency & ratio of f and transfer frequency f*
    double *xi, *fonfs;
    // Polarization basis tensors, convenient quantities
    double **eplus, **ecross, **dplus, **dcross;    
    double *x, *y, *z;  // Spacecraft position and separation vectors
    double **xv, **yv, **zv;
    double *r12, *r13, *r21, *r23, *r31, *r32;
    // Time varying quantities (Re & Im) broken up into convenient segments
    double **TR, **TI;
    // Fourier coefficients before FFT and after convolution:
    double *data12, *data13, *data21, *data23, *data31, *data32;    
    // Fourier coefficients of slowly evolving terms (numerical)
    double *a12, *a13, *a21, *a23, *a31, *a32;                      
    double *b; // Fourier coefficients of rapidly evolving terms (analytical)
    double *an, *bn; // MV: Fourier transforms for convolve_fft
    // Fourier coefficients of entire response (convolution)
    double *c12, *c13, *c21, *c23, *c31, *c32; 
    double *ReA, *ImA, *ReB, *ImB, *ReC, *ImC;
    double *X, *Y, *Z;

  double L = 2.5e9;   // Assuming analytical orbit.
  double fstar = clight/(2.0*pi*L);
  double ec = L/(2.0*AU*sqrt(3.0));
  
    fftw_complex *in, *out;
    fftw_plan plan_forward, plan_backward;

    void spacecraft(double t);

    void convolve(double *a, double *b, double *cn, int method);

    void XYZ(double f0, long q, double *XLS, double *XSL,
	     double *YLS, double *YSL, double *ZLS, double *ZSL);

  public:
  FastBinary(long Nreq,double Treq,double dtreq);
  ~FastBinary();
  void response(double f0,double fdot,double theta,double phi,
		double A,double iota,double psi,double phio,
		double *XLS, double *XSL,
		double *YLS, double *YSL,
		double *ZLS, double *ZSL, long len,
		int method);

  void setL(double l);
};

#define FASTBINARY_H
#endif
