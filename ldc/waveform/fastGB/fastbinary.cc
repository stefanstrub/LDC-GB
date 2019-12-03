#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fastbinary.hpp"

extern "C" {
    #include <fftw3.h>
    #include "arrays.h"
}

// LISA -- e = L / (2 AU sqrt(3)); fstar = c / (2 pi L)
//double L = 5.0e9;             // armlength
double L = 2.5e9;             // armlength
//double fstar = 0.00954269032; // transfer frequency
//double ec = 0.009648370435;   // eccentricity



// Assuming analytical orbit.
double fstar = clight/(2.0*pi*L);
double ec = L/(2.0*AU*sqrt(3.0));

void setL(double l) {
    L = l;

    ec = L / (2.0 * AU * sq3);
    fstar = clight / (2.0 * pi * L);
}
FastBinary::FastBinary(){} //nullary constructor for cython

FastBinary::FastBinary(long Nreq,double Treq,double dtreq) : N(Nreq), M(Nreq), T(Treq), dt(dtreq) {
    u = dvector(1,3); v = dvector(1,3); k = dvector(1,3);

    kdotx = dvector(1,3); kdotr = dmatrix(1,3,1,3);

    xi = dvector(1,3); fonfs = dvector(1,3);

    eplus  = dmatrix(1,3,1,3); ecross = dmatrix(1,3,1,3);
    dplus  = dmatrix(1,3,1,3); dcross = dmatrix(1,3,1,3);

    x = dvector(1,3); y = dvector(1,3); z = dvector(1,3);

    r12 = dvector(1,3); r21 = dvector(1,3); r31 = dvector(1,3);
    r13 = dvector(1,3); r23 = dvector(1,3); r32 = dvector(1,3);

    TR = dmatrix(1,3,1,3); TI = dmatrix(1,3,1,3);

    // these are the long vectors...

    data12 = dvector(1,2*N); data21 = dvector(1,2*N); data31 = dvector(1,2*N);
    data13 = dvector(1,2*N); data23 = dvector(1,2*N); data32 = dvector(1,2*N);

    b = dvector(1,2*M);

    c12 = dvector(1,2*M); c21 = dvector(1,2*M); c31 = dvector(1,2*M);
    c13 = dvector(1,2*M); c23 = dvector(1,2*M); c32 = dvector(1,2*M);

    in  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);

    plan_forward  = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD,  FFTW_ESTIMATE);
    plan_backward = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    X = dvector(1,2*M);  Y = dvector(1,2*M);  Z = dvector(1,2*M);
};

FastBinary::~FastBinary() {
    fftw_destroy_plan(plan_backward);
    fftw_destroy_plan(plan_forward);

    fftw_free(in); fftw_free(out);

    free_dvector(u,1,3); free_dvector(v,1,3); free_dvector(k,1,3);

    free_dvector(kdotx,1,3); free_dmatrix(kdotr,1,3,1,3);

    free_dvector(xi,1,3);

    free_dvector(fonfs,1,3);

    free_dmatrix(eplus,1,3,1,3); free_dmatrix(ecross,1,3,1,3);

    free_dmatrix(dplus,1,3,1,3); free_dmatrix(dcross,1,3,1,3);

    free_dvector(x,1,3); free_dvector(y,1,3); free_dvector(z,1,3);

    free_dvector(r12,1,3); free_dvector(r21,1,3); free_dvector(r31,1,3);
    free_dvector(r13,1,3); free_dvector(r23,1,3); free_dvector(r32,1,3);

    free_dmatrix(TR,1,3,1,3); free_dmatrix(TI,1,3,1,3);

    free_dvector(data12,1,2*N); free_dvector(data21,1,2*N); free_dvector(data31,1,2*N);
    free_dvector(data13,1,2*N); free_dvector(data23,1,2*N); free_dvector(data32,1,2*N);

    free_dvector(b,1,2*M+2);

    free_dvector(c12,1,2*M+2); free_dvector(c21,1,2*M+2); free_dvector(c31,1,2*M+2);
    free_dvector(c13,1,2*M+2); free_dvector(c23,1,2*M+2); free_dvector(c32,1,2*M+2);

    free_dvector(X,1,2*M);  free_dvector(Y,1,2*M);  free_dvector(Z,1,2*M);
};

// note that this works only if M == N, appropriate for algorithm=='legacy'/method==0
void FastBinary::convolve(double *a,double *bn,double *cn,int method) {
    // multiply, inverting frequencies in a

    in[0][0] = a[2*0+1]*bn[2*0+1] - a[2*0+2]*bn[2*0+2];
    in[0][1] = a[2*0+1]*bn[2*0+2] + a[2*0+2]*bn[2*0+1];

    for(int n=1; n<N; n++) {
        in[n][0] = a[2*(N-n)+1]*bn[2*n+1] - a[2*(N-n)+2]*bn[2*n+2];
        in[n][1] = a[2*(N-n)+1]*bn[2*n+2] + a[2*(N-n)+2]*bn[2*n+1];
    }

    // run the inverse FFT
    fftw_execute(plan_backward);

    for(int n=0; n<N; n++) {
        cn[2*n+1] = out[n][0] / N;
        cn[2*n+2] = out[n][1] / N;
    }
}

void FastBinary::response(double f,double fdot,double theta,double phi,double A,
			  double iota,double psi,double phio,
			  double *XLS,double *XSL,
			  double *YLS,double *YSL,
			  double *ZLS,double *ZSL,long len,
			  int method=0)
{
  // Calculate cos and sin of sky position, inclination, polarization
    double costh = cos(theta);   double sinth = sin(theta);
    double cosph = cos(phi);     double sinph = sin(phi);
    double cosps = cos(2.*psi);  double sinps = sin(2.*psi);
    double cosi  = cos(iota);

    // Note: We don't bother to compute the amplitude evolution as it is immeasurably small
    double Aplus  = A*(1.+cosi*cosi);
    double Across = -2.0*A*cosi;

    // with the new computing method, shift the central frequency to its half-evolved value
    double f0 = (method == 0) ? f : (f + 0.5 * fdot * T);

    // Calculate carrier frequency bin
    long q = (long)(f0*T);

    // Calculate constant pieces of transfer functions
    double DPr =  Aplus*cosps;
    double DPi = -Across*sinps;
    double DCr = -Aplus*sinps;
    double DCi = -Across*cosps;

    /* Tensor stuff */
    u[1] =  costh*cosph;  u[2] =  costh*sinph;    u[3] = -sinth;
    v[1] =  sinph;        v[2] = -cosph;          v[3] =  0.0;
    k[1] = -sinth*cosph;  k[2] = -sinth*sinph;    k[3] = -costh;

    for(int i=1; i<=3; i++) {
        for(int j=1; j<=3; j++) {
            eplus[i][j]  = u[i]*u[j] - v[i]*v[j];
            ecross[i][j] = u[i]*v[j] + v[i]*u[j];
        }
    }

    /* Analytical Fourier transform of high frequency components.
       TO DO: replace this and the FFT with the time-domain evaluation of the "slow" exponential. */
    for(int i=1; i<=M; i++) {
        int m = q + i-1 - M/2;
        double xm = pi*(f0*T - (double)m);
        double sinc = (xm != 0.0) ? (sin(xm) / xm) : 1.0;

        b[2*i-1] = cos(xm) * sinc;
        b[2*i]   = sin(xm) * sinc;
    }

    // this B transform needs to be done only once
    for(int n=0; n<N; n++) { in[n][0] = b[2*n+1]; in[n][1] = b[2*n+2]; }
    fftw_execute(plan_forward);
    for(int n=0; n<N; n++) { b[2*n+1] = out[n][0]; b[2*n+2] = out[n][1]; }

    // TO DO: cache the spacecraft position vectors
    for(int n=1; n<=N; n++) {
        double t = T * (double)(n-1) / (double)N;

        // Calculate position of each spacecraft at time t
        spacecraft(t);

        for(int i=1; i<=3; i++)  {
            kdotx[i] = (x[i]*k[1]+y[i]*k[2]+z[i]*k[3])/clight;
            xi[i]    = t - kdotx[i];
            fonfs[i] = (f + fdot*xi[i]) / fstar;    // works with both methods...
        }

        // Unit separation vector from spacecrafts i to j
        r12[1] = (x[2] - x[1])/L;   r13[1] = (x[3] - x[1])/L;   r23[1] = (x[3] - x[2])/L;
        r12[2] = (y[2] - y[1])/L;   r13[2] = (y[3] - y[1])/L;   r23[2] = (y[3] - y[2])/L;
        r12[3] = (z[2] - z[1])/L;   r13[3] = (z[3] - z[1])/L;   r23[3] = (z[3] - z[2])/L;

        // Make use of symmetry
        for(int i=1; i<=3; i++) {
            r21[i] = -r12[i];
            r31[i] = -r13[i];
            r32[i] = -r23[i];
        }

        // TO DO: is this needed?
        dplus[1][2]  = dplus[1][3]  = dplus[2][1]  = dplus[2][3]  = dplus[3][1]  = dplus[3][2]  = 0.0;
        dcross[1][2] = dcross[1][3] = dcross[2][1] = dcross[2][3] = dcross[3][1] = dcross[3][2] = 0.0;

        // Convenient quantities d+ & dx
        for(int i=1; i<=3; i++) {
            for(int j=1; j<=3; j++) {
                dplus[1][2]  += r12[i]*r12[j]*eplus[i][j];   dcross[1][2] += r12[i]*r12[j]*ecross[i][j];
                dplus[2][3]  += r23[i]*r23[j]*eplus[i][j];   dcross[2][3] += r23[i]*r23[j]*ecross[i][j];
                dplus[1][3]  += r13[i]*r13[j]*eplus[i][j];   dcross[1][3] += r13[i]*r13[j]*ecross[i][j];
            }
        }

        dplus[2][1] = dplus[1][2];  dcross[2][1] = dcross[1][2];
        dplus[3][2] = dplus[2][3];  dcross[3][2] = dcross[2][3];
        dplus[3][1] = dplus[1][3];  dcross[3][1] = dcross[1][3];

        kdotr[1][2] = kdotr[1][3] = kdotr[2][1] = kdotr[2][3] = kdotr[3][1] = kdotr[3][2] = 0.0;

        for(int i=1; i<=3; i++) {
            kdotr[1][2] += k[i]*r12[i];   kdotr[1][3] += k[i]*r13[i];   kdotr[2][3] += k[i]*r23[i];
        }

        kdotr[2][1] = -kdotr[1][2]; kdotr[3][1] = -kdotr[1][3]; kdotr[3][2] = -kdotr[2][3];

        // Calculating Transfer function
        for(int i=1; i<=3; i++) {
            for(int j=1; j<=3; j++) {
                if(i!=j) {
                    double arg1 = 0.5*fonfs[i]*(1 - kdotr[i][j]);

                    // first term is phase evolution from fdot; second term is initial phase; third term is Doppler
                    double arg2 = pi*fdot*xi[i]*xi[i] + phio - 2.*pi*kdotx[i]*f0;
                    arg2 += (method == 0) ? 0 : -pi*T*fdot*xi[i];   // additional contribution from recentering f0

                    double sinc = (arg1 != 0.0) ? (sin(arg1) / arg1) : 1.0;

                    double tran1r = dplus[i][j]*DPr + dcross[i][j]*DCr;
                    double tran1i = dplus[i][j]*DPi + dcross[i][j]*DCi;

                    double tran2r = cos(arg1 + arg2);
                    double tran2i = sin(arg1 + arg2);

                    // Real & Imaginary part of the slowly evolving signal
                    TR[i][j] = 0.25*sinc*(tran1r*tran2r - tran1i*tran2i);
                    TI[i][j] = 0.25*sinc*(tran1r*tran2i + tran1i*tran2r);
                }
            }
        }

        data12[2*n-1] = TR[1][2];   data21[2*n-1] = TR[2][1];   data31[2*n-1] = TR[3][1];
        data12[2*n]   = TI[1][2];   data21[2*n]   = TI[2][1];   data31[2*n]   = TI[3][1];
        data13[2*n-1] = TR[1][3];   data23[2*n-1] = TR[2][3];   data32[2*n-1] = TR[3][2];
        data13[2*n]   = TI[1][3];   data23[2*n]   = TI[2][3];   data32[2*n]   = TI[3][2];
    }

    // TO DO: these could be multithreaded, but we'll need more FFTW plans
    convolve(data12,b,c12,method);  convolve(data21,b,c21,method);  convolve(data31,b,c31,method);
    convolve(data13,b,c13,method);  convolve(data23,b,c23,method);  convolve(data32,b,c32,method);

    // Renormalize so that the time series is real (could do it within convolve I suppose)
    for(int i=1; i<=2*M; i++) {
        c12[i] = 0.5*c12[i];  c21[i] = 0.5*c21[i];  c31[i] = 0.5*c31[i];
        c13[i] = 0.5*c13[i];  c23[i] = 0.5*c23[i];  c32[i] = 0.5*c32[i];
    }

    /* Call subroutines for synthesizing different TDI data channels */

    /* X Y Z-Channel: note the pointer change since we're passed 0-based arrays */
    XYZ(f0, q, --XLS, --XSL, --YLS, --YSL, --ZLS, --ZSL);

    return;
}

/* Rigid approximation position of each LISA spacecraft */
/* Could be cached... */
void FastBinary::spacecraft(double t) {
    double alpha = 2.0*pi*fm*t + kappa;

    double beta1 = 0.0 + lambda;
    double beta2 = 2.0*pi/3. + lambda;
    double beta3 = 4.0*pi/3. + lambda;

    double sa = sin(alpha);
    double ca = cos(alpha);

    double sb, cb;

    sb = sin(beta1);
    cb = cos(beta1);

    x[1] = AU*ca + AU*ec*(sa*ca*sb - (1. + sa*sa)*cb);
    y[1] = AU*sa + AU*ec*(sa*ca*cb - (1. + ca*ca)*sb);
    z[1] = -sq3*AU*ec*(ca*cb + sa*sb);

    sb = sin(beta2);
    cb = cos(beta2);

    x[2] = AU*ca + AU*ec*(sa*ca*sb - (1. + sa*sa)*cb);
    y[2] = AU*sa + AU*ec*(sa*ca*cb - (1. + ca*ca)*sb);
    z[2] = -sq3*AU*ec*(ca*cb + sa*sb);

    sb = sin(beta3);
    cb = cos(beta3);

    x[3] = AU*ca + AU*ec*(sa*ca*sb - (1. + sa*sa)*cb);
    y[3] = AU*sa + AU*ec*(sa*ca*cb - (1. + ca*ca)*sb);
    z[3] = -sq3*AU*ec*(ca*cb + sa*sb);
}

void FastBinary::XYZ(double f0, long q, double *XLS, double *XSL,
		     double *YLS, double *YSL, double *ZLS, double *ZSL) {
    double phiLS = 2.0*pi*f0*(dt/2.0-L/clight);
    double cLS = cos(phiLS); double sLS = sin(phiLS);

    double phiSL = pi/2.0-2.0*pi*f0*(L/clight);
    double cSL = cos(phiSL); double sSL = sin(phiSL);

    for(int i=1; i<=M; i++) {
        double f = ((double)(q + i-1 - M/2))/T;
        double fonfs = f/fstar;

        double c3 = cos(3.0*fonfs);  double c2 = cos(2.0*fonfs);  double c1 = cos(1.0*fonfs);
        double s3 = sin(3.0*fonfs);  double s2 = sin(2.0*fonfs);  double s1 = sin(1.0*fonfs);

        X[2*i-1] = (c12[2*i-1] - c13[2*i-1])*c3 + (c12[2*i]   - c13[2*i]) * s3 +
                   (c21[2*i-1] - c31[2*i-1])*c2 + (c21[2*i]   - c31[2*i]) * s2 +
                   (c13[2*i-1] - c12[2*i-1])*c1 + (c13[2*i]   - c12[2*i]) * s1 +
                   (c31[2*i-1] - c21[2*i-1]);

        X[2*i]   = (c12[2*i]   - c13[2*i]) * c3 - (c12[2*i-1] - c13[2*i-1])*s3 +
                   (c21[2*i]   - c31[2*i]) * c2 - (c21[2*i-1] - c31[2*i-1])*s2 +
                   (c13[2*i]   - c12[2*i]) * c1 - (c13[2*i-1] - c12[2*i-1])*s1 +
                   (c31[2*i]   - c21[2*i]);


        Y[2*i-1] = (c23[2*i-1] - c21[2*i-1])*c3 + (c23[2*i]   - c21[2*i]) * s3 +
                   (c32[2*i-1] - c12[2*i-1])*c2 + (c32[2*i]   - c12[2*i]) * s2 +
                   (c21[2*i-1] - c23[2*i-1])*c1 + (c21[2*i]   - c23[2*i]) * s1 +
                   (c12[2*i-1] - c32[2*i-1]);

        Y[2*i]   = (c23[2*i]   - c21[2*i]) * c3 - (c23[2*i-1] - c21[2*i-1])*s3 +
                   (c32[2*i]   - c12[2*i]) * c2 - (c32[2*i-1] - c12[2*i-1])*s2 +
                   (c21[2*i]   - c23[2*i]) * c1 - (c21[2*i-1] - c23[2*i-1])*s1 +
                   (c12[2*i]   - c32[2*i]);


        Z[2*i-1] = (c31[2*i-1] - c32[2*i-1])*c3 + (c31[2*i]   - c32[2*i]) * s3 +
                   (c13[2*i-1] - c23[2*i-1])*c2 + (c13[2*i]   - c23[2*i]) * s2 +
                   (c32[2*i-1] - c31[2*i-1])*c1 + (c32[2*i]   - c31[2*i]) * s1 +
                   (c23[2*i-1] - c13[2*i-1]);

        Z[2*i]   = (c31[2*i]   - c32[2*i]) * c3 - (c31[2*i-1] - c32[2*i-1])*s3 +
                   (c13[2*i]   - c23[2*i]) * c2 - (c13[2*i-1] - c23[2*i-1])*s2 +
                   (c32[2*i]   - c31[2*i]) * c1 - (c32[2*i-1] - c31[2*i-1])*s1 +
                   (c23[2*i]   - c13[2*i]);

        XSL[2*i-1] = 2.0*fonfs*(X[2*i-1]*cSL-X[2*i]*sSL);
        XSL[2*i]   = 2.0*fonfs*(X[2*i-1]*sSL+X[2*i]*cSL);
        YSL[2*i-1] = 2.0*fonfs*(Y[2*i-1]*cSL-Y[2*i]*sSL);
        YSL[2*i]   = 2.0*fonfs*(Y[2*i-1]*sSL+Y[2*i]*cSL);
        ZSL[2*i-1] = 2.0*fonfs*(Z[2*i-1]*cSL-Z[2*i]*sSL);
        ZSL[2*i]   = 2.0*fonfs*(Z[2*i-1]*sSL+Z[2*i]*cSL);

        XLS[2*i-1] = (X[2*i-1]*cLS - X[2*i]*sLS);
        XLS[2*i]   = (X[2*i-1]*sLS + X[2*i]*cLS);
        YLS[2*i-1] = (Y[2*i-1]*cLS - Y[2*i]*sLS);
        YLS[2*i]   = (Y[2*i-1]*sLS + Y[2*i]*cLS);
        ZLS[2*i-1] = (Z[2*i-1]*cLS - Z[2*i]*sLS);
        ZLS[2*i]   = (Z[2*i-1]*sLS + Z[2*i]*cLS);
    }
}
