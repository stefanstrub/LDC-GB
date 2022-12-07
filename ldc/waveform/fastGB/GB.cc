#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#include "LISA.h"
#include "GB.h"
#include "c_wrapper.h"

void Fast_GB(double *params, long N, double Tobs, double dt,
	     double *XLS, double *YLS, double *ZLS,
	     double* XSL, double* YSL, double* ZSL, int NP)
{
  double* orbit_params = (double*) malloc(3*sizeof(double));
  orbit_params[0] = LARM;
  orbit_params[1] = LAMBDA;
  orbit_params[2] = KAPPA;

  Fast_GB_with_orbits(params, N, Tobs,dt, orbit_params, XLS, YLS, ZLS, XSL, YSL, ZSL, NP);
}

void Fast_GB_with_orbits(double *params, long N, double Tobs, double dt,
			 double *orbit_params,
			 double *XLS, double *YLS, double *ZLS,
			 double* XSL, double* YSL, double* ZSL, int NP)
{

	long n;     // iterator
	double t;	// time

	struct AnalyticOrbits* lisa = newAnalyticOrbits(orbit_params[0], orbit_params[1],
							orbit_params[2]);
        double Larm = AnalyticOrbits_get_armlength(lisa);
	double fstar = (C/Larm)/(2*PI); //0.01908538063694777;

	// waveform struct to hold pieces for calculation
	struct Waveform *wfm = (struct Waveform*) malloc(sizeof(struct Waveform));

	wfm->N  = N; // set number of samples
	wfm->T  = Tobs; // set observation period
	wfm->NP = NP; // inform model of number of parameters being used
	alloc_waveform(wfm);	 // allocate memory to hold pieces of waveform
	copy_params(wfm, params);    // copy parameters to waveform structure


	get_basis_tensors(wfm);      //  Tensor construction for building slowly evolving LISA response

	for(n=0; n<N; n++)
	{
	  t = wfm->T*(double)(n)/(double)N; // First time sample must be at t=0 for phasing
	  calc_xi_f(wfm, lisa, fstar, t);  // calc frequency and time variables
	  calc_sep_vecs(wfm, Larm);       // calculate the S/C separation vectors
	  calc_d_matrices(wfm);     // calculate pieces of waveform
	  calc_kdotr(wfm);		  // calculate dot product
	  get_transfer(wfm, t);     // Calculating Transfer function
	  fill_time_series(wfm, n); // Fill  time series data arrays with slowly evolving signal.
	}

	fft_data(wfm);     // Numerical Fourier transform of slowly evolving signal
	unpack_data(wfm);  // Unpack arrays from FFT and normalize


	XYZ(wfm->d, wfm->params[0]/wfm->T, wfm->q, N, dt, Tobs, Larm, fstar, XLS, YLS, ZLS, XSL, YSL, ZSL);

	free_waveform(wfm);  // Deallocate memory
	free(wfm);
	deleteAnalyticOrbits(lisa);

	return;
}



void calc_xi_f(struct Waveform *wfm, struct AnalyticOrbits* lisa, double fstar, double t)
{
	long i;

	double f0;
	double dfdt_0 = 0;
	double d2fdt2_0 = 0;

	f0       = wfm->params[0]/wfm->T;
	if (wfm->NP > 7) dfdt_0   = wfm->params[7]/wfm->T/wfm->T;
	if (wfm->NP > 8) d2fdt2_0 = wfm->params[8]/wfm->T/wfm->T/wfm->T;

	// Calculate position of each spacecraft at time t
	//spacecraft(t, wfm->x, wfm->y, wfm->z);
	for(i=0; i<3; i++){
	  wfm->x[i] = AnalyticOrbits_get_position_x(lisa, i+1, t);
	  wfm->y[i] = AnalyticOrbits_get_position_y(lisa, i+1, t);
	  wfm->z[i] = AnalyticOrbits_get_position_z(lisa, i+1, t);
	}

	for(i=0; i<3; i++)
	{
		wfm->kdotx[i] = (wfm->x[i]*wfm->k[0] + wfm->y[i]*wfm->k[1] + wfm->z[i]*wfm->k[2])/C;
		//Wave arrival time at spacecraft i
		wfm->xi[i]    = t - wfm->kdotx[i];
		//FIXME
		//wfm->xi[i]    = t + wfm->kdotx[i];
		//First order approximation to frequency at spacecraft i
		wfm->f[i]     = f0;
		if (wfm->NP > 7) wfm->f[i] += dfdt_0*wfm->xi[i];
		if (wfm->NP > 8) wfm->f[i] += 0.5*d2fdt2_0*wfm->xi[i]*wfm->xi[i];

		//Ratio of true frequency to transfer frequency
		wfm->fonfs[i] = wfm->f[i]/fstar;
	}

	return;
}

void copy_params(struct Waveform *wfm, double *params)
{
	long i;

	int NP = wfm->NP;

	wfm->params = (double*) malloc(NP*sizeof(double));

	for (i=0; i<NP; i++) wfm->params[i] = params[i];

	wfm->q  = (long)(params[0]); //Calculate carrier frequency bin

	return;
}

void fill_time_series(struct Waveform *wfm, int n)
{
	wfm->data12[2*n]   = wfm->TR[0][1];
	wfm->data21[2*n]   = wfm->TR[1][0];
	wfm->data31[2*n]   = wfm->TR[2][0];
	wfm->data12[2*n+1] = wfm->TI[0][1];
	wfm->data21[2*n+1] = wfm->TI[1][0];
	wfm->data31[2*n+1] = wfm->TI[2][0];
	wfm->data13[2*n]   = wfm->TR[0][2];
	wfm->data23[2*n]   = wfm->TR[1][2];
	wfm->data32[2*n]   = wfm->TR[2][1];
	wfm->data13[2*n+1] = wfm->TI[0][2];
	wfm->data23[2*n+1] = wfm->TI[1][2];
	wfm->data32[2*n+1] = wfm->TI[2][1];

	return;
}

void unpack_data(struct Waveform *wfm)
{
	long i;
	long N = wfm->N;

	for(i=0; i<N; i++)
	{	// populate from most negative (Nyquist) to most positive (Nyquist-1)
		wfm->a12[i]   = wfm->data12[N+i]/(double)N;
		wfm->a21[i]   = wfm->data21[N+i]/(double)N;
		wfm->a31[i]   = wfm->data31[N+i]/(double)N;
		wfm->a12[i+N] = wfm->data12[i]/(double)N;
		wfm->a21[i+N] = wfm->data21[i]/(double)N;
		wfm->a31[i+N] = wfm->data31[i]/(double)N;
		wfm->a13[i]   = wfm->data13[N+i]/(double)N;
		wfm->a23[i]   = wfm->data23[N+i]/(double)N;
		wfm->a32[i]   = wfm->data32[N+i]/(double)N;
		wfm->a13[i+N] = wfm->data13[i]/(double)N;
		wfm->a23[i+N] = wfm->data23[i]/(double)N;
		wfm->a32[i+N] = wfm->data32[i]/(double)N;
	}

	//   Renormalize so that the resulting time series is real
	for(i=0; i<2*N; i++)
	{
		wfm->d[0][1][i] = 0.5*wfm->a12[i];
		wfm->d[1][0][i] = 0.5*wfm->a21[i];
		wfm->d[2][0][i] = 0.5*wfm->a31[i];
		wfm->d[0][2][i] = 0.5*wfm->a13[i];
		wfm->d[1][2][i] = 0.5*wfm->a23[i];
		wfm->d[2][1][i] = 0.5*wfm->a32[i];
	}

	return;
}

void fft_data(struct Waveform *wfm)
{
	long N = wfm->N;
	gsl_fft_complex_radix2_forward(wfm->data12, 1, N);
	gsl_fft_complex_radix2_forward(wfm->data21, 1, N);
	gsl_fft_complex_radix2_forward(wfm->data31, 1, N);
	gsl_fft_complex_radix2_forward(wfm->data13, 1, N);
	gsl_fft_complex_radix2_forward(wfm->data23, 1, N);
	gsl_fft_complex_radix2_forward(wfm->data32, 1, N);

	return;
}

void alloc_waveform(struct Waveform *wfm)
{
	long i, j, n;
	long N;

	N = wfm->N;

	wfm->k = (double*) malloc(3*sizeof(double));

	wfm->kdotx = (double*) malloc(3*sizeof(double));
	wfm->kdotr = (double**) malloc(3*sizeof(double *));
	for (i=0; i<3; i++) wfm->kdotr[i] = (double*) malloc(3*sizeof(double));

	for (i=0; i<3; i++)
	{
	  for (j=0; j<3; j++) wfm->kdotr[i][j] = 0.;
	  wfm->kdotx[i] = 0.;
	}

	wfm->xi    = (double*) malloc(3*sizeof(double));
	wfm->f     = (double*) malloc(3*sizeof(double));
	wfm->fonfs = (double*) malloc(3*sizeof(double));
	for (i=0; i<3; i++)
	{
		wfm->xi[i]    = 0.;
		wfm->f[i]     = 0.;
		wfm->fonfs[i] = 0.;
	}

	// Polarization basis tensors
	wfm->eplus = (double**) malloc(3*sizeof(double*));
	for (i=0; i<3; i++) wfm->eplus[i] = (double*) malloc(3*sizeof(double));
	wfm->ecross = (double**) malloc(3*sizeof(double*));
	for (i=0; i<3; i++) wfm->ecross[i] = (double*) malloc(3*sizeof(double));

	wfm->dplus = (double**) malloc(3*sizeof(double*));
	for (i=0; i<3; i++) wfm->dplus[i] = (double*) malloc(3*sizeof(double));
	wfm->dcross = (double**) malloc(3*sizeof(double*));
	for (i=0; i<3; i++) wfm->dcross[i] = (double*) malloc(3*sizeof(double));

	wfm->r12 = (double*) malloc(3*sizeof(double));
	wfm->r21 = (double*) malloc(3*sizeof(double));
	wfm->r31 = (double*) malloc(3*sizeof(double));
	wfm->r13 = (double*) malloc(3*sizeof(double));
	wfm->r23 = (double*) malloc(3*sizeof(double));
	wfm->r32 = (double*) malloc(3*sizeof(double));


	wfm->data12 = (double*) malloc(2*N*sizeof(double));
	wfm->data21 = (double*) malloc(2*N*sizeof(double));
	wfm->data31 = (double*) malloc(2*N*sizeof(double));
	wfm->data13 = (double*) malloc(2*N*sizeof(double));
	wfm->data23 = (double*) malloc(2*N*sizeof(double));
	wfm->data32 = (double*) malloc(2*N*sizeof(double));
	for (i=0; i<2*N; i++)
	{
		wfm->data12[i] = 0.;
		wfm->data21[i] = 0.;
		wfm->data31[i] = 0.;
		wfm->data13[i] = 0.;
		wfm->data23[i] = 0.;
		wfm->data32[i] = 0.;
	}

	wfm->a12 = (double*) malloc(2*N*sizeof(double));
	wfm->a21 = (double*) malloc(2*N*sizeof(double));
	wfm->a31 = (double*) malloc(2*N*sizeof(double));
	wfm->a13 = (double*) malloc(2*N*sizeof(double));
	wfm->a23 = (double*) malloc(2*N*sizeof(double));
	wfm->a32 = (double*) malloc(2*N*sizeof(double));
	for (i=0; i<2*N; i++)
	{
		wfm->a12[i] = 0.;
		wfm->a21[i] = 0.;
		wfm->a31[i] = 0.;
		wfm->a13[i] = 0.;
		wfm->a23[i] = 0.;
		wfm->a32[i] = 0.;
	}

	wfm->TR = (double**) malloc(3*sizeof(double*));
	for (i=0;i<3; i++) wfm->TR[i] = (double*) malloc(3*sizeof(double));
	wfm->TI = (double**) malloc(3*sizeof(double*));
	for (i=0;i<3; i++) wfm->TI[i] = (double*) malloc(3*sizeof(double));

	wfm->x = (double*) malloc(3*sizeof(double));
	wfm->y = (double*) malloc(3*sizeof(double));
	wfm->z = (double*) malloc(3*sizeof(double));

	for (i=0; i<3; i++)
	{
		for(j=0; j<3; j++)
		{
			wfm->eplus[i][j]  = 0.;
			wfm->ecross[i][j] = 0.;
			wfm->dplus[i][j]  = 0.;
			wfm->dcross[i][j] = 0.;
			wfm->TR[i][j]     = 0.;
			wfm->TI[i][j]     = 0.;
		}
		wfm->x[i]   = 0.;
		wfm->y[i]   = 0.;
		wfm->z[i]   = 0.;
		wfm->r12[i] = 0.;
		wfm->r21[i] = 0.;
		wfm->r31[i] = 0.;
		wfm->r13[i] = 0.;
		wfm->r23[i] = 0.;
		wfm->r32[i] = 0.;
	}

	wfm->d = (double***) malloc(3*sizeof(double**));
	for (i=0; i<3; i++) wfm->d[i] = (double**) malloc(3*sizeof(double*));

	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
		{
			wfm->d[i][j] =(double*) malloc(2*N*sizeof(double));
		}
	}

	for (i=0; i<3; i++)
	{
		for(j=0; j<3; j++)
		{
			for (n=0; n<2*N; n++)
			{
				wfm->d[i][j][n] = 0.;
			}
		}
	}

	return;
}

void free_waveform(struct Waveform *wfm)
{
	long i, j;

	free(wfm->k);
	free(wfm->kdotx);
	for (i=0; i<3; i++) free(wfm->kdotr[i]);
	free(wfm->kdotr);

	free(wfm->xi);
	free(wfm->f);
	free(wfm->fonfs);

	for (i=0; i<3; i++){
	  free(wfm->eplus[i]);
	  free(wfm->ecross[i]);
	  free(wfm->dplus[i]);
	  free(wfm->dcross[i]);
	}
	free(wfm->eplus);
	free(wfm->ecross);
	free(wfm->dplus);
	free(wfm->dcross);

	free(wfm->r12);
	free(wfm->r21);
	free(wfm->r31);
	free(wfm->r13);
	free(wfm->r23);
	free(wfm->r32);

	free(wfm->data12);
	free(wfm->data21);
	free(wfm->data31);
	free(wfm->data13);
	free(wfm->data23);
	free(wfm->data32);

	free(wfm->a12);
	free(wfm->a21);
	free(wfm->a31);
	free(wfm->a13);
	free(wfm->a23);
	free(wfm->a32);

	free(wfm->x);
	free(wfm->y);
	free(wfm->z);

	for (i=0; i<3; i++){
	  free(wfm->TR[i]);
	  free(wfm->TI[i]);
	}
	free(wfm->TR);
	free(wfm->TI);

	for (i=0; i<3; i++)
	  for(j=0; j<3; j++)
	    free(wfm->d[i][j]);

	for (i=0; i<3; i++)
	  free(wfm->d[i]);
	free(wfm->d);

	free(wfm->params);

	return;
}

void calc_d_matrices(struct Waveform *wfm)
{
	long i, j;

	//Zero arrays to be summed
	wfm->dplus [0][1] = wfm->dplus [0][2] = wfm->dplus [1][0] = 0.;
	wfm->dplus [1][2] = wfm->dplus [2][0] = wfm->dplus [2][1] = 0.;
	wfm->dcross[0][1] = wfm->dcross[0][2] = wfm->dcross[1][0] = 0.;
	wfm->dcross[1][2] = wfm->dcross[2][0] = wfm->dcross[2][1] = 0.;

	//Convenient quantities d+ & dx
	for(i=0; i<3; i++)
	{
		for(j=0; j<3; j++)
		{
			wfm->dplus [0][1] += wfm->r12[i]*wfm->r12[j]*wfm->eplus[i][j];
			wfm->dcross[0][1] += wfm->r12[i]*wfm->r12[j]*wfm->ecross[i][j];
			wfm->dplus [1][2] += wfm->r23[i]*wfm->r23[j]*wfm->eplus[i][j];
			wfm->dcross[1][2] += wfm->r23[i]*wfm->r23[j]*wfm->ecross[i][j];
			wfm->dplus [0][2] += wfm->r13[i]*wfm->r13[j]*wfm->eplus[i][j];
			wfm->dcross[0][2] += wfm->r13[i]*wfm->r13[j]*wfm->ecross[i][j];
		}
	}
	//Makng use of symmetry
	wfm->dplus[1][0] = wfm->dplus[0][1];  wfm->dcross[1][0] = wfm->dcross[0][1];
	wfm->dplus[2][1] = wfm->dplus[1][2];  wfm->dcross[2][1] = wfm->dcross[1][2];
	wfm->dplus[2][0] = wfm->dplus[0][2];  wfm->dcross[2][0] = wfm->dcross[0][2];

	return;
}

void calc_sep_vecs(struct Waveform *wfm, double Larm)
{
	long i;

	//Unit separation vector from spacecrafts i to j
	wfm->r12[0] = (wfm->x[1] - wfm->x[0])/Larm;
	wfm->r13[0] = (wfm->x[2] - wfm->x[0])/Larm;
	wfm->r23[0] = (wfm->x[2] - wfm->x[1])/Larm;
	wfm->r12[1] = (wfm->y[1] - wfm->y[0])/Larm;
	wfm->r13[1] = (wfm->y[2] - wfm->y[0])/Larm;
	wfm->r23[1] = (wfm->y[2] - wfm->y[1])/Larm;
	wfm->r12[2] = (wfm->z[1] - wfm->z[0])/Larm;
	wfm->r13[2] = (wfm->z[2] - wfm->z[0])/Larm;
	wfm->r23[2] = (wfm->z[2] - wfm->z[1])/Larm;

	//Make use of symmetry
	for(i=0; i<3; i++)
	{
		wfm->r21[i] = -wfm->r12[i];
		wfm->r31[i] = -wfm->r13[i];
		wfm->r32[i] = -wfm->r23[i];
	}
	return;
}

void calc_kdotr(struct Waveform *wfm)
{
	long i;

	//Zero arrays to be summed
	wfm->kdotr[0][1] = wfm->kdotr[0][2] = wfm->kdotr[1][0] = 0.;
	wfm->kdotr[1][2] = wfm->kdotr[2][0] = wfm->kdotr[2][1] = 0.;

	for(i=0; i<3; i++)
	{
		wfm->kdotr[0][1] += wfm->k[i]*wfm->r12[i];
		wfm->kdotr[0][2] += wfm->k[i]*wfm->r13[i];
		wfm->kdotr[1][2] += wfm->k[i]*wfm->r23[i];
	}

	//Making use of antisymmetry
	wfm->kdotr[1][0] = -wfm->kdotr[0][1];
	wfm->kdotr[2][0] = -wfm->kdotr[0][2];
	wfm->kdotr[2][1] = -wfm->kdotr[1][2];

	return;
}

void set_const_trans(struct Waveform *wfm)
{
	double amp, cosiota;
	double Aplus, Across;
	double psi;
	double sinps, cosps;

	amp      = exp(wfm->params[3]);
	cosiota  = wfm->params[4];
	psi      = wfm->params[5];

	//Calculate GW polarization amplitudes
	Aplus  = amp*(1. + cosiota*cosiota);
	// Aplus  = -amp*(1. + cosiota*cosiota);
	Across = -2.0*amp*cosiota;
	//Across = 2.0*amp*cosiota;

	//Calculate cos and sin of polarization
	cosps = cos(2.*psi);
	sinps = sin(2.*psi);

	//Calculate constant pieces of transfer functions
	wfm->DPr    =  Aplus*cosps;
	wfm->DPi    = -Across*sinps;
	wfm->DCr    = -Aplus*sinps;
	wfm->DCi    = -Across*cosps;

	return;
}

void get_basis_tensors(struct Waveform *wfm)
{
	long i, j;

	double *u, *v;		  // GW basis vectors

	u = (double*) malloc(3*sizeof(double));
	v = (double*) malloc(3*sizeof(double));

	set_const_trans(wfm);  // set the constant pieces of transfer function

	get_basis_vecs(wfm->params, u, v, wfm->k); //Gravitational Wave source basis vectors

	//GW polarization basis tensors
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			//wfm->eplus[i][j]  = u[i]*u[j] - v[i]*v[j];
			wfm->eplus[i][j]  = v[i]*v[j] - u[i]*u[j];
			wfm->ecross[i][j] = u[i]*v[j] + v[i]*u[j];
			//wfm->ecross[i][j] = -u[i]*v[j] - v[i]*u[j];
		}
	}

	free(u);
	free(v);

	return;
}

void get_basis_vecs(double *params, double *u, double *v, double *k)
{
	long i;

	double phi;
	double costh, sinth, cosph, sinph;

	for (i=0; i<3; i++)
	{
		u[i] = 0.;
		v[i] = 0.;
		k[i] = 0.;
	}

	phi	  = params[2];
	costh = params[1];

	sinth = sqrt(1.0-costh*costh);

	cosph = cos(phi);
	sinph = sin(phi);

	u[0] =  costh*cosph;  u[1] =  costh*sinph;  u[2] = -sinth;
	v[0] =  sinph;        v[1] = -cosph;        v[2] =  0.;
	k[0] = -sinth*cosph;  k[1] = -sinth*sinph;  k[2] = -costh;

	return;
}

void get_transfer(struct Waveform *wfm, double t)
{
	long i, j;
	long q;

	double tran1r, tran1i;
	double tran2r, tran2i;
	double aevol;			// amplitude evolution factor
	double arg1, arg2, sinc;
	double f0, dfdt_0, d2fdt2_0;
	double df, phi0;

	f0       = wfm->params[0]/wfm->T;
	phi0     = wfm->params[6];

	if (wfm->NP > 7) dfdt_0   = wfm->params[7]/wfm->T/wfm->T;
 	if (wfm->NP > 8) d2fdt2_0 = wfm->params[8]/wfm->T/wfm->T/wfm->T;

	q  = wfm->q;
	df = PI*2*(((double)q)/wfm->T);

	for(i=0; i<3; i++)
	{
		for(j=0; j<3; j++)
		{
			if(i!=j)
			{
				//Argument of transfer function
				// FIXME
				//arg1 = 0.5*wfm->fonfs[i]*(1. - wfm->kdotr[i][j]);
				arg1 = 0.5*wfm->fonfs[i]*(1. + wfm->kdotr[i][j]);

				//Argument of complex exponentials
				arg2 =  PI*2*f0*wfm->xi[i] + phi0 - df*t;

				if (wfm->NP > 7) arg2 += PI*dfdt_0*wfm->xi[i]*wfm->xi[i];
				if (wfm->NP > 8) arg2 += PI*d2fdt2_0*wfm->xi[i]*wfm->xi[i]*wfm->xi[i]/3.0 ;

				//Transfer function
				sinc = 0.25*sin(arg1)/arg1;
				if(arg1==0)  sinc = 0.25;
				  

				//Evolution of amplitude
				aevol = 1.0;
				if (wfm->NP > 7) aevol += 0.66666666666666666666*dfdt_0/f0*wfm->xi[i];

				///Real and imaginary pieces of time series (no complex exponential)
				tran1r = aevol*(wfm->dplus[i][j]*wfm->DPr + wfm->dcross[i][j]*wfm->DCr);
				tran1i = aevol*(wfm->dplus[i][j]*wfm->DPi + wfm->dcross[i][j]*wfm->DCi);

				//Real and imaginry components of complex exponential
				tran2r = cos(arg1 + arg2);
				tran2i = sin(arg1 + arg2);

				//Real & Imaginary part of the slowly evolving signal
				wfm->TR[i][j] = sinc*(tran1r*tran2r - tran1i*tran2i);
				wfm->TI[i][j] = sinc*(tran1r*tran2i + tran1i*tran2r);
			}
		}
	}

	return;
}

void XYZ(double ***d, double f0, long q, long M, double dt, double Tobs, double Larm, double fstar,
	 double *XLS, double *YLS, double *ZLS,
	 double* XSL, double* YSL, double* ZSL)
{
	int i;
	double fonfs;
	double c3, s3, c2, s2, c1, s1;
	double f;
	double *X, *Y, *Z;
	double phiLS, cLS, sLS, phiSL, cSL, sSL;

	X   = (double*) malloc(2*M*sizeof(double));
	Y   = (double*) malloc(2*M*sizeof(double));
	Z   = (double*) malloc(2*M*sizeof(double));

	// YLS = malloc(2*M*sizeof(double));
	// ZLS = malloc(2*M*sizeof(double));

	phiLS = 2*PI*f0*(dt/2.0-Larm/C);

	cLS = cos(phiLS);
	sLS = sin(phiLS);

	//double phiLS = 2.0*pi*f0*(dt/2.0-L/clight);
	//double cLS = cos(phiLS); double sLS = sin(phiLS);

	phiSL = PI/2.0-2.0*PI*f0*(Larm/C);
	cSL = cos(phiSL);
	sSL = sin(phiSL);

  //printf("Stas, q=%ld, f0=%f, check: %f, %f \n", q, f0, q/Tobs, Tobs);
	for(i=0; i<M; i++)
	{
		f = ((double)(q + i - M/2))/Tobs;
		fonfs = f/fstar;
		//printf("Stas fonfs = %f, %f, %f, %f \n", fonfs, f, fstar, Tobs);
		c3 = cos(3.*fonfs);  c2 = cos(2.*fonfs);  c1 = cos(1.*fonfs);
		s3 = sin(3.*fonfs);  s2 = sin(2.*fonfs);  s1 = sin(1.*fonfs);

		X[2*i]   = (d[0][1][2*i]-d[0][2][2*i])*c3 + (d[0][1][2*i+1]-d[0][2][2*i+1])*s3 +
		           (d[1][0][2*i]-d[2][0][2*i])*c2 + (d[1][0][2*i+1]-d[2][0][2*i+1])*s2 +
		           (d[0][2][2*i]-d[0][1][2*i])*c1 + (d[0][2][2*i+1]-d[0][1][2*i+1])*s1 +
		           (d[2][0][2*i]-d[1][0][2*i]);

		X[2*i+1] = (d[0][1][2*i+1]-d[0][2][2*i+1])*c3 - (d[0][1][2*i]-d[0][2][2*i])*s3 +
		           (d[1][0][2*i+1]-d[2][0][2*i+1])*c2 - (d[1][0][2*i]-d[2][0][2*i])*s2 +
		           (d[0][2][2*i+1]-d[0][1][2*i+1])*c1 - (d[0][2][2*i]-d[0][1][2*i])*s1 +
		           (d[2][0][2*i+1]-d[1][0][2*i+1]);

		Y[2*i]   = (d[1][2][2*i]-d[1][0][2*i])*c3 + (d[1][2][2*i+1]-d[1][0][2*i+1])*s3 +
		           (d[2][1][2*i]-d[0][1][2*i])*c2 + (d[2][1][2*i+1]-d[0][1][2*i+1])*s2+
		           (d[1][0][2*i]-d[1][2][2*i])*c1 + (d[1][0][2*i+1]-d[1][2][2*i+1])*s1+
		           (d[0][1][2*i]-d[2][1][2*i]);

		Y[2*i+1] = (d[1][2][2*i+1]-d[1][0][2*i+1])*c3 - (d[1][2][2*i]-d[1][0][2*i])*s3+
		           (d[2][1][2*i+1]-d[0][1][2*i+1])*c2 - (d[2][1][2*i]-d[0][1][2*i])*s2+
		           (d[1][0][2*i+1]-d[1][2][2*i+1])*c1 - (d[1][0][2*i]-d[1][2][2*i])*s1+
		           (d[0][1][2*i+1]-d[2][1][2*i+1]);

		Z[2*i]   = (d[2][0][2*i]-d[2][1][2*i])*c3 + (d[2][0][2*i+1]-d[2][1][2*i+1])*s3+
		           (d[0][2][2*i]-d[1][2][2*i])*c2 + (d[0][2][2*i+1]-d[1][2][2*i+1])*s2+
		           (d[2][1][2*i]-d[2][0][2*i])*c1 + (d[2][1][2*i+1]-d[2][0][2*i+1])*s1+
		           (d[1][2][2*i]-d[0][2][2*i]);

		Z[2*i+1] = (d[2][0][2*i+1]-d[2][1][2*i+1])*c3 - (d[2][0][2*i]-d[2][1][2*i])*s3+
		           (d[0][2][2*i+1]-d[1][2][2*i+1])*c2 - (d[0][2][2*i]-d[1][2][2*i])*s2+
		           (d[2][1][2*i+1]-d[2][0][2*i+1])*c1 - (d[2][1][2*i]-d[2][0][2*i])*s1+
		           (d[1][2][2*i+1]-d[0][2][2*i+1]);

		// XLS[2*i]   =  (X[2*i]*cLS - X[2*i+1]*sLS);
		// XLS[2*i+1] = -(X[2*i]*sLS + X[2*i+1]*cLS);
		// YLS[2*i]   =  (Y[2*i]*cLS - Y[2*i+1]*sLS);
		// YLS[2*i+1] = -(Y[2*i]*sLS + Y[2*i+1]*cLS);
		// ZLS[2*i]   =  (Z[2*i]*cLS - Z[2*i+1]*sLS);
		// ZLS[2*i+1] = -(Z[2*i]*sLS + Z[2*i+1]*cLS);
    //
		// XSL[2*i]   =  2.0*fonfs*(X[2*i]*cSL - X[2*i+1]*sSL);
		// XSL[2*i+1] = -2.0*fonfs*(X[2*i]*sSL + X[2*i+1]*cSL);
		// YSL[2*i]   =  2.0*fonfs*(Y[2*i]*cSL - Y[2*i+1]*sSL);
		// YSL[2*i+1] = -2.0*fonfs*(Y[2*i]*sSL + Y[2*i+1]*cSL);
		// ZSL[2*i]   =  2.0*fonfs*(Z[2*i]*cSL - Z[2*i+1]*sSL);
		// ZSL[2*i+1] = -2.0*fonfs*(Z[2*i]*sSL + Z[2*i+1]*cSL);

		// Alternative polarization definition
		XLS[2*i]   =  (X[2*i]*cLS - X[2*i+1]*sLS);
		XLS[2*i+1] =  (X[2*i]*sLS + X[2*i+1]*cLS);
		YLS[2*i]   =  (Y[2*i]*cLS - Y[2*i+1]*sLS);
		YLS[2*i+1] =  (Y[2*i]*sLS + Y[2*i+1]*cLS);
		ZLS[2*i]   =  (Z[2*i]*cLS - Z[2*i+1]*sLS);
		ZLS[2*i+1] =  (Z[2*i]*sLS + Z[2*i+1]*cLS);

		XSL[2*i]   =  2.0*fonfs*(X[2*i]*cSL - X[2*i+1]*sSL);
		XSL[2*i+1] =  2.0*fonfs*(X[2*i]*sSL + X[2*i+1]*cSL);
		YSL[2*i]   =  2.0*fonfs*(Y[2*i]*cSL - Y[2*i+1]*sSL);
		YSL[2*i+1] =  2.0*fonfs*(Y[2*i]*sSL + Y[2*i+1]*cSL);
		ZSL[2*i]   =  2.0*fonfs*(Z[2*i]*cSL - Z[2*i+1]*sSL);
		ZSL[2*i+1] =  2.0*fonfs*(Z[2*i]*sSL + Z[2*i+1]*cSL);

	}

	// for(i=0; i<2*M; i++)
	// {
	// 	// A channel
	// 	ALS[i] = (2.0*XLS[i] - YLS[i] - ZLS[i])/3.0;
	// 	// E channel
	// 	ELS[i] = (ZLS[i]-YLS[i])/SQ3;
	// }

	free(X);
	free(Y);
	free(Z);

	//free(YLS);
	//free(ZLS);

	return;
}

long get_N(double *params, double Tobs, double fstar)
{
	// This determines the number of samples to take of the slowly evolving bit
	// of the GB waveform. Right now only instrument noise is used in the estimate

	long mult, N, M;

	double amp, f0, fonfs;
	double SnX, SnAE;
	double Acut, Sm;

	f0  = params[0]/Tobs;
	amp = exp(params[3]);

	mult = 8;
	if((Tobs/YEAR) <= 8.0) mult = 8;
	if((Tobs/YEAR) <= 4.0) mult = 4;
	if((Tobs/YEAR) <= 2.0) mult = 2;
	if((Tobs/YEAR) <= 1.0) mult = 1;

	N = 32*mult;
	if(f0 > 0.001) N = 64*mult;
	if(f0 > 0.01)  N = 256*mult;
	if(f0 > 0.03)  N = 512*mult;
	if(f0 > 0.1)   N = 1024*mult;

	fonfs = f0/fstar;

	instrument_noise(f0, &SnAE, &SnX);

	//  calculate michelson noise
	Sm = SnX/(4.0*sin(fonfs)*sin(fonfs));

	Acut = amp*sqrt(Tobs/Sm);

	M = (long)(pow(2.0,(rint(log(Acut)/log(2.0))+1.0)));

	if(M < N)    M = N;
	if(N < M)    N = M;
	if(M > 8192) M = 8192;

	N = M;

	return N;
}
