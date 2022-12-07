#ifndef FREQAKRAHH
#define FREQAKRAHH

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include "EMRItemplate.h"
#include "orbits.hpp"
#include "lisaconstants.hpp" // LISA constants

/*** Class FreqAK_RA
 *
 * Generator of analytic kludge waveforms described in
 * Barack & Cutler gr-qc/0310125. it generates waveform in the freq. domain using rigid adiabatic approximation
 *
 * @author  Jon Gair, Stas Babak 2010
 */

//namespace LISAWP{
//
void EccentricLISAMotion(float kappa0, float lambda0, double t, double* R, double** q, double** n);
void EccentricLISAMotion2(AnalyticOrbits lisa, double t, double** q, double** n);

class FreqAK_RA{

  public:

    FreqAK_RA();

   /*** Constructor
    * @param timestep cadence for the time domain
    * @param maxDur maximum duration of the observation, stop waveform if plunge didn't occure
    * @param dtPhase cadence (internal) for evolving orbit and LISA motion, (usually 2048 sec).
    */
    FreqAK_RA(double timestep, double maxDur, double dtPhase, int Extended, double* orbit_params);

    /***  destructor */
    ~FreqAK_RA();

    /*** Evolves orbit and responce (RA) with large time step dtPhase
     * Note! that integration starting from any point is not yet implemented it is hardcoded
     * that template starts at 0 at up to maxDuration long....
     * @param  S  EMRItemplate object with parameters etc.
     * @param t0  moment of time where IC are specified
     * @param dur duration of the template
     * @param tStart time beginning of template: tStart <= t0 <= tStart+dur
     */

    void PhaseEv_RA(EMRItemplate& S, double t0, double dur, double tStart);
    void NewPhaseEv_RA(EMRItemplate& S, double t_ini, double dur, double tStart);

    // double ComputeHAmpl(double ecc, int n);
    void ComputeHAmpl(double ecc, int n, double* XaXb, double* Xa_Xb, double* Xc);

    void GetEvolution(double* t_phase, double* nus, double* eccs, double* phis, double* alps, double* gams, double* fhigh);

    /** The same as above but for long-wavelength limit response */
    //void PhaseEv_LW(EMRItemplate& S, double t0, double dur, double tStart);

    /** Constructing waveform in time domain using LW response */
    //void ConstructWaveTime_LW(EMRItemplate& S, Matrix<double>& X, Matrix<double>& Y, Matrix<double>& Z);

    /** Constructing waveform in time domain using RA response */
    //void ConstructWaveTime_RA(EMRItemplate& S, Matrix<double>& X, Matrix<double>& Y, Matrix<double>& Z);

    /**   Constructing waveform in time domain using LW response */
    //void ConstructWaveTime_LW(EMRItemplate& S, Matrix<double>& X, Matrix<double>& Z_Y);

   /** Constructing waveform in freq domain using RA response */
    void ConstructWaveFreq_RA(EMRItemplate& S, double* Xf_r, double* Xf_im,  double* Yf_r, double* Yf_im, double* Zf_r, double* Zf_im);
    /**
    @param useBessel if >0  use bessel functions for amplitude evaluation (requires gsl), otherwise use e<<1 Taylor expansion
    **/

    /** Returns harmonics and phases separately */
    // void GetHarmonicsXYZ(EMRItemplate& S, double** X_r, double** X_i, double** Y_r, double** Y_i, double** Z_r, double** Z_i, double* phi, double* gam, double* alph, double* tim, double* frqs);

    void GetHarmonicsXYZ(EMRItemplate& S,  double** X_r, double** X_i, double** Y_r, double** Y_i, double** Z_r,
                double** Z_i, double* phi, double* gam, double* alph, double* tim, double** frqs, double** phase);

   /** Constructing waveform in freq domain using LW response */
    //void ConstructWaveFreq_LW(EMRItemplate& S, Matrix<double>& Xf, Matrix<double>& Z_Yf);

    void ComputeResponse(double om, double* kn, double* kq, double* u, double *v, double c2psi, double s2psi,
    		std::complex<double>* Xplus, std::complex<double>* Xcross, std::complex<double>* Yplus,
            std::complex<double>* Ycross, std::complex<double>* Zplus, std::complex<double>* Zcross);

    void ComputeOrbPrecAmpls(EMRItemplate& S, std::complex<double>* Ap, std::complex<double>* Ac,
    						std::complex<double>*Bp, std::complex<double>* Bc );

    void InterpolateWaveform(int nn, int ll, int mm, int ci, EMRItemplate& S, double xi_in, double xi_fin,
                            std::complex<double>* Ap, std::complex<double>* Ac, double* Xf_r, double* Xf_im,
                            double* Yf_r, double* Yf_im, double* Zf_r, double* Zf_im);

    void ConstructWaveFreq_RAv2(EMRItemplate& S, double* Xf_r, double* Xf_im,  double* Yf_r, double* Yf_im,
                                                    double* Zf_r, double* Zf_im);

    void FreeMemory();

  private:

     double* tm_ph;
     double* phiev;
     double* alpev;
     double* gamev;
     double* eccev;
     double* nuev;
     double* gamdotev;
     double* alpdotev;
     double* phiddotev;
     double* gamddotev;
     double* alpddotev;
     double* Ampev;
     int Ext;

     std::complex<double>** Xp;
     std::complex<double>** Xc;
     std::complex<double>** Yp;
     std::complex<double>** Yc;
     std::complex<double>** Zp;
     std::complex<double>** Zc;

     std::complex<double>** XXp;
     std::complex<double>** XXc;
     std::complex<double>** YYp;
     std::complex<double>** YYc;
     std::complex<double>** ZZp;
     std::complex<double>** ZZc;

     int Nph;     // number of points for phase and response evaluation
     int Nps;    // number of points for waveform evaluation

     double dt_w;
     double dt_ph;
     double df;    // freq resolution df = 1/Tobs
     double Tobs;  // duration of observation in sec
     int imax;  // maximum index of the orbital evolution
     double arm; // armlength
     AnalyticOrbits lisa_orbits;
     //double ArcTan(double up, double down);

     double Sinc(double y);


};
//} // end of the namespace

#endif
