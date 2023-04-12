/****

@author  Jon Gair, Stas Babak 2010

***/

#include "FreqAK_RAv2.h"
#include <gsl/gsl_sf_bessel.h>



void EccentricLISAMotion2(AnalyticOrbits lisa, double t, double** q, double** n){

  for(int i=0; i<3; i++){ // loop over s/c
    lisa.position_x(i+1, &t, &q[i][0], 1); 
    lisa.position_y(i+1, &t, &q[i][1], 1);
    lisa.position_z(i+1, &t, &q[i][2], 1);

    double alpha = lisa.alpha(t);
    q[i][0] -= LISAWP_AU_SI*cos(alpha);
    q[i][1] -= LISAWP_AU_SI*sin(alpha);
    q[i][0] /= lisa.spacecraft_separation;
    q[i][1] /= lisa.spacecraft_separation;
    q[i][2] /= lisa.spacecraft_separation;
  }
  
  for(int i=0; i<3; i++){ // links: 1st index - s/c, 2nd index - coordinates
    n[0][i] = q[1][i] - q[2][i];
    n[1][i] = q[2][i] - q[0][i];
    n[2][i] = q[0][i] - q[1][i];
  }
}

void EccentricLISAMotion(float kappa0, float lambda0, double t, double* R, \
	                         double** q, double** n){
	//	double Omega = LISAWP_TWOPI/year;

  double Omega = LISAWP_TWOPI/LISAWP_YRSID_SI; //*3.1709791983764586e-8;
  double alpha = Omega*t + kappa0;
  double xi[3];
  xi[0] = lambda0;             //////////  MOVE THOSE CALCULATIONS OUT TO THE CONSTRUCTOR
  xi[1] = lambda0 + 2.0*LISAWP_PI/3.0;
  xi[2] = lambda0 + 4.0*LISAWP_PI/3.0;
  double RAU = LISAWP_AU_SI/LISAWP_C_SI;
  
  R[0] = RAU*cos(alpha);
  R[1] = RAU*sin(alpha);
  R[2] = 0;
  
  for(int i=0; i<3; i++){ // loop over s/c
    q[i][0] = ( sin(alpha)*cos(alpha)*sin(xi[i]) -			\
		(1.+ sin(alpha)*sin(alpha))*cos(xi[i])  )/(2.0*sqrt(3.));
    q[i][1] = ( sin(alpha)*cos(alpha)*cos(xi[i]) -			\
		(1.+ cos(alpha)*cos(alpha))*sin(xi[i])  )/(2.0*sqrt(3.));
    q[i][2] = -0.5*cos(alpha - xi[i]);
  }
  for(int i=0; i<3; i++){ // links: 1st index - s/c, 2nd index - coordinates
    n[0][i] = q[1][i] - q[2][i];
    n[1][i] = q[2][i] - q[0][i];
    n[2][i] = q[0][i] - q[1][i];
  }
}


FreqAK_RA::FreqAK_RA(){
		// std::cout << "null" << std::endl;
}

FreqAK_RA::FreqAK_RA(double timestep, double maxDur, double dtPhase, int Extended,
		     double *orbit_params){
  
  lisa_orbits = AnalyticOrbits(orbit_params[0], orbit_params[1],
			       orbit_params[2]);
  double LISA_arm = lisa_orbits.spacecraft_separation;
  
  dt_w = timestep;
  Tobs = maxDur;
  Nps = (int)floor(Tobs/dt_w);
  df = 1./Tobs;
  dt_ph = dtPhase;
  Nph = (int)floor(Tobs/dt_ph);
  arm =  LISA_arm/ LISAWP_C_SI; // LISA's arm in sec
  Ext = Extended;

    phiev = NULL;
    alpev = NULL;
    gamev = NULL;
    eccev = NULL;
    nuev = NULL;
    gamdotev =  NULL;
    alpdotev = NULL;
    phiddotev = NULL;
    gamddotev = NULL;
    alpddotev = NULL;
    Ampev = NULL;
    tm_ph = NULL;

    Xp = NULL;
    Xc = NULL;
    Yp = NULL;
    Yc = NULL;
    Zp = NULL;
    Zc = NULL;
    XXp = NULL;
    XXc = NULL;
    YYp = NULL;
    YYc = NULL;
    ZZp = NULL;
    ZZc = NULL;
    
    phiev = new double[Nph];
    tm_ph = new double[Nph];
    alpev = new double[Nph];
    gamev = new double[Nph];
    eccev = new double[Nph];
    nuev = new double[Nph];
    gamdotev = new double[Nph];
    alpdotev = new double[Nph];
    phiddotev = new double[Nph];
    gamddotev = new double[Nph];
    alpddotev = new double[Nph];
    Ampev = new double[Nph];


    Xp = new std::complex<double>*[25];
    Xc = new std::complex<double>*[25];
    Yp = new std::complex<double>*[25];
    Yc = new std::complex<double>*[25];
    Zp = new std::complex<double>*[25];
    Zc = new std::complex<double>*[25];

	if (Extended>0){

		XXp = new std::complex<double>*[25];
	    XXc = new std::complex<double>*[25];
	    YYp = new std::complex<double>*[25];
	    YYc = new std::complex<double>*[25];
	    ZZp = new std::complex<double>*[25];
	    ZZc = new std::complex<double>*[25];

	}
		//std::cout << "here 3" << std::endl;
    for (int i=0; i<25; i++){
			Xp[i] = new std::complex<double>[Nph];
            Xc[i] = new std::complex<double>[Nph];
            Yp[i] = new std::complex<double>[Nph];
            Yc[i] = new std::complex<double>[Nph];
            Zp[i] = new std::complex<double>[Nph];
            Zc[i] = new std::complex<double>[Nph];
			if (Extended >0){
				XXp[i] = new std::complex<double>[Nph];
	            XXc[i] = new std::complex<double>[Nph];
	            YYp[i] = new std::complex<double>[Nph];
	            YYc[i] = new std::complex<double>[Nph];
	            ZZp[i] = new std::complex<double>[Nph];
	            ZZc[i] = new std::complex<double>[Nph];

			}
    }
			// std::cout << "here 4" << std::endl;
		for (int i=0; i<Nph; i++){
					phiev[i] = 0.0;
					tm_ph[i] = (double)i*dt_ph;
					alpev[i] = 0.0;
					gamev[i] = 0.0;
					eccev[i] = 0.0;
					nuev[i] = 0.0;
					gamdotev[i] = 0.0;
					alpdotev[i] = 0.0;
					phiddotev[i] = 0.0;
					gamddotev[i] = 0.0;
					alpddotev[i] = 0.0;
					Ampev[i] = 0.0;
		}
		std::complex<double> zr(0.0, 0.0);
			for (int i=0; i<25; i++){
					for (int ii=0; ii<Nph; ii++){
						Xp[i][ii] = zr;
            Xc[i][ii] = zr;
            Yp[i][ii] = zr;
            Yc[i][ii] = zr;
            Zp[i][ii] = zr;
            Zc[i][ii] = zr;
					}
      }


}

FreqAK_RA::~FreqAK_RA(){
	// std::cout << "destructor" << std::endl;
}

void FreqAK_RA::FreeMemory(){

// 	std::cout << "here d1" <<std::endl;
       if (phiev != NULL){
		//   std::cout << "here d1.2" <<std::endl;
          delete [] phiev;
			phiev=NULL;
		}
// 	  std::cout << "here d1.5" <<std::endl;
		if (alpev != NULL){
            delete [] alpev;
			alpev = NULL;
		}
		if	(gamev != NULL){
                delete gamev;
				gamev = NULL;
		}
		if (eccev != NULL){
		   delete [] eccev;
			 eccev = NULL;
		 }
		if (nuev != NULL){
		   delete [] nuev;
			 nuev = NULL;
		 }
		if (gamdotev !=  NULL){
		   delete [] gamdotev;
			 gamdotev = NULL;
		 }
		if (alpdotev != NULL){
		   delete [] alpdotev;
			 alpdotev = NULL;
		 }
		if (phiddotev != NULL){
		   delete [] phiddotev;
			 phiddotev = NULL;
		 }
		if (gamddotev != NULL){
		   delete [] gamddotev;
			 gamddotev = NULL;
		 }
		if (alpddotev != NULL){
		   delete [] alpddotev;
			 alpddotev = NULL;
		 }
		if (Ampev != NULL){
		   delete [] Ampev;
			 Ampev = NULL;
		 }
		 if (tm_ph != NULL){
		   delete [] tm_ph;
			 tm_ph = NULL;
		 }

		// std::cout << "here d2" <<std::endl;

		if (Xp != NULL){
		  for (int i=0; i<25; i++){
			     delete [] Xp[i];
			     delete [] Xc[i];
			     delete [] Yp[i];
			     delete [] Yc[i];
			     delete [] Zp[i];
			     delete [] Zc[i];
		  }
				  delete [] Xp;
				  delete [] Xc;
				  delete [] Yp;
				  delete [] Yc;
				  delete [] Zp;
				  delete [] Zc;
		}
		if (XXp != NULL){
 		 for (int i=0; i<25; i++){
 				delete [] XXp[i];
 				delete [] XXc[i];
 				delete [] YYp[i];
 				delete [] YYc[i];
 				delete [] ZZp[i];
 				delete [] ZZc[i];
 		 }
 				 delete [] XXp;
 				 delete [] XXc;
 				 delete [] YYp;
 				 delete [] YYc;
 				 delete [] ZZp;
 				 delete [] ZZc;
 	   }
// 	std::cout << "here d3" <<std::endl;
}


void FreqAK_RA::PhaseEv_RA(EMRItemplate& S, double t_ini, double dur, double tStart){
	/**
	@param S EMRItemplate object storing the parameters
	@param t0 moment of time for which intitial conditions are given
	@param maximum dur durarion of the signal
	@param tStart - how much back in time we should go.
	**/

	double k[3];
	double uhat[3];
	double vhat[3];
	double u[3];
	double v[3];
	double kn[3];
	double kq[3];

   double up = (S.ctS*S.stK*cos(S.phS - S.phK) - S.ctK*S.stS);
   double dw = (S.stK*sin(S.phS-S.phK));
	//  std::cout << up << "  " << dw << "  " << S.ctS << "  " << S.stS << "   " <<  S.phS << std::endl;
   double psi;
   if (dw != 0.0) {
      // psi = atan2(up, dw);
      // psi = 0.5*LISAWP_PI -atan2(up, dw);
      psi = -atan2(up, dw);
      // psi = atan2(up, dw);
			// std::cout << "psi = " << psi << std::endl;
   }else {
      psi = 0.5*LISAWP_PI;
   }

   double c2psi=cos(2.*psi);
   double s2psi=sin(2.*psi);


   // note that k = -n, where k is propagation vector and n is sky location
   k[0] =  -S.stS*S.cpS;
   k[1] = -S.stS*S.spS;
   k[2] = -S.ctS;

   uhat[0] = S.ctS*S.cpS;
   uhat[1] = S.ctS*S.spS;
   uhat[2] = -S.stS;

   vhat[0] = S.spS;
   vhat[1] = -S.cpS;
   vhat[2] = 0.0;
   double nU, nV;

   S.t0 = 0.0;
   // double e0 = S.e0;
   // double nu0 = S.nu0;
   // double ph0 = S.Phi0;
   // double al0 =  S.alpha0;
   // double gam0 = S.gamma0;

   double* R;
   double** q;
   double** n;
   R = new double[3];
   q = new double*[3];
   n = new double*[3];
   for (int i=0; i<3; i++){
       q[i] = new double[3];
       n[i] = new double[3];
   }

   double AUsec =  LISAWP_AU_SI/LISAWP_C_SI;
   double clam=cos(S.lam);
   double slam=sin(S.lam);



   double M = S.Mt;
   double mu = S.mt;
   double e = S.e0;
   double e2=e*e;
   double Sp = S.a;
   double nu = S.nu0;
   double phi = S.Phi0;
   double gam = S.gamma0;
   double alp = S.alpha0;
   double Y, Z;
   double edotm, nudotm, phidotm, alpdotm, gamdotm;
   double edot, nudot, phidot, alpdot, gamdot;
   double dalpdnu, dalpde, dgamdnu, dgamde, alpddot, gamddot, phiddot;
   double de, dnu, dphi, dgam, dalp, rhs;
   double T=0.0;
   int imax, ind;
   double ampfct= S.Ampl;

   double nn, mm, om;
   std::complex<double> chi0;
   std::complex<double> chi1;
   std::complex<double> chi2;
   std::complex<double> img(0.0, 1.0);
   double x, x2;

   // given the courseneess of the steps we need first to interpolate to the nearest dt*i
   int i_ini = floor(t_ini/dt_ph);
   double t_i_ini = i_ini*dt_ph;
   double dt_ini = t_i_ini - t_ini;

   Y=1./(1.-e2);
   Z=pow(LISAWP_TWOPI*M*nu,1./3.);
   edot = e*mu/M/M*(-1./15.*pow(Y,3.5)*pow(Z,8.)*((304.+121.*e2)/Y+Z*Z*(70648.-231960.*e2-56101.*e2*e2)/56.)+\
			Sp*clam*pow(Z,11.)*pow(Y,4.)*(8184.+10064.*e2+789.*e2*e2)/30.);
   nudot = 96./(10.*M_PI)*mu/pow(M,3)*(pow(Z,11.)*pow(Y,4.5)*((96.+292.*e2+37.*e2*e2)/Y/96.\
				+Z*Z*(20368.-61464.*e2-163170.*e2*e2-13147.*e2*e2*e2)/5376.) - \
			   pow(Z,14.)*pow(Y,5.)*Sp*clam*(1168.+9688.*e2+6286.*e2*e2 +195.*e2*e2*e2)/192.);
   phidot = 2.*LISAWP_PI*nu;
   alpdot = 8.*LISAWP_PI*LISAWP_PI*nu*nu*Sp*M*pow(Y,1.5);
   gamdot = 6.*LISAWP_PI*nu*Z*Z*Y*(1.+.25*Z*Z*Y*(26.-15.*e2))-3.*clam*alpdot;

	// change the initial conditions to t_i_ini.
   S.e0 =  e+ edot*dt_ini;
   S.Phi0 = phi + phidot*dt_ini;
   S.gamma0 =  gam + gamdot*dt_ini;
   S.alpha0 = alp + alpdot*dt_ini;
   S.nu0 = nu + nudot*dt_ini;
   S.t0 = t_i_ini;

   T = t_ini + dt_ini;

   // calculate how much we need forward and how much back

   double t_back = t_i_ini - tStart;
   if (t_back < 0.0){
	   std::cout << "tstart before t_ini  " << tStart << "  " << t_ini << "  " << t_i_ini << std::endl;
	   exit(1);
   }

   double t_forward = dur + tStart - t_ini;
   if (t_forward < 0.0){
	   std::cout << "duration is not compatible with t_ini and tStart  " << dur << "   " << tStart << "   " << t_ini << std::endl;
   }
   double t_end = dur + tStart; // end time

   if (dur > Tobs){
	   std::cout << "duration is larger than observation time " << dur << "   " << Tobs << std::endl;
   }

   // start forward loop

   // std::cout << "t_ini = " << t_ini << "  t_i_ini = " << t_i_ini << "  T = " << T << " tStart = " << tStart << \
   // 			"  t_end = " << t_end << "  dur = " << dur << "  " << t_i_ini - dur << std::endl;

	e = S.e0;
    e2=e*e;
    Sp = S.a;
    nu = S.nu0;
    phi = S.Phi0;
    gam = S.gamma0;
    alp = S.alpha0;

   int i = 0;
   while (T <= t_end){
	   e2=e*e;
	   Y=1./(1.-e2);
	   Z=pow(LISAWP_TWOPI*M*nu,1./3.);
	   if (i != 0){
		  edotm=edot;
		  nudotm=nudot;
		  phidotm=phidot;
		  alpdotm=alpdot;
		  gamdotm=gamdot;
	   }
	   edot = e*mu/M/M*(-1./15.*pow(Y,3.5)*pow(Z,8.)*((304.+121.*e2)/Y+Z*Z*(70648.-231960.*e2-56101.*e2*e2)/56.)+\
				Sp*clam*pow(Z,11.)*pow(Y,4.)*(8184.+10064.*e2+789.*e2*e2)/30.);
	   nudot = 96./(10.*M_PI)*mu/pow(M,3)*(pow(Z,11.)*pow(Y,4.5)*((96.+292.*e2+37.*e2*e2)/Y/96.\
					+Z*Z*(20368.-61464.*e2-163170.*e2*e2-13147.*e2*e2*e2)/5376.) - \
				   pow(Z,14.)*pow(Y,5.)*Sp*clam*(1168.+9688.*e2+6286.*e2*e2 +195.*e2*e2*e2)/192.);
	   phidot = 2.*LISAWP_PI*nu;
	   alpdot = 8.*LISAWP_PI*LISAWP_PI*nu*nu*Sp*M*pow(Y,1.5);
	   gamdot = 6.*LISAWP_PI*nu*Z*Z*Y*(1.+.25*Z*Z*Y*(26.-15.*e2))-3.*clam*alpdot;

	   dalpdnu=16.*Sp*LISAWP_PI*nu*M*pow(Y,1.5);
	   dalpde=12.*Sp*LISAWP_PI*nu*nu*M*sqrt(Y)*2.*e*Y*Y;
	   dgamdnu=6.*Z*Z*Y*(1.+.25*Z*Z*Y*(26.-15.*e2))-3.*clam*dalpdnu+12.*Z*Y*(1.+.5*Z*Z*Y*(25.-15.*e2))*Z/(3.);
	   dgamde=(6.*nu*Z*Z*(1.+.5*Z*Z*Y*(26.-15.*e2)))*2.*e*Y*Y-45.*nu*Z*Z*Z*Z*Y*Y*e-3.*clam*dalpde;

	   alpddot=LISAWP_PI*(dalpdnu*nudot+dalpde*edot);
	   gamddot=LISAWP_PI*(dgamdnu*nudot+dgamde*edot);
	   phiddot=LISAWP_TWOPI*nudot;

	   if (i == 0) {
		  edotm=edot;
		  nudotm=nudot;
		  alpdotm=alpdot;
		  gamdotm=gamdot;
		  phidotm=phidot;
	   }

	   de=(1.5*edot-.5*edotm)*dt_ph;
	   dnu=(1.5*nudot-.5*nudotm)*dt_ph;
	   dphi=(1.5*phidot-.5*phidotm)*dt_ph;
	   dgam=(1.5*gamdot-.5*gamdotm)*dt_ph;
	   dalp=(1.5*alpdot-.5*alpdotm)*dt_ph;

	   phiev[i_ini+i]=phi;
	   alpev[i_ini+i]=alp;
	   gamev[i_ini+i]=gam;
	   eccev[i_ini+i]=e;
	   nuev[i_ini+i]=nu;
	   gamdotev[i_ini+i]=gamdot;
	   alpdotev[i_ini+i]=alpdot;
	   Ampev[i_ini+i]=ampfct*pow(LISAWP_TWOPI*S.Mt*nu, 2./3.);
	   phiddotev[i_ini+i]=phiddot;
	   gamddotev[i_ini+i]=gamddot;
	   alpddotev[i_ini+i]=alpddot;
	   tm_ph[i_ini+i] = T;
	   e+=de;
	   phi+=dphi;
	   gam+=dgam;
	   alp+=dalp;
	   e2=e*e;
	   nu+=dnu;
	   T+= dt_ph;

	   rhs = pow( (1.0-e2)/(6.0+2.0*e), 1.5 )/(LISAWP_TWOPI * M);
	   if(rhs - nu <= 0.0){
		  std::cout << "*** we reached plunge at t = " << T << std::endl;
		  // std::cout << " i = " << i << std::endl;
		  imax = i;
		  S.tPl = T;
		  S.e_pl = e;
		  S.nu_pl = nu;
		  S.alpha_pl = alp;
		  S.gamma_pl = gam;
		  S.Phi_pl = phi;
			 break;
	   }
	   // LISA's motion
	   //  EccentricLISAMotion(float kappa0, float lambda0, double t, double* R, \
						  double** q, double** n)
	   EccentricLISAMotion(0.0, 0.0, T, R, q, n);
	   for(int j =0; j<3; j++){
		   kn[j] = 0.0;
			kq[j] = 0.0;
			nU = 0.0;
			nV = 0.0;
			for(int ii=0; ii<3; ii++){
			  kn[j] += k[ii]*n[j][ii];
				  kq[j] += k[ii]*q[j][ii];
				  nU += uhat[ii]*n[j][ii];
				  nV += vhat[ii]*n[j][ii];
			}
			u[j] = 0.5*(nU*nU - nV*nV);
			v[j] = nU*nV;
			// u[j] = 0.5*(nU*nU - nV*nV);
			// v[j] = -nU*nV;
	  }

	   ind = 0;
	   for (int ii=0; ii<5; ii++){ // nn harmonic
		   nn = (double)ii+1.;
		   for (int jj=0; jj<5; jj++){
			  mm = (double)jj-2.;
			  om = nn*LISAWP_TWOPI*nu + 2.*gamdot + mm*alpdot;
			  x = om*arm;
			  x2 = 0.5*x;
			  chi1 = -x*sin(x)*( Sinc(x2*(1.-kn[1]))*exp(-img*x) \
								 + Sinc(x2*(1.+kn[1])) )*exp(-img*x2*(3.0 + kq[0] + kq[2]));
			  chi2 = x*sin(x)*( Sinc(x2*(1.-kn[2])) + exp(-img*x)*\
										 Sinc(x2*(1.+kn[2])) )*exp(-img*x2*(3.0 + kq[1] + kq[0]));

			  Xp[ind][i_ini+i] = (u[1]*c2psi - v[1]*s2psi)*chi1 + (u[2]*c2psi - v[2]*s2psi)*chi2;
			  Xc[ind][i_ini+i] = (v[1]*c2psi + u[1]*s2psi)*chi1 + (v[2]*c2psi + u[2]*s2psi)*chi2;


			  chi2 = -x*sin(x)*( Sinc(x2*(1.-kn[2]))*exp(-img*x) \
								   + Sinc(x2*(1.+kn[2])) )*exp(-img*x2*(3.0 + kq[1] + kq[0]));
			  chi0 = x*sin(x)*( Sinc(x2*(1.-kn[0])) + exp(-img*x)*\
										 Sinc(x2*(1.+kn[0])) )*exp(-img*x2*(3.0 + kq[2] + kq[1]));

			  Yp[ind][i_ini+i] = (u[2]*c2psi - v[2]*s2psi)*chi2 + (u[0]*c2psi - v[0]*s2psi)*chi0;
			  Yc[ind][i_ini+i] = (v[2]*c2psi + u[2]*s2psi)*chi2 + (v[0]*c2psi + u[0]*s2psi)*chi0;


			  chi0 = -x*sin(x)*( Sinc(x2*(1.-kn[0]))*exp(-img*x) \
									 + Sinc(x2*(1.+kn[0])) )*exp(-img*x2*(3.0 + kq[2] + kq[1]));
			  chi1 = x*sin(x)*( Sinc(x2*(1.-kn[1])) + exp(-img*x)*\
											 Sinc(x2*(1.+kn[1])) )*exp(-img*x2*(3.0 + kq[0] + kq[2]));

			  Zp[ind][i_ini+i] = (u[0]*c2psi - v[0]*s2psi)*chi0 + (u[1]*c2psi - v[1]*s2psi)*chi1;
			  Zc[ind][i_ini+i] = (v[0]*c2psi + u[0]*s2psi)*chi0 + (v[1]*c2psi + u[1]*s2psi)*chi1;


	//          fout03  << T << spr << om << spr  << tp.real() << spr << tp.imag() << std::endl; //<< "   " << Xp[ind][i].imag() << "   " << \
						  Xc[ind][i].real() << "   " << Xc[ind][i].imag() << "   ";
			  ind++;

		   }
	   } // end of double loop over harmonics
	   i += 1;
  }
	// restore the initial conditions
  i = 0;
  e = S.e0;
  e2=e*e;
  Sp = S.a;
  nu = S.nu0;
  phi = S.Phi0;
  gam = S.gamma0;
  alp = S.alpha0;
  T = t_i_ini;

//    std::cout << "starting the loop backwards \n";
   //for (int i=0; i<Nph; i++){
   while ((T>=tStart) && (i_ini>=i)){

          e2=e*e;
          Y=1./(1.-e2);
          Z=pow(LISAWP_TWOPI*M*nu,1./3.);
          if (i != 0){
             edotm=edot;
             nudotm=nudot;
             phidotm=phidot;
             alpdotm=alpdot;
             gamdotm=gamdot;
          }
          edot = e*mu/M/M*(-1./15.*pow(Y,3.5)*pow(Z,8.)*((304.+121.*e2)/Y+Z*Z*(70648.-231960.*e2-56101.*e2*e2)/56.)+\
                   Sp*clam*pow(Z,11.)*pow(Y,4.)*(8184.+10064.*e2+789.*e2*e2)/30.);
          nudot = 96./(10.*M_PI)*mu/pow(M,3)*(pow(Z,11.)*pow(Y,4.5)*((96.+292.*e2+37.*e2*e2)/Y/96.\
       			       +Z*Z*(20368.-61464.*e2-163170.*e2*e2-13147.*e2*e2*e2)/5376.) - \
       			      pow(Z,14.)*pow(Y,5.)*Sp*clam*(1168.+9688.*e2+6286.*e2*e2 +195.*e2*e2*e2)/192.);
          phidot = 2.*LISAWP_PI*nu;
          alpdot = 8.*LISAWP_PI*LISAWP_PI*nu*nu*Sp*M*pow(Y,1.5);
          gamdot = 6.*LISAWP_PI*nu*Z*Z*Y*(1.+.25*Z*Z*Y*(26.-15.*e2))-3.*clam*alpdot;

          dalpdnu=16.*Sp*LISAWP_PI*nu*M*pow(Y,1.5);
          dalpde=12.*Sp*LISAWP_PI*nu*nu*M*sqrt(Y)*2.*e*Y*Y;
          dgamdnu=6.*Z*Z*Y*(1.+.25*Z*Z*Y*(26.-15.*e2))-3.*clam*dalpdnu+12.*Z*Y*(1.+.5*Z*Z*Y*(25.-15.*e2))*Z/(3.);
          dgamde=(6.*nu*Z*Z*(1.+.5*Z*Z*Y*(26.-15.*e2)))*2.*e*Y*Y-45.*nu*Z*Z*Z*Z*Y*Y*e-3.*clam*dalpde;

          alpddot=LISAWP_PI*(dalpdnu*nudot+dalpde*edot);
          gamddot=LISAWP_PI*(dgamdnu*nudot+dgamde*edot);
          phiddot=LISAWP_TWOPI*nudot;

          if (i == 0) {
             edotm=edot;
             nudotm=nudot;
             alpdotm=alpdot;
             gamdotm=gamdot;
             phidotm=phidot;
          }

          de=(1.5*edot-.5*edotm)*dt_ph;
          dnu=(1.5*nudot-.5*nudotm)*dt_ph;
          dphi=(1.5*phidot-.5*phidotm)*dt_ph;
          dgam=(1.5*gamdot-.5*gamdotm)*dt_ph;
          dalp=(1.5*alpdot-.5*alpdotm)*dt_ph;

          phiev[i_ini-i]=phi;
          alpev[i_ini-i]=alp;
          gamev[i_ini-i]=gam;
          eccev[i_ini-i]=e;
          nuev[i_ini-i]=nu;
          gamdotev[i_ini-i]=gamdot;
          alpdotev[i_ini-i]=alpdot;
          Ampev[i_ini-i]=ampfct*pow(LISAWP_TWOPI*S.Mt*nu, 2./3.);
          phiddotev[i_ini-i]=phiddot;
          gamddotev[i_ini-i]=gamddot;
          alpddotev[i_ini-i]=alpddot;
          tm_ph[i_ini-i] = T;
          e-=de;
          phi-=dphi;
          gam-=dgam;
          alp-=dalp;
          e2=e*e;
          nu-=dnu;
          T-= dt_ph;

		  // use additional stopping conditions
		  if (e > 0.9999)
		  		break;
		  if (isnan(e) || isnan(nu))
		  		break;


          // LISA's motion
          //  EccentricLISAMotion(float kappa0, float lambda0, double t, double* R, \
	                         double** q, double** n)
          EccentricLISAMotion(0.0, 0.0, T, R, q, n);
          for(int j =0; j<3; j++){
       	      kn[j] = 0.0;
               kq[j] = 0.0;
       		   nU = 0.0;
       		   nV = 0.0;
       		   for(int ii=0; ii<3; ii++){
       		   	 kn[j] += k[ii]*n[j][ii];
	       		 	 kq[j] += k[ii]*q[j][ii];
	       			 nU += uhat[ii]*n[j][ii];
	       			 nV += vhat[ii]*n[j][ii];
       		   }
       		   u[j] = 0.5*(nU*nU - nV*nV);
       		   v[j] = nU*nV;
       		   // u[j] = 0.5*(nU*nU - nV*nV);
       		   // v[j] = -nU*nV;
       	 }

          ind = 0;
          for (int ii=0; ii<5; ii++){ // nn harmonic
              nn = (double)ii+1.;
              for (int jj=0; jj<5; jj++){
                 mm = (double)jj-2.;
                 om = nn*LISAWP_TWOPI*nu + 2.*gamdot + mm*alpdot;
                 x = om*arm;
                 x2 = 0.5*x;
                 chi1 = -x*sin(x)*( Sinc(x2*(1.-kn[1]))*exp(-img*x) \
                                    + Sinc(x2*(1.+kn[1])) )*exp(-img*x2*(3.0 + kq[0] + kq[2]));
                 chi2 = x*sin(x)*( Sinc(x2*(1.-kn[2])) + exp(-img*x)*\
                 							Sinc(x2*(1.+kn[2])) )*exp(-img*x2*(3.0 + kq[1] + kq[0]));

                 Xp[ind][i_ini-i] = (u[1]*c2psi - v[1]*s2psi)*chi1 + (u[2]*c2psi - v[2]*s2psi)*chi2;
                 Xc[ind][i_ini-i] = (v[1]*c2psi + u[1]*s2psi)*chi1 + (v[2]*c2psi + u[2]*s2psi)*chi2;


                 chi2 = -x*sin(x)*( Sinc(x2*(1.-kn[2]))*exp(-img*x) \
                                      + Sinc(x2*(1.+kn[2])) )*exp(-img*x2*(3.0 + kq[1] + kq[0]));
                 chi0 = x*sin(x)*( Sinc(x2*(1.-kn[0])) + exp(-img*x)*\
                   							Sinc(x2*(1.+kn[0])) )*exp(-img*x2*(3.0 + kq[2] + kq[1]));

                 Yp[ind][i_ini-i] = (u[2]*c2psi - v[2]*s2psi)*chi2 + (u[0]*c2psi - v[0]*s2psi)*chi0;
                 Yc[ind][i_ini-i] = (v[2]*c2psi + u[2]*s2psi)*chi2 + (v[0]*c2psi + u[0]*s2psi)*chi0;


                 chi0 = -x*sin(x)*( Sinc(x2*(1.-kn[0]))*exp(-img*x) \
                                        + Sinc(x2*(1.+kn[0])) )*exp(-img*x2*(3.0 + kq[2] + kq[1]));
                 chi1 = x*sin(x)*( Sinc(x2*(1.-kn[1])) + exp(-img*x)*\
                     							Sinc(x2*(1.+kn[1])) )*exp(-img*x2*(3.0 + kq[0] + kq[2]));

                 Zp[ind][i_ini-i] = (u[0]*c2psi - v[0]*s2psi)*chi0 + (u[1]*c2psi - v[1]*s2psi)*chi1;
                 Zc[ind][i_ini-i] = (v[0]*c2psi + u[0]*s2psi)*chi0 + (v[1]*c2psi + u[1]*s2psi)*chi1;


       //          fout03  << T << spr << om << spr  << tp.real() << spr << tp.imag() << std::endl; //<< "   " << Xp[ind][i].imag() << "   " << \
                             Xc[ind][i].real() << "   " << Xc[ind][i].imag() << "   ";
                 ind++;

              }
          }
	 i += 1;
     }
	 // std::cout << "i_ini = " << i_ini << " check " << i_ini*dt_ph << " i = " << i << "  T = " << T << std::endl;
   delete [] R;
   for (int i=0; i<3; i++){
       delete [] q[i];
       delete [] n[i];
   }
   delete [] q;
   delete [] n;

}

void FreqAK_RA::NewPhaseEv_RA(EMRItemplate& S, double t_ini, double dur, double tStart){
	/**
	@param S EMRItemplate object storing the parameters
	@param t0 moment of time for which intitial conditions are given
	@param maximum dur durarion of the signal
	@param tStart - how much back in time we should go.
	**/

	double k[3];
	double uhat[3];
	double vhat[3];
	double u[3];
	double v[3];
	double kn[3];
	double kq[3];

   double up = (S.ctS*S.stK*cos(S.phS - S.phK) - S.ctK*S.stS);
   double dw = (S.stK*sin(S.phS-S.phK));
	//  std::cout << up << "  " << dw << "  " << S.ctS << "  " << S.stS << "   " <<  S.phS << std::endl;
   double psi;
   if (dw != 0.0) {
      // psi = atan2(up, dw);
      // psi = 0.5*LISAWP_PI -atan2(up, dw);
      psi = -atan2(up, dw);
      // psi = atan2(up, dw);
			// std::cout << "psi = " << psi << std::endl;
   }else {
      psi = 0.5*LISAWP_PI;
   }

   double c2psi=cos(2.*psi);
   double s2psi=sin(2.*psi);


   // note that k = -n, where k is propagation vector and n is sky location
   k[0] =  -S.stS*S.cpS;
   k[1] = -S.stS*S.spS;
   k[2] = -S.ctS;

   uhat[0] = S.ctS*S.cpS;
   uhat[1] = S.ctS*S.spS;
   uhat[2] = -S.stS;

   vhat[0] = S.spS;
   vhat[1] = -S.cpS;
   vhat[2] = 0.0;
   double nU, nV;

   S.t0 = 0.0;
   // double e0 = S.e0;
   // double nu0 = S.nu0;
   // double ph0 = S.Phi0;
   // double al0 =  S.alpha0;
   // double gam0 = S.gamma0;

   double* R;
   double** q;
   double** n;
   R = new double[3];
   q = new double*[3];
   n = new double*[3];
   for (int i=0; i<3; i++){
       q[i] = new double[3];
       n[i] = new double[3];
   }


   double clam=cos(S.lam);
   double slam=sin(S.lam);



   double M = S.Mt;
   double mu = S.mt;
   double e = S.e0;
   double e2=e*e;
   double Sp = S.a;
   double nu = S.nu0;
   double phi = S.Phi0;
   double gam = S.gamma0;
   double alp = S.alpha0;
   double quadmom = S.qm;
   double Y, Z;
   double edotm, nudotm, phidotm, alpdotm, gamdotm;
   double edot, nudot, phidot, alpdot, gamdot;
   double dalpdnu, dalpde, dgamdnu, dgamde, alpddot, gamddot, phiddot;
   double de, dnu, dphi, dgam, dalp, rhs;
   double T=0.0;
   int imax, ind;
   double ampfct= S.Ampl;

   double nn, mm, om;
   std::complex<double> chi0;
   std::complex<double> chi1;
   std::complex<double> chi2;
   std::complex<double> img(0.0, 1.0);
   // double x, x2;

   // given the courseneess of the steps we need first to interpolate to the nearest dt*i
   int i_ini = floor(t_ini/dt_ph);
   double t_i_ini = i_ini*dt_ph;
   double dt_ini = t_i_ini - t_ini;

   Y=1./(1.-e2);
   Z=pow(LISAWP_TWOPI*M*nu,1./3.);
   edot = e*mu/M/M*(-1./15.*pow(Y,3.5)*pow(Z,8.)*((304.+121.*e2)/Y+Z*Z*(70648.-231960.*e2-56101.*e2*e2)/56.)+\
			Sp*clam*pow(Z,11.)*pow(Y,4.)*(8184.+10064.*e2+789.*e2*e2)/30.);
   nudot = 96./(10.*M_PI)*mu/pow(M,3)*(pow(Z,11.)*pow(Y,4.5)*((96.+292.*e2+37.*e2*e2)/Y/96.\
							   +Z*Z*(20368.-61464.*e2-163170.*e2*e2-13147.*e2*e2*e2)/5376.)
				    -pow(Z,14.)*pow(Y,5.)*Sp*clam*(1168.+9688.*e2+6286.*e2*e2 +195.*e2*e2*e2)/192.\
				    -pow(Z, 15.)*pow(Y,5.5)*quadmom*(33./16.+359.*e2/32.-527.*(1.-clam*clam)/96.));
//TODO Done !

   phidot = 2.*LISAWP_PI*nu;
   alpdot = 8.*LISAWP_PI*LISAWP_PI*nu*nu*Sp*M*pow(Y,1.5);
   gamdot = 6.*LISAWP_PI*nu*Z*Z*Y*(1.+.25*Z*Z*Y*(26.-15.*e2))-3.*clam*alpdot;
   alpdot = alpdot + 3.*M_PI*nu*quadmom*Z*Z*Z*Z*Y*Y*clam;
   gamdot = gamdot - 1.5*M_PI*nu*quadmom*Z*Z*Z*Z*Y*Y*(5.*clam-1.);
//TODO Done !

	// change the initial conditions to t_i_ini.
   S.e0 =  e + edot*dt_ini;
   S.Phi0 = phi + phidot*dt_ini;
   S.gamma0 =  gam + gamdot*dt_ini;
   S.alpha0 = alp + alpdot*dt_ini;
   S.nu0 = nu + nudot*dt_ini;
   S.t0 = t_i_ini;
   T = t_ini + dt_ini;

   // calculate how much we need forward and how much back

   double t_back = t_i_ini - tStart;
   if (t_back < 0.0){
	   std::cout << "tstart before t_ini  " << tStart << "  " << t_ini << "  " << t_i_ini << std::endl;
	   exit(1);
   }

   double t_forward = dur + tStart - t_ini;
   if (t_forward < 0.0){
	   std::cout << "duration is not compatible with t_ini and tStart  " << dur << "   " << tStart << "   " << t_ini << std::endl;
   }
   double t_end = dur + tStart; // end time

   if (dur > Tobs){
	   std::cout << "new duration is larger than observation time " << dur << "   " << Tobs << std::endl;
   }

   // start forward loop

   // std::cout << "New: t_ini = " << t_ini << "  t_i_ini = " << t_i_ini << "  T = " << T << " tStart = " << tStart << \
   // 			"  t_end = " << t_end << "  dur = " << dur << "  " << t_i_ini - dur << std::endl;

 	e = S.e0;
	e2=e*e;
	Sp = S.a;
	nu = S.nu0;
	phi = S.Phi0;
	gam = S.gamma0;
	alp = S.alpha0;

	// std::cout << "IC: e0 = " << e << " nu= " << nu << " spin = " << Sp << std::endl;

   int i = 0;
   while (T <= t_end){
	   e2=e*e;
	   Y=1./(1.-e2);
	   Z=pow(LISAWP_TWOPI*M*nu,1./3.);
	   if (i != 0){
		  edotm=edot;
		  nudotm=nudot;
		  phidotm=phidot;
		  alpdotm=alpdot;
		  gamdotm=gamdot;
	   }
	   edot = e*mu/M/M*(-1./15.*pow(Y,3.5)*pow(Z,8.)*((304.+121.*e2)/Y+Z*Z*(70648.-231960.*e2-56101.*e2*e2)/56.)+\
				Sp*clam*pow(Z,11.)*pow(Y,4.)*(8184.+10064.*e2+789.*e2*e2)/30.);
	   nudot = 96./(10.*M_PI)*mu/pow(M,3)*(pow(Z,11.)*pow(Y,4.5)*((96.+292.*e2+37.*e2*e2)/Y/96.\
							   +Z*Z*(20368.-61464.*e2-163170.*e2*e2-13147.*e2*e2*e2)/5376.)\
				    -pow(Z,14.)*pow(Y,5.)*Sp*clam*(1168.+9688.*e2+6286.*e2*e2 +195.*e2*e2*e2)/192.\
				     -pow(Z, 15.)*pow(Y,5.5)*quadmom*(33./16.+359.*e2/32.-527.*(1.-clam*clam)/96.));
//TODO Done !
	   phidot = 2.*LISAWP_PI*nu;
	   alpdot = 8.*LISAWP_PI*LISAWP_PI*nu*nu*Sp*M*pow(Y,1.5);
	   gamdot = 6.*LISAWP_PI*nu*Z*Z*Y*(1.+.25*Z*Z*Y*(26.-15.*e2))-3.*clam*alpdot;
	   alpdot = alpdot + 3.*M_PI*nu*quadmom*Z*Z*Z*Z*Y*Y*clam;
	   gamdot = gamdot - 1.5*M_PI*nu*quadmom*Z*Z*Z*Z*Y*Y*(5.*clam-1.);
//TODO Done !

	   dalpdnu=16.*Sp*LISAWP_PI*nu*M*pow(Y,1.5);
	   dalpde=12.*Sp*LISAWP_PI*nu*nu*M*sqrt(Y)*2.*e*Y*Y;
	   dgamdnu=6.*Z*Z*Y*(1.+.25*Z*Z*Y*(26.-15.*e2))-3.*clam*dalpdnu+12.*Z*Y*(1.+.5*Z*Z*Y*(25.-15.*e2))*Z/(3.);
	   dgamde=(6.*nu*Z*Z*(1.+.5*Z*Z*Y*(26.-15.*e2)))*2.*e*Y*Y-45.*nu*Z*Z*Z*Z*Y*Y*e-3.*clam*dalpde;

	   alpddot=LISAWP_PI*(dalpdnu*nudot+dalpde*edot);
	   gamddot=LISAWP_PI*(dgamdnu*nudot+dgamde*edot);
	   phiddot=LISAWP_TWOPI*nudot;

	   if (i == 0) {
		  edotm=edot;
		  nudotm=nudot;
		  alpdotm=alpdot;
		  gamdotm=gamdot;
		  phidotm=phidot;
	   }

	   de=(1.5*edot-.5*edotm)*dt_ph;
	   dnu=(1.5*nudot-.5*nudotm)*dt_ph;
	   dphi=(1.5*phidot-.5*phidotm)*dt_ph;
	   dgam=(1.5*gamdot-.5*gamdotm)*dt_ph;
	   dalp=(1.5*alpdot-.5*alpdotm)*dt_ph;

	   phiev[i_ini+i]=phi;
	   alpev[i_ini+i]=alp;
	   gamev[i_ini+i]=gam;
	   eccev[i_ini+i]=e;
	   nuev[i_ini+i]=nu;
	   gamdotev[i_ini+i]=gamdot;
	   alpdotev[i_ini+i]=alpdot;
	   Ampev[i_ini+i]=ampfct*pow(LISAWP_TWOPI*S.Mt*nu, 2./3.);
	   phiddotev[i_ini+i]=phiddot;
	   gamddotev[i_ini+i]=gamddot;
	   alpddotev[i_ini+i]=alpddot;
	   tm_ph[i_ini+i] = T;
	   e+=de;
	   phi+=dphi;
	   gam+=dgam;
	   alp+=dalp;
	   e2=e*e;
	   nu+=dnu;
	   T+= dt_ph;

	   rhs = pow( (1.0-e2)/(6.0+2.0*e), 1.5 )/(LISAWP_TWOPI * M);
	   if(rhs - nu <= 0.0){
		 //  std::cout << "*** we reached plunge at t = " << T << std::endl;
		  // std::cout << " i = " << i << std::endl;
		  imax = i;
		  S.tPl = T;
		  S.e_pl = e;
		  S.nu_pl = nu;
		  S.alpha_pl = alp;
		  S.gamma_pl = gam;
		  S.Phi_pl = phi;
			 break;
	   }
	   // LISA's motion
	   //  EccentricLISAMotion(float kappa0, float lambda0, double t, double* R, \
	   //					  double** q, double** n)
	   EccentricLISAMotion2(lisa_orbits, T, q, n);
	   for(int j =0; j<3; j++){
		   kn[j] = 0.0;
			kq[j] = 0.0;
			nU = 0.0;
			nV = 0.0;
			for(int ii=0; ii<3; ii++){
			  kn[j] += k[ii]*n[j][ii];
				  kq[j] += k[ii]*q[j][ii];
				  nU += uhat[ii]*n[j][ii];
				  nV += vhat[ii]*n[j][ii];
			}
			u[j] = 0.5*(nU*nU - nV*nV);
			v[j] = nU*nV;
			// u[j] = 0.5*(nU*nU - nV*nV);
			// v[j] = -nU*nV;
	  }

	   ind = 0;
	   for (int ii=0; ii<5; ii++){ // nn harmonic
		   nn = (double)ii+1.;
		   for (int jj=0; jj<5; jj++){
			  mm = (double)jj-2.;
			  om = nn*LISAWP_TWOPI*nu + 2.*gamdot + mm*alpdot;
			  ComputeResponse(om, kn, kq, u, v, c2psi, s2psi, &Xp[ind][i_ini+i], &Xc[ind][i_ini+i], &Yp[ind][i_ini+i],\
				  &Yc[ind][i_ini+i], &Zp[ind][i_ini+i], &Zc[ind][i_ini+i]);
			  if (Ext >0){
			  	// next subdominant harmonic
			  	om = nn*LISAWP_TWOPI*nu  + mm*alpdot;
			  	ComputeResponse(om, kn, kq, u, v, c2psi, s2psi, &XXp[ind][i_ini+i], &XXc[ind][i_ini+i], &YYp[ind][i_ini+i],\
					  	&YYc[ind][i_ini+i], &ZZp[ind][i_ini+i], &ZZc[ind][i_ini+i]);
			  }

			  ind++;

		   }
	   } // end of double loop over harmonics
	   i += 1;
  }
	// restore the initial conditions
  i = 0;
  e = S.e0;
  e2=e*e;
  Sp = S.a;
  nu = S.nu0;
  phi = S.Phi0;
  gam = S.gamma0;
  alp = S.alpha0;
  T = t_i_ini;

//    std::cout << "starting the loop backwards \n";
   //for (int i=0; i<Nph; i++){
   while ((T>=tStart) && (i_ini>=i)){

          e2=e*e;
          Y=1./(1.-e2);
          Z=pow(LISAWP_TWOPI*M*nu,1./3.);
          if (i != 0){
             edotm=edot;
             nudotm=nudot;
             phidotm=phidot;
             alpdotm=alpdot;
             gamdotm=gamdot;
          }
          edot = e*mu/M/M*(-1./15.*pow(Y,3.5)*pow(Z,8.)*((304.+121.*e2)/Y+Z*Z*(70648.-231960.*e2-56101.*e2*e2)/56.)+\
                   Sp*clam*pow(Z,11.)*pow(Y,4.)*(8184.+10064.*e2+789.*e2*e2)/30.);
          nudot = 96./(10.*M_PI)*mu/pow(M,3)*(pow(Z,11.)*pow(Y,4.5)*((96.+292.*e2+37.*e2*e2)/Y/96.\
							   +Z*Z*(20368.-61464.*e2-163170.*e2*e2-13147.*e2*e2*e2)/5376.)\
				    -pow(Z,14.)*pow(Y,5.)*Sp*clam*(1168.+9688.*e2+6286.*e2*e2 +195.*e2*e2*e2)/192.\
				     -pow(Z, 15.)*pow(Y,5.5)*quadmom*(33./16.+359.*e2/32.-527.*(1.-clam*clam)/96.));
//TODO Done !
          phidot = 2.*LISAWP_PI*nu;
          alpdot = 8.*LISAWP_PI*LISAWP_PI*nu*nu*Sp*M*pow(Y,1.5);
          gamdot = 6.*LISAWP_PI*nu*Z*Z*Y*(1.+.25*Z*Z*Y*(26.-15.*e2))-3.*clam*alpdot;
          alpdot = alpdot + 3.*M_PI*nu*quadmom*Z*Z*Z*Z*Y*Y*clam;
          gamdot = gamdot - 1.5*M_PI*nu*quadmom*Z*Z*Z*Z*Y*Y*(5.*clam-1.);
//TODO Done !

          dalpdnu=16.*Sp*LISAWP_PI*nu*M*pow(Y,1.5);
          dalpde=12.*Sp*LISAWP_PI*nu*nu*M*sqrt(Y)*2.*e*Y*Y;
          dgamdnu=6.*Z*Z*Y*(1.+.25*Z*Z*Y*(26.-15.*e2))-3.*clam*dalpdnu+12.*Z*Y*(1.+.5*Z*Z*Y*(25.-15.*e2))*Z/(3.);
          dgamde=(6.*nu*Z*Z*(1.+.5*Z*Z*Y*(26.-15.*e2)))*2.*e*Y*Y-45.*nu*Z*Z*Z*Z*Y*Y*e-3.*clam*dalpde;

          alpddot=LISAWP_PI*(dalpdnu*nudot+dalpde*edot);
          gamddot=LISAWP_PI*(dgamdnu*nudot+dgamde*edot);
          phiddot=LISAWP_TWOPI*nudot;

          if (i == 0) {
             edotm=edot;
             nudotm=nudot;
             alpdotm=alpdot;
             gamdotm=gamdot;
             phidotm=phidot;
          }

          de=(1.5*edot-.5*edotm)*dt_ph;
          dnu=(1.5*nudot-.5*nudotm)*dt_ph;
          dphi=(1.5*phidot-.5*phidotm)*dt_ph;
          dgam=(1.5*gamdot-.5*gamdotm)*dt_ph;
          dalp=(1.5*alpdot-.5*alpdotm)*dt_ph;

          phiev[i_ini-i]=phi;
          alpev[i_ini-i]=alp;
          gamev[i_ini-i]=gam;
          eccev[i_ini-i]=e;
          nuev[i_ini-i]=nu;
          gamdotev[i_ini-i]=gamdot;
          alpdotev[i_ini-i]=alpdot;
          Ampev[i_ini-i]=ampfct*pow(LISAWP_TWOPI*S.Mt*nu, 2./3.);
          phiddotev[i_ini-i]=phiddot;
          gamddotev[i_ini-i]=gamddot;
          alpddotev[i_ini-i]=alpddot;
          tm_ph[i_ini-i] = T;
          e-=de;
          phi-=dphi;
          gam-=dgam;
          alp-=dalp;
          e2=e*e;
          nu-=dnu;
          T-= dt_ph;

		  // use additional stopping conditions
		  if (e > 0.9999)
		  		break;
		  if (::isnan(e) || ::isnan(nu))
		  		break;


          // LISA's motion
          //  EccentricLISAMotion(float kappa0, float lambda0, double t, double* R, \
		  //double** q, double** n)
	  EccentricLISAMotion2(lisa_orbits, T, q, n);
          for(int j =0; j<3; j++){
       	      kn[j] = 0.0;
               kq[j] = 0.0;
       		   nU = 0.0;
       		   nV = 0.0;
       		   for(int ii=0; ii<3; ii++){
       		   	 kn[j] += k[ii]*n[j][ii];
	       		 	 kq[j] += k[ii]*q[j][ii];
	       			 nU += uhat[ii]*n[j][ii];
	       			 nV += vhat[ii]*n[j][ii];
       		   }
       		   u[j] = 0.5*(nU*nU - nV*nV);
       		   v[j] = nU*nV;
       		   // u[j] = 0.5*(nU*nU - nV*nV);
       		   // v[j] = -nU*nV;
       	 }

          ind = 0;
          for (int ii=0; ii<5; ii++){ // nn harmonic
              nn = (double)ii+1.;
              for (int jj=0; jj<5; jj++){
                 mm = (double)jj-2.;
                 om = nn*LISAWP_TWOPI*nu + 2.*gamdot + mm*alpdot;
				 ComputeResponse(om, kn, kq, u, v, c2psi, s2psi, &Xp[ind][i_ini-i], &Xc[ind][i_ini-i], &Yp[ind][i_ini-i],
					 &Yc[ind][i_ini-i], &Zp[ind][i_ini-i], &Zc[ind][i_ini-i]);

				 if (Ext >0){
				 	// next subdominant harmonic
				 	om = nn*LISAWP_TWOPI*nu  + mm*alpdot;
				 	ComputeResponse(om, kn, kq, u, v, c2psi, s2psi, &XXp[ind][i_ini-i], &XXc[ind][i_ini-i],
						&YYp[ind][i_ini-i], &YYc[ind][i_ini-i], &ZZp[ind][i_ini-i], &ZZc[ind][i_ini-i]);
				 }
                 ind++;

              }
          }
	 i += 1;
     }
	 // std::cout << "i_ini = " << i_ini << " check " << i_ini*dt_ph << " i = " << i << "  T = " << T << std::endl;
   delete [] R;
   for (int i=0; i<3; i++){
       delete [] q[i];
       delete [] n[i];
   }
   delete [] q;
   delete [] n;

}

void FreqAK_RA::ComputeResponse(double om, double* kn, double* kq, double* u, double *v, double c2psi, double s2psi,
	std::complex<double>* Xplus, std::complex<double>* Xcross, std::complex<double>* Yplus,
	std::complex<double>* Ycross, std::complex<double>* Zplus, std::complex<double>* Zcross){

	double x = om*arm;
	double x2 = 0.5*x;
	std::complex<double> chi0;
    std::complex<double> chi1;
    std::complex<double> chi2;
    std::complex<double> img(0.0, 1.0);

	chi1 = -x*sin(x)*( Sinc(x2*(1.-kn[1]))*exp(-img*x) \
					   + Sinc(x2*(1.+kn[1])) )*exp(-img*x2*(3.0 + kq[0] + kq[2]));
	chi2 = x*sin(x)*( Sinc(x2*(1.-kn[2])) + exp(-img*x)*\
							   Sinc(x2*(1.+kn[2])) )*exp(-img*x2*(3.0 + kq[1] + kq[0]));

	*Xplus = (u[1]*c2psi - v[1]*s2psi)*chi1 + (u[2]*c2psi - v[2]*s2psi)*chi2;
	*Xcross = (v[1]*c2psi + u[1]*s2psi)*chi1 + (v[2]*c2psi + u[2]*s2psi)*chi2;


	chi2 = -x*sin(x)*( Sinc(x2*(1.-kn[2]))*exp(-img*x) \
						 + Sinc(x2*(1.+kn[2])) )*exp(-img*x2*(3.0 + kq[1] + kq[0]));
	chi0 = x*sin(x)*( Sinc(x2*(1.-kn[0])) + exp(-img*x)*\
							   Sinc(x2*(1.+kn[0])) )*exp(-img*x2*(3.0 + kq[2] + kq[1]));

	*Yplus = (u[2]*c2psi - v[2]*s2psi)*chi2 + (u[0]*c2psi - v[0]*s2psi)*chi0;
	*Ycross = (v[2]*c2psi + u[2]*s2psi)*chi2 + (v[0]*c2psi + u[0]*s2psi)*chi0;


	chi0 = -x*sin(x)*( Sinc(x2*(1.-kn[0]))*exp(-img*x) \
						   + Sinc(x2*(1.+kn[0])) )*exp(-img*x2*(3.0 + kq[2] + kq[1]));
	chi1 = x*sin(x)*( Sinc(x2*(1.-kn[1])) + exp(-img*x)*\
								   Sinc(x2*(1.+kn[1])) )*exp(-img*x2*(3.0 + kq[0] + kq[2]));

	*Zplus = (u[0]*c2psi - v[0]*s2psi)*chi0 + (u[1]*c2psi - v[1]*s2psi)*chi1;
	*Zcross = (v[0]*c2psi + u[0]*s2psi)*chi0 + (v[1]*c2psi + u[1]*s2psi)*chi1;

}

void FreqAK_RA::GetEvolution(double* t_phase, double* nus, double* eccs, double* phis, double* alps, double* gams, double* fhigh){

	// double StasBessel[8];
	// int state;
	for (int i=0; i<Nph; i++){
			t_phase[i] = tm_ph[i];
			nus[i] = nuev[i];
			eccs[i] = eccev[i];
			// StasBessel = gsl_sf_bessel_J0(eccs[i]);
			// state = gsl_sf_bessel_Jn_array(0, 7, eccs[i], StasBessel);
			// std::cout << "stas bessel test " << StasBessel[0] << "   " << StasBessel[7] << "  "<< eccs[i] << std::endl;
			phis[i] = phiev[i];
			alps[i] = alpev[i];
			gams[i] = gamev[i];
			fhigh[i]=  5.0*LISAWP_TWOPI*nuev[i] + 2.*gamdotev[i]  + 2.0*alpdotev[i];
	}

}

void FreqAK_RA::ComputeHAmpl(double ec, int n, double* XaXb, double* Xa_Xb, double* Xc){

	double output;

	// Compute Bessel functions
	double Bss[5];
	double Xa, Xb;
	// for (int i=1; i<=nH; i++){
	if (ec <= 0.45){ // Use taylor expansion
		double ec2=ec*ec;
		double ec3=ec2*ec;
		double ec4=ec3*ec;
		double ec5=ec4*ec;
		double ec6=ec5*ec;
		double ec7=ec6*ec;
		*Xa_Xb = 0.0;
	   switch (n) {
		  case 1:
			  *XaXb =3.*ec-1.625*ec3;
			  *Xc = ec - ec3/8. + ec5/192.;
			  break;
		  case 2:
			  *XaXb =-4.+10.*ec2-5.75*ec4;
			  *Xc = ec2 - ec4/3.+ ec6/24.;
			  break;
		  case 3:
			  *XaXb =-9.*ec+21.375*ec3-15.046875*ec5;
			  *Xc = 9.*ec3/8. - 81.*ec5/128.;
			  break;
		  case 4:
			  *XaXb =-16.*ec2+40.*ec4-101.*ec6/3.;
			  *Xc = 4.*ec4/3. - 16.*ec6/15.;
			  break;
		  case 5:
			  *XaXb =-625.*ec3/24.+26875.*ec5/384.-210625.*ec7/3072.;
			  *Xc = 625.*ec5/384. - 15625.*ec7/9216.;
			  break;
		}
	}else{
		if (n==1){
			Bss[0] = gsl_sf_bessel_Jn(n-2, n*ec);
			Bss[1] = gsl_sf_bessel_Jn(n-1, n*ec);
			Bss[2] = gsl_sf_bessel_Jn(n, n*ec);
			Bss[3] = gsl_sf_bessel_Jn(n+1, n*ec);
			Bss[4] = gsl_sf_bessel_Jn(n+2, n*ec);
		}
		else{
			int state = gsl_sf_bessel_Jn_array(n-2, n+2, n*ec, Bss);
		}
		Xa = Bss[0] - Bss[4] - 2.0*ec*(Bss[1] - Bss[3]) + 2.0/n*Bss[2];
		Xb = sqrt(1. - ec*ec)*( Bss[0] + Bss[4] - 2.0*Bss[2] );
		*Xc = 2.0*Bss[2];
		*XaXb = -n*(Xa+Xb); // dominant harmonic
		*Xa_Xb = -n*(Xa - Xb);
	}

}

// TODO interface is different now, we do not us MAtrix anymore!!!!
void FreqAK_RA::ConstructWaveFreq_RA(EMRItemplate& S, double* Xf_r, double* Xf_im,
				     double* Yf_r, double* Yf_im,
				     double* Zf_r, double* Zf_im)

{

     // NOTE!!!! Xf, Yf, Zf must have proper size and be zero

	 int debug =0;
     double T = 0.0;
     double ec, ec2, ec3, ec4, ec5, ec6, ec7, hamp;
     int ind_low, ind_up;
     double delta, eps, amp;
     double  xi, Apc, Aps, Acc, Acs;
     int harmms[]={-2,-1,0,1,2};
     double fact;
     std::complex<double> hplus, hcross;
     std::complex<double> cFp;
     std::complex<double> cFc;
     std::complex<double> x, xp, xc;
     std::complex<double> img(0.0, 1.0);
     double om_in;
     double om_fin;
     double dOm = df*LISAWP_TWOPI;
     double Om, dom, delom;
	 double Om_Nyq = LISAWP_TWOPI*0.5/dt_w; //we limit freq. content of the waveform to this freq.
     double xi_in, xi_fin, dxi;
     double orbOm = LISAWP_TWOPI/LISAWP_YRSID_SI;
     double AUsec =  LISAWP_AU_SI/LISAWP_C_SI;
     double DM = AUsec*S.stS;
     double faza, sinph, cosph;
     int ind;
     double nn,mm;
     std::string spr = "    ";
	 double XXc, Xb_a;

     double cS= S.ctS;
     double sS= S.stS;
     double cK= S.ctK;
     double sK= S.stK;
     double cSp= S.cpS;
     double sSp= S.spS;
     double cKp= S.cpK;
     double sKp= S.spK;
     double clam=cos(S.lam);
     double slam=sin(S.lam);
     double Sn = cS*cK + sS*sK*cos(S.phS  - S. phK);
     double cX = Sn;

     double sX = sqrt( sS*sS*cK*cK - 2.*sS*sSp*cK*cS*sK*sKp + \
                      cS*cS*sK*sK - 2.*cS*sK*cKp*sS*cSp*cK + \
                      sS*sS*cSp*cSp*sK*sK*sKp*sKp - 2.*sS*sS*cSp*sK*sK*sKp*sSp*cKp +\
                      sS*sS*sSp*sSp*sK*sK*cKp*cKp);

     double Apc1, Aps1, Apcn1, Apc2, Aps2, Apcn2, Aqc1, Aqs1, Aqcn1, Aqc2, Aqs2, Aqcn2;

       Apc1 = ( -cK*cKp*sS*cSp - cK*sKp*sS*sSp + sK*cS )/(sX);
       Aps1 = ( sKp*sS*cSp - cKp*sS*sSp )/(sX);
       Apcn1 = ( cK*cS + sK*cKp*sS*cSp + sK*sKp*sS*sSp - cX)*clam/(sX*slam);

       Apc2 = (sS*cSp*sKp - sS*sSp*cKp )*clam/(sX);
       Aps2 = ( cK*cKp*sS*cSp + cK*sKp*sS*sSp - cS*sK )*clam/(sX);
       Apcn2 = 0.0;


       Aqc1 = ( sS*cSp*sKp - sS*sSp*cKp  )*cX/(sX);
       Aqs1 = ( cK*cKp*sS*cSp + cK*sKp*sS*sSp - cS*sK )*cX/(sX);
       Aqcn1 = 0.0;


       Aqc2 = cX*clam*( cK*cKp*sS*cSp + cK*sKp*sS*sSp - sK*cS)/(sX);
       Aqs2 = -cX*clam*sS*( sKp*cSp - cKp*sSp )/(sX);
       Aqcn2 = -( cX*clam*clam*( cK*cS + sK*cKp*sS*cSp + sK*sKp*sS*sSp ) + \
                                1.- cX*cX - clam*clam )/(sX*slam);


     double Bp1c1 = 2.0*(Apc1*Apcn1 - Aqc1*Aqcn1 + Aqc2*Aqcn2 - Apc2*Apcn2);
     double Bp1c2 =  0.5*(Aps2*Aps2 - Aqc1*Aqc1  + Apc1*Apc1  - Aps1*Aps1 + \
                                Aqc2*Aqc2 + Aqs1*Aqs1 - Apc2*Apc2 - Aqs2*Aqs2);
     double Bp1s1 = 2.0*(Aqs2*Aqcn2 - Aps2*Apcn2 - Aqs1*Aqcn1 + Aps1*Apcn1);
     double Bp1s2 = (Apc1*Aps1 + Aqc2*Aqs2 - Apc2*Aps2 - Aqc1*Aqs1);
     double Bp1cn = 0.5*(Apc1*Apc1 + Aps1*Aps1 - Aqc1*Aqc1 - Aqs1*Aqs1 - Apc2*Apc2 \
                                + Aqc2*Aqc2 + Aqs2*Aqs2 - Aps2*Aps2) + Aqcn2*Aqcn2 - Aqcn1*Aqcn1 \
                                + Apcn1*Apcn1 - Apcn2*Apcn2;

     double Bp2c1 = (Apcn1*Apc2 + Apc1*Apcn2 - Aqcn1*Aqc2 - Aqc1*Aqcn2);
     double Bp2c2 = 0.5*(Aqs1*Aqs2 - Aps1*Aps2 + Apc1*Apc2 - Aqc1*Aqc2);
     double Bp2s1 = (Aps1*Apcn2 + Apcn1*Aps2 - Aqcn1*Aqs2 - Aqs1*Aqcn2);
     double Bp2s2 = 0.5*( Apc1*Aps2 - Aqc1*Aqs2 + Aps1*Apc2 - Aqs1*Aqc2);
     double Bp2cn = 0.5*(Aps1*Aps2 - Aqs1*Aqs2 - Aqc1*Aqc2 + Apc1*Apc2) -Aqcn1*Aqcn2 + Apcn1*Apcn2;

     double Bc1c1 = (-Apc2*Aqcn2 - Apcn2*Aqc2 + Apc1*Aqcn1 + Apcn1*Aqc1);
     double Bc1c2 = 0.5*( Apc1*Aqc1 - Aps1*Aqs1 - Apc2*Aqc2 + Aps2*Aqs2);
     double Bc1s1 = (Apcn1*Aqs1 - Aps2*Aqcn2 + Aps1*Aqcn1 - Apcn2*Aqs2);
     double Bc1s2 = 0.5*(-Apc2*Aqs2 + Apc1*Aqs1 - Aps2*Aqc2 + Aps1*Aqc1);
     double Bc1cn = -Apcn2*Aqcn2 + Apcn1*Aqcn1 + 0.5*(Apc1*Aqc1 - Aps2*Aqs2 + Aps1*Aqs1 - Apc2*Aqc2);

     double Bc2c1 = (Aqc1*Apcn2 + Aqcn1*Apc2 + Apc1*Aqcn2 + Apcn1*Aqc2);
     double Bc2c2 = 0.5*( Apc1*Aqc2 - Aps1*Aqs2 + Aqc1*Apc2 - Aqs1*Aps2);
     double Bc2s1 = (Apcn1*Aqs2 + Aqs1*Apcn2 + Aps1*Aqcn2 + Aqcn1*Aps2);
     double Bc2s2 = 0.5*(Aqc1*Aps2 + Apc1*Aqs2 + Aqs1*Apc2 + Aps1*Aqc2);
     double Bc2cn = Aqcn1*Apcn2 + Apcn1*Aqcn2 + 0.5*(Apc1*Aqc2 + Aqs1*Aps2 +Aps1*Aqs2 + Aqc1*Apc2);

     double AApcos[5],AApsin[5],AAccos[5],AAcsin[5];
       AApcos[0]=0.5*(Bp1c2+Bp2s2);
       AApsin[0]=0.5*(Bp2c2-Bp1s2);
       AAccos[0]=0.5*(Bc1c2+Bc2s2);
       AAcsin[0]=0.5*(Bc2c2-Bc1s2);
       AApcos[1]=0.5*(Bp1c1+Bp2s1);
       AApsin[1]=0.5*(Bp2c1-Bp1s1);
       AAccos[1]=0.5*(Bc1c1+Bc2s1);
       AAcsin[1]=0.5*(Bc2c1-Bc1s1);
       AApcos[2]=Bp1cn;
       AApsin[2]=Bp2cn;
       AAccos[2]=Bc1cn;
       AAcsin[2]=Bc2cn;
       AApcos[3]=0.5*(Bp1c1-Bp2s1);
       AApsin[3]=0.5*(Bp2c1+Bp1s1);
       AAccos[3]=0.5*(Bc1c1-Bc2s1);
       AAcsin[3]=0.5*(Bc2c1+Bc1s1);
       AApcos[4]=0.5*(Bp1c2-Bp2s2);
       AApsin[4]=0.5*(Bp2c2+Bp1s2);
       AAccos[4]=0.5*(Bc1c2-Bc2s2);
       AAcsin[4]=0.5*(Bc2c2+Bc1s2);
       int ki;

	 // int useBessel = 1;
	 /*if (useBessel){
		 std::cout << " testing negative order " << std::endl;
		 std::cout << " n=-1 " << gsl_sf_bessel_Jn(-1, 0.999184) << std::endl;
		 double Bss[5];
	 	 double Xa, Xb, Xc;
		 double Se = 0.999184;
	 	 // for (int i=1; i<=nH; i++){
	 	 // int state = gsl_sf_bessel_Jn_array(nn-2, nn+2, nn*Se, Bss);
		 int state;
		 for (int nn=1; nn<=5; nn++){
			 // state = gsl_sf_bessel_Jn_array(nn-2, nn+2, nn*Se, Bss);
			 Bss[0] = gsl_sf_bessel_Jn(nn-2, nn*Se);
			 Bss[1] = gsl_sf_bessel_Jn(nn-1, nn*Se);
			 Bss[2] = gsl_sf_bessel_Jn(nn, nn*Se);
			 Bss[3] = gsl_sf_bessel_Jn(nn+1, nn*Se);
			 Bss[4] = gsl_sf_bessel_Jn(nn+2, nn*Se);

			 std::cout << nn << "  " << nn*Se << "  " << Bss[0] << "  " << Bss[1] << "  " << Bss[2] << \
			 			"  " << Bss[3] << "  " << Bss[4] << std::endl;
		 	// std::cout << " nn = " << nn << "    " << ComputeHAmpl(0.999184, nn) << std::endl;
		 }

	 }*/

       // std::complex<double> test;
		// std::ofstream fout09("Data/Evolution.dat");
		// if (debug){
	  //   std::ofstream fout09("Data/Evolution.dat");
		// }
     for (int i=1; i<Nph; i++){
			 /* std::cout << i << "  out of " << Nph << "  time  " << tm_ph[i] <<  std::endl;
			  std::cout << "orbit: " << tm_ph[i] << "   " << eccev[i] << "   " << nuev[i] << "    " << phiev[i] <<
							 "    " << alpev[i] << "   " << gamev[i] << "   " << gamdotev[i] << "   " <<  alpdotev[i] << "   " <<
							 phiddotev[i] << "    "  << gamddotev[i] << "    " << alpddotev[i] << "    " << Ampev[i] << std::endl;
			 if (debug){
				  fout09  << std::setprecision(10) << tm_ph[i] << "   " << eccev[i] << "   " << nuev[i] << "    " << phiev[i] <<
									"    " << alpev[i] << "   " << gamev[i] << std::endl;
				}*/
        xi_in = tm_ph[i-1] - DM*cos(orbOm*tm_ph[i-1]-S.phS);
        xi_fin = tm_ph[i] - DM*cos(orbOm*tm_ph[i]-S.phS);
        dxi = xi_fin - xi_in;
        ind = 0;
        // loops over harmonics
        for (int j=0;j<5;j++) {
        	  nn=(double)(j+1);
        	  for(int jj=0; jj<5; jj++){
        	     mm = (double)jj-2.;
        	     om_in = nn*LISAWP_TWOPI*nuev[i-1] + 2.*gamdotev[i-1]  + mm*alpdotev[i-1];
        	     om_fin = nn*LISAWP_TWOPI*nuev[i] + 2.*gamdotev[i]  + mm*alpdotev[i];
              delom = om_fin - om_in;
              ind_low = (int) ceil(om_in/dOm);
              ind_up = (int) floor(om_fin/dOm);
              if (om_fin == (double)ind_up*dOm && om_fin != 0.){
                  // std::cout << "Fourier freq = harm freq\n" << om_fin << spr << (double)ind_up*dOm << std::endl;
                  ind_up = ind_up-1; // otherwise we count twice this bin
              }
              if (om_in != 0. && om_fin != 0.){
                 // loop over fourier bins between two values of harmonic
								 // std::cout << "harm  " << j << "  " << jj << "  ind  " << ind_low << "   " << ind_up << std::endl;
                 for (int ii=ind_low; ii<=ind_up; ii++){
                     	 Om = (double)ii * dOm;
						 if ( Om < Om_Nyq ){
		                     delta = (Om - om_in)/delom;
		                     eps = 1.-delta;
		                     T =  tm_ph[i]*delta + tm_ph[i-1]*eps;
		                     xi = T - DM*cos(orbOm*T - S.phS);
		                     delta = (xi - xi_in)/dxi;
		                     eps = 1.-delta;
		                     dom = (nn*phiddotev[i-1] + 2.*gamddotev[i-1]  + mm*alpddotev[i-1])*eps +
		                           (nn*phiddotev[i] + 2.*gamddotev[i]  + mm*alpddotev[i])*delta;
		                     faza = (nn*phiev[i-1] + 2.*gamev[i-1]  + mm*alpev[i-1])*eps +
		                               (nn*phiev[i] + 2.*gamev[i]  + mm*alpev[i])*delta;
		                     amp = Ampev[i-1]*eps + Ampev[i]*delta;
		                     ec = eccev[i-1]*eps +  eccev[i]*delta;

							 if (ec>=0.45){
								 // std::cout << "call for Bessels:  " << ec << "   " << j+1 << std::endl;
								 // hamp = -(j+1)*ComputeHAmpl(ec, j+1);
								 ComputeHAmpl(ec, j+1, &hamp, &Xb_a, &XXc);

							 }else{
								 ec2=ec*ec;
			                     ec3=ec2*ec;
			                 	 ec4=ec3*ec;
			                     ec5=ec4*ec;
			                     ec6=ec5*ec;
			                     ec7=ec6*ec;
		                     	switch (j+1) {
		               			   case 1:
		               				   hamp =3.*ec-1.625*ec3;
		               				   break;
		               			   case 2:
		               				   hamp =-4.+10.*ec2-5.75*ec4;
		               				   break;
		               			   case 3:
		               				   hamp =-9.*ec+21.375*ec3-15.046875*ec5;
		               				   break;
		               			   case 4:
		               				   hamp =-16.*ec2+40.*ec4-101.*ec6/3.;
		               				   break;
		               			   case 5:
		               				   hamp =-625.*ec3/24.+26875.*ec5/384.-210625.*ec7/3072.;
		               				   break;
		               			 }
							 }
		               		 // delta = (T - tm_ph[i-1])/dt_ph;
							 delta = (Om - om_in)/delom; // restoring original definition
		                     eps = 1.-delta;
		                     cFp = Xp[ind][i-1]*eps + Xp[ind][i]*delta;
		                     cFc = Xc[ind][i-1]*eps + Xc[ind][i]*delta;

		                     Apc = hamp*AApcos[harmms[jj]+2];
		         		     Aps = hamp*AApsin[harmms[jj]+2]; // should be "-" for strain
		         		     Acc = 2.*hamp*AAccos[harmms[jj]+2];
		         		     Acs = 2.*hamp*AAcsin[harmms[jj]+2];

		         		     sinph = sin(faza - Om*xi + LISAWP_PI*0.25);
		                     cosph = cos(faza - Om*xi + LISAWP_PI*0.25);

		                     fact = 0.5*amp*sqrt(LISAWP_TWOPI/dom);
		                     //hplus = Apc - img*Aps;
		                     //hcross = Acc - img*Acs;
		                    // Xf[ii] = Xf[ii] + fact*(cFp*hplus + cFc*hcross)*(cosph + img*sinph);
		                     ki = Nps - ii;
		                     //xp = fact*hplus*(cosph + img*sinph);
		                     //xc = fact*hcross*(cosph + img*sinph);
		                     xp = fact*(Apc*cosph + Aps*sinph + img*(Apc*sinph - Aps*cosph));
		                     xc = fact*(Acc*cosph + Acs*sinph + img*(Acc*sinph - Acs*cosph));
		                     // TODO output was changed
		                     // -------------------------------------------------------------
		                     //Xf(ii) += (cFp*xp + cFc*xc).real();
		                     //Xf(ki) += (cFp*xp + cFc*xc).imag();
		                     Xf_r[ii] += (cFp*xp + cFc*xc).real();
		                     Xf_im[ii] += (cFp*xp + cFc*xc).imag();
		                     cFp = Yp[ind][i-1]*eps + Yp[ind][i]*delta;
		                     cFc = Yc[ind][i-1]*eps + Yc[ind][i]*delta;
		                     Yf_r[ii] += (cFp*xp + cFc*xc).real();
		                     Yf_im[ii] += (cFp*xp + cFc*xc).imag();
		                     //Yf(ii) += (cFp*xp + cFc*xc).real();
		                     //Yf(ki) += (cFp*xp + cFc*xc).imag();
		                     cFp = Zp[ind][i-1]*eps + Zp[ind][i]*delta;
		                     cFc = Zc[ind][i-1]*eps + Zc[ind][i]*delta;
		                     Zf_r[ii] += (cFp*xp + cFc*xc).real();
		                     Zf_im[ii] += (cFp*xp + cFc*xc).imag();
		                     //Zf(ii) += (cFp*xp + cFc*xc).real();
		                     //Zf(ki) += (cFp*xp + cFc*xc).imag();
						}// in Nyquist
                  } // ii-loop
              } // another if
              ind ++;
        	  }// jj loop
     	  } //j loop
     }

     // fout09.close();
}

void FreqAK_RA::ConstructWaveFreq_RAv2(EMRItemplate& S, double* Xf_r, double* Xf_im,  double* Yf_r, double* Yf_im,
                        double* Zf_r, double* Zf_im)

{
	int debug =0;
	double AUsec =  LISAWP_AU_SI/LISAWP_C_SI;
	double orbOm = LISAWP_TWOPI/LISAWP_YRSID_SI;
	double DM = AUsec*S.stS;
	double xi_in, xi_fin;
	int nn, ll, mm;

	if (debug>0){
		double cS= S.ctS;
        double sS= S.stS;
        double cK= S.ctK;
        double sK= S.stK;
        double cSp= S.cpS;
        double sSp= S.spS;
        double cKp= S.cpK;
        double sKp= S.spK;
        double clam=cos(S.lam);
        double slam=sin(S.lam);
        double Sn = cS*cK + sS*sK*cos(S.phS  - S. phK);
        double cX = Sn;

        double sX = sqrt( sS*sS*cK*cK - 2.*sS*sSp*cK*cS*sK*sKp + \
                         cS*cS*sK*sK - 2.*cS*sK*cKp*sS*cSp*cK + \
                         sS*sS*cSp*cSp*sK*sK*sKp*sKp - 2.*sS*sS*cSp*sK*sK*sKp*sSp*cKp +\
                         sS*sS*sSp*sSp*sK*sK*cKp*cKp);

        double Apc1, Aps1, Apcn1, Apc2, Aps2, Apcn2, Aqc1, Aqs1, Aqcn1, Aqc2, Aqs2, Aqcn2;

          Apc1 = ( -cK*cKp*sS*cSp - cK*sKp*sS*sSp + sK*cS )/(sX);
          Aps1 = ( sKp*sS*cSp - cKp*sS*sSp )/(sX);
          Apcn1 = ( cK*cS + sK*cKp*sS*cSp + sK*sKp*sS*sSp - cX)*clam/(sX*slam);

          Apc2 = (sS*cSp*sKp - sS*sSp*cKp )*clam/(sX);
          Aps2 = ( cK*cKp*sS*cSp + cK*sKp*sS*sSp - cS*sK )*clam/(sX);
          Apcn2 = 0.0;


          Aqc1 = ( sS*cSp*sKp - sS*sSp*cKp  )*cX/(sX);
          Aqs1 = ( cK*cKp*sS*cSp + cK*sKp*sS*sSp - cS*sK )*cX/(sX);
          Aqcn1 = 0.0;


          Aqc2 = cX*clam*( cK*cKp*sS*cSp + cK*sKp*sS*sSp - sK*cS)/(sX);
          Aqs2 = -cX*clam*sS*( sKp*cSp - cKp*sSp )/(sX);
          Aqcn2 = -( cX*clam*clam*( cK*cS + sK*cKp*sS*cSp + sK*sKp*sS*sSp ) + \
                                   1.- cX*cX - clam*clam )/(sX*slam);


        double Bp1c1 = 2.0*(Apc1*Apcn1 - Aqc1*Aqcn1 + Aqc2*Aqcn2 - Apc2*Apcn2);
        double Bp1c2 =  0.5*(Aps2*Aps2 - Aqc1*Aqc1  + Apc1*Apc1  - Aps1*Aps1 + \
                                   Aqc2*Aqc2 + Aqs1*Aqs1 - Apc2*Apc2 - Aqs2*Aqs2);
        double Bp1s1 = 2.0*(Aqs2*Aqcn2 - Aps2*Apcn2 - Aqs1*Aqcn1 + Aps1*Apcn1);
        double Bp1s2 = (Apc1*Aps1 + Aqc2*Aqs2 - Apc2*Aps2 - Aqc1*Aqs1);
        double Bp1cn = 0.5*(Apc1*Apc1 + Aps1*Aps1 - Aqc1*Aqc1 - Aqs1*Aqs1 - Apc2*Apc2 \
                                   + Aqc2*Aqc2 + Aqs2*Aqs2 - Aps2*Aps2) + Aqcn2*Aqcn2 - Aqcn1*Aqcn1 \
                                   + Apcn1*Apcn1 - Apcn2*Apcn2;

        double Bp2c1 = (Apcn1*Apc2 + Apc1*Apcn2 - Aqcn1*Aqc2 - Aqc1*Aqcn2);
        double Bp2c2 = 0.5*(Aqs1*Aqs2 - Aps1*Aps2 + Apc1*Apc2 - Aqc1*Aqc2);
        double Bp2s1 = (Aps1*Apcn2 + Apcn1*Aps2 - Aqcn1*Aqs2 - Aqs1*Aqcn2);
        double Bp2s2 = 0.5*( Apc1*Aps2 - Aqc1*Aqs2 + Aps1*Apc2 - Aqs1*Aqc2);
        double Bp2cn = 0.5*(Aps1*Aps2 - Aqs1*Aqs2 - Aqc1*Aqc2 + Apc1*Apc2) -Aqcn1*Aqcn2 + Apcn1*Apcn2;

        double Bc1c1 = (-Apc2*Aqcn2 - Apcn2*Aqc2 + Apc1*Aqcn1 + Apcn1*Aqc1);
        double Bc1c2 = 0.5*( Apc1*Aqc1 - Aps1*Aqs1 - Apc2*Aqc2 + Aps2*Aqs2);
        double Bc1s1 = (Apcn1*Aqs1 - Aps2*Aqcn2 + Aps1*Aqcn1 - Apcn2*Aqs2);
        double Bc1s2 = 0.5*(-Apc2*Aqs2 + Apc1*Aqs1 - Aps2*Aqc2 + Aps1*Aqc1);
        double Bc1cn = -Apcn2*Aqcn2 + Apcn1*Aqcn1 + 0.5*(Apc1*Aqc1 - Aps2*Aqs2 + Aps1*Aqs1 - Apc2*Aqc2);

        double Bc2c1 = (Aqc1*Apcn2 + Aqcn1*Apc2 + Apc1*Aqcn2 + Apcn1*Aqc2);
        double Bc2c2 = 0.5*( Apc1*Aqc2 - Aps1*Aqs2 + Aqc1*Apc2 - Aqs1*Aps2);
        double Bc2s1 = (Apcn1*Aqs2 + Aqs1*Apcn2 + Aps1*Aqcn2 + Aqcn1*Aps2);
        double Bc2s2 = 0.5*(Aqc1*Aps2 + Apc1*Aqs2 + Aqs1*Apc2 + Aps1*Aqc2);
        double Bc2cn = Aqcn1*Apcn2 + Apcn1*Aqcn2 + 0.5*(Apc1*Aqc2 + Aqs1*Aps2 +Aps1*Aqs2 + Aqc1*Apc2);

        double AApcos[5],AApsin[5],AAccos[5],AAcsin[5];
         AApcos[0]=0.5*(Bp1c2+Bp2s2);
          AApsin[0]=0.5*(Bp2c2-Bp1s2);
          AAccos[0]=0.5*(Bc1c2+Bc2s2);
          AAcsin[0]=0.5*(Bc2c2-Bc1s2);
          AApcos[1]=0.5*(Bp1c1+Bp2s1);
          AApsin[1]=0.5*(Bp2c1-Bp1s1);
          AAccos[1]=0.5*(Bc1c1+Bc2s1);
          AAcsin[1]=0.5*(Bc2c1-Bc1s1);
          AApcos[2]=Bp1cn;
          AApsin[2]=Bp2cn;
          AAccos[2]=Bc1cn;
          AAcsin[2]=Bc2cn;
          AApcos[3]=0.5*(Bp1c1-Bp2s1);
          AApsin[3]=0.5*(Bp2c1+Bp1s1);
          AAccos[3]=0.5*(Bc1c1-Bc2s1);
          AAcsin[3]=0.5*(Bc2c1+Bc1s1);
          AApcos[4]=0.5*(Bp1c2-Bp2s2);
          AApsin[4]=0.5*(Bp2c2+Bp1s2);
          AAccos[4]=0.5*(Bc1c2-Bc2s2);
          AAcsin[4]=0.5*(Bc2c2+Bc1s2);

		  for (int i=0; i<5; i++){
			  std::cout << i << "   " << AApcos[i] << "   " << AApsin[i] << "   " <<
			   AAccos[i] << "   " << AAcsin[i] << std::endl;
		  }
	}

	std::complex<double> Ap[5];
	std::complex<double> Ac[5];
	std::complex<double> Bp[5];
	std::complex<double> Bc[5];

	ComputeOrbPrecAmpls(S, Ap, Ac, Bp, Bc);
	if (debug){
		for (int i=0; i<5; i++){
			std::cout << i<<" to compare with " << Ap[i] << "  " << Ac[i] << " ,  " << Bp[i] << "   " << Bc[i] << std::endl;
		}
	}

	for (int i=1; i<Nph; i++){
		// time in the LISA' frame
		xi_in = tm_ph[i-1] - DM*cos(orbOm*tm_ph[i-1]-S.phS);
        xi_fin = tm_ph[i] - DM*cos(orbOm*tm_ph[i]-S.phS);


		for (int j=0;j<5;j++) {
        	  nn=(double)(j+1);
        	  for(int jj=0; jj<5; jj++){
        	     mm = (double)jj-2.;
				 // dominant modes
				 ll = 2;
				 InterpolateWaveform(nn, ll, mm, i, S, xi_in, xi_fin, Ap, Ac, Xf_r, Xf_im,
						Yf_r, Yf_im, Zf_r, Zf_im);
				 if (Ext >0){
					 ll = 0;
					 InterpolateWaveform(nn, ll, mm, i, S, xi_in, xi_fin, Bp, Bc, Xf_r, Xf_im,
							Yf_r, Yf_im, Zf_r, Zf_im);
				 }
			 }// end of jj-loop
		 }// end of j loop
	 } // end of orb evol loop
// double* Apcos, double* Accos, double* Apsin, double* Acsin
}


void FreqAK_RA::ComputeOrbPrecAmpls(EMRItemplate& S, std::complex<double>* Ap, std::complex<double>* Ac,
						std::complex<double>* Bp, std::complex<double>* Bc )
{

	// double* AApcos, double* AApsin, double* AAccos, double* AAcsin){

	// def pe1pe1():
	double lam = S.lam;
	double th = S.thS;
	double phi = S.phS;
	double thK = S.thK;
	double phiK = S.phK;
	double c_chi = cos(th)*cos(thK) + sin(th)*sin(thK)*cos(phi  - phiK);
	double chi = acos(c_chi);
	double pe1pe1[5];
	double pe2pe2[5];
	double qe1qe1[5];
	double qe2qe2[5];
	double pe1pe2[5];
	double qe1qe2[5];
	double pe1qe1[5];
	double pe2qe2[5];
	double pe1qe2[5];
	double qe1pe2[5];
	std::complex<double> img(0.0, 1.0);
	std::complex<double> B1[5];
	std::complex<double> B2[5];
	std::complex<double> B3[5];
	std::complex<double> B4[5];
	std::complex<double> B5[5];
	std::complex<double> B6[5];



	pe1pe1[0] =   pow(1./tan(lam), 2)  * pow( 1/tan(chi) - 1./sin(chi) * \
	        ( cos(th) * cos(thK) + cos(phi - phiK) * sin(th)* sin(thK) ) , 2) +\
	        0.5 * pow(1./sin(chi),2) * ( pow(cos(thK) * cos(phi - phiK) * sin(th) - \
	        cos(th) * sin(thK), 2)  + pow(sin(th) * sin(phi-phiK),2)  );

	pe1pe1[1] = -2. * (1./tan(lam)) * pow(1./sin(chi),2) * ( cos(thK) * cos( phi - phiK ) * \
	        sin(th) - cos(th)*sin(thK) ) * ( cos(th) * cos(thK) - cos(chi) + \
	        cos(phi-phiK) * sin(th)* sin(thK) );

	pe1pe1[2] = 0.5 * pow(1/sin(chi),2) * ( pow(cos(thK)*cos(phi - phiK)*sin(th) - \
	        cos(th)*sin(thK), 2)  - pow(sin(th) *sin(phi - phiK), 2)  );

	pe1pe1[3] = -2./tan(lam) * pow(1./sin(chi),2) * sin(th) * ( cos(th)*cos(thK) - \
	        cos(chi) + cos(phi-phiK)*sin(th)*sin(thK) ) * sin(phi-phiK);

	pe1pe1[4] = 0.5*pow(1./sin(chi),2)  * ( -sin(2.*th) * sin(thK) * sin(phi-phiK) +\
	        cos(thK) * pow(sin(th),2)  * sin(2.*(phi-phiK)) );


	// def pe2pe2():

	pe2pe2[0] = 0.5*pow(cos(lam)*(1./sin(chi)),2)  * ( pow( cos(thK) * cos(phi-phiK) * sin(th) - \
	        cos(th)*sin(thK) ,2) + pow(sin(th) * sin(phi-phiK), 2)  );

	pe2pe2[2] = 0.5*pow(cos(lam)*(1./sin(chi)),2) * ( -pow(cos(thK) * cos(phi-phiK) * sin(th) - \
	        cos(th)*sin(thK) ,2)  + pow(sin(th) * sin(phi-phiK), 2)  );

	pe2pe2[1] = 0.0;

	pe2pe2[3] = 0.0;

	pe2pe2[4] = 0.5*pow(cos(lam) * (1./sin(chi)), 2)  * ( sin(2*th)*sin(thK)*sin(phi-phiK) -\
	        cos(thK)*pow(sin(th),2) * sin(2.*(phi-phiK)) );


	// def qe1qe1():

	qe1qe1[0] = 0.5*pow(1./tan(chi), 2) * ( pow( cos(thK) * cos(phi-phiK) * sin(th) - \
	        cos(th) * sin(thK) ,2) + pow(sin(th) * sin(phi-phiK), 2) );

	qe1qe1[2] = 0.5*pow(1./tan(chi), 2) * ( -pow( cos(thK) * cos(phi-phiK) * sin(th) - \
	        cos(th) * sin(thK) ,2) + pow(sin(th) * sin(phi-phiK), 2) );

	qe1qe1[1] = 0.0;
	qe1qe1[3] = 0.0;

	qe1qe1[4] = 0.5*pow(1/tan(chi),2) * ( sin(2*th) * sin(thK) * sin(phi-phiK) - \
	        cos(thK) * pow(sin(th),2) * sin(2*(phi-phiK)) );


	// def qe2qe2():

	qe2qe2[0] = 0.5*pow(cos(lam)*(1./tan(chi)), 2) * ( pow(cos(thK) * cos(phi-phiK)*sin(th) -\
	        cos(th) * sin(thK) ,2) + pow(sin(th) * sin(phi-phiK), 2)  ) \
	        + pow( cos(lam)*(1./tan(lam)) * ( cos(th)*cos(thK)*(1./tan(chi)) - \
	        1./sin(chi) + cos(phi-phiK) * (1./tan(chi)) *sin(th)*sin(thK) ) + \
	        (1./sin(lam))*sin(chi) ,2);

	qe2qe2[2] = 0.5*pow(cos(lam)*(1./tan(chi)), 2) * ( pow(cos(thK) * cos(phi-phiK)*sin(th) -\
	        cos(th) * sin(thK) ,2) - pow(sin(th) * sin(phi-phiK),2) );

	qe2qe2[1] = -2./tan(lam) * 1./tan(chi) * 1./sin(chi) * ( cos(thK)*cos(phi-phiK)*sin(th) -\
	        cos(th)*sin(thK) ) * ( pow(cos(lam),2)  * ( -1 + cos(th)*cos(thK)*cos(chi) + \
	        cos(phi-phiK)*cos(chi)*sin(th)*sin(thK) ) + pow(sin(chi),2)  );

	qe2qe2[3] = -2./tan(lam) * 1./tan(chi) * 1./sin(chi) * sin(th)*sin(phi-phiK) * \
	        ( pow(cos(lam),2) * ( -1 + cos(th)*cos(thK)*cos(chi) + \
	        cos(phi-phiK)* cos(chi)*sin(th)*sin(thK) ) + pow(sin(chi),2) );

	qe2qe2[4] =  pow(cos(lam) * (1./tan(chi)), 2) * sin(th)*( cos(thK)*cos(phi-phiK)*sin(th) -\
	        cos(th)*sin(thK) )*sin(phi-phiK);


	// def pe1pe2():

	pe1pe2[0] = 0.0;

	pe1pe2[2] = cos(lam) * pow(1./sin(chi),2) * sin(th)*( -cos(thK)*cos(phi-phiK)*sin(th) +\
	        cos(th)*sin(thK) ) * sin(phi-phiK);

	pe1pe2[1] = cos(lam) * 1./tan(lam) * pow(1./sin(chi),2) * sin(th) * ( cos(th)*cos(thK) - \
	        cos(chi) + cos(phi-phiK)*sin(th)*sin(thK) ) * sin(phi-phiK);

	pe1pe2[3] = cos(lam) * 1./tan(lam) * pow(1./sin(chi),2) * ( cos(thK)*cos(phi-phiK)*sin(th) -\
	        cos(th)*sin(thK) ) * ( -cos(th)*cos(thK) + cos(chi) - cos(phi-phiK)*sin(th)*sin(thK) );

	pe1pe2[4] = 1./8. * cos(lam) * pow(1./sin(chi),2) * ( (3. + cos(2.*thK)) * cos(2.*(phi-phiK))*pow(sin(th),2) +\
	        (1. + 3.*cos(2.*th)) * pow(sin(thK),2) - 2.*cos(phi-phiK)*sin(2.*th)*sin(2.*thK) );


	// def qe1qe2():

	qe1qe2[0] = 0.0;

	qe1qe2[2] = cos(lam)*pow(1./tan(chi),2) * sin(th)*( cos(thK)*cos(phi-phiK)*sin(th) -\
	        cos(th)*sin(thK) ) * sin(phi-phiK);

	qe1qe2[1] = -1./tan(chi) * 1./sin(lam) * 1./sin(chi) * sin(th) * sin(phi-phiK) *\
	        ( pow(cos(lam),2) * ( -1 + cos(th)*cos(thK)*cos(chi) + \
	        cos(phi-phiK)*cos(chi)*sin(th)*sin(thK) ) + pow(sin(chi),2) );

	qe1qe2[3] = -1./tan(chi) * 1./sin(lam) * 1./sin(chi) * ( -cos(thK)*cos(phi-phiK)*sin(th) +\
	        cos(th)*sin(thK) ) * ( pow(cos(lam),2) * (-1. + cos(th)*cos(thK)*cos(chi) + \
	        cos(phi-phiK)*cos(chi)*sin(th)*sin(thK) ) + pow(sin(chi),2)  );

	qe1qe2[4] =  -1./8. * cos(lam) * pow(1./tan(chi),2) * ( (3. + cos(2.*thK)) * cos(2.*(phi-phiK))*pow(sin(th),2) +\
	        (1. + 3.*cos(2.*th)) * pow(sin(thK),2) - 2.*cos(phi-phiK)*sin(2.*th)*sin(2.*thK) );


	// def pe1qe1():

	pe1qe1[0] = 0.0;

	pe1qe1[2] = 1./tan(chi) * 1./sin(chi) * sin(th) * ( cos(thK)*cos(phi-phiK)*sin(th) -\
	        cos(th)*sin(thK) ) * sin(phi-phiK);

	pe1qe1[1] = 1./tan(lam) * 1./tan(chi) * 1./sin(chi) * sin(th) *\
	        ( -cos(th)*cos(thK) + cos(chi) - cos(phi-phiK)*sin(th)*sin(thK) )* sin(phi-phiK);

	pe1qe1[3] = 1/tan(lam) * 1./tan(chi) * ( cos(thK)*cos(phi-phiK)*sin(th) - cos(th)*sin(thK) ) *\
	        ( -1./tan(chi) + 1./sin(chi) * ( cos(th)*cos(thK) + cos(phi-phiK)*sin(th)*sin(thK) ) );

	pe1qe1[4] = -1./8. * 1/tan(chi) * 1./sin(chi) * ( (3. + cos(2.*thK)) * cos(2.*(phi-phiK))*pow(sin(th),2) +\
	        (1. + 3.*cos(2.*th)) * pow(sin(thK),2) - 2.*cos(phi-phiK)*sin(2.*th)*sin(2.*thK) );


	// def pe2qe2():

	pe2qe2[0] = 0.0;

	pe2qe2[2] = pow(cos(lam),2) * 1./tan(chi) * 1./sin(chi) * sin(th) * ( -cos(thK)*cos(phi-phiK)*sin(th) +\
	        cos(th)*sin(thK) ) * sin(phi-phiK);

	pe2qe2[1] = 1./tan(lam) * pow(1./sin(chi),2)  * sin(th)*sin(phi-phiK)* \
	    ( pow(cos(lam),2) * ( -1. + cos(th)*cos(thK)*cos(chi) + \
	    cos(phi-phiK)*cos(chi)*sin(th)* sin(thK) ) + pow(sin(chi),2)  );

	pe2qe2[3] = 1./tan(lam) * pow(1./sin(chi),2) * ( -cos(thK)*cos(phi-phiK)*sin(th) + cos(th)*sin(thK) ) *\
	        ( pow(cos(lam),2) * ( -1. + cos(th)*cos(thK)*cos(chi) + \
	        cos(phi-phiK)*cos(chi)*sin(th)* sin(thK) ) + pow(sin(chi),2)  );

	pe2qe2[4] = 1./8. * pow(cos(lam),2) * 1/tan(chi) * 1./sin(chi) * ( (3. + cos(2.*thK)) * cos(2.*(phi-phiK))*pow(sin(th),2) +\
	        (1. + 3.*cos(2.*th)) * pow(sin(thK),2) - 2.*cos(phi-phiK)*sin(2.*th)*sin(2.*thK) );


	// def pe1qe2():

	pe1qe2[0] =  0.5/sin(chi) * ( cos(lam) * 1./tan(chi) * pow( cos(thK)*cos(phi-phiK)*sin(th) -\
	        cos(th)*sin(thK) ,2) + cos(lam)*1./tan(chi)*pow(sin(th) * sin(phi-phiK),2) +\
	        2./tan(lam) * 1./sin(lam) * 1./sin(chi) * ( -cos(th)*cos(thK) + cos(chi) - \
	        cos(phi-phiK)*sin(th)*sin(thK) ) * ( -pow(cos(lam),2) * ( -1. + cos(th)*cos(thK)*cos(chi) +\
	        cos(phi-phiK)*cos(chi)*sin(th)*sin(thK) )  - pow(sin(chi),2) ) );

	pe1qe2[2] =  0.5 * cos(lam)*1./tan(chi) * 1./sin(chi) * ( pow(cos(thK)*cos(phi-phiK)*sin(th) -\
	        cos(th) * sin(thK) ,2) - pow(sin(th) * sin(phi-phiK),2) );

	pe1qe2[1] = -1./sin(lam) * pow(1./sin(chi),2) * ( cos(thK)*cos(phi-phiK)*sin(th) -\
	        cos(th)*sin(thK) ) * ( -pow(cos(lam),2)  * ( 1. + pow(cos(chi),2) - 2.*cos(chi) *\
	         ( cos(th)*cos(thK) + cos(phi -phiK)*sin(th)*sin(thK) ) )  + pow(sin(chi),2)  );

	pe1qe2[3] = -1./sin(lam) * pow(1./sin(chi),2) * sin(th) * sin(phi-phiK) *\
	        ( -pow(cos(lam),2) * ( 1. + pow(cos(chi),2) - 2.*cos(chi)*( cos(th)*cos(thK) +\
	         cos(phi-phiK)*sin(th)*sin(thK) ) ) + pow(sin(chi),2) );

	pe1qe2[4] = 0.5*cos(lam) * 1./tan(chi) * 1./sin(chi) *( -sin(2*th)*sin(thK)*sin(phi-phiK) + \
	        cos(thK) * pow(sin(th),2) * sin(2*(phi-phiK)) );


	// def qe1pe2():

	qe1pe2[0] = 0.5 * cos(lam) * 1./tan(chi) * 1./sin(chi) * ( -pow( cos(thK)*cos(phi-phiK)*sin(th) - \
	        cos(th)*sin(thK) ,2) - pow(sin(th) * sin(phi-phiK),2)  );

	qe1pe2[2] = 0.5*cos(lam) * 1./tan(chi) * 1./sin(chi) * ( pow( cos(thK)*cos(phi-phiK)*sin(th) - \
	        cos(th)*sin(thK) ,2) - pow(sin(th) * sin(phi-phiK),2)  );

	qe1pe2[1] = 0.0;

	qe1pe2[3] = 0.0;

	qe1pe2[4] = 0.5*cos(lam) * 1./tan(chi) * 1./sin(chi)* ( -sin(2.*th) * sin(thK) * sin(phi-phiK) +\
	        cos(thK) * pow(sin(th),2)  * sin(2.*(phi-phiK)) );


	// had to do it (historically)
	double e1pe1p[5], e1qe1q[5], e2pe2p[5], e2qe2q[5];
	double e1pe2p[5], e1qe2q[5], e1pe1q[5], e2pe2q[5];
	double e1pe2q[5], e1qe2p[5];

	for (int i=0; i<5; i++){
		e1pe1p[i] = pe1pe1[i];
		e1qe1q[i] = qe1qe1[i];
		e2pe2p[i] = pe2pe2[i];
		e2qe2q[i] = qe2qe2[i];
		e1pe2p[i] = pe1pe2[i];
		e1qe2q[i] = qe1qe2[i];
		e1pe1q[i] = pe1qe1[i];
		e2pe2q[i] = pe2qe2[i];
		e1pe2q[i] = pe1qe2[i];
		e1qe2p[i] = qe1pe2[i];
		// std::cout << "e1pe1p  " << e1pe1p[i] << std::endl;
		// std::cout << "e2pe2p  " << e2pe2p[i] << std::endl;
	}

	B1[2] = e1pe1p[0] - e1qe1q[0] - e2pe2p[0] + e2qe2q[0];
	B1[1] = 0.5*( e1pe1p[1] - e1qe1q[1] - e2pe2p[1] + e2qe2q[1]) - \
		        0.5*img* (e1pe1p[3] - e1qe1q[3] - e2pe2p[3] + e2qe2q[3]);
	B1[0] =  0.5*( e1pe1p[2] - e1qe1q[2] - e2pe2p[2] + e2qe2q[2]) - \
		        0.5*img* (e1pe1p[4] - e1qe1q[4] - e2pe2p[4] + e2qe2q[4]);
	B1[3] = 0.5*( e1pe1p[1] - e1qe1q[1] - e2pe2p[1] + e2qe2q[1]) + \
		        0.5*img* (e1pe1p[3] - e1qe1q[3] - e2pe2p[3] + e2qe2q[3]);
	B1[4] =  0.5*( e1pe1p[2] - e1qe1q[2] - e2pe2p[2] + e2qe2q[2]) + \
		        0.5*img* (e1pe1p[4] - e1qe1q[4] - e2pe2p[4] + e2qe2q[4]);


	B2[2] = e1pe1p[0] - e1qe1q[0] + e2pe2p[0] - e2qe2q[0];
	B2[1] = 0.5*(e1pe1p[1] - e1qe1q[1] + e2pe2p[1] - e2qe2q[1]) -\
		        0.5*img*(e1pe1p[3] - e1qe1q[3] + e2pe2p[3] - e2qe2q[3]);
	B2[0] = 0.5*(e1pe1p[2] - e1qe1q[2] + e2pe2p[2] - e2qe2q[2]) -\
		        0.5*img*(e1pe1p[4] - e1qe1q[4] + e2pe2p[4] - e2qe2q[4]);
	B2[3] = 0.5*(e1pe1p[1] - e1qe1q[1] + e2pe2p[1] - e2qe2q[1]) +\
		        0.5*img*(e1pe1p[3] - e1qe1q[3] + e2pe2p[3] - e2qe2q[3]);
	B2[4] = 0.5*(e1pe1p[2] - e1qe1q[2] + e2pe2p[2] - e2qe2q[2]) +\
		        0.5*img*(e1pe1p[4] - e1qe1q[4] + e2pe2p[4] - e2qe2q[4]);


	B3[2] = e1pe2p[0] - e1qe2q[0];
	B3[1] = 0.5*(e1pe2p[1] - e1qe2q[1]) -0.5*img*(e1pe2p[3] - e1qe2q[3]);
	B3[0] = 0.5*(e1pe2p[2] - e1qe2q[2]) -0.5*img*(e1pe2p[4] - e1qe2q[4]);
	B3[3] = 0.5*(e1pe2p[1] - e1qe2q[1]) +0.5*img*(e1pe2p[3] - e1qe2q[3]);
	B3[4] = 0.5*(e1pe2p[2] - e1qe2q[2]) +0.5*img*(e1pe2p[4] - e1qe2q[4]);

	B4[2] = e1pe1q[0] - e2pe2q[0];
	B4[1] = 0.5*(e1pe1q[1] - e2pe2q[1]) - 0.5*img*(e1pe1q[3] - e2pe2q[3]);
	B4[0] = 0.5*(e1pe1q[2] - e2pe2q[2]) - 0.5*img*(e1pe1q[4] - e2pe2q[4]);
	B4[3] = 0.5*(e1pe1q[1] - e2pe2q[1]) + 0.5*img*(e1pe1q[3] - e2pe2q[3]);
	B4[4] = 0.5*(e1pe1q[2] - e2pe2q[2]) + 0.5*img*(e1pe1q[4] - e2pe2q[4]);

	B5[2] = e1pe1q[0] + e2pe2q[0];
	B5[1] = 0.5*(e1pe1q[1] + e2pe2q[1]) -  0.5*img*(e1pe1q[3] + e2pe2q[3]);
	B5[0] = 0.5*(e1pe1q[2] + e2pe2q[2]) -  0.5*img*(e1pe1q[4] + e2pe2q[4]);
	B5[3] = 0.5*(e1pe1q[1] + e2pe2q[1]) +  0.5*img*(e1pe1q[3] + e2pe2q[3]);
	B5[4] = 0.5*(e1pe1q[2] + e2pe2q[2]) +  0.5*img*(e1pe1q[4] + e2pe2q[4]);

	B6[2] = e1pe2q[0] + e1qe2p[0];
	B6[1] = 0.5*(e1pe2q[1] + e1qe2p[1]) - 0.5*img*(e1pe2q[3] + e1qe2p[3]);
	B6[0] = 0.5*(e1pe2q[2] + e1qe2p[2]) - 0.5*img*(e1pe2q[4] + e1qe2p[4]);
	B6[3] = 0.5*(e1pe2q[1] + e1qe2p[1]) + 0.5*img*(e1pe2q[3] + e1qe2p[3]);
	B6[4] = 0.5*(e1pe2q[2] + e1qe2p[2]) + 0.5*img*(e1pe2q[4] + e1qe2p[4]);


	// Ap =  np.conjugate(Bmy[0, :]) + 1.j*np.conjugate(Bmy[2,:])
	for (int i=0; i<5; i++){
		Ap[i] = conj(B1[i]) + img*conj(B3[i]);
		Ac[i] = conj(B4[i]) + img*conj(B6[i]);
		Bp[i] = conj(B2[i]);
		Bc[i] = conj(B5[i]);
	}
	// Ac =  np.conjugate(Bmy[3, :] + 1.j*Bmy[5,:]))
}

// void FreqAK_RA::InterpolateWaveform(int nn, int ll, int mm, EMRItemplate& S, double xi_in, double xi_fin,
// 				double* Apcos, double* Accos, double* Apsin, double* Acsin, double* Xf_r, double* Xf_im,
// 				double* Yf_r, double* Yf_im, double* Zf_r, double* Zf_im){

void FreqAK_RA::InterpolateWaveform(int nn, int ll, int mm, int ci, EMRItemplate& S, double xi_in, double xi_fin,
				std::complex<double>* Ap, std::complex<double>* Ac, double* Xf_r, double* Xf_im,
				double* Yf_r, double* Yf_im, double* Zf_r, double* Zf_im)
{

	double  xi;
	// , Apc, Aps, Acc, Acs;
	// int harmms[]={-2,-1,0,1,2};
	double fact;
	std::complex<double> cFp;
	std::complex<double> cFc;
	std::complex<double> x, xp, xc;
	std::complex<double> img(0.0, 1.0);
	double om_in;
	double om_fin;
	double dOm = df*LISAWP_TWOPI;
	double Om, dom, delom;
	double Om_Nyq = LISAWP_TWOPI*0.5/dt_w; //we limit freq. content of the waveform to this freq.
	double dxi;
	double orbOm = LISAWP_TWOPI/LISAWP_YRSID_SI;
	double AUsec =  LISAWP_AU_SI/LISAWP_C_SI;
	double DM = AUsec*S.stS;
	double delta, eps;
	double T, ec;
	double amp, faza, sinph, cosph, hamp;
	double XaXb, Xa_Xb, Xic;
	int ind, ind_up, ind_low;


	// I assume that m starts from -2 -> 2, n: 1->5
	ind = (nn-1)*5 + (mm+2);
	dxi = xi_fin - xi_in;
	om_in = nn*LISAWP_TWOPI*nuev[ci-1] + ll*gamdotev[ci-1]  + mm*alpdotev[ci-1];
	om_fin = nn*LISAWP_TWOPI*nuev[ci] + ll*gamdotev[ci]  + mm*alpdotev[ci];
	delom = om_fin - om_in;
	ind_low = (int) ceil(om_in/dOm);
	ind_up = (int) floor(om_fin/dOm);
	if (om_fin == (double)ind_up*dOm && om_fin != 0.){
		ind_up = ind_up-1; // otherwise we count twice this bin
	}
	if (om_in != 0. && om_fin != 0.){
	   // loop over fourier bins between two values of harmonic
					   // std::cout << "harm  " << j << "  " << jj << "  ind  " << ind_low << "   " << ind_up << std::endl;
	   for (int ii=ind_low; ii<=ind_up; ii++){
		   Om = (double)ii * dOm;
		   if ( Om < Om_Nyq ){
			   delta = (Om - om_in)/delom;
			   eps = 1.-delta;
			   T =  tm_ph[ci]*delta + tm_ph[ci-1]*eps;
			   xi = T - DM*cos(orbOm*T - S.phS);
			   delta = (xi - xi_in)/dxi;
			   eps = 1.-delta;
			   dom = (nn*phiddotev[ci-1] + ll*gamddotev[ci-1]  + mm*alpddotev[ci-1])*eps +
					 (nn*phiddotev[ci] + ll*gamddotev[ci]  + mm*alpddotev[ci])*delta;
			   faza = (nn*phiev[ci-1] + ll*gamev[ci-1]  + mm*alpev[ci-1])*eps +
						 (nn*phiev[ci] + ll*gamev[ci]  + mm*alpev[ci])*delta;
			   amp = Ampev[ci-1]*eps + Ampev[ci]*delta;
			   ec = eccev[ci-1]*eps +  eccev[ci]*delta;


			   ComputeHAmpl(ec, nn, &XaXb, &Xa_Xb, &Xic);
			   if (ll == 2){
				   hamp = XaXb;
			   }
			   if (ll == 0){
				   hamp = Xic;
				   // std::cout << "ll = 0  " << "nn = " << nn << "   ecc  " << ec << "   " << Om << "  " << hamp << std::endl;
			   }
			   if (ll == -2){
				   hamp = Xa_Xb;
			   }

			   delta = (Om - om_in)/delom; // restoring original definition
			   eps = 1.-delta;
			   if (ll == 2){
			   		cFp = Xp[ind][ci-1]*eps + Xp[ind][ci]*delta;
			   		cFc = Xc[ind][ci-1]*eps + Xc[ind][ci]*delta;
			   }else{
				   cFp = XXp[ind][ci-1]*eps + XXp[ind][ci]*delta;
				   cFc = XXc[ind][ci-1]*eps + XXc[ind][ci]*delta;
			   }

			   // Apc = hamp*Apcos[mm+2];
			   // Aps = hamp*Apsin[mm+2]; // should be "-" for strain
			   // Acc = 2.*hamp*Accos[mm+2];
			   // Acs = 2.*hamp*Acsin[mm+2];

			   sinph = sin(faza - Om*xi + LISAWP_PI*0.25);
			   cosph = cos(faza - Om*xi + LISAWP_PI*0.25);

			   fact = 0.5*amp*sqrt(LISAWP_TWOPI/dom);
			   //hplus = Apc - img*Aps;
			   //hcross = Acc - img*Acs;
			  // Xf[ii] = Xf[ii] + fact*(cFp*hplus + cFc*hcross)*(cosph + img*sinph);
			   // ki = Nps - ii;
			   //xp = fact*hplus*(cosph + img*sinph);
			   //xc = fact*hcross*(cosph + img*sinph);

			   // xp = fact*(Apc*cosph + Aps*sinph + img*(Apc*sinph - Aps*cosph));
			   // xc = fact*(Acc*cosph + Acs*sinph + img*(Acc*sinph - Acs*cosph));

			   xp = fact*hamp*Ap[mm+2]*(cosph + img*sinph);
			   xc = 2*fact*hamp*Ac[mm+2]*(cosph + img*sinph);

			   // TODO output was changed
			   // -------------------------------------------------------------
			   //Xf(ii) += (cFp*xp + cFc*xc).real();
			   //Xf(ki) += (cFp*xp + cFc*xc).imag();
			   // Ampl = hamp*fact*( cFp*Ap[mm+2] + 2.0*cFc*Ac[mm+2] )


			   Xf_r[ii] += (cFp*xp + cFc*xc).real();
			   Xf_im[ii] += (cFp*xp + cFc*xc).imag();
			   if (ll == 2){
				   cFp = Yp[ind][ci-1]*eps + Yp[ind][ci]*delta;
				   cFc = Yc[ind][ci-1]*eps + Yc[ind][ci]*delta;
			   }else{
				   cFp = YYp[ind][ci-1]*eps + YYp[ind][ci]*delta;
				   cFc = YYc[ind][ci-1]*eps + YYc[ind][ci]*delta;
			   }

			   Yf_r[ii] += (cFp*xp + cFc*xc).real();
			   Yf_im[ii] += (cFp*xp + cFc*xc).imag();
			   //Yf(ii) += (cFp*xp + cFc*xc).real();
			   //Yf(ki) += (cFp*xp + cFc*xc).imag();
			   if (ll == 2){
				   cFp = Zp[ind][ci-1]*eps + Zp[ind][ci]*delta;
				   cFc = Zc[ind][ci-1]*eps + Zc[ind][ci]*delta;
			   }else{
				   cFp = ZZp[ind][ci-1]*eps + ZZp[ind][ci]*delta;
				   cFc = ZZc[ind][ci-1]*eps + ZZc[ind][ci]*delta;
			   }

			   Zf_r[ii] += (cFp*xp + cFc*xc).real();
			   Zf_im[ii] += (cFp*xp + cFc*xc).imag();


		   }// end Nyquist if

	   }// end of for-loop
   }// end of if

}



void FreqAK_RA::GetHarmonicsXYZ(EMRItemplate& S,  double** X_r, double** X_i, double** Y_r, double** Y_i, double** Z_r,
            double** Z_i, double* phi, double* gam, double* alph, double* tim, double** frqs, double** phase){


     // Currently I do not interpolate between the points for which the orbital evolution was computed

     double T = 0.0;
     double ec, ec2, ec3, ec4, ec5, ec6, ec7, hamp;
     int ind_low, ind_up;
     double delta, eps, amp;
     double  xi, Apc, Aps, Acc, Acs;
     int harmms[]={-2,-1,0,1,2};
     double fact;
     std::complex<double> hplus, hcross;
     std::complex<double> cFp;
     std::complex<double> cFc;
     std::complex<double> x, xp, xc;
     std::complex<double> img(0.0, 1.0);
     double om;
     double dOm = df*LISAWP_TWOPI;
     double Om, dom, delom;
     double xi_in, xi_fin, dxi;
     double orbOm = LISAWP_TWOPI/LISAWP_YRSID_SI;
     double AUsec =  LISAWP_AU_SI/LISAWP_C_SI;
     double DM = AUsec*S.stS;
     double faza, sinph, cosph;
     int ind;
     double nn,mm;
     std::string spr = "    ";

     double cS= S.ctS;
     double sS= S.stS;
     double cK= S.ctK;
     double sK= S.stK;
     double cSp= S.cpS;
     double sSp= S.spS;
     double cKp= S.cpK;
     double sKp= S.spK;
     double clam=cos(S.lam);
     double slam=sin(S.lam);
     double Sn = cS*cK + sS*sK*cos(S.phS  - S. phK);
     double cX = Sn;
     std::complex<double> Ampl;


     double sX = sqrt( sS*sS*cK*cK - 2.*sS*sSp*cK*cS*sK*sKp + \
                      cS*cS*sK*sK - 2.*cS*sK*cKp*sS*cSp*cK + \
                      sS*sS*cSp*cSp*sK*sK*sKp*sKp - 2.*sS*sS*cSp*sK*sK*sKp*sSp*cKp +\
                      sS*sS*sSp*sSp*sK*sK*cKp*cKp);

     double Apc1, Aps1, Apcn1, Apc2, Aps2, Apcn2, Aqc1, Aqs1, Aqcn1, Aqc2, Aqs2, Aqcn2;

       Apc1 = ( -cK*cKp*sS*cSp - cK*sKp*sS*sSp + sK*cS )/(sX);
       Aps1 = ( sKp*sS*cSp - cKp*sS*sSp )/(sX);
       Apcn1 = ( cK*cS + sK*cKp*sS*cSp + sK*sKp*sS*sSp - cX)*clam/(sX*slam);

       Apc2 = (sS*cSp*sKp - sS*sSp*cKp )*clam/(sX);
       Aps2 = ( cK*cKp*sS*cSp + cK*sKp*sS*sSp - cS*sK )*clam/(sX);
       Apcn2 = 0.0;


       Aqc1 = ( sS*cSp*sKp - sS*sSp*cKp  )*cX/(sX);
       Aqs1 = ( cK*cKp*sS*cSp + cK*sKp*sS*sSp - cS*sK )*cX/(sX);
       Aqcn1 = 0.0;


       Aqc2 = cX*clam*( cK*cKp*sS*cSp + cK*sKp*sS*sSp - sK*cS)/(sX);
       Aqs2 = -cX*clam*sS*( sKp*cSp - cKp*sSp )/(sX);
       Aqcn2 = -( cX*clam*clam*( cK*cS + sK*cKp*sS*cSp + sK*sKp*sS*sSp ) + \
                                1.- cX*cX - clam*clam )/(sX*slam);


     double Bp1c1 = 2.0*(Apc1*Apcn1 - Aqc1*Aqcn1 + Aqc2*Aqcn2 - Apc2*Apcn2);
     double Bp1c2 =  0.5*(Aps2*Aps2 - Aqc1*Aqc1  + Apc1*Apc1  - Aps1*Aps1 + \
                                Aqc2*Aqc2 + Aqs1*Aqs1 - Apc2*Apc2 - Aqs2*Aqs2);
     double Bp1s1 = 2.0*(Aqs2*Aqcn2 - Aps2*Apcn2 - Aqs1*Aqcn1 + Aps1*Apcn1);
     double Bp1s2 = (Apc1*Aps1 + Aqc2*Aqs2 - Apc2*Aps2 - Aqc1*Aqs1);
     double Bp1cn = 0.5*(Apc1*Apc1 + Aps1*Aps1 - Aqc1*Aqc1 - Aqs1*Aqs1 - Apc2*Apc2 \
                                + Aqc2*Aqc2 + Aqs2*Aqs2 - Aps2*Aps2) + Aqcn2*Aqcn2 - Aqcn1*Aqcn1 \
                                + Apcn1*Apcn1 - Apcn2*Apcn2;

     double Bp2c1 = (Apcn1*Apc2 + Apc1*Apcn2 - Aqcn1*Aqc2 - Aqc1*Aqcn2);
     double Bp2c2 = 0.5*(Aqs1*Aqs2 - Aps1*Aps2 + Apc1*Apc2 - Aqc1*Aqc2);
     double Bp2s1 = (Aps1*Apcn2 + Apcn1*Aps2 - Aqcn1*Aqs2 - Aqs1*Aqcn2);
     double Bp2s2 = 0.5*( Apc1*Aps2 - Aqc1*Aqs2 + Aps1*Apc2 - Aqs1*Aqc2);
     double Bp2cn = 0.5*(Aps1*Aps2 - Aqs1*Aqs2 - Aqc1*Aqc2 + Apc1*Apc2) -Aqcn1*Aqcn2 + Apcn1*Apcn2;

     double Bc1c1 = (-Apc2*Aqcn2 - Apcn2*Aqc2 + Apc1*Aqcn1 + Apcn1*Aqc1);
     double Bc1c2 = 0.5*( Apc1*Aqc1 - Aps1*Aqs1 - Apc2*Aqc2 + Aps2*Aqs2);
     double Bc1s1 = (Apcn1*Aqs1 - Aps2*Aqcn2 + Aps1*Aqcn1 - Apcn2*Aqs2);
     double Bc1s2 = 0.5*(-Apc2*Aqs2 + Apc1*Aqs1 - Aps2*Aqc2 + Aps1*Aqc1);
     double Bc1cn = -Apcn2*Aqcn2 + Apcn1*Aqcn1 + 0.5*(Apc1*Aqc1 - Aps2*Aqs2 + Aps1*Aqs1 - Apc2*Aqc2);

     double Bc2c1 = (Aqc1*Apcn2 + Aqcn1*Apc2 + Apc1*Aqcn2 + Apcn1*Aqc2);
     double Bc2c2 = 0.5*( Apc1*Aqc2 - Aps1*Aqs2 + Aqc1*Apc2 - Aqs1*Aps2);
     double Bc2s1 = (Apcn1*Aqs2 + Aqs1*Apcn2 + Aps1*Aqcn2 + Aqcn1*Aps2);
     double Bc2s2 = 0.5*(Aqc1*Aps2 + Apc1*Aqs2 + Aqs1*Apc2 + Aps1*Aqc2);
     double Bc2cn = Aqcn1*Apcn2 + Apcn1*Aqcn2 + 0.5*(Apc1*Aqc2 + Aqs1*Aps2 +Aps1*Aqs2 + Aqc1*Apc2);

     double AApcos[5],AApsin[5],AAccos[5],AAcsin[5];
       AApcos[0]=0.5*(Bp1c2+Bp2s2);
       AApsin[0]=0.5*(Bp2c2-Bp1s2);
       AAccos[0]=0.5*(Bc1c2+Bc2s2);
       AAcsin[0]=0.5*(Bc2c2-Bc1s2);
       AApcos[1]=0.5*(Bp1c1+Bp2s1);
       AApsin[1]=0.5*(Bp2c1-Bp1s1);
       AAccos[1]=0.5*(Bc1c1+Bc2s1);
       AAcsin[1]=0.5*(Bc2c1-Bc1s1);
       AApcos[2]=Bp1cn;
       AApsin[2]=Bp2cn;
       AAccos[2]=Bc1cn;
       AAcsin[2]=Bc2cn;
       AApcos[3]=0.5*(Bp1c1-Bp2s1);
       AApsin[3]=0.5*(Bp2c1+Bp1s1);
       AAccos[3]=0.5*(Bc1c1-Bc2s1);
       AAcsin[3]=0.5*(Bc2c1+Bc1s1);
       AApcos[4]=0.5*(Bp1c2-Bp2s2);
       AApsin[4]=0.5*(Bp2c2+Bp1s2);
       AAccos[4]=0.5*(Bc1c2-Bc2s2);
       AAcsin[4]=0.5*(Bc2c2+Bc1s2);
       int ki;

     // std::complex<double> test;
      //std::ofstream fout09("Data/Scratch.dat");
     for (int i=0; i<Nph; i++){
        T = tm_ph[i];
        tim[i] = T;
        xi = tm_ph[i] - DM*cos(orbOm*tm_ph[i]-S.phS);
        ind = 0;
        phi[i] = phiev[i];
        gam[i] = gamev[i];
        alph[i] = alpev[i];
        // loops over harmonics
        for (int j=0;j<5;j++) {
        	  nn=(double)(j+1);
        	  for(int jj=0; jj<5; jj++){
        	     mm = (double)jj-2.;
                     dom = nn*phiddotev[i] + 2.*gamddotev[i]  + mm*alpddotev[i];
        	     om = nn*LISAWP_TWOPI*nuev[i] + 2.*gamdotev[i]  + mm*alpdotev[i];
                     faza = (nn*phiev[i] + 2.*gamev[i]  + mm*alpev[i]);
                     ec = eccev[i];
                     ec2=ec*ec;
                     ec3=ec2*ec;
                     ec4=ec3*ec;
                     ec5=ec4*ec;
                     ec6=ec5*ec;
                     ec7=ec6*ec;
                     amp = Ampev[i];
                     switch (j+1) {
               			   case 1:
               				   hamp =3.*ec-1.625*ec3;
               				   break;
               			   case 2:
               				   hamp =-4.+10.*ec2-5.75*ec4;
               				   break;
               			   case 3:
               				   hamp =-9.*ec+21.375*ec3-15.046875*ec5;
               				   break;
               			   case 4:
               				   hamp =-16.*ec2+40.*ec4-101.*ec6/3.;
               				   break;
               			   case 5:
               				   hamp =-625.*ec3/24.+26875.*ec5/384.-210625.*ec7/3072.;
               				   break;
               	     }
                     cFp = Xp[ind][i];
                     cFc = Xc[ind][i];

                     Apc = hamp*AApcos[harmms[jj]+2];
         	     Aps = hamp*AApsin[harmms[jj]+2]; // should be "-" for strain
         	     Acc = 2.*hamp*AAccos[harmms[jj]+2];
         	     Acs = 2.*hamp*AAcsin[harmms[jj]+2];

                     phase[ind][i] =  faza - om*xi + LISAWP_PI*0.25; // NOTE that this is SPA phase
                     frqs[ind][i] = om;

                     fact = 0.5*amp*sqrt(LISAWP_TWOPI/dom);

                     Ampl = fact*( cFp*(Apc - img*Aps) + cFc*(Acc - img*Acs) );
                     X_r[ind][i] = Ampl.real();
                     X_i[ind][i] = Ampl.imag();


                     cFp = Yp[ind][i];
                     cFc = Yc[ind][i];
                     Ampl = fact*( cFp*(Apc - img*Aps) + cFc*(Acc - img*Acs) );
                     Y_r[ind][i] = Ampl.real();
                     Y_i[ind][i] = Ampl.imag();

                     cFp = Zp[ind][i];
                     cFc = Zc[ind][i];
                     Ampl = fact*( cFp*(Apc - img*Aps) + cFc*(Acc - img*Acs) );
                     Z_r[ind][i] = Ampl.real();
                     Z_i[ind][i] = Ampl.imag();

                     ind += 1;


                  }// jj loop
              } //j loop
     } // loop over time

}


double FreqAK_RA::Sinc(double y){
	double z;
	if(y != 0.0){
		z = sin(y)/y;
	}else{
		z = 1.0;
	}
	return(z);
}
