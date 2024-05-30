/*
Copyright (C) 2005  S. Babak

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/


/******************  CVS info ************************
#define CVSTAG "$Name:  $"
#define CVSHEADER "$Header: /afs/aeiw/cvsroot/waves/people/stas/LISAWP/base/src/EMRItemplate.cc,v 1.2 2008/04/08 11:20:54 stba Exp $"
*/

#include "EMRItemplate.h"

EMRItemplate::EMRItemplate() {}

EMRItemplate::EMRItemplate(double Mass, double mu, double spin, double lambda, double quadmom){

    M = Mass;
    m = mu;
    a = spin;
    lam = lambda;
    qm = quadmom;

    Mt = M*LISAWP_MTSUN_SI;
    mt = m*LISAWP_MTSUN_SI;

    SkySet = false;

}

void EMRItemplate::SetPosition(double thetaS, double phiS, double thetaK, double phiK, double D){

    thS = thetaS;
    phS = phiS;
    thK = thetaK;
    phK = phiK;

    ctS = cos(thetaS);
    stS = sin(thetaS);
    ctK = cos(thetaK);
    stK = sin(thetaK);

    cpS = cos(phiS);
    spS = sin(phiS);
    cpK = cos(phiK);
    spK = sin(phiK);
     Mt = M*LISAWP_MTSUN_SI;
     mt = m*LISAWP_MTSUN_SI;

    dist = D*(LISAWP_PC_SI/LISAWP_C_SI); // distance in seconds

    Ampl = mt/dist; // dimensionless amplitude mu/D

    SkySet = true;



}

EMRItemplate::EMRItemplate(const EMRItemplate& p){

    M = p.M;
    Mt = p.Mt;
    m = p.m;
    mt = p.mt;
    a = p.a;
    lam = p.lam;
    qm = p.qm;

    e0 = p.e0;
    nu0 = p.nu0;
    Phi0 = p.Phi0;
    gamma0 = p.gamma0;
    alpha0 = p.alpha0;
    t0 = p.t0;
    fgam0 = p.fgam0;
    falph0 = p.falph0;

    tPl = p.tPl;
    e_pl =  p.e_pl;
    nu_pl = p.nu_pl;
    Phi_pl = p.Phi_pl;
    alpha_pl = p.alpha_pl;
    gamma_pl = p.gamma_pl;

    SetPosition(p.thS, p.phS, p.thK, p.phK, p.dist);

    dist = p.dist;
    Ampl = p.Ampl;

    SNR = p.SNR;
    LogL = p.LogL;

//

}

EMRItemplate EMRItemplate::operator=(const EMRItemplate& p){

    this->M = p.M;
    this->Mt = p.Mt;
    this->m = p.m;
    this->mt = p.mt;
    this->a = p.a;
    this->lam = p.lam;
    this->qm = p.qm;

    this->e0 = p.e0;
    this->nu0 = p.nu0;
    this->Phi0 = p.Phi0;
    this->gamma0 = p.gamma0;
    this->alpha0 = p.alpha0;
    this->t0 = p.t0;
    this->fgam0 = p.fgam0;
    this->falph0 = p.falph0;

    this->tPl = p.tPl;
    this->e_pl =  p.e_pl;
    this->nu_pl = p.nu_pl;
    this->Phi_pl = p.Phi_pl;
    this->alpha_pl = p.alpha_pl;
    this->gamma_pl = p.gamma_pl;

    this->SetPosition(p.thS, p.phS, p.thK, p.phK, p.dist);

    this->dist = p.dist;
    this->Ampl = p.Ampl;


    this->SNR = p.SNR;
    this->LogL = p.LogL;


    return *this;

}


//EMRItemplate EMRItemplate::operator=(const EMRItemplate&){
//}
