//
// 
// This code is generated from CPPINTS, a C++ program to generate the 
// analytical integrals based on Gaussian form primitive functions. 
// Copyright (C) 2015 The State University of New York at Buffalo 
// This software uses the MIT license as below: 
// 
// Permission is hereby granted, free of charge, to any person obtaining 
// a copy of this software and associated documentation files (the "Software"), 
// to deal in the Software without restriction, including without limitation 
// the rights to use, copy, modify, merge, publish, distribute, sublicense, 
// and/or sell copies of the Software, and to permit persons to whom the Software 
// is furnished to do so, subject to the following conditions: 
// 
// The above copyright notice and this permission notice shall be included in all 
// copies or substantial portions of the Software. 
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
// PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
// FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER  
// DEALINGS IN THE SOFTWARE. 
// 
// 

#include "constants.h"
#include <cstddef>
#include <math.h>

typedef int             Int;
typedef size_t          UInt;
#define THRESHOLD_MATH  0.00000000000001
#ifdef WITH_SINGLE_PRECISION
typedef float           Double;
#else
typedef double          Double;
#endif

//
//  here below is a list of variables used in the program
//
//  alpha is the bra1's exponent
//  beta  is the bra2's exponent
//  gamma is the ket1's exponent
//  delta is the ket2's exponent
//  A is the nuclear center for bra1
//  B is the nuclear center for bra2
//  C is the nuclear center for ket1
//  D is the nuclear center for ket2
//  P is the new center after bra1 combined with bra2
//  Q is the new center after ket1 combined with ket2
//  W is the new center after P combined with Q
//  pMax is maximum value of corresponding density matrix block(or value pair), used for ERI
//  omega is the exponential factor used for operator in form of erf(omega*r12)/r12
//  also omega could be the exponential factor used for operator in form of e^(-omega*r12^2)
//
//  variables:
//
//  zeta      = alpha + beta
//  eta       = gamma + delta
//  oned2z    = 1/(2*zeta)
//  oned2e    = 1/(2*eta)
//  onedz     = 1/zeta
//  onede     = 1/eta
//  kappa     = zeta + eta
//  onedk     = 1/kappa
//  oned2zeta = 1/(2*(alpha+beta+gamma))
//  xi        = alpha*beta*onedz
//  twoxi     = 2*alpha*beta*onedz
//  rho       = zeta*eta*onedk
//  rhod2zsq  = rho/(2*zeta*zeta)
//  rhod2esq  = rho/(2*eta*eta)
//  odorho    = omega/(rho+omega)
//  rhodorho  = rho/(rho+omega)
//  orhod2z2  = (rho/(2*zeta*zeta))*(omega/(rho+omega))
//  orhod2e2  = (rho/(2*eta*eta))*(omega/(rho+omega))
//  od2k      = (1/(2*kappa))*(omega/(rho+omega))
//  adz       = alpha*onedz
//  bdz       = beta*onedz
//  gde       = gamma*onede
//  gde       = delta*onede
//
//  input parameters based on primitive functions pair:
//
//  bra side shell pair is index as i
//  inp2  is the number of primitive pairs
//  iexp  is the array of 1/(alpha+beta)
//  icoe  is the array of ic_bra1*ic_bra2
//  ifac  is the array of pre-factor on bra side 
//  for (SS|SS)^{m} etc. type of integrals
//  ket side shell pair is index as j
//  jnp2  is the number of primitive pairs
//  jexp  is the array of 1/(gamma+delta)
//  jcoe  is the array of jc_ket1*jc_ket2
//  jfac  is the array of pre-factor on ket side 
//  for (SS|SS)^{m} etc. type of integrals
//

void hgp_os_twobodyoverlap_f_d(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_TWOBODYOVERLAP_H5x_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H4xy_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H4xz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H3x2y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H3xyz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H3x2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x3y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x2yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H2xy2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx4y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx3yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx2y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Hxy3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H5y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H4yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H3y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H2y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Hy4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_H5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G4x_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xy_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G2xyz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx2yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Gxy2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G4y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G3yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G2y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Gy3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_G4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_F3x_S = 0.0E0;
  Double I_TWOBODYOVERLAP_F2xy_S = 0.0E0;
  Double I_TWOBODYOVERLAP_F2xz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Fx2y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Fxyz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Fx2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_F3y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_F2yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Fy2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_F3z_S = 0.0E0;

  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double ic2   = icoe[ip2];
    Double fbra  = ifac[ip2];
    Double onedz = iexp[ip2];
    Double oned2z= 0.5E0*onedz;
    UInt offsetP = 3*ip2;
    Double PX    = P[offsetP  ];
    Double PY    = P[offsetP+1];
    Double PZ    = P[offsetP+2];
    Double PAX   = PX - A[0];
    Double PAY   = PY - A[1];
    Double PAZ   = PZ - A[2];
    Double I_TWOBODYOVERLAP_S_S_vrr = ic2*fbra;
    if(fabs(I_TWOBODYOVERLAP_S_S_vrr)<THRESHOLD_MATH) continue;


    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_P_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_Px_S_vrr = PAX*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Py_S_vrr = PAY*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Pz_S_vrr = PAZ*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 2 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_S_vrr = PAX*I_TWOBODYOVERLAP_Px_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_S_vrr = PAY*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_S_vrr = PAY*I_TWOBODYOVERLAP_Py_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_S_vrr = PAZ*I_TWOBODYOVERLAP_Pz_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_S_vrr = PAX*I_TWOBODYOVERLAP_D2x_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_S_vrr = PAY*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_S_vrr = PAZ*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_S_vrr = PAX*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fxyz_S_vrr = PAZ*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_S_vrr = PAX*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_S_vrr = PAY*I_TWOBODYOVERLAP_D2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fy2z_S_vrr = PAY*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_S_vrr = PAZ*I_TWOBODYOVERLAP_D2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_S_vrr = PAX*I_TWOBODYOVERLAP_F3x_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_S_vrr = PAY*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_S_vrr = PAZ*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_S_vrr = PAY*I_TWOBODYOVERLAP_F2xy_S_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_S_vrr = PAZ*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_S_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_S_vrr = PAX*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_S_vrr = PAY*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_S_vrr = PAX*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_S_vrr = PAY*I_TWOBODYOVERLAP_F3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_F2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_S_vrr = PAY*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_S_vrr = PAZ*I_TWOBODYOVERLAP_F3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_H5x_S_vrr = PAX*I_TWOBODYOVERLAP_G4x_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xy_S_vrr = PAY*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xz_S_vrr = PAZ*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_S_vrr = PAY*I_TWOBODYOVERLAP_G3xy_S_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_S_vrr = PAZ*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_G3xz_S_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_S_vrr = PAX*I_TWOBODYOVERLAP_Gx3y_S_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_S_vrr = PAY*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_S_vrr = PAX*I_TWOBODYOVERLAP_Gx3z_S_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_S_vrr = PAX*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_S_vrr = PAX*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_S_vrr = PAY*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_S_vrr = PAX*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5y_S_vrr = PAY*I_TWOBODYOVERLAP_G4y_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_H4yz_S_vrr = PAZ*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_G3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_S_vrr = PAY*I_TWOBODYOVERLAP_Gy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_S_vrr = PAY*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5z_S_vrr = PAZ*I_TWOBODYOVERLAP_G4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_S
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_TWOBODYOVERLAP_H5x_S += I_TWOBODYOVERLAP_H5x_S_vrr;
    I_TWOBODYOVERLAP_H4xy_S += I_TWOBODYOVERLAP_H4xy_S_vrr;
    I_TWOBODYOVERLAP_H4xz_S += I_TWOBODYOVERLAP_H4xz_S_vrr;
    I_TWOBODYOVERLAP_H3x2y_S += I_TWOBODYOVERLAP_H3x2y_S_vrr;
    I_TWOBODYOVERLAP_H3xyz_S += I_TWOBODYOVERLAP_H3xyz_S_vrr;
    I_TWOBODYOVERLAP_H3x2z_S += I_TWOBODYOVERLAP_H3x2z_S_vrr;
    I_TWOBODYOVERLAP_H2x3y_S += I_TWOBODYOVERLAP_H2x3y_S_vrr;
    I_TWOBODYOVERLAP_H2x2yz_S += I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    I_TWOBODYOVERLAP_H2xy2z_S += I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    I_TWOBODYOVERLAP_H2x3z_S += I_TWOBODYOVERLAP_H2x3z_S_vrr;
    I_TWOBODYOVERLAP_Hx4y_S += I_TWOBODYOVERLAP_Hx4y_S_vrr;
    I_TWOBODYOVERLAP_Hx3yz_S += I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    I_TWOBODYOVERLAP_Hx2y2z_S += I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    I_TWOBODYOVERLAP_Hxy3z_S += I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    I_TWOBODYOVERLAP_Hx4z_S += I_TWOBODYOVERLAP_Hx4z_S_vrr;
    I_TWOBODYOVERLAP_H5y_S += I_TWOBODYOVERLAP_H5y_S_vrr;
    I_TWOBODYOVERLAP_H4yz_S += I_TWOBODYOVERLAP_H4yz_S_vrr;
    I_TWOBODYOVERLAP_H3y2z_S += I_TWOBODYOVERLAP_H3y2z_S_vrr;
    I_TWOBODYOVERLAP_H2y3z_S += I_TWOBODYOVERLAP_H2y3z_S_vrr;
    I_TWOBODYOVERLAP_Hy4z_S += I_TWOBODYOVERLAP_Hy4z_S_vrr;
    I_TWOBODYOVERLAP_H5z_S += I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_TWOBODYOVERLAP_G4x_S += I_TWOBODYOVERLAP_G4x_S_vrr;
    I_TWOBODYOVERLAP_G3xy_S += I_TWOBODYOVERLAP_G3xy_S_vrr;
    I_TWOBODYOVERLAP_G3xz_S += I_TWOBODYOVERLAP_G3xz_S_vrr;
    I_TWOBODYOVERLAP_G2x2y_S += I_TWOBODYOVERLAP_G2x2y_S_vrr;
    I_TWOBODYOVERLAP_G2xyz_S += I_TWOBODYOVERLAP_G2xyz_S_vrr;
    I_TWOBODYOVERLAP_G2x2z_S += I_TWOBODYOVERLAP_G2x2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3y_S += I_TWOBODYOVERLAP_Gx3y_S_vrr;
    I_TWOBODYOVERLAP_Gx2yz_S += I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    I_TWOBODYOVERLAP_Gxy2z_S += I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3z_S += I_TWOBODYOVERLAP_Gx3z_S_vrr;
    I_TWOBODYOVERLAP_G4y_S += I_TWOBODYOVERLAP_G4y_S_vrr;
    I_TWOBODYOVERLAP_G3yz_S += I_TWOBODYOVERLAP_G3yz_S_vrr;
    I_TWOBODYOVERLAP_G2y2z_S += I_TWOBODYOVERLAP_G2y2z_S_vrr;
    I_TWOBODYOVERLAP_Gy3z_S += I_TWOBODYOVERLAP_Gy3z_S_vrr;
    I_TWOBODYOVERLAP_G4z_S += I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_TWOBODYOVERLAP_F3x_S += I_TWOBODYOVERLAP_F3x_S_vrr;
    I_TWOBODYOVERLAP_F2xy_S += I_TWOBODYOVERLAP_F2xy_S_vrr;
    I_TWOBODYOVERLAP_F2xz_S += I_TWOBODYOVERLAP_F2xz_S_vrr;
    I_TWOBODYOVERLAP_Fx2y_S += I_TWOBODYOVERLAP_Fx2y_S_vrr;
    I_TWOBODYOVERLAP_Fxyz_S += I_TWOBODYOVERLAP_Fxyz_S_vrr;
    I_TWOBODYOVERLAP_Fx2z_S += I_TWOBODYOVERLAP_Fx2z_S_vrr;
    I_TWOBODYOVERLAP_F3y_S += I_TWOBODYOVERLAP_F3y_S_vrr;
    I_TWOBODYOVERLAP_F2yz_S += I_TWOBODYOVERLAP_F2yz_S_vrr;
    I_TWOBODYOVERLAP_Fy2z_S += I_TWOBODYOVERLAP_Fy2z_S_vrr;
    I_TWOBODYOVERLAP_F3z_S += I_TWOBODYOVERLAP_F3z_S_vrr;
  }

  /************************************************************
   * declare the HRR1 result shell quartets in array form
   ************************************************************/

  /************************************************************
   * initilize the HRR steps : build the AB/CD variables
   ************************************************************/
  Double ABX = A[0] - B[0];
  Double ABY = A[1] - B[1];
  Double ABZ = A[2] - B[2];

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
   ************************************************************/
  Double I_TWOBODYOVERLAP_F3x_Px = I_TWOBODYOVERLAP_G4x_S+ABX*I_TWOBODYOVERLAP_F3x_S;
  Double I_TWOBODYOVERLAP_F2xy_Px = I_TWOBODYOVERLAP_G3xy_S+ABX*I_TWOBODYOVERLAP_F2xy_S;
  Double I_TWOBODYOVERLAP_F2xz_Px = I_TWOBODYOVERLAP_G3xz_S+ABX*I_TWOBODYOVERLAP_F2xz_S;
  Double I_TWOBODYOVERLAP_Fx2y_Px = I_TWOBODYOVERLAP_G2x2y_S+ABX*I_TWOBODYOVERLAP_Fx2y_S;
  Double I_TWOBODYOVERLAP_Fxyz_Px = I_TWOBODYOVERLAP_G2xyz_S+ABX*I_TWOBODYOVERLAP_Fxyz_S;
  Double I_TWOBODYOVERLAP_Fx2z_Px = I_TWOBODYOVERLAP_G2x2z_S+ABX*I_TWOBODYOVERLAP_Fx2z_S;
  Double I_TWOBODYOVERLAP_F3y_Px = I_TWOBODYOVERLAP_Gx3y_S+ABX*I_TWOBODYOVERLAP_F3y_S;
  Double I_TWOBODYOVERLAP_F2yz_Px = I_TWOBODYOVERLAP_Gx2yz_S+ABX*I_TWOBODYOVERLAP_F2yz_S;
  Double I_TWOBODYOVERLAP_Fy2z_Px = I_TWOBODYOVERLAP_Gxy2z_S+ABX*I_TWOBODYOVERLAP_Fy2z_S;
  Double I_TWOBODYOVERLAP_F3z_Px = I_TWOBODYOVERLAP_Gx3z_S+ABX*I_TWOBODYOVERLAP_F3z_S;
  Double I_TWOBODYOVERLAP_F3x_Py = I_TWOBODYOVERLAP_G3xy_S+ABY*I_TWOBODYOVERLAP_F3x_S;
  Double I_TWOBODYOVERLAP_F2xy_Py = I_TWOBODYOVERLAP_G2x2y_S+ABY*I_TWOBODYOVERLAP_F2xy_S;
  Double I_TWOBODYOVERLAP_F2xz_Py = I_TWOBODYOVERLAP_G2xyz_S+ABY*I_TWOBODYOVERLAP_F2xz_S;
  Double I_TWOBODYOVERLAP_Fx2y_Py = I_TWOBODYOVERLAP_Gx3y_S+ABY*I_TWOBODYOVERLAP_Fx2y_S;
  Double I_TWOBODYOVERLAP_Fxyz_Py = I_TWOBODYOVERLAP_Gx2yz_S+ABY*I_TWOBODYOVERLAP_Fxyz_S;
  Double I_TWOBODYOVERLAP_Fx2z_Py = I_TWOBODYOVERLAP_Gxy2z_S+ABY*I_TWOBODYOVERLAP_Fx2z_S;
  Double I_TWOBODYOVERLAP_F3y_Py = I_TWOBODYOVERLAP_G4y_S+ABY*I_TWOBODYOVERLAP_F3y_S;
  Double I_TWOBODYOVERLAP_F2yz_Py = I_TWOBODYOVERLAP_G3yz_S+ABY*I_TWOBODYOVERLAP_F2yz_S;
  Double I_TWOBODYOVERLAP_Fy2z_Py = I_TWOBODYOVERLAP_G2y2z_S+ABY*I_TWOBODYOVERLAP_Fy2z_S;
  Double I_TWOBODYOVERLAP_F3z_Py = I_TWOBODYOVERLAP_Gy3z_S+ABY*I_TWOBODYOVERLAP_F3z_S;
  Double I_TWOBODYOVERLAP_F3x_Pz = I_TWOBODYOVERLAP_G3xz_S+ABZ*I_TWOBODYOVERLAP_F3x_S;
  Double I_TWOBODYOVERLAP_F2xy_Pz = I_TWOBODYOVERLAP_G2xyz_S+ABZ*I_TWOBODYOVERLAP_F2xy_S;
  Double I_TWOBODYOVERLAP_F2xz_Pz = I_TWOBODYOVERLAP_G2x2z_S+ABZ*I_TWOBODYOVERLAP_F2xz_S;
  Double I_TWOBODYOVERLAP_Fx2y_Pz = I_TWOBODYOVERLAP_Gx2yz_S+ABZ*I_TWOBODYOVERLAP_Fx2y_S;
  Double I_TWOBODYOVERLAP_Fxyz_Pz = I_TWOBODYOVERLAP_Gxy2z_S+ABZ*I_TWOBODYOVERLAP_Fxyz_S;
  Double I_TWOBODYOVERLAP_Fx2z_Pz = I_TWOBODYOVERLAP_Gx3z_S+ABZ*I_TWOBODYOVERLAP_Fx2z_S;
  Double I_TWOBODYOVERLAP_F3y_Pz = I_TWOBODYOVERLAP_G3yz_S+ABZ*I_TWOBODYOVERLAP_F3y_S;
  Double I_TWOBODYOVERLAP_F2yz_Pz = I_TWOBODYOVERLAP_G2y2z_S+ABZ*I_TWOBODYOVERLAP_F2yz_S;
  Double I_TWOBODYOVERLAP_Fy2z_Pz = I_TWOBODYOVERLAP_Gy3z_S+ABZ*I_TWOBODYOVERLAP_Fy2z_S;
  Double I_TWOBODYOVERLAP_F3z_Pz = I_TWOBODYOVERLAP_G4z_S+ABZ*I_TWOBODYOVERLAP_F3z_S;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 6 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
   ************************************************************/
  Double I_TWOBODYOVERLAP_G4x_Px = I_TWOBODYOVERLAP_H5x_S+ABX*I_TWOBODYOVERLAP_G4x_S;
  Double I_TWOBODYOVERLAP_G3xy_Px = I_TWOBODYOVERLAP_H4xy_S+ABX*I_TWOBODYOVERLAP_G3xy_S;
  Double I_TWOBODYOVERLAP_G3xz_Px = I_TWOBODYOVERLAP_H4xz_S+ABX*I_TWOBODYOVERLAP_G3xz_S;
  Double I_TWOBODYOVERLAP_G2x2y_Px = I_TWOBODYOVERLAP_H3x2y_S+ABX*I_TWOBODYOVERLAP_G2x2y_S;
  Double I_TWOBODYOVERLAP_G2xyz_Px = I_TWOBODYOVERLAP_H3xyz_S+ABX*I_TWOBODYOVERLAP_G2xyz_S;
  Double I_TWOBODYOVERLAP_G2x2z_Px = I_TWOBODYOVERLAP_H3x2z_S+ABX*I_TWOBODYOVERLAP_G2x2z_S;
  Double I_TWOBODYOVERLAP_Gx3y_Px = I_TWOBODYOVERLAP_H2x3y_S+ABX*I_TWOBODYOVERLAP_Gx3y_S;
  Double I_TWOBODYOVERLAP_Gx2yz_Px = I_TWOBODYOVERLAP_H2x2yz_S+ABX*I_TWOBODYOVERLAP_Gx2yz_S;
  Double I_TWOBODYOVERLAP_Gxy2z_Px = I_TWOBODYOVERLAP_H2xy2z_S+ABX*I_TWOBODYOVERLAP_Gxy2z_S;
  Double I_TWOBODYOVERLAP_Gx3z_Px = I_TWOBODYOVERLAP_H2x3z_S+ABX*I_TWOBODYOVERLAP_Gx3z_S;
  Double I_TWOBODYOVERLAP_G4y_Px = I_TWOBODYOVERLAP_Hx4y_S+ABX*I_TWOBODYOVERLAP_G4y_S;
  Double I_TWOBODYOVERLAP_G3yz_Px = I_TWOBODYOVERLAP_Hx3yz_S+ABX*I_TWOBODYOVERLAP_G3yz_S;
  Double I_TWOBODYOVERLAP_G2y2z_Px = I_TWOBODYOVERLAP_Hx2y2z_S+ABX*I_TWOBODYOVERLAP_G2y2z_S;
  Double I_TWOBODYOVERLAP_Gy3z_Px = I_TWOBODYOVERLAP_Hxy3z_S+ABX*I_TWOBODYOVERLAP_Gy3z_S;
  Double I_TWOBODYOVERLAP_G4z_Px = I_TWOBODYOVERLAP_Hx4z_S+ABX*I_TWOBODYOVERLAP_G4z_S;
  Double I_TWOBODYOVERLAP_G3xy_Py = I_TWOBODYOVERLAP_H3x2y_S+ABY*I_TWOBODYOVERLAP_G3xy_S;
  Double I_TWOBODYOVERLAP_G3xz_Py = I_TWOBODYOVERLAP_H3xyz_S+ABY*I_TWOBODYOVERLAP_G3xz_S;
  Double I_TWOBODYOVERLAP_G2x2y_Py = I_TWOBODYOVERLAP_H2x3y_S+ABY*I_TWOBODYOVERLAP_G2x2y_S;
  Double I_TWOBODYOVERLAP_G2xyz_Py = I_TWOBODYOVERLAP_H2x2yz_S+ABY*I_TWOBODYOVERLAP_G2xyz_S;
  Double I_TWOBODYOVERLAP_G2x2z_Py = I_TWOBODYOVERLAP_H2xy2z_S+ABY*I_TWOBODYOVERLAP_G2x2z_S;
  Double I_TWOBODYOVERLAP_Gx3y_Py = I_TWOBODYOVERLAP_Hx4y_S+ABY*I_TWOBODYOVERLAP_Gx3y_S;
  Double I_TWOBODYOVERLAP_Gx2yz_Py = I_TWOBODYOVERLAP_Hx3yz_S+ABY*I_TWOBODYOVERLAP_Gx2yz_S;
  Double I_TWOBODYOVERLAP_Gxy2z_Py = I_TWOBODYOVERLAP_Hx2y2z_S+ABY*I_TWOBODYOVERLAP_Gxy2z_S;
  Double I_TWOBODYOVERLAP_Gx3z_Py = I_TWOBODYOVERLAP_Hxy3z_S+ABY*I_TWOBODYOVERLAP_Gx3z_S;
  Double I_TWOBODYOVERLAP_G4y_Py = I_TWOBODYOVERLAP_H5y_S+ABY*I_TWOBODYOVERLAP_G4y_S;
  Double I_TWOBODYOVERLAP_G3yz_Py = I_TWOBODYOVERLAP_H4yz_S+ABY*I_TWOBODYOVERLAP_G3yz_S;
  Double I_TWOBODYOVERLAP_G2y2z_Py = I_TWOBODYOVERLAP_H3y2z_S+ABY*I_TWOBODYOVERLAP_G2y2z_S;
  Double I_TWOBODYOVERLAP_Gy3z_Py = I_TWOBODYOVERLAP_H2y3z_S+ABY*I_TWOBODYOVERLAP_Gy3z_S;
  Double I_TWOBODYOVERLAP_G4z_Py = I_TWOBODYOVERLAP_Hy4z_S+ABY*I_TWOBODYOVERLAP_G4z_S;
  Double I_TWOBODYOVERLAP_G3xz_Pz = I_TWOBODYOVERLAP_H3x2z_S+ABZ*I_TWOBODYOVERLAP_G3xz_S;
  Double I_TWOBODYOVERLAP_G2xyz_Pz = I_TWOBODYOVERLAP_H2xy2z_S+ABZ*I_TWOBODYOVERLAP_G2xyz_S;
  Double I_TWOBODYOVERLAP_G2x2z_Pz = I_TWOBODYOVERLAP_H2x3z_S+ABZ*I_TWOBODYOVERLAP_G2x2z_S;
  Double I_TWOBODYOVERLAP_Gx2yz_Pz = I_TWOBODYOVERLAP_Hx2y2z_S+ABZ*I_TWOBODYOVERLAP_Gx2yz_S;
  Double I_TWOBODYOVERLAP_Gxy2z_Pz = I_TWOBODYOVERLAP_Hxy3z_S+ABZ*I_TWOBODYOVERLAP_Gxy2z_S;
  Double I_TWOBODYOVERLAP_Gx3z_Pz = I_TWOBODYOVERLAP_Hx4z_S+ABZ*I_TWOBODYOVERLAP_Gx3z_S;
  Double I_TWOBODYOVERLAP_G3yz_Pz = I_TWOBODYOVERLAP_H3y2z_S+ABZ*I_TWOBODYOVERLAP_G3yz_S;
  Double I_TWOBODYOVERLAP_G2y2z_Pz = I_TWOBODYOVERLAP_H2y3z_S+ABZ*I_TWOBODYOVERLAP_G2y2z_S;
  Double I_TWOBODYOVERLAP_Gy3z_Pz = I_TWOBODYOVERLAP_Hy4z_S+ABZ*I_TWOBODYOVERLAP_Gy3z_S;
  Double I_TWOBODYOVERLAP_G4z_Pz = I_TWOBODYOVERLAP_H5z_S+ABZ*I_TWOBODYOVERLAP_G4z_S;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_D
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
   ************************************************************/
  abcd[0] = I_TWOBODYOVERLAP_G4x_Px+ABX*I_TWOBODYOVERLAP_F3x_Px;
  abcd[1] = I_TWOBODYOVERLAP_G3xy_Px+ABX*I_TWOBODYOVERLAP_F2xy_Px;
  abcd[2] = I_TWOBODYOVERLAP_G3xz_Px+ABX*I_TWOBODYOVERLAP_F2xz_Px;
  abcd[3] = I_TWOBODYOVERLAP_G2x2y_Px+ABX*I_TWOBODYOVERLAP_Fx2y_Px;
  abcd[4] = I_TWOBODYOVERLAP_G2xyz_Px+ABX*I_TWOBODYOVERLAP_Fxyz_Px;
  abcd[5] = I_TWOBODYOVERLAP_G2x2z_Px+ABX*I_TWOBODYOVERLAP_Fx2z_Px;
  abcd[6] = I_TWOBODYOVERLAP_Gx3y_Px+ABX*I_TWOBODYOVERLAP_F3y_Px;
  abcd[7] = I_TWOBODYOVERLAP_Gx2yz_Px+ABX*I_TWOBODYOVERLAP_F2yz_Px;
  abcd[8] = I_TWOBODYOVERLAP_Gxy2z_Px+ABX*I_TWOBODYOVERLAP_Fy2z_Px;
  abcd[9] = I_TWOBODYOVERLAP_Gx3z_Px+ABX*I_TWOBODYOVERLAP_F3z_Px;
  abcd[10] = I_TWOBODYOVERLAP_G3xy_Px+ABY*I_TWOBODYOVERLAP_F3x_Px;
  abcd[11] = I_TWOBODYOVERLAP_G2x2y_Px+ABY*I_TWOBODYOVERLAP_F2xy_Px;
  abcd[12] = I_TWOBODYOVERLAP_G2xyz_Px+ABY*I_TWOBODYOVERLAP_F2xz_Px;
  abcd[13] = I_TWOBODYOVERLAP_Gx3y_Px+ABY*I_TWOBODYOVERLAP_Fx2y_Px;
  abcd[14] = I_TWOBODYOVERLAP_Gx2yz_Px+ABY*I_TWOBODYOVERLAP_Fxyz_Px;
  abcd[15] = I_TWOBODYOVERLAP_Gxy2z_Px+ABY*I_TWOBODYOVERLAP_Fx2z_Px;
  abcd[16] = I_TWOBODYOVERLAP_G4y_Px+ABY*I_TWOBODYOVERLAP_F3y_Px;
  abcd[17] = I_TWOBODYOVERLAP_G3yz_Px+ABY*I_TWOBODYOVERLAP_F2yz_Px;
  abcd[18] = I_TWOBODYOVERLAP_G2y2z_Px+ABY*I_TWOBODYOVERLAP_Fy2z_Px;
  abcd[19] = I_TWOBODYOVERLAP_Gy3z_Px+ABY*I_TWOBODYOVERLAP_F3z_Px;
  abcd[20] = I_TWOBODYOVERLAP_G3xz_Px+ABZ*I_TWOBODYOVERLAP_F3x_Px;
  abcd[21] = I_TWOBODYOVERLAP_G2xyz_Px+ABZ*I_TWOBODYOVERLAP_F2xy_Px;
  abcd[22] = I_TWOBODYOVERLAP_G2x2z_Px+ABZ*I_TWOBODYOVERLAP_F2xz_Px;
  abcd[23] = I_TWOBODYOVERLAP_Gx2yz_Px+ABZ*I_TWOBODYOVERLAP_Fx2y_Px;
  abcd[24] = I_TWOBODYOVERLAP_Gxy2z_Px+ABZ*I_TWOBODYOVERLAP_Fxyz_Px;
  abcd[25] = I_TWOBODYOVERLAP_Gx3z_Px+ABZ*I_TWOBODYOVERLAP_Fx2z_Px;
  abcd[26] = I_TWOBODYOVERLAP_G3yz_Px+ABZ*I_TWOBODYOVERLAP_F3y_Px;
  abcd[27] = I_TWOBODYOVERLAP_G2y2z_Px+ABZ*I_TWOBODYOVERLAP_F2yz_Px;
  abcd[28] = I_TWOBODYOVERLAP_Gy3z_Px+ABZ*I_TWOBODYOVERLAP_Fy2z_Px;
  abcd[29] = I_TWOBODYOVERLAP_G4z_Px+ABZ*I_TWOBODYOVERLAP_F3z_Px;
  abcd[30] = I_TWOBODYOVERLAP_G3xy_Py+ABY*I_TWOBODYOVERLAP_F3x_Py;
  abcd[31] = I_TWOBODYOVERLAP_G2x2y_Py+ABY*I_TWOBODYOVERLAP_F2xy_Py;
  abcd[32] = I_TWOBODYOVERLAP_G2xyz_Py+ABY*I_TWOBODYOVERLAP_F2xz_Py;
  abcd[33] = I_TWOBODYOVERLAP_Gx3y_Py+ABY*I_TWOBODYOVERLAP_Fx2y_Py;
  abcd[34] = I_TWOBODYOVERLAP_Gx2yz_Py+ABY*I_TWOBODYOVERLAP_Fxyz_Py;
  abcd[35] = I_TWOBODYOVERLAP_Gxy2z_Py+ABY*I_TWOBODYOVERLAP_Fx2z_Py;
  abcd[36] = I_TWOBODYOVERLAP_G4y_Py+ABY*I_TWOBODYOVERLAP_F3y_Py;
  abcd[37] = I_TWOBODYOVERLAP_G3yz_Py+ABY*I_TWOBODYOVERLAP_F2yz_Py;
  abcd[38] = I_TWOBODYOVERLAP_G2y2z_Py+ABY*I_TWOBODYOVERLAP_Fy2z_Py;
  abcd[39] = I_TWOBODYOVERLAP_Gy3z_Py+ABY*I_TWOBODYOVERLAP_F3z_Py;
  abcd[40] = I_TWOBODYOVERLAP_G3xz_Py+ABZ*I_TWOBODYOVERLAP_F3x_Py;
  abcd[41] = I_TWOBODYOVERLAP_G2xyz_Py+ABZ*I_TWOBODYOVERLAP_F2xy_Py;
  abcd[42] = I_TWOBODYOVERLAP_G2x2z_Py+ABZ*I_TWOBODYOVERLAP_F2xz_Py;
  abcd[43] = I_TWOBODYOVERLAP_Gx2yz_Py+ABZ*I_TWOBODYOVERLAP_Fx2y_Py;
  abcd[44] = I_TWOBODYOVERLAP_Gxy2z_Py+ABZ*I_TWOBODYOVERLAP_Fxyz_Py;
  abcd[45] = I_TWOBODYOVERLAP_Gx3z_Py+ABZ*I_TWOBODYOVERLAP_Fx2z_Py;
  abcd[46] = I_TWOBODYOVERLAP_G3yz_Py+ABZ*I_TWOBODYOVERLAP_F3y_Py;
  abcd[47] = I_TWOBODYOVERLAP_G2y2z_Py+ABZ*I_TWOBODYOVERLAP_F2yz_Py;
  abcd[48] = I_TWOBODYOVERLAP_Gy3z_Py+ABZ*I_TWOBODYOVERLAP_Fy2z_Py;
  abcd[49] = I_TWOBODYOVERLAP_G4z_Py+ABZ*I_TWOBODYOVERLAP_F3z_Py;
  abcd[50] = I_TWOBODYOVERLAP_G3xz_Pz+ABZ*I_TWOBODYOVERLAP_F3x_Pz;
  abcd[51] = I_TWOBODYOVERLAP_G2xyz_Pz+ABZ*I_TWOBODYOVERLAP_F2xy_Pz;
  abcd[52] = I_TWOBODYOVERLAP_G2x2z_Pz+ABZ*I_TWOBODYOVERLAP_F2xz_Pz;
  abcd[53] = I_TWOBODYOVERLAP_Gx2yz_Pz+ABZ*I_TWOBODYOVERLAP_Fx2y_Pz;
  abcd[54] = I_TWOBODYOVERLAP_Gxy2z_Pz+ABZ*I_TWOBODYOVERLAP_Fxyz_Pz;
  abcd[55] = I_TWOBODYOVERLAP_Gx3z_Pz+ABZ*I_TWOBODYOVERLAP_Fx2z_Pz;
  abcd[56] = I_TWOBODYOVERLAP_G3yz_Pz+ABZ*I_TWOBODYOVERLAP_F3y_Pz;
  abcd[57] = I_TWOBODYOVERLAP_G2y2z_Pz+ABZ*I_TWOBODYOVERLAP_F2yz_Pz;
  abcd[58] = I_TWOBODYOVERLAP_Gy3z_Pz+ABZ*I_TWOBODYOVERLAP_Fy2z_Pz;
  abcd[59] = I_TWOBODYOVERLAP_G4z_Pz+ABZ*I_TWOBODYOVERLAP_F3z_Pz;
}
