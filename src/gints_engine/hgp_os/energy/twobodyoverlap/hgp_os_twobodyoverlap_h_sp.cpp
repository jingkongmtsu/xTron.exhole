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

void hgp_os_twobodyoverlap_h_sp(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_TWOBODYOVERLAP_I6x_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xy_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_I4xyz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x2yz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_I3xy2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x3yz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x2y2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_I2xy3z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix4yz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix3y2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix2y3z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Ixy4z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_I6y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_I5yz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_I4y2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_I3y3z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_I2y4z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Iy5z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_I6z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H5x_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H4xy_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H4xz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H3x2y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H3xyz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H3x2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x3y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x2yz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H2xy2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x3z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx4y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx3yz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx2y2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Hxy3z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx4z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H5y_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H4yz_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H3y2z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H2y3z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_Hy4z_S_C1005 = 0.0E0;
  Double I_TWOBODYOVERLAP_H5z_S_C1005 = 0.0E0;

  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double ic2   = icoe[ip2];
    Double ic2_1 = icoe[ip2+1*inp2];
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
    Double I_TWOBODYOVERLAP_S_S_vrr = fbra;
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
     * totally 3 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_S_vrr = PAX*I_TWOBODYOVERLAP_Px_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_S_vrr = PAY*I_TWOBODYOVERLAP_Py_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_D2z_S_vrr = PAZ*I_TWOBODYOVERLAP_Pz_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 4 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_S_vrr = PAX*I_TWOBODYOVERLAP_D2x_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_S_vrr = PAY*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_S_vrr = PAZ*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_S_vrr = PAY*I_TWOBODYOVERLAP_D2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_F3z_S_vrr = PAZ*I_TWOBODYOVERLAP_D2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 3 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_S_vrr = PAX*I_TWOBODYOVERLAP_F3x_S_vrr+3*oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_S_vrr = PAY*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_S_vrr = PAZ*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_S_vrr = PAY*I_TWOBODYOVERLAP_F2xy_S_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_F2xz_S_vrr+oned2z*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_S_vrr = PAX*I_TWOBODYOVERLAP_F3y_S_vrr;
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
     * shell quartet name: SQ_TWOBODYOVERLAP_I_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_I6x_S_vrr = PAX*I_TWOBODYOVERLAP_H5x_S_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xy_S_vrr = PAY*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xz_S_vrr = PAZ*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_S_vrr = PAY*I_TWOBODYOVERLAP_H4xy_S_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_S_vrr = PAZ*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_H4xz_S_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_S_vrr = PAY*I_TWOBODYOVERLAP_H3x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_S_vrr = PAY*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_S_vrr = PAZ*I_TWOBODYOVERLAP_H3x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_S_vrr = PAX*I_TWOBODYOVERLAP_Hx4y_S_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_H2x2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_S_vrr = PAY*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_S_vrr = PAX*I_TWOBODYOVERLAP_Hx4z_S_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_S_vrr = PAX*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_S_vrr = PAZ*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_S_vrr = PAX*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_S_vrr = PAX*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_S_vrr = PAY*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_S_vrr = PAX*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6y_S_vrr = PAY*I_TWOBODYOVERLAP_H5y_S_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_I5yz_S_vrr = PAZ*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_H4yz_S_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_S_vrr = PAZ*I_TWOBODYOVERLAP_H3y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_S_vrr = PAY*I_TWOBODYOVERLAP_Hy4z_S_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_S_vrr = PAY*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6z_S_vrr = PAZ*I_TWOBODYOVERLAP_H5z_S_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_S_C5
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_H_S_C5_coefs = ic2;
    abcd[0] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_H5x_S_vrr;
    abcd[1] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_H4xy_S_vrr;
    abcd[2] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_H4xz_S_vrr;
    abcd[3] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    abcd[4] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    abcd[5] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    abcd[6] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    abcd[7] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    abcd[8] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    abcd[9] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    abcd[10] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    abcd[11] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    abcd[12] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    abcd[13] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    abcd[14] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    abcd[15] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_H5y_S_vrr;
    abcd[16] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_H4yz_S_vrr;
    abcd[17] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    abcd[18] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    abcd[19] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    abcd[20] += SQ_TWOBODYOVERLAP_H_S_C5_coefs*I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_I_S_C1005
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_I_S_C1005_coefs = ic2_1;
    I_TWOBODYOVERLAP_I6x_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I6x_S_vrr;
    I_TWOBODYOVERLAP_I5xy_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I5xy_S_vrr;
    I_TWOBODYOVERLAP_I5xz_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I5xz_S_vrr;
    I_TWOBODYOVERLAP_I4x2y_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    I_TWOBODYOVERLAP_I4xyz_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I4xyz_S_vrr;
    I_TWOBODYOVERLAP_I4x2z_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3y_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    I_TWOBODYOVERLAP_I3x2yz_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    I_TWOBODYOVERLAP_I3xy2z_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I3xy2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3z_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4y_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    I_TWOBODYOVERLAP_I2x3yz_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    I_TWOBODYOVERLAP_I2x2y2z_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I2x2y2z_S_vrr;
    I_TWOBODYOVERLAP_I2xy3z_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4z_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5y_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    I_TWOBODYOVERLAP_Ix4yz_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_Ix4yz_S_vrr;
    I_TWOBODYOVERLAP_Ix3y2z_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_Ix3y2z_S_vrr;
    I_TWOBODYOVERLAP_Ix2y3z_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_Ix2y3z_S_vrr;
    I_TWOBODYOVERLAP_Ixy4z_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_Ixy4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5z_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    I_TWOBODYOVERLAP_I6y_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I6y_S_vrr;
    I_TWOBODYOVERLAP_I5yz_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I5yz_S_vrr;
    I_TWOBODYOVERLAP_I4y2z_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    I_TWOBODYOVERLAP_I3y3z_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    I_TWOBODYOVERLAP_I2y4z_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    I_TWOBODYOVERLAP_Iy5z_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    I_TWOBODYOVERLAP_I6z_S_C1005 += SQ_TWOBODYOVERLAP_I_S_C1005_coefs*I_TWOBODYOVERLAP_I6z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_S_C1005
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_H_S_C1005_coefs = ic2_1;
    I_TWOBODYOVERLAP_H5x_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H5x_S_vrr;
    I_TWOBODYOVERLAP_H4xy_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H4xy_S_vrr;
    I_TWOBODYOVERLAP_H4xz_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H4xz_S_vrr;
    I_TWOBODYOVERLAP_H3x2y_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    I_TWOBODYOVERLAP_H3xyz_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    I_TWOBODYOVERLAP_H3x2z_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    I_TWOBODYOVERLAP_H2x3y_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    I_TWOBODYOVERLAP_H2x2yz_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    I_TWOBODYOVERLAP_H2xy2z_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    I_TWOBODYOVERLAP_H2x3z_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    I_TWOBODYOVERLAP_Hx4y_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    I_TWOBODYOVERLAP_Hx3yz_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    I_TWOBODYOVERLAP_Hx2y2z_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    I_TWOBODYOVERLAP_Hxy3z_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    I_TWOBODYOVERLAP_Hx4z_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    I_TWOBODYOVERLAP_H5y_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H5y_S_vrr;
    I_TWOBODYOVERLAP_H4yz_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H4yz_S_vrr;
    I_TWOBODYOVERLAP_H3y2z_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    I_TWOBODYOVERLAP_H2y3z_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    I_TWOBODYOVERLAP_Hy4z_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    I_TWOBODYOVERLAP_H5z_S_C1005 += SQ_TWOBODYOVERLAP_H_S_C1005_coefs*I_TWOBODYOVERLAP_H5z_S_vrr;
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
   * shell quartet name: SQ_TWOBODYOVERLAP_H_P_C1005
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S_C1005
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_C1005
   ************************************************************/
  abcd[21] = I_TWOBODYOVERLAP_I6x_S_C1005+ABX*I_TWOBODYOVERLAP_H5x_S_C1005;
  abcd[22] = I_TWOBODYOVERLAP_I5xy_S_C1005+ABX*I_TWOBODYOVERLAP_H4xy_S_C1005;
  abcd[23] = I_TWOBODYOVERLAP_I5xz_S_C1005+ABX*I_TWOBODYOVERLAP_H4xz_S_C1005;
  abcd[24] = I_TWOBODYOVERLAP_I4x2y_S_C1005+ABX*I_TWOBODYOVERLAP_H3x2y_S_C1005;
  abcd[25] = I_TWOBODYOVERLAP_I4xyz_S_C1005+ABX*I_TWOBODYOVERLAP_H3xyz_S_C1005;
  abcd[26] = I_TWOBODYOVERLAP_I4x2z_S_C1005+ABX*I_TWOBODYOVERLAP_H3x2z_S_C1005;
  abcd[27] = I_TWOBODYOVERLAP_I3x3y_S_C1005+ABX*I_TWOBODYOVERLAP_H2x3y_S_C1005;
  abcd[28] = I_TWOBODYOVERLAP_I3x2yz_S_C1005+ABX*I_TWOBODYOVERLAP_H2x2yz_S_C1005;
  abcd[29] = I_TWOBODYOVERLAP_I3xy2z_S_C1005+ABX*I_TWOBODYOVERLAP_H2xy2z_S_C1005;
  abcd[30] = I_TWOBODYOVERLAP_I3x3z_S_C1005+ABX*I_TWOBODYOVERLAP_H2x3z_S_C1005;
  abcd[31] = I_TWOBODYOVERLAP_I2x4y_S_C1005+ABX*I_TWOBODYOVERLAP_Hx4y_S_C1005;
  abcd[32] = I_TWOBODYOVERLAP_I2x3yz_S_C1005+ABX*I_TWOBODYOVERLAP_Hx3yz_S_C1005;
  abcd[33] = I_TWOBODYOVERLAP_I2x2y2z_S_C1005+ABX*I_TWOBODYOVERLAP_Hx2y2z_S_C1005;
  abcd[34] = I_TWOBODYOVERLAP_I2xy3z_S_C1005+ABX*I_TWOBODYOVERLAP_Hxy3z_S_C1005;
  abcd[35] = I_TWOBODYOVERLAP_I2x4z_S_C1005+ABX*I_TWOBODYOVERLAP_Hx4z_S_C1005;
  abcd[36] = I_TWOBODYOVERLAP_Ix5y_S_C1005+ABX*I_TWOBODYOVERLAP_H5y_S_C1005;
  abcd[37] = I_TWOBODYOVERLAP_Ix4yz_S_C1005+ABX*I_TWOBODYOVERLAP_H4yz_S_C1005;
  abcd[38] = I_TWOBODYOVERLAP_Ix3y2z_S_C1005+ABX*I_TWOBODYOVERLAP_H3y2z_S_C1005;
  abcd[39] = I_TWOBODYOVERLAP_Ix2y3z_S_C1005+ABX*I_TWOBODYOVERLAP_H2y3z_S_C1005;
  abcd[40] = I_TWOBODYOVERLAP_Ixy4z_S_C1005+ABX*I_TWOBODYOVERLAP_Hy4z_S_C1005;
  abcd[41] = I_TWOBODYOVERLAP_Ix5z_S_C1005+ABX*I_TWOBODYOVERLAP_H5z_S_C1005;
  abcd[42] = I_TWOBODYOVERLAP_I5xy_S_C1005+ABY*I_TWOBODYOVERLAP_H5x_S_C1005;
  abcd[43] = I_TWOBODYOVERLAP_I4x2y_S_C1005+ABY*I_TWOBODYOVERLAP_H4xy_S_C1005;
  abcd[44] = I_TWOBODYOVERLAP_I4xyz_S_C1005+ABY*I_TWOBODYOVERLAP_H4xz_S_C1005;
  abcd[45] = I_TWOBODYOVERLAP_I3x3y_S_C1005+ABY*I_TWOBODYOVERLAP_H3x2y_S_C1005;
  abcd[46] = I_TWOBODYOVERLAP_I3x2yz_S_C1005+ABY*I_TWOBODYOVERLAP_H3xyz_S_C1005;
  abcd[47] = I_TWOBODYOVERLAP_I3xy2z_S_C1005+ABY*I_TWOBODYOVERLAP_H3x2z_S_C1005;
  abcd[48] = I_TWOBODYOVERLAP_I2x4y_S_C1005+ABY*I_TWOBODYOVERLAP_H2x3y_S_C1005;
  abcd[49] = I_TWOBODYOVERLAP_I2x3yz_S_C1005+ABY*I_TWOBODYOVERLAP_H2x2yz_S_C1005;
  abcd[50] = I_TWOBODYOVERLAP_I2x2y2z_S_C1005+ABY*I_TWOBODYOVERLAP_H2xy2z_S_C1005;
  abcd[51] = I_TWOBODYOVERLAP_I2xy3z_S_C1005+ABY*I_TWOBODYOVERLAP_H2x3z_S_C1005;
  abcd[52] = I_TWOBODYOVERLAP_Ix5y_S_C1005+ABY*I_TWOBODYOVERLAP_Hx4y_S_C1005;
  abcd[53] = I_TWOBODYOVERLAP_Ix4yz_S_C1005+ABY*I_TWOBODYOVERLAP_Hx3yz_S_C1005;
  abcd[54] = I_TWOBODYOVERLAP_Ix3y2z_S_C1005+ABY*I_TWOBODYOVERLAP_Hx2y2z_S_C1005;
  abcd[55] = I_TWOBODYOVERLAP_Ix2y3z_S_C1005+ABY*I_TWOBODYOVERLAP_Hxy3z_S_C1005;
  abcd[56] = I_TWOBODYOVERLAP_Ixy4z_S_C1005+ABY*I_TWOBODYOVERLAP_Hx4z_S_C1005;
  abcd[57] = I_TWOBODYOVERLAP_I6y_S_C1005+ABY*I_TWOBODYOVERLAP_H5y_S_C1005;
  abcd[58] = I_TWOBODYOVERLAP_I5yz_S_C1005+ABY*I_TWOBODYOVERLAP_H4yz_S_C1005;
  abcd[59] = I_TWOBODYOVERLAP_I4y2z_S_C1005+ABY*I_TWOBODYOVERLAP_H3y2z_S_C1005;
  abcd[60] = I_TWOBODYOVERLAP_I3y3z_S_C1005+ABY*I_TWOBODYOVERLAP_H2y3z_S_C1005;
  abcd[61] = I_TWOBODYOVERLAP_I2y4z_S_C1005+ABY*I_TWOBODYOVERLAP_Hy4z_S_C1005;
  abcd[62] = I_TWOBODYOVERLAP_Iy5z_S_C1005+ABY*I_TWOBODYOVERLAP_H5z_S_C1005;
  abcd[63] = I_TWOBODYOVERLAP_I5xz_S_C1005+ABZ*I_TWOBODYOVERLAP_H5x_S_C1005;
  abcd[64] = I_TWOBODYOVERLAP_I4xyz_S_C1005+ABZ*I_TWOBODYOVERLAP_H4xy_S_C1005;
  abcd[65] = I_TWOBODYOVERLAP_I4x2z_S_C1005+ABZ*I_TWOBODYOVERLAP_H4xz_S_C1005;
  abcd[66] = I_TWOBODYOVERLAP_I3x2yz_S_C1005+ABZ*I_TWOBODYOVERLAP_H3x2y_S_C1005;
  abcd[67] = I_TWOBODYOVERLAP_I3xy2z_S_C1005+ABZ*I_TWOBODYOVERLAP_H3xyz_S_C1005;
  abcd[68] = I_TWOBODYOVERLAP_I3x3z_S_C1005+ABZ*I_TWOBODYOVERLAP_H3x2z_S_C1005;
  abcd[69] = I_TWOBODYOVERLAP_I2x3yz_S_C1005+ABZ*I_TWOBODYOVERLAP_H2x3y_S_C1005;
  abcd[70] = I_TWOBODYOVERLAP_I2x2y2z_S_C1005+ABZ*I_TWOBODYOVERLAP_H2x2yz_S_C1005;
  abcd[71] = I_TWOBODYOVERLAP_I2xy3z_S_C1005+ABZ*I_TWOBODYOVERLAP_H2xy2z_S_C1005;
  abcd[72] = I_TWOBODYOVERLAP_I2x4z_S_C1005+ABZ*I_TWOBODYOVERLAP_H2x3z_S_C1005;
  abcd[73] = I_TWOBODYOVERLAP_Ix4yz_S_C1005+ABZ*I_TWOBODYOVERLAP_Hx4y_S_C1005;
  abcd[74] = I_TWOBODYOVERLAP_Ix3y2z_S_C1005+ABZ*I_TWOBODYOVERLAP_Hx3yz_S_C1005;
  abcd[75] = I_TWOBODYOVERLAP_Ix2y3z_S_C1005+ABZ*I_TWOBODYOVERLAP_Hx2y2z_S_C1005;
  abcd[76] = I_TWOBODYOVERLAP_Ixy4z_S_C1005+ABZ*I_TWOBODYOVERLAP_Hxy3z_S_C1005;
  abcd[77] = I_TWOBODYOVERLAP_Ix5z_S_C1005+ABZ*I_TWOBODYOVERLAP_Hx4z_S_C1005;
  abcd[78] = I_TWOBODYOVERLAP_I5yz_S_C1005+ABZ*I_TWOBODYOVERLAP_H5y_S_C1005;
  abcd[79] = I_TWOBODYOVERLAP_I4y2z_S_C1005+ABZ*I_TWOBODYOVERLAP_H4yz_S_C1005;
  abcd[80] = I_TWOBODYOVERLAP_I3y3z_S_C1005+ABZ*I_TWOBODYOVERLAP_H3y2z_S_C1005;
  abcd[81] = I_TWOBODYOVERLAP_I2y4z_S_C1005+ABZ*I_TWOBODYOVERLAP_H2y3z_S_C1005;
  abcd[82] = I_TWOBODYOVERLAP_Iy5z_S_C1005+ABZ*I_TWOBODYOVERLAP_Hy4z_S_C1005;
  abcd[83] = I_TWOBODYOVERLAP_I6z_S_C1005+ABZ*I_TWOBODYOVERLAP_H5z_S_C1005;
}
