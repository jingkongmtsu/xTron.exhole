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

//
// print out the information regarding of derivatives 
// here below we count on all of RHS integrals(including the repeat ones)
// this is used to simulate the FLOPS counting
// BRA1 as redundant position, total RHS integrals evaluated as: 4581
// BRA2 as redundant position, total RHS integrals evaluated as: 3078
// KET1 as redundant position, total RHS integrals evaluated as: 0
// KET2 as redundant position, total RHS integrals evaluated as: 0
// the redundant position is: BRA2
//

//
// @@@@ derivative position-direction information
// BRA1
// X
// Y
// Z
// ####

void hgp_os_kinetic_h_p_d1(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_KINETIC_I6x_Px_a = 0.0E0;
  Double I_KINETIC_I5xy_Px_a = 0.0E0;
  Double I_KINETIC_I5xz_Px_a = 0.0E0;
  Double I_KINETIC_I4x2y_Px_a = 0.0E0;
  Double I_KINETIC_I4xyz_Px_a = 0.0E0;
  Double I_KINETIC_I4x2z_Px_a = 0.0E0;
  Double I_KINETIC_I3x3y_Px_a = 0.0E0;
  Double I_KINETIC_I3x2yz_Px_a = 0.0E0;
  Double I_KINETIC_I3xy2z_Px_a = 0.0E0;
  Double I_KINETIC_I3x3z_Px_a = 0.0E0;
  Double I_KINETIC_I2x4y_Px_a = 0.0E0;
  Double I_KINETIC_I2x3yz_Px_a = 0.0E0;
  Double I_KINETIC_I2x2y2z_Px_a = 0.0E0;
  Double I_KINETIC_I2xy3z_Px_a = 0.0E0;
  Double I_KINETIC_I2x4z_Px_a = 0.0E0;
  Double I_KINETIC_Ix5y_Px_a = 0.0E0;
  Double I_KINETIC_Ix4yz_Px_a = 0.0E0;
  Double I_KINETIC_Ix3y2z_Px_a = 0.0E0;
  Double I_KINETIC_Ix2y3z_Px_a = 0.0E0;
  Double I_KINETIC_Ixy4z_Px_a = 0.0E0;
  Double I_KINETIC_Ix5z_Px_a = 0.0E0;
  Double I_KINETIC_I6y_Px_a = 0.0E0;
  Double I_KINETIC_I5yz_Px_a = 0.0E0;
  Double I_KINETIC_I4y2z_Px_a = 0.0E0;
  Double I_KINETIC_I3y3z_Px_a = 0.0E0;
  Double I_KINETIC_I2y4z_Px_a = 0.0E0;
  Double I_KINETIC_Iy5z_Px_a = 0.0E0;
  Double I_KINETIC_I6z_Px_a = 0.0E0;
  Double I_KINETIC_I6x_Py_a = 0.0E0;
  Double I_KINETIC_I5xy_Py_a = 0.0E0;
  Double I_KINETIC_I5xz_Py_a = 0.0E0;
  Double I_KINETIC_I4x2y_Py_a = 0.0E0;
  Double I_KINETIC_I4xyz_Py_a = 0.0E0;
  Double I_KINETIC_I4x2z_Py_a = 0.0E0;
  Double I_KINETIC_I3x3y_Py_a = 0.0E0;
  Double I_KINETIC_I3x2yz_Py_a = 0.0E0;
  Double I_KINETIC_I3xy2z_Py_a = 0.0E0;
  Double I_KINETIC_I3x3z_Py_a = 0.0E0;
  Double I_KINETIC_I2x4y_Py_a = 0.0E0;
  Double I_KINETIC_I2x3yz_Py_a = 0.0E0;
  Double I_KINETIC_I2x2y2z_Py_a = 0.0E0;
  Double I_KINETIC_I2xy3z_Py_a = 0.0E0;
  Double I_KINETIC_I2x4z_Py_a = 0.0E0;
  Double I_KINETIC_Ix5y_Py_a = 0.0E0;
  Double I_KINETIC_Ix4yz_Py_a = 0.0E0;
  Double I_KINETIC_Ix3y2z_Py_a = 0.0E0;
  Double I_KINETIC_Ix2y3z_Py_a = 0.0E0;
  Double I_KINETIC_Ixy4z_Py_a = 0.0E0;
  Double I_KINETIC_Ix5z_Py_a = 0.0E0;
  Double I_KINETIC_I6y_Py_a = 0.0E0;
  Double I_KINETIC_I5yz_Py_a = 0.0E0;
  Double I_KINETIC_I4y2z_Py_a = 0.0E0;
  Double I_KINETIC_I3y3z_Py_a = 0.0E0;
  Double I_KINETIC_I2y4z_Py_a = 0.0E0;
  Double I_KINETIC_Iy5z_Py_a = 0.0E0;
  Double I_KINETIC_I6z_Py_a = 0.0E0;
  Double I_KINETIC_I6x_Pz_a = 0.0E0;
  Double I_KINETIC_I5xy_Pz_a = 0.0E0;
  Double I_KINETIC_I5xz_Pz_a = 0.0E0;
  Double I_KINETIC_I4x2y_Pz_a = 0.0E0;
  Double I_KINETIC_I4xyz_Pz_a = 0.0E0;
  Double I_KINETIC_I4x2z_Pz_a = 0.0E0;
  Double I_KINETIC_I3x3y_Pz_a = 0.0E0;
  Double I_KINETIC_I3x2yz_Pz_a = 0.0E0;
  Double I_KINETIC_I3xy2z_Pz_a = 0.0E0;
  Double I_KINETIC_I3x3z_Pz_a = 0.0E0;
  Double I_KINETIC_I2x4y_Pz_a = 0.0E0;
  Double I_KINETIC_I2x3yz_Pz_a = 0.0E0;
  Double I_KINETIC_I2x2y2z_Pz_a = 0.0E0;
  Double I_KINETIC_I2xy3z_Pz_a = 0.0E0;
  Double I_KINETIC_I2x4z_Pz_a = 0.0E0;
  Double I_KINETIC_Ix5y_Pz_a = 0.0E0;
  Double I_KINETIC_Ix4yz_Pz_a = 0.0E0;
  Double I_KINETIC_Ix3y2z_Pz_a = 0.0E0;
  Double I_KINETIC_Ix2y3z_Pz_a = 0.0E0;
  Double I_KINETIC_Ixy4z_Pz_a = 0.0E0;
  Double I_KINETIC_Ix5z_Pz_a = 0.0E0;
  Double I_KINETIC_I6y_Pz_a = 0.0E0;
  Double I_KINETIC_I5yz_Pz_a = 0.0E0;
  Double I_KINETIC_I4y2z_Pz_a = 0.0E0;
  Double I_KINETIC_I3y3z_Pz_a = 0.0E0;
  Double I_KINETIC_I2y4z_Pz_a = 0.0E0;
  Double I_KINETIC_Iy5z_Pz_a = 0.0E0;
  Double I_KINETIC_I6z_Pz_a = 0.0E0;
  Double I_KINETIC_G4x_Px = 0.0E0;
  Double I_KINETIC_G3xy_Px = 0.0E0;
  Double I_KINETIC_G3xz_Px = 0.0E0;
  Double I_KINETIC_G2x2y_Px = 0.0E0;
  Double I_KINETIC_G2xyz_Px = 0.0E0;
  Double I_KINETIC_G2x2z_Px = 0.0E0;
  Double I_KINETIC_Gx3y_Px = 0.0E0;
  Double I_KINETIC_Gx2yz_Px = 0.0E0;
  Double I_KINETIC_Gxy2z_Px = 0.0E0;
  Double I_KINETIC_Gx3z_Px = 0.0E0;
  Double I_KINETIC_G4y_Px = 0.0E0;
  Double I_KINETIC_G3yz_Px = 0.0E0;
  Double I_KINETIC_G2y2z_Px = 0.0E0;
  Double I_KINETIC_Gy3z_Px = 0.0E0;
  Double I_KINETIC_G4z_Px = 0.0E0;
  Double I_KINETIC_G4x_Py = 0.0E0;
  Double I_KINETIC_G3xy_Py = 0.0E0;
  Double I_KINETIC_G3xz_Py = 0.0E0;
  Double I_KINETIC_G2x2y_Py = 0.0E0;
  Double I_KINETIC_G2xyz_Py = 0.0E0;
  Double I_KINETIC_G2x2z_Py = 0.0E0;
  Double I_KINETIC_Gx3y_Py = 0.0E0;
  Double I_KINETIC_Gx2yz_Py = 0.0E0;
  Double I_KINETIC_Gxy2z_Py = 0.0E0;
  Double I_KINETIC_Gx3z_Py = 0.0E0;
  Double I_KINETIC_G4y_Py = 0.0E0;
  Double I_KINETIC_G3yz_Py = 0.0E0;
  Double I_KINETIC_G2y2z_Py = 0.0E0;
  Double I_KINETIC_Gy3z_Py = 0.0E0;
  Double I_KINETIC_G4z_Py = 0.0E0;
  Double I_KINETIC_G4x_Pz = 0.0E0;
  Double I_KINETIC_G3xy_Pz = 0.0E0;
  Double I_KINETIC_G3xz_Pz = 0.0E0;
  Double I_KINETIC_G2x2y_Pz = 0.0E0;
  Double I_KINETIC_G2xyz_Pz = 0.0E0;
  Double I_KINETIC_G2x2z_Pz = 0.0E0;
  Double I_KINETIC_Gx3y_Pz = 0.0E0;
  Double I_KINETIC_Gx2yz_Pz = 0.0E0;
  Double I_KINETIC_Gxy2z_Pz = 0.0E0;
  Double I_KINETIC_Gx3z_Pz = 0.0E0;
  Double I_KINETIC_G4y_Pz = 0.0E0;
  Double I_KINETIC_G3yz_Pz = 0.0E0;
  Double I_KINETIC_G2y2z_Pz = 0.0E0;
  Double I_KINETIC_Gy3z_Pz = 0.0E0;
  Double I_KINETIC_G4z_Pz = 0.0E0;

  Double AB2 = (A[0]-B[0])*(A[0]-B[0])+(A[1]-B[1])*(A[1]-B[1])+(A[2]-B[2])*(A[2]-B[2]);
  for(UInt ip2=0; ip2<inp2; ip2++) {
    Double ic2   = icoe[ip2];
    Double onedz = iexp[ip2];
    Double zeta  = 1.0E0/onedz;
    Double zdiff = iexpdiff[ip2];
    Double alpha = 0.5E0*(zeta+zdiff);
    Double beta  = 0.5E0*(zeta-zdiff);
    Double xi    = alpha*beta*onedz;
    Double twoxi = 2.0E0*xi;
    Double fbra  = ifac[ip2];
    Double oned2z= 0.5E0*onedz;
    Double adz   = alpha*onedz;
    Double bdz   = beta*onedz;
    UInt offsetP = 3*ip2;
    Double PX    = P[offsetP  ];
    Double PY    = P[offsetP+1];
    Double PZ    = P[offsetP+2];
    Double PAX   = PX - A[0];
    Double PAY   = PY - A[1];
    Double PAZ   = PZ - A[2];
    Double PBX   = PX - B[0];
    Double PBY   = PY - B[1];
    Double PBZ   = PZ - B[2];
    Double I_KINETIC_S_S_vrr = ic2*fbra*xi*(3.0E0-twoxi*AB2);
    Double I_TWOBODYOVERLAP_S_S_vrr = ic2*fbra;
    if(fabs(I_KINETIC_S_S_vrr)<THRESHOLD_MATH) continue;


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
     * shell quartet name: SQ_TWOBODYOVERLAP_G_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_G4x_Px_vrr = PBX*I_TWOBODYOVERLAP_G4x_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Px_vrr = PBX*I_TWOBODYOVERLAP_G3xy_S_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Px_vrr = PBX*I_TWOBODYOVERLAP_G3xz_S_vrr+3*oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Px_vrr = PBX*I_TWOBODYOVERLAP_G2x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Px_vrr = PBX*I_TWOBODYOVERLAP_G2xyz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Px_vrr = PBX*I_TWOBODYOVERLAP_G2x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Px_vrr = PBX*I_TWOBODYOVERLAP_Gx3y_S_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Px_vrr = PBX*I_TWOBODYOVERLAP_Gx2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Px_vrr = PBX*I_TWOBODYOVERLAP_Gxy2z_S_vrr+oned2z*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Px_vrr = PBX*I_TWOBODYOVERLAP_Gx3z_S_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_Px_vrr = PBX*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Px_vrr = PBX*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Px_vrr = PBX*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Px_vrr = PBX*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_Px_vrr = PBX*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_G4x_Py_vrr = PBY*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Py_vrr = PBY*I_TWOBODYOVERLAP_G3xy_S_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Py_vrr = PBY*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Py_vrr = PBY*I_TWOBODYOVERLAP_G2x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Py_vrr = PBY*I_TWOBODYOVERLAP_G2xyz_S_vrr+oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Py_vrr = PBY*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Py_vrr = PBY*I_TWOBODYOVERLAP_Gx3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Py_vrr = PBY*I_TWOBODYOVERLAP_Gx2yz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Py_vrr = PBY*I_TWOBODYOVERLAP_Gxy2z_S_vrr+oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Py_vrr = PBY*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_Py_vrr = PBY*I_TWOBODYOVERLAP_G4y_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Py_vrr = PBY*I_TWOBODYOVERLAP_G3yz_S_vrr+3*oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Py_vrr = PBY*I_TWOBODYOVERLAP_G2y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Py_vrr = PBY*I_TWOBODYOVERLAP_Gy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_Py_vrr = PBY*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_G4x_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_G3xy_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_G3xz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G3xz_S_vrr+oned2z*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_G2xyz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G2xyz_S_vrr+oned2z*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_TWOBODYOVERLAP_G2x2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G2x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Gx2yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Gx2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_TWOBODYOVERLAP_Gxy2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Gxy2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_TWOBODYOVERLAP_Gx3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Gx3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_TWOBODYOVERLAP_G4y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_G3yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_TWOBODYOVERLAP_G2y2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G2y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Gy3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Gy3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_TWOBODYOVERLAP_G4z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_G4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_F3z_S_vrr;

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
     * shell quartet name: SQ_TWOBODYOVERLAP_I_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_I6x_Px_vrr = PBX*I_TWOBODYOVERLAP_I6x_S_vrr+6*oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xy_Px_vrr = PBX*I_TWOBODYOVERLAP_I5xy_S_vrr+5*oned2z*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_I5xz_Px_vrr = PBX*I_TWOBODYOVERLAP_I5xz_S_vrr+5*oned2z*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_Px_vrr = PBX*I_TWOBODYOVERLAP_I4x2y_S_vrr+4*oned2z*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_Px_vrr = PBX*I_TWOBODYOVERLAP_I4xyz_S_vrr+4*oned2z*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_Px_vrr = PBX*I_TWOBODYOVERLAP_I4x2z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_Px_vrr = PBX*I_TWOBODYOVERLAP_I3x3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_Px_vrr = PBX*I_TWOBODYOVERLAP_I3x2yz_S_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_Px_vrr = PBX*I_TWOBODYOVERLAP_I3xy2z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_Px_vrr = PBX*I_TWOBODYOVERLAP_I3x3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_Px_vrr = PBX*I_TWOBODYOVERLAP_I2x4y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_Px_vrr = PBX*I_TWOBODYOVERLAP_I2x3yz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_Px_vrr = PBX*I_TWOBODYOVERLAP_I2x2y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_Px_vrr = PBX*I_TWOBODYOVERLAP_I2xy3z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_Px_vrr = PBX*I_TWOBODYOVERLAP_I2x4z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_Px_vrr = PBX*I_TWOBODYOVERLAP_Ix5y_S_vrr+oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_Px_vrr = PBX*I_TWOBODYOVERLAP_Ix4yz_S_vrr+oned2z*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_Px_vrr = PBX*I_TWOBODYOVERLAP_Ix3y2z_S_vrr+oned2z*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_Px_vrr = PBX*I_TWOBODYOVERLAP_Ix2y3z_S_vrr+oned2z*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_Px_vrr = PBX*I_TWOBODYOVERLAP_Ixy4z_S_vrr+oned2z*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_Px_vrr = PBX*I_TWOBODYOVERLAP_Ix5z_S_vrr+oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6y_Px_vrr = PBX*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_I5yz_Px_vrr = PBX*I_TWOBODYOVERLAP_I5yz_S_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_Px_vrr = PBX*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_Px_vrr = PBX*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_Px_vrr = PBX*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_Px_vrr = PBX*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6z_Px_vrr = PBX*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_I6x_Py_vrr = PBY*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xy_Py_vrr = PBY*I_TWOBODYOVERLAP_I5xy_S_vrr+oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xz_Py_vrr = PBY*I_TWOBODYOVERLAP_I5xz_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_Py_vrr = PBY*I_TWOBODYOVERLAP_I4x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_Py_vrr = PBY*I_TWOBODYOVERLAP_I4xyz_S_vrr+oned2z*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_Py_vrr = PBY*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_Py_vrr = PBY*I_TWOBODYOVERLAP_I3x3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_Py_vrr = PBY*I_TWOBODYOVERLAP_I3x2yz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_Py_vrr = PBY*I_TWOBODYOVERLAP_I3xy2z_S_vrr+oned2z*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_Py_vrr = PBY*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_Py_vrr = PBY*I_TWOBODYOVERLAP_I2x4y_S_vrr+4*oned2z*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_Py_vrr = PBY*I_TWOBODYOVERLAP_I2x3yz_S_vrr+3*oned2z*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_Py_vrr = PBY*I_TWOBODYOVERLAP_I2x2y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_Py_vrr = PBY*I_TWOBODYOVERLAP_I2xy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_Py_vrr = PBY*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_Py_vrr = PBY*I_TWOBODYOVERLAP_Ix5y_S_vrr+5*oned2z*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_Py_vrr = PBY*I_TWOBODYOVERLAP_Ix4yz_S_vrr+4*oned2z*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_Py_vrr = PBY*I_TWOBODYOVERLAP_Ix3y2z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_Py_vrr = PBY*I_TWOBODYOVERLAP_Ix2y3z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_Py_vrr = PBY*I_TWOBODYOVERLAP_Ixy4z_S_vrr+oned2z*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_Py_vrr = PBY*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6y_Py_vrr = PBY*I_TWOBODYOVERLAP_I6y_S_vrr+6*oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_I5yz_Py_vrr = PBY*I_TWOBODYOVERLAP_I5yz_S_vrr+5*oned2z*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_Py_vrr = PBY*I_TWOBODYOVERLAP_I4y2z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_Py_vrr = PBY*I_TWOBODYOVERLAP_I3y3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_Py_vrr = PBY*I_TWOBODYOVERLAP_I2y4z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_Py_vrr = PBY*I_TWOBODYOVERLAP_Iy5z_S_vrr+oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_I6z_Py_vrr = PBY*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_I6x_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_I5xy_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_TWOBODYOVERLAP_I5xz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I5xz_S_vrr+oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I4xyz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I4xyz_S_vrr+oned2z*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_I4x2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I4x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    Double I_TWOBODYOVERLAP_I3x2yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I3x2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_I3xy2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I3xy2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    Double I_TWOBODYOVERLAP_I3x3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I3x3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    Double I_TWOBODYOVERLAP_I2x3yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I2x3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_I2x2y2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I2x2y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_TWOBODYOVERLAP_I2xy3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I2xy3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    Double I_TWOBODYOVERLAP_I2x4z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I2x4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_TWOBODYOVERLAP_Ix4yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Ix4yz_S_vrr+oned2z*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_Ix3y2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Ix3y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    Double I_TWOBODYOVERLAP_Ix2y3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Ix2y3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Ixy4z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Ixy4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    Double I_TWOBODYOVERLAP_Ix5z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Ix5z_S_vrr+5*oned2z*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_I6y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_I5yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I5yz_S_vrr+oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_I4y2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I4y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_TWOBODYOVERLAP_I3y3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I3y3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_I2y4z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I2y4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Iy5z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Iy5z_S_vrr+5*oned2z*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_TWOBODYOVERLAP_I6z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_I6z_S_vrr+6*oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_P_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_S_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_KINETIC_Px_S_vrr = PAX*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_KINETIC_Py_S_vrr = PAY*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_KINETIC_Pz_S_vrr = PAZ*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_D_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 2 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_KINETIC_S_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_KINETIC_D2x_S_vrr = PAX*I_KINETIC_Px_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_S_vrr-bdz*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_KINETIC_Dxy_S_vrr = PAY*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_Dxy_S_vrr;
    Double I_KINETIC_D2y_S_vrr = PAY*I_KINETIC_Py_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_S_vrr-bdz*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_KINETIC_D2z_S_vrr = PAZ*I_KINETIC_Pz_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_S_vrr-bdz*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_D_S
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_KINETIC_F3x_S_vrr = PAX*I_KINETIC_D2x_S_vrr+2*oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_S_vrr-2*bdz*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_KINETIC_F2xy_S_vrr = PAY*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_KINETIC_F2xz_S_vrr = PAZ*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_KINETIC_Fx2y_S_vrr = PAX*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_KINETIC_Fxyz_S_vrr = PAZ*I_KINETIC_Dxy_S_vrr+twoxi*I_TWOBODYOVERLAP_Fxyz_S_vrr;
    Double I_KINETIC_Fx2z_S_vrr = PAX*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_KINETIC_F3y_S_vrr = PAY*I_KINETIC_D2y_S_vrr+2*oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_S_vrr-2*bdz*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_KINETIC_F2yz_S_vrr = PAZ*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_S_vrr;
    Double I_KINETIC_Fy2z_S_vrr = PAY*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fy2z_S_vrr;
    Double I_KINETIC_F3z_S_vrr = PAZ*I_KINETIC_D2z_S_vrr+2*oned2z*I_KINETIC_Pz_S_vrr+twoxi*I_TWOBODYOVERLAP_F3z_S_vrr-2*bdz*I_TWOBODYOVERLAP_Pz_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_KINETIC_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     ************************************************************/
    Double I_KINETIC_G4x_S_vrr = PAX*I_KINETIC_F3x_S_vrr+3*oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_S_vrr-3*bdz*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_KINETIC_G3xy_S_vrr = PAY*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_KINETIC_G3xz_S_vrr = PAZ*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_KINETIC_G2x2y_S_vrr = PAY*I_KINETIC_F2xy_S_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_S_vrr-bdz*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_KINETIC_G2xyz_S_vrr = PAZ*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_KINETIC_G2x2z_S_vrr = PAZ*I_KINETIC_F2xz_S_vrr+oned2z*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_S_vrr-bdz*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_KINETIC_Gx3y_S_vrr = PAX*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_KINETIC_Gx2yz_S_vrr = PAZ*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_KINETIC_Gxy2z_S_vrr = PAY*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    Double I_KINETIC_Gx3z_S_vrr = PAX*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_KINETIC_G4y_S_vrr = PAY*I_KINETIC_F3y_S_vrr+3*oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_S_vrr-3*bdz*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_KINETIC_G3yz_S_vrr = PAZ*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_KINETIC_G2y2z_S_vrr = PAZ*I_KINETIC_F2yz_S_vrr+oned2z*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_S_vrr-bdz*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_KINETIC_Gy3z_S_vrr = PAY*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_KINETIC_G4z_S_vrr = PAZ*I_KINETIC_F3z_S_vrr+3*oned2z*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_S_vrr-3*bdz*I_TWOBODYOVERLAP_D2z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_G_S
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
     ************************************************************/
    Double I_KINETIC_G4x_Px_vrr = PBX*I_KINETIC_G4x_S_vrr+4*oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Px_vrr;
    Double I_KINETIC_G3xy_Px_vrr = PBX*I_KINETIC_G3xy_S_vrr+3*oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Px_vrr;
    Double I_KINETIC_G3xz_Px_vrr = PBX*I_KINETIC_G3xz_S_vrr+3*oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Px_vrr;
    Double I_KINETIC_G2x2y_Px_vrr = PBX*I_KINETIC_G2x2y_S_vrr+2*oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Px_vrr;
    Double I_KINETIC_G2xyz_Px_vrr = PBX*I_KINETIC_G2xyz_S_vrr+2*oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Px_vrr;
    Double I_KINETIC_G2x2z_Px_vrr = PBX*I_KINETIC_G2x2z_S_vrr+2*oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Px_vrr;
    Double I_KINETIC_Gx3y_Px_vrr = PBX*I_KINETIC_Gx3y_S_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Px_vrr;
    Double I_KINETIC_Gx2yz_Px_vrr = PBX*I_KINETIC_Gx2yz_S_vrr+oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Px_vrr;
    Double I_KINETIC_Gxy2z_Px_vrr = PBX*I_KINETIC_Gxy2z_S_vrr+oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Px_vrr;
    Double I_KINETIC_Gx3z_Px_vrr = PBX*I_KINETIC_Gx3z_S_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Px_vrr;
    Double I_KINETIC_G4y_Px_vrr = PBX*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Px_vrr;
    Double I_KINETIC_G3yz_Px_vrr = PBX*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Px_vrr;
    Double I_KINETIC_G2y2z_Px_vrr = PBX*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Px_vrr;
    Double I_KINETIC_Gy3z_Px_vrr = PBX*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Px_vrr;
    Double I_KINETIC_G4z_Px_vrr = PBX*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Px_vrr;
    Double I_KINETIC_G4x_Py_vrr = PBY*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Py_vrr;
    Double I_KINETIC_G3xy_Py_vrr = PBY*I_KINETIC_G3xy_S_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Py_vrr;
    Double I_KINETIC_G3xz_Py_vrr = PBY*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Py_vrr;
    Double I_KINETIC_G2x2y_Py_vrr = PBY*I_KINETIC_G2x2y_S_vrr+2*oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Py_vrr;
    Double I_KINETIC_G2xyz_Py_vrr = PBY*I_KINETIC_G2xyz_S_vrr+oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Py_vrr;
    Double I_KINETIC_G2x2z_Py_vrr = PBY*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Py_vrr;
    Double I_KINETIC_Gx3y_Py_vrr = PBY*I_KINETIC_Gx3y_S_vrr+3*oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Py_vrr;
    Double I_KINETIC_Gx2yz_Py_vrr = PBY*I_KINETIC_Gx2yz_S_vrr+2*oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Py_vrr;
    Double I_KINETIC_Gxy2z_Py_vrr = PBY*I_KINETIC_Gxy2z_S_vrr+oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Py_vrr;
    Double I_KINETIC_Gx3z_Py_vrr = PBY*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Py_vrr;
    Double I_KINETIC_G4y_Py_vrr = PBY*I_KINETIC_G4y_S_vrr+4*oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Py_vrr;
    Double I_KINETIC_G3yz_Py_vrr = PBY*I_KINETIC_G3yz_S_vrr+3*oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Py_vrr;
    Double I_KINETIC_G2y2z_Py_vrr = PBY*I_KINETIC_G2y2z_S_vrr+2*oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Py_vrr;
    Double I_KINETIC_Gy3z_Py_vrr = PBY*I_KINETIC_Gy3z_S_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Py_vrr;
    Double I_KINETIC_G4z_Py_vrr = PBY*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Py_vrr;
    Double I_KINETIC_G4x_Pz_vrr = PBZ*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_G4x_Pz_vrr;
    Double I_KINETIC_G3xy_Pz_vrr = PBZ*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xy_Pz_vrr;
    Double I_KINETIC_G3xz_Pz_vrr = PBZ*I_KINETIC_G3xz_S_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_G3xz_Pz_vrr;
    Double I_KINETIC_G2x2y_Pz_vrr = PBZ*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2y_Pz_vrr;
    Double I_KINETIC_G2xyz_Pz_vrr = PBZ*I_KINETIC_G2xyz_S_vrr+oned2z*I_KINETIC_F2xy_S_vrr+twoxi*I_TWOBODYOVERLAP_G2xyz_Pz_vrr;
    Double I_KINETIC_G2x2z_Pz_vrr = PBZ*I_KINETIC_G2x2z_S_vrr+2*oned2z*I_KINETIC_F2xz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2x2z_Pz_vrr;
    Double I_KINETIC_Gx3y_Pz_vrr = PBZ*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3y_Pz_vrr;
    Double I_KINETIC_Gx2yz_Pz_vrr = PBZ*I_KINETIC_Gx2yz_S_vrr+oned2z*I_KINETIC_Fx2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx2yz_Pz_vrr;
    Double I_KINETIC_Gxy2z_Pz_vrr = PBZ*I_KINETIC_Gxy2z_S_vrr+2*oned2z*I_KINETIC_Fxyz_S_vrr+twoxi*I_TWOBODYOVERLAP_Gxy2z_Pz_vrr;
    Double I_KINETIC_Gx3z_Pz_vrr = PBZ*I_KINETIC_Gx3z_S_vrr+3*oned2z*I_KINETIC_Fx2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gx3z_Pz_vrr;
    Double I_KINETIC_G4y_Pz_vrr = PBZ*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_G4y_Pz_vrr;
    Double I_KINETIC_G3yz_Pz_vrr = PBZ*I_KINETIC_G3yz_S_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_G3yz_Pz_vrr;
    Double I_KINETIC_G2y2z_Pz_vrr = PBZ*I_KINETIC_G2y2z_S_vrr+2*oned2z*I_KINETIC_F2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_G2y2z_Pz_vrr;
    Double I_KINETIC_Gy3z_Pz_vrr = PBZ*I_KINETIC_Gy3z_S_vrr+3*oned2z*I_KINETIC_Fy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Gy3z_Pz_vrr;
    Double I_KINETIC_G4z_Pz_vrr = PBZ*I_KINETIC_G4z_S_vrr+4*oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_G4z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_G_S
     * RHS shell quartet name: SQ_KINETIC_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     ************************************************************/
    Double I_KINETIC_H5x_S_vrr = PAX*I_KINETIC_G4x_S_vrr+4*oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_S_vrr-4*bdz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_H4xy_S_vrr = PAY*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_KINETIC_H4xz_S_vrr = PAZ*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_KINETIC_H3x2y_S_vrr = PAY*I_KINETIC_G3xy_S_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_S_vrr-bdz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_H3xyz_S_vrr = PAZ*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    Double I_KINETIC_H3x2z_S_vrr = PAZ*I_KINETIC_G3xz_S_vrr+oned2z*I_KINETIC_F3x_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_S_vrr-bdz*I_TWOBODYOVERLAP_F3x_S_vrr;
    Double I_KINETIC_H2x3y_S_vrr = PAX*I_KINETIC_Gx3y_S_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_S_vrr-bdz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_H2x2yz_S_vrr = PAZ*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    Double I_KINETIC_H2xy2z_S_vrr = PAY*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    Double I_KINETIC_H2x3z_S_vrr = PAX*I_KINETIC_Gx3z_S_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_S_vrr-bdz*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_KINETIC_Hx4y_S_vrr = PAX*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_KINETIC_Hx3yz_S_vrr = PAZ*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    Double I_KINETIC_Hx2y2z_S_vrr = PAX*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    Double I_KINETIC_Hxy3z_S_vrr = PAY*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    Double I_KINETIC_Hx4z_S_vrr = PAX*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_KINETIC_H5y_S_vrr = PAY*I_KINETIC_G4y_S_vrr+4*oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_S_vrr-4*bdz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_H4yz_S_vrr = PAZ*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_KINETIC_H3y2z_S_vrr = PAZ*I_KINETIC_G3yz_S_vrr+oned2z*I_KINETIC_F3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_S_vrr-bdz*I_TWOBODYOVERLAP_F3y_S_vrr;
    Double I_KINETIC_H2y3z_S_vrr = PAY*I_KINETIC_Gy3z_S_vrr+oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_S_vrr-bdz*I_TWOBODYOVERLAP_F3z_S_vrr;
    Double I_KINETIC_Hy4z_S_vrr = PAY*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_KINETIC_H5z_S_vrr = PAZ*I_KINETIC_G4z_S_vrr+4*oned2z*I_KINETIC_F3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_S_vrr-4*bdz*I_TWOBODYOVERLAP_F3z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_I_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_H_S
     * RHS shell quartet name: SQ_KINETIC_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     ************************************************************/
    Double I_KINETIC_I6x_S_vrr = PAX*I_KINETIC_H5x_S_vrr+5*oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_I6x_S_vrr-5*bdz*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_KINETIC_I5xy_S_vrr = PAY*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_KINETIC_I5xz_S_vrr = PAZ*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_S_vrr;
    Double I_KINETIC_I4x2y_S_vrr = PAY*I_KINETIC_H4xy_S_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_S_vrr-bdz*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_KINETIC_I4xyz_S_vrr = PAZ*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_S_vrr;
    Double I_KINETIC_I4x2z_S_vrr = PAZ*I_KINETIC_H4xz_S_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_S_vrr-bdz*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_KINETIC_I3x3y_S_vrr = PAY*I_KINETIC_H3x2y_S_vrr+2*oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_S_vrr-2*bdz*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_KINETIC_I3x2yz_S_vrr = PAZ*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    Double I_KINETIC_I3xy2z_S_vrr = PAY*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_S_vrr;
    Double I_KINETIC_I3x3z_S_vrr = PAZ*I_KINETIC_H3x2z_S_vrr+2*oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_S_vrr-2*bdz*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_KINETIC_I2x4y_S_vrr = PAX*I_KINETIC_Hx4y_S_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_S_vrr-bdz*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_KINETIC_I2x3yz_S_vrr = PAZ*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    Double I_KINETIC_I2x2y2z_S_vrr = PAZ*I_KINETIC_H2x2yz_S_vrr+oned2z*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_S_vrr-bdz*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_KINETIC_I2xy3z_S_vrr = PAY*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    Double I_KINETIC_I2x4z_S_vrr = PAX*I_KINETIC_Hx4z_S_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_S_vrr-bdz*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_KINETIC_Ix5y_S_vrr = PAX*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_KINETIC_Ix4yz_S_vrr = PAZ*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_S_vrr;
    Double I_KINETIC_Ix3y2z_S_vrr = PAX*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_S_vrr;
    Double I_KINETIC_Ix2y3z_S_vrr = PAX*I_KINETIC_H2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_S_vrr;
    Double I_KINETIC_Ixy4z_S_vrr = PAY*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_S_vrr;
    Double I_KINETIC_Ix5z_S_vrr = PAX*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_KINETIC_I6y_S_vrr = PAY*I_KINETIC_H5y_S_vrr+5*oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_I6y_S_vrr-5*bdz*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_KINETIC_I5yz_S_vrr = PAZ*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_S_vrr;
    Double I_KINETIC_I4y2z_S_vrr = PAZ*I_KINETIC_H4yz_S_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_S_vrr-bdz*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_KINETIC_I3y3z_S_vrr = PAZ*I_KINETIC_H3y2z_S_vrr+2*oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_S_vrr-2*bdz*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_KINETIC_I2y4z_S_vrr = PAY*I_KINETIC_Hy4z_S_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_S_vrr-bdz*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_KINETIC_Iy5z_S_vrr = PAY*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    Double I_KINETIC_I6z_S_vrr = PAZ*I_KINETIC_H5z_S_vrr+5*oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I6z_S_vrr-5*bdz*I_TWOBODYOVERLAP_G4z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_I_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_I_S
     * RHS shell quartet name: SQ_KINETIC_H_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_P
     ************************************************************/
    Double I_KINETIC_I6x_Px_vrr = PBX*I_KINETIC_I6x_S_vrr+6*oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_I6x_Px_vrr;
    Double I_KINETIC_I5xy_Px_vrr = PBX*I_KINETIC_I5xy_S_vrr+5*oned2z*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_Px_vrr;
    Double I_KINETIC_I5xz_Px_vrr = PBX*I_KINETIC_I5xz_S_vrr+5*oned2z*I_KINETIC_H4xz_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_Px_vrr;
    Double I_KINETIC_I4x2y_Px_vrr = PBX*I_KINETIC_I4x2y_S_vrr+4*oned2z*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_Px_vrr;
    Double I_KINETIC_I4xyz_Px_vrr = PBX*I_KINETIC_I4xyz_S_vrr+4*oned2z*I_KINETIC_H3xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_Px_vrr;
    Double I_KINETIC_I4x2z_Px_vrr = PBX*I_KINETIC_I4x2z_S_vrr+4*oned2z*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_Px_vrr;
    Double I_KINETIC_I3x3y_Px_vrr = PBX*I_KINETIC_I3x3y_S_vrr+3*oned2z*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_Px_vrr;
    Double I_KINETIC_I3x2yz_Px_vrr = PBX*I_KINETIC_I3x2yz_S_vrr+3*oned2z*I_KINETIC_H2x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_Px_vrr;
    Double I_KINETIC_I3xy2z_Px_vrr = PBX*I_KINETIC_I3xy2z_S_vrr+3*oned2z*I_KINETIC_H2xy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_Px_vrr;
    Double I_KINETIC_I3x3z_Px_vrr = PBX*I_KINETIC_I3x3z_S_vrr+3*oned2z*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_Px_vrr;
    Double I_KINETIC_I2x4y_Px_vrr = PBX*I_KINETIC_I2x4y_S_vrr+2*oned2z*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_Px_vrr;
    Double I_KINETIC_I2x3yz_Px_vrr = PBX*I_KINETIC_I2x3yz_S_vrr+2*oned2z*I_KINETIC_Hx3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_Px_vrr;
    Double I_KINETIC_I2x2y2z_Px_vrr = PBX*I_KINETIC_I2x2y2z_S_vrr+2*oned2z*I_KINETIC_Hx2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_Px_vrr;
    Double I_KINETIC_I2xy3z_Px_vrr = PBX*I_KINETIC_I2xy3z_S_vrr+2*oned2z*I_KINETIC_Hxy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_Px_vrr;
    Double I_KINETIC_I2x4z_Px_vrr = PBX*I_KINETIC_I2x4z_S_vrr+2*oned2z*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_Px_vrr;
    Double I_KINETIC_Ix5y_Px_vrr = PBX*I_KINETIC_Ix5y_S_vrr+oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_Px_vrr;
    Double I_KINETIC_Ix4yz_Px_vrr = PBX*I_KINETIC_Ix4yz_S_vrr+oned2z*I_KINETIC_H4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_Px_vrr;
    Double I_KINETIC_Ix3y2z_Px_vrr = PBX*I_KINETIC_Ix3y2z_S_vrr+oned2z*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_Px_vrr;
    Double I_KINETIC_Ix2y3z_Px_vrr = PBX*I_KINETIC_Ix2y3z_S_vrr+oned2z*I_KINETIC_H2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_Px_vrr;
    Double I_KINETIC_Ixy4z_Px_vrr = PBX*I_KINETIC_Ixy4z_S_vrr+oned2z*I_KINETIC_Hy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_Px_vrr;
    Double I_KINETIC_Ix5z_Px_vrr = PBX*I_KINETIC_Ix5z_S_vrr+oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_Px_vrr;
    Double I_KINETIC_I6y_Px_vrr = PBX*I_KINETIC_I6y_S_vrr+twoxi*I_TWOBODYOVERLAP_I6y_Px_vrr;
    Double I_KINETIC_I5yz_Px_vrr = PBX*I_KINETIC_I5yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_Px_vrr;
    Double I_KINETIC_I4y2z_Px_vrr = PBX*I_KINETIC_I4y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_Px_vrr;
    Double I_KINETIC_I3y3z_Px_vrr = PBX*I_KINETIC_I3y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_Px_vrr;
    Double I_KINETIC_I2y4z_Px_vrr = PBX*I_KINETIC_I2y4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_Px_vrr;
    Double I_KINETIC_Iy5z_Px_vrr = PBX*I_KINETIC_Iy5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_Px_vrr;
    Double I_KINETIC_I6z_Px_vrr = PBX*I_KINETIC_I6z_S_vrr+twoxi*I_TWOBODYOVERLAP_I6z_Px_vrr;
    Double I_KINETIC_I6x_Py_vrr = PBY*I_KINETIC_I6x_S_vrr+twoxi*I_TWOBODYOVERLAP_I6x_Py_vrr;
    Double I_KINETIC_I5xy_Py_vrr = PBY*I_KINETIC_I5xy_S_vrr+oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_Py_vrr;
    Double I_KINETIC_I5xz_Py_vrr = PBY*I_KINETIC_I5xz_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_Py_vrr;
    Double I_KINETIC_I4x2y_Py_vrr = PBY*I_KINETIC_I4x2y_S_vrr+2*oned2z*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_Py_vrr;
    Double I_KINETIC_I4xyz_Py_vrr = PBY*I_KINETIC_I4xyz_S_vrr+oned2z*I_KINETIC_H4xz_S_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_Py_vrr;
    Double I_KINETIC_I4x2z_Py_vrr = PBY*I_KINETIC_I4x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_Py_vrr;
    Double I_KINETIC_I3x3y_Py_vrr = PBY*I_KINETIC_I3x3y_S_vrr+3*oned2z*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_Py_vrr;
    Double I_KINETIC_I3x2yz_Py_vrr = PBY*I_KINETIC_I3x2yz_S_vrr+2*oned2z*I_KINETIC_H3xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_Py_vrr;
    Double I_KINETIC_I3xy2z_Py_vrr = PBY*I_KINETIC_I3xy2z_S_vrr+oned2z*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_Py_vrr;
    Double I_KINETIC_I3x3z_Py_vrr = PBY*I_KINETIC_I3x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_Py_vrr;
    Double I_KINETIC_I2x4y_Py_vrr = PBY*I_KINETIC_I2x4y_S_vrr+4*oned2z*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_Py_vrr;
    Double I_KINETIC_I2x3yz_Py_vrr = PBY*I_KINETIC_I2x3yz_S_vrr+3*oned2z*I_KINETIC_H2x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_Py_vrr;
    Double I_KINETIC_I2x2y2z_Py_vrr = PBY*I_KINETIC_I2x2y2z_S_vrr+2*oned2z*I_KINETIC_H2xy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_Py_vrr;
    Double I_KINETIC_I2xy3z_Py_vrr = PBY*I_KINETIC_I2xy3z_S_vrr+oned2z*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_Py_vrr;
    Double I_KINETIC_I2x4z_Py_vrr = PBY*I_KINETIC_I2x4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_Py_vrr;
    Double I_KINETIC_Ix5y_Py_vrr = PBY*I_KINETIC_Ix5y_S_vrr+5*oned2z*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_Py_vrr;
    Double I_KINETIC_Ix4yz_Py_vrr = PBY*I_KINETIC_Ix4yz_S_vrr+4*oned2z*I_KINETIC_Hx3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_Py_vrr;
    Double I_KINETIC_Ix3y2z_Py_vrr = PBY*I_KINETIC_Ix3y2z_S_vrr+3*oned2z*I_KINETIC_Hx2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_Py_vrr;
    Double I_KINETIC_Ix2y3z_Py_vrr = PBY*I_KINETIC_Ix2y3z_S_vrr+2*oned2z*I_KINETIC_Hxy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_Py_vrr;
    Double I_KINETIC_Ixy4z_Py_vrr = PBY*I_KINETIC_Ixy4z_S_vrr+oned2z*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_Py_vrr;
    Double I_KINETIC_Ix5z_Py_vrr = PBY*I_KINETIC_Ix5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_Py_vrr;
    Double I_KINETIC_I6y_Py_vrr = PBY*I_KINETIC_I6y_S_vrr+6*oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_I6y_Py_vrr;
    Double I_KINETIC_I5yz_Py_vrr = PBY*I_KINETIC_I5yz_S_vrr+5*oned2z*I_KINETIC_H4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_Py_vrr;
    Double I_KINETIC_I4y2z_Py_vrr = PBY*I_KINETIC_I4y2z_S_vrr+4*oned2z*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_Py_vrr;
    Double I_KINETIC_I3y3z_Py_vrr = PBY*I_KINETIC_I3y3z_S_vrr+3*oned2z*I_KINETIC_H2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_Py_vrr;
    Double I_KINETIC_I2y4z_Py_vrr = PBY*I_KINETIC_I2y4z_S_vrr+2*oned2z*I_KINETIC_Hy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_Py_vrr;
    Double I_KINETIC_Iy5z_Py_vrr = PBY*I_KINETIC_Iy5z_S_vrr+oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_Py_vrr;
    Double I_KINETIC_I6z_Py_vrr = PBY*I_KINETIC_I6z_S_vrr+twoxi*I_TWOBODYOVERLAP_I6z_Py_vrr;
    Double I_KINETIC_I6x_Pz_vrr = PBZ*I_KINETIC_I6x_S_vrr+twoxi*I_TWOBODYOVERLAP_I6x_Pz_vrr;
    Double I_KINETIC_I5xy_Pz_vrr = PBZ*I_KINETIC_I5xy_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xy_Pz_vrr;
    Double I_KINETIC_I5xz_Pz_vrr = PBZ*I_KINETIC_I5xz_S_vrr+oned2z*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_I5xz_Pz_vrr;
    Double I_KINETIC_I4x2y_Pz_vrr = PBZ*I_KINETIC_I4x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2y_Pz_vrr;
    Double I_KINETIC_I4xyz_Pz_vrr = PBZ*I_KINETIC_I4xyz_S_vrr+oned2z*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_I4xyz_Pz_vrr;
    Double I_KINETIC_I4x2z_Pz_vrr = PBZ*I_KINETIC_I4x2z_S_vrr+2*oned2z*I_KINETIC_H4xz_S_vrr+twoxi*I_TWOBODYOVERLAP_I4x2z_Pz_vrr;
    Double I_KINETIC_I3x3y_Pz_vrr = PBZ*I_KINETIC_I3x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3y_Pz_vrr;
    Double I_KINETIC_I3x2yz_Pz_vrr = PBZ*I_KINETIC_I3x2yz_S_vrr+oned2z*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x2yz_Pz_vrr;
    Double I_KINETIC_I3xy2z_Pz_vrr = PBZ*I_KINETIC_I3xy2z_S_vrr+2*oned2z*I_KINETIC_H3xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_I3xy2z_Pz_vrr;
    Double I_KINETIC_I3x3z_Pz_vrr = PBZ*I_KINETIC_I3x3z_S_vrr+3*oned2z*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3x3z_Pz_vrr;
    Double I_KINETIC_I2x4y_Pz_vrr = PBZ*I_KINETIC_I2x4y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4y_Pz_vrr;
    Double I_KINETIC_I2x3yz_Pz_vrr = PBZ*I_KINETIC_I2x3yz_S_vrr+oned2z*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x3yz_Pz_vrr;
    Double I_KINETIC_I2x2y2z_Pz_vrr = PBZ*I_KINETIC_I2x2y2z_S_vrr+2*oned2z*I_KINETIC_H2x2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x2y2z_Pz_vrr;
    Double I_KINETIC_I2xy3z_Pz_vrr = PBZ*I_KINETIC_I2xy3z_S_vrr+3*oned2z*I_KINETIC_H2xy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2xy3z_Pz_vrr;
    Double I_KINETIC_I2x4z_Pz_vrr = PBZ*I_KINETIC_I2x4z_S_vrr+4*oned2z*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2x4z_Pz_vrr;
    Double I_KINETIC_Ix5y_Pz_vrr = PBZ*I_KINETIC_Ix5y_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5y_Pz_vrr;
    Double I_KINETIC_Ix4yz_Pz_vrr = PBZ*I_KINETIC_Ix4yz_S_vrr+oned2z*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix4yz_Pz_vrr;
    Double I_KINETIC_Ix3y2z_Pz_vrr = PBZ*I_KINETIC_Ix3y2z_S_vrr+2*oned2z*I_KINETIC_Hx3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix3y2z_Pz_vrr;
    Double I_KINETIC_Ix2y3z_Pz_vrr = PBZ*I_KINETIC_Ix2y3z_S_vrr+3*oned2z*I_KINETIC_Hx2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix2y3z_Pz_vrr;
    Double I_KINETIC_Ixy4z_Pz_vrr = PBZ*I_KINETIC_Ixy4z_S_vrr+4*oned2z*I_KINETIC_Hxy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ixy4z_Pz_vrr;
    Double I_KINETIC_Ix5z_Pz_vrr = PBZ*I_KINETIC_Ix5z_S_vrr+5*oned2z*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Ix5z_Pz_vrr;
    Double I_KINETIC_I6y_Pz_vrr = PBZ*I_KINETIC_I6y_S_vrr+twoxi*I_TWOBODYOVERLAP_I6y_Pz_vrr;
    Double I_KINETIC_I5yz_Pz_vrr = PBZ*I_KINETIC_I5yz_S_vrr+oned2z*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_I5yz_Pz_vrr;
    Double I_KINETIC_I4y2z_Pz_vrr = PBZ*I_KINETIC_I4y2z_S_vrr+2*oned2z*I_KINETIC_H4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_I4y2z_Pz_vrr;
    Double I_KINETIC_I3y3z_Pz_vrr = PBZ*I_KINETIC_I3y3z_S_vrr+3*oned2z*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_I3y3z_Pz_vrr;
    Double I_KINETIC_I2y4z_Pz_vrr = PBZ*I_KINETIC_I2y4z_S_vrr+4*oned2z*I_KINETIC_H2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_I2y4z_Pz_vrr;
    Double I_KINETIC_Iy5z_Pz_vrr = PBZ*I_KINETIC_Iy5z_S_vrr+5*oned2z*I_KINETIC_Hy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Iy5z_Pz_vrr;
    Double I_KINETIC_I6z_Pz_vrr = PBZ*I_KINETIC_I6z_S_vrr+6*oned2z*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_I6z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_I_P_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_KINETIC_I_P_a_coefs = alpha;
    I_KINETIC_I6x_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I6x_Px_vrr;
    I_KINETIC_I5xy_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I5xy_Px_vrr;
    I_KINETIC_I5xz_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I5xz_Px_vrr;
    I_KINETIC_I4x2y_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I4x2y_Px_vrr;
    I_KINETIC_I4xyz_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I4xyz_Px_vrr;
    I_KINETIC_I4x2z_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I4x2z_Px_vrr;
    I_KINETIC_I3x3y_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I3x3y_Px_vrr;
    I_KINETIC_I3x2yz_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I3x2yz_Px_vrr;
    I_KINETIC_I3xy2z_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I3xy2z_Px_vrr;
    I_KINETIC_I3x3z_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I3x3z_Px_vrr;
    I_KINETIC_I2x4y_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I2x4y_Px_vrr;
    I_KINETIC_I2x3yz_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I2x3yz_Px_vrr;
    I_KINETIC_I2x2y2z_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I2x2y2z_Px_vrr;
    I_KINETIC_I2xy3z_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I2xy3z_Px_vrr;
    I_KINETIC_I2x4z_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I2x4z_Px_vrr;
    I_KINETIC_Ix5y_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Ix5y_Px_vrr;
    I_KINETIC_Ix4yz_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Ix4yz_Px_vrr;
    I_KINETIC_Ix3y2z_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Ix3y2z_Px_vrr;
    I_KINETIC_Ix2y3z_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Ix2y3z_Px_vrr;
    I_KINETIC_Ixy4z_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Ixy4z_Px_vrr;
    I_KINETIC_Ix5z_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Ix5z_Px_vrr;
    I_KINETIC_I6y_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I6y_Px_vrr;
    I_KINETIC_I5yz_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I5yz_Px_vrr;
    I_KINETIC_I4y2z_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I4y2z_Px_vrr;
    I_KINETIC_I3y3z_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I3y3z_Px_vrr;
    I_KINETIC_I2y4z_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I2y4z_Px_vrr;
    I_KINETIC_Iy5z_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Iy5z_Px_vrr;
    I_KINETIC_I6z_Px_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I6z_Px_vrr;
    I_KINETIC_I6x_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I6x_Py_vrr;
    I_KINETIC_I5xy_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I5xy_Py_vrr;
    I_KINETIC_I5xz_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I5xz_Py_vrr;
    I_KINETIC_I4x2y_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I4x2y_Py_vrr;
    I_KINETIC_I4xyz_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I4xyz_Py_vrr;
    I_KINETIC_I4x2z_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I4x2z_Py_vrr;
    I_KINETIC_I3x3y_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I3x3y_Py_vrr;
    I_KINETIC_I3x2yz_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I3x2yz_Py_vrr;
    I_KINETIC_I3xy2z_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I3xy2z_Py_vrr;
    I_KINETIC_I3x3z_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I3x3z_Py_vrr;
    I_KINETIC_I2x4y_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I2x4y_Py_vrr;
    I_KINETIC_I2x3yz_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I2x3yz_Py_vrr;
    I_KINETIC_I2x2y2z_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I2x2y2z_Py_vrr;
    I_KINETIC_I2xy3z_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I2xy3z_Py_vrr;
    I_KINETIC_I2x4z_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I2x4z_Py_vrr;
    I_KINETIC_Ix5y_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Ix5y_Py_vrr;
    I_KINETIC_Ix4yz_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Ix4yz_Py_vrr;
    I_KINETIC_Ix3y2z_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Ix3y2z_Py_vrr;
    I_KINETIC_Ix2y3z_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Ix2y3z_Py_vrr;
    I_KINETIC_Ixy4z_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Ixy4z_Py_vrr;
    I_KINETIC_Ix5z_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Ix5z_Py_vrr;
    I_KINETIC_I6y_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I6y_Py_vrr;
    I_KINETIC_I5yz_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I5yz_Py_vrr;
    I_KINETIC_I4y2z_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I4y2z_Py_vrr;
    I_KINETIC_I3y3z_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I3y3z_Py_vrr;
    I_KINETIC_I2y4z_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I2y4z_Py_vrr;
    I_KINETIC_Iy5z_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Iy5z_Py_vrr;
    I_KINETIC_I6z_Py_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I6z_Py_vrr;
    I_KINETIC_I6x_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I6x_Pz_vrr;
    I_KINETIC_I5xy_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I5xy_Pz_vrr;
    I_KINETIC_I5xz_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I5xz_Pz_vrr;
    I_KINETIC_I4x2y_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I4x2y_Pz_vrr;
    I_KINETIC_I4xyz_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I4xyz_Pz_vrr;
    I_KINETIC_I4x2z_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I4x2z_Pz_vrr;
    I_KINETIC_I3x3y_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I3x3y_Pz_vrr;
    I_KINETIC_I3x2yz_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I3x2yz_Pz_vrr;
    I_KINETIC_I3xy2z_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I3xy2z_Pz_vrr;
    I_KINETIC_I3x3z_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I3x3z_Pz_vrr;
    I_KINETIC_I2x4y_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I2x4y_Pz_vrr;
    I_KINETIC_I2x3yz_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I2x3yz_Pz_vrr;
    I_KINETIC_I2x2y2z_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I2x2y2z_Pz_vrr;
    I_KINETIC_I2xy3z_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I2xy3z_Pz_vrr;
    I_KINETIC_I2x4z_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I2x4z_Pz_vrr;
    I_KINETIC_Ix5y_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Ix5y_Pz_vrr;
    I_KINETIC_Ix4yz_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Ix4yz_Pz_vrr;
    I_KINETIC_Ix3y2z_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Ix3y2z_Pz_vrr;
    I_KINETIC_Ix2y3z_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Ix2y3z_Pz_vrr;
    I_KINETIC_Ixy4z_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Ixy4z_Pz_vrr;
    I_KINETIC_Ix5z_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Ix5z_Pz_vrr;
    I_KINETIC_I6y_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I6y_Pz_vrr;
    I_KINETIC_I5yz_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I5yz_Pz_vrr;
    I_KINETIC_I4y2z_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I4y2z_Pz_vrr;
    I_KINETIC_I3y3z_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I3y3z_Pz_vrr;
    I_KINETIC_I2y4z_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I2y4z_Pz_vrr;
    I_KINETIC_Iy5z_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_Iy5z_Pz_vrr;
    I_KINETIC_I6z_Pz_a += SQ_KINETIC_I_P_a_coefs*I_KINETIC_I6z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_G_P
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_KINETIC_G4x_Px += I_KINETIC_G4x_Px_vrr;
    I_KINETIC_G3xy_Px += I_KINETIC_G3xy_Px_vrr;
    I_KINETIC_G3xz_Px += I_KINETIC_G3xz_Px_vrr;
    I_KINETIC_G2x2y_Px += I_KINETIC_G2x2y_Px_vrr;
    I_KINETIC_G2xyz_Px += I_KINETIC_G2xyz_Px_vrr;
    I_KINETIC_G2x2z_Px += I_KINETIC_G2x2z_Px_vrr;
    I_KINETIC_Gx3y_Px += I_KINETIC_Gx3y_Px_vrr;
    I_KINETIC_Gx2yz_Px += I_KINETIC_Gx2yz_Px_vrr;
    I_KINETIC_Gxy2z_Px += I_KINETIC_Gxy2z_Px_vrr;
    I_KINETIC_Gx3z_Px += I_KINETIC_Gx3z_Px_vrr;
    I_KINETIC_G4y_Px += I_KINETIC_G4y_Px_vrr;
    I_KINETIC_G3yz_Px += I_KINETIC_G3yz_Px_vrr;
    I_KINETIC_G2y2z_Px += I_KINETIC_G2y2z_Px_vrr;
    I_KINETIC_Gy3z_Px += I_KINETIC_Gy3z_Px_vrr;
    I_KINETIC_G4z_Px += I_KINETIC_G4z_Px_vrr;
    I_KINETIC_G4x_Py += I_KINETIC_G4x_Py_vrr;
    I_KINETIC_G3xy_Py += I_KINETIC_G3xy_Py_vrr;
    I_KINETIC_G3xz_Py += I_KINETIC_G3xz_Py_vrr;
    I_KINETIC_G2x2y_Py += I_KINETIC_G2x2y_Py_vrr;
    I_KINETIC_G2xyz_Py += I_KINETIC_G2xyz_Py_vrr;
    I_KINETIC_G2x2z_Py += I_KINETIC_G2x2z_Py_vrr;
    I_KINETIC_Gx3y_Py += I_KINETIC_Gx3y_Py_vrr;
    I_KINETIC_Gx2yz_Py += I_KINETIC_Gx2yz_Py_vrr;
    I_KINETIC_Gxy2z_Py += I_KINETIC_Gxy2z_Py_vrr;
    I_KINETIC_Gx3z_Py += I_KINETIC_Gx3z_Py_vrr;
    I_KINETIC_G4y_Py += I_KINETIC_G4y_Py_vrr;
    I_KINETIC_G3yz_Py += I_KINETIC_G3yz_Py_vrr;
    I_KINETIC_G2y2z_Py += I_KINETIC_G2y2z_Py_vrr;
    I_KINETIC_Gy3z_Py += I_KINETIC_Gy3z_Py_vrr;
    I_KINETIC_G4z_Py += I_KINETIC_G4z_Py_vrr;
    I_KINETIC_G4x_Pz += I_KINETIC_G4x_Pz_vrr;
    I_KINETIC_G3xy_Pz += I_KINETIC_G3xy_Pz_vrr;
    I_KINETIC_G3xz_Pz += I_KINETIC_G3xz_Pz_vrr;
    I_KINETIC_G2x2y_Pz += I_KINETIC_G2x2y_Pz_vrr;
    I_KINETIC_G2xyz_Pz += I_KINETIC_G2xyz_Pz_vrr;
    I_KINETIC_G2x2z_Pz += I_KINETIC_G2x2z_Pz_vrr;
    I_KINETIC_Gx3y_Pz += I_KINETIC_Gx3y_Pz_vrr;
    I_KINETIC_Gx2yz_Pz += I_KINETIC_Gx2yz_Pz_vrr;
    I_KINETIC_Gxy2z_Pz += I_KINETIC_Gxy2z_Pz_vrr;
    I_KINETIC_Gx3z_Pz += I_KINETIC_Gx3z_Pz_vrr;
    I_KINETIC_G4y_Pz += I_KINETIC_G4y_Pz_vrr;
    I_KINETIC_G3yz_Pz += I_KINETIC_G3yz_Pz_vrr;
    I_KINETIC_G2y2z_Pz += I_KINETIC_G2y2z_Pz_vrr;
    I_KINETIC_Gy3z_Pz += I_KINETIC_Gy3z_Pz_vrr;
    I_KINETIC_G4z_Pz += I_KINETIC_G4z_Pz_vrr;
  }

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_P_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_I_P_a
   * RHS shell quartet name: SQ_KINETIC_G_P
   ************************************************************/
  abcd[0] = 2.0E0*I_KINETIC_I6x_Px_a-5*I_KINETIC_G4x_Px;
  abcd[1] = 2.0E0*I_KINETIC_I5xy_Px_a-4*I_KINETIC_G3xy_Px;
  abcd[2] = 2.0E0*I_KINETIC_I5xz_Px_a-4*I_KINETIC_G3xz_Px;
  abcd[3] = 2.0E0*I_KINETIC_I4x2y_Px_a-3*I_KINETIC_G2x2y_Px;
  abcd[4] = 2.0E0*I_KINETIC_I4xyz_Px_a-3*I_KINETIC_G2xyz_Px;
  abcd[5] = 2.0E0*I_KINETIC_I4x2z_Px_a-3*I_KINETIC_G2x2z_Px;
  abcd[6] = 2.0E0*I_KINETIC_I3x3y_Px_a-2*I_KINETIC_Gx3y_Px;
  abcd[7] = 2.0E0*I_KINETIC_I3x2yz_Px_a-2*I_KINETIC_Gx2yz_Px;
  abcd[8] = 2.0E0*I_KINETIC_I3xy2z_Px_a-2*I_KINETIC_Gxy2z_Px;
  abcd[9] = 2.0E0*I_KINETIC_I3x3z_Px_a-2*I_KINETIC_Gx3z_Px;
  abcd[10] = 2.0E0*I_KINETIC_I2x4y_Px_a-1*I_KINETIC_G4y_Px;
  abcd[11] = 2.0E0*I_KINETIC_I2x3yz_Px_a-1*I_KINETIC_G3yz_Px;
  abcd[12] = 2.0E0*I_KINETIC_I2x2y2z_Px_a-1*I_KINETIC_G2y2z_Px;
  abcd[13] = 2.0E0*I_KINETIC_I2xy3z_Px_a-1*I_KINETIC_Gy3z_Px;
  abcd[14] = 2.0E0*I_KINETIC_I2x4z_Px_a-1*I_KINETIC_G4z_Px;
  abcd[15] = 2.0E0*I_KINETIC_Ix5y_Px_a;
  abcd[16] = 2.0E0*I_KINETIC_Ix4yz_Px_a;
  abcd[17] = 2.0E0*I_KINETIC_Ix3y2z_Px_a;
  abcd[18] = 2.0E0*I_KINETIC_Ix2y3z_Px_a;
  abcd[19] = 2.0E0*I_KINETIC_Ixy4z_Px_a;
  abcd[20] = 2.0E0*I_KINETIC_Ix5z_Px_a;
  abcd[21] = 2.0E0*I_KINETIC_I6x_Py_a-5*I_KINETIC_G4x_Py;
  abcd[22] = 2.0E0*I_KINETIC_I5xy_Py_a-4*I_KINETIC_G3xy_Py;
  abcd[23] = 2.0E0*I_KINETIC_I5xz_Py_a-4*I_KINETIC_G3xz_Py;
  abcd[24] = 2.0E0*I_KINETIC_I4x2y_Py_a-3*I_KINETIC_G2x2y_Py;
  abcd[25] = 2.0E0*I_KINETIC_I4xyz_Py_a-3*I_KINETIC_G2xyz_Py;
  abcd[26] = 2.0E0*I_KINETIC_I4x2z_Py_a-3*I_KINETIC_G2x2z_Py;
  abcd[27] = 2.0E0*I_KINETIC_I3x3y_Py_a-2*I_KINETIC_Gx3y_Py;
  abcd[28] = 2.0E0*I_KINETIC_I3x2yz_Py_a-2*I_KINETIC_Gx2yz_Py;
  abcd[29] = 2.0E0*I_KINETIC_I3xy2z_Py_a-2*I_KINETIC_Gxy2z_Py;
  abcd[30] = 2.0E0*I_KINETIC_I3x3z_Py_a-2*I_KINETIC_Gx3z_Py;
  abcd[31] = 2.0E0*I_KINETIC_I2x4y_Py_a-1*I_KINETIC_G4y_Py;
  abcd[32] = 2.0E0*I_KINETIC_I2x3yz_Py_a-1*I_KINETIC_G3yz_Py;
  abcd[33] = 2.0E0*I_KINETIC_I2x2y2z_Py_a-1*I_KINETIC_G2y2z_Py;
  abcd[34] = 2.0E0*I_KINETIC_I2xy3z_Py_a-1*I_KINETIC_Gy3z_Py;
  abcd[35] = 2.0E0*I_KINETIC_I2x4z_Py_a-1*I_KINETIC_G4z_Py;
  abcd[36] = 2.0E0*I_KINETIC_Ix5y_Py_a;
  abcd[37] = 2.0E0*I_KINETIC_Ix4yz_Py_a;
  abcd[38] = 2.0E0*I_KINETIC_Ix3y2z_Py_a;
  abcd[39] = 2.0E0*I_KINETIC_Ix2y3z_Py_a;
  abcd[40] = 2.0E0*I_KINETIC_Ixy4z_Py_a;
  abcd[41] = 2.0E0*I_KINETIC_Ix5z_Py_a;
  abcd[42] = 2.0E0*I_KINETIC_I6x_Pz_a-5*I_KINETIC_G4x_Pz;
  abcd[43] = 2.0E0*I_KINETIC_I5xy_Pz_a-4*I_KINETIC_G3xy_Pz;
  abcd[44] = 2.0E0*I_KINETIC_I5xz_Pz_a-4*I_KINETIC_G3xz_Pz;
  abcd[45] = 2.0E0*I_KINETIC_I4x2y_Pz_a-3*I_KINETIC_G2x2y_Pz;
  abcd[46] = 2.0E0*I_KINETIC_I4xyz_Pz_a-3*I_KINETIC_G2xyz_Pz;
  abcd[47] = 2.0E0*I_KINETIC_I4x2z_Pz_a-3*I_KINETIC_G2x2z_Pz;
  abcd[48] = 2.0E0*I_KINETIC_I3x3y_Pz_a-2*I_KINETIC_Gx3y_Pz;
  abcd[49] = 2.0E0*I_KINETIC_I3x2yz_Pz_a-2*I_KINETIC_Gx2yz_Pz;
  abcd[50] = 2.0E0*I_KINETIC_I3xy2z_Pz_a-2*I_KINETIC_Gxy2z_Pz;
  abcd[51] = 2.0E0*I_KINETIC_I3x3z_Pz_a-2*I_KINETIC_Gx3z_Pz;
  abcd[52] = 2.0E0*I_KINETIC_I2x4y_Pz_a-1*I_KINETIC_G4y_Pz;
  abcd[53] = 2.0E0*I_KINETIC_I2x3yz_Pz_a-1*I_KINETIC_G3yz_Pz;
  abcd[54] = 2.0E0*I_KINETIC_I2x2y2z_Pz_a-1*I_KINETIC_G2y2z_Pz;
  abcd[55] = 2.0E0*I_KINETIC_I2xy3z_Pz_a-1*I_KINETIC_Gy3z_Pz;
  abcd[56] = 2.0E0*I_KINETIC_I2x4z_Pz_a-1*I_KINETIC_G4z_Pz;
  abcd[57] = 2.0E0*I_KINETIC_Ix5y_Pz_a;
  abcd[58] = 2.0E0*I_KINETIC_Ix4yz_Pz_a;
  abcd[59] = 2.0E0*I_KINETIC_Ix3y2z_Pz_a;
  abcd[60] = 2.0E0*I_KINETIC_Ix2y3z_Pz_a;
  abcd[61] = 2.0E0*I_KINETIC_Ixy4z_Pz_a;
  abcd[62] = 2.0E0*I_KINETIC_Ix5z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_P_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_I_P_a
   * RHS shell quartet name: SQ_KINETIC_G_P
   ************************************************************/
  abcd[63] = 2.0E0*I_KINETIC_I5xy_Px_a;
  abcd[64] = 2.0E0*I_KINETIC_I4x2y_Px_a-1*I_KINETIC_G4x_Px;
  abcd[65] = 2.0E0*I_KINETIC_I4xyz_Px_a;
  abcd[66] = 2.0E0*I_KINETIC_I3x3y_Px_a-2*I_KINETIC_G3xy_Px;
  abcd[67] = 2.0E0*I_KINETIC_I3x2yz_Px_a-1*I_KINETIC_G3xz_Px;
  abcd[68] = 2.0E0*I_KINETIC_I3xy2z_Px_a;
  abcd[69] = 2.0E0*I_KINETIC_I2x4y_Px_a-3*I_KINETIC_G2x2y_Px;
  abcd[70] = 2.0E0*I_KINETIC_I2x3yz_Px_a-2*I_KINETIC_G2xyz_Px;
  abcd[71] = 2.0E0*I_KINETIC_I2x2y2z_Px_a-1*I_KINETIC_G2x2z_Px;
  abcd[72] = 2.0E0*I_KINETIC_I2xy3z_Px_a;
  abcd[73] = 2.0E0*I_KINETIC_Ix5y_Px_a-4*I_KINETIC_Gx3y_Px;
  abcd[74] = 2.0E0*I_KINETIC_Ix4yz_Px_a-3*I_KINETIC_Gx2yz_Px;
  abcd[75] = 2.0E0*I_KINETIC_Ix3y2z_Px_a-2*I_KINETIC_Gxy2z_Px;
  abcd[76] = 2.0E0*I_KINETIC_Ix2y3z_Px_a-1*I_KINETIC_Gx3z_Px;
  abcd[77] = 2.0E0*I_KINETIC_Ixy4z_Px_a;
  abcd[78] = 2.0E0*I_KINETIC_I6y_Px_a-5*I_KINETIC_G4y_Px;
  abcd[79] = 2.0E0*I_KINETIC_I5yz_Px_a-4*I_KINETIC_G3yz_Px;
  abcd[80] = 2.0E0*I_KINETIC_I4y2z_Px_a-3*I_KINETIC_G2y2z_Px;
  abcd[81] = 2.0E0*I_KINETIC_I3y3z_Px_a-2*I_KINETIC_Gy3z_Px;
  abcd[82] = 2.0E0*I_KINETIC_I2y4z_Px_a-1*I_KINETIC_G4z_Px;
  abcd[83] = 2.0E0*I_KINETIC_Iy5z_Px_a;
  abcd[84] = 2.0E0*I_KINETIC_I5xy_Py_a;
  abcd[85] = 2.0E0*I_KINETIC_I4x2y_Py_a-1*I_KINETIC_G4x_Py;
  abcd[86] = 2.0E0*I_KINETIC_I4xyz_Py_a;
  abcd[87] = 2.0E0*I_KINETIC_I3x3y_Py_a-2*I_KINETIC_G3xy_Py;
  abcd[88] = 2.0E0*I_KINETIC_I3x2yz_Py_a-1*I_KINETIC_G3xz_Py;
  abcd[89] = 2.0E0*I_KINETIC_I3xy2z_Py_a;
  abcd[90] = 2.0E0*I_KINETIC_I2x4y_Py_a-3*I_KINETIC_G2x2y_Py;
  abcd[91] = 2.0E0*I_KINETIC_I2x3yz_Py_a-2*I_KINETIC_G2xyz_Py;
  abcd[92] = 2.0E0*I_KINETIC_I2x2y2z_Py_a-1*I_KINETIC_G2x2z_Py;
  abcd[93] = 2.0E0*I_KINETIC_I2xy3z_Py_a;
  abcd[94] = 2.0E0*I_KINETIC_Ix5y_Py_a-4*I_KINETIC_Gx3y_Py;
  abcd[95] = 2.0E0*I_KINETIC_Ix4yz_Py_a-3*I_KINETIC_Gx2yz_Py;
  abcd[96] = 2.0E0*I_KINETIC_Ix3y2z_Py_a-2*I_KINETIC_Gxy2z_Py;
  abcd[97] = 2.0E0*I_KINETIC_Ix2y3z_Py_a-1*I_KINETIC_Gx3z_Py;
  abcd[98] = 2.0E0*I_KINETIC_Ixy4z_Py_a;
  abcd[99] = 2.0E0*I_KINETIC_I6y_Py_a-5*I_KINETIC_G4y_Py;
  abcd[100] = 2.0E0*I_KINETIC_I5yz_Py_a-4*I_KINETIC_G3yz_Py;
  abcd[101] = 2.0E0*I_KINETIC_I4y2z_Py_a-3*I_KINETIC_G2y2z_Py;
  abcd[102] = 2.0E0*I_KINETIC_I3y3z_Py_a-2*I_KINETIC_Gy3z_Py;
  abcd[103] = 2.0E0*I_KINETIC_I2y4z_Py_a-1*I_KINETIC_G4z_Py;
  abcd[104] = 2.0E0*I_KINETIC_Iy5z_Py_a;
  abcd[105] = 2.0E0*I_KINETIC_I5xy_Pz_a;
  abcd[106] = 2.0E0*I_KINETIC_I4x2y_Pz_a-1*I_KINETIC_G4x_Pz;
  abcd[107] = 2.0E0*I_KINETIC_I4xyz_Pz_a;
  abcd[108] = 2.0E0*I_KINETIC_I3x3y_Pz_a-2*I_KINETIC_G3xy_Pz;
  abcd[109] = 2.0E0*I_KINETIC_I3x2yz_Pz_a-1*I_KINETIC_G3xz_Pz;
  abcd[110] = 2.0E0*I_KINETIC_I3xy2z_Pz_a;
  abcd[111] = 2.0E0*I_KINETIC_I2x4y_Pz_a-3*I_KINETIC_G2x2y_Pz;
  abcd[112] = 2.0E0*I_KINETIC_I2x3yz_Pz_a-2*I_KINETIC_G2xyz_Pz;
  abcd[113] = 2.0E0*I_KINETIC_I2x2y2z_Pz_a-1*I_KINETIC_G2x2z_Pz;
  abcd[114] = 2.0E0*I_KINETIC_I2xy3z_Pz_a;
  abcd[115] = 2.0E0*I_KINETIC_Ix5y_Pz_a-4*I_KINETIC_Gx3y_Pz;
  abcd[116] = 2.0E0*I_KINETIC_Ix4yz_Pz_a-3*I_KINETIC_Gx2yz_Pz;
  abcd[117] = 2.0E0*I_KINETIC_Ix3y2z_Pz_a-2*I_KINETIC_Gxy2z_Pz;
  abcd[118] = 2.0E0*I_KINETIC_Ix2y3z_Pz_a-1*I_KINETIC_Gx3z_Pz;
  abcd[119] = 2.0E0*I_KINETIC_Ixy4z_Pz_a;
  abcd[120] = 2.0E0*I_KINETIC_I6y_Pz_a-5*I_KINETIC_G4y_Pz;
  abcd[121] = 2.0E0*I_KINETIC_I5yz_Pz_a-4*I_KINETIC_G3yz_Pz;
  abcd[122] = 2.0E0*I_KINETIC_I4y2z_Pz_a-3*I_KINETIC_G2y2z_Pz;
  abcd[123] = 2.0E0*I_KINETIC_I3y3z_Pz_a-2*I_KINETIC_Gy3z_Pz;
  abcd[124] = 2.0E0*I_KINETIC_I2y4z_Pz_a-1*I_KINETIC_G4z_Pz;
  abcd[125] = 2.0E0*I_KINETIC_Iy5z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_KINETIC_H_P_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_KINETIC_I_P_a
   * RHS shell quartet name: SQ_KINETIC_G_P
   ************************************************************/
  abcd[126] = 2.0E0*I_KINETIC_I5xz_Px_a;
  abcd[127] = 2.0E0*I_KINETIC_I4xyz_Px_a;
  abcd[128] = 2.0E0*I_KINETIC_I4x2z_Px_a-1*I_KINETIC_G4x_Px;
  abcd[129] = 2.0E0*I_KINETIC_I3x2yz_Px_a;
  abcd[130] = 2.0E0*I_KINETIC_I3xy2z_Px_a-1*I_KINETIC_G3xy_Px;
  abcd[131] = 2.0E0*I_KINETIC_I3x3z_Px_a-2*I_KINETIC_G3xz_Px;
  abcd[132] = 2.0E0*I_KINETIC_I2x3yz_Px_a;
  abcd[133] = 2.0E0*I_KINETIC_I2x2y2z_Px_a-1*I_KINETIC_G2x2y_Px;
  abcd[134] = 2.0E0*I_KINETIC_I2xy3z_Px_a-2*I_KINETIC_G2xyz_Px;
  abcd[135] = 2.0E0*I_KINETIC_I2x4z_Px_a-3*I_KINETIC_G2x2z_Px;
  abcd[136] = 2.0E0*I_KINETIC_Ix4yz_Px_a;
  abcd[137] = 2.0E0*I_KINETIC_Ix3y2z_Px_a-1*I_KINETIC_Gx3y_Px;
  abcd[138] = 2.0E0*I_KINETIC_Ix2y3z_Px_a-2*I_KINETIC_Gx2yz_Px;
  abcd[139] = 2.0E0*I_KINETIC_Ixy4z_Px_a-3*I_KINETIC_Gxy2z_Px;
  abcd[140] = 2.0E0*I_KINETIC_Ix5z_Px_a-4*I_KINETIC_Gx3z_Px;
  abcd[141] = 2.0E0*I_KINETIC_I5yz_Px_a;
  abcd[142] = 2.0E0*I_KINETIC_I4y2z_Px_a-1*I_KINETIC_G4y_Px;
  abcd[143] = 2.0E0*I_KINETIC_I3y3z_Px_a-2*I_KINETIC_G3yz_Px;
  abcd[144] = 2.0E0*I_KINETIC_I2y4z_Px_a-3*I_KINETIC_G2y2z_Px;
  abcd[145] = 2.0E0*I_KINETIC_Iy5z_Px_a-4*I_KINETIC_Gy3z_Px;
  abcd[146] = 2.0E0*I_KINETIC_I6z_Px_a-5*I_KINETIC_G4z_Px;
  abcd[147] = 2.0E0*I_KINETIC_I5xz_Py_a;
  abcd[148] = 2.0E0*I_KINETIC_I4xyz_Py_a;
  abcd[149] = 2.0E0*I_KINETIC_I4x2z_Py_a-1*I_KINETIC_G4x_Py;
  abcd[150] = 2.0E0*I_KINETIC_I3x2yz_Py_a;
  abcd[151] = 2.0E0*I_KINETIC_I3xy2z_Py_a-1*I_KINETIC_G3xy_Py;
  abcd[152] = 2.0E0*I_KINETIC_I3x3z_Py_a-2*I_KINETIC_G3xz_Py;
  abcd[153] = 2.0E0*I_KINETIC_I2x3yz_Py_a;
  abcd[154] = 2.0E0*I_KINETIC_I2x2y2z_Py_a-1*I_KINETIC_G2x2y_Py;
  abcd[155] = 2.0E0*I_KINETIC_I2xy3z_Py_a-2*I_KINETIC_G2xyz_Py;
  abcd[156] = 2.0E0*I_KINETIC_I2x4z_Py_a-3*I_KINETIC_G2x2z_Py;
  abcd[157] = 2.0E0*I_KINETIC_Ix4yz_Py_a;
  abcd[158] = 2.0E0*I_KINETIC_Ix3y2z_Py_a-1*I_KINETIC_Gx3y_Py;
  abcd[159] = 2.0E0*I_KINETIC_Ix2y3z_Py_a-2*I_KINETIC_Gx2yz_Py;
  abcd[160] = 2.0E0*I_KINETIC_Ixy4z_Py_a-3*I_KINETIC_Gxy2z_Py;
  abcd[161] = 2.0E0*I_KINETIC_Ix5z_Py_a-4*I_KINETIC_Gx3z_Py;
  abcd[162] = 2.0E0*I_KINETIC_I5yz_Py_a;
  abcd[163] = 2.0E0*I_KINETIC_I4y2z_Py_a-1*I_KINETIC_G4y_Py;
  abcd[164] = 2.0E0*I_KINETIC_I3y3z_Py_a-2*I_KINETIC_G3yz_Py;
  abcd[165] = 2.0E0*I_KINETIC_I2y4z_Py_a-3*I_KINETIC_G2y2z_Py;
  abcd[166] = 2.0E0*I_KINETIC_Iy5z_Py_a-4*I_KINETIC_Gy3z_Py;
  abcd[167] = 2.0E0*I_KINETIC_I6z_Py_a-5*I_KINETIC_G4z_Py;
  abcd[168] = 2.0E0*I_KINETIC_I5xz_Pz_a;
  abcd[169] = 2.0E0*I_KINETIC_I4xyz_Pz_a;
  abcd[170] = 2.0E0*I_KINETIC_I4x2z_Pz_a-1*I_KINETIC_G4x_Pz;
  abcd[171] = 2.0E0*I_KINETIC_I3x2yz_Pz_a;
  abcd[172] = 2.0E0*I_KINETIC_I3xy2z_Pz_a-1*I_KINETIC_G3xy_Pz;
  abcd[173] = 2.0E0*I_KINETIC_I3x3z_Pz_a-2*I_KINETIC_G3xz_Pz;
  abcd[174] = 2.0E0*I_KINETIC_I2x3yz_Pz_a;
  abcd[175] = 2.0E0*I_KINETIC_I2x2y2z_Pz_a-1*I_KINETIC_G2x2y_Pz;
  abcd[176] = 2.0E0*I_KINETIC_I2xy3z_Pz_a-2*I_KINETIC_G2xyz_Pz;
  abcd[177] = 2.0E0*I_KINETIC_I2x4z_Pz_a-3*I_KINETIC_G2x2z_Pz;
  abcd[178] = 2.0E0*I_KINETIC_Ix4yz_Pz_a;
  abcd[179] = 2.0E0*I_KINETIC_Ix3y2z_Pz_a-1*I_KINETIC_Gx3y_Pz;
  abcd[180] = 2.0E0*I_KINETIC_Ix2y3z_Pz_a-2*I_KINETIC_Gx2yz_Pz;
  abcd[181] = 2.0E0*I_KINETIC_Ixy4z_Pz_a-3*I_KINETIC_Gxy2z_Pz;
  abcd[182] = 2.0E0*I_KINETIC_Ix5z_Pz_a-4*I_KINETIC_Gx3z_Pz;
  abcd[183] = 2.0E0*I_KINETIC_I5yz_Pz_a;
  abcd[184] = 2.0E0*I_KINETIC_I4y2z_Pz_a-1*I_KINETIC_G4y_Pz;
  abcd[185] = 2.0E0*I_KINETIC_I3y3z_Pz_a-2*I_KINETIC_G3yz_Pz;
  abcd[186] = 2.0E0*I_KINETIC_I2y4z_Pz_a-3*I_KINETIC_G2y2z_Pz;
  abcd[187] = 2.0E0*I_KINETIC_Iy5z_Pz_a-4*I_KINETIC_Gy3z_Pz;
  abcd[188] = 2.0E0*I_KINETIC_I6z_Pz_a-5*I_KINETIC_G4z_Pz;
}
