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
// BRA1 as redundant position, total RHS integrals evaluated as: 0
// BRA2 as redundant position, total RHS integrals evaluated as: 0
// KET1 as redundant position, total RHS integrals evaluated as: 0
// KET2 as redundant position, total RHS integrals evaluated as: 0
// the redundant position is: NOT AVIALABLE
//

//
// @@@@ derivative position-direction information
// BRA1
// X
// Y
// Z
// ####

void hgp_os_twobodyoverlap_g_g_d1(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_TWOBODYOVERLAP_M9x_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M8xy_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M8xz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M7x2y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M7xyz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M7x2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M6x3y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M6x2yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M6xy2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M6x3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M5x4y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M5x3yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M5x2y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M5xy3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M5x4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M4x5y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M4x4yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M4x3y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M4x2y3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M4xy4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M4x5z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M3x6y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M3x5yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M3x4y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M3x3y3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M3x2y4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M3xy5z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M3x6z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x7y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x6yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x5y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x4y3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x3y4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x2y5z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M2xy6z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x7z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx8y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx7yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx6y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx5y3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx4y4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx3y5z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx2y6z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Mxy7z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx8z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M9y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M8yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M7y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M6y3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M5y4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M4y5z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M3y6z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M2y7z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_My8z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_M9z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L8x_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L7xy_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L7xz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L6x2y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L6xyz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L6x2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L5x3y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L5x2yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L5xy2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L5x3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L4x4y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L4x3yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L4x2y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L4xy3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L4x4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x5y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x4yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x3y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x2y3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L3xy4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x5z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x6y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x5yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x4y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x3y3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x2y4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L2xy5z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x6z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx7y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx6yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx5y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx4y3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx3y4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx2y5z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Lxy6z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx7z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L8y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L7yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L6y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L5y3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L4y4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L3y5z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L2y6z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ly7z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_L8z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K7x_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K6xy_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K6xz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K5x2y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K5xyz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K5x2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x3y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x2yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K4xy2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x4y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x3yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x2y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K3xy3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x5y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x4yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x3y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x2y3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K2xy4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x5z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx6y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx5yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx4y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx3y3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx2y4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Kxy5z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx6z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K7y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K6yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K5y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K4y3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K3y4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K2y5z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ky6z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K7z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I6x_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xy_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I4xyz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x2yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3xy2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x3yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x2y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2xy3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix4yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix3y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix2y3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ixy4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I6y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I5yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I4y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I3y3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I2y4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Iy5z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_I6z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H5x_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H4xy_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H4xz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H3x2y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H3xyz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H3x2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x3y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x2yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H2xy2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H2x3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx4y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx3yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx2y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hxy3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hx4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H5y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H4yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H3y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H2y3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Hy4z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_H5z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_K7x_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K6xy_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K6xz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K5x2y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K5xyz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K5x2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x3y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x2yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K4xy2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x4y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x3yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x2y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K3xy3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x5y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x4yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x3y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x2y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K2xy4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx6y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx5yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx4y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx3y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx2y4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Kxy5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx6z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K7y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K6yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K5y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K4y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K3y4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K2y5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Ky6z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_K7z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I6x_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xy_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I4xyz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x2yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I3xy2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x3yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x2y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I2xy3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix4yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix3y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix2y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Ixy4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I6y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I5yz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I4y2z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I3y3z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I2y4z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Iy5z_S = 0.0E0;
  Double I_TWOBODYOVERLAP_I6z_S = 0.0E0;
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
    Double zeta  = 1.0E0/onedz;
    Double zdiff = iexpdiff[ip2];
    Double alpha = 0.5E0*(zeta+zdiff);
    Double beta  = 0.5E0*(zeta-zdiff);
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
     * shell quartet name: SQ_TWOBODYOVERLAP_K_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_K7x_S_vrr = PAX*I_TWOBODYOVERLAP_I6x_S_vrr+6*oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_K6xy_S_vrr = PAY*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_K6xz_S_vrr = PAZ*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_K5x2y_S_vrr = PAY*I_TWOBODYOVERLAP_I5xy_S_vrr+oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_K5xyz_S_vrr = PAZ*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_TWOBODYOVERLAP_K5x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_I5xz_S_vrr+oned2z*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_K4x3y_S_vrr = PAY*I_TWOBODYOVERLAP_I4x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_K4x2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_TWOBODYOVERLAP_K4xy2z_S_vrr = PAY*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    Double I_TWOBODYOVERLAP_K4x3z_S_vrr = PAZ*I_TWOBODYOVERLAP_I4x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_TWOBODYOVERLAP_K3x4y_S_vrr = PAX*I_TWOBODYOVERLAP_I2x4y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_K3x3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    Double I_TWOBODYOVERLAP_K3x2y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_I3x2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_K3xy3z_S_vrr = PAY*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    Double I_TWOBODYOVERLAP_K3x4z_S_vrr = PAX*I_TWOBODYOVERLAP_I2x4z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_K2x5y_S_vrr = PAX*I_TWOBODYOVERLAP_Ix5y_S_vrr+oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_K2x4yz_S_vrr = PAZ*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    Double I_TWOBODYOVERLAP_K2x3y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_I2x3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_K2x2y3z_S_vrr = PAY*I_TWOBODYOVERLAP_I2xy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_K2xy4z_S_vrr = PAY*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    Double I_TWOBODYOVERLAP_K2x5z_S_vrr = PAX*I_TWOBODYOVERLAP_Ix5z_S_vrr+oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx6y_S_vrr = PAX*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_Kx5yz_S_vrr = PAZ*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_TWOBODYOVERLAP_Kx4y2z_S_vrr = PAX*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx3y3z_S_vrr = PAX*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx2y4z_S_vrr = PAX*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    Double I_TWOBODYOVERLAP_Kxy5z_S_vrr = PAY*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_TWOBODYOVERLAP_Kx6z_S_vrr = PAX*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_K7y_S_vrr = PAY*I_TWOBODYOVERLAP_I6y_S_vrr+6*oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_K6yz_S_vrr = PAZ*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_K5y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_I5yz_S_vrr+oned2z*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_K4y3z_S_vrr = PAZ*I_TWOBODYOVERLAP_I4y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_TWOBODYOVERLAP_K3y4z_S_vrr = PAY*I_TWOBODYOVERLAP_I2y4z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_TWOBODYOVERLAP_K2y5z_S_vrr = PAY*I_TWOBODYOVERLAP_Iy5z_S_vrr+oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_Ky6z_S_vrr = PAY*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_K7z_S_vrr = PAZ*I_TWOBODYOVERLAP_I6z_S_vrr+6*oned2z*I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_L_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_L8x_S_vrr = PAX*I_TWOBODYOVERLAP_K7x_S_vrr+7*oned2z*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_L7xy_S_vrr = PAY*I_TWOBODYOVERLAP_K7x_S_vrr;
    Double I_TWOBODYOVERLAP_L7xz_S_vrr = PAZ*I_TWOBODYOVERLAP_K7x_S_vrr;
    Double I_TWOBODYOVERLAP_L6x2y_S_vrr = PAY*I_TWOBODYOVERLAP_K6xy_S_vrr+oned2z*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_L6xyz_S_vrr = PAZ*I_TWOBODYOVERLAP_K6xy_S_vrr;
    Double I_TWOBODYOVERLAP_L6x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_K6xz_S_vrr+oned2z*I_TWOBODYOVERLAP_I6x_S_vrr;
    Double I_TWOBODYOVERLAP_L5x3y_S_vrr = PAY*I_TWOBODYOVERLAP_K5x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I5xy_S_vrr;
    Double I_TWOBODYOVERLAP_L5x2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_K5x2y_S_vrr;
    Double I_TWOBODYOVERLAP_L5xy2z_S_vrr = PAY*I_TWOBODYOVERLAP_K5x2z_S_vrr;
    Double I_TWOBODYOVERLAP_L5x3z_S_vrr = PAZ*I_TWOBODYOVERLAP_K5x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I5xz_S_vrr;
    Double I_TWOBODYOVERLAP_L4x4y_S_vrr = PAY*I_TWOBODYOVERLAP_K4x3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_TWOBODYOVERLAP_L4x3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_K4x3y_S_vrr;
    Double I_TWOBODYOVERLAP_L4x2y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_K4x2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    Double I_TWOBODYOVERLAP_L4xy3z_S_vrr = PAY*I_TWOBODYOVERLAP_K4x3z_S_vrr;
    Double I_TWOBODYOVERLAP_L4x4z_S_vrr = PAZ*I_TWOBODYOVERLAP_K4x3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    Double I_TWOBODYOVERLAP_L3x5y_S_vrr = PAX*I_TWOBODYOVERLAP_K2x5y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    Double I_TWOBODYOVERLAP_L3x4yz_S_vrr = PAZ*I_TWOBODYOVERLAP_K3x4y_S_vrr;
    Double I_TWOBODYOVERLAP_L3x3y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_K3x3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    Double I_TWOBODYOVERLAP_L3x2y3z_S_vrr = PAY*I_TWOBODYOVERLAP_K3xy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    Double I_TWOBODYOVERLAP_L3xy4z_S_vrr = PAY*I_TWOBODYOVERLAP_K3x4z_S_vrr;
    Double I_TWOBODYOVERLAP_L3x5z_S_vrr = PAX*I_TWOBODYOVERLAP_K2x5z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    Double I_TWOBODYOVERLAP_L2x6y_S_vrr = PAX*I_TWOBODYOVERLAP_Kx6y_S_vrr+oned2z*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_L2x5yz_S_vrr = PAZ*I_TWOBODYOVERLAP_K2x5y_S_vrr;
    Double I_TWOBODYOVERLAP_L2x4y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_K2x4yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    Double I_TWOBODYOVERLAP_L2x3y3z_S_vrr = PAX*I_TWOBODYOVERLAP_Kx3y3z_S_vrr+oned2z*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    Double I_TWOBODYOVERLAP_L2x2y4z_S_vrr = PAY*I_TWOBODYOVERLAP_K2xy4z_S_vrr+oned2z*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    Double I_TWOBODYOVERLAP_L2xy5z_S_vrr = PAY*I_TWOBODYOVERLAP_K2x5z_S_vrr;
    Double I_TWOBODYOVERLAP_L2x6z_S_vrr = PAX*I_TWOBODYOVERLAP_Kx6z_S_vrr+oned2z*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_Lx7y_S_vrr = PAX*I_TWOBODYOVERLAP_K7y_S_vrr;
    Double I_TWOBODYOVERLAP_Lx6yz_S_vrr = PAZ*I_TWOBODYOVERLAP_Kx6y_S_vrr;
    Double I_TWOBODYOVERLAP_Lx5y2z_S_vrr = PAX*I_TWOBODYOVERLAP_K5y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Lx4y3z_S_vrr = PAX*I_TWOBODYOVERLAP_K4y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Lx3y4z_S_vrr = PAX*I_TWOBODYOVERLAP_K3y4z_S_vrr;
    Double I_TWOBODYOVERLAP_Lx2y5z_S_vrr = PAX*I_TWOBODYOVERLAP_K2y5z_S_vrr;
    Double I_TWOBODYOVERLAP_Lxy6z_S_vrr = PAY*I_TWOBODYOVERLAP_Kx6z_S_vrr;
    Double I_TWOBODYOVERLAP_Lx7z_S_vrr = PAX*I_TWOBODYOVERLAP_K7z_S_vrr;
    Double I_TWOBODYOVERLAP_L8y_S_vrr = PAY*I_TWOBODYOVERLAP_K7y_S_vrr+7*oned2z*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_L7yz_S_vrr = PAZ*I_TWOBODYOVERLAP_K7y_S_vrr;
    Double I_TWOBODYOVERLAP_L6y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_K6yz_S_vrr+oned2z*I_TWOBODYOVERLAP_I6y_S_vrr;
    Double I_TWOBODYOVERLAP_L5y3z_S_vrr = PAZ*I_TWOBODYOVERLAP_K5y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_I5yz_S_vrr;
    Double I_TWOBODYOVERLAP_L4y4z_S_vrr = PAZ*I_TWOBODYOVERLAP_K4y3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    Double I_TWOBODYOVERLAP_L3y5z_S_vrr = PAY*I_TWOBODYOVERLAP_K2y5z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    Double I_TWOBODYOVERLAP_L2y6z_S_vrr = PAY*I_TWOBODYOVERLAP_Ky6z_S_vrr+oned2z*I_TWOBODYOVERLAP_I6z_S_vrr;
    Double I_TWOBODYOVERLAP_Ly7z_S_vrr = PAY*I_TWOBODYOVERLAP_K7z_S_vrr;
    Double I_TWOBODYOVERLAP_L8z_S_vrr = PAZ*I_TWOBODYOVERLAP_K7z_S_vrr+7*oned2z*I_TWOBODYOVERLAP_I6z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_M_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_L_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_M9x_S_vrr = PAX*I_TWOBODYOVERLAP_L8x_S_vrr+8*oned2z*I_TWOBODYOVERLAP_K7x_S_vrr;
    Double I_TWOBODYOVERLAP_M8xy_S_vrr = PAY*I_TWOBODYOVERLAP_L8x_S_vrr;
    Double I_TWOBODYOVERLAP_M8xz_S_vrr = PAZ*I_TWOBODYOVERLAP_L8x_S_vrr;
    Double I_TWOBODYOVERLAP_M7x2y_S_vrr = PAY*I_TWOBODYOVERLAP_L7xy_S_vrr+oned2z*I_TWOBODYOVERLAP_K7x_S_vrr;
    Double I_TWOBODYOVERLAP_M7xyz_S_vrr = PAZ*I_TWOBODYOVERLAP_L7xy_S_vrr;
    Double I_TWOBODYOVERLAP_M7x2z_S_vrr = PAZ*I_TWOBODYOVERLAP_L7xz_S_vrr+oned2z*I_TWOBODYOVERLAP_K7x_S_vrr;
    Double I_TWOBODYOVERLAP_M6x3y_S_vrr = PAY*I_TWOBODYOVERLAP_L6x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_K6xy_S_vrr;
    Double I_TWOBODYOVERLAP_M6x2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_L6x2y_S_vrr;
    Double I_TWOBODYOVERLAP_M6xy2z_S_vrr = PAY*I_TWOBODYOVERLAP_L6x2z_S_vrr;
    Double I_TWOBODYOVERLAP_M6x3z_S_vrr = PAZ*I_TWOBODYOVERLAP_L6x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_K6xz_S_vrr;
    Double I_TWOBODYOVERLAP_M5x4y_S_vrr = PAY*I_TWOBODYOVERLAP_L5x3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_K5x2y_S_vrr;
    Double I_TWOBODYOVERLAP_M5x3yz_S_vrr = PAZ*I_TWOBODYOVERLAP_L5x3y_S_vrr;
    Double I_TWOBODYOVERLAP_M5x2y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_L5x2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_K5x2y_S_vrr;
    Double I_TWOBODYOVERLAP_M5xy3z_S_vrr = PAY*I_TWOBODYOVERLAP_L5x3z_S_vrr;
    Double I_TWOBODYOVERLAP_M5x4z_S_vrr = PAZ*I_TWOBODYOVERLAP_L5x3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_K5x2z_S_vrr;
    Double I_TWOBODYOVERLAP_M4x5y_S_vrr = PAX*I_TWOBODYOVERLAP_L3x5y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_K2x5y_S_vrr;
    Double I_TWOBODYOVERLAP_M4x4yz_S_vrr = PAZ*I_TWOBODYOVERLAP_L4x4y_S_vrr;
    Double I_TWOBODYOVERLAP_M4x3y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_L4x3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_K4x3y_S_vrr;
    Double I_TWOBODYOVERLAP_M4x2y3z_S_vrr = PAY*I_TWOBODYOVERLAP_L4xy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_K4x3z_S_vrr;
    Double I_TWOBODYOVERLAP_M4xy4z_S_vrr = PAY*I_TWOBODYOVERLAP_L4x4z_S_vrr;
    Double I_TWOBODYOVERLAP_M4x5z_S_vrr = PAX*I_TWOBODYOVERLAP_L3x5z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_K2x5z_S_vrr;
    Double I_TWOBODYOVERLAP_M3x6y_S_vrr = PAX*I_TWOBODYOVERLAP_L2x6y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Kx6y_S_vrr;
    Double I_TWOBODYOVERLAP_M3x5yz_S_vrr = PAZ*I_TWOBODYOVERLAP_L3x5y_S_vrr;
    Double I_TWOBODYOVERLAP_M3x4y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_L3x4yz_S_vrr+oned2z*I_TWOBODYOVERLAP_K3x4y_S_vrr;
    Double I_TWOBODYOVERLAP_M3x3y3z_S_vrr = PAZ*I_TWOBODYOVERLAP_L3x3y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_K3x3yz_S_vrr;
    Double I_TWOBODYOVERLAP_M3x2y4z_S_vrr = PAY*I_TWOBODYOVERLAP_L3xy4z_S_vrr+oned2z*I_TWOBODYOVERLAP_K3x4z_S_vrr;
    Double I_TWOBODYOVERLAP_M3xy5z_S_vrr = PAY*I_TWOBODYOVERLAP_L3x5z_S_vrr;
    Double I_TWOBODYOVERLAP_M3x6z_S_vrr = PAX*I_TWOBODYOVERLAP_L2x6z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Kx6z_S_vrr;
    Double I_TWOBODYOVERLAP_M2x7y_S_vrr = PAX*I_TWOBODYOVERLAP_Lx7y_S_vrr+oned2z*I_TWOBODYOVERLAP_K7y_S_vrr;
    Double I_TWOBODYOVERLAP_M2x6yz_S_vrr = PAZ*I_TWOBODYOVERLAP_L2x6y_S_vrr;
    Double I_TWOBODYOVERLAP_M2x5y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_L2x5yz_S_vrr+oned2z*I_TWOBODYOVERLAP_K2x5y_S_vrr;
    Double I_TWOBODYOVERLAP_M2x4y3z_S_vrr = PAX*I_TWOBODYOVERLAP_Lx4y3z_S_vrr+oned2z*I_TWOBODYOVERLAP_K4y3z_S_vrr;
    Double I_TWOBODYOVERLAP_M2x3y4z_S_vrr = PAX*I_TWOBODYOVERLAP_Lx3y4z_S_vrr+oned2z*I_TWOBODYOVERLAP_K3y4z_S_vrr;
    Double I_TWOBODYOVERLAP_M2x2y5z_S_vrr = PAY*I_TWOBODYOVERLAP_L2xy5z_S_vrr+oned2z*I_TWOBODYOVERLAP_K2x5z_S_vrr;
    Double I_TWOBODYOVERLAP_M2xy6z_S_vrr = PAY*I_TWOBODYOVERLAP_L2x6z_S_vrr;
    Double I_TWOBODYOVERLAP_M2x7z_S_vrr = PAX*I_TWOBODYOVERLAP_Lx7z_S_vrr+oned2z*I_TWOBODYOVERLAP_K7z_S_vrr;
    Double I_TWOBODYOVERLAP_Mx8y_S_vrr = PAX*I_TWOBODYOVERLAP_L8y_S_vrr;
    Double I_TWOBODYOVERLAP_Mx7yz_S_vrr = PAZ*I_TWOBODYOVERLAP_Lx7y_S_vrr;
    Double I_TWOBODYOVERLAP_Mx6y2z_S_vrr = PAX*I_TWOBODYOVERLAP_L6y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Mx5y3z_S_vrr = PAX*I_TWOBODYOVERLAP_L5y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Mx4y4z_S_vrr = PAX*I_TWOBODYOVERLAP_L4y4z_S_vrr;
    Double I_TWOBODYOVERLAP_Mx3y5z_S_vrr = PAX*I_TWOBODYOVERLAP_L3y5z_S_vrr;
    Double I_TWOBODYOVERLAP_Mx2y6z_S_vrr = PAX*I_TWOBODYOVERLAP_L2y6z_S_vrr;
    Double I_TWOBODYOVERLAP_Mxy7z_S_vrr = PAY*I_TWOBODYOVERLAP_Lx7z_S_vrr;
    Double I_TWOBODYOVERLAP_Mx8z_S_vrr = PAX*I_TWOBODYOVERLAP_L8z_S_vrr;
    Double I_TWOBODYOVERLAP_M9y_S_vrr = PAY*I_TWOBODYOVERLAP_L8y_S_vrr+8*oned2z*I_TWOBODYOVERLAP_K7y_S_vrr;
    Double I_TWOBODYOVERLAP_M8yz_S_vrr = PAZ*I_TWOBODYOVERLAP_L8y_S_vrr;
    Double I_TWOBODYOVERLAP_M7y2z_S_vrr = PAZ*I_TWOBODYOVERLAP_L7yz_S_vrr+oned2z*I_TWOBODYOVERLAP_K7y_S_vrr;
    Double I_TWOBODYOVERLAP_M6y3z_S_vrr = PAZ*I_TWOBODYOVERLAP_L6y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_K6yz_S_vrr;
    Double I_TWOBODYOVERLAP_M5y4z_S_vrr = PAZ*I_TWOBODYOVERLAP_L5y3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_K5y2z_S_vrr;
    Double I_TWOBODYOVERLAP_M4y5z_S_vrr = PAY*I_TWOBODYOVERLAP_L3y5z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_K2y5z_S_vrr;
    Double I_TWOBODYOVERLAP_M3y6z_S_vrr = PAY*I_TWOBODYOVERLAP_L2y6z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Ky6z_S_vrr;
    Double I_TWOBODYOVERLAP_M2y7z_S_vrr = PAY*I_TWOBODYOVERLAP_Ly7z_S_vrr+oned2z*I_TWOBODYOVERLAP_K7z_S_vrr;
    Double I_TWOBODYOVERLAP_My8z_S_vrr = PAY*I_TWOBODYOVERLAP_L8z_S_vrr;
    Double I_TWOBODYOVERLAP_M9z_S_vrr = PAZ*I_TWOBODYOVERLAP_L8z_S_vrr+8*oned2z*I_TWOBODYOVERLAP_K7z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_M_S_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_M_S_a_coefs = alpha;
    I_TWOBODYOVERLAP_M9x_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M9x_S_vrr;
    I_TWOBODYOVERLAP_M8xy_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M8xy_S_vrr;
    I_TWOBODYOVERLAP_M8xz_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M8xz_S_vrr;
    I_TWOBODYOVERLAP_M7x2y_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M7x2y_S_vrr;
    I_TWOBODYOVERLAP_M7xyz_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M7xyz_S_vrr;
    I_TWOBODYOVERLAP_M7x2z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M7x2z_S_vrr;
    I_TWOBODYOVERLAP_M6x3y_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M6x3y_S_vrr;
    I_TWOBODYOVERLAP_M6x2yz_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M6x2yz_S_vrr;
    I_TWOBODYOVERLAP_M6xy2z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M6xy2z_S_vrr;
    I_TWOBODYOVERLAP_M6x3z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M6x3z_S_vrr;
    I_TWOBODYOVERLAP_M5x4y_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M5x4y_S_vrr;
    I_TWOBODYOVERLAP_M5x3yz_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M5x3yz_S_vrr;
    I_TWOBODYOVERLAP_M5x2y2z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M5x2y2z_S_vrr;
    I_TWOBODYOVERLAP_M5xy3z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M5xy3z_S_vrr;
    I_TWOBODYOVERLAP_M5x4z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M5x4z_S_vrr;
    I_TWOBODYOVERLAP_M4x5y_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M4x5y_S_vrr;
    I_TWOBODYOVERLAP_M4x4yz_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M4x4yz_S_vrr;
    I_TWOBODYOVERLAP_M4x3y2z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M4x3y2z_S_vrr;
    I_TWOBODYOVERLAP_M4x2y3z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M4x2y3z_S_vrr;
    I_TWOBODYOVERLAP_M4xy4z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M4xy4z_S_vrr;
    I_TWOBODYOVERLAP_M4x5z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M4x5z_S_vrr;
    I_TWOBODYOVERLAP_M3x6y_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M3x6y_S_vrr;
    I_TWOBODYOVERLAP_M3x5yz_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M3x5yz_S_vrr;
    I_TWOBODYOVERLAP_M3x4y2z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M3x4y2z_S_vrr;
    I_TWOBODYOVERLAP_M3x3y3z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M3x3y3z_S_vrr;
    I_TWOBODYOVERLAP_M3x2y4z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M3x2y4z_S_vrr;
    I_TWOBODYOVERLAP_M3xy5z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M3xy5z_S_vrr;
    I_TWOBODYOVERLAP_M3x6z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M3x6z_S_vrr;
    I_TWOBODYOVERLAP_M2x7y_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M2x7y_S_vrr;
    I_TWOBODYOVERLAP_M2x6yz_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M2x6yz_S_vrr;
    I_TWOBODYOVERLAP_M2x5y2z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M2x5y2z_S_vrr;
    I_TWOBODYOVERLAP_M2x4y3z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M2x4y3z_S_vrr;
    I_TWOBODYOVERLAP_M2x3y4z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M2x3y4z_S_vrr;
    I_TWOBODYOVERLAP_M2x2y5z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M2x2y5z_S_vrr;
    I_TWOBODYOVERLAP_M2xy6z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M2xy6z_S_vrr;
    I_TWOBODYOVERLAP_M2x7z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M2x7z_S_vrr;
    I_TWOBODYOVERLAP_Mx8y_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_Mx8y_S_vrr;
    I_TWOBODYOVERLAP_Mx7yz_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_Mx7yz_S_vrr;
    I_TWOBODYOVERLAP_Mx6y2z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_Mx6y2z_S_vrr;
    I_TWOBODYOVERLAP_Mx5y3z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_Mx5y3z_S_vrr;
    I_TWOBODYOVERLAP_Mx4y4z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_Mx4y4z_S_vrr;
    I_TWOBODYOVERLAP_Mx3y5z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_Mx3y5z_S_vrr;
    I_TWOBODYOVERLAP_Mx2y6z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_Mx2y6z_S_vrr;
    I_TWOBODYOVERLAP_Mxy7z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_Mxy7z_S_vrr;
    I_TWOBODYOVERLAP_Mx8z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_Mx8z_S_vrr;
    I_TWOBODYOVERLAP_M9y_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M9y_S_vrr;
    I_TWOBODYOVERLAP_M8yz_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M8yz_S_vrr;
    I_TWOBODYOVERLAP_M7y2z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M7y2z_S_vrr;
    I_TWOBODYOVERLAP_M6y3z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M6y3z_S_vrr;
    I_TWOBODYOVERLAP_M5y4z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M5y4z_S_vrr;
    I_TWOBODYOVERLAP_M4y5z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M4y5z_S_vrr;
    I_TWOBODYOVERLAP_M3y6z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M3y6z_S_vrr;
    I_TWOBODYOVERLAP_M2y7z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M2y7z_S_vrr;
    I_TWOBODYOVERLAP_My8z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_My8z_S_vrr;
    I_TWOBODYOVERLAP_M9z_S_a += SQ_TWOBODYOVERLAP_M_S_a_coefs*I_TWOBODYOVERLAP_M9z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_L_S_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_L_S_a_coefs = alpha;
    I_TWOBODYOVERLAP_L8x_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L8x_S_vrr;
    I_TWOBODYOVERLAP_L7xy_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L7xy_S_vrr;
    I_TWOBODYOVERLAP_L7xz_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L7xz_S_vrr;
    I_TWOBODYOVERLAP_L6x2y_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L6x2y_S_vrr;
    I_TWOBODYOVERLAP_L6xyz_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L6xyz_S_vrr;
    I_TWOBODYOVERLAP_L6x2z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L6x2z_S_vrr;
    I_TWOBODYOVERLAP_L5x3y_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L5x3y_S_vrr;
    I_TWOBODYOVERLAP_L5x2yz_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L5x2yz_S_vrr;
    I_TWOBODYOVERLAP_L5xy2z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L5xy2z_S_vrr;
    I_TWOBODYOVERLAP_L5x3z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L5x3z_S_vrr;
    I_TWOBODYOVERLAP_L4x4y_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L4x4y_S_vrr;
    I_TWOBODYOVERLAP_L4x3yz_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L4x3yz_S_vrr;
    I_TWOBODYOVERLAP_L4x2y2z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L4x2y2z_S_vrr;
    I_TWOBODYOVERLAP_L4xy3z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L4xy3z_S_vrr;
    I_TWOBODYOVERLAP_L4x4z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L4x4z_S_vrr;
    I_TWOBODYOVERLAP_L3x5y_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L3x5y_S_vrr;
    I_TWOBODYOVERLAP_L3x4yz_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L3x4yz_S_vrr;
    I_TWOBODYOVERLAP_L3x3y2z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L3x3y2z_S_vrr;
    I_TWOBODYOVERLAP_L3x2y3z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L3x2y3z_S_vrr;
    I_TWOBODYOVERLAP_L3xy4z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L3xy4z_S_vrr;
    I_TWOBODYOVERLAP_L3x5z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L3x5z_S_vrr;
    I_TWOBODYOVERLAP_L2x6y_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L2x6y_S_vrr;
    I_TWOBODYOVERLAP_L2x5yz_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L2x5yz_S_vrr;
    I_TWOBODYOVERLAP_L2x4y2z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L2x4y2z_S_vrr;
    I_TWOBODYOVERLAP_L2x3y3z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L2x3y3z_S_vrr;
    I_TWOBODYOVERLAP_L2x2y4z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L2x2y4z_S_vrr;
    I_TWOBODYOVERLAP_L2xy5z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L2xy5z_S_vrr;
    I_TWOBODYOVERLAP_L2x6z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L2x6z_S_vrr;
    I_TWOBODYOVERLAP_Lx7y_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_Lx7y_S_vrr;
    I_TWOBODYOVERLAP_Lx6yz_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_Lx6yz_S_vrr;
    I_TWOBODYOVERLAP_Lx5y2z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_Lx5y2z_S_vrr;
    I_TWOBODYOVERLAP_Lx4y3z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_Lx4y3z_S_vrr;
    I_TWOBODYOVERLAP_Lx3y4z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_Lx3y4z_S_vrr;
    I_TWOBODYOVERLAP_Lx2y5z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_Lx2y5z_S_vrr;
    I_TWOBODYOVERLAP_Lxy6z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_Lxy6z_S_vrr;
    I_TWOBODYOVERLAP_Lx7z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_Lx7z_S_vrr;
    I_TWOBODYOVERLAP_L8y_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L8y_S_vrr;
    I_TWOBODYOVERLAP_L7yz_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L7yz_S_vrr;
    I_TWOBODYOVERLAP_L6y2z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L6y2z_S_vrr;
    I_TWOBODYOVERLAP_L5y3z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L5y3z_S_vrr;
    I_TWOBODYOVERLAP_L4y4z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L4y4z_S_vrr;
    I_TWOBODYOVERLAP_L3y5z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L3y5z_S_vrr;
    I_TWOBODYOVERLAP_L2y6z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L2y6z_S_vrr;
    I_TWOBODYOVERLAP_Ly7z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_Ly7z_S_vrr;
    I_TWOBODYOVERLAP_L8z_S_a += SQ_TWOBODYOVERLAP_L_S_a_coefs*I_TWOBODYOVERLAP_L8z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_K_S_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_K_S_a_coefs = alpha;
    I_TWOBODYOVERLAP_K7x_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K7x_S_vrr;
    I_TWOBODYOVERLAP_K6xy_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K6xy_S_vrr;
    I_TWOBODYOVERLAP_K6xz_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K6xz_S_vrr;
    I_TWOBODYOVERLAP_K5x2y_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K5x2y_S_vrr;
    I_TWOBODYOVERLAP_K5xyz_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K5xyz_S_vrr;
    I_TWOBODYOVERLAP_K5x2z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K5x2z_S_vrr;
    I_TWOBODYOVERLAP_K4x3y_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K4x3y_S_vrr;
    I_TWOBODYOVERLAP_K4x2yz_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K4x2yz_S_vrr;
    I_TWOBODYOVERLAP_K4xy2z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K4xy2z_S_vrr;
    I_TWOBODYOVERLAP_K4x3z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K4x3z_S_vrr;
    I_TWOBODYOVERLAP_K3x4y_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K3x4y_S_vrr;
    I_TWOBODYOVERLAP_K3x3yz_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K3x3yz_S_vrr;
    I_TWOBODYOVERLAP_K3x2y2z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K3x2y2z_S_vrr;
    I_TWOBODYOVERLAP_K3xy3z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K3xy3z_S_vrr;
    I_TWOBODYOVERLAP_K3x4z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K3x4z_S_vrr;
    I_TWOBODYOVERLAP_K2x5y_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K2x5y_S_vrr;
    I_TWOBODYOVERLAP_K2x4yz_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K2x4yz_S_vrr;
    I_TWOBODYOVERLAP_K2x3y2z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K2x3y2z_S_vrr;
    I_TWOBODYOVERLAP_K2x2y3z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K2x2y3z_S_vrr;
    I_TWOBODYOVERLAP_K2xy4z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K2xy4z_S_vrr;
    I_TWOBODYOVERLAP_K2x5z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K2x5z_S_vrr;
    I_TWOBODYOVERLAP_Kx6y_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_Kx6y_S_vrr;
    I_TWOBODYOVERLAP_Kx5yz_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_Kx5yz_S_vrr;
    I_TWOBODYOVERLAP_Kx4y2z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_Kx4y2z_S_vrr;
    I_TWOBODYOVERLAP_Kx3y3z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_Kx3y3z_S_vrr;
    I_TWOBODYOVERLAP_Kx2y4z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_Kx2y4z_S_vrr;
    I_TWOBODYOVERLAP_Kxy5z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_Kxy5z_S_vrr;
    I_TWOBODYOVERLAP_Kx6z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_Kx6z_S_vrr;
    I_TWOBODYOVERLAP_K7y_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K7y_S_vrr;
    I_TWOBODYOVERLAP_K6yz_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K6yz_S_vrr;
    I_TWOBODYOVERLAP_K5y2z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K5y2z_S_vrr;
    I_TWOBODYOVERLAP_K4y3z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K4y3z_S_vrr;
    I_TWOBODYOVERLAP_K3y4z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K3y4z_S_vrr;
    I_TWOBODYOVERLAP_K2y5z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K2y5z_S_vrr;
    I_TWOBODYOVERLAP_Ky6z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_Ky6z_S_vrr;
    I_TWOBODYOVERLAP_K7z_S_a += SQ_TWOBODYOVERLAP_K_S_a_coefs*I_TWOBODYOVERLAP_K7z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_I_S_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_I_S_a_coefs = alpha;
    I_TWOBODYOVERLAP_I6x_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I6x_S_vrr;
    I_TWOBODYOVERLAP_I5xy_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I5xy_S_vrr;
    I_TWOBODYOVERLAP_I5xz_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I5xz_S_vrr;
    I_TWOBODYOVERLAP_I4x2y_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    I_TWOBODYOVERLAP_I4xyz_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I4xyz_S_vrr;
    I_TWOBODYOVERLAP_I4x2z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3y_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    I_TWOBODYOVERLAP_I3x2yz_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    I_TWOBODYOVERLAP_I3xy2z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I3xy2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4y_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    I_TWOBODYOVERLAP_I2x3yz_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    I_TWOBODYOVERLAP_I2x2y2z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I2x2y2z_S_vrr;
    I_TWOBODYOVERLAP_I2xy3z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5y_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    I_TWOBODYOVERLAP_Ix4yz_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_Ix4yz_S_vrr;
    I_TWOBODYOVERLAP_Ix3y2z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_Ix3y2z_S_vrr;
    I_TWOBODYOVERLAP_Ix2y3z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_Ix2y3z_S_vrr;
    I_TWOBODYOVERLAP_Ixy4z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_Ixy4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    I_TWOBODYOVERLAP_I6y_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I6y_S_vrr;
    I_TWOBODYOVERLAP_I5yz_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I5yz_S_vrr;
    I_TWOBODYOVERLAP_I4y2z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    I_TWOBODYOVERLAP_I3y3z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    I_TWOBODYOVERLAP_I2y4z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    I_TWOBODYOVERLAP_Iy5z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    I_TWOBODYOVERLAP_I6z_S_a += SQ_TWOBODYOVERLAP_I_S_a_coefs*I_TWOBODYOVERLAP_I6z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_H_S_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_H_S_a_coefs = alpha;
    I_TWOBODYOVERLAP_H5x_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H5x_S_vrr;
    I_TWOBODYOVERLAP_H4xy_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H4xy_S_vrr;
    I_TWOBODYOVERLAP_H4xz_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H4xz_S_vrr;
    I_TWOBODYOVERLAP_H3x2y_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    I_TWOBODYOVERLAP_H3xyz_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H3xyz_S_vrr;
    I_TWOBODYOVERLAP_H3x2z_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    I_TWOBODYOVERLAP_H2x3y_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    I_TWOBODYOVERLAP_H2x2yz_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H2x2yz_S_vrr;
    I_TWOBODYOVERLAP_H2xy2z_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H2xy2z_S_vrr;
    I_TWOBODYOVERLAP_H2x3z_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    I_TWOBODYOVERLAP_Hx4y_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    I_TWOBODYOVERLAP_Hx3yz_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_Hx3yz_S_vrr;
    I_TWOBODYOVERLAP_Hx2y2z_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_Hx2y2z_S_vrr;
    I_TWOBODYOVERLAP_Hxy3z_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_Hxy3z_S_vrr;
    I_TWOBODYOVERLAP_Hx4z_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    I_TWOBODYOVERLAP_H5y_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H5y_S_vrr;
    I_TWOBODYOVERLAP_H4yz_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H4yz_S_vrr;
    I_TWOBODYOVERLAP_H3y2z_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    I_TWOBODYOVERLAP_H2y3z_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    I_TWOBODYOVERLAP_Hy4z_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    I_TWOBODYOVERLAP_H5z_S_a += SQ_TWOBODYOVERLAP_H_S_a_coefs*I_TWOBODYOVERLAP_H5z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_K_S
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_TWOBODYOVERLAP_K7x_S += I_TWOBODYOVERLAP_K7x_S_vrr;
    I_TWOBODYOVERLAP_K6xy_S += I_TWOBODYOVERLAP_K6xy_S_vrr;
    I_TWOBODYOVERLAP_K6xz_S += I_TWOBODYOVERLAP_K6xz_S_vrr;
    I_TWOBODYOVERLAP_K5x2y_S += I_TWOBODYOVERLAP_K5x2y_S_vrr;
    I_TWOBODYOVERLAP_K5xyz_S += I_TWOBODYOVERLAP_K5xyz_S_vrr;
    I_TWOBODYOVERLAP_K5x2z_S += I_TWOBODYOVERLAP_K5x2z_S_vrr;
    I_TWOBODYOVERLAP_K4x3y_S += I_TWOBODYOVERLAP_K4x3y_S_vrr;
    I_TWOBODYOVERLAP_K4x2yz_S += I_TWOBODYOVERLAP_K4x2yz_S_vrr;
    I_TWOBODYOVERLAP_K4xy2z_S += I_TWOBODYOVERLAP_K4xy2z_S_vrr;
    I_TWOBODYOVERLAP_K4x3z_S += I_TWOBODYOVERLAP_K4x3z_S_vrr;
    I_TWOBODYOVERLAP_K3x4y_S += I_TWOBODYOVERLAP_K3x4y_S_vrr;
    I_TWOBODYOVERLAP_K3x3yz_S += I_TWOBODYOVERLAP_K3x3yz_S_vrr;
    I_TWOBODYOVERLAP_K3x2y2z_S += I_TWOBODYOVERLAP_K3x2y2z_S_vrr;
    I_TWOBODYOVERLAP_K3xy3z_S += I_TWOBODYOVERLAP_K3xy3z_S_vrr;
    I_TWOBODYOVERLAP_K3x4z_S += I_TWOBODYOVERLAP_K3x4z_S_vrr;
    I_TWOBODYOVERLAP_K2x5y_S += I_TWOBODYOVERLAP_K2x5y_S_vrr;
    I_TWOBODYOVERLAP_K2x4yz_S += I_TWOBODYOVERLAP_K2x4yz_S_vrr;
    I_TWOBODYOVERLAP_K2x3y2z_S += I_TWOBODYOVERLAP_K2x3y2z_S_vrr;
    I_TWOBODYOVERLAP_K2x2y3z_S += I_TWOBODYOVERLAP_K2x2y3z_S_vrr;
    I_TWOBODYOVERLAP_K2xy4z_S += I_TWOBODYOVERLAP_K2xy4z_S_vrr;
    I_TWOBODYOVERLAP_K2x5z_S += I_TWOBODYOVERLAP_K2x5z_S_vrr;
    I_TWOBODYOVERLAP_Kx6y_S += I_TWOBODYOVERLAP_Kx6y_S_vrr;
    I_TWOBODYOVERLAP_Kx5yz_S += I_TWOBODYOVERLAP_Kx5yz_S_vrr;
    I_TWOBODYOVERLAP_Kx4y2z_S += I_TWOBODYOVERLAP_Kx4y2z_S_vrr;
    I_TWOBODYOVERLAP_Kx3y3z_S += I_TWOBODYOVERLAP_Kx3y3z_S_vrr;
    I_TWOBODYOVERLAP_Kx2y4z_S += I_TWOBODYOVERLAP_Kx2y4z_S_vrr;
    I_TWOBODYOVERLAP_Kxy5z_S += I_TWOBODYOVERLAP_Kxy5z_S_vrr;
    I_TWOBODYOVERLAP_Kx6z_S += I_TWOBODYOVERLAP_Kx6z_S_vrr;
    I_TWOBODYOVERLAP_K7y_S += I_TWOBODYOVERLAP_K7y_S_vrr;
    I_TWOBODYOVERLAP_K6yz_S += I_TWOBODYOVERLAP_K6yz_S_vrr;
    I_TWOBODYOVERLAP_K5y2z_S += I_TWOBODYOVERLAP_K5y2z_S_vrr;
    I_TWOBODYOVERLAP_K4y3z_S += I_TWOBODYOVERLAP_K4y3z_S_vrr;
    I_TWOBODYOVERLAP_K3y4z_S += I_TWOBODYOVERLAP_K3y4z_S_vrr;
    I_TWOBODYOVERLAP_K2y5z_S += I_TWOBODYOVERLAP_K2y5z_S_vrr;
    I_TWOBODYOVERLAP_Ky6z_S += I_TWOBODYOVERLAP_Ky6z_S_vrr;
    I_TWOBODYOVERLAP_K7z_S += I_TWOBODYOVERLAP_K7z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_I_S
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_TWOBODYOVERLAP_I6x_S += I_TWOBODYOVERLAP_I6x_S_vrr;
    I_TWOBODYOVERLAP_I5xy_S += I_TWOBODYOVERLAP_I5xy_S_vrr;
    I_TWOBODYOVERLAP_I5xz_S += I_TWOBODYOVERLAP_I5xz_S_vrr;
    I_TWOBODYOVERLAP_I4x2y_S += I_TWOBODYOVERLAP_I4x2y_S_vrr;
    I_TWOBODYOVERLAP_I4xyz_S += I_TWOBODYOVERLAP_I4xyz_S_vrr;
    I_TWOBODYOVERLAP_I4x2z_S += I_TWOBODYOVERLAP_I4x2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3y_S += I_TWOBODYOVERLAP_I3x3y_S_vrr;
    I_TWOBODYOVERLAP_I3x2yz_S += I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    I_TWOBODYOVERLAP_I3xy2z_S += I_TWOBODYOVERLAP_I3xy2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3z_S += I_TWOBODYOVERLAP_I3x3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4y_S += I_TWOBODYOVERLAP_I2x4y_S_vrr;
    I_TWOBODYOVERLAP_I2x3yz_S += I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    I_TWOBODYOVERLAP_I2x2y2z_S += I_TWOBODYOVERLAP_I2x2y2z_S_vrr;
    I_TWOBODYOVERLAP_I2xy3z_S += I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4z_S += I_TWOBODYOVERLAP_I2x4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5y_S += I_TWOBODYOVERLAP_Ix5y_S_vrr;
    I_TWOBODYOVERLAP_Ix4yz_S += I_TWOBODYOVERLAP_Ix4yz_S_vrr;
    I_TWOBODYOVERLAP_Ix3y2z_S += I_TWOBODYOVERLAP_Ix3y2z_S_vrr;
    I_TWOBODYOVERLAP_Ix2y3z_S += I_TWOBODYOVERLAP_Ix2y3z_S_vrr;
    I_TWOBODYOVERLAP_Ixy4z_S += I_TWOBODYOVERLAP_Ixy4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5z_S += I_TWOBODYOVERLAP_Ix5z_S_vrr;
    I_TWOBODYOVERLAP_I6y_S += I_TWOBODYOVERLAP_I6y_S_vrr;
    I_TWOBODYOVERLAP_I5yz_S += I_TWOBODYOVERLAP_I5yz_S_vrr;
    I_TWOBODYOVERLAP_I4y2z_S += I_TWOBODYOVERLAP_I4y2z_S_vrr;
    I_TWOBODYOVERLAP_I3y3z_S += I_TWOBODYOVERLAP_I3y3z_S_vrr;
    I_TWOBODYOVERLAP_I2y4z_S += I_TWOBODYOVERLAP_I2y4z_S_vrr;
    I_TWOBODYOVERLAP_Iy5z_S += I_TWOBODYOVERLAP_Iy5z_S_vrr;
    I_TWOBODYOVERLAP_I6z_S += I_TWOBODYOVERLAP_I6z_S_vrr;

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
   * totally 0 integrals are omitted 
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
  Double I_TWOBODYOVERLAP_G4x_Py = I_TWOBODYOVERLAP_H4xy_S+ABY*I_TWOBODYOVERLAP_G4x_S;
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
  Double I_TWOBODYOVERLAP_G4x_Pz = I_TWOBODYOVERLAP_H4xz_S+ABZ*I_TWOBODYOVERLAP_G4x_S;
  Double I_TWOBODYOVERLAP_G3xy_Pz = I_TWOBODYOVERLAP_H3xyz_S+ABZ*I_TWOBODYOVERLAP_G3xy_S;
  Double I_TWOBODYOVERLAP_G3xz_Pz = I_TWOBODYOVERLAP_H3x2z_S+ABZ*I_TWOBODYOVERLAP_G3xz_S;
  Double I_TWOBODYOVERLAP_G2x2y_Pz = I_TWOBODYOVERLAP_H2x2yz_S+ABZ*I_TWOBODYOVERLAP_G2x2y_S;
  Double I_TWOBODYOVERLAP_G2xyz_Pz = I_TWOBODYOVERLAP_H2xy2z_S+ABZ*I_TWOBODYOVERLAP_G2xyz_S;
  Double I_TWOBODYOVERLAP_G2x2z_Pz = I_TWOBODYOVERLAP_H2x3z_S+ABZ*I_TWOBODYOVERLAP_G2x2z_S;
  Double I_TWOBODYOVERLAP_Gx3y_Pz = I_TWOBODYOVERLAP_Hx3yz_S+ABZ*I_TWOBODYOVERLAP_Gx3y_S;
  Double I_TWOBODYOVERLAP_Gx2yz_Pz = I_TWOBODYOVERLAP_Hx2y2z_S+ABZ*I_TWOBODYOVERLAP_Gx2yz_S;
  Double I_TWOBODYOVERLAP_Gxy2z_Pz = I_TWOBODYOVERLAP_Hxy3z_S+ABZ*I_TWOBODYOVERLAP_Gxy2z_S;
  Double I_TWOBODYOVERLAP_Gx3z_Pz = I_TWOBODYOVERLAP_Hx4z_S+ABZ*I_TWOBODYOVERLAP_Gx3z_S;
  Double I_TWOBODYOVERLAP_G4y_Pz = I_TWOBODYOVERLAP_H4yz_S+ABZ*I_TWOBODYOVERLAP_G4y_S;
  Double I_TWOBODYOVERLAP_G3yz_Pz = I_TWOBODYOVERLAP_H3y2z_S+ABZ*I_TWOBODYOVERLAP_G3yz_S;
  Double I_TWOBODYOVERLAP_G2y2z_Pz = I_TWOBODYOVERLAP_H2y3z_S+ABZ*I_TWOBODYOVERLAP_G2y2z_S;
  Double I_TWOBODYOVERLAP_Gy3z_Pz = I_TWOBODYOVERLAP_Hy4z_S+ABZ*I_TWOBODYOVERLAP_Gy3z_S;
  Double I_TWOBODYOVERLAP_G4z_Pz = I_TWOBODYOVERLAP_H5z_S+ABZ*I_TWOBODYOVERLAP_G4z_S;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_D
   * expanding position: BRA2
   * code section is: HRR
   * totally 30 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
   ************************************************************/
  Double I_TWOBODYOVERLAP_F3x_D2x = I_TWOBODYOVERLAP_G4x_Px+ABX*I_TWOBODYOVERLAP_F3x_Px;
  Double I_TWOBODYOVERLAP_F2xy_D2x = I_TWOBODYOVERLAP_G3xy_Px+ABX*I_TWOBODYOVERLAP_F2xy_Px;
  Double I_TWOBODYOVERLAP_F2xz_D2x = I_TWOBODYOVERLAP_G3xz_Px+ABX*I_TWOBODYOVERLAP_F2xz_Px;
  Double I_TWOBODYOVERLAP_Fx2y_D2x = I_TWOBODYOVERLAP_G2x2y_Px+ABX*I_TWOBODYOVERLAP_Fx2y_Px;
  Double I_TWOBODYOVERLAP_Fxyz_D2x = I_TWOBODYOVERLAP_G2xyz_Px+ABX*I_TWOBODYOVERLAP_Fxyz_Px;
  Double I_TWOBODYOVERLAP_Fx2z_D2x = I_TWOBODYOVERLAP_G2x2z_Px+ABX*I_TWOBODYOVERLAP_Fx2z_Px;
  Double I_TWOBODYOVERLAP_F3y_D2x = I_TWOBODYOVERLAP_Gx3y_Px+ABX*I_TWOBODYOVERLAP_F3y_Px;
  Double I_TWOBODYOVERLAP_F2yz_D2x = I_TWOBODYOVERLAP_Gx2yz_Px+ABX*I_TWOBODYOVERLAP_F2yz_Px;
  Double I_TWOBODYOVERLAP_Fy2z_D2x = I_TWOBODYOVERLAP_Gxy2z_Px+ABX*I_TWOBODYOVERLAP_Fy2z_Px;
  Double I_TWOBODYOVERLAP_F3z_D2x = I_TWOBODYOVERLAP_Gx3z_Px+ABX*I_TWOBODYOVERLAP_F3z_Px;
  Double I_TWOBODYOVERLAP_F3x_D2y = I_TWOBODYOVERLAP_G3xy_Py+ABY*I_TWOBODYOVERLAP_F3x_Py;
  Double I_TWOBODYOVERLAP_F2xy_D2y = I_TWOBODYOVERLAP_G2x2y_Py+ABY*I_TWOBODYOVERLAP_F2xy_Py;
  Double I_TWOBODYOVERLAP_F2xz_D2y = I_TWOBODYOVERLAP_G2xyz_Py+ABY*I_TWOBODYOVERLAP_F2xz_Py;
  Double I_TWOBODYOVERLAP_Fx2y_D2y = I_TWOBODYOVERLAP_Gx3y_Py+ABY*I_TWOBODYOVERLAP_Fx2y_Py;
  Double I_TWOBODYOVERLAP_Fxyz_D2y = I_TWOBODYOVERLAP_Gx2yz_Py+ABY*I_TWOBODYOVERLAP_Fxyz_Py;
  Double I_TWOBODYOVERLAP_Fx2z_D2y = I_TWOBODYOVERLAP_Gxy2z_Py+ABY*I_TWOBODYOVERLAP_Fx2z_Py;
  Double I_TWOBODYOVERLAP_F3y_D2y = I_TWOBODYOVERLAP_G4y_Py+ABY*I_TWOBODYOVERLAP_F3y_Py;
  Double I_TWOBODYOVERLAP_F2yz_D2y = I_TWOBODYOVERLAP_G3yz_Py+ABY*I_TWOBODYOVERLAP_F2yz_Py;
  Double I_TWOBODYOVERLAP_Fy2z_D2y = I_TWOBODYOVERLAP_G2y2z_Py+ABY*I_TWOBODYOVERLAP_Fy2z_Py;
  Double I_TWOBODYOVERLAP_F3z_D2y = I_TWOBODYOVERLAP_Gy3z_Py+ABY*I_TWOBODYOVERLAP_F3z_Py;
  Double I_TWOBODYOVERLAP_F3x_D2z = I_TWOBODYOVERLAP_G3xz_Pz+ABZ*I_TWOBODYOVERLAP_F3x_Pz;
  Double I_TWOBODYOVERLAP_F2xy_D2z = I_TWOBODYOVERLAP_G2xyz_Pz+ABZ*I_TWOBODYOVERLAP_F2xy_Pz;
  Double I_TWOBODYOVERLAP_F2xz_D2z = I_TWOBODYOVERLAP_G2x2z_Pz+ABZ*I_TWOBODYOVERLAP_F2xz_Pz;
  Double I_TWOBODYOVERLAP_Fx2y_D2z = I_TWOBODYOVERLAP_Gx2yz_Pz+ABZ*I_TWOBODYOVERLAP_Fx2y_Pz;
  Double I_TWOBODYOVERLAP_Fxyz_D2z = I_TWOBODYOVERLAP_Gxy2z_Pz+ABZ*I_TWOBODYOVERLAP_Fxyz_Pz;
  Double I_TWOBODYOVERLAP_Fx2z_D2z = I_TWOBODYOVERLAP_Gx3z_Pz+ABZ*I_TWOBODYOVERLAP_Fx2z_Pz;
  Double I_TWOBODYOVERLAP_F3y_D2z = I_TWOBODYOVERLAP_G3yz_Pz+ABZ*I_TWOBODYOVERLAP_F3y_Pz;
  Double I_TWOBODYOVERLAP_F2yz_D2z = I_TWOBODYOVERLAP_G2y2z_Pz+ABZ*I_TWOBODYOVERLAP_F2yz_Pz;
  Double I_TWOBODYOVERLAP_Fy2z_D2z = I_TWOBODYOVERLAP_Gy3z_Pz+ABZ*I_TWOBODYOVERLAP_Fy2z_Pz;
  Double I_TWOBODYOVERLAP_F3z_D2z = I_TWOBODYOVERLAP_G4z_Pz+ABZ*I_TWOBODYOVERLAP_F3z_Pz;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 3 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
   ************************************************************/
  Double I_TWOBODYOVERLAP_H5x_Px = I_TWOBODYOVERLAP_I6x_S+ABX*I_TWOBODYOVERLAP_H5x_S;
  Double I_TWOBODYOVERLAP_H4xy_Px = I_TWOBODYOVERLAP_I5xy_S+ABX*I_TWOBODYOVERLAP_H4xy_S;
  Double I_TWOBODYOVERLAP_H4xz_Px = I_TWOBODYOVERLAP_I5xz_S+ABX*I_TWOBODYOVERLAP_H4xz_S;
  Double I_TWOBODYOVERLAP_H3x2y_Px = I_TWOBODYOVERLAP_I4x2y_S+ABX*I_TWOBODYOVERLAP_H3x2y_S;
  Double I_TWOBODYOVERLAP_H3xyz_Px = I_TWOBODYOVERLAP_I4xyz_S+ABX*I_TWOBODYOVERLAP_H3xyz_S;
  Double I_TWOBODYOVERLAP_H3x2z_Px = I_TWOBODYOVERLAP_I4x2z_S+ABX*I_TWOBODYOVERLAP_H3x2z_S;
  Double I_TWOBODYOVERLAP_H2x3y_Px = I_TWOBODYOVERLAP_I3x3y_S+ABX*I_TWOBODYOVERLAP_H2x3y_S;
  Double I_TWOBODYOVERLAP_H2x2yz_Px = I_TWOBODYOVERLAP_I3x2yz_S+ABX*I_TWOBODYOVERLAP_H2x2yz_S;
  Double I_TWOBODYOVERLAP_H2xy2z_Px = I_TWOBODYOVERLAP_I3xy2z_S+ABX*I_TWOBODYOVERLAP_H2xy2z_S;
  Double I_TWOBODYOVERLAP_H2x3z_Px = I_TWOBODYOVERLAP_I3x3z_S+ABX*I_TWOBODYOVERLAP_H2x3z_S;
  Double I_TWOBODYOVERLAP_Hx4y_Px = I_TWOBODYOVERLAP_I2x4y_S+ABX*I_TWOBODYOVERLAP_Hx4y_S;
  Double I_TWOBODYOVERLAP_Hx3yz_Px = I_TWOBODYOVERLAP_I2x3yz_S+ABX*I_TWOBODYOVERLAP_Hx3yz_S;
  Double I_TWOBODYOVERLAP_Hx2y2z_Px = I_TWOBODYOVERLAP_I2x2y2z_S+ABX*I_TWOBODYOVERLAP_Hx2y2z_S;
  Double I_TWOBODYOVERLAP_Hxy3z_Px = I_TWOBODYOVERLAP_I2xy3z_S+ABX*I_TWOBODYOVERLAP_Hxy3z_S;
  Double I_TWOBODYOVERLAP_Hx4z_Px = I_TWOBODYOVERLAP_I2x4z_S+ABX*I_TWOBODYOVERLAP_Hx4z_S;
  Double I_TWOBODYOVERLAP_H5y_Px = I_TWOBODYOVERLAP_Ix5y_S+ABX*I_TWOBODYOVERLAP_H5y_S;
  Double I_TWOBODYOVERLAP_H4yz_Px = I_TWOBODYOVERLAP_Ix4yz_S+ABX*I_TWOBODYOVERLAP_H4yz_S;
  Double I_TWOBODYOVERLAP_H3y2z_Px = I_TWOBODYOVERLAP_Ix3y2z_S+ABX*I_TWOBODYOVERLAP_H3y2z_S;
  Double I_TWOBODYOVERLAP_H2y3z_Px = I_TWOBODYOVERLAP_Ix2y3z_S+ABX*I_TWOBODYOVERLAP_H2y3z_S;
  Double I_TWOBODYOVERLAP_Hy4z_Px = I_TWOBODYOVERLAP_Ixy4z_S+ABX*I_TWOBODYOVERLAP_Hy4z_S;
  Double I_TWOBODYOVERLAP_H5z_Px = I_TWOBODYOVERLAP_Ix5z_S+ABX*I_TWOBODYOVERLAP_H5z_S;
  Double I_TWOBODYOVERLAP_H4xy_Py = I_TWOBODYOVERLAP_I4x2y_S+ABY*I_TWOBODYOVERLAP_H4xy_S;
  Double I_TWOBODYOVERLAP_H4xz_Py = I_TWOBODYOVERLAP_I4xyz_S+ABY*I_TWOBODYOVERLAP_H4xz_S;
  Double I_TWOBODYOVERLAP_H3x2y_Py = I_TWOBODYOVERLAP_I3x3y_S+ABY*I_TWOBODYOVERLAP_H3x2y_S;
  Double I_TWOBODYOVERLAP_H3xyz_Py = I_TWOBODYOVERLAP_I3x2yz_S+ABY*I_TWOBODYOVERLAP_H3xyz_S;
  Double I_TWOBODYOVERLAP_H3x2z_Py = I_TWOBODYOVERLAP_I3xy2z_S+ABY*I_TWOBODYOVERLAP_H3x2z_S;
  Double I_TWOBODYOVERLAP_H2x3y_Py = I_TWOBODYOVERLAP_I2x4y_S+ABY*I_TWOBODYOVERLAP_H2x3y_S;
  Double I_TWOBODYOVERLAP_H2x2yz_Py = I_TWOBODYOVERLAP_I2x3yz_S+ABY*I_TWOBODYOVERLAP_H2x2yz_S;
  Double I_TWOBODYOVERLAP_H2xy2z_Py = I_TWOBODYOVERLAP_I2x2y2z_S+ABY*I_TWOBODYOVERLAP_H2xy2z_S;
  Double I_TWOBODYOVERLAP_H2x3z_Py = I_TWOBODYOVERLAP_I2xy3z_S+ABY*I_TWOBODYOVERLAP_H2x3z_S;
  Double I_TWOBODYOVERLAP_Hx4y_Py = I_TWOBODYOVERLAP_Ix5y_S+ABY*I_TWOBODYOVERLAP_Hx4y_S;
  Double I_TWOBODYOVERLAP_Hx3yz_Py = I_TWOBODYOVERLAP_Ix4yz_S+ABY*I_TWOBODYOVERLAP_Hx3yz_S;
  Double I_TWOBODYOVERLAP_Hx2y2z_Py = I_TWOBODYOVERLAP_Ix3y2z_S+ABY*I_TWOBODYOVERLAP_Hx2y2z_S;
  Double I_TWOBODYOVERLAP_Hxy3z_Py = I_TWOBODYOVERLAP_Ix2y3z_S+ABY*I_TWOBODYOVERLAP_Hxy3z_S;
  Double I_TWOBODYOVERLAP_Hx4z_Py = I_TWOBODYOVERLAP_Ixy4z_S+ABY*I_TWOBODYOVERLAP_Hx4z_S;
  Double I_TWOBODYOVERLAP_H5y_Py = I_TWOBODYOVERLAP_I6y_S+ABY*I_TWOBODYOVERLAP_H5y_S;
  Double I_TWOBODYOVERLAP_H4yz_Py = I_TWOBODYOVERLAP_I5yz_S+ABY*I_TWOBODYOVERLAP_H4yz_S;
  Double I_TWOBODYOVERLAP_H3y2z_Py = I_TWOBODYOVERLAP_I4y2z_S+ABY*I_TWOBODYOVERLAP_H3y2z_S;
  Double I_TWOBODYOVERLAP_H2y3z_Py = I_TWOBODYOVERLAP_I3y3z_S+ABY*I_TWOBODYOVERLAP_H2y3z_S;
  Double I_TWOBODYOVERLAP_Hy4z_Py = I_TWOBODYOVERLAP_I2y4z_S+ABY*I_TWOBODYOVERLAP_Hy4z_S;
  Double I_TWOBODYOVERLAP_H5z_Py = I_TWOBODYOVERLAP_Iy5z_S+ABY*I_TWOBODYOVERLAP_H5z_S;
  Double I_TWOBODYOVERLAP_H4xy_Pz = I_TWOBODYOVERLAP_I4xyz_S+ABZ*I_TWOBODYOVERLAP_H4xy_S;
  Double I_TWOBODYOVERLAP_H4xz_Pz = I_TWOBODYOVERLAP_I4x2z_S+ABZ*I_TWOBODYOVERLAP_H4xz_S;
  Double I_TWOBODYOVERLAP_H3x2y_Pz = I_TWOBODYOVERLAP_I3x2yz_S+ABZ*I_TWOBODYOVERLAP_H3x2y_S;
  Double I_TWOBODYOVERLAP_H3xyz_Pz = I_TWOBODYOVERLAP_I3xy2z_S+ABZ*I_TWOBODYOVERLAP_H3xyz_S;
  Double I_TWOBODYOVERLAP_H3x2z_Pz = I_TWOBODYOVERLAP_I3x3z_S+ABZ*I_TWOBODYOVERLAP_H3x2z_S;
  Double I_TWOBODYOVERLAP_H2x3y_Pz = I_TWOBODYOVERLAP_I2x3yz_S+ABZ*I_TWOBODYOVERLAP_H2x3y_S;
  Double I_TWOBODYOVERLAP_H2x2yz_Pz = I_TWOBODYOVERLAP_I2x2y2z_S+ABZ*I_TWOBODYOVERLAP_H2x2yz_S;
  Double I_TWOBODYOVERLAP_H2xy2z_Pz = I_TWOBODYOVERLAP_I2xy3z_S+ABZ*I_TWOBODYOVERLAP_H2xy2z_S;
  Double I_TWOBODYOVERLAP_H2x3z_Pz = I_TWOBODYOVERLAP_I2x4z_S+ABZ*I_TWOBODYOVERLAP_H2x3z_S;
  Double I_TWOBODYOVERLAP_Hx4y_Pz = I_TWOBODYOVERLAP_Ix4yz_S+ABZ*I_TWOBODYOVERLAP_Hx4y_S;
  Double I_TWOBODYOVERLAP_Hx3yz_Pz = I_TWOBODYOVERLAP_Ix3y2z_S+ABZ*I_TWOBODYOVERLAP_Hx3yz_S;
  Double I_TWOBODYOVERLAP_Hx2y2z_Pz = I_TWOBODYOVERLAP_Ix2y3z_S+ABZ*I_TWOBODYOVERLAP_Hx2y2z_S;
  Double I_TWOBODYOVERLAP_Hxy3z_Pz = I_TWOBODYOVERLAP_Ixy4z_S+ABZ*I_TWOBODYOVERLAP_Hxy3z_S;
  Double I_TWOBODYOVERLAP_Hx4z_Pz = I_TWOBODYOVERLAP_Ix5z_S+ABZ*I_TWOBODYOVERLAP_Hx4z_S;
  Double I_TWOBODYOVERLAP_H4yz_Pz = I_TWOBODYOVERLAP_I4y2z_S+ABZ*I_TWOBODYOVERLAP_H4yz_S;
  Double I_TWOBODYOVERLAP_H3y2z_Pz = I_TWOBODYOVERLAP_I3y3z_S+ABZ*I_TWOBODYOVERLAP_H3y2z_S;
  Double I_TWOBODYOVERLAP_H2y3z_Pz = I_TWOBODYOVERLAP_I2y4z_S+ABZ*I_TWOBODYOVERLAP_H2y3z_S;
  Double I_TWOBODYOVERLAP_Hy4z_Pz = I_TWOBODYOVERLAP_Iy5z_S+ABZ*I_TWOBODYOVERLAP_Hy4z_S;
  Double I_TWOBODYOVERLAP_H5z_Pz = I_TWOBODYOVERLAP_I6z_S+ABZ*I_TWOBODYOVERLAP_H5z_S;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_D
   * expanding position: BRA2
   * code section is: HRR
   * totally 45 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P
   ************************************************************/
  Double I_TWOBODYOVERLAP_G4x_D2x = I_TWOBODYOVERLAP_H5x_Px+ABX*I_TWOBODYOVERLAP_G4x_Px;
  Double I_TWOBODYOVERLAP_G3xy_D2x = I_TWOBODYOVERLAP_H4xy_Px+ABX*I_TWOBODYOVERLAP_G3xy_Px;
  Double I_TWOBODYOVERLAP_G3xz_D2x = I_TWOBODYOVERLAP_H4xz_Px+ABX*I_TWOBODYOVERLAP_G3xz_Px;
  Double I_TWOBODYOVERLAP_G2x2y_D2x = I_TWOBODYOVERLAP_H3x2y_Px+ABX*I_TWOBODYOVERLAP_G2x2y_Px;
  Double I_TWOBODYOVERLAP_G2xyz_D2x = I_TWOBODYOVERLAP_H3xyz_Px+ABX*I_TWOBODYOVERLAP_G2xyz_Px;
  Double I_TWOBODYOVERLAP_G2x2z_D2x = I_TWOBODYOVERLAP_H3x2z_Px+ABX*I_TWOBODYOVERLAP_G2x2z_Px;
  Double I_TWOBODYOVERLAP_Gx3y_D2x = I_TWOBODYOVERLAP_H2x3y_Px+ABX*I_TWOBODYOVERLAP_Gx3y_Px;
  Double I_TWOBODYOVERLAP_Gx2yz_D2x = I_TWOBODYOVERLAP_H2x2yz_Px+ABX*I_TWOBODYOVERLAP_Gx2yz_Px;
  Double I_TWOBODYOVERLAP_Gxy2z_D2x = I_TWOBODYOVERLAP_H2xy2z_Px+ABX*I_TWOBODYOVERLAP_Gxy2z_Px;
  Double I_TWOBODYOVERLAP_Gx3z_D2x = I_TWOBODYOVERLAP_H2x3z_Px+ABX*I_TWOBODYOVERLAP_Gx3z_Px;
  Double I_TWOBODYOVERLAP_G4y_D2x = I_TWOBODYOVERLAP_Hx4y_Px+ABX*I_TWOBODYOVERLAP_G4y_Px;
  Double I_TWOBODYOVERLAP_G3yz_D2x = I_TWOBODYOVERLAP_Hx3yz_Px+ABX*I_TWOBODYOVERLAP_G3yz_Px;
  Double I_TWOBODYOVERLAP_G2y2z_D2x = I_TWOBODYOVERLAP_Hx2y2z_Px+ABX*I_TWOBODYOVERLAP_G2y2z_Px;
  Double I_TWOBODYOVERLAP_Gy3z_D2x = I_TWOBODYOVERLAP_Hxy3z_Px+ABX*I_TWOBODYOVERLAP_Gy3z_Px;
  Double I_TWOBODYOVERLAP_G4z_D2x = I_TWOBODYOVERLAP_Hx4z_Px+ABX*I_TWOBODYOVERLAP_G4z_Px;
  Double I_TWOBODYOVERLAP_G4x_D2y = I_TWOBODYOVERLAP_H4xy_Py+ABY*I_TWOBODYOVERLAP_G4x_Py;
  Double I_TWOBODYOVERLAP_G3xy_D2y = I_TWOBODYOVERLAP_H3x2y_Py+ABY*I_TWOBODYOVERLAP_G3xy_Py;
  Double I_TWOBODYOVERLAP_G3xz_D2y = I_TWOBODYOVERLAP_H3xyz_Py+ABY*I_TWOBODYOVERLAP_G3xz_Py;
  Double I_TWOBODYOVERLAP_G2x2y_D2y = I_TWOBODYOVERLAP_H2x3y_Py+ABY*I_TWOBODYOVERLAP_G2x2y_Py;
  Double I_TWOBODYOVERLAP_G2xyz_D2y = I_TWOBODYOVERLAP_H2x2yz_Py+ABY*I_TWOBODYOVERLAP_G2xyz_Py;
  Double I_TWOBODYOVERLAP_G2x2z_D2y = I_TWOBODYOVERLAP_H2xy2z_Py+ABY*I_TWOBODYOVERLAP_G2x2z_Py;
  Double I_TWOBODYOVERLAP_Gx3y_D2y = I_TWOBODYOVERLAP_Hx4y_Py+ABY*I_TWOBODYOVERLAP_Gx3y_Py;
  Double I_TWOBODYOVERLAP_Gx2yz_D2y = I_TWOBODYOVERLAP_Hx3yz_Py+ABY*I_TWOBODYOVERLAP_Gx2yz_Py;
  Double I_TWOBODYOVERLAP_Gxy2z_D2y = I_TWOBODYOVERLAP_Hx2y2z_Py+ABY*I_TWOBODYOVERLAP_Gxy2z_Py;
  Double I_TWOBODYOVERLAP_Gx3z_D2y = I_TWOBODYOVERLAP_Hxy3z_Py+ABY*I_TWOBODYOVERLAP_Gx3z_Py;
  Double I_TWOBODYOVERLAP_G4y_D2y = I_TWOBODYOVERLAP_H5y_Py+ABY*I_TWOBODYOVERLAP_G4y_Py;
  Double I_TWOBODYOVERLAP_G3yz_D2y = I_TWOBODYOVERLAP_H4yz_Py+ABY*I_TWOBODYOVERLAP_G3yz_Py;
  Double I_TWOBODYOVERLAP_G2y2z_D2y = I_TWOBODYOVERLAP_H3y2z_Py+ABY*I_TWOBODYOVERLAP_G2y2z_Py;
  Double I_TWOBODYOVERLAP_Gy3z_D2y = I_TWOBODYOVERLAP_H2y3z_Py+ABY*I_TWOBODYOVERLAP_Gy3z_Py;
  Double I_TWOBODYOVERLAP_G4z_D2y = I_TWOBODYOVERLAP_Hy4z_Py+ABY*I_TWOBODYOVERLAP_G4z_Py;
  Double I_TWOBODYOVERLAP_G4x_D2z = I_TWOBODYOVERLAP_H4xz_Pz+ABZ*I_TWOBODYOVERLAP_G4x_Pz;
  Double I_TWOBODYOVERLAP_G3xy_D2z = I_TWOBODYOVERLAP_H3xyz_Pz+ABZ*I_TWOBODYOVERLAP_G3xy_Pz;
  Double I_TWOBODYOVERLAP_G3xz_D2z = I_TWOBODYOVERLAP_H3x2z_Pz+ABZ*I_TWOBODYOVERLAP_G3xz_Pz;
  Double I_TWOBODYOVERLAP_G2x2y_D2z = I_TWOBODYOVERLAP_H2x2yz_Pz+ABZ*I_TWOBODYOVERLAP_G2x2y_Pz;
  Double I_TWOBODYOVERLAP_G2xyz_D2z = I_TWOBODYOVERLAP_H2xy2z_Pz+ABZ*I_TWOBODYOVERLAP_G2xyz_Pz;
  Double I_TWOBODYOVERLAP_G2x2z_D2z = I_TWOBODYOVERLAP_H2x3z_Pz+ABZ*I_TWOBODYOVERLAP_G2x2z_Pz;
  Double I_TWOBODYOVERLAP_Gx3y_D2z = I_TWOBODYOVERLAP_Hx3yz_Pz+ABZ*I_TWOBODYOVERLAP_Gx3y_Pz;
  Double I_TWOBODYOVERLAP_Gx2yz_D2z = I_TWOBODYOVERLAP_Hx2y2z_Pz+ABZ*I_TWOBODYOVERLAP_Gx2yz_Pz;
  Double I_TWOBODYOVERLAP_Gxy2z_D2z = I_TWOBODYOVERLAP_Hxy3z_Pz+ABZ*I_TWOBODYOVERLAP_Gxy2z_Pz;
  Double I_TWOBODYOVERLAP_Gx3z_D2z = I_TWOBODYOVERLAP_Hx4z_Pz+ABZ*I_TWOBODYOVERLAP_Gx3z_Pz;
  Double I_TWOBODYOVERLAP_G4y_D2z = I_TWOBODYOVERLAP_H4yz_Pz+ABZ*I_TWOBODYOVERLAP_G4y_Pz;
  Double I_TWOBODYOVERLAP_G3yz_D2z = I_TWOBODYOVERLAP_H3y2z_Pz+ABZ*I_TWOBODYOVERLAP_G3yz_Pz;
  Double I_TWOBODYOVERLAP_G2y2z_D2z = I_TWOBODYOVERLAP_H2y3z_Pz+ABZ*I_TWOBODYOVERLAP_G2y2z_Pz;
  Double I_TWOBODYOVERLAP_Gy3z_D2z = I_TWOBODYOVERLAP_Hy4z_Pz+ABZ*I_TWOBODYOVERLAP_Gy3z_Pz;
  Double I_TWOBODYOVERLAP_G4z_D2z = I_TWOBODYOVERLAP_H5z_Pz+ABZ*I_TWOBODYOVERLAP_G4z_Pz;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_F
   * expanding position: BRA2
   * code section is: HRR
   * totally 20 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_D
   ************************************************************/
  Double I_TWOBODYOVERLAP_F3x_F3x = I_TWOBODYOVERLAP_G4x_D2x+ABX*I_TWOBODYOVERLAP_F3x_D2x;
  Double I_TWOBODYOVERLAP_F2xy_F3x = I_TWOBODYOVERLAP_G3xy_D2x+ABX*I_TWOBODYOVERLAP_F2xy_D2x;
  Double I_TWOBODYOVERLAP_F2xz_F3x = I_TWOBODYOVERLAP_G3xz_D2x+ABX*I_TWOBODYOVERLAP_F2xz_D2x;
  Double I_TWOBODYOVERLAP_Fx2y_F3x = I_TWOBODYOVERLAP_G2x2y_D2x+ABX*I_TWOBODYOVERLAP_Fx2y_D2x;
  Double I_TWOBODYOVERLAP_Fxyz_F3x = I_TWOBODYOVERLAP_G2xyz_D2x+ABX*I_TWOBODYOVERLAP_Fxyz_D2x;
  Double I_TWOBODYOVERLAP_Fx2z_F3x = I_TWOBODYOVERLAP_G2x2z_D2x+ABX*I_TWOBODYOVERLAP_Fx2z_D2x;
  Double I_TWOBODYOVERLAP_F3y_F3x = I_TWOBODYOVERLAP_Gx3y_D2x+ABX*I_TWOBODYOVERLAP_F3y_D2x;
  Double I_TWOBODYOVERLAP_F2yz_F3x = I_TWOBODYOVERLAP_Gx2yz_D2x+ABX*I_TWOBODYOVERLAP_F2yz_D2x;
  Double I_TWOBODYOVERLAP_Fy2z_F3x = I_TWOBODYOVERLAP_Gxy2z_D2x+ABX*I_TWOBODYOVERLAP_Fy2z_D2x;
  Double I_TWOBODYOVERLAP_F3z_F3x = I_TWOBODYOVERLAP_Gx3z_D2x+ABX*I_TWOBODYOVERLAP_F3z_D2x;
  Double I_TWOBODYOVERLAP_F3x_F2xy = I_TWOBODYOVERLAP_G3xy_D2x+ABY*I_TWOBODYOVERLAP_F3x_D2x;
  Double I_TWOBODYOVERLAP_F2xy_F2xy = I_TWOBODYOVERLAP_G2x2y_D2x+ABY*I_TWOBODYOVERLAP_F2xy_D2x;
  Double I_TWOBODYOVERLAP_F2xz_F2xy = I_TWOBODYOVERLAP_G2xyz_D2x+ABY*I_TWOBODYOVERLAP_F2xz_D2x;
  Double I_TWOBODYOVERLAP_Fx2y_F2xy = I_TWOBODYOVERLAP_Gx3y_D2x+ABY*I_TWOBODYOVERLAP_Fx2y_D2x;
  Double I_TWOBODYOVERLAP_Fxyz_F2xy = I_TWOBODYOVERLAP_Gx2yz_D2x+ABY*I_TWOBODYOVERLAP_Fxyz_D2x;
  Double I_TWOBODYOVERLAP_Fx2z_F2xy = I_TWOBODYOVERLAP_Gxy2z_D2x+ABY*I_TWOBODYOVERLAP_Fx2z_D2x;
  Double I_TWOBODYOVERLAP_F3y_F2xy = I_TWOBODYOVERLAP_G4y_D2x+ABY*I_TWOBODYOVERLAP_F3y_D2x;
  Double I_TWOBODYOVERLAP_F2yz_F2xy = I_TWOBODYOVERLAP_G3yz_D2x+ABY*I_TWOBODYOVERLAP_F2yz_D2x;
  Double I_TWOBODYOVERLAP_Fy2z_F2xy = I_TWOBODYOVERLAP_G2y2z_D2x+ABY*I_TWOBODYOVERLAP_Fy2z_D2x;
  Double I_TWOBODYOVERLAP_F3z_F2xy = I_TWOBODYOVERLAP_Gy3z_D2x+ABY*I_TWOBODYOVERLAP_F3z_D2x;
  Double I_TWOBODYOVERLAP_F3x_F2xz = I_TWOBODYOVERLAP_G3xz_D2x+ABZ*I_TWOBODYOVERLAP_F3x_D2x;
  Double I_TWOBODYOVERLAP_F2xy_F2xz = I_TWOBODYOVERLAP_G2xyz_D2x+ABZ*I_TWOBODYOVERLAP_F2xy_D2x;
  Double I_TWOBODYOVERLAP_F2xz_F2xz = I_TWOBODYOVERLAP_G2x2z_D2x+ABZ*I_TWOBODYOVERLAP_F2xz_D2x;
  Double I_TWOBODYOVERLAP_Fx2y_F2xz = I_TWOBODYOVERLAP_Gx2yz_D2x+ABZ*I_TWOBODYOVERLAP_Fx2y_D2x;
  Double I_TWOBODYOVERLAP_Fxyz_F2xz = I_TWOBODYOVERLAP_Gxy2z_D2x+ABZ*I_TWOBODYOVERLAP_Fxyz_D2x;
  Double I_TWOBODYOVERLAP_Fx2z_F2xz = I_TWOBODYOVERLAP_Gx3z_D2x+ABZ*I_TWOBODYOVERLAP_Fx2z_D2x;
  Double I_TWOBODYOVERLAP_F3y_F2xz = I_TWOBODYOVERLAP_G3yz_D2x+ABZ*I_TWOBODYOVERLAP_F3y_D2x;
  Double I_TWOBODYOVERLAP_F2yz_F2xz = I_TWOBODYOVERLAP_G2y2z_D2x+ABZ*I_TWOBODYOVERLAP_F2yz_D2x;
  Double I_TWOBODYOVERLAP_Fy2z_F2xz = I_TWOBODYOVERLAP_Gy3z_D2x+ABZ*I_TWOBODYOVERLAP_Fy2z_D2x;
  Double I_TWOBODYOVERLAP_F3z_F2xz = I_TWOBODYOVERLAP_G4z_D2x+ABZ*I_TWOBODYOVERLAP_F3z_D2x;
  Double I_TWOBODYOVERLAP_F3x_Fx2y = I_TWOBODYOVERLAP_G4x_D2y+ABX*I_TWOBODYOVERLAP_F3x_D2y;
  Double I_TWOBODYOVERLAP_F2xy_Fx2y = I_TWOBODYOVERLAP_G3xy_D2y+ABX*I_TWOBODYOVERLAP_F2xy_D2y;
  Double I_TWOBODYOVERLAP_F2xz_Fx2y = I_TWOBODYOVERLAP_G3xz_D2y+ABX*I_TWOBODYOVERLAP_F2xz_D2y;
  Double I_TWOBODYOVERLAP_Fx2y_Fx2y = I_TWOBODYOVERLAP_G2x2y_D2y+ABX*I_TWOBODYOVERLAP_Fx2y_D2y;
  Double I_TWOBODYOVERLAP_Fxyz_Fx2y = I_TWOBODYOVERLAP_G2xyz_D2y+ABX*I_TWOBODYOVERLAP_Fxyz_D2y;
  Double I_TWOBODYOVERLAP_Fx2z_Fx2y = I_TWOBODYOVERLAP_G2x2z_D2y+ABX*I_TWOBODYOVERLAP_Fx2z_D2y;
  Double I_TWOBODYOVERLAP_F3y_Fx2y = I_TWOBODYOVERLAP_Gx3y_D2y+ABX*I_TWOBODYOVERLAP_F3y_D2y;
  Double I_TWOBODYOVERLAP_F2yz_Fx2y = I_TWOBODYOVERLAP_Gx2yz_D2y+ABX*I_TWOBODYOVERLAP_F2yz_D2y;
  Double I_TWOBODYOVERLAP_Fy2z_Fx2y = I_TWOBODYOVERLAP_Gxy2z_D2y+ABX*I_TWOBODYOVERLAP_Fy2z_D2y;
  Double I_TWOBODYOVERLAP_F3z_Fx2y = I_TWOBODYOVERLAP_Gx3z_D2y+ABX*I_TWOBODYOVERLAP_F3z_D2y;
  Double I_TWOBODYOVERLAP_F3x_Fx2z = I_TWOBODYOVERLAP_G4x_D2z+ABX*I_TWOBODYOVERLAP_F3x_D2z;
  Double I_TWOBODYOVERLAP_F2xy_Fx2z = I_TWOBODYOVERLAP_G3xy_D2z+ABX*I_TWOBODYOVERLAP_F2xy_D2z;
  Double I_TWOBODYOVERLAP_F2xz_Fx2z = I_TWOBODYOVERLAP_G3xz_D2z+ABX*I_TWOBODYOVERLAP_F2xz_D2z;
  Double I_TWOBODYOVERLAP_Fx2y_Fx2z = I_TWOBODYOVERLAP_G2x2y_D2z+ABX*I_TWOBODYOVERLAP_Fx2y_D2z;
  Double I_TWOBODYOVERLAP_Fxyz_Fx2z = I_TWOBODYOVERLAP_G2xyz_D2z+ABX*I_TWOBODYOVERLAP_Fxyz_D2z;
  Double I_TWOBODYOVERLAP_Fx2z_Fx2z = I_TWOBODYOVERLAP_G2x2z_D2z+ABX*I_TWOBODYOVERLAP_Fx2z_D2z;
  Double I_TWOBODYOVERLAP_F3y_Fx2z = I_TWOBODYOVERLAP_Gx3y_D2z+ABX*I_TWOBODYOVERLAP_F3y_D2z;
  Double I_TWOBODYOVERLAP_F2yz_Fx2z = I_TWOBODYOVERLAP_Gx2yz_D2z+ABX*I_TWOBODYOVERLAP_F2yz_D2z;
  Double I_TWOBODYOVERLAP_Fy2z_Fx2z = I_TWOBODYOVERLAP_Gxy2z_D2z+ABX*I_TWOBODYOVERLAP_Fy2z_D2z;
  Double I_TWOBODYOVERLAP_F3z_Fx2z = I_TWOBODYOVERLAP_Gx3z_D2z+ABX*I_TWOBODYOVERLAP_F3z_D2z;
  Double I_TWOBODYOVERLAP_F3x_F3y = I_TWOBODYOVERLAP_G3xy_D2y+ABY*I_TWOBODYOVERLAP_F3x_D2y;
  Double I_TWOBODYOVERLAP_F2xy_F3y = I_TWOBODYOVERLAP_G2x2y_D2y+ABY*I_TWOBODYOVERLAP_F2xy_D2y;
  Double I_TWOBODYOVERLAP_F2xz_F3y = I_TWOBODYOVERLAP_G2xyz_D2y+ABY*I_TWOBODYOVERLAP_F2xz_D2y;
  Double I_TWOBODYOVERLAP_Fx2y_F3y = I_TWOBODYOVERLAP_Gx3y_D2y+ABY*I_TWOBODYOVERLAP_Fx2y_D2y;
  Double I_TWOBODYOVERLAP_Fxyz_F3y = I_TWOBODYOVERLAP_Gx2yz_D2y+ABY*I_TWOBODYOVERLAP_Fxyz_D2y;
  Double I_TWOBODYOVERLAP_Fx2z_F3y = I_TWOBODYOVERLAP_Gxy2z_D2y+ABY*I_TWOBODYOVERLAP_Fx2z_D2y;
  Double I_TWOBODYOVERLAP_F3y_F3y = I_TWOBODYOVERLAP_G4y_D2y+ABY*I_TWOBODYOVERLAP_F3y_D2y;
  Double I_TWOBODYOVERLAP_F2yz_F3y = I_TWOBODYOVERLAP_G3yz_D2y+ABY*I_TWOBODYOVERLAP_F2yz_D2y;
  Double I_TWOBODYOVERLAP_Fy2z_F3y = I_TWOBODYOVERLAP_G2y2z_D2y+ABY*I_TWOBODYOVERLAP_Fy2z_D2y;
  Double I_TWOBODYOVERLAP_F3z_F3y = I_TWOBODYOVERLAP_Gy3z_D2y+ABY*I_TWOBODYOVERLAP_F3z_D2y;
  Double I_TWOBODYOVERLAP_F3x_F2yz = I_TWOBODYOVERLAP_G3xz_D2y+ABZ*I_TWOBODYOVERLAP_F3x_D2y;
  Double I_TWOBODYOVERLAP_F2xy_F2yz = I_TWOBODYOVERLAP_G2xyz_D2y+ABZ*I_TWOBODYOVERLAP_F2xy_D2y;
  Double I_TWOBODYOVERLAP_F2xz_F2yz = I_TWOBODYOVERLAP_G2x2z_D2y+ABZ*I_TWOBODYOVERLAP_F2xz_D2y;
  Double I_TWOBODYOVERLAP_Fx2y_F2yz = I_TWOBODYOVERLAP_Gx2yz_D2y+ABZ*I_TWOBODYOVERLAP_Fx2y_D2y;
  Double I_TWOBODYOVERLAP_Fxyz_F2yz = I_TWOBODYOVERLAP_Gxy2z_D2y+ABZ*I_TWOBODYOVERLAP_Fxyz_D2y;
  Double I_TWOBODYOVERLAP_Fx2z_F2yz = I_TWOBODYOVERLAP_Gx3z_D2y+ABZ*I_TWOBODYOVERLAP_Fx2z_D2y;
  Double I_TWOBODYOVERLAP_F3y_F2yz = I_TWOBODYOVERLAP_G3yz_D2y+ABZ*I_TWOBODYOVERLAP_F3y_D2y;
  Double I_TWOBODYOVERLAP_F2yz_F2yz = I_TWOBODYOVERLAP_G2y2z_D2y+ABZ*I_TWOBODYOVERLAP_F2yz_D2y;
  Double I_TWOBODYOVERLAP_Fy2z_F2yz = I_TWOBODYOVERLAP_Gy3z_D2y+ABZ*I_TWOBODYOVERLAP_Fy2z_D2y;
  Double I_TWOBODYOVERLAP_F3z_F2yz = I_TWOBODYOVERLAP_G4z_D2y+ABZ*I_TWOBODYOVERLAP_F3z_D2y;
  Double I_TWOBODYOVERLAP_F3x_F3z = I_TWOBODYOVERLAP_G3xz_D2z+ABZ*I_TWOBODYOVERLAP_F3x_D2z;
  Double I_TWOBODYOVERLAP_F2xy_F3z = I_TWOBODYOVERLAP_G2xyz_D2z+ABZ*I_TWOBODYOVERLAP_F2xy_D2z;
  Double I_TWOBODYOVERLAP_F2xz_F3z = I_TWOBODYOVERLAP_G2x2z_D2z+ABZ*I_TWOBODYOVERLAP_F2xz_D2z;
  Double I_TWOBODYOVERLAP_Fx2y_F3z = I_TWOBODYOVERLAP_Gx2yz_D2z+ABZ*I_TWOBODYOVERLAP_Fx2y_D2z;
  Double I_TWOBODYOVERLAP_Fxyz_F3z = I_TWOBODYOVERLAP_Gxy2z_D2z+ABZ*I_TWOBODYOVERLAP_Fxyz_D2z;
  Double I_TWOBODYOVERLAP_Fx2z_F3z = I_TWOBODYOVERLAP_Gx3z_D2z+ABZ*I_TWOBODYOVERLAP_Fx2z_D2z;
  Double I_TWOBODYOVERLAP_F3y_F3z = I_TWOBODYOVERLAP_G3yz_D2z+ABZ*I_TWOBODYOVERLAP_F3y_D2z;
  Double I_TWOBODYOVERLAP_F2yz_F3z = I_TWOBODYOVERLAP_G2y2z_D2z+ABZ*I_TWOBODYOVERLAP_F2yz_D2z;
  Double I_TWOBODYOVERLAP_Fy2z_F3z = I_TWOBODYOVERLAP_Gy3z_D2z+ABZ*I_TWOBODYOVERLAP_Fy2z_D2z;
  Double I_TWOBODYOVERLAP_F3z_F3z = I_TWOBODYOVERLAP_G4z_D2z+ABZ*I_TWOBODYOVERLAP_F3z_D2z;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_I_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 24 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S
   ************************************************************/
  Double I_TWOBODYOVERLAP_I6x_Px = I_TWOBODYOVERLAP_K7x_S+ABX*I_TWOBODYOVERLAP_I6x_S;
  Double I_TWOBODYOVERLAP_I5xy_Px = I_TWOBODYOVERLAP_K6xy_S+ABX*I_TWOBODYOVERLAP_I5xy_S;
  Double I_TWOBODYOVERLAP_I5xz_Px = I_TWOBODYOVERLAP_K6xz_S+ABX*I_TWOBODYOVERLAP_I5xz_S;
  Double I_TWOBODYOVERLAP_I4x2y_Px = I_TWOBODYOVERLAP_K5x2y_S+ABX*I_TWOBODYOVERLAP_I4x2y_S;
  Double I_TWOBODYOVERLAP_I4xyz_Px = I_TWOBODYOVERLAP_K5xyz_S+ABX*I_TWOBODYOVERLAP_I4xyz_S;
  Double I_TWOBODYOVERLAP_I4x2z_Px = I_TWOBODYOVERLAP_K5x2z_S+ABX*I_TWOBODYOVERLAP_I4x2z_S;
  Double I_TWOBODYOVERLAP_I3x3y_Px = I_TWOBODYOVERLAP_K4x3y_S+ABX*I_TWOBODYOVERLAP_I3x3y_S;
  Double I_TWOBODYOVERLAP_I3x2yz_Px = I_TWOBODYOVERLAP_K4x2yz_S+ABX*I_TWOBODYOVERLAP_I3x2yz_S;
  Double I_TWOBODYOVERLAP_I3xy2z_Px = I_TWOBODYOVERLAP_K4xy2z_S+ABX*I_TWOBODYOVERLAP_I3xy2z_S;
  Double I_TWOBODYOVERLAP_I3x3z_Px = I_TWOBODYOVERLAP_K4x3z_S+ABX*I_TWOBODYOVERLAP_I3x3z_S;
  Double I_TWOBODYOVERLAP_I2x4y_Px = I_TWOBODYOVERLAP_K3x4y_S+ABX*I_TWOBODYOVERLAP_I2x4y_S;
  Double I_TWOBODYOVERLAP_I2x3yz_Px = I_TWOBODYOVERLAP_K3x3yz_S+ABX*I_TWOBODYOVERLAP_I2x3yz_S;
  Double I_TWOBODYOVERLAP_I2x2y2z_Px = I_TWOBODYOVERLAP_K3x2y2z_S+ABX*I_TWOBODYOVERLAP_I2x2y2z_S;
  Double I_TWOBODYOVERLAP_I2xy3z_Px = I_TWOBODYOVERLAP_K3xy3z_S+ABX*I_TWOBODYOVERLAP_I2xy3z_S;
  Double I_TWOBODYOVERLAP_I2x4z_Px = I_TWOBODYOVERLAP_K3x4z_S+ABX*I_TWOBODYOVERLAP_I2x4z_S;
  Double I_TWOBODYOVERLAP_Ix5y_Px = I_TWOBODYOVERLAP_K2x5y_S+ABX*I_TWOBODYOVERLAP_Ix5y_S;
  Double I_TWOBODYOVERLAP_Ix4yz_Px = I_TWOBODYOVERLAP_K2x4yz_S+ABX*I_TWOBODYOVERLAP_Ix4yz_S;
  Double I_TWOBODYOVERLAP_Ix3y2z_Px = I_TWOBODYOVERLAP_K2x3y2z_S+ABX*I_TWOBODYOVERLAP_Ix3y2z_S;
  Double I_TWOBODYOVERLAP_Ix2y3z_Px = I_TWOBODYOVERLAP_K2x2y3z_S+ABX*I_TWOBODYOVERLAP_Ix2y3z_S;
  Double I_TWOBODYOVERLAP_Ixy4z_Px = I_TWOBODYOVERLAP_K2xy4z_S+ABX*I_TWOBODYOVERLAP_Ixy4z_S;
  Double I_TWOBODYOVERLAP_Ix5z_Px = I_TWOBODYOVERLAP_K2x5z_S+ABX*I_TWOBODYOVERLAP_Ix5z_S;
  Double I_TWOBODYOVERLAP_I4x2y_Py = I_TWOBODYOVERLAP_K4x3y_S+ABY*I_TWOBODYOVERLAP_I4x2y_S;
  Double I_TWOBODYOVERLAP_I4xyz_Py = I_TWOBODYOVERLAP_K4x2yz_S+ABY*I_TWOBODYOVERLAP_I4xyz_S;
  Double I_TWOBODYOVERLAP_I3x3y_Py = I_TWOBODYOVERLAP_K3x4y_S+ABY*I_TWOBODYOVERLAP_I3x3y_S;
  Double I_TWOBODYOVERLAP_I3x2yz_Py = I_TWOBODYOVERLAP_K3x3yz_S+ABY*I_TWOBODYOVERLAP_I3x2yz_S;
  Double I_TWOBODYOVERLAP_I3xy2z_Py = I_TWOBODYOVERLAP_K3x2y2z_S+ABY*I_TWOBODYOVERLAP_I3xy2z_S;
  Double I_TWOBODYOVERLAP_I2x4y_Py = I_TWOBODYOVERLAP_K2x5y_S+ABY*I_TWOBODYOVERLAP_I2x4y_S;
  Double I_TWOBODYOVERLAP_I2x3yz_Py = I_TWOBODYOVERLAP_K2x4yz_S+ABY*I_TWOBODYOVERLAP_I2x3yz_S;
  Double I_TWOBODYOVERLAP_I2x2y2z_Py = I_TWOBODYOVERLAP_K2x3y2z_S+ABY*I_TWOBODYOVERLAP_I2x2y2z_S;
  Double I_TWOBODYOVERLAP_I2xy3z_Py = I_TWOBODYOVERLAP_K2x2y3z_S+ABY*I_TWOBODYOVERLAP_I2xy3z_S;
  Double I_TWOBODYOVERLAP_Ix5y_Py = I_TWOBODYOVERLAP_Kx6y_S+ABY*I_TWOBODYOVERLAP_Ix5y_S;
  Double I_TWOBODYOVERLAP_Ix4yz_Py = I_TWOBODYOVERLAP_Kx5yz_S+ABY*I_TWOBODYOVERLAP_Ix4yz_S;
  Double I_TWOBODYOVERLAP_Ix3y2z_Py = I_TWOBODYOVERLAP_Kx4y2z_S+ABY*I_TWOBODYOVERLAP_Ix3y2z_S;
  Double I_TWOBODYOVERLAP_Ix2y3z_Py = I_TWOBODYOVERLAP_Kx3y3z_S+ABY*I_TWOBODYOVERLAP_Ix2y3z_S;
  Double I_TWOBODYOVERLAP_Ixy4z_Py = I_TWOBODYOVERLAP_Kx2y4z_S+ABY*I_TWOBODYOVERLAP_Ixy4z_S;
  Double I_TWOBODYOVERLAP_I6y_Py = I_TWOBODYOVERLAP_K7y_S+ABY*I_TWOBODYOVERLAP_I6y_S;
  Double I_TWOBODYOVERLAP_I5yz_Py = I_TWOBODYOVERLAP_K6yz_S+ABY*I_TWOBODYOVERLAP_I5yz_S;
  Double I_TWOBODYOVERLAP_I4y2z_Py = I_TWOBODYOVERLAP_K5y2z_S+ABY*I_TWOBODYOVERLAP_I4y2z_S;
  Double I_TWOBODYOVERLAP_I3y3z_Py = I_TWOBODYOVERLAP_K4y3z_S+ABY*I_TWOBODYOVERLAP_I3y3z_S;
  Double I_TWOBODYOVERLAP_I2y4z_Py = I_TWOBODYOVERLAP_K3y4z_S+ABY*I_TWOBODYOVERLAP_I2y4z_S;
  Double I_TWOBODYOVERLAP_Iy5z_Py = I_TWOBODYOVERLAP_K2y5z_S+ABY*I_TWOBODYOVERLAP_Iy5z_S;
  Double I_TWOBODYOVERLAP_I4xyz_Pz = I_TWOBODYOVERLAP_K4xy2z_S+ABZ*I_TWOBODYOVERLAP_I4xyz_S;
  Double I_TWOBODYOVERLAP_I4x2z_Pz = I_TWOBODYOVERLAP_K4x3z_S+ABZ*I_TWOBODYOVERLAP_I4x2z_S;
  Double I_TWOBODYOVERLAP_I3x2yz_Pz = I_TWOBODYOVERLAP_K3x2y2z_S+ABZ*I_TWOBODYOVERLAP_I3x2yz_S;
  Double I_TWOBODYOVERLAP_I3xy2z_Pz = I_TWOBODYOVERLAP_K3xy3z_S+ABZ*I_TWOBODYOVERLAP_I3xy2z_S;
  Double I_TWOBODYOVERLAP_I3x3z_Pz = I_TWOBODYOVERLAP_K3x4z_S+ABZ*I_TWOBODYOVERLAP_I3x3z_S;
  Double I_TWOBODYOVERLAP_I2x3yz_Pz = I_TWOBODYOVERLAP_K2x3y2z_S+ABZ*I_TWOBODYOVERLAP_I2x3yz_S;
  Double I_TWOBODYOVERLAP_I2x2y2z_Pz = I_TWOBODYOVERLAP_K2x2y3z_S+ABZ*I_TWOBODYOVERLAP_I2x2y2z_S;
  Double I_TWOBODYOVERLAP_I2xy3z_Pz = I_TWOBODYOVERLAP_K2xy4z_S+ABZ*I_TWOBODYOVERLAP_I2xy3z_S;
  Double I_TWOBODYOVERLAP_I2x4z_Pz = I_TWOBODYOVERLAP_K2x5z_S+ABZ*I_TWOBODYOVERLAP_I2x4z_S;
  Double I_TWOBODYOVERLAP_Ix4yz_Pz = I_TWOBODYOVERLAP_Kx4y2z_S+ABZ*I_TWOBODYOVERLAP_Ix4yz_S;
  Double I_TWOBODYOVERLAP_Ix3y2z_Pz = I_TWOBODYOVERLAP_Kx3y3z_S+ABZ*I_TWOBODYOVERLAP_Ix3y2z_S;
  Double I_TWOBODYOVERLAP_Ix2y3z_Pz = I_TWOBODYOVERLAP_Kx2y4z_S+ABZ*I_TWOBODYOVERLAP_Ix2y3z_S;
  Double I_TWOBODYOVERLAP_Ixy4z_Pz = I_TWOBODYOVERLAP_Kxy5z_S+ABZ*I_TWOBODYOVERLAP_Ixy4z_S;
  Double I_TWOBODYOVERLAP_Ix5z_Pz = I_TWOBODYOVERLAP_Kx6z_S+ABZ*I_TWOBODYOVERLAP_Ix5z_S;
  Double I_TWOBODYOVERLAP_I4y2z_Pz = I_TWOBODYOVERLAP_K4y3z_S+ABZ*I_TWOBODYOVERLAP_I4y2z_S;
  Double I_TWOBODYOVERLAP_I3y3z_Pz = I_TWOBODYOVERLAP_K3y4z_S+ABZ*I_TWOBODYOVERLAP_I3y3z_S;
  Double I_TWOBODYOVERLAP_I2y4z_Pz = I_TWOBODYOVERLAP_K2y5z_S+ABZ*I_TWOBODYOVERLAP_I2y4z_S;
  Double I_TWOBODYOVERLAP_Iy5z_Pz = I_TWOBODYOVERLAP_Ky6z_S+ABZ*I_TWOBODYOVERLAP_Iy5z_S;
  Double I_TWOBODYOVERLAP_I6z_Pz = I_TWOBODYOVERLAP_K7z_S+ABZ*I_TWOBODYOVERLAP_I6z_S;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_D
   * expanding position: BRA2
   * code section is: HRR
   * totally 66 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_P
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P
   ************************************************************/
  Double I_TWOBODYOVERLAP_H5x_D2x = I_TWOBODYOVERLAP_I6x_Px+ABX*I_TWOBODYOVERLAP_H5x_Px;
  Double I_TWOBODYOVERLAP_H4xy_D2x = I_TWOBODYOVERLAP_I5xy_Px+ABX*I_TWOBODYOVERLAP_H4xy_Px;
  Double I_TWOBODYOVERLAP_H4xz_D2x = I_TWOBODYOVERLAP_I5xz_Px+ABX*I_TWOBODYOVERLAP_H4xz_Px;
  Double I_TWOBODYOVERLAP_H3x2y_D2x = I_TWOBODYOVERLAP_I4x2y_Px+ABX*I_TWOBODYOVERLAP_H3x2y_Px;
  Double I_TWOBODYOVERLAP_H3xyz_D2x = I_TWOBODYOVERLAP_I4xyz_Px+ABX*I_TWOBODYOVERLAP_H3xyz_Px;
  Double I_TWOBODYOVERLAP_H3x2z_D2x = I_TWOBODYOVERLAP_I4x2z_Px+ABX*I_TWOBODYOVERLAP_H3x2z_Px;
  Double I_TWOBODYOVERLAP_H2x3y_D2x = I_TWOBODYOVERLAP_I3x3y_Px+ABX*I_TWOBODYOVERLAP_H2x3y_Px;
  Double I_TWOBODYOVERLAP_H2x2yz_D2x = I_TWOBODYOVERLAP_I3x2yz_Px+ABX*I_TWOBODYOVERLAP_H2x2yz_Px;
  Double I_TWOBODYOVERLAP_H2xy2z_D2x = I_TWOBODYOVERLAP_I3xy2z_Px+ABX*I_TWOBODYOVERLAP_H2xy2z_Px;
  Double I_TWOBODYOVERLAP_H2x3z_D2x = I_TWOBODYOVERLAP_I3x3z_Px+ABX*I_TWOBODYOVERLAP_H2x3z_Px;
  Double I_TWOBODYOVERLAP_Hx4y_D2x = I_TWOBODYOVERLAP_I2x4y_Px+ABX*I_TWOBODYOVERLAP_Hx4y_Px;
  Double I_TWOBODYOVERLAP_Hx3yz_D2x = I_TWOBODYOVERLAP_I2x3yz_Px+ABX*I_TWOBODYOVERLAP_Hx3yz_Px;
  Double I_TWOBODYOVERLAP_Hx2y2z_D2x = I_TWOBODYOVERLAP_I2x2y2z_Px+ABX*I_TWOBODYOVERLAP_Hx2y2z_Px;
  Double I_TWOBODYOVERLAP_Hxy3z_D2x = I_TWOBODYOVERLAP_I2xy3z_Px+ABX*I_TWOBODYOVERLAP_Hxy3z_Px;
  Double I_TWOBODYOVERLAP_Hx4z_D2x = I_TWOBODYOVERLAP_I2x4z_Px+ABX*I_TWOBODYOVERLAP_Hx4z_Px;
  Double I_TWOBODYOVERLAP_H5y_D2x = I_TWOBODYOVERLAP_Ix5y_Px+ABX*I_TWOBODYOVERLAP_H5y_Px;
  Double I_TWOBODYOVERLAP_H4yz_D2x = I_TWOBODYOVERLAP_Ix4yz_Px+ABX*I_TWOBODYOVERLAP_H4yz_Px;
  Double I_TWOBODYOVERLAP_H3y2z_D2x = I_TWOBODYOVERLAP_Ix3y2z_Px+ABX*I_TWOBODYOVERLAP_H3y2z_Px;
  Double I_TWOBODYOVERLAP_H2y3z_D2x = I_TWOBODYOVERLAP_Ix2y3z_Px+ABX*I_TWOBODYOVERLAP_H2y3z_Px;
  Double I_TWOBODYOVERLAP_Hy4z_D2x = I_TWOBODYOVERLAP_Ixy4z_Px+ABX*I_TWOBODYOVERLAP_Hy4z_Px;
  Double I_TWOBODYOVERLAP_H5z_D2x = I_TWOBODYOVERLAP_Ix5z_Px+ABX*I_TWOBODYOVERLAP_H5z_Px;
  Double I_TWOBODYOVERLAP_H4xy_D2y = I_TWOBODYOVERLAP_I4x2y_Py+ABY*I_TWOBODYOVERLAP_H4xy_Py;
  Double I_TWOBODYOVERLAP_H4xz_D2y = I_TWOBODYOVERLAP_I4xyz_Py+ABY*I_TWOBODYOVERLAP_H4xz_Py;
  Double I_TWOBODYOVERLAP_H3x2y_D2y = I_TWOBODYOVERLAP_I3x3y_Py+ABY*I_TWOBODYOVERLAP_H3x2y_Py;
  Double I_TWOBODYOVERLAP_H3xyz_D2y = I_TWOBODYOVERLAP_I3x2yz_Py+ABY*I_TWOBODYOVERLAP_H3xyz_Py;
  Double I_TWOBODYOVERLAP_H3x2z_D2y = I_TWOBODYOVERLAP_I3xy2z_Py+ABY*I_TWOBODYOVERLAP_H3x2z_Py;
  Double I_TWOBODYOVERLAP_H2x3y_D2y = I_TWOBODYOVERLAP_I2x4y_Py+ABY*I_TWOBODYOVERLAP_H2x3y_Py;
  Double I_TWOBODYOVERLAP_H2x2yz_D2y = I_TWOBODYOVERLAP_I2x3yz_Py+ABY*I_TWOBODYOVERLAP_H2x2yz_Py;
  Double I_TWOBODYOVERLAP_H2xy2z_D2y = I_TWOBODYOVERLAP_I2x2y2z_Py+ABY*I_TWOBODYOVERLAP_H2xy2z_Py;
  Double I_TWOBODYOVERLAP_H2x3z_D2y = I_TWOBODYOVERLAP_I2xy3z_Py+ABY*I_TWOBODYOVERLAP_H2x3z_Py;
  Double I_TWOBODYOVERLAP_Hx4y_D2y = I_TWOBODYOVERLAP_Ix5y_Py+ABY*I_TWOBODYOVERLAP_Hx4y_Py;
  Double I_TWOBODYOVERLAP_Hx3yz_D2y = I_TWOBODYOVERLAP_Ix4yz_Py+ABY*I_TWOBODYOVERLAP_Hx3yz_Py;
  Double I_TWOBODYOVERLAP_Hx2y2z_D2y = I_TWOBODYOVERLAP_Ix3y2z_Py+ABY*I_TWOBODYOVERLAP_Hx2y2z_Py;
  Double I_TWOBODYOVERLAP_Hxy3z_D2y = I_TWOBODYOVERLAP_Ix2y3z_Py+ABY*I_TWOBODYOVERLAP_Hxy3z_Py;
  Double I_TWOBODYOVERLAP_Hx4z_D2y = I_TWOBODYOVERLAP_Ixy4z_Py+ABY*I_TWOBODYOVERLAP_Hx4z_Py;
  Double I_TWOBODYOVERLAP_H5y_D2y = I_TWOBODYOVERLAP_I6y_Py+ABY*I_TWOBODYOVERLAP_H5y_Py;
  Double I_TWOBODYOVERLAP_H4yz_D2y = I_TWOBODYOVERLAP_I5yz_Py+ABY*I_TWOBODYOVERLAP_H4yz_Py;
  Double I_TWOBODYOVERLAP_H3y2z_D2y = I_TWOBODYOVERLAP_I4y2z_Py+ABY*I_TWOBODYOVERLAP_H3y2z_Py;
  Double I_TWOBODYOVERLAP_H2y3z_D2y = I_TWOBODYOVERLAP_I3y3z_Py+ABY*I_TWOBODYOVERLAP_H2y3z_Py;
  Double I_TWOBODYOVERLAP_Hy4z_D2y = I_TWOBODYOVERLAP_I2y4z_Py+ABY*I_TWOBODYOVERLAP_Hy4z_Py;
  Double I_TWOBODYOVERLAP_H5z_D2y = I_TWOBODYOVERLAP_Iy5z_Py+ABY*I_TWOBODYOVERLAP_H5z_Py;
  Double I_TWOBODYOVERLAP_H4xy_D2z = I_TWOBODYOVERLAP_I4xyz_Pz+ABZ*I_TWOBODYOVERLAP_H4xy_Pz;
  Double I_TWOBODYOVERLAP_H4xz_D2z = I_TWOBODYOVERLAP_I4x2z_Pz+ABZ*I_TWOBODYOVERLAP_H4xz_Pz;
  Double I_TWOBODYOVERLAP_H3x2y_D2z = I_TWOBODYOVERLAP_I3x2yz_Pz+ABZ*I_TWOBODYOVERLAP_H3x2y_Pz;
  Double I_TWOBODYOVERLAP_H3xyz_D2z = I_TWOBODYOVERLAP_I3xy2z_Pz+ABZ*I_TWOBODYOVERLAP_H3xyz_Pz;
  Double I_TWOBODYOVERLAP_H3x2z_D2z = I_TWOBODYOVERLAP_I3x3z_Pz+ABZ*I_TWOBODYOVERLAP_H3x2z_Pz;
  Double I_TWOBODYOVERLAP_H2x3y_D2z = I_TWOBODYOVERLAP_I2x3yz_Pz+ABZ*I_TWOBODYOVERLAP_H2x3y_Pz;
  Double I_TWOBODYOVERLAP_H2x2yz_D2z = I_TWOBODYOVERLAP_I2x2y2z_Pz+ABZ*I_TWOBODYOVERLAP_H2x2yz_Pz;
  Double I_TWOBODYOVERLAP_H2xy2z_D2z = I_TWOBODYOVERLAP_I2xy3z_Pz+ABZ*I_TWOBODYOVERLAP_H2xy2z_Pz;
  Double I_TWOBODYOVERLAP_H2x3z_D2z = I_TWOBODYOVERLAP_I2x4z_Pz+ABZ*I_TWOBODYOVERLAP_H2x3z_Pz;
  Double I_TWOBODYOVERLAP_Hx4y_D2z = I_TWOBODYOVERLAP_Ix4yz_Pz+ABZ*I_TWOBODYOVERLAP_Hx4y_Pz;
  Double I_TWOBODYOVERLAP_Hx3yz_D2z = I_TWOBODYOVERLAP_Ix3y2z_Pz+ABZ*I_TWOBODYOVERLAP_Hx3yz_Pz;
  Double I_TWOBODYOVERLAP_Hx2y2z_D2z = I_TWOBODYOVERLAP_Ix2y3z_Pz+ABZ*I_TWOBODYOVERLAP_Hx2y2z_Pz;
  Double I_TWOBODYOVERLAP_Hxy3z_D2z = I_TWOBODYOVERLAP_Ixy4z_Pz+ABZ*I_TWOBODYOVERLAP_Hxy3z_Pz;
  Double I_TWOBODYOVERLAP_Hx4z_D2z = I_TWOBODYOVERLAP_Ix5z_Pz+ABZ*I_TWOBODYOVERLAP_Hx4z_Pz;
  Double I_TWOBODYOVERLAP_H4yz_D2z = I_TWOBODYOVERLAP_I4y2z_Pz+ABZ*I_TWOBODYOVERLAP_H4yz_Pz;
  Double I_TWOBODYOVERLAP_H3y2z_D2z = I_TWOBODYOVERLAP_I3y3z_Pz+ABZ*I_TWOBODYOVERLAP_H3y2z_Pz;
  Double I_TWOBODYOVERLAP_H2y3z_D2z = I_TWOBODYOVERLAP_I2y4z_Pz+ABZ*I_TWOBODYOVERLAP_H2y3z_Pz;
  Double I_TWOBODYOVERLAP_Hy4z_D2z = I_TWOBODYOVERLAP_Iy5z_Pz+ABZ*I_TWOBODYOVERLAP_Hy4z_Pz;
  Double I_TWOBODYOVERLAP_H5z_D2z = I_TWOBODYOVERLAP_I6z_Pz+ABZ*I_TWOBODYOVERLAP_H5z_Pz;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_F
   * expanding position: BRA2
   * code section is: HRR
   * totally 51 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_D
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D
   ************************************************************/
  Double I_TWOBODYOVERLAP_G4x_F3x = I_TWOBODYOVERLAP_H5x_D2x+ABX*I_TWOBODYOVERLAP_G4x_D2x;
  Double I_TWOBODYOVERLAP_G3xy_F3x = I_TWOBODYOVERLAP_H4xy_D2x+ABX*I_TWOBODYOVERLAP_G3xy_D2x;
  Double I_TWOBODYOVERLAP_G3xz_F3x = I_TWOBODYOVERLAP_H4xz_D2x+ABX*I_TWOBODYOVERLAP_G3xz_D2x;
  Double I_TWOBODYOVERLAP_G2x2y_F3x = I_TWOBODYOVERLAP_H3x2y_D2x+ABX*I_TWOBODYOVERLAP_G2x2y_D2x;
  Double I_TWOBODYOVERLAP_G2xyz_F3x = I_TWOBODYOVERLAP_H3xyz_D2x+ABX*I_TWOBODYOVERLAP_G2xyz_D2x;
  Double I_TWOBODYOVERLAP_G2x2z_F3x = I_TWOBODYOVERLAP_H3x2z_D2x+ABX*I_TWOBODYOVERLAP_G2x2z_D2x;
  Double I_TWOBODYOVERLAP_Gx3y_F3x = I_TWOBODYOVERLAP_H2x3y_D2x+ABX*I_TWOBODYOVERLAP_Gx3y_D2x;
  Double I_TWOBODYOVERLAP_Gx2yz_F3x = I_TWOBODYOVERLAP_H2x2yz_D2x+ABX*I_TWOBODYOVERLAP_Gx2yz_D2x;
  Double I_TWOBODYOVERLAP_Gxy2z_F3x = I_TWOBODYOVERLAP_H2xy2z_D2x+ABX*I_TWOBODYOVERLAP_Gxy2z_D2x;
  Double I_TWOBODYOVERLAP_Gx3z_F3x = I_TWOBODYOVERLAP_H2x3z_D2x+ABX*I_TWOBODYOVERLAP_Gx3z_D2x;
  Double I_TWOBODYOVERLAP_G4y_F3x = I_TWOBODYOVERLAP_Hx4y_D2x+ABX*I_TWOBODYOVERLAP_G4y_D2x;
  Double I_TWOBODYOVERLAP_G3yz_F3x = I_TWOBODYOVERLAP_Hx3yz_D2x+ABX*I_TWOBODYOVERLAP_G3yz_D2x;
  Double I_TWOBODYOVERLAP_G2y2z_F3x = I_TWOBODYOVERLAP_Hx2y2z_D2x+ABX*I_TWOBODYOVERLAP_G2y2z_D2x;
  Double I_TWOBODYOVERLAP_Gy3z_F3x = I_TWOBODYOVERLAP_Hxy3z_D2x+ABX*I_TWOBODYOVERLAP_Gy3z_D2x;
  Double I_TWOBODYOVERLAP_G4z_F3x = I_TWOBODYOVERLAP_Hx4z_D2x+ABX*I_TWOBODYOVERLAP_G4z_D2x;
  Double I_TWOBODYOVERLAP_G3xy_F2xy = I_TWOBODYOVERLAP_H3x2y_D2x+ABY*I_TWOBODYOVERLAP_G3xy_D2x;
  Double I_TWOBODYOVERLAP_G3xz_F2xy = I_TWOBODYOVERLAP_H3xyz_D2x+ABY*I_TWOBODYOVERLAP_G3xz_D2x;
  Double I_TWOBODYOVERLAP_G2x2y_F2xy = I_TWOBODYOVERLAP_H2x3y_D2x+ABY*I_TWOBODYOVERLAP_G2x2y_D2x;
  Double I_TWOBODYOVERLAP_G2xyz_F2xy = I_TWOBODYOVERLAP_H2x2yz_D2x+ABY*I_TWOBODYOVERLAP_G2xyz_D2x;
  Double I_TWOBODYOVERLAP_G2x2z_F2xy = I_TWOBODYOVERLAP_H2xy2z_D2x+ABY*I_TWOBODYOVERLAP_G2x2z_D2x;
  Double I_TWOBODYOVERLAP_Gx3y_F2xy = I_TWOBODYOVERLAP_Hx4y_D2x+ABY*I_TWOBODYOVERLAP_Gx3y_D2x;
  Double I_TWOBODYOVERLAP_Gx2yz_F2xy = I_TWOBODYOVERLAP_Hx3yz_D2x+ABY*I_TWOBODYOVERLAP_Gx2yz_D2x;
  Double I_TWOBODYOVERLAP_Gxy2z_F2xy = I_TWOBODYOVERLAP_Hx2y2z_D2x+ABY*I_TWOBODYOVERLAP_Gxy2z_D2x;
  Double I_TWOBODYOVERLAP_Gx3z_F2xy = I_TWOBODYOVERLAP_Hxy3z_D2x+ABY*I_TWOBODYOVERLAP_Gx3z_D2x;
  Double I_TWOBODYOVERLAP_G4y_F2xy = I_TWOBODYOVERLAP_H5y_D2x+ABY*I_TWOBODYOVERLAP_G4y_D2x;
  Double I_TWOBODYOVERLAP_G3yz_F2xy = I_TWOBODYOVERLAP_H4yz_D2x+ABY*I_TWOBODYOVERLAP_G3yz_D2x;
  Double I_TWOBODYOVERLAP_G2y2z_F2xy = I_TWOBODYOVERLAP_H3y2z_D2x+ABY*I_TWOBODYOVERLAP_G2y2z_D2x;
  Double I_TWOBODYOVERLAP_Gy3z_F2xy = I_TWOBODYOVERLAP_H2y3z_D2x+ABY*I_TWOBODYOVERLAP_Gy3z_D2x;
  Double I_TWOBODYOVERLAP_G4z_F2xy = I_TWOBODYOVERLAP_Hy4z_D2x+ABY*I_TWOBODYOVERLAP_G4z_D2x;
  Double I_TWOBODYOVERLAP_G3xz_F2xz = I_TWOBODYOVERLAP_H3x2z_D2x+ABZ*I_TWOBODYOVERLAP_G3xz_D2x;
  Double I_TWOBODYOVERLAP_G2xyz_F2xz = I_TWOBODYOVERLAP_H2xy2z_D2x+ABZ*I_TWOBODYOVERLAP_G2xyz_D2x;
  Double I_TWOBODYOVERLAP_G2x2z_F2xz = I_TWOBODYOVERLAP_H2x3z_D2x+ABZ*I_TWOBODYOVERLAP_G2x2z_D2x;
  Double I_TWOBODYOVERLAP_Gx2yz_F2xz = I_TWOBODYOVERLAP_Hx2y2z_D2x+ABZ*I_TWOBODYOVERLAP_Gx2yz_D2x;
  Double I_TWOBODYOVERLAP_Gxy2z_F2xz = I_TWOBODYOVERLAP_Hxy3z_D2x+ABZ*I_TWOBODYOVERLAP_Gxy2z_D2x;
  Double I_TWOBODYOVERLAP_Gx3z_F2xz = I_TWOBODYOVERLAP_Hx4z_D2x+ABZ*I_TWOBODYOVERLAP_Gx3z_D2x;
  Double I_TWOBODYOVERLAP_G3yz_F2xz = I_TWOBODYOVERLAP_H3y2z_D2x+ABZ*I_TWOBODYOVERLAP_G3yz_D2x;
  Double I_TWOBODYOVERLAP_G2y2z_F2xz = I_TWOBODYOVERLAP_H2y3z_D2x+ABZ*I_TWOBODYOVERLAP_G2y2z_D2x;
  Double I_TWOBODYOVERLAP_Gy3z_F2xz = I_TWOBODYOVERLAP_Hy4z_D2x+ABZ*I_TWOBODYOVERLAP_Gy3z_D2x;
  Double I_TWOBODYOVERLAP_G4z_F2xz = I_TWOBODYOVERLAP_H5z_D2x+ABZ*I_TWOBODYOVERLAP_G4z_D2x;
  Double I_TWOBODYOVERLAP_G3xz_Fx2y = I_TWOBODYOVERLAP_H4xz_D2y+ABX*I_TWOBODYOVERLAP_G3xz_D2y;
  Double I_TWOBODYOVERLAP_G2xyz_Fx2y = I_TWOBODYOVERLAP_H3xyz_D2y+ABX*I_TWOBODYOVERLAP_G2xyz_D2y;
  Double I_TWOBODYOVERLAP_G2x2z_Fx2y = I_TWOBODYOVERLAP_H3x2z_D2y+ABX*I_TWOBODYOVERLAP_G2x2z_D2y;
  Double I_TWOBODYOVERLAP_Gx2yz_Fx2y = I_TWOBODYOVERLAP_H2x2yz_D2y+ABX*I_TWOBODYOVERLAP_Gx2yz_D2y;
  Double I_TWOBODYOVERLAP_Gxy2z_Fx2y = I_TWOBODYOVERLAP_H2xy2z_D2y+ABX*I_TWOBODYOVERLAP_Gxy2z_D2y;
  Double I_TWOBODYOVERLAP_Gx3z_Fx2y = I_TWOBODYOVERLAP_H2x3z_D2y+ABX*I_TWOBODYOVERLAP_Gx3z_D2y;
  Double I_TWOBODYOVERLAP_G3yz_Fx2y = I_TWOBODYOVERLAP_Hx3yz_D2y+ABX*I_TWOBODYOVERLAP_G3yz_D2y;
  Double I_TWOBODYOVERLAP_G2y2z_Fx2y = I_TWOBODYOVERLAP_Hx2y2z_D2y+ABX*I_TWOBODYOVERLAP_G2y2z_D2y;
  Double I_TWOBODYOVERLAP_Gy3z_Fx2y = I_TWOBODYOVERLAP_Hxy3z_D2y+ABX*I_TWOBODYOVERLAP_Gy3z_D2y;
  Double I_TWOBODYOVERLAP_G4z_Fx2y = I_TWOBODYOVERLAP_Hx4z_D2y+ABX*I_TWOBODYOVERLAP_G4z_D2y;
  Double I_TWOBODYOVERLAP_G3xy_Fx2z = I_TWOBODYOVERLAP_H4xy_D2z+ABX*I_TWOBODYOVERLAP_G3xy_D2z;
  Double I_TWOBODYOVERLAP_G2x2y_Fx2z = I_TWOBODYOVERLAP_H3x2y_D2z+ABX*I_TWOBODYOVERLAP_G2x2y_D2z;
  Double I_TWOBODYOVERLAP_G2xyz_Fx2z = I_TWOBODYOVERLAP_H3xyz_D2z+ABX*I_TWOBODYOVERLAP_G2xyz_D2z;
  Double I_TWOBODYOVERLAP_Gx3y_Fx2z = I_TWOBODYOVERLAP_H2x3y_D2z+ABX*I_TWOBODYOVERLAP_Gx3y_D2z;
  Double I_TWOBODYOVERLAP_Gx2yz_Fx2z = I_TWOBODYOVERLAP_H2x2yz_D2z+ABX*I_TWOBODYOVERLAP_Gx2yz_D2z;
  Double I_TWOBODYOVERLAP_Gxy2z_Fx2z = I_TWOBODYOVERLAP_H2xy2z_D2z+ABX*I_TWOBODYOVERLAP_Gxy2z_D2z;
  Double I_TWOBODYOVERLAP_G4y_Fx2z = I_TWOBODYOVERLAP_Hx4y_D2z+ABX*I_TWOBODYOVERLAP_G4y_D2z;
  Double I_TWOBODYOVERLAP_G3yz_Fx2z = I_TWOBODYOVERLAP_Hx3yz_D2z+ABX*I_TWOBODYOVERLAP_G3yz_D2z;
  Double I_TWOBODYOVERLAP_G2y2z_Fx2z = I_TWOBODYOVERLAP_Hx2y2z_D2z+ABX*I_TWOBODYOVERLAP_G2y2z_D2z;
  Double I_TWOBODYOVERLAP_Gy3z_Fx2z = I_TWOBODYOVERLAP_Hxy3z_D2z+ABX*I_TWOBODYOVERLAP_Gy3z_D2z;
  Double I_TWOBODYOVERLAP_G4x_F3y = I_TWOBODYOVERLAP_H4xy_D2y+ABY*I_TWOBODYOVERLAP_G4x_D2y;
  Double I_TWOBODYOVERLAP_G3xy_F3y = I_TWOBODYOVERLAP_H3x2y_D2y+ABY*I_TWOBODYOVERLAP_G3xy_D2y;
  Double I_TWOBODYOVERLAP_G3xz_F3y = I_TWOBODYOVERLAP_H3xyz_D2y+ABY*I_TWOBODYOVERLAP_G3xz_D2y;
  Double I_TWOBODYOVERLAP_G2x2y_F3y = I_TWOBODYOVERLAP_H2x3y_D2y+ABY*I_TWOBODYOVERLAP_G2x2y_D2y;
  Double I_TWOBODYOVERLAP_G2xyz_F3y = I_TWOBODYOVERLAP_H2x2yz_D2y+ABY*I_TWOBODYOVERLAP_G2xyz_D2y;
  Double I_TWOBODYOVERLAP_G2x2z_F3y = I_TWOBODYOVERLAP_H2xy2z_D2y+ABY*I_TWOBODYOVERLAP_G2x2z_D2y;
  Double I_TWOBODYOVERLAP_Gx3y_F3y = I_TWOBODYOVERLAP_Hx4y_D2y+ABY*I_TWOBODYOVERLAP_Gx3y_D2y;
  Double I_TWOBODYOVERLAP_Gx2yz_F3y = I_TWOBODYOVERLAP_Hx3yz_D2y+ABY*I_TWOBODYOVERLAP_Gx2yz_D2y;
  Double I_TWOBODYOVERLAP_Gxy2z_F3y = I_TWOBODYOVERLAP_Hx2y2z_D2y+ABY*I_TWOBODYOVERLAP_Gxy2z_D2y;
  Double I_TWOBODYOVERLAP_Gx3z_F3y = I_TWOBODYOVERLAP_Hxy3z_D2y+ABY*I_TWOBODYOVERLAP_Gx3z_D2y;
  Double I_TWOBODYOVERLAP_G4y_F3y = I_TWOBODYOVERLAP_H5y_D2y+ABY*I_TWOBODYOVERLAP_G4y_D2y;
  Double I_TWOBODYOVERLAP_G3yz_F3y = I_TWOBODYOVERLAP_H4yz_D2y+ABY*I_TWOBODYOVERLAP_G3yz_D2y;
  Double I_TWOBODYOVERLAP_G2y2z_F3y = I_TWOBODYOVERLAP_H3y2z_D2y+ABY*I_TWOBODYOVERLAP_G2y2z_D2y;
  Double I_TWOBODYOVERLAP_Gy3z_F3y = I_TWOBODYOVERLAP_H2y3z_D2y+ABY*I_TWOBODYOVERLAP_Gy3z_D2y;
  Double I_TWOBODYOVERLAP_G4z_F3y = I_TWOBODYOVERLAP_Hy4z_D2y+ABY*I_TWOBODYOVERLAP_G4z_D2y;
  Double I_TWOBODYOVERLAP_G3xz_F2yz = I_TWOBODYOVERLAP_H3x2z_D2y+ABZ*I_TWOBODYOVERLAP_G3xz_D2y;
  Double I_TWOBODYOVERLAP_G2xyz_F2yz = I_TWOBODYOVERLAP_H2xy2z_D2y+ABZ*I_TWOBODYOVERLAP_G2xyz_D2y;
  Double I_TWOBODYOVERLAP_G2x2z_F2yz = I_TWOBODYOVERLAP_H2x3z_D2y+ABZ*I_TWOBODYOVERLAP_G2x2z_D2y;
  Double I_TWOBODYOVERLAP_Gx2yz_F2yz = I_TWOBODYOVERLAP_Hx2y2z_D2y+ABZ*I_TWOBODYOVERLAP_Gx2yz_D2y;
  Double I_TWOBODYOVERLAP_Gxy2z_F2yz = I_TWOBODYOVERLAP_Hxy3z_D2y+ABZ*I_TWOBODYOVERLAP_Gxy2z_D2y;
  Double I_TWOBODYOVERLAP_Gx3z_F2yz = I_TWOBODYOVERLAP_Hx4z_D2y+ABZ*I_TWOBODYOVERLAP_Gx3z_D2y;
  Double I_TWOBODYOVERLAP_G3yz_F2yz = I_TWOBODYOVERLAP_H3y2z_D2y+ABZ*I_TWOBODYOVERLAP_G3yz_D2y;
  Double I_TWOBODYOVERLAP_G2y2z_F2yz = I_TWOBODYOVERLAP_H2y3z_D2y+ABZ*I_TWOBODYOVERLAP_G2y2z_D2y;
  Double I_TWOBODYOVERLAP_Gy3z_F2yz = I_TWOBODYOVERLAP_Hy4z_D2y+ABZ*I_TWOBODYOVERLAP_Gy3z_D2y;
  Double I_TWOBODYOVERLAP_G4z_F2yz = I_TWOBODYOVERLAP_H5z_D2y+ABZ*I_TWOBODYOVERLAP_G4z_D2y;
  Double I_TWOBODYOVERLAP_G4x_F3z = I_TWOBODYOVERLAP_H4xz_D2z+ABZ*I_TWOBODYOVERLAP_G4x_D2z;
  Double I_TWOBODYOVERLAP_G3xy_F3z = I_TWOBODYOVERLAP_H3xyz_D2z+ABZ*I_TWOBODYOVERLAP_G3xy_D2z;
  Double I_TWOBODYOVERLAP_G3xz_F3z = I_TWOBODYOVERLAP_H3x2z_D2z+ABZ*I_TWOBODYOVERLAP_G3xz_D2z;
  Double I_TWOBODYOVERLAP_G2x2y_F3z = I_TWOBODYOVERLAP_H2x2yz_D2z+ABZ*I_TWOBODYOVERLAP_G2x2y_D2z;
  Double I_TWOBODYOVERLAP_G2xyz_F3z = I_TWOBODYOVERLAP_H2xy2z_D2z+ABZ*I_TWOBODYOVERLAP_G2xyz_D2z;
  Double I_TWOBODYOVERLAP_G2x2z_F3z = I_TWOBODYOVERLAP_H2x3z_D2z+ABZ*I_TWOBODYOVERLAP_G2x2z_D2z;
  Double I_TWOBODYOVERLAP_Gx3y_F3z = I_TWOBODYOVERLAP_Hx3yz_D2z+ABZ*I_TWOBODYOVERLAP_Gx3y_D2z;
  Double I_TWOBODYOVERLAP_Gx2yz_F3z = I_TWOBODYOVERLAP_Hx2y2z_D2z+ABZ*I_TWOBODYOVERLAP_Gx2yz_D2z;
  Double I_TWOBODYOVERLAP_Gxy2z_F3z = I_TWOBODYOVERLAP_Hxy3z_D2z+ABZ*I_TWOBODYOVERLAP_Gxy2z_D2z;
  Double I_TWOBODYOVERLAP_Gx3z_F3z = I_TWOBODYOVERLAP_Hx4z_D2z+ABZ*I_TWOBODYOVERLAP_Gx3z_D2z;
  Double I_TWOBODYOVERLAP_G4y_F3z = I_TWOBODYOVERLAP_H4yz_D2z+ABZ*I_TWOBODYOVERLAP_G4y_D2z;
  Double I_TWOBODYOVERLAP_G3yz_F3z = I_TWOBODYOVERLAP_H3y2z_D2z+ABZ*I_TWOBODYOVERLAP_G3yz_D2z;
  Double I_TWOBODYOVERLAP_G2y2z_F3z = I_TWOBODYOVERLAP_H2y3z_D2z+ABZ*I_TWOBODYOVERLAP_G2y2z_D2z;
  Double I_TWOBODYOVERLAP_Gy3z_F3z = I_TWOBODYOVERLAP_Hy4z_D2z+ABZ*I_TWOBODYOVERLAP_Gy3z_D2z;
  Double I_TWOBODYOVERLAP_G4z_F3z = I_TWOBODYOVERLAP_H5z_D2z+ABZ*I_TWOBODYOVERLAP_G4z_D2z;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_F_G
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_F
   ************************************************************/
  Double I_TWOBODYOVERLAP_F3x_G4x = I_TWOBODYOVERLAP_G4x_F3x+ABX*I_TWOBODYOVERLAP_F3x_F3x;
  Double I_TWOBODYOVERLAP_F2xy_G4x = I_TWOBODYOVERLAP_G3xy_F3x+ABX*I_TWOBODYOVERLAP_F2xy_F3x;
  Double I_TWOBODYOVERLAP_F2xz_G4x = I_TWOBODYOVERLAP_G3xz_F3x+ABX*I_TWOBODYOVERLAP_F2xz_F3x;
  Double I_TWOBODYOVERLAP_Fx2y_G4x = I_TWOBODYOVERLAP_G2x2y_F3x+ABX*I_TWOBODYOVERLAP_Fx2y_F3x;
  Double I_TWOBODYOVERLAP_Fxyz_G4x = I_TWOBODYOVERLAP_G2xyz_F3x+ABX*I_TWOBODYOVERLAP_Fxyz_F3x;
  Double I_TWOBODYOVERLAP_Fx2z_G4x = I_TWOBODYOVERLAP_G2x2z_F3x+ABX*I_TWOBODYOVERLAP_Fx2z_F3x;
  Double I_TWOBODYOVERLAP_F3y_G4x = I_TWOBODYOVERLAP_Gx3y_F3x+ABX*I_TWOBODYOVERLAP_F3y_F3x;
  Double I_TWOBODYOVERLAP_F2yz_G4x = I_TWOBODYOVERLAP_Gx2yz_F3x+ABX*I_TWOBODYOVERLAP_F2yz_F3x;
  Double I_TWOBODYOVERLAP_Fy2z_G4x = I_TWOBODYOVERLAP_Gxy2z_F3x+ABX*I_TWOBODYOVERLAP_Fy2z_F3x;
  Double I_TWOBODYOVERLAP_F3z_G4x = I_TWOBODYOVERLAP_Gx3z_F3x+ABX*I_TWOBODYOVERLAP_F3z_F3x;
  Double I_TWOBODYOVERLAP_F3x_G3xy = I_TWOBODYOVERLAP_G3xy_F3x+ABY*I_TWOBODYOVERLAP_F3x_F3x;
  Double I_TWOBODYOVERLAP_F2xy_G3xy = I_TWOBODYOVERLAP_G2x2y_F3x+ABY*I_TWOBODYOVERLAP_F2xy_F3x;
  Double I_TWOBODYOVERLAP_F2xz_G3xy = I_TWOBODYOVERLAP_G2xyz_F3x+ABY*I_TWOBODYOVERLAP_F2xz_F3x;
  Double I_TWOBODYOVERLAP_Fx2y_G3xy = I_TWOBODYOVERLAP_Gx3y_F3x+ABY*I_TWOBODYOVERLAP_Fx2y_F3x;
  Double I_TWOBODYOVERLAP_Fxyz_G3xy = I_TWOBODYOVERLAP_Gx2yz_F3x+ABY*I_TWOBODYOVERLAP_Fxyz_F3x;
  Double I_TWOBODYOVERLAP_Fx2z_G3xy = I_TWOBODYOVERLAP_Gxy2z_F3x+ABY*I_TWOBODYOVERLAP_Fx2z_F3x;
  Double I_TWOBODYOVERLAP_F3y_G3xy = I_TWOBODYOVERLAP_G4y_F3x+ABY*I_TWOBODYOVERLAP_F3y_F3x;
  Double I_TWOBODYOVERLAP_F2yz_G3xy = I_TWOBODYOVERLAP_G3yz_F3x+ABY*I_TWOBODYOVERLAP_F2yz_F3x;
  Double I_TWOBODYOVERLAP_Fy2z_G3xy = I_TWOBODYOVERLAP_G2y2z_F3x+ABY*I_TWOBODYOVERLAP_Fy2z_F3x;
  Double I_TWOBODYOVERLAP_F3z_G3xy = I_TWOBODYOVERLAP_Gy3z_F3x+ABY*I_TWOBODYOVERLAP_F3z_F3x;
  Double I_TWOBODYOVERLAP_F3x_G3xz = I_TWOBODYOVERLAP_G3xz_F3x+ABZ*I_TWOBODYOVERLAP_F3x_F3x;
  Double I_TWOBODYOVERLAP_F2xy_G3xz = I_TWOBODYOVERLAP_G2xyz_F3x+ABZ*I_TWOBODYOVERLAP_F2xy_F3x;
  Double I_TWOBODYOVERLAP_F2xz_G3xz = I_TWOBODYOVERLAP_G2x2z_F3x+ABZ*I_TWOBODYOVERLAP_F2xz_F3x;
  Double I_TWOBODYOVERLAP_Fx2y_G3xz = I_TWOBODYOVERLAP_Gx2yz_F3x+ABZ*I_TWOBODYOVERLAP_Fx2y_F3x;
  Double I_TWOBODYOVERLAP_Fxyz_G3xz = I_TWOBODYOVERLAP_Gxy2z_F3x+ABZ*I_TWOBODYOVERLAP_Fxyz_F3x;
  Double I_TWOBODYOVERLAP_Fx2z_G3xz = I_TWOBODYOVERLAP_Gx3z_F3x+ABZ*I_TWOBODYOVERLAP_Fx2z_F3x;
  Double I_TWOBODYOVERLAP_F3y_G3xz = I_TWOBODYOVERLAP_G3yz_F3x+ABZ*I_TWOBODYOVERLAP_F3y_F3x;
  Double I_TWOBODYOVERLAP_F2yz_G3xz = I_TWOBODYOVERLAP_G2y2z_F3x+ABZ*I_TWOBODYOVERLAP_F2yz_F3x;
  Double I_TWOBODYOVERLAP_Fy2z_G3xz = I_TWOBODYOVERLAP_Gy3z_F3x+ABZ*I_TWOBODYOVERLAP_Fy2z_F3x;
  Double I_TWOBODYOVERLAP_F3z_G3xz = I_TWOBODYOVERLAP_G4z_F3x+ABZ*I_TWOBODYOVERLAP_F3z_F3x;
  Double I_TWOBODYOVERLAP_F3x_G2x2y = I_TWOBODYOVERLAP_G3xy_F2xy+ABY*I_TWOBODYOVERLAP_F3x_F2xy;
  Double I_TWOBODYOVERLAP_F2xy_G2x2y = I_TWOBODYOVERLAP_G2x2y_F2xy+ABY*I_TWOBODYOVERLAP_F2xy_F2xy;
  Double I_TWOBODYOVERLAP_F2xz_G2x2y = I_TWOBODYOVERLAP_G2xyz_F2xy+ABY*I_TWOBODYOVERLAP_F2xz_F2xy;
  Double I_TWOBODYOVERLAP_Fx2y_G2x2y = I_TWOBODYOVERLAP_Gx3y_F2xy+ABY*I_TWOBODYOVERLAP_Fx2y_F2xy;
  Double I_TWOBODYOVERLAP_Fxyz_G2x2y = I_TWOBODYOVERLAP_Gx2yz_F2xy+ABY*I_TWOBODYOVERLAP_Fxyz_F2xy;
  Double I_TWOBODYOVERLAP_Fx2z_G2x2y = I_TWOBODYOVERLAP_Gxy2z_F2xy+ABY*I_TWOBODYOVERLAP_Fx2z_F2xy;
  Double I_TWOBODYOVERLAP_F3y_G2x2y = I_TWOBODYOVERLAP_G4y_F2xy+ABY*I_TWOBODYOVERLAP_F3y_F2xy;
  Double I_TWOBODYOVERLAP_F2yz_G2x2y = I_TWOBODYOVERLAP_G3yz_F2xy+ABY*I_TWOBODYOVERLAP_F2yz_F2xy;
  Double I_TWOBODYOVERLAP_Fy2z_G2x2y = I_TWOBODYOVERLAP_G2y2z_F2xy+ABY*I_TWOBODYOVERLAP_Fy2z_F2xy;
  Double I_TWOBODYOVERLAP_F3z_G2x2y = I_TWOBODYOVERLAP_Gy3z_F2xy+ABY*I_TWOBODYOVERLAP_F3z_F2xy;
  Double I_TWOBODYOVERLAP_F3x_G2xyz = I_TWOBODYOVERLAP_G3xz_F2xy+ABZ*I_TWOBODYOVERLAP_F3x_F2xy;
  Double I_TWOBODYOVERLAP_F2xy_G2xyz = I_TWOBODYOVERLAP_G2xyz_F2xy+ABZ*I_TWOBODYOVERLAP_F2xy_F2xy;
  Double I_TWOBODYOVERLAP_F2xz_G2xyz = I_TWOBODYOVERLAP_G2x2z_F2xy+ABZ*I_TWOBODYOVERLAP_F2xz_F2xy;
  Double I_TWOBODYOVERLAP_Fx2y_G2xyz = I_TWOBODYOVERLAP_Gx2yz_F2xy+ABZ*I_TWOBODYOVERLAP_Fx2y_F2xy;
  Double I_TWOBODYOVERLAP_Fxyz_G2xyz = I_TWOBODYOVERLAP_Gxy2z_F2xy+ABZ*I_TWOBODYOVERLAP_Fxyz_F2xy;
  Double I_TWOBODYOVERLAP_Fx2z_G2xyz = I_TWOBODYOVERLAP_Gx3z_F2xy+ABZ*I_TWOBODYOVERLAP_Fx2z_F2xy;
  Double I_TWOBODYOVERLAP_F3y_G2xyz = I_TWOBODYOVERLAP_G3yz_F2xy+ABZ*I_TWOBODYOVERLAP_F3y_F2xy;
  Double I_TWOBODYOVERLAP_F2yz_G2xyz = I_TWOBODYOVERLAP_G2y2z_F2xy+ABZ*I_TWOBODYOVERLAP_F2yz_F2xy;
  Double I_TWOBODYOVERLAP_Fy2z_G2xyz = I_TWOBODYOVERLAP_Gy3z_F2xy+ABZ*I_TWOBODYOVERLAP_Fy2z_F2xy;
  Double I_TWOBODYOVERLAP_F3z_G2xyz = I_TWOBODYOVERLAP_G4z_F2xy+ABZ*I_TWOBODYOVERLAP_F3z_F2xy;
  Double I_TWOBODYOVERLAP_F3x_G2x2z = I_TWOBODYOVERLAP_G3xz_F2xz+ABZ*I_TWOBODYOVERLAP_F3x_F2xz;
  Double I_TWOBODYOVERLAP_F2xy_G2x2z = I_TWOBODYOVERLAP_G2xyz_F2xz+ABZ*I_TWOBODYOVERLAP_F2xy_F2xz;
  Double I_TWOBODYOVERLAP_F2xz_G2x2z = I_TWOBODYOVERLAP_G2x2z_F2xz+ABZ*I_TWOBODYOVERLAP_F2xz_F2xz;
  Double I_TWOBODYOVERLAP_Fx2y_G2x2z = I_TWOBODYOVERLAP_Gx2yz_F2xz+ABZ*I_TWOBODYOVERLAP_Fx2y_F2xz;
  Double I_TWOBODYOVERLAP_Fxyz_G2x2z = I_TWOBODYOVERLAP_Gxy2z_F2xz+ABZ*I_TWOBODYOVERLAP_Fxyz_F2xz;
  Double I_TWOBODYOVERLAP_Fx2z_G2x2z = I_TWOBODYOVERLAP_Gx3z_F2xz+ABZ*I_TWOBODYOVERLAP_Fx2z_F2xz;
  Double I_TWOBODYOVERLAP_F3y_G2x2z = I_TWOBODYOVERLAP_G3yz_F2xz+ABZ*I_TWOBODYOVERLAP_F3y_F2xz;
  Double I_TWOBODYOVERLAP_F2yz_G2x2z = I_TWOBODYOVERLAP_G2y2z_F2xz+ABZ*I_TWOBODYOVERLAP_F2yz_F2xz;
  Double I_TWOBODYOVERLAP_Fy2z_G2x2z = I_TWOBODYOVERLAP_Gy3z_F2xz+ABZ*I_TWOBODYOVERLAP_Fy2z_F2xz;
  Double I_TWOBODYOVERLAP_F3z_G2x2z = I_TWOBODYOVERLAP_G4z_F2xz+ABZ*I_TWOBODYOVERLAP_F3z_F2xz;
  Double I_TWOBODYOVERLAP_F3x_Gx3y = I_TWOBODYOVERLAP_G4x_F3y+ABX*I_TWOBODYOVERLAP_F3x_F3y;
  Double I_TWOBODYOVERLAP_F2xy_Gx3y = I_TWOBODYOVERLAP_G3xy_F3y+ABX*I_TWOBODYOVERLAP_F2xy_F3y;
  Double I_TWOBODYOVERLAP_F2xz_Gx3y = I_TWOBODYOVERLAP_G3xz_F3y+ABX*I_TWOBODYOVERLAP_F2xz_F3y;
  Double I_TWOBODYOVERLAP_Fx2y_Gx3y = I_TWOBODYOVERLAP_G2x2y_F3y+ABX*I_TWOBODYOVERLAP_Fx2y_F3y;
  Double I_TWOBODYOVERLAP_Fxyz_Gx3y = I_TWOBODYOVERLAP_G2xyz_F3y+ABX*I_TWOBODYOVERLAP_Fxyz_F3y;
  Double I_TWOBODYOVERLAP_Fx2z_Gx3y = I_TWOBODYOVERLAP_G2x2z_F3y+ABX*I_TWOBODYOVERLAP_Fx2z_F3y;
  Double I_TWOBODYOVERLAP_F3y_Gx3y = I_TWOBODYOVERLAP_Gx3y_F3y+ABX*I_TWOBODYOVERLAP_F3y_F3y;
  Double I_TWOBODYOVERLAP_F2yz_Gx3y = I_TWOBODYOVERLAP_Gx2yz_F3y+ABX*I_TWOBODYOVERLAP_F2yz_F3y;
  Double I_TWOBODYOVERLAP_Fy2z_Gx3y = I_TWOBODYOVERLAP_Gxy2z_F3y+ABX*I_TWOBODYOVERLAP_Fy2z_F3y;
  Double I_TWOBODYOVERLAP_F3z_Gx3y = I_TWOBODYOVERLAP_Gx3z_F3y+ABX*I_TWOBODYOVERLAP_F3z_F3y;
  Double I_TWOBODYOVERLAP_F3x_Gx2yz = I_TWOBODYOVERLAP_G3xz_Fx2y+ABZ*I_TWOBODYOVERLAP_F3x_Fx2y;
  Double I_TWOBODYOVERLAP_F2xy_Gx2yz = I_TWOBODYOVERLAP_G2xyz_Fx2y+ABZ*I_TWOBODYOVERLAP_F2xy_Fx2y;
  Double I_TWOBODYOVERLAP_F2xz_Gx2yz = I_TWOBODYOVERLAP_G2x2z_Fx2y+ABZ*I_TWOBODYOVERLAP_F2xz_Fx2y;
  Double I_TWOBODYOVERLAP_Fx2y_Gx2yz = I_TWOBODYOVERLAP_Gx2yz_Fx2y+ABZ*I_TWOBODYOVERLAP_Fx2y_Fx2y;
  Double I_TWOBODYOVERLAP_Fxyz_Gx2yz = I_TWOBODYOVERLAP_Gxy2z_Fx2y+ABZ*I_TWOBODYOVERLAP_Fxyz_Fx2y;
  Double I_TWOBODYOVERLAP_Fx2z_Gx2yz = I_TWOBODYOVERLAP_Gx3z_Fx2y+ABZ*I_TWOBODYOVERLAP_Fx2z_Fx2y;
  Double I_TWOBODYOVERLAP_F3y_Gx2yz = I_TWOBODYOVERLAP_G3yz_Fx2y+ABZ*I_TWOBODYOVERLAP_F3y_Fx2y;
  Double I_TWOBODYOVERLAP_F2yz_Gx2yz = I_TWOBODYOVERLAP_G2y2z_Fx2y+ABZ*I_TWOBODYOVERLAP_F2yz_Fx2y;
  Double I_TWOBODYOVERLAP_Fy2z_Gx2yz = I_TWOBODYOVERLAP_Gy3z_Fx2y+ABZ*I_TWOBODYOVERLAP_Fy2z_Fx2y;
  Double I_TWOBODYOVERLAP_F3z_Gx2yz = I_TWOBODYOVERLAP_G4z_Fx2y+ABZ*I_TWOBODYOVERLAP_F3z_Fx2y;
  Double I_TWOBODYOVERLAP_F3x_Gxy2z = I_TWOBODYOVERLAP_G3xy_Fx2z+ABY*I_TWOBODYOVERLAP_F3x_Fx2z;
  Double I_TWOBODYOVERLAP_F2xy_Gxy2z = I_TWOBODYOVERLAP_G2x2y_Fx2z+ABY*I_TWOBODYOVERLAP_F2xy_Fx2z;
  Double I_TWOBODYOVERLAP_F2xz_Gxy2z = I_TWOBODYOVERLAP_G2xyz_Fx2z+ABY*I_TWOBODYOVERLAP_F2xz_Fx2z;
  Double I_TWOBODYOVERLAP_Fx2y_Gxy2z = I_TWOBODYOVERLAP_Gx3y_Fx2z+ABY*I_TWOBODYOVERLAP_Fx2y_Fx2z;
  Double I_TWOBODYOVERLAP_Fxyz_Gxy2z = I_TWOBODYOVERLAP_Gx2yz_Fx2z+ABY*I_TWOBODYOVERLAP_Fxyz_Fx2z;
  Double I_TWOBODYOVERLAP_Fx2z_Gxy2z = I_TWOBODYOVERLAP_Gxy2z_Fx2z+ABY*I_TWOBODYOVERLAP_Fx2z_Fx2z;
  Double I_TWOBODYOVERLAP_F3y_Gxy2z = I_TWOBODYOVERLAP_G4y_Fx2z+ABY*I_TWOBODYOVERLAP_F3y_Fx2z;
  Double I_TWOBODYOVERLAP_F2yz_Gxy2z = I_TWOBODYOVERLAP_G3yz_Fx2z+ABY*I_TWOBODYOVERLAP_F2yz_Fx2z;
  Double I_TWOBODYOVERLAP_Fy2z_Gxy2z = I_TWOBODYOVERLAP_G2y2z_Fx2z+ABY*I_TWOBODYOVERLAP_Fy2z_Fx2z;
  Double I_TWOBODYOVERLAP_F3z_Gxy2z = I_TWOBODYOVERLAP_Gy3z_Fx2z+ABY*I_TWOBODYOVERLAP_F3z_Fx2z;
  Double I_TWOBODYOVERLAP_F3x_Gx3z = I_TWOBODYOVERLAP_G4x_F3z+ABX*I_TWOBODYOVERLAP_F3x_F3z;
  Double I_TWOBODYOVERLAP_F2xy_Gx3z = I_TWOBODYOVERLAP_G3xy_F3z+ABX*I_TWOBODYOVERLAP_F2xy_F3z;
  Double I_TWOBODYOVERLAP_F2xz_Gx3z = I_TWOBODYOVERLAP_G3xz_F3z+ABX*I_TWOBODYOVERLAP_F2xz_F3z;
  Double I_TWOBODYOVERLAP_Fx2y_Gx3z = I_TWOBODYOVERLAP_G2x2y_F3z+ABX*I_TWOBODYOVERLAP_Fx2y_F3z;
  Double I_TWOBODYOVERLAP_Fxyz_Gx3z = I_TWOBODYOVERLAP_G2xyz_F3z+ABX*I_TWOBODYOVERLAP_Fxyz_F3z;
  Double I_TWOBODYOVERLAP_Fx2z_Gx3z = I_TWOBODYOVERLAP_G2x2z_F3z+ABX*I_TWOBODYOVERLAP_Fx2z_F3z;
  Double I_TWOBODYOVERLAP_F3y_Gx3z = I_TWOBODYOVERLAP_Gx3y_F3z+ABX*I_TWOBODYOVERLAP_F3y_F3z;
  Double I_TWOBODYOVERLAP_F2yz_Gx3z = I_TWOBODYOVERLAP_Gx2yz_F3z+ABX*I_TWOBODYOVERLAP_F2yz_F3z;
  Double I_TWOBODYOVERLAP_Fy2z_Gx3z = I_TWOBODYOVERLAP_Gxy2z_F3z+ABX*I_TWOBODYOVERLAP_Fy2z_F3z;
  Double I_TWOBODYOVERLAP_F3z_Gx3z = I_TWOBODYOVERLAP_Gx3z_F3z+ABX*I_TWOBODYOVERLAP_F3z_F3z;
  Double I_TWOBODYOVERLAP_F3x_G4y = I_TWOBODYOVERLAP_G3xy_F3y+ABY*I_TWOBODYOVERLAP_F3x_F3y;
  Double I_TWOBODYOVERLAP_F2xy_G4y = I_TWOBODYOVERLAP_G2x2y_F3y+ABY*I_TWOBODYOVERLAP_F2xy_F3y;
  Double I_TWOBODYOVERLAP_F2xz_G4y = I_TWOBODYOVERLAP_G2xyz_F3y+ABY*I_TWOBODYOVERLAP_F2xz_F3y;
  Double I_TWOBODYOVERLAP_Fx2y_G4y = I_TWOBODYOVERLAP_Gx3y_F3y+ABY*I_TWOBODYOVERLAP_Fx2y_F3y;
  Double I_TWOBODYOVERLAP_Fxyz_G4y = I_TWOBODYOVERLAP_Gx2yz_F3y+ABY*I_TWOBODYOVERLAP_Fxyz_F3y;
  Double I_TWOBODYOVERLAP_Fx2z_G4y = I_TWOBODYOVERLAP_Gxy2z_F3y+ABY*I_TWOBODYOVERLAP_Fx2z_F3y;
  Double I_TWOBODYOVERLAP_F3y_G4y = I_TWOBODYOVERLAP_G4y_F3y+ABY*I_TWOBODYOVERLAP_F3y_F3y;
  Double I_TWOBODYOVERLAP_F2yz_G4y = I_TWOBODYOVERLAP_G3yz_F3y+ABY*I_TWOBODYOVERLAP_F2yz_F3y;
  Double I_TWOBODYOVERLAP_Fy2z_G4y = I_TWOBODYOVERLAP_G2y2z_F3y+ABY*I_TWOBODYOVERLAP_Fy2z_F3y;
  Double I_TWOBODYOVERLAP_F3z_G4y = I_TWOBODYOVERLAP_Gy3z_F3y+ABY*I_TWOBODYOVERLAP_F3z_F3y;
  Double I_TWOBODYOVERLAP_F3x_G3yz = I_TWOBODYOVERLAP_G3xz_F3y+ABZ*I_TWOBODYOVERLAP_F3x_F3y;
  Double I_TWOBODYOVERLAP_F2xy_G3yz = I_TWOBODYOVERLAP_G2xyz_F3y+ABZ*I_TWOBODYOVERLAP_F2xy_F3y;
  Double I_TWOBODYOVERLAP_F2xz_G3yz = I_TWOBODYOVERLAP_G2x2z_F3y+ABZ*I_TWOBODYOVERLAP_F2xz_F3y;
  Double I_TWOBODYOVERLAP_Fx2y_G3yz = I_TWOBODYOVERLAP_Gx2yz_F3y+ABZ*I_TWOBODYOVERLAP_Fx2y_F3y;
  Double I_TWOBODYOVERLAP_Fxyz_G3yz = I_TWOBODYOVERLAP_Gxy2z_F3y+ABZ*I_TWOBODYOVERLAP_Fxyz_F3y;
  Double I_TWOBODYOVERLAP_Fx2z_G3yz = I_TWOBODYOVERLAP_Gx3z_F3y+ABZ*I_TWOBODYOVERLAP_Fx2z_F3y;
  Double I_TWOBODYOVERLAP_F3y_G3yz = I_TWOBODYOVERLAP_G3yz_F3y+ABZ*I_TWOBODYOVERLAP_F3y_F3y;
  Double I_TWOBODYOVERLAP_F2yz_G3yz = I_TWOBODYOVERLAP_G2y2z_F3y+ABZ*I_TWOBODYOVERLAP_F2yz_F3y;
  Double I_TWOBODYOVERLAP_Fy2z_G3yz = I_TWOBODYOVERLAP_Gy3z_F3y+ABZ*I_TWOBODYOVERLAP_Fy2z_F3y;
  Double I_TWOBODYOVERLAP_F3z_G3yz = I_TWOBODYOVERLAP_G4z_F3y+ABZ*I_TWOBODYOVERLAP_F3z_F3y;
  Double I_TWOBODYOVERLAP_F3x_G2y2z = I_TWOBODYOVERLAP_G3xz_F2yz+ABZ*I_TWOBODYOVERLAP_F3x_F2yz;
  Double I_TWOBODYOVERLAP_F2xy_G2y2z = I_TWOBODYOVERLAP_G2xyz_F2yz+ABZ*I_TWOBODYOVERLAP_F2xy_F2yz;
  Double I_TWOBODYOVERLAP_F2xz_G2y2z = I_TWOBODYOVERLAP_G2x2z_F2yz+ABZ*I_TWOBODYOVERLAP_F2xz_F2yz;
  Double I_TWOBODYOVERLAP_Fx2y_G2y2z = I_TWOBODYOVERLAP_Gx2yz_F2yz+ABZ*I_TWOBODYOVERLAP_Fx2y_F2yz;
  Double I_TWOBODYOVERLAP_Fxyz_G2y2z = I_TWOBODYOVERLAP_Gxy2z_F2yz+ABZ*I_TWOBODYOVERLAP_Fxyz_F2yz;
  Double I_TWOBODYOVERLAP_Fx2z_G2y2z = I_TWOBODYOVERLAP_Gx3z_F2yz+ABZ*I_TWOBODYOVERLAP_Fx2z_F2yz;
  Double I_TWOBODYOVERLAP_F3y_G2y2z = I_TWOBODYOVERLAP_G3yz_F2yz+ABZ*I_TWOBODYOVERLAP_F3y_F2yz;
  Double I_TWOBODYOVERLAP_F2yz_G2y2z = I_TWOBODYOVERLAP_G2y2z_F2yz+ABZ*I_TWOBODYOVERLAP_F2yz_F2yz;
  Double I_TWOBODYOVERLAP_Fy2z_G2y2z = I_TWOBODYOVERLAP_Gy3z_F2yz+ABZ*I_TWOBODYOVERLAP_Fy2z_F2yz;
  Double I_TWOBODYOVERLAP_F3z_G2y2z = I_TWOBODYOVERLAP_G4z_F2yz+ABZ*I_TWOBODYOVERLAP_F3z_F2yz;
  Double I_TWOBODYOVERLAP_F3x_Gy3z = I_TWOBODYOVERLAP_G3xy_F3z+ABY*I_TWOBODYOVERLAP_F3x_F3z;
  Double I_TWOBODYOVERLAP_F2xy_Gy3z = I_TWOBODYOVERLAP_G2x2y_F3z+ABY*I_TWOBODYOVERLAP_F2xy_F3z;
  Double I_TWOBODYOVERLAP_F2xz_Gy3z = I_TWOBODYOVERLAP_G2xyz_F3z+ABY*I_TWOBODYOVERLAP_F2xz_F3z;
  Double I_TWOBODYOVERLAP_Fx2y_Gy3z = I_TWOBODYOVERLAP_Gx3y_F3z+ABY*I_TWOBODYOVERLAP_Fx2y_F3z;
  Double I_TWOBODYOVERLAP_Fxyz_Gy3z = I_TWOBODYOVERLAP_Gx2yz_F3z+ABY*I_TWOBODYOVERLAP_Fxyz_F3z;
  Double I_TWOBODYOVERLAP_Fx2z_Gy3z = I_TWOBODYOVERLAP_Gxy2z_F3z+ABY*I_TWOBODYOVERLAP_Fx2z_F3z;
  Double I_TWOBODYOVERLAP_F3y_Gy3z = I_TWOBODYOVERLAP_G4y_F3z+ABY*I_TWOBODYOVERLAP_F3y_F3z;
  Double I_TWOBODYOVERLAP_F2yz_Gy3z = I_TWOBODYOVERLAP_G3yz_F3z+ABY*I_TWOBODYOVERLAP_F2yz_F3z;
  Double I_TWOBODYOVERLAP_Fy2z_Gy3z = I_TWOBODYOVERLAP_G2y2z_F3z+ABY*I_TWOBODYOVERLAP_Fy2z_F3z;
  Double I_TWOBODYOVERLAP_F3z_Gy3z = I_TWOBODYOVERLAP_Gy3z_F3z+ABY*I_TWOBODYOVERLAP_F3z_F3z;
  Double I_TWOBODYOVERLAP_F3x_G4z = I_TWOBODYOVERLAP_G3xz_F3z+ABZ*I_TWOBODYOVERLAP_F3x_F3z;
  Double I_TWOBODYOVERLAP_F2xy_G4z = I_TWOBODYOVERLAP_G2xyz_F3z+ABZ*I_TWOBODYOVERLAP_F2xy_F3z;
  Double I_TWOBODYOVERLAP_F2xz_G4z = I_TWOBODYOVERLAP_G2x2z_F3z+ABZ*I_TWOBODYOVERLAP_F2xz_F3z;
  Double I_TWOBODYOVERLAP_Fx2y_G4z = I_TWOBODYOVERLAP_Gx2yz_F3z+ABZ*I_TWOBODYOVERLAP_Fx2y_F3z;
  Double I_TWOBODYOVERLAP_Fxyz_G4z = I_TWOBODYOVERLAP_Gxy2z_F3z+ABZ*I_TWOBODYOVERLAP_Fxyz_F3z;
  Double I_TWOBODYOVERLAP_Fx2z_G4z = I_TWOBODYOVERLAP_Gx3z_F3z+ABZ*I_TWOBODYOVERLAP_Fx2z_F3z;
  Double I_TWOBODYOVERLAP_F3y_G4z = I_TWOBODYOVERLAP_G3yz_F3z+ABZ*I_TWOBODYOVERLAP_F3y_F3z;
  Double I_TWOBODYOVERLAP_F2yz_G4z = I_TWOBODYOVERLAP_G2y2z_F3z+ABZ*I_TWOBODYOVERLAP_F2yz_F3z;
  Double I_TWOBODYOVERLAP_Fy2z_G4z = I_TWOBODYOVERLAP_Gy3z_F3z+ABZ*I_TWOBODYOVERLAP_Fy2z_F3z;
  Double I_TWOBODYOVERLAP_F3z_G4z = I_TWOBODYOVERLAP_G4z_F3z+ABZ*I_TWOBODYOVERLAP_F3z_F3z;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_H5x_Px_a = I_TWOBODYOVERLAP_I6x_S_a+ABX*I_TWOBODYOVERLAP_H5x_S_a;
  Double I_TWOBODYOVERLAP_H4xy_Px_a = I_TWOBODYOVERLAP_I5xy_S_a+ABX*I_TWOBODYOVERLAP_H4xy_S_a;
  Double I_TWOBODYOVERLAP_H4xz_Px_a = I_TWOBODYOVERLAP_I5xz_S_a+ABX*I_TWOBODYOVERLAP_H4xz_S_a;
  Double I_TWOBODYOVERLAP_H3x2y_Px_a = I_TWOBODYOVERLAP_I4x2y_S_a+ABX*I_TWOBODYOVERLAP_H3x2y_S_a;
  Double I_TWOBODYOVERLAP_H3xyz_Px_a = I_TWOBODYOVERLAP_I4xyz_S_a+ABX*I_TWOBODYOVERLAP_H3xyz_S_a;
  Double I_TWOBODYOVERLAP_H3x2z_Px_a = I_TWOBODYOVERLAP_I4x2z_S_a+ABX*I_TWOBODYOVERLAP_H3x2z_S_a;
  Double I_TWOBODYOVERLAP_H2x3y_Px_a = I_TWOBODYOVERLAP_I3x3y_S_a+ABX*I_TWOBODYOVERLAP_H2x3y_S_a;
  Double I_TWOBODYOVERLAP_H2x2yz_Px_a = I_TWOBODYOVERLAP_I3x2yz_S_a+ABX*I_TWOBODYOVERLAP_H2x2yz_S_a;
  Double I_TWOBODYOVERLAP_H2xy2z_Px_a = I_TWOBODYOVERLAP_I3xy2z_S_a+ABX*I_TWOBODYOVERLAP_H2xy2z_S_a;
  Double I_TWOBODYOVERLAP_H2x3z_Px_a = I_TWOBODYOVERLAP_I3x3z_S_a+ABX*I_TWOBODYOVERLAP_H2x3z_S_a;
  Double I_TWOBODYOVERLAP_Hx4y_Px_a = I_TWOBODYOVERLAP_I2x4y_S_a+ABX*I_TWOBODYOVERLAP_Hx4y_S_a;
  Double I_TWOBODYOVERLAP_Hx3yz_Px_a = I_TWOBODYOVERLAP_I2x3yz_S_a+ABX*I_TWOBODYOVERLAP_Hx3yz_S_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_Px_a = I_TWOBODYOVERLAP_I2x2y2z_S_a+ABX*I_TWOBODYOVERLAP_Hx2y2z_S_a;
  Double I_TWOBODYOVERLAP_Hxy3z_Px_a = I_TWOBODYOVERLAP_I2xy3z_S_a+ABX*I_TWOBODYOVERLAP_Hxy3z_S_a;
  Double I_TWOBODYOVERLAP_Hx4z_Px_a = I_TWOBODYOVERLAP_I2x4z_S_a+ABX*I_TWOBODYOVERLAP_Hx4z_S_a;
  Double I_TWOBODYOVERLAP_H5y_Px_a = I_TWOBODYOVERLAP_Ix5y_S_a+ABX*I_TWOBODYOVERLAP_H5y_S_a;
  Double I_TWOBODYOVERLAP_H4yz_Px_a = I_TWOBODYOVERLAP_Ix4yz_S_a+ABX*I_TWOBODYOVERLAP_H4yz_S_a;
  Double I_TWOBODYOVERLAP_H3y2z_Px_a = I_TWOBODYOVERLAP_Ix3y2z_S_a+ABX*I_TWOBODYOVERLAP_H3y2z_S_a;
  Double I_TWOBODYOVERLAP_H2y3z_Px_a = I_TWOBODYOVERLAP_Ix2y3z_S_a+ABX*I_TWOBODYOVERLAP_H2y3z_S_a;
  Double I_TWOBODYOVERLAP_Hy4z_Px_a = I_TWOBODYOVERLAP_Ixy4z_S_a+ABX*I_TWOBODYOVERLAP_Hy4z_S_a;
  Double I_TWOBODYOVERLAP_H5z_Px_a = I_TWOBODYOVERLAP_Ix5z_S_a+ABX*I_TWOBODYOVERLAP_H5z_S_a;
  Double I_TWOBODYOVERLAP_H5x_Py_a = I_TWOBODYOVERLAP_I5xy_S_a+ABY*I_TWOBODYOVERLAP_H5x_S_a;
  Double I_TWOBODYOVERLAP_H4xy_Py_a = I_TWOBODYOVERLAP_I4x2y_S_a+ABY*I_TWOBODYOVERLAP_H4xy_S_a;
  Double I_TWOBODYOVERLAP_H4xz_Py_a = I_TWOBODYOVERLAP_I4xyz_S_a+ABY*I_TWOBODYOVERLAP_H4xz_S_a;
  Double I_TWOBODYOVERLAP_H3x2y_Py_a = I_TWOBODYOVERLAP_I3x3y_S_a+ABY*I_TWOBODYOVERLAP_H3x2y_S_a;
  Double I_TWOBODYOVERLAP_H3xyz_Py_a = I_TWOBODYOVERLAP_I3x2yz_S_a+ABY*I_TWOBODYOVERLAP_H3xyz_S_a;
  Double I_TWOBODYOVERLAP_H3x2z_Py_a = I_TWOBODYOVERLAP_I3xy2z_S_a+ABY*I_TWOBODYOVERLAP_H3x2z_S_a;
  Double I_TWOBODYOVERLAP_H2x3y_Py_a = I_TWOBODYOVERLAP_I2x4y_S_a+ABY*I_TWOBODYOVERLAP_H2x3y_S_a;
  Double I_TWOBODYOVERLAP_H2x2yz_Py_a = I_TWOBODYOVERLAP_I2x3yz_S_a+ABY*I_TWOBODYOVERLAP_H2x2yz_S_a;
  Double I_TWOBODYOVERLAP_H2xy2z_Py_a = I_TWOBODYOVERLAP_I2x2y2z_S_a+ABY*I_TWOBODYOVERLAP_H2xy2z_S_a;
  Double I_TWOBODYOVERLAP_H2x3z_Py_a = I_TWOBODYOVERLAP_I2xy3z_S_a+ABY*I_TWOBODYOVERLAP_H2x3z_S_a;
  Double I_TWOBODYOVERLAP_Hx4y_Py_a = I_TWOBODYOVERLAP_Ix5y_S_a+ABY*I_TWOBODYOVERLAP_Hx4y_S_a;
  Double I_TWOBODYOVERLAP_Hx3yz_Py_a = I_TWOBODYOVERLAP_Ix4yz_S_a+ABY*I_TWOBODYOVERLAP_Hx3yz_S_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_Py_a = I_TWOBODYOVERLAP_Ix3y2z_S_a+ABY*I_TWOBODYOVERLAP_Hx2y2z_S_a;
  Double I_TWOBODYOVERLAP_Hxy3z_Py_a = I_TWOBODYOVERLAP_Ix2y3z_S_a+ABY*I_TWOBODYOVERLAP_Hxy3z_S_a;
  Double I_TWOBODYOVERLAP_Hx4z_Py_a = I_TWOBODYOVERLAP_Ixy4z_S_a+ABY*I_TWOBODYOVERLAP_Hx4z_S_a;
  Double I_TWOBODYOVERLAP_H5y_Py_a = I_TWOBODYOVERLAP_I6y_S_a+ABY*I_TWOBODYOVERLAP_H5y_S_a;
  Double I_TWOBODYOVERLAP_H4yz_Py_a = I_TWOBODYOVERLAP_I5yz_S_a+ABY*I_TWOBODYOVERLAP_H4yz_S_a;
  Double I_TWOBODYOVERLAP_H3y2z_Py_a = I_TWOBODYOVERLAP_I4y2z_S_a+ABY*I_TWOBODYOVERLAP_H3y2z_S_a;
  Double I_TWOBODYOVERLAP_H2y3z_Py_a = I_TWOBODYOVERLAP_I3y3z_S_a+ABY*I_TWOBODYOVERLAP_H2y3z_S_a;
  Double I_TWOBODYOVERLAP_Hy4z_Py_a = I_TWOBODYOVERLAP_I2y4z_S_a+ABY*I_TWOBODYOVERLAP_Hy4z_S_a;
  Double I_TWOBODYOVERLAP_H5z_Py_a = I_TWOBODYOVERLAP_Iy5z_S_a+ABY*I_TWOBODYOVERLAP_H5z_S_a;
  Double I_TWOBODYOVERLAP_H5x_Pz_a = I_TWOBODYOVERLAP_I5xz_S_a+ABZ*I_TWOBODYOVERLAP_H5x_S_a;
  Double I_TWOBODYOVERLAP_H4xy_Pz_a = I_TWOBODYOVERLAP_I4xyz_S_a+ABZ*I_TWOBODYOVERLAP_H4xy_S_a;
  Double I_TWOBODYOVERLAP_H4xz_Pz_a = I_TWOBODYOVERLAP_I4x2z_S_a+ABZ*I_TWOBODYOVERLAP_H4xz_S_a;
  Double I_TWOBODYOVERLAP_H3x2y_Pz_a = I_TWOBODYOVERLAP_I3x2yz_S_a+ABZ*I_TWOBODYOVERLAP_H3x2y_S_a;
  Double I_TWOBODYOVERLAP_H3xyz_Pz_a = I_TWOBODYOVERLAP_I3xy2z_S_a+ABZ*I_TWOBODYOVERLAP_H3xyz_S_a;
  Double I_TWOBODYOVERLAP_H3x2z_Pz_a = I_TWOBODYOVERLAP_I3x3z_S_a+ABZ*I_TWOBODYOVERLAP_H3x2z_S_a;
  Double I_TWOBODYOVERLAP_H2x3y_Pz_a = I_TWOBODYOVERLAP_I2x3yz_S_a+ABZ*I_TWOBODYOVERLAP_H2x3y_S_a;
  Double I_TWOBODYOVERLAP_H2x2yz_Pz_a = I_TWOBODYOVERLAP_I2x2y2z_S_a+ABZ*I_TWOBODYOVERLAP_H2x2yz_S_a;
  Double I_TWOBODYOVERLAP_H2xy2z_Pz_a = I_TWOBODYOVERLAP_I2xy3z_S_a+ABZ*I_TWOBODYOVERLAP_H2xy2z_S_a;
  Double I_TWOBODYOVERLAP_H2x3z_Pz_a = I_TWOBODYOVERLAP_I2x4z_S_a+ABZ*I_TWOBODYOVERLAP_H2x3z_S_a;
  Double I_TWOBODYOVERLAP_Hx4y_Pz_a = I_TWOBODYOVERLAP_Ix4yz_S_a+ABZ*I_TWOBODYOVERLAP_Hx4y_S_a;
  Double I_TWOBODYOVERLAP_Hx3yz_Pz_a = I_TWOBODYOVERLAP_Ix3y2z_S_a+ABZ*I_TWOBODYOVERLAP_Hx3yz_S_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_Pz_a = I_TWOBODYOVERLAP_Ix2y3z_S_a+ABZ*I_TWOBODYOVERLAP_Hx2y2z_S_a;
  Double I_TWOBODYOVERLAP_Hxy3z_Pz_a = I_TWOBODYOVERLAP_Ixy4z_S_a+ABZ*I_TWOBODYOVERLAP_Hxy3z_S_a;
  Double I_TWOBODYOVERLAP_Hx4z_Pz_a = I_TWOBODYOVERLAP_Ix5z_S_a+ABZ*I_TWOBODYOVERLAP_Hx4z_S_a;
  Double I_TWOBODYOVERLAP_H5y_Pz_a = I_TWOBODYOVERLAP_I5yz_S_a+ABZ*I_TWOBODYOVERLAP_H5y_S_a;
  Double I_TWOBODYOVERLAP_H4yz_Pz_a = I_TWOBODYOVERLAP_I4y2z_S_a+ABZ*I_TWOBODYOVERLAP_H4yz_S_a;
  Double I_TWOBODYOVERLAP_H3y2z_Pz_a = I_TWOBODYOVERLAP_I3y3z_S_a+ABZ*I_TWOBODYOVERLAP_H3y2z_S_a;
  Double I_TWOBODYOVERLAP_H2y3z_Pz_a = I_TWOBODYOVERLAP_I2y4z_S_a+ABZ*I_TWOBODYOVERLAP_H2y3z_S_a;
  Double I_TWOBODYOVERLAP_Hy4z_Pz_a = I_TWOBODYOVERLAP_Iy5z_S_a+ABZ*I_TWOBODYOVERLAP_Hy4z_S_a;
  Double I_TWOBODYOVERLAP_H5z_Pz_a = I_TWOBODYOVERLAP_I6z_S_a+ABZ*I_TWOBODYOVERLAP_H5z_S_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_I_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_I6x_Px_a = I_TWOBODYOVERLAP_K7x_S_a+ABX*I_TWOBODYOVERLAP_I6x_S_a;
  Double I_TWOBODYOVERLAP_I5xy_Px_a = I_TWOBODYOVERLAP_K6xy_S_a+ABX*I_TWOBODYOVERLAP_I5xy_S_a;
  Double I_TWOBODYOVERLAP_I5xz_Px_a = I_TWOBODYOVERLAP_K6xz_S_a+ABX*I_TWOBODYOVERLAP_I5xz_S_a;
  Double I_TWOBODYOVERLAP_I4x2y_Px_a = I_TWOBODYOVERLAP_K5x2y_S_a+ABX*I_TWOBODYOVERLAP_I4x2y_S_a;
  Double I_TWOBODYOVERLAP_I4xyz_Px_a = I_TWOBODYOVERLAP_K5xyz_S_a+ABX*I_TWOBODYOVERLAP_I4xyz_S_a;
  Double I_TWOBODYOVERLAP_I4x2z_Px_a = I_TWOBODYOVERLAP_K5x2z_S_a+ABX*I_TWOBODYOVERLAP_I4x2z_S_a;
  Double I_TWOBODYOVERLAP_I3x3y_Px_a = I_TWOBODYOVERLAP_K4x3y_S_a+ABX*I_TWOBODYOVERLAP_I3x3y_S_a;
  Double I_TWOBODYOVERLAP_I3x2yz_Px_a = I_TWOBODYOVERLAP_K4x2yz_S_a+ABX*I_TWOBODYOVERLAP_I3x2yz_S_a;
  Double I_TWOBODYOVERLAP_I3xy2z_Px_a = I_TWOBODYOVERLAP_K4xy2z_S_a+ABX*I_TWOBODYOVERLAP_I3xy2z_S_a;
  Double I_TWOBODYOVERLAP_I3x3z_Px_a = I_TWOBODYOVERLAP_K4x3z_S_a+ABX*I_TWOBODYOVERLAP_I3x3z_S_a;
  Double I_TWOBODYOVERLAP_I2x4y_Px_a = I_TWOBODYOVERLAP_K3x4y_S_a+ABX*I_TWOBODYOVERLAP_I2x4y_S_a;
  Double I_TWOBODYOVERLAP_I2x3yz_Px_a = I_TWOBODYOVERLAP_K3x3yz_S_a+ABX*I_TWOBODYOVERLAP_I2x3yz_S_a;
  Double I_TWOBODYOVERLAP_I2x2y2z_Px_a = I_TWOBODYOVERLAP_K3x2y2z_S_a+ABX*I_TWOBODYOVERLAP_I2x2y2z_S_a;
  Double I_TWOBODYOVERLAP_I2xy3z_Px_a = I_TWOBODYOVERLAP_K3xy3z_S_a+ABX*I_TWOBODYOVERLAP_I2xy3z_S_a;
  Double I_TWOBODYOVERLAP_I2x4z_Px_a = I_TWOBODYOVERLAP_K3x4z_S_a+ABX*I_TWOBODYOVERLAP_I2x4z_S_a;
  Double I_TWOBODYOVERLAP_Ix5y_Px_a = I_TWOBODYOVERLAP_K2x5y_S_a+ABX*I_TWOBODYOVERLAP_Ix5y_S_a;
  Double I_TWOBODYOVERLAP_Ix4yz_Px_a = I_TWOBODYOVERLAP_K2x4yz_S_a+ABX*I_TWOBODYOVERLAP_Ix4yz_S_a;
  Double I_TWOBODYOVERLAP_Ix3y2z_Px_a = I_TWOBODYOVERLAP_K2x3y2z_S_a+ABX*I_TWOBODYOVERLAP_Ix3y2z_S_a;
  Double I_TWOBODYOVERLAP_Ix2y3z_Px_a = I_TWOBODYOVERLAP_K2x2y3z_S_a+ABX*I_TWOBODYOVERLAP_Ix2y3z_S_a;
  Double I_TWOBODYOVERLAP_Ixy4z_Px_a = I_TWOBODYOVERLAP_K2xy4z_S_a+ABX*I_TWOBODYOVERLAP_Ixy4z_S_a;
  Double I_TWOBODYOVERLAP_Ix5z_Px_a = I_TWOBODYOVERLAP_K2x5z_S_a+ABX*I_TWOBODYOVERLAP_Ix5z_S_a;
  Double I_TWOBODYOVERLAP_I6y_Px_a = I_TWOBODYOVERLAP_Kx6y_S_a+ABX*I_TWOBODYOVERLAP_I6y_S_a;
  Double I_TWOBODYOVERLAP_I5yz_Px_a = I_TWOBODYOVERLAP_Kx5yz_S_a+ABX*I_TWOBODYOVERLAP_I5yz_S_a;
  Double I_TWOBODYOVERLAP_I4y2z_Px_a = I_TWOBODYOVERLAP_Kx4y2z_S_a+ABX*I_TWOBODYOVERLAP_I4y2z_S_a;
  Double I_TWOBODYOVERLAP_I3y3z_Px_a = I_TWOBODYOVERLAP_Kx3y3z_S_a+ABX*I_TWOBODYOVERLAP_I3y3z_S_a;
  Double I_TWOBODYOVERLAP_I2y4z_Px_a = I_TWOBODYOVERLAP_Kx2y4z_S_a+ABX*I_TWOBODYOVERLAP_I2y4z_S_a;
  Double I_TWOBODYOVERLAP_Iy5z_Px_a = I_TWOBODYOVERLAP_Kxy5z_S_a+ABX*I_TWOBODYOVERLAP_Iy5z_S_a;
  Double I_TWOBODYOVERLAP_I6z_Px_a = I_TWOBODYOVERLAP_Kx6z_S_a+ABX*I_TWOBODYOVERLAP_I6z_S_a;
  Double I_TWOBODYOVERLAP_I6x_Py_a = I_TWOBODYOVERLAP_K6xy_S_a+ABY*I_TWOBODYOVERLAP_I6x_S_a;
  Double I_TWOBODYOVERLAP_I5xy_Py_a = I_TWOBODYOVERLAP_K5x2y_S_a+ABY*I_TWOBODYOVERLAP_I5xy_S_a;
  Double I_TWOBODYOVERLAP_I5xz_Py_a = I_TWOBODYOVERLAP_K5xyz_S_a+ABY*I_TWOBODYOVERLAP_I5xz_S_a;
  Double I_TWOBODYOVERLAP_I4x2y_Py_a = I_TWOBODYOVERLAP_K4x3y_S_a+ABY*I_TWOBODYOVERLAP_I4x2y_S_a;
  Double I_TWOBODYOVERLAP_I4xyz_Py_a = I_TWOBODYOVERLAP_K4x2yz_S_a+ABY*I_TWOBODYOVERLAP_I4xyz_S_a;
  Double I_TWOBODYOVERLAP_I4x2z_Py_a = I_TWOBODYOVERLAP_K4xy2z_S_a+ABY*I_TWOBODYOVERLAP_I4x2z_S_a;
  Double I_TWOBODYOVERLAP_I3x3y_Py_a = I_TWOBODYOVERLAP_K3x4y_S_a+ABY*I_TWOBODYOVERLAP_I3x3y_S_a;
  Double I_TWOBODYOVERLAP_I3x2yz_Py_a = I_TWOBODYOVERLAP_K3x3yz_S_a+ABY*I_TWOBODYOVERLAP_I3x2yz_S_a;
  Double I_TWOBODYOVERLAP_I3xy2z_Py_a = I_TWOBODYOVERLAP_K3x2y2z_S_a+ABY*I_TWOBODYOVERLAP_I3xy2z_S_a;
  Double I_TWOBODYOVERLAP_I3x3z_Py_a = I_TWOBODYOVERLAP_K3xy3z_S_a+ABY*I_TWOBODYOVERLAP_I3x3z_S_a;
  Double I_TWOBODYOVERLAP_I2x4y_Py_a = I_TWOBODYOVERLAP_K2x5y_S_a+ABY*I_TWOBODYOVERLAP_I2x4y_S_a;
  Double I_TWOBODYOVERLAP_I2x3yz_Py_a = I_TWOBODYOVERLAP_K2x4yz_S_a+ABY*I_TWOBODYOVERLAP_I2x3yz_S_a;
  Double I_TWOBODYOVERLAP_I2x2y2z_Py_a = I_TWOBODYOVERLAP_K2x3y2z_S_a+ABY*I_TWOBODYOVERLAP_I2x2y2z_S_a;
  Double I_TWOBODYOVERLAP_I2xy3z_Py_a = I_TWOBODYOVERLAP_K2x2y3z_S_a+ABY*I_TWOBODYOVERLAP_I2xy3z_S_a;
  Double I_TWOBODYOVERLAP_I2x4z_Py_a = I_TWOBODYOVERLAP_K2xy4z_S_a+ABY*I_TWOBODYOVERLAP_I2x4z_S_a;
  Double I_TWOBODYOVERLAP_Ix5y_Py_a = I_TWOBODYOVERLAP_Kx6y_S_a+ABY*I_TWOBODYOVERLAP_Ix5y_S_a;
  Double I_TWOBODYOVERLAP_Ix4yz_Py_a = I_TWOBODYOVERLAP_Kx5yz_S_a+ABY*I_TWOBODYOVERLAP_Ix4yz_S_a;
  Double I_TWOBODYOVERLAP_Ix3y2z_Py_a = I_TWOBODYOVERLAP_Kx4y2z_S_a+ABY*I_TWOBODYOVERLAP_Ix3y2z_S_a;
  Double I_TWOBODYOVERLAP_Ix2y3z_Py_a = I_TWOBODYOVERLAP_Kx3y3z_S_a+ABY*I_TWOBODYOVERLAP_Ix2y3z_S_a;
  Double I_TWOBODYOVERLAP_Ixy4z_Py_a = I_TWOBODYOVERLAP_Kx2y4z_S_a+ABY*I_TWOBODYOVERLAP_Ixy4z_S_a;
  Double I_TWOBODYOVERLAP_Ix5z_Py_a = I_TWOBODYOVERLAP_Kxy5z_S_a+ABY*I_TWOBODYOVERLAP_Ix5z_S_a;
  Double I_TWOBODYOVERLAP_I6y_Py_a = I_TWOBODYOVERLAP_K7y_S_a+ABY*I_TWOBODYOVERLAP_I6y_S_a;
  Double I_TWOBODYOVERLAP_I5yz_Py_a = I_TWOBODYOVERLAP_K6yz_S_a+ABY*I_TWOBODYOVERLAP_I5yz_S_a;
  Double I_TWOBODYOVERLAP_I4y2z_Py_a = I_TWOBODYOVERLAP_K5y2z_S_a+ABY*I_TWOBODYOVERLAP_I4y2z_S_a;
  Double I_TWOBODYOVERLAP_I3y3z_Py_a = I_TWOBODYOVERLAP_K4y3z_S_a+ABY*I_TWOBODYOVERLAP_I3y3z_S_a;
  Double I_TWOBODYOVERLAP_I2y4z_Py_a = I_TWOBODYOVERLAP_K3y4z_S_a+ABY*I_TWOBODYOVERLAP_I2y4z_S_a;
  Double I_TWOBODYOVERLAP_Iy5z_Py_a = I_TWOBODYOVERLAP_K2y5z_S_a+ABY*I_TWOBODYOVERLAP_Iy5z_S_a;
  Double I_TWOBODYOVERLAP_I6z_Py_a = I_TWOBODYOVERLAP_Ky6z_S_a+ABY*I_TWOBODYOVERLAP_I6z_S_a;
  Double I_TWOBODYOVERLAP_I6x_Pz_a = I_TWOBODYOVERLAP_K6xz_S_a+ABZ*I_TWOBODYOVERLAP_I6x_S_a;
  Double I_TWOBODYOVERLAP_I5xy_Pz_a = I_TWOBODYOVERLAP_K5xyz_S_a+ABZ*I_TWOBODYOVERLAP_I5xy_S_a;
  Double I_TWOBODYOVERLAP_I5xz_Pz_a = I_TWOBODYOVERLAP_K5x2z_S_a+ABZ*I_TWOBODYOVERLAP_I5xz_S_a;
  Double I_TWOBODYOVERLAP_I4x2y_Pz_a = I_TWOBODYOVERLAP_K4x2yz_S_a+ABZ*I_TWOBODYOVERLAP_I4x2y_S_a;
  Double I_TWOBODYOVERLAP_I4xyz_Pz_a = I_TWOBODYOVERLAP_K4xy2z_S_a+ABZ*I_TWOBODYOVERLAP_I4xyz_S_a;
  Double I_TWOBODYOVERLAP_I4x2z_Pz_a = I_TWOBODYOVERLAP_K4x3z_S_a+ABZ*I_TWOBODYOVERLAP_I4x2z_S_a;
  Double I_TWOBODYOVERLAP_I3x3y_Pz_a = I_TWOBODYOVERLAP_K3x3yz_S_a+ABZ*I_TWOBODYOVERLAP_I3x3y_S_a;
  Double I_TWOBODYOVERLAP_I3x2yz_Pz_a = I_TWOBODYOVERLAP_K3x2y2z_S_a+ABZ*I_TWOBODYOVERLAP_I3x2yz_S_a;
  Double I_TWOBODYOVERLAP_I3xy2z_Pz_a = I_TWOBODYOVERLAP_K3xy3z_S_a+ABZ*I_TWOBODYOVERLAP_I3xy2z_S_a;
  Double I_TWOBODYOVERLAP_I3x3z_Pz_a = I_TWOBODYOVERLAP_K3x4z_S_a+ABZ*I_TWOBODYOVERLAP_I3x3z_S_a;
  Double I_TWOBODYOVERLAP_I2x4y_Pz_a = I_TWOBODYOVERLAP_K2x4yz_S_a+ABZ*I_TWOBODYOVERLAP_I2x4y_S_a;
  Double I_TWOBODYOVERLAP_I2x3yz_Pz_a = I_TWOBODYOVERLAP_K2x3y2z_S_a+ABZ*I_TWOBODYOVERLAP_I2x3yz_S_a;
  Double I_TWOBODYOVERLAP_I2x2y2z_Pz_a = I_TWOBODYOVERLAP_K2x2y3z_S_a+ABZ*I_TWOBODYOVERLAP_I2x2y2z_S_a;
  Double I_TWOBODYOVERLAP_I2xy3z_Pz_a = I_TWOBODYOVERLAP_K2xy4z_S_a+ABZ*I_TWOBODYOVERLAP_I2xy3z_S_a;
  Double I_TWOBODYOVERLAP_I2x4z_Pz_a = I_TWOBODYOVERLAP_K2x5z_S_a+ABZ*I_TWOBODYOVERLAP_I2x4z_S_a;
  Double I_TWOBODYOVERLAP_Ix5y_Pz_a = I_TWOBODYOVERLAP_Kx5yz_S_a+ABZ*I_TWOBODYOVERLAP_Ix5y_S_a;
  Double I_TWOBODYOVERLAP_Ix4yz_Pz_a = I_TWOBODYOVERLAP_Kx4y2z_S_a+ABZ*I_TWOBODYOVERLAP_Ix4yz_S_a;
  Double I_TWOBODYOVERLAP_Ix3y2z_Pz_a = I_TWOBODYOVERLAP_Kx3y3z_S_a+ABZ*I_TWOBODYOVERLAP_Ix3y2z_S_a;
  Double I_TWOBODYOVERLAP_Ix2y3z_Pz_a = I_TWOBODYOVERLAP_Kx2y4z_S_a+ABZ*I_TWOBODYOVERLAP_Ix2y3z_S_a;
  Double I_TWOBODYOVERLAP_Ixy4z_Pz_a = I_TWOBODYOVERLAP_Kxy5z_S_a+ABZ*I_TWOBODYOVERLAP_Ixy4z_S_a;
  Double I_TWOBODYOVERLAP_Ix5z_Pz_a = I_TWOBODYOVERLAP_Kx6z_S_a+ABZ*I_TWOBODYOVERLAP_Ix5z_S_a;
  Double I_TWOBODYOVERLAP_I6y_Pz_a = I_TWOBODYOVERLAP_K6yz_S_a+ABZ*I_TWOBODYOVERLAP_I6y_S_a;
  Double I_TWOBODYOVERLAP_I5yz_Pz_a = I_TWOBODYOVERLAP_K5y2z_S_a+ABZ*I_TWOBODYOVERLAP_I5yz_S_a;
  Double I_TWOBODYOVERLAP_I4y2z_Pz_a = I_TWOBODYOVERLAP_K4y3z_S_a+ABZ*I_TWOBODYOVERLAP_I4y2z_S_a;
  Double I_TWOBODYOVERLAP_I3y3z_Pz_a = I_TWOBODYOVERLAP_K3y4z_S_a+ABZ*I_TWOBODYOVERLAP_I3y3z_S_a;
  Double I_TWOBODYOVERLAP_I2y4z_Pz_a = I_TWOBODYOVERLAP_K2y5z_S_a+ABZ*I_TWOBODYOVERLAP_I2y4z_S_a;
  Double I_TWOBODYOVERLAP_Iy5z_Pz_a = I_TWOBODYOVERLAP_Ky6z_S_a+ABZ*I_TWOBODYOVERLAP_Iy5z_S_a;
  Double I_TWOBODYOVERLAP_I6z_Pz_a = I_TWOBODYOVERLAP_K7z_S_a+ABZ*I_TWOBODYOVERLAP_I6z_S_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_D_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 63 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_P_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_H5x_D2x_a = I_TWOBODYOVERLAP_I6x_Px_a+ABX*I_TWOBODYOVERLAP_H5x_Px_a;
  Double I_TWOBODYOVERLAP_H4xy_D2x_a = I_TWOBODYOVERLAP_I5xy_Px_a+ABX*I_TWOBODYOVERLAP_H4xy_Px_a;
  Double I_TWOBODYOVERLAP_H4xz_D2x_a = I_TWOBODYOVERLAP_I5xz_Px_a+ABX*I_TWOBODYOVERLAP_H4xz_Px_a;
  Double I_TWOBODYOVERLAP_H3x2y_D2x_a = I_TWOBODYOVERLAP_I4x2y_Px_a+ABX*I_TWOBODYOVERLAP_H3x2y_Px_a;
  Double I_TWOBODYOVERLAP_H3xyz_D2x_a = I_TWOBODYOVERLAP_I4xyz_Px_a+ABX*I_TWOBODYOVERLAP_H3xyz_Px_a;
  Double I_TWOBODYOVERLAP_H3x2z_D2x_a = I_TWOBODYOVERLAP_I4x2z_Px_a+ABX*I_TWOBODYOVERLAP_H3x2z_Px_a;
  Double I_TWOBODYOVERLAP_H2x3y_D2x_a = I_TWOBODYOVERLAP_I3x3y_Px_a+ABX*I_TWOBODYOVERLAP_H2x3y_Px_a;
  Double I_TWOBODYOVERLAP_H2x2yz_D2x_a = I_TWOBODYOVERLAP_I3x2yz_Px_a+ABX*I_TWOBODYOVERLAP_H2x2yz_Px_a;
  Double I_TWOBODYOVERLAP_H2xy2z_D2x_a = I_TWOBODYOVERLAP_I3xy2z_Px_a+ABX*I_TWOBODYOVERLAP_H2xy2z_Px_a;
  Double I_TWOBODYOVERLAP_H2x3z_D2x_a = I_TWOBODYOVERLAP_I3x3z_Px_a+ABX*I_TWOBODYOVERLAP_H2x3z_Px_a;
  Double I_TWOBODYOVERLAP_Hx4y_D2x_a = I_TWOBODYOVERLAP_I2x4y_Px_a+ABX*I_TWOBODYOVERLAP_Hx4y_Px_a;
  Double I_TWOBODYOVERLAP_Hx3yz_D2x_a = I_TWOBODYOVERLAP_I2x3yz_Px_a+ABX*I_TWOBODYOVERLAP_Hx3yz_Px_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_D2x_a = I_TWOBODYOVERLAP_I2x2y2z_Px_a+ABX*I_TWOBODYOVERLAP_Hx2y2z_Px_a;
  Double I_TWOBODYOVERLAP_Hxy3z_D2x_a = I_TWOBODYOVERLAP_I2xy3z_Px_a+ABX*I_TWOBODYOVERLAP_Hxy3z_Px_a;
  Double I_TWOBODYOVERLAP_Hx4z_D2x_a = I_TWOBODYOVERLAP_I2x4z_Px_a+ABX*I_TWOBODYOVERLAP_Hx4z_Px_a;
  Double I_TWOBODYOVERLAP_H5y_D2x_a = I_TWOBODYOVERLAP_Ix5y_Px_a+ABX*I_TWOBODYOVERLAP_H5y_Px_a;
  Double I_TWOBODYOVERLAP_H4yz_D2x_a = I_TWOBODYOVERLAP_Ix4yz_Px_a+ABX*I_TWOBODYOVERLAP_H4yz_Px_a;
  Double I_TWOBODYOVERLAP_H3y2z_D2x_a = I_TWOBODYOVERLAP_Ix3y2z_Px_a+ABX*I_TWOBODYOVERLAP_H3y2z_Px_a;
  Double I_TWOBODYOVERLAP_H2y3z_D2x_a = I_TWOBODYOVERLAP_Ix2y3z_Px_a+ABX*I_TWOBODYOVERLAP_H2y3z_Px_a;
  Double I_TWOBODYOVERLAP_Hy4z_D2x_a = I_TWOBODYOVERLAP_Ixy4z_Px_a+ABX*I_TWOBODYOVERLAP_Hy4z_Px_a;
  Double I_TWOBODYOVERLAP_H5z_D2x_a = I_TWOBODYOVERLAP_Ix5z_Px_a+ABX*I_TWOBODYOVERLAP_H5z_Px_a;
  Double I_TWOBODYOVERLAP_H5x_D2y_a = I_TWOBODYOVERLAP_I5xy_Py_a+ABY*I_TWOBODYOVERLAP_H5x_Py_a;
  Double I_TWOBODYOVERLAP_H4xy_D2y_a = I_TWOBODYOVERLAP_I4x2y_Py_a+ABY*I_TWOBODYOVERLAP_H4xy_Py_a;
  Double I_TWOBODYOVERLAP_H4xz_D2y_a = I_TWOBODYOVERLAP_I4xyz_Py_a+ABY*I_TWOBODYOVERLAP_H4xz_Py_a;
  Double I_TWOBODYOVERLAP_H3x2y_D2y_a = I_TWOBODYOVERLAP_I3x3y_Py_a+ABY*I_TWOBODYOVERLAP_H3x2y_Py_a;
  Double I_TWOBODYOVERLAP_H3xyz_D2y_a = I_TWOBODYOVERLAP_I3x2yz_Py_a+ABY*I_TWOBODYOVERLAP_H3xyz_Py_a;
  Double I_TWOBODYOVERLAP_H3x2z_D2y_a = I_TWOBODYOVERLAP_I3xy2z_Py_a+ABY*I_TWOBODYOVERLAP_H3x2z_Py_a;
  Double I_TWOBODYOVERLAP_H2x3y_D2y_a = I_TWOBODYOVERLAP_I2x4y_Py_a+ABY*I_TWOBODYOVERLAP_H2x3y_Py_a;
  Double I_TWOBODYOVERLAP_H2x2yz_D2y_a = I_TWOBODYOVERLAP_I2x3yz_Py_a+ABY*I_TWOBODYOVERLAP_H2x2yz_Py_a;
  Double I_TWOBODYOVERLAP_H2xy2z_D2y_a = I_TWOBODYOVERLAP_I2x2y2z_Py_a+ABY*I_TWOBODYOVERLAP_H2xy2z_Py_a;
  Double I_TWOBODYOVERLAP_H2x3z_D2y_a = I_TWOBODYOVERLAP_I2xy3z_Py_a+ABY*I_TWOBODYOVERLAP_H2x3z_Py_a;
  Double I_TWOBODYOVERLAP_Hx4y_D2y_a = I_TWOBODYOVERLAP_Ix5y_Py_a+ABY*I_TWOBODYOVERLAP_Hx4y_Py_a;
  Double I_TWOBODYOVERLAP_Hx3yz_D2y_a = I_TWOBODYOVERLAP_Ix4yz_Py_a+ABY*I_TWOBODYOVERLAP_Hx3yz_Py_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_D2y_a = I_TWOBODYOVERLAP_Ix3y2z_Py_a+ABY*I_TWOBODYOVERLAP_Hx2y2z_Py_a;
  Double I_TWOBODYOVERLAP_Hxy3z_D2y_a = I_TWOBODYOVERLAP_Ix2y3z_Py_a+ABY*I_TWOBODYOVERLAP_Hxy3z_Py_a;
  Double I_TWOBODYOVERLAP_Hx4z_D2y_a = I_TWOBODYOVERLAP_Ixy4z_Py_a+ABY*I_TWOBODYOVERLAP_Hx4z_Py_a;
  Double I_TWOBODYOVERLAP_H5y_D2y_a = I_TWOBODYOVERLAP_I6y_Py_a+ABY*I_TWOBODYOVERLAP_H5y_Py_a;
  Double I_TWOBODYOVERLAP_H4yz_D2y_a = I_TWOBODYOVERLAP_I5yz_Py_a+ABY*I_TWOBODYOVERLAP_H4yz_Py_a;
  Double I_TWOBODYOVERLAP_H3y2z_D2y_a = I_TWOBODYOVERLAP_I4y2z_Py_a+ABY*I_TWOBODYOVERLAP_H3y2z_Py_a;
  Double I_TWOBODYOVERLAP_H2y3z_D2y_a = I_TWOBODYOVERLAP_I3y3z_Py_a+ABY*I_TWOBODYOVERLAP_H2y3z_Py_a;
  Double I_TWOBODYOVERLAP_Hy4z_D2y_a = I_TWOBODYOVERLAP_I2y4z_Py_a+ABY*I_TWOBODYOVERLAP_Hy4z_Py_a;
  Double I_TWOBODYOVERLAP_H5z_D2y_a = I_TWOBODYOVERLAP_Iy5z_Py_a+ABY*I_TWOBODYOVERLAP_H5z_Py_a;
  Double I_TWOBODYOVERLAP_H5x_D2z_a = I_TWOBODYOVERLAP_I5xz_Pz_a+ABZ*I_TWOBODYOVERLAP_H5x_Pz_a;
  Double I_TWOBODYOVERLAP_H4xy_D2z_a = I_TWOBODYOVERLAP_I4xyz_Pz_a+ABZ*I_TWOBODYOVERLAP_H4xy_Pz_a;
  Double I_TWOBODYOVERLAP_H4xz_D2z_a = I_TWOBODYOVERLAP_I4x2z_Pz_a+ABZ*I_TWOBODYOVERLAP_H4xz_Pz_a;
  Double I_TWOBODYOVERLAP_H3x2y_D2z_a = I_TWOBODYOVERLAP_I3x2yz_Pz_a+ABZ*I_TWOBODYOVERLAP_H3x2y_Pz_a;
  Double I_TWOBODYOVERLAP_H3xyz_D2z_a = I_TWOBODYOVERLAP_I3xy2z_Pz_a+ABZ*I_TWOBODYOVERLAP_H3xyz_Pz_a;
  Double I_TWOBODYOVERLAP_H3x2z_D2z_a = I_TWOBODYOVERLAP_I3x3z_Pz_a+ABZ*I_TWOBODYOVERLAP_H3x2z_Pz_a;
  Double I_TWOBODYOVERLAP_H2x3y_D2z_a = I_TWOBODYOVERLAP_I2x3yz_Pz_a+ABZ*I_TWOBODYOVERLAP_H2x3y_Pz_a;
  Double I_TWOBODYOVERLAP_H2x2yz_D2z_a = I_TWOBODYOVERLAP_I2x2y2z_Pz_a+ABZ*I_TWOBODYOVERLAP_H2x2yz_Pz_a;
  Double I_TWOBODYOVERLAP_H2xy2z_D2z_a = I_TWOBODYOVERLAP_I2xy3z_Pz_a+ABZ*I_TWOBODYOVERLAP_H2xy2z_Pz_a;
  Double I_TWOBODYOVERLAP_H2x3z_D2z_a = I_TWOBODYOVERLAP_I2x4z_Pz_a+ABZ*I_TWOBODYOVERLAP_H2x3z_Pz_a;
  Double I_TWOBODYOVERLAP_Hx4y_D2z_a = I_TWOBODYOVERLAP_Ix4yz_Pz_a+ABZ*I_TWOBODYOVERLAP_Hx4y_Pz_a;
  Double I_TWOBODYOVERLAP_Hx3yz_D2z_a = I_TWOBODYOVERLAP_Ix3y2z_Pz_a+ABZ*I_TWOBODYOVERLAP_Hx3yz_Pz_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_D2z_a = I_TWOBODYOVERLAP_Ix2y3z_Pz_a+ABZ*I_TWOBODYOVERLAP_Hx2y2z_Pz_a;
  Double I_TWOBODYOVERLAP_Hxy3z_D2z_a = I_TWOBODYOVERLAP_Ixy4z_Pz_a+ABZ*I_TWOBODYOVERLAP_Hxy3z_Pz_a;
  Double I_TWOBODYOVERLAP_Hx4z_D2z_a = I_TWOBODYOVERLAP_Ix5z_Pz_a+ABZ*I_TWOBODYOVERLAP_Hx4z_Pz_a;
  Double I_TWOBODYOVERLAP_H5y_D2z_a = I_TWOBODYOVERLAP_I5yz_Pz_a+ABZ*I_TWOBODYOVERLAP_H5y_Pz_a;
  Double I_TWOBODYOVERLAP_H4yz_D2z_a = I_TWOBODYOVERLAP_I4y2z_Pz_a+ABZ*I_TWOBODYOVERLAP_H4yz_Pz_a;
  Double I_TWOBODYOVERLAP_H3y2z_D2z_a = I_TWOBODYOVERLAP_I3y3z_Pz_a+ABZ*I_TWOBODYOVERLAP_H3y2z_Pz_a;
  Double I_TWOBODYOVERLAP_H2y3z_D2z_a = I_TWOBODYOVERLAP_I2y4z_Pz_a+ABZ*I_TWOBODYOVERLAP_H2y3z_Pz_a;
  Double I_TWOBODYOVERLAP_Hy4z_D2z_a = I_TWOBODYOVERLAP_Iy5z_Pz_a+ABZ*I_TWOBODYOVERLAP_Hy4z_Pz_a;
  Double I_TWOBODYOVERLAP_H5z_D2z_a = I_TWOBODYOVERLAP_I6z_Pz_a+ABZ*I_TWOBODYOVERLAP_H5z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_K_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 3 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_L_S_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_K7x_Px_a = I_TWOBODYOVERLAP_L8x_S_a+ABX*I_TWOBODYOVERLAP_K7x_S_a;
  Double I_TWOBODYOVERLAP_K6xy_Px_a = I_TWOBODYOVERLAP_L7xy_S_a+ABX*I_TWOBODYOVERLAP_K6xy_S_a;
  Double I_TWOBODYOVERLAP_K6xz_Px_a = I_TWOBODYOVERLAP_L7xz_S_a+ABX*I_TWOBODYOVERLAP_K6xz_S_a;
  Double I_TWOBODYOVERLAP_K5x2y_Px_a = I_TWOBODYOVERLAP_L6x2y_S_a+ABX*I_TWOBODYOVERLAP_K5x2y_S_a;
  Double I_TWOBODYOVERLAP_K5xyz_Px_a = I_TWOBODYOVERLAP_L6xyz_S_a+ABX*I_TWOBODYOVERLAP_K5xyz_S_a;
  Double I_TWOBODYOVERLAP_K5x2z_Px_a = I_TWOBODYOVERLAP_L6x2z_S_a+ABX*I_TWOBODYOVERLAP_K5x2z_S_a;
  Double I_TWOBODYOVERLAP_K4x3y_Px_a = I_TWOBODYOVERLAP_L5x3y_S_a+ABX*I_TWOBODYOVERLAP_K4x3y_S_a;
  Double I_TWOBODYOVERLAP_K4x2yz_Px_a = I_TWOBODYOVERLAP_L5x2yz_S_a+ABX*I_TWOBODYOVERLAP_K4x2yz_S_a;
  Double I_TWOBODYOVERLAP_K4xy2z_Px_a = I_TWOBODYOVERLAP_L5xy2z_S_a+ABX*I_TWOBODYOVERLAP_K4xy2z_S_a;
  Double I_TWOBODYOVERLAP_K4x3z_Px_a = I_TWOBODYOVERLAP_L5x3z_S_a+ABX*I_TWOBODYOVERLAP_K4x3z_S_a;
  Double I_TWOBODYOVERLAP_K3x4y_Px_a = I_TWOBODYOVERLAP_L4x4y_S_a+ABX*I_TWOBODYOVERLAP_K3x4y_S_a;
  Double I_TWOBODYOVERLAP_K3x3yz_Px_a = I_TWOBODYOVERLAP_L4x3yz_S_a+ABX*I_TWOBODYOVERLAP_K3x3yz_S_a;
  Double I_TWOBODYOVERLAP_K3x2y2z_Px_a = I_TWOBODYOVERLAP_L4x2y2z_S_a+ABX*I_TWOBODYOVERLAP_K3x2y2z_S_a;
  Double I_TWOBODYOVERLAP_K3xy3z_Px_a = I_TWOBODYOVERLAP_L4xy3z_S_a+ABX*I_TWOBODYOVERLAP_K3xy3z_S_a;
  Double I_TWOBODYOVERLAP_K3x4z_Px_a = I_TWOBODYOVERLAP_L4x4z_S_a+ABX*I_TWOBODYOVERLAP_K3x4z_S_a;
  Double I_TWOBODYOVERLAP_K2x5y_Px_a = I_TWOBODYOVERLAP_L3x5y_S_a+ABX*I_TWOBODYOVERLAP_K2x5y_S_a;
  Double I_TWOBODYOVERLAP_K2x4yz_Px_a = I_TWOBODYOVERLAP_L3x4yz_S_a+ABX*I_TWOBODYOVERLAP_K2x4yz_S_a;
  Double I_TWOBODYOVERLAP_K2x3y2z_Px_a = I_TWOBODYOVERLAP_L3x3y2z_S_a+ABX*I_TWOBODYOVERLAP_K2x3y2z_S_a;
  Double I_TWOBODYOVERLAP_K2x2y3z_Px_a = I_TWOBODYOVERLAP_L3x2y3z_S_a+ABX*I_TWOBODYOVERLAP_K2x2y3z_S_a;
  Double I_TWOBODYOVERLAP_K2xy4z_Px_a = I_TWOBODYOVERLAP_L3xy4z_S_a+ABX*I_TWOBODYOVERLAP_K2xy4z_S_a;
  Double I_TWOBODYOVERLAP_K2x5z_Px_a = I_TWOBODYOVERLAP_L3x5z_S_a+ABX*I_TWOBODYOVERLAP_K2x5z_S_a;
  Double I_TWOBODYOVERLAP_Kx6y_Px_a = I_TWOBODYOVERLAP_L2x6y_S_a+ABX*I_TWOBODYOVERLAP_Kx6y_S_a;
  Double I_TWOBODYOVERLAP_Kx5yz_Px_a = I_TWOBODYOVERLAP_L2x5yz_S_a+ABX*I_TWOBODYOVERLAP_Kx5yz_S_a;
  Double I_TWOBODYOVERLAP_Kx4y2z_Px_a = I_TWOBODYOVERLAP_L2x4y2z_S_a+ABX*I_TWOBODYOVERLAP_Kx4y2z_S_a;
  Double I_TWOBODYOVERLAP_Kx3y3z_Px_a = I_TWOBODYOVERLAP_L2x3y3z_S_a+ABX*I_TWOBODYOVERLAP_Kx3y3z_S_a;
  Double I_TWOBODYOVERLAP_Kx2y4z_Px_a = I_TWOBODYOVERLAP_L2x2y4z_S_a+ABX*I_TWOBODYOVERLAP_Kx2y4z_S_a;
  Double I_TWOBODYOVERLAP_Kxy5z_Px_a = I_TWOBODYOVERLAP_L2xy5z_S_a+ABX*I_TWOBODYOVERLAP_Kxy5z_S_a;
  Double I_TWOBODYOVERLAP_Kx6z_Px_a = I_TWOBODYOVERLAP_L2x6z_S_a+ABX*I_TWOBODYOVERLAP_Kx6z_S_a;
  Double I_TWOBODYOVERLAP_K7y_Px_a = I_TWOBODYOVERLAP_Lx7y_S_a+ABX*I_TWOBODYOVERLAP_K7y_S_a;
  Double I_TWOBODYOVERLAP_K6yz_Px_a = I_TWOBODYOVERLAP_Lx6yz_S_a+ABX*I_TWOBODYOVERLAP_K6yz_S_a;
  Double I_TWOBODYOVERLAP_K5y2z_Px_a = I_TWOBODYOVERLAP_Lx5y2z_S_a+ABX*I_TWOBODYOVERLAP_K5y2z_S_a;
  Double I_TWOBODYOVERLAP_K4y3z_Px_a = I_TWOBODYOVERLAP_Lx4y3z_S_a+ABX*I_TWOBODYOVERLAP_K4y3z_S_a;
  Double I_TWOBODYOVERLAP_K3y4z_Px_a = I_TWOBODYOVERLAP_Lx3y4z_S_a+ABX*I_TWOBODYOVERLAP_K3y4z_S_a;
  Double I_TWOBODYOVERLAP_K2y5z_Px_a = I_TWOBODYOVERLAP_Lx2y5z_S_a+ABX*I_TWOBODYOVERLAP_K2y5z_S_a;
  Double I_TWOBODYOVERLAP_Ky6z_Px_a = I_TWOBODYOVERLAP_Lxy6z_S_a+ABX*I_TWOBODYOVERLAP_Ky6z_S_a;
  Double I_TWOBODYOVERLAP_K7z_Px_a = I_TWOBODYOVERLAP_Lx7z_S_a+ABX*I_TWOBODYOVERLAP_K7z_S_a;
  Double I_TWOBODYOVERLAP_K6xy_Py_a = I_TWOBODYOVERLAP_L6x2y_S_a+ABY*I_TWOBODYOVERLAP_K6xy_S_a;
  Double I_TWOBODYOVERLAP_K6xz_Py_a = I_TWOBODYOVERLAP_L6xyz_S_a+ABY*I_TWOBODYOVERLAP_K6xz_S_a;
  Double I_TWOBODYOVERLAP_K5x2y_Py_a = I_TWOBODYOVERLAP_L5x3y_S_a+ABY*I_TWOBODYOVERLAP_K5x2y_S_a;
  Double I_TWOBODYOVERLAP_K5xyz_Py_a = I_TWOBODYOVERLAP_L5x2yz_S_a+ABY*I_TWOBODYOVERLAP_K5xyz_S_a;
  Double I_TWOBODYOVERLAP_K5x2z_Py_a = I_TWOBODYOVERLAP_L5xy2z_S_a+ABY*I_TWOBODYOVERLAP_K5x2z_S_a;
  Double I_TWOBODYOVERLAP_K4x3y_Py_a = I_TWOBODYOVERLAP_L4x4y_S_a+ABY*I_TWOBODYOVERLAP_K4x3y_S_a;
  Double I_TWOBODYOVERLAP_K4x2yz_Py_a = I_TWOBODYOVERLAP_L4x3yz_S_a+ABY*I_TWOBODYOVERLAP_K4x2yz_S_a;
  Double I_TWOBODYOVERLAP_K4xy2z_Py_a = I_TWOBODYOVERLAP_L4x2y2z_S_a+ABY*I_TWOBODYOVERLAP_K4xy2z_S_a;
  Double I_TWOBODYOVERLAP_K4x3z_Py_a = I_TWOBODYOVERLAP_L4xy3z_S_a+ABY*I_TWOBODYOVERLAP_K4x3z_S_a;
  Double I_TWOBODYOVERLAP_K3x4y_Py_a = I_TWOBODYOVERLAP_L3x5y_S_a+ABY*I_TWOBODYOVERLAP_K3x4y_S_a;
  Double I_TWOBODYOVERLAP_K3x3yz_Py_a = I_TWOBODYOVERLAP_L3x4yz_S_a+ABY*I_TWOBODYOVERLAP_K3x3yz_S_a;
  Double I_TWOBODYOVERLAP_K3x2y2z_Py_a = I_TWOBODYOVERLAP_L3x3y2z_S_a+ABY*I_TWOBODYOVERLAP_K3x2y2z_S_a;
  Double I_TWOBODYOVERLAP_K3xy3z_Py_a = I_TWOBODYOVERLAP_L3x2y3z_S_a+ABY*I_TWOBODYOVERLAP_K3xy3z_S_a;
  Double I_TWOBODYOVERLAP_K3x4z_Py_a = I_TWOBODYOVERLAP_L3xy4z_S_a+ABY*I_TWOBODYOVERLAP_K3x4z_S_a;
  Double I_TWOBODYOVERLAP_K2x5y_Py_a = I_TWOBODYOVERLAP_L2x6y_S_a+ABY*I_TWOBODYOVERLAP_K2x5y_S_a;
  Double I_TWOBODYOVERLAP_K2x4yz_Py_a = I_TWOBODYOVERLAP_L2x5yz_S_a+ABY*I_TWOBODYOVERLAP_K2x4yz_S_a;
  Double I_TWOBODYOVERLAP_K2x3y2z_Py_a = I_TWOBODYOVERLAP_L2x4y2z_S_a+ABY*I_TWOBODYOVERLAP_K2x3y2z_S_a;
  Double I_TWOBODYOVERLAP_K2x2y3z_Py_a = I_TWOBODYOVERLAP_L2x3y3z_S_a+ABY*I_TWOBODYOVERLAP_K2x2y3z_S_a;
  Double I_TWOBODYOVERLAP_K2xy4z_Py_a = I_TWOBODYOVERLAP_L2x2y4z_S_a+ABY*I_TWOBODYOVERLAP_K2xy4z_S_a;
  Double I_TWOBODYOVERLAP_K2x5z_Py_a = I_TWOBODYOVERLAP_L2xy5z_S_a+ABY*I_TWOBODYOVERLAP_K2x5z_S_a;
  Double I_TWOBODYOVERLAP_Kx6y_Py_a = I_TWOBODYOVERLAP_Lx7y_S_a+ABY*I_TWOBODYOVERLAP_Kx6y_S_a;
  Double I_TWOBODYOVERLAP_Kx5yz_Py_a = I_TWOBODYOVERLAP_Lx6yz_S_a+ABY*I_TWOBODYOVERLAP_Kx5yz_S_a;
  Double I_TWOBODYOVERLAP_Kx4y2z_Py_a = I_TWOBODYOVERLAP_Lx5y2z_S_a+ABY*I_TWOBODYOVERLAP_Kx4y2z_S_a;
  Double I_TWOBODYOVERLAP_Kx3y3z_Py_a = I_TWOBODYOVERLAP_Lx4y3z_S_a+ABY*I_TWOBODYOVERLAP_Kx3y3z_S_a;
  Double I_TWOBODYOVERLAP_Kx2y4z_Py_a = I_TWOBODYOVERLAP_Lx3y4z_S_a+ABY*I_TWOBODYOVERLAP_Kx2y4z_S_a;
  Double I_TWOBODYOVERLAP_Kxy5z_Py_a = I_TWOBODYOVERLAP_Lx2y5z_S_a+ABY*I_TWOBODYOVERLAP_Kxy5z_S_a;
  Double I_TWOBODYOVERLAP_Kx6z_Py_a = I_TWOBODYOVERLAP_Lxy6z_S_a+ABY*I_TWOBODYOVERLAP_Kx6z_S_a;
  Double I_TWOBODYOVERLAP_K7y_Py_a = I_TWOBODYOVERLAP_L8y_S_a+ABY*I_TWOBODYOVERLAP_K7y_S_a;
  Double I_TWOBODYOVERLAP_K6yz_Py_a = I_TWOBODYOVERLAP_L7yz_S_a+ABY*I_TWOBODYOVERLAP_K6yz_S_a;
  Double I_TWOBODYOVERLAP_K5y2z_Py_a = I_TWOBODYOVERLAP_L6y2z_S_a+ABY*I_TWOBODYOVERLAP_K5y2z_S_a;
  Double I_TWOBODYOVERLAP_K4y3z_Py_a = I_TWOBODYOVERLAP_L5y3z_S_a+ABY*I_TWOBODYOVERLAP_K4y3z_S_a;
  Double I_TWOBODYOVERLAP_K3y4z_Py_a = I_TWOBODYOVERLAP_L4y4z_S_a+ABY*I_TWOBODYOVERLAP_K3y4z_S_a;
  Double I_TWOBODYOVERLAP_K2y5z_Py_a = I_TWOBODYOVERLAP_L3y5z_S_a+ABY*I_TWOBODYOVERLAP_K2y5z_S_a;
  Double I_TWOBODYOVERLAP_Ky6z_Py_a = I_TWOBODYOVERLAP_L2y6z_S_a+ABY*I_TWOBODYOVERLAP_Ky6z_S_a;
  Double I_TWOBODYOVERLAP_K7z_Py_a = I_TWOBODYOVERLAP_Ly7z_S_a+ABY*I_TWOBODYOVERLAP_K7z_S_a;
  Double I_TWOBODYOVERLAP_K6xy_Pz_a = I_TWOBODYOVERLAP_L6xyz_S_a+ABZ*I_TWOBODYOVERLAP_K6xy_S_a;
  Double I_TWOBODYOVERLAP_K6xz_Pz_a = I_TWOBODYOVERLAP_L6x2z_S_a+ABZ*I_TWOBODYOVERLAP_K6xz_S_a;
  Double I_TWOBODYOVERLAP_K5x2y_Pz_a = I_TWOBODYOVERLAP_L5x2yz_S_a+ABZ*I_TWOBODYOVERLAP_K5x2y_S_a;
  Double I_TWOBODYOVERLAP_K5xyz_Pz_a = I_TWOBODYOVERLAP_L5xy2z_S_a+ABZ*I_TWOBODYOVERLAP_K5xyz_S_a;
  Double I_TWOBODYOVERLAP_K5x2z_Pz_a = I_TWOBODYOVERLAP_L5x3z_S_a+ABZ*I_TWOBODYOVERLAP_K5x2z_S_a;
  Double I_TWOBODYOVERLAP_K4x3y_Pz_a = I_TWOBODYOVERLAP_L4x3yz_S_a+ABZ*I_TWOBODYOVERLAP_K4x3y_S_a;
  Double I_TWOBODYOVERLAP_K4x2yz_Pz_a = I_TWOBODYOVERLAP_L4x2y2z_S_a+ABZ*I_TWOBODYOVERLAP_K4x2yz_S_a;
  Double I_TWOBODYOVERLAP_K4xy2z_Pz_a = I_TWOBODYOVERLAP_L4xy3z_S_a+ABZ*I_TWOBODYOVERLAP_K4xy2z_S_a;
  Double I_TWOBODYOVERLAP_K4x3z_Pz_a = I_TWOBODYOVERLAP_L4x4z_S_a+ABZ*I_TWOBODYOVERLAP_K4x3z_S_a;
  Double I_TWOBODYOVERLAP_K3x4y_Pz_a = I_TWOBODYOVERLAP_L3x4yz_S_a+ABZ*I_TWOBODYOVERLAP_K3x4y_S_a;
  Double I_TWOBODYOVERLAP_K3x3yz_Pz_a = I_TWOBODYOVERLAP_L3x3y2z_S_a+ABZ*I_TWOBODYOVERLAP_K3x3yz_S_a;
  Double I_TWOBODYOVERLAP_K3x2y2z_Pz_a = I_TWOBODYOVERLAP_L3x2y3z_S_a+ABZ*I_TWOBODYOVERLAP_K3x2y2z_S_a;
  Double I_TWOBODYOVERLAP_K3xy3z_Pz_a = I_TWOBODYOVERLAP_L3xy4z_S_a+ABZ*I_TWOBODYOVERLAP_K3xy3z_S_a;
  Double I_TWOBODYOVERLAP_K3x4z_Pz_a = I_TWOBODYOVERLAP_L3x5z_S_a+ABZ*I_TWOBODYOVERLAP_K3x4z_S_a;
  Double I_TWOBODYOVERLAP_K2x5y_Pz_a = I_TWOBODYOVERLAP_L2x5yz_S_a+ABZ*I_TWOBODYOVERLAP_K2x5y_S_a;
  Double I_TWOBODYOVERLAP_K2x4yz_Pz_a = I_TWOBODYOVERLAP_L2x4y2z_S_a+ABZ*I_TWOBODYOVERLAP_K2x4yz_S_a;
  Double I_TWOBODYOVERLAP_K2x3y2z_Pz_a = I_TWOBODYOVERLAP_L2x3y3z_S_a+ABZ*I_TWOBODYOVERLAP_K2x3y2z_S_a;
  Double I_TWOBODYOVERLAP_K2x2y3z_Pz_a = I_TWOBODYOVERLAP_L2x2y4z_S_a+ABZ*I_TWOBODYOVERLAP_K2x2y3z_S_a;
  Double I_TWOBODYOVERLAP_K2xy4z_Pz_a = I_TWOBODYOVERLAP_L2xy5z_S_a+ABZ*I_TWOBODYOVERLAP_K2xy4z_S_a;
  Double I_TWOBODYOVERLAP_K2x5z_Pz_a = I_TWOBODYOVERLAP_L2x6z_S_a+ABZ*I_TWOBODYOVERLAP_K2x5z_S_a;
  Double I_TWOBODYOVERLAP_Kx6y_Pz_a = I_TWOBODYOVERLAP_Lx6yz_S_a+ABZ*I_TWOBODYOVERLAP_Kx6y_S_a;
  Double I_TWOBODYOVERLAP_Kx5yz_Pz_a = I_TWOBODYOVERLAP_Lx5y2z_S_a+ABZ*I_TWOBODYOVERLAP_Kx5yz_S_a;
  Double I_TWOBODYOVERLAP_Kx4y2z_Pz_a = I_TWOBODYOVERLAP_Lx4y3z_S_a+ABZ*I_TWOBODYOVERLAP_Kx4y2z_S_a;
  Double I_TWOBODYOVERLAP_Kx3y3z_Pz_a = I_TWOBODYOVERLAP_Lx3y4z_S_a+ABZ*I_TWOBODYOVERLAP_Kx3y3z_S_a;
  Double I_TWOBODYOVERLAP_Kx2y4z_Pz_a = I_TWOBODYOVERLAP_Lx2y5z_S_a+ABZ*I_TWOBODYOVERLAP_Kx2y4z_S_a;
  Double I_TWOBODYOVERLAP_Kxy5z_Pz_a = I_TWOBODYOVERLAP_Lxy6z_S_a+ABZ*I_TWOBODYOVERLAP_Kxy5z_S_a;
  Double I_TWOBODYOVERLAP_Kx6z_Pz_a = I_TWOBODYOVERLAP_Lx7z_S_a+ABZ*I_TWOBODYOVERLAP_Kx6z_S_a;
  Double I_TWOBODYOVERLAP_K6yz_Pz_a = I_TWOBODYOVERLAP_L6y2z_S_a+ABZ*I_TWOBODYOVERLAP_K6yz_S_a;
  Double I_TWOBODYOVERLAP_K5y2z_Pz_a = I_TWOBODYOVERLAP_L5y3z_S_a+ABZ*I_TWOBODYOVERLAP_K5y2z_S_a;
  Double I_TWOBODYOVERLAP_K4y3z_Pz_a = I_TWOBODYOVERLAP_L4y4z_S_a+ABZ*I_TWOBODYOVERLAP_K4y3z_S_a;
  Double I_TWOBODYOVERLAP_K3y4z_Pz_a = I_TWOBODYOVERLAP_L3y5z_S_a+ABZ*I_TWOBODYOVERLAP_K3y4z_S_a;
  Double I_TWOBODYOVERLAP_K2y5z_Pz_a = I_TWOBODYOVERLAP_L2y6z_S_a+ABZ*I_TWOBODYOVERLAP_K2y5z_S_a;
  Double I_TWOBODYOVERLAP_Ky6z_Pz_a = I_TWOBODYOVERLAP_Ly7z_S_a+ABZ*I_TWOBODYOVERLAP_Ky6z_S_a;
  Double I_TWOBODYOVERLAP_K7z_Pz_a = I_TWOBODYOVERLAP_L8z_S_a+ABZ*I_TWOBODYOVERLAP_K7z_S_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_I_D_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 84 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_P_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_P_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_I6x_D2x_a = I_TWOBODYOVERLAP_K7x_Px_a+ABX*I_TWOBODYOVERLAP_I6x_Px_a;
  Double I_TWOBODYOVERLAP_I5xy_D2x_a = I_TWOBODYOVERLAP_K6xy_Px_a+ABX*I_TWOBODYOVERLAP_I5xy_Px_a;
  Double I_TWOBODYOVERLAP_I5xz_D2x_a = I_TWOBODYOVERLAP_K6xz_Px_a+ABX*I_TWOBODYOVERLAP_I5xz_Px_a;
  Double I_TWOBODYOVERLAP_I4x2y_D2x_a = I_TWOBODYOVERLAP_K5x2y_Px_a+ABX*I_TWOBODYOVERLAP_I4x2y_Px_a;
  Double I_TWOBODYOVERLAP_I4xyz_D2x_a = I_TWOBODYOVERLAP_K5xyz_Px_a+ABX*I_TWOBODYOVERLAP_I4xyz_Px_a;
  Double I_TWOBODYOVERLAP_I4x2z_D2x_a = I_TWOBODYOVERLAP_K5x2z_Px_a+ABX*I_TWOBODYOVERLAP_I4x2z_Px_a;
  Double I_TWOBODYOVERLAP_I3x3y_D2x_a = I_TWOBODYOVERLAP_K4x3y_Px_a+ABX*I_TWOBODYOVERLAP_I3x3y_Px_a;
  Double I_TWOBODYOVERLAP_I3x2yz_D2x_a = I_TWOBODYOVERLAP_K4x2yz_Px_a+ABX*I_TWOBODYOVERLAP_I3x2yz_Px_a;
  Double I_TWOBODYOVERLAP_I3xy2z_D2x_a = I_TWOBODYOVERLAP_K4xy2z_Px_a+ABX*I_TWOBODYOVERLAP_I3xy2z_Px_a;
  Double I_TWOBODYOVERLAP_I3x3z_D2x_a = I_TWOBODYOVERLAP_K4x3z_Px_a+ABX*I_TWOBODYOVERLAP_I3x3z_Px_a;
  Double I_TWOBODYOVERLAP_I2x4y_D2x_a = I_TWOBODYOVERLAP_K3x4y_Px_a+ABX*I_TWOBODYOVERLAP_I2x4y_Px_a;
  Double I_TWOBODYOVERLAP_I2x3yz_D2x_a = I_TWOBODYOVERLAP_K3x3yz_Px_a+ABX*I_TWOBODYOVERLAP_I2x3yz_Px_a;
  Double I_TWOBODYOVERLAP_I2x2y2z_D2x_a = I_TWOBODYOVERLAP_K3x2y2z_Px_a+ABX*I_TWOBODYOVERLAP_I2x2y2z_Px_a;
  Double I_TWOBODYOVERLAP_I2xy3z_D2x_a = I_TWOBODYOVERLAP_K3xy3z_Px_a+ABX*I_TWOBODYOVERLAP_I2xy3z_Px_a;
  Double I_TWOBODYOVERLAP_I2x4z_D2x_a = I_TWOBODYOVERLAP_K3x4z_Px_a+ABX*I_TWOBODYOVERLAP_I2x4z_Px_a;
  Double I_TWOBODYOVERLAP_Ix5y_D2x_a = I_TWOBODYOVERLAP_K2x5y_Px_a+ABX*I_TWOBODYOVERLAP_Ix5y_Px_a;
  Double I_TWOBODYOVERLAP_Ix4yz_D2x_a = I_TWOBODYOVERLAP_K2x4yz_Px_a+ABX*I_TWOBODYOVERLAP_Ix4yz_Px_a;
  Double I_TWOBODYOVERLAP_Ix3y2z_D2x_a = I_TWOBODYOVERLAP_K2x3y2z_Px_a+ABX*I_TWOBODYOVERLAP_Ix3y2z_Px_a;
  Double I_TWOBODYOVERLAP_Ix2y3z_D2x_a = I_TWOBODYOVERLAP_K2x2y3z_Px_a+ABX*I_TWOBODYOVERLAP_Ix2y3z_Px_a;
  Double I_TWOBODYOVERLAP_Ixy4z_D2x_a = I_TWOBODYOVERLAP_K2xy4z_Px_a+ABX*I_TWOBODYOVERLAP_Ixy4z_Px_a;
  Double I_TWOBODYOVERLAP_Ix5z_D2x_a = I_TWOBODYOVERLAP_K2x5z_Px_a+ABX*I_TWOBODYOVERLAP_Ix5z_Px_a;
  Double I_TWOBODYOVERLAP_I6y_D2x_a = I_TWOBODYOVERLAP_Kx6y_Px_a+ABX*I_TWOBODYOVERLAP_I6y_Px_a;
  Double I_TWOBODYOVERLAP_I5yz_D2x_a = I_TWOBODYOVERLAP_Kx5yz_Px_a+ABX*I_TWOBODYOVERLAP_I5yz_Px_a;
  Double I_TWOBODYOVERLAP_I4y2z_D2x_a = I_TWOBODYOVERLAP_Kx4y2z_Px_a+ABX*I_TWOBODYOVERLAP_I4y2z_Px_a;
  Double I_TWOBODYOVERLAP_I3y3z_D2x_a = I_TWOBODYOVERLAP_Kx3y3z_Px_a+ABX*I_TWOBODYOVERLAP_I3y3z_Px_a;
  Double I_TWOBODYOVERLAP_I2y4z_D2x_a = I_TWOBODYOVERLAP_Kx2y4z_Px_a+ABX*I_TWOBODYOVERLAP_I2y4z_Px_a;
  Double I_TWOBODYOVERLAP_Iy5z_D2x_a = I_TWOBODYOVERLAP_Kxy5z_Px_a+ABX*I_TWOBODYOVERLAP_Iy5z_Px_a;
  Double I_TWOBODYOVERLAP_I6z_D2x_a = I_TWOBODYOVERLAP_Kx6z_Px_a+ABX*I_TWOBODYOVERLAP_I6z_Px_a;
  Double I_TWOBODYOVERLAP_I6x_D2y_a = I_TWOBODYOVERLAP_K6xy_Py_a+ABY*I_TWOBODYOVERLAP_I6x_Py_a;
  Double I_TWOBODYOVERLAP_I5xy_D2y_a = I_TWOBODYOVERLAP_K5x2y_Py_a+ABY*I_TWOBODYOVERLAP_I5xy_Py_a;
  Double I_TWOBODYOVERLAP_I5xz_D2y_a = I_TWOBODYOVERLAP_K5xyz_Py_a+ABY*I_TWOBODYOVERLAP_I5xz_Py_a;
  Double I_TWOBODYOVERLAP_I4x2y_D2y_a = I_TWOBODYOVERLAP_K4x3y_Py_a+ABY*I_TWOBODYOVERLAP_I4x2y_Py_a;
  Double I_TWOBODYOVERLAP_I4xyz_D2y_a = I_TWOBODYOVERLAP_K4x2yz_Py_a+ABY*I_TWOBODYOVERLAP_I4xyz_Py_a;
  Double I_TWOBODYOVERLAP_I4x2z_D2y_a = I_TWOBODYOVERLAP_K4xy2z_Py_a+ABY*I_TWOBODYOVERLAP_I4x2z_Py_a;
  Double I_TWOBODYOVERLAP_I3x3y_D2y_a = I_TWOBODYOVERLAP_K3x4y_Py_a+ABY*I_TWOBODYOVERLAP_I3x3y_Py_a;
  Double I_TWOBODYOVERLAP_I3x2yz_D2y_a = I_TWOBODYOVERLAP_K3x3yz_Py_a+ABY*I_TWOBODYOVERLAP_I3x2yz_Py_a;
  Double I_TWOBODYOVERLAP_I3xy2z_D2y_a = I_TWOBODYOVERLAP_K3x2y2z_Py_a+ABY*I_TWOBODYOVERLAP_I3xy2z_Py_a;
  Double I_TWOBODYOVERLAP_I3x3z_D2y_a = I_TWOBODYOVERLAP_K3xy3z_Py_a+ABY*I_TWOBODYOVERLAP_I3x3z_Py_a;
  Double I_TWOBODYOVERLAP_I2x4y_D2y_a = I_TWOBODYOVERLAP_K2x5y_Py_a+ABY*I_TWOBODYOVERLAP_I2x4y_Py_a;
  Double I_TWOBODYOVERLAP_I2x3yz_D2y_a = I_TWOBODYOVERLAP_K2x4yz_Py_a+ABY*I_TWOBODYOVERLAP_I2x3yz_Py_a;
  Double I_TWOBODYOVERLAP_I2x2y2z_D2y_a = I_TWOBODYOVERLAP_K2x3y2z_Py_a+ABY*I_TWOBODYOVERLAP_I2x2y2z_Py_a;
  Double I_TWOBODYOVERLAP_I2xy3z_D2y_a = I_TWOBODYOVERLAP_K2x2y3z_Py_a+ABY*I_TWOBODYOVERLAP_I2xy3z_Py_a;
  Double I_TWOBODYOVERLAP_I2x4z_D2y_a = I_TWOBODYOVERLAP_K2xy4z_Py_a+ABY*I_TWOBODYOVERLAP_I2x4z_Py_a;
  Double I_TWOBODYOVERLAP_Ix5y_D2y_a = I_TWOBODYOVERLAP_Kx6y_Py_a+ABY*I_TWOBODYOVERLAP_Ix5y_Py_a;
  Double I_TWOBODYOVERLAP_Ix4yz_D2y_a = I_TWOBODYOVERLAP_Kx5yz_Py_a+ABY*I_TWOBODYOVERLAP_Ix4yz_Py_a;
  Double I_TWOBODYOVERLAP_Ix3y2z_D2y_a = I_TWOBODYOVERLAP_Kx4y2z_Py_a+ABY*I_TWOBODYOVERLAP_Ix3y2z_Py_a;
  Double I_TWOBODYOVERLAP_Ix2y3z_D2y_a = I_TWOBODYOVERLAP_Kx3y3z_Py_a+ABY*I_TWOBODYOVERLAP_Ix2y3z_Py_a;
  Double I_TWOBODYOVERLAP_Ixy4z_D2y_a = I_TWOBODYOVERLAP_Kx2y4z_Py_a+ABY*I_TWOBODYOVERLAP_Ixy4z_Py_a;
  Double I_TWOBODYOVERLAP_Ix5z_D2y_a = I_TWOBODYOVERLAP_Kxy5z_Py_a+ABY*I_TWOBODYOVERLAP_Ix5z_Py_a;
  Double I_TWOBODYOVERLAP_I6y_D2y_a = I_TWOBODYOVERLAP_K7y_Py_a+ABY*I_TWOBODYOVERLAP_I6y_Py_a;
  Double I_TWOBODYOVERLAP_I5yz_D2y_a = I_TWOBODYOVERLAP_K6yz_Py_a+ABY*I_TWOBODYOVERLAP_I5yz_Py_a;
  Double I_TWOBODYOVERLAP_I4y2z_D2y_a = I_TWOBODYOVERLAP_K5y2z_Py_a+ABY*I_TWOBODYOVERLAP_I4y2z_Py_a;
  Double I_TWOBODYOVERLAP_I3y3z_D2y_a = I_TWOBODYOVERLAP_K4y3z_Py_a+ABY*I_TWOBODYOVERLAP_I3y3z_Py_a;
  Double I_TWOBODYOVERLAP_I2y4z_D2y_a = I_TWOBODYOVERLAP_K3y4z_Py_a+ABY*I_TWOBODYOVERLAP_I2y4z_Py_a;
  Double I_TWOBODYOVERLAP_Iy5z_D2y_a = I_TWOBODYOVERLAP_K2y5z_Py_a+ABY*I_TWOBODYOVERLAP_Iy5z_Py_a;
  Double I_TWOBODYOVERLAP_I6z_D2y_a = I_TWOBODYOVERLAP_Ky6z_Py_a+ABY*I_TWOBODYOVERLAP_I6z_Py_a;
  Double I_TWOBODYOVERLAP_I6x_D2z_a = I_TWOBODYOVERLAP_K6xz_Pz_a+ABZ*I_TWOBODYOVERLAP_I6x_Pz_a;
  Double I_TWOBODYOVERLAP_I5xy_D2z_a = I_TWOBODYOVERLAP_K5xyz_Pz_a+ABZ*I_TWOBODYOVERLAP_I5xy_Pz_a;
  Double I_TWOBODYOVERLAP_I5xz_D2z_a = I_TWOBODYOVERLAP_K5x2z_Pz_a+ABZ*I_TWOBODYOVERLAP_I5xz_Pz_a;
  Double I_TWOBODYOVERLAP_I4x2y_D2z_a = I_TWOBODYOVERLAP_K4x2yz_Pz_a+ABZ*I_TWOBODYOVERLAP_I4x2y_Pz_a;
  Double I_TWOBODYOVERLAP_I4xyz_D2z_a = I_TWOBODYOVERLAP_K4xy2z_Pz_a+ABZ*I_TWOBODYOVERLAP_I4xyz_Pz_a;
  Double I_TWOBODYOVERLAP_I4x2z_D2z_a = I_TWOBODYOVERLAP_K4x3z_Pz_a+ABZ*I_TWOBODYOVERLAP_I4x2z_Pz_a;
  Double I_TWOBODYOVERLAP_I3x3y_D2z_a = I_TWOBODYOVERLAP_K3x3yz_Pz_a+ABZ*I_TWOBODYOVERLAP_I3x3y_Pz_a;
  Double I_TWOBODYOVERLAP_I3x2yz_D2z_a = I_TWOBODYOVERLAP_K3x2y2z_Pz_a+ABZ*I_TWOBODYOVERLAP_I3x2yz_Pz_a;
  Double I_TWOBODYOVERLAP_I3xy2z_D2z_a = I_TWOBODYOVERLAP_K3xy3z_Pz_a+ABZ*I_TWOBODYOVERLAP_I3xy2z_Pz_a;
  Double I_TWOBODYOVERLAP_I3x3z_D2z_a = I_TWOBODYOVERLAP_K3x4z_Pz_a+ABZ*I_TWOBODYOVERLAP_I3x3z_Pz_a;
  Double I_TWOBODYOVERLAP_I2x4y_D2z_a = I_TWOBODYOVERLAP_K2x4yz_Pz_a+ABZ*I_TWOBODYOVERLAP_I2x4y_Pz_a;
  Double I_TWOBODYOVERLAP_I2x3yz_D2z_a = I_TWOBODYOVERLAP_K2x3y2z_Pz_a+ABZ*I_TWOBODYOVERLAP_I2x3yz_Pz_a;
  Double I_TWOBODYOVERLAP_I2x2y2z_D2z_a = I_TWOBODYOVERLAP_K2x2y3z_Pz_a+ABZ*I_TWOBODYOVERLAP_I2x2y2z_Pz_a;
  Double I_TWOBODYOVERLAP_I2xy3z_D2z_a = I_TWOBODYOVERLAP_K2xy4z_Pz_a+ABZ*I_TWOBODYOVERLAP_I2xy3z_Pz_a;
  Double I_TWOBODYOVERLAP_I2x4z_D2z_a = I_TWOBODYOVERLAP_K2x5z_Pz_a+ABZ*I_TWOBODYOVERLAP_I2x4z_Pz_a;
  Double I_TWOBODYOVERLAP_Ix5y_D2z_a = I_TWOBODYOVERLAP_Kx5yz_Pz_a+ABZ*I_TWOBODYOVERLAP_Ix5y_Pz_a;
  Double I_TWOBODYOVERLAP_Ix4yz_D2z_a = I_TWOBODYOVERLAP_Kx4y2z_Pz_a+ABZ*I_TWOBODYOVERLAP_Ix4yz_Pz_a;
  Double I_TWOBODYOVERLAP_Ix3y2z_D2z_a = I_TWOBODYOVERLAP_Kx3y3z_Pz_a+ABZ*I_TWOBODYOVERLAP_Ix3y2z_Pz_a;
  Double I_TWOBODYOVERLAP_Ix2y3z_D2z_a = I_TWOBODYOVERLAP_Kx2y4z_Pz_a+ABZ*I_TWOBODYOVERLAP_Ix2y3z_Pz_a;
  Double I_TWOBODYOVERLAP_Ixy4z_D2z_a = I_TWOBODYOVERLAP_Kxy5z_Pz_a+ABZ*I_TWOBODYOVERLAP_Ixy4z_Pz_a;
  Double I_TWOBODYOVERLAP_Ix5z_D2z_a = I_TWOBODYOVERLAP_Kx6z_Pz_a+ABZ*I_TWOBODYOVERLAP_Ix5z_Pz_a;
  Double I_TWOBODYOVERLAP_I6y_D2z_a = I_TWOBODYOVERLAP_K6yz_Pz_a+ABZ*I_TWOBODYOVERLAP_I6y_Pz_a;
  Double I_TWOBODYOVERLAP_I5yz_D2z_a = I_TWOBODYOVERLAP_K5y2z_Pz_a+ABZ*I_TWOBODYOVERLAP_I5yz_Pz_a;
  Double I_TWOBODYOVERLAP_I4y2z_D2z_a = I_TWOBODYOVERLAP_K4y3z_Pz_a+ABZ*I_TWOBODYOVERLAP_I4y2z_Pz_a;
  Double I_TWOBODYOVERLAP_I3y3z_D2z_a = I_TWOBODYOVERLAP_K3y4z_Pz_a+ABZ*I_TWOBODYOVERLAP_I3y3z_Pz_a;
  Double I_TWOBODYOVERLAP_I2y4z_D2z_a = I_TWOBODYOVERLAP_K2y5z_Pz_a+ABZ*I_TWOBODYOVERLAP_I2y4z_Pz_a;
  Double I_TWOBODYOVERLAP_Iy5z_D2z_a = I_TWOBODYOVERLAP_Ky6z_Pz_a+ABZ*I_TWOBODYOVERLAP_Iy5z_Pz_a;
  Double I_TWOBODYOVERLAP_I6z_D2z_a = I_TWOBODYOVERLAP_K7z_Pz_a+ABZ*I_TWOBODYOVERLAP_I6z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_F_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 42 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_D_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_D_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_H5x_F3x_a = I_TWOBODYOVERLAP_I6x_D2x_a+ABX*I_TWOBODYOVERLAP_H5x_D2x_a;
  Double I_TWOBODYOVERLAP_H4xy_F3x_a = I_TWOBODYOVERLAP_I5xy_D2x_a+ABX*I_TWOBODYOVERLAP_H4xy_D2x_a;
  Double I_TWOBODYOVERLAP_H4xz_F3x_a = I_TWOBODYOVERLAP_I5xz_D2x_a+ABX*I_TWOBODYOVERLAP_H4xz_D2x_a;
  Double I_TWOBODYOVERLAP_H3x2y_F3x_a = I_TWOBODYOVERLAP_I4x2y_D2x_a+ABX*I_TWOBODYOVERLAP_H3x2y_D2x_a;
  Double I_TWOBODYOVERLAP_H3xyz_F3x_a = I_TWOBODYOVERLAP_I4xyz_D2x_a+ABX*I_TWOBODYOVERLAP_H3xyz_D2x_a;
  Double I_TWOBODYOVERLAP_H3x2z_F3x_a = I_TWOBODYOVERLAP_I4x2z_D2x_a+ABX*I_TWOBODYOVERLAP_H3x2z_D2x_a;
  Double I_TWOBODYOVERLAP_H2x3y_F3x_a = I_TWOBODYOVERLAP_I3x3y_D2x_a+ABX*I_TWOBODYOVERLAP_H2x3y_D2x_a;
  Double I_TWOBODYOVERLAP_H2x2yz_F3x_a = I_TWOBODYOVERLAP_I3x2yz_D2x_a+ABX*I_TWOBODYOVERLAP_H2x2yz_D2x_a;
  Double I_TWOBODYOVERLAP_H2xy2z_F3x_a = I_TWOBODYOVERLAP_I3xy2z_D2x_a+ABX*I_TWOBODYOVERLAP_H2xy2z_D2x_a;
  Double I_TWOBODYOVERLAP_H2x3z_F3x_a = I_TWOBODYOVERLAP_I3x3z_D2x_a+ABX*I_TWOBODYOVERLAP_H2x3z_D2x_a;
  Double I_TWOBODYOVERLAP_Hx4y_F3x_a = I_TWOBODYOVERLAP_I2x4y_D2x_a+ABX*I_TWOBODYOVERLAP_Hx4y_D2x_a;
  Double I_TWOBODYOVERLAP_Hx3yz_F3x_a = I_TWOBODYOVERLAP_I2x3yz_D2x_a+ABX*I_TWOBODYOVERLAP_Hx3yz_D2x_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_F3x_a = I_TWOBODYOVERLAP_I2x2y2z_D2x_a+ABX*I_TWOBODYOVERLAP_Hx2y2z_D2x_a;
  Double I_TWOBODYOVERLAP_Hxy3z_F3x_a = I_TWOBODYOVERLAP_I2xy3z_D2x_a+ABX*I_TWOBODYOVERLAP_Hxy3z_D2x_a;
  Double I_TWOBODYOVERLAP_Hx4z_F3x_a = I_TWOBODYOVERLAP_I2x4z_D2x_a+ABX*I_TWOBODYOVERLAP_Hx4z_D2x_a;
  Double I_TWOBODYOVERLAP_H5y_F3x_a = I_TWOBODYOVERLAP_Ix5y_D2x_a+ABX*I_TWOBODYOVERLAP_H5y_D2x_a;
  Double I_TWOBODYOVERLAP_H4yz_F3x_a = I_TWOBODYOVERLAP_Ix4yz_D2x_a+ABX*I_TWOBODYOVERLAP_H4yz_D2x_a;
  Double I_TWOBODYOVERLAP_H3y2z_F3x_a = I_TWOBODYOVERLAP_Ix3y2z_D2x_a+ABX*I_TWOBODYOVERLAP_H3y2z_D2x_a;
  Double I_TWOBODYOVERLAP_H2y3z_F3x_a = I_TWOBODYOVERLAP_Ix2y3z_D2x_a+ABX*I_TWOBODYOVERLAP_H2y3z_D2x_a;
  Double I_TWOBODYOVERLAP_Hy4z_F3x_a = I_TWOBODYOVERLAP_Ixy4z_D2x_a+ABX*I_TWOBODYOVERLAP_Hy4z_D2x_a;
  Double I_TWOBODYOVERLAP_H5z_F3x_a = I_TWOBODYOVERLAP_Ix5z_D2x_a+ABX*I_TWOBODYOVERLAP_H5z_D2x_a;
  Double I_TWOBODYOVERLAP_H5x_F2xy_a = I_TWOBODYOVERLAP_I5xy_D2x_a+ABY*I_TWOBODYOVERLAP_H5x_D2x_a;
  Double I_TWOBODYOVERLAP_H4xy_F2xy_a = I_TWOBODYOVERLAP_I4x2y_D2x_a+ABY*I_TWOBODYOVERLAP_H4xy_D2x_a;
  Double I_TWOBODYOVERLAP_H4xz_F2xy_a = I_TWOBODYOVERLAP_I4xyz_D2x_a+ABY*I_TWOBODYOVERLAP_H4xz_D2x_a;
  Double I_TWOBODYOVERLAP_H3x2y_F2xy_a = I_TWOBODYOVERLAP_I3x3y_D2x_a+ABY*I_TWOBODYOVERLAP_H3x2y_D2x_a;
  Double I_TWOBODYOVERLAP_H3xyz_F2xy_a = I_TWOBODYOVERLAP_I3x2yz_D2x_a+ABY*I_TWOBODYOVERLAP_H3xyz_D2x_a;
  Double I_TWOBODYOVERLAP_H3x2z_F2xy_a = I_TWOBODYOVERLAP_I3xy2z_D2x_a+ABY*I_TWOBODYOVERLAP_H3x2z_D2x_a;
  Double I_TWOBODYOVERLAP_H2x3y_F2xy_a = I_TWOBODYOVERLAP_I2x4y_D2x_a+ABY*I_TWOBODYOVERLAP_H2x3y_D2x_a;
  Double I_TWOBODYOVERLAP_H2x2yz_F2xy_a = I_TWOBODYOVERLAP_I2x3yz_D2x_a+ABY*I_TWOBODYOVERLAP_H2x2yz_D2x_a;
  Double I_TWOBODYOVERLAP_H2xy2z_F2xy_a = I_TWOBODYOVERLAP_I2x2y2z_D2x_a+ABY*I_TWOBODYOVERLAP_H2xy2z_D2x_a;
  Double I_TWOBODYOVERLAP_H2x3z_F2xy_a = I_TWOBODYOVERLAP_I2xy3z_D2x_a+ABY*I_TWOBODYOVERLAP_H2x3z_D2x_a;
  Double I_TWOBODYOVERLAP_Hx4y_F2xy_a = I_TWOBODYOVERLAP_Ix5y_D2x_a+ABY*I_TWOBODYOVERLAP_Hx4y_D2x_a;
  Double I_TWOBODYOVERLAP_Hx3yz_F2xy_a = I_TWOBODYOVERLAP_Ix4yz_D2x_a+ABY*I_TWOBODYOVERLAP_Hx3yz_D2x_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_F2xy_a = I_TWOBODYOVERLAP_Ix3y2z_D2x_a+ABY*I_TWOBODYOVERLAP_Hx2y2z_D2x_a;
  Double I_TWOBODYOVERLAP_Hxy3z_F2xy_a = I_TWOBODYOVERLAP_Ix2y3z_D2x_a+ABY*I_TWOBODYOVERLAP_Hxy3z_D2x_a;
  Double I_TWOBODYOVERLAP_Hx4z_F2xy_a = I_TWOBODYOVERLAP_Ixy4z_D2x_a+ABY*I_TWOBODYOVERLAP_Hx4z_D2x_a;
  Double I_TWOBODYOVERLAP_H5y_F2xy_a = I_TWOBODYOVERLAP_I6y_D2x_a+ABY*I_TWOBODYOVERLAP_H5y_D2x_a;
  Double I_TWOBODYOVERLAP_H4yz_F2xy_a = I_TWOBODYOVERLAP_I5yz_D2x_a+ABY*I_TWOBODYOVERLAP_H4yz_D2x_a;
  Double I_TWOBODYOVERLAP_H3y2z_F2xy_a = I_TWOBODYOVERLAP_I4y2z_D2x_a+ABY*I_TWOBODYOVERLAP_H3y2z_D2x_a;
  Double I_TWOBODYOVERLAP_H2y3z_F2xy_a = I_TWOBODYOVERLAP_I3y3z_D2x_a+ABY*I_TWOBODYOVERLAP_H2y3z_D2x_a;
  Double I_TWOBODYOVERLAP_Hy4z_F2xy_a = I_TWOBODYOVERLAP_I2y4z_D2x_a+ABY*I_TWOBODYOVERLAP_Hy4z_D2x_a;
  Double I_TWOBODYOVERLAP_H5z_F2xy_a = I_TWOBODYOVERLAP_Iy5z_D2x_a+ABY*I_TWOBODYOVERLAP_H5z_D2x_a;
  Double I_TWOBODYOVERLAP_H5x_F2xz_a = I_TWOBODYOVERLAP_I5xz_D2x_a+ABZ*I_TWOBODYOVERLAP_H5x_D2x_a;
  Double I_TWOBODYOVERLAP_H4xy_F2xz_a = I_TWOBODYOVERLAP_I4xyz_D2x_a+ABZ*I_TWOBODYOVERLAP_H4xy_D2x_a;
  Double I_TWOBODYOVERLAP_H4xz_F2xz_a = I_TWOBODYOVERLAP_I4x2z_D2x_a+ABZ*I_TWOBODYOVERLAP_H4xz_D2x_a;
  Double I_TWOBODYOVERLAP_H3x2y_F2xz_a = I_TWOBODYOVERLAP_I3x2yz_D2x_a+ABZ*I_TWOBODYOVERLAP_H3x2y_D2x_a;
  Double I_TWOBODYOVERLAP_H3xyz_F2xz_a = I_TWOBODYOVERLAP_I3xy2z_D2x_a+ABZ*I_TWOBODYOVERLAP_H3xyz_D2x_a;
  Double I_TWOBODYOVERLAP_H3x2z_F2xz_a = I_TWOBODYOVERLAP_I3x3z_D2x_a+ABZ*I_TWOBODYOVERLAP_H3x2z_D2x_a;
  Double I_TWOBODYOVERLAP_H2x3y_F2xz_a = I_TWOBODYOVERLAP_I2x3yz_D2x_a+ABZ*I_TWOBODYOVERLAP_H2x3y_D2x_a;
  Double I_TWOBODYOVERLAP_H2x2yz_F2xz_a = I_TWOBODYOVERLAP_I2x2y2z_D2x_a+ABZ*I_TWOBODYOVERLAP_H2x2yz_D2x_a;
  Double I_TWOBODYOVERLAP_H2xy2z_F2xz_a = I_TWOBODYOVERLAP_I2xy3z_D2x_a+ABZ*I_TWOBODYOVERLAP_H2xy2z_D2x_a;
  Double I_TWOBODYOVERLAP_H2x3z_F2xz_a = I_TWOBODYOVERLAP_I2x4z_D2x_a+ABZ*I_TWOBODYOVERLAP_H2x3z_D2x_a;
  Double I_TWOBODYOVERLAP_Hx4y_F2xz_a = I_TWOBODYOVERLAP_Ix4yz_D2x_a+ABZ*I_TWOBODYOVERLAP_Hx4y_D2x_a;
  Double I_TWOBODYOVERLAP_Hx3yz_F2xz_a = I_TWOBODYOVERLAP_Ix3y2z_D2x_a+ABZ*I_TWOBODYOVERLAP_Hx3yz_D2x_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_F2xz_a = I_TWOBODYOVERLAP_Ix2y3z_D2x_a+ABZ*I_TWOBODYOVERLAP_Hx2y2z_D2x_a;
  Double I_TWOBODYOVERLAP_Hxy3z_F2xz_a = I_TWOBODYOVERLAP_Ixy4z_D2x_a+ABZ*I_TWOBODYOVERLAP_Hxy3z_D2x_a;
  Double I_TWOBODYOVERLAP_Hx4z_F2xz_a = I_TWOBODYOVERLAP_Ix5z_D2x_a+ABZ*I_TWOBODYOVERLAP_Hx4z_D2x_a;
  Double I_TWOBODYOVERLAP_H5y_F2xz_a = I_TWOBODYOVERLAP_I5yz_D2x_a+ABZ*I_TWOBODYOVERLAP_H5y_D2x_a;
  Double I_TWOBODYOVERLAP_H4yz_F2xz_a = I_TWOBODYOVERLAP_I4y2z_D2x_a+ABZ*I_TWOBODYOVERLAP_H4yz_D2x_a;
  Double I_TWOBODYOVERLAP_H3y2z_F2xz_a = I_TWOBODYOVERLAP_I3y3z_D2x_a+ABZ*I_TWOBODYOVERLAP_H3y2z_D2x_a;
  Double I_TWOBODYOVERLAP_H2y3z_F2xz_a = I_TWOBODYOVERLAP_I2y4z_D2x_a+ABZ*I_TWOBODYOVERLAP_H2y3z_D2x_a;
  Double I_TWOBODYOVERLAP_Hy4z_F2xz_a = I_TWOBODYOVERLAP_Iy5z_D2x_a+ABZ*I_TWOBODYOVERLAP_Hy4z_D2x_a;
  Double I_TWOBODYOVERLAP_H5z_F2xz_a = I_TWOBODYOVERLAP_I6z_D2x_a+ABZ*I_TWOBODYOVERLAP_H5z_D2x_a;
  Double I_TWOBODYOVERLAP_H5x_Fx2y_a = I_TWOBODYOVERLAP_I6x_D2y_a+ABX*I_TWOBODYOVERLAP_H5x_D2y_a;
  Double I_TWOBODYOVERLAP_H4xy_Fx2y_a = I_TWOBODYOVERLAP_I5xy_D2y_a+ABX*I_TWOBODYOVERLAP_H4xy_D2y_a;
  Double I_TWOBODYOVERLAP_H4xz_Fx2y_a = I_TWOBODYOVERLAP_I5xz_D2y_a+ABX*I_TWOBODYOVERLAP_H4xz_D2y_a;
  Double I_TWOBODYOVERLAP_H3x2y_Fx2y_a = I_TWOBODYOVERLAP_I4x2y_D2y_a+ABX*I_TWOBODYOVERLAP_H3x2y_D2y_a;
  Double I_TWOBODYOVERLAP_H3xyz_Fx2y_a = I_TWOBODYOVERLAP_I4xyz_D2y_a+ABX*I_TWOBODYOVERLAP_H3xyz_D2y_a;
  Double I_TWOBODYOVERLAP_H3x2z_Fx2y_a = I_TWOBODYOVERLAP_I4x2z_D2y_a+ABX*I_TWOBODYOVERLAP_H3x2z_D2y_a;
  Double I_TWOBODYOVERLAP_H2x3y_Fx2y_a = I_TWOBODYOVERLAP_I3x3y_D2y_a+ABX*I_TWOBODYOVERLAP_H2x3y_D2y_a;
  Double I_TWOBODYOVERLAP_H2x2yz_Fx2y_a = I_TWOBODYOVERLAP_I3x2yz_D2y_a+ABX*I_TWOBODYOVERLAP_H2x2yz_D2y_a;
  Double I_TWOBODYOVERLAP_H2xy2z_Fx2y_a = I_TWOBODYOVERLAP_I3xy2z_D2y_a+ABX*I_TWOBODYOVERLAP_H2xy2z_D2y_a;
  Double I_TWOBODYOVERLAP_H2x3z_Fx2y_a = I_TWOBODYOVERLAP_I3x3z_D2y_a+ABX*I_TWOBODYOVERLAP_H2x3z_D2y_a;
  Double I_TWOBODYOVERLAP_Hx4y_Fx2y_a = I_TWOBODYOVERLAP_I2x4y_D2y_a+ABX*I_TWOBODYOVERLAP_Hx4y_D2y_a;
  Double I_TWOBODYOVERLAP_Hx3yz_Fx2y_a = I_TWOBODYOVERLAP_I2x3yz_D2y_a+ABX*I_TWOBODYOVERLAP_Hx3yz_D2y_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_Fx2y_a = I_TWOBODYOVERLAP_I2x2y2z_D2y_a+ABX*I_TWOBODYOVERLAP_Hx2y2z_D2y_a;
  Double I_TWOBODYOVERLAP_Hxy3z_Fx2y_a = I_TWOBODYOVERLAP_I2xy3z_D2y_a+ABX*I_TWOBODYOVERLAP_Hxy3z_D2y_a;
  Double I_TWOBODYOVERLAP_Hx4z_Fx2y_a = I_TWOBODYOVERLAP_I2x4z_D2y_a+ABX*I_TWOBODYOVERLAP_Hx4z_D2y_a;
  Double I_TWOBODYOVERLAP_H5y_Fx2y_a = I_TWOBODYOVERLAP_Ix5y_D2y_a+ABX*I_TWOBODYOVERLAP_H5y_D2y_a;
  Double I_TWOBODYOVERLAP_H4yz_Fx2y_a = I_TWOBODYOVERLAP_Ix4yz_D2y_a+ABX*I_TWOBODYOVERLAP_H4yz_D2y_a;
  Double I_TWOBODYOVERLAP_H3y2z_Fx2y_a = I_TWOBODYOVERLAP_Ix3y2z_D2y_a+ABX*I_TWOBODYOVERLAP_H3y2z_D2y_a;
  Double I_TWOBODYOVERLAP_H2y3z_Fx2y_a = I_TWOBODYOVERLAP_Ix2y3z_D2y_a+ABX*I_TWOBODYOVERLAP_H2y3z_D2y_a;
  Double I_TWOBODYOVERLAP_Hy4z_Fx2y_a = I_TWOBODYOVERLAP_Ixy4z_D2y_a+ABX*I_TWOBODYOVERLAP_Hy4z_D2y_a;
  Double I_TWOBODYOVERLAP_H5z_Fx2y_a = I_TWOBODYOVERLAP_Ix5z_D2y_a+ABX*I_TWOBODYOVERLAP_H5z_D2y_a;
  Double I_TWOBODYOVERLAP_H5x_Fx2z_a = I_TWOBODYOVERLAP_I6x_D2z_a+ABX*I_TWOBODYOVERLAP_H5x_D2z_a;
  Double I_TWOBODYOVERLAP_H4xy_Fx2z_a = I_TWOBODYOVERLAP_I5xy_D2z_a+ABX*I_TWOBODYOVERLAP_H4xy_D2z_a;
  Double I_TWOBODYOVERLAP_H4xz_Fx2z_a = I_TWOBODYOVERLAP_I5xz_D2z_a+ABX*I_TWOBODYOVERLAP_H4xz_D2z_a;
  Double I_TWOBODYOVERLAP_H3x2y_Fx2z_a = I_TWOBODYOVERLAP_I4x2y_D2z_a+ABX*I_TWOBODYOVERLAP_H3x2y_D2z_a;
  Double I_TWOBODYOVERLAP_H3xyz_Fx2z_a = I_TWOBODYOVERLAP_I4xyz_D2z_a+ABX*I_TWOBODYOVERLAP_H3xyz_D2z_a;
  Double I_TWOBODYOVERLAP_H3x2z_Fx2z_a = I_TWOBODYOVERLAP_I4x2z_D2z_a+ABX*I_TWOBODYOVERLAP_H3x2z_D2z_a;
  Double I_TWOBODYOVERLAP_H2x3y_Fx2z_a = I_TWOBODYOVERLAP_I3x3y_D2z_a+ABX*I_TWOBODYOVERLAP_H2x3y_D2z_a;
  Double I_TWOBODYOVERLAP_H2x2yz_Fx2z_a = I_TWOBODYOVERLAP_I3x2yz_D2z_a+ABX*I_TWOBODYOVERLAP_H2x2yz_D2z_a;
  Double I_TWOBODYOVERLAP_H2xy2z_Fx2z_a = I_TWOBODYOVERLAP_I3xy2z_D2z_a+ABX*I_TWOBODYOVERLAP_H2xy2z_D2z_a;
  Double I_TWOBODYOVERLAP_H2x3z_Fx2z_a = I_TWOBODYOVERLAP_I3x3z_D2z_a+ABX*I_TWOBODYOVERLAP_H2x3z_D2z_a;
  Double I_TWOBODYOVERLAP_Hx4y_Fx2z_a = I_TWOBODYOVERLAP_I2x4y_D2z_a+ABX*I_TWOBODYOVERLAP_Hx4y_D2z_a;
  Double I_TWOBODYOVERLAP_Hx3yz_Fx2z_a = I_TWOBODYOVERLAP_I2x3yz_D2z_a+ABX*I_TWOBODYOVERLAP_Hx3yz_D2z_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_Fx2z_a = I_TWOBODYOVERLAP_I2x2y2z_D2z_a+ABX*I_TWOBODYOVERLAP_Hx2y2z_D2z_a;
  Double I_TWOBODYOVERLAP_Hxy3z_Fx2z_a = I_TWOBODYOVERLAP_I2xy3z_D2z_a+ABX*I_TWOBODYOVERLAP_Hxy3z_D2z_a;
  Double I_TWOBODYOVERLAP_Hx4z_Fx2z_a = I_TWOBODYOVERLAP_I2x4z_D2z_a+ABX*I_TWOBODYOVERLAP_Hx4z_D2z_a;
  Double I_TWOBODYOVERLAP_H5y_Fx2z_a = I_TWOBODYOVERLAP_Ix5y_D2z_a+ABX*I_TWOBODYOVERLAP_H5y_D2z_a;
  Double I_TWOBODYOVERLAP_H4yz_Fx2z_a = I_TWOBODYOVERLAP_Ix4yz_D2z_a+ABX*I_TWOBODYOVERLAP_H4yz_D2z_a;
  Double I_TWOBODYOVERLAP_H3y2z_Fx2z_a = I_TWOBODYOVERLAP_Ix3y2z_D2z_a+ABX*I_TWOBODYOVERLAP_H3y2z_D2z_a;
  Double I_TWOBODYOVERLAP_H2y3z_Fx2z_a = I_TWOBODYOVERLAP_Ix2y3z_D2z_a+ABX*I_TWOBODYOVERLAP_H2y3z_D2z_a;
  Double I_TWOBODYOVERLAP_Hy4z_Fx2z_a = I_TWOBODYOVERLAP_Ixy4z_D2z_a+ABX*I_TWOBODYOVERLAP_Hy4z_D2z_a;
  Double I_TWOBODYOVERLAP_H5z_Fx2z_a = I_TWOBODYOVERLAP_Ix5z_D2z_a+ABX*I_TWOBODYOVERLAP_H5z_D2z_a;
  Double I_TWOBODYOVERLAP_H5x_F3y_a = I_TWOBODYOVERLAP_I5xy_D2y_a+ABY*I_TWOBODYOVERLAP_H5x_D2y_a;
  Double I_TWOBODYOVERLAP_H4xy_F3y_a = I_TWOBODYOVERLAP_I4x2y_D2y_a+ABY*I_TWOBODYOVERLAP_H4xy_D2y_a;
  Double I_TWOBODYOVERLAP_H4xz_F3y_a = I_TWOBODYOVERLAP_I4xyz_D2y_a+ABY*I_TWOBODYOVERLAP_H4xz_D2y_a;
  Double I_TWOBODYOVERLAP_H3x2y_F3y_a = I_TWOBODYOVERLAP_I3x3y_D2y_a+ABY*I_TWOBODYOVERLAP_H3x2y_D2y_a;
  Double I_TWOBODYOVERLAP_H3xyz_F3y_a = I_TWOBODYOVERLAP_I3x2yz_D2y_a+ABY*I_TWOBODYOVERLAP_H3xyz_D2y_a;
  Double I_TWOBODYOVERLAP_H3x2z_F3y_a = I_TWOBODYOVERLAP_I3xy2z_D2y_a+ABY*I_TWOBODYOVERLAP_H3x2z_D2y_a;
  Double I_TWOBODYOVERLAP_H2x3y_F3y_a = I_TWOBODYOVERLAP_I2x4y_D2y_a+ABY*I_TWOBODYOVERLAP_H2x3y_D2y_a;
  Double I_TWOBODYOVERLAP_H2x2yz_F3y_a = I_TWOBODYOVERLAP_I2x3yz_D2y_a+ABY*I_TWOBODYOVERLAP_H2x2yz_D2y_a;
  Double I_TWOBODYOVERLAP_H2xy2z_F3y_a = I_TWOBODYOVERLAP_I2x2y2z_D2y_a+ABY*I_TWOBODYOVERLAP_H2xy2z_D2y_a;
  Double I_TWOBODYOVERLAP_H2x3z_F3y_a = I_TWOBODYOVERLAP_I2xy3z_D2y_a+ABY*I_TWOBODYOVERLAP_H2x3z_D2y_a;
  Double I_TWOBODYOVERLAP_Hx4y_F3y_a = I_TWOBODYOVERLAP_Ix5y_D2y_a+ABY*I_TWOBODYOVERLAP_Hx4y_D2y_a;
  Double I_TWOBODYOVERLAP_Hx3yz_F3y_a = I_TWOBODYOVERLAP_Ix4yz_D2y_a+ABY*I_TWOBODYOVERLAP_Hx3yz_D2y_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_F3y_a = I_TWOBODYOVERLAP_Ix3y2z_D2y_a+ABY*I_TWOBODYOVERLAP_Hx2y2z_D2y_a;
  Double I_TWOBODYOVERLAP_Hxy3z_F3y_a = I_TWOBODYOVERLAP_Ix2y3z_D2y_a+ABY*I_TWOBODYOVERLAP_Hxy3z_D2y_a;
  Double I_TWOBODYOVERLAP_Hx4z_F3y_a = I_TWOBODYOVERLAP_Ixy4z_D2y_a+ABY*I_TWOBODYOVERLAP_Hx4z_D2y_a;
  Double I_TWOBODYOVERLAP_H5y_F3y_a = I_TWOBODYOVERLAP_I6y_D2y_a+ABY*I_TWOBODYOVERLAP_H5y_D2y_a;
  Double I_TWOBODYOVERLAP_H4yz_F3y_a = I_TWOBODYOVERLAP_I5yz_D2y_a+ABY*I_TWOBODYOVERLAP_H4yz_D2y_a;
  Double I_TWOBODYOVERLAP_H3y2z_F3y_a = I_TWOBODYOVERLAP_I4y2z_D2y_a+ABY*I_TWOBODYOVERLAP_H3y2z_D2y_a;
  Double I_TWOBODYOVERLAP_H2y3z_F3y_a = I_TWOBODYOVERLAP_I3y3z_D2y_a+ABY*I_TWOBODYOVERLAP_H2y3z_D2y_a;
  Double I_TWOBODYOVERLAP_Hy4z_F3y_a = I_TWOBODYOVERLAP_I2y4z_D2y_a+ABY*I_TWOBODYOVERLAP_Hy4z_D2y_a;
  Double I_TWOBODYOVERLAP_H5z_F3y_a = I_TWOBODYOVERLAP_Iy5z_D2y_a+ABY*I_TWOBODYOVERLAP_H5z_D2y_a;
  Double I_TWOBODYOVERLAP_H5x_F2yz_a = I_TWOBODYOVERLAP_I5xz_D2y_a+ABZ*I_TWOBODYOVERLAP_H5x_D2y_a;
  Double I_TWOBODYOVERLAP_H4xy_F2yz_a = I_TWOBODYOVERLAP_I4xyz_D2y_a+ABZ*I_TWOBODYOVERLAP_H4xy_D2y_a;
  Double I_TWOBODYOVERLAP_H4xz_F2yz_a = I_TWOBODYOVERLAP_I4x2z_D2y_a+ABZ*I_TWOBODYOVERLAP_H4xz_D2y_a;
  Double I_TWOBODYOVERLAP_H3x2y_F2yz_a = I_TWOBODYOVERLAP_I3x2yz_D2y_a+ABZ*I_TWOBODYOVERLAP_H3x2y_D2y_a;
  Double I_TWOBODYOVERLAP_H3xyz_F2yz_a = I_TWOBODYOVERLAP_I3xy2z_D2y_a+ABZ*I_TWOBODYOVERLAP_H3xyz_D2y_a;
  Double I_TWOBODYOVERLAP_H3x2z_F2yz_a = I_TWOBODYOVERLAP_I3x3z_D2y_a+ABZ*I_TWOBODYOVERLAP_H3x2z_D2y_a;
  Double I_TWOBODYOVERLAP_H2x3y_F2yz_a = I_TWOBODYOVERLAP_I2x3yz_D2y_a+ABZ*I_TWOBODYOVERLAP_H2x3y_D2y_a;
  Double I_TWOBODYOVERLAP_H2x2yz_F2yz_a = I_TWOBODYOVERLAP_I2x2y2z_D2y_a+ABZ*I_TWOBODYOVERLAP_H2x2yz_D2y_a;
  Double I_TWOBODYOVERLAP_H2xy2z_F2yz_a = I_TWOBODYOVERLAP_I2xy3z_D2y_a+ABZ*I_TWOBODYOVERLAP_H2xy2z_D2y_a;
  Double I_TWOBODYOVERLAP_H2x3z_F2yz_a = I_TWOBODYOVERLAP_I2x4z_D2y_a+ABZ*I_TWOBODYOVERLAP_H2x3z_D2y_a;
  Double I_TWOBODYOVERLAP_Hx4y_F2yz_a = I_TWOBODYOVERLAP_Ix4yz_D2y_a+ABZ*I_TWOBODYOVERLAP_Hx4y_D2y_a;
  Double I_TWOBODYOVERLAP_Hx3yz_F2yz_a = I_TWOBODYOVERLAP_Ix3y2z_D2y_a+ABZ*I_TWOBODYOVERLAP_Hx3yz_D2y_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_F2yz_a = I_TWOBODYOVERLAP_Ix2y3z_D2y_a+ABZ*I_TWOBODYOVERLAP_Hx2y2z_D2y_a;
  Double I_TWOBODYOVERLAP_Hxy3z_F2yz_a = I_TWOBODYOVERLAP_Ixy4z_D2y_a+ABZ*I_TWOBODYOVERLAP_Hxy3z_D2y_a;
  Double I_TWOBODYOVERLAP_Hx4z_F2yz_a = I_TWOBODYOVERLAP_Ix5z_D2y_a+ABZ*I_TWOBODYOVERLAP_Hx4z_D2y_a;
  Double I_TWOBODYOVERLAP_H5y_F2yz_a = I_TWOBODYOVERLAP_I5yz_D2y_a+ABZ*I_TWOBODYOVERLAP_H5y_D2y_a;
  Double I_TWOBODYOVERLAP_H4yz_F2yz_a = I_TWOBODYOVERLAP_I4y2z_D2y_a+ABZ*I_TWOBODYOVERLAP_H4yz_D2y_a;
  Double I_TWOBODYOVERLAP_H3y2z_F2yz_a = I_TWOBODYOVERLAP_I3y3z_D2y_a+ABZ*I_TWOBODYOVERLAP_H3y2z_D2y_a;
  Double I_TWOBODYOVERLAP_H2y3z_F2yz_a = I_TWOBODYOVERLAP_I2y4z_D2y_a+ABZ*I_TWOBODYOVERLAP_H2y3z_D2y_a;
  Double I_TWOBODYOVERLAP_Hy4z_F2yz_a = I_TWOBODYOVERLAP_Iy5z_D2y_a+ABZ*I_TWOBODYOVERLAP_Hy4z_D2y_a;
  Double I_TWOBODYOVERLAP_H5z_F2yz_a = I_TWOBODYOVERLAP_I6z_D2y_a+ABZ*I_TWOBODYOVERLAP_H5z_D2y_a;
  Double I_TWOBODYOVERLAP_H5x_F3z_a = I_TWOBODYOVERLAP_I5xz_D2z_a+ABZ*I_TWOBODYOVERLAP_H5x_D2z_a;
  Double I_TWOBODYOVERLAP_H4xy_F3z_a = I_TWOBODYOVERLAP_I4xyz_D2z_a+ABZ*I_TWOBODYOVERLAP_H4xy_D2z_a;
  Double I_TWOBODYOVERLAP_H4xz_F3z_a = I_TWOBODYOVERLAP_I4x2z_D2z_a+ABZ*I_TWOBODYOVERLAP_H4xz_D2z_a;
  Double I_TWOBODYOVERLAP_H3x2y_F3z_a = I_TWOBODYOVERLAP_I3x2yz_D2z_a+ABZ*I_TWOBODYOVERLAP_H3x2y_D2z_a;
  Double I_TWOBODYOVERLAP_H3xyz_F3z_a = I_TWOBODYOVERLAP_I3xy2z_D2z_a+ABZ*I_TWOBODYOVERLAP_H3xyz_D2z_a;
  Double I_TWOBODYOVERLAP_H3x2z_F3z_a = I_TWOBODYOVERLAP_I3x3z_D2z_a+ABZ*I_TWOBODYOVERLAP_H3x2z_D2z_a;
  Double I_TWOBODYOVERLAP_H2x3y_F3z_a = I_TWOBODYOVERLAP_I2x3yz_D2z_a+ABZ*I_TWOBODYOVERLAP_H2x3y_D2z_a;
  Double I_TWOBODYOVERLAP_H2x2yz_F3z_a = I_TWOBODYOVERLAP_I2x2y2z_D2z_a+ABZ*I_TWOBODYOVERLAP_H2x2yz_D2z_a;
  Double I_TWOBODYOVERLAP_H2xy2z_F3z_a = I_TWOBODYOVERLAP_I2xy3z_D2z_a+ABZ*I_TWOBODYOVERLAP_H2xy2z_D2z_a;
  Double I_TWOBODYOVERLAP_H2x3z_F3z_a = I_TWOBODYOVERLAP_I2x4z_D2z_a+ABZ*I_TWOBODYOVERLAP_H2x3z_D2z_a;
  Double I_TWOBODYOVERLAP_Hx4y_F3z_a = I_TWOBODYOVERLAP_Ix4yz_D2z_a+ABZ*I_TWOBODYOVERLAP_Hx4y_D2z_a;
  Double I_TWOBODYOVERLAP_Hx3yz_F3z_a = I_TWOBODYOVERLAP_Ix3y2z_D2z_a+ABZ*I_TWOBODYOVERLAP_Hx3yz_D2z_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_F3z_a = I_TWOBODYOVERLAP_Ix2y3z_D2z_a+ABZ*I_TWOBODYOVERLAP_Hx2y2z_D2z_a;
  Double I_TWOBODYOVERLAP_Hxy3z_F3z_a = I_TWOBODYOVERLAP_Ixy4z_D2z_a+ABZ*I_TWOBODYOVERLAP_Hxy3z_D2z_a;
  Double I_TWOBODYOVERLAP_Hx4z_F3z_a = I_TWOBODYOVERLAP_Ix5z_D2z_a+ABZ*I_TWOBODYOVERLAP_Hx4z_D2z_a;
  Double I_TWOBODYOVERLAP_H5y_F3z_a = I_TWOBODYOVERLAP_I5yz_D2z_a+ABZ*I_TWOBODYOVERLAP_H5y_D2z_a;
  Double I_TWOBODYOVERLAP_H4yz_F3z_a = I_TWOBODYOVERLAP_I4y2z_D2z_a+ABZ*I_TWOBODYOVERLAP_H4yz_D2z_a;
  Double I_TWOBODYOVERLAP_H3y2z_F3z_a = I_TWOBODYOVERLAP_I3y3z_D2z_a+ABZ*I_TWOBODYOVERLAP_H3y2z_D2z_a;
  Double I_TWOBODYOVERLAP_H2y3z_F3z_a = I_TWOBODYOVERLAP_I2y4z_D2z_a+ABZ*I_TWOBODYOVERLAP_H2y3z_D2z_a;
  Double I_TWOBODYOVERLAP_Hy4z_F3z_a = I_TWOBODYOVERLAP_Iy5z_D2z_a+ABZ*I_TWOBODYOVERLAP_Hy4z_D2z_a;
  Double I_TWOBODYOVERLAP_H5z_F3z_a = I_TWOBODYOVERLAP_I6z_D2z_a+ABZ*I_TWOBODYOVERLAP_H5z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_L_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 30 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_M_S_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_L_S_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_L8x_Px_a = I_TWOBODYOVERLAP_M9x_S_a+ABX*I_TWOBODYOVERLAP_L8x_S_a;
  Double I_TWOBODYOVERLAP_L7xy_Px_a = I_TWOBODYOVERLAP_M8xy_S_a+ABX*I_TWOBODYOVERLAP_L7xy_S_a;
  Double I_TWOBODYOVERLAP_L7xz_Px_a = I_TWOBODYOVERLAP_M8xz_S_a+ABX*I_TWOBODYOVERLAP_L7xz_S_a;
  Double I_TWOBODYOVERLAP_L6x2y_Px_a = I_TWOBODYOVERLAP_M7x2y_S_a+ABX*I_TWOBODYOVERLAP_L6x2y_S_a;
  Double I_TWOBODYOVERLAP_L6xyz_Px_a = I_TWOBODYOVERLAP_M7xyz_S_a+ABX*I_TWOBODYOVERLAP_L6xyz_S_a;
  Double I_TWOBODYOVERLAP_L6x2z_Px_a = I_TWOBODYOVERLAP_M7x2z_S_a+ABX*I_TWOBODYOVERLAP_L6x2z_S_a;
  Double I_TWOBODYOVERLAP_L5x3y_Px_a = I_TWOBODYOVERLAP_M6x3y_S_a+ABX*I_TWOBODYOVERLAP_L5x3y_S_a;
  Double I_TWOBODYOVERLAP_L5x2yz_Px_a = I_TWOBODYOVERLAP_M6x2yz_S_a+ABX*I_TWOBODYOVERLAP_L5x2yz_S_a;
  Double I_TWOBODYOVERLAP_L5xy2z_Px_a = I_TWOBODYOVERLAP_M6xy2z_S_a+ABX*I_TWOBODYOVERLAP_L5xy2z_S_a;
  Double I_TWOBODYOVERLAP_L5x3z_Px_a = I_TWOBODYOVERLAP_M6x3z_S_a+ABX*I_TWOBODYOVERLAP_L5x3z_S_a;
  Double I_TWOBODYOVERLAP_L4x4y_Px_a = I_TWOBODYOVERLAP_M5x4y_S_a+ABX*I_TWOBODYOVERLAP_L4x4y_S_a;
  Double I_TWOBODYOVERLAP_L4x3yz_Px_a = I_TWOBODYOVERLAP_M5x3yz_S_a+ABX*I_TWOBODYOVERLAP_L4x3yz_S_a;
  Double I_TWOBODYOVERLAP_L4x2y2z_Px_a = I_TWOBODYOVERLAP_M5x2y2z_S_a+ABX*I_TWOBODYOVERLAP_L4x2y2z_S_a;
  Double I_TWOBODYOVERLAP_L4xy3z_Px_a = I_TWOBODYOVERLAP_M5xy3z_S_a+ABX*I_TWOBODYOVERLAP_L4xy3z_S_a;
  Double I_TWOBODYOVERLAP_L4x4z_Px_a = I_TWOBODYOVERLAP_M5x4z_S_a+ABX*I_TWOBODYOVERLAP_L4x4z_S_a;
  Double I_TWOBODYOVERLAP_L3x5y_Px_a = I_TWOBODYOVERLAP_M4x5y_S_a+ABX*I_TWOBODYOVERLAP_L3x5y_S_a;
  Double I_TWOBODYOVERLAP_L3x4yz_Px_a = I_TWOBODYOVERLAP_M4x4yz_S_a+ABX*I_TWOBODYOVERLAP_L3x4yz_S_a;
  Double I_TWOBODYOVERLAP_L3x3y2z_Px_a = I_TWOBODYOVERLAP_M4x3y2z_S_a+ABX*I_TWOBODYOVERLAP_L3x3y2z_S_a;
  Double I_TWOBODYOVERLAP_L3x2y3z_Px_a = I_TWOBODYOVERLAP_M4x2y3z_S_a+ABX*I_TWOBODYOVERLAP_L3x2y3z_S_a;
  Double I_TWOBODYOVERLAP_L3xy4z_Px_a = I_TWOBODYOVERLAP_M4xy4z_S_a+ABX*I_TWOBODYOVERLAP_L3xy4z_S_a;
  Double I_TWOBODYOVERLAP_L3x5z_Px_a = I_TWOBODYOVERLAP_M4x5z_S_a+ABX*I_TWOBODYOVERLAP_L3x5z_S_a;
  Double I_TWOBODYOVERLAP_L2x6y_Px_a = I_TWOBODYOVERLAP_M3x6y_S_a+ABX*I_TWOBODYOVERLAP_L2x6y_S_a;
  Double I_TWOBODYOVERLAP_L2x5yz_Px_a = I_TWOBODYOVERLAP_M3x5yz_S_a+ABX*I_TWOBODYOVERLAP_L2x5yz_S_a;
  Double I_TWOBODYOVERLAP_L2x4y2z_Px_a = I_TWOBODYOVERLAP_M3x4y2z_S_a+ABX*I_TWOBODYOVERLAP_L2x4y2z_S_a;
  Double I_TWOBODYOVERLAP_L2x3y3z_Px_a = I_TWOBODYOVERLAP_M3x3y3z_S_a+ABX*I_TWOBODYOVERLAP_L2x3y3z_S_a;
  Double I_TWOBODYOVERLAP_L2x2y4z_Px_a = I_TWOBODYOVERLAP_M3x2y4z_S_a+ABX*I_TWOBODYOVERLAP_L2x2y4z_S_a;
  Double I_TWOBODYOVERLAP_L2xy5z_Px_a = I_TWOBODYOVERLAP_M3xy5z_S_a+ABX*I_TWOBODYOVERLAP_L2xy5z_S_a;
  Double I_TWOBODYOVERLAP_L2x6z_Px_a = I_TWOBODYOVERLAP_M3x6z_S_a+ABX*I_TWOBODYOVERLAP_L2x6z_S_a;
  Double I_TWOBODYOVERLAP_Lx7y_Px_a = I_TWOBODYOVERLAP_M2x7y_S_a+ABX*I_TWOBODYOVERLAP_Lx7y_S_a;
  Double I_TWOBODYOVERLAP_Lx6yz_Px_a = I_TWOBODYOVERLAP_M2x6yz_S_a+ABX*I_TWOBODYOVERLAP_Lx6yz_S_a;
  Double I_TWOBODYOVERLAP_Lx5y2z_Px_a = I_TWOBODYOVERLAP_M2x5y2z_S_a+ABX*I_TWOBODYOVERLAP_Lx5y2z_S_a;
  Double I_TWOBODYOVERLAP_Lx4y3z_Px_a = I_TWOBODYOVERLAP_M2x4y3z_S_a+ABX*I_TWOBODYOVERLAP_Lx4y3z_S_a;
  Double I_TWOBODYOVERLAP_Lx3y4z_Px_a = I_TWOBODYOVERLAP_M2x3y4z_S_a+ABX*I_TWOBODYOVERLAP_Lx3y4z_S_a;
  Double I_TWOBODYOVERLAP_Lx2y5z_Px_a = I_TWOBODYOVERLAP_M2x2y5z_S_a+ABX*I_TWOBODYOVERLAP_Lx2y5z_S_a;
  Double I_TWOBODYOVERLAP_Lxy6z_Px_a = I_TWOBODYOVERLAP_M2xy6z_S_a+ABX*I_TWOBODYOVERLAP_Lxy6z_S_a;
  Double I_TWOBODYOVERLAP_Lx7z_Px_a = I_TWOBODYOVERLAP_M2x7z_S_a+ABX*I_TWOBODYOVERLAP_Lx7z_S_a;
  Double I_TWOBODYOVERLAP_L6x2y_Py_a = I_TWOBODYOVERLAP_M6x3y_S_a+ABY*I_TWOBODYOVERLAP_L6x2y_S_a;
  Double I_TWOBODYOVERLAP_L6xyz_Py_a = I_TWOBODYOVERLAP_M6x2yz_S_a+ABY*I_TWOBODYOVERLAP_L6xyz_S_a;
  Double I_TWOBODYOVERLAP_L5x3y_Py_a = I_TWOBODYOVERLAP_M5x4y_S_a+ABY*I_TWOBODYOVERLAP_L5x3y_S_a;
  Double I_TWOBODYOVERLAP_L5x2yz_Py_a = I_TWOBODYOVERLAP_M5x3yz_S_a+ABY*I_TWOBODYOVERLAP_L5x2yz_S_a;
  Double I_TWOBODYOVERLAP_L5xy2z_Py_a = I_TWOBODYOVERLAP_M5x2y2z_S_a+ABY*I_TWOBODYOVERLAP_L5xy2z_S_a;
  Double I_TWOBODYOVERLAP_L4x4y_Py_a = I_TWOBODYOVERLAP_M4x5y_S_a+ABY*I_TWOBODYOVERLAP_L4x4y_S_a;
  Double I_TWOBODYOVERLAP_L4x3yz_Py_a = I_TWOBODYOVERLAP_M4x4yz_S_a+ABY*I_TWOBODYOVERLAP_L4x3yz_S_a;
  Double I_TWOBODYOVERLAP_L4x2y2z_Py_a = I_TWOBODYOVERLAP_M4x3y2z_S_a+ABY*I_TWOBODYOVERLAP_L4x2y2z_S_a;
  Double I_TWOBODYOVERLAP_L4xy3z_Py_a = I_TWOBODYOVERLAP_M4x2y3z_S_a+ABY*I_TWOBODYOVERLAP_L4xy3z_S_a;
  Double I_TWOBODYOVERLAP_L3x5y_Py_a = I_TWOBODYOVERLAP_M3x6y_S_a+ABY*I_TWOBODYOVERLAP_L3x5y_S_a;
  Double I_TWOBODYOVERLAP_L3x4yz_Py_a = I_TWOBODYOVERLAP_M3x5yz_S_a+ABY*I_TWOBODYOVERLAP_L3x4yz_S_a;
  Double I_TWOBODYOVERLAP_L3x3y2z_Py_a = I_TWOBODYOVERLAP_M3x4y2z_S_a+ABY*I_TWOBODYOVERLAP_L3x3y2z_S_a;
  Double I_TWOBODYOVERLAP_L3x2y3z_Py_a = I_TWOBODYOVERLAP_M3x3y3z_S_a+ABY*I_TWOBODYOVERLAP_L3x2y3z_S_a;
  Double I_TWOBODYOVERLAP_L3xy4z_Py_a = I_TWOBODYOVERLAP_M3x2y4z_S_a+ABY*I_TWOBODYOVERLAP_L3xy4z_S_a;
  Double I_TWOBODYOVERLAP_L2x6y_Py_a = I_TWOBODYOVERLAP_M2x7y_S_a+ABY*I_TWOBODYOVERLAP_L2x6y_S_a;
  Double I_TWOBODYOVERLAP_L2x5yz_Py_a = I_TWOBODYOVERLAP_M2x6yz_S_a+ABY*I_TWOBODYOVERLAP_L2x5yz_S_a;
  Double I_TWOBODYOVERLAP_L2x4y2z_Py_a = I_TWOBODYOVERLAP_M2x5y2z_S_a+ABY*I_TWOBODYOVERLAP_L2x4y2z_S_a;
  Double I_TWOBODYOVERLAP_L2x3y3z_Py_a = I_TWOBODYOVERLAP_M2x4y3z_S_a+ABY*I_TWOBODYOVERLAP_L2x3y3z_S_a;
  Double I_TWOBODYOVERLAP_L2x2y4z_Py_a = I_TWOBODYOVERLAP_M2x3y4z_S_a+ABY*I_TWOBODYOVERLAP_L2x2y4z_S_a;
  Double I_TWOBODYOVERLAP_L2xy5z_Py_a = I_TWOBODYOVERLAP_M2x2y5z_S_a+ABY*I_TWOBODYOVERLAP_L2xy5z_S_a;
  Double I_TWOBODYOVERLAP_Lx7y_Py_a = I_TWOBODYOVERLAP_Mx8y_S_a+ABY*I_TWOBODYOVERLAP_Lx7y_S_a;
  Double I_TWOBODYOVERLAP_Lx6yz_Py_a = I_TWOBODYOVERLAP_Mx7yz_S_a+ABY*I_TWOBODYOVERLAP_Lx6yz_S_a;
  Double I_TWOBODYOVERLAP_Lx5y2z_Py_a = I_TWOBODYOVERLAP_Mx6y2z_S_a+ABY*I_TWOBODYOVERLAP_Lx5y2z_S_a;
  Double I_TWOBODYOVERLAP_Lx4y3z_Py_a = I_TWOBODYOVERLAP_Mx5y3z_S_a+ABY*I_TWOBODYOVERLAP_Lx4y3z_S_a;
  Double I_TWOBODYOVERLAP_Lx3y4z_Py_a = I_TWOBODYOVERLAP_Mx4y4z_S_a+ABY*I_TWOBODYOVERLAP_Lx3y4z_S_a;
  Double I_TWOBODYOVERLAP_Lx2y5z_Py_a = I_TWOBODYOVERLAP_Mx3y5z_S_a+ABY*I_TWOBODYOVERLAP_Lx2y5z_S_a;
  Double I_TWOBODYOVERLAP_Lxy6z_Py_a = I_TWOBODYOVERLAP_Mx2y6z_S_a+ABY*I_TWOBODYOVERLAP_Lxy6z_S_a;
  Double I_TWOBODYOVERLAP_L8y_Py_a = I_TWOBODYOVERLAP_M9y_S_a+ABY*I_TWOBODYOVERLAP_L8y_S_a;
  Double I_TWOBODYOVERLAP_L7yz_Py_a = I_TWOBODYOVERLAP_M8yz_S_a+ABY*I_TWOBODYOVERLAP_L7yz_S_a;
  Double I_TWOBODYOVERLAP_L6y2z_Py_a = I_TWOBODYOVERLAP_M7y2z_S_a+ABY*I_TWOBODYOVERLAP_L6y2z_S_a;
  Double I_TWOBODYOVERLAP_L5y3z_Py_a = I_TWOBODYOVERLAP_M6y3z_S_a+ABY*I_TWOBODYOVERLAP_L5y3z_S_a;
  Double I_TWOBODYOVERLAP_L4y4z_Py_a = I_TWOBODYOVERLAP_M5y4z_S_a+ABY*I_TWOBODYOVERLAP_L4y4z_S_a;
  Double I_TWOBODYOVERLAP_L3y5z_Py_a = I_TWOBODYOVERLAP_M4y5z_S_a+ABY*I_TWOBODYOVERLAP_L3y5z_S_a;
  Double I_TWOBODYOVERLAP_L2y6z_Py_a = I_TWOBODYOVERLAP_M3y6z_S_a+ABY*I_TWOBODYOVERLAP_L2y6z_S_a;
  Double I_TWOBODYOVERLAP_Ly7z_Py_a = I_TWOBODYOVERLAP_M2y7z_S_a+ABY*I_TWOBODYOVERLAP_Ly7z_S_a;
  Double I_TWOBODYOVERLAP_L6xyz_Pz_a = I_TWOBODYOVERLAP_M6xy2z_S_a+ABZ*I_TWOBODYOVERLAP_L6xyz_S_a;
  Double I_TWOBODYOVERLAP_L6x2z_Pz_a = I_TWOBODYOVERLAP_M6x3z_S_a+ABZ*I_TWOBODYOVERLAP_L6x2z_S_a;
  Double I_TWOBODYOVERLAP_L5x2yz_Pz_a = I_TWOBODYOVERLAP_M5x2y2z_S_a+ABZ*I_TWOBODYOVERLAP_L5x2yz_S_a;
  Double I_TWOBODYOVERLAP_L5xy2z_Pz_a = I_TWOBODYOVERLAP_M5xy3z_S_a+ABZ*I_TWOBODYOVERLAP_L5xy2z_S_a;
  Double I_TWOBODYOVERLAP_L5x3z_Pz_a = I_TWOBODYOVERLAP_M5x4z_S_a+ABZ*I_TWOBODYOVERLAP_L5x3z_S_a;
  Double I_TWOBODYOVERLAP_L4x3yz_Pz_a = I_TWOBODYOVERLAP_M4x3y2z_S_a+ABZ*I_TWOBODYOVERLAP_L4x3yz_S_a;
  Double I_TWOBODYOVERLAP_L4x2y2z_Pz_a = I_TWOBODYOVERLAP_M4x2y3z_S_a+ABZ*I_TWOBODYOVERLAP_L4x2y2z_S_a;
  Double I_TWOBODYOVERLAP_L4xy3z_Pz_a = I_TWOBODYOVERLAP_M4xy4z_S_a+ABZ*I_TWOBODYOVERLAP_L4xy3z_S_a;
  Double I_TWOBODYOVERLAP_L4x4z_Pz_a = I_TWOBODYOVERLAP_M4x5z_S_a+ABZ*I_TWOBODYOVERLAP_L4x4z_S_a;
  Double I_TWOBODYOVERLAP_L3x4yz_Pz_a = I_TWOBODYOVERLAP_M3x4y2z_S_a+ABZ*I_TWOBODYOVERLAP_L3x4yz_S_a;
  Double I_TWOBODYOVERLAP_L3x3y2z_Pz_a = I_TWOBODYOVERLAP_M3x3y3z_S_a+ABZ*I_TWOBODYOVERLAP_L3x3y2z_S_a;
  Double I_TWOBODYOVERLAP_L3x2y3z_Pz_a = I_TWOBODYOVERLAP_M3x2y4z_S_a+ABZ*I_TWOBODYOVERLAP_L3x2y3z_S_a;
  Double I_TWOBODYOVERLAP_L3xy4z_Pz_a = I_TWOBODYOVERLAP_M3xy5z_S_a+ABZ*I_TWOBODYOVERLAP_L3xy4z_S_a;
  Double I_TWOBODYOVERLAP_L3x5z_Pz_a = I_TWOBODYOVERLAP_M3x6z_S_a+ABZ*I_TWOBODYOVERLAP_L3x5z_S_a;
  Double I_TWOBODYOVERLAP_L2x5yz_Pz_a = I_TWOBODYOVERLAP_M2x5y2z_S_a+ABZ*I_TWOBODYOVERLAP_L2x5yz_S_a;
  Double I_TWOBODYOVERLAP_L2x4y2z_Pz_a = I_TWOBODYOVERLAP_M2x4y3z_S_a+ABZ*I_TWOBODYOVERLAP_L2x4y2z_S_a;
  Double I_TWOBODYOVERLAP_L2x3y3z_Pz_a = I_TWOBODYOVERLAP_M2x3y4z_S_a+ABZ*I_TWOBODYOVERLAP_L2x3y3z_S_a;
  Double I_TWOBODYOVERLAP_L2x2y4z_Pz_a = I_TWOBODYOVERLAP_M2x2y5z_S_a+ABZ*I_TWOBODYOVERLAP_L2x2y4z_S_a;
  Double I_TWOBODYOVERLAP_L2xy5z_Pz_a = I_TWOBODYOVERLAP_M2xy6z_S_a+ABZ*I_TWOBODYOVERLAP_L2xy5z_S_a;
  Double I_TWOBODYOVERLAP_L2x6z_Pz_a = I_TWOBODYOVERLAP_M2x7z_S_a+ABZ*I_TWOBODYOVERLAP_L2x6z_S_a;
  Double I_TWOBODYOVERLAP_Lx6yz_Pz_a = I_TWOBODYOVERLAP_Mx6y2z_S_a+ABZ*I_TWOBODYOVERLAP_Lx6yz_S_a;
  Double I_TWOBODYOVERLAP_Lx5y2z_Pz_a = I_TWOBODYOVERLAP_Mx5y3z_S_a+ABZ*I_TWOBODYOVERLAP_Lx5y2z_S_a;
  Double I_TWOBODYOVERLAP_Lx4y3z_Pz_a = I_TWOBODYOVERLAP_Mx4y4z_S_a+ABZ*I_TWOBODYOVERLAP_Lx4y3z_S_a;
  Double I_TWOBODYOVERLAP_Lx3y4z_Pz_a = I_TWOBODYOVERLAP_Mx3y5z_S_a+ABZ*I_TWOBODYOVERLAP_Lx3y4z_S_a;
  Double I_TWOBODYOVERLAP_Lx2y5z_Pz_a = I_TWOBODYOVERLAP_Mx2y6z_S_a+ABZ*I_TWOBODYOVERLAP_Lx2y5z_S_a;
  Double I_TWOBODYOVERLAP_Lxy6z_Pz_a = I_TWOBODYOVERLAP_Mxy7z_S_a+ABZ*I_TWOBODYOVERLAP_Lxy6z_S_a;
  Double I_TWOBODYOVERLAP_Lx7z_Pz_a = I_TWOBODYOVERLAP_Mx8z_S_a+ABZ*I_TWOBODYOVERLAP_Lx7z_S_a;
  Double I_TWOBODYOVERLAP_L6y2z_Pz_a = I_TWOBODYOVERLAP_M6y3z_S_a+ABZ*I_TWOBODYOVERLAP_L6y2z_S_a;
  Double I_TWOBODYOVERLAP_L5y3z_Pz_a = I_TWOBODYOVERLAP_M5y4z_S_a+ABZ*I_TWOBODYOVERLAP_L5y3z_S_a;
  Double I_TWOBODYOVERLAP_L4y4z_Pz_a = I_TWOBODYOVERLAP_M4y5z_S_a+ABZ*I_TWOBODYOVERLAP_L4y4z_S_a;
  Double I_TWOBODYOVERLAP_L3y5z_Pz_a = I_TWOBODYOVERLAP_M3y6z_S_a+ABZ*I_TWOBODYOVERLAP_L3y5z_S_a;
  Double I_TWOBODYOVERLAP_L2y6z_Pz_a = I_TWOBODYOVERLAP_M2y7z_S_a+ABZ*I_TWOBODYOVERLAP_L2y6z_S_a;
  Double I_TWOBODYOVERLAP_Ly7z_Pz_a = I_TWOBODYOVERLAP_My8z_S_a+ABZ*I_TWOBODYOVERLAP_Ly7z_S_a;
  Double I_TWOBODYOVERLAP_L8z_Pz_a = I_TWOBODYOVERLAP_M9z_S_a+ABZ*I_TWOBODYOVERLAP_L8z_S_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_K_D_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 111 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_L_P_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_P_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_K7x_D2x_a = I_TWOBODYOVERLAP_L8x_Px_a+ABX*I_TWOBODYOVERLAP_K7x_Px_a;
  Double I_TWOBODYOVERLAP_K6xy_D2x_a = I_TWOBODYOVERLAP_L7xy_Px_a+ABX*I_TWOBODYOVERLAP_K6xy_Px_a;
  Double I_TWOBODYOVERLAP_K6xz_D2x_a = I_TWOBODYOVERLAP_L7xz_Px_a+ABX*I_TWOBODYOVERLAP_K6xz_Px_a;
  Double I_TWOBODYOVERLAP_K5x2y_D2x_a = I_TWOBODYOVERLAP_L6x2y_Px_a+ABX*I_TWOBODYOVERLAP_K5x2y_Px_a;
  Double I_TWOBODYOVERLAP_K5xyz_D2x_a = I_TWOBODYOVERLAP_L6xyz_Px_a+ABX*I_TWOBODYOVERLAP_K5xyz_Px_a;
  Double I_TWOBODYOVERLAP_K5x2z_D2x_a = I_TWOBODYOVERLAP_L6x2z_Px_a+ABX*I_TWOBODYOVERLAP_K5x2z_Px_a;
  Double I_TWOBODYOVERLAP_K4x3y_D2x_a = I_TWOBODYOVERLAP_L5x3y_Px_a+ABX*I_TWOBODYOVERLAP_K4x3y_Px_a;
  Double I_TWOBODYOVERLAP_K4x2yz_D2x_a = I_TWOBODYOVERLAP_L5x2yz_Px_a+ABX*I_TWOBODYOVERLAP_K4x2yz_Px_a;
  Double I_TWOBODYOVERLAP_K4xy2z_D2x_a = I_TWOBODYOVERLAP_L5xy2z_Px_a+ABX*I_TWOBODYOVERLAP_K4xy2z_Px_a;
  Double I_TWOBODYOVERLAP_K4x3z_D2x_a = I_TWOBODYOVERLAP_L5x3z_Px_a+ABX*I_TWOBODYOVERLAP_K4x3z_Px_a;
  Double I_TWOBODYOVERLAP_K3x4y_D2x_a = I_TWOBODYOVERLAP_L4x4y_Px_a+ABX*I_TWOBODYOVERLAP_K3x4y_Px_a;
  Double I_TWOBODYOVERLAP_K3x3yz_D2x_a = I_TWOBODYOVERLAP_L4x3yz_Px_a+ABX*I_TWOBODYOVERLAP_K3x3yz_Px_a;
  Double I_TWOBODYOVERLAP_K3x2y2z_D2x_a = I_TWOBODYOVERLAP_L4x2y2z_Px_a+ABX*I_TWOBODYOVERLAP_K3x2y2z_Px_a;
  Double I_TWOBODYOVERLAP_K3xy3z_D2x_a = I_TWOBODYOVERLAP_L4xy3z_Px_a+ABX*I_TWOBODYOVERLAP_K3xy3z_Px_a;
  Double I_TWOBODYOVERLAP_K3x4z_D2x_a = I_TWOBODYOVERLAP_L4x4z_Px_a+ABX*I_TWOBODYOVERLAP_K3x4z_Px_a;
  Double I_TWOBODYOVERLAP_K2x5y_D2x_a = I_TWOBODYOVERLAP_L3x5y_Px_a+ABX*I_TWOBODYOVERLAP_K2x5y_Px_a;
  Double I_TWOBODYOVERLAP_K2x4yz_D2x_a = I_TWOBODYOVERLAP_L3x4yz_Px_a+ABX*I_TWOBODYOVERLAP_K2x4yz_Px_a;
  Double I_TWOBODYOVERLAP_K2x3y2z_D2x_a = I_TWOBODYOVERLAP_L3x3y2z_Px_a+ABX*I_TWOBODYOVERLAP_K2x3y2z_Px_a;
  Double I_TWOBODYOVERLAP_K2x2y3z_D2x_a = I_TWOBODYOVERLAP_L3x2y3z_Px_a+ABX*I_TWOBODYOVERLAP_K2x2y3z_Px_a;
  Double I_TWOBODYOVERLAP_K2xy4z_D2x_a = I_TWOBODYOVERLAP_L3xy4z_Px_a+ABX*I_TWOBODYOVERLAP_K2xy4z_Px_a;
  Double I_TWOBODYOVERLAP_K2x5z_D2x_a = I_TWOBODYOVERLAP_L3x5z_Px_a+ABX*I_TWOBODYOVERLAP_K2x5z_Px_a;
  Double I_TWOBODYOVERLAP_Kx6y_D2x_a = I_TWOBODYOVERLAP_L2x6y_Px_a+ABX*I_TWOBODYOVERLAP_Kx6y_Px_a;
  Double I_TWOBODYOVERLAP_Kx5yz_D2x_a = I_TWOBODYOVERLAP_L2x5yz_Px_a+ABX*I_TWOBODYOVERLAP_Kx5yz_Px_a;
  Double I_TWOBODYOVERLAP_Kx4y2z_D2x_a = I_TWOBODYOVERLAP_L2x4y2z_Px_a+ABX*I_TWOBODYOVERLAP_Kx4y2z_Px_a;
  Double I_TWOBODYOVERLAP_Kx3y3z_D2x_a = I_TWOBODYOVERLAP_L2x3y3z_Px_a+ABX*I_TWOBODYOVERLAP_Kx3y3z_Px_a;
  Double I_TWOBODYOVERLAP_Kx2y4z_D2x_a = I_TWOBODYOVERLAP_L2x2y4z_Px_a+ABX*I_TWOBODYOVERLAP_Kx2y4z_Px_a;
  Double I_TWOBODYOVERLAP_Kxy5z_D2x_a = I_TWOBODYOVERLAP_L2xy5z_Px_a+ABX*I_TWOBODYOVERLAP_Kxy5z_Px_a;
  Double I_TWOBODYOVERLAP_Kx6z_D2x_a = I_TWOBODYOVERLAP_L2x6z_Px_a+ABX*I_TWOBODYOVERLAP_Kx6z_Px_a;
  Double I_TWOBODYOVERLAP_K7y_D2x_a = I_TWOBODYOVERLAP_Lx7y_Px_a+ABX*I_TWOBODYOVERLAP_K7y_Px_a;
  Double I_TWOBODYOVERLAP_K6yz_D2x_a = I_TWOBODYOVERLAP_Lx6yz_Px_a+ABX*I_TWOBODYOVERLAP_K6yz_Px_a;
  Double I_TWOBODYOVERLAP_K5y2z_D2x_a = I_TWOBODYOVERLAP_Lx5y2z_Px_a+ABX*I_TWOBODYOVERLAP_K5y2z_Px_a;
  Double I_TWOBODYOVERLAP_K4y3z_D2x_a = I_TWOBODYOVERLAP_Lx4y3z_Px_a+ABX*I_TWOBODYOVERLAP_K4y3z_Px_a;
  Double I_TWOBODYOVERLAP_K3y4z_D2x_a = I_TWOBODYOVERLAP_Lx3y4z_Px_a+ABX*I_TWOBODYOVERLAP_K3y4z_Px_a;
  Double I_TWOBODYOVERLAP_K2y5z_D2x_a = I_TWOBODYOVERLAP_Lx2y5z_Px_a+ABX*I_TWOBODYOVERLAP_K2y5z_Px_a;
  Double I_TWOBODYOVERLAP_Ky6z_D2x_a = I_TWOBODYOVERLAP_Lxy6z_Px_a+ABX*I_TWOBODYOVERLAP_Ky6z_Px_a;
  Double I_TWOBODYOVERLAP_K7z_D2x_a = I_TWOBODYOVERLAP_Lx7z_Px_a+ABX*I_TWOBODYOVERLAP_K7z_Px_a;
  Double I_TWOBODYOVERLAP_K6xy_D2y_a = I_TWOBODYOVERLAP_L6x2y_Py_a+ABY*I_TWOBODYOVERLAP_K6xy_Py_a;
  Double I_TWOBODYOVERLAP_K6xz_D2y_a = I_TWOBODYOVERLAP_L6xyz_Py_a+ABY*I_TWOBODYOVERLAP_K6xz_Py_a;
  Double I_TWOBODYOVERLAP_K5x2y_D2y_a = I_TWOBODYOVERLAP_L5x3y_Py_a+ABY*I_TWOBODYOVERLAP_K5x2y_Py_a;
  Double I_TWOBODYOVERLAP_K5xyz_D2y_a = I_TWOBODYOVERLAP_L5x2yz_Py_a+ABY*I_TWOBODYOVERLAP_K5xyz_Py_a;
  Double I_TWOBODYOVERLAP_K5x2z_D2y_a = I_TWOBODYOVERLAP_L5xy2z_Py_a+ABY*I_TWOBODYOVERLAP_K5x2z_Py_a;
  Double I_TWOBODYOVERLAP_K4x3y_D2y_a = I_TWOBODYOVERLAP_L4x4y_Py_a+ABY*I_TWOBODYOVERLAP_K4x3y_Py_a;
  Double I_TWOBODYOVERLAP_K4x2yz_D2y_a = I_TWOBODYOVERLAP_L4x3yz_Py_a+ABY*I_TWOBODYOVERLAP_K4x2yz_Py_a;
  Double I_TWOBODYOVERLAP_K4xy2z_D2y_a = I_TWOBODYOVERLAP_L4x2y2z_Py_a+ABY*I_TWOBODYOVERLAP_K4xy2z_Py_a;
  Double I_TWOBODYOVERLAP_K4x3z_D2y_a = I_TWOBODYOVERLAP_L4xy3z_Py_a+ABY*I_TWOBODYOVERLAP_K4x3z_Py_a;
  Double I_TWOBODYOVERLAP_K3x4y_D2y_a = I_TWOBODYOVERLAP_L3x5y_Py_a+ABY*I_TWOBODYOVERLAP_K3x4y_Py_a;
  Double I_TWOBODYOVERLAP_K3x3yz_D2y_a = I_TWOBODYOVERLAP_L3x4yz_Py_a+ABY*I_TWOBODYOVERLAP_K3x3yz_Py_a;
  Double I_TWOBODYOVERLAP_K3x2y2z_D2y_a = I_TWOBODYOVERLAP_L3x3y2z_Py_a+ABY*I_TWOBODYOVERLAP_K3x2y2z_Py_a;
  Double I_TWOBODYOVERLAP_K3xy3z_D2y_a = I_TWOBODYOVERLAP_L3x2y3z_Py_a+ABY*I_TWOBODYOVERLAP_K3xy3z_Py_a;
  Double I_TWOBODYOVERLAP_K3x4z_D2y_a = I_TWOBODYOVERLAP_L3xy4z_Py_a+ABY*I_TWOBODYOVERLAP_K3x4z_Py_a;
  Double I_TWOBODYOVERLAP_K2x5y_D2y_a = I_TWOBODYOVERLAP_L2x6y_Py_a+ABY*I_TWOBODYOVERLAP_K2x5y_Py_a;
  Double I_TWOBODYOVERLAP_K2x4yz_D2y_a = I_TWOBODYOVERLAP_L2x5yz_Py_a+ABY*I_TWOBODYOVERLAP_K2x4yz_Py_a;
  Double I_TWOBODYOVERLAP_K2x3y2z_D2y_a = I_TWOBODYOVERLAP_L2x4y2z_Py_a+ABY*I_TWOBODYOVERLAP_K2x3y2z_Py_a;
  Double I_TWOBODYOVERLAP_K2x2y3z_D2y_a = I_TWOBODYOVERLAP_L2x3y3z_Py_a+ABY*I_TWOBODYOVERLAP_K2x2y3z_Py_a;
  Double I_TWOBODYOVERLAP_K2xy4z_D2y_a = I_TWOBODYOVERLAP_L2x2y4z_Py_a+ABY*I_TWOBODYOVERLAP_K2xy4z_Py_a;
  Double I_TWOBODYOVERLAP_K2x5z_D2y_a = I_TWOBODYOVERLAP_L2xy5z_Py_a+ABY*I_TWOBODYOVERLAP_K2x5z_Py_a;
  Double I_TWOBODYOVERLAP_Kx6y_D2y_a = I_TWOBODYOVERLAP_Lx7y_Py_a+ABY*I_TWOBODYOVERLAP_Kx6y_Py_a;
  Double I_TWOBODYOVERLAP_Kx5yz_D2y_a = I_TWOBODYOVERLAP_Lx6yz_Py_a+ABY*I_TWOBODYOVERLAP_Kx5yz_Py_a;
  Double I_TWOBODYOVERLAP_Kx4y2z_D2y_a = I_TWOBODYOVERLAP_Lx5y2z_Py_a+ABY*I_TWOBODYOVERLAP_Kx4y2z_Py_a;
  Double I_TWOBODYOVERLAP_Kx3y3z_D2y_a = I_TWOBODYOVERLAP_Lx4y3z_Py_a+ABY*I_TWOBODYOVERLAP_Kx3y3z_Py_a;
  Double I_TWOBODYOVERLAP_Kx2y4z_D2y_a = I_TWOBODYOVERLAP_Lx3y4z_Py_a+ABY*I_TWOBODYOVERLAP_Kx2y4z_Py_a;
  Double I_TWOBODYOVERLAP_Kxy5z_D2y_a = I_TWOBODYOVERLAP_Lx2y5z_Py_a+ABY*I_TWOBODYOVERLAP_Kxy5z_Py_a;
  Double I_TWOBODYOVERLAP_Kx6z_D2y_a = I_TWOBODYOVERLAP_Lxy6z_Py_a+ABY*I_TWOBODYOVERLAP_Kx6z_Py_a;
  Double I_TWOBODYOVERLAP_K7y_D2y_a = I_TWOBODYOVERLAP_L8y_Py_a+ABY*I_TWOBODYOVERLAP_K7y_Py_a;
  Double I_TWOBODYOVERLAP_K6yz_D2y_a = I_TWOBODYOVERLAP_L7yz_Py_a+ABY*I_TWOBODYOVERLAP_K6yz_Py_a;
  Double I_TWOBODYOVERLAP_K5y2z_D2y_a = I_TWOBODYOVERLAP_L6y2z_Py_a+ABY*I_TWOBODYOVERLAP_K5y2z_Py_a;
  Double I_TWOBODYOVERLAP_K4y3z_D2y_a = I_TWOBODYOVERLAP_L5y3z_Py_a+ABY*I_TWOBODYOVERLAP_K4y3z_Py_a;
  Double I_TWOBODYOVERLAP_K3y4z_D2y_a = I_TWOBODYOVERLAP_L4y4z_Py_a+ABY*I_TWOBODYOVERLAP_K3y4z_Py_a;
  Double I_TWOBODYOVERLAP_K2y5z_D2y_a = I_TWOBODYOVERLAP_L3y5z_Py_a+ABY*I_TWOBODYOVERLAP_K2y5z_Py_a;
  Double I_TWOBODYOVERLAP_Ky6z_D2y_a = I_TWOBODYOVERLAP_L2y6z_Py_a+ABY*I_TWOBODYOVERLAP_Ky6z_Py_a;
  Double I_TWOBODYOVERLAP_K7z_D2y_a = I_TWOBODYOVERLAP_Ly7z_Py_a+ABY*I_TWOBODYOVERLAP_K7z_Py_a;
  Double I_TWOBODYOVERLAP_K6xy_D2z_a = I_TWOBODYOVERLAP_L6xyz_Pz_a+ABZ*I_TWOBODYOVERLAP_K6xy_Pz_a;
  Double I_TWOBODYOVERLAP_K6xz_D2z_a = I_TWOBODYOVERLAP_L6x2z_Pz_a+ABZ*I_TWOBODYOVERLAP_K6xz_Pz_a;
  Double I_TWOBODYOVERLAP_K5x2y_D2z_a = I_TWOBODYOVERLAP_L5x2yz_Pz_a+ABZ*I_TWOBODYOVERLAP_K5x2y_Pz_a;
  Double I_TWOBODYOVERLAP_K5xyz_D2z_a = I_TWOBODYOVERLAP_L5xy2z_Pz_a+ABZ*I_TWOBODYOVERLAP_K5xyz_Pz_a;
  Double I_TWOBODYOVERLAP_K5x2z_D2z_a = I_TWOBODYOVERLAP_L5x3z_Pz_a+ABZ*I_TWOBODYOVERLAP_K5x2z_Pz_a;
  Double I_TWOBODYOVERLAP_K4x3y_D2z_a = I_TWOBODYOVERLAP_L4x3yz_Pz_a+ABZ*I_TWOBODYOVERLAP_K4x3y_Pz_a;
  Double I_TWOBODYOVERLAP_K4x2yz_D2z_a = I_TWOBODYOVERLAP_L4x2y2z_Pz_a+ABZ*I_TWOBODYOVERLAP_K4x2yz_Pz_a;
  Double I_TWOBODYOVERLAP_K4xy2z_D2z_a = I_TWOBODYOVERLAP_L4xy3z_Pz_a+ABZ*I_TWOBODYOVERLAP_K4xy2z_Pz_a;
  Double I_TWOBODYOVERLAP_K4x3z_D2z_a = I_TWOBODYOVERLAP_L4x4z_Pz_a+ABZ*I_TWOBODYOVERLAP_K4x3z_Pz_a;
  Double I_TWOBODYOVERLAP_K3x4y_D2z_a = I_TWOBODYOVERLAP_L3x4yz_Pz_a+ABZ*I_TWOBODYOVERLAP_K3x4y_Pz_a;
  Double I_TWOBODYOVERLAP_K3x3yz_D2z_a = I_TWOBODYOVERLAP_L3x3y2z_Pz_a+ABZ*I_TWOBODYOVERLAP_K3x3yz_Pz_a;
  Double I_TWOBODYOVERLAP_K3x2y2z_D2z_a = I_TWOBODYOVERLAP_L3x2y3z_Pz_a+ABZ*I_TWOBODYOVERLAP_K3x2y2z_Pz_a;
  Double I_TWOBODYOVERLAP_K3xy3z_D2z_a = I_TWOBODYOVERLAP_L3xy4z_Pz_a+ABZ*I_TWOBODYOVERLAP_K3xy3z_Pz_a;
  Double I_TWOBODYOVERLAP_K3x4z_D2z_a = I_TWOBODYOVERLAP_L3x5z_Pz_a+ABZ*I_TWOBODYOVERLAP_K3x4z_Pz_a;
  Double I_TWOBODYOVERLAP_K2x5y_D2z_a = I_TWOBODYOVERLAP_L2x5yz_Pz_a+ABZ*I_TWOBODYOVERLAP_K2x5y_Pz_a;
  Double I_TWOBODYOVERLAP_K2x4yz_D2z_a = I_TWOBODYOVERLAP_L2x4y2z_Pz_a+ABZ*I_TWOBODYOVERLAP_K2x4yz_Pz_a;
  Double I_TWOBODYOVERLAP_K2x3y2z_D2z_a = I_TWOBODYOVERLAP_L2x3y3z_Pz_a+ABZ*I_TWOBODYOVERLAP_K2x3y2z_Pz_a;
  Double I_TWOBODYOVERLAP_K2x2y3z_D2z_a = I_TWOBODYOVERLAP_L2x2y4z_Pz_a+ABZ*I_TWOBODYOVERLAP_K2x2y3z_Pz_a;
  Double I_TWOBODYOVERLAP_K2xy4z_D2z_a = I_TWOBODYOVERLAP_L2xy5z_Pz_a+ABZ*I_TWOBODYOVERLAP_K2xy4z_Pz_a;
  Double I_TWOBODYOVERLAP_K2x5z_D2z_a = I_TWOBODYOVERLAP_L2x6z_Pz_a+ABZ*I_TWOBODYOVERLAP_K2x5z_Pz_a;
  Double I_TWOBODYOVERLAP_Kx6y_D2z_a = I_TWOBODYOVERLAP_Lx6yz_Pz_a+ABZ*I_TWOBODYOVERLAP_Kx6y_Pz_a;
  Double I_TWOBODYOVERLAP_Kx5yz_D2z_a = I_TWOBODYOVERLAP_Lx5y2z_Pz_a+ABZ*I_TWOBODYOVERLAP_Kx5yz_Pz_a;
  Double I_TWOBODYOVERLAP_Kx4y2z_D2z_a = I_TWOBODYOVERLAP_Lx4y3z_Pz_a+ABZ*I_TWOBODYOVERLAP_Kx4y2z_Pz_a;
  Double I_TWOBODYOVERLAP_Kx3y3z_D2z_a = I_TWOBODYOVERLAP_Lx3y4z_Pz_a+ABZ*I_TWOBODYOVERLAP_Kx3y3z_Pz_a;
  Double I_TWOBODYOVERLAP_Kx2y4z_D2z_a = I_TWOBODYOVERLAP_Lx2y5z_Pz_a+ABZ*I_TWOBODYOVERLAP_Kx2y4z_Pz_a;
  Double I_TWOBODYOVERLAP_Kxy5z_D2z_a = I_TWOBODYOVERLAP_Lxy6z_Pz_a+ABZ*I_TWOBODYOVERLAP_Kxy5z_Pz_a;
  Double I_TWOBODYOVERLAP_Kx6z_D2z_a = I_TWOBODYOVERLAP_Lx7z_Pz_a+ABZ*I_TWOBODYOVERLAP_Kx6z_Pz_a;
  Double I_TWOBODYOVERLAP_K6yz_D2z_a = I_TWOBODYOVERLAP_L6y2z_Pz_a+ABZ*I_TWOBODYOVERLAP_K6yz_Pz_a;
  Double I_TWOBODYOVERLAP_K5y2z_D2z_a = I_TWOBODYOVERLAP_L5y3z_Pz_a+ABZ*I_TWOBODYOVERLAP_K5y2z_Pz_a;
  Double I_TWOBODYOVERLAP_K4y3z_D2z_a = I_TWOBODYOVERLAP_L4y4z_Pz_a+ABZ*I_TWOBODYOVERLAP_K4y3z_Pz_a;
  Double I_TWOBODYOVERLAP_K3y4z_D2z_a = I_TWOBODYOVERLAP_L3y5z_Pz_a+ABZ*I_TWOBODYOVERLAP_K3y4z_Pz_a;
  Double I_TWOBODYOVERLAP_K2y5z_D2z_a = I_TWOBODYOVERLAP_L2y6z_Pz_a+ABZ*I_TWOBODYOVERLAP_K2y5z_Pz_a;
  Double I_TWOBODYOVERLAP_Ky6z_D2z_a = I_TWOBODYOVERLAP_Ly7z_Pz_a+ABZ*I_TWOBODYOVERLAP_Ky6z_Pz_a;
  Double I_TWOBODYOVERLAP_K7z_D2z_a = I_TWOBODYOVERLAP_L8z_Pz_a+ABZ*I_TWOBODYOVERLAP_K7z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_I_F_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 85 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_D_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_D_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_I6x_F3x_a = I_TWOBODYOVERLAP_K7x_D2x_a+ABX*I_TWOBODYOVERLAP_I6x_D2x_a;
  Double I_TWOBODYOVERLAP_I5xy_F3x_a = I_TWOBODYOVERLAP_K6xy_D2x_a+ABX*I_TWOBODYOVERLAP_I5xy_D2x_a;
  Double I_TWOBODYOVERLAP_I5xz_F3x_a = I_TWOBODYOVERLAP_K6xz_D2x_a+ABX*I_TWOBODYOVERLAP_I5xz_D2x_a;
  Double I_TWOBODYOVERLAP_I4x2y_F3x_a = I_TWOBODYOVERLAP_K5x2y_D2x_a+ABX*I_TWOBODYOVERLAP_I4x2y_D2x_a;
  Double I_TWOBODYOVERLAP_I4xyz_F3x_a = I_TWOBODYOVERLAP_K5xyz_D2x_a+ABX*I_TWOBODYOVERLAP_I4xyz_D2x_a;
  Double I_TWOBODYOVERLAP_I4x2z_F3x_a = I_TWOBODYOVERLAP_K5x2z_D2x_a+ABX*I_TWOBODYOVERLAP_I4x2z_D2x_a;
  Double I_TWOBODYOVERLAP_I3x3y_F3x_a = I_TWOBODYOVERLAP_K4x3y_D2x_a+ABX*I_TWOBODYOVERLAP_I3x3y_D2x_a;
  Double I_TWOBODYOVERLAP_I3x2yz_F3x_a = I_TWOBODYOVERLAP_K4x2yz_D2x_a+ABX*I_TWOBODYOVERLAP_I3x2yz_D2x_a;
  Double I_TWOBODYOVERLAP_I3xy2z_F3x_a = I_TWOBODYOVERLAP_K4xy2z_D2x_a+ABX*I_TWOBODYOVERLAP_I3xy2z_D2x_a;
  Double I_TWOBODYOVERLAP_I3x3z_F3x_a = I_TWOBODYOVERLAP_K4x3z_D2x_a+ABX*I_TWOBODYOVERLAP_I3x3z_D2x_a;
  Double I_TWOBODYOVERLAP_I2x4y_F3x_a = I_TWOBODYOVERLAP_K3x4y_D2x_a+ABX*I_TWOBODYOVERLAP_I2x4y_D2x_a;
  Double I_TWOBODYOVERLAP_I2x3yz_F3x_a = I_TWOBODYOVERLAP_K3x3yz_D2x_a+ABX*I_TWOBODYOVERLAP_I2x3yz_D2x_a;
  Double I_TWOBODYOVERLAP_I2x2y2z_F3x_a = I_TWOBODYOVERLAP_K3x2y2z_D2x_a+ABX*I_TWOBODYOVERLAP_I2x2y2z_D2x_a;
  Double I_TWOBODYOVERLAP_I2xy3z_F3x_a = I_TWOBODYOVERLAP_K3xy3z_D2x_a+ABX*I_TWOBODYOVERLAP_I2xy3z_D2x_a;
  Double I_TWOBODYOVERLAP_I2x4z_F3x_a = I_TWOBODYOVERLAP_K3x4z_D2x_a+ABX*I_TWOBODYOVERLAP_I2x4z_D2x_a;
  Double I_TWOBODYOVERLAP_Ix5y_F3x_a = I_TWOBODYOVERLAP_K2x5y_D2x_a+ABX*I_TWOBODYOVERLAP_Ix5y_D2x_a;
  Double I_TWOBODYOVERLAP_Ix4yz_F3x_a = I_TWOBODYOVERLAP_K2x4yz_D2x_a+ABX*I_TWOBODYOVERLAP_Ix4yz_D2x_a;
  Double I_TWOBODYOVERLAP_Ix3y2z_F3x_a = I_TWOBODYOVERLAP_K2x3y2z_D2x_a+ABX*I_TWOBODYOVERLAP_Ix3y2z_D2x_a;
  Double I_TWOBODYOVERLAP_Ix2y3z_F3x_a = I_TWOBODYOVERLAP_K2x2y3z_D2x_a+ABX*I_TWOBODYOVERLAP_Ix2y3z_D2x_a;
  Double I_TWOBODYOVERLAP_Ixy4z_F3x_a = I_TWOBODYOVERLAP_K2xy4z_D2x_a+ABX*I_TWOBODYOVERLAP_Ixy4z_D2x_a;
  Double I_TWOBODYOVERLAP_Ix5z_F3x_a = I_TWOBODYOVERLAP_K2x5z_D2x_a+ABX*I_TWOBODYOVERLAP_Ix5z_D2x_a;
  Double I_TWOBODYOVERLAP_I6y_F3x_a = I_TWOBODYOVERLAP_Kx6y_D2x_a+ABX*I_TWOBODYOVERLAP_I6y_D2x_a;
  Double I_TWOBODYOVERLAP_I5yz_F3x_a = I_TWOBODYOVERLAP_Kx5yz_D2x_a+ABX*I_TWOBODYOVERLAP_I5yz_D2x_a;
  Double I_TWOBODYOVERLAP_I4y2z_F3x_a = I_TWOBODYOVERLAP_Kx4y2z_D2x_a+ABX*I_TWOBODYOVERLAP_I4y2z_D2x_a;
  Double I_TWOBODYOVERLAP_I3y3z_F3x_a = I_TWOBODYOVERLAP_Kx3y3z_D2x_a+ABX*I_TWOBODYOVERLAP_I3y3z_D2x_a;
  Double I_TWOBODYOVERLAP_I2y4z_F3x_a = I_TWOBODYOVERLAP_Kx2y4z_D2x_a+ABX*I_TWOBODYOVERLAP_I2y4z_D2x_a;
  Double I_TWOBODYOVERLAP_Iy5z_F3x_a = I_TWOBODYOVERLAP_Kxy5z_D2x_a+ABX*I_TWOBODYOVERLAP_Iy5z_D2x_a;
  Double I_TWOBODYOVERLAP_I6z_F3x_a = I_TWOBODYOVERLAP_Kx6z_D2x_a+ABX*I_TWOBODYOVERLAP_I6z_D2x_a;
  Double I_TWOBODYOVERLAP_I5xy_F2xy_a = I_TWOBODYOVERLAP_K5x2y_D2x_a+ABY*I_TWOBODYOVERLAP_I5xy_D2x_a;
  Double I_TWOBODYOVERLAP_I5xz_F2xy_a = I_TWOBODYOVERLAP_K5xyz_D2x_a+ABY*I_TWOBODYOVERLAP_I5xz_D2x_a;
  Double I_TWOBODYOVERLAP_I4x2y_F2xy_a = I_TWOBODYOVERLAP_K4x3y_D2x_a+ABY*I_TWOBODYOVERLAP_I4x2y_D2x_a;
  Double I_TWOBODYOVERLAP_I4xyz_F2xy_a = I_TWOBODYOVERLAP_K4x2yz_D2x_a+ABY*I_TWOBODYOVERLAP_I4xyz_D2x_a;
  Double I_TWOBODYOVERLAP_I4x2z_F2xy_a = I_TWOBODYOVERLAP_K4xy2z_D2x_a+ABY*I_TWOBODYOVERLAP_I4x2z_D2x_a;
  Double I_TWOBODYOVERLAP_I3x3y_F2xy_a = I_TWOBODYOVERLAP_K3x4y_D2x_a+ABY*I_TWOBODYOVERLAP_I3x3y_D2x_a;
  Double I_TWOBODYOVERLAP_I3x2yz_F2xy_a = I_TWOBODYOVERLAP_K3x3yz_D2x_a+ABY*I_TWOBODYOVERLAP_I3x2yz_D2x_a;
  Double I_TWOBODYOVERLAP_I3xy2z_F2xy_a = I_TWOBODYOVERLAP_K3x2y2z_D2x_a+ABY*I_TWOBODYOVERLAP_I3xy2z_D2x_a;
  Double I_TWOBODYOVERLAP_I3x3z_F2xy_a = I_TWOBODYOVERLAP_K3xy3z_D2x_a+ABY*I_TWOBODYOVERLAP_I3x3z_D2x_a;
  Double I_TWOBODYOVERLAP_I2x4y_F2xy_a = I_TWOBODYOVERLAP_K2x5y_D2x_a+ABY*I_TWOBODYOVERLAP_I2x4y_D2x_a;
  Double I_TWOBODYOVERLAP_I2x3yz_F2xy_a = I_TWOBODYOVERLAP_K2x4yz_D2x_a+ABY*I_TWOBODYOVERLAP_I2x3yz_D2x_a;
  Double I_TWOBODYOVERLAP_I2x2y2z_F2xy_a = I_TWOBODYOVERLAP_K2x3y2z_D2x_a+ABY*I_TWOBODYOVERLAP_I2x2y2z_D2x_a;
  Double I_TWOBODYOVERLAP_I2xy3z_F2xy_a = I_TWOBODYOVERLAP_K2x2y3z_D2x_a+ABY*I_TWOBODYOVERLAP_I2xy3z_D2x_a;
  Double I_TWOBODYOVERLAP_I2x4z_F2xy_a = I_TWOBODYOVERLAP_K2xy4z_D2x_a+ABY*I_TWOBODYOVERLAP_I2x4z_D2x_a;
  Double I_TWOBODYOVERLAP_Ix5y_F2xy_a = I_TWOBODYOVERLAP_Kx6y_D2x_a+ABY*I_TWOBODYOVERLAP_Ix5y_D2x_a;
  Double I_TWOBODYOVERLAP_Ix4yz_F2xy_a = I_TWOBODYOVERLAP_Kx5yz_D2x_a+ABY*I_TWOBODYOVERLAP_Ix4yz_D2x_a;
  Double I_TWOBODYOVERLAP_Ix3y2z_F2xy_a = I_TWOBODYOVERLAP_Kx4y2z_D2x_a+ABY*I_TWOBODYOVERLAP_Ix3y2z_D2x_a;
  Double I_TWOBODYOVERLAP_Ix2y3z_F2xy_a = I_TWOBODYOVERLAP_Kx3y3z_D2x_a+ABY*I_TWOBODYOVERLAP_Ix2y3z_D2x_a;
  Double I_TWOBODYOVERLAP_Ixy4z_F2xy_a = I_TWOBODYOVERLAP_Kx2y4z_D2x_a+ABY*I_TWOBODYOVERLAP_Ixy4z_D2x_a;
  Double I_TWOBODYOVERLAP_Ix5z_F2xy_a = I_TWOBODYOVERLAP_Kxy5z_D2x_a+ABY*I_TWOBODYOVERLAP_Ix5z_D2x_a;
  Double I_TWOBODYOVERLAP_I6y_F2xy_a = I_TWOBODYOVERLAP_K7y_D2x_a+ABY*I_TWOBODYOVERLAP_I6y_D2x_a;
  Double I_TWOBODYOVERLAP_I5yz_F2xy_a = I_TWOBODYOVERLAP_K6yz_D2x_a+ABY*I_TWOBODYOVERLAP_I5yz_D2x_a;
  Double I_TWOBODYOVERLAP_I4y2z_F2xy_a = I_TWOBODYOVERLAP_K5y2z_D2x_a+ABY*I_TWOBODYOVERLAP_I4y2z_D2x_a;
  Double I_TWOBODYOVERLAP_I3y3z_F2xy_a = I_TWOBODYOVERLAP_K4y3z_D2x_a+ABY*I_TWOBODYOVERLAP_I3y3z_D2x_a;
  Double I_TWOBODYOVERLAP_I2y4z_F2xy_a = I_TWOBODYOVERLAP_K3y4z_D2x_a+ABY*I_TWOBODYOVERLAP_I2y4z_D2x_a;
  Double I_TWOBODYOVERLAP_Iy5z_F2xy_a = I_TWOBODYOVERLAP_K2y5z_D2x_a+ABY*I_TWOBODYOVERLAP_Iy5z_D2x_a;
  Double I_TWOBODYOVERLAP_I6z_F2xy_a = I_TWOBODYOVERLAP_Ky6z_D2x_a+ABY*I_TWOBODYOVERLAP_I6z_D2x_a;
  Double I_TWOBODYOVERLAP_I5xz_F2xz_a = I_TWOBODYOVERLAP_K5x2z_D2x_a+ABZ*I_TWOBODYOVERLAP_I5xz_D2x_a;
  Double I_TWOBODYOVERLAP_I4xyz_F2xz_a = I_TWOBODYOVERLAP_K4xy2z_D2x_a+ABZ*I_TWOBODYOVERLAP_I4xyz_D2x_a;
  Double I_TWOBODYOVERLAP_I4x2z_F2xz_a = I_TWOBODYOVERLAP_K4x3z_D2x_a+ABZ*I_TWOBODYOVERLAP_I4x2z_D2x_a;
  Double I_TWOBODYOVERLAP_I3x2yz_F2xz_a = I_TWOBODYOVERLAP_K3x2y2z_D2x_a+ABZ*I_TWOBODYOVERLAP_I3x2yz_D2x_a;
  Double I_TWOBODYOVERLAP_I3xy2z_F2xz_a = I_TWOBODYOVERLAP_K3xy3z_D2x_a+ABZ*I_TWOBODYOVERLAP_I3xy2z_D2x_a;
  Double I_TWOBODYOVERLAP_I3x3z_F2xz_a = I_TWOBODYOVERLAP_K3x4z_D2x_a+ABZ*I_TWOBODYOVERLAP_I3x3z_D2x_a;
  Double I_TWOBODYOVERLAP_I2x3yz_F2xz_a = I_TWOBODYOVERLAP_K2x3y2z_D2x_a+ABZ*I_TWOBODYOVERLAP_I2x3yz_D2x_a;
  Double I_TWOBODYOVERLAP_I2x2y2z_F2xz_a = I_TWOBODYOVERLAP_K2x2y3z_D2x_a+ABZ*I_TWOBODYOVERLAP_I2x2y2z_D2x_a;
  Double I_TWOBODYOVERLAP_I2xy3z_F2xz_a = I_TWOBODYOVERLAP_K2xy4z_D2x_a+ABZ*I_TWOBODYOVERLAP_I2xy3z_D2x_a;
  Double I_TWOBODYOVERLAP_I2x4z_F2xz_a = I_TWOBODYOVERLAP_K2x5z_D2x_a+ABZ*I_TWOBODYOVERLAP_I2x4z_D2x_a;
  Double I_TWOBODYOVERLAP_Ix4yz_F2xz_a = I_TWOBODYOVERLAP_Kx4y2z_D2x_a+ABZ*I_TWOBODYOVERLAP_Ix4yz_D2x_a;
  Double I_TWOBODYOVERLAP_Ix3y2z_F2xz_a = I_TWOBODYOVERLAP_Kx3y3z_D2x_a+ABZ*I_TWOBODYOVERLAP_Ix3y2z_D2x_a;
  Double I_TWOBODYOVERLAP_Ix2y3z_F2xz_a = I_TWOBODYOVERLAP_Kx2y4z_D2x_a+ABZ*I_TWOBODYOVERLAP_Ix2y3z_D2x_a;
  Double I_TWOBODYOVERLAP_Ixy4z_F2xz_a = I_TWOBODYOVERLAP_Kxy5z_D2x_a+ABZ*I_TWOBODYOVERLAP_Ixy4z_D2x_a;
  Double I_TWOBODYOVERLAP_Ix5z_F2xz_a = I_TWOBODYOVERLAP_Kx6z_D2x_a+ABZ*I_TWOBODYOVERLAP_Ix5z_D2x_a;
  Double I_TWOBODYOVERLAP_I5yz_F2xz_a = I_TWOBODYOVERLAP_K5y2z_D2x_a+ABZ*I_TWOBODYOVERLAP_I5yz_D2x_a;
  Double I_TWOBODYOVERLAP_I4y2z_F2xz_a = I_TWOBODYOVERLAP_K4y3z_D2x_a+ABZ*I_TWOBODYOVERLAP_I4y2z_D2x_a;
  Double I_TWOBODYOVERLAP_I3y3z_F2xz_a = I_TWOBODYOVERLAP_K3y4z_D2x_a+ABZ*I_TWOBODYOVERLAP_I3y3z_D2x_a;
  Double I_TWOBODYOVERLAP_I2y4z_F2xz_a = I_TWOBODYOVERLAP_K2y5z_D2x_a+ABZ*I_TWOBODYOVERLAP_I2y4z_D2x_a;
  Double I_TWOBODYOVERLAP_Iy5z_F2xz_a = I_TWOBODYOVERLAP_Ky6z_D2x_a+ABZ*I_TWOBODYOVERLAP_Iy5z_D2x_a;
  Double I_TWOBODYOVERLAP_I6z_F2xz_a = I_TWOBODYOVERLAP_K7z_D2x_a+ABZ*I_TWOBODYOVERLAP_I6z_D2x_a;
  Double I_TWOBODYOVERLAP_I5xz_Fx2y_a = I_TWOBODYOVERLAP_K6xz_D2y_a+ABX*I_TWOBODYOVERLAP_I5xz_D2y_a;
  Double I_TWOBODYOVERLAP_I4xyz_Fx2y_a = I_TWOBODYOVERLAP_K5xyz_D2y_a+ABX*I_TWOBODYOVERLAP_I4xyz_D2y_a;
  Double I_TWOBODYOVERLAP_I4x2z_Fx2y_a = I_TWOBODYOVERLAP_K5x2z_D2y_a+ABX*I_TWOBODYOVERLAP_I4x2z_D2y_a;
  Double I_TWOBODYOVERLAP_I3x2yz_Fx2y_a = I_TWOBODYOVERLAP_K4x2yz_D2y_a+ABX*I_TWOBODYOVERLAP_I3x2yz_D2y_a;
  Double I_TWOBODYOVERLAP_I3xy2z_Fx2y_a = I_TWOBODYOVERLAP_K4xy2z_D2y_a+ABX*I_TWOBODYOVERLAP_I3xy2z_D2y_a;
  Double I_TWOBODYOVERLAP_I3x3z_Fx2y_a = I_TWOBODYOVERLAP_K4x3z_D2y_a+ABX*I_TWOBODYOVERLAP_I3x3z_D2y_a;
  Double I_TWOBODYOVERLAP_I2x3yz_Fx2y_a = I_TWOBODYOVERLAP_K3x3yz_D2y_a+ABX*I_TWOBODYOVERLAP_I2x3yz_D2y_a;
  Double I_TWOBODYOVERLAP_I2x2y2z_Fx2y_a = I_TWOBODYOVERLAP_K3x2y2z_D2y_a+ABX*I_TWOBODYOVERLAP_I2x2y2z_D2y_a;
  Double I_TWOBODYOVERLAP_I2xy3z_Fx2y_a = I_TWOBODYOVERLAP_K3xy3z_D2y_a+ABX*I_TWOBODYOVERLAP_I2xy3z_D2y_a;
  Double I_TWOBODYOVERLAP_I2x4z_Fx2y_a = I_TWOBODYOVERLAP_K3x4z_D2y_a+ABX*I_TWOBODYOVERLAP_I2x4z_D2y_a;
  Double I_TWOBODYOVERLAP_Ix4yz_Fx2y_a = I_TWOBODYOVERLAP_K2x4yz_D2y_a+ABX*I_TWOBODYOVERLAP_Ix4yz_D2y_a;
  Double I_TWOBODYOVERLAP_Ix3y2z_Fx2y_a = I_TWOBODYOVERLAP_K2x3y2z_D2y_a+ABX*I_TWOBODYOVERLAP_Ix3y2z_D2y_a;
  Double I_TWOBODYOVERLAP_Ix2y3z_Fx2y_a = I_TWOBODYOVERLAP_K2x2y3z_D2y_a+ABX*I_TWOBODYOVERLAP_Ix2y3z_D2y_a;
  Double I_TWOBODYOVERLAP_Ixy4z_Fx2y_a = I_TWOBODYOVERLAP_K2xy4z_D2y_a+ABX*I_TWOBODYOVERLAP_Ixy4z_D2y_a;
  Double I_TWOBODYOVERLAP_Ix5z_Fx2y_a = I_TWOBODYOVERLAP_K2x5z_D2y_a+ABX*I_TWOBODYOVERLAP_Ix5z_D2y_a;
  Double I_TWOBODYOVERLAP_I5yz_Fx2y_a = I_TWOBODYOVERLAP_Kx5yz_D2y_a+ABX*I_TWOBODYOVERLAP_I5yz_D2y_a;
  Double I_TWOBODYOVERLAP_I4y2z_Fx2y_a = I_TWOBODYOVERLAP_Kx4y2z_D2y_a+ABX*I_TWOBODYOVERLAP_I4y2z_D2y_a;
  Double I_TWOBODYOVERLAP_I3y3z_Fx2y_a = I_TWOBODYOVERLAP_Kx3y3z_D2y_a+ABX*I_TWOBODYOVERLAP_I3y3z_D2y_a;
  Double I_TWOBODYOVERLAP_I2y4z_Fx2y_a = I_TWOBODYOVERLAP_Kx2y4z_D2y_a+ABX*I_TWOBODYOVERLAP_I2y4z_D2y_a;
  Double I_TWOBODYOVERLAP_Iy5z_Fx2y_a = I_TWOBODYOVERLAP_Kxy5z_D2y_a+ABX*I_TWOBODYOVERLAP_Iy5z_D2y_a;
  Double I_TWOBODYOVERLAP_I6z_Fx2y_a = I_TWOBODYOVERLAP_Kx6z_D2y_a+ABX*I_TWOBODYOVERLAP_I6z_D2y_a;
  Double I_TWOBODYOVERLAP_I5xy_Fx2z_a = I_TWOBODYOVERLAP_K6xy_D2z_a+ABX*I_TWOBODYOVERLAP_I5xy_D2z_a;
  Double I_TWOBODYOVERLAP_I4x2y_Fx2z_a = I_TWOBODYOVERLAP_K5x2y_D2z_a+ABX*I_TWOBODYOVERLAP_I4x2y_D2z_a;
  Double I_TWOBODYOVERLAP_I4xyz_Fx2z_a = I_TWOBODYOVERLAP_K5xyz_D2z_a+ABX*I_TWOBODYOVERLAP_I4xyz_D2z_a;
  Double I_TWOBODYOVERLAP_I3x3y_Fx2z_a = I_TWOBODYOVERLAP_K4x3y_D2z_a+ABX*I_TWOBODYOVERLAP_I3x3y_D2z_a;
  Double I_TWOBODYOVERLAP_I3x2yz_Fx2z_a = I_TWOBODYOVERLAP_K4x2yz_D2z_a+ABX*I_TWOBODYOVERLAP_I3x2yz_D2z_a;
  Double I_TWOBODYOVERLAP_I3xy2z_Fx2z_a = I_TWOBODYOVERLAP_K4xy2z_D2z_a+ABX*I_TWOBODYOVERLAP_I3xy2z_D2z_a;
  Double I_TWOBODYOVERLAP_I2x4y_Fx2z_a = I_TWOBODYOVERLAP_K3x4y_D2z_a+ABX*I_TWOBODYOVERLAP_I2x4y_D2z_a;
  Double I_TWOBODYOVERLAP_I2x3yz_Fx2z_a = I_TWOBODYOVERLAP_K3x3yz_D2z_a+ABX*I_TWOBODYOVERLAP_I2x3yz_D2z_a;
  Double I_TWOBODYOVERLAP_I2x2y2z_Fx2z_a = I_TWOBODYOVERLAP_K3x2y2z_D2z_a+ABX*I_TWOBODYOVERLAP_I2x2y2z_D2z_a;
  Double I_TWOBODYOVERLAP_I2xy3z_Fx2z_a = I_TWOBODYOVERLAP_K3xy3z_D2z_a+ABX*I_TWOBODYOVERLAP_I2xy3z_D2z_a;
  Double I_TWOBODYOVERLAP_Ix5y_Fx2z_a = I_TWOBODYOVERLAP_K2x5y_D2z_a+ABX*I_TWOBODYOVERLAP_Ix5y_D2z_a;
  Double I_TWOBODYOVERLAP_Ix4yz_Fx2z_a = I_TWOBODYOVERLAP_K2x4yz_D2z_a+ABX*I_TWOBODYOVERLAP_Ix4yz_D2z_a;
  Double I_TWOBODYOVERLAP_Ix3y2z_Fx2z_a = I_TWOBODYOVERLAP_K2x3y2z_D2z_a+ABX*I_TWOBODYOVERLAP_Ix3y2z_D2z_a;
  Double I_TWOBODYOVERLAP_Ix2y3z_Fx2z_a = I_TWOBODYOVERLAP_K2x2y3z_D2z_a+ABX*I_TWOBODYOVERLAP_Ix2y3z_D2z_a;
  Double I_TWOBODYOVERLAP_Ixy4z_Fx2z_a = I_TWOBODYOVERLAP_K2xy4z_D2z_a+ABX*I_TWOBODYOVERLAP_Ixy4z_D2z_a;
  Double I_TWOBODYOVERLAP_I6y_Fx2z_a = I_TWOBODYOVERLAP_Kx6y_D2z_a+ABX*I_TWOBODYOVERLAP_I6y_D2z_a;
  Double I_TWOBODYOVERLAP_I5yz_Fx2z_a = I_TWOBODYOVERLAP_Kx5yz_D2z_a+ABX*I_TWOBODYOVERLAP_I5yz_D2z_a;
  Double I_TWOBODYOVERLAP_I4y2z_Fx2z_a = I_TWOBODYOVERLAP_Kx4y2z_D2z_a+ABX*I_TWOBODYOVERLAP_I4y2z_D2z_a;
  Double I_TWOBODYOVERLAP_I3y3z_Fx2z_a = I_TWOBODYOVERLAP_Kx3y3z_D2z_a+ABX*I_TWOBODYOVERLAP_I3y3z_D2z_a;
  Double I_TWOBODYOVERLAP_I2y4z_Fx2z_a = I_TWOBODYOVERLAP_Kx2y4z_D2z_a+ABX*I_TWOBODYOVERLAP_I2y4z_D2z_a;
  Double I_TWOBODYOVERLAP_Iy5z_Fx2z_a = I_TWOBODYOVERLAP_Kxy5z_D2z_a+ABX*I_TWOBODYOVERLAP_Iy5z_D2z_a;
  Double I_TWOBODYOVERLAP_I6x_F3y_a = I_TWOBODYOVERLAP_K6xy_D2y_a+ABY*I_TWOBODYOVERLAP_I6x_D2y_a;
  Double I_TWOBODYOVERLAP_I5xy_F3y_a = I_TWOBODYOVERLAP_K5x2y_D2y_a+ABY*I_TWOBODYOVERLAP_I5xy_D2y_a;
  Double I_TWOBODYOVERLAP_I5xz_F3y_a = I_TWOBODYOVERLAP_K5xyz_D2y_a+ABY*I_TWOBODYOVERLAP_I5xz_D2y_a;
  Double I_TWOBODYOVERLAP_I4x2y_F3y_a = I_TWOBODYOVERLAP_K4x3y_D2y_a+ABY*I_TWOBODYOVERLAP_I4x2y_D2y_a;
  Double I_TWOBODYOVERLAP_I4xyz_F3y_a = I_TWOBODYOVERLAP_K4x2yz_D2y_a+ABY*I_TWOBODYOVERLAP_I4xyz_D2y_a;
  Double I_TWOBODYOVERLAP_I4x2z_F3y_a = I_TWOBODYOVERLAP_K4xy2z_D2y_a+ABY*I_TWOBODYOVERLAP_I4x2z_D2y_a;
  Double I_TWOBODYOVERLAP_I3x3y_F3y_a = I_TWOBODYOVERLAP_K3x4y_D2y_a+ABY*I_TWOBODYOVERLAP_I3x3y_D2y_a;
  Double I_TWOBODYOVERLAP_I3x2yz_F3y_a = I_TWOBODYOVERLAP_K3x3yz_D2y_a+ABY*I_TWOBODYOVERLAP_I3x2yz_D2y_a;
  Double I_TWOBODYOVERLAP_I3xy2z_F3y_a = I_TWOBODYOVERLAP_K3x2y2z_D2y_a+ABY*I_TWOBODYOVERLAP_I3xy2z_D2y_a;
  Double I_TWOBODYOVERLAP_I3x3z_F3y_a = I_TWOBODYOVERLAP_K3xy3z_D2y_a+ABY*I_TWOBODYOVERLAP_I3x3z_D2y_a;
  Double I_TWOBODYOVERLAP_I2x4y_F3y_a = I_TWOBODYOVERLAP_K2x5y_D2y_a+ABY*I_TWOBODYOVERLAP_I2x4y_D2y_a;
  Double I_TWOBODYOVERLAP_I2x3yz_F3y_a = I_TWOBODYOVERLAP_K2x4yz_D2y_a+ABY*I_TWOBODYOVERLAP_I2x3yz_D2y_a;
  Double I_TWOBODYOVERLAP_I2x2y2z_F3y_a = I_TWOBODYOVERLAP_K2x3y2z_D2y_a+ABY*I_TWOBODYOVERLAP_I2x2y2z_D2y_a;
  Double I_TWOBODYOVERLAP_I2xy3z_F3y_a = I_TWOBODYOVERLAP_K2x2y3z_D2y_a+ABY*I_TWOBODYOVERLAP_I2xy3z_D2y_a;
  Double I_TWOBODYOVERLAP_I2x4z_F3y_a = I_TWOBODYOVERLAP_K2xy4z_D2y_a+ABY*I_TWOBODYOVERLAP_I2x4z_D2y_a;
  Double I_TWOBODYOVERLAP_Ix5y_F3y_a = I_TWOBODYOVERLAP_Kx6y_D2y_a+ABY*I_TWOBODYOVERLAP_Ix5y_D2y_a;
  Double I_TWOBODYOVERLAP_Ix4yz_F3y_a = I_TWOBODYOVERLAP_Kx5yz_D2y_a+ABY*I_TWOBODYOVERLAP_Ix4yz_D2y_a;
  Double I_TWOBODYOVERLAP_Ix3y2z_F3y_a = I_TWOBODYOVERLAP_Kx4y2z_D2y_a+ABY*I_TWOBODYOVERLAP_Ix3y2z_D2y_a;
  Double I_TWOBODYOVERLAP_Ix2y3z_F3y_a = I_TWOBODYOVERLAP_Kx3y3z_D2y_a+ABY*I_TWOBODYOVERLAP_Ix2y3z_D2y_a;
  Double I_TWOBODYOVERLAP_Ixy4z_F3y_a = I_TWOBODYOVERLAP_Kx2y4z_D2y_a+ABY*I_TWOBODYOVERLAP_Ixy4z_D2y_a;
  Double I_TWOBODYOVERLAP_Ix5z_F3y_a = I_TWOBODYOVERLAP_Kxy5z_D2y_a+ABY*I_TWOBODYOVERLAP_Ix5z_D2y_a;
  Double I_TWOBODYOVERLAP_I6y_F3y_a = I_TWOBODYOVERLAP_K7y_D2y_a+ABY*I_TWOBODYOVERLAP_I6y_D2y_a;
  Double I_TWOBODYOVERLAP_I5yz_F3y_a = I_TWOBODYOVERLAP_K6yz_D2y_a+ABY*I_TWOBODYOVERLAP_I5yz_D2y_a;
  Double I_TWOBODYOVERLAP_I4y2z_F3y_a = I_TWOBODYOVERLAP_K5y2z_D2y_a+ABY*I_TWOBODYOVERLAP_I4y2z_D2y_a;
  Double I_TWOBODYOVERLAP_I3y3z_F3y_a = I_TWOBODYOVERLAP_K4y3z_D2y_a+ABY*I_TWOBODYOVERLAP_I3y3z_D2y_a;
  Double I_TWOBODYOVERLAP_I2y4z_F3y_a = I_TWOBODYOVERLAP_K3y4z_D2y_a+ABY*I_TWOBODYOVERLAP_I2y4z_D2y_a;
  Double I_TWOBODYOVERLAP_Iy5z_F3y_a = I_TWOBODYOVERLAP_K2y5z_D2y_a+ABY*I_TWOBODYOVERLAP_Iy5z_D2y_a;
  Double I_TWOBODYOVERLAP_I6z_F3y_a = I_TWOBODYOVERLAP_Ky6z_D2y_a+ABY*I_TWOBODYOVERLAP_I6z_D2y_a;
  Double I_TWOBODYOVERLAP_I5xz_F2yz_a = I_TWOBODYOVERLAP_K5x2z_D2y_a+ABZ*I_TWOBODYOVERLAP_I5xz_D2y_a;
  Double I_TWOBODYOVERLAP_I4xyz_F2yz_a = I_TWOBODYOVERLAP_K4xy2z_D2y_a+ABZ*I_TWOBODYOVERLAP_I4xyz_D2y_a;
  Double I_TWOBODYOVERLAP_I4x2z_F2yz_a = I_TWOBODYOVERLAP_K4x3z_D2y_a+ABZ*I_TWOBODYOVERLAP_I4x2z_D2y_a;
  Double I_TWOBODYOVERLAP_I3x2yz_F2yz_a = I_TWOBODYOVERLAP_K3x2y2z_D2y_a+ABZ*I_TWOBODYOVERLAP_I3x2yz_D2y_a;
  Double I_TWOBODYOVERLAP_I3xy2z_F2yz_a = I_TWOBODYOVERLAP_K3xy3z_D2y_a+ABZ*I_TWOBODYOVERLAP_I3xy2z_D2y_a;
  Double I_TWOBODYOVERLAP_I3x3z_F2yz_a = I_TWOBODYOVERLAP_K3x4z_D2y_a+ABZ*I_TWOBODYOVERLAP_I3x3z_D2y_a;
  Double I_TWOBODYOVERLAP_I2x3yz_F2yz_a = I_TWOBODYOVERLAP_K2x3y2z_D2y_a+ABZ*I_TWOBODYOVERLAP_I2x3yz_D2y_a;
  Double I_TWOBODYOVERLAP_I2x2y2z_F2yz_a = I_TWOBODYOVERLAP_K2x2y3z_D2y_a+ABZ*I_TWOBODYOVERLAP_I2x2y2z_D2y_a;
  Double I_TWOBODYOVERLAP_I2xy3z_F2yz_a = I_TWOBODYOVERLAP_K2xy4z_D2y_a+ABZ*I_TWOBODYOVERLAP_I2xy3z_D2y_a;
  Double I_TWOBODYOVERLAP_I2x4z_F2yz_a = I_TWOBODYOVERLAP_K2x5z_D2y_a+ABZ*I_TWOBODYOVERLAP_I2x4z_D2y_a;
  Double I_TWOBODYOVERLAP_Ix4yz_F2yz_a = I_TWOBODYOVERLAP_Kx4y2z_D2y_a+ABZ*I_TWOBODYOVERLAP_Ix4yz_D2y_a;
  Double I_TWOBODYOVERLAP_Ix3y2z_F2yz_a = I_TWOBODYOVERLAP_Kx3y3z_D2y_a+ABZ*I_TWOBODYOVERLAP_Ix3y2z_D2y_a;
  Double I_TWOBODYOVERLAP_Ix2y3z_F2yz_a = I_TWOBODYOVERLAP_Kx2y4z_D2y_a+ABZ*I_TWOBODYOVERLAP_Ix2y3z_D2y_a;
  Double I_TWOBODYOVERLAP_Ixy4z_F2yz_a = I_TWOBODYOVERLAP_Kxy5z_D2y_a+ABZ*I_TWOBODYOVERLAP_Ixy4z_D2y_a;
  Double I_TWOBODYOVERLAP_Ix5z_F2yz_a = I_TWOBODYOVERLAP_Kx6z_D2y_a+ABZ*I_TWOBODYOVERLAP_Ix5z_D2y_a;
  Double I_TWOBODYOVERLAP_I5yz_F2yz_a = I_TWOBODYOVERLAP_K5y2z_D2y_a+ABZ*I_TWOBODYOVERLAP_I5yz_D2y_a;
  Double I_TWOBODYOVERLAP_I4y2z_F2yz_a = I_TWOBODYOVERLAP_K4y3z_D2y_a+ABZ*I_TWOBODYOVERLAP_I4y2z_D2y_a;
  Double I_TWOBODYOVERLAP_I3y3z_F2yz_a = I_TWOBODYOVERLAP_K3y4z_D2y_a+ABZ*I_TWOBODYOVERLAP_I3y3z_D2y_a;
  Double I_TWOBODYOVERLAP_I2y4z_F2yz_a = I_TWOBODYOVERLAP_K2y5z_D2y_a+ABZ*I_TWOBODYOVERLAP_I2y4z_D2y_a;
  Double I_TWOBODYOVERLAP_Iy5z_F2yz_a = I_TWOBODYOVERLAP_Ky6z_D2y_a+ABZ*I_TWOBODYOVERLAP_Iy5z_D2y_a;
  Double I_TWOBODYOVERLAP_I6z_F2yz_a = I_TWOBODYOVERLAP_K7z_D2y_a+ABZ*I_TWOBODYOVERLAP_I6z_D2y_a;
  Double I_TWOBODYOVERLAP_I6x_F3z_a = I_TWOBODYOVERLAP_K6xz_D2z_a+ABZ*I_TWOBODYOVERLAP_I6x_D2z_a;
  Double I_TWOBODYOVERLAP_I5xy_F3z_a = I_TWOBODYOVERLAP_K5xyz_D2z_a+ABZ*I_TWOBODYOVERLAP_I5xy_D2z_a;
  Double I_TWOBODYOVERLAP_I5xz_F3z_a = I_TWOBODYOVERLAP_K5x2z_D2z_a+ABZ*I_TWOBODYOVERLAP_I5xz_D2z_a;
  Double I_TWOBODYOVERLAP_I4x2y_F3z_a = I_TWOBODYOVERLAP_K4x2yz_D2z_a+ABZ*I_TWOBODYOVERLAP_I4x2y_D2z_a;
  Double I_TWOBODYOVERLAP_I4xyz_F3z_a = I_TWOBODYOVERLAP_K4xy2z_D2z_a+ABZ*I_TWOBODYOVERLAP_I4xyz_D2z_a;
  Double I_TWOBODYOVERLAP_I4x2z_F3z_a = I_TWOBODYOVERLAP_K4x3z_D2z_a+ABZ*I_TWOBODYOVERLAP_I4x2z_D2z_a;
  Double I_TWOBODYOVERLAP_I3x3y_F3z_a = I_TWOBODYOVERLAP_K3x3yz_D2z_a+ABZ*I_TWOBODYOVERLAP_I3x3y_D2z_a;
  Double I_TWOBODYOVERLAP_I3x2yz_F3z_a = I_TWOBODYOVERLAP_K3x2y2z_D2z_a+ABZ*I_TWOBODYOVERLAP_I3x2yz_D2z_a;
  Double I_TWOBODYOVERLAP_I3xy2z_F3z_a = I_TWOBODYOVERLAP_K3xy3z_D2z_a+ABZ*I_TWOBODYOVERLAP_I3xy2z_D2z_a;
  Double I_TWOBODYOVERLAP_I3x3z_F3z_a = I_TWOBODYOVERLAP_K3x4z_D2z_a+ABZ*I_TWOBODYOVERLAP_I3x3z_D2z_a;
  Double I_TWOBODYOVERLAP_I2x4y_F3z_a = I_TWOBODYOVERLAP_K2x4yz_D2z_a+ABZ*I_TWOBODYOVERLAP_I2x4y_D2z_a;
  Double I_TWOBODYOVERLAP_I2x3yz_F3z_a = I_TWOBODYOVERLAP_K2x3y2z_D2z_a+ABZ*I_TWOBODYOVERLAP_I2x3yz_D2z_a;
  Double I_TWOBODYOVERLAP_I2x2y2z_F3z_a = I_TWOBODYOVERLAP_K2x2y3z_D2z_a+ABZ*I_TWOBODYOVERLAP_I2x2y2z_D2z_a;
  Double I_TWOBODYOVERLAP_I2xy3z_F3z_a = I_TWOBODYOVERLAP_K2xy4z_D2z_a+ABZ*I_TWOBODYOVERLAP_I2xy3z_D2z_a;
  Double I_TWOBODYOVERLAP_I2x4z_F3z_a = I_TWOBODYOVERLAP_K2x5z_D2z_a+ABZ*I_TWOBODYOVERLAP_I2x4z_D2z_a;
  Double I_TWOBODYOVERLAP_Ix5y_F3z_a = I_TWOBODYOVERLAP_Kx5yz_D2z_a+ABZ*I_TWOBODYOVERLAP_Ix5y_D2z_a;
  Double I_TWOBODYOVERLAP_Ix4yz_F3z_a = I_TWOBODYOVERLAP_Kx4y2z_D2z_a+ABZ*I_TWOBODYOVERLAP_Ix4yz_D2z_a;
  Double I_TWOBODYOVERLAP_Ix3y2z_F3z_a = I_TWOBODYOVERLAP_Kx3y3z_D2z_a+ABZ*I_TWOBODYOVERLAP_Ix3y2z_D2z_a;
  Double I_TWOBODYOVERLAP_Ix2y3z_F3z_a = I_TWOBODYOVERLAP_Kx2y4z_D2z_a+ABZ*I_TWOBODYOVERLAP_Ix2y3z_D2z_a;
  Double I_TWOBODYOVERLAP_Ixy4z_F3z_a = I_TWOBODYOVERLAP_Kxy5z_D2z_a+ABZ*I_TWOBODYOVERLAP_Ixy4z_D2z_a;
  Double I_TWOBODYOVERLAP_Ix5z_F3z_a = I_TWOBODYOVERLAP_Kx6z_D2z_a+ABZ*I_TWOBODYOVERLAP_Ix5z_D2z_a;
  Double I_TWOBODYOVERLAP_I6y_F3z_a = I_TWOBODYOVERLAP_K6yz_D2z_a+ABZ*I_TWOBODYOVERLAP_I6y_D2z_a;
  Double I_TWOBODYOVERLAP_I5yz_F3z_a = I_TWOBODYOVERLAP_K5y2z_D2z_a+ABZ*I_TWOBODYOVERLAP_I5yz_D2z_a;
  Double I_TWOBODYOVERLAP_I4y2z_F3z_a = I_TWOBODYOVERLAP_K4y3z_D2z_a+ABZ*I_TWOBODYOVERLAP_I4y2z_D2z_a;
  Double I_TWOBODYOVERLAP_I3y3z_F3z_a = I_TWOBODYOVERLAP_K3y4z_D2z_a+ABZ*I_TWOBODYOVERLAP_I3y3z_D2z_a;
  Double I_TWOBODYOVERLAP_I2y4z_F3z_a = I_TWOBODYOVERLAP_K2y5z_D2z_a+ABZ*I_TWOBODYOVERLAP_I2y4z_D2z_a;
  Double I_TWOBODYOVERLAP_Iy5z_F3z_a = I_TWOBODYOVERLAP_Ky6z_D2z_a+ABZ*I_TWOBODYOVERLAP_Iy5z_D2z_a;
  Double I_TWOBODYOVERLAP_I6z_F3z_a = I_TWOBODYOVERLAP_K7z_D2z_a+ABZ*I_TWOBODYOVERLAP_I6z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_H_G_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_F_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_F_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_H5x_G4x_a = I_TWOBODYOVERLAP_I6x_F3x_a+ABX*I_TWOBODYOVERLAP_H5x_F3x_a;
  Double I_TWOBODYOVERLAP_H4xy_G4x_a = I_TWOBODYOVERLAP_I5xy_F3x_a+ABX*I_TWOBODYOVERLAP_H4xy_F3x_a;
  Double I_TWOBODYOVERLAP_H4xz_G4x_a = I_TWOBODYOVERLAP_I5xz_F3x_a+ABX*I_TWOBODYOVERLAP_H4xz_F3x_a;
  Double I_TWOBODYOVERLAP_H3x2y_G4x_a = I_TWOBODYOVERLAP_I4x2y_F3x_a+ABX*I_TWOBODYOVERLAP_H3x2y_F3x_a;
  Double I_TWOBODYOVERLAP_H3xyz_G4x_a = I_TWOBODYOVERLAP_I4xyz_F3x_a+ABX*I_TWOBODYOVERLAP_H3xyz_F3x_a;
  Double I_TWOBODYOVERLAP_H3x2z_G4x_a = I_TWOBODYOVERLAP_I4x2z_F3x_a+ABX*I_TWOBODYOVERLAP_H3x2z_F3x_a;
  Double I_TWOBODYOVERLAP_H2x3y_G4x_a = I_TWOBODYOVERLAP_I3x3y_F3x_a+ABX*I_TWOBODYOVERLAP_H2x3y_F3x_a;
  Double I_TWOBODYOVERLAP_H2x2yz_G4x_a = I_TWOBODYOVERLAP_I3x2yz_F3x_a+ABX*I_TWOBODYOVERLAP_H2x2yz_F3x_a;
  Double I_TWOBODYOVERLAP_H2xy2z_G4x_a = I_TWOBODYOVERLAP_I3xy2z_F3x_a+ABX*I_TWOBODYOVERLAP_H2xy2z_F3x_a;
  Double I_TWOBODYOVERLAP_H2x3z_G4x_a = I_TWOBODYOVERLAP_I3x3z_F3x_a+ABX*I_TWOBODYOVERLAP_H2x3z_F3x_a;
  Double I_TWOBODYOVERLAP_Hx4y_G4x_a = I_TWOBODYOVERLAP_I2x4y_F3x_a+ABX*I_TWOBODYOVERLAP_Hx4y_F3x_a;
  Double I_TWOBODYOVERLAP_Hx3yz_G4x_a = I_TWOBODYOVERLAP_I2x3yz_F3x_a+ABX*I_TWOBODYOVERLAP_Hx3yz_F3x_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_G4x_a = I_TWOBODYOVERLAP_I2x2y2z_F3x_a+ABX*I_TWOBODYOVERLAP_Hx2y2z_F3x_a;
  Double I_TWOBODYOVERLAP_Hxy3z_G4x_a = I_TWOBODYOVERLAP_I2xy3z_F3x_a+ABX*I_TWOBODYOVERLAP_Hxy3z_F3x_a;
  Double I_TWOBODYOVERLAP_Hx4z_G4x_a = I_TWOBODYOVERLAP_I2x4z_F3x_a+ABX*I_TWOBODYOVERLAP_Hx4z_F3x_a;
  Double I_TWOBODYOVERLAP_H5y_G4x_a = I_TWOBODYOVERLAP_Ix5y_F3x_a+ABX*I_TWOBODYOVERLAP_H5y_F3x_a;
  Double I_TWOBODYOVERLAP_H4yz_G4x_a = I_TWOBODYOVERLAP_Ix4yz_F3x_a+ABX*I_TWOBODYOVERLAP_H4yz_F3x_a;
  Double I_TWOBODYOVERLAP_H3y2z_G4x_a = I_TWOBODYOVERLAP_Ix3y2z_F3x_a+ABX*I_TWOBODYOVERLAP_H3y2z_F3x_a;
  Double I_TWOBODYOVERLAP_H2y3z_G4x_a = I_TWOBODYOVERLAP_Ix2y3z_F3x_a+ABX*I_TWOBODYOVERLAP_H2y3z_F3x_a;
  Double I_TWOBODYOVERLAP_Hy4z_G4x_a = I_TWOBODYOVERLAP_Ixy4z_F3x_a+ABX*I_TWOBODYOVERLAP_Hy4z_F3x_a;
  Double I_TWOBODYOVERLAP_H5z_G4x_a = I_TWOBODYOVERLAP_Ix5z_F3x_a+ABX*I_TWOBODYOVERLAP_H5z_F3x_a;
  Double I_TWOBODYOVERLAP_H5x_G3xy_a = I_TWOBODYOVERLAP_I5xy_F3x_a+ABY*I_TWOBODYOVERLAP_H5x_F3x_a;
  Double I_TWOBODYOVERLAP_H4xy_G3xy_a = I_TWOBODYOVERLAP_I4x2y_F3x_a+ABY*I_TWOBODYOVERLAP_H4xy_F3x_a;
  Double I_TWOBODYOVERLAP_H4xz_G3xy_a = I_TWOBODYOVERLAP_I4xyz_F3x_a+ABY*I_TWOBODYOVERLAP_H4xz_F3x_a;
  Double I_TWOBODYOVERLAP_H3x2y_G3xy_a = I_TWOBODYOVERLAP_I3x3y_F3x_a+ABY*I_TWOBODYOVERLAP_H3x2y_F3x_a;
  Double I_TWOBODYOVERLAP_H3xyz_G3xy_a = I_TWOBODYOVERLAP_I3x2yz_F3x_a+ABY*I_TWOBODYOVERLAP_H3xyz_F3x_a;
  Double I_TWOBODYOVERLAP_H3x2z_G3xy_a = I_TWOBODYOVERLAP_I3xy2z_F3x_a+ABY*I_TWOBODYOVERLAP_H3x2z_F3x_a;
  Double I_TWOBODYOVERLAP_H2x3y_G3xy_a = I_TWOBODYOVERLAP_I2x4y_F3x_a+ABY*I_TWOBODYOVERLAP_H2x3y_F3x_a;
  Double I_TWOBODYOVERLAP_H2x2yz_G3xy_a = I_TWOBODYOVERLAP_I2x3yz_F3x_a+ABY*I_TWOBODYOVERLAP_H2x2yz_F3x_a;
  Double I_TWOBODYOVERLAP_H2xy2z_G3xy_a = I_TWOBODYOVERLAP_I2x2y2z_F3x_a+ABY*I_TWOBODYOVERLAP_H2xy2z_F3x_a;
  Double I_TWOBODYOVERLAP_H2x3z_G3xy_a = I_TWOBODYOVERLAP_I2xy3z_F3x_a+ABY*I_TWOBODYOVERLAP_H2x3z_F3x_a;
  Double I_TWOBODYOVERLAP_Hx4y_G3xy_a = I_TWOBODYOVERLAP_Ix5y_F3x_a+ABY*I_TWOBODYOVERLAP_Hx4y_F3x_a;
  Double I_TWOBODYOVERLAP_Hx3yz_G3xy_a = I_TWOBODYOVERLAP_Ix4yz_F3x_a+ABY*I_TWOBODYOVERLAP_Hx3yz_F3x_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_G3xy_a = I_TWOBODYOVERLAP_Ix3y2z_F3x_a+ABY*I_TWOBODYOVERLAP_Hx2y2z_F3x_a;
  Double I_TWOBODYOVERLAP_Hxy3z_G3xy_a = I_TWOBODYOVERLAP_Ix2y3z_F3x_a+ABY*I_TWOBODYOVERLAP_Hxy3z_F3x_a;
  Double I_TWOBODYOVERLAP_Hx4z_G3xy_a = I_TWOBODYOVERLAP_Ixy4z_F3x_a+ABY*I_TWOBODYOVERLAP_Hx4z_F3x_a;
  Double I_TWOBODYOVERLAP_H5y_G3xy_a = I_TWOBODYOVERLAP_I6y_F3x_a+ABY*I_TWOBODYOVERLAP_H5y_F3x_a;
  Double I_TWOBODYOVERLAP_H4yz_G3xy_a = I_TWOBODYOVERLAP_I5yz_F3x_a+ABY*I_TWOBODYOVERLAP_H4yz_F3x_a;
  Double I_TWOBODYOVERLAP_H3y2z_G3xy_a = I_TWOBODYOVERLAP_I4y2z_F3x_a+ABY*I_TWOBODYOVERLAP_H3y2z_F3x_a;
  Double I_TWOBODYOVERLAP_H2y3z_G3xy_a = I_TWOBODYOVERLAP_I3y3z_F3x_a+ABY*I_TWOBODYOVERLAP_H2y3z_F3x_a;
  Double I_TWOBODYOVERLAP_Hy4z_G3xy_a = I_TWOBODYOVERLAP_I2y4z_F3x_a+ABY*I_TWOBODYOVERLAP_Hy4z_F3x_a;
  Double I_TWOBODYOVERLAP_H5z_G3xy_a = I_TWOBODYOVERLAP_Iy5z_F3x_a+ABY*I_TWOBODYOVERLAP_H5z_F3x_a;
  Double I_TWOBODYOVERLAP_H5x_G3xz_a = I_TWOBODYOVERLAP_I5xz_F3x_a+ABZ*I_TWOBODYOVERLAP_H5x_F3x_a;
  Double I_TWOBODYOVERLAP_H4xy_G3xz_a = I_TWOBODYOVERLAP_I4xyz_F3x_a+ABZ*I_TWOBODYOVERLAP_H4xy_F3x_a;
  Double I_TWOBODYOVERLAP_H4xz_G3xz_a = I_TWOBODYOVERLAP_I4x2z_F3x_a+ABZ*I_TWOBODYOVERLAP_H4xz_F3x_a;
  Double I_TWOBODYOVERLAP_H3x2y_G3xz_a = I_TWOBODYOVERLAP_I3x2yz_F3x_a+ABZ*I_TWOBODYOVERLAP_H3x2y_F3x_a;
  Double I_TWOBODYOVERLAP_H3xyz_G3xz_a = I_TWOBODYOVERLAP_I3xy2z_F3x_a+ABZ*I_TWOBODYOVERLAP_H3xyz_F3x_a;
  Double I_TWOBODYOVERLAP_H3x2z_G3xz_a = I_TWOBODYOVERLAP_I3x3z_F3x_a+ABZ*I_TWOBODYOVERLAP_H3x2z_F3x_a;
  Double I_TWOBODYOVERLAP_H2x3y_G3xz_a = I_TWOBODYOVERLAP_I2x3yz_F3x_a+ABZ*I_TWOBODYOVERLAP_H2x3y_F3x_a;
  Double I_TWOBODYOVERLAP_H2x2yz_G3xz_a = I_TWOBODYOVERLAP_I2x2y2z_F3x_a+ABZ*I_TWOBODYOVERLAP_H2x2yz_F3x_a;
  Double I_TWOBODYOVERLAP_H2xy2z_G3xz_a = I_TWOBODYOVERLAP_I2xy3z_F3x_a+ABZ*I_TWOBODYOVERLAP_H2xy2z_F3x_a;
  Double I_TWOBODYOVERLAP_H2x3z_G3xz_a = I_TWOBODYOVERLAP_I2x4z_F3x_a+ABZ*I_TWOBODYOVERLAP_H2x3z_F3x_a;
  Double I_TWOBODYOVERLAP_Hx4y_G3xz_a = I_TWOBODYOVERLAP_Ix4yz_F3x_a+ABZ*I_TWOBODYOVERLAP_Hx4y_F3x_a;
  Double I_TWOBODYOVERLAP_Hx3yz_G3xz_a = I_TWOBODYOVERLAP_Ix3y2z_F3x_a+ABZ*I_TWOBODYOVERLAP_Hx3yz_F3x_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_G3xz_a = I_TWOBODYOVERLAP_Ix2y3z_F3x_a+ABZ*I_TWOBODYOVERLAP_Hx2y2z_F3x_a;
  Double I_TWOBODYOVERLAP_Hxy3z_G3xz_a = I_TWOBODYOVERLAP_Ixy4z_F3x_a+ABZ*I_TWOBODYOVERLAP_Hxy3z_F3x_a;
  Double I_TWOBODYOVERLAP_Hx4z_G3xz_a = I_TWOBODYOVERLAP_Ix5z_F3x_a+ABZ*I_TWOBODYOVERLAP_Hx4z_F3x_a;
  Double I_TWOBODYOVERLAP_H5y_G3xz_a = I_TWOBODYOVERLAP_I5yz_F3x_a+ABZ*I_TWOBODYOVERLAP_H5y_F3x_a;
  Double I_TWOBODYOVERLAP_H4yz_G3xz_a = I_TWOBODYOVERLAP_I4y2z_F3x_a+ABZ*I_TWOBODYOVERLAP_H4yz_F3x_a;
  Double I_TWOBODYOVERLAP_H3y2z_G3xz_a = I_TWOBODYOVERLAP_I3y3z_F3x_a+ABZ*I_TWOBODYOVERLAP_H3y2z_F3x_a;
  Double I_TWOBODYOVERLAP_H2y3z_G3xz_a = I_TWOBODYOVERLAP_I2y4z_F3x_a+ABZ*I_TWOBODYOVERLAP_H2y3z_F3x_a;
  Double I_TWOBODYOVERLAP_Hy4z_G3xz_a = I_TWOBODYOVERLAP_Iy5z_F3x_a+ABZ*I_TWOBODYOVERLAP_Hy4z_F3x_a;
  Double I_TWOBODYOVERLAP_H5z_G3xz_a = I_TWOBODYOVERLAP_I6z_F3x_a+ABZ*I_TWOBODYOVERLAP_H5z_F3x_a;
  Double I_TWOBODYOVERLAP_H5x_G2x2y_a = I_TWOBODYOVERLAP_I5xy_F2xy_a+ABY*I_TWOBODYOVERLAP_H5x_F2xy_a;
  Double I_TWOBODYOVERLAP_H4xy_G2x2y_a = I_TWOBODYOVERLAP_I4x2y_F2xy_a+ABY*I_TWOBODYOVERLAP_H4xy_F2xy_a;
  Double I_TWOBODYOVERLAP_H4xz_G2x2y_a = I_TWOBODYOVERLAP_I4xyz_F2xy_a+ABY*I_TWOBODYOVERLAP_H4xz_F2xy_a;
  Double I_TWOBODYOVERLAP_H3x2y_G2x2y_a = I_TWOBODYOVERLAP_I3x3y_F2xy_a+ABY*I_TWOBODYOVERLAP_H3x2y_F2xy_a;
  Double I_TWOBODYOVERLAP_H3xyz_G2x2y_a = I_TWOBODYOVERLAP_I3x2yz_F2xy_a+ABY*I_TWOBODYOVERLAP_H3xyz_F2xy_a;
  Double I_TWOBODYOVERLAP_H3x2z_G2x2y_a = I_TWOBODYOVERLAP_I3xy2z_F2xy_a+ABY*I_TWOBODYOVERLAP_H3x2z_F2xy_a;
  Double I_TWOBODYOVERLAP_H2x3y_G2x2y_a = I_TWOBODYOVERLAP_I2x4y_F2xy_a+ABY*I_TWOBODYOVERLAP_H2x3y_F2xy_a;
  Double I_TWOBODYOVERLAP_H2x2yz_G2x2y_a = I_TWOBODYOVERLAP_I2x3yz_F2xy_a+ABY*I_TWOBODYOVERLAP_H2x2yz_F2xy_a;
  Double I_TWOBODYOVERLAP_H2xy2z_G2x2y_a = I_TWOBODYOVERLAP_I2x2y2z_F2xy_a+ABY*I_TWOBODYOVERLAP_H2xy2z_F2xy_a;
  Double I_TWOBODYOVERLAP_H2x3z_G2x2y_a = I_TWOBODYOVERLAP_I2xy3z_F2xy_a+ABY*I_TWOBODYOVERLAP_H2x3z_F2xy_a;
  Double I_TWOBODYOVERLAP_Hx4y_G2x2y_a = I_TWOBODYOVERLAP_Ix5y_F2xy_a+ABY*I_TWOBODYOVERLAP_Hx4y_F2xy_a;
  Double I_TWOBODYOVERLAP_Hx3yz_G2x2y_a = I_TWOBODYOVERLAP_Ix4yz_F2xy_a+ABY*I_TWOBODYOVERLAP_Hx3yz_F2xy_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_G2x2y_a = I_TWOBODYOVERLAP_Ix3y2z_F2xy_a+ABY*I_TWOBODYOVERLAP_Hx2y2z_F2xy_a;
  Double I_TWOBODYOVERLAP_Hxy3z_G2x2y_a = I_TWOBODYOVERLAP_Ix2y3z_F2xy_a+ABY*I_TWOBODYOVERLAP_Hxy3z_F2xy_a;
  Double I_TWOBODYOVERLAP_Hx4z_G2x2y_a = I_TWOBODYOVERLAP_Ixy4z_F2xy_a+ABY*I_TWOBODYOVERLAP_Hx4z_F2xy_a;
  Double I_TWOBODYOVERLAP_H5y_G2x2y_a = I_TWOBODYOVERLAP_I6y_F2xy_a+ABY*I_TWOBODYOVERLAP_H5y_F2xy_a;
  Double I_TWOBODYOVERLAP_H4yz_G2x2y_a = I_TWOBODYOVERLAP_I5yz_F2xy_a+ABY*I_TWOBODYOVERLAP_H4yz_F2xy_a;
  Double I_TWOBODYOVERLAP_H3y2z_G2x2y_a = I_TWOBODYOVERLAP_I4y2z_F2xy_a+ABY*I_TWOBODYOVERLAP_H3y2z_F2xy_a;
  Double I_TWOBODYOVERLAP_H2y3z_G2x2y_a = I_TWOBODYOVERLAP_I3y3z_F2xy_a+ABY*I_TWOBODYOVERLAP_H2y3z_F2xy_a;
  Double I_TWOBODYOVERLAP_Hy4z_G2x2y_a = I_TWOBODYOVERLAP_I2y4z_F2xy_a+ABY*I_TWOBODYOVERLAP_Hy4z_F2xy_a;
  Double I_TWOBODYOVERLAP_H5z_G2x2y_a = I_TWOBODYOVERLAP_Iy5z_F2xy_a+ABY*I_TWOBODYOVERLAP_H5z_F2xy_a;
  Double I_TWOBODYOVERLAP_H5x_G2xyz_a = I_TWOBODYOVERLAP_I5xz_F2xy_a+ABZ*I_TWOBODYOVERLAP_H5x_F2xy_a;
  Double I_TWOBODYOVERLAP_H4xy_G2xyz_a = I_TWOBODYOVERLAP_I4xyz_F2xy_a+ABZ*I_TWOBODYOVERLAP_H4xy_F2xy_a;
  Double I_TWOBODYOVERLAP_H4xz_G2xyz_a = I_TWOBODYOVERLAP_I4x2z_F2xy_a+ABZ*I_TWOBODYOVERLAP_H4xz_F2xy_a;
  Double I_TWOBODYOVERLAP_H3x2y_G2xyz_a = I_TWOBODYOVERLAP_I3x2yz_F2xy_a+ABZ*I_TWOBODYOVERLAP_H3x2y_F2xy_a;
  Double I_TWOBODYOVERLAP_H3xyz_G2xyz_a = I_TWOBODYOVERLAP_I3xy2z_F2xy_a+ABZ*I_TWOBODYOVERLAP_H3xyz_F2xy_a;
  Double I_TWOBODYOVERLAP_H3x2z_G2xyz_a = I_TWOBODYOVERLAP_I3x3z_F2xy_a+ABZ*I_TWOBODYOVERLAP_H3x2z_F2xy_a;
  Double I_TWOBODYOVERLAP_H2x3y_G2xyz_a = I_TWOBODYOVERLAP_I2x3yz_F2xy_a+ABZ*I_TWOBODYOVERLAP_H2x3y_F2xy_a;
  Double I_TWOBODYOVERLAP_H2x2yz_G2xyz_a = I_TWOBODYOVERLAP_I2x2y2z_F2xy_a+ABZ*I_TWOBODYOVERLAP_H2x2yz_F2xy_a;
  Double I_TWOBODYOVERLAP_H2xy2z_G2xyz_a = I_TWOBODYOVERLAP_I2xy3z_F2xy_a+ABZ*I_TWOBODYOVERLAP_H2xy2z_F2xy_a;
  Double I_TWOBODYOVERLAP_H2x3z_G2xyz_a = I_TWOBODYOVERLAP_I2x4z_F2xy_a+ABZ*I_TWOBODYOVERLAP_H2x3z_F2xy_a;
  Double I_TWOBODYOVERLAP_Hx4y_G2xyz_a = I_TWOBODYOVERLAP_Ix4yz_F2xy_a+ABZ*I_TWOBODYOVERLAP_Hx4y_F2xy_a;
  Double I_TWOBODYOVERLAP_Hx3yz_G2xyz_a = I_TWOBODYOVERLAP_Ix3y2z_F2xy_a+ABZ*I_TWOBODYOVERLAP_Hx3yz_F2xy_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_G2xyz_a = I_TWOBODYOVERLAP_Ix2y3z_F2xy_a+ABZ*I_TWOBODYOVERLAP_Hx2y2z_F2xy_a;
  Double I_TWOBODYOVERLAP_Hxy3z_G2xyz_a = I_TWOBODYOVERLAP_Ixy4z_F2xy_a+ABZ*I_TWOBODYOVERLAP_Hxy3z_F2xy_a;
  Double I_TWOBODYOVERLAP_Hx4z_G2xyz_a = I_TWOBODYOVERLAP_Ix5z_F2xy_a+ABZ*I_TWOBODYOVERLAP_Hx4z_F2xy_a;
  Double I_TWOBODYOVERLAP_H5y_G2xyz_a = I_TWOBODYOVERLAP_I5yz_F2xy_a+ABZ*I_TWOBODYOVERLAP_H5y_F2xy_a;
  Double I_TWOBODYOVERLAP_H4yz_G2xyz_a = I_TWOBODYOVERLAP_I4y2z_F2xy_a+ABZ*I_TWOBODYOVERLAP_H4yz_F2xy_a;
  Double I_TWOBODYOVERLAP_H3y2z_G2xyz_a = I_TWOBODYOVERLAP_I3y3z_F2xy_a+ABZ*I_TWOBODYOVERLAP_H3y2z_F2xy_a;
  Double I_TWOBODYOVERLAP_H2y3z_G2xyz_a = I_TWOBODYOVERLAP_I2y4z_F2xy_a+ABZ*I_TWOBODYOVERLAP_H2y3z_F2xy_a;
  Double I_TWOBODYOVERLAP_Hy4z_G2xyz_a = I_TWOBODYOVERLAP_Iy5z_F2xy_a+ABZ*I_TWOBODYOVERLAP_Hy4z_F2xy_a;
  Double I_TWOBODYOVERLAP_H5z_G2xyz_a = I_TWOBODYOVERLAP_I6z_F2xy_a+ABZ*I_TWOBODYOVERLAP_H5z_F2xy_a;
  Double I_TWOBODYOVERLAP_H5x_G2x2z_a = I_TWOBODYOVERLAP_I5xz_F2xz_a+ABZ*I_TWOBODYOVERLAP_H5x_F2xz_a;
  Double I_TWOBODYOVERLAP_H4xy_G2x2z_a = I_TWOBODYOVERLAP_I4xyz_F2xz_a+ABZ*I_TWOBODYOVERLAP_H4xy_F2xz_a;
  Double I_TWOBODYOVERLAP_H4xz_G2x2z_a = I_TWOBODYOVERLAP_I4x2z_F2xz_a+ABZ*I_TWOBODYOVERLAP_H4xz_F2xz_a;
  Double I_TWOBODYOVERLAP_H3x2y_G2x2z_a = I_TWOBODYOVERLAP_I3x2yz_F2xz_a+ABZ*I_TWOBODYOVERLAP_H3x2y_F2xz_a;
  Double I_TWOBODYOVERLAP_H3xyz_G2x2z_a = I_TWOBODYOVERLAP_I3xy2z_F2xz_a+ABZ*I_TWOBODYOVERLAP_H3xyz_F2xz_a;
  Double I_TWOBODYOVERLAP_H3x2z_G2x2z_a = I_TWOBODYOVERLAP_I3x3z_F2xz_a+ABZ*I_TWOBODYOVERLAP_H3x2z_F2xz_a;
  Double I_TWOBODYOVERLAP_H2x3y_G2x2z_a = I_TWOBODYOVERLAP_I2x3yz_F2xz_a+ABZ*I_TWOBODYOVERLAP_H2x3y_F2xz_a;
  Double I_TWOBODYOVERLAP_H2x2yz_G2x2z_a = I_TWOBODYOVERLAP_I2x2y2z_F2xz_a+ABZ*I_TWOBODYOVERLAP_H2x2yz_F2xz_a;
  Double I_TWOBODYOVERLAP_H2xy2z_G2x2z_a = I_TWOBODYOVERLAP_I2xy3z_F2xz_a+ABZ*I_TWOBODYOVERLAP_H2xy2z_F2xz_a;
  Double I_TWOBODYOVERLAP_H2x3z_G2x2z_a = I_TWOBODYOVERLAP_I2x4z_F2xz_a+ABZ*I_TWOBODYOVERLAP_H2x3z_F2xz_a;
  Double I_TWOBODYOVERLAP_Hx4y_G2x2z_a = I_TWOBODYOVERLAP_Ix4yz_F2xz_a+ABZ*I_TWOBODYOVERLAP_Hx4y_F2xz_a;
  Double I_TWOBODYOVERLAP_Hx3yz_G2x2z_a = I_TWOBODYOVERLAP_Ix3y2z_F2xz_a+ABZ*I_TWOBODYOVERLAP_Hx3yz_F2xz_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_G2x2z_a = I_TWOBODYOVERLAP_Ix2y3z_F2xz_a+ABZ*I_TWOBODYOVERLAP_Hx2y2z_F2xz_a;
  Double I_TWOBODYOVERLAP_Hxy3z_G2x2z_a = I_TWOBODYOVERLAP_Ixy4z_F2xz_a+ABZ*I_TWOBODYOVERLAP_Hxy3z_F2xz_a;
  Double I_TWOBODYOVERLAP_Hx4z_G2x2z_a = I_TWOBODYOVERLAP_Ix5z_F2xz_a+ABZ*I_TWOBODYOVERLAP_Hx4z_F2xz_a;
  Double I_TWOBODYOVERLAP_H5y_G2x2z_a = I_TWOBODYOVERLAP_I5yz_F2xz_a+ABZ*I_TWOBODYOVERLAP_H5y_F2xz_a;
  Double I_TWOBODYOVERLAP_H4yz_G2x2z_a = I_TWOBODYOVERLAP_I4y2z_F2xz_a+ABZ*I_TWOBODYOVERLAP_H4yz_F2xz_a;
  Double I_TWOBODYOVERLAP_H3y2z_G2x2z_a = I_TWOBODYOVERLAP_I3y3z_F2xz_a+ABZ*I_TWOBODYOVERLAP_H3y2z_F2xz_a;
  Double I_TWOBODYOVERLAP_H2y3z_G2x2z_a = I_TWOBODYOVERLAP_I2y4z_F2xz_a+ABZ*I_TWOBODYOVERLAP_H2y3z_F2xz_a;
  Double I_TWOBODYOVERLAP_Hy4z_G2x2z_a = I_TWOBODYOVERLAP_Iy5z_F2xz_a+ABZ*I_TWOBODYOVERLAP_Hy4z_F2xz_a;
  Double I_TWOBODYOVERLAP_H5z_G2x2z_a = I_TWOBODYOVERLAP_I6z_F2xz_a+ABZ*I_TWOBODYOVERLAP_H5z_F2xz_a;
  Double I_TWOBODYOVERLAP_H5x_Gx3y_a = I_TWOBODYOVERLAP_I6x_F3y_a+ABX*I_TWOBODYOVERLAP_H5x_F3y_a;
  Double I_TWOBODYOVERLAP_H4xy_Gx3y_a = I_TWOBODYOVERLAP_I5xy_F3y_a+ABX*I_TWOBODYOVERLAP_H4xy_F3y_a;
  Double I_TWOBODYOVERLAP_H4xz_Gx3y_a = I_TWOBODYOVERLAP_I5xz_F3y_a+ABX*I_TWOBODYOVERLAP_H4xz_F3y_a;
  Double I_TWOBODYOVERLAP_H3x2y_Gx3y_a = I_TWOBODYOVERLAP_I4x2y_F3y_a+ABX*I_TWOBODYOVERLAP_H3x2y_F3y_a;
  Double I_TWOBODYOVERLAP_H3xyz_Gx3y_a = I_TWOBODYOVERLAP_I4xyz_F3y_a+ABX*I_TWOBODYOVERLAP_H3xyz_F3y_a;
  Double I_TWOBODYOVERLAP_H3x2z_Gx3y_a = I_TWOBODYOVERLAP_I4x2z_F3y_a+ABX*I_TWOBODYOVERLAP_H3x2z_F3y_a;
  Double I_TWOBODYOVERLAP_H2x3y_Gx3y_a = I_TWOBODYOVERLAP_I3x3y_F3y_a+ABX*I_TWOBODYOVERLAP_H2x3y_F3y_a;
  Double I_TWOBODYOVERLAP_H2x2yz_Gx3y_a = I_TWOBODYOVERLAP_I3x2yz_F3y_a+ABX*I_TWOBODYOVERLAP_H2x2yz_F3y_a;
  Double I_TWOBODYOVERLAP_H2xy2z_Gx3y_a = I_TWOBODYOVERLAP_I3xy2z_F3y_a+ABX*I_TWOBODYOVERLAP_H2xy2z_F3y_a;
  Double I_TWOBODYOVERLAP_H2x3z_Gx3y_a = I_TWOBODYOVERLAP_I3x3z_F3y_a+ABX*I_TWOBODYOVERLAP_H2x3z_F3y_a;
  Double I_TWOBODYOVERLAP_Hx4y_Gx3y_a = I_TWOBODYOVERLAP_I2x4y_F3y_a+ABX*I_TWOBODYOVERLAP_Hx4y_F3y_a;
  Double I_TWOBODYOVERLAP_Hx3yz_Gx3y_a = I_TWOBODYOVERLAP_I2x3yz_F3y_a+ABX*I_TWOBODYOVERLAP_Hx3yz_F3y_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_Gx3y_a = I_TWOBODYOVERLAP_I2x2y2z_F3y_a+ABX*I_TWOBODYOVERLAP_Hx2y2z_F3y_a;
  Double I_TWOBODYOVERLAP_Hxy3z_Gx3y_a = I_TWOBODYOVERLAP_I2xy3z_F3y_a+ABX*I_TWOBODYOVERLAP_Hxy3z_F3y_a;
  Double I_TWOBODYOVERLAP_Hx4z_Gx3y_a = I_TWOBODYOVERLAP_I2x4z_F3y_a+ABX*I_TWOBODYOVERLAP_Hx4z_F3y_a;
  Double I_TWOBODYOVERLAP_H5y_Gx3y_a = I_TWOBODYOVERLAP_Ix5y_F3y_a+ABX*I_TWOBODYOVERLAP_H5y_F3y_a;
  Double I_TWOBODYOVERLAP_H4yz_Gx3y_a = I_TWOBODYOVERLAP_Ix4yz_F3y_a+ABX*I_TWOBODYOVERLAP_H4yz_F3y_a;
  Double I_TWOBODYOVERLAP_H3y2z_Gx3y_a = I_TWOBODYOVERLAP_Ix3y2z_F3y_a+ABX*I_TWOBODYOVERLAP_H3y2z_F3y_a;
  Double I_TWOBODYOVERLAP_H2y3z_Gx3y_a = I_TWOBODYOVERLAP_Ix2y3z_F3y_a+ABX*I_TWOBODYOVERLAP_H2y3z_F3y_a;
  Double I_TWOBODYOVERLAP_Hy4z_Gx3y_a = I_TWOBODYOVERLAP_Ixy4z_F3y_a+ABX*I_TWOBODYOVERLAP_Hy4z_F3y_a;
  Double I_TWOBODYOVERLAP_H5z_Gx3y_a = I_TWOBODYOVERLAP_Ix5z_F3y_a+ABX*I_TWOBODYOVERLAP_H5z_F3y_a;
  Double I_TWOBODYOVERLAP_H5x_Gx2yz_a = I_TWOBODYOVERLAP_I5xz_Fx2y_a+ABZ*I_TWOBODYOVERLAP_H5x_Fx2y_a;
  Double I_TWOBODYOVERLAP_H4xy_Gx2yz_a = I_TWOBODYOVERLAP_I4xyz_Fx2y_a+ABZ*I_TWOBODYOVERLAP_H4xy_Fx2y_a;
  Double I_TWOBODYOVERLAP_H4xz_Gx2yz_a = I_TWOBODYOVERLAP_I4x2z_Fx2y_a+ABZ*I_TWOBODYOVERLAP_H4xz_Fx2y_a;
  Double I_TWOBODYOVERLAP_H3x2y_Gx2yz_a = I_TWOBODYOVERLAP_I3x2yz_Fx2y_a+ABZ*I_TWOBODYOVERLAP_H3x2y_Fx2y_a;
  Double I_TWOBODYOVERLAP_H3xyz_Gx2yz_a = I_TWOBODYOVERLAP_I3xy2z_Fx2y_a+ABZ*I_TWOBODYOVERLAP_H3xyz_Fx2y_a;
  Double I_TWOBODYOVERLAP_H3x2z_Gx2yz_a = I_TWOBODYOVERLAP_I3x3z_Fx2y_a+ABZ*I_TWOBODYOVERLAP_H3x2z_Fx2y_a;
  Double I_TWOBODYOVERLAP_H2x3y_Gx2yz_a = I_TWOBODYOVERLAP_I2x3yz_Fx2y_a+ABZ*I_TWOBODYOVERLAP_H2x3y_Fx2y_a;
  Double I_TWOBODYOVERLAP_H2x2yz_Gx2yz_a = I_TWOBODYOVERLAP_I2x2y2z_Fx2y_a+ABZ*I_TWOBODYOVERLAP_H2x2yz_Fx2y_a;
  Double I_TWOBODYOVERLAP_H2xy2z_Gx2yz_a = I_TWOBODYOVERLAP_I2xy3z_Fx2y_a+ABZ*I_TWOBODYOVERLAP_H2xy2z_Fx2y_a;
  Double I_TWOBODYOVERLAP_H2x3z_Gx2yz_a = I_TWOBODYOVERLAP_I2x4z_Fx2y_a+ABZ*I_TWOBODYOVERLAP_H2x3z_Fx2y_a;
  Double I_TWOBODYOVERLAP_Hx4y_Gx2yz_a = I_TWOBODYOVERLAP_Ix4yz_Fx2y_a+ABZ*I_TWOBODYOVERLAP_Hx4y_Fx2y_a;
  Double I_TWOBODYOVERLAP_Hx3yz_Gx2yz_a = I_TWOBODYOVERLAP_Ix3y2z_Fx2y_a+ABZ*I_TWOBODYOVERLAP_Hx3yz_Fx2y_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_Gx2yz_a = I_TWOBODYOVERLAP_Ix2y3z_Fx2y_a+ABZ*I_TWOBODYOVERLAP_Hx2y2z_Fx2y_a;
  Double I_TWOBODYOVERLAP_Hxy3z_Gx2yz_a = I_TWOBODYOVERLAP_Ixy4z_Fx2y_a+ABZ*I_TWOBODYOVERLAP_Hxy3z_Fx2y_a;
  Double I_TWOBODYOVERLAP_Hx4z_Gx2yz_a = I_TWOBODYOVERLAP_Ix5z_Fx2y_a+ABZ*I_TWOBODYOVERLAP_Hx4z_Fx2y_a;
  Double I_TWOBODYOVERLAP_H5y_Gx2yz_a = I_TWOBODYOVERLAP_I5yz_Fx2y_a+ABZ*I_TWOBODYOVERLAP_H5y_Fx2y_a;
  Double I_TWOBODYOVERLAP_H4yz_Gx2yz_a = I_TWOBODYOVERLAP_I4y2z_Fx2y_a+ABZ*I_TWOBODYOVERLAP_H4yz_Fx2y_a;
  Double I_TWOBODYOVERLAP_H3y2z_Gx2yz_a = I_TWOBODYOVERLAP_I3y3z_Fx2y_a+ABZ*I_TWOBODYOVERLAP_H3y2z_Fx2y_a;
  Double I_TWOBODYOVERLAP_H2y3z_Gx2yz_a = I_TWOBODYOVERLAP_I2y4z_Fx2y_a+ABZ*I_TWOBODYOVERLAP_H2y3z_Fx2y_a;
  Double I_TWOBODYOVERLAP_Hy4z_Gx2yz_a = I_TWOBODYOVERLAP_Iy5z_Fx2y_a+ABZ*I_TWOBODYOVERLAP_Hy4z_Fx2y_a;
  Double I_TWOBODYOVERLAP_H5z_Gx2yz_a = I_TWOBODYOVERLAP_I6z_Fx2y_a+ABZ*I_TWOBODYOVERLAP_H5z_Fx2y_a;
  Double I_TWOBODYOVERLAP_H5x_Gxy2z_a = I_TWOBODYOVERLAP_I5xy_Fx2z_a+ABY*I_TWOBODYOVERLAP_H5x_Fx2z_a;
  Double I_TWOBODYOVERLAP_H4xy_Gxy2z_a = I_TWOBODYOVERLAP_I4x2y_Fx2z_a+ABY*I_TWOBODYOVERLAP_H4xy_Fx2z_a;
  Double I_TWOBODYOVERLAP_H4xz_Gxy2z_a = I_TWOBODYOVERLAP_I4xyz_Fx2z_a+ABY*I_TWOBODYOVERLAP_H4xz_Fx2z_a;
  Double I_TWOBODYOVERLAP_H3x2y_Gxy2z_a = I_TWOBODYOVERLAP_I3x3y_Fx2z_a+ABY*I_TWOBODYOVERLAP_H3x2y_Fx2z_a;
  Double I_TWOBODYOVERLAP_H3xyz_Gxy2z_a = I_TWOBODYOVERLAP_I3x2yz_Fx2z_a+ABY*I_TWOBODYOVERLAP_H3xyz_Fx2z_a;
  Double I_TWOBODYOVERLAP_H3x2z_Gxy2z_a = I_TWOBODYOVERLAP_I3xy2z_Fx2z_a+ABY*I_TWOBODYOVERLAP_H3x2z_Fx2z_a;
  Double I_TWOBODYOVERLAP_H2x3y_Gxy2z_a = I_TWOBODYOVERLAP_I2x4y_Fx2z_a+ABY*I_TWOBODYOVERLAP_H2x3y_Fx2z_a;
  Double I_TWOBODYOVERLAP_H2x2yz_Gxy2z_a = I_TWOBODYOVERLAP_I2x3yz_Fx2z_a+ABY*I_TWOBODYOVERLAP_H2x2yz_Fx2z_a;
  Double I_TWOBODYOVERLAP_H2xy2z_Gxy2z_a = I_TWOBODYOVERLAP_I2x2y2z_Fx2z_a+ABY*I_TWOBODYOVERLAP_H2xy2z_Fx2z_a;
  Double I_TWOBODYOVERLAP_H2x3z_Gxy2z_a = I_TWOBODYOVERLAP_I2xy3z_Fx2z_a+ABY*I_TWOBODYOVERLAP_H2x3z_Fx2z_a;
  Double I_TWOBODYOVERLAP_Hx4y_Gxy2z_a = I_TWOBODYOVERLAP_Ix5y_Fx2z_a+ABY*I_TWOBODYOVERLAP_Hx4y_Fx2z_a;
  Double I_TWOBODYOVERLAP_Hx3yz_Gxy2z_a = I_TWOBODYOVERLAP_Ix4yz_Fx2z_a+ABY*I_TWOBODYOVERLAP_Hx3yz_Fx2z_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_Gxy2z_a = I_TWOBODYOVERLAP_Ix3y2z_Fx2z_a+ABY*I_TWOBODYOVERLAP_Hx2y2z_Fx2z_a;
  Double I_TWOBODYOVERLAP_Hxy3z_Gxy2z_a = I_TWOBODYOVERLAP_Ix2y3z_Fx2z_a+ABY*I_TWOBODYOVERLAP_Hxy3z_Fx2z_a;
  Double I_TWOBODYOVERLAP_Hx4z_Gxy2z_a = I_TWOBODYOVERLAP_Ixy4z_Fx2z_a+ABY*I_TWOBODYOVERLAP_Hx4z_Fx2z_a;
  Double I_TWOBODYOVERLAP_H5y_Gxy2z_a = I_TWOBODYOVERLAP_I6y_Fx2z_a+ABY*I_TWOBODYOVERLAP_H5y_Fx2z_a;
  Double I_TWOBODYOVERLAP_H4yz_Gxy2z_a = I_TWOBODYOVERLAP_I5yz_Fx2z_a+ABY*I_TWOBODYOVERLAP_H4yz_Fx2z_a;
  Double I_TWOBODYOVERLAP_H3y2z_Gxy2z_a = I_TWOBODYOVERLAP_I4y2z_Fx2z_a+ABY*I_TWOBODYOVERLAP_H3y2z_Fx2z_a;
  Double I_TWOBODYOVERLAP_H2y3z_Gxy2z_a = I_TWOBODYOVERLAP_I3y3z_Fx2z_a+ABY*I_TWOBODYOVERLAP_H2y3z_Fx2z_a;
  Double I_TWOBODYOVERLAP_Hy4z_Gxy2z_a = I_TWOBODYOVERLAP_I2y4z_Fx2z_a+ABY*I_TWOBODYOVERLAP_Hy4z_Fx2z_a;
  Double I_TWOBODYOVERLAP_H5z_Gxy2z_a = I_TWOBODYOVERLAP_Iy5z_Fx2z_a+ABY*I_TWOBODYOVERLAP_H5z_Fx2z_a;
  Double I_TWOBODYOVERLAP_H5x_Gx3z_a = I_TWOBODYOVERLAP_I6x_F3z_a+ABX*I_TWOBODYOVERLAP_H5x_F3z_a;
  Double I_TWOBODYOVERLAP_H4xy_Gx3z_a = I_TWOBODYOVERLAP_I5xy_F3z_a+ABX*I_TWOBODYOVERLAP_H4xy_F3z_a;
  Double I_TWOBODYOVERLAP_H4xz_Gx3z_a = I_TWOBODYOVERLAP_I5xz_F3z_a+ABX*I_TWOBODYOVERLAP_H4xz_F3z_a;
  Double I_TWOBODYOVERLAP_H3x2y_Gx3z_a = I_TWOBODYOVERLAP_I4x2y_F3z_a+ABX*I_TWOBODYOVERLAP_H3x2y_F3z_a;
  Double I_TWOBODYOVERLAP_H3xyz_Gx3z_a = I_TWOBODYOVERLAP_I4xyz_F3z_a+ABX*I_TWOBODYOVERLAP_H3xyz_F3z_a;
  Double I_TWOBODYOVERLAP_H3x2z_Gx3z_a = I_TWOBODYOVERLAP_I4x2z_F3z_a+ABX*I_TWOBODYOVERLAP_H3x2z_F3z_a;
  Double I_TWOBODYOVERLAP_H2x3y_Gx3z_a = I_TWOBODYOVERLAP_I3x3y_F3z_a+ABX*I_TWOBODYOVERLAP_H2x3y_F3z_a;
  Double I_TWOBODYOVERLAP_H2x2yz_Gx3z_a = I_TWOBODYOVERLAP_I3x2yz_F3z_a+ABX*I_TWOBODYOVERLAP_H2x2yz_F3z_a;
  Double I_TWOBODYOVERLAP_H2xy2z_Gx3z_a = I_TWOBODYOVERLAP_I3xy2z_F3z_a+ABX*I_TWOBODYOVERLAP_H2xy2z_F3z_a;
  Double I_TWOBODYOVERLAP_H2x3z_Gx3z_a = I_TWOBODYOVERLAP_I3x3z_F3z_a+ABX*I_TWOBODYOVERLAP_H2x3z_F3z_a;
  Double I_TWOBODYOVERLAP_Hx4y_Gx3z_a = I_TWOBODYOVERLAP_I2x4y_F3z_a+ABX*I_TWOBODYOVERLAP_Hx4y_F3z_a;
  Double I_TWOBODYOVERLAP_Hx3yz_Gx3z_a = I_TWOBODYOVERLAP_I2x3yz_F3z_a+ABX*I_TWOBODYOVERLAP_Hx3yz_F3z_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_Gx3z_a = I_TWOBODYOVERLAP_I2x2y2z_F3z_a+ABX*I_TWOBODYOVERLAP_Hx2y2z_F3z_a;
  Double I_TWOBODYOVERLAP_Hxy3z_Gx3z_a = I_TWOBODYOVERLAP_I2xy3z_F3z_a+ABX*I_TWOBODYOVERLAP_Hxy3z_F3z_a;
  Double I_TWOBODYOVERLAP_Hx4z_Gx3z_a = I_TWOBODYOVERLAP_I2x4z_F3z_a+ABX*I_TWOBODYOVERLAP_Hx4z_F3z_a;
  Double I_TWOBODYOVERLAP_H5y_Gx3z_a = I_TWOBODYOVERLAP_Ix5y_F3z_a+ABX*I_TWOBODYOVERLAP_H5y_F3z_a;
  Double I_TWOBODYOVERLAP_H4yz_Gx3z_a = I_TWOBODYOVERLAP_Ix4yz_F3z_a+ABX*I_TWOBODYOVERLAP_H4yz_F3z_a;
  Double I_TWOBODYOVERLAP_H3y2z_Gx3z_a = I_TWOBODYOVERLAP_Ix3y2z_F3z_a+ABX*I_TWOBODYOVERLAP_H3y2z_F3z_a;
  Double I_TWOBODYOVERLAP_H2y3z_Gx3z_a = I_TWOBODYOVERLAP_Ix2y3z_F3z_a+ABX*I_TWOBODYOVERLAP_H2y3z_F3z_a;
  Double I_TWOBODYOVERLAP_Hy4z_Gx3z_a = I_TWOBODYOVERLAP_Ixy4z_F3z_a+ABX*I_TWOBODYOVERLAP_Hy4z_F3z_a;
  Double I_TWOBODYOVERLAP_H5z_Gx3z_a = I_TWOBODYOVERLAP_Ix5z_F3z_a+ABX*I_TWOBODYOVERLAP_H5z_F3z_a;
  Double I_TWOBODYOVERLAP_H5x_G4y_a = I_TWOBODYOVERLAP_I5xy_F3y_a+ABY*I_TWOBODYOVERLAP_H5x_F3y_a;
  Double I_TWOBODYOVERLAP_H4xy_G4y_a = I_TWOBODYOVERLAP_I4x2y_F3y_a+ABY*I_TWOBODYOVERLAP_H4xy_F3y_a;
  Double I_TWOBODYOVERLAP_H4xz_G4y_a = I_TWOBODYOVERLAP_I4xyz_F3y_a+ABY*I_TWOBODYOVERLAP_H4xz_F3y_a;
  Double I_TWOBODYOVERLAP_H3x2y_G4y_a = I_TWOBODYOVERLAP_I3x3y_F3y_a+ABY*I_TWOBODYOVERLAP_H3x2y_F3y_a;
  Double I_TWOBODYOVERLAP_H3xyz_G4y_a = I_TWOBODYOVERLAP_I3x2yz_F3y_a+ABY*I_TWOBODYOVERLAP_H3xyz_F3y_a;
  Double I_TWOBODYOVERLAP_H3x2z_G4y_a = I_TWOBODYOVERLAP_I3xy2z_F3y_a+ABY*I_TWOBODYOVERLAP_H3x2z_F3y_a;
  Double I_TWOBODYOVERLAP_H2x3y_G4y_a = I_TWOBODYOVERLAP_I2x4y_F3y_a+ABY*I_TWOBODYOVERLAP_H2x3y_F3y_a;
  Double I_TWOBODYOVERLAP_H2x2yz_G4y_a = I_TWOBODYOVERLAP_I2x3yz_F3y_a+ABY*I_TWOBODYOVERLAP_H2x2yz_F3y_a;
  Double I_TWOBODYOVERLAP_H2xy2z_G4y_a = I_TWOBODYOVERLAP_I2x2y2z_F3y_a+ABY*I_TWOBODYOVERLAP_H2xy2z_F3y_a;
  Double I_TWOBODYOVERLAP_H2x3z_G4y_a = I_TWOBODYOVERLAP_I2xy3z_F3y_a+ABY*I_TWOBODYOVERLAP_H2x3z_F3y_a;
  Double I_TWOBODYOVERLAP_Hx4y_G4y_a = I_TWOBODYOVERLAP_Ix5y_F3y_a+ABY*I_TWOBODYOVERLAP_Hx4y_F3y_a;
  Double I_TWOBODYOVERLAP_Hx3yz_G4y_a = I_TWOBODYOVERLAP_Ix4yz_F3y_a+ABY*I_TWOBODYOVERLAP_Hx3yz_F3y_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_G4y_a = I_TWOBODYOVERLAP_Ix3y2z_F3y_a+ABY*I_TWOBODYOVERLAP_Hx2y2z_F3y_a;
  Double I_TWOBODYOVERLAP_Hxy3z_G4y_a = I_TWOBODYOVERLAP_Ix2y3z_F3y_a+ABY*I_TWOBODYOVERLAP_Hxy3z_F3y_a;
  Double I_TWOBODYOVERLAP_Hx4z_G4y_a = I_TWOBODYOVERLAP_Ixy4z_F3y_a+ABY*I_TWOBODYOVERLAP_Hx4z_F3y_a;
  Double I_TWOBODYOVERLAP_H5y_G4y_a = I_TWOBODYOVERLAP_I6y_F3y_a+ABY*I_TWOBODYOVERLAP_H5y_F3y_a;
  Double I_TWOBODYOVERLAP_H4yz_G4y_a = I_TWOBODYOVERLAP_I5yz_F3y_a+ABY*I_TWOBODYOVERLAP_H4yz_F3y_a;
  Double I_TWOBODYOVERLAP_H3y2z_G4y_a = I_TWOBODYOVERLAP_I4y2z_F3y_a+ABY*I_TWOBODYOVERLAP_H3y2z_F3y_a;
  Double I_TWOBODYOVERLAP_H2y3z_G4y_a = I_TWOBODYOVERLAP_I3y3z_F3y_a+ABY*I_TWOBODYOVERLAP_H2y3z_F3y_a;
  Double I_TWOBODYOVERLAP_Hy4z_G4y_a = I_TWOBODYOVERLAP_I2y4z_F3y_a+ABY*I_TWOBODYOVERLAP_Hy4z_F3y_a;
  Double I_TWOBODYOVERLAP_H5z_G4y_a = I_TWOBODYOVERLAP_Iy5z_F3y_a+ABY*I_TWOBODYOVERLAP_H5z_F3y_a;
  Double I_TWOBODYOVERLAP_H5x_G3yz_a = I_TWOBODYOVERLAP_I5xz_F3y_a+ABZ*I_TWOBODYOVERLAP_H5x_F3y_a;
  Double I_TWOBODYOVERLAP_H4xy_G3yz_a = I_TWOBODYOVERLAP_I4xyz_F3y_a+ABZ*I_TWOBODYOVERLAP_H4xy_F3y_a;
  Double I_TWOBODYOVERLAP_H4xz_G3yz_a = I_TWOBODYOVERLAP_I4x2z_F3y_a+ABZ*I_TWOBODYOVERLAP_H4xz_F3y_a;
  Double I_TWOBODYOVERLAP_H3x2y_G3yz_a = I_TWOBODYOVERLAP_I3x2yz_F3y_a+ABZ*I_TWOBODYOVERLAP_H3x2y_F3y_a;
  Double I_TWOBODYOVERLAP_H3xyz_G3yz_a = I_TWOBODYOVERLAP_I3xy2z_F3y_a+ABZ*I_TWOBODYOVERLAP_H3xyz_F3y_a;
  Double I_TWOBODYOVERLAP_H3x2z_G3yz_a = I_TWOBODYOVERLAP_I3x3z_F3y_a+ABZ*I_TWOBODYOVERLAP_H3x2z_F3y_a;
  Double I_TWOBODYOVERLAP_H2x3y_G3yz_a = I_TWOBODYOVERLAP_I2x3yz_F3y_a+ABZ*I_TWOBODYOVERLAP_H2x3y_F3y_a;
  Double I_TWOBODYOVERLAP_H2x2yz_G3yz_a = I_TWOBODYOVERLAP_I2x2y2z_F3y_a+ABZ*I_TWOBODYOVERLAP_H2x2yz_F3y_a;
  Double I_TWOBODYOVERLAP_H2xy2z_G3yz_a = I_TWOBODYOVERLAP_I2xy3z_F3y_a+ABZ*I_TWOBODYOVERLAP_H2xy2z_F3y_a;
  Double I_TWOBODYOVERLAP_H2x3z_G3yz_a = I_TWOBODYOVERLAP_I2x4z_F3y_a+ABZ*I_TWOBODYOVERLAP_H2x3z_F3y_a;
  Double I_TWOBODYOVERLAP_Hx4y_G3yz_a = I_TWOBODYOVERLAP_Ix4yz_F3y_a+ABZ*I_TWOBODYOVERLAP_Hx4y_F3y_a;
  Double I_TWOBODYOVERLAP_Hx3yz_G3yz_a = I_TWOBODYOVERLAP_Ix3y2z_F3y_a+ABZ*I_TWOBODYOVERLAP_Hx3yz_F3y_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_G3yz_a = I_TWOBODYOVERLAP_Ix2y3z_F3y_a+ABZ*I_TWOBODYOVERLAP_Hx2y2z_F3y_a;
  Double I_TWOBODYOVERLAP_Hxy3z_G3yz_a = I_TWOBODYOVERLAP_Ixy4z_F3y_a+ABZ*I_TWOBODYOVERLAP_Hxy3z_F3y_a;
  Double I_TWOBODYOVERLAP_Hx4z_G3yz_a = I_TWOBODYOVERLAP_Ix5z_F3y_a+ABZ*I_TWOBODYOVERLAP_Hx4z_F3y_a;
  Double I_TWOBODYOVERLAP_H5y_G3yz_a = I_TWOBODYOVERLAP_I5yz_F3y_a+ABZ*I_TWOBODYOVERLAP_H5y_F3y_a;
  Double I_TWOBODYOVERLAP_H4yz_G3yz_a = I_TWOBODYOVERLAP_I4y2z_F3y_a+ABZ*I_TWOBODYOVERLAP_H4yz_F3y_a;
  Double I_TWOBODYOVERLAP_H3y2z_G3yz_a = I_TWOBODYOVERLAP_I3y3z_F3y_a+ABZ*I_TWOBODYOVERLAP_H3y2z_F3y_a;
  Double I_TWOBODYOVERLAP_H2y3z_G3yz_a = I_TWOBODYOVERLAP_I2y4z_F3y_a+ABZ*I_TWOBODYOVERLAP_H2y3z_F3y_a;
  Double I_TWOBODYOVERLAP_Hy4z_G3yz_a = I_TWOBODYOVERLAP_Iy5z_F3y_a+ABZ*I_TWOBODYOVERLAP_Hy4z_F3y_a;
  Double I_TWOBODYOVERLAP_H5z_G3yz_a = I_TWOBODYOVERLAP_I6z_F3y_a+ABZ*I_TWOBODYOVERLAP_H5z_F3y_a;
  Double I_TWOBODYOVERLAP_H5x_G2y2z_a = I_TWOBODYOVERLAP_I5xz_F2yz_a+ABZ*I_TWOBODYOVERLAP_H5x_F2yz_a;
  Double I_TWOBODYOVERLAP_H4xy_G2y2z_a = I_TWOBODYOVERLAP_I4xyz_F2yz_a+ABZ*I_TWOBODYOVERLAP_H4xy_F2yz_a;
  Double I_TWOBODYOVERLAP_H4xz_G2y2z_a = I_TWOBODYOVERLAP_I4x2z_F2yz_a+ABZ*I_TWOBODYOVERLAP_H4xz_F2yz_a;
  Double I_TWOBODYOVERLAP_H3x2y_G2y2z_a = I_TWOBODYOVERLAP_I3x2yz_F2yz_a+ABZ*I_TWOBODYOVERLAP_H3x2y_F2yz_a;
  Double I_TWOBODYOVERLAP_H3xyz_G2y2z_a = I_TWOBODYOVERLAP_I3xy2z_F2yz_a+ABZ*I_TWOBODYOVERLAP_H3xyz_F2yz_a;
  Double I_TWOBODYOVERLAP_H3x2z_G2y2z_a = I_TWOBODYOVERLAP_I3x3z_F2yz_a+ABZ*I_TWOBODYOVERLAP_H3x2z_F2yz_a;
  Double I_TWOBODYOVERLAP_H2x3y_G2y2z_a = I_TWOBODYOVERLAP_I2x3yz_F2yz_a+ABZ*I_TWOBODYOVERLAP_H2x3y_F2yz_a;
  Double I_TWOBODYOVERLAP_H2x2yz_G2y2z_a = I_TWOBODYOVERLAP_I2x2y2z_F2yz_a+ABZ*I_TWOBODYOVERLAP_H2x2yz_F2yz_a;
  Double I_TWOBODYOVERLAP_H2xy2z_G2y2z_a = I_TWOBODYOVERLAP_I2xy3z_F2yz_a+ABZ*I_TWOBODYOVERLAP_H2xy2z_F2yz_a;
  Double I_TWOBODYOVERLAP_H2x3z_G2y2z_a = I_TWOBODYOVERLAP_I2x4z_F2yz_a+ABZ*I_TWOBODYOVERLAP_H2x3z_F2yz_a;
  Double I_TWOBODYOVERLAP_Hx4y_G2y2z_a = I_TWOBODYOVERLAP_Ix4yz_F2yz_a+ABZ*I_TWOBODYOVERLAP_Hx4y_F2yz_a;
  Double I_TWOBODYOVERLAP_Hx3yz_G2y2z_a = I_TWOBODYOVERLAP_Ix3y2z_F2yz_a+ABZ*I_TWOBODYOVERLAP_Hx3yz_F2yz_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_G2y2z_a = I_TWOBODYOVERLAP_Ix2y3z_F2yz_a+ABZ*I_TWOBODYOVERLAP_Hx2y2z_F2yz_a;
  Double I_TWOBODYOVERLAP_Hxy3z_G2y2z_a = I_TWOBODYOVERLAP_Ixy4z_F2yz_a+ABZ*I_TWOBODYOVERLAP_Hxy3z_F2yz_a;
  Double I_TWOBODYOVERLAP_Hx4z_G2y2z_a = I_TWOBODYOVERLAP_Ix5z_F2yz_a+ABZ*I_TWOBODYOVERLAP_Hx4z_F2yz_a;
  Double I_TWOBODYOVERLAP_H5y_G2y2z_a = I_TWOBODYOVERLAP_I5yz_F2yz_a+ABZ*I_TWOBODYOVERLAP_H5y_F2yz_a;
  Double I_TWOBODYOVERLAP_H4yz_G2y2z_a = I_TWOBODYOVERLAP_I4y2z_F2yz_a+ABZ*I_TWOBODYOVERLAP_H4yz_F2yz_a;
  Double I_TWOBODYOVERLAP_H3y2z_G2y2z_a = I_TWOBODYOVERLAP_I3y3z_F2yz_a+ABZ*I_TWOBODYOVERLAP_H3y2z_F2yz_a;
  Double I_TWOBODYOVERLAP_H2y3z_G2y2z_a = I_TWOBODYOVERLAP_I2y4z_F2yz_a+ABZ*I_TWOBODYOVERLAP_H2y3z_F2yz_a;
  Double I_TWOBODYOVERLAP_Hy4z_G2y2z_a = I_TWOBODYOVERLAP_Iy5z_F2yz_a+ABZ*I_TWOBODYOVERLAP_Hy4z_F2yz_a;
  Double I_TWOBODYOVERLAP_H5z_G2y2z_a = I_TWOBODYOVERLAP_I6z_F2yz_a+ABZ*I_TWOBODYOVERLAP_H5z_F2yz_a;
  Double I_TWOBODYOVERLAP_H5x_Gy3z_a = I_TWOBODYOVERLAP_I5xy_F3z_a+ABY*I_TWOBODYOVERLAP_H5x_F3z_a;
  Double I_TWOBODYOVERLAP_H4xy_Gy3z_a = I_TWOBODYOVERLAP_I4x2y_F3z_a+ABY*I_TWOBODYOVERLAP_H4xy_F3z_a;
  Double I_TWOBODYOVERLAP_H4xz_Gy3z_a = I_TWOBODYOVERLAP_I4xyz_F3z_a+ABY*I_TWOBODYOVERLAP_H4xz_F3z_a;
  Double I_TWOBODYOVERLAP_H3x2y_Gy3z_a = I_TWOBODYOVERLAP_I3x3y_F3z_a+ABY*I_TWOBODYOVERLAP_H3x2y_F3z_a;
  Double I_TWOBODYOVERLAP_H3xyz_Gy3z_a = I_TWOBODYOVERLAP_I3x2yz_F3z_a+ABY*I_TWOBODYOVERLAP_H3xyz_F3z_a;
  Double I_TWOBODYOVERLAP_H3x2z_Gy3z_a = I_TWOBODYOVERLAP_I3xy2z_F3z_a+ABY*I_TWOBODYOVERLAP_H3x2z_F3z_a;
  Double I_TWOBODYOVERLAP_H2x3y_Gy3z_a = I_TWOBODYOVERLAP_I2x4y_F3z_a+ABY*I_TWOBODYOVERLAP_H2x3y_F3z_a;
  Double I_TWOBODYOVERLAP_H2x2yz_Gy3z_a = I_TWOBODYOVERLAP_I2x3yz_F3z_a+ABY*I_TWOBODYOVERLAP_H2x2yz_F3z_a;
  Double I_TWOBODYOVERLAP_H2xy2z_Gy3z_a = I_TWOBODYOVERLAP_I2x2y2z_F3z_a+ABY*I_TWOBODYOVERLAP_H2xy2z_F3z_a;
  Double I_TWOBODYOVERLAP_H2x3z_Gy3z_a = I_TWOBODYOVERLAP_I2xy3z_F3z_a+ABY*I_TWOBODYOVERLAP_H2x3z_F3z_a;
  Double I_TWOBODYOVERLAP_Hx4y_Gy3z_a = I_TWOBODYOVERLAP_Ix5y_F3z_a+ABY*I_TWOBODYOVERLAP_Hx4y_F3z_a;
  Double I_TWOBODYOVERLAP_Hx3yz_Gy3z_a = I_TWOBODYOVERLAP_Ix4yz_F3z_a+ABY*I_TWOBODYOVERLAP_Hx3yz_F3z_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_Gy3z_a = I_TWOBODYOVERLAP_Ix3y2z_F3z_a+ABY*I_TWOBODYOVERLAP_Hx2y2z_F3z_a;
  Double I_TWOBODYOVERLAP_Hxy3z_Gy3z_a = I_TWOBODYOVERLAP_Ix2y3z_F3z_a+ABY*I_TWOBODYOVERLAP_Hxy3z_F3z_a;
  Double I_TWOBODYOVERLAP_Hx4z_Gy3z_a = I_TWOBODYOVERLAP_Ixy4z_F3z_a+ABY*I_TWOBODYOVERLAP_Hx4z_F3z_a;
  Double I_TWOBODYOVERLAP_H5y_Gy3z_a = I_TWOBODYOVERLAP_I6y_F3z_a+ABY*I_TWOBODYOVERLAP_H5y_F3z_a;
  Double I_TWOBODYOVERLAP_H4yz_Gy3z_a = I_TWOBODYOVERLAP_I5yz_F3z_a+ABY*I_TWOBODYOVERLAP_H4yz_F3z_a;
  Double I_TWOBODYOVERLAP_H3y2z_Gy3z_a = I_TWOBODYOVERLAP_I4y2z_F3z_a+ABY*I_TWOBODYOVERLAP_H3y2z_F3z_a;
  Double I_TWOBODYOVERLAP_H2y3z_Gy3z_a = I_TWOBODYOVERLAP_I3y3z_F3z_a+ABY*I_TWOBODYOVERLAP_H2y3z_F3z_a;
  Double I_TWOBODYOVERLAP_Hy4z_Gy3z_a = I_TWOBODYOVERLAP_I2y4z_F3z_a+ABY*I_TWOBODYOVERLAP_Hy4z_F3z_a;
  Double I_TWOBODYOVERLAP_H5z_Gy3z_a = I_TWOBODYOVERLAP_Iy5z_F3z_a+ABY*I_TWOBODYOVERLAP_H5z_F3z_a;
  Double I_TWOBODYOVERLAP_H5x_G4z_a = I_TWOBODYOVERLAP_I5xz_F3z_a+ABZ*I_TWOBODYOVERLAP_H5x_F3z_a;
  Double I_TWOBODYOVERLAP_H4xy_G4z_a = I_TWOBODYOVERLAP_I4xyz_F3z_a+ABZ*I_TWOBODYOVERLAP_H4xy_F3z_a;
  Double I_TWOBODYOVERLAP_H4xz_G4z_a = I_TWOBODYOVERLAP_I4x2z_F3z_a+ABZ*I_TWOBODYOVERLAP_H4xz_F3z_a;
  Double I_TWOBODYOVERLAP_H3x2y_G4z_a = I_TWOBODYOVERLAP_I3x2yz_F3z_a+ABZ*I_TWOBODYOVERLAP_H3x2y_F3z_a;
  Double I_TWOBODYOVERLAP_H3xyz_G4z_a = I_TWOBODYOVERLAP_I3xy2z_F3z_a+ABZ*I_TWOBODYOVERLAP_H3xyz_F3z_a;
  Double I_TWOBODYOVERLAP_H3x2z_G4z_a = I_TWOBODYOVERLAP_I3x3z_F3z_a+ABZ*I_TWOBODYOVERLAP_H3x2z_F3z_a;
  Double I_TWOBODYOVERLAP_H2x3y_G4z_a = I_TWOBODYOVERLAP_I2x3yz_F3z_a+ABZ*I_TWOBODYOVERLAP_H2x3y_F3z_a;
  Double I_TWOBODYOVERLAP_H2x2yz_G4z_a = I_TWOBODYOVERLAP_I2x2y2z_F3z_a+ABZ*I_TWOBODYOVERLAP_H2x2yz_F3z_a;
  Double I_TWOBODYOVERLAP_H2xy2z_G4z_a = I_TWOBODYOVERLAP_I2xy3z_F3z_a+ABZ*I_TWOBODYOVERLAP_H2xy2z_F3z_a;
  Double I_TWOBODYOVERLAP_H2x3z_G4z_a = I_TWOBODYOVERLAP_I2x4z_F3z_a+ABZ*I_TWOBODYOVERLAP_H2x3z_F3z_a;
  Double I_TWOBODYOVERLAP_Hx4y_G4z_a = I_TWOBODYOVERLAP_Ix4yz_F3z_a+ABZ*I_TWOBODYOVERLAP_Hx4y_F3z_a;
  Double I_TWOBODYOVERLAP_Hx3yz_G4z_a = I_TWOBODYOVERLAP_Ix3y2z_F3z_a+ABZ*I_TWOBODYOVERLAP_Hx3yz_F3z_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_G4z_a = I_TWOBODYOVERLAP_Ix2y3z_F3z_a+ABZ*I_TWOBODYOVERLAP_Hx2y2z_F3z_a;
  Double I_TWOBODYOVERLAP_Hxy3z_G4z_a = I_TWOBODYOVERLAP_Ixy4z_F3z_a+ABZ*I_TWOBODYOVERLAP_Hxy3z_F3z_a;
  Double I_TWOBODYOVERLAP_Hx4z_G4z_a = I_TWOBODYOVERLAP_Ix5z_F3z_a+ABZ*I_TWOBODYOVERLAP_Hx4z_F3z_a;
  Double I_TWOBODYOVERLAP_H5y_G4z_a = I_TWOBODYOVERLAP_I5yz_F3z_a+ABZ*I_TWOBODYOVERLAP_H5y_F3z_a;
  Double I_TWOBODYOVERLAP_H4yz_G4z_a = I_TWOBODYOVERLAP_I4y2z_F3z_a+ABZ*I_TWOBODYOVERLAP_H4yz_F3z_a;
  Double I_TWOBODYOVERLAP_H3y2z_G4z_a = I_TWOBODYOVERLAP_I3y3z_F3z_a+ABZ*I_TWOBODYOVERLAP_H3y2z_F3z_a;
  Double I_TWOBODYOVERLAP_H2y3z_G4z_a = I_TWOBODYOVERLAP_I2y4z_F3z_a+ABZ*I_TWOBODYOVERLAP_H2y3z_F3z_a;
  Double I_TWOBODYOVERLAP_Hy4z_G4z_a = I_TWOBODYOVERLAP_Iy5z_F3z_a+ABZ*I_TWOBODYOVERLAP_Hy4z_F3z_a;
  Double I_TWOBODYOVERLAP_H5z_G4z_a = I_TWOBODYOVERLAP_I6z_F3z_a+ABZ*I_TWOBODYOVERLAP_H5z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_G_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_G_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_G
   ************************************************************/
  abcd[0] = 2.0E0*I_TWOBODYOVERLAP_H5x_G4x_a-4*I_TWOBODYOVERLAP_F3x_G4x;
  abcd[1] = 2.0E0*I_TWOBODYOVERLAP_H4xy_G4x_a-3*I_TWOBODYOVERLAP_F2xy_G4x;
  abcd[2] = 2.0E0*I_TWOBODYOVERLAP_H4xz_G4x_a-3*I_TWOBODYOVERLAP_F2xz_G4x;
  abcd[3] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_G4x_a-2*I_TWOBODYOVERLAP_Fx2y_G4x;
  abcd[4] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G4x_a-2*I_TWOBODYOVERLAP_Fxyz_G4x;
  abcd[5] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_G4x_a-2*I_TWOBODYOVERLAP_Fx2z_G4x;
  abcd[6] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_G4x_a-1*I_TWOBODYOVERLAP_F3y_G4x;
  abcd[7] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G4x_a-1*I_TWOBODYOVERLAP_F2yz_G4x;
  abcd[8] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G4x_a-1*I_TWOBODYOVERLAP_Fy2z_G4x;
  abcd[9] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_G4x_a-1*I_TWOBODYOVERLAP_F3z_G4x;
  abcd[10] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_G4x_a;
  abcd[11] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G4x_a;
  abcd[12] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G4x_a;
  abcd[13] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G4x_a;
  abcd[14] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_G4x_a;
  abcd[15] = 2.0E0*I_TWOBODYOVERLAP_H5x_G3xy_a-4*I_TWOBODYOVERLAP_F3x_G3xy;
  abcd[16] = 2.0E0*I_TWOBODYOVERLAP_H4xy_G3xy_a-3*I_TWOBODYOVERLAP_F2xy_G3xy;
  abcd[17] = 2.0E0*I_TWOBODYOVERLAP_H4xz_G3xy_a-3*I_TWOBODYOVERLAP_F2xz_G3xy;
  abcd[18] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_G3xy_a-2*I_TWOBODYOVERLAP_Fx2y_G3xy;
  abcd[19] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G3xy_a-2*I_TWOBODYOVERLAP_Fxyz_G3xy;
  abcd[20] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_G3xy_a-2*I_TWOBODYOVERLAP_Fx2z_G3xy;
  abcd[21] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_G3xy_a-1*I_TWOBODYOVERLAP_F3y_G3xy;
  abcd[22] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G3xy_a-1*I_TWOBODYOVERLAP_F2yz_G3xy;
  abcd[23] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G3xy_a-1*I_TWOBODYOVERLAP_Fy2z_G3xy;
  abcd[24] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_G3xy_a-1*I_TWOBODYOVERLAP_F3z_G3xy;
  abcd[25] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_G3xy_a;
  abcd[26] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G3xy_a;
  abcd[27] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G3xy_a;
  abcd[28] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G3xy_a;
  abcd[29] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_G3xy_a;
  abcd[30] = 2.0E0*I_TWOBODYOVERLAP_H5x_G3xz_a-4*I_TWOBODYOVERLAP_F3x_G3xz;
  abcd[31] = 2.0E0*I_TWOBODYOVERLAP_H4xy_G3xz_a-3*I_TWOBODYOVERLAP_F2xy_G3xz;
  abcd[32] = 2.0E0*I_TWOBODYOVERLAP_H4xz_G3xz_a-3*I_TWOBODYOVERLAP_F2xz_G3xz;
  abcd[33] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_G3xz_a-2*I_TWOBODYOVERLAP_Fx2y_G3xz;
  abcd[34] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G3xz_a-2*I_TWOBODYOVERLAP_Fxyz_G3xz;
  abcd[35] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_G3xz_a-2*I_TWOBODYOVERLAP_Fx2z_G3xz;
  abcd[36] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_G3xz_a-1*I_TWOBODYOVERLAP_F3y_G3xz;
  abcd[37] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G3xz_a-1*I_TWOBODYOVERLAP_F2yz_G3xz;
  abcd[38] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G3xz_a-1*I_TWOBODYOVERLAP_Fy2z_G3xz;
  abcd[39] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_G3xz_a-1*I_TWOBODYOVERLAP_F3z_G3xz;
  abcd[40] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_G3xz_a;
  abcd[41] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G3xz_a;
  abcd[42] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G3xz_a;
  abcd[43] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G3xz_a;
  abcd[44] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_G3xz_a;
  abcd[45] = 2.0E0*I_TWOBODYOVERLAP_H5x_G2x2y_a-4*I_TWOBODYOVERLAP_F3x_G2x2y;
  abcd[46] = 2.0E0*I_TWOBODYOVERLAP_H4xy_G2x2y_a-3*I_TWOBODYOVERLAP_F2xy_G2x2y;
  abcd[47] = 2.0E0*I_TWOBODYOVERLAP_H4xz_G2x2y_a-3*I_TWOBODYOVERLAP_F2xz_G2x2y;
  abcd[48] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_G2x2y_a-2*I_TWOBODYOVERLAP_Fx2y_G2x2y;
  abcd[49] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G2x2y_a-2*I_TWOBODYOVERLAP_Fxyz_G2x2y;
  abcd[50] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_G2x2y_a-2*I_TWOBODYOVERLAP_Fx2z_G2x2y;
  abcd[51] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_G2x2y_a-1*I_TWOBODYOVERLAP_F3y_G2x2y;
  abcd[52] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G2x2y_a-1*I_TWOBODYOVERLAP_F2yz_G2x2y;
  abcd[53] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G2x2y_a-1*I_TWOBODYOVERLAP_Fy2z_G2x2y;
  abcd[54] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_G2x2y_a-1*I_TWOBODYOVERLAP_F3z_G2x2y;
  abcd[55] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_G2x2y_a;
  abcd[56] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G2x2y_a;
  abcd[57] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G2x2y_a;
  abcd[58] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G2x2y_a;
  abcd[59] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_G2x2y_a;
  abcd[60] = 2.0E0*I_TWOBODYOVERLAP_H5x_G2xyz_a-4*I_TWOBODYOVERLAP_F3x_G2xyz;
  abcd[61] = 2.0E0*I_TWOBODYOVERLAP_H4xy_G2xyz_a-3*I_TWOBODYOVERLAP_F2xy_G2xyz;
  abcd[62] = 2.0E0*I_TWOBODYOVERLAP_H4xz_G2xyz_a-3*I_TWOBODYOVERLAP_F2xz_G2xyz;
  abcd[63] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_G2xyz_a-2*I_TWOBODYOVERLAP_Fx2y_G2xyz;
  abcd[64] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G2xyz_a-2*I_TWOBODYOVERLAP_Fxyz_G2xyz;
  abcd[65] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_G2xyz_a-2*I_TWOBODYOVERLAP_Fx2z_G2xyz;
  abcd[66] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_G2xyz_a-1*I_TWOBODYOVERLAP_F3y_G2xyz;
  abcd[67] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G2xyz_a-1*I_TWOBODYOVERLAP_F2yz_G2xyz;
  abcd[68] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G2xyz_a-1*I_TWOBODYOVERLAP_Fy2z_G2xyz;
  abcd[69] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_G2xyz_a-1*I_TWOBODYOVERLAP_F3z_G2xyz;
  abcd[70] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_G2xyz_a;
  abcd[71] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G2xyz_a;
  abcd[72] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G2xyz_a;
  abcd[73] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G2xyz_a;
  abcd[74] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_G2xyz_a;
  abcd[75] = 2.0E0*I_TWOBODYOVERLAP_H5x_G2x2z_a-4*I_TWOBODYOVERLAP_F3x_G2x2z;
  abcd[76] = 2.0E0*I_TWOBODYOVERLAP_H4xy_G2x2z_a-3*I_TWOBODYOVERLAP_F2xy_G2x2z;
  abcd[77] = 2.0E0*I_TWOBODYOVERLAP_H4xz_G2x2z_a-3*I_TWOBODYOVERLAP_F2xz_G2x2z;
  abcd[78] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_G2x2z_a-2*I_TWOBODYOVERLAP_Fx2y_G2x2z;
  abcd[79] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G2x2z_a-2*I_TWOBODYOVERLAP_Fxyz_G2x2z;
  abcd[80] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_G2x2z_a-2*I_TWOBODYOVERLAP_Fx2z_G2x2z;
  abcd[81] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_G2x2z_a-1*I_TWOBODYOVERLAP_F3y_G2x2z;
  abcd[82] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G2x2z_a-1*I_TWOBODYOVERLAP_F2yz_G2x2z;
  abcd[83] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G2x2z_a-1*I_TWOBODYOVERLAP_Fy2z_G2x2z;
  abcd[84] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_G2x2z_a-1*I_TWOBODYOVERLAP_F3z_G2x2z;
  abcd[85] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_G2x2z_a;
  abcd[86] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G2x2z_a;
  abcd[87] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G2x2z_a;
  abcd[88] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G2x2z_a;
  abcd[89] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_G2x2z_a;
  abcd[90] = 2.0E0*I_TWOBODYOVERLAP_H5x_Gx3y_a-4*I_TWOBODYOVERLAP_F3x_Gx3y;
  abcd[91] = 2.0E0*I_TWOBODYOVERLAP_H4xy_Gx3y_a-3*I_TWOBODYOVERLAP_F2xy_Gx3y;
  abcd[92] = 2.0E0*I_TWOBODYOVERLAP_H4xz_Gx3y_a-3*I_TWOBODYOVERLAP_F2xz_Gx3y;
  abcd[93] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_Gx3y_a-2*I_TWOBODYOVERLAP_Fx2y_Gx3y;
  abcd[94] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_Gx3y_a-2*I_TWOBODYOVERLAP_Fxyz_Gx3y;
  abcd[95] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_Gx3y_a-2*I_TWOBODYOVERLAP_Fx2z_Gx3y;
  abcd[96] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_Gx3y_a-1*I_TWOBODYOVERLAP_F3y_Gx3y;
  abcd[97] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_Gx3y_a-1*I_TWOBODYOVERLAP_F2yz_Gx3y;
  abcd[98] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_Gx3y_a-1*I_TWOBODYOVERLAP_Fy2z_Gx3y;
  abcd[99] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_Gx3y_a-1*I_TWOBODYOVERLAP_F3z_Gx3y;
  abcd[100] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_Gx3y_a;
  abcd[101] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_Gx3y_a;
  abcd[102] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_Gx3y_a;
  abcd[103] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_Gx3y_a;
  abcd[104] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_Gx3y_a;
  abcd[105] = 2.0E0*I_TWOBODYOVERLAP_H5x_Gx2yz_a-4*I_TWOBODYOVERLAP_F3x_Gx2yz;
  abcd[106] = 2.0E0*I_TWOBODYOVERLAP_H4xy_Gx2yz_a-3*I_TWOBODYOVERLAP_F2xy_Gx2yz;
  abcd[107] = 2.0E0*I_TWOBODYOVERLAP_H4xz_Gx2yz_a-3*I_TWOBODYOVERLAP_F2xz_Gx2yz;
  abcd[108] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_Gx2yz_a-2*I_TWOBODYOVERLAP_Fx2y_Gx2yz;
  abcd[109] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_Gx2yz_a-2*I_TWOBODYOVERLAP_Fxyz_Gx2yz;
  abcd[110] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_Gx2yz_a-2*I_TWOBODYOVERLAP_Fx2z_Gx2yz;
  abcd[111] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_Gx2yz_a-1*I_TWOBODYOVERLAP_F3y_Gx2yz;
  abcd[112] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_Gx2yz_a-1*I_TWOBODYOVERLAP_F2yz_Gx2yz;
  abcd[113] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_Gx2yz_a-1*I_TWOBODYOVERLAP_Fy2z_Gx2yz;
  abcd[114] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_Gx2yz_a-1*I_TWOBODYOVERLAP_F3z_Gx2yz;
  abcd[115] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_Gx2yz_a;
  abcd[116] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_Gx2yz_a;
  abcd[117] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_Gx2yz_a;
  abcd[118] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_Gx2yz_a;
  abcd[119] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_Gx2yz_a;
  abcd[120] = 2.0E0*I_TWOBODYOVERLAP_H5x_Gxy2z_a-4*I_TWOBODYOVERLAP_F3x_Gxy2z;
  abcd[121] = 2.0E0*I_TWOBODYOVERLAP_H4xy_Gxy2z_a-3*I_TWOBODYOVERLAP_F2xy_Gxy2z;
  abcd[122] = 2.0E0*I_TWOBODYOVERLAP_H4xz_Gxy2z_a-3*I_TWOBODYOVERLAP_F2xz_Gxy2z;
  abcd[123] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_Gxy2z_a-2*I_TWOBODYOVERLAP_Fx2y_Gxy2z;
  abcd[124] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_Gxy2z_a-2*I_TWOBODYOVERLAP_Fxyz_Gxy2z;
  abcd[125] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_Gxy2z_a-2*I_TWOBODYOVERLAP_Fx2z_Gxy2z;
  abcd[126] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_Gxy2z_a-1*I_TWOBODYOVERLAP_F3y_Gxy2z;
  abcd[127] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_Gxy2z_a-1*I_TWOBODYOVERLAP_F2yz_Gxy2z;
  abcd[128] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_Gxy2z_a-1*I_TWOBODYOVERLAP_Fy2z_Gxy2z;
  abcd[129] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_Gxy2z_a-1*I_TWOBODYOVERLAP_F3z_Gxy2z;
  abcd[130] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_Gxy2z_a;
  abcd[131] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_Gxy2z_a;
  abcd[132] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_Gxy2z_a;
  abcd[133] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_Gxy2z_a;
  abcd[134] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_Gxy2z_a;
  abcd[135] = 2.0E0*I_TWOBODYOVERLAP_H5x_Gx3z_a-4*I_TWOBODYOVERLAP_F3x_Gx3z;
  abcd[136] = 2.0E0*I_TWOBODYOVERLAP_H4xy_Gx3z_a-3*I_TWOBODYOVERLAP_F2xy_Gx3z;
  abcd[137] = 2.0E0*I_TWOBODYOVERLAP_H4xz_Gx3z_a-3*I_TWOBODYOVERLAP_F2xz_Gx3z;
  abcd[138] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_Gx3z_a-2*I_TWOBODYOVERLAP_Fx2y_Gx3z;
  abcd[139] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_Gx3z_a-2*I_TWOBODYOVERLAP_Fxyz_Gx3z;
  abcd[140] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_Gx3z_a-2*I_TWOBODYOVERLAP_Fx2z_Gx3z;
  abcd[141] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_Gx3z_a-1*I_TWOBODYOVERLAP_F3y_Gx3z;
  abcd[142] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_Gx3z_a-1*I_TWOBODYOVERLAP_F2yz_Gx3z;
  abcd[143] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_Gx3z_a-1*I_TWOBODYOVERLAP_Fy2z_Gx3z;
  abcd[144] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_Gx3z_a-1*I_TWOBODYOVERLAP_F3z_Gx3z;
  abcd[145] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_Gx3z_a;
  abcd[146] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_Gx3z_a;
  abcd[147] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_Gx3z_a;
  abcd[148] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_Gx3z_a;
  abcd[149] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_Gx3z_a;
  abcd[150] = 2.0E0*I_TWOBODYOVERLAP_H5x_G4y_a-4*I_TWOBODYOVERLAP_F3x_G4y;
  abcd[151] = 2.0E0*I_TWOBODYOVERLAP_H4xy_G4y_a-3*I_TWOBODYOVERLAP_F2xy_G4y;
  abcd[152] = 2.0E0*I_TWOBODYOVERLAP_H4xz_G4y_a-3*I_TWOBODYOVERLAP_F2xz_G4y;
  abcd[153] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_G4y_a-2*I_TWOBODYOVERLAP_Fx2y_G4y;
  abcd[154] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G4y_a-2*I_TWOBODYOVERLAP_Fxyz_G4y;
  abcd[155] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_G4y_a-2*I_TWOBODYOVERLAP_Fx2z_G4y;
  abcd[156] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_G4y_a-1*I_TWOBODYOVERLAP_F3y_G4y;
  abcd[157] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G4y_a-1*I_TWOBODYOVERLAP_F2yz_G4y;
  abcd[158] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G4y_a-1*I_TWOBODYOVERLAP_Fy2z_G4y;
  abcd[159] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_G4y_a-1*I_TWOBODYOVERLAP_F3z_G4y;
  abcd[160] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_G4y_a;
  abcd[161] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G4y_a;
  abcd[162] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G4y_a;
  abcd[163] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G4y_a;
  abcd[164] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_G4y_a;
  abcd[165] = 2.0E0*I_TWOBODYOVERLAP_H5x_G3yz_a-4*I_TWOBODYOVERLAP_F3x_G3yz;
  abcd[166] = 2.0E0*I_TWOBODYOVERLAP_H4xy_G3yz_a-3*I_TWOBODYOVERLAP_F2xy_G3yz;
  abcd[167] = 2.0E0*I_TWOBODYOVERLAP_H4xz_G3yz_a-3*I_TWOBODYOVERLAP_F2xz_G3yz;
  abcd[168] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_G3yz_a-2*I_TWOBODYOVERLAP_Fx2y_G3yz;
  abcd[169] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G3yz_a-2*I_TWOBODYOVERLAP_Fxyz_G3yz;
  abcd[170] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_G3yz_a-2*I_TWOBODYOVERLAP_Fx2z_G3yz;
  abcd[171] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_G3yz_a-1*I_TWOBODYOVERLAP_F3y_G3yz;
  abcd[172] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G3yz_a-1*I_TWOBODYOVERLAP_F2yz_G3yz;
  abcd[173] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G3yz_a-1*I_TWOBODYOVERLAP_Fy2z_G3yz;
  abcd[174] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_G3yz_a-1*I_TWOBODYOVERLAP_F3z_G3yz;
  abcd[175] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_G3yz_a;
  abcd[176] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G3yz_a;
  abcd[177] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G3yz_a;
  abcd[178] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G3yz_a;
  abcd[179] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_G3yz_a;
  abcd[180] = 2.0E0*I_TWOBODYOVERLAP_H5x_G2y2z_a-4*I_TWOBODYOVERLAP_F3x_G2y2z;
  abcd[181] = 2.0E0*I_TWOBODYOVERLAP_H4xy_G2y2z_a-3*I_TWOBODYOVERLAP_F2xy_G2y2z;
  abcd[182] = 2.0E0*I_TWOBODYOVERLAP_H4xz_G2y2z_a-3*I_TWOBODYOVERLAP_F2xz_G2y2z;
  abcd[183] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_G2y2z_a-2*I_TWOBODYOVERLAP_Fx2y_G2y2z;
  abcd[184] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G2y2z_a-2*I_TWOBODYOVERLAP_Fxyz_G2y2z;
  abcd[185] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_G2y2z_a-2*I_TWOBODYOVERLAP_Fx2z_G2y2z;
  abcd[186] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_G2y2z_a-1*I_TWOBODYOVERLAP_F3y_G2y2z;
  abcd[187] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G2y2z_a-1*I_TWOBODYOVERLAP_F2yz_G2y2z;
  abcd[188] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G2y2z_a-1*I_TWOBODYOVERLAP_Fy2z_G2y2z;
  abcd[189] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_G2y2z_a-1*I_TWOBODYOVERLAP_F3z_G2y2z;
  abcd[190] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_G2y2z_a;
  abcd[191] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G2y2z_a;
  abcd[192] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G2y2z_a;
  abcd[193] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G2y2z_a;
  abcd[194] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_G2y2z_a;
  abcd[195] = 2.0E0*I_TWOBODYOVERLAP_H5x_Gy3z_a-4*I_TWOBODYOVERLAP_F3x_Gy3z;
  abcd[196] = 2.0E0*I_TWOBODYOVERLAP_H4xy_Gy3z_a-3*I_TWOBODYOVERLAP_F2xy_Gy3z;
  abcd[197] = 2.0E0*I_TWOBODYOVERLAP_H4xz_Gy3z_a-3*I_TWOBODYOVERLAP_F2xz_Gy3z;
  abcd[198] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_Gy3z_a-2*I_TWOBODYOVERLAP_Fx2y_Gy3z;
  abcd[199] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_Gy3z_a-2*I_TWOBODYOVERLAP_Fxyz_Gy3z;
  abcd[200] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_Gy3z_a-2*I_TWOBODYOVERLAP_Fx2z_Gy3z;
  abcd[201] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_Gy3z_a-1*I_TWOBODYOVERLAP_F3y_Gy3z;
  abcd[202] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_Gy3z_a-1*I_TWOBODYOVERLAP_F2yz_Gy3z;
  abcd[203] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_Gy3z_a-1*I_TWOBODYOVERLAP_Fy2z_Gy3z;
  abcd[204] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_Gy3z_a-1*I_TWOBODYOVERLAP_F3z_Gy3z;
  abcd[205] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_Gy3z_a;
  abcd[206] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_Gy3z_a;
  abcd[207] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_Gy3z_a;
  abcd[208] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_Gy3z_a;
  abcd[209] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_Gy3z_a;
  abcd[210] = 2.0E0*I_TWOBODYOVERLAP_H5x_G4z_a-4*I_TWOBODYOVERLAP_F3x_G4z;
  abcd[211] = 2.0E0*I_TWOBODYOVERLAP_H4xy_G4z_a-3*I_TWOBODYOVERLAP_F2xy_G4z;
  abcd[212] = 2.0E0*I_TWOBODYOVERLAP_H4xz_G4z_a-3*I_TWOBODYOVERLAP_F2xz_G4z;
  abcd[213] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_G4z_a-2*I_TWOBODYOVERLAP_Fx2y_G4z;
  abcd[214] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G4z_a-2*I_TWOBODYOVERLAP_Fxyz_G4z;
  abcd[215] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_G4z_a-2*I_TWOBODYOVERLAP_Fx2z_G4z;
  abcd[216] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_G4z_a-1*I_TWOBODYOVERLAP_F3y_G4z;
  abcd[217] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G4z_a-1*I_TWOBODYOVERLAP_F2yz_G4z;
  abcd[218] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G4z_a-1*I_TWOBODYOVERLAP_Fy2z_G4z;
  abcd[219] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_G4z_a-1*I_TWOBODYOVERLAP_F3z_G4z;
  abcd[220] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_G4z_a;
  abcd[221] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G4z_a;
  abcd[222] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G4z_a;
  abcd[223] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G4z_a;
  abcd[224] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_G4z_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_G_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_G_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_G
   ************************************************************/
  abcd[225] = 2.0E0*I_TWOBODYOVERLAP_H4xy_G4x_a;
  abcd[226] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_G4x_a-1*I_TWOBODYOVERLAP_F3x_G4x;
  abcd[227] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G4x_a;
  abcd[228] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_G4x_a-2*I_TWOBODYOVERLAP_F2xy_G4x;
  abcd[229] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G4x_a-1*I_TWOBODYOVERLAP_F2xz_G4x;
  abcd[230] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G4x_a;
  abcd[231] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_G4x_a-3*I_TWOBODYOVERLAP_Fx2y_G4x;
  abcd[232] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G4x_a-2*I_TWOBODYOVERLAP_Fxyz_G4x;
  abcd[233] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G4x_a-1*I_TWOBODYOVERLAP_Fx2z_G4x;
  abcd[234] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G4x_a;
  abcd[235] = 2.0E0*I_TWOBODYOVERLAP_H5y_G4x_a-4*I_TWOBODYOVERLAP_F3y_G4x;
  abcd[236] = 2.0E0*I_TWOBODYOVERLAP_H4yz_G4x_a-3*I_TWOBODYOVERLAP_F2yz_G4x;
  abcd[237] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_G4x_a-2*I_TWOBODYOVERLAP_Fy2z_G4x;
  abcd[238] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_G4x_a-1*I_TWOBODYOVERLAP_F3z_G4x;
  abcd[239] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_G4x_a;
  abcd[240] = 2.0E0*I_TWOBODYOVERLAP_H4xy_G3xy_a;
  abcd[241] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_G3xy_a-1*I_TWOBODYOVERLAP_F3x_G3xy;
  abcd[242] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G3xy_a;
  abcd[243] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_G3xy_a-2*I_TWOBODYOVERLAP_F2xy_G3xy;
  abcd[244] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G3xy_a-1*I_TWOBODYOVERLAP_F2xz_G3xy;
  abcd[245] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G3xy_a;
  abcd[246] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_G3xy_a-3*I_TWOBODYOVERLAP_Fx2y_G3xy;
  abcd[247] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G3xy_a-2*I_TWOBODYOVERLAP_Fxyz_G3xy;
  abcd[248] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G3xy_a-1*I_TWOBODYOVERLAP_Fx2z_G3xy;
  abcd[249] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G3xy_a;
  abcd[250] = 2.0E0*I_TWOBODYOVERLAP_H5y_G3xy_a-4*I_TWOBODYOVERLAP_F3y_G3xy;
  abcd[251] = 2.0E0*I_TWOBODYOVERLAP_H4yz_G3xy_a-3*I_TWOBODYOVERLAP_F2yz_G3xy;
  abcd[252] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_G3xy_a-2*I_TWOBODYOVERLAP_Fy2z_G3xy;
  abcd[253] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_G3xy_a-1*I_TWOBODYOVERLAP_F3z_G3xy;
  abcd[254] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_G3xy_a;
  abcd[255] = 2.0E0*I_TWOBODYOVERLAP_H4xy_G3xz_a;
  abcd[256] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_G3xz_a-1*I_TWOBODYOVERLAP_F3x_G3xz;
  abcd[257] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G3xz_a;
  abcd[258] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_G3xz_a-2*I_TWOBODYOVERLAP_F2xy_G3xz;
  abcd[259] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G3xz_a-1*I_TWOBODYOVERLAP_F2xz_G3xz;
  abcd[260] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G3xz_a;
  abcd[261] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_G3xz_a-3*I_TWOBODYOVERLAP_Fx2y_G3xz;
  abcd[262] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G3xz_a-2*I_TWOBODYOVERLAP_Fxyz_G3xz;
  abcd[263] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G3xz_a-1*I_TWOBODYOVERLAP_Fx2z_G3xz;
  abcd[264] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G3xz_a;
  abcd[265] = 2.0E0*I_TWOBODYOVERLAP_H5y_G3xz_a-4*I_TWOBODYOVERLAP_F3y_G3xz;
  abcd[266] = 2.0E0*I_TWOBODYOVERLAP_H4yz_G3xz_a-3*I_TWOBODYOVERLAP_F2yz_G3xz;
  abcd[267] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_G3xz_a-2*I_TWOBODYOVERLAP_Fy2z_G3xz;
  abcd[268] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_G3xz_a-1*I_TWOBODYOVERLAP_F3z_G3xz;
  abcd[269] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_G3xz_a;
  abcd[270] = 2.0E0*I_TWOBODYOVERLAP_H4xy_G2x2y_a;
  abcd[271] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_G2x2y_a-1*I_TWOBODYOVERLAP_F3x_G2x2y;
  abcd[272] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G2x2y_a;
  abcd[273] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_G2x2y_a-2*I_TWOBODYOVERLAP_F2xy_G2x2y;
  abcd[274] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G2x2y_a-1*I_TWOBODYOVERLAP_F2xz_G2x2y;
  abcd[275] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G2x2y_a;
  abcd[276] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_G2x2y_a-3*I_TWOBODYOVERLAP_Fx2y_G2x2y;
  abcd[277] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G2x2y_a-2*I_TWOBODYOVERLAP_Fxyz_G2x2y;
  abcd[278] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G2x2y_a-1*I_TWOBODYOVERLAP_Fx2z_G2x2y;
  abcd[279] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G2x2y_a;
  abcd[280] = 2.0E0*I_TWOBODYOVERLAP_H5y_G2x2y_a-4*I_TWOBODYOVERLAP_F3y_G2x2y;
  abcd[281] = 2.0E0*I_TWOBODYOVERLAP_H4yz_G2x2y_a-3*I_TWOBODYOVERLAP_F2yz_G2x2y;
  abcd[282] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_G2x2y_a-2*I_TWOBODYOVERLAP_Fy2z_G2x2y;
  abcd[283] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_G2x2y_a-1*I_TWOBODYOVERLAP_F3z_G2x2y;
  abcd[284] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_G2x2y_a;
  abcd[285] = 2.0E0*I_TWOBODYOVERLAP_H4xy_G2xyz_a;
  abcd[286] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_G2xyz_a-1*I_TWOBODYOVERLAP_F3x_G2xyz;
  abcd[287] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G2xyz_a;
  abcd[288] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_G2xyz_a-2*I_TWOBODYOVERLAP_F2xy_G2xyz;
  abcd[289] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G2xyz_a-1*I_TWOBODYOVERLAP_F2xz_G2xyz;
  abcd[290] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G2xyz_a;
  abcd[291] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_G2xyz_a-3*I_TWOBODYOVERLAP_Fx2y_G2xyz;
  abcd[292] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G2xyz_a-2*I_TWOBODYOVERLAP_Fxyz_G2xyz;
  abcd[293] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G2xyz_a-1*I_TWOBODYOVERLAP_Fx2z_G2xyz;
  abcd[294] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G2xyz_a;
  abcd[295] = 2.0E0*I_TWOBODYOVERLAP_H5y_G2xyz_a-4*I_TWOBODYOVERLAP_F3y_G2xyz;
  abcd[296] = 2.0E0*I_TWOBODYOVERLAP_H4yz_G2xyz_a-3*I_TWOBODYOVERLAP_F2yz_G2xyz;
  abcd[297] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_G2xyz_a-2*I_TWOBODYOVERLAP_Fy2z_G2xyz;
  abcd[298] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_G2xyz_a-1*I_TWOBODYOVERLAP_F3z_G2xyz;
  abcd[299] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_G2xyz_a;
  abcd[300] = 2.0E0*I_TWOBODYOVERLAP_H4xy_G2x2z_a;
  abcd[301] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_G2x2z_a-1*I_TWOBODYOVERLAP_F3x_G2x2z;
  abcd[302] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G2x2z_a;
  abcd[303] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_G2x2z_a-2*I_TWOBODYOVERLAP_F2xy_G2x2z;
  abcd[304] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G2x2z_a-1*I_TWOBODYOVERLAP_F2xz_G2x2z;
  abcd[305] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G2x2z_a;
  abcd[306] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_G2x2z_a-3*I_TWOBODYOVERLAP_Fx2y_G2x2z;
  abcd[307] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G2x2z_a-2*I_TWOBODYOVERLAP_Fxyz_G2x2z;
  abcd[308] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G2x2z_a-1*I_TWOBODYOVERLAP_Fx2z_G2x2z;
  abcd[309] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G2x2z_a;
  abcd[310] = 2.0E0*I_TWOBODYOVERLAP_H5y_G2x2z_a-4*I_TWOBODYOVERLAP_F3y_G2x2z;
  abcd[311] = 2.0E0*I_TWOBODYOVERLAP_H4yz_G2x2z_a-3*I_TWOBODYOVERLAP_F2yz_G2x2z;
  abcd[312] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_G2x2z_a-2*I_TWOBODYOVERLAP_Fy2z_G2x2z;
  abcd[313] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_G2x2z_a-1*I_TWOBODYOVERLAP_F3z_G2x2z;
  abcd[314] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_G2x2z_a;
  abcd[315] = 2.0E0*I_TWOBODYOVERLAP_H4xy_Gx3y_a;
  abcd[316] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_Gx3y_a-1*I_TWOBODYOVERLAP_F3x_Gx3y;
  abcd[317] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_Gx3y_a;
  abcd[318] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_Gx3y_a-2*I_TWOBODYOVERLAP_F2xy_Gx3y;
  abcd[319] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_Gx3y_a-1*I_TWOBODYOVERLAP_F2xz_Gx3y;
  abcd[320] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_Gx3y_a;
  abcd[321] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_Gx3y_a-3*I_TWOBODYOVERLAP_Fx2y_Gx3y;
  abcd[322] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_Gx3y_a-2*I_TWOBODYOVERLAP_Fxyz_Gx3y;
  abcd[323] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_Gx3y_a-1*I_TWOBODYOVERLAP_Fx2z_Gx3y;
  abcd[324] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_Gx3y_a;
  abcd[325] = 2.0E0*I_TWOBODYOVERLAP_H5y_Gx3y_a-4*I_TWOBODYOVERLAP_F3y_Gx3y;
  abcd[326] = 2.0E0*I_TWOBODYOVERLAP_H4yz_Gx3y_a-3*I_TWOBODYOVERLAP_F2yz_Gx3y;
  abcd[327] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_Gx3y_a-2*I_TWOBODYOVERLAP_Fy2z_Gx3y;
  abcd[328] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_Gx3y_a-1*I_TWOBODYOVERLAP_F3z_Gx3y;
  abcd[329] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_Gx3y_a;
  abcd[330] = 2.0E0*I_TWOBODYOVERLAP_H4xy_Gx2yz_a;
  abcd[331] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_Gx2yz_a-1*I_TWOBODYOVERLAP_F3x_Gx2yz;
  abcd[332] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_Gx2yz_a;
  abcd[333] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_Gx2yz_a-2*I_TWOBODYOVERLAP_F2xy_Gx2yz;
  abcd[334] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_Gx2yz_a-1*I_TWOBODYOVERLAP_F2xz_Gx2yz;
  abcd[335] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_Gx2yz_a;
  abcd[336] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_Gx2yz_a-3*I_TWOBODYOVERLAP_Fx2y_Gx2yz;
  abcd[337] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_Gx2yz_a-2*I_TWOBODYOVERLAP_Fxyz_Gx2yz;
  abcd[338] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_Gx2yz_a-1*I_TWOBODYOVERLAP_Fx2z_Gx2yz;
  abcd[339] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_Gx2yz_a;
  abcd[340] = 2.0E0*I_TWOBODYOVERLAP_H5y_Gx2yz_a-4*I_TWOBODYOVERLAP_F3y_Gx2yz;
  abcd[341] = 2.0E0*I_TWOBODYOVERLAP_H4yz_Gx2yz_a-3*I_TWOBODYOVERLAP_F2yz_Gx2yz;
  abcd[342] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_Gx2yz_a-2*I_TWOBODYOVERLAP_Fy2z_Gx2yz;
  abcd[343] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_Gx2yz_a-1*I_TWOBODYOVERLAP_F3z_Gx2yz;
  abcd[344] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_Gx2yz_a;
  abcd[345] = 2.0E0*I_TWOBODYOVERLAP_H4xy_Gxy2z_a;
  abcd[346] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_Gxy2z_a-1*I_TWOBODYOVERLAP_F3x_Gxy2z;
  abcd[347] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_Gxy2z_a;
  abcd[348] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_Gxy2z_a-2*I_TWOBODYOVERLAP_F2xy_Gxy2z;
  abcd[349] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_Gxy2z_a-1*I_TWOBODYOVERLAP_F2xz_Gxy2z;
  abcd[350] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_Gxy2z_a;
  abcd[351] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_Gxy2z_a-3*I_TWOBODYOVERLAP_Fx2y_Gxy2z;
  abcd[352] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_Gxy2z_a-2*I_TWOBODYOVERLAP_Fxyz_Gxy2z;
  abcd[353] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_Gxy2z_a-1*I_TWOBODYOVERLAP_Fx2z_Gxy2z;
  abcd[354] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_Gxy2z_a;
  abcd[355] = 2.0E0*I_TWOBODYOVERLAP_H5y_Gxy2z_a-4*I_TWOBODYOVERLAP_F3y_Gxy2z;
  abcd[356] = 2.0E0*I_TWOBODYOVERLAP_H4yz_Gxy2z_a-3*I_TWOBODYOVERLAP_F2yz_Gxy2z;
  abcd[357] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_Gxy2z_a-2*I_TWOBODYOVERLAP_Fy2z_Gxy2z;
  abcd[358] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_Gxy2z_a-1*I_TWOBODYOVERLAP_F3z_Gxy2z;
  abcd[359] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_Gxy2z_a;
  abcd[360] = 2.0E0*I_TWOBODYOVERLAP_H4xy_Gx3z_a;
  abcd[361] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_Gx3z_a-1*I_TWOBODYOVERLAP_F3x_Gx3z;
  abcd[362] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_Gx3z_a;
  abcd[363] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_Gx3z_a-2*I_TWOBODYOVERLAP_F2xy_Gx3z;
  abcd[364] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_Gx3z_a-1*I_TWOBODYOVERLAP_F2xz_Gx3z;
  abcd[365] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_Gx3z_a;
  abcd[366] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_Gx3z_a-3*I_TWOBODYOVERLAP_Fx2y_Gx3z;
  abcd[367] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_Gx3z_a-2*I_TWOBODYOVERLAP_Fxyz_Gx3z;
  abcd[368] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_Gx3z_a-1*I_TWOBODYOVERLAP_Fx2z_Gx3z;
  abcd[369] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_Gx3z_a;
  abcd[370] = 2.0E0*I_TWOBODYOVERLAP_H5y_Gx3z_a-4*I_TWOBODYOVERLAP_F3y_Gx3z;
  abcd[371] = 2.0E0*I_TWOBODYOVERLAP_H4yz_Gx3z_a-3*I_TWOBODYOVERLAP_F2yz_Gx3z;
  abcd[372] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_Gx3z_a-2*I_TWOBODYOVERLAP_Fy2z_Gx3z;
  abcd[373] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_Gx3z_a-1*I_TWOBODYOVERLAP_F3z_Gx3z;
  abcd[374] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_Gx3z_a;
  abcd[375] = 2.0E0*I_TWOBODYOVERLAP_H4xy_G4y_a;
  abcd[376] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_G4y_a-1*I_TWOBODYOVERLAP_F3x_G4y;
  abcd[377] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G4y_a;
  abcd[378] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_G4y_a-2*I_TWOBODYOVERLAP_F2xy_G4y;
  abcd[379] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G4y_a-1*I_TWOBODYOVERLAP_F2xz_G4y;
  abcd[380] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G4y_a;
  abcd[381] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_G4y_a-3*I_TWOBODYOVERLAP_Fx2y_G4y;
  abcd[382] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G4y_a-2*I_TWOBODYOVERLAP_Fxyz_G4y;
  abcd[383] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G4y_a-1*I_TWOBODYOVERLAP_Fx2z_G4y;
  abcd[384] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G4y_a;
  abcd[385] = 2.0E0*I_TWOBODYOVERLAP_H5y_G4y_a-4*I_TWOBODYOVERLAP_F3y_G4y;
  abcd[386] = 2.0E0*I_TWOBODYOVERLAP_H4yz_G4y_a-3*I_TWOBODYOVERLAP_F2yz_G4y;
  abcd[387] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_G4y_a-2*I_TWOBODYOVERLAP_Fy2z_G4y;
  abcd[388] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_G4y_a-1*I_TWOBODYOVERLAP_F3z_G4y;
  abcd[389] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_G4y_a;
  abcd[390] = 2.0E0*I_TWOBODYOVERLAP_H4xy_G3yz_a;
  abcd[391] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_G3yz_a-1*I_TWOBODYOVERLAP_F3x_G3yz;
  abcd[392] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G3yz_a;
  abcd[393] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_G3yz_a-2*I_TWOBODYOVERLAP_F2xy_G3yz;
  abcd[394] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G3yz_a-1*I_TWOBODYOVERLAP_F2xz_G3yz;
  abcd[395] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G3yz_a;
  abcd[396] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_G3yz_a-3*I_TWOBODYOVERLAP_Fx2y_G3yz;
  abcd[397] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G3yz_a-2*I_TWOBODYOVERLAP_Fxyz_G3yz;
  abcd[398] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G3yz_a-1*I_TWOBODYOVERLAP_Fx2z_G3yz;
  abcd[399] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G3yz_a;
  abcd[400] = 2.0E0*I_TWOBODYOVERLAP_H5y_G3yz_a-4*I_TWOBODYOVERLAP_F3y_G3yz;
  abcd[401] = 2.0E0*I_TWOBODYOVERLAP_H4yz_G3yz_a-3*I_TWOBODYOVERLAP_F2yz_G3yz;
  abcd[402] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_G3yz_a-2*I_TWOBODYOVERLAP_Fy2z_G3yz;
  abcd[403] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_G3yz_a-1*I_TWOBODYOVERLAP_F3z_G3yz;
  abcd[404] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_G3yz_a;
  abcd[405] = 2.0E0*I_TWOBODYOVERLAP_H4xy_G2y2z_a;
  abcd[406] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_G2y2z_a-1*I_TWOBODYOVERLAP_F3x_G2y2z;
  abcd[407] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G2y2z_a;
  abcd[408] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_G2y2z_a-2*I_TWOBODYOVERLAP_F2xy_G2y2z;
  abcd[409] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G2y2z_a-1*I_TWOBODYOVERLAP_F2xz_G2y2z;
  abcd[410] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G2y2z_a;
  abcd[411] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_G2y2z_a-3*I_TWOBODYOVERLAP_Fx2y_G2y2z;
  abcd[412] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G2y2z_a-2*I_TWOBODYOVERLAP_Fxyz_G2y2z;
  abcd[413] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G2y2z_a-1*I_TWOBODYOVERLAP_Fx2z_G2y2z;
  abcd[414] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G2y2z_a;
  abcd[415] = 2.0E0*I_TWOBODYOVERLAP_H5y_G2y2z_a-4*I_TWOBODYOVERLAP_F3y_G2y2z;
  abcd[416] = 2.0E0*I_TWOBODYOVERLAP_H4yz_G2y2z_a-3*I_TWOBODYOVERLAP_F2yz_G2y2z;
  abcd[417] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_G2y2z_a-2*I_TWOBODYOVERLAP_Fy2z_G2y2z;
  abcd[418] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_G2y2z_a-1*I_TWOBODYOVERLAP_F3z_G2y2z;
  abcd[419] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_G2y2z_a;
  abcd[420] = 2.0E0*I_TWOBODYOVERLAP_H4xy_Gy3z_a;
  abcd[421] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_Gy3z_a-1*I_TWOBODYOVERLAP_F3x_Gy3z;
  abcd[422] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_Gy3z_a;
  abcd[423] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_Gy3z_a-2*I_TWOBODYOVERLAP_F2xy_Gy3z;
  abcd[424] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_Gy3z_a-1*I_TWOBODYOVERLAP_F2xz_Gy3z;
  abcd[425] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_Gy3z_a;
  abcd[426] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_Gy3z_a-3*I_TWOBODYOVERLAP_Fx2y_Gy3z;
  abcd[427] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_Gy3z_a-2*I_TWOBODYOVERLAP_Fxyz_Gy3z;
  abcd[428] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_Gy3z_a-1*I_TWOBODYOVERLAP_Fx2z_Gy3z;
  abcd[429] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_Gy3z_a;
  abcd[430] = 2.0E0*I_TWOBODYOVERLAP_H5y_Gy3z_a-4*I_TWOBODYOVERLAP_F3y_Gy3z;
  abcd[431] = 2.0E0*I_TWOBODYOVERLAP_H4yz_Gy3z_a-3*I_TWOBODYOVERLAP_F2yz_Gy3z;
  abcd[432] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_Gy3z_a-2*I_TWOBODYOVERLAP_Fy2z_Gy3z;
  abcd[433] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_Gy3z_a-1*I_TWOBODYOVERLAP_F3z_Gy3z;
  abcd[434] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_Gy3z_a;
  abcd[435] = 2.0E0*I_TWOBODYOVERLAP_H4xy_G4z_a;
  abcd[436] = 2.0E0*I_TWOBODYOVERLAP_H3x2y_G4z_a-1*I_TWOBODYOVERLAP_F3x_G4z;
  abcd[437] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G4z_a;
  abcd[438] = 2.0E0*I_TWOBODYOVERLAP_H2x3y_G4z_a-2*I_TWOBODYOVERLAP_F2xy_G4z;
  abcd[439] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G4z_a-1*I_TWOBODYOVERLAP_F2xz_G4z;
  abcd[440] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G4z_a;
  abcd[441] = 2.0E0*I_TWOBODYOVERLAP_Hx4y_G4z_a-3*I_TWOBODYOVERLAP_Fx2y_G4z;
  abcd[442] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G4z_a-2*I_TWOBODYOVERLAP_Fxyz_G4z;
  abcd[443] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G4z_a-1*I_TWOBODYOVERLAP_Fx2z_G4z;
  abcd[444] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G4z_a;
  abcd[445] = 2.0E0*I_TWOBODYOVERLAP_H5y_G4z_a-4*I_TWOBODYOVERLAP_F3y_G4z;
  abcd[446] = 2.0E0*I_TWOBODYOVERLAP_H4yz_G4z_a-3*I_TWOBODYOVERLAP_F2yz_G4z;
  abcd[447] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_G4z_a-2*I_TWOBODYOVERLAP_Fy2z_G4z;
  abcd[448] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_G4z_a-1*I_TWOBODYOVERLAP_F3z_G4z;
  abcd[449] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_G4z_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_G_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_G_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_G
   ************************************************************/
  abcd[450] = 2.0E0*I_TWOBODYOVERLAP_H4xz_G4x_a;
  abcd[451] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G4x_a;
  abcd[452] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_G4x_a-1*I_TWOBODYOVERLAP_F3x_G4x;
  abcd[453] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G4x_a;
  abcd[454] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G4x_a-1*I_TWOBODYOVERLAP_F2xy_G4x;
  abcd[455] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_G4x_a-2*I_TWOBODYOVERLAP_F2xz_G4x;
  abcd[456] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G4x_a;
  abcd[457] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G4x_a-1*I_TWOBODYOVERLAP_Fx2y_G4x;
  abcd[458] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G4x_a-2*I_TWOBODYOVERLAP_Fxyz_G4x;
  abcd[459] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_G4x_a-3*I_TWOBODYOVERLAP_Fx2z_G4x;
  abcd[460] = 2.0E0*I_TWOBODYOVERLAP_H4yz_G4x_a;
  abcd[461] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_G4x_a-1*I_TWOBODYOVERLAP_F3y_G4x;
  abcd[462] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_G4x_a-2*I_TWOBODYOVERLAP_F2yz_G4x;
  abcd[463] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_G4x_a-3*I_TWOBODYOVERLAP_Fy2z_G4x;
  abcd[464] = 2.0E0*I_TWOBODYOVERLAP_H5z_G4x_a-4*I_TWOBODYOVERLAP_F3z_G4x;
  abcd[465] = 2.0E0*I_TWOBODYOVERLAP_H4xz_G3xy_a;
  abcd[466] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G3xy_a;
  abcd[467] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_G3xy_a-1*I_TWOBODYOVERLAP_F3x_G3xy;
  abcd[468] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G3xy_a;
  abcd[469] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G3xy_a-1*I_TWOBODYOVERLAP_F2xy_G3xy;
  abcd[470] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_G3xy_a-2*I_TWOBODYOVERLAP_F2xz_G3xy;
  abcd[471] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G3xy_a;
  abcd[472] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G3xy_a-1*I_TWOBODYOVERLAP_Fx2y_G3xy;
  abcd[473] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G3xy_a-2*I_TWOBODYOVERLAP_Fxyz_G3xy;
  abcd[474] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_G3xy_a-3*I_TWOBODYOVERLAP_Fx2z_G3xy;
  abcd[475] = 2.0E0*I_TWOBODYOVERLAP_H4yz_G3xy_a;
  abcd[476] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_G3xy_a-1*I_TWOBODYOVERLAP_F3y_G3xy;
  abcd[477] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_G3xy_a-2*I_TWOBODYOVERLAP_F2yz_G3xy;
  abcd[478] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_G3xy_a-3*I_TWOBODYOVERLAP_Fy2z_G3xy;
  abcd[479] = 2.0E0*I_TWOBODYOVERLAP_H5z_G3xy_a-4*I_TWOBODYOVERLAP_F3z_G3xy;
  abcd[480] = 2.0E0*I_TWOBODYOVERLAP_H4xz_G3xz_a;
  abcd[481] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G3xz_a;
  abcd[482] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_G3xz_a-1*I_TWOBODYOVERLAP_F3x_G3xz;
  abcd[483] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G3xz_a;
  abcd[484] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G3xz_a-1*I_TWOBODYOVERLAP_F2xy_G3xz;
  abcd[485] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_G3xz_a-2*I_TWOBODYOVERLAP_F2xz_G3xz;
  abcd[486] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G3xz_a;
  abcd[487] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G3xz_a-1*I_TWOBODYOVERLAP_Fx2y_G3xz;
  abcd[488] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G3xz_a-2*I_TWOBODYOVERLAP_Fxyz_G3xz;
  abcd[489] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_G3xz_a-3*I_TWOBODYOVERLAP_Fx2z_G3xz;
  abcd[490] = 2.0E0*I_TWOBODYOVERLAP_H4yz_G3xz_a;
  abcd[491] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_G3xz_a-1*I_TWOBODYOVERLAP_F3y_G3xz;
  abcd[492] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_G3xz_a-2*I_TWOBODYOVERLAP_F2yz_G3xz;
  abcd[493] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_G3xz_a-3*I_TWOBODYOVERLAP_Fy2z_G3xz;
  abcd[494] = 2.0E0*I_TWOBODYOVERLAP_H5z_G3xz_a-4*I_TWOBODYOVERLAP_F3z_G3xz;
  abcd[495] = 2.0E0*I_TWOBODYOVERLAP_H4xz_G2x2y_a;
  abcd[496] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G2x2y_a;
  abcd[497] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_G2x2y_a-1*I_TWOBODYOVERLAP_F3x_G2x2y;
  abcd[498] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G2x2y_a;
  abcd[499] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G2x2y_a-1*I_TWOBODYOVERLAP_F2xy_G2x2y;
  abcd[500] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_G2x2y_a-2*I_TWOBODYOVERLAP_F2xz_G2x2y;
  abcd[501] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G2x2y_a;
  abcd[502] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G2x2y_a-1*I_TWOBODYOVERLAP_Fx2y_G2x2y;
  abcd[503] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G2x2y_a-2*I_TWOBODYOVERLAP_Fxyz_G2x2y;
  abcd[504] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_G2x2y_a-3*I_TWOBODYOVERLAP_Fx2z_G2x2y;
  abcd[505] = 2.0E0*I_TWOBODYOVERLAP_H4yz_G2x2y_a;
  abcd[506] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_G2x2y_a-1*I_TWOBODYOVERLAP_F3y_G2x2y;
  abcd[507] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_G2x2y_a-2*I_TWOBODYOVERLAP_F2yz_G2x2y;
  abcd[508] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_G2x2y_a-3*I_TWOBODYOVERLAP_Fy2z_G2x2y;
  abcd[509] = 2.0E0*I_TWOBODYOVERLAP_H5z_G2x2y_a-4*I_TWOBODYOVERLAP_F3z_G2x2y;
  abcd[510] = 2.0E0*I_TWOBODYOVERLAP_H4xz_G2xyz_a;
  abcd[511] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G2xyz_a;
  abcd[512] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_G2xyz_a-1*I_TWOBODYOVERLAP_F3x_G2xyz;
  abcd[513] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G2xyz_a;
  abcd[514] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G2xyz_a-1*I_TWOBODYOVERLAP_F2xy_G2xyz;
  abcd[515] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_G2xyz_a-2*I_TWOBODYOVERLAP_F2xz_G2xyz;
  abcd[516] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G2xyz_a;
  abcd[517] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G2xyz_a-1*I_TWOBODYOVERLAP_Fx2y_G2xyz;
  abcd[518] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G2xyz_a-2*I_TWOBODYOVERLAP_Fxyz_G2xyz;
  abcd[519] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_G2xyz_a-3*I_TWOBODYOVERLAP_Fx2z_G2xyz;
  abcd[520] = 2.0E0*I_TWOBODYOVERLAP_H4yz_G2xyz_a;
  abcd[521] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_G2xyz_a-1*I_TWOBODYOVERLAP_F3y_G2xyz;
  abcd[522] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_G2xyz_a-2*I_TWOBODYOVERLAP_F2yz_G2xyz;
  abcd[523] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_G2xyz_a-3*I_TWOBODYOVERLAP_Fy2z_G2xyz;
  abcd[524] = 2.0E0*I_TWOBODYOVERLAP_H5z_G2xyz_a-4*I_TWOBODYOVERLAP_F3z_G2xyz;
  abcd[525] = 2.0E0*I_TWOBODYOVERLAP_H4xz_G2x2z_a;
  abcd[526] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G2x2z_a;
  abcd[527] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_G2x2z_a-1*I_TWOBODYOVERLAP_F3x_G2x2z;
  abcd[528] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G2x2z_a;
  abcd[529] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G2x2z_a-1*I_TWOBODYOVERLAP_F2xy_G2x2z;
  abcd[530] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_G2x2z_a-2*I_TWOBODYOVERLAP_F2xz_G2x2z;
  abcd[531] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G2x2z_a;
  abcd[532] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G2x2z_a-1*I_TWOBODYOVERLAP_Fx2y_G2x2z;
  abcd[533] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G2x2z_a-2*I_TWOBODYOVERLAP_Fxyz_G2x2z;
  abcd[534] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_G2x2z_a-3*I_TWOBODYOVERLAP_Fx2z_G2x2z;
  abcd[535] = 2.0E0*I_TWOBODYOVERLAP_H4yz_G2x2z_a;
  abcd[536] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_G2x2z_a-1*I_TWOBODYOVERLAP_F3y_G2x2z;
  abcd[537] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_G2x2z_a-2*I_TWOBODYOVERLAP_F2yz_G2x2z;
  abcd[538] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_G2x2z_a-3*I_TWOBODYOVERLAP_Fy2z_G2x2z;
  abcd[539] = 2.0E0*I_TWOBODYOVERLAP_H5z_G2x2z_a-4*I_TWOBODYOVERLAP_F3z_G2x2z;
  abcd[540] = 2.0E0*I_TWOBODYOVERLAP_H4xz_Gx3y_a;
  abcd[541] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_Gx3y_a;
  abcd[542] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_Gx3y_a-1*I_TWOBODYOVERLAP_F3x_Gx3y;
  abcd[543] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_Gx3y_a;
  abcd[544] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_Gx3y_a-1*I_TWOBODYOVERLAP_F2xy_Gx3y;
  abcd[545] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_Gx3y_a-2*I_TWOBODYOVERLAP_F2xz_Gx3y;
  abcd[546] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_Gx3y_a;
  abcd[547] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_Gx3y_a-1*I_TWOBODYOVERLAP_Fx2y_Gx3y;
  abcd[548] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_Gx3y_a-2*I_TWOBODYOVERLAP_Fxyz_Gx3y;
  abcd[549] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_Gx3y_a-3*I_TWOBODYOVERLAP_Fx2z_Gx3y;
  abcd[550] = 2.0E0*I_TWOBODYOVERLAP_H4yz_Gx3y_a;
  abcd[551] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_Gx3y_a-1*I_TWOBODYOVERLAP_F3y_Gx3y;
  abcd[552] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_Gx3y_a-2*I_TWOBODYOVERLAP_F2yz_Gx3y;
  abcd[553] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_Gx3y_a-3*I_TWOBODYOVERLAP_Fy2z_Gx3y;
  abcd[554] = 2.0E0*I_TWOBODYOVERLAP_H5z_Gx3y_a-4*I_TWOBODYOVERLAP_F3z_Gx3y;
  abcd[555] = 2.0E0*I_TWOBODYOVERLAP_H4xz_Gx2yz_a;
  abcd[556] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_Gx2yz_a;
  abcd[557] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_Gx2yz_a-1*I_TWOBODYOVERLAP_F3x_Gx2yz;
  abcd[558] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_Gx2yz_a;
  abcd[559] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_Gx2yz_a-1*I_TWOBODYOVERLAP_F2xy_Gx2yz;
  abcd[560] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_Gx2yz_a-2*I_TWOBODYOVERLAP_F2xz_Gx2yz;
  abcd[561] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_Gx2yz_a;
  abcd[562] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_Gx2yz_a-1*I_TWOBODYOVERLAP_Fx2y_Gx2yz;
  abcd[563] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_Gx2yz_a-2*I_TWOBODYOVERLAP_Fxyz_Gx2yz;
  abcd[564] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_Gx2yz_a-3*I_TWOBODYOVERLAP_Fx2z_Gx2yz;
  abcd[565] = 2.0E0*I_TWOBODYOVERLAP_H4yz_Gx2yz_a;
  abcd[566] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_Gx2yz_a-1*I_TWOBODYOVERLAP_F3y_Gx2yz;
  abcd[567] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_Gx2yz_a-2*I_TWOBODYOVERLAP_F2yz_Gx2yz;
  abcd[568] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_Gx2yz_a-3*I_TWOBODYOVERLAP_Fy2z_Gx2yz;
  abcd[569] = 2.0E0*I_TWOBODYOVERLAP_H5z_Gx2yz_a-4*I_TWOBODYOVERLAP_F3z_Gx2yz;
  abcd[570] = 2.0E0*I_TWOBODYOVERLAP_H4xz_Gxy2z_a;
  abcd[571] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_Gxy2z_a;
  abcd[572] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_Gxy2z_a-1*I_TWOBODYOVERLAP_F3x_Gxy2z;
  abcd[573] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_Gxy2z_a;
  abcd[574] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_Gxy2z_a-1*I_TWOBODYOVERLAP_F2xy_Gxy2z;
  abcd[575] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_Gxy2z_a-2*I_TWOBODYOVERLAP_F2xz_Gxy2z;
  abcd[576] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_Gxy2z_a;
  abcd[577] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_Gxy2z_a-1*I_TWOBODYOVERLAP_Fx2y_Gxy2z;
  abcd[578] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_Gxy2z_a-2*I_TWOBODYOVERLAP_Fxyz_Gxy2z;
  abcd[579] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_Gxy2z_a-3*I_TWOBODYOVERLAP_Fx2z_Gxy2z;
  abcd[580] = 2.0E0*I_TWOBODYOVERLAP_H4yz_Gxy2z_a;
  abcd[581] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_Gxy2z_a-1*I_TWOBODYOVERLAP_F3y_Gxy2z;
  abcd[582] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_Gxy2z_a-2*I_TWOBODYOVERLAP_F2yz_Gxy2z;
  abcd[583] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_Gxy2z_a-3*I_TWOBODYOVERLAP_Fy2z_Gxy2z;
  abcd[584] = 2.0E0*I_TWOBODYOVERLAP_H5z_Gxy2z_a-4*I_TWOBODYOVERLAP_F3z_Gxy2z;
  abcd[585] = 2.0E0*I_TWOBODYOVERLAP_H4xz_Gx3z_a;
  abcd[586] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_Gx3z_a;
  abcd[587] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_Gx3z_a-1*I_TWOBODYOVERLAP_F3x_Gx3z;
  abcd[588] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_Gx3z_a;
  abcd[589] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_Gx3z_a-1*I_TWOBODYOVERLAP_F2xy_Gx3z;
  abcd[590] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_Gx3z_a-2*I_TWOBODYOVERLAP_F2xz_Gx3z;
  abcd[591] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_Gx3z_a;
  abcd[592] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_Gx3z_a-1*I_TWOBODYOVERLAP_Fx2y_Gx3z;
  abcd[593] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_Gx3z_a-2*I_TWOBODYOVERLAP_Fxyz_Gx3z;
  abcd[594] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_Gx3z_a-3*I_TWOBODYOVERLAP_Fx2z_Gx3z;
  abcd[595] = 2.0E0*I_TWOBODYOVERLAP_H4yz_Gx3z_a;
  abcd[596] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_Gx3z_a-1*I_TWOBODYOVERLAP_F3y_Gx3z;
  abcd[597] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_Gx3z_a-2*I_TWOBODYOVERLAP_F2yz_Gx3z;
  abcd[598] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_Gx3z_a-3*I_TWOBODYOVERLAP_Fy2z_Gx3z;
  abcd[599] = 2.0E0*I_TWOBODYOVERLAP_H5z_Gx3z_a-4*I_TWOBODYOVERLAP_F3z_Gx3z;
  abcd[600] = 2.0E0*I_TWOBODYOVERLAP_H4xz_G4y_a;
  abcd[601] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G4y_a;
  abcd[602] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_G4y_a-1*I_TWOBODYOVERLAP_F3x_G4y;
  abcd[603] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G4y_a;
  abcd[604] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G4y_a-1*I_TWOBODYOVERLAP_F2xy_G4y;
  abcd[605] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_G4y_a-2*I_TWOBODYOVERLAP_F2xz_G4y;
  abcd[606] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G4y_a;
  abcd[607] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G4y_a-1*I_TWOBODYOVERLAP_Fx2y_G4y;
  abcd[608] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G4y_a-2*I_TWOBODYOVERLAP_Fxyz_G4y;
  abcd[609] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_G4y_a-3*I_TWOBODYOVERLAP_Fx2z_G4y;
  abcd[610] = 2.0E0*I_TWOBODYOVERLAP_H4yz_G4y_a;
  abcd[611] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_G4y_a-1*I_TWOBODYOVERLAP_F3y_G4y;
  abcd[612] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_G4y_a-2*I_TWOBODYOVERLAP_F2yz_G4y;
  abcd[613] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_G4y_a-3*I_TWOBODYOVERLAP_Fy2z_G4y;
  abcd[614] = 2.0E0*I_TWOBODYOVERLAP_H5z_G4y_a-4*I_TWOBODYOVERLAP_F3z_G4y;
  abcd[615] = 2.0E0*I_TWOBODYOVERLAP_H4xz_G3yz_a;
  abcd[616] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G3yz_a;
  abcd[617] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_G3yz_a-1*I_TWOBODYOVERLAP_F3x_G3yz;
  abcd[618] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G3yz_a;
  abcd[619] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G3yz_a-1*I_TWOBODYOVERLAP_F2xy_G3yz;
  abcd[620] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_G3yz_a-2*I_TWOBODYOVERLAP_F2xz_G3yz;
  abcd[621] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G3yz_a;
  abcd[622] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G3yz_a-1*I_TWOBODYOVERLAP_Fx2y_G3yz;
  abcd[623] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G3yz_a-2*I_TWOBODYOVERLAP_Fxyz_G3yz;
  abcd[624] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_G3yz_a-3*I_TWOBODYOVERLAP_Fx2z_G3yz;
  abcd[625] = 2.0E0*I_TWOBODYOVERLAP_H4yz_G3yz_a;
  abcd[626] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_G3yz_a-1*I_TWOBODYOVERLAP_F3y_G3yz;
  abcd[627] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_G3yz_a-2*I_TWOBODYOVERLAP_F2yz_G3yz;
  abcd[628] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_G3yz_a-3*I_TWOBODYOVERLAP_Fy2z_G3yz;
  abcd[629] = 2.0E0*I_TWOBODYOVERLAP_H5z_G3yz_a-4*I_TWOBODYOVERLAP_F3z_G3yz;
  abcd[630] = 2.0E0*I_TWOBODYOVERLAP_H4xz_G2y2z_a;
  abcd[631] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G2y2z_a;
  abcd[632] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_G2y2z_a-1*I_TWOBODYOVERLAP_F3x_G2y2z;
  abcd[633] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G2y2z_a;
  abcd[634] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G2y2z_a-1*I_TWOBODYOVERLAP_F2xy_G2y2z;
  abcd[635] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_G2y2z_a-2*I_TWOBODYOVERLAP_F2xz_G2y2z;
  abcd[636] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G2y2z_a;
  abcd[637] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G2y2z_a-1*I_TWOBODYOVERLAP_Fx2y_G2y2z;
  abcd[638] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G2y2z_a-2*I_TWOBODYOVERLAP_Fxyz_G2y2z;
  abcd[639] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_G2y2z_a-3*I_TWOBODYOVERLAP_Fx2z_G2y2z;
  abcd[640] = 2.0E0*I_TWOBODYOVERLAP_H4yz_G2y2z_a;
  abcd[641] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_G2y2z_a-1*I_TWOBODYOVERLAP_F3y_G2y2z;
  abcd[642] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_G2y2z_a-2*I_TWOBODYOVERLAP_F2yz_G2y2z;
  abcd[643] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_G2y2z_a-3*I_TWOBODYOVERLAP_Fy2z_G2y2z;
  abcd[644] = 2.0E0*I_TWOBODYOVERLAP_H5z_G2y2z_a-4*I_TWOBODYOVERLAP_F3z_G2y2z;
  abcd[645] = 2.0E0*I_TWOBODYOVERLAP_H4xz_Gy3z_a;
  abcd[646] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_Gy3z_a;
  abcd[647] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_Gy3z_a-1*I_TWOBODYOVERLAP_F3x_Gy3z;
  abcd[648] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_Gy3z_a;
  abcd[649] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_Gy3z_a-1*I_TWOBODYOVERLAP_F2xy_Gy3z;
  abcd[650] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_Gy3z_a-2*I_TWOBODYOVERLAP_F2xz_Gy3z;
  abcd[651] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_Gy3z_a;
  abcd[652] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_Gy3z_a-1*I_TWOBODYOVERLAP_Fx2y_Gy3z;
  abcd[653] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_Gy3z_a-2*I_TWOBODYOVERLAP_Fxyz_Gy3z;
  abcd[654] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_Gy3z_a-3*I_TWOBODYOVERLAP_Fx2z_Gy3z;
  abcd[655] = 2.0E0*I_TWOBODYOVERLAP_H4yz_Gy3z_a;
  abcd[656] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_Gy3z_a-1*I_TWOBODYOVERLAP_F3y_Gy3z;
  abcd[657] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_Gy3z_a-2*I_TWOBODYOVERLAP_F2yz_Gy3z;
  abcd[658] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_Gy3z_a-3*I_TWOBODYOVERLAP_Fy2z_Gy3z;
  abcd[659] = 2.0E0*I_TWOBODYOVERLAP_H5z_Gy3z_a-4*I_TWOBODYOVERLAP_F3z_Gy3z;
  abcd[660] = 2.0E0*I_TWOBODYOVERLAP_H4xz_G4z_a;
  abcd[661] = 2.0E0*I_TWOBODYOVERLAP_H3xyz_G4z_a;
  abcd[662] = 2.0E0*I_TWOBODYOVERLAP_H3x2z_G4z_a-1*I_TWOBODYOVERLAP_F3x_G4z;
  abcd[663] = 2.0E0*I_TWOBODYOVERLAP_H2x2yz_G4z_a;
  abcd[664] = 2.0E0*I_TWOBODYOVERLAP_H2xy2z_G4z_a-1*I_TWOBODYOVERLAP_F2xy_G4z;
  abcd[665] = 2.0E0*I_TWOBODYOVERLAP_H2x3z_G4z_a-2*I_TWOBODYOVERLAP_F2xz_G4z;
  abcd[666] = 2.0E0*I_TWOBODYOVERLAP_Hx3yz_G4z_a;
  abcd[667] = 2.0E0*I_TWOBODYOVERLAP_Hx2y2z_G4z_a-1*I_TWOBODYOVERLAP_Fx2y_G4z;
  abcd[668] = 2.0E0*I_TWOBODYOVERLAP_Hxy3z_G4z_a-2*I_TWOBODYOVERLAP_Fxyz_G4z;
  abcd[669] = 2.0E0*I_TWOBODYOVERLAP_Hx4z_G4z_a-3*I_TWOBODYOVERLAP_Fx2z_G4z;
  abcd[670] = 2.0E0*I_TWOBODYOVERLAP_H4yz_G4z_a;
  abcd[671] = 2.0E0*I_TWOBODYOVERLAP_H3y2z_G4z_a-1*I_TWOBODYOVERLAP_F3y_G4z;
  abcd[672] = 2.0E0*I_TWOBODYOVERLAP_H2y3z_G4z_a-2*I_TWOBODYOVERLAP_F2yz_G4z;
  abcd[673] = 2.0E0*I_TWOBODYOVERLAP_Hy4z_G4z_a-3*I_TWOBODYOVERLAP_Fy2z_G4z;
  abcd[674] = 2.0E0*I_TWOBODYOVERLAP_H5z_G4z_a-4*I_TWOBODYOVERLAP_F3z_G4z;
}
