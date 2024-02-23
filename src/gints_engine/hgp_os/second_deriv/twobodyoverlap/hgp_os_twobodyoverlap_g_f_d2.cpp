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
// BRA1 as redundant position, total RHS integrals evaluated as: 7126
// BRA2 as redundant position, total RHS integrals evaluated as: 5936
// KET1 as redundant position, total RHS integrals evaluated as: 0
// KET2 as redundant position, total RHS integrals evaluated as: 0
// the redundant position is: BRA2
//

//
// @@@@ derivative position-direction information
// BRA1  BRA1
// X  X
// X  Y
// X  Z
// Y  Y
// Y  Z
// Z  Z
// ####

void hgp_os_twobodyoverlap_g_f_d2(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
  //
  // declare the variables as result of VRR process
  //
  Double I_TWOBODYOVERLAP_M9x_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M8xy_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M8xz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M7x2y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M7xyz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M7x2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M6x3y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M6x2yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M6xy2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M6x3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M5x4y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M5x3yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M5x2y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M5xy3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M5x4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M4x5y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M4x4yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M4x3y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M4x2y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M4xy4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M4x5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M3x6y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M3x5yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M3x4y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M3x3y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M3x2y4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M3xy5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M3x6z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x7y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x6yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x5y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x4y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x3y4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x2y5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M2xy6z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M2x7z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx8y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx7yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx6y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx5y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx4y4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx3y5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx2y6z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Mxy7z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Mx8z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M9y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M8yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M7y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M6y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M5y4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M4y5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M3y6z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M2y7z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_My8z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_M9z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L8x_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L7xy_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L7xz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L6x2y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L6xyz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L6x2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L5x3y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L5x2yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L5xy2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L5x3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L4x4y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L4x3yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L4x2y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L4xy3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L4x4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x5y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x4yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x3y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x2y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3xy4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3x5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x6y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x5yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x4y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x3y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x2y4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2xy5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2x6z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx7y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx6yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx5y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx4y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx3y4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx2y5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lxy6z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Lx7z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L8y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L7yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L6y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L5y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L4y4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L3y5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L2y6z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ly7z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_L8z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K7x_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K6xy_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K6xz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K5x2y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K5xyz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K5x2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x3y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x2yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K4xy2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K4x3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x4y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x3yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x2y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3xy3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3x4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x5y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x4yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x3y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x2y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2xy4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2x5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx6y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx5yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx4y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx3y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx2y4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kxy5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Kx6z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K7y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K6yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K5y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K4y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K3y4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K2y5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ky6z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_K7z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I6x_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xy_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I5xz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I4xyz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I4x2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x2yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I3xy2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I3x3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x3yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x2y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I2xy3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I2x4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix4yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix3y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix2y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ixy4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Ix5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I6y_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I5yz_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I4y2z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I3y3z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I2y4z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_Iy5z_S_aa = 0.0E0;
  Double I_TWOBODYOVERLAP_I6z_S_aa = 0.0E0;
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
  Double I_TWOBODYOVERLAP_G4x_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xy_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G3xz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G2xyz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G2x2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx2yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Gxy2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Gx3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G4y_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G3yz_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G2y2z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_Gy3z_S_a = 0.0E0;
  Double I_TWOBODYOVERLAP_G4z_S_a = 0.0E0;
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
  Double I_TWOBODYOVERLAP_D2x_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Dxy_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Dxz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_D2y_S = 0.0E0;
  Double I_TWOBODYOVERLAP_Dyz_S = 0.0E0;
  Double I_TWOBODYOVERLAP_D2z_S = 0.0E0;

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
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_D2x_S_vrr = PAX*I_TWOBODYOVERLAP_Px_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Dxy_S_vrr = PAY*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_Dxz_S_vrr = PAZ*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_D2y_S_vrr = PAY*I_TWOBODYOVERLAP_Py_S_vrr+oned2z*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_TWOBODYOVERLAP_Dyz_S_vrr = PAZ*I_TWOBODYOVERLAP_Py_S_vrr;
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
     * shell quartet name: SQ_TWOBODYOVERLAP_M_S_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_M_S_aa_coefs = alpha*alpha;
    I_TWOBODYOVERLAP_M9x_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M9x_S_vrr;
    I_TWOBODYOVERLAP_M8xy_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M8xy_S_vrr;
    I_TWOBODYOVERLAP_M8xz_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M8xz_S_vrr;
    I_TWOBODYOVERLAP_M7x2y_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M7x2y_S_vrr;
    I_TWOBODYOVERLAP_M7xyz_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M7xyz_S_vrr;
    I_TWOBODYOVERLAP_M7x2z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M7x2z_S_vrr;
    I_TWOBODYOVERLAP_M6x3y_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M6x3y_S_vrr;
    I_TWOBODYOVERLAP_M6x2yz_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M6x2yz_S_vrr;
    I_TWOBODYOVERLAP_M6xy2z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M6xy2z_S_vrr;
    I_TWOBODYOVERLAP_M6x3z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M6x3z_S_vrr;
    I_TWOBODYOVERLAP_M5x4y_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M5x4y_S_vrr;
    I_TWOBODYOVERLAP_M5x3yz_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M5x3yz_S_vrr;
    I_TWOBODYOVERLAP_M5x2y2z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M5x2y2z_S_vrr;
    I_TWOBODYOVERLAP_M5xy3z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M5xy3z_S_vrr;
    I_TWOBODYOVERLAP_M5x4z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M5x4z_S_vrr;
    I_TWOBODYOVERLAP_M4x5y_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M4x5y_S_vrr;
    I_TWOBODYOVERLAP_M4x4yz_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M4x4yz_S_vrr;
    I_TWOBODYOVERLAP_M4x3y2z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M4x3y2z_S_vrr;
    I_TWOBODYOVERLAP_M4x2y3z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M4x2y3z_S_vrr;
    I_TWOBODYOVERLAP_M4xy4z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M4xy4z_S_vrr;
    I_TWOBODYOVERLAP_M4x5z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M4x5z_S_vrr;
    I_TWOBODYOVERLAP_M3x6y_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M3x6y_S_vrr;
    I_TWOBODYOVERLAP_M3x5yz_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M3x5yz_S_vrr;
    I_TWOBODYOVERLAP_M3x4y2z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M3x4y2z_S_vrr;
    I_TWOBODYOVERLAP_M3x3y3z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M3x3y3z_S_vrr;
    I_TWOBODYOVERLAP_M3x2y4z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M3x2y4z_S_vrr;
    I_TWOBODYOVERLAP_M3xy5z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M3xy5z_S_vrr;
    I_TWOBODYOVERLAP_M3x6z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M3x6z_S_vrr;
    I_TWOBODYOVERLAP_M2x7y_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M2x7y_S_vrr;
    I_TWOBODYOVERLAP_M2x6yz_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M2x6yz_S_vrr;
    I_TWOBODYOVERLAP_M2x5y2z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M2x5y2z_S_vrr;
    I_TWOBODYOVERLAP_M2x4y3z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M2x4y3z_S_vrr;
    I_TWOBODYOVERLAP_M2x3y4z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M2x3y4z_S_vrr;
    I_TWOBODYOVERLAP_M2x2y5z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M2x2y5z_S_vrr;
    I_TWOBODYOVERLAP_M2xy6z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M2xy6z_S_vrr;
    I_TWOBODYOVERLAP_M2x7z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M2x7z_S_vrr;
    I_TWOBODYOVERLAP_Mx8y_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_Mx8y_S_vrr;
    I_TWOBODYOVERLAP_Mx7yz_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_Mx7yz_S_vrr;
    I_TWOBODYOVERLAP_Mx6y2z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_Mx6y2z_S_vrr;
    I_TWOBODYOVERLAP_Mx5y3z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_Mx5y3z_S_vrr;
    I_TWOBODYOVERLAP_Mx4y4z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_Mx4y4z_S_vrr;
    I_TWOBODYOVERLAP_Mx3y5z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_Mx3y5z_S_vrr;
    I_TWOBODYOVERLAP_Mx2y6z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_Mx2y6z_S_vrr;
    I_TWOBODYOVERLAP_Mxy7z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_Mxy7z_S_vrr;
    I_TWOBODYOVERLAP_Mx8z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_Mx8z_S_vrr;
    I_TWOBODYOVERLAP_M9y_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M9y_S_vrr;
    I_TWOBODYOVERLAP_M8yz_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M8yz_S_vrr;
    I_TWOBODYOVERLAP_M7y2z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M7y2z_S_vrr;
    I_TWOBODYOVERLAP_M6y3z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M6y3z_S_vrr;
    I_TWOBODYOVERLAP_M5y4z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M5y4z_S_vrr;
    I_TWOBODYOVERLAP_M4y5z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M4y5z_S_vrr;
    I_TWOBODYOVERLAP_M3y6z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M3y6z_S_vrr;
    I_TWOBODYOVERLAP_M2y7z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M2y7z_S_vrr;
    I_TWOBODYOVERLAP_My8z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_My8z_S_vrr;
    I_TWOBODYOVERLAP_M9z_S_aa += SQ_TWOBODYOVERLAP_M_S_aa_coefs*I_TWOBODYOVERLAP_M9z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_L_S_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_L_S_aa_coefs = alpha*alpha;
    I_TWOBODYOVERLAP_L8x_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L8x_S_vrr;
    I_TWOBODYOVERLAP_L7xy_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L7xy_S_vrr;
    I_TWOBODYOVERLAP_L7xz_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L7xz_S_vrr;
    I_TWOBODYOVERLAP_L6x2y_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L6x2y_S_vrr;
    I_TWOBODYOVERLAP_L6xyz_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L6xyz_S_vrr;
    I_TWOBODYOVERLAP_L6x2z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L6x2z_S_vrr;
    I_TWOBODYOVERLAP_L5x3y_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L5x3y_S_vrr;
    I_TWOBODYOVERLAP_L5x2yz_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L5x2yz_S_vrr;
    I_TWOBODYOVERLAP_L5xy2z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L5xy2z_S_vrr;
    I_TWOBODYOVERLAP_L5x3z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L5x3z_S_vrr;
    I_TWOBODYOVERLAP_L4x4y_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L4x4y_S_vrr;
    I_TWOBODYOVERLAP_L4x3yz_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L4x3yz_S_vrr;
    I_TWOBODYOVERLAP_L4x2y2z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L4x2y2z_S_vrr;
    I_TWOBODYOVERLAP_L4xy3z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L4xy3z_S_vrr;
    I_TWOBODYOVERLAP_L4x4z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L4x4z_S_vrr;
    I_TWOBODYOVERLAP_L3x5y_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L3x5y_S_vrr;
    I_TWOBODYOVERLAP_L3x4yz_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L3x4yz_S_vrr;
    I_TWOBODYOVERLAP_L3x3y2z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L3x3y2z_S_vrr;
    I_TWOBODYOVERLAP_L3x2y3z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L3x2y3z_S_vrr;
    I_TWOBODYOVERLAP_L3xy4z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L3xy4z_S_vrr;
    I_TWOBODYOVERLAP_L3x5z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L3x5z_S_vrr;
    I_TWOBODYOVERLAP_L2x6y_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L2x6y_S_vrr;
    I_TWOBODYOVERLAP_L2x5yz_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L2x5yz_S_vrr;
    I_TWOBODYOVERLAP_L2x4y2z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L2x4y2z_S_vrr;
    I_TWOBODYOVERLAP_L2x3y3z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L2x3y3z_S_vrr;
    I_TWOBODYOVERLAP_L2x2y4z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L2x2y4z_S_vrr;
    I_TWOBODYOVERLAP_L2xy5z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L2xy5z_S_vrr;
    I_TWOBODYOVERLAP_L2x6z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L2x6z_S_vrr;
    I_TWOBODYOVERLAP_Lx7y_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_Lx7y_S_vrr;
    I_TWOBODYOVERLAP_Lx6yz_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_Lx6yz_S_vrr;
    I_TWOBODYOVERLAP_Lx5y2z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_Lx5y2z_S_vrr;
    I_TWOBODYOVERLAP_Lx4y3z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_Lx4y3z_S_vrr;
    I_TWOBODYOVERLAP_Lx3y4z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_Lx3y4z_S_vrr;
    I_TWOBODYOVERLAP_Lx2y5z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_Lx2y5z_S_vrr;
    I_TWOBODYOVERLAP_Lxy6z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_Lxy6z_S_vrr;
    I_TWOBODYOVERLAP_Lx7z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_Lx7z_S_vrr;
    I_TWOBODYOVERLAP_L8y_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L8y_S_vrr;
    I_TWOBODYOVERLAP_L7yz_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L7yz_S_vrr;
    I_TWOBODYOVERLAP_L6y2z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L6y2z_S_vrr;
    I_TWOBODYOVERLAP_L5y3z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L5y3z_S_vrr;
    I_TWOBODYOVERLAP_L4y4z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L4y4z_S_vrr;
    I_TWOBODYOVERLAP_L3y5z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L3y5z_S_vrr;
    I_TWOBODYOVERLAP_L2y6z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L2y6z_S_vrr;
    I_TWOBODYOVERLAP_Ly7z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_Ly7z_S_vrr;
    I_TWOBODYOVERLAP_L8z_S_aa += SQ_TWOBODYOVERLAP_L_S_aa_coefs*I_TWOBODYOVERLAP_L8z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_K_S_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_K_S_aa_coefs = alpha*alpha;
    I_TWOBODYOVERLAP_K7x_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K7x_S_vrr;
    I_TWOBODYOVERLAP_K6xy_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K6xy_S_vrr;
    I_TWOBODYOVERLAP_K6xz_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K6xz_S_vrr;
    I_TWOBODYOVERLAP_K5x2y_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K5x2y_S_vrr;
    I_TWOBODYOVERLAP_K5xyz_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K5xyz_S_vrr;
    I_TWOBODYOVERLAP_K5x2z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K5x2z_S_vrr;
    I_TWOBODYOVERLAP_K4x3y_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K4x3y_S_vrr;
    I_TWOBODYOVERLAP_K4x2yz_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K4x2yz_S_vrr;
    I_TWOBODYOVERLAP_K4xy2z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K4xy2z_S_vrr;
    I_TWOBODYOVERLAP_K4x3z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K4x3z_S_vrr;
    I_TWOBODYOVERLAP_K3x4y_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K3x4y_S_vrr;
    I_TWOBODYOVERLAP_K3x3yz_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K3x3yz_S_vrr;
    I_TWOBODYOVERLAP_K3x2y2z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K3x2y2z_S_vrr;
    I_TWOBODYOVERLAP_K3xy3z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K3xy3z_S_vrr;
    I_TWOBODYOVERLAP_K3x4z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K3x4z_S_vrr;
    I_TWOBODYOVERLAP_K2x5y_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K2x5y_S_vrr;
    I_TWOBODYOVERLAP_K2x4yz_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K2x4yz_S_vrr;
    I_TWOBODYOVERLAP_K2x3y2z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K2x3y2z_S_vrr;
    I_TWOBODYOVERLAP_K2x2y3z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K2x2y3z_S_vrr;
    I_TWOBODYOVERLAP_K2xy4z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K2xy4z_S_vrr;
    I_TWOBODYOVERLAP_K2x5z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K2x5z_S_vrr;
    I_TWOBODYOVERLAP_Kx6y_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_Kx6y_S_vrr;
    I_TWOBODYOVERLAP_Kx5yz_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_Kx5yz_S_vrr;
    I_TWOBODYOVERLAP_Kx4y2z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_Kx4y2z_S_vrr;
    I_TWOBODYOVERLAP_Kx3y3z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_Kx3y3z_S_vrr;
    I_TWOBODYOVERLAP_Kx2y4z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_Kx2y4z_S_vrr;
    I_TWOBODYOVERLAP_Kxy5z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_Kxy5z_S_vrr;
    I_TWOBODYOVERLAP_Kx6z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_Kx6z_S_vrr;
    I_TWOBODYOVERLAP_K7y_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K7y_S_vrr;
    I_TWOBODYOVERLAP_K6yz_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K6yz_S_vrr;
    I_TWOBODYOVERLAP_K5y2z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K5y2z_S_vrr;
    I_TWOBODYOVERLAP_K4y3z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K4y3z_S_vrr;
    I_TWOBODYOVERLAP_K3y4z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K3y4z_S_vrr;
    I_TWOBODYOVERLAP_K2y5z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K2y5z_S_vrr;
    I_TWOBODYOVERLAP_Ky6z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_Ky6z_S_vrr;
    I_TWOBODYOVERLAP_K7z_S_aa += SQ_TWOBODYOVERLAP_K_S_aa_coefs*I_TWOBODYOVERLAP_K7z_S_vrr;

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_I_S_aa
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_I_S_aa_coefs = alpha*alpha;
    I_TWOBODYOVERLAP_I6x_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I6x_S_vrr;
    I_TWOBODYOVERLAP_I5xy_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I5xy_S_vrr;
    I_TWOBODYOVERLAP_I5xz_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I5xz_S_vrr;
    I_TWOBODYOVERLAP_I4x2y_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I4x2y_S_vrr;
    I_TWOBODYOVERLAP_I4xyz_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I4xyz_S_vrr;
    I_TWOBODYOVERLAP_I4x2z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I4x2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3y_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I3x3y_S_vrr;
    I_TWOBODYOVERLAP_I3x2yz_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I3x2yz_S_vrr;
    I_TWOBODYOVERLAP_I3xy2z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I3xy2z_S_vrr;
    I_TWOBODYOVERLAP_I3x3z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I3x3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4y_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I2x4y_S_vrr;
    I_TWOBODYOVERLAP_I2x3yz_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I2x3yz_S_vrr;
    I_TWOBODYOVERLAP_I2x2y2z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I2x2y2z_S_vrr;
    I_TWOBODYOVERLAP_I2xy3z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I2xy3z_S_vrr;
    I_TWOBODYOVERLAP_I2x4z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I2x4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5y_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_Ix5y_S_vrr;
    I_TWOBODYOVERLAP_Ix4yz_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_Ix4yz_S_vrr;
    I_TWOBODYOVERLAP_Ix3y2z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_Ix3y2z_S_vrr;
    I_TWOBODYOVERLAP_Ix2y3z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_Ix2y3z_S_vrr;
    I_TWOBODYOVERLAP_Ixy4z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_Ixy4z_S_vrr;
    I_TWOBODYOVERLAP_Ix5z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_Ix5z_S_vrr;
    I_TWOBODYOVERLAP_I6y_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I6y_S_vrr;
    I_TWOBODYOVERLAP_I5yz_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I5yz_S_vrr;
    I_TWOBODYOVERLAP_I4y2z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I4y2z_S_vrr;
    I_TWOBODYOVERLAP_I3y3z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I3y3z_S_vrr;
    I_TWOBODYOVERLAP_I2y4z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I2y4z_S_vrr;
    I_TWOBODYOVERLAP_Iy5z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_Iy5z_S_vrr;
    I_TWOBODYOVERLAP_I6z_S_aa += SQ_TWOBODYOVERLAP_I_S_aa_coefs*I_TWOBODYOVERLAP_I6z_S_vrr;

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
     * shell quartet name: SQ_TWOBODYOVERLAP_G_S_a
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    Double SQ_TWOBODYOVERLAP_G_S_a_coefs = alpha;
    I_TWOBODYOVERLAP_G4x_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_G4x_S_vrr;
    I_TWOBODYOVERLAP_G3xy_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_G3xy_S_vrr;
    I_TWOBODYOVERLAP_G3xz_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_G3xz_S_vrr;
    I_TWOBODYOVERLAP_G2x2y_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    I_TWOBODYOVERLAP_G2xyz_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    I_TWOBODYOVERLAP_G2x2z_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3y_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    I_TWOBODYOVERLAP_Gx2yz_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    I_TWOBODYOVERLAP_Gxy2z_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    I_TWOBODYOVERLAP_Gx3z_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    I_TWOBODYOVERLAP_G4y_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_G4y_S_vrr;
    I_TWOBODYOVERLAP_G3yz_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_G3yz_S_vrr;
    I_TWOBODYOVERLAP_G2y2z_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    I_TWOBODYOVERLAP_Gy3z_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    I_TWOBODYOVERLAP_G4z_S_a += SQ_TWOBODYOVERLAP_G_S_a_coefs*I_TWOBODYOVERLAP_G4z_S_vrr;

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

    /************************************************************
     * shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    I_TWOBODYOVERLAP_D2x_S += I_TWOBODYOVERLAP_D2x_S_vrr;
    I_TWOBODYOVERLAP_Dxy_S += I_TWOBODYOVERLAP_Dxy_S_vrr;
    I_TWOBODYOVERLAP_Dxz_S += I_TWOBODYOVERLAP_Dxz_S_vrr;
    I_TWOBODYOVERLAP_D2y_S += I_TWOBODYOVERLAP_D2y_S_vrr;
    I_TWOBODYOVERLAP_Dyz_S += I_TWOBODYOVERLAP_Dyz_S_vrr;
    I_TWOBODYOVERLAP_D2z_S += I_TWOBODYOVERLAP_D2z_S_vrr;
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
   * shell quartet name: SQ_TWOBODYOVERLAP_D_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
   ************************************************************/
  Double I_TWOBODYOVERLAP_D2x_Px = I_TWOBODYOVERLAP_F3x_S+ABX*I_TWOBODYOVERLAP_D2x_S;
  Double I_TWOBODYOVERLAP_Dxy_Px = I_TWOBODYOVERLAP_F2xy_S+ABX*I_TWOBODYOVERLAP_Dxy_S;
  Double I_TWOBODYOVERLAP_Dxz_Px = I_TWOBODYOVERLAP_F2xz_S+ABX*I_TWOBODYOVERLAP_Dxz_S;
  Double I_TWOBODYOVERLAP_D2y_Px = I_TWOBODYOVERLAP_Fx2y_S+ABX*I_TWOBODYOVERLAP_D2y_S;
  Double I_TWOBODYOVERLAP_Dyz_Px = I_TWOBODYOVERLAP_Fxyz_S+ABX*I_TWOBODYOVERLAP_Dyz_S;
  Double I_TWOBODYOVERLAP_D2z_Px = I_TWOBODYOVERLAP_Fx2z_S+ABX*I_TWOBODYOVERLAP_D2z_S;
  Double I_TWOBODYOVERLAP_D2x_Py = I_TWOBODYOVERLAP_F2xy_S+ABY*I_TWOBODYOVERLAP_D2x_S;
  Double I_TWOBODYOVERLAP_Dxy_Py = I_TWOBODYOVERLAP_Fx2y_S+ABY*I_TWOBODYOVERLAP_Dxy_S;
  Double I_TWOBODYOVERLAP_Dxz_Py = I_TWOBODYOVERLAP_Fxyz_S+ABY*I_TWOBODYOVERLAP_Dxz_S;
  Double I_TWOBODYOVERLAP_D2y_Py = I_TWOBODYOVERLAP_F3y_S+ABY*I_TWOBODYOVERLAP_D2y_S;
  Double I_TWOBODYOVERLAP_Dyz_Py = I_TWOBODYOVERLAP_F2yz_S+ABY*I_TWOBODYOVERLAP_Dyz_S;
  Double I_TWOBODYOVERLAP_D2z_Py = I_TWOBODYOVERLAP_Fy2z_S+ABY*I_TWOBODYOVERLAP_D2z_S;
  Double I_TWOBODYOVERLAP_D2x_Pz = I_TWOBODYOVERLAP_F2xz_S+ABZ*I_TWOBODYOVERLAP_D2x_S;
  Double I_TWOBODYOVERLAP_Dxy_Pz = I_TWOBODYOVERLAP_Fxyz_S+ABZ*I_TWOBODYOVERLAP_Dxy_S;
  Double I_TWOBODYOVERLAP_Dxz_Pz = I_TWOBODYOVERLAP_Fx2z_S+ABZ*I_TWOBODYOVERLAP_Dxz_S;
  Double I_TWOBODYOVERLAP_D2y_Pz = I_TWOBODYOVERLAP_F2yz_S+ABZ*I_TWOBODYOVERLAP_D2y_S;
  Double I_TWOBODYOVERLAP_Dyz_Pz = I_TWOBODYOVERLAP_Fy2z_S+ABZ*I_TWOBODYOVERLAP_Dyz_S;
  Double I_TWOBODYOVERLAP_D2z_Pz = I_TWOBODYOVERLAP_F3z_S+ABZ*I_TWOBODYOVERLAP_D2z_S;

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
   * shell quartet name: SQ_TWOBODYOVERLAP_D_D
   * expanding position: BRA2
   * code section is: HRR
   * totally 12 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_P
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_P
   ************************************************************/
  Double I_TWOBODYOVERLAP_D2x_D2x = I_TWOBODYOVERLAP_F3x_Px+ABX*I_TWOBODYOVERLAP_D2x_Px;
  Double I_TWOBODYOVERLAP_Dxy_D2x = I_TWOBODYOVERLAP_F2xy_Px+ABX*I_TWOBODYOVERLAP_Dxy_Px;
  Double I_TWOBODYOVERLAP_Dxz_D2x = I_TWOBODYOVERLAP_F2xz_Px+ABX*I_TWOBODYOVERLAP_Dxz_Px;
  Double I_TWOBODYOVERLAP_D2y_D2x = I_TWOBODYOVERLAP_Fx2y_Px+ABX*I_TWOBODYOVERLAP_D2y_Px;
  Double I_TWOBODYOVERLAP_Dyz_D2x = I_TWOBODYOVERLAP_Fxyz_Px+ABX*I_TWOBODYOVERLAP_Dyz_Px;
  Double I_TWOBODYOVERLAP_D2z_D2x = I_TWOBODYOVERLAP_Fx2z_Px+ABX*I_TWOBODYOVERLAP_D2z_Px;
  Double I_TWOBODYOVERLAP_D2x_Dxy = I_TWOBODYOVERLAP_F2xy_Px+ABY*I_TWOBODYOVERLAP_D2x_Px;
  Double I_TWOBODYOVERLAP_Dxy_Dxy = I_TWOBODYOVERLAP_Fx2y_Px+ABY*I_TWOBODYOVERLAP_Dxy_Px;
  Double I_TWOBODYOVERLAP_Dxz_Dxy = I_TWOBODYOVERLAP_Fxyz_Px+ABY*I_TWOBODYOVERLAP_Dxz_Px;
  Double I_TWOBODYOVERLAP_D2y_Dxy = I_TWOBODYOVERLAP_F3y_Px+ABY*I_TWOBODYOVERLAP_D2y_Px;
  Double I_TWOBODYOVERLAP_Dyz_Dxy = I_TWOBODYOVERLAP_F2yz_Px+ABY*I_TWOBODYOVERLAP_Dyz_Px;
  Double I_TWOBODYOVERLAP_D2z_Dxy = I_TWOBODYOVERLAP_Fy2z_Px+ABY*I_TWOBODYOVERLAP_D2z_Px;
  Double I_TWOBODYOVERLAP_D2x_D2y = I_TWOBODYOVERLAP_F2xy_Py+ABY*I_TWOBODYOVERLAP_D2x_Py;
  Double I_TWOBODYOVERLAP_Dxy_D2y = I_TWOBODYOVERLAP_Fx2y_Py+ABY*I_TWOBODYOVERLAP_Dxy_Py;
  Double I_TWOBODYOVERLAP_Dxz_D2y = I_TWOBODYOVERLAP_Fxyz_Py+ABY*I_TWOBODYOVERLAP_Dxz_Py;
  Double I_TWOBODYOVERLAP_D2y_D2y = I_TWOBODYOVERLAP_F3y_Py+ABY*I_TWOBODYOVERLAP_D2y_Py;
  Double I_TWOBODYOVERLAP_Dyz_D2y = I_TWOBODYOVERLAP_F2yz_Py+ABY*I_TWOBODYOVERLAP_Dyz_Py;
  Double I_TWOBODYOVERLAP_D2z_D2y = I_TWOBODYOVERLAP_Fy2z_Py+ABY*I_TWOBODYOVERLAP_D2z_Py;
  Double I_TWOBODYOVERLAP_D2x_D2z = I_TWOBODYOVERLAP_F2xz_Pz+ABZ*I_TWOBODYOVERLAP_D2x_Pz;
  Double I_TWOBODYOVERLAP_Dxy_D2z = I_TWOBODYOVERLAP_Fxyz_Pz+ABZ*I_TWOBODYOVERLAP_Dxy_Pz;
  Double I_TWOBODYOVERLAP_Dxz_D2z = I_TWOBODYOVERLAP_Fx2z_Pz+ABZ*I_TWOBODYOVERLAP_Dxz_Pz;
  Double I_TWOBODYOVERLAP_D2y_D2z = I_TWOBODYOVERLAP_F2yz_Pz+ABZ*I_TWOBODYOVERLAP_D2y_Pz;
  Double I_TWOBODYOVERLAP_Dyz_D2z = I_TWOBODYOVERLAP_Fy2z_Pz+ABZ*I_TWOBODYOVERLAP_Dyz_Pz;
  Double I_TWOBODYOVERLAP_D2z_D2z = I_TWOBODYOVERLAP_F3z_Pz+ABZ*I_TWOBODYOVERLAP_D2z_Pz;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_P
   * expanding position: BRA2
   * code section is: HRR
   * totally 12 integrals are omitted 
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
  Double I_TWOBODYOVERLAP_G3yz_Px = I_TWOBODYOVERLAP_Hx3yz_S+ABX*I_TWOBODYOVERLAP_G3yz_S;
  Double I_TWOBODYOVERLAP_G2y2z_Px = I_TWOBODYOVERLAP_Hx2y2z_S+ABX*I_TWOBODYOVERLAP_G2y2z_S;
  Double I_TWOBODYOVERLAP_Gy3z_Px = I_TWOBODYOVERLAP_Hxy3z_S+ABX*I_TWOBODYOVERLAP_Gy3z_S;
  Double I_TWOBODYOVERLAP_G3xy_Py = I_TWOBODYOVERLAP_H3x2y_S+ABY*I_TWOBODYOVERLAP_G3xy_S;
  Double I_TWOBODYOVERLAP_G2x2y_Py = I_TWOBODYOVERLAP_H2x3y_S+ABY*I_TWOBODYOVERLAP_G2x2y_S;
  Double I_TWOBODYOVERLAP_G2xyz_Py = I_TWOBODYOVERLAP_H2x2yz_S+ABY*I_TWOBODYOVERLAP_G2xyz_S;
  Double I_TWOBODYOVERLAP_Gx3y_Py = I_TWOBODYOVERLAP_Hx4y_S+ABY*I_TWOBODYOVERLAP_Gx3y_S;
  Double I_TWOBODYOVERLAP_Gx2yz_Py = I_TWOBODYOVERLAP_Hx3yz_S+ABY*I_TWOBODYOVERLAP_Gx2yz_S;
  Double I_TWOBODYOVERLAP_Gxy2z_Py = I_TWOBODYOVERLAP_Hx2y2z_S+ABY*I_TWOBODYOVERLAP_Gxy2z_S;
  Double I_TWOBODYOVERLAP_G4y_Py = I_TWOBODYOVERLAP_H5y_S+ABY*I_TWOBODYOVERLAP_G4y_S;
  Double I_TWOBODYOVERLAP_G3yz_Py = I_TWOBODYOVERLAP_H4yz_S+ABY*I_TWOBODYOVERLAP_G3yz_S;
  Double I_TWOBODYOVERLAP_G2y2z_Py = I_TWOBODYOVERLAP_H3y2z_S+ABY*I_TWOBODYOVERLAP_G2y2z_S;
  Double I_TWOBODYOVERLAP_Gy3z_Py = I_TWOBODYOVERLAP_H2y3z_S+ABY*I_TWOBODYOVERLAP_Gy3z_S;
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
   * totally 24 integrals are omitted 
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
  Double I_TWOBODYOVERLAP_F2xz_Dxy = I_TWOBODYOVERLAP_G2xyz_Px+ABY*I_TWOBODYOVERLAP_F2xz_Px;
  Double I_TWOBODYOVERLAP_Fxyz_Dxy = I_TWOBODYOVERLAP_Gx2yz_Px+ABY*I_TWOBODYOVERLAP_Fxyz_Px;
  Double I_TWOBODYOVERLAP_Fx2z_Dxy = I_TWOBODYOVERLAP_Gxy2z_Px+ABY*I_TWOBODYOVERLAP_Fx2z_Px;
  Double I_TWOBODYOVERLAP_F2yz_Dxy = I_TWOBODYOVERLAP_G3yz_Px+ABY*I_TWOBODYOVERLAP_F2yz_Px;
  Double I_TWOBODYOVERLAP_Fy2z_Dxy = I_TWOBODYOVERLAP_G2y2z_Px+ABY*I_TWOBODYOVERLAP_Fy2z_Px;
  Double I_TWOBODYOVERLAP_F3z_Dxy = I_TWOBODYOVERLAP_Gy3z_Px+ABY*I_TWOBODYOVERLAP_F3z_Px;
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
   * shell quartet name: SQ_TWOBODYOVERLAP_D_F
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_D
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_D
   ************************************************************/
  Double I_TWOBODYOVERLAP_D2x_F3x = I_TWOBODYOVERLAP_F3x_D2x+ABX*I_TWOBODYOVERLAP_D2x_D2x;
  Double I_TWOBODYOVERLAP_Dxy_F3x = I_TWOBODYOVERLAP_F2xy_D2x+ABX*I_TWOBODYOVERLAP_Dxy_D2x;
  Double I_TWOBODYOVERLAP_Dxz_F3x = I_TWOBODYOVERLAP_F2xz_D2x+ABX*I_TWOBODYOVERLAP_Dxz_D2x;
  Double I_TWOBODYOVERLAP_D2y_F3x = I_TWOBODYOVERLAP_Fx2y_D2x+ABX*I_TWOBODYOVERLAP_D2y_D2x;
  Double I_TWOBODYOVERLAP_Dyz_F3x = I_TWOBODYOVERLAP_Fxyz_D2x+ABX*I_TWOBODYOVERLAP_Dyz_D2x;
  Double I_TWOBODYOVERLAP_D2z_F3x = I_TWOBODYOVERLAP_Fx2z_D2x+ABX*I_TWOBODYOVERLAP_D2z_D2x;
  Double I_TWOBODYOVERLAP_D2x_F2xy = I_TWOBODYOVERLAP_F2xy_D2x+ABY*I_TWOBODYOVERLAP_D2x_D2x;
  Double I_TWOBODYOVERLAP_Dxy_F2xy = I_TWOBODYOVERLAP_Fx2y_D2x+ABY*I_TWOBODYOVERLAP_Dxy_D2x;
  Double I_TWOBODYOVERLAP_Dxz_F2xy = I_TWOBODYOVERLAP_Fxyz_D2x+ABY*I_TWOBODYOVERLAP_Dxz_D2x;
  Double I_TWOBODYOVERLAP_D2y_F2xy = I_TWOBODYOVERLAP_F3y_D2x+ABY*I_TWOBODYOVERLAP_D2y_D2x;
  Double I_TWOBODYOVERLAP_Dyz_F2xy = I_TWOBODYOVERLAP_F2yz_D2x+ABY*I_TWOBODYOVERLAP_Dyz_D2x;
  Double I_TWOBODYOVERLAP_D2z_F2xy = I_TWOBODYOVERLAP_Fy2z_D2x+ABY*I_TWOBODYOVERLAP_D2z_D2x;
  Double I_TWOBODYOVERLAP_D2x_F2xz = I_TWOBODYOVERLAP_F2xz_D2x+ABZ*I_TWOBODYOVERLAP_D2x_D2x;
  Double I_TWOBODYOVERLAP_Dxy_F2xz = I_TWOBODYOVERLAP_Fxyz_D2x+ABZ*I_TWOBODYOVERLAP_Dxy_D2x;
  Double I_TWOBODYOVERLAP_Dxz_F2xz = I_TWOBODYOVERLAP_Fx2z_D2x+ABZ*I_TWOBODYOVERLAP_Dxz_D2x;
  Double I_TWOBODYOVERLAP_D2y_F2xz = I_TWOBODYOVERLAP_F2yz_D2x+ABZ*I_TWOBODYOVERLAP_D2y_D2x;
  Double I_TWOBODYOVERLAP_Dyz_F2xz = I_TWOBODYOVERLAP_Fy2z_D2x+ABZ*I_TWOBODYOVERLAP_Dyz_D2x;
  Double I_TWOBODYOVERLAP_D2z_F2xz = I_TWOBODYOVERLAP_F3z_D2x+ABZ*I_TWOBODYOVERLAP_D2z_D2x;
  Double I_TWOBODYOVERLAP_D2x_Fx2y = I_TWOBODYOVERLAP_F3x_D2y+ABX*I_TWOBODYOVERLAP_D2x_D2y;
  Double I_TWOBODYOVERLAP_Dxy_Fx2y = I_TWOBODYOVERLAP_F2xy_D2y+ABX*I_TWOBODYOVERLAP_Dxy_D2y;
  Double I_TWOBODYOVERLAP_Dxz_Fx2y = I_TWOBODYOVERLAP_F2xz_D2y+ABX*I_TWOBODYOVERLAP_Dxz_D2y;
  Double I_TWOBODYOVERLAP_D2y_Fx2y = I_TWOBODYOVERLAP_Fx2y_D2y+ABX*I_TWOBODYOVERLAP_D2y_D2y;
  Double I_TWOBODYOVERLAP_Dyz_Fx2y = I_TWOBODYOVERLAP_Fxyz_D2y+ABX*I_TWOBODYOVERLAP_Dyz_D2y;
  Double I_TWOBODYOVERLAP_D2z_Fx2y = I_TWOBODYOVERLAP_Fx2z_D2y+ABX*I_TWOBODYOVERLAP_D2z_D2y;
  Double I_TWOBODYOVERLAP_D2x_Fxyz = I_TWOBODYOVERLAP_F2xz_Dxy+ABZ*I_TWOBODYOVERLAP_D2x_Dxy;
  Double I_TWOBODYOVERLAP_Dxy_Fxyz = I_TWOBODYOVERLAP_Fxyz_Dxy+ABZ*I_TWOBODYOVERLAP_Dxy_Dxy;
  Double I_TWOBODYOVERLAP_Dxz_Fxyz = I_TWOBODYOVERLAP_Fx2z_Dxy+ABZ*I_TWOBODYOVERLAP_Dxz_Dxy;
  Double I_TWOBODYOVERLAP_D2y_Fxyz = I_TWOBODYOVERLAP_F2yz_Dxy+ABZ*I_TWOBODYOVERLAP_D2y_Dxy;
  Double I_TWOBODYOVERLAP_Dyz_Fxyz = I_TWOBODYOVERLAP_Fy2z_Dxy+ABZ*I_TWOBODYOVERLAP_Dyz_Dxy;
  Double I_TWOBODYOVERLAP_D2z_Fxyz = I_TWOBODYOVERLAP_F3z_Dxy+ABZ*I_TWOBODYOVERLAP_D2z_Dxy;
  Double I_TWOBODYOVERLAP_D2x_Fx2z = I_TWOBODYOVERLAP_F3x_D2z+ABX*I_TWOBODYOVERLAP_D2x_D2z;
  Double I_TWOBODYOVERLAP_Dxy_Fx2z = I_TWOBODYOVERLAP_F2xy_D2z+ABX*I_TWOBODYOVERLAP_Dxy_D2z;
  Double I_TWOBODYOVERLAP_Dxz_Fx2z = I_TWOBODYOVERLAP_F2xz_D2z+ABX*I_TWOBODYOVERLAP_Dxz_D2z;
  Double I_TWOBODYOVERLAP_D2y_Fx2z = I_TWOBODYOVERLAP_Fx2y_D2z+ABX*I_TWOBODYOVERLAP_D2y_D2z;
  Double I_TWOBODYOVERLAP_Dyz_Fx2z = I_TWOBODYOVERLAP_Fxyz_D2z+ABX*I_TWOBODYOVERLAP_Dyz_D2z;
  Double I_TWOBODYOVERLAP_D2z_Fx2z = I_TWOBODYOVERLAP_Fx2z_D2z+ABX*I_TWOBODYOVERLAP_D2z_D2z;
  Double I_TWOBODYOVERLAP_D2x_F3y = I_TWOBODYOVERLAP_F2xy_D2y+ABY*I_TWOBODYOVERLAP_D2x_D2y;
  Double I_TWOBODYOVERLAP_Dxy_F3y = I_TWOBODYOVERLAP_Fx2y_D2y+ABY*I_TWOBODYOVERLAP_Dxy_D2y;
  Double I_TWOBODYOVERLAP_Dxz_F3y = I_TWOBODYOVERLAP_Fxyz_D2y+ABY*I_TWOBODYOVERLAP_Dxz_D2y;
  Double I_TWOBODYOVERLAP_D2y_F3y = I_TWOBODYOVERLAP_F3y_D2y+ABY*I_TWOBODYOVERLAP_D2y_D2y;
  Double I_TWOBODYOVERLAP_Dyz_F3y = I_TWOBODYOVERLAP_F2yz_D2y+ABY*I_TWOBODYOVERLAP_Dyz_D2y;
  Double I_TWOBODYOVERLAP_D2z_F3y = I_TWOBODYOVERLAP_Fy2z_D2y+ABY*I_TWOBODYOVERLAP_D2z_D2y;
  Double I_TWOBODYOVERLAP_D2x_F2yz = I_TWOBODYOVERLAP_F2xz_D2y+ABZ*I_TWOBODYOVERLAP_D2x_D2y;
  Double I_TWOBODYOVERLAP_Dxy_F2yz = I_TWOBODYOVERLAP_Fxyz_D2y+ABZ*I_TWOBODYOVERLAP_Dxy_D2y;
  Double I_TWOBODYOVERLAP_Dxz_F2yz = I_TWOBODYOVERLAP_Fx2z_D2y+ABZ*I_TWOBODYOVERLAP_Dxz_D2y;
  Double I_TWOBODYOVERLAP_D2y_F2yz = I_TWOBODYOVERLAP_F2yz_D2y+ABZ*I_TWOBODYOVERLAP_D2y_D2y;
  Double I_TWOBODYOVERLAP_Dyz_F2yz = I_TWOBODYOVERLAP_Fy2z_D2y+ABZ*I_TWOBODYOVERLAP_Dyz_D2y;
  Double I_TWOBODYOVERLAP_D2z_F2yz = I_TWOBODYOVERLAP_F3z_D2y+ABZ*I_TWOBODYOVERLAP_D2z_D2y;
  Double I_TWOBODYOVERLAP_D2x_Fy2z = I_TWOBODYOVERLAP_F2xy_D2z+ABY*I_TWOBODYOVERLAP_D2x_D2z;
  Double I_TWOBODYOVERLAP_Dxy_Fy2z = I_TWOBODYOVERLAP_Fx2y_D2z+ABY*I_TWOBODYOVERLAP_Dxy_D2z;
  Double I_TWOBODYOVERLAP_Dxz_Fy2z = I_TWOBODYOVERLAP_Fxyz_D2z+ABY*I_TWOBODYOVERLAP_Dxz_D2z;
  Double I_TWOBODYOVERLAP_D2y_Fy2z = I_TWOBODYOVERLAP_F3y_D2z+ABY*I_TWOBODYOVERLAP_D2y_D2z;
  Double I_TWOBODYOVERLAP_Dyz_Fy2z = I_TWOBODYOVERLAP_F2yz_D2z+ABY*I_TWOBODYOVERLAP_Dyz_D2z;
  Double I_TWOBODYOVERLAP_D2z_Fy2z = I_TWOBODYOVERLAP_Fy2z_D2z+ABY*I_TWOBODYOVERLAP_D2z_D2z;
  Double I_TWOBODYOVERLAP_D2x_F3z = I_TWOBODYOVERLAP_F2xz_D2z+ABZ*I_TWOBODYOVERLAP_D2x_D2z;
  Double I_TWOBODYOVERLAP_Dxy_F3z = I_TWOBODYOVERLAP_Fxyz_D2z+ABZ*I_TWOBODYOVERLAP_Dxy_D2z;
  Double I_TWOBODYOVERLAP_Dxz_F3z = I_TWOBODYOVERLAP_Fx2z_D2z+ABZ*I_TWOBODYOVERLAP_Dxz_D2z;
  Double I_TWOBODYOVERLAP_D2y_F3z = I_TWOBODYOVERLAP_F2yz_D2z+ABZ*I_TWOBODYOVERLAP_D2y_D2z;
  Double I_TWOBODYOVERLAP_Dyz_F3z = I_TWOBODYOVERLAP_Fy2z_D2z+ABZ*I_TWOBODYOVERLAP_Dyz_D2z;
  Double I_TWOBODYOVERLAP_D2z_F3z = I_TWOBODYOVERLAP_F3z_D2z+ABZ*I_TWOBODYOVERLAP_D2z_D2z;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_G4x_Px_a = I_TWOBODYOVERLAP_H5x_S_a+ABX*I_TWOBODYOVERLAP_G4x_S_a;
  Double I_TWOBODYOVERLAP_G3xy_Px_a = I_TWOBODYOVERLAP_H4xy_S_a+ABX*I_TWOBODYOVERLAP_G3xy_S_a;
  Double I_TWOBODYOVERLAP_G3xz_Px_a = I_TWOBODYOVERLAP_H4xz_S_a+ABX*I_TWOBODYOVERLAP_G3xz_S_a;
  Double I_TWOBODYOVERLAP_G2x2y_Px_a = I_TWOBODYOVERLAP_H3x2y_S_a+ABX*I_TWOBODYOVERLAP_G2x2y_S_a;
  Double I_TWOBODYOVERLAP_G2xyz_Px_a = I_TWOBODYOVERLAP_H3xyz_S_a+ABX*I_TWOBODYOVERLAP_G2xyz_S_a;
  Double I_TWOBODYOVERLAP_G2x2z_Px_a = I_TWOBODYOVERLAP_H3x2z_S_a+ABX*I_TWOBODYOVERLAP_G2x2z_S_a;
  Double I_TWOBODYOVERLAP_Gx3y_Px_a = I_TWOBODYOVERLAP_H2x3y_S_a+ABX*I_TWOBODYOVERLAP_Gx3y_S_a;
  Double I_TWOBODYOVERLAP_Gx2yz_Px_a = I_TWOBODYOVERLAP_H2x2yz_S_a+ABX*I_TWOBODYOVERLAP_Gx2yz_S_a;
  Double I_TWOBODYOVERLAP_Gxy2z_Px_a = I_TWOBODYOVERLAP_H2xy2z_S_a+ABX*I_TWOBODYOVERLAP_Gxy2z_S_a;
  Double I_TWOBODYOVERLAP_Gx3z_Px_a = I_TWOBODYOVERLAP_H2x3z_S_a+ABX*I_TWOBODYOVERLAP_Gx3z_S_a;
  Double I_TWOBODYOVERLAP_G4y_Px_a = I_TWOBODYOVERLAP_Hx4y_S_a+ABX*I_TWOBODYOVERLAP_G4y_S_a;
  Double I_TWOBODYOVERLAP_G3yz_Px_a = I_TWOBODYOVERLAP_Hx3yz_S_a+ABX*I_TWOBODYOVERLAP_G3yz_S_a;
  Double I_TWOBODYOVERLAP_G2y2z_Px_a = I_TWOBODYOVERLAP_Hx2y2z_S_a+ABX*I_TWOBODYOVERLAP_G2y2z_S_a;
  Double I_TWOBODYOVERLAP_Gy3z_Px_a = I_TWOBODYOVERLAP_Hxy3z_S_a+ABX*I_TWOBODYOVERLAP_Gy3z_S_a;
  Double I_TWOBODYOVERLAP_G4z_Px_a = I_TWOBODYOVERLAP_Hx4z_S_a+ABX*I_TWOBODYOVERLAP_G4z_S_a;
  Double I_TWOBODYOVERLAP_G4x_Py_a = I_TWOBODYOVERLAP_H4xy_S_a+ABY*I_TWOBODYOVERLAP_G4x_S_a;
  Double I_TWOBODYOVERLAP_G3xy_Py_a = I_TWOBODYOVERLAP_H3x2y_S_a+ABY*I_TWOBODYOVERLAP_G3xy_S_a;
  Double I_TWOBODYOVERLAP_G3xz_Py_a = I_TWOBODYOVERLAP_H3xyz_S_a+ABY*I_TWOBODYOVERLAP_G3xz_S_a;
  Double I_TWOBODYOVERLAP_G2x2y_Py_a = I_TWOBODYOVERLAP_H2x3y_S_a+ABY*I_TWOBODYOVERLAP_G2x2y_S_a;
  Double I_TWOBODYOVERLAP_G2xyz_Py_a = I_TWOBODYOVERLAP_H2x2yz_S_a+ABY*I_TWOBODYOVERLAP_G2xyz_S_a;
  Double I_TWOBODYOVERLAP_G2x2z_Py_a = I_TWOBODYOVERLAP_H2xy2z_S_a+ABY*I_TWOBODYOVERLAP_G2x2z_S_a;
  Double I_TWOBODYOVERLAP_Gx3y_Py_a = I_TWOBODYOVERLAP_Hx4y_S_a+ABY*I_TWOBODYOVERLAP_Gx3y_S_a;
  Double I_TWOBODYOVERLAP_Gx2yz_Py_a = I_TWOBODYOVERLAP_Hx3yz_S_a+ABY*I_TWOBODYOVERLAP_Gx2yz_S_a;
  Double I_TWOBODYOVERLAP_Gxy2z_Py_a = I_TWOBODYOVERLAP_Hx2y2z_S_a+ABY*I_TWOBODYOVERLAP_Gxy2z_S_a;
  Double I_TWOBODYOVERLAP_Gx3z_Py_a = I_TWOBODYOVERLAP_Hxy3z_S_a+ABY*I_TWOBODYOVERLAP_Gx3z_S_a;
  Double I_TWOBODYOVERLAP_G4y_Py_a = I_TWOBODYOVERLAP_H5y_S_a+ABY*I_TWOBODYOVERLAP_G4y_S_a;
  Double I_TWOBODYOVERLAP_G3yz_Py_a = I_TWOBODYOVERLAP_H4yz_S_a+ABY*I_TWOBODYOVERLAP_G3yz_S_a;
  Double I_TWOBODYOVERLAP_G2y2z_Py_a = I_TWOBODYOVERLAP_H3y2z_S_a+ABY*I_TWOBODYOVERLAP_G2y2z_S_a;
  Double I_TWOBODYOVERLAP_Gy3z_Py_a = I_TWOBODYOVERLAP_H2y3z_S_a+ABY*I_TWOBODYOVERLAP_Gy3z_S_a;
  Double I_TWOBODYOVERLAP_G4z_Py_a = I_TWOBODYOVERLAP_Hy4z_S_a+ABY*I_TWOBODYOVERLAP_G4z_S_a;
  Double I_TWOBODYOVERLAP_G4x_Pz_a = I_TWOBODYOVERLAP_H4xz_S_a+ABZ*I_TWOBODYOVERLAP_G4x_S_a;
  Double I_TWOBODYOVERLAP_G3xy_Pz_a = I_TWOBODYOVERLAP_H3xyz_S_a+ABZ*I_TWOBODYOVERLAP_G3xy_S_a;
  Double I_TWOBODYOVERLAP_G3xz_Pz_a = I_TWOBODYOVERLAP_H3x2z_S_a+ABZ*I_TWOBODYOVERLAP_G3xz_S_a;
  Double I_TWOBODYOVERLAP_G2x2y_Pz_a = I_TWOBODYOVERLAP_H2x2yz_S_a+ABZ*I_TWOBODYOVERLAP_G2x2y_S_a;
  Double I_TWOBODYOVERLAP_G2xyz_Pz_a = I_TWOBODYOVERLAP_H2xy2z_S_a+ABZ*I_TWOBODYOVERLAP_G2xyz_S_a;
  Double I_TWOBODYOVERLAP_G2x2z_Pz_a = I_TWOBODYOVERLAP_H2x3z_S_a+ABZ*I_TWOBODYOVERLAP_G2x2z_S_a;
  Double I_TWOBODYOVERLAP_Gx3y_Pz_a = I_TWOBODYOVERLAP_Hx3yz_S_a+ABZ*I_TWOBODYOVERLAP_Gx3y_S_a;
  Double I_TWOBODYOVERLAP_Gx2yz_Pz_a = I_TWOBODYOVERLAP_Hx2y2z_S_a+ABZ*I_TWOBODYOVERLAP_Gx2yz_S_a;
  Double I_TWOBODYOVERLAP_Gxy2z_Pz_a = I_TWOBODYOVERLAP_Hxy3z_S_a+ABZ*I_TWOBODYOVERLAP_Gxy2z_S_a;
  Double I_TWOBODYOVERLAP_Gx3z_Pz_a = I_TWOBODYOVERLAP_Hx4z_S_a+ABZ*I_TWOBODYOVERLAP_Gx3z_S_a;
  Double I_TWOBODYOVERLAP_G4y_Pz_a = I_TWOBODYOVERLAP_H4yz_S_a+ABZ*I_TWOBODYOVERLAP_G4y_S_a;
  Double I_TWOBODYOVERLAP_G3yz_Pz_a = I_TWOBODYOVERLAP_H3y2z_S_a+ABZ*I_TWOBODYOVERLAP_G3yz_S_a;
  Double I_TWOBODYOVERLAP_G2y2z_Pz_a = I_TWOBODYOVERLAP_H2y3z_S_a+ABZ*I_TWOBODYOVERLAP_G2y2z_S_a;
  Double I_TWOBODYOVERLAP_Gy3z_Pz_a = I_TWOBODYOVERLAP_Hy4z_S_a+ABZ*I_TWOBODYOVERLAP_Gy3z_S_a;
  Double I_TWOBODYOVERLAP_G4z_Pz_a = I_TWOBODYOVERLAP_H5z_S_a+ABZ*I_TWOBODYOVERLAP_G4z_S_a;

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
   * shell quartet name: SQ_TWOBODYOVERLAP_G_D_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 30 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_P_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_G4x_D2x_a = I_TWOBODYOVERLAP_H5x_Px_a+ABX*I_TWOBODYOVERLAP_G4x_Px_a;
  Double I_TWOBODYOVERLAP_G3xy_D2x_a = I_TWOBODYOVERLAP_H4xy_Px_a+ABX*I_TWOBODYOVERLAP_G3xy_Px_a;
  Double I_TWOBODYOVERLAP_G3xz_D2x_a = I_TWOBODYOVERLAP_H4xz_Px_a+ABX*I_TWOBODYOVERLAP_G3xz_Px_a;
  Double I_TWOBODYOVERLAP_G2x2y_D2x_a = I_TWOBODYOVERLAP_H3x2y_Px_a+ABX*I_TWOBODYOVERLAP_G2x2y_Px_a;
  Double I_TWOBODYOVERLAP_G2xyz_D2x_a = I_TWOBODYOVERLAP_H3xyz_Px_a+ABX*I_TWOBODYOVERLAP_G2xyz_Px_a;
  Double I_TWOBODYOVERLAP_G2x2z_D2x_a = I_TWOBODYOVERLAP_H3x2z_Px_a+ABX*I_TWOBODYOVERLAP_G2x2z_Px_a;
  Double I_TWOBODYOVERLAP_Gx3y_D2x_a = I_TWOBODYOVERLAP_H2x3y_Px_a+ABX*I_TWOBODYOVERLAP_Gx3y_Px_a;
  Double I_TWOBODYOVERLAP_Gx2yz_D2x_a = I_TWOBODYOVERLAP_H2x2yz_Px_a+ABX*I_TWOBODYOVERLAP_Gx2yz_Px_a;
  Double I_TWOBODYOVERLAP_Gxy2z_D2x_a = I_TWOBODYOVERLAP_H2xy2z_Px_a+ABX*I_TWOBODYOVERLAP_Gxy2z_Px_a;
  Double I_TWOBODYOVERLAP_Gx3z_D2x_a = I_TWOBODYOVERLAP_H2x3z_Px_a+ABX*I_TWOBODYOVERLAP_Gx3z_Px_a;
  Double I_TWOBODYOVERLAP_G4y_D2x_a = I_TWOBODYOVERLAP_Hx4y_Px_a+ABX*I_TWOBODYOVERLAP_G4y_Px_a;
  Double I_TWOBODYOVERLAP_G3yz_D2x_a = I_TWOBODYOVERLAP_Hx3yz_Px_a+ABX*I_TWOBODYOVERLAP_G3yz_Px_a;
  Double I_TWOBODYOVERLAP_G2y2z_D2x_a = I_TWOBODYOVERLAP_Hx2y2z_Px_a+ABX*I_TWOBODYOVERLAP_G2y2z_Px_a;
  Double I_TWOBODYOVERLAP_Gy3z_D2x_a = I_TWOBODYOVERLAP_Hxy3z_Px_a+ABX*I_TWOBODYOVERLAP_Gy3z_Px_a;
  Double I_TWOBODYOVERLAP_G4z_D2x_a = I_TWOBODYOVERLAP_Hx4z_Px_a+ABX*I_TWOBODYOVERLAP_G4z_Px_a;
  Double I_TWOBODYOVERLAP_G4x_Dxy_a = I_TWOBODYOVERLAP_H4xy_Px_a+ABY*I_TWOBODYOVERLAP_G4x_Px_a;
  Double I_TWOBODYOVERLAP_G3xy_Dxy_a = I_TWOBODYOVERLAP_H3x2y_Px_a+ABY*I_TWOBODYOVERLAP_G3xy_Px_a;
  Double I_TWOBODYOVERLAP_G3xz_Dxy_a = I_TWOBODYOVERLAP_H3xyz_Px_a+ABY*I_TWOBODYOVERLAP_G3xz_Px_a;
  Double I_TWOBODYOVERLAP_G2x2y_Dxy_a = I_TWOBODYOVERLAP_H2x3y_Px_a+ABY*I_TWOBODYOVERLAP_G2x2y_Px_a;
  Double I_TWOBODYOVERLAP_G2xyz_Dxy_a = I_TWOBODYOVERLAP_H2x2yz_Px_a+ABY*I_TWOBODYOVERLAP_G2xyz_Px_a;
  Double I_TWOBODYOVERLAP_G2x2z_Dxy_a = I_TWOBODYOVERLAP_H2xy2z_Px_a+ABY*I_TWOBODYOVERLAP_G2x2z_Px_a;
  Double I_TWOBODYOVERLAP_Gx3y_Dxy_a = I_TWOBODYOVERLAP_Hx4y_Px_a+ABY*I_TWOBODYOVERLAP_Gx3y_Px_a;
  Double I_TWOBODYOVERLAP_Gx2yz_Dxy_a = I_TWOBODYOVERLAP_Hx3yz_Px_a+ABY*I_TWOBODYOVERLAP_Gx2yz_Px_a;
  Double I_TWOBODYOVERLAP_Gxy2z_Dxy_a = I_TWOBODYOVERLAP_Hx2y2z_Px_a+ABY*I_TWOBODYOVERLAP_Gxy2z_Px_a;
  Double I_TWOBODYOVERLAP_Gx3z_Dxy_a = I_TWOBODYOVERLAP_Hxy3z_Px_a+ABY*I_TWOBODYOVERLAP_Gx3z_Px_a;
  Double I_TWOBODYOVERLAP_G4y_Dxy_a = I_TWOBODYOVERLAP_H5y_Px_a+ABY*I_TWOBODYOVERLAP_G4y_Px_a;
  Double I_TWOBODYOVERLAP_G3yz_Dxy_a = I_TWOBODYOVERLAP_H4yz_Px_a+ABY*I_TWOBODYOVERLAP_G3yz_Px_a;
  Double I_TWOBODYOVERLAP_G2y2z_Dxy_a = I_TWOBODYOVERLAP_H3y2z_Px_a+ABY*I_TWOBODYOVERLAP_G2y2z_Px_a;
  Double I_TWOBODYOVERLAP_Gy3z_Dxy_a = I_TWOBODYOVERLAP_H2y3z_Px_a+ABY*I_TWOBODYOVERLAP_Gy3z_Px_a;
  Double I_TWOBODYOVERLAP_G4z_Dxy_a = I_TWOBODYOVERLAP_Hy4z_Px_a+ABY*I_TWOBODYOVERLAP_G4z_Px_a;
  Double I_TWOBODYOVERLAP_G4x_D2y_a = I_TWOBODYOVERLAP_H4xy_Py_a+ABY*I_TWOBODYOVERLAP_G4x_Py_a;
  Double I_TWOBODYOVERLAP_G3xy_D2y_a = I_TWOBODYOVERLAP_H3x2y_Py_a+ABY*I_TWOBODYOVERLAP_G3xy_Py_a;
  Double I_TWOBODYOVERLAP_G3xz_D2y_a = I_TWOBODYOVERLAP_H3xyz_Py_a+ABY*I_TWOBODYOVERLAP_G3xz_Py_a;
  Double I_TWOBODYOVERLAP_G2x2y_D2y_a = I_TWOBODYOVERLAP_H2x3y_Py_a+ABY*I_TWOBODYOVERLAP_G2x2y_Py_a;
  Double I_TWOBODYOVERLAP_G2xyz_D2y_a = I_TWOBODYOVERLAP_H2x2yz_Py_a+ABY*I_TWOBODYOVERLAP_G2xyz_Py_a;
  Double I_TWOBODYOVERLAP_G2x2z_D2y_a = I_TWOBODYOVERLAP_H2xy2z_Py_a+ABY*I_TWOBODYOVERLAP_G2x2z_Py_a;
  Double I_TWOBODYOVERLAP_Gx3y_D2y_a = I_TWOBODYOVERLAP_Hx4y_Py_a+ABY*I_TWOBODYOVERLAP_Gx3y_Py_a;
  Double I_TWOBODYOVERLAP_Gx2yz_D2y_a = I_TWOBODYOVERLAP_Hx3yz_Py_a+ABY*I_TWOBODYOVERLAP_Gx2yz_Py_a;
  Double I_TWOBODYOVERLAP_Gxy2z_D2y_a = I_TWOBODYOVERLAP_Hx2y2z_Py_a+ABY*I_TWOBODYOVERLAP_Gxy2z_Py_a;
  Double I_TWOBODYOVERLAP_Gx3z_D2y_a = I_TWOBODYOVERLAP_Hxy3z_Py_a+ABY*I_TWOBODYOVERLAP_Gx3z_Py_a;
  Double I_TWOBODYOVERLAP_G4y_D2y_a = I_TWOBODYOVERLAP_H5y_Py_a+ABY*I_TWOBODYOVERLAP_G4y_Py_a;
  Double I_TWOBODYOVERLAP_G3yz_D2y_a = I_TWOBODYOVERLAP_H4yz_Py_a+ABY*I_TWOBODYOVERLAP_G3yz_Py_a;
  Double I_TWOBODYOVERLAP_G2y2z_D2y_a = I_TWOBODYOVERLAP_H3y2z_Py_a+ABY*I_TWOBODYOVERLAP_G2y2z_Py_a;
  Double I_TWOBODYOVERLAP_Gy3z_D2y_a = I_TWOBODYOVERLAP_H2y3z_Py_a+ABY*I_TWOBODYOVERLAP_Gy3z_Py_a;
  Double I_TWOBODYOVERLAP_G4z_D2y_a = I_TWOBODYOVERLAP_Hy4z_Py_a+ABY*I_TWOBODYOVERLAP_G4z_Py_a;
  Double I_TWOBODYOVERLAP_G4x_D2z_a = I_TWOBODYOVERLAP_H4xz_Pz_a+ABZ*I_TWOBODYOVERLAP_G4x_Pz_a;
  Double I_TWOBODYOVERLAP_G3xy_D2z_a = I_TWOBODYOVERLAP_H3xyz_Pz_a+ABZ*I_TWOBODYOVERLAP_G3xy_Pz_a;
  Double I_TWOBODYOVERLAP_G3xz_D2z_a = I_TWOBODYOVERLAP_H3x2z_Pz_a+ABZ*I_TWOBODYOVERLAP_G3xz_Pz_a;
  Double I_TWOBODYOVERLAP_G2x2y_D2z_a = I_TWOBODYOVERLAP_H2x2yz_Pz_a+ABZ*I_TWOBODYOVERLAP_G2x2y_Pz_a;
  Double I_TWOBODYOVERLAP_G2xyz_D2z_a = I_TWOBODYOVERLAP_H2xy2z_Pz_a+ABZ*I_TWOBODYOVERLAP_G2xyz_Pz_a;
  Double I_TWOBODYOVERLAP_G2x2z_D2z_a = I_TWOBODYOVERLAP_H2x3z_Pz_a+ABZ*I_TWOBODYOVERLAP_G2x2z_Pz_a;
  Double I_TWOBODYOVERLAP_Gx3y_D2z_a = I_TWOBODYOVERLAP_Hx3yz_Pz_a+ABZ*I_TWOBODYOVERLAP_Gx3y_Pz_a;
  Double I_TWOBODYOVERLAP_Gx2yz_D2z_a = I_TWOBODYOVERLAP_Hx2y2z_Pz_a+ABZ*I_TWOBODYOVERLAP_Gx2yz_Pz_a;
  Double I_TWOBODYOVERLAP_Gxy2z_D2z_a = I_TWOBODYOVERLAP_Hxy3z_Pz_a+ABZ*I_TWOBODYOVERLAP_Gxy2z_Pz_a;
  Double I_TWOBODYOVERLAP_Gx3z_D2z_a = I_TWOBODYOVERLAP_Hx4z_Pz_a+ABZ*I_TWOBODYOVERLAP_Gx3z_Pz_a;
  Double I_TWOBODYOVERLAP_G4y_D2z_a = I_TWOBODYOVERLAP_H4yz_Pz_a+ABZ*I_TWOBODYOVERLAP_G4y_Pz_a;
  Double I_TWOBODYOVERLAP_G3yz_D2z_a = I_TWOBODYOVERLAP_H3y2z_Pz_a+ABZ*I_TWOBODYOVERLAP_G3yz_Pz_a;
  Double I_TWOBODYOVERLAP_G2y2z_D2z_a = I_TWOBODYOVERLAP_H2y3z_Pz_a+ABZ*I_TWOBODYOVERLAP_G2y2z_Pz_a;
  Double I_TWOBODYOVERLAP_Gy3z_D2z_a = I_TWOBODYOVERLAP_Hy4z_Pz_a+ABZ*I_TWOBODYOVERLAP_Gy3z_Pz_a;
  Double I_TWOBODYOVERLAP_G4z_D2z_a = I_TWOBODYOVERLAP_H5z_Pz_a+ABZ*I_TWOBODYOVERLAP_G4z_Pz_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_I_P_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 16 integrals are omitted 
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
  Double I_TWOBODYOVERLAP_I5yz_Px_a = I_TWOBODYOVERLAP_Kx5yz_S_a+ABX*I_TWOBODYOVERLAP_I5yz_S_a;
  Double I_TWOBODYOVERLAP_I4y2z_Px_a = I_TWOBODYOVERLAP_Kx4y2z_S_a+ABX*I_TWOBODYOVERLAP_I4y2z_S_a;
  Double I_TWOBODYOVERLAP_I3y3z_Px_a = I_TWOBODYOVERLAP_Kx3y3z_S_a+ABX*I_TWOBODYOVERLAP_I3y3z_S_a;
  Double I_TWOBODYOVERLAP_I2y4z_Px_a = I_TWOBODYOVERLAP_Kx2y4z_S_a+ABX*I_TWOBODYOVERLAP_I2y4z_S_a;
  Double I_TWOBODYOVERLAP_Iy5z_Px_a = I_TWOBODYOVERLAP_Kxy5z_S_a+ABX*I_TWOBODYOVERLAP_Iy5z_S_a;
  Double I_TWOBODYOVERLAP_I5xy_Py_a = I_TWOBODYOVERLAP_K5x2y_S_a+ABY*I_TWOBODYOVERLAP_I5xy_S_a;
  Double I_TWOBODYOVERLAP_I4x2y_Py_a = I_TWOBODYOVERLAP_K4x3y_S_a+ABY*I_TWOBODYOVERLAP_I4x2y_S_a;
  Double I_TWOBODYOVERLAP_I4xyz_Py_a = I_TWOBODYOVERLAP_K4x2yz_S_a+ABY*I_TWOBODYOVERLAP_I4xyz_S_a;
  Double I_TWOBODYOVERLAP_I3x3y_Py_a = I_TWOBODYOVERLAP_K3x4y_S_a+ABY*I_TWOBODYOVERLAP_I3x3y_S_a;
  Double I_TWOBODYOVERLAP_I3x2yz_Py_a = I_TWOBODYOVERLAP_K3x3yz_S_a+ABY*I_TWOBODYOVERLAP_I3x2yz_S_a;
  Double I_TWOBODYOVERLAP_I3xy2z_Py_a = I_TWOBODYOVERLAP_K3x2y2z_S_a+ABY*I_TWOBODYOVERLAP_I3xy2z_S_a;
  Double I_TWOBODYOVERLAP_I2x4y_Py_a = I_TWOBODYOVERLAP_K2x5y_S_a+ABY*I_TWOBODYOVERLAP_I2x4y_S_a;
  Double I_TWOBODYOVERLAP_I2x3yz_Py_a = I_TWOBODYOVERLAP_K2x4yz_S_a+ABY*I_TWOBODYOVERLAP_I2x3yz_S_a;
  Double I_TWOBODYOVERLAP_I2x2y2z_Py_a = I_TWOBODYOVERLAP_K2x3y2z_S_a+ABY*I_TWOBODYOVERLAP_I2x2y2z_S_a;
  Double I_TWOBODYOVERLAP_I2xy3z_Py_a = I_TWOBODYOVERLAP_K2x2y3z_S_a+ABY*I_TWOBODYOVERLAP_I2xy3z_S_a;
  Double I_TWOBODYOVERLAP_Ix5y_Py_a = I_TWOBODYOVERLAP_Kx6y_S_a+ABY*I_TWOBODYOVERLAP_Ix5y_S_a;
  Double I_TWOBODYOVERLAP_Ix4yz_Py_a = I_TWOBODYOVERLAP_Kx5yz_S_a+ABY*I_TWOBODYOVERLAP_Ix4yz_S_a;
  Double I_TWOBODYOVERLAP_Ix3y2z_Py_a = I_TWOBODYOVERLAP_Kx4y2z_S_a+ABY*I_TWOBODYOVERLAP_Ix3y2z_S_a;
  Double I_TWOBODYOVERLAP_Ix2y3z_Py_a = I_TWOBODYOVERLAP_Kx3y3z_S_a+ABY*I_TWOBODYOVERLAP_Ix2y3z_S_a;
  Double I_TWOBODYOVERLAP_Ixy4z_Py_a = I_TWOBODYOVERLAP_Kx2y4z_S_a+ABY*I_TWOBODYOVERLAP_Ixy4z_S_a;
  Double I_TWOBODYOVERLAP_I6y_Py_a = I_TWOBODYOVERLAP_K7y_S_a+ABY*I_TWOBODYOVERLAP_I6y_S_a;
  Double I_TWOBODYOVERLAP_I5yz_Py_a = I_TWOBODYOVERLAP_K6yz_S_a+ABY*I_TWOBODYOVERLAP_I5yz_S_a;
  Double I_TWOBODYOVERLAP_I4y2z_Py_a = I_TWOBODYOVERLAP_K5y2z_S_a+ABY*I_TWOBODYOVERLAP_I4y2z_S_a;
  Double I_TWOBODYOVERLAP_I3y3z_Py_a = I_TWOBODYOVERLAP_K4y3z_S_a+ABY*I_TWOBODYOVERLAP_I3y3z_S_a;
  Double I_TWOBODYOVERLAP_I2y4z_Py_a = I_TWOBODYOVERLAP_K3y4z_S_a+ABY*I_TWOBODYOVERLAP_I2y4z_S_a;
  Double I_TWOBODYOVERLAP_Iy5z_Py_a = I_TWOBODYOVERLAP_K2y5z_S_a+ABY*I_TWOBODYOVERLAP_Iy5z_S_a;
  Double I_TWOBODYOVERLAP_I5xz_Pz_a = I_TWOBODYOVERLAP_K5x2z_S_a+ABZ*I_TWOBODYOVERLAP_I5xz_S_a;
  Double I_TWOBODYOVERLAP_I4xyz_Pz_a = I_TWOBODYOVERLAP_K4xy2z_S_a+ABZ*I_TWOBODYOVERLAP_I4xyz_S_a;
  Double I_TWOBODYOVERLAP_I4x2z_Pz_a = I_TWOBODYOVERLAP_K4x3z_S_a+ABZ*I_TWOBODYOVERLAP_I4x2z_S_a;
  Double I_TWOBODYOVERLAP_I3x2yz_Pz_a = I_TWOBODYOVERLAP_K3x2y2z_S_a+ABZ*I_TWOBODYOVERLAP_I3x2yz_S_a;
  Double I_TWOBODYOVERLAP_I3xy2z_Pz_a = I_TWOBODYOVERLAP_K3xy3z_S_a+ABZ*I_TWOBODYOVERLAP_I3xy2z_S_a;
  Double I_TWOBODYOVERLAP_I3x3z_Pz_a = I_TWOBODYOVERLAP_K3x4z_S_a+ABZ*I_TWOBODYOVERLAP_I3x3z_S_a;
  Double I_TWOBODYOVERLAP_I2x3yz_Pz_a = I_TWOBODYOVERLAP_K2x3y2z_S_a+ABZ*I_TWOBODYOVERLAP_I2x3yz_S_a;
  Double I_TWOBODYOVERLAP_I2x2y2z_Pz_a = I_TWOBODYOVERLAP_K2x2y3z_S_a+ABZ*I_TWOBODYOVERLAP_I2x2y2z_S_a;
  Double I_TWOBODYOVERLAP_I2xy3z_Pz_a = I_TWOBODYOVERLAP_K2xy4z_S_a+ABZ*I_TWOBODYOVERLAP_I2xy3z_S_a;
  Double I_TWOBODYOVERLAP_I2x4z_Pz_a = I_TWOBODYOVERLAP_K2x5z_S_a+ABZ*I_TWOBODYOVERLAP_I2x4z_S_a;
  Double I_TWOBODYOVERLAP_Ix4yz_Pz_a = I_TWOBODYOVERLAP_Kx4y2z_S_a+ABZ*I_TWOBODYOVERLAP_Ix4yz_S_a;
  Double I_TWOBODYOVERLAP_Ix3y2z_Pz_a = I_TWOBODYOVERLAP_Kx3y3z_S_a+ABZ*I_TWOBODYOVERLAP_Ix3y2z_S_a;
  Double I_TWOBODYOVERLAP_Ix2y3z_Pz_a = I_TWOBODYOVERLAP_Kx2y4z_S_a+ABZ*I_TWOBODYOVERLAP_Ix2y3z_S_a;
  Double I_TWOBODYOVERLAP_Ixy4z_Pz_a = I_TWOBODYOVERLAP_Kxy5z_S_a+ABZ*I_TWOBODYOVERLAP_Ixy4z_S_a;
  Double I_TWOBODYOVERLAP_Ix5z_Pz_a = I_TWOBODYOVERLAP_Kx6z_S_a+ABZ*I_TWOBODYOVERLAP_Ix5z_S_a;
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
   * totally 48 integrals are omitted 
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
  Double I_TWOBODYOVERLAP_H4xz_Dxy_a = I_TWOBODYOVERLAP_I4xyz_Px_a+ABY*I_TWOBODYOVERLAP_H4xz_Px_a;
  Double I_TWOBODYOVERLAP_H3xyz_Dxy_a = I_TWOBODYOVERLAP_I3x2yz_Px_a+ABY*I_TWOBODYOVERLAP_H3xyz_Px_a;
  Double I_TWOBODYOVERLAP_H3x2z_Dxy_a = I_TWOBODYOVERLAP_I3xy2z_Px_a+ABY*I_TWOBODYOVERLAP_H3x2z_Px_a;
  Double I_TWOBODYOVERLAP_H2x2yz_Dxy_a = I_TWOBODYOVERLAP_I2x3yz_Px_a+ABY*I_TWOBODYOVERLAP_H2x2yz_Px_a;
  Double I_TWOBODYOVERLAP_H2xy2z_Dxy_a = I_TWOBODYOVERLAP_I2x2y2z_Px_a+ABY*I_TWOBODYOVERLAP_H2xy2z_Px_a;
  Double I_TWOBODYOVERLAP_H2x3z_Dxy_a = I_TWOBODYOVERLAP_I2xy3z_Px_a+ABY*I_TWOBODYOVERLAP_H2x3z_Px_a;
  Double I_TWOBODYOVERLAP_Hx3yz_Dxy_a = I_TWOBODYOVERLAP_Ix4yz_Px_a+ABY*I_TWOBODYOVERLAP_Hx3yz_Px_a;
  Double I_TWOBODYOVERLAP_Hx2y2z_Dxy_a = I_TWOBODYOVERLAP_Ix3y2z_Px_a+ABY*I_TWOBODYOVERLAP_Hx2y2z_Px_a;
  Double I_TWOBODYOVERLAP_Hxy3z_Dxy_a = I_TWOBODYOVERLAP_Ix2y3z_Px_a+ABY*I_TWOBODYOVERLAP_Hxy3z_Px_a;
  Double I_TWOBODYOVERLAP_Hx4z_Dxy_a = I_TWOBODYOVERLAP_Ixy4z_Px_a+ABY*I_TWOBODYOVERLAP_Hx4z_Px_a;
  Double I_TWOBODYOVERLAP_H4yz_Dxy_a = I_TWOBODYOVERLAP_I5yz_Px_a+ABY*I_TWOBODYOVERLAP_H4yz_Px_a;
  Double I_TWOBODYOVERLAP_H3y2z_Dxy_a = I_TWOBODYOVERLAP_I4y2z_Px_a+ABY*I_TWOBODYOVERLAP_H3y2z_Px_a;
  Double I_TWOBODYOVERLAP_H2y3z_Dxy_a = I_TWOBODYOVERLAP_I3y3z_Px_a+ABY*I_TWOBODYOVERLAP_H2y3z_Px_a;
  Double I_TWOBODYOVERLAP_Hy4z_Dxy_a = I_TWOBODYOVERLAP_I2y4z_Px_a+ABY*I_TWOBODYOVERLAP_Hy4z_Px_a;
  Double I_TWOBODYOVERLAP_H5z_Dxy_a = I_TWOBODYOVERLAP_Iy5z_Px_a+ABY*I_TWOBODYOVERLAP_H5z_Px_a;
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
   * shell quartet name: SQ_TWOBODYOVERLAP_G_F_a
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_D_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_D_a
   ************************************************************/
  Double I_TWOBODYOVERLAP_G4x_F3x_a = I_TWOBODYOVERLAP_H5x_D2x_a+ABX*I_TWOBODYOVERLAP_G4x_D2x_a;
  Double I_TWOBODYOVERLAP_G3xy_F3x_a = I_TWOBODYOVERLAP_H4xy_D2x_a+ABX*I_TWOBODYOVERLAP_G3xy_D2x_a;
  Double I_TWOBODYOVERLAP_G3xz_F3x_a = I_TWOBODYOVERLAP_H4xz_D2x_a+ABX*I_TWOBODYOVERLAP_G3xz_D2x_a;
  Double I_TWOBODYOVERLAP_G2x2y_F3x_a = I_TWOBODYOVERLAP_H3x2y_D2x_a+ABX*I_TWOBODYOVERLAP_G2x2y_D2x_a;
  Double I_TWOBODYOVERLAP_G2xyz_F3x_a = I_TWOBODYOVERLAP_H3xyz_D2x_a+ABX*I_TWOBODYOVERLAP_G2xyz_D2x_a;
  Double I_TWOBODYOVERLAP_G2x2z_F3x_a = I_TWOBODYOVERLAP_H3x2z_D2x_a+ABX*I_TWOBODYOVERLAP_G2x2z_D2x_a;
  Double I_TWOBODYOVERLAP_Gx3y_F3x_a = I_TWOBODYOVERLAP_H2x3y_D2x_a+ABX*I_TWOBODYOVERLAP_Gx3y_D2x_a;
  Double I_TWOBODYOVERLAP_Gx2yz_F3x_a = I_TWOBODYOVERLAP_H2x2yz_D2x_a+ABX*I_TWOBODYOVERLAP_Gx2yz_D2x_a;
  Double I_TWOBODYOVERLAP_Gxy2z_F3x_a = I_TWOBODYOVERLAP_H2xy2z_D2x_a+ABX*I_TWOBODYOVERLAP_Gxy2z_D2x_a;
  Double I_TWOBODYOVERLAP_Gx3z_F3x_a = I_TWOBODYOVERLAP_H2x3z_D2x_a+ABX*I_TWOBODYOVERLAP_Gx3z_D2x_a;
  Double I_TWOBODYOVERLAP_G4y_F3x_a = I_TWOBODYOVERLAP_Hx4y_D2x_a+ABX*I_TWOBODYOVERLAP_G4y_D2x_a;
  Double I_TWOBODYOVERLAP_G3yz_F3x_a = I_TWOBODYOVERLAP_Hx3yz_D2x_a+ABX*I_TWOBODYOVERLAP_G3yz_D2x_a;
  Double I_TWOBODYOVERLAP_G2y2z_F3x_a = I_TWOBODYOVERLAP_Hx2y2z_D2x_a+ABX*I_TWOBODYOVERLAP_G2y2z_D2x_a;
  Double I_TWOBODYOVERLAP_Gy3z_F3x_a = I_TWOBODYOVERLAP_Hxy3z_D2x_a+ABX*I_TWOBODYOVERLAP_Gy3z_D2x_a;
  Double I_TWOBODYOVERLAP_G4z_F3x_a = I_TWOBODYOVERLAP_Hx4z_D2x_a+ABX*I_TWOBODYOVERLAP_G4z_D2x_a;
  Double I_TWOBODYOVERLAP_G4x_F2xy_a = I_TWOBODYOVERLAP_H4xy_D2x_a+ABY*I_TWOBODYOVERLAP_G4x_D2x_a;
  Double I_TWOBODYOVERLAP_G3xy_F2xy_a = I_TWOBODYOVERLAP_H3x2y_D2x_a+ABY*I_TWOBODYOVERLAP_G3xy_D2x_a;
  Double I_TWOBODYOVERLAP_G3xz_F2xy_a = I_TWOBODYOVERLAP_H3xyz_D2x_a+ABY*I_TWOBODYOVERLAP_G3xz_D2x_a;
  Double I_TWOBODYOVERLAP_G2x2y_F2xy_a = I_TWOBODYOVERLAP_H2x3y_D2x_a+ABY*I_TWOBODYOVERLAP_G2x2y_D2x_a;
  Double I_TWOBODYOVERLAP_G2xyz_F2xy_a = I_TWOBODYOVERLAP_H2x2yz_D2x_a+ABY*I_TWOBODYOVERLAP_G2xyz_D2x_a;
  Double I_TWOBODYOVERLAP_G2x2z_F2xy_a = I_TWOBODYOVERLAP_H2xy2z_D2x_a+ABY*I_TWOBODYOVERLAP_G2x2z_D2x_a;
  Double I_TWOBODYOVERLAP_Gx3y_F2xy_a = I_TWOBODYOVERLAP_Hx4y_D2x_a+ABY*I_TWOBODYOVERLAP_Gx3y_D2x_a;
  Double I_TWOBODYOVERLAP_Gx2yz_F2xy_a = I_TWOBODYOVERLAP_Hx3yz_D2x_a+ABY*I_TWOBODYOVERLAP_Gx2yz_D2x_a;
  Double I_TWOBODYOVERLAP_Gxy2z_F2xy_a = I_TWOBODYOVERLAP_Hx2y2z_D2x_a+ABY*I_TWOBODYOVERLAP_Gxy2z_D2x_a;
  Double I_TWOBODYOVERLAP_Gx3z_F2xy_a = I_TWOBODYOVERLAP_Hxy3z_D2x_a+ABY*I_TWOBODYOVERLAP_Gx3z_D2x_a;
  Double I_TWOBODYOVERLAP_G4y_F2xy_a = I_TWOBODYOVERLAP_H5y_D2x_a+ABY*I_TWOBODYOVERLAP_G4y_D2x_a;
  Double I_TWOBODYOVERLAP_G3yz_F2xy_a = I_TWOBODYOVERLAP_H4yz_D2x_a+ABY*I_TWOBODYOVERLAP_G3yz_D2x_a;
  Double I_TWOBODYOVERLAP_G2y2z_F2xy_a = I_TWOBODYOVERLAP_H3y2z_D2x_a+ABY*I_TWOBODYOVERLAP_G2y2z_D2x_a;
  Double I_TWOBODYOVERLAP_Gy3z_F2xy_a = I_TWOBODYOVERLAP_H2y3z_D2x_a+ABY*I_TWOBODYOVERLAP_Gy3z_D2x_a;
  Double I_TWOBODYOVERLAP_G4z_F2xy_a = I_TWOBODYOVERLAP_Hy4z_D2x_a+ABY*I_TWOBODYOVERLAP_G4z_D2x_a;
  Double I_TWOBODYOVERLAP_G4x_F2xz_a = I_TWOBODYOVERLAP_H4xz_D2x_a+ABZ*I_TWOBODYOVERLAP_G4x_D2x_a;
  Double I_TWOBODYOVERLAP_G3xy_F2xz_a = I_TWOBODYOVERLAP_H3xyz_D2x_a+ABZ*I_TWOBODYOVERLAP_G3xy_D2x_a;
  Double I_TWOBODYOVERLAP_G3xz_F2xz_a = I_TWOBODYOVERLAP_H3x2z_D2x_a+ABZ*I_TWOBODYOVERLAP_G3xz_D2x_a;
  Double I_TWOBODYOVERLAP_G2x2y_F2xz_a = I_TWOBODYOVERLAP_H2x2yz_D2x_a+ABZ*I_TWOBODYOVERLAP_G2x2y_D2x_a;
  Double I_TWOBODYOVERLAP_G2xyz_F2xz_a = I_TWOBODYOVERLAP_H2xy2z_D2x_a+ABZ*I_TWOBODYOVERLAP_G2xyz_D2x_a;
  Double I_TWOBODYOVERLAP_G2x2z_F2xz_a = I_TWOBODYOVERLAP_H2x3z_D2x_a+ABZ*I_TWOBODYOVERLAP_G2x2z_D2x_a;
  Double I_TWOBODYOVERLAP_Gx3y_F2xz_a = I_TWOBODYOVERLAP_Hx3yz_D2x_a+ABZ*I_TWOBODYOVERLAP_Gx3y_D2x_a;
  Double I_TWOBODYOVERLAP_Gx2yz_F2xz_a = I_TWOBODYOVERLAP_Hx2y2z_D2x_a+ABZ*I_TWOBODYOVERLAP_Gx2yz_D2x_a;
  Double I_TWOBODYOVERLAP_Gxy2z_F2xz_a = I_TWOBODYOVERLAP_Hxy3z_D2x_a+ABZ*I_TWOBODYOVERLAP_Gxy2z_D2x_a;
  Double I_TWOBODYOVERLAP_Gx3z_F2xz_a = I_TWOBODYOVERLAP_Hx4z_D2x_a+ABZ*I_TWOBODYOVERLAP_Gx3z_D2x_a;
  Double I_TWOBODYOVERLAP_G4y_F2xz_a = I_TWOBODYOVERLAP_H4yz_D2x_a+ABZ*I_TWOBODYOVERLAP_G4y_D2x_a;
  Double I_TWOBODYOVERLAP_G3yz_F2xz_a = I_TWOBODYOVERLAP_H3y2z_D2x_a+ABZ*I_TWOBODYOVERLAP_G3yz_D2x_a;
  Double I_TWOBODYOVERLAP_G2y2z_F2xz_a = I_TWOBODYOVERLAP_H2y3z_D2x_a+ABZ*I_TWOBODYOVERLAP_G2y2z_D2x_a;
  Double I_TWOBODYOVERLAP_Gy3z_F2xz_a = I_TWOBODYOVERLAP_Hy4z_D2x_a+ABZ*I_TWOBODYOVERLAP_Gy3z_D2x_a;
  Double I_TWOBODYOVERLAP_G4z_F2xz_a = I_TWOBODYOVERLAP_H5z_D2x_a+ABZ*I_TWOBODYOVERLAP_G4z_D2x_a;
  Double I_TWOBODYOVERLAP_G4x_Fx2y_a = I_TWOBODYOVERLAP_H5x_D2y_a+ABX*I_TWOBODYOVERLAP_G4x_D2y_a;
  Double I_TWOBODYOVERLAP_G3xy_Fx2y_a = I_TWOBODYOVERLAP_H4xy_D2y_a+ABX*I_TWOBODYOVERLAP_G3xy_D2y_a;
  Double I_TWOBODYOVERLAP_G3xz_Fx2y_a = I_TWOBODYOVERLAP_H4xz_D2y_a+ABX*I_TWOBODYOVERLAP_G3xz_D2y_a;
  Double I_TWOBODYOVERLAP_G2x2y_Fx2y_a = I_TWOBODYOVERLAP_H3x2y_D2y_a+ABX*I_TWOBODYOVERLAP_G2x2y_D2y_a;
  Double I_TWOBODYOVERLAP_G2xyz_Fx2y_a = I_TWOBODYOVERLAP_H3xyz_D2y_a+ABX*I_TWOBODYOVERLAP_G2xyz_D2y_a;
  Double I_TWOBODYOVERLAP_G2x2z_Fx2y_a = I_TWOBODYOVERLAP_H3x2z_D2y_a+ABX*I_TWOBODYOVERLAP_G2x2z_D2y_a;
  Double I_TWOBODYOVERLAP_Gx3y_Fx2y_a = I_TWOBODYOVERLAP_H2x3y_D2y_a+ABX*I_TWOBODYOVERLAP_Gx3y_D2y_a;
  Double I_TWOBODYOVERLAP_Gx2yz_Fx2y_a = I_TWOBODYOVERLAP_H2x2yz_D2y_a+ABX*I_TWOBODYOVERLAP_Gx2yz_D2y_a;
  Double I_TWOBODYOVERLAP_Gxy2z_Fx2y_a = I_TWOBODYOVERLAP_H2xy2z_D2y_a+ABX*I_TWOBODYOVERLAP_Gxy2z_D2y_a;
  Double I_TWOBODYOVERLAP_Gx3z_Fx2y_a = I_TWOBODYOVERLAP_H2x3z_D2y_a+ABX*I_TWOBODYOVERLAP_Gx3z_D2y_a;
  Double I_TWOBODYOVERLAP_G4y_Fx2y_a = I_TWOBODYOVERLAP_Hx4y_D2y_a+ABX*I_TWOBODYOVERLAP_G4y_D2y_a;
  Double I_TWOBODYOVERLAP_G3yz_Fx2y_a = I_TWOBODYOVERLAP_Hx3yz_D2y_a+ABX*I_TWOBODYOVERLAP_G3yz_D2y_a;
  Double I_TWOBODYOVERLAP_G2y2z_Fx2y_a = I_TWOBODYOVERLAP_Hx2y2z_D2y_a+ABX*I_TWOBODYOVERLAP_G2y2z_D2y_a;
  Double I_TWOBODYOVERLAP_Gy3z_Fx2y_a = I_TWOBODYOVERLAP_Hxy3z_D2y_a+ABX*I_TWOBODYOVERLAP_Gy3z_D2y_a;
  Double I_TWOBODYOVERLAP_G4z_Fx2y_a = I_TWOBODYOVERLAP_Hx4z_D2y_a+ABX*I_TWOBODYOVERLAP_G4z_D2y_a;
  Double I_TWOBODYOVERLAP_G4x_Fxyz_a = I_TWOBODYOVERLAP_H4xz_Dxy_a+ABZ*I_TWOBODYOVERLAP_G4x_Dxy_a;
  Double I_TWOBODYOVERLAP_G3xy_Fxyz_a = I_TWOBODYOVERLAP_H3xyz_Dxy_a+ABZ*I_TWOBODYOVERLAP_G3xy_Dxy_a;
  Double I_TWOBODYOVERLAP_G3xz_Fxyz_a = I_TWOBODYOVERLAP_H3x2z_Dxy_a+ABZ*I_TWOBODYOVERLAP_G3xz_Dxy_a;
  Double I_TWOBODYOVERLAP_G2x2y_Fxyz_a = I_TWOBODYOVERLAP_H2x2yz_Dxy_a+ABZ*I_TWOBODYOVERLAP_G2x2y_Dxy_a;
  Double I_TWOBODYOVERLAP_G2xyz_Fxyz_a = I_TWOBODYOVERLAP_H2xy2z_Dxy_a+ABZ*I_TWOBODYOVERLAP_G2xyz_Dxy_a;
  Double I_TWOBODYOVERLAP_G2x2z_Fxyz_a = I_TWOBODYOVERLAP_H2x3z_Dxy_a+ABZ*I_TWOBODYOVERLAP_G2x2z_Dxy_a;
  Double I_TWOBODYOVERLAP_Gx3y_Fxyz_a = I_TWOBODYOVERLAP_Hx3yz_Dxy_a+ABZ*I_TWOBODYOVERLAP_Gx3y_Dxy_a;
  Double I_TWOBODYOVERLAP_Gx2yz_Fxyz_a = I_TWOBODYOVERLAP_Hx2y2z_Dxy_a+ABZ*I_TWOBODYOVERLAP_Gx2yz_Dxy_a;
  Double I_TWOBODYOVERLAP_Gxy2z_Fxyz_a = I_TWOBODYOVERLAP_Hxy3z_Dxy_a+ABZ*I_TWOBODYOVERLAP_Gxy2z_Dxy_a;
  Double I_TWOBODYOVERLAP_Gx3z_Fxyz_a = I_TWOBODYOVERLAP_Hx4z_Dxy_a+ABZ*I_TWOBODYOVERLAP_Gx3z_Dxy_a;
  Double I_TWOBODYOVERLAP_G4y_Fxyz_a = I_TWOBODYOVERLAP_H4yz_Dxy_a+ABZ*I_TWOBODYOVERLAP_G4y_Dxy_a;
  Double I_TWOBODYOVERLAP_G3yz_Fxyz_a = I_TWOBODYOVERLAP_H3y2z_Dxy_a+ABZ*I_TWOBODYOVERLAP_G3yz_Dxy_a;
  Double I_TWOBODYOVERLAP_G2y2z_Fxyz_a = I_TWOBODYOVERLAP_H2y3z_Dxy_a+ABZ*I_TWOBODYOVERLAP_G2y2z_Dxy_a;
  Double I_TWOBODYOVERLAP_Gy3z_Fxyz_a = I_TWOBODYOVERLAP_Hy4z_Dxy_a+ABZ*I_TWOBODYOVERLAP_Gy3z_Dxy_a;
  Double I_TWOBODYOVERLAP_G4z_Fxyz_a = I_TWOBODYOVERLAP_H5z_Dxy_a+ABZ*I_TWOBODYOVERLAP_G4z_Dxy_a;
  Double I_TWOBODYOVERLAP_G4x_Fx2z_a = I_TWOBODYOVERLAP_H5x_D2z_a+ABX*I_TWOBODYOVERLAP_G4x_D2z_a;
  Double I_TWOBODYOVERLAP_G3xy_Fx2z_a = I_TWOBODYOVERLAP_H4xy_D2z_a+ABX*I_TWOBODYOVERLAP_G3xy_D2z_a;
  Double I_TWOBODYOVERLAP_G3xz_Fx2z_a = I_TWOBODYOVERLAP_H4xz_D2z_a+ABX*I_TWOBODYOVERLAP_G3xz_D2z_a;
  Double I_TWOBODYOVERLAP_G2x2y_Fx2z_a = I_TWOBODYOVERLAP_H3x2y_D2z_a+ABX*I_TWOBODYOVERLAP_G2x2y_D2z_a;
  Double I_TWOBODYOVERLAP_G2xyz_Fx2z_a = I_TWOBODYOVERLAP_H3xyz_D2z_a+ABX*I_TWOBODYOVERLAP_G2xyz_D2z_a;
  Double I_TWOBODYOVERLAP_G2x2z_Fx2z_a = I_TWOBODYOVERLAP_H3x2z_D2z_a+ABX*I_TWOBODYOVERLAP_G2x2z_D2z_a;
  Double I_TWOBODYOVERLAP_Gx3y_Fx2z_a = I_TWOBODYOVERLAP_H2x3y_D2z_a+ABX*I_TWOBODYOVERLAP_Gx3y_D2z_a;
  Double I_TWOBODYOVERLAP_Gx2yz_Fx2z_a = I_TWOBODYOVERLAP_H2x2yz_D2z_a+ABX*I_TWOBODYOVERLAP_Gx2yz_D2z_a;
  Double I_TWOBODYOVERLAP_Gxy2z_Fx2z_a = I_TWOBODYOVERLAP_H2xy2z_D2z_a+ABX*I_TWOBODYOVERLAP_Gxy2z_D2z_a;
  Double I_TWOBODYOVERLAP_Gx3z_Fx2z_a = I_TWOBODYOVERLAP_H2x3z_D2z_a+ABX*I_TWOBODYOVERLAP_Gx3z_D2z_a;
  Double I_TWOBODYOVERLAP_G4y_Fx2z_a = I_TWOBODYOVERLAP_Hx4y_D2z_a+ABX*I_TWOBODYOVERLAP_G4y_D2z_a;
  Double I_TWOBODYOVERLAP_G3yz_Fx2z_a = I_TWOBODYOVERLAP_Hx3yz_D2z_a+ABX*I_TWOBODYOVERLAP_G3yz_D2z_a;
  Double I_TWOBODYOVERLAP_G2y2z_Fx2z_a = I_TWOBODYOVERLAP_Hx2y2z_D2z_a+ABX*I_TWOBODYOVERLAP_G2y2z_D2z_a;
  Double I_TWOBODYOVERLAP_Gy3z_Fx2z_a = I_TWOBODYOVERLAP_Hxy3z_D2z_a+ABX*I_TWOBODYOVERLAP_Gy3z_D2z_a;
  Double I_TWOBODYOVERLAP_G4z_Fx2z_a = I_TWOBODYOVERLAP_Hx4z_D2z_a+ABX*I_TWOBODYOVERLAP_G4z_D2z_a;
  Double I_TWOBODYOVERLAP_G4x_F3y_a = I_TWOBODYOVERLAP_H4xy_D2y_a+ABY*I_TWOBODYOVERLAP_G4x_D2y_a;
  Double I_TWOBODYOVERLAP_G3xy_F3y_a = I_TWOBODYOVERLAP_H3x2y_D2y_a+ABY*I_TWOBODYOVERLAP_G3xy_D2y_a;
  Double I_TWOBODYOVERLAP_G3xz_F3y_a = I_TWOBODYOVERLAP_H3xyz_D2y_a+ABY*I_TWOBODYOVERLAP_G3xz_D2y_a;
  Double I_TWOBODYOVERLAP_G2x2y_F3y_a = I_TWOBODYOVERLAP_H2x3y_D2y_a+ABY*I_TWOBODYOVERLAP_G2x2y_D2y_a;
  Double I_TWOBODYOVERLAP_G2xyz_F3y_a = I_TWOBODYOVERLAP_H2x2yz_D2y_a+ABY*I_TWOBODYOVERLAP_G2xyz_D2y_a;
  Double I_TWOBODYOVERLAP_G2x2z_F3y_a = I_TWOBODYOVERLAP_H2xy2z_D2y_a+ABY*I_TWOBODYOVERLAP_G2x2z_D2y_a;
  Double I_TWOBODYOVERLAP_Gx3y_F3y_a = I_TWOBODYOVERLAP_Hx4y_D2y_a+ABY*I_TWOBODYOVERLAP_Gx3y_D2y_a;
  Double I_TWOBODYOVERLAP_Gx2yz_F3y_a = I_TWOBODYOVERLAP_Hx3yz_D2y_a+ABY*I_TWOBODYOVERLAP_Gx2yz_D2y_a;
  Double I_TWOBODYOVERLAP_Gxy2z_F3y_a = I_TWOBODYOVERLAP_Hx2y2z_D2y_a+ABY*I_TWOBODYOVERLAP_Gxy2z_D2y_a;
  Double I_TWOBODYOVERLAP_Gx3z_F3y_a = I_TWOBODYOVERLAP_Hxy3z_D2y_a+ABY*I_TWOBODYOVERLAP_Gx3z_D2y_a;
  Double I_TWOBODYOVERLAP_G4y_F3y_a = I_TWOBODYOVERLAP_H5y_D2y_a+ABY*I_TWOBODYOVERLAP_G4y_D2y_a;
  Double I_TWOBODYOVERLAP_G3yz_F3y_a = I_TWOBODYOVERLAP_H4yz_D2y_a+ABY*I_TWOBODYOVERLAP_G3yz_D2y_a;
  Double I_TWOBODYOVERLAP_G2y2z_F3y_a = I_TWOBODYOVERLAP_H3y2z_D2y_a+ABY*I_TWOBODYOVERLAP_G2y2z_D2y_a;
  Double I_TWOBODYOVERLAP_Gy3z_F3y_a = I_TWOBODYOVERLAP_H2y3z_D2y_a+ABY*I_TWOBODYOVERLAP_Gy3z_D2y_a;
  Double I_TWOBODYOVERLAP_G4z_F3y_a = I_TWOBODYOVERLAP_Hy4z_D2y_a+ABY*I_TWOBODYOVERLAP_G4z_D2y_a;
  Double I_TWOBODYOVERLAP_G4x_F2yz_a = I_TWOBODYOVERLAP_H4xz_D2y_a+ABZ*I_TWOBODYOVERLAP_G4x_D2y_a;
  Double I_TWOBODYOVERLAP_G3xy_F2yz_a = I_TWOBODYOVERLAP_H3xyz_D2y_a+ABZ*I_TWOBODYOVERLAP_G3xy_D2y_a;
  Double I_TWOBODYOVERLAP_G3xz_F2yz_a = I_TWOBODYOVERLAP_H3x2z_D2y_a+ABZ*I_TWOBODYOVERLAP_G3xz_D2y_a;
  Double I_TWOBODYOVERLAP_G2x2y_F2yz_a = I_TWOBODYOVERLAP_H2x2yz_D2y_a+ABZ*I_TWOBODYOVERLAP_G2x2y_D2y_a;
  Double I_TWOBODYOVERLAP_G2xyz_F2yz_a = I_TWOBODYOVERLAP_H2xy2z_D2y_a+ABZ*I_TWOBODYOVERLAP_G2xyz_D2y_a;
  Double I_TWOBODYOVERLAP_G2x2z_F2yz_a = I_TWOBODYOVERLAP_H2x3z_D2y_a+ABZ*I_TWOBODYOVERLAP_G2x2z_D2y_a;
  Double I_TWOBODYOVERLAP_Gx3y_F2yz_a = I_TWOBODYOVERLAP_Hx3yz_D2y_a+ABZ*I_TWOBODYOVERLAP_Gx3y_D2y_a;
  Double I_TWOBODYOVERLAP_Gx2yz_F2yz_a = I_TWOBODYOVERLAP_Hx2y2z_D2y_a+ABZ*I_TWOBODYOVERLAP_Gx2yz_D2y_a;
  Double I_TWOBODYOVERLAP_Gxy2z_F2yz_a = I_TWOBODYOVERLAP_Hxy3z_D2y_a+ABZ*I_TWOBODYOVERLAP_Gxy2z_D2y_a;
  Double I_TWOBODYOVERLAP_Gx3z_F2yz_a = I_TWOBODYOVERLAP_Hx4z_D2y_a+ABZ*I_TWOBODYOVERLAP_Gx3z_D2y_a;
  Double I_TWOBODYOVERLAP_G4y_F2yz_a = I_TWOBODYOVERLAP_H4yz_D2y_a+ABZ*I_TWOBODYOVERLAP_G4y_D2y_a;
  Double I_TWOBODYOVERLAP_G3yz_F2yz_a = I_TWOBODYOVERLAP_H3y2z_D2y_a+ABZ*I_TWOBODYOVERLAP_G3yz_D2y_a;
  Double I_TWOBODYOVERLAP_G2y2z_F2yz_a = I_TWOBODYOVERLAP_H2y3z_D2y_a+ABZ*I_TWOBODYOVERLAP_G2y2z_D2y_a;
  Double I_TWOBODYOVERLAP_Gy3z_F2yz_a = I_TWOBODYOVERLAP_Hy4z_D2y_a+ABZ*I_TWOBODYOVERLAP_Gy3z_D2y_a;
  Double I_TWOBODYOVERLAP_G4z_F2yz_a = I_TWOBODYOVERLAP_H5z_D2y_a+ABZ*I_TWOBODYOVERLAP_G4z_D2y_a;
  Double I_TWOBODYOVERLAP_G4x_Fy2z_a = I_TWOBODYOVERLAP_H4xy_D2z_a+ABY*I_TWOBODYOVERLAP_G4x_D2z_a;
  Double I_TWOBODYOVERLAP_G3xy_Fy2z_a = I_TWOBODYOVERLAP_H3x2y_D2z_a+ABY*I_TWOBODYOVERLAP_G3xy_D2z_a;
  Double I_TWOBODYOVERLAP_G3xz_Fy2z_a = I_TWOBODYOVERLAP_H3xyz_D2z_a+ABY*I_TWOBODYOVERLAP_G3xz_D2z_a;
  Double I_TWOBODYOVERLAP_G2x2y_Fy2z_a = I_TWOBODYOVERLAP_H2x3y_D2z_a+ABY*I_TWOBODYOVERLAP_G2x2y_D2z_a;
  Double I_TWOBODYOVERLAP_G2xyz_Fy2z_a = I_TWOBODYOVERLAP_H2x2yz_D2z_a+ABY*I_TWOBODYOVERLAP_G2xyz_D2z_a;
  Double I_TWOBODYOVERLAP_G2x2z_Fy2z_a = I_TWOBODYOVERLAP_H2xy2z_D2z_a+ABY*I_TWOBODYOVERLAP_G2x2z_D2z_a;
  Double I_TWOBODYOVERLAP_Gx3y_Fy2z_a = I_TWOBODYOVERLAP_Hx4y_D2z_a+ABY*I_TWOBODYOVERLAP_Gx3y_D2z_a;
  Double I_TWOBODYOVERLAP_Gx2yz_Fy2z_a = I_TWOBODYOVERLAP_Hx3yz_D2z_a+ABY*I_TWOBODYOVERLAP_Gx2yz_D2z_a;
  Double I_TWOBODYOVERLAP_Gxy2z_Fy2z_a = I_TWOBODYOVERLAP_Hx2y2z_D2z_a+ABY*I_TWOBODYOVERLAP_Gxy2z_D2z_a;
  Double I_TWOBODYOVERLAP_Gx3z_Fy2z_a = I_TWOBODYOVERLAP_Hxy3z_D2z_a+ABY*I_TWOBODYOVERLAP_Gx3z_D2z_a;
  Double I_TWOBODYOVERLAP_G4y_Fy2z_a = I_TWOBODYOVERLAP_H5y_D2z_a+ABY*I_TWOBODYOVERLAP_G4y_D2z_a;
  Double I_TWOBODYOVERLAP_G3yz_Fy2z_a = I_TWOBODYOVERLAP_H4yz_D2z_a+ABY*I_TWOBODYOVERLAP_G3yz_D2z_a;
  Double I_TWOBODYOVERLAP_G2y2z_Fy2z_a = I_TWOBODYOVERLAP_H3y2z_D2z_a+ABY*I_TWOBODYOVERLAP_G2y2z_D2z_a;
  Double I_TWOBODYOVERLAP_Gy3z_Fy2z_a = I_TWOBODYOVERLAP_H2y3z_D2z_a+ABY*I_TWOBODYOVERLAP_Gy3z_D2z_a;
  Double I_TWOBODYOVERLAP_G4z_Fy2z_a = I_TWOBODYOVERLAP_Hy4z_D2z_a+ABY*I_TWOBODYOVERLAP_G4z_D2z_a;
  Double I_TWOBODYOVERLAP_G4x_F3z_a = I_TWOBODYOVERLAP_H4xz_D2z_a+ABZ*I_TWOBODYOVERLAP_G4x_D2z_a;
  Double I_TWOBODYOVERLAP_G3xy_F3z_a = I_TWOBODYOVERLAP_H3xyz_D2z_a+ABZ*I_TWOBODYOVERLAP_G3xy_D2z_a;
  Double I_TWOBODYOVERLAP_G3xz_F3z_a = I_TWOBODYOVERLAP_H3x2z_D2z_a+ABZ*I_TWOBODYOVERLAP_G3xz_D2z_a;
  Double I_TWOBODYOVERLAP_G2x2y_F3z_a = I_TWOBODYOVERLAP_H2x2yz_D2z_a+ABZ*I_TWOBODYOVERLAP_G2x2y_D2z_a;
  Double I_TWOBODYOVERLAP_G2xyz_F3z_a = I_TWOBODYOVERLAP_H2xy2z_D2z_a+ABZ*I_TWOBODYOVERLAP_G2xyz_D2z_a;
  Double I_TWOBODYOVERLAP_G2x2z_F3z_a = I_TWOBODYOVERLAP_H2x3z_D2z_a+ABZ*I_TWOBODYOVERLAP_G2x2z_D2z_a;
  Double I_TWOBODYOVERLAP_Gx3y_F3z_a = I_TWOBODYOVERLAP_Hx3yz_D2z_a+ABZ*I_TWOBODYOVERLAP_Gx3y_D2z_a;
  Double I_TWOBODYOVERLAP_Gx2yz_F3z_a = I_TWOBODYOVERLAP_Hx2y2z_D2z_a+ABZ*I_TWOBODYOVERLAP_Gx2yz_D2z_a;
  Double I_TWOBODYOVERLAP_Gxy2z_F3z_a = I_TWOBODYOVERLAP_Hxy3z_D2z_a+ABZ*I_TWOBODYOVERLAP_Gxy2z_D2z_a;
  Double I_TWOBODYOVERLAP_Gx3z_F3z_a = I_TWOBODYOVERLAP_Hx4z_D2z_a+ABZ*I_TWOBODYOVERLAP_Gx3z_D2z_a;
  Double I_TWOBODYOVERLAP_G4y_F3z_a = I_TWOBODYOVERLAP_H4yz_D2z_a+ABZ*I_TWOBODYOVERLAP_G4y_D2z_a;
  Double I_TWOBODYOVERLAP_G3yz_F3z_a = I_TWOBODYOVERLAP_H3y2z_D2z_a+ABZ*I_TWOBODYOVERLAP_G3yz_D2z_a;
  Double I_TWOBODYOVERLAP_G2y2z_F3z_a = I_TWOBODYOVERLAP_H2y3z_D2z_a+ABZ*I_TWOBODYOVERLAP_G2y2z_D2z_a;
  Double I_TWOBODYOVERLAP_Gy3z_F3z_a = I_TWOBODYOVERLAP_Hy4z_D2z_a+ABZ*I_TWOBODYOVERLAP_Gy3z_D2z_a;
  Double I_TWOBODYOVERLAP_G4z_F3z_a = I_TWOBODYOVERLAP_H5z_D2z_a+ABZ*I_TWOBODYOVERLAP_G4z_D2z_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_I_P_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_S_aa
   ************************************************************/
  Double I_TWOBODYOVERLAP_I6x_Px_aa = I_TWOBODYOVERLAP_K7x_S_aa+ABX*I_TWOBODYOVERLAP_I6x_S_aa;
  Double I_TWOBODYOVERLAP_I5xy_Px_aa = I_TWOBODYOVERLAP_K6xy_S_aa+ABX*I_TWOBODYOVERLAP_I5xy_S_aa;
  Double I_TWOBODYOVERLAP_I5xz_Px_aa = I_TWOBODYOVERLAP_K6xz_S_aa+ABX*I_TWOBODYOVERLAP_I5xz_S_aa;
  Double I_TWOBODYOVERLAP_I4x2y_Px_aa = I_TWOBODYOVERLAP_K5x2y_S_aa+ABX*I_TWOBODYOVERLAP_I4x2y_S_aa;
  Double I_TWOBODYOVERLAP_I4xyz_Px_aa = I_TWOBODYOVERLAP_K5xyz_S_aa+ABX*I_TWOBODYOVERLAP_I4xyz_S_aa;
  Double I_TWOBODYOVERLAP_I4x2z_Px_aa = I_TWOBODYOVERLAP_K5x2z_S_aa+ABX*I_TWOBODYOVERLAP_I4x2z_S_aa;
  Double I_TWOBODYOVERLAP_I3x3y_Px_aa = I_TWOBODYOVERLAP_K4x3y_S_aa+ABX*I_TWOBODYOVERLAP_I3x3y_S_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_Px_aa = I_TWOBODYOVERLAP_K4x2yz_S_aa+ABX*I_TWOBODYOVERLAP_I3x2yz_S_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_Px_aa = I_TWOBODYOVERLAP_K4xy2z_S_aa+ABX*I_TWOBODYOVERLAP_I3xy2z_S_aa;
  Double I_TWOBODYOVERLAP_I3x3z_Px_aa = I_TWOBODYOVERLAP_K4x3z_S_aa+ABX*I_TWOBODYOVERLAP_I3x3z_S_aa;
  Double I_TWOBODYOVERLAP_I2x4y_Px_aa = I_TWOBODYOVERLAP_K3x4y_S_aa+ABX*I_TWOBODYOVERLAP_I2x4y_S_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_Px_aa = I_TWOBODYOVERLAP_K3x3yz_S_aa+ABX*I_TWOBODYOVERLAP_I2x3yz_S_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_Px_aa = I_TWOBODYOVERLAP_K3x2y2z_S_aa+ABX*I_TWOBODYOVERLAP_I2x2y2z_S_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_Px_aa = I_TWOBODYOVERLAP_K3xy3z_S_aa+ABX*I_TWOBODYOVERLAP_I2xy3z_S_aa;
  Double I_TWOBODYOVERLAP_I2x4z_Px_aa = I_TWOBODYOVERLAP_K3x4z_S_aa+ABX*I_TWOBODYOVERLAP_I2x4z_S_aa;
  Double I_TWOBODYOVERLAP_Ix5y_Px_aa = I_TWOBODYOVERLAP_K2x5y_S_aa+ABX*I_TWOBODYOVERLAP_Ix5y_S_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_Px_aa = I_TWOBODYOVERLAP_K2x4yz_S_aa+ABX*I_TWOBODYOVERLAP_Ix4yz_S_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_Px_aa = I_TWOBODYOVERLAP_K2x3y2z_S_aa+ABX*I_TWOBODYOVERLAP_Ix3y2z_S_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_Px_aa = I_TWOBODYOVERLAP_K2x2y3z_S_aa+ABX*I_TWOBODYOVERLAP_Ix2y3z_S_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_Px_aa = I_TWOBODYOVERLAP_K2xy4z_S_aa+ABX*I_TWOBODYOVERLAP_Ixy4z_S_aa;
  Double I_TWOBODYOVERLAP_Ix5z_Px_aa = I_TWOBODYOVERLAP_K2x5z_S_aa+ABX*I_TWOBODYOVERLAP_Ix5z_S_aa;
  Double I_TWOBODYOVERLAP_I6y_Px_aa = I_TWOBODYOVERLAP_Kx6y_S_aa+ABX*I_TWOBODYOVERLAP_I6y_S_aa;
  Double I_TWOBODYOVERLAP_I5yz_Px_aa = I_TWOBODYOVERLAP_Kx5yz_S_aa+ABX*I_TWOBODYOVERLAP_I5yz_S_aa;
  Double I_TWOBODYOVERLAP_I4y2z_Px_aa = I_TWOBODYOVERLAP_Kx4y2z_S_aa+ABX*I_TWOBODYOVERLAP_I4y2z_S_aa;
  Double I_TWOBODYOVERLAP_I3y3z_Px_aa = I_TWOBODYOVERLAP_Kx3y3z_S_aa+ABX*I_TWOBODYOVERLAP_I3y3z_S_aa;
  Double I_TWOBODYOVERLAP_I2y4z_Px_aa = I_TWOBODYOVERLAP_Kx2y4z_S_aa+ABX*I_TWOBODYOVERLAP_I2y4z_S_aa;
  Double I_TWOBODYOVERLAP_Iy5z_Px_aa = I_TWOBODYOVERLAP_Kxy5z_S_aa+ABX*I_TWOBODYOVERLAP_Iy5z_S_aa;
  Double I_TWOBODYOVERLAP_I6z_Px_aa = I_TWOBODYOVERLAP_Kx6z_S_aa+ABX*I_TWOBODYOVERLAP_I6z_S_aa;
  Double I_TWOBODYOVERLAP_I6x_Py_aa = I_TWOBODYOVERLAP_K6xy_S_aa+ABY*I_TWOBODYOVERLAP_I6x_S_aa;
  Double I_TWOBODYOVERLAP_I5xy_Py_aa = I_TWOBODYOVERLAP_K5x2y_S_aa+ABY*I_TWOBODYOVERLAP_I5xy_S_aa;
  Double I_TWOBODYOVERLAP_I5xz_Py_aa = I_TWOBODYOVERLAP_K5xyz_S_aa+ABY*I_TWOBODYOVERLAP_I5xz_S_aa;
  Double I_TWOBODYOVERLAP_I4x2y_Py_aa = I_TWOBODYOVERLAP_K4x3y_S_aa+ABY*I_TWOBODYOVERLAP_I4x2y_S_aa;
  Double I_TWOBODYOVERLAP_I4xyz_Py_aa = I_TWOBODYOVERLAP_K4x2yz_S_aa+ABY*I_TWOBODYOVERLAP_I4xyz_S_aa;
  Double I_TWOBODYOVERLAP_I4x2z_Py_aa = I_TWOBODYOVERLAP_K4xy2z_S_aa+ABY*I_TWOBODYOVERLAP_I4x2z_S_aa;
  Double I_TWOBODYOVERLAP_I3x3y_Py_aa = I_TWOBODYOVERLAP_K3x4y_S_aa+ABY*I_TWOBODYOVERLAP_I3x3y_S_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_Py_aa = I_TWOBODYOVERLAP_K3x3yz_S_aa+ABY*I_TWOBODYOVERLAP_I3x2yz_S_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_Py_aa = I_TWOBODYOVERLAP_K3x2y2z_S_aa+ABY*I_TWOBODYOVERLAP_I3xy2z_S_aa;
  Double I_TWOBODYOVERLAP_I3x3z_Py_aa = I_TWOBODYOVERLAP_K3xy3z_S_aa+ABY*I_TWOBODYOVERLAP_I3x3z_S_aa;
  Double I_TWOBODYOVERLAP_I2x4y_Py_aa = I_TWOBODYOVERLAP_K2x5y_S_aa+ABY*I_TWOBODYOVERLAP_I2x4y_S_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_Py_aa = I_TWOBODYOVERLAP_K2x4yz_S_aa+ABY*I_TWOBODYOVERLAP_I2x3yz_S_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_Py_aa = I_TWOBODYOVERLAP_K2x3y2z_S_aa+ABY*I_TWOBODYOVERLAP_I2x2y2z_S_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_Py_aa = I_TWOBODYOVERLAP_K2x2y3z_S_aa+ABY*I_TWOBODYOVERLAP_I2xy3z_S_aa;
  Double I_TWOBODYOVERLAP_I2x4z_Py_aa = I_TWOBODYOVERLAP_K2xy4z_S_aa+ABY*I_TWOBODYOVERLAP_I2x4z_S_aa;
  Double I_TWOBODYOVERLAP_Ix5y_Py_aa = I_TWOBODYOVERLAP_Kx6y_S_aa+ABY*I_TWOBODYOVERLAP_Ix5y_S_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_Py_aa = I_TWOBODYOVERLAP_Kx5yz_S_aa+ABY*I_TWOBODYOVERLAP_Ix4yz_S_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_Py_aa = I_TWOBODYOVERLAP_Kx4y2z_S_aa+ABY*I_TWOBODYOVERLAP_Ix3y2z_S_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_Py_aa = I_TWOBODYOVERLAP_Kx3y3z_S_aa+ABY*I_TWOBODYOVERLAP_Ix2y3z_S_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_Py_aa = I_TWOBODYOVERLAP_Kx2y4z_S_aa+ABY*I_TWOBODYOVERLAP_Ixy4z_S_aa;
  Double I_TWOBODYOVERLAP_Ix5z_Py_aa = I_TWOBODYOVERLAP_Kxy5z_S_aa+ABY*I_TWOBODYOVERLAP_Ix5z_S_aa;
  Double I_TWOBODYOVERLAP_I6y_Py_aa = I_TWOBODYOVERLAP_K7y_S_aa+ABY*I_TWOBODYOVERLAP_I6y_S_aa;
  Double I_TWOBODYOVERLAP_I5yz_Py_aa = I_TWOBODYOVERLAP_K6yz_S_aa+ABY*I_TWOBODYOVERLAP_I5yz_S_aa;
  Double I_TWOBODYOVERLAP_I4y2z_Py_aa = I_TWOBODYOVERLAP_K5y2z_S_aa+ABY*I_TWOBODYOVERLAP_I4y2z_S_aa;
  Double I_TWOBODYOVERLAP_I3y3z_Py_aa = I_TWOBODYOVERLAP_K4y3z_S_aa+ABY*I_TWOBODYOVERLAP_I3y3z_S_aa;
  Double I_TWOBODYOVERLAP_I2y4z_Py_aa = I_TWOBODYOVERLAP_K3y4z_S_aa+ABY*I_TWOBODYOVERLAP_I2y4z_S_aa;
  Double I_TWOBODYOVERLAP_Iy5z_Py_aa = I_TWOBODYOVERLAP_K2y5z_S_aa+ABY*I_TWOBODYOVERLAP_Iy5z_S_aa;
  Double I_TWOBODYOVERLAP_I6z_Py_aa = I_TWOBODYOVERLAP_Ky6z_S_aa+ABY*I_TWOBODYOVERLAP_I6z_S_aa;
  Double I_TWOBODYOVERLAP_I6x_Pz_aa = I_TWOBODYOVERLAP_K6xz_S_aa+ABZ*I_TWOBODYOVERLAP_I6x_S_aa;
  Double I_TWOBODYOVERLAP_I5xy_Pz_aa = I_TWOBODYOVERLAP_K5xyz_S_aa+ABZ*I_TWOBODYOVERLAP_I5xy_S_aa;
  Double I_TWOBODYOVERLAP_I5xz_Pz_aa = I_TWOBODYOVERLAP_K5x2z_S_aa+ABZ*I_TWOBODYOVERLAP_I5xz_S_aa;
  Double I_TWOBODYOVERLAP_I4x2y_Pz_aa = I_TWOBODYOVERLAP_K4x2yz_S_aa+ABZ*I_TWOBODYOVERLAP_I4x2y_S_aa;
  Double I_TWOBODYOVERLAP_I4xyz_Pz_aa = I_TWOBODYOVERLAP_K4xy2z_S_aa+ABZ*I_TWOBODYOVERLAP_I4xyz_S_aa;
  Double I_TWOBODYOVERLAP_I4x2z_Pz_aa = I_TWOBODYOVERLAP_K4x3z_S_aa+ABZ*I_TWOBODYOVERLAP_I4x2z_S_aa;
  Double I_TWOBODYOVERLAP_I3x3y_Pz_aa = I_TWOBODYOVERLAP_K3x3yz_S_aa+ABZ*I_TWOBODYOVERLAP_I3x3y_S_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_Pz_aa = I_TWOBODYOVERLAP_K3x2y2z_S_aa+ABZ*I_TWOBODYOVERLAP_I3x2yz_S_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_Pz_aa = I_TWOBODYOVERLAP_K3xy3z_S_aa+ABZ*I_TWOBODYOVERLAP_I3xy2z_S_aa;
  Double I_TWOBODYOVERLAP_I3x3z_Pz_aa = I_TWOBODYOVERLAP_K3x4z_S_aa+ABZ*I_TWOBODYOVERLAP_I3x3z_S_aa;
  Double I_TWOBODYOVERLAP_I2x4y_Pz_aa = I_TWOBODYOVERLAP_K2x4yz_S_aa+ABZ*I_TWOBODYOVERLAP_I2x4y_S_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_Pz_aa = I_TWOBODYOVERLAP_K2x3y2z_S_aa+ABZ*I_TWOBODYOVERLAP_I2x3yz_S_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_Pz_aa = I_TWOBODYOVERLAP_K2x2y3z_S_aa+ABZ*I_TWOBODYOVERLAP_I2x2y2z_S_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_Pz_aa = I_TWOBODYOVERLAP_K2xy4z_S_aa+ABZ*I_TWOBODYOVERLAP_I2xy3z_S_aa;
  Double I_TWOBODYOVERLAP_I2x4z_Pz_aa = I_TWOBODYOVERLAP_K2x5z_S_aa+ABZ*I_TWOBODYOVERLAP_I2x4z_S_aa;
  Double I_TWOBODYOVERLAP_Ix5y_Pz_aa = I_TWOBODYOVERLAP_Kx5yz_S_aa+ABZ*I_TWOBODYOVERLAP_Ix5y_S_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_Pz_aa = I_TWOBODYOVERLAP_Kx4y2z_S_aa+ABZ*I_TWOBODYOVERLAP_Ix4yz_S_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_Pz_aa = I_TWOBODYOVERLAP_Kx3y3z_S_aa+ABZ*I_TWOBODYOVERLAP_Ix3y2z_S_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_Pz_aa = I_TWOBODYOVERLAP_Kx2y4z_S_aa+ABZ*I_TWOBODYOVERLAP_Ix2y3z_S_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_Pz_aa = I_TWOBODYOVERLAP_Kxy5z_S_aa+ABZ*I_TWOBODYOVERLAP_Ixy4z_S_aa;
  Double I_TWOBODYOVERLAP_Ix5z_Pz_aa = I_TWOBODYOVERLAP_Kx6z_S_aa+ABZ*I_TWOBODYOVERLAP_Ix5z_S_aa;
  Double I_TWOBODYOVERLAP_I6y_Pz_aa = I_TWOBODYOVERLAP_K6yz_S_aa+ABZ*I_TWOBODYOVERLAP_I6y_S_aa;
  Double I_TWOBODYOVERLAP_I5yz_Pz_aa = I_TWOBODYOVERLAP_K5y2z_S_aa+ABZ*I_TWOBODYOVERLAP_I5yz_S_aa;
  Double I_TWOBODYOVERLAP_I4y2z_Pz_aa = I_TWOBODYOVERLAP_K4y3z_S_aa+ABZ*I_TWOBODYOVERLAP_I4y2z_S_aa;
  Double I_TWOBODYOVERLAP_I3y3z_Pz_aa = I_TWOBODYOVERLAP_K3y4z_S_aa+ABZ*I_TWOBODYOVERLAP_I3y3z_S_aa;
  Double I_TWOBODYOVERLAP_I2y4z_Pz_aa = I_TWOBODYOVERLAP_K2y5z_S_aa+ABZ*I_TWOBODYOVERLAP_I2y4z_S_aa;
  Double I_TWOBODYOVERLAP_Iy5z_Pz_aa = I_TWOBODYOVERLAP_Ky6z_S_aa+ABZ*I_TWOBODYOVERLAP_Iy5z_S_aa;
  Double I_TWOBODYOVERLAP_I6z_Pz_aa = I_TWOBODYOVERLAP_K7z_S_aa+ABZ*I_TWOBODYOVERLAP_I6z_S_aa;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_K_P_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_L_S_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_S_aa
   ************************************************************/
  Double I_TWOBODYOVERLAP_K7x_Px_aa = I_TWOBODYOVERLAP_L8x_S_aa+ABX*I_TWOBODYOVERLAP_K7x_S_aa;
  Double I_TWOBODYOVERLAP_K6xy_Px_aa = I_TWOBODYOVERLAP_L7xy_S_aa+ABX*I_TWOBODYOVERLAP_K6xy_S_aa;
  Double I_TWOBODYOVERLAP_K6xz_Px_aa = I_TWOBODYOVERLAP_L7xz_S_aa+ABX*I_TWOBODYOVERLAP_K6xz_S_aa;
  Double I_TWOBODYOVERLAP_K5x2y_Px_aa = I_TWOBODYOVERLAP_L6x2y_S_aa+ABX*I_TWOBODYOVERLAP_K5x2y_S_aa;
  Double I_TWOBODYOVERLAP_K5xyz_Px_aa = I_TWOBODYOVERLAP_L6xyz_S_aa+ABX*I_TWOBODYOVERLAP_K5xyz_S_aa;
  Double I_TWOBODYOVERLAP_K5x2z_Px_aa = I_TWOBODYOVERLAP_L6x2z_S_aa+ABX*I_TWOBODYOVERLAP_K5x2z_S_aa;
  Double I_TWOBODYOVERLAP_K4x3y_Px_aa = I_TWOBODYOVERLAP_L5x3y_S_aa+ABX*I_TWOBODYOVERLAP_K4x3y_S_aa;
  Double I_TWOBODYOVERLAP_K4x2yz_Px_aa = I_TWOBODYOVERLAP_L5x2yz_S_aa+ABX*I_TWOBODYOVERLAP_K4x2yz_S_aa;
  Double I_TWOBODYOVERLAP_K4xy2z_Px_aa = I_TWOBODYOVERLAP_L5xy2z_S_aa+ABX*I_TWOBODYOVERLAP_K4xy2z_S_aa;
  Double I_TWOBODYOVERLAP_K4x3z_Px_aa = I_TWOBODYOVERLAP_L5x3z_S_aa+ABX*I_TWOBODYOVERLAP_K4x3z_S_aa;
  Double I_TWOBODYOVERLAP_K3x4y_Px_aa = I_TWOBODYOVERLAP_L4x4y_S_aa+ABX*I_TWOBODYOVERLAP_K3x4y_S_aa;
  Double I_TWOBODYOVERLAP_K3x3yz_Px_aa = I_TWOBODYOVERLAP_L4x3yz_S_aa+ABX*I_TWOBODYOVERLAP_K3x3yz_S_aa;
  Double I_TWOBODYOVERLAP_K3x2y2z_Px_aa = I_TWOBODYOVERLAP_L4x2y2z_S_aa+ABX*I_TWOBODYOVERLAP_K3x2y2z_S_aa;
  Double I_TWOBODYOVERLAP_K3xy3z_Px_aa = I_TWOBODYOVERLAP_L4xy3z_S_aa+ABX*I_TWOBODYOVERLAP_K3xy3z_S_aa;
  Double I_TWOBODYOVERLAP_K3x4z_Px_aa = I_TWOBODYOVERLAP_L4x4z_S_aa+ABX*I_TWOBODYOVERLAP_K3x4z_S_aa;
  Double I_TWOBODYOVERLAP_K2x5y_Px_aa = I_TWOBODYOVERLAP_L3x5y_S_aa+ABX*I_TWOBODYOVERLAP_K2x5y_S_aa;
  Double I_TWOBODYOVERLAP_K2x4yz_Px_aa = I_TWOBODYOVERLAP_L3x4yz_S_aa+ABX*I_TWOBODYOVERLAP_K2x4yz_S_aa;
  Double I_TWOBODYOVERLAP_K2x3y2z_Px_aa = I_TWOBODYOVERLAP_L3x3y2z_S_aa+ABX*I_TWOBODYOVERLAP_K2x3y2z_S_aa;
  Double I_TWOBODYOVERLAP_K2x2y3z_Px_aa = I_TWOBODYOVERLAP_L3x2y3z_S_aa+ABX*I_TWOBODYOVERLAP_K2x2y3z_S_aa;
  Double I_TWOBODYOVERLAP_K2xy4z_Px_aa = I_TWOBODYOVERLAP_L3xy4z_S_aa+ABX*I_TWOBODYOVERLAP_K2xy4z_S_aa;
  Double I_TWOBODYOVERLAP_K2x5z_Px_aa = I_TWOBODYOVERLAP_L3x5z_S_aa+ABX*I_TWOBODYOVERLAP_K2x5z_S_aa;
  Double I_TWOBODYOVERLAP_Kx6y_Px_aa = I_TWOBODYOVERLAP_L2x6y_S_aa+ABX*I_TWOBODYOVERLAP_Kx6y_S_aa;
  Double I_TWOBODYOVERLAP_Kx5yz_Px_aa = I_TWOBODYOVERLAP_L2x5yz_S_aa+ABX*I_TWOBODYOVERLAP_Kx5yz_S_aa;
  Double I_TWOBODYOVERLAP_Kx4y2z_Px_aa = I_TWOBODYOVERLAP_L2x4y2z_S_aa+ABX*I_TWOBODYOVERLAP_Kx4y2z_S_aa;
  Double I_TWOBODYOVERLAP_Kx3y3z_Px_aa = I_TWOBODYOVERLAP_L2x3y3z_S_aa+ABX*I_TWOBODYOVERLAP_Kx3y3z_S_aa;
  Double I_TWOBODYOVERLAP_Kx2y4z_Px_aa = I_TWOBODYOVERLAP_L2x2y4z_S_aa+ABX*I_TWOBODYOVERLAP_Kx2y4z_S_aa;
  Double I_TWOBODYOVERLAP_Kxy5z_Px_aa = I_TWOBODYOVERLAP_L2xy5z_S_aa+ABX*I_TWOBODYOVERLAP_Kxy5z_S_aa;
  Double I_TWOBODYOVERLAP_Kx6z_Px_aa = I_TWOBODYOVERLAP_L2x6z_S_aa+ABX*I_TWOBODYOVERLAP_Kx6z_S_aa;
  Double I_TWOBODYOVERLAP_K7y_Px_aa = I_TWOBODYOVERLAP_Lx7y_S_aa+ABX*I_TWOBODYOVERLAP_K7y_S_aa;
  Double I_TWOBODYOVERLAP_K6yz_Px_aa = I_TWOBODYOVERLAP_Lx6yz_S_aa+ABX*I_TWOBODYOVERLAP_K6yz_S_aa;
  Double I_TWOBODYOVERLAP_K5y2z_Px_aa = I_TWOBODYOVERLAP_Lx5y2z_S_aa+ABX*I_TWOBODYOVERLAP_K5y2z_S_aa;
  Double I_TWOBODYOVERLAP_K4y3z_Px_aa = I_TWOBODYOVERLAP_Lx4y3z_S_aa+ABX*I_TWOBODYOVERLAP_K4y3z_S_aa;
  Double I_TWOBODYOVERLAP_K3y4z_Px_aa = I_TWOBODYOVERLAP_Lx3y4z_S_aa+ABX*I_TWOBODYOVERLAP_K3y4z_S_aa;
  Double I_TWOBODYOVERLAP_K2y5z_Px_aa = I_TWOBODYOVERLAP_Lx2y5z_S_aa+ABX*I_TWOBODYOVERLAP_K2y5z_S_aa;
  Double I_TWOBODYOVERLAP_Ky6z_Px_aa = I_TWOBODYOVERLAP_Lxy6z_S_aa+ABX*I_TWOBODYOVERLAP_Ky6z_S_aa;
  Double I_TWOBODYOVERLAP_K7z_Px_aa = I_TWOBODYOVERLAP_Lx7z_S_aa+ABX*I_TWOBODYOVERLAP_K7z_S_aa;
  Double I_TWOBODYOVERLAP_K7x_Py_aa = I_TWOBODYOVERLAP_L7xy_S_aa+ABY*I_TWOBODYOVERLAP_K7x_S_aa;
  Double I_TWOBODYOVERLAP_K6xy_Py_aa = I_TWOBODYOVERLAP_L6x2y_S_aa+ABY*I_TWOBODYOVERLAP_K6xy_S_aa;
  Double I_TWOBODYOVERLAP_K6xz_Py_aa = I_TWOBODYOVERLAP_L6xyz_S_aa+ABY*I_TWOBODYOVERLAP_K6xz_S_aa;
  Double I_TWOBODYOVERLAP_K5x2y_Py_aa = I_TWOBODYOVERLAP_L5x3y_S_aa+ABY*I_TWOBODYOVERLAP_K5x2y_S_aa;
  Double I_TWOBODYOVERLAP_K5xyz_Py_aa = I_TWOBODYOVERLAP_L5x2yz_S_aa+ABY*I_TWOBODYOVERLAP_K5xyz_S_aa;
  Double I_TWOBODYOVERLAP_K5x2z_Py_aa = I_TWOBODYOVERLAP_L5xy2z_S_aa+ABY*I_TWOBODYOVERLAP_K5x2z_S_aa;
  Double I_TWOBODYOVERLAP_K4x3y_Py_aa = I_TWOBODYOVERLAP_L4x4y_S_aa+ABY*I_TWOBODYOVERLAP_K4x3y_S_aa;
  Double I_TWOBODYOVERLAP_K4x2yz_Py_aa = I_TWOBODYOVERLAP_L4x3yz_S_aa+ABY*I_TWOBODYOVERLAP_K4x2yz_S_aa;
  Double I_TWOBODYOVERLAP_K4xy2z_Py_aa = I_TWOBODYOVERLAP_L4x2y2z_S_aa+ABY*I_TWOBODYOVERLAP_K4xy2z_S_aa;
  Double I_TWOBODYOVERLAP_K4x3z_Py_aa = I_TWOBODYOVERLAP_L4xy3z_S_aa+ABY*I_TWOBODYOVERLAP_K4x3z_S_aa;
  Double I_TWOBODYOVERLAP_K3x4y_Py_aa = I_TWOBODYOVERLAP_L3x5y_S_aa+ABY*I_TWOBODYOVERLAP_K3x4y_S_aa;
  Double I_TWOBODYOVERLAP_K3x3yz_Py_aa = I_TWOBODYOVERLAP_L3x4yz_S_aa+ABY*I_TWOBODYOVERLAP_K3x3yz_S_aa;
  Double I_TWOBODYOVERLAP_K3x2y2z_Py_aa = I_TWOBODYOVERLAP_L3x3y2z_S_aa+ABY*I_TWOBODYOVERLAP_K3x2y2z_S_aa;
  Double I_TWOBODYOVERLAP_K3xy3z_Py_aa = I_TWOBODYOVERLAP_L3x2y3z_S_aa+ABY*I_TWOBODYOVERLAP_K3xy3z_S_aa;
  Double I_TWOBODYOVERLAP_K3x4z_Py_aa = I_TWOBODYOVERLAP_L3xy4z_S_aa+ABY*I_TWOBODYOVERLAP_K3x4z_S_aa;
  Double I_TWOBODYOVERLAP_K2x5y_Py_aa = I_TWOBODYOVERLAP_L2x6y_S_aa+ABY*I_TWOBODYOVERLAP_K2x5y_S_aa;
  Double I_TWOBODYOVERLAP_K2x4yz_Py_aa = I_TWOBODYOVERLAP_L2x5yz_S_aa+ABY*I_TWOBODYOVERLAP_K2x4yz_S_aa;
  Double I_TWOBODYOVERLAP_K2x3y2z_Py_aa = I_TWOBODYOVERLAP_L2x4y2z_S_aa+ABY*I_TWOBODYOVERLAP_K2x3y2z_S_aa;
  Double I_TWOBODYOVERLAP_K2x2y3z_Py_aa = I_TWOBODYOVERLAP_L2x3y3z_S_aa+ABY*I_TWOBODYOVERLAP_K2x2y3z_S_aa;
  Double I_TWOBODYOVERLAP_K2xy4z_Py_aa = I_TWOBODYOVERLAP_L2x2y4z_S_aa+ABY*I_TWOBODYOVERLAP_K2xy4z_S_aa;
  Double I_TWOBODYOVERLAP_K2x5z_Py_aa = I_TWOBODYOVERLAP_L2xy5z_S_aa+ABY*I_TWOBODYOVERLAP_K2x5z_S_aa;
  Double I_TWOBODYOVERLAP_Kx6y_Py_aa = I_TWOBODYOVERLAP_Lx7y_S_aa+ABY*I_TWOBODYOVERLAP_Kx6y_S_aa;
  Double I_TWOBODYOVERLAP_Kx5yz_Py_aa = I_TWOBODYOVERLAP_Lx6yz_S_aa+ABY*I_TWOBODYOVERLAP_Kx5yz_S_aa;
  Double I_TWOBODYOVERLAP_Kx4y2z_Py_aa = I_TWOBODYOVERLAP_Lx5y2z_S_aa+ABY*I_TWOBODYOVERLAP_Kx4y2z_S_aa;
  Double I_TWOBODYOVERLAP_Kx3y3z_Py_aa = I_TWOBODYOVERLAP_Lx4y3z_S_aa+ABY*I_TWOBODYOVERLAP_Kx3y3z_S_aa;
  Double I_TWOBODYOVERLAP_Kx2y4z_Py_aa = I_TWOBODYOVERLAP_Lx3y4z_S_aa+ABY*I_TWOBODYOVERLAP_Kx2y4z_S_aa;
  Double I_TWOBODYOVERLAP_Kxy5z_Py_aa = I_TWOBODYOVERLAP_Lx2y5z_S_aa+ABY*I_TWOBODYOVERLAP_Kxy5z_S_aa;
  Double I_TWOBODYOVERLAP_Kx6z_Py_aa = I_TWOBODYOVERLAP_Lxy6z_S_aa+ABY*I_TWOBODYOVERLAP_Kx6z_S_aa;
  Double I_TWOBODYOVERLAP_K7y_Py_aa = I_TWOBODYOVERLAP_L8y_S_aa+ABY*I_TWOBODYOVERLAP_K7y_S_aa;
  Double I_TWOBODYOVERLAP_K6yz_Py_aa = I_TWOBODYOVERLAP_L7yz_S_aa+ABY*I_TWOBODYOVERLAP_K6yz_S_aa;
  Double I_TWOBODYOVERLAP_K5y2z_Py_aa = I_TWOBODYOVERLAP_L6y2z_S_aa+ABY*I_TWOBODYOVERLAP_K5y2z_S_aa;
  Double I_TWOBODYOVERLAP_K4y3z_Py_aa = I_TWOBODYOVERLAP_L5y3z_S_aa+ABY*I_TWOBODYOVERLAP_K4y3z_S_aa;
  Double I_TWOBODYOVERLAP_K3y4z_Py_aa = I_TWOBODYOVERLAP_L4y4z_S_aa+ABY*I_TWOBODYOVERLAP_K3y4z_S_aa;
  Double I_TWOBODYOVERLAP_K2y5z_Py_aa = I_TWOBODYOVERLAP_L3y5z_S_aa+ABY*I_TWOBODYOVERLAP_K2y5z_S_aa;
  Double I_TWOBODYOVERLAP_Ky6z_Py_aa = I_TWOBODYOVERLAP_L2y6z_S_aa+ABY*I_TWOBODYOVERLAP_Ky6z_S_aa;
  Double I_TWOBODYOVERLAP_K7z_Py_aa = I_TWOBODYOVERLAP_Ly7z_S_aa+ABY*I_TWOBODYOVERLAP_K7z_S_aa;
  Double I_TWOBODYOVERLAP_K7x_Pz_aa = I_TWOBODYOVERLAP_L7xz_S_aa+ABZ*I_TWOBODYOVERLAP_K7x_S_aa;
  Double I_TWOBODYOVERLAP_K6xy_Pz_aa = I_TWOBODYOVERLAP_L6xyz_S_aa+ABZ*I_TWOBODYOVERLAP_K6xy_S_aa;
  Double I_TWOBODYOVERLAP_K6xz_Pz_aa = I_TWOBODYOVERLAP_L6x2z_S_aa+ABZ*I_TWOBODYOVERLAP_K6xz_S_aa;
  Double I_TWOBODYOVERLAP_K5x2y_Pz_aa = I_TWOBODYOVERLAP_L5x2yz_S_aa+ABZ*I_TWOBODYOVERLAP_K5x2y_S_aa;
  Double I_TWOBODYOVERLAP_K5xyz_Pz_aa = I_TWOBODYOVERLAP_L5xy2z_S_aa+ABZ*I_TWOBODYOVERLAP_K5xyz_S_aa;
  Double I_TWOBODYOVERLAP_K5x2z_Pz_aa = I_TWOBODYOVERLAP_L5x3z_S_aa+ABZ*I_TWOBODYOVERLAP_K5x2z_S_aa;
  Double I_TWOBODYOVERLAP_K4x3y_Pz_aa = I_TWOBODYOVERLAP_L4x3yz_S_aa+ABZ*I_TWOBODYOVERLAP_K4x3y_S_aa;
  Double I_TWOBODYOVERLAP_K4x2yz_Pz_aa = I_TWOBODYOVERLAP_L4x2y2z_S_aa+ABZ*I_TWOBODYOVERLAP_K4x2yz_S_aa;
  Double I_TWOBODYOVERLAP_K4xy2z_Pz_aa = I_TWOBODYOVERLAP_L4xy3z_S_aa+ABZ*I_TWOBODYOVERLAP_K4xy2z_S_aa;
  Double I_TWOBODYOVERLAP_K4x3z_Pz_aa = I_TWOBODYOVERLAP_L4x4z_S_aa+ABZ*I_TWOBODYOVERLAP_K4x3z_S_aa;
  Double I_TWOBODYOVERLAP_K3x4y_Pz_aa = I_TWOBODYOVERLAP_L3x4yz_S_aa+ABZ*I_TWOBODYOVERLAP_K3x4y_S_aa;
  Double I_TWOBODYOVERLAP_K3x3yz_Pz_aa = I_TWOBODYOVERLAP_L3x3y2z_S_aa+ABZ*I_TWOBODYOVERLAP_K3x3yz_S_aa;
  Double I_TWOBODYOVERLAP_K3x2y2z_Pz_aa = I_TWOBODYOVERLAP_L3x2y3z_S_aa+ABZ*I_TWOBODYOVERLAP_K3x2y2z_S_aa;
  Double I_TWOBODYOVERLAP_K3xy3z_Pz_aa = I_TWOBODYOVERLAP_L3xy4z_S_aa+ABZ*I_TWOBODYOVERLAP_K3xy3z_S_aa;
  Double I_TWOBODYOVERLAP_K3x4z_Pz_aa = I_TWOBODYOVERLAP_L3x5z_S_aa+ABZ*I_TWOBODYOVERLAP_K3x4z_S_aa;
  Double I_TWOBODYOVERLAP_K2x5y_Pz_aa = I_TWOBODYOVERLAP_L2x5yz_S_aa+ABZ*I_TWOBODYOVERLAP_K2x5y_S_aa;
  Double I_TWOBODYOVERLAP_K2x4yz_Pz_aa = I_TWOBODYOVERLAP_L2x4y2z_S_aa+ABZ*I_TWOBODYOVERLAP_K2x4yz_S_aa;
  Double I_TWOBODYOVERLAP_K2x3y2z_Pz_aa = I_TWOBODYOVERLAP_L2x3y3z_S_aa+ABZ*I_TWOBODYOVERLAP_K2x3y2z_S_aa;
  Double I_TWOBODYOVERLAP_K2x2y3z_Pz_aa = I_TWOBODYOVERLAP_L2x2y4z_S_aa+ABZ*I_TWOBODYOVERLAP_K2x2y3z_S_aa;
  Double I_TWOBODYOVERLAP_K2xy4z_Pz_aa = I_TWOBODYOVERLAP_L2xy5z_S_aa+ABZ*I_TWOBODYOVERLAP_K2xy4z_S_aa;
  Double I_TWOBODYOVERLAP_K2x5z_Pz_aa = I_TWOBODYOVERLAP_L2x6z_S_aa+ABZ*I_TWOBODYOVERLAP_K2x5z_S_aa;
  Double I_TWOBODYOVERLAP_Kx6y_Pz_aa = I_TWOBODYOVERLAP_Lx6yz_S_aa+ABZ*I_TWOBODYOVERLAP_Kx6y_S_aa;
  Double I_TWOBODYOVERLAP_Kx5yz_Pz_aa = I_TWOBODYOVERLAP_Lx5y2z_S_aa+ABZ*I_TWOBODYOVERLAP_Kx5yz_S_aa;
  Double I_TWOBODYOVERLAP_Kx4y2z_Pz_aa = I_TWOBODYOVERLAP_Lx4y3z_S_aa+ABZ*I_TWOBODYOVERLAP_Kx4y2z_S_aa;
  Double I_TWOBODYOVERLAP_Kx3y3z_Pz_aa = I_TWOBODYOVERLAP_Lx3y4z_S_aa+ABZ*I_TWOBODYOVERLAP_Kx3y3z_S_aa;
  Double I_TWOBODYOVERLAP_Kx2y4z_Pz_aa = I_TWOBODYOVERLAP_Lx2y5z_S_aa+ABZ*I_TWOBODYOVERLAP_Kx2y4z_S_aa;
  Double I_TWOBODYOVERLAP_Kxy5z_Pz_aa = I_TWOBODYOVERLAP_Lxy6z_S_aa+ABZ*I_TWOBODYOVERLAP_Kxy5z_S_aa;
  Double I_TWOBODYOVERLAP_Kx6z_Pz_aa = I_TWOBODYOVERLAP_Lx7z_S_aa+ABZ*I_TWOBODYOVERLAP_Kx6z_S_aa;
  Double I_TWOBODYOVERLAP_K7y_Pz_aa = I_TWOBODYOVERLAP_L7yz_S_aa+ABZ*I_TWOBODYOVERLAP_K7y_S_aa;
  Double I_TWOBODYOVERLAP_K6yz_Pz_aa = I_TWOBODYOVERLAP_L6y2z_S_aa+ABZ*I_TWOBODYOVERLAP_K6yz_S_aa;
  Double I_TWOBODYOVERLAP_K5y2z_Pz_aa = I_TWOBODYOVERLAP_L5y3z_S_aa+ABZ*I_TWOBODYOVERLAP_K5y2z_S_aa;
  Double I_TWOBODYOVERLAP_K4y3z_Pz_aa = I_TWOBODYOVERLAP_L4y4z_S_aa+ABZ*I_TWOBODYOVERLAP_K4y3z_S_aa;
  Double I_TWOBODYOVERLAP_K3y4z_Pz_aa = I_TWOBODYOVERLAP_L3y5z_S_aa+ABZ*I_TWOBODYOVERLAP_K3y4z_S_aa;
  Double I_TWOBODYOVERLAP_K2y5z_Pz_aa = I_TWOBODYOVERLAP_L2y6z_S_aa+ABZ*I_TWOBODYOVERLAP_K2y5z_S_aa;
  Double I_TWOBODYOVERLAP_Ky6z_Pz_aa = I_TWOBODYOVERLAP_Ly7z_S_aa+ABZ*I_TWOBODYOVERLAP_Ky6z_S_aa;
  Double I_TWOBODYOVERLAP_K7z_Pz_aa = I_TWOBODYOVERLAP_L8z_S_aa+ABZ*I_TWOBODYOVERLAP_K7z_S_aa;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_I_D_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 56 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_P_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_P_aa
   ************************************************************/
  Double I_TWOBODYOVERLAP_I6x_D2x_aa = I_TWOBODYOVERLAP_K7x_Px_aa+ABX*I_TWOBODYOVERLAP_I6x_Px_aa;
  Double I_TWOBODYOVERLAP_I5xy_D2x_aa = I_TWOBODYOVERLAP_K6xy_Px_aa+ABX*I_TWOBODYOVERLAP_I5xy_Px_aa;
  Double I_TWOBODYOVERLAP_I5xz_D2x_aa = I_TWOBODYOVERLAP_K6xz_Px_aa+ABX*I_TWOBODYOVERLAP_I5xz_Px_aa;
  Double I_TWOBODYOVERLAP_I4x2y_D2x_aa = I_TWOBODYOVERLAP_K5x2y_Px_aa+ABX*I_TWOBODYOVERLAP_I4x2y_Px_aa;
  Double I_TWOBODYOVERLAP_I4xyz_D2x_aa = I_TWOBODYOVERLAP_K5xyz_Px_aa+ABX*I_TWOBODYOVERLAP_I4xyz_Px_aa;
  Double I_TWOBODYOVERLAP_I4x2z_D2x_aa = I_TWOBODYOVERLAP_K5x2z_Px_aa+ABX*I_TWOBODYOVERLAP_I4x2z_Px_aa;
  Double I_TWOBODYOVERLAP_I3x3y_D2x_aa = I_TWOBODYOVERLAP_K4x3y_Px_aa+ABX*I_TWOBODYOVERLAP_I3x3y_Px_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_D2x_aa = I_TWOBODYOVERLAP_K4x2yz_Px_aa+ABX*I_TWOBODYOVERLAP_I3x2yz_Px_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_D2x_aa = I_TWOBODYOVERLAP_K4xy2z_Px_aa+ABX*I_TWOBODYOVERLAP_I3xy2z_Px_aa;
  Double I_TWOBODYOVERLAP_I3x3z_D2x_aa = I_TWOBODYOVERLAP_K4x3z_Px_aa+ABX*I_TWOBODYOVERLAP_I3x3z_Px_aa;
  Double I_TWOBODYOVERLAP_I2x4y_D2x_aa = I_TWOBODYOVERLAP_K3x4y_Px_aa+ABX*I_TWOBODYOVERLAP_I2x4y_Px_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_D2x_aa = I_TWOBODYOVERLAP_K3x3yz_Px_aa+ABX*I_TWOBODYOVERLAP_I2x3yz_Px_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_D2x_aa = I_TWOBODYOVERLAP_K3x2y2z_Px_aa+ABX*I_TWOBODYOVERLAP_I2x2y2z_Px_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_D2x_aa = I_TWOBODYOVERLAP_K3xy3z_Px_aa+ABX*I_TWOBODYOVERLAP_I2xy3z_Px_aa;
  Double I_TWOBODYOVERLAP_I2x4z_D2x_aa = I_TWOBODYOVERLAP_K3x4z_Px_aa+ABX*I_TWOBODYOVERLAP_I2x4z_Px_aa;
  Double I_TWOBODYOVERLAP_Ix5y_D2x_aa = I_TWOBODYOVERLAP_K2x5y_Px_aa+ABX*I_TWOBODYOVERLAP_Ix5y_Px_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_D2x_aa = I_TWOBODYOVERLAP_K2x4yz_Px_aa+ABX*I_TWOBODYOVERLAP_Ix4yz_Px_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_D2x_aa = I_TWOBODYOVERLAP_K2x3y2z_Px_aa+ABX*I_TWOBODYOVERLAP_Ix3y2z_Px_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_D2x_aa = I_TWOBODYOVERLAP_K2x2y3z_Px_aa+ABX*I_TWOBODYOVERLAP_Ix2y3z_Px_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_D2x_aa = I_TWOBODYOVERLAP_K2xy4z_Px_aa+ABX*I_TWOBODYOVERLAP_Ixy4z_Px_aa;
  Double I_TWOBODYOVERLAP_Ix5z_D2x_aa = I_TWOBODYOVERLAP_K2x5z_Px_aa+ABX*I_TWOBODYOVERLAP_Ix5z_Px_aa;
  Double I_TWOBODYOVERLAP_I6y_D2x_aa = I_TWOBODYOVERLAP_Kx6y_Px_aa+ABX*I_TWOBODYOVERLAP_I6y_Px_aa;
  Double I_TWOBODYOVERLAP_I5yz_D2x_aa = I_TWOBODYOVERLAP_Kx5yz_Px_aa+ABX*I_TWOBODYOVERLAP_I5yz_Px_aa;
  Double I_TWOBODYOVERLAP_I4y2z_D2x_aa = I_TWOBODYOVERLAP_Kx4y2z_Px_aa+ABX*I_TWOBODYOVERLAP_I4y2z_Px_aa;
  Double I_TWOBODYOVERLAP_I3y3z_D2x_aa = I_TWOBODYOVERLAP_Kx3y3z_Px_aa+ABX*I_TWOBODYOVERLAP_I3y3z_Px_aa;
  Double I_TWOBODYOVERLAP_I2y4z_D2x_aa = I_TWOBODYOVERLAP_Kx2y4z_Px_aa+ABX*I_TWOBODYOVERLAP_I2y4z_Px_aa;
  Double I_TWOBODYOVERLAP_Iy5z_D2x_aa = I_TWOBODYOVERLAP_Kxy5z_Px_aa+ABX*I_TWOBODYOVERLAP_Iy5z_Px_aa;
  Double I_TWOBODYOVERLAP_I6z_D2x_aa = I_TWOBODYOVERLAP_Kx6z_Px_aa+ABX*I_TWOBODYOVERLAP_I6z_Px_aa;
  Double I_TWOBODYOVERLAP_I6x_Dxy_aa = I_TWOBODYOVERLAP_K6xy_Px_aa+ABY*I_TWOBODYOVERLAP_I6x_Px_aa;
  Double I_TWOBODYOVERLAP_I5xy_Dxy_aa = I_TWOBODYOVERLAP_K5x2y_Px_aa+ABY*I_TWOBODYOVERLAP_I5xy_Px_aa;
  Double I_TWOBODYOVERLAP_I5xz_Dxy_aa = I_TWOBODYOVERLAP_K5xyz_Px_aa+ABY*I_TWOBODYOVERLAP_I5xz_Px_aa;
  Double I_TWOBODYOVERLAP_I4x2y_Dxy_aa = I_TWOBODYOVERLAP_K4x3y_Px_aa+ABY*I_TWOBODYOVERLAP_I4x2y_Px_aa;
  Double I_TWOBODYOVERLAP_I4xyz_Dxy_aa = I_TWOBODYOVERLAP_K4x2yz_Px_aa+ABY*I_TWOBODYOVERLAP_I4xyz_Px_aa;
  Double I_TWOBODYOVERLAP_I4x2z_Dxy_aa = I_TWOBODYOVERLAP_K4xy2z_Px_aa+ABY*I_TWOBODYOVERLAP_I4x2z_Px_aa;
  Double I_TWOBODYOVERLAP_I3x3y_Dxy_aa = I_TWOBODYOVERLAP_K3x4y_Px_aa+ABY*I_TWOBODYOVERLAP_I3x3y_Px_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_Dxy_aa = I_TWOBODYOVERLAP_K3x3yz_Px_aa+ABY*I_TWOBODYOVERLAP_I3x2yz_Px_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_Dxy_aa = I_TWOBODYOVERLAP_K3x2y2z_Px_aa+ABY*I_TWOBODYOVERLAP_I3xy2z_Px_aa;
  Double I_TWOBODYOVERLAP_I3x3z_Dxy_aa = I_TWOBODYOVERLAP_K3xy3z_Px_aa+ABY*I_TWOBODYOVERLAP_I3x3z_Px_aa;
  Double I_TWOBODYOVERLAP_I2x4y_Dxy_aa = I_TWOBODYOVERLAP_K2x5y_Px_aa+ABY*I_TWOBODYOVERLAP_I2x4y_Px_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_Dxy_aa = I_TWOBODYOVERLAP_K2x4yz_Px_aa+ABY*I_TWOBODYOVERLAP_I2x3yz_Px_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_Dxy_aa = I_TWOBODYOVERLAP_K2x3y2z_Px_aa+ABY*I_TWOBODYOVERLAP_I2x2y2z_Px_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_Dxy_aa = I_TWOBODYOVERLAP_K2x2y3z_Px_aa+ABY*I_TWOBODYOVERLAP_I2xy3z_Px_aa;
  Double I_TWOBODYOVERLAP_I2x4z_Dxy_aa = I_TWOBODYOVERLAP_K2xy4z_Px_aa+ABY*I_TWOBODYOVERLAP_I2x4z_Px_aa;
  Double I_TWOBODYOVERLAP_Ix5y_Dxy_aa = I_TWOBODYOVERLAP_Kx6y_Px_aa+ABY*I_TWOBODYOVERLAP_Ix5y_Px_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_Dxy_aa = I_TWOBODYOVERLAP_Kx5yz_Px_aa+ABY*I_TWOBODYOVERLAP_Ix4yz_Px_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_Dxy_aa = I_TWOBODYOVERLAP_Kx4y2z_Px_aa+ABY*I_TWOBODYOVERLAP_Ix3y2z_Px_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_Dxy_aa = I_TWOBODYOVERLAP_Kx3y3z_Px_aa+ABY*I_TWOBODYOVERLAP_Ix2y3z_Px_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_Dxy_aa = I_TWOBODYOVERLAP_Kx2y4z_Px_aa+ABY*I_TWOBODYOVERLAP_Ixy4z_Px_aa;
  Double I_TWOBODYOVERLAP_Ix5z_Dxy_aa = I_TWOBODYOVERLAP_Kxy5z_Px_aa+ABY*I_TWOBODYOVERLAP_Ix5z_Px_aa;
  Double I_TWOBODYOVERLAP_I6y_Dxy_aa = I_TWOBODYOVERLAP_K7y_Px_aa+ABY*I_TWOBODYOVERLAP_I6y_Px_aa;
  Double I_TWOBODYOVERLAP_I5yz_Dxy_aa = I_TWOBODYOVERLAP_K6yz_Px_aa+ABY*I_TWOBODYOVERLAP_I5yz_Px_aa;
  Double I_TWOBODYOVERLAP_I4y2z_Dxy_aa = I_TWOBODYOVERLAP_K5y2z_Px_aa+ABY*I_TWOBODYOVERLAP_I4y2z_Px_aa;
  Double I_TWOBODYOVERLAP_I3y3z_Dxy_aa = I_TWOBODYOVERLAP_K4y3z_Px_aa+ABY*I_TWOBODYOVERLAP_I3y3z_Px_aa;
  Double I_TWOBODYOVERLAP_I2y4z_Dxy_aa = I_TWOBODYOVERLAP_K3y4z_Px_aa+ABY*I_TWOBODYOVERLAP_I2y4z_Px_aa;
  Double I_TWOBODYOVERLAP_Iy5z_Dxy_aa = I_TWOBODYOVERLAP_K2y5z_Px_aa+ABY*I_TWOBODYOVERLAP_Iy5z_Px_aa;
  Double I_TWOBODYOVERLAP_I6z_Dxy_aa = I_TWOBODYOVERLAP_Ky6z_Px_aa+ABY*I_TWOBODYOVERLAP_I6z_Px_aa;
  Double I_TWOBODYOVERLAP_I6x_D2y_aa = I_TWOBODYOVERLAP_K6xy_Py_aa+ABY*I_TWOBODYOVERLAP_I6x_Py_aa;
  Double I_TWOBODYOVERLAP_I5xy_D2y_aa = I_TWOBODYOVERLAP_K5x2y_Py_aa+ABY*I_TWOBODYOVERLAP_I5xy_Py_aa;
  Double I_TWOBODYOVERLAP_I5xz_D2y_aa = I_TWOBODYOVERLAP_K5xyz_Py_aa+ABY*I_TWOBODYOVERLAP_I5xz_Py_aa;
  Double I_TWOBODYOVERLAP_I4x2y_D2y_aa = I_TWOBODYOVERLAP_K4x3y_Py_aa+ABY*I_TWOBODYOVERLAP_I4x2y_Py_aa;
  Double I_TWOBODYOVERLAP_I4xyz_D2y_aa = I_TWOBODYOVERLAP_K4x2yz_Py_aa+ABY*I_TWOBODYOVERLAP_I4xyz_Py_aa;
  Double I_TWOBODYOVERLAP_I4x2z_D2y_aa = I_TWOBODYOVERLAP_K4xy2z_Py_aa+ABY*I_TWOBODYOVERLAP_I4x2z_Py_aa;
  Double I_TWOBODYOVERLAP_I3x3y_D2y_aa = I_TWOBODYOVERLAP_K3x4y_Py_aa+ABY*I_TWOBODYOVERLAP_I3x3y_Py_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_D2y_aa = I_TWOBODYOVERLAP_K3x3yz_Py_aa+ABY*I_TWOBODYOVERLAP_I3x2yz_Py_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_D2y_aa = I_TWOBODYOVERLAP_K3x2y2z_Py_aa+ABY*I_TWOBODYOVERLAP_I3xy2z_Py_aa;
  Double I_TWOBODYOVERLAP_I3x3z_D2y_aa = I_TWOBODYOVERLAP_K3xy3z_Py_aa+ABY*I_TWOBODYOVERLAP_I3x3z_Py_aa;
  Double I_TWOBODYOVERLAP_I2x4y_D2y_aa = I_TWOBODYOVERLAP_K2x5y_Py_aa+ABY*I_TWOBODYOVERLAP_I2x4y_Py_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_D2y_aa = I_TWOBODYOVERLAP_K2x4yz_Py_aa+ABY*I_TWOBODYOVERLAP_I2x3yz_Py_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_D2y_aa = I_TWOBODYOVERLAP_K2x3y2z_Py_aa+ABY*I_TWOBODYOVERLAP_I2x2y2z_Py_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_D2y_aa = I_TWOBODYOVERLAP_K2x2y3z_Py_aa+ABY*I_TWOBODYOVERLAP_I2xy3z_Py_aa;
  Double I_TWOBODYOVERLAP_I2x4z_D2y_aa = I_TWOBODYOVERLAP_K2xy4z_Py_aa+ABY*I_TWOBODYOVERLAP_I2x4z_Py_aa;
  Double I_TWOBODYOVERLAP_Ix5y_D2y_aa = I_TWOBODYOVERLAP_Kx6y_Py_aa+ABY*I_TWOBODYOVERLAP_Ix5y_Py_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_D2y_aa = I_TWOBODYOVERLAP_Kx5yz_Py_aa+ABY*I_TWOBODYOVERLAP_Ix4yz_Py_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_D2y_aa = I_TWOBODYOVERLAP_Kx4y2z_Py_aa+ABY*I_TWOBODYOVERLAP_Ix3y2z_Py_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_D2y_aa = I_TWOBODYOVERLAP_Kx3y3z_Py_aa+ABY*I_TWOBODYOVERLAP_Ix2y3z_Py_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_D2y_aa = I_TWOBODYOVERLAP_Kx2y4z_Py_aa+ABY*I_TWOBODYOVERLAP_Ixy4z_Py_aa;
  Double I_TWOBODYOVERLAP_Ix5z_D2y_aa = I_TWOBODYOVERLAP_Kxy5z_Py_aa+ABY*I_TWOBODYOVERLAP_Ix5z_Py_aa;
  Double I_TWOBODYOVERLAP_I6y_D2y_aa = I_TWOBODYOVERLAP_K7y_Py_aa+ABY*I_TWOBODYOVERLAP_I6y_Py_aa;
  Double I_TWOBODYOVERLAP_I5yz_D2y_aa = I_TWOBODYOVERLAP_K6yz_Py_aa+ABY*I_TWOBODYOVERLAP_I5yz_Py_aa;
  Double I_TWOBODYOVERLAP_I4y2z_D2y_aa = I_TWOBODYOVERLAP_K5y2z_Py_aa+ABY*I_TWOBODYOVERLAP_I4y2z_Py_aa;
  Double I_TWOBODYOVERLAP_I3y3z_D2y_aa = I_TWOBODYOVERLAP_K4y3z_Py_aa+ABY*I_TWOBODYOVERLAP_I3y3z_Py_aa;
  Double I_TWOBODYOVERLAP_I2y4z_D2y_aa = I_TWOBODYOVERLAP_K3y4z_Py_aa+ABY*I_TWOBODYOVERLAP_I2y4z_Py_aa;
  Double I_TWOBODYOVERLAP_Iy5z_D2y_aa = I_TWOBODYOVERLAP_K2y5z_Py_aa+ABY*I_TWOBODYOVERLAP_Iy5z_Py_aa;
  Double I_TWOBODYOVERLAP_I6z_D2y_aa = I_TWOBODYOVERLAP_Ky6z_Py_aa+ABY*I_TWOBODYOVERLAP_I6z_Py_aa;
  Double I_TWOBODYOVERLAP_I6x_D2z_aa = I_TWOBODYOVERLAP_K6xz_Pz_aa+ABZ*I_TWOBODYOVERLAP_I6x_Pz_aa;
  Double I_TWOBODYOVERLAP_I5xy_D2z_aa = I_TWOBODYOVERLAP_K5xyz_Pz_aa+ABZ*I_TWOBODYOVERLAP_I5xy_Pz_aa;
  Double I_TWOBODYOVERLAP_I5xz_D2z_aa = I_TWOBODYOVERLAP_K5x2z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I5xz_Pz_aa;
  Double I_TWOBODYOVERLAP_I4x2y_D2z_aa = I_TWOBODYOVERLAP_K4x2yz_Pz_aa+ABZ*I_TWOBODYOVERLAP_I4x2y_Pz_aa;
  Double I_TWOBODYOVERLAP_I4xyz_D2z_aa = I_TWOBODYOVERLAP_K4xy2z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I4xyz_Pz_aa;
  Double I_TWOBODYOVERLAP_I4x2z_D2z_aa = I_TWOBODYOVERLAP_K4x3z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I4x2z_Pz_aa;
  Double I_TWOBODYOVERLAP_I3x3y_D2z_aa = I_TWOBODYOVERLAP_K3x3yz_Pz_aa+ABZ*I_TWOBODYOVERLAP_I3x3y_Pz_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_D2z_aa = I_TWOBODYOVERLAP_K3x2y2z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I3x2yz_Pz_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_D2z_aa = I_TWOBODYOVERLAP_K3xy3z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I3xy2z_Pz_aa;
  Double I_TWOBODYOVERLAP_I3x3z_D2z_aa = I_TWOBODYOVERLAP_K3x4z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I3x3z_Pz_aa;
  Double I_TWOBODYOVERLAP_I2x4y_D2z_aa = I_TWOBODYOVERLAP_K2x4yz_Pz_aa+ABZ*I_TWOBODYOVERLAP_I2x4y_Pz_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_D2z_aa = I_TWOBODYOVERLAP_K2x3y2z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I2x3yz_Pz_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_D2z_aa = I_TWOBODYOVERLAP_K2x2y3z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I2x2y2z_Pz_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_D2z_aa = I_TWOBODYOVERLAP_K2xy4z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I2xy3z_Pz_aa;
  Double I_TWOBODYOVERLAP_I2x4z_D2z_aa = I_TWOBODYOVERLAP_K2x5z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I2x4z_Pz_aa;
  Double I_TWOBODYOVERLAP_Ix5y_D2z_aa = I_TWOBODYOVERLAP_Kx5yz_Pz_aa+ABZ*I_TWOBODYOVERLAP_Ix5y_Pz_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_D2z_aa = I_TWOBODYOVERLAP_Kx4y2z_Pz_aa+ABZ*I_TWOBODYOVERLAP_Ix4yz_Pz_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_D2z_aa = I_TWOBODYOVERLAP_Kx3y3z_Pz_aa+ABZ*I_TWOBODYOVERLAP_Ix3y2z_Pz_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_D2z_aa = I_TWOBODYOVERLAP_Kx2y4z_Pz_aa+ABZ*I_TWOBODYOVERLAP_Ix2y3z_Pz_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_D2z_aa = I_TWOBODYOVERLAP_Kxy5z_Pz_aa+ABZ*I_TWOBODYOVERLAP_Ixy4z_Pz_aa;
  Double I_TWOBODYOVERLAP_Ix5z_D2z_aa = I_TWOBODYOVERLAP_Kx6z_Pz_aa+ABZ*I_TWOBODYOVERLAP_Ix5z_Pz_aa;
  Double I_TWOBODYOVERLAP_I6y_D2z_aa = I_TWOBODYOVERLAP_K6yz_Pz_aa+ABZ*I_TWOBODYOVERLAP_I6y_Pz_aa;
  Double I_TWOBODYOVERLAP_I5yz_D2z_aa = I_TWOBODYOVERLAP_K5y2z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I5yz_Pz_aa;
  Double I_TWOBODYOVERLAP_I4y2z_D2z_aa = I_TWOBODYOVERLAP_K4y3z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I4y2z_Pz_aa;
  Double I_TWOBODYOVERLAP_I3y3z_D2z_aa = I_TWOBODYOVERLAP_K3y4z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I3y3z_Pz_aa;
  Double I_TWOBODYOVERLAP_I2y4z_D2z_aa = I_TWOBODYOVERLAP_K2y5z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I2y4z_Pz_aa;
  Double I_TWOBODYOVERLAP_Iy5z_D2z_aa = I_TWOBODYOVERLAP_Ky6z_Pz_aa+ABZ*I_TWOBODYOVERLAP_Iy5z_Pz_aa;
  Double I_TWOBODYOVERLAP_I6z_D2z_aa = I_TWOBODYOVERLAP_K7z_Pz_aa+ABZ*I_TWOBODYOVERLAP_I6z_Pz_aa;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_L_P_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 20 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_M_S_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_L_S_aa
   ************************************************************/
  Double I_TWOBODYOVERLAP_L8x_Px_aa = I_TWOBODYOVERLAP_M9x_S_aa+ABX*I_TWOBODYOVERLAP_L8x_S_aa;
  Double I_TWOBODYOVERLAP_L7xy_Px_aa = I_TWOBODYOVERLAP_M8xy_S_aa+ABX*I_TWOBODYOVERLAP_L7xy_S_aa;
  Double I_TWOBODYOVERLAP_L7xz_Px_aa = I_TWOBODYOVERLAP_M8xz_S_aa+ABX*I_TWOBODYOVERLAP_L7xz_S_aa;
  Double I_TWOBODYOVERLAP_L6x2y_Px_aa = I_TWOBODYOVERLAP_M7x2y_S_aa+ABX*I_TWOBODYOVERLAP_L6x2y_S_aa;
  Double I_TWOBODYOVERLAP_L6xyz_Px_aa = I_TWOBODYOVERLAP_M7xyz_S_aa+ABX*I_TWOBODYOVERLAP_L6xyz_S_aa;
  Double I_TWOBODYOVERLAP_L6x2z_Px_aa = I_TWOBODYOVERLAP_M7x2z_S_aa+ABX*I_TWOBODYOVERLAP_L6x2z_S_aa;
  Double I_TWOBODYOVERLAP_L5x3y_Px_aa = I_TWOBODYOVERLAP_M6x3y_S_aa+ABX*I_TWOBODYOVERLAP_L5x3y_S_aa;
  Double I_TWOBODYOVERLAP_L5x2yz_Px_aa = I_TWOBODYOVERLAP_M6x2yz_S_aa+ABX*I_TWOBODYOVERLAP_L5x2yz_S_aa;
  Double I_TWOBODYOVERLAP_L5xy2z_Px_aa = I_TWOBODYOVERLAP_M6xy2z_S_aa+ABX*I_TWOBODYOVERLAP_L5xy2z_S_aa;
  Double I_TWOBODYOVERLAP_L5x3z_Px_aa = I_TWOBODYOVERLAP_M6x3z_S_aa+ABX*I_TWOBODYOVERLAP_L5x3z_S_aa;
  Double I_TWOBODYOVERLAP_L4x4y_Px_aa = I_TWOBODYOVERLAP_M5x4y_S_aa+ABX*I_TWOBODYOVERLAP_L4x4y_S_aa;
  Double I_TWOBODYOVERLAP_L4x3yz_Px_aa = I_TWOBODYOVERLAP_M5x3yz_S_aa+ABX*I_TWOBODYOVERLAP_L4x3yz_S_aa;
  Double I_TWOBODYOVERLAP_L4x2y2z_Px_aa = I_TWOBODYOVERLAP_M5x2y2z_S_aa+ABX*I_TWOBODYOVERLAP_L4x2y2z_S_aa;
  Double I_TWOBODYOVERLAP_L4xy3z_Px_aa = I_TWOBODYOVERLAP_M5xy3z_S_aa+ABX*I_TWOBODYOVERLAP_L4xy3z_S_aa;
  Double I_TWOBODYOVERLAP_L4x4z_Px_aa = I_TWOBODYOVERLAP_M5x4z_S_aa+ABX*I_TWOBODYOVERLAP_L4x4z_S_aa;
  Double I_TWOBODYOVERLAP_L3x5y_Px_aa = I_TWOBODYOVERLAP_M4x5y_S_aa+ABX*I_TWOBODYOVERLAP_L3x5y_S_aa;
  Double I_TWOBODYOVERLAP_L3x4yz_Px_aa = I_TWOBODYOVERLAP_M4x4yz_S_aa+ABX*I_TWOBODYOVERLAP_L3x4yz_S_aa;
  Double I_TWOBODYOVERLAP_L3x3y2z_Px_aa = I_TWOBODYOVERLAP_M4x3y2z_S_aa+ABX*I_TWOBODYOVERLAP_L3x3y2z_S_aa;
  Double I_TWOBODYOVERLAP_L3x2y3z_Px_aa = I_TWOBODYOVERLAP_M4x2y3z_S_aa+ABX*I_TWOBODYOVERLAP_L3x2y3z_S_aa;
  Double I_TWOBODYOVERLAP_L3xy4z_Px_aa = I_TWOBODYOVERLAP_M4xy4z_S_aa+ABX*I_TWOBODYOVERLAP_L3xy4z_S_aa;
  Double I_TWOBODYOVERLAP_L3x5z_Px_aa = I_TWOBODYOVERLAP_M4x5z_S_aa+ABX*I_TWOBODYOVERLAP_L3x5z_S_aa;
  Double I_TWOBODYOVERLAP_L2x6y_Px_aa = I_TWOBODYOVERLAP_M3x6y_S_aa+ABX*I_TWOBODYOVERLAP_L2x6y_S_aa;
  Double I_TWOBODYOVERLAP_L2x5yz_Px_aa = I_TWOBODYOVERLAP_M3x5yz_S_aa+ABX*I_TWOBODYOVERLAP_L2x5yz_S_aa;
  Double I_TWOBODYOVERLAP_L2x4y2z_Px_aa = I_TWOBODYOVERLAP_M3x4y2z_S_aa+ABX*I_TWOBODYOVERLAP_L2x4y2z_S_aa;
  Double I_TWOBODYOVERLAP_L2x3y3z_Px_aa = I_TWOBODYOVERLAP_M3x3y3z_S_aa+ABX*I_TWOBODYOVERLAP_L2x3y3z_S_aa;
  Double I_TWOBODYOVERLAP_L2x2y4z_Px_aa = I_TWOBODYOVERLAP_M3x2y4z_S_aa+ABX*I_TWOBODYOVERLAP_L2x2y4z_S_aa;
  Double I_TWOBODYOVERLAP_L2xy5z_Px_aa = I_TWOBODYOVERLAP_M3xy5z_S_aa+ABX*I_TWOBODYOVERLAP_L2xy5z_S_aa;
  Double I_TWOBODYOVERLAP_L2x6z_Px_aa = I_TWOBODYOVERLAP_M3x6z_S_aa+ABX*I_TWOBODYOVERLAP_L2x6z_S_aa;
  Double I_TWOBODYOVERLAP_Lx7y_Px_aa = I_TWOBODYOVERLAP_M2x7y_S_aa+ABX*I_TWOBODYOVERLAP_Lx7y_S_aa;
  Double I_TWOBODYOVERLAP_Lx6yz_Px_aa = I_TWOBODYOVERLAP_M2x6yz_S_aa+ABX*I_TWOBODYOVERLAP_Lx6yz_S_aa;
  Double I_TWOBODYOVERLAP_Lx5y2z_Px_aa = I_TWOBODYOVERLAP_M2x5y2z_S_aa+ABX*I_TWOBODYOVERLAP_Lx5y2z_S_aa;
  Double I_TWOBODYOVERLAP_Lx4y3z_Px_aa = I_TWOBODYOVERLAP_M2x4y3z_S_aa+ABX*I_TWOBODYOVERLAP_Lx4y3z_S_aa;
  Double I_TWOBODYOVERLAP_Lx3y4z_Px_aa = I_TWOBODYOVERLAP_M2x3y4z_S_aa+ABX*I_TWOBODYOVERLAP_Lx3y4z_S_aa;
  Double I_TWOBODYOVERLAP_Lx2y5z_Px_aa = I_TWOBODYOVERLAP_M2x2y5z_S_aa+ABX*I_TWOBODYOVERLAP_Lx2y5z_S_aa;
  Double I_TWOBODYOVERLAP_Lxy6z_Px_aa = I_TWOBODYOVERLAP_M2xy6z_S_aa+ABX*I_TWOBODYOVERLAP_Lxy6z_S_aa;
  Double I_TWOBODYOVERLAP_Lx7z_Px_aa = I_TWOBODYOVERLAP_M2x7z_S_aa+ABX*I_TWOBODYOVERLAP_Lx7z_S_aa;
  Double I_TWOBODYOVERLAP_L7yz_Px_aa = I_TWOBODYOVERLAP_Mx7yz_S_aa+ABX*I_TWOBODYOVERLAP_L7yz_S_aa;
  Double I_TWOBODYOVERLAP_L6y2z_Px_aa = I_TWOBODYOVERLAP_Mx6y2z_S_aa+ABX*I_TWOBODYOVERLAP_L6y2z_S_aa;
  Double I_TWOBODYOVERLAP_L5y3z_Px_aa = I_TWOBODYOVERLAP_Mx5y3z_S_aa+ABX*I_TWOBODYOVERLAP_L5y3z_S_aa;
  Double I_TWOBODYOVERLAP_L4y4z_Px_aa = I_TWOBODYOVERLAP_Mx4y4z_S_aa+ABX*I_TWOBODYOVERLAP_L4y4z_S_aa;
  Double I_TWOBODYOVERLAP_L3y5z_Px_aa = I_TWOBODYOVERLAP_Mx3y5z_S_aa+ABX*I_TWOBODYOVERLAP_L3y5z_S_aa;
  Double I_TWOBODYOVERLAP_L2y6z_Px_aa = I_TWOBODYOVERLAP_Mx2y6z_S_aa+ABX*I_TWOBODYOVERLAP_L2y6z_S_aa;
  Double I_TWOBODYOVERLAP_Ly7z_Px_aa = I_TWOBODYOVERLAP_Mxy7z_S_aa+ABX*I_TWOBODYOVERLAP_Ly7z_S_aa;
  Double I_TWOBODYOVERLAP_L7xy_Py_aa = I_TWOBODYOVERLAP_M7x2y_S_aa+ABY*I_TWOBODYOVERLAP_L7xy_S_aa;
  Double I_TWOBODYOVERLAP_L6x2y_Py_aa = I_TWOBODYOVERLAP_M6x3y_S_aa+ABY*I_TWOBODYOVERLAP_L6x2y_S_aa;
  Double I_TWOBODYOVERLAP_L6xyz_Py_aa = I_TWOBODYOVERLAP_M6x2yz_S_aa+ABY*I_TWOBODYOVERLAP_L6xyz_S_aa;
  Double I_TWOBODYOVERLAP_L5x3y_Py_aa = I_TWOBODYOVERLAP_M5x4y_S_aa+ABY*I_TWOBODYOVERLAP_L5x3y_S_aa;
  Double I_TWOBODYOVERLAP_L5x2yz_Py_aa = I_TWOBODYOVERLAP_M5x3yz_S_aa+ABY*I_TWOBODYOVERLAP_L5x2yz_S_aa;
  Double I_TWOBODYOVERLAP_L5xy2z_Py_aa = I_TWOBODYOVERLAP_M5x2y2z_S_aa+ABY*I_TWOBODYOVERLAP_L5xy2z_S_aa;
  Double I_TWOBODYOVERLAP_L4x4y_Py_aa = I_TWOBODYOVERLAP_M4x5y_S_aa+ABY*I_TWOBODYOVERLAP_L4x4y_S_aa;
  Double I_TWOBODYOVERLAP_L4x3yz_Py_aa = I_TWOBODYOVERLAP_M4x4yz_S_aa+ABY*I_TWOBODYOVERLAP_L4x3yz_S_aa;
  Double I_TWOBODYOVERLAP_L4x2y2z_Py_aa = I_TWOBODYOVERLAP_M4x3y2z_S_aa+ABY*I_TWOBODYOVERLAP_L4x2y2z_S_aa;
  Double I_TWOBODYOVERLAP_L4xy3z_Py_aa = I_TWOBODYOVERLAP_M4x2y3z_S_aa+ABY*I_TWOBODYOVERLAP_L4xy3z_S_aa;
  Double I_TWOBODYOVERLAP_L3x5y_Py_aa = I_TWOBODYOVERLAP_M3x6y_S_aa+ABY*I_TWOBODYOVERLAP_L3x5y_S_aa;
  Double I_TWOBODYOVERLAP_L3x4yz_Py_aa = I_TWOBODYOVERLAP_M3x5yz_S_aa+ABY*I_TWOBODYOVERLAP_L3x4yz_S_aa;
  Double I_TWOBODYOVERLAP_L3x3y2z_Py_aa = I_TWOBODYOVERLAP_M3x4y2z_S_aa+ABY*I_TWOBODYOVERLAP_L3x3y2z_S_aa;
  Double I_TWOBODYOVERLAP_L3x2y3z_Py_aa = I_TWOBODYOVERLAP_M3x3y3z_S_aa+ABY*I_TWOBODYOVERLAP_L3x2y3z_S_aa;
  Double I_TWOBODYOVERLAP_L3xy4z_Py_aa = I_TWOBODYOVERLAP_M3x2y4z_S_aa+ABY*I_TWOBODYOVERLAP_L3xy4z_S_aa;
  Double I_TWOBODYOVERLAP_L2x6y_Py_aa = I_TWOBODYOVERLAP_M2x7y_S_aa+ABY*I_TWOBODYOVERLAP_L2x6y_S_aa;
  Double I_TWOBODYOVERLAP_L2x5yz_Py_aa = I_TWOBODYOVERLAP_M2x6yz_S_aa+ABY*I_TWOBODYOVERLAP_L2x5yz_S_aa;
  Double I_TWOBODYOVERLAP_L2x4y2z_Py_aa = I_TWOBODYOVERLAP_M2x5y2z_S_aa+ABY*I_TWOBODYOVERLAP_L2x4y2z_S_aa;
  Double I_TWOBODYOVERLAP_L2x3y3z_Py_aa = I_TWOBODYOVERLAP_M2x4y3z_S_aa+ABY*I_TWOBODYOVERLAP_L2x3y3z_S_aa;
  Double I_TWOBODYOVERLAP_L2x2y4z_Py_aa = I_TWOBODYOVERLAP_M2x3y4z_S_aa+ABY*I_TWOBODYOVERLAP_L2x2y4z_S_aa;
  Double I_TWOBODYOVERLAP_L2xy5z_Py_aa = I_TWOBODYOVERLAP_M2x2y5z_S_aa+ABY*I_TWOBODYOVERLAP_L2xy5z_S_aa;
  Double I_TWOBODYOVERLAP_Lx7y_Py_aa = I_TWOBODYOVERLAP_Mx8y_S_aa+ABY*I_TWOBODYOVERLAP_Lx7y_S_aa;
  Double I_TWOBODYOVERLAP_Lx6yz_Py_aa = I_TWOBODYOVERLAP_Mx7yz_S_aa+ABY*I_TWOBODYOVERLAP_Lx6yz_S_aa;
  Double I_TWOBODYOVERLAP_Lx5y2z_Py_aa = I_TWOBODYOVERLAP_Mx6y2z_S_aa+ABY*I_TWOBODYOVERLAP_Lx5y2z_S_aa;
  Double I_TWOBODYOVERLAP_Lx4y3z_Py_aa = I_TWOBODYOVERLAP_Mx5y3z_S_aa+ABY*I_TWOBODYOVERLAP_Lx4y3z_S_aa;
  Double I_TWOBODYOVERLAP_Lx3y4z_Py_aa = I_TWOBODYOVERLAP_Mx4y4z_S_aa+ABY*I_TWOBODYOVERLAP_Lx3y4z_S_aa;
  Double I_TWOBODYOVERLAP_Lx2y5z_Py_aa = I_TWOBODYOVERLAP_Mx3y5z_S_aa+ABY*I_TWOBODYOVERLAP_Lx2y5z_S_aa;
  Double I_TWOBODYOVERLAP_Lxy6z_Py_aa = I_TWOBODYOVERLAP_Mx2y6z_S_aa+ABY*I_TWOBODYOVERLAP_Lxy6z_S_aa;
  Double I_TWOBODYOVERLAP_L8y_Py_aa = I_TWOBODYOVERLAP_M9y_S_aa+ABY*I_TWOBODYOVERLAP_L8y_S_aa;
  Double I_TWOBODYOVERLAP_L7yz_Py_aa = I_TWOBODYOVERLAP_M8yz_S_aa+ABY*I_TWOBODYOVERLAP_L7yz_S_aa;
  Double I_TWOBODYOVERLAP_L6y2z_Py_aa = I_TWOBODYOVERLAP_M7y2z_S_aa+ABY*I_TWOBODYOVERLAP_L6y2z_S_aa;
  Double I_TWOBODYOVERLAP_L5y3z_Py_aa = I_TWOBODYOVERLAP_M6y3z_S_aa+ABY*I_TWOBODYOVERLAP_L5y3z_S_aa;
  Double I_TWOBODYOVERLAP_L4y4z_Py_aa = I_TWOBODYOVERLAP_M5y4z_S_aa+ABY*I_TWOBODYOVERLAP_L4y4z_S_aa;
  Double I_TWOBODYOVERLAP_L3y5z_Py_aa = I_TWOBODYOVERLAP_M4y5z_S_aa+ABY*I_TWOBODYOVERLAP_L3y5z_S_aa;
  Double I_TWOBODYOVERLAP_L2y6z_Py_aa = I_TWOBODYOVERLAP_M3y6z_S_aa+ABY*I_TWOBODYOVERLAP_L2y6z_S_aa;
  Double I_TWOBODYOVERLAP_Ly7z_Py_aa = I_TWOBODYOVERLAP_M2y7z_S_aa+ABY*I_TWOBODYOVERLAP_Ly7z_S_aa;
  Double I_TWOBODYOVERLAP_L7xz_Pz_aa = I_TWOBODYOVERLAP_M7x2z_S_aa+ABZ*I_TWOBODYOVERLAP_L7xz_S_aa;
  Double I_TWOBODYOVERLAP_L6xyz_Pz_aa = I_TWOBODYOVERLAP_M6xy2z_S_aa+ABZ*I_TWOBODYOVERLAP_L6xyz_S_aa;
  Double I_TWOBODYOVERLAP_L6x2z_Pz_aa = I_TWOBODYOVERLAP_M6x3z_S_aa+ABZ*I_TWOBODYOVERLAP_L6x2z_S_aa;
  Double I_TWOBODYOVERLAP_L5x2yz_Pz_aa = I_TWOBODYOVERLAP_M5x2y2z_S_aa+ABZ*I_TWOBODYOVERLAP_L5x2yz_S_aa;
  Double I_TWOBODYOVERLAP_L5xy2z_Pz_aa = I_TWOBODYOVERLAP_M5xy3z_S_aa+ABZ*I_TWOBODYOVERLAP_L5xy2z_S_aa;
  Double I_TWOBODYOVERLAP_L5x3z_Pz_aa = I_TWOBODYOVERLAP_M5x4z_S_aa+ABZ*I_TWOBODYOVERLAP_L5x3z_S_aa;
  Double I_TWOBODYOVERLAP_L4x3yz_Pz_aa = I_TWOBODYOVERLAP_M4x3y2z_S_aa+ABZ*I_TWOBODYOVERLAP_L4x3yz_S_aa;
  Double I_TWOBODYOVERLAP_L4x2y2z_Pz_aa = I_TWOBODYOVERLAP_M4x2y3z_S_aa+ABZ*I_TWOBODYOVERLAP_L4x2y2z_S_aa;
  Double I_TWOBODYOVERLAP_L4xy3z_Pz_aa = I_TWOBODYOVERLAP_M4xy4z_S_aa+ABZ*I_TWOBODYOVERLAP_L4xy3z_S_aa;
  Double I_TWOBODYOVERLAP_L4x4z_Pz_aa = I_TWOBODYOVERLAP_M4x5z_S_aa+ABZ*I_TWOBODYOVERLAP_L4x4z_S_aa;
  Double I_TWOBODYOVERLAP_L3x4yz_Pz_aa = I_TWOBODYOVERLAP_M3x4y2z_S_aa+ABZ*I_TWOBODYOVERLAP_L3x4yz_S_aa;
  Double I_TWOBODYOVERLAP_L3x3y2z_Pz_aa = I_TWOBODYOVERLAP_M3x3y3z_S_aa+ABZ*I_TWOBODYOVERLAP_L3x3y2z_S_aa;
  Double I_TWOBODYOVERLAP_L3x2y3z_Pz_aa = I_TWOBODYOVERLAP_M3x2y4z_S_aa+ABZ*I_TWOBODYOVERLAP_L3x2y3z_S_aa;
  Double I_TWOBODYOVERLAP_L3xy4z_Pz_aa = I_TWOBODYOVERLAP_M3xy5z_S_aa+ABZ*I_TWOBODYOVERLAP_L3xy4z_S_aa;
  Double I_TWOBODYOVERLAP_L3x5z_Pz_aa = I_TWOBODYOVERLAP_M3x6z_S_aa+ABZ*I_TWOBODYOVERLAP_L3x5z_S_aa;
  Double I_TWOBODYOVERLAP_L2x5yz_Pz_aa = I_TWOBODYOVERLAP_M2x5y2z_S_aa+ABZ*I_TWOBODYOVERLAP_L2x5yz_S_aa;
  Double I_TWOBODYOVERLAP_L2x4y2z_Pz_aa = I_TWOBODYOVERLAP_M2x4y3z_S_aa+ABZ*I_TWOBODYOVERLAP_L2x4y2z_S_aa;
  Double I_TWOBODYOVERLAP_L2x3y3z_Pz_aa = I_TWOBODYOVERLAP_M2x3y4z_S_aa+ABZ*I_TWOBODYOVERLAP_L2x3y3z_S_aa;
  Double I_TWOBODYOVERLAP_L2x2y4z_Pz_aa = I_TWOBODYOVERLAP_M2x2y5z_S_aa+ABZ*I_TWOBODYOVERLAP_L2x2y4z_S_aa;
  Double I_TWOBODYOVERLAP_L2xy5z_Pz_aa = I_TWOBODYOVERLAP_M2xy6z_S_aa+ABZ*I_TWOBODYOVERLAP_L2xy5z_S_aa;
  Double I_TWOBODYOVERLAP_L2x6z_Pz_aa = I_TWOBODYOVERLAP_M2x7z_S_aa+ABZ*I_TWOBODYOVERLAP_L2x6z_S_aa;
  Double I_TWOBODYOVERLAP_Lx6yz_Pz_aa = I_TWOBODYOVERLAP_Mx6y2z_S_aa+ABZ*I_TWOBODYOVERLAP_Lx6yz_S_aa;
  Double I_TWOBODYOVERLAP_Lx5y2z_Pz_aa = I_TWOBODYOVERLAP_Mx5y3z_S_aa+ABZ*I_TWOBODYOVERLAP_Lx5y2z_S_aa;
  Double I_TWOBODYOVERLAP_Lx4y3z_Pz_aa = I_TWOBODYOVERLAP_Mx4y4z_S_aa+ABZ*I_TWOBODYOVERLAP_Lx4y3z_S_aa;
  Double I_TWOBODYOVERLAP_Lx3y4z_Pz_aa = I_TWOBODYOVERLAP_Mx3y5z_S_aa+ABZ*I_TWOBODYOVERLAP_Lx3y4z_S_aa;
  Double I_TWOBODYOVERLAP_Lx2y5z_Pz_aa = I_TWOBODYOVERLAP_Mx2y6z_S_aa+ABZ*I_TWOBODYOVERLAP_Lx2y5z_S_aa;
  Double I_TWOBODYOVERLAP_Lxy6z_Pz_aa = I_TWOBODYOVERLAP_Mxy7z_S_aa+ABZ*I_TWOBODYOVERLAP_Lxy6z_S_aa;
  Double I_TWOBODYOVERLAP_Lx7z_Pz_aa = I_TWOBODYOVERLAP_Mx8z_S_aa+ABZ*I_TWOBODYOVERLAP_Lx7z_S_aa;
  Double I_TWOBODYOVERLAP_L7yz_Pz_aa = I_TWOBODYOVERLAP_M7y2z_S_aa+ABZ*I_TWOBODYOVERLAP_L7yz_S_aa;
  Double I_TWOBODYOVERLAP_L6y2z_Pz_aa = I_TWOBODYOVERLAP_M6y3z_S_aa+ABZ*I_TWOBODYOVERLAP_L6y2z_S_aa;
  Double I_TWOBODYOVERLAP_L5y3z_Pz_aa = I_TWOBODYOVERLAP_M5y4z_S_aa+ABZ*I_TWOBODYOVERLAP_L5y3z_S_aa;
  Double I_TWOBODYOVERLAP_L4y4z_Pz_aa = I_TWOBODYOVERLAP_M4y5z_S_aa+ABZ*I_TWOBODYOVERLAP_L4y4z_S_aa;
  Double I_TWOBODYOVERLAP_L3y5z_Pz_aa = I_TWOBODYOVERLAP_M3y6z_S_aa+ABZ*I_TWOBODYOVERLAP_L3y5z_S_aa;
  Double I_TWOBODYOVERLAP_L2y6z_Pz_aa = I_TWOBODYOVERLAP_M2y7z_S_aa+ABZ*I_TWOBODYOVERLAP_L2y6z_S_aa;
  Double I_TWOBODYOVERLAP_Ly7z_Pz_aa = I_TWOBODYOVERLAP_My8z_S_aa+ABZ*I_TWOBODYOVERLAP_Ly7z_S_aa;
  Double I_TWOBODYOVERLAP_L8z_Pz_aa = I_TWOBODYOVERLAP_M9z_S_aa+ABZ*I_TWOBODYOVERLAP_L8z_S_aa;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_K_D_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 80 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_L_P_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_P_aa
   ************************************************************/
  Double I_TWOBODYOVERLAP_K7x_D2x_aa = I_TWOBODYOVERLAP_L8x_Px_aa+ABX*I_TWOBODYOVERLAP_K7x_Px_aa;
  Double I_TWOBODYOVERLAP_K6xy_D2x_aa = I_TWOBODYOVERLAP_L7xy_Px_aa+ABX*I_TWOBODYOVERLAP_K6xy_Px_aa;
  Double I_TWOBODYOVERLAP_K6xz_D2x_aa = I_TWOBODYOVERLAP_L7xz_Px_aa+ABX*I_TWOBODYOVERLAP_K6xz_Px_aa;
  Double I_TWOBODYOVERLAP_K5x2y_D2x_aa = I_TWOBODYOVERLAP_L6x2y_Px_aa+ABX*I_TWOBODYOVERLAP_K5x2y_Px_aa;
  Double I_TWOBODYOVERLAP_K5xyz_D2x_aa = I_TWOBODYOVERLAP_L6xyz_Px_aa+ABX*I_TWOBODYOVERLAP_K5xyz_Px_aa;
  Double I_TWOBODYOVERLAP_K5x2z_D2x_aa = I_TWOBODYOVERLAP_L6x2z_Px_aa+ABX*I_TWOBODYOVERLAP_K5x2z_Px_aa;
  Double I_TWOBODYOVERLAP_K4x3y_D2x_aa = I_TWOBODYOVERLAP_L5x3y_Px_aa+ABX*I_TWOBODYOVERLAP_K4x3y_Px_aa;
  Double I_TWOBODYOVERLAP_K4x2yz_D2x_aa = I_TWOBODYOVERLAP_L5x2yz_Px_aa+ABX*I_TWOBODYOVERLAP_K4x2yz_Px_aa;
  Double I_TWOBODYOVERLAP_K4xy2z_D2x_aa = I_TWOBODYOVERLAP_L5xy2z_Px_aa+ABX*I_TWOBODYOVERLAP_K4xy2z_Px_aa;
  Double I_TWOBODYOVERLAP_K4x3z_D2x_aa = I_TWOBODYOVERLAP_L5x3z_Px_aa+ABX*I_TWOBODYOVERLAP_K4x3z_Px_aa;
  Double I_TWOBODYOVERLAP_K3x4y_D2x_aa = I_TWOBODYOVERLAP_L4x4y_Px_aa+ABX*I_TWOBODYOVERLAP_K3x4y_Px_aa;
  Double I_TWOBODYOVERLAP_K3x3yz_D2x_aa = I_TWOBODYOVERLAP_L4x3yz_Px_aa+ABX*I_TWOBODYOVERLAP_K3x3yz_Px_aa;
  Double I_TWOBODYOVERLAP_K3x2y2z_D2x_aa = I_TWOBODYOVERLAP_L4x2y2z_Px_aa+ABX*I_TWOBODYOVERLAP_K3x2y2z_Px_aa;
  Double I_TWOBODYOVERLAP_K3xy3z_D2x_aa = I_TWOBODYOVERLAP_L4xy3z_Px_aa+ABX*I_TWOBODYOVERLAP_K3xy3z_Px_aa;
  Double I_TWOBODYOVERLAP_K3x4z_D2x_aa = I_TWOBODYOVERLAP_L4x4z_Px_aa+ABX*I_TWOBODYOVERLAP_K3x4z_Px_aa;
  Double I_TWOBODYOVERLAP_K2x5y_D2x_aa = I_TWOBODYOVERLAP_L3x5y_Px_aa+ABX*I_TWOBODYOVERLAP_K2x5y_Px_aa;
  Double I_TWOBODYOVERLAP_K2x4yz_D2x_aa = I_TWOBODYOVERLAP_L3x4yz_Px_aa+ABX*I_TWOBODYOVERLAP_K2x4yz_Px_aa;
  Double I_TWOBODYOVERLAP_K2x3y2z_D2x_aa = I_TWOBODYOVERLAP_L3x3y2z_Px_aa+ABX*I_TWOBODYOVERLAP_K2x3y2z_Px_aa;
  Double I_TWOBODYOVERLAP_K2x2y3z_D2x_aa = I_TWOBODYOVERLAP_L3x2y3z_Px_aa+ABX*I_TWOBODYOVERLAP_K2x2y3z_Px_aa;
  Double I_TWOBODYOVERLAP_K2xy4z_D2x_aa = I_TWOBODYOVERLAP_L3xy4z_Px_aa+ABX*I_TWOBODYOVERLAP_K2xy4z_Px_aa;
  Double I_TWOBODYOVERLAP_K2x5z_D2x_aa = I_TWOBODYOVERLAP_L3x5z_Px_aa+ABX*I_TWOBODYOVERLAP_K2x5z_Px_aa;
  Double I_TWOBODYOVERLAP_Kx6y_D2x_aa = I_TWOBODYOVERLAP_L2x6y_Px_aa+ABX*I_TWOBODYOVERLAP_Kx6y_Px_aa;
  Double I_TWOBODYOVERLAP_Kx5yz_D2x_aa = I_TWOBODYOVERLAP_L2x5yz_Px_aa+ABX*I_TWOBODYOVERLAP_Kx5yz_Px_aa;
  Double I_TWOBODYOVERLAP_Kx4y2z_D2x_aa = I_TWOBODYOVERLAP_L2x4y2z_Px_aa+ABX*I_TWOBODYOVERLAP_Kx4y2z_Px_aa;
  Double I_TWOBODYOVERLAP_Kx3y3z_D2x_aa = I_TWOBODYOVERLAP_L2x3y3z_Px_aa+ABX*I_TWOBODYOVERLAP_Kx3y3z_Px_aa;
  Double I_TWOBODYOVERLAP_Kx2y4z_D2x_aa = I_TWOBODYOVERLAP_L2x2y4z_Px_aa+ABX*I_TWOBODYOVERLAP_Kx2y4z_Px_aa;
  Double I_TWOBODYOVERLAP_Kxy5z_D2x_aa = I_TWOBODYOVERLAP_L2xy5z_Px_aa+ABX*I_TWOBODYOVERLAP_Kxy5z_Px_aa;
  Double I_TWOBODYOVERLAP_Kx6z_D2x_aa = I_TWOBODYOVERLAP_L2x6z_Px_aa+ABX*I_TWOBODYOVERLAP_Kx6z_Px_aa;
  Double I_TWOBODYOVERLAP_K7y_D2x_aa = I_TWOBODYOVERLAP_Lx7y_Px_aa+ABX*I_TWOBODYOVERLAP_K7y_Px_aa;
  Double I_TWOBODYOVERLAP_K6yz_D2x_aa = I_TWOBODYOVERLAP_Lx6yz_Px_aa+ABX*I_TWOBODYOVERLAP_K6yz_Px_aa;
  Double I_TWOBODYOVERLAP_K5y2z_D2x_aa = I_TWOBODYOVERLAP_Lx5y2z_Px_aa+ABX*I_TWOBODYOVERLAP_K5y2z_Px_aa;
  Double I_TWOBODYOVERLAP_K4y3z_D2x_aa = I_TWOBODYOVERLAP_Lx4y3z_Px_aa+ABX*I_TWOBODYOVERLAP_K4y3z_Px_aa;
  Double I_TWOBODYOVERLAP_K3y4z_D2x_aa = I_TWOBODYOVERLAP_Lx3y4z_Px_aa+ABX*I_TWOBODYOVERLAP_K3y4z_Px_aa;
  Double I_TWOBODYOVERLAP_K2y5z_D2x_aa = I_TWOBODYOVERLAP_Lx2y5z_Px_aa+ABX*I_TWOBODYOVERLAP_K2y5z_Px_aa;
  Double I_TWOBODYOVERLAP_Ky6z_D2x_aa = I_TWOBODYOVERLAP_Lxy6z_Px_aa+ABX*I_TWOBODYOVERLAP_Ky6z_Px_aa;
  Double I_TWOBODYOVERLAP_K7z_D2x_aa = I_TWOBODYOVERLAP_Lx7z_Px_aa+ABX*I_TWOBODYOVERLAP_K7z_Px_aa;
  Double I_TWOBODYOVERLAP_K6xz_Dxy_aa = I_TWOBODYOVERLAP_L6xyz_Px_aa+ABY*I_TWOBODYOVERLAP_K6xz_Px_aa;
  Double I_TWOBODYOVERLAP_K5xyz_Dxy_aa = I_TWOBODYOVERLAP_L5x2yz_Px_aa+ABY*I_TWOBODYOVERLAP_K5xyz_Px_aa;
  Double I_TWOBODYOVERLAP_K5x2z_Dxy_aa = I_TWOBODYOVERLAP_L5xy2z_Px_aa+ABY*I_TWOBODYOVERLAP_K5x2z_Px_aa;
  Double I_TWOBODYOVERLAP_K4x2yz_Dxy_aa = I_TWOBODYOVERLAP_L4x3yz_Px_aa+ABY*I_TWOBODYOVERLAP_K4x2yz_Px_aa;
  Double I_TWOBODYOVERLAP_K4xy2z_Dxy_aa = I_TWOBODYOVERLAP_L4x2y2z_Px_aa+ABY*I_TWOBODYOVERLAP_K4xy2z_Px_aa;
  Double I_TWOBODYOVERLAP_K4x3z_Dxy_aa = I_TWOBODYOVERLAP_L4xy3z_Px_aa+ABY*I_TWOBODYOVERLAP_K4x3z_Px_aa;
  Double I_TWOBODYOVERLAP_K3x3yz_Dxy_aa = I_TWOBODYOVERLAP_L3x4yz_Px_aa+ABY*I_TWOBODYOVERLAP_K3x3yz_Px_aa;
  Double I_TWOBODYOVERLAP_K3x2y2z_Dxy_aa = I_TWOBODYOVERLAP_L3x3y2z_Px_aa+ABY*I_TWOBODYOVERLAP_K3x2y2z_Px_aa;
  Double I_TWOBODYOVERLAP_K3xy3z_Dxy_aa = I_TWOBODYOVERLAP_L3x2y3z_Px_aa+ABY*I_TWOBODYOVERLAP_K3xy3z_Px_aa;
  Double I_TWOBODYOVERLAP_K3x4z_Dxy_aa = I_TWOBODYOVERLAP_L3xy4z_Px_aa+ABY*I_TWOBODYOVERLAP_K3x4z_Px_aa;
  Double I_TWOBODYOVERLAP_K2x4yz_Dxy_aa = I_TWOBODYOVERLAP_L2x5yz_Px_aa+ABY*I_TWOBODYOVERLAP_K2x4yz_Px_aa;
  Double I_TWOBODYOVERLAP_K2x3y2z_Dxy_aa = I_TWOBODYOVERLAP_L2x4y2z_Px_aa+ABY*I_TWOBODYOVERLAP_K2x3y2z_Px_aa;
  Double I_TWOBODYOVERLAP_K2x2y3z_Dxy_aa = I_TWOBODYOVERLAP_L2x3y3z_Px_aa+ABY*I_TWOBODYOVERLAP_K2x2y3z_Px_aa;
  Double I_TWOBODYOVERLAP_K2xy4z_Dxy_aa = I_TWOBODYOVERLAP_L2x2y4z_Px_aa+ABY*I_TWOBODYOVERLAP_K2xy4z_Px_aa;
  Double I_TWOBODYOVERLAP_K2x5z_Dxy_aa = I_TWOBODYOVERLAP_L2xy5z_Px_aa+ABY*I_TWOBODYOVERLAP_K2x5z_Px_aa;
  Double I_TWOBODYOVERLAP_Kx5yz_Dxy_aa = I_TWOBODYOVERLAP_Lx6yz_Px_aa+ABY*I_TWOBODYOVERLAP_Kx5yz_Px_aa;
  Double I_TWOBODYOVERLAP_Kx4y2z_Dxy_aa = I_TWOBODYOVERLAP_Lx5y2z_Px_aa+ABY*I_TWOBODYOVERLAP_Kx4y2z_Px_aa;
  Double I_TWOBODYOVERLAP_Kx3y3z_Dxy_aa = I_TWOBODYOVERLAP_Lx4y3z_Px_aa+ABY*I_TWOBODYOVERLAP_Kx3y3z_Px_aa;
  Double I_TWOBODYOVERLAP_Kx2y4z_Dxy_aa = I_TWOBODYOVERLAP_Lx3y4z_Px_aa+ABY*I_TWOBODYOVERLAP_Kx2y4z_Px_aa;
  Double I_TWOBODYOVERLAP_Kxy5z_Dxy_aa = I_TWOBODYOVERLAP_Lx2y5z_Px_aa+ABY*I_TWOBODYOVERLAP_Kxy5z_Px_aa;
  Double I_TWOBODYOVERLAP_Kx6z_Dxy_aa = I_TWOBODYOVERLAP_Lxy6z_Px_aa+ABY*I_TWOBODYOVERLAP_Kx6z_Px_aa;
  Double I_TWOBODYOVERLAP_K6yz_Dxy_aa = I_TWOBODYOVERLAP_L7yz_Px_aa+ABY*I_TWOBODYOVERLAP_K6yz_Px_aa;
  Double I_TWOBODYOVERLAP_K5y2z_Dxy_aa = I_TWOBODYOVERLAP_L6y2z_Px_aa+ABY*I_TWOBODYOVERLAP_K5y2z_Px_aa;
  Double I_TWOBODYOVERLAP_K4y3z_Dxy_aa = I_TWOBODYOVERLAP_L5y3z_Px_aa+ABY*I_TWOBODYOVERLAP_K4y3z_Px_aa;
  Double I_TWOBODYOVERLAP_K3y4z_Dxy_aa = I_TWOBODYOVERLAP_L4y4z_Px_aa+ABY*I_TWOBODYOVERLAP_K3y4z_Px_aa;
  Double I_TWOBODYOVERLAP_K2y5z_Dxy_aa = I_TWOBODYOVERLAP_L3y5z_Px_aa+ABY*I_TWOBODYOVERLAP_K2y5z_Px_aa;
  Double I_TWOBODYOVERLAP_Ky6z_Dxy_aa = I_TWOBODYOVERLAP_L2y6z_Px_aa+ABY*I_TWOBODYOVERLAP_Ky6z_Px_aa;
  Double I_TWOBODYOVERLAP_K7z_Dxy_aa = I_TWOBODYOVERLAP_Ly7z_Px_aa+ABY*I_TWOBODYOVERLAP_K7z_Px_aa;
  Double I_TWOBODYOVERLAP_K7x_D2y_aa = I_TWOBODYOVERLAP_L7xy_Py_aa+ABY*I_TWOBODYOVERLAP_K7x_Py_aa;
  Double I_TWOBODYOVERLAP_K6xy_D2y_aa = I_TWOBODYOVERLAP_L6x2y_Py_aa+ABY*I_TWOBODYOVERLAP_K6xy_Py_aa;
  Double I_TWOBODYOVERLAP_K6xz_D2y_aa = I_TWOBODYOVERLAP_L6xyz_Py_aa+ABY*I_TWOBODYOVERLAP_K6xz_Py_aa;
  Double I_TWOBODYOVERLAP_K5x2y_D2y_aa = I_TWOBODYOVERLAP_L5x3y_Py_aa+ABY*I_TWOBODYOVERLAP_K5x2y_Py_aa;
  Double I_TWOBODYOVERLAP_K5xyz_D2y_aa = I_TWOBODYOVERLAP_L5x2yz_Py_aa+ABY*I_TWOBODYOVERLAP_K5xyz_Py_aa;
  Double I_TWOBODYOVERLAP_K5x2z_D2y_aa = I_TWOBODYOVERLAP_L5xy2z_Py_aa+ABY*I_TWOBODYOVERLAP_K5x2z_Py_aa;
  Double I_TWOBODYOVERLAP_K4x3y_D2y_aa = I_TWOBODYOVERLAP_L4x4y_Py_aa+ABY*I_TWOBODYOVERLAP_K4x3y_Py_aa;
  Double I_TWOBODYOVERLAP_K4x2yz_D2y_aa = I_TWOBODYOVERLAP_L4x3yz_Py_aa+ABY*I_TWOBODYOVERLAP_K4x2yz_Py_aa;
  Double I_TWOBODYOVERLAP_K4xy2z_D2y_aa = I_TWOBODYOVERLAP_L4x2y2z_Py_aa+ABY*I_TWOBODYOVERLAP_K4xy2z_Py_aa;
  Double I_TWOBODYOVERLAP_K4x3z_D2y_aa = I_TWOBODYOVERLAP_L4xy3z_Py_aa+ABY*I_TWOBODYOVERLAP_K4x3z_Py_aa;
  Double I_TWOBODYOVERLAP_K3x4y_D2y_aa = I_TWOBODYOVERLAP_L3x5y_Py_aa+ABY*I_TWOBODYOVERLAP_K3x4y_Py_aa;
  Double I_TWOBODYOVERLAP_K3x3yz_D2y_aa = I_TWOBODYOVERLAP_L3x4yz_Py_aa+ABY*I_TWOBODYOVERLAP_K3x3yz_Py_aa;
  Double I_TWOBODYOVERLAP_K3x2y2z_D2y_aa = I_TWOBODYOVERLAP_L3x3y2z_Py_aa+ABY*I_TWOBODYOVERLAP_K3x2y2z_Py_aa;
  Double I_TWOBODYOVERLAP_K3xy3z_D2y_aa = I_TWOBODYOVERLAP_L3x2y3z_Py_aa+ABY*I_TWOBODYOVERLAP_K3xy3z_Py_aa;
  Double I_TWOBODYOVERLAP_K3x4z_D2y_aa = I_TWOBODYOVERLAP_L3xy4z_Py_aa+ABY*I_TWOBODYOVERLAP_K3x4z_Py_aa;
  Double I_TWOBODYOVERLAP_K2x5y_D2y_aa = I_TWOBODYOVERLAP_L2x6y_Py_aa+ABY*I_TWOBODYOVERLAP_K2x5y_Py_aa;
  Double I_TWOBODYOVERLAP_K2x4yz_D2y_aa = I_TWOBODYOVERLAP_L2x5yz_Py_aa+ABY*I_TWOBODYOVERLAP_K2x4yz_Py_aa;
  Double I_TWOBODYOVERLAP_K2x3y2z_D2y_aa = I_TWOBODYOVERLAP_L2x4y2z_Py_aa+ABY*I_TWOBODYOVERLAP_K2x3y2z_Py_aa;
  Double I_TWOBODYOVERLAP_K2x2y3z_D2y_aa = I_TWOBODYOVERLAP_L2x3y3z_Py_aa+ABY*I_TWOBODYOVERLAP_K2x2y3z_Py_aa;
  Double I_TWOBODYOVERLAP_K2xy4z_D2y_aa = I_TWOBODYOVERLAP_L2x2y4z_Py_aa+ABY*I_TWOBODYOVERLAP_K2xy4z_Py_aa;
  Double I_TWOBODYOVERLAP_K2x5z_D2y_aa = I_TWOBODYOVERLAP_L2xy5z_Py_aa+ABY*I_TWOBODYOVERLAP_K2x5z_Py_aa;
  Double I_TWOBODYOVERLAP_Kx6y_D2y_aa = I_TWOBODYOVERLAP_Lx7y_Py_aa+ABY*I_TWOBODYOVERLAP_Kx6y_Py_aa;
  Double I_TWOBODYOVERLAP_Kx5yz_D2y_aa = I_TWOBODYOVERLAP_Lx6yz_Py_aa+ABY*I_TWOBODYOVERLAP_Kx5yz_Py_aa;
  Double I_TWOBODYOVERLAP_Kx4y2z_D2y_aa = I_TWOBODYOVERLAP_Lx5y2z_Py_aa+ABY*I_TWOBODYOVERLAP_Kx4y2z_Py_aa;
  Double I_TWOBODYOVERLAP_Kx3y3z_D2y_aa = I_TWOBODYOVERLAP_Lx4y3z_Py_aa+ABY*I_TWOBODYOVERLAP_Kx3y3z_Py_aa;
  Double I_TWOBODYOVERLAP_Kx2y4z_D2y_aa = I_TWOBODYOVERLAP_Lx3y4z_Py_aa+ABY*I_TWOBODYOVERLAP_Kx2y4z_Py_aa;
  Double I_TWOBODYOVERLAP_Kxy5z_D2y_aa = I_TWOBODYOVERLAP_Lx2y5z_Py_aa+ABY*I_TWOBODYOVERLAP_Kxy5z_Py_aa;
  Double I_TWOBODYOVERLAP_Kx6z_D2y_aa = I_TWOBODYOVERLAP_Lxy6z_Py_aa+ABY*I_TWOBODYOVERLAP_Kx6z_Py_aa;
  Double I_TWOBODYOVERLAP_K7y_D2y_aa = I_TWOBODYOVERLAP_L8y_Py_aa+ABY*I_TWOBODYOVERLAP_K7y_Py_aa;
  Double I_TWOBODYOVERLAP_K6yz_D2y_aa = I_TWOBODYOVERLAP_L7yz_Py_aa+ABY*I_TWOBODYOVERLAP_K6yz_Py_aa;
  Double I_TWOBODYOVERLAP_K5y2z_D2y_aa = I_TWOBODYOVERLAP_L6y2z_Py_aa+ABY*I_TWOBODYOVERLAP_K5y2z_Py_aa;
  Double I_TWOBODYOVERLAP_K4y3z_D2y_aa = I_TWOBODYOVERLAP_L5y3z_Py_aa+ABY*I_TWOBODYOVERLAP_K4y3z_Py_aa;
  Double I_TWOBODYOVERLAP_K3y4z_D2y_aa = I_TWOBODYOVERLAP_L4y4z_Py_aa+ABY*I_TWOBODYOVERLAP_K3y4z_Py_aa;
  Double I_TWOBODYOVERLAP_K2y5z_D2y_aa = I_TWOBODYOVERLAP_L3y5z_Py_aa+ABY*I_TWOBODYOVERLAP_K2y5z_Py_aa;
  Double I_TWOBODYOVERLAP_Ky6z_D2y_aa = I_TWOBODYOVERLAP_L2y6z_Py_aa+ABY*I_TWOBODYOVERLAP_Ky6z_Py_aa;
  Double I_TWOBODYOVERLAP_K7z_D2y_aa = I_TWOBODYOVERLAP_Ly7z_Py_aa+ABY*I_TWOBODYOVERLAP_K7z_Py_aa;
  Double I_TWOBODYOVERLAP_K7x_D2z_aa = I_TWOBODYOVERLAP_L7xz_Pz_aa+ABZ*I_TWOBODYOVERLAP_K7x_Pz_aa;
  Double I_TWOBODYOVERLAP_K6xy_D2z_aa = I_TWOBODYOVERLAP_L6xyz_Pz_aa+ABZ*I_TWOBODYOVERLAP_K6xy_Pz_aa;
  Double I_TWOBODYOVERLAP_K6xz_D2z_aa = I_TWOBODYOVERLAP_L6x2z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K6xz_Pz_aa;
  Double I_TWOBODYOVERLAP_K5x2y_D2z_aa = I_TWOBODYOVERLAP_L5x2yz_Pz_aa+ABZ*I_TWOBODYOVERLAP_K5x2y_Pz_aa;
  Double I_TWOBODYOVERLAP_K5xyz_D2z_aa = I_TWOBODYOVERLAP_L5xy2z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K5xyz_Pz_aa;
  Double I_TWOBODYOVERLAP_K5x2z_D2z_aa = I_TWOBODYOVERLAP_L5x3z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K5x2z_Pz_aa;
  Double I_TWOBODYOVERLAP_K4x3y_D2z_aa = I_TWOBODYOVERLAP_L4x3yz_Pz_aa+ABZ*I_TWOBODYOVERLAP_K4x3y_Pz_aa;
  Double I_TWOBODYOVERLAP_K4x2yz_D2z_aa = I_TWOBODYOVERLAP_L4x2y2z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K4x2yz_Pz_aa;
  Double I_TWOBODYOVERLAP_K4xy2z_D2z_aa = I_TWOBODYOVERLAP_L4xy3z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K4xy2z_Pz_aa;
  Double I_TWOBODYOVERLAP_K4x3z_D2z_aa = I_TWOBODYOVERLAP_L4x4z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K4x3z_Pz_aa;
  Double I_TWOBODYOVERLAP_K3x4y_D2z_aa = I_TWOBODYOVERLAP_L3x4yz_Pz_aa+ABZ*I_TWOBODYOVERLAP_K3x4y_Pz_aa;
  Double I_TWOBODYOVERLAP_K3x3yz_D2z_aa = I_TWOBODYOVERLAP_L3x3y2z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K3x3yz_Pz_aa;
  Double I_TWOBODYOVERLAP_K3x2y2z_D2z_aa = I_TWOBODYOVERLAP_L3x2y3z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K3x2y2z_Pz_aa;
  Double I_TWOBODYOVERLAP_K3xy3z_D2z_aa = I_TWOBODYOVERLAP_L3xy4z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K3xy3z_Pz_aa;
  Double I_TWOBODYOVERLAP_K3x4z_D2z_aa = I_TWOBODYOVERLAP_L3x5z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K3x4z_Pz_aa;
  Double I_TWOBODYOVERLAP_K2x5y_D2z_aa = I_TWOBODYOVERLAP_L2x5yz_Pz_aa+ABZ*I_TWOBODYOVERLAP_K2x5y_Pz_aa;
  Double I_TWOBODYOVERLAP_K2x4yz_D2z_aa = I_TWOBODYOVERLAP_L2x4y2z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K2x4yz_Pz_aa;
  Double I_TWOBODYOVERLAP_K2x3y2z_D2z_aa = I_TWOBODYOVERLAP_L2x3y3z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K2x3y2z_Pz_aa;
  Double I_TWOBODYOVERLAP_K2x2y3z_D2z_aa = I_TWOBODYOVERLAP_L2x2y4z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K2x2y3z_Pz_aa;
  Double I_TWOBODYOVERLAP_K2xy4z_D2z_aa = I_TWOBODYOVERLAP_L2xy5z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K2xy4z_Pz_aa;
  Double I_TWOBODYOVERLAP_K2x5z_D2z_aa = I_TWOBODYOVERLAP_L2x6z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K2x5z_Pz_aa;
  Double I_TWOBODYOVERLAP_Kx6y_D2z_aa = I_TWOBODYOVERLAP_Lx6yz_Pz_aa+ABZ*I_TWOBODYOVERLAP_Kx6y_Pz_aa;
  Double I_TWOBODYOVERLAP_Kx5yz_D2z_aa = I_TWOBODYOVERLAP_Lx5y2z_Pz_aa+ABZ*I_TWOBODYOVERLAP_Kx5yz_Pz_aa;
  Double I_TWOBODYOVERLAP_Kx4y2z_D2z_aa = I_TWOBODYOVERLAP_Lx4y3z_Pz_aa+ABZ*I_TWOBODYOVERLAP_Kx4y2z_Pz_aa;
  Double I_TWOBODYOVERLAP_Kx3y3z_D2z_aa = I_TWOBODYOVERLAP_Lx3y4z_Pz_aa+ABZ*I_TWOBODYOVERLAP_Kx3y3z_Pz_aa;
  Double I_TWOBODYOVERLAP_Kx2y4z_D2z_aa = I_TWOBODYOVERLAP_Lx2y5z_Pz_aa+ABZ*I_TWOBODYOVERLAP_Kx2y4z_Pz_aa;
  Double I_TWOBODYOVERLAP_Kxy5z_D2z_aa = I_TWOBODYOVERLAP_Lxy6z_Pz_aa+ABZ*I_TWOBODYOVERLAP_Kxy5z_Pz_aa;
  Double I_TWOBODYOVERLAP_Kx6z_D2z_aa = I_TWOBODYOVERLAP_Lx7z_Pz_aa+ABZ*I_TWOBODYOVERLAP_Kx6z_Pz_aa;
  Double I_TWOBODYOVERLAP_K7y_D2z_aa = I_TWOBODYOVERLAP_L7yz_Pz_aa+ABZ*I_TWOBODYOVERLAP_K7y_Pz_aa;
  Double I_TWOBODYOVERLAP_K6yz_D2z_aa = I_TWOBODYOVERLAP_L6y2z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K6yz_Pz_aa;
  Double I_TWOBODYOVERLAP_K5y2z_D2z_aa = I_TWOBODYOVERLAP_L5y3z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K5y2z_Pz_aa;
  Double I_TWOBODYOVERLAP_K4y3z_D2z_aa = I_TWOBODYOVERLAP_L4y4z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K4y3z_Pz_aa;
  Double I_TWOBODYOVERLAP_K3y4z_D2z_aa = I_TWOBODYOVERLAP_L3y5z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K3y4z_Pz_aa;
  Double I_TWOBODYOVERLAP_K2y5z_D2z_aa = I_TWOBODYOVERLAP_L2y6z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K2y5z_Pz_aa;
  Double I_TWOBODYOVERLAP_Ky6z_D2z_aa = I_TWOBODYOVERLAP_Ly7z_Pz_aa+ABZ*I_TWOBODYOVERLAP_Ky6z_Pz_aa;
  Double I_TWOBODYOVERLAP_K7z_D2z_aa = I_TWOBODYOVERLAP_L8z_Pz_aa+ABZ*I_TWOBODYOVERLAP_K7z_Pz_aa;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_I_F_aa
   * expanding position: BRA2
   * code section is: HRR
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_K_D_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_D_aa
   ************************************************************/
  Double I_TWOBODYOVERLAP_I6x_F3x_aa = I_TWOBODYOVERLAP_K7x_D2x_aa+ABX*I_TWOBODYOVERLAP_I6x_D2x_aa;
  Double I_TWOBODYOVERLAP_I5xy_F3x_aa = I_TWOBODYOVERLAP_K6xy_D2x_aa+ABX*I_TWOBODYOVERLAP_I5xy_D2x_aa;
  Double I_TWOBODYOVERLAP_I5xz_F3x_aa = I_TWOBODYOVERLAP_K6xz_D2x_aa+ABX*I_TWOBODYOVERLAP_I5xz_D2x_aa;
  Double I_TWOBODYOVERLAP_I4x2y_F3x_aa = I_TWOBODYOVERLAP_K5x2y_D2x_aa+ABX*I_TWOBODYOVERLAP_I4x2y_D2x_aa;
  Double I_TWOBODYOVERLAP_I4xyz_F3x_aa = I_TWOBODYOVERLAP_K5xyz_D2x_aa+ABX*I_TWOBODYOVERLAP_I4xyz_D2x_aa;
  Double I_TWOBODYOVERLAP_I4x2z_F3x_aa = I_TWOBODYOVERLAP_K5x2z_D2x_aa+ABX*I_TWOBODYOVERLAP_I4x2z_D2x_aa;
  Double I_TWOBODYOVERLAP_I3x3y_F3x_aa = I_TWOBODYOVERLAP_K4x3y_D2x_aa+ABX*I_TWOBODYOVERLAP_I3x3y_D2x_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_F3x_aa = I_TWOBODYOVERLAP_K4x2yz_D2x_aa+ABX*I_TWOBODYOVERLAP_I3x2yz_D2x_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_F3x_aa = I_TWOBODYOVERLAP_K4xy2z_D2x_aa+ABX*I_TWOBODYOVERLAP_I3xy2z_D2x_aa;
  Double I_TWOBODYOVERLAP_I3x3z_F3x_aa = I_TWOBODYOVERLAP_K4x3z_D2x_aa+ABX*I_TWOBODYOVERLAP_I3x3z_D2x_aa;
  Double I_TWOBODYOVERLAP_I2x4y_F3x_aa = I_TWOBODYOVERLAP_K3x4y_D2x_aa+ABX*I_TWOBODYOVERLAP_I2x4y_D2x_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_F3x_aa = I_TWOBODYOVERLAP_K3x3yz_D2x_aa+ABX*I_TWOBODYOVERLAP_I2x3yz_D2x_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_F3x_aa = I_TWOBODYOVERLAP_K3x2y2z_D2x_aa+ABX*I_TWOBODYOVERLAP_I2x2y2z_D2x_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_F3x_aa = I_TWOBODYOVERLAP_K3xy3z_D2x_aa+ABX*I_TWOBODYOVERLAP_I2xy3z_D2x_aa;
  Double I_TWOBODYOVERLAP_I2x4z_F3x_aa = I_TWOBODYOVERLAP_K3x4z_D2x_aa+ABX*I_TWOBODYOVERLAP_I2x4z_D2x_aa;
  Double I_TWOBODYOVERLAP_Ix5y_F3x_aa = I_TWOBODYOVERLAP_K2x5y_D2x_aa+ABX*I_TWOBODYOVERLAP_Ix5y_D2x_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_F3x_aa = I_TWOBODYOVERLAP_K2x4yz_D2x_aa+ABX*I_TWOBODYOVERLAP_Ix4yz_D2x_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_F3x_aa = I_TWOBODYOVERLAP_K2x3y2z_D2x_aa+ABX*I_TWOBODYOVERLAP_Ix3y2z_D2x_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_F3x_aa = I_TWOBODYOVERLAP_K2x2y3z_D2x_aa+ABX*I_TWOBODYOVERLAP_Ix2y3z_D2x_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_F3x_aa = I_TWOBODYOVERLAP_K2xy4z_D2x_aa+ABX*I_TWOBODYOVERLAP_Ixy4z_D2x_aa;
  Double I_TWOBODYOVERLAP_Ix5z_F3x_aa = I_TWOBODYOVERLAP_K2x5z_D2x_aa+ABX*I_TWOBODYOVERLAP_Ix5z_D2x_aa;
  Double I_TWOBODYOVERLAP_I6y_F3x_aa = I_TWOBODYOVERLAP_Kx6y_D2x_aa+ABX*I_TWOBODYOVERLAP_I6y_D2x_aa;
  Double I_TWOBODYOVERLAP_I5yz_F3x_aa = I_TWOBODYOVERLAP_Kx5yz_D2x_aa+ABX*I_TWOBODYOVERLAP_I5yz_D2x_aa;
  Double I_TWOBODYOVERLAP_I4y2z_F3x_aa = I_TWOBODYOVERLAP_Kx4y2z_D2x_aa+ABX*I_TWOBODYOVERLAP_I4y2z_D2x_aa;
  Double I_TWOBODYOVERLAP_I3y3z_F3x_aa = I_TWOBODYOVERLAP_Kx3y3z_D2x_aa+ABX*I_TWOBODYOVERLAP_I3y3z_D2x_aa;
  Double I_TWOBODYOVERLAP_I2y4z_F3x_aa = I_TWOBODYOVERLAP_Kx2y4z_D2x_aa+ABX*I_TWOBODYOVERLAP_I2y4z_D2x_aa;
  Double I_TWOBODYOVERLAP_Iy5z_F3x_aa = I_TWOBODYOVERLAP_Kxy5z_D2x_aa+ABX*I_TWOBODYOVERLAP_Iy5z_D2x_aa;
  Double I_TWOBODYOVERLAP_I6z_F3x_aa = I_TWOBODYOVERLAP_Kx6z_D2x_aa+ABX*I_TWOBODYOVERLAP_I6z_D2x_aa;
  Double I_TWOBODYOVERLAP_I6x_F2xy_aa = I_TWOBODYOVERLAP_K6xy_D2x_aa+ABY*I_TWOBODYOVERLAP_I6x_D2x_aa;
  Double I_TWOBODYOVERLAP_I5xy_F2xy_aa = I_TWOBODYOVERLAP_K5x2y_D2x_aa+ABY*I_TWOBODYOVERLAP_I5xy_D2x_aa;
  Double I_TWOBODYOVERLAP_I5xz_F2xy_aa = I_TWOBODYOVERLAP_K5xyz_D2x_aa+ABY*I_TWOBODYOVERLAP_I5xz_D2x_aa;
  Double I_TWOBODYOVERLAP_I4x2y_F2xy_aa = I_TWOBODYOVERLAP_K4x3y_D2x_aa+ABY*I_TWOBODYOVERLAP_I4x2y_D2x_aa;
  Double I_TWOBODYOVERLAP_I4xyz_F2xy_aa = I_TWOBODYOVERLAP_K4x2yz_D2x_aa+ABY*I_TWOBODYOVERLAP_I4xyz_D2x_aa;
  Double I_TWOBODYOVERLAP_I4x2z_F2xy_aa = I_TWOBODYOVERLAP_K4xy2z_D2x_aa+ABY*I_TWOBODYOVERLAP_I4x2z_D2x_aa;
  Double I_TWOBODYOVERLAP_I3x3y_F2xy_aa = I_TWOBODYOVERLAP_K3x4y_D2x_aa+ABY*I_TWOBODYOVERLAP_I3x3y_D2x_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_F2xy_aa = I_TWOBODYOVERLAP_K3x3yz_D2x_aa+ABY*I_TWOBODYOVERLAP_I3x2yz_D2x_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_F2xy_aa = I_TWOBODYOVERLAP_K3x2y2z_D2x_aa+ABY*I_TWOBODYOVERLAP_I3xy2z_D2x_aa;
  Double I_TWOBODYOVERLAP_I3x3z_F2xy_aa = I_TWOBODYOVERLAP_K3xy3z_D2x_aa+ABY*I_TWOBODYOVERLAP_I3x3z_D2x_aa;
  Double I_TWOBODYOVERLAP_I2x4y_F2xy_aa = I_TWOBODYOVERLAP_K2x5y_D2x_aa+ABY*I_TWOBODYOVERLAP_I2x4y_D2x_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_F2xy_aa = I_TWOBODYOVERLAP_K2x4yz_D2x_aa+ABY*I_TWOBODYOVERLAP_I2x3yz_D2x_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_F2xy_aa = I_TWOBODYOVERLAP_K2x3y2z_D2x_aa+ABY*I_TWOBODYOVERLAP_I2x2y2z_D2x_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_F2xy_aa = I_TWOBODYOVERLAP_K2x2y3z_D2x_aa+ABY*I_TWOBODYOVERLAP_I2xy3z_D2x_aa;
  Double I_TWOBODYOVERLAP_I2x4z_F2xy_aa = I_TWOBODYOVERLAP_K2xy4z_D2x_aa+ABY*I_TWOBODYOVERLAP_I2x4z_D2x_aa;
  Double I_TWOBODYOVERLAP_Ix5y_F2xy_aa = I_TWOBODYOVERLAP_Kx6y_D2x_aa+ABY*I_TWOBODYOVERLAP_Ix5y_D2x_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_F2xy_aa = I_TWOBODYOVERLAP_Kx5yz_D2x_aa+ABY*I_TWOBODYOVERLAP_Ix4yz_D2x_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_F2xy_aa = I_TWOBODYOVERLAP_Kx4y2z_D2x_aa+ABY*I_TWOBODYOVERLAP_Ix3y2z_D2x_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_F2xy_aa = I_TWOBODYOVERLAP_Kx3y3z_D2x_aa+ABY*I_TWOBODYOVERLAP_Ix2y3z_D2x_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_F2xy_aa = I_TWOBODYOVERLAP_Kx2y4z_D2x_aa+ABY*I_TWOBODYOVERLAP_Ixy4z_D2x_aa;
  Double I_TWOBODYOVERLAP_Ix5z_F2xy_aa = I_TWOBODYOVERLAP_Kxy5z_D2x_aa+ABY*I_TWOBODYOVERLAP_Ix5z_D2x_aa;
  Double I_TWOBODYOVERLAP_I6y_F2xy_aa = I_TWOBODYOVERLAP_K7y_D2x_aa+ABY*I_TWOBODYOVERLAP_I6y_D2x_aa;
  Double I_TWOBODYOVERLAP_I5yz_F2xy_aa = I_TWOBODYOVERLAP_K6yz_D2x_aa+ABY*I_TWOBODYOVERLAP_I5yz_D2x_aa;
  Double I_TWOBODYOVERLAP_I4y2z_F2xy_aa = I_TWOBODYOVERLAP_K5y2z_D2x_aa+ABY*I_TWOBODYOVERLAP_I4y2z_D2x_aa;
  Double I_TWOBODYOVERLAP_I3y3z_F2xy_aa = I_TWOBODYOVERLAP_K4y3z_D2x_aa+ABY*I_TWOBODYOVERLAP_I3y3z_D2x_aa;
  Double I_TWOBODYOVERLAP_I2y4z_F2xy_aa = I_TWOBODYOVERLAP_K3y4z_D2x_aa+ABY*I_TWOBODYOVERLAP_I2y4z_D2x_aa;
  Double I_TWOBODYOVERLAP_Iy5z_F2xy_aa = I_TWOBODYOVERLAP_K2y5z_D2x_aa+ABY*I_TWOBODYOVERLAP_Iy5z_D2x_aa;
  Double I_TWOBODYOVERLAP_I6z_F2xy_aa = I_TWOBODYOVERLAP_Ky6z_D2x_aa+ABY*I_TWOBODYOVERLAP_I6z_D2x_aa;
  Double I_TWOBODYOVERLAP_I6x_F2xz_aa = I_TWOBODYOVERLAP_K6xz_D2x_aa+ABZ*I_TWOBODYOVERLAP_I6x_D2x_aa;
  Double I_TWOBODYOVERLAP_I5xy_F2xz_aa = I_TWOBODYOVERLAP_K5xyz_D2x_aa+ABZ*I_TWOBODYOVERLAP_I5xy_D2x_aa;
  Double I_TWOBODYOVERLAP_I5xz_F2xz_aa = I_TWOBODYOVERLAP_K5x2z_D2x_aa+ABZ*I_TWOBODYOVERLAP_I5xz_D2x_aa;
  Double I_TWOBODYOVERLAP_I4x2y_F2xz_aa = I_TWOBODYOVERLAP_K4x2yz_D2x_aa+ABZ*I_TWOBODYOVERLAP_I4x2y_D2x_aa;
  Double I_TWOBODYOVERLAP_I4xyz_F2xz_aa = I_TWOBODYOVERLAP_K4xy2z_D2x_aa+ABZ*I_TWOBODYOVERLAP_I4xyz_D2x_aa;
  Double I_TWOBODYOVERLAP_I4x2z_F2xz_aa = I_TWOBODYOVERLAP_K4x3z_D2x_aa+ABZ*I_TWOBODYOVERLAP_I4x2z_D2x_aa;
  Double I_TWOBODYOVERLAP_I3x3y_F2xz_aa = I_TWOBODYOVERLAP_K3x3yz_D2x_aa+ABZ*I_TWOBODYOVERLAP_I3x3y_D2x_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_F2xz_aa = I_TWOBODYOVERLAP_K3x2y2z_D2x_aa+ABZ*I_TWOBODYOVERLAP_I3x2yz_D2x_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_F2xz_aa = I_TWOBODYOVERLAP_K3xy3z_D2x_aa+ABZ*I_TWOBODYOVERLAP_I3xy2z_D2x_aa;
  Double I_TWOBODYOVERLAP_I3x3z_F2xz_aa = I_TWOBODYOVERLAP_K3x4z_D2x_aa+ABZ*I_TWOBODYOVERLAP_I3x3z_D2x_aa;
  Double I_TWOBODYOVERLAP_I2x4y_F2xz_aa = I_TWOBODYOVERLAP_K2x4yz_D2x_aa+ABZ*I_TWOBODYOVERLAP_I2x4y_D2x_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_F2xz_aa = I_TWOBODYOVERLAP_K2x3y2z_D2x_aa+ABZ*I_TWOBODYOVERLAP_I2x3yz_D2x_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_F2xz_aa = I_TWOBODYOVERLAP_K2x2y3z_D2x_aa+ABZ*I_TWOBODYOVERLAP_I2x2y2z_D2x_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_F2xz_aa = I_TWOBODYOVERLAP_K2xy4z_D2x_aa+ABZ*I_TWOBODYOVERLAP_I2xy3z_D2x_aa;
  Double I_TWOBODYOVERLAP_I2x4z_F2xz_aa = I_TWOBODYOVERLAP_K2x5z_D2x_aa+ABZ*I_TWOBODYOVERLAP_I2x4z_D2x_aa;
  Double I_TWOBODYOVERLAP_Ix5y_F2xz_aa = I_TWOBODYOVERLAP_Kx5yz_D2x_aa+ABZ*I_TWOBODYOVERLAP_Ix5y_D2x_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_F2xz_aa = I_TWOBODYOVERLAP_Kx4y2z_D2x_aa+ABZ*I_TWOBODYOVERLAP_Ix4yz_D2x_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_F2xz_aa = I_TWOBODYOVERLAP_Kx3y3z_D2x_aa+ABZ*I_TWOBODYOVERLAP_Ix3y2z_D2x_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_F2xz_aa = I_TWOBODYOVERLAP_Kx2y4z_D2x_aa+ABZ*I_TWOBODYOVERLAP_Ix2y3z_D2x_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_F2xz_aa = I_TWOBODYOVERLAP_Kxy5z_D2x_aa+ABZ*I_TWOBODYOVERLAP_Ixy4z_D2x_aa;
  Double I_TWOBODYOVERLAP_Ix5z_F2xz_aa = I_TWOBODYOVERLAP_Kx6z_D2x_aa+ABZ*I_TWOBODYOVERLAP_Ix5z_D2x_aa;
  Double I_TWOBODYOVERLAP_I6y_F2xz_aa = I_TWOBODYOVERLAP_K6yz_D2x_aa+ABZ*I_TWOBODYOVERLAP_I6y_D2x_aa;
  Double I_TWOBODYOVERLAP_I5yz_F2xz_aa = I_TWOBODYOVERLAP_K5y2z_D2x_aa+ABZ*I_TWOBODYOVERLAP_I5yz_D2x_aa;
  Double I_TWOBODYOVERLAP_I4y2z_F2xz_aa = I_TWOBODYOVERLAP_K4y3z_D2x_aa+ABZ*I_TWOBODYOVERLAP_I4y2z_D2x_aa;
  Double I_TWOBODYOVERLAP_I3y3z_F2xz_aa = I_TWOBODYOVERLAP_K3y4z_D2x_aa+ABZ*I_TWOBODYOVERLAP_I3y3z_D2x_aa;
  Double I_TWOBODYOVERLAP_I2y4z_F2xz_aa = I_TWOBODYOVERLAP_K2y5z_D2x_aa+ABZ*I_TWOBODYOVERLAP_I2y4z_D2x_aa;
  Double I_TWOBODYOVERLAP_Iy5z_F2xz_aa = I_TWOBODYOVERLAP_Ky6z_D2x_aa+ABZ*I_TWOBODYOVERLAP_Iy5z_D2x_aa;
  Double I_TWOBODYOVERLAP_I6z_F2xz_aa = I_TWOBODYOVERLAP_K7z_D2x_aa+ABZ*I_TWOBODYOVERLAP_I6z_D2x_aa;
  Double I_TWOBODYOVERLAP_I6x_Fx2y_aa = I_TWOBODYOVERLAP_K7x_D2y_aa+ABX*I_TWOBODYOVERLAP_I6x_D2y_aa;
  Double I_TWOBODYOVERLAP_I5xy_Fx2y_aa = I_TWOBODYOVERLAP_K6xy_D2y_aa+ABX*I_TWOBODYOVERLAP_I5xy_D2y_aa;
  Double I_TWOBODYOVERLAP_I5xz_Fx2y_aa = I_TWOBODYOVERLAP_K6xz_D2y_aa+ABX*I_TWOBODYOVERLAP_I5xz_D2y_aa;
  Double I_TWOBODYOVERLAP_I4x2y_Fx2y_aa = I_TWOBODYOVERLAP_K5x2y_D2y_aa+ABX*I_TWOBODYOVERLAP_I4x2y_D2y_aa;
  Double I_TWOBODYOVERLAP_I4xyz_Fx2y_aa = I_TWOBODYOVERLAP_K5xyz_D2y_aa+ABX*I_TWOBODYOVERLAP_I4xyz_D2y_aa;
  Double I_TWOBODYOVERLAP_I4x2z_Fx2y_aa = I_TWOBODYOVERLAP_K5x2z_D2y_aa+ABX*I_TWOBODYOVERLAP_I4x2z_D2y_aa;
  Double I_TWOBODYOVERLAP_I3x3y_Fx2y_aa = I_TWOBODYOVERLAP_K4x3y_D2y_aa+ABX*I_TWOBODYOVERLAP_I3x3y_D2y_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_Fx2y_aa = I_TWOBODYOVERLAP_K4x2yz_D2y_aa+ABX*I_TWOBODYOVERLAP_I3x2yz_D2y_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_Fx2y_aa = I_TWOBODYOVERLAP_K4xy2z_D2y_aa+ABX*I_TWOBODYOVERLAP_I3xy2z_D2y_aa;
  Double I_TWOBODYOVERLAP_I3x3z_Fx2y_aa = I_TWOBODYOVERLAP_K4x3z_D2y_aa+ABX*I_TWOBODYOVERLAP_I3x3z_D2y_aa;
  Double I_TWOBODYOVERLAP_I2x4y_Fx2y_aa = I_TWOBODYOVERLAP_K3x4y_D2y_aa+ABX*I_TWOBODYOVERLAP_I2x4y_D2y_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_Fx2y_aa = I_TWOBODYOVERLAP_K3x3yz_D2y_aa+ABX*I_TWOBODYOVERLAP_I2x3yz_D2y_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_Fx2y_aa = I_TWOBODYOVERLAP_K3x2y2z_D2y_aa+ABX*I_TWOBODYOVERLAP_I2x2y2z_D2y_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_Fx2y_aa = I_TWOBODYOVERLAP_K3xy3z_D2y_aa+ABX*I_TWOBODYOVERLAP_I2xy3z_D2y_aa;
  Double I_TWOBODYOVERLAP_I2x4z_Fx2y_aa = I_TWOBODYOVERLAP_K3x4z_D2y_aa+ABX*I_TWOBODYOVERLAP_I2x4z_D2y_aa;
  Double I_TWOBODYOVERLAP_Ix5y_Fx2y_aa = I_TWOBODYOVERLAP_K2x5y_D2y_aa+ABX*I_TWOBODYOVERLAP_Ix5y_D2y_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_Fx2y_aa = I_TWOBODYOVERLAP_K2x4yz_D2y_aa+ABX*I_TWOBODYOVERLAP_Ix4yz_D2y_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_Fx2y_aa = I_TWOBODYOVERLAP_K2x3y2z_D2y_aa+ABX*I_TWOBODYOVERLAP_Ix3y2z_D2y_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_Fx2y_aa = I_TWOBODYOVERLAP_K2x2y3z_D2y_aa+ABX*I_TWOBODYOVERLAP_Ix2y3z_D2y_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_Fx2y_aa = I_TWOBODYOVERLAP_K2xy4z_D2y_aa+ABX*I_TWOBODYOVERLAP_Ixy4z_D2y_aa;
  Double I_TWOBODYOVERLAP_Ix5z_Fx2y_aa = I_TWOBODYOVERLAP_K2x5z_D2y_aa+ABX*I_TWOBODYOVERLAP_Ix5z_D2y_aa;
  Double I_TWOBODYOVERLAP_I6y_Fx2y_aa = I_TWOBODYOVERLAP_Kx6y_D2y_aa+ABX*I_TWOBODYOVERLAP_I6y_D2y_aa;
  Double I_TWOBODYOVERLAP_I5yz_Fx2y_aa = I_TWOBODYOVERLAP_Kx5yz_D2y_aa+ABX*I_TWOBODYOVERLAP_I5yz_D2y_aa;
  Double I_TWOBODYOVERLAP_I4y2z_Fx2y_aa = I_TWOBODYOVERLAP_Kx4y2z_D2y_aa+ABX*I_TWOBODYOVERLAP_I4y2z_D2y_aa;
  Double I_TWOBODYOVERLAP_I3y3z_Fx2y_aa = I_TWOBODYOVERLAP_Kx3y3z_D2y_aa+ABX*I_TWOBODYOVERLAP_I3y3z_D2y_aa;
  Double I_TWOBODYOVERLAP_I2y4z_Fx2y_aa = I_TWOBODYOVERLAP_Kx2y4z_D2y_aa+ABX*I_TWOBODYOVERLAP_I2y4z_D2y_aa;
  Double I_TWOBODYOVERLAP_Iy5z_Fx2y_aa = I_TWOBODYOVERLAP_Kxy5z_D2y_aa+ABX*I_TWOBODYOVERLAP_Iy5z_D2y_aa;
  Double I_TWOBODYOVERLAP_I6z_Fx2y_aa = I_TWOBODYOVERLAP_Kx6z_D2y_aa+ABX*I_TWOBODYOVERLAP_I6z_D2y_aa;
  Double I_TWOBODYOVERLAP_I6x_Fxyz_aa = I_TWOBODYOVERLAP_K6xz_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I6x_Dxy_aa;
  Double I_TWOBODYOVERLAP_I5xy_Fxyz_aa = I_TWOBODYOVERLAP_K5xyz_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I5xy_Dxy_aa;
  Double I_TWOBODYOVERLAP_I5xz_Fxyz_aa = I_TWOBODYOVERLAP_K5x2z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I5xz_Dxy_aa;
  Double I_TWOBODYOVERLAP_I4x2y_Fxyz_aa = I_TWOBODYOVERLAP_K4x2yz_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I4x2y_Dxy_aa;
  Double I_TWOBODYOVERLAP_I4xyz_Fxyz_aa = I_TWOBODYOVERLAP_K4xy2z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I4xyz_Dxy_aa;
  Double I_TWOBODYOVERLAP_I4x2z_Fxyz_aa = I_TWOBODYOVERLAP_K4x3z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I4x2z_Dxy_aa;
  Double I_TWOBODYOVERLAP_I3x3y_Fxyz_aa = I_TWOBODYOVERLAP_K3x3yz_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I3x3y_Dxy_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_Fxyz_aa = I_TWOBODYOVERLAP_K3x2y2z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I3x2yz_Dxy_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_Fxyz_aa = I_TWOBODYOVERLAP_K3xy3z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I3xy2z_Dxy_aa;
  Double I_TWOBODYOVERLAP_I3x3z_Fxyz_aa = I_TWOBODYOVERLAP_K3x4z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I3x3z_Dxy_aa;
  Double I_TWOBODYOVERLAP_I2x4y_Fxyz_aa = I_TWOBODYOVERLAP_K2x4yz_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I2x4y_Dxy_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_Fxyz_aa = I_TWOBODYOVERLAP_K2x3y2z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I2x3yz_Dxy_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_Fxyz_aa = I_TWOBODYOVERLAP_K2x2y3z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I2x2y2z_Dxy_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_Fxyz_aa = I_TWOBODYOVERLAP_K2xy4z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I2xy3z_Dxy_aa;
  Double I_TWOBODYOVERLAP_I2x4z_Fxyz_aa = I_TWOBODYOVERLAP_K2x5z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I2x4z_Dxy_aa;
  Double I_TWOBODYOVERLAP_Ix5y_Fxyz_aa = I_TWOBODYOVERLAP_Kx5yz_Dxy_aa+ABZ*I_TWOBODYOVERLAP_Ix5y_Dxy_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_Fxyz_aa = I_TWOBODYOVERLAP_Kx4y2z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_Ix4yz_Dxy_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_Fxyz_aa = I_TWOBODYOVERLAP_Kx3y3z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_Ix3y2z_Dxy_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_Fxyz_aa = I_TWOBODYOVERLAP_Kx2y4z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_Ix2y3z_Dxy_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_Fxyz_aa = I_TWOBODYOVERLAP_Kxy5z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_Ixy4z_Dxy_aa;
  Double I_TWOBODYOVERLAP_Ix5z_Fxyz_aa = I_TWOBODYOVERLAP_Kx6z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_Ix5z_Dxy_aa;
  Double I_TWOBODYOVERLAP_I6y_Fxyz_aa = I_TWOBODYOVERLAP_K6yz_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I6y_Dxy_aa;
  Double I_TWOBODYOVERLAP_I5yz_Fxyz_aa = I_TWOBODYOVERLAP_K5y2z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I5yz_Dxy_aa;
  Double I_TWOBODYOVERLAP_I4y2z_Fxyz_aa = I_TWOBODYOVERLAP_K4y3z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I4y2z_Dxy_aa;
  Double I_TWOBODYOVERLAP_I3y3z_Fxyz_aa = I_TWOBODYOVERLAP_K3y4z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I3y3z_Dxy_aa;
  Double I_TWOBODYOVERLAP_I2y4z_Fxyz_aa = I_TWOBODYOVERLAP_K2y5z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I2y4z_Dxy_aa;
  Double I_TWOBODYOVERLAP_Iy5z_Fxyz_aa = I_TWOBODYOVERLAP_Ky6z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_Iy5z_Dxy_aa;
  Double I_TWOBODYOVERLAP_I6z_Fxyz_aa = I_TWOBODYOVERLAP_K7z_Dxy_aa+ABZ*I_TWOBODYOVERLAP_I6z_Dxy_aa;
  Double I_TWOBODYOVERLAP_I6x_Fx2z_aa = I_TWOBODYOVERLAP_K7x_D2z_aa+ABX*I_TWOBODYOVERLAP_I6x_D2z_aa;
  Double I_TWOBODYOVERLAP_I5xy_Fx2z_aa = I_TWOBODYOVERLAP_K6xy_D2z_aa+ABX*I_TWOBODYOVERLAP_I5xy_D2z_aa;
  Double I_TWOBODYOVERLAP_I5xz_Fx2z_aa = I_TWOBODYOVERLAP_K6xz_D2z_aa+ABX*I_TWOBODYOVERLAP_I5xz_D2z_aa;
  Double I_TWOBODYOVERLAP_I4x2y_Fx2z_aa = I_TWOBODYOVERLAP_K5x2y_D2z_aa+ABX*I_TWOBODYOVERLAP_I4x2y_D2z_aa;
  Double I_TWOBODYOVERLAP_I4xyz_Fx2z_aa = I_TWOBODYOVERLAP_K5xyz_D2z_aa+ABX*I_TWOBODYOVERLAP_I4xyz_D2z_aa;
  Double I_TWOBODYOVERLAP_I4x2z_Fx2z_aa = I_TWOBODYOVERLAP_K5x2z_D2z_aa+ABX*I_TWOBODYOVERLAP_I4x2z_D2z_aa;
  Double I_TWOBODYOVERLAP_I3x3y_Fx2z_aa = I_TWOBODYOVERLAP_K4x3y_D2z_aa+ABX*I_TWOBODYOVERLAP_I3x3y_D2z_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_Fx2z_aa = I_TWOBODYOVERLAP_K4x2yz_D2z_aa+ABX*I_TWOBODYOVERLAP_I3x2yz_D2z_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_Fx2z_aa = I_TWOBODYOVERLAP_K4xy2z_D2z_aa+ABX*I_TWOBODYOVERLAP_I3xy2z_D2z_aa;
  Double I_TWOBODYOVERLAP_I3x3z_Fx2z_aa = I_TWOBODYOVERLAP_K4x3z_D2z_aa+ABX*I_TWOBODYOVERLAP_I3x3z_D2z_aa;
  Double I_TWOBODYOVERLAP_I2x4y_Fx2z_aa = I_TWOBODYOVERLAP_K3x4y_D2z_aa+ABX*I_TWOBODYOVERLAP_I2x4y_D2z_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_Fx2z_aa = I_TWOBODYOVERLAP_K3x3yz_D2z_aa+ABX*I_TWOBODYOVERLAP_I2x3yz_D2z_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_Fx2z_aa = I_TWOBODYOVERLAP_K3x2y2z_D2z_aa+ABX*I_TWOBODYOVERLAP_I2x2y2z_D2z_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_Fx2z_aa = I_TWOBODYOVERLAP_K3xy3z_D2z_aa+ABX*I_TWOBODYOVERLAP_I2xy3z_D2z_aa;
  Double I_TWOBODYOVERLAP_I2x4z_Fx2z_aa = I_TWOBODYOVERLAP_K3x4z_D2z_aa+ABX*I_TWOBODYOVERLAP_I2x4z_D2z_aa;
  Double I_TWOBODYOVERLAP_Ix5y_Fx2z_aa = I_TWOBODYOVERLAP_K2x5y_D2z_aa+ABX*I_TWOBODYOVERLAP_Ix5y_D2z_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_Fx2z_aa = I_TWOBODYOVERLAP_K2x4yz_D2z_aa+ABX*I_TWOBODYOVERLAP_Ix4yz_D2z_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_Fx2z_aa = I_TWOBODYOVERLAP_K2x3y2z_D2z_aa+ABX*I_TWOBODYOVERLAP_Ix3y2z_D2z_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_Fx2z_aa = I_TWOBODYOVERLAP_K2x2y3z_D2z_aa+ABX*I_TWOBODYOVERLAP_Ix2y3z_D2z_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_Fx2z_aa = I_TWOBODYOVERLAP_K2xy4z_D2z_aa+ABX*I_TWOBODYOVERLAP_Ixy4z_D2z_aa;
  Double I_TWOBODYOVERLAP_Ix5z_Fx2z_aa = I_TWOBODYOVERLAP_K2x5z_D2z_aa+ABX*I_TWOBODYOVERLAP_Ix5z_D2z_aa;
  Double I_TWOBODYOVERLAP_I6y_Fx2z_aa = I_TWOBODYOVERLAP_Kx6y_D2z_aa+ABX*I_TWOBODYOVERLAP_I6y_D2z_aa;
  Double I_TWOBODYOVERLAP_I5yz_Fx2z_aa = I_TWOBODYOVERLAP_Kx5yz_D2z_aa+ABX*I_TWOBODYOVERLAP_I5yz_D2z_aa;
  Double I_TWOBODYOVERLAP_I4y2z_Fx2z_aa = I_TWOBODYOVERLAP_Kx4y2z_D2z_aa+ABX*I_TWOBODYOVERLAP_I4y2z_D2z_aa;
  Double I_TWOBODYOVERLAP_I3y3z_Fx2z_aa = I_TWOBODYOVERLAP_Kx3y3z_D2z_aa+ABX*I_TWOBODYOVERLAP_I3y3z_D2z_aa;
  Double I_TWOBODYOVERLAP_I2y4z_Fx2z_aa = I_TWOBODYOVERLAP_Kx2y4z_D2z_aa+ABX*I_TWOBODYOVERLAP_I2y4z_D2z_aa;
  Double I_TWOBODYOVERLAP_Iy5z_Fx2z_aa = I_TWOBODYOVERLAP_Kxy5z_D2z_aa+ABX*I_TWOBODYOVERLAP_Iy5z_D2z_aa;
  Double I_TWOBODYOVERLAP_I6z_Fx2z_aa = I_TWOBODYOVERLAP_Kx6z_D2z_aa+ABX*I_TWOBODYOVERLAP_I6z_D2z_aa;
  Double I_TWOBODYOVERLAP_I6x_F3y_aa = I_TWOBODYOVERLAP_K6xy_D2y_aa+ABY*I_TWOBODYOVERLAP_I6x_D2y_aa;
  Double I_TWOBODYOVERLAP_I5xy_F3y_aa = I_TWOBODYOVERLAP_K5x2y_D2y_aa+ABY*I_TWOBODYOVERLAP_I5xy_D2y_aa;
  Double I_TWOBODYOVERLAP_I5xz_F3y_aa = I_TWOBODYOVERLAP_K5xyz_D2y_aa+ABY*I_TWOBODYOVERLAP_I5xz_D2y_aa;
  Double I_TWOBODYOVERLAP_I4x2y_F3y_aa = I_TWOBODYOVERLAP_K4x3y_D2y_aa+ABY*I_TWOBODYOVERLAP_I4x2y_D2y_aa;
  Double I_TWOBODYOVERLAP_I4xyz_F3y_aa = I_TWOBODYOVERLAP_K4x2yz_D2y_aa+ABY*I_TWOBODYOVERLAP_I4xyz_D2y_aa;
  Double I_TWOBODYOVERLAP_I4x2z_F3y_aa = I_TWOBODYOVERLAP_K4xy2z_D2y_aa+ABY*I_TWOBODYOVERLAP_I4x2z_D2y_aa;
  Double I_TWOBODYOVERLAP_I3x3y_F3y_aa = I_TWOBODYOVERLAP_K3x4y_D2y_aa+ABY*I_TWOBODYOVERLAP_I3x3y_D2y_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_F3y_aa = I_TWOBODYOVERLAP_K3x3yz_D2y_aa+ABY*I_TWOBODYOVERLAP_I3x2yz_D2y_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_F3y_aa = I_TWOBODYOVERLAP_K3x2y2z_D2y_aa+ABY*I_TWOBODYOVERLAP_I3xy2z_D2y_aa;
  Double I_TWOBODYOVERLAP_I3x3z_F3y_aa = I_TWOBODYOVERLAP_K3xy3z_D2y_aa+ABY*I_TWOBODYOVERLAP_I3x3z_D2y_aa;
  Double I_TWOBODYOVERLAP_I2x4y_F3y_aa = I_TWOBODYOVERLAP_K2x5y_D2y_aa+ABY*I_TWOBODYOVERLAP_I2x4y_D2y_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_F3y_aa = I_TWOBODYOVERLAP_K2x4yz_D2y_aa+ABY*I_TWOBODYOVERLAP_I2x3yz_D2y_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_F3y_aa = I_TWOBODYOVERLAP_K2x3y2z_D2y_aa+ABY*I_TWOBODYOVERLAP_I2x2y2z_D2y_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_F3y_aa = I_TWOBODYOVERLAP_K2x2y3z_D2y_aa+ABY*I_TWOBODYOVERLAP_I2xy3z_D2y_aa;
  Double I_TWOBODYOVERLAP_I2x4z_F3y_aa = I_TWOBODYOVERLAP_K2xy4z_D2y_aa+ABY*I_TWOBODYOVERLAP_I2x4z_D2y_aa;
  Double I_TWOBODYOVERLAP_Ix5y_F3y_aa = I_TWOBODYOVERLAP_Kx6y_D2y_aa+ABY*I_TWOBODYOVERLAP_Ix5y_D2y_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_F3y_aa = I_TWOBODYOVERLAP_Kx5yz_D2y_aa+ABY*I_TWOBODYOVERLAP_Ix4yz_D2y_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_F3y_aa = I_TWOBODYOVERLAP_Kx4y2z_D2y_aa+ABY*I_TWOBODYOVERLAP_Ix3y2z_D2y_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_F3y_aa = I_TWOBODYOVERLAP_Kx3y3z_D2y_aa+ABY*I_TWOBODYOVERLAP_Ix2y3z_D2y_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_F3y_aa = I_TWOBODYOVERLAP_Kx2y4z_D2y_aa+ABY*I_TWOBODYOVERLAP_Ixy4z_D2y_aa;
  Double I_TWOBODYOVERLAP_Ix5z_F3y_aa = I_TWOBODYOVERLAP_Kxy5z_D2y_aa+ABY*I_TWOBODYOVERLAP_Ix5z_D2y_aa;
  Double I_TWOBODYOVERLAP_I6y_F3y_aa = I_TWOBODYOVERLAP_K7y_D2y_aa+ABY*I_TWOBODYOVERLAP_I6y_D2y_aa;
  Double I_TWOBODYOVERLAP_I5yz_F3y_aa = I_TWOBODYOVERLAP_K6yz_D2y_aa+ABY*I_TWOBODYOVERLAP_I5yz_D2y_aa;
  Double I_TWOBODYOVERLAP_I4y2z_F3y_aa = I_TWOBODYOVERLAP_K5y2z_D2y_aa+ABY*I_TWOBODYOVERLAP_I4y2z_D2y_aa;
  Double I_TWOBODYOVERLAP_I3y3z_F3y_aa = I_TWOBODYOVERLAP_K4y3z_D2y_aa+ABY*I_TWOBODYOVERLAP_I3y3z_D2y_aa;
  Double I_TWOBODYOVERLAP_I2y4z_F3y_aa = I_TWOBODYOVERLAP_K3y4z_D2y_aa+ABY*I_TWOBODYOVERLAP_I2y4z_D2y_aa;
  Double I_TWOBODYOVERLAP_Iy5z_F3y_aa = I_TWOBODYOVERLAP_K2y5z_D2y_aa+ABY*I_TWOBODYOVERLAP_Iy5z_D2y_aa;
  Double I_TWOBODYOVERLAP_I6z_F3y_aa = I_TWOBODYOVERLAP_Ky6z_D2y_aa+ABY*I_TWOBODYOVERLAP_I6z_D2y_aa;
  Double I_TWOBODYOVERLAP_I6x_F2yz_aa = I_TWOBODYOVERLAP_K6xz_D2y_aa+ABZ*I_TWOBODYOVERLAP_I6x_D2y_aa;
  Double I_TWOBODYOVERLAP_I5xy_F2yz_aa = I_TWOBODYOVERLAP_K5xyz_D2y_aa+ABZ*I_TWOBODYOVERLAP_I5xy_D2y_aa;
  Double I_TWOBODYOVERLAP_I5xz_F2yz_aa = I_TWOBODYOVERLAP_K5x2z_D2y_aa+ABZ*I_TWOBODYOVERLAP_I5xz_D2y_aa;
  Double I_TWOBODYOVERLAP_I4x2y_F2yz_aa = I_TWOBODYOVERLAP_K4x2yz_D2y_aa+ABZ*I_TWOBODYOVERLAP_I4x2y_D2y_aa;
  Double I_TWOBODYOVERLAP_I4xyz_F2yz_aa = I_TWOBODYOVERLAP_K4xy2z_D2y_aa+ABZ*I_TWOBODYOVERLAP_I4xyz_D2y_aa;
  Double I_TWOBODYOVERLAP_I4x2z_F2yz_aa = I_TWOBODYOVERLAP_K4x3z_D2y_aa+ABZ*I_TWOBODYOVERLAP_I4x2z_D2y_aa;
  Double I_TWOBODYOVERLAP_I3x3y_F2yz_aa = I_TWOBODYOVERLAP_K3x3yz_D2y_aa+ABZ*I_TWOBODYOVERLAP_I3x3y_D2y_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_F2yz_aa = I_TWOBODYOVERLAP_K3x2y2z_D2y_aa+ABZ*I_TWOBODYOVERLAP_I3x2yz_D2y_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_F2yz_aa = I_TWOBODYOVERLAP_K3xy3z_D2y_aa+ABZ*I_TWOBODYOVERLAP_I3xy2z_D2y_aa;
  Double I_TWOBODYOVERLAP_I3x3z_F2yz_aa = I_TWOBODYOVERLAP_K3x4z_D2y_aa+ABZ*I_TWOBODYOVERLAP_I3x3z_D2y_aa;
  Double I_TWOBODYOVERLAP_I2x4y_F2yz_aa = I_TWOBODYOVERLAP_K2x4yz_D2y_aa+ABZ*I_TWOBODYOVERLAP_I2x4y_D2y_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_F2yz_aa = I_TWOBODYOVERLAP_K2x3y2z_D2y_aa+ABZ*I_TWOBODYOVERLAP_I2x3yz_D2y_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_F2yz_aa = I_TWOBODYOVERLAP_K2x2y3z_D2y_aa+ABZ*I_TWOBODYOVERLAP_I2x2y2z_D2y_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_F2yz_aa = I_TWOBODYOVERLAP_K2xy4z_D2y_aa+ABZ*I_TWOBODYOVERLAP_I2xy3z_D2y_aa;
  Double I_TWOBODYOVERLAP_I2x4z_F2yz_aa = I_TWOBODYOVERLAP_K2x5z_D2y_aa+ABZ*I_TWOBODYOVERLAP_I2x4z_D2y_aa;
  Double I_TWOBODYOVERLAP_Ix5y_F2yz_aa = I_TWOBODYOVERLAP_Kx5yz_D2y_aa+ABZ*I_TWOBODYOVERLAP_Ix5y_D2y_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_F2yz_aa = I_TWOBODYOVERLAP_Kx4y2z_D2y_aa+ABZ*I_TWOBODYOVERLAP_Ix4yz_D2y_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_F2yz_aa = I_TWOBODYOVERLAP_Kx3y3z_D2y_aa+ABZ*I_TWOBODYOVERLAP_Ix3y2z_D2y_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_F2yz_aa = I_TWOBODYOVERLAP_Kx2y4z_D2y_aa+ABZ*I_TWOBODYOVERLAP_Ix2y3z_D2y_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_F2yz_aa = I_TWOBODYOVERLAP_Kxy5z_D2y_aa+ABZ*I_TWOBODYOVERLAP_Ixy4z_D2y_aa;
  Double I_TWOBODYOVERLAP_Ix5z_F2yz_aa = I_TWOBODYOVERLAP_Kx6z_D2y_aa+ABZ*I_TWOBODYOVERLAP_Ix5z_D2y_aa;
  Double I_TWOBODYOVERLAP_I6y_F2yz_aa = I_TWOBODYOVERLAP_K6yz_D2y_aa+ABZ*I_TWOBODYOVERLAP_I6y_D2y_aa;
  Double I_TWOBODYOVERLAP_I5yz_F2yz_aa = I_TWOBODYOVERLAP_K5y2z_D2y_aa+ABZ*I_TWOBODYOVERLAP_I5yz_D2y_aa;
  Double I_TWOBODYOVERLAP_I4y2z_F2yz_aa = I_TWOBODYOVERLAP_K4y3z_D2y_aa+ABZ*I_TWOBODYOVERLAP_I4y2z_D2y_aa;
  Double I_TWOBODYOVERLAP_I3y3z_F2yz_aa = I_TWOBODYOVERLAP_K3y4z_D2y_aa+ABZ*I_TWOBODYOVERLAP_I3y3z_D2y_aa;
  Double I_TWOBODYOVERLAP_I2y4z_F2yz_aa = I_TWOBODYOVERLAP_K2y5z_D2y_aa+ABZ*I_TWOBODYOVERLAP_I2y4z_D2y_aa;
  Double I_TWOBODYOVERLAP_Iy5z_F2yz_aa = I_TWOBODYOVERLAP_Ky6z_D2y_aa+ABZ*I_TWOBODYOVERLAP_Iy5z_D2y_aa;
  Double I_TWOBODYOVERLAP_I6z_F2yz_aa = I_TWOBODYOVERLAP_K7z_D2y_aa+ABZ*I_TWOBODYOVERLAP_I6z_D2y_aa;
  Double I_TWOBODYOVERLAP_I6x_Fy2z_aa = I_TWOBODYOVERLAP_K6xy_D2z_aa+ABY*I_TWOBODYOVERLAP_I6x_D2z_aa;
  Double I_TWOBODYOVERLAP_I5xy_Fy2z_aa = I_TWOBODYOVERLAP_K5x2y_D2z_aa+ABY*I_TWOBODYOVERLAP_I5xy_D2z_aa;
  Double I_TWOBODYOVERLAP_I5xz_Fy2z_aa = I_TWOBODYOVERLAP_K5xyz_D2z_aa+ABY*I_TWOBODYOVERLAP_I5xz_D2z_aa;
  Double I_TWOBODYOVERLAP_I4x2y_Fy2z_aa = I_TWOBODYOVERLAP_K4x3y_D2z_aa+ABY*I_TWOBODYOVERLAP_I4x2y_D2z_aa;
  Double I_TWOBODYOVERLAP_I4xyz_Fy2z_aa = I_TWOBODYOVERLAP_K4x2yz_D2z_aa+ABY*I_TWOBODYOVERLAP_I4xyz_D2z_aa;
  Double I_TWOBODYOVERLAP_I4x2z_Fy2z_aa = I_TWOBODYOVERLAP_K4xy2z_D2z_aa+ABY*I_TWOBODYOVERLAP_I4x2z_D2z_aa;
  Double I_TWOBODYOVERLAP_I3x3y_Fy2z_aa = I_TWOBODYOVERLAP_K3x4y_D2z_aa+ABY*I_TWOBODYOVERLAP_I3x3y_D2z_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_Fy2z_aa = I_TWOBODYOVERLAP_K3x3yz_D2z_aa+ABY*I_TWOBODYOVERLAP_I3x2yz_D2z_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_Fy2z_aa = I_TWOBODYOVERLAP_K3x2y2z_D2z_aa+ABY*I_TWOBODYOVERLAP_I3xy2z_D2z_aa;
  Double I_TWOBODYOVERLAP_I3x3z_Fy2z_aa = I_TWOBODYOVERLAP_K3xy3z_D2z_aa+ABY*I_TWOBODYOVERLAP_I3x3z_D2z_aa;
  Double I_TWOBODYOVERLAP_I2x4y_Fy2z_aa = I_TWOBODYOVERLAP_K2x5y_D2z_aa+ABY*I_TWOBODYOVERLAP_I2x4y_D2z_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_Fy2z_aa = I_TWOBODYOVERLAP_K2x4yz_D2z_aa+ABY*I_TWOBODYOVERLAP_I2x3yz_D2z_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_Fy2z_aa = I_TWOBODYOVERLAP_K2x3y2z_D2z_aa+ABY*I_TWOBODYOVERLAP_I2x2y2z_D2z_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_Fy2z_aa = I_TWOBODYOVERLAP_K2x2y3z_D2z_aa+ABY*I_TWOBODYOVERLAP_I2xy3z_D2z_aa;
  Double I_TWOBODYOVERLAP_I2x4z_Fy2z_aa = I_TWOBODYOVERLAP_K2xy4z_D2z_aa+ABY*I_TWOBODYOVERLAP_I2x4z_D2z_aa;
  Double I_TWOBODYOVERLAP_Ix5y_Fy2z_aa = I_TWOBODYOVERLAP_Kx6y_D2z_aa+ABY*I_TWOBODYOVERLAP_Ix5y_D2z_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_Fy2z_aa = I_TWOBODYOVERLAP_Kx5yz_D2z_aa+ABY*I_TWOBODYOVERLAP_Ix4yz_D2z_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_Fy2z_aa = I_TWOBODYOVERLAP_Kx4y2z_D2z_aa+ABY*I_TWOBODYOVERLAP_Ix3y2z_D2z_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_Fy2z_aa = I_TWOBODYOVERLAP_Kx3y3z_D2z_aa+ABY*I_TWOBODYOVERLAP_Ix2y3z_D2z_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_Fy2z_aa = I_TWOBODYOVERLAP_Kx2y4z_D2z_aa+ABY*I_TWOBODYOVERLAP_Ixy4z_D2z_aa;
  Double I_TWOBODYOVERLAP_Ix5z_Fy2z_aa = I_TWOBODYOVERLAP_Kxy5z_D2z_aa+ABY*I_TWOBODYOVERLAP_Ix5z_D2z_aa;
  Double I_TWOBODYOVERLAP_I6y_Fy2z_aa = I_TWOBODYOVERLAP_K7y_D2z_aa+ABY*I_TWOBODYOVERLAP_I6y_D2z_aa;
  Double I_TWOBODYOVERLAP_I5yz_Fy2z_aa = I_TWOBODYOVERLAP_K6yz_D2z_aa+ABY*I_TWOBODYOVERLAP_I5yz_D2z_aa;
  Double I_TWOBODYOVERLAP_I4y2z_Fy2z_aa = I_TWOBODYOVERLAP_K5y2z_D2z_aa+ABY*I_TWOBODYOVERLAP_I4y2z_D2z_aa;
  Double I_TWOBODYOVERLAP_I3y3z_Fy2z_aa = I_TWOBODYOVERLAP_K4y3z_D2z_aa+ABY*I_TWOBODYOVERLAP_I3y3z_D2z_aa;
  Double I_TWOBODYOVERLAP_I2y4z_Fy2z_aa = I_TWOBODYOVERLAP_K3y4z_D2z_aa+ABY*I_TWOBODYOVERLAP_I2y4z_D2z_aa;
  Double I_TWOBODYOVERLAP_Iy5z_Fy2z_aa = I_TWOBODYOVERLAP_K2y5z_D2z_aa+ABY*I_TWOBODYOVERLAP_Iy5z_D2z_aa;
  Double I_TWOBODYOVERLAP_I6z_Fy2z_aa = I_TWOBODYOVERLAP_Ky6z_D2z_aa+ABY*I_TWOBODYOVERLAP_I6z_D2z_aa;
  Double I_TWOBODYOVERLAP_I6x_F3z_aa = I_TWOBODYOVERLAP_K6xz_D2z_aa+ABZ*I_TWOBODYOVERLAP_I6x_D2z_aa;
  Double I_TWOBODYOVERLAP_I5xy_F3z_aa = I_TWOBODYOVERLAP_K5xyz_D2z_aa+ABZ*I_TWOBODYOVERLAP_I5xy_D2z_aa;
  Double I_TWOBODYOVERLAP_I5xz_F3z_aa = I_TWOBODYOVERLAP_K5x2z_D2z_aa+ABZ*I_TWOBODYOVERLAP_I5xz_D2z_aa;
  Double I_TWOBODYOVERLAP_I4x2y_F3z_aa = I_TWOBODYOVERLAP_K4x2yz_D2z_aa+ABZ*I_TWOBODYOVERLAP_I4x2y_D2z_aa;
  Double I_TWOBODYOVERLAP_I4xyz_F3z_aa = I_TWOBODYOVERLAP_K4xy2z_D2z_aa+ABZ*I_TWOBODYOVERLAP_I4xyz_D2z_aa;
  Double I_TWOBODYOVERLAP_I4x2z_F3z_aa = I_TWOBODYOVERLAP_K4x3z_D2z_aa+ABZ*I_TWOBODYOVERLAP_I4x2z_D2z_aa;
  Double I_TWOBODYOVERLAP_I3x3y_F3z_aa = I_TWOBODYOVERLAP_K3x3yz_D2z_aa+ABZ*I_TWOBODYOVERLAP_I3x3y_D2z_aa;
  Double I_TWOBODYOVERLAP_I3x2yz_F3z_aa = I_TWOBODYOVERLAP_K3x2y2z_D2z_aa+ABZ*I_TWOBODYOVERLAP_I3x2yz_D2z_aa;
  Double I_TWOBODYOVERLAP_I3xy2z_F3z_aa = I_TWOBODYOVERLAP_K3xy3z_D2z_aa+ABZ*I_TWOBODYOVERLAP_I3xy2z_D2z_aa;
  Double I_TWOBODYOVERLAP_I3x3z_F3z_aa = I_TWOBODYOVERLAP_K3x4z_D2z_aa+ABZ*I_TWOBODYOVERLAP_I3x3z_D2z_aa;
  Double I_TWOBODYOVERLAP_I2x4y_F3z_aa = I_TWOBODYOVERLAP_K2x4yz_D2z_aa+ABZ*I_TWOBODYOVERLAP_I2x4y_D2z_aa;
  Double I_TWOBODYOVERLAP_I2x3yz_F3z_aa = I_TWOBODYOVERLAP_K2x3y2z_D2z_aa+ABZ*I_TWOBODYOVERLAP_I2x3yz_D2z_aa;
  Double I_TWOBODYOVERLAP_I2x2y2z_F3z_aa = I_TWOBODYOVERLAP_K2x2y3z_D2z_aa+ABZ*I_TWOBODYOVERLAP_I2x2y2z_D2z_aa;
  Double I_TWOBODYOVERLAP_I2xy3z_F3z_aa = I_TWOBODYOVERLAP_K2xy4z_D2z_aa+ABZ*I_TWOBODYOVERLAP_I2xy3z_D2z_aa;
  Double I_TWOBODYOVERLAP_I2x4z_F3z_aa = I_TWOBODYOVERLAP_K2x5z_D2z_aa+ABZ*I_TWOBODYOVERLAP_I2x4z_D2z_aa;
  Double I_TWOBODYOVERLAP_Ix5y_F3z_aa = I_TWOBODYOVERLAP_Kx5yz_D2z_aa+ABZ*I_TWOBODYOVERLAP_Ix5y_D2z_aa;
  Double I_TWOBODYOVERLAP_Ix4yz_F3z_aa = I_TWOBODYOVERLAP_Kx4y2z_D2z_aa+ABZ*I_TWOBODYOVERLAP_Ix4yz_D2z_aa;
  Double I_TWOBODYOVERLAP_Ix3y2z_F3z_aa = I_TWOBODYOVERLAP_Kx3y3z_D2z_aa+ABZ*I_TWOBODYOVERLAP_Ix3y2z_D2z_aa;
  Double I_TWOBODYOVERLAP_Ix2y3z_F3z_aa = I_TWOBODYOVERLAP_Kx2y4z_D2z_aa+ABZ*I_TWOBODYOVERLAP_Ix2y3z_D2z_aa;
  Double I_TWOBODYOVERLAP_Ixy4z_F3z_aa = I_TWOBODYOVERLAP_Kxy5z_D2z_aa+ABZ*I_TWOBODYOVERLAP_Ixy4z_D2z_aa;
  Double I_TWOBODYOVERLAP_Ix5z_F3z_aa = I_TWOBODYOVERLAP_Kx6z_D2z_aa+ABZ*I_TWOBODYOVERLAP_Ix5z_D2z_aa;
  Double I_TWOBODYOVERLAP_I6y_F3z_aa = I_TWOBODYOVERLAP_K6yz_D2z_aa+ABZ*I_TWOBODYOVERLAP_I6y_D2z_aa;
  Double I_TWOBODYOVERLAP_I5yz_F3z_aa = I_TWOBODYOVERLAP_K5y2z_D2z_aa+ABZ*I_TWOBODYOVERLAP_I5yz_D2z_aa;
  Double I_TWOBODYOVERLAP_I4y2z_F3z_aa = I_TWOBODYOVERLAP_K4y3z_D2z_aa+ABZ*I_TWOBODYOVERLAP_I4y2z_D2z_aa;
  Double I_TWOBODYOVERLAP_I3y3z_F3z_aa = I_TWOBODYOVERLAP_K3y4z_D2z_aa+ABZ*I_TWOBODYOVERLAP_I3y3z_D2z_aa;
  Double I_TWOBODYOVERLAP_I2y4z_F3z_aa = I_TWOBODYOVERLAP_K2y5z_D2z_aa+ABZ*I_TWOBODYOVERLAP_I2y4z_D2z_aa;
  Double I_TWOBODYOVERLAP_Iy5z_F3z_aa = I_TWOBODYOVERLAP_Ky6z_D2z_aa+ABZ*I_TWOBODYOVERLAP_Iy5z_D2z_aa;
  Double I_TWOBODYOVERLAP_I6z_F3z_aa = I_TWOBODYOVERLAP_K7z_D2z_aa+ABZ*I_TWOBODYOVERLAP_I6z_D2z_aa;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_F_dax_dax
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_F_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_F
   ************************************************************/
  abcd[0] = 4.0E0*I_TWOBODYOVERLAP_I6x_F3x_aa-2.0E0*4*I_TWOBODYOVERLAP_G4x_F3x_a-2.0E0*5*I_TWOBODYOVERLAP_G4x_F3x_a+4*3*I_TWOBODYOVERLAP_D2x_F3x;
  abcd[1] = 4.0E0*I_TWOBODYOVERLAP_I5xy_F3x_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xy_F3x_a-2.0E0*4*I_TWOBODYOVERLAP_G3xy_F3x_a+3*2*I_TWOBODYOVERLAP_Dxy_F3x;
  abcd[2] = 4.0E0*I_TWOBODYOVERLAP_I5xz_F3x_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xz_F3x_a-2.0E0*4*I_TWOBODYOVERLAP_G3xz_F3x_a+3*2*I_TWOBODYOVERLAP_Dxz_F3x;
  abcd[3] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_F3x_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F3x_a+2*1*I_TWOBODYOVERLAP_D2y_F3x;
  abcd[4] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3x_a-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_F3x_a+2*1*I_TWOBODYOVERLAP_Dyz_F3x;
  abcd[5] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_F3x_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F3x_a+2*1*I_TWOBODYOVERLAP_D2z_F3x;
  abcd[6] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F3x_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_F3x_a;
  abcd[7] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_F3x_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3x_a;
  abcd[8] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_F3x_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3x_a;
  abcd[9] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F3x_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_F3x_a;
  abcd[10] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_F3x_a;
  abcd[11] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F3x_a;
  abcd[12] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_F3x_a;
  abcd[13] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F3x_a;
  abcd[14] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_F3x_a;
  abcd[15] = 4.0E0*I_TWOBODYOVERLAP_I6x_F2xy_aa-2.0E0*4*I_TWOBODYOVERLAP_G4x_F2xy_a-2.0E0*5*I_TWOBODYOVERLAP_G4x_F2xy_a+4*3*I_TWOBODYOVERLAP_D2x_F2xy;
  abcd[16] = 4.0E0*I_TWOBODYOVERLAP_I5xy_F2xy_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xy_F2xy_a-2.0E0*4*I_TWOBODYOVERLAP_G3xy_F2xy_a+3*2*I_TWOBODYOVERLAP_Dxy_F2xy;
  abcd[17] = 4.0E0*I_TWOBODYOVERLAP_I5xz_F2xy_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xz_F2xy_a-2.0E0*4*I_TWOBODYOVERLAP_G3xz_F2xy_a+3*2*I_TWOBODYOVERLAP_Dxz_F2xy;
  abcd[18] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_F2xy_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F2xy_a+2*1*I_TWOBODYOVERLAP_D2y_F2xy;
  abcd[19] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2xy_a-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_F2xy_a+2*1*I_TWOBODYOVERLAP_Dyz_F2xy;
  abcd[20] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_F2xy_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F2xy_a+2*1*I_TWOBODYOVERLAP_D2z_F2xy;
  abcd[21] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F2xy_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_F2xy_a;
  abcd[22] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_F2xy_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2xy_a;
  abcd[23] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_F2xy_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2xy_a;
  abcd[24] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F2xy_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_F2xy_a;
  abcd[25] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_F2xy_a;
  abcd[26] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F2xy_a;
  abcd[27] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_F2xy_a;
  abcd[28] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F2xy_a;
  abcd[29] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_F2xy_a;
  abcd[30] = 4.0E0*I_TWOBODYOVERLAP_I6x_F2xz_aa-2.0E0*4*I_TWOBODYOVERLAP_G4x_F2xz_a-2.0E0*5*I_TWOBODYOVERLAP_G4x_F2xz_a+4*3*I_TWOBODYOVERLAP_D2x_F2xz;
  abcd[31] = 4.0E0*I_TWOBODYOVERLAP_I5xy_F2xz_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xy_F2xz_a-2.0E0*4*I_TWOBODYOVERLAP_G3xy_F2xz_a+3*2*I_TWOBODYOVERLAP_Dxy_F2xz;
  abcd[32] = 4.0E0*I_TWOBODYOVERLAP_I5xz_F2xz_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xz_F2xz_a-2.0E0*4*I_TWOBODYOVERLAP_G3xz_F2xz_a+3*2*I_TWOBODYOVERLAP_Dxz_F2xz;
  abcd[33] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_F2xz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F2xz_a+2*1*I_TWOBODYOVERLAP_D2y_F2xz;
  abcd[34] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2xz_a-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_F2xz_a+2*1*I_TWOBODYOVERLAP_Dyz_F2xz;
  abcd[35] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_F2xz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F2xz_a+2*1*I_TWOBODYOVERLAP_D2z_F2xz;
  abcd[36] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F2xz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_F2xz_a;
  abcd[37] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_F2xz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2xz_a;
  abcd[38] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_F2xz_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2xz_a;
  abcd[39] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F2xz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_F2xz_a;
  abcd[40] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_F2xz_a;
  abcd[41] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F2xz_a;
  abcd[42] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_F2xz_a;
  abcd[43] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F2xz_a;
  abcd[44] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_F2xz_a;
  abcd[45] = 4.0E0*I_TWOBODYOVERLAP_I6x_Fx2y_aa-2.0E0*4*I_TWOBODYOVERLAP_G4x_Fx2y_a-2.0E0*5*I_TWOBODYOVERLAP_G4x_Fx2y_a+4*3*I_TWOBODYOVERLAP_D2x_Fx2y;
  abcd[46] = 4.0E0*I_TWOBODYOVERLAP_I5xy_Fx2y_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xy_Fx2y_a-2.0E0*4*I_TWOBODYOVERLAP_G3xy_Fx2y_a+3*2*I_TWOBODYOVERLAP_Dxy_Fx2y;
  abcd[47] = 4.0E0*I_TWOBODYOVERLAP_I5xz_Fx2y_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xz_Fx2y_a-2.0E0*4*I_TWOBODYOVERLAP_G3xz_Fx2y_a+3*2*I_TWOBODYOVERLAP_Dxz_Fx2y;
  abcd[48] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_Fx2y_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Fx2y_a+2*1*I_TWOBODYOVERLAP_D2y_Fx2y;
  abcd[49] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fx2y_a-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Fx2y_a+2*1*I_TWOBODYOVERLAP_Dyz_Fx2y;
  abcd[50] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_Fx2y_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Fx2y_a+2*1*I_TWOBODYOVERLAP_D2z_Fx2y;
  abcd[51] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Fx2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_Fx2y_a;
  abcd[52] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_Fx2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fx2y_a;
  abcd[53] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_Fx2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fx2y_a;
  abcd[54] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Fx2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_Fx2y_a;
  abcd[55] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Fx2y_a;
  abcd[56] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Fx2y_a;
  abcd[57] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Fx2y_a;
  abcd[58] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Fx2y_a;
  abcd[59] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_Fx2y_a;
  abcd[60] = 4.0E0*I_TWOBODYOVERLAP_I6x_Fxyz_aa-2.0E0*4*I_TWOBODYOVERLAP_G4x_Fxyz_a-2.0E0*5*I_TWOBODYOVERLAP_G4x_Fxyz_a+4*3*I_TWOBODYOVERLAP_D2x_Fxyz;
  abcd[61] = 4.0E0*I_TWOBODYOVERLAP_I5xy_Fxyz_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xy_Fxyz_a-2.0E0*4*I_TWOBODYOVERLAP_G3xy_Fxyz_a+3*2*I_TWOBODYOVERLAP_Dxy_Fxyz;
  abcd[62] = 4.0E0*I_TWOBODYOVERLAP_I5xz_Fxyz_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xz_Fxyz_a-2.0E0*4*I_TWOBODYOVERLAP_G3xz_Fxyz_a+3*2*I_TWOBODYOVERLAP_Dxz_Fxyz;
  abcd[63] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_Fxyz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Fxyz_a+2*1*I_TWOBODYOVERLAP_D2y_Fxyz;
  abcd[64] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fxyz_a-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Fxyz_a+2*1*I_TWOBODYOVERLAP_Dyz_Fxyz;
  abcd[65] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_Fxyz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Fxyz_a+2*1*I_TWOBODYOVERLAP_D2z_Fxyz;
  abcd[66] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Fxyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_Fxyz_a;
  abcd[67] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_Fxyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fxyz_a;
  abcd[68] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_Fxyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fxyz_a;
  abcd[69] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Fxyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_Fxyz_a;
  abcd[70] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Fxyz_a;
  abcd[71] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Fxyz_a;
  abcd[72] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Fxyz_a;
  abcd[73] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Fxyz_a;
  abcd[74] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_Fxyz_a;
  abcd[75] = 4.0E0*I_TWOBODYOVERLAP_I6x_Fx2z_aa-2.0E0*4*I_TWOBODYOVERLAP_G4x_Fx2z_a-2.0E0*5*I_TWOBODYOVERLAP_G4x_Fx2z_a+4*3*I_TWOBODYOVERLAP_D2x_Fx2z;
  abcd[76] = 4.0E0*I_TWOBODYOVERLAP_I5xy_Fx2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xy_Fx2z_a-2.0E0*4*I_TWOBODYOVERLAP_G3xy_Fx2z_a+3*2*I_TWOBODYOVERLAP_Dxy_Fx2z;
  abcd[77] = 4.0E0*I_TWOBODYOVERLAP_I5xz_Fx2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xz_Fx2z_a-2.0E0*4*I_TWOBODYOVERLAP_G3xz_Fx2z_a+3*2*I_TWOBODYOVERLAP_Dxz_Fx2z;
  abcd[78] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_Fx2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Fx2z_a+2*1*I_TWOBODYOVERLAP_D2y_Fx2z;
  abcd[79] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fx2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Fx2z_a+2*1*I_TWOBODYOVERLAP_Dyz_Fx2z;
  abcd[80] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_Fx2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Fx2z_a+2*1*I_TWOBODYOVERLAP_D2z_Fx2z;
  abcd[81] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Fx2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_Fx2z_a;
  abcd[82] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_Fx2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fx2z_a;
  abcd[83] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_Fx2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fx2z_a;
  abcd[84] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Fx2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_Fx2z_a;
  abcd[85] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Fx2z_a;
  abcd[86] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Fx2z_a;
  abcd[87] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Fx2z_a;
  abcd[88] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Fx2z_a;
  abcd[89] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_Fx2z_a;
  abcd[90] = 4.0E0*I_TWOBODYOVERLAP_I6x_F3y_aa-2.0E0*4*I_TWOBODYOVERLAP_G4x_F3y_a-2.0E0*5*I_TWOBODYOVERLAP_G4x_F3y_a+4*3*I_TWOBODYOVERLAP_D2x_F3y;
  abcd[91] = 4.0E0*I_TWOBODYOVERLAP_I5xy_F3y_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xy_F3y_a-2.0E0*4*I_TWOBODYOVERLAP_G3xy_F3y_a+3*2*I_TWOBODYOVERLAP_Dxy_F3y;
  abcd[92] = 4.0E0*I_TWOBODYOVERLAP_I5xz_F3y_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xz_F3y_a-2.0E0*4*I_TWOBODYOVERLAP_G3xz_F3y_a+3*2*I_TWOBODYOVERLAP_Dxz_F3y;
  abcd[93] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_F3y_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F3y_a+2*1*I_TWOBODYOVERLAP_D2y_F3y;
  abcd[94] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3y_a-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_F3y_a+2*1*I_TWOBODYOVERLAP_Dyz_F3y;
  abcd[95] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_F3y_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F3y_a+2*1*I_TWOBODYOVERLAP_D2z_F3y;
  abcd[96] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F3y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_F3y_a;
  abcd[97] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_F3y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3y_a;
  abcd[98] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_F3y_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3y_a;
  abcd[99] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F3y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_F3y_a;
  abcd[100] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_F3y_a;
  abcd[101] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F3y_a;
  abcd[102] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_F3y_a;
  abcd[103] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F3y_a;
  abcd[104] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_F3y_a;
  abcd[105] = 4.0E0*I_TWOBODYOVERLAP_I6x_F2yz_aa-2.0E0*4*I_TWOBODYOVERLAP_G4x_F2yz_a-2.0E0*5*I_TWOBODYOVERLAP_G4x_F2yz_a+4*3*I_TWOBODYOVERLAP_D2x_F2yz;
  abcd[106] = 4.0E0*I_TWOBODYOVERLAP_I5xy_F2yz_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xy_F2yz_a-2.0E0*4*I_TWOBODYOVERLAP_G3xy_F2yz_a+3*2*I_TWOBODYOVERLAP_Dxy_F2yz;
  abcd[107] = 4.0E0*I_TWOBODYOVERLAP_I5xz_F2yz_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xz_F2yz_a-2.0E0*4*I_TWOBODYOVERLAP_G3xz_F2yz_a+3*2*I_TWOBODYOVERLAP_Dxz_F2yz;
  abcd[108] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_F2yz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F2yz_a+2*1*I_TWOBODYOVERLAP_D2y_F2yz;
  abcd[109] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2yz_a-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_F2yz_a+2*1*I_TWOBODYOVERLAP_Dyz_F2yz;
  abcd[110] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_F2yz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F2yz_a+2*1*I_TWOBODYOVERLAP_D2z_F2yz;
  abcd[111] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F2yz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_F2yz_a;
  abcd[112] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_F2yz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2yz_a;
  abcd[113] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_F2yz_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2yz_a;
  abcd[114] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F2yz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_F2yz_a;
  abcd[115] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_F2yz_a;
  abcd[116] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F2yz_a;
  abcd[117] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_F2yz_a;
  abcd[118] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F2yz_a;
  abcd[119] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_F2yz_a;
  abcd[120] = 4.0E0*I_TWOBODYOVERLAP_I6x_Fy2z_aa-2.0E0*4*I_TWOBODYOVERLAP_G4x_Fy2z_a-2.0E0*5*I_TWOBODYOVERLAP_G4x_Fy2z_a+4*3*I_TWOBODYOVERLAP_D2x_Fy2z;
  abcd[121] = 4.0E0*I_TWOBODYOVERLAP_I5xy_Fy2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xy_Fy2z_a-2.0E0*4*I_TWOBODYOVERLAP_G3xy_Fy2z_a+3*2*I_TWOBODYOVERLAP_Dxy_Fy2z;
  abcd[122] = 4.0E0*I_TWOBODYOVERLAP_I5xz_Fy2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xz_Fy2z_a-2.0E0*4*I_TWOBODYOVERLAP_G3xz_Fy2z_a+3*2*I_TWOBODYOVERLAP_Dxz_Fy2z;
  abcd[123] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_Fy2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Fy2z_a+2*1*I_TWOBODYOVERLAP_D2y_Fy2z;
  abcd[124] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fy2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Fy2z_a+2*1*I_TWOBODYOVERLAP_Dyz_Fy2z;
  abcd[125] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_Fy2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Fy2z_a+2*1*I_TWOBODYOVERLAP_D2z_Fy2z;
  abcd[126] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Fy2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_Fy2z_a;
  abcd[127] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_Fy2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fy2z_a;
  abcd[128] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_Fy2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fy2z_a;
  abcd[129] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Fy2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_Fy2z_a;
  abcd[130] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Fy2z_a;
  abcd[131] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Fy2z_a;
  abcd[132] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Fy2z_a;
  abcd[133] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Fy2z_a;
  abcd[134] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_Fy2z_a;
  abcd[135] = 4.0E0*I_TWOBODYOVERLAP_I6x_F3z_aa-2.0E0*4*I_TWOBODYOVERLAP_G4x_F3z_a-2.0E0*5*I_TWOBODYOVERLAP_G4x_F3z_a+4*3*I_TWOBODYOVERLAP_D2x_F3z;
  abcd[136] = 4.0E0*I_TWOBODYOVERLAP_I5xy_F3z_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xy_F3z_a-2.0E0*4*I_TWOBODYOVERLAP_G3xy_F3z_a+3*2*I_TWOBODYOVERLAP_Dxy_F3z;
  abcd[137] = 4.0E0*I_TWOBODYOVERLAP_I5xz_F3z_aa-2.0E0*3*I_TWOBODYOVERLAP_G3xz_F3z_a-2.0E0*4*I_TWOBODYOVERLAP_G3xz_F3z_a+3*2*I_TWOBODYOVERLAP_Dxz_F3z;
  abcd[138] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_F3z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F3z_a+2*1*I_TWOBODYOVERLAP_D2y_F3z;
  abcd[139] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3z_a-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_F3z_a+2*1*I_TWOBODYOVERLAP_Dyz_F3z;
  abcd[140] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_F3z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F3z_a+2*1*I_TWOBODYOVERLAP_D2z_F3z;
  abcd[141] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F3z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_F3z_a;
  abcd[142] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_F3z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3z_a;
  abcd[143] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_F3z_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3z_a;
  abcd[144] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F3z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_F3z_a;
  abcd[145] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_F3z_a;
  abcd[146] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F3z_a;
  abcd[147] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_F3z_a;
  abcd[148] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F3z_a;
  abcd[149] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_F_dax_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_F_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_F
   ************************************************************/
  abcd[150] = 4.0E0*I_TWOBODYOVERLAP_I5xy_F3x_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xy_F3x_a;
  abcd[151] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F3x_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F3x_a+3*1*I_TWOBODYOVERLAP_D2x_F3x;
  abcd[152] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F3x_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_F3x_a;
  abcd[153] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xy_F3x_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_F3x_a+2*2*I_TWOBODYOVERLAP_Dxy_F3x;
  abcd[154] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F3x_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3x_a+2*1*I_TWOBODYOVERLAP_Dxz_F3x;
  abcd[155] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3x_a;
  abcd[156] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_F3x_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F3x_a-2.0E0*1*I_TWOBODYOVERLAP_G4y_F3x_a+3*I_TWOBODYOVERLAP_D2y_F3x;
  abcd[157] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3x_a-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F3x_a+2*I_TWOBODYOVERLAP_Dyz_F3x;
  abcd[158] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_F3x_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_F3x_a+1*I_TWOBODYOVERLAP_D2z_F3x;
  abcd[159] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F3x_a;
  abcd[160] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_F3x_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_F3x_a;
  abcd[161] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F3x_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_F3x_a;
  abcd[162] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3x_a;
  abcd[163] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F3x_a;
  abcd[164] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F3x_aa;
  abcd[165] = 4.0E0*I_TWOBODYOVERLAP_I5xy_F2xy_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xy_F2xy_a;
  abcd[166] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F2xy_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F2xy_a+3*1*I_TWOBODYOVERLAP_D2x_F2xy;
  abcd[167] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F2xy_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_F2xy_a;
  abcd[168] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xy_F2xy_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_F2xy_a+2*2*I_TWOBODYOVERLAP_Dxy_F2xy;
  abcd[169] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F2xy_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2xy_a+2*1*I_TWOBODYOVERLAP_Dxz_F2xy;
  abcd[170] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2xy_a;
  abcd[171] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_F2xy_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F2xy_a-2.0E0*1*I_TWOBODYOVERLAP_G4y_F2xy_a+3*I_TWOBODYOVERLAP_D2y_F2xy;
  abcd[172] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2xy_a-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F2xy_a+2*I_TWOBODYOVERLAP_Dyz_F2xy;
  abcd[173] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_F2xy_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_F2xy_a+1*I_TWOBODYOVERLAP_D2z_F2xy;
  abcd[174] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F2xy_a;
  abcd[175] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_F2xy_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_F2xy_a;
  abcd[176] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F2xy_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_F2xy_a;
  abcd[177] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2xy_a;
  abcd[178] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F2xy_a;
  abcd[179] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F2xy_aa;
  abcd[180] = 4.0E0*I_TWOBODYOVERLAP_I5xy_F2xz_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xy_F2xz_a;
  abcd[181] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F2xz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F2xz_a+3*1*I_TWOBODYOVERLAP_D2x_F2xz;
  abcd[182] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F2xz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_F2xz_a;
  abcd[183] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xy_F2xz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_F2xz_a+2*2*I_TWOBODYOVERLAP_Dxy_F2xz;
  abcd[184] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F2xz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2xz_a+2*1*I_TWOBODYOVERLAP_Dxz_F2xz;
  abcd[185] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2xz_a;
  abcd[186] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_F2xz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F2xz_a-2.0E0*1*I_TWOBODYOVERLAP_G4y_F2xz_a+3*I_TWOBODYOVERLAP_D2y_F2xz;
  abcd[187] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2xz_a-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F2xz_a+2*I_TWOBODYOVERLAP_Dyz_F2xz;
  abcd[188] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_F2xz_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_F2xz_a+1*I_TWOBODYOVERLAP_D2z_F2xz;
  abcd[189] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F2xz_a;
  abcd[190] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_F2xz_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_F2xz_a;
  abcd[191] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F2xz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_F2xz_a;
  abcd[192] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2xz_a;
  abcd[193] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F2xz_a;
  abcd[194] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F2xz_aa;
  abcd[195] = 4.0E0*I_TWOBODYOVERLAP_I5xy_Fx2y_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xy_Fx2y_a;
  abcd[196] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Fx2y_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Fx2y_a+3*1*I_TWOBODYOVERLAP_D2x_Fx2y;
  abcd[197] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Fx2y_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Fx2y_a;
  abcd[198] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xy_Fx2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_Fx2y_a+2*2*I_TWOBODYOVERLAP_Dxy_Fx2y;
  abcd[199] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Fx2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fx2y_a+2*1*I_TWOBODYOVERLAP_Dxz_Fx2y;
  abcd[200] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fx2y_a;
  abcd[201] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Fx2y_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Fx2y_a-2.0E0*1*I_TWOBODYOVERLAP_G4y_Fx2y_a+3*I_TWOBODYOVERLAP_D2y_Fx2y;
  abcd[202] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fx2y_a-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Fx2y_a+2*I_TWOBODYOVERLAP_Dyz_Fx2y;
  abcd[203] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Fx2y_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Fx2y_a+1*I_TWOBODYOVERLAP_D2z_Fx2y;
  abcd[204] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Fx2y_a;
  abcd[205] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_Fx2y_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_Fx2y_a;
  abcd[206] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Fx2y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Fx2y_a;
  abcd[207] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fx2y_a;
  abcd[208] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Fx2y_a;
  abcd[209] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Fx2y_aa;
  abcd[210] = 4.0E0*I_TWOBODYOVERLAP_I5xy_Fxyz_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xy_Fxyz_a;
  abcd[211] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Fxyz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Fxyz_a+3*1*I_TWOBODYOVERLAP_D2x_Fxyz;
  abcd[212] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Fxyz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Fxyz_a;
  abcd[213] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xy_Fxyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_Fxyz_a+2*2*I_TWOBODYOVERLAP_Dxy_Fxyz;
  abcd[214] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Fxyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fxyz_a+2*1*I_TWOBODYOVERLAP_Dxz_Fxyz;
  abcd[215] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fxyz_a;
  abcd[216] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Fxyz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Fxyz_a-2.0E0*1*I_TWOBODYOVERLAP_G4y_Fxyz_a+3*I_TWOBODYOVERLAP_D2y_Fxyz;
  abcd[217] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fxyz_a-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Fxyz_a+2*I_TWOBODYOVERLAP_Dyz_Fxyz;
  abcd[218] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Fxyz_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Fxyz_a+1*I_TWOBODYOVERLAP_D2z_Fxyz;
  abcd[219] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Fxyz_a;
  abcd[220] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_Fxyz_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_Fxyz_a;
  abcd[221] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Fxyz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Fxyz_a;
  abcd[222] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fxyz_a;
  abcd[223] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Fxyz_a;
  abcd[224] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Fxyz_aa;
  abcd[225] = 4.0E0*I_TWOBODYOVERLAP_I5xy_Fx2z_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xy_Fx2z_a;
  abcd[226] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Fx2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Fx2z_a+3*1*I_TWOBODYOVERLAP_D2x_Fx2z;
  abcd[227] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Fx2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Fx2z_a;
  abcd[228] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xy_Fx2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_Fx2z_a+2*2*I_TWOBODYOVERLAP_Dxy_Fx2z;
  abcd[229] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Fx2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fx2z_a+2*1*I_TWOBODYOVERLAP_Dxz_Fx2z;
  abcd[230] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fx2z_a;
  abcd[231] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Fx2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Fx2z_a-2.0E0*1*I_TWOBODYOVERLAP_G4y_Fx2z_a+3*I_TWOBODYOVERLAP_D2y_Fx2z;
  abcd[232] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fx2z_a-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Fx2z_a+2*I_TWOBODYOVERLAP_Dyz_Fx2z;
  abcd[233] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Fx2z_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Fx2z_a+1*I_TWOBODYOVERLAP_D2z_Fx2z;
  abcd[234] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Fx2z_a;
  abcd[235] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_Fx2z_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_Fx2z_a;
  abcd[236] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Fx2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Fx2z_a;
  abcd[237] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fx2z_a;
  abcd[238] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Fx2z_a;
  abcd[239] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Fx2z_aa;
  abcd[240] = 4.0E0*I_TWOBODYOVERLAP_I5xy_F3y_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xy_F3y_a;
  abcd[241] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F3y_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F3y_a+3*1*I_TWOBODYOVERLAP_D2x_F3y;
  abcd[242] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F3y_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_F3y_a;
  abcd[243] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xy_F3y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_F3y_a+2*2*I_TWOBODYOVERLAP_Dxy_F3y;
  abcd[244] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F3y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3y_a+2*1*I_TWOBODYOVERLAP_Dxz_F3y;
  abcd[245] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3y_a;
  abcd[246] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_F3y_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F3y_a-2.0E0*1*I_TWOBODYOVERLAP_G4y_F3y_a+3*I_TWOBODYOVERLAP_D2y_F3y;
  abcd[247] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3y_a-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F3y_a+2*I_TWOBODYOVERLAP_Dyz_F3y;
  abcd[248] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_F3y_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_F3y_a+1*I_TWOBODYOVERLAP_D2z_F3y;
  abcd[249] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F3y_a;
  abcd[250] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_F3y_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_F3y_a;
  abcd[251] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F3y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_F3y_a;
  abcd[252] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3y_a;
  abcd[253] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F3y_a;
  abcd[254] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F3y_aa;
  abcd[255] = 4.0E0*I_TWOBODYOVERLAP_I5xy_F2yz_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xy_F2yz_a;
  abcd[256] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F2yz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F2yz_a+3*1*I_TWOBODYOVERLAP_D2x_F2yz;
  abcd[257] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F2yz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_F2yz_a;
  abcd[258] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xy_F2yz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_F2yz_a+2*2*I_TWOBODYOVERLAP_Dxy_F2yz;
  abcd[259] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F2yz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2yz_a+2*1*I_TWOBODYOVERLAP_Dxz_F2yz;
  abcd[260] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2yz_a;
  abcd[261] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_F2yz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F2yz_a-2.0E0*1*I_TWOBODYOVERLAP_G4y_F2yz_a+3*I_TWOBODYOVERLAP_D2y_F2yz;
  abcd[262] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2yz_a-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F2yz_a+2*I_TWOBODYOVERLAP_Dyz_F2yz;
  abcd[263] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_F2yz_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_F2yz_a+1*I_TWOBODYOVERLAP_D2z_F2yz;
  abcd[264] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F2yz_a;
  abcd[265] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_F2yz_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_F2yz_a;
  abcd[266] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F2yz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_F2yz_a;
  abcd[267] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2yz_a;
  abcd[268] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F2yz_a;
  abcd[269] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F2yz_aa;
  abcd[270] = 4.0E0*I_TWOBODYOVERLAP_I5xy_Fy2z_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xy_Fy2z_a;
  abcd[271] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Fy2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Fy2z_a+3*1*I_TWOBODYOVERLAP_D2x_Fy2z;
  abcd[272] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Fy2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Fy2z_a;
  abcd[273] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xy_Fy2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_Fy2z_a+2*2*I_TWOBODYOVERLAP_Dxy_Fy2z;
  abcd[274] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Fy2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fy2z_a+2*1*I_TWOBODYOVERLAP_Dxz_Fy2z;
  abcd[275] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fy2z_a;
  abcd[276] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Fy2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Fy2z_a-2.0E0*1*I_TWOBODYOVERLAP_G4y_Fy2z_a+3*I_TWOBODYOVERLAP_D2y_Fy2z;
  abcd[277] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fy2z_a-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Fy2z_a+2*I_TWOBODYOVERLAP_Dyz_Fy2z;
  abcd[278] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Fy2z_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Fy2z_a+1*I_TWOBODYOVERLAP_D2z_Fy2z;
  abcd[279] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Fy2z_a;
  abcd[280] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_Fy2z_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_Fy2z_a;
  abcd[281] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Fy2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Fy2z_a;
  abcd[282] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fy2z_a;
  abcd[283] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Fy2z_a;
  abcd[284] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Fy2z_aa;
  abcd[285] = 4.0E0*I_TWOBODYOVERLAP_I5xy_F3z_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xy_F3z_a;
  abcd[286] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F3z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F3z_a+3*1*I_TWOBODYOVERLAP_D2x_F3z;
  abcd[287] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F3z_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_F3z_a;
  abcd[288] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xy_F3z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3y_F3z_a+2*2*I_TWOBODYOVERLAP_Dxy_F3z;
  abcd[289] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F3z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3z_a+2*1*I_TWOBODYOVERLAP_Dxz_F3z;
  abcd[290] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3z_a;
  abcd[291] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_F3z_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F3z_a-2.0E0*1*I_TWOBODYOVERLAP_G4y_F3z_a+3*I_TWOBODYOVERLAP_D2y_F3z;
  abcd[292] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3z_a-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F3z_a+2*I_TWOBODYOVERLAP_Dyz_F3z;
  abcd[293] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_F3z_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_F3z_a+1*I_TWOBODYOVERLAP_D2z_F3z;
  abcd[294] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F3z_a;
  abcd[295] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_F3z_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_F3z_a;
  abcd[296] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F3z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_F3z_a;
  abcd[297] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3z_a;
  abcd[298] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F3z_a;
  abcd[299] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F3z_aa;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_F_dax_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_F_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_F
   ************************************************************/
  abcd[300] = 4.0E0*I_TWOBODYOVERLAP_I5xz_F3x_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xz_F3x_a;
  abcd[301] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F3x_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_F3x_a;
  abcd[302] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F3x_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F3x_a+3*1*I_TWOBODYOVERLAP_D2x_F3x;
  abcd[303] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3x_a;
  abcd[304] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F3x_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3x_a+2*1*I_TWOBODYOVERLAP_Dxy_F3x;
  abcd[305] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xz_F3x_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_F3x_a+2*2*I_TWOBODYOVERLAP_Dxz_F3x;
  abcd[306] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F3x_a;
  abcd[307] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_F3x_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_F3x_a+1*I_TWOBODYOVERLAP_D2y_F3x;
  abcd[308] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3x_a-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F3x_a+2*I_TWOBODYOVERLAP_Dyz_F3x;
  abcd[309] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_F3x_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F3x_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_F3x_a+3*I_TWOBODYOVERLAP_D2z_F3x;
  abcd[310] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F3x_aa;
  abcd[311] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F3x_a;
  abcd[312] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3x_a;
  abcd[313] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F3x_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_F3x_a;
  abcd[314] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_F3x_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_F3x_a;
  abcd[315] = 4.0E0*I_TWOBODYOVERLAP_I5xz_F2xy_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xz_F2xy_a;
  abcd[316] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F2xy_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_F2xy_a;
  abcd[317] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F2xy_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F2xy_a+3*1*I_TWOBODYOVERLAP_D2x_F2xy;
  abcd[318] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2xy_a;
  abcd[319] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F2xy_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2xy_a+2*1*I_TWOBODYOVERLAP_Dxy_F2xy;
  abcd[320] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xz_F2xy_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_F2xy_a+2*2*I_TWOBODYOVERLAP_Dxz_F2xy;
  abcd[321] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F2xy_a;
  abcd[322] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_F2xy_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_F2xy_a+1*I_TWOBODYOVERLAP_D2y_F2xy;
  abcd[323] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2xy_a-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F2xy_a+2*I_TWOBODYOVERLAP_Dyz_F2xy;
  abcd[324] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_F2xy_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F2xy_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_F2xy_a+3*I_TWOBODYOVERLAP_D2z_F2xy;
  abcd[325] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F2xy_aa;
  abcd[326] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F2xy_a;
  abcd[327] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2xy_a;
  abcd[328] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F2xy_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_F2xy_a;
  abcd[329] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_F2xy_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_F2xy_a;
  abcd[330] = 4.0E0*I_TWOBODYOVERLAP_I5xz_F2xz_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xz_F2xz_a;
  abcd[331] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F2xz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_F2xz_a;
  abcd[332] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F2xz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F2xz_a+3*1*I_TWOBODYOVERLAP_D2x_F2xz;
  abcd[333] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2xz_a;
  abcd[334] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F2xz_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2xz_a+2*1*I_TWOBODYOVERLAP_Dxy_F2xz;
  abcd[335] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xz_F2xz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_F2xz_a+2*2*I_TWOBODYOVERLAP_Dxz_F2xz;
  abcd[336] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F2xz_a;
  abcd[337] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_F2xz_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_F2xz_a+1*I_TWOBODYOVERLAP_D2y_F2xz;
  abcd[338] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2xz_a-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F2xz_a+2*I_TWOBODYOVERLAP_Dyz_F2xz;
  abcd[339] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_F2xz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F2xz_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_F2xz_a+3*I_TWOBODYOVERLAP_D2z_F2xz;
  abcd[340] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F2xz_aa;
  abcd[341] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F2xz_a;
  abcd[342] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2xz_a;
  abcd[343] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F2xz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_F2xz_a;
  abcd[344] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_F2xz_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_F2xz_a;
  abcd[345] = 4.0E0*I_TWOBODYOVERLAP_I5xz_Fx2y_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xz_Fx2y_a;
  abcd[346] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Fx2y_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Fx2y_a;
  abcd[347] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Fx2y_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Fx2y_a+3*1*I_TWOBODYOVERLAP_D2x_Fx2y;
  abcd[348] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fx2y_a;
  abcd[349] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Fx2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fx2y_a+2*1*I_TWOBODYOVERLAP_Dxy_Fx2y;
  abcd[350] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xz_Fx2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_Fx2y_a+2*2*I_TWOBODYOVERLAP_Dxz_Fx2y;
  abcd[351] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Fx2y_a;
  abcd[352] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Fx2y_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Fx2y_a+1*I_TWOBODYOVERLAP_D2y_Fx2y;
  abcd[353] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fx2y_a-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Fx2y_a+2*I_TWOBODYOVERLAP_Dyz_Fx2y;
  abcd[354] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Fx2y_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Fx2y_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_Fx2y_a+3*I_TWOBODYOVERLAP_D2z_Fx2y;
  abcd[355] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Fx2y_aa;
  abcd[356] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Fx2y_a;
  abcd[357] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fx2y_a;
  abcd[358] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Fx2y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Fx2y_a;
  abcd[359] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_Fx2y_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_Fx2y_a;
  abcd[360] = 4.0E0*I_TWOBODYOVERLAP_I5xz_Fxyz_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xz_Fxyz_a;
  abcd[361] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Fxyz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Fxyz_a;
  abcd[362] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Fxyz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Fxyz_a+3*1*I_TWOBODYOVERLAP_D2x_Fxyz;
  abcd[363] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fxyz_a;
  abcd[364] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Fxyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fxyz_a+2*1*I_TWOBODYOVERLAP_Dxy_Fxyz;
  abcd[365] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xz_Fxyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_Fxyz_a+2*2*I_TWOBODYOVERLAP_Dxz_Fxyz;
  abcd[366] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Fxyz_a;
  abcd[367] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Fxyz_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Fxyz_a+1*I_TWOBODYOVERLAP_D2y_Fxyz;
  abcd[368] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fxyz_a-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Fxyz_a+2*I_TWOBODYOVERLAP_Dyz_Fxyz;
  abcd[369] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Fxyz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Fxyz_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_Fxyz_a+3*I_TWOBODYOVERLAP_D2z_Fxyz;
  abcd[370] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Fxyz_aa;
  abcd[371] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Fxyz_a;
  abcd[372] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fxyz_a;
  abcd[373] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Fxyz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Fxyz_a;
  abcd[374] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_Fxyz_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_Fxyz_a;
  abcd[375] = 4.0E0*I_TWOBODYOVERLAP_I5xz_Fx2z_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xz_Fx2z_a;
  abcd[376] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Fx2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Fx2z_a;
  abcd[377] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Fx2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Fx2z_a+3*1*I_TWOBODYOVERLAP_D2x_Fx2z;
  abcd[378] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fx2z_a;
  abcd[379] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Fx2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fx2z_a+2*1*I_TWOBODYOVERLAP_Dxy_Fx2z;
  abcd[380] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xz_Fx2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_Fx2z_a+2*2*I_TWOBODYOVERLAP_Dxz_Fx2z;
  abcd[381] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Fx2z_a;
  abcd[382] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Fx2z_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Fx2z_a+1*I_TWOBODYOVERLAP_D2y_Fx2z;
  abcd[383] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fx2z_a-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Fx2z_a+2*I_TWOBODYOVERLAP_Dyz_Fx2z;
  abcd[384] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Fx2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Fx2z_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_Fx2z_a+3*I_TWOBODYOVERLAP_D2z_Fx2z;
  abcd[385] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Fx2z_aa;
  abcd[386] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Fx2z_a;
  abcd[387] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fx2z_a;
  abcd[388] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Fx2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Fx2z_a;
  abcd[389] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_Fx2z_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_Fx2z_a;
  abcd[390] = 4.0E0*I_TWOBODYOVERLAP_I5xz_F3y_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xz_F3y_a;
  abcd[391] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F3y_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_F3y_a;
  abcd[392] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F3y_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F3y_a+3*1*I_TWOBODYOVERLAP_D2x_F3y;
  abcd[393] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3y_a;
  abcd[394] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F3y_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3y_a+2*1*I_TWOBODYOVERLAP_Dxy_F3y;
  abcd[395] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xz_F3y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_F3y_a+2*2*I_TWOBODYOVERLAP_Dxz_F3y;
  abcd[396] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F3y_a;
  abcd[397] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_F3y_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_F3y_a+1*I_TWOBODYOVERLAP_D2y_F3y;
  abcd[398] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3y_a-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F3y_a+2*I_TWOBODYOVERLAP_Dyz_F3y;
  abcd[399] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_F3y_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F3y_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_F3y_a+3*I_TWOBODYOVERLAP_D2z_F3y;
  abcd[400] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F3y_aa;
  abcd[401] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F3y_a;
  abcd[402] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3y_a;
  abcd[403] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F3y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_F3y_a;
  abcd[404] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_F3y_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_F3y_a;
  abcd[405] = 4.0E0*I_TWOBODYOVERLAP_I5xz_F2yz_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xz_F2yz_a;
  abcd[406] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F2yz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_F2yz_a;
  abcd[407] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F2yz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F2yz_a+3*1*I_TWOBODYOVERLAP_D2x_F2yz;
  abcd[408] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2yz_a;
  abcd[409] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F2yz_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2yz_a+2*1*I_TWOBODYOVERLAP_Dxy_F2yz;
  abcd[410] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xz_F2yz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_F2yz_a+2*2*I_TWOBODYOVERLAP_Dxz_F2yz;
  abcd[411] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F2yz_a;
  abcd[412] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_F2yz_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_F2yz_a+1*I_TWOBODYOVERLAP_D2y_F2yz;
  abcd[413] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2yz_a-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F2yz_a+2*I_TWOBODYOVERLAP_Dyz_F2yz;
  abcd[414] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_F2yz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F2yz_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_F2yz_a+3*I_TWOBODYOVERLAP_D2z_F2yz;
  abcd[415] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F2yz_aa;
  abcd[416] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F2yz_a;
  abcd[417] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2yz_a;
  abcd[418] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F2yz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_F2yz_a;
  abcd[419] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_F2yz_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_F2yz_a;
  abcd[420] = 4.0E0*I_TWOBODYOVERLAP_I5xz_Fy2z_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xz_Fy2z_a;
  abcd[421] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Fy2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_Fy2z_a;
  abcd[422] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Fy2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Fy2z_a+3*1*I_TWOBODYOVERLAP_D2x_Fy2z;
  abcd[423] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fy2z_a;
  abcd[424] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Fy2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fy2z_a+2*1*I_TWOBODYOVERLAP_Dxy_Fy2z;
  abcd[425] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xz_Fy2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_Fy2z_a+2*2*I_TWOBODYOVERLAP_Dxz_Fy2z;
  abcd[426] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Fy2z_a;
  abcd[427] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Fy2z_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_Fy2z_a+1*I_TWOBODYOVERLAP_D2y_Fy2z;
  abcd[428] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fy2z_a-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Fy2z_a+2*I_TWOBODYOVERLAP_Dyz_Fy2z;
  abcd[429] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Fy2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Fy2z_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_Fy2z_a+3*I_TWOBODYOVERLAP_D2z_Fy2z;
  abcd[430] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Fy2z_aa;
  abcd[431] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Fy2z_a;
  abcd[432] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fy2z_a;
  abcd[433] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Fy2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Fy2z_a;
  abcd[434] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_Fy2z_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_Fy2z_a;
  abcd[435] = 4.0E0*I_TWOBODYOVERLAP_I5xz_F3z_aa-2.0E0*4*I_TWOBODYOVERLAP_G3xz_F3z_a;
  abcd[436] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F3z_aa-2.0E0*3*I_TWOBODYOVERLAP_G2xyz_F3z_a;
  abcd[437] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F3z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F3z_a+3*1*I_TWOBODYOVERLAP_D2x_F3z;
  abcd[438] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3z_a;
  abcd[439] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F3z_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3z_a+2*1*I_TWOBODYOVERLAP_Dxy_F3z;
  abcd[440] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_G3xz_F3z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx3z_F3z_a+2*2*I_TWOBODYOVERLAP_Dxz_F3z;
  abcd[441] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F3z_a;
  abcd[442] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_F3z_a-2.0E0*1*I_TWOBODYOVERLAP_G2y2z_F3z_a+1*I_TWOBODYOVERLAP_D2y_F3z;
  abcd[443] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3z_a-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F3z_a+2*I_TWOBODYOVERLAP_Dyz_F3z;
  abcd[444] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_F3z_aa-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F3z_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_F3z_a+3*I_TWOBODYOVERLAP_D2z_F3z;
  abcd[445] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F3z_aa;
  abcd[446] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F3z_a;
  abcd[447] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3z_a;
  abcd[448] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F3z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_F3z_a;
  abcd[449] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_F3z_aa-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_F_day_day
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_F_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_F
   ************************************************************/
  abcd[450] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F3x_a;
  abcd[451] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F3x_a-2.0E0*2*I_TWOBODYOVERLAP_G3xy_F3x_a;
  abcd[452] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F3x_a;
  abcd[453] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_F3x_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F3x_a+2*1*I_TWOBODYOVERLAP_D2x_F3x;
  abcd[454] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_F3x_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3x_a;
  abcd[455] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_F3x_a;
  abcd[456] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_F3x_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3y_F3x_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_F3x_a+3*2*I_TWOBODYOVERLAP_Dxy_F3x;
  abcd[457] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3x_a-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_F3x_a+2*1*I_TWOBODYOVERLAP_Dxz_F3x;
  abcd[458] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_F3x_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3x_a;
  abcd[459] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F3x_a;
  abcd[460] = 4.0E0*I_TWOBODYOVERLAP_I6y_F3x_aa-2.0E0*4*I_TWOBODYOVERLAP_G4y_F3x_a-2.0E0*5*I_TWOBODYOVERLAP_G4y_F3x_a+4*3*I_TWOBODYOVERLAP_D2y_F3x;
  abcd[461] = 4.0E0*I_TWOBODYOVERLAP_I5yz_F3x_aa-2.0E0*3*I_TWOBODYOVERLAP_G3yz_F3x_a-2.0E0*4*I_TWOBODYOVERLAP_G3yz_F3x_a+3*2*I_TWOBODYOVERLAP_Dyz_F3x;
  abcd[462] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_F3x_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F3x_a+2*1*I_TWOBODYOVERLAP_D2z_F3x;
  abcd[463] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F3x_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_F3x_a;
  abcd[464] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_F3x_a;
  abcd[465] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F2xy_a;
  abcd[466] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F2xy_a-2.0E0*2*I_TWOBODYOVERLAP_G3xy_F2xy_a;
  abcd[467] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F2xy_a;
  abcd[468] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_F2xy_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F2xy_a+2*1*I_TWOBODYOVERLAP_D2x_F2xy;
  abcd[469] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_F2xy_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2xy_a;
  abcd[470] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_F2xy_a;
  abcd[471] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_F2xy_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3y_F2xy_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_F2xy_a+3*2*I_TWOBODYOVERLAP_Dxy_F2xy;
  abcd[472] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2xy_a-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_F2xy_a+2*1*I_TWOBODYOVERLAP_Dxz_F2xy;
  abcd[473] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_F2xy_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2xy_a;
  abcd[474] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F2xy_a;
  abcd[475] = 4.0E0*I_TWOBODYOVERLAP_I6y_F2xy_aa-2.0E0*4*I_TWOBODYOVERLAP_G4y_F2xy_a-2.0E0*5*I_TWOBODYOVERLAP_G4y_F2xy_a+4*3*I_TWOBODYOVERLAP_D2y_F2xy;
  abcd[476] = 4.0E0*I_TWOBODYOVERLAP_I5yz_F2xy_aa-2.0E0*3*I_TWOBODYOVERLAP_G3yz_F2xy_a-2.0E0*4*I_TWOBODYOVERLAP_G3yz_F2xy_a+3*2*I_TWOBODYOVERLAP_Dyz_F2xy;
  abcd[477] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_F2xy_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F2xy_a+2*1*I_TWOBODYOVERLAP_D2z_F2xy;
  abcd[478] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F2xy_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_F2xy_a;
  abcd[479] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_F2xy_a;
  abcd[480] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F2xz_a;
  abcd[481] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F2xz_a-2.0E0*2*I_TWOBODYOVERLAP_G3xy_F2xz_a;
  abcd[482] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F2xz_a;
  abcd[483] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_F2xz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F2xz_a+2*1*I_TWOBODYOVERLAP_D2x_F2xz;
  abcd[484] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_F2xz_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2xz_a;
  abcd[485] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_F2xz_a;
  abcd[486] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_F2xz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3y_F2xz_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_F2xz_a+3*2*I_TWOBODYOVERLAP_Dxy_F2xz;
  abcd[487] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2xz_a-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_F2xz_a+2*1*I_TWOBODYOVERLAP_Dxz_F2xz;
  abcd[488] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_F2xz_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2xz_a;
  abcd[489] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F2xz_a;
  abcd[490] = 4.0E0*I_TWOBODYOVERLAP_I6y_F2xz_aa-2.0E0*4*I_TWOBODYOVERLAP_G4y_F2xz_a-2.0E0*5*I_TWOBODYOVERLAP_G4y_F2xz_a+4*3*I_TWOBODYOVERLAP_D2y_F2xz;
  abcd[491] = 4.0E0*I_TWOBODYOVERLAP_I5yz_F2xz_aa-2.0E0*3*I_TWOBODYOVERLAP_G3yz_F2xz_a-2.0E0*4*I_TWOBODYOVERLAP_G3yz_F2xz_a+3*2*I_TWOBODYOVERLAP_Dyz_F2xz;
  abcd[492] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_F2xz_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F2xz_a+2*1*I_TWOBODYOVERLAP_D2z_F2xz;
  abcd[493] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F2xz_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_F2xz_a;
  abcd[494] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_F2xz_a;
  abcd[495] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Fx2y_a;
  abcd[496] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Fx2y_a-2.0E0*2*I_TWOBODYOVERLAP_G3xy_Fx2y_a;
  abcd[497] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Fx2y_a;
  abcd[498] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_Fx2y_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Fx2y_a+2*1*I_TWOBODYOVERLAP_D2x_Fx2y;
  abcd[499] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_Fx2y_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fx2y_a;
  abcd[500] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Fx2y_a;
  abcd[501] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_Fx2y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3y_Fx2y_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_Fx2y_a+3*2*I_TWOBODYOVERLAP_Dxy_Fx2y;
  abcd[502] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fx2y_a-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Fx2y_a+2*1*I_TWOBODYOVERLAP_Dxz_Fx2y;
  abcd[503] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_Fx2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fx2y_a;
  abcd[504] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Fx2y_a;
  abcd[505] = 4.0E0*I_TWOBODYOVERLAP_I6y_Fx2y_aa-2.0E0*4*I_TWOBODYOVERLAP_G4y_Fx2y_a-2.0E0*5*I_TWOBODYOVERLAP_G4y_Fx2y_a+4*3*I_TWOBODYOVERLAP_D2y_Fx2y;
  abcd[506] = 4.0E0*I_TWOBODYOVERLAP_I5yz_Fx2y_aa-2.0E0*3*I_TWOBODYOVERLAP_G3yz_Fx2y_a-2.0E0*4*I_TWOBODYOVERLAP_G3yz_Fx2y_a+3*2*I_TWOBODYOVERLAP_Dyz_Fx2y;
  abcd[507] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_Fx2y_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Fx2y_a+2*1*I_TWOBODYOVERLAP_D2z_Fx2y;
  abcd[508] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Fx2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_Fx2y_a;
  abcd[509] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_Fx2y_a;
  abcd[510] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Fxyz_a;
  abcd[511] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Fxyz_a-2.0E0*2*I_TWOBODYOVERLAP_G3xy_Fxyz_a;
  abcd[512] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Fxyz_a;
  abcd[513] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_Fxyz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Fxyz_a+2*1*I_TWOBODYOVERLAP_D2x_Fxyz;
  abcd[514] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_Fxyz_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fxyz_a;
  abcd[515] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Fxyz_a;
  abcd[516] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_Fxyz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3y_Fxyz_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_Fxyz_a+3*2*I_TWOBODYOVERLAP_Dxy_Fxyz;
  abcd[517] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fxyz_a-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Fxyz_a+2*1*I_TWOBODYOVERLAP_Dxz_Fxyz;
  abcd[518] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_Fxyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fxyz_a;
  abcd[519] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Fxyz_a;
  abcd[520] = 4.0E0*I_TWOBODYOVERLAP_I6y_Fxyz_aa-2.0E0*4*I_TWOBODYOVERLAP_G4y_Fxyz_a-2.0E0*5*I_TWOBODYOVERLAP_G4y_Fxyz_a+4*3*I_TWOBODYOVERLAP_D2y_Fxyz;
  abcd[521] = 4.0E0*I_TWOBODYOVERLAP_I5yz_Fxyz_aa-2.0E0*3*I_TWOBODYOVERLAP_G3yz_Fxyz_a-2.0E0*4*I_TWOBODYOVERLAP_G3yz_Fxyz_a+3*2*I_TWOBODYOVERLAP_Dyz_Fxyz;
  abcd[522] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_Fxyz_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Fxyz_a+2*1*I_TWOBODYOVERLAP_D2z_Fxyz;
  abcd[523] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Fxyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_Fxyz_a;
  abcd[524] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_Fxyz_a;
  abcd[525] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Fx2z_a;
  abcd[526] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Fx2z_a-2.0E0*2*I_TWOBODYOVERLAP_G3xy_Fx2z_a;
  abcd[527] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Fx2z_a;
  abcd[528] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_Fx2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Fx2z_a+2*1*I_TWOBODYOVERLAP_D2x_Fx2z;
  abcd[529] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_Fx2z_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fx2z_a;
  abcd[530] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Fx2z_a;
  abcd[531] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_Fx2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3y_Fx2z_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_Fx2z_a+3*2*I_TWOBODYOVERLAP_Dxy_Fx2z;
  abcd[532] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fx2z_a-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Fx2z_a+2*1*I_TWOBODYOVERLAP_Dxz_Fx2z;
  abcd[533] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_Fx2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fx2z_a;
  abcd[534] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Fx2z_a;
  abcd[535] = 4.0E0*I_TWOBODYOVERLAP_I6y_Fx2z_aa-2.0E0*4*I_TWOBODYOVERLAP_G4y_Fx2z_a-2.0E0*5*I_TWOBODYOVERLAP_G4y_Fx2z_a+4*3*I_TWOBODYOVERLAP_D2y_Fx2z;
  abcd[536] = 4.0E0*I_TWOBODYOVERLAP_I5yz_Fx2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G3yz_Fx2z_a-2.0E0*4*I_TWOBODYOVERLAP_G3yz_Fx2z_a+3*2*I_TWOBODYOVERLAP_Dyz_Fx2z;
  abcd[537] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_Fx2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Fx2z_a+2*1*I_TWOBODYOVERLAP_D2z_Fx2z;
  abcd[538] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Fx2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_Fx2z_a;
  abcd[539] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_Fx2z_a;
  abcd[540] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F3y_a;
  abcd[541] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F3y_a-2.0E0*2*I_TWOBODYOVERLAP_G3xy_F3y_a;
  abcd[542] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F3y_a;
  abcd[543] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_F3y_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F3y_a+2*1*I_TWOBODYOVERLAP_D2x_F3y;
  abcd[544] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_F3y_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3y_a;
  abcd[545] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_F3y_a;
  abcd[546] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_F3y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3y_F3y_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_F3y_a+3*2*I_TWOBODYOVERLAP_Dxy_F3y;
  abcd[547] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3y_a-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_F3y_a+2*1*I_TWOBODYOVERLAP_Dxz_F3y;
  abcd[548] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_F3y_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3y_a;
  abcd[549] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F3y_a;
  abcd[550] = 4.0E0*I_TWOBODYOVERLAP_I6y_F3y_aa-2.0E0*4*I_TWOBODYOVERLAP_G4y_F3y_a-2.0E0*5*I_TWOBODYOVERLAP_G4y_F3y_a+4*3*I_TWOBODYOVERLAP_D2y_F3y;
  abcd[551] = 4.0E0*I_TWOBODYOVERLAP_I5yz_F3y_aa-2.0E0*3*I_TWOBODYOVERLAP_G3yz_F3y_a-2.0E0*4*I_TWOBODYOVERLAP_G3yz_F3y_a+3*2*I_TWOBODYOVERLAP_Dyz_F3y;
  abcd[552] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_F3y_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F3y_a+2*1*I_TWOBODYOVERLAP_D2z_F3y;
  abcd[553] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F3y_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_F3y_a;
  abcd[554] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_F3y_a;
  abcd[555] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F2yz_a;
  abcd[556] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F2yz_a-2.0E0*2*I_TWOBODYOVERLAP_G3xy_F2yz_a;
  abcd[557] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F2yz_a;
  abcd[558] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_F2yz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F2yz_a+2*1*I_TWOBODYOVERLAP_D2x_F2yz;
  abcd[559] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_F2yz_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2yz_a;
  abcd[560] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_F2yz_a;
  abcd[561] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_F2yz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3y_F2yz_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_F2yz_a+3*2*I_TWOBODYOVERLAP_Dxy_F2yz;
  abcd[562] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2yz_a-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_F2yz_a+2*1*I_TWOBODYOVERLAP_Dxz_F2yz;
  abcd[563] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_F2yz_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2yz_a;
  abcd[564] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F2yz_a;
  abcd[565] = 4.0E0*I_TWOBODYOVERLAP_I6y_F2yz_aa-2.0E0*4*I_TWOBODYOVERLAP_G4y_F2yz_a-2.0E0*5*I_TWOBODYOVERLAP_G4y_F2yz_a+4*3*I_TWOBODYOVERLAP_D2y_F2yz;
  abcd[566] = 4.0E0*I_TWOBODYOVERLAP_I5yz_F2yz_aa-2.0E0*3*I_TWOBODYOVERLAP_G3yz_F2yz_a-2.0E0*4*I_TWOBODYOVERLAP_G3yz_F2yz_a+3*2*I_TWOBODYOVERLAP_Dyz_F2yz;
  abcd[567] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_F2yz_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F2yz_a+2*1*I_TWOBODYOVERLAP_D2z_F2yz;
  abcd[568] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F2yz_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_F2yz_a;
  abcd[569] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_F2yz_a;
  abcd[570] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Fy2z_a;
  abcd[571] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Fy2z_a-2.0E0*2*I_TWOBODYOVERLAP_G3xy_Fy2z_a;
  abcd[572] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Fy2z_a;
  abcd[573] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_Fy2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_Fy2z_a+2*1*I_TWOBODYOVERLAP_D2x_Fy2z;
  abcd[574] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_Fy2z_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fy2z_a;
  abcd[575] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Fy2z_a;
  abcd[576] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_Fy2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3y_Fy2z_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_Fy2z_a+3*2*I_TWOBODYOVERLAP_Dxy_Fy2z;
  abcd[577] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fy2z_a-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Fy2z_a+2*1*I_TWOBODYOVERLAP_Dxz_Fy2z;
  abcd[578] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_Fy2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fy2z_a;
  abcd[579] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Fy2z_a;
  abcd[580] = 4.0E0*I_TWOBODYOVERLAP_I6y_Fy2z_aa-2.0E0*4*I_TWOBODYOVERLAP_G4y_Fy2z_a-2.0E0*5*I_TWOBODYOVERLAP_G4y_Fy2z_a+4*3*I_TWOBODYOVERLAP_D2y_Fy2z;
  abcd[581] = 4.0E0*I_TWOBODYOVERLAP_I5yz_Fy2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G3yz_Fy2z_a-2.0E0*4*I_TWOBODYOVERLAP_G3yz_Fy2z_a+3*2*I_TWOBODYOVERLAP_Dyz_Fy2z;
  abcd[582] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_Fy2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Fy2z_a+2*1*I_TWOBODYOVERLAP_D2z_Fy2z;
  abcd[583] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_Fy2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_Fy2z_a;
  abcd[584] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_Fy2z_a;
  abcd[585] = 4.0E0*I_TWOBODYOVERLAP_I4x2y_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F3z_a;
  abcd[586] = 4.0E0*I_TWOBODYOVERLAP_I3x3y_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F3z_a-2.0E0*2*I_TWOBODYOVERLAP_G3xy_F3z_a;
  abcd[587] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F3z_a;
  abcd[588] = 4.0E0*I_TWOBODYOVERLAP_I2x4y_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2y_F3z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2y_F3z_a+2*1*I_TWOBODYOVERLAP_D2x_F3z;
  abcd[589] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_F3z_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3z_a;
  abcd[590] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_F3z_a;
  abcd[591] = 4.0E0*I_TWOBODYOVERLAP_Ix5y_F3z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3y_F3z_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3y_F3z_a+3*2*I_TWOBODYOVERLAP_Dxy_F3z;
  abcd[592] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3z_a-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_F3z_a+2*1*I_TWOBODYOVERLAP_Dxz_F3z;
  abcd[593] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gxy2z_F3z_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3z_a;
  abcd[594] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F3z_a;
  abcd[595] = 4.0E0*I_TWOBODYOVERLAP_I6y_F3z_aa-2.0E0*4*I_TWOBODYOVERLAP_G4y_F3z_a-2.0E0*5*I_TWOBODYOVERLAP_G4y_F3z_a+4*3*I_TWOBODYOVERLAP_D2y_F3z;
  abcd[596] = 4.0E0*I_TWOBODYOVERLAP_I5yz_F3z_aa-2.0E0*3*I_TWOBODYOVERLAP_G3yz_F3z_a-2.0E0*4*I_TWOBODYOVERLAP_G3yz_F3z_a+3*2*I_TWOBODYOVERLAP_Dyz_F3z;
  abcd[597] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_F3z_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F3z_a+2*1*I_TWOBODYOVERLAP_D2z_F3z;
  abcd[598] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gy3z_F3z_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_F3z_a;
  abcd[599] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_F_day_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_F_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_F
   ************************************************************/
  abcd[600] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F3x_aa;
  abcd[601] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F3x_a;
  abcd[602] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F3x_a;
  abcd[603] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3x_a;
  abcd[604] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_F3x_a-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_F3x_a+1*I_TWOBODYOVERLAP_D2x_F3x;
  abcd[605] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3x_a;
  abcd[606] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F3x_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_F3x_a;
  abcd[607] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F3x_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3x_a+2*1*I_TWOBODYOVERLAP_Dxy_F3x;
  abcd[608] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3x_a-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F3x_a+2*I_TWOBODYOVERLAP_Dxz_F3x;
  abcd[609] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F3x_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_F3x_a;
  abcd[610] = 4.0E0*I_TWOBODYOVERLAP_I5yz_F3x_aa-2.0E0*4*I_TWOBODYOVERLAP_G3yz_F3x_a;
  abcd[611] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_F3x_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F3x_a+3*1*I_TWOBODYOVERLAP_D2y_F3x;
  abcd[612] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_G3yz_F3x_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_F3x_a+2*2*I_TWOBODYOVERLAP_Dyz_F3x;
  abcd[613] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_F3x_aa-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F3x_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_F3x_a+3*I_TWOBODYOVERLAP_D2z_F3x;
  abcd[614] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_F3x_aa-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_F3x_a;
  abcd[615] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F2xy_aa;
  abcd[616] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F2xy_a;
  abcd[617] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F2xy_a;
  abcd[618] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2xy_a;
  abcd[619] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_F2xy_a-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_F2xy_a+1*I_TWOBODYOVERLAP_D2x_F2xy;
  abcd[620] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2xy_a;
  abcd[621] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F2xy_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_F2xy_a;
  abcd[622] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F2xy_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2xy_a+2*1*I_TWOBODYOVERLAP_Dxy_F2xy;
  abcd[623] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2xy_a-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F2xy_a+2*I_TWOBODYOVERLAP_Dxz_F2xy;
  abcd[624] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F2xy_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_F2xy_a;
  abcd[625] = 4.0E0*I_TWOBODYOVERLAP_I5yz_F2xy_aa-2.0E0*4*I_TWOBODYOVERLAP_G3yz_F2xy_a;
  abcd[626] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_F2xy_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F2xy_a+3*1*I_TWOBODYOVERLAP_D2y_F2xy;
  abcd[627] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_G3yz_F2xy_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_F2xy_a+2*2*I_TWOBODYOVERLAP_Dyz_F2xy;
  abcd[628] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_F2xy_aa-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F2xy_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_F2xy_a+3*I_TWOBODYOVERLAP_D2z_F2xy;
  abcd[629] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_F2xy_aa-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_F2xy_a;
  abcd[630] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F2xz_aa;
  abcd[631] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F2xz_a;
  abcd[632] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F2xz_a;
  abcd[633] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2xz_a;
  abcd[634] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_F2xz_a-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_F2xz_a+1*I_TWOBODYOVERLAP_D2x_F2xz;
  abcd[635] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2xz_a;
  abcd[636] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F2xz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_F2xz_a;
  abcd[637] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F2xz_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2xz_a+2*1*I_TWOBODYOVERLAP_Dxy_F2xz;
  abcd[638] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2xz_a-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F2xz_a+2*I_TWOBODYOVERLAP_Dxz_F2xz;
  abcd[639] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F2xz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_F2xz_a;
  abcd[640] = 4.0E0*I_TWOBODYOVERLAP_I5yz_F2xz_aa-2.0E0*4*I_TWOBODYOVERLAP_G3yz_F2xz_a;
  abcd[641] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_F2xz_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F2xz_a+3*1*I_TWOBODYOVERLAP_D2y_F2xz;
  abcd[642] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_G3yz_F2xz_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_F2xz_a+2*2*I_TWOBODYOVERLAP_Dyz_F2xz;
  abcd[643] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_F2xz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F2xz_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_F2xz_a+3*I_TWOBODYOVERLAP_D2z_F2xz;
  abcd[644] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_F2xz_aa-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_F2xz_a;
  abcd[645] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Fx2y_aa;
  abcd[646] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Fx2y_a;
  abcd[647] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Fx2y_a;
  abcd[648] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fx2y_a;
  abcd[649] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Fx2y_a-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Fx2y_a+1*I_TWOBODYOVERLAP_D2x_Fx2y;
  abcd[650] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fx2y_a;
  abcd[651] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Fx2y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Fx2y_a;
  abcd[652] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Fx2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fx2y_a+2*1*I_TWOBODYOVERLAP_Dxy_Fx2y;
  abcd[653] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fx2y_a-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Fx2y_a+2*I_TWOBODYOVERLAP_Dxz_Fx2y;
  abcd[654] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Fx2y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Fx2y_a;
  abcd[655] = 4.0E0*I_TWOBODYOVERLAP_I5yz_Fx2y_aa-2.0E0*4*I_TWOBODYOVERLAP_G3yz_Fx2y_a;
  abcd[656] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Fx2y_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Fx2y_a+3*1*I_TWOBODYOVERLAP_D2y_Fx2y;
  abcd[657] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G3yz_Fx2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_Fx2y_a+2*2*I_TWOBODYOVERLAP_Dyz_Fx2y;
  abcd[658] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Fx2y_aa-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Fx2y_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_Fx2y_a+3*I_TWOBODYOVERLAP_D2z_Fx2y;
  abcd[659] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_Fx2y_aa-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_Fx2y_a;
  abcd[660] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Fxyz_aa;
  abcd[661] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Fxyz_a;
  abcd[662] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Fxyz_a;
  abcd[663] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fxyz_a;
  abcd[664] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Fxyz_a-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Fxyz_a+1*I_TWOBODYOVERLAP_D2x_Fxyz;
  abcd[665] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fxyz_a;
  abcd[666] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Fxyz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Fxyz_a;
  abcd[667] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Fxyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fxyz_a+2*1*I_TWOBODYOVERLAP_Dxy_Fxyz;
  abcd[668] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fxyz_a-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Fxyz_a+2*I_TWOBODYOVERLAP_Dxz_Fxyz;
  abcd[669] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Fxyz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Fxyz_a;
  abcd[670] = 4.0E0*I_TWOBODYOVERLAP_I5yz_Fxyz_aa-2.0E0*4*I_TWOBODYOVERLAP_G3yz_Fxyz_a;
  abcd[671] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Fxyz_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Fxyz_a+3*1*I_TWOBODYOVERLAP_D2y_Fxyz;
  abcd[672] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G3yz_Fxyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_Fxyz_a+2*2*I_TWOBODYOVERLAP_Dyz_Fxyz;
  abcd[673] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Fxyz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Fxyz_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_Fxyz_a+3*I_TWOBODYOVERLAP_D2z_Fxyz;
  abcd[674] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_Fxyz_aa-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_Fxyz_a;
  abcd[675] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Fx2z_aa;
  abcd[676] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Fx2z_a;
  abcd[677] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Fx2z_a;
  abcd[678] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fx2z_a;
  abcd[679] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Fx2z_a-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Fx2z_a+1*I_TWOBODYOVERLAP_D2x_Fx2z;
  abcd[680] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fx2z_a;
  abcd[681] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Fx2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Fx2z_a;
  abcd[682] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Fx2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fx2z_a+2*1*I_TWOBODYOVERLAP_Dxy_Fx2z;
  abcd[683] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fx2z_a-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Fx2z_a+2*I_TWOBODYOVERLAP_Dxz_Fx2z;
  abcd[684] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Fx2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Fx2z_a;
  abcd[685] = 4.0E0*I_TWOBODYOVERLAP_I5yz_Fx2z_aa-2.0E0*4*I_TWOBODYOVERLAP_G3yz_Fx2z_a;
  abcd[686] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Fx2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Fx2z_a+3*1*I_TWOBODYOVERLAP_D2y_Fx2z;
  abcd[687] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G3yz_Fx2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_Fx2z_a+2*2*I_TWOBODYOVERLAP_Dyz_Fx2z;
  abcd[688] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Fx2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Fx2z_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_Fx2z_a+3*I_TWOBODYOVERLAP_D2z_Fx2z;
  abcd[689] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_Fx2z_aa-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_Fx2z_a;
  abcd[690] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F3y_aa;
  abcd[691] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F3y_a;
  abcd[692] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F3y_a;
  abcd[693] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3y_a;
  abcd[694] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_F3y_a-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_F3y_a+1*I_TWOBODYOVERLAP_D2x_F3y;
  abcd[695] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3y_a;
  abcd[696] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F3y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_F3y_a;
  abcd[697] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F3y_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3y_a+2*1*I_TWOBODYOVERLAP_Dxy_F3y;
  abcd[698] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3y_a-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F3y_a+2*I_TWOBODYOVERLAP_Dxz_F3y;
  abcd[699] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F3y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_F3y_a;
  abcd[700] = 4.0E0*I_TWOBODYOVERLAP_I5yz_F3y_aa-2.0E0*4*I_TWOBODYOVERLAP_G3yz_F3y_a;
  abcd[701] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_F3y_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F3y_a+3*1*I_TWOBODYOVERLAP_D2y_F3y;
  abcd[702] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_G3yz_F3y_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_F3y_a+2*2*I_TWOBODYOVERLAP_Dyz_F3y;
  abcd[703] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_F3y_aa-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F3y_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_F3y_a+3*I_TWOBODYOVERLAP_D2z_F3y;
  abcd[704] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_F3y_aa-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_F3y_a;
  abcd[705] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F2yz_aa;
  abcd[706] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F2yz_a;
  abcd[707] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F2yz_a;
  abcd[708] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2yz_a;
  abcd[709] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_F2yz_a-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_F2yz_a+1*I_TWOBODYOVERLAP_D2x_F2yz;
  abcd[710] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2yz_a;
  abcd[711] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F2yz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_F2yz_a;
  abcd[712] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F2yz_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2yz_a+2*1*I_TWOBODYOVERLAP_Dxy_F2yz;
  abcd[713] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2yz_a-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F2yz_a+2*I_TWOBODYOVERLAP_Dxz_F2yz;
  abcd[714] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F2yz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_F2yz_a;
  abcd[715] = 4.0E0*I_TWOBODYOVERLAP_I5yz_F2yz_aa-2.0E0*4*I_TWOBODYOVERLAP_G3yz_F2yz_a;
  abcd[716] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_F2yz_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F2yz_a+3*1*I_TWOBODYOVERLAP_D2y_F2yz;
  abcd[717] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_G3yz_F2yz_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_F2yz_a+2*2*I_TWOBODYOVERLAP_Dyz_F2yz;
  abcd[718] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_F2yz_aa-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F2yz_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_F2yz_a+3*I_TWOBODYOVERLAP_D2z_F2yz;
  abcd[719] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_F2yz_aa-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_F2yz_a;
  abcd[720] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_Fy2z_aa;
  abcd[721] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Fy2z_a;
  abcd[722] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Fy2z_a;
  abcd[723] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fy2z_a;
  abcd[724] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Fy2z_a-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_Fy2z_a+1*I_TWOBODYOVERLAP_D2x_Fy2z;
  abcd[725] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fy2z_a;
  abcd[726] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_Fy2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_Fy2z_a;
  abcd[727] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Fy2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fy2z_a+2*1*I_TWOBODYOVERLAP_Dxy_Fy2z;
  abcd[728] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fy2z_a-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_Fy2z_a+2*I_TWOBODYOVERLAP_Dxz_Fy2z;
  abcd[729] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Fy2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Fy2z_a;
  abcd[730] = 4.0E0*I_TWOBODYOVERLAP_I5yz_Fy2z_aa-2.0E0*4*I_TWOBODYOVERLAP_G3yz_Fy2z_a;
  abcd[731] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Fy2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Fy2z_a+3*1*I_TWOBODYOVERLAP_D2y_Fy2z;
  abcd[732] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G3yz_Fy2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_Fy2z_a+2*2*I_TWOBODYOVERLAP_Dyz_Fy2z;
  abcd[733] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Fy2z_aa-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Fy2z_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_Fy2z_a+3*I_TWOBODYOVERLAP_D2z_Fy2z;
  abcd[734] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_Fy2z_aa-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_Fy2z_a;
  abcd[735] = 4.0E0*I_TWOBODYOVERLAP_I4xyz_F3z_aa;
  abcd[736] = 4.0E0*I_TWOBODYOVERLAP_I3x2yz_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F3z_a;
  abcd[737] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F3z_a;
  abcd[738] = 4.0E0*I_TWOBODYOVERLAP_I2x3yz_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3z_a;
  abcd[739] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_F3z_a-2.0E0*1*I_TWOBODYOVERLAP_G2x2z_F3z_a+1*I_TWOBODYOVERLAP_D2x_F3z;
  abcd[740] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3z_a;
  abcd[741] = 4.0E0*I_TWOBODYOVERLAP_Ix4yz_F3z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx2yz_F3z_a;
  abcd[742] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F3z_a-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3z_a+2*1*I_TWOBODYOVERLAP_Dxy_F3z;
  abcd[743] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3z_a-2.0E0*1*I_TWOBODYOVERLAP_Gx3z_F3z_a+2*I_TWOBODYOVERLAP_Dxz_F3z;
  abcd[744] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F3z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_F3z_a;
  abcd[745] = 4.0E0*I_TWOBODYOVERLAP_I5yz_F3z_aa-2.0E0*4*I_TWOBODYOVERLAP_G3yz_F3z_a;
  abcd[746] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_F3z_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F3z_a+3*1*I_TWOBODYOVERLAP_D2y_F3z;
  abcd[747] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_G3yz_F3z_a-2.0E0*2*I_TWOBODYOVERLAP_Gy3z_F3z_a+2*2*I_TWOBODYOVERLAP_Dyz_F3z;
  abcd[748] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_F3z_aa-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F3z_a-2.0E0*1*I_TWOBODYOVERLAP_G4z_F3z_a+3*I_TWOBODYOVERLAP_D2z_F3z;
  abcd[749] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_F3z_aa-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_F3z_a;

  /************************************************************
   * shell quartet name: SQ_TWOBODYOVERLAP_G_F_daz_daz
   * code section is: DERIV
   * totally 0 integrals are omitted 
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_I_F_aa
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_F_a
   * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_F
   ************************************************************/
  abcd[750] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F3x_a;
  abcd[751] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F3x_a;
  abcd[752] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F3x_a-2.0E0*2*I_TWOBODYOVERLAP_G3xz_F3x_a;
  abcd[753] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_F3x_a;
  abcd[754] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_F3x_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3x_a;
  abcd[755] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_F3x_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F3x_a+2*1*I_TWOBODYOVERLAP_D2x_F3x;
  abcd[756] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F3x_a;
  abcd[757] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_F3x_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3x_a;
  abcd[758] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3x_a-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_F3x_a+2*1*I_TWOBODYOVERLAP_Dxy_F3x;
  abcd[759] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_F3x_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3z_F3x_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_F3x_a+3*2*I_TWOBODYOVERLAP_Dxz_F3x;
  abcd[760] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_F3x_a;
  abcd[761] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_F3x_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F3x_a-2.0E0*2*I_TWOBODYOVERLAP_G3yz_F3x_a;
  abcd[762] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_F3x_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_F3x_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F3x_a+2*1*I_TWOBODYOVERLAP_D2y_F3x;
  abcd[763] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_F3x_aa-2.0E0*3*I_TWOBODYOVERLAP_Gy3z_F3x_a-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_F3x_a+3*2*I_TWOBODYOVERLAP_Dyz_F3x;
  abcd[764] = 4.0E0*I_TWOBODYOVERLAP_I6z_F3x_aa-2.0E0*4*I_TWOBODYOVERLAP_G4z_F3x_a-2.0E0*5*I_TWOBODYOVERLAP_G4z_F3x_a+4*3*I_TWOBODYOVERLAP_D2z_F3x;
  abcd[765] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F2xy_a;
  abcd[766] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F2xy_a;
  abcd[767] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F2xy_a-2.0E0*2*I_TWOBODYOVERLAP_G3xz_F2xy_a;
  abcd[768] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_F2xy_a;
  abcd[769] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_F2xy_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2xy_a;
  abcd[770] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_F2xy_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F2xy_a+2*1*I_TWOBODYOVERLAP_D2x_F2xy;
  abcd[771] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F2xy_a;
  abcd[772] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_F2xy_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2xy_a;
  abcd[773] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2xy_a-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_F2xy_a+2*1*I_TWOBODYOVERLAP_Dxy_F2xy;
  abcd[774] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_F2xy_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3z_F2xy_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_F2xy_a+3*2*I_TWOBODYOVERLAP_Dxz_F2xy;
  abcd[775] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_F2xy_a;
  abcd[776] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_F2xy_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F2xy_a-2.0E0*2*I_TWOBODYOVERLAP_G3yz_F2xy_a;
  abcd[777] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_F2xy_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_F2xy_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F2xy_a+2*1*I_TWOBODYOVERLAP_D2y_F2xy;
  abcd[778] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_F2xy_aa-2.0E0*3*I_TWOBODYOVERLAP_Gy3z_F2xy_a-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_F2xy_a+3*2*I_TWOBODYOVERLAP_Dyz_F2xy;
  abcd[779] = 4.0E0*I_TWOBODYOVERLAP_I6z_F2xy_aa-2.0E0*4*I_TWOBODYOVERLAP_G4z_F2xy_a-2.0E0*5*I_TWOBODYOVERLAP_G4z_F2xy_a+4*3*I_TWOBODYOVERLAP_D2z_F2xy;
  abcd[780] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F2xz_a;
  abcd[781] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F2xz_a;
  abcd[782] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F2xz_a-2.0E0*2*I_TWOBODYOVERLAP_G3xz_F2xz_a;
  abcd[783] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_F2xz_a;
  abcd[784] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_F2xz_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2xz_a;
  abcd[785] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_F2xz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F2xz_a+2*1*I_TWOBODYOVERLAP_D2x_F2xz;
  abcd[786] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F2xz_a;
  abcd[787] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_F2xz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2xz_a;
  abcd[788] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2xz_a-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_F2xz_a+2*1*I_TWOBODYOVERLAP_Dxy_F2xz;
  abcd[789] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_F2xz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3z_F2xz_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_F2xz_a+3*2*I_TWOBODYOVERLAP_Dxz_F2xz;
  abcd[790] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_F2xz_a;
  abcd[791] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_F2xz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F2xz_a-2.0E0*2*I_TWOBODYOVERLAP_G3yz_F2xz_a;
  abcd[792] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_F2xz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_F2xz_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F2xz_a+2*1*I_TWOBODYOVERLAP_D2y_F2xz;
  abcd[793] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_F2xz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gy3z_F2xz_a-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_F2xz_a+3*2*I_TWOBODYOVERLAP_Dyz_F2xz;
  abcd[794] = 4.0E0*I_TWOBODYOVERLAP_I6z_F2xz_aa-2.0E0*4*I_TWOBODYOVERLAP_G4z_F2xz_a-2.0E0*5*I_TWOBODYOVERLAP_G4z_F2xz_a+4*3*I_TWOBODYOVERLAP_D2z_F2xz;
  abcd[795] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Fx2y_a;
  abcd[796] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Fx2y_a;
  abcd[797] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Fx2y_a-2.0E0*2*I_TWOBODYOVERLAP_G3xz_Fx2y_a;
  abcd[798] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Fx2y_a;
  abcd[799] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_Fx2y_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fx2y_a;
  abcd[800] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_Fx2y_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Fx2y_a+2*1*I_TWOBODYOVERLAP_D2x_Fx2y;
  abcd[801] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Fx2y_a;
  abcd[802] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_Fx2y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fx2y_a;
  abcd[803] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fx2y_a-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Fx2y_a+2*1*I_TWOBODYOVERLAP_Dxy_Fx2y;
  abcd[804] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_Fx2y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3z_Fx2y_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_Fx2y_a+3*2*I_TWOBODYOVERLAP_Dxz_Fx2y;
  abcd[805] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Fx2y_a;
  abcd[806] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Fx2y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Fx2y_a-2.0E0*2*I_TWOBODYOVERLAP_G3yz_Fx2y_a;
  abcd[807] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Fx2y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_Fx2y_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Fx2y_a+2*1*I_TWOBODYOVERLAP_D2y_Fx2y;
  abcd[808] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_Fx2y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gy3z_Fx2y_a-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_Fx2y_a+3*2*I_TWOBODYOVERLAP_Dyz_Fx2y;
  abcd[809] = 4.0E0*I_TWOBODYOVERLAP_I6z_Fx2y_aa-2.0E0*4*I_TWOBODYOVERLAP_G4z_Fx2y_a-2.0E0*5*I_TWOBODYOVERLAP_G4z_Fx2y_a+4*3*I_TWOBODYOVERLAP_D2z_Fx2y;
  abcd[810] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Fxyz_a;
  abcd[811] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Fxyz_a;
  abcd[812] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Fxyz_a-2.0E0*2*I_TWOBODYOVERLAP_G3xz_Fxyz_a;
  abcd[813] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Fxyz_a;
  abcd[814] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_Fxyz_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fxyz_a;
  abcd[815] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_Fxyz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Fxyz_a+2*1*I_TWOBODYOVERLAP_D2x_Fxyz;
  abcd[816] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Fxyz_a;
  abcd[817] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_Fxyz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fxyz_a;
  abcd[818] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fxyz_a-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Fxyz_a+2*1*I_TWOBODYOVERLAP_Dxy_Fxyz;
  abcd[819] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_Fxyz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3z_Fxyz_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_Fxyz_a+3*2*I_TWOBODYOVERLAP_Dxz_Fxyz;
  abcd[820] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Fxyz_a;
  abcd[821] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Fxyz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Fxyz_a-2.0E0*2*I_TWOBODYOVERLAP_G3yz_Fxyz_a;
  abcd[822] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Fxyz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_Fxyz_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Fxyz_a+2*1*I_TWOBODYOVERLAP_D2y_Fxyz;
  abcd[823] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_Fxyz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gy3z_Fxyz_a-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_Fxyz_a+3*2*I_TWOBODYOVERLAP_Dyz_Fxyz;
  abcd[824] = 4.0E0*I_TWOBODYOVERLAP_I6z_Fxyz_aa-2.0E0*4*I_TWOBODYOVERLAP_G4z_Fxyz_a-2.0E0*5*I_TWOBODYOVERLAP_G4z_Fxyz_a+4*3*I_TWOBODYOVERLAP_D2z_Fxyz;
  abcd[825] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Fx2z_a;
  abcd[826] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Fx2z_a;
  abcd[827] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Fx2z_a-2.0E0*2*I_TWOBODYOVERLAP_G3xz_Fx2z_a;
  abcd[828] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Fx2z_a;
  abcd[829] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_Fx2z_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fx2z_a;
  abcd[830] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_Fx2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Fx2z_a+2*1*I_TWOBODYOVERLAP_D2x_Fx2z;
  abcd[831] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Fx2z_a;
  abcd[832] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_Fx2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fx2z_a;
  abcd[833] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fx2z_a-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Fx2z_a+2*1*I_TWOBODYOVERLAP_Dxy_Fx2z;
  abcd[834] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_Fx2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3z_Fx2z_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_Fx2z_a+3*2*I_TWOBODYOVERLAP_Dxz_Fx2z;
  abcd[835] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Fx2z_a;
  abcd[836] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Fx2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Fx2z_a-2.0E0*2*I_TWOBODYOVERLAP_G3yz_Fx2z_a;
  abcd[837] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Fx2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_Fx2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Fx2z_a+2*1*I_TWOBODYOVERLAP_D2y_Fx2z;
  abcd[838] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_Fx2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gy3z_Fx2z_a-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_Fx2z_a+3*2*I_TWOBODYOVERLAP_Dyz_Fx2z;
  abcd[839] = 4.0E0*I_TWOBODYOVERLAP_I6z_Fx2z_aa-2.0E0*4*I_TWOBODYOVERLAP_G4z_Fx2z_a-2.0E0*5*I_TWOBODYOVERLAP_G4z_Fx2z_a+4*3*I_TWOBODYOVERLAP_D2z_Fx2z;
  abcd[840] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F3y_a;
  abcd[841] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F3y_a;
  abcd[842] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F3y_a-2.0E0*2*I_TWOBODYOVERLAP_G3xz_F3y_a;
  abcd[843] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_F3y_a;
  abcd[844] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_F3y_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3y_a;
  abcd[845] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_F3y_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F3y_a+2*1*I_TWOBODYOVERLAP_D2x_F3y;
  abcd[846] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F3y_a;
  abcd[847] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_F3y_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3y_a;
  abcd[848] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3y_a-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_F3y_a+2*1*I_TWOBODYOVERLAP_Dxy_F3y;
  abcd[849] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_F3y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3z_F3y_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_F3y_a+3*2*I_TWOBODYOVERLAP_Dxz_F3y;
  abcd[850] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_F3y_a;
  abcd[851] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_F3y_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F3y_a-2.0E0*2*I_TWOBODYOVERLAP_G3yz_F3y_a;
  abcd[852] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_F3y_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_F3y_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F3y_a+2*1*I_TWOBODYOVERLAP_D2y_F3y;
  abcd[853] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_F3y_aa-2.0E0*3*I_TWOBODYOVERLAP_Gy3z_F3y_a-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_F3y_a+3*2*I_TWOBODYOVERLAP_Dyz_F3y;
  abcd[854] = 4.0E0*I_TWOBODYOVERLAP_I6z_F3y_aa-2.0E0*4*I_TWOBODYOVERLAP_G4z_F3y_a-2.0E0*5*I_TWOBODYOVERLAP_G4z_F3y_a+4*3*I_TWOBODYOVERLAP_D2z_F3y;
  abcd[855] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F2yz_a;
  abcd[856] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F2yz_a;
  abcd[857] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F2yz_a-2.0E0*2*I_TWOBODYOVERLAP_G3xz_F2yz_a;
  abcd[858] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_F2yz_a;
  abcd[859] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_F2yz_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F2yz_a;
  abcd[860] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_F2yz_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F2yz_a+2*1*I_TWOBODYOVERLAP_D2x_F2yz;
  abcd[861] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F2yz_a;
  abcd[862] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_F2yz_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F2yz_a;
  abcd[863] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F2yz_a-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_F2yz_a+2*1*I_TWOBODYOVERLAP_Dxy_F2yz;
  abcd[864] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_F2yz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3z_F2yz_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_F2yz_a+3*2*I_TWOBODYOVERLAP_Dxz_F2yz;
  abcd[865] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_F2yz_a;
  abcd[866] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_F2yz_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F2yz_a-2.0E0*2*I_TWOBODYOVERLAP_G3yz_F2yz_a;
  abcd[867] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_F2yz_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_F2yz_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F2yz_a+2*1*I_TWOBODYOVERLAP_D2y_F2yz;
  abcd[868] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_F2yz_aa-2.0E0*3*I_TWOBODYOVERLAP_Gy3z_F2yz_a-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_F2yz_a+3*2*I_TWOBODYOVERLAP_Dyz_F2yz;
  abcd[869] = 4.0E0*I_TWOBODYOVERLAP_I6z_F2yz_aa-2.0E0*4*I_TWOBODYOVERLAP_G4z_F2yz_a-2.0E0*5*I_TWOBODYOVERLAP_G4z_F2yz_a+4*3*I_TWOBODYOVERLAP_D2z_F2yz;
  abcd[870] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_Fy2z_a;
  abcd[871] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_Fy2z_a;
  abcd[872] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_Fy2z_a-2.0E0*2*I_TWOBODYOVERLAP_G3xz_Fy2z_a;
  abcd[873] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_Fy2z_a;
  abcd[874] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_Fy2z_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_Fy2z_a;
  abcd[875] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_Fy2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_Fy2z_a+2*1*I_TWOBODYOVERLAP_D2x_Fy2z;
  abcd[876] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_Fy2z_a;
  abcd[877] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_Fy2z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_Fy2z_a;
  abcd[878] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_Fy2z_a-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_Fy2z_a+2*1*I_TWOBODYOVERLAP_Dxy_Fy2z;
  abcd[879] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_Fy2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3z_Fy2z_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_Fy2z_a+3*2*I_TWOBODYOVERLAP_Dxz_Fy2z;
  abcd[880] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_Fy2z_a;
  abcd[881] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_Fy2z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_Fy2z_a-2.0E0*2*I_TWOBODYOVERLAP_G3yz_Fy2z_a;
  abcd[882] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_Fy2z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_Fy2z_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_Fy2z_a+2*1*I_TWOBODYOVERLAP_D2y_Fy2z;
  abcd[883] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_Fy2z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gy3z_Fy2z_a-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_Fy2z_a+3*2*I_TWOBODYOVERLAP_Dyz_Fy2z;
  abcd[884] = 4.0E0*I_TWOBODYOVERLAP_I6z_Fy2z_aa-2.0E0*4*I_TWOBODYOVERLAP_G4z_Fy2z_a-2.0E0*5*I_TWOBODYOVERLAP_G4z_Fy2z_a+4*3*I_TWOBODYOVERLAP_D2z_Fy2z;
  abcd[885] = 4.0E0*I_TWOBODYOVERLAP_I4x2z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4x_F3z_a;
  abcd[886] = 4.0E0*I_TWOBODYOVERLAP_I3xy2z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xy_F3z_a;
  abcd[887] = 4.0E0*I_TWOBODYOVERLAP_I3x3z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3xz_F3z_a-2.0E0*2*I_TWOBODYOVERLAP_G3xz_F3z_a;
  abcd[888] = 4.0E0*I_TWOBODYOVERLAP_I2x2y2z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2x2y_F3z_a;
  abcd[889] = 4.0E0*I_TWOBODYOVERLAP_I2xy3z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G2xyz_F3z_a-2.0E0*2*I_TWOBODYOVERLAP_G2xyz_F3z_a;
  abcd[890] = 4.0E0*I_TWOBODYOVERLAP_I2x4z_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2x2z_F3z_a-2.0E0*3*I_TWOBODYOVERLAP_G2x2z_F3z_a+2*1*I_TWOBODYOVERLAP_D2x_F3z;
  abcd[891] = 4.0E0*I_TWOBODYOVERLAP_Ix3y2z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx3y_F3z_a;
  abcd[892] = 4.0E0*I_TWOBODYOVERLAP_Ix2y3z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_Gx2yz_F3z_a-2.0E0*2*I_TWOBODYOVERLAP_Gx2yz_F3z_a;
  abcd[893] = 4.0E0*I_TWOBODYOVERLAP_Ixy4z_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_Gxy2z_F3z_a-2.0E0*3*I_TWOBODYOVERLAP_Gxy2z_F3z_a+2*1*I_TWOBODYOVERLAP_Dxy_F3z;
  abcd[894] = 4.0E0*I_TWOBODYOVERLAP_Ix5z_F3z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gx3z_F3z_a-2.0E0*4*I_TWOBODYOVERLAP_Gx3z_F3z_a+3*2*I_TWOBODYOVERLAP_Dxz_F3z;
  abcd[895] = 4.0E0*I_TWOBODYOVERLAP_I4y2z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G4y_F3z_a;
  abcd[896] = 4.0E0*I_TWOBODYOVERLAP_I3y3z_F3z_aa-2.0E0*1*I_TWOBODYOVERLAP_G3yz_F3z_a-2.0E0*2*I_TWOBODYOVERLAP_G3yz_F3z_a;
  abcd[897] = 4.0E0*I_TWOBODYOVERLAP_I2y4z_F3z_aa-2.0E0*2*I_TWOBODYOVERLAP_G2y2z_F3z_a-2.0E0*3*I_TWOBODYOVERLAP_G2y2z_F3z_a+2*1*I_TWOBODYOVERLAP_D2y_F3z;
  abcd[898] = 4.0E0*I_TWOBODYOVERLAP_Iy5z_F3z_aa-2.0E0*3*I_TWOBODYOVERLAP_Gy3z_F3z_a-2.0E0*4*I_TWOBODYOVERLAP_Gy3z_F3z_a+3*2*I_TWOBODYOVERLAP_Dyz_F3z;
  abcd[899] = 4.0E0*I_TWOBODYOVERLAP_I6z_F3z_aa-2.0E0*4*I_TWOBODYOVERLAP_G4z_F3z_a-2.0E0*5*I_TWOBODYOVERLAP_G4z_F3z_a+4*3*I_TWOBODYOVERLAP_D2z_F3z;
}
