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

void hgp_os_kinetic_h_p(const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd)
{
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
     * totally 2 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_F3x_S_vrr = PAX*I_TWOBODYOVERLAP_D2x_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_TWOBODYOVERLAP_F2xy_S_vrr = PAY*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_F2xz_S_vrr = PAZ*I_TWOBODYOVERLAP_D2x_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2y_S_vrr = PAX*I_TWOBODYOVERLAP_D2y_S_vrr;
    Double I_TWOBODYOVERLAP_Fx2z_S_vrr = PAX*I_TWOBODYOVERLAP_D2z_S_vrr;
    Double I_TWOBODYOVERLAP_F3y_S_vrr = PAY*I_TWOBODYOVERLAP_D2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_TWOBODYOVERLAP_F2yz_S_vrr = PAZ*I_TWOBODYOVERLAP_D2y_S_vrr;
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
     * shell quartet name: SQ_TWOBODYOVERLAP_H_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_G_S
     ************************************************************/
    Double I_TWOBODYOVERLAP_H5x_Px_vrr = PBX*I_TWOBODYOVERLAP_H5x_S_vrr+5*oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Px_vrr = PBX*I_TWOBODYOVERLAP_H4xy_S_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Px_vrr = PBX*I_TWOBODYOVERLAP_H4xz_S_vrr+4*oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Px_vrr = PBX*I_TWOBODYOVERLAP_H3x2y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Px_vrr = PBX*I_TWOBODYOVERLAP_H3xyz_S_vrr+3*oned2z*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Px_vrr = PBX*I_TWOBODYOVERLAP_H3x2z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Px_vrr = PBX*I_TWOBODYOVERLAP_H2x3y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Px_vrr = PBX*I_TWOBODYOVERLAP_H2x2yz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Px_vrr = PBX*I_TWOBODYOVERLAP_H2xy2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Px_vrr = PBX*I_TWOBODYOVERLAP_H2x3z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Px_vrr = PBX*I_TWOBODYOVERLAP_Hx4y_S_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Px_vrr = PBX*I_TWOBODYOVERLAP_Hx3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Px_vrr = PBX*I_TWOBODYOVERLAP_Hx2y2z_S_vrr+oned2z*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Px_vrr = PBX*I_TWOBODYOVERLAP_Hxy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Px_vrr = PBX*I_TWOBODYOVERLAP_Hx4z_S_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5y_Px_vrr = PBX*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Px_vrr = PBX*I_TWOBODYOVERLAP_H4yz_S_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Px_vrr = PBX*I_TWOBODYOVERLAP_H3y2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Px_vrr = PBX*I_TWOBODYOVERLAP_H2y3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Px_vrr = PBX*I_TWOBODYOVERLAP_Hy4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5z_Px_vrr = PBX*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_H5x_Py_vrr = PBY*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Py_vrr = PBY*I_TWOBODYOVERLAP_H4xy_S_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Py_vrr = PBY*I_TWOBODYOVERLAP_H4xz_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Py_vrr = PBY*I_TWOBODYOVERLAP_H3x2y_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Py_vrr = PBY*I_TWOBODYOVERLAP_H3xyz_S_vrr+oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Py_vrr = PBY*I_TWOBODYOVERLAP_H3x2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Py_vrr = PBY*I_TWOBODYOVERLAP_H2x3y_S_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Py_vrr = PBY*I_TWOBODYOVERLAP_H2x2yz_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Py_vrr = PBY*I_TWOBODYOVERLAP_H2xy2z_S_vrr+oned2z*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Py_vrr = PBY*I_TWOBODYOVERLAP_H2x3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Py_vrr = PBY*I_TWOBODYOVERLAP_Hx4y_S_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Py_vrr = PBY*I_TWOBODYOVERLAP_Hx3yz_S_vrr+3*oned2z*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Py_vrr = PBY*I_TWOBODYOVERLAP_Hx2y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Py_vrr = PBY*I_TWOBODYOVERLAP_Hxy3z_S_vrr+oned2z*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Py_vrr = PBY*I_TWOBODYOVERLAP_Hx4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5y_Py_vrr = PBY*I_TWOBODYOVERLAP_H5y_S_vrr+5*oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Py_vrr = PBY*I_TWOBODYOVERLAP_H4yz_S_vrr+4*oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Py_vrr = PBY*I_TWOBODYOVERLAP_H3y2z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Py_vrr = PBY*I_TWOBODYOVERLAP_H2y3z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Py_vrr = PBY*I_TWOBODYOVERLAP_Hy4z_S_vrr+oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;
    Double I_TWOBODYOVERLAP_H5z_Py_vrr = PBY*I_TWOBODYOVERLAP_H5z_S_vrr;
    Double I_TWOBODYOVERLAP_H5x_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H5x_S_vrr;
    Double I_TWOBODYOVERLAP_H4xy_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H4xy_S_vrr;
    Double I_TWOBODYOVERLAP_H4xz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H4xz_S_vrr+oned2z*I_TWOBODYOVERLAP_G4x_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H3xyz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H3xyz_S_vrr+oned2z*I_TWOBODYOVERLAP_G3xy_S_vrr;
    Double I_TWOBODYOVERLAP_H3x2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H3x2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G3xz_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3y_S_vrr;
    Double I_TWOBODYOVERLAP_H2x2yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H2x2yz_S_vrr+oned2z*I_TWOBODYOVERLAP_G2x2y_S_vrr;
    Double I_TWOBODYOVERLAP_H2xy2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H2xy2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G2xyz_S_vrr;
    Double I_TWOBODYOVERLAP_H2x3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H2x3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_G2x2z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx3yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Hx3yz_S_vrr+oned2z*I_TWOBODYOVERLAP_Gx3y_S_vrr;
    Double I_TWOBODYOVERLAP_Hx2y2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Hx2y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_Gx2yz_S_vrr;
    Double I_TWOBODYOVERLAP_Hxy3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Hxy3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_Gxy2z_S_vrr;
    Double I_TWOBODYOVERLAP_Hx4z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Hx4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_Gx3z_S_vrr;
    Double I_TWOBODYOVERLAP_H5y_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H5y_S_vrr;
    Double I_TWOBODYOVERLAP_H4yz_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H4yz_S_vrr+oned2z*I_TWOBODYOVERLAP_G4y_S_vrr;
    Double I_TWOBODYOVERLAP_H3y2z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H3y2z_S_vrr+2*oned2z*I_TWOBODYOVERLAP_G3yz_S_vrr;
    Double I_TWOBODYOVERLAP_H2y3z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H2y3z_S_vrr+3*oned2z*I_TWOBODYOVERLAP_G2y2z_S_vrr;
    Double I_TWOBODYOVERLAP_Hy4z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_Hy4z_S_vrr+4*oned2z*I_TWOBODYOVERLAP_Gy3z_S_vrr;
    Double I_TWOBODYOVERLAP_H5z_Pz_vrr = PBZ*I_TWOBODYOVERLAP_H5z_S_vrr+5*oned2z*I_TWOBODYOVERLAP_G4z_S_vrr;

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
     * totally 3 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_KINETIC_S_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_D_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_S_S
     ************************************************************/
    Double I_KINETIC_D2x_S_vrr = PAX*I_KINETIC_Px_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_D2x_S_vrr-bdz*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_KINETIC_D2y_S_vrr = PAY*I_KINETIC_Py_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_D2y_S_vrr-bdz*I_TWOBODYOVERLAP_S_S_vrr;
    Double I_KINETIC_D2z_S_vrr = PAZ*I_KINETIC_Pz_S_vrr+oned2z*I_KINETIC_S_S_vrr+twoxi*I_TWOBODYOVERLAP_D2z_S_vrr-bdz*I_TWOBODYOVERLAP_S_S_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_F_S
     * expanding position: BRA1
     * code section is: VRR
     * totally 2 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_D_S
     * RHS shell quartet name: SQ_KINETIC_P_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_F_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_P_S
     ************************************************************/
    Double I_KINETIC_F3x_S_vrr = PAX*I_KINETIC_D2x_S_vrr+2*oned2z*I_KINETIC_Px_S_vrr+twoxi*I_TWOBODYOVERLAP_F3x_S_vrr-2*bdz*I_TWOBODYOVERLAP_Px_S_vrr;
    Double I_KINETIC_F2xy_S_vrr = PAY*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xy_S_vrr;
    Double I_KINETIC_F2xz_S_vrr = PAZ*I_KINETIC_D2x_S_vrr+twoxi*I_TWOBODYOVERLAP_F2xz_S_vrr;
    Double I_KINETIC_Fx2y_S_vrr = PAX*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2y_S_vrr;
    Double I_KINETIC_Fx2z_S_vrr = PAX*I_KINETIC_D2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Fx2z_S_vrr;
    Double I_KINETIC_F3y_S_vrr = PAY*I_KINETIC_D2y_S_vrr+2*oned2z*I_KINETIC_Py_S_vrr+twoxi*I_TWOBODYOVERLAP_F3y_S_vrr-2*bdz*I_TWOBODYOVERLAP_Py_S_vrr;
    Double I_KINETIC_F2yz_S_vrr = PAZ*I_KINETIC_D2y_S_vrr+twoxi*I_TWOBODYOVERLAP_F2yz_S_vrr;
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
     * shell quartet name: SQ_KINETIC_H_P
     * expanding position: BRA2
     * code section is: VRR
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_KINETIC_H_S
     * RHS shell quartet name: SQ_KINETIC_G_S
     * RHS shell quartet name: SQ_TWOBODYOVERLAP_H_P
     ************************************************************/
    Double I_KINETIC_H5x_Px_vrr = PBX*I_KINETIC_H5x_S_vrr+5*oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Px_vrr;
    Double I_KINETIC_H4xy_Px_vrr = PBX*I_KINETIC_H4xy_S_vrr+4*oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Px_vrr;
    Double I_KINETIC_H4xz_Px_vrr = PBX*I_KINETIC_H4xz_S_vrr+4*oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Px_vrr;
    Double I_KINETIC_H3x2y_Px_vrr = PBX*I_KINETIC_H3x2y_S_vrr+3*oned2z*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Px_vrr;
    Double I_KINETIC_H3xyz_Px_vrr = PBX*I_KINETIC_H3xyz_S_vrr+3*oned2z*I_KINETIC_G2xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Px_vrr;
    Double I_KINETIC_H3x2z_Px_vrr = PBX*I_KINETIC_H3x2z_S_vrr+3*oned2z*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Px_vrr;
    Double I_KINETIC_H2x3y_Px_vrr = PBX*I_KINETIC_H2x3y_S_vrr+2*oned2z*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Px_vrr;
    Double I_KINETIC_H2x2yz_Px_vrr = PBX*I_KINETIC_H2x2yz_S_vrr+2*oned2z*I_KINETIC_Gx2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Px_vrr;
    Double I_KINETIC_H2xy2z_Px_vrr = PBX*I_KINETIC_H2xy2z_S_vrr+2*oned2z*I_KINETIC_Gxy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Px_vrr;
    Double I_KINETIC_H2x3z_Px_vrr = PBX*I_KINETIC_H2x3z_S_vrr+2*oned2z*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Px_vrr;
    Double I_KINETIC_Hx4y_Px_vrr = PBX*I_KINETIC_Hx4y_S_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Px_vrr;
    Double I_KINETIC_Hx3yz_Px_vrr = PBX*I_KINETIC_Hx3yz_S_vrr+oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Px_vrr;
    Double I_KINETIC_Hx2y2z_Px_vrr = PBX*I_KINETIC_Hx2y2z_S_vrr+oned2z*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Px_vrr;
    Double I_KINETIC_Hxy3z_Px_vrr = PBX*I_KINETIC_Hxy3z_S_vrr+oned2z*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Px_vrr;
    Double I_KINETIC_Hx4z_Px_vrr = PBX*I_KINETIC_Hx4z_S_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Px_vrr;
    Double I_KINETIC_H5y_Px_vrr = PBX*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Px_vrr;
    Double I_KINETIC_H4yz_Px_vrr = PBX*I_KINETIC_H4yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Px_vrr;
    Double I_KINETIC_H3y2z_Px_vrr = PBX*I_KINETIC_H3y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Px_vrr;
    Double I_KINETIC_H2y3z_Px_vrr = PBX*I_KINETIC_H2y3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Px_vrr;
    Double I_KINETIC_Hy4z_Px_vrr = PBX*I_KINETIC_Hy4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Px_vrr;
    Double I_KINETIC_H5z_Px_vrr = PBX*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Px_vrr;
    Double I_KINETIC_H5x_Py_vrr = PBY*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Py_vrr;
    Double I_KINETIC_H4xy_Py_vrr = PBY*I_KINETIC_H4xy_S_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Py_vrr;
    Double I_KINETIC_H4xz_Py_vrr = PBY*I_KINETIC_H4xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Py_vrr;
    Double I_KINETIC_H3x2y_Py_vrr = PBY*I_KINETIC_H3x2y_S_vrr+2*oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Py_vrr;
    Double I_KINETIC_H3xyz_Py_vrr = PBY*I_KINETIC_H3xyz_S_vrr+oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Py_vrr;
    Double I_KINETIC_H3x2z_Py_vrr = PBY*I_KINETIC_H3x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Py_vrr;
    Double I_KINETIC_H2x3y_Py_vrr = PBY*I_KINETIC_H2x3y_S_vrr+3*oned2z*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Py_vrr;
    Double I_KINETIC_H2x2yz_Py_vrr = PBY*I_KINETIC_H2x2yz_S_vrr+2*oned2z*I_KINETIC_G2xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Py_vrr;
    Double I_KINETIC_H2xy2z_Py_vrr = PBY*I_KINETIC_H2xy2z_S_vrr+oned2z*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Py_vrr;
    Double I_KINETIC_H2x3z_Py_vrr = PBY*I_KINETIC_H2x3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Py_vrr;
    Double I_KINETIC_Hx4y_Py_vrr = PBY*I_KINETIC_Hx4y_S_vrr+4*oned2z*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Py_vrr;
    Double I_KINETIC_Hx3yz_Py_vrr = PBY*I_KINETIC_Hx3yz_S_vrr+3*oned2z*I_KINETIC_Gx2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Py_vrr;
    Double I_KINETIC_Hx2y2z_Py_vrr = PBY*I_KINETIC_Hx2y2z_S_vrr+2*oned2z*I_KINETIC_Gxy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Py_vrr;
    Double I_KINETIC_Hxy3z_Py_vrr = PBY*I_KINETIC_Hxy3z_S_vrr+oned2z*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Py_vrr;
    Double I_KINETIC_Hx4z_Py_vrr = PBY*I_KINETIC_Hx4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Py_vrr;
    Double I_KINETIC_H5y_Py_vrr = PBY*I_KINETIC_H5y_S_vrr+5*oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Py_vrr;
    Double I_KINETIC_H4yz_Py_vrr = PBY*I_KINETIC_H4yz_S_vrr+4*oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Py_vrr;
    Double I_KINETIC_H3y2z_Py_vrr = PBY*I_KINETIC_H3y2z_S_vrr+3*oned2z*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Py_vrr;
    Double I_KINETIC_H2y3z_Py_vrr = PBY*I_KINETIC_H2y3z_S_vrr+2*oned2z*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Py_vrr;
    Double I_KINETIC_Hy4z_Py_vrr = PBY*I_KINETIC_Hy4z_S_vrr+oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Py_vrr;
    Double I_KINETIC_H5z_Py_vrr = PBY*I_KINETIC_H5z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Py_vrr;
    Double I_KINETIC_H5x_Pz_vrr = PBZ*I_KINETIC_H5x_S_vrr+twoxi*I_TWOBODYOVERLAP_H5x_Pz_vrr;
    Double I_KINETIC_H4xy_Pz_vrr = PBZ*I_KINETIC_H4xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xy_Pz_vrr;
    Double I_KINETIC_H4xz_Pz_vrr = PBZ*I_KINETIC_H4xz_S_vrr+oned2z*I_KINETIC_G4x_S_vrr+twoxi*I_TWOBODYOVERLAP_H4xz_Pz_vrr;
    Double I_KINETIC_H3x2y_Pz_vrr = PBZ*I_KINETIC_H3x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2y_Pz_vrr;
    Double I_KINETIC_H3xyz_Pz_vrr = PBZ*I_KINETIC_H3xyz_S_vrr+oned2z*I_KINETIC_G3xy_S_vrr+twoxi*I_TWOBODYOVERLAP_H3xyz_Pz_vrr;
    Double I_KINETIC_H3x2z_Pz_vrr = PBZ*I_KINETIC_H3x2z_S_vrr+2*oned2z*I_KINETIC_G3xz_S_vrr+twoxi*I_TWOBODYOVERLAP_H3x2z_Pz_vrr;
    Double I_KINETIC_H2x3y_Pz_vrr = PBZ*I_KINETIC_H2x3y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3y_Pz_vrr;
    Double I_KINETIC_H2x2yz_Pz_vrr = PBZ*I_KINETIC_H2x2yz_S_vrr+oned2z*I_KINETIC_G2x2y_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x2yz_Pz_vrr;
    Double I_KINETIC_H2xy2z_Pz_vrr = PBZ*I_KINETIC_H2xy2z_S_vrr+2*oned2z*I_KINETIC_G2xyz_S_vrr+twoxi*I_TWOBODYOVERLAP_H2xy2z_Pz_vrr;
    Double I_KINETIC_H2x3z_Pz_vrr = PBZ*I_KINETIC_H2x3z_S_vrr+3*oned2z*I_KINETIC_G2x2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2x3z_Pz_vrr;
    Double I_KINETIC_Hx4y_Pz_vrr = PBZ*I_KINETIC_Hx4y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4y_Pz_vrr;
    Double I_KINETIC_Hx3yz_Pz_vrr = PBZ*I_KINETIC_Hx3yz_S_vrr+oned2z*I_KINETIC_Gx3y_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx3yz_Pz_vrr;
    Double I_KINETIC_Hx2y2z_Pz_vrr = PBZ*I_KINETIC_Hx2y2z_S_vrr+2*oned2z*I_KINETIC_Gx2yz_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx2y2z_Pz_vrr;
    Double I_KINETIC_Hxy3z_Pz_vrr = PBZ*I_KINETIC_Hxy3z_S_vrr+3*oned2z*I_KINETIC_Gxy2z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hxy3z_Pz_vrr;
    Double I_KINETIC_Hx4z_Pz_vrr = PBZ*I_KINETIC_Hx4z_S_vrr+4*oned2z*I_KINETIC_Gx3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hx4z_Pz_vrr;
    Double I_KINETIC_H5y_Pz_vrr = PBZ*I_KINETIC_H5y_S_vrr+twoxi*I_TWOBODYOVERLAP_H5y_Pz_vrr;
    Double I_KINETIC_H4yz_Pz_vrr = PBZ*I_KINETIC_H4yz_S_vrr+oned2z*I_KINETIC_G4y_S_vrr+twoxi*I_TWOBODYOVERLAP_H4yz_Pz_vrr;
    Double I_KINETIC_H3y2z_Pz_vrr = PBZ*I_KINETIC_H3y2z_S_vrr+2*oned2z*I_KINETIC_G3yz_S_vrr+twoxi*I_TWOBODYOVERLAP_H3y2z_Pz_vrr;
    Double I_KINETIC_H2y3z_Pz_vrr = PBZ*I_KINETIC_H2y3z_S_vrr+3*oned2z*I_KINETIC_G2y2z_S_vrr+twoxi*I_TWOBODYOVERLAP_H2y3z_Pz_vrr;
    Double I_KINETIC_Hy4z_Pz_vrr = PBZ*I_KINETIC_Hy4z_S_vrr+4*oned2z*I_KINETIC_Gy3z_S_vrr+twoxi*I_TWOBODYOVERLAP_Hy4z_Pz_vrr;
    Double I_KINETIC_H5z_Pz_vrr = PBZ*I_KINETIC_H5z_S_vrr+5*oned2z*I_KINETIC_G4z_S_vrr+twoxi*I_TWOBODYOVERLAP_H5z_Pz_vrr;

    /************************************************************
     * shell quartet name: SQ_KINETIC_H_P
     * doing contraction work for VRR part 
     * totally 0 integrals are omitted 
     ************************************************************/
    abcd[0] += I_KINETIC_H5x_Px_vrr;
    abcd[1] += I_KINETIC_H4xy_Px_vrr;
    abcd[2] += I_KINETIC_H4xz_Px_vrr;
    abcd[3] += I_KINETIC_H3x2y_Px_vrr;
    abcd[4] += I_KINETIC_H3xyz_Px_vrr;
    abcd[5] += I_KINETIC_H3x2z_Px_vrr;
    abcd[6] += I_KINETIC_H2x3y_Px_vrr;
    abcd[7] += I_KINETIC_H2x2yz_Px_vrr;
    abcd[8] += I_KINETIC_H2xy2z_Px_vrr;
    abcd[9] += I_KINETIC_H2x3z_Px_vrr;
    abcd[10] += I_KINETIC_Hx4y_Px_vrr;
    abcd[11] += I_KINETIC_Hx3yz_Px_vrr;
    abcd[12] += I_KINETIC_Hx2y2z_Px_vrr;
    abcd[13] += I_KINETIC_Hxy3z_Px_vrr;
    abcd[14] += I_KINETIC_Hx4z_Px_vrr;
    abcd[15] += I_KINETIC_H5y_Px_vrr;
    abcd[16] += I_KINETIC_H4yz_Px_vrr;
    abcd[17] += I_KINETIC_H3y2z_Px_vrr;
    abcd[18] += I_KINETIC_H2y3z_Px_vrr;
    abcd[19] += I_KINETIC_Hy4z_Px_vrr;
    abcd[20] += I_KINETIC_H5z_Px_vrr;
    abcd[21] += I_KINETIC_H5x_Py_vrr;
    abcd[22] += I_KINETIC_H4xy_Py_vrr;
    abcd[23] += I_KINETIC_H4xz_Py_vrr;
    abcd[24] += I_KINETIC_H3x2y_Py_vrr;
    abcd[25] += I_KINETIC_H3xyz_Py_vrr;
    abcd[26] += I_KINETIC_H3x2z_Py_vrr;
    abcd[27] += I_KINETIC_H2x3y_Py_vrr;
    abcd[28] += I_KINETIC_H2x2yz_Py_vrr;
    abcd[29] += I_KINETIC_H2xy2z_Py_vrr;
    abcd[30] += I_KINETIC_H2x3z_Py_vrr;
    abcd[31] += I_KINETIC_Hx4y_Py_vrr;
    abcd[32] += I_KINETIC_Hx3yz_Py_vrr;
    abcd[33] += I_KINETIC_Hx2y2z_Py_vrr;
    abcd[34] += I_KINETIC_Hxy3z_Py_vrr;
    abcd[35] += I_KINETIC_Hx4z_Py_vrr;
    abcd[36] += I_KINETIC_H5y_Py_vrr;
    abcd[37] += I_KINETIC_H4yz_Py_vrr;
    abcd[38] += I_KINETIC_H3y2z_Py_vrr;
    abcd[39] += I_KINETIC_H2y3z_Py_vrr;
    abcd[40] += I_KINETIC_Hy4z_Py_vrr;
    abcd[41] += I_KINETIC_H5z_Py_vrr;
    abcd[42] += I_KINETIC_H5x_Pz_vrr;
    abcd[43] += I_KINETIC_H4xy_Pz_vrr;
    abcd[44] += I_KINETIC_H4xz_Pz_vrr;
    abcd[45] += I_KINETIC_H3x2y_Pz_vrr;
    abcd[46] += I_KINETIC_H3xyz_Pz_vrr;
    abcd[47] += I_KINETIC_H3x2z_Pz_vrr;
    abcd[48] += I_KINETIC_H2x3y_Pz_vrr;
    abcd[49] += I_KINETIC_H2x2yz_Pz_vrr;
    abcd[50] += I_KINETIC_H2xy2z_Pz_vrr;
    abcd[51] += I_KINETIC_H2x3z_Pz_vrr;
    abcd[52] += I_KINETIC_Hx4y_Pz_vrr;
    abcd[53] += I_KINETIC_Hx3yz_Pz_vrr;
    abcd[54] += I_KINETIC_Hx2y2z_Pz_vrr;
    abcd[55] += I_KINETIC_Hxy3z_Pz_vrr;
    abcd[56] += I_KINETIC_Hx4z_Pz_vrr;
    abcd[57] += I_KINETIC_H5y_Pz_vrr;
    abcd[58] += I_KINETIC_H4yz_Pz_vrr;
    abcd[59] += I_KINETIC_H3y2z_Pz_vrr;
    abcd[60] += I_KINETIC_H2y3z_Pz_vrr;
    abcd[61] += I_KINETIC_Hy4z_Pz_vrr;
    abcd[62] += I_KINETIC_H5z_Pz_vrr;
  }
}
