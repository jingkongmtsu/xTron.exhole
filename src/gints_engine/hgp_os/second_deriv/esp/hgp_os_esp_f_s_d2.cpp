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
// BRA1 as redundant position, total RHS integrals evaluated as: 1206
// BRA2 as redundant position, total RHS integrals evaluated as: 990
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

void hgp_os_esp_f_s_d2(const UInt& inp2, const UInt& nGrids, const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* R, Double* abcd)
{
  // loop over grid points 
  for(UInt iGrid=0; iGrid<nGrids; iGrid++) {

    //
    // declare the variables as result of VRR process
    //
    Double I_ESP_H5x_S_aa = 0.0E0;
    Double I_ESP_H4xy_S_aa = 0.0E0;
    Double I_ESP_H4xz_S_aa = 0.0E0;
    Double I_ESP_H3x2y_S_aa = 0.0E0;
    Double I_ESP_H3xyz_S_aa = 0.0E0;
    Double I_ESP_H3x2z_S_aa = 0.0E0;
    Double I_ESP_H2x3y_S_aa = 0.0E0;
    Double I_ESP_H2x2yz_S_aa = 0.0E0;
    Double I_ESP_H2xy2z_S_aa = 0.0E0;
    Double I_ESP_H2x3z_S_aa = 0.0E0;
    Double I_ESP_Hx4y_S_aa = 0.0E0;
    Double I_ESP_Hx3yz_S_aa = 0.0E0;
    Double I_ESP_Hx2y2z_S_aa = 0.0E0;
    Double I_ESP_Hxy3z_S_aa = 0.0E0;
    Double I_ESP_Hx4z_S_aa = 0.0E0;
    Double I_ESP_H5y_S_aa = 0.0E0;
    Double I_ESP_H4yz_S_aa = 0.0E0;
    Double I_ESP_H3y2z_S_aa = 0.0E0;
    Double I_ESP_H2y3z_S_aa = 0.0E0;
    Double I_ESP_Hy4z_S_aa = 0.0E0;
    Double I_ESP_H5z_S_aa = 0.0E0;
    Double I_ESP_F3x_S_a = 0.0E0;
    Double I_ESP_F2xy_S_a = 0.0E0;
    Double I_ESP_F2xz_S_a = 0.0E0;
    Double I_ESP_Fx2y_S_a = 0.0E0;
    Double I_ESP_Fxyz_S_a = 0.0E0;
    Double I_ESP_Fx2z_S_a = 0.0E0;
    Double I_ESP_F3y_S_a = 0.0E0;
    Double I_ESP_F2yz_S_a = 0.0E0;
    Double I_ESP_Fy2z_S_a = 0.0E0;
    Double I_ESP_F3z_S_a = 0.0E0;
    Double I_ESP_Px_S = 0.0E0;
    Double I_ESP_Py_S = 0.0E0;
    Double I_ESP_Pz_S = 0.0E0;

    for(UInt ip2=0; ip2<inp2; ip2++) {
      Double ic2   = icoe[ip2];
      Double onedz = iexp[ip2];
      Double rho   = 1.0E0/onedz;
      Double zeta  = rho;
      Double zdiff = iexpdiff[ip2];
      Double alpha = 0.5E0*(zeta+zdiff);
      Double beta  = 0.5E0*(zeta-zdiff);
      Double sqrho = sqrt(rho);
      Double fbra  = ifac[ip2];
      Double oned2z= 0.5E0*onedz;
      UInt offsetP = 3*ip2;
      Double PX    = P[offsetP  ];
      Double PY    = P[offsetP+1];
      Double PZ    = P[offsetP+2];
      Double PAX   = PX - A[0];
      Double PAY   = PY - A[1];
      Double PAZ   = PZ - A[2];
      Double PRX   = PX - R[iGrid*3  ];
      Double PRY   = PY - R[iGrid*3+1];
      Double PRZ   = PZ - R[iGrid*3+2];
      Double PR2   = PRX*PRX+PRY*PRY+PRZ*PRZ;
      Double u     = rho*PR2;
      Double squ   = sqrt(u);
      Double prefactor = ic2*fbra;

      //
      //
      // now here for maxM>0 to compute the infamous incomplete Gamma function f_{m}(u)
      // the implementation is divided in two situations:
      // 1  if u <=1.8; use power series to get f_{Mmax}(u), then use down recursive
      //    relation to get the rest of incomplete Gamma functions;
      // 2  for u >1.8 and M <= 10 we calculate erf(u), then use up recursive
      //    relation to calculate the rest of results
      // 3  for u> 1.8 and M >  10 we calculate f_{Mmax}(u) then use down 
      //    recursive relation to get rest of incomplete Gamma functions 
      // The above procedure is tested for u between 0 to 40 with step length 1.0E-6
      // (or 1.0E-5 for float double data), for up recursive relation it shows the error
      // within 1.0E-12 (for M_limit = 12 or error within 1.0E-6 for float type of data
      // For the polynomial expansion and down recursive procedure the error is within 
      // 1.0E-14. All of the testing details please refer to the fmt_test folder
      // 
      // There's one thing need to note for up recursive process. We found that the up
      // recursive procedure is only stable for maxM<=10 and u>1.8 with double
      // precision data, single precision data will lose accuracy quickly so the result
      // for single precision calculation is not doable. Therefore if the "WITH_SINGLE_PRECISION"
      // is defined, then for erf function calculation as well as up recursive
      // process we will use the double type of data
      // 
      //

      Double I_ESP_S_S_vrr  = 0.0E0;
      Double I_ESP_S_S_M1_vrr  = 0.0E0;
      Double I_ESP_S_S_M2_vrr  = 0.0E0;
      Double I_ESP_S_S_M3_vrr  = 0.0E0;
      Double I_ESP_S_S_M4_vrr  = 0.0E0;
      Double I_ESP_S_S_M5_vrr  = 0.0E0;

      // declare double type of results to store the accuracy
#ifdef WITH_SINGLE_PRECISION
      double I_ESP_S_S_vrr_d  = 0.0E0;
      double I_ESP_S_S_M1_vrr_d  = 0.0E0;
      double I_ESP_S_S_M2_vrr_d  = 0.0E0;
      double I_ESP_S_S_M3_vrr_d  = 0.0E0;
      double I_ESP_S_S_M4_vrr_d  = 0.0E0;
      double I_ESP_S_S_M5_vrr_d  = 0.0E0;
#endif

      if (u<=1.8E0) {

        // calculate (SS|SS)^{Mmax}
        // use 18 terms power series to expand the (SS|SS)^{Mmax}
        Double u2 = 2.0E0*u;
        I_ESP_S_S_M5_vrr = 1.0E0+u2*ONEOVER45;
        I_ESP_S_S_M5_vrr = 1.0E0+u2*ONEOVER43*I_ESP_S_S_M5_vrr;
        I_ESP_S_S_M5_vrr = 1.0E0+u2*ONEOVER41*I_ESP_S_S_M5_vrr;
        I_ESP_S_S_M5_vrr = 1.0E0+u2*ONEOVER39*I_ESP_S_S_M5_vrr;
        I_ESP_S_S_M5_vrr = 1.0E0+u2*ONEOVER37*I_ESP_S_S_M5_vrr;
        I_ESP_S_S_M5_vrr = 1.0E0+u2*ONEOVER35*I_ESP_S_S_M5_vrr;
        I_ESP_S_S_M5_vrr = 1.0E0+u2*ONEOVER33*I_ESP_S_S_M5_vrr;
        I_ESP_S_S_M5_vrr = 1.0E0+u2*ONEOVER31*I_ESP_S_S_M5_vrr;
        I_ESP_S_S_M5_vrr = 1.0E0+u2*ONEOVER29*I_ESP_S_S_M5_vrr;
        I_ESP_S_S_M5_vrr = 1.0E0+u2*ONEOVER27*I_ESP_S_S_M5_vrr;
        I_ESP_S_S_M5_vrr = 1.0E0+u2*ONEOVER25*I_ESP_S_S_M5_vrr;
        I_ESP_S_S_M5_vrr = 1.0E0+u2*ONEOVER23*I_ESP_S_S_M5_vrr;
        I_ESP_S_S_M5_vrr = 1.0E0+u2*ONEOVER21*I_ESP_S_S_M5_vrr;
        I_ESP_S_S_M5_vrr = 1.0E0+u2*ONEOVER19*I_ESP_S_S_M5_vrr;
        I_ESP_S_S_M5_vrr = 1.0E0+u2*ONEOVER17*I_ESP_S_S_M5_vrr;
        I_ESP_S_S_M5_vrr = 1.0E0+u2*ONEOVER15*I_ESP_S_S_M5_vrr;
        I_ESP_S_S_M5_vrr = 1.0E0+u2*ONEOVER13*I_ESP_S_S_M5_vrr;
        I_ESP_S_S_M5_vrr = ONEOVER11*I_ESP_S_S_M5_vrr;
        Double eu = exp(-u);
        Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ESP_S_S_M5_vrr  = f*I_ESP_S_S_M5_vrr;

        // now use down recursive relation to get
        // rest of (SS|SS)^{m}
        I_ESP_S_S_M4_vrr  = ONEOVER9*(u2*I_ESP_S_S_M5_vrr+f);
        I_ESP_S_S_M3_vrr  = ONEOVER7*(u2*I_ESP_S_S_M4_vrr+f);
        I_ESP_S_S_M2_vrr  = ONEOVER5*(u2*I_ESP_S_S_M3_vrr+f);
        I_ESP_S_S_M1_vrr  = ONEOVER3*(u2*I_ESP_S_S_M2_vrr+f);
        I_ESP_S_S_vrr  = ONEOVER1*(u2*I_ESP_S_S_M1_vrr+f);

      }else{
#ifdef WITH_SINGLE_PRECISION

        // recompute the variable in terms of double accuracy
        double u_d     = u;
        double rho_d   = rho;
        double fac_d   = prefactor;
        double sqrho_d = sqrt(rho_d);
        double squ_d   = sqrt(u_d);

        // use erf function to get (SS|SS)^{0}
        if (fabs(u_d)<THRESHOLD_MATH) {
          I_ESP_S_S_vrr_d = fac_d*sqrho_d*TWOOVERSQRTPI;
        }else{
          I_ESP_S_S_vrr_d = (fac_d*sqrho_d/squ_d)*erf(squ_d);
        }

        // now use up recursive relation to get
        // rest of (SS|SS)^{m}
        double oneO2u  = 0.5E0/u_d;
        double eu      = exp(-u_d);
        double f       = TWOOVERSQRTPI*fac_d*sqrho_d*eu;
        I_ESP_S_S_M1_vrr_d = oneO2u*(1.0E0*I_ESP_S_S_vrr_d-f);
        I_ESP_S_S_M2_vrr_d = oneO2u*(3.0E0*I_ESP_S_S_M1_vrr_d-f);
        I_ESP_S_S_M3_vrr_d = oneO2u*(5.0E0*I_ESP_S_S_M2_vrr_d-f);
        I_ESP_S_S_M4_vrr_d = oneO2u*(7.0E0*I_ESP_S_S_M3_vrr_d-f);
        I_ESP_S_S_M5_vrr_d = oneO2u*(9.0E0*I_ESP_S_S_M4_vrr_d-f);

        // write the double result back to the float var
        I_ESP_S_S_vrr = static_cast<Double>(I_ESP_S_S_vrr_d);
        I_ESP_S_S_M1_vrr = static_cast<Double>(I_ESP_S_S_M1_vrr_d);
        I_ESP_S_S_M2_vrr = static_cast<Double>(I_ESP_S_S_M2_vrr_d);
        I_ESP_S_S_M3_vrr = static_cast<Double>(I_ESP_S_S_M3_vrr_d);
        I_ESP_S_S_M4_vrr = static_cast<Double>(I_ESP_S_S_M4_vrr_d);
        I_ESP_S_S_M5_vrr = static_cast<Double>(I_ESP_S_S_M5_vrr_d);

#else

        // use erf function to get (SS|SS)^{0}
        if (fabs(u)<THRESHOLD_MATH) {
          I_ESP_S_S_vrr = prefactor*sqrho*TWOOVERSQRTPI;
        }else{
          I_ESP_S_S_vrr = (prefactor*sqrho/squ)*erf(squ);
        }

        // now use up recursive relation to get
        // rest of (SS|SS)^{m}
        Double oneO2u = 0.5E0/u;
        Double eu     = exp(-u);
        Double f      = TWOOVERSQRTPI*prefactor*sqrho*eu;
        I_ESP_S_S_M1_vrr = oneO2u*(1.0E0*I_ESP_S_S_vrr-f);
        I_ESP_S_S_M2_vrr = oneO2u*(3.0E0*I_ESP_S_S_M1_vrr-f);
        I_ESP_S_S_M3_vrr = oneO2u*(5.0E0*I_ESP_S_S_M2_vrr-f);
        I_ESP_S_S_M4_vrr = oneO2u*(7.0E0*I_ESP_S_S_M3_vrr-f);
        I_ESP_S_S_M5_vrr = oneO2u*(9.0E0*I_ESP_S_S_M4_vrr-f);

#endif

      }


        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M4
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M4
         * RHS shell quartet name: SQ_ESP_S_S_M5
         ************************************************************/
        Double I_ESP_Px_S_M4_vrr = PAX*I_ESP_S_S_M4_vrr-PRX*I_ESP_S_S_M5_vrr;
        Double I_ESP_Py_S_M4_vrr = PAY*I_ESP_S_S_M4_vrr-PRY*I_ESP_S_S_M5_vrr;
        Double I_ESP_Pz_S_M4_vrr = PAZ*I_ESP_S_S_M4_vrr-PRZ*I_ESP_S_S_M5_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M3
         * RHS shell quartet name: SQ_ESP_S_S_M4
         ************************************************************/
        Double I_ESP_Px_S_M3_vrr = PAX*I_ESP_S_S_M3_vrr-PRX*I_ESP_S_S_M4_vrr;
        Double I_ESP_Py_S_M3_vrr = PAY*I_ESP_S_S_M3_vrr-PRY*I_ESP_S_S_M4_vrr;
        Double I_ESP_Pz_S_M3_vrr = PAZ*I_ESP_S_S_M3_vrr-PRZ*I_ESP_S_S_M4_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M3
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M3
         * RHS shell quartet name: SQ_ESP_P_S_M4
         * RHS shell quartet name: SQ_ESP_S_S_M3
         * RHS shell quartet name: SQ_ESP_S_S_M4
         ************************************************************/
        Double I_ESP_D2x_S_M3_vrr = PAX*I_ESP_Px_S_M3_vrr-PRX*I_ESP_Px_S_M4_vrr+oned2z*I_ESP_S_S_M3_vrr-oned2z*I_ESP_S_S_M4_vrr;
        Double I_ESP_D2y_S_M3_vrr = PAY*I_ESP_Py_S_M3_vrr-PRY*I_ESP_Py_S_M4_vrr+oned2z*I_ESP_S_S_M3_vrr-oned2z*I_ESP_S_S_M4_vrr;
        Double I_ESP_D2z_S_M3_vrr = PAZ*I_ESP_Pz_S_M3_vrr-PRZ*I_ESP_Pz_S_M4_vrr+oned2z*I_ESP_S_S_M3_vrr-oned2z*I_ESP_S_S_M4_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M2
         * RHS shell quartet name: SQ_ESP_S_S_M3
         ************************************************************/
        Double I_ESP_Px_S_M2_vrr = PAX*I_ESP_S_S_M2_vrr-PRX*I_ESP_S_S_M3_vrr;
        Double I_ESP_Py_S_M2_vrr = PAY*I_ESP_S_S_M2_vrr-PRY*I_ESP_S_S_M3_vrr;
        Double I_ESP_Pz_S_M2_vrr = PAZ*I_ESP_S_S_M2_vrr-PRZ*I_ESP_S_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M2
         * RHS shell quartet name: SQ_ESP_P_S_M3
         * RHS shell quartet name: SQ_ESP_S_S_M2
         * RHS shell quartet name: SQ_ESP_S_S_M3
         ************************************************************/
        Double I_ESP_D2x_S_M2_vrr = PAX*I_ESP_Px_S_M2_vrr-PRX*I_ESP_Px_S_M3_vrr+oned2z*I_ESP_S_S_M2_vrr-oned2z*I_ESP_S_S_M3_vrr;
        Double I_ESP_D2y_S_M2_vrr = PAY*I_ESP_Py_S_M2_vrr-PRY*I_ESP_Py_S_M3_vrr+oned2z*I_ESP_S_S_M2_vrr-oned2z*I_ESP_S_S_M3_vrr;
        Double I_ESP_D2z_S_M2_vrr = PAZ*I_ESP_Pz_S_M2_vrr-PRZ*I_ESP_Pz_S_M3_vrr+oned2z*I_ESP_S_S_M2_vrr-oned2z*I_ESP_S_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_M2
         * expanding position: BRA1
         * code section is: VRR
         * totally 4 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M2
         * RHS shell quartet name: SQ_ESP_D_S_M3
         * RHS shell quartet name: SQ_ESP_P_S_M2
         * RHS shell quartet name: SQ_ESP_P_S_M3
         ************************************************************/
        Double I_ESP_F3x_S_M2_vrr = PAX*I_ESP_D2x_S_M2_vrr-PRX*I_ESP_D2x_S_M3_vrr+2*oned2z*I_ESP_Px_S_M2_vrr-2*oned2z*I_ESP_Px_S_M3_vrr;
        Double I_ESP_F2xy_S_M2_vrr = PAY*I_ESP_D2x_S_M2_vrr-PRY*I_ESP_D2x_S_M3_vrr;
        Double I_ESP_F2xz_S_M2_vrr = PAZ*I_ESP_D2x_S_M2_vrr-PRZ*I_ESP_D2x_S_M3_vrr;
        Double I_ESP_F3y_S_M2_vrr = PAY*I_ESP_D2y_S_M2_vrr-PRY*I_ESP_D2y_S_M3_vrr+2*oned2z*I_ESP_Py_S_M2_vrr-2*oned2z*I_ESP_Py_S_M3_vrr;
        Double I_ESP_F2yz_S_M2_vrr = PAZ*I_ESP_D2y_S_M2_vrr-PRZ*I_ESP_D2y_S_M3_vrr;
        Double I_ESP_F3z_S_M2_vrr = PAZ*I_ESP_D2z_S_M2_vrr-PRZ*I_ESP_D2z_S_M3_vrr+2*oned2z*I_ESP_Pz_S_M2_vrr-2*oned2z*I_ESP_Pz_S_M3_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S_M1
         * RHS shell quartet name: SQ_ESP_S_S_M2
         ************************************************************/
        Double I_ESP_Px_S_M1_vrr = PAX*I_ESP_S_S_M1_vrr-PRX*I_ESP_S_S_M2_vrr;
        Double I_ESP_Py_S_M1_vrr = PAY*I_ESP_S_S_M1_vrr-PRY*I_ESP_S_S_M2_vrr;
        Double I_ESP_Pz_S_M1_vrr = PAZ*I_ESP_S_S_M1_vrr-PRZ*I_ESP_S_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 2 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_P_S_M2
         * RHS shell quartet name: SQ_ESP_S_S_M1
         * RHS shell quartet name: SQ_ESP_S_S_M2
         ************************************************************/
        Double I_ESP_D2x_S_M1_vrr = PAX*I_ESP_Px_S_M1_vrr-PRX*I_ESP_Px_S_M2_vrr+oned2z*I_ESP_S_S_M1_vrr-oned2z*I_ESP_S_S_M2_vrr;
        Double I_ESP_Dxy_S_M1_vrr = PAY*I_ESP_Px_S_M1_vrr-PRY*I_ESP_Px_S_M2_vrr;
        Double I_ESP_D2y_S_M1_vrr = PAY*I_ESP_Py_S_M1_vrr-PRY*I_ESP_Py_S_M2_vrr+oned2z*I_ESP_S_S_M1_vrr-oned2z*I_ESP_S_S_M2_vrr;
        Double I_ESP_D2z_S_M1_vrr = PAZ*I_ESP_Pz_S_M1_vrr-PRZ*I_ESP_Pz_S_M2_vrr+oned2z*I_ESP_S_S_M1_vrr-oned2z*I_ESP_S_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 4 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S_M1
         * RHS shell quartet name: SQ_ESP_D_S_M2
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_P_S_M2
         ************************************************************/
        Double I_ESP_F3x_S_M1_vrr = PAX*I_ESP_D2x_S_M1_vrr-PRX*I_ESP_D2x_S_M2_vrr+2*oned2z*I_ESP_Px_S_M1_vrr-2*oned2z*I_ESP_Px_S_M2_vrr;
        Double I_ESP_F2xy_S_M1_vrr = PAY*I_ESP_D2x_S_M1_vrr-PRY*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_F2xz_S_M1_vrr = PAZ*I_ESP_D2x_S_M1_vrr-PRZ*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_F3y_S_M1_vrr = PAY*I_ESP_D2y_S_M1_vrr-PRY*I_ESP_D2y_S_M2_vrr+2*oned2z*I_ESP_Py_S_M1_vrr-2*oned2z*I_ESP_Py_S_M2_vrr;
        Double I_ESP_F2yz_S_M1_vrr = PAZ*I_ESP_D2y_S_M1_vrr-PRZ*I_ESP_D2y_S_M2_vrr;
        Double I_ESP_F3z_S_M1_vrr = PAZ*I_ESP_D2z_S_M1_vrr-PRZ*I_ESP_D2z_S_M2_vrr+2*oned2z*I_ESP_Pz_S_M1_vrr-2*oned2z*I_ESP_Pz_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S_M1
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S_M1
         * RHS shell quartet name: SQ_ESP_F_S_M2
         * RHS shell quartet name: SQ_ESP_D_S_M1
         * RHS shell quartet name: SQ_ESP_D_S_M2
         ************************************************************/
        Double I_ESP_G4x_S_M1_vrr = PAX*I_ESP_F3x_S_M1_vrr-PRX*I_ESP_F3x_S_M2_vrr+3*oned2z*I_ESP_D2x_S_M1_vrr-3*oned2z*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_G3xy_S_M1_vrr = PAY*I_ESP_F3x_S_M1_vrr-PRY*I_ESP_F3x_S_M2_vrr;
        Double I_ESP_G3xz_S_M1_vrr = PAZ*I_ESP_F3x_S_M1_vrr-PRZ*I_ESP_F3x_S_M2_vrr;
        Double I_ESP_G2x2y_S_M1_vrr = PAY*I_ESP_F2xy_S_M1_vrr-PRY*I_ESP_F2xy_S_M2_vrr+oned2z*I_ESP_D2x_S_M1_vrr-oned2z*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_G2x2z_S_M1_vrr = PAZ*I_ESP_F2xz_S_M1_vrr-PRZ*I_ESP_F2xz_S_M2_vrr+oned2z*I_ESP_D2x_S_M1_vrr-oned2z*I_ESP_D2x_S_M2_vrr;
        Double I_ESP_Gx3y_S_M1_vrr = PAX*I_ESP_F3y_S_M1_vrr-PRX*I_ESP_F3y_S_M2_vrr;
        Double I_ESP_Gx3z_S_M1_vrr = PAX*I_ESP_F3z_S_M1_vrr-PRX*I_ESP_F3z_S_M2_vrr;
        Double I_ESP_G4y_S_M1_vrr = PAY*I_ESP_F3y_S_M1_vrr-PRY*I_ESP_F3y_S_M2_vrr+3*oned2z*I_ESP_D2y_S_M1_vrr-3*oned2z*I_ESP_D2y_S_M2_vrr;
        Double I_ESP_G3yz_S_M1_vrr = PAZ*I_ESP_F3y_S_M1_vrr-PRZ*I_ESP_F3y_S_M2_vrr;
        Double I_ESP_G2y2z_S_M1_vrr = PAZ*I_ESP_F2yz_S_M1_vrr-PRZ*I_ESP_F2yz_S_M2_vrr+oned2z*I_ESP_D2y_S_M1_vrr-oned2z*I_ESP_D2y_S_M2_vrr;
        Double I_ESP_Gy3z_S_M1_vrr = PAY*I_ESP_F3z_S_M1_vrr-PRY*I_ESP_F3z_S_M2_vrr;
        Double I_ESP_G4z_S_M1_vrr = PAZ*I_ESP_F3z_S_M1_vrr-PRZ*I_ESP_F3z_S_M2_vrr+3*oned2z*I_ESP_D2z_S_M1_vrr-3*oned2z*I_ESP_D2z_S_M2_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_S_S
         * RHS shell quartet name: SQ_ESP_S_S_M1
         ************************************************************/
        Double I_ESP_Px_S_vrr = PAX*I_ESP_S_S_vrr-PRX*I_ESP_S_S_M1_vrr;
        Double I_ESP_Py_S_vrr = PAY*I_ESP_S_S_vrr-PRY*I_ESP_S_S_M1_vrr;
        Double I_ESP_Pz_S_vrr = PAZ*I_ESP_S_S_vrr-PRZ*I_ESP_S_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_D_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 2 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_P_S
         * RHS shell quartet name: SQ_ESP_P_S_M1
         * RHS shell quartet name: SQ_ESP_S_S
         * RHS shell quartet name: SQ_ESP_S_S_M1
         ************************************************************/
        Double I_ESP_D2x_S_vrr = PAX*I_ESP_Px_S_vrr-PRX*I_ESP_Px_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;
        Double I_ESP_Dxy_S_vrr = PAY*I_ESP_Px_S_vrr-PRY*I_ESP_Px_S_M1_vrr;
        Double I_ESP_D2y_S_vrr = PAY*I_ESP_Py_S_vrr-PRY*I_ESP_Py_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;
        Double I_ESP_D2z_S_vrr = PAZ*I_ESP_Pz_S_vrr-PRZ*I_ESP_Pz_S_M1_vrr+oned2z*I_ESP_S_S_vrr-oned2z*I_ESP_S_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_D_S
         * RHS shell quartet name: SQ_ESP_D_S_M1
         * RHS shell quartet name: SQ_ESP_P_S
         * RHS shell quartet name: SQ_ESP_P_S_M1
         ************************************************************/
        Double I_ESP_F3x_S_vrr = PAX*I_ESP_D2x_S_vrr-PRX*I_ESP_D2x_S_M1_vrr+2*oned2z*I_ESP_Px_S_vrr-2*oned2z*I_ESP_Px_S_M1_vrr;
        Double I_ESP_F2xy_S_vrr = PAY*I_ESP_D2x_S_vrr-PRY*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_F2xz_S_vrr = PAZ*I_ESP_D2x_S_vrr-PRZ*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_Fx2y_S_vrr = PAX*I_ESP_D2y_S_vrr-PRX*I_ESP_D2y_S_M1_vrr;
        Double I_ESP_Fxyz_S_vrr = PAZ*I_ESP_Dxy_S_vrr-PRZ*I_ESP_Dxy_S_M1_vrr;
        Double I_ESP_Fx2z_S_vrr = PAX*I_ESP_D2z_S_vrr-PRX*I_ESP_D2z_S_M1_vrr;
        Double I_ESP_F3y_S_vrr = PAY*I_ESP_D2y_S_vrr-PRY*I_ESP_D2y_S_M1_vrr+2*oned2z*I_ESP_Py_S_vrr-2*oned2z*I_ESP_Py_S_M1_vrr;
        Double I_ESP_F2yz_S_vrr = PAZ*I_ESP_D2y_S_vrr-PRZ*I_ESP_D2y_S_M1_vrr;
        Double I_ESP_Fy2z_S_vrr = PAY*I_ESP_D2z_S_vrr-PRY*I_ESP_D2z_S_M1_vrr;
        Double I_ESP_F3z_S_vrr = PAZ*I_ESP_D2z_S_vrr-PRZ*I_ESP_D2z_S_M1_vrr+2*oned2z*I_ESP_Pz_S_vrr-2*oned2z*I_ESP_Pz_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_G_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 3 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_F_S
         * RHS shell quartet name: SQ_ESP_F_S_M1
         * RHS shell quartet name: SQ_ESP_D_S
         * RHS shell quartet name: SQ_ESP_D_S_M1
         ************************************************************/
        Double I_ESP_G4x_S_vrr = PAX*I_ESP_F3x_S_vrr-PRX*I_ESP_F3x_S_M1_vrr+3*oned2z*I_ESP_D2x_S_vrr-3*oned2z*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_G3xy_S_vrr = PAY*I_ESP_F3x_S_vrr-PRY*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_G3xz_S_vrr = PAZ*I_ESP_F3x_S_vrr-PRZ*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_G2x2y_S_vrr = PAY*I_ESP_F2xy_S_vrr-PRY*I_ESP_F2xy_S_M1_vrr+oned2z*I_ESP_D2x_S_vrr-oned2z*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_G2x2z_S_vrr = PAZ*I_ESP_F2xz_S_vrr-PRZ*I_ESP_F2xz_S_M1_vrr+oned2z*I_ESP_D2x_S_vrr-oned2z*I_ESP_D2x_S_M1_vrr;
        Double I_ESP_Gx3y_S_vrr = PAX*I_ESP_F3y_S_vrr-PRX*I_ESP_F3y_S_M1_vrr;
        Double I_ESP_Gx3z_S_vrr = PAX*I_ESP_F3z_S_vrr-PRX*I_ESP_F3z_S_M1_vrr;
        Double I_ESP_G4y_S_vrr = PAY*I_ESP_F3y_S_vrr-PRY*I_ESP_F3y_S_M1_vrr+3*oned2z*I_ESP_D2y_S_vrr-3*oned2z*I_ESP_D2y_S_M1_vrr;
        Double I_ESP_G3yz_S_vrr = PAZ*I_ESP_F3y_S_vrr-PRZ*I_ESP_F3y_S_M1_vrr;
        Double I_ESP_G2y2z_S_vrr = PAZ*I_ESP_F2yz_S_vrr-PRZ*I_ESP_F2yz_S_M1_vrr+oned2z*I_ESP_D2y_S_vrr-oned2z*I_ESP_D2y_S_M1_vrr;
        Double I_ESP_Gy3z_S_vrr = PAY*I_ESP_F3z_S_vrr-PRY*I_ESP_F3z_S_M1_vrr;
        Double I_ESP_G4z_S_vrr = PAZ*I_ESP_F3z_S_vrr-PRZ*I_ESP_F3z_S_M1_vrr+3*oned2z*I_ESP_D2z_S_vrr-3*oned2z*I_ESP_D2z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S
         * expanding position: BRA1
         * code section is: VRR
         * totally 0 integrals are omitted 
         * RHS shell quartet name: SQ_ESP_G_S
         * RHS shell quartet name: SQ_ESP_G_S_M1
         * RHS shell quartet name: SQ_ESP_F_S
         * RHS shell quartet name: SQ_ESP_F_S_M1
         ************************************************************/
        Double I_ESP_H5x_S_vrr = PAX*I_ESP_G4x_S_vrr-PRX*I_ESP_G4x_S_M1_vrr+4*oned2z*I_ESP_F3x_S_vrr-4*oned2z*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_H4xy_S_vrr = PAY*I_ESP_G4x_S_vrr-PRY*I_ESP_G4x_S_M1_vrr;
        Double I_ESP_H4xz_S_vrr = PAZ*I_ESP_G4x_S_vrr-PRZ*I_ESP_G4x_S_M1_vrr;
        Double I_ESP_H3x2y_S_vrr = PAY*I_ESP_G3xy_S_vrr-PRY*I_ESP_G3xy_S_M1_vrr+oned2z*I_ESP_F3x_S_vrr-oned2z*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_H3xyz_S_vrr = PAZ*I_ESP_G3xy_S_vrr-PRZ*I_ESP_G3xy_S_M1_vrr;
        Double I_ESP_H3x2z_S_vrr = PAZ*I_ESP_G3xz_S_vrr-PRZ*I_ESP_G3xz_S_M1_vrr+oned2z*I_ESP_F3x_S_vrr-oned2z*I_ESP_F3x_S_M1_vrr;
        Double I_ESP_H2x3y_S_vrr = PAX*I_ESP_Gx3y_S_vrr-PRX*I_ESP_Gx3y_S_M1_vrr+oned2z*I_ESP_F3y_S_vrr-oned2z*I_ESP_F3y_S_M1_vrr;
        Double I_ESP_H2x2yz_S_vrr = PAZ*I_ESP_G2x2y_S_vrr-PRZ*I_ESP_G2x2y_S_M1_vrr;
        Double I_ESP_H2xy2z_S_vrr = PAY*I_ESP_G2x2z_S_vrr-PRY*I_ESP_G2x2z_S_M1_vrr;
        Double I_ESP_H2x3z_S_vrr = PAX*I_ESP_Gx3z_S_vrr-PRX*I_ESP_Gx3z_S_M1_vrr+oned2z*I_ESP_F3z_S_vrr-oned2z*I_ESP_F3z_S_M1_vrr;
        Double I_ESP_Hx4y_S_vrr = PAX*I_ESP_G4y_S_vrr-PRX*I_ESP_G4y_S_M1_vrr;
        Double I_ESP_Hx3yz_S_vrr = PAZ*I_ESP_Gx3y_S_vrr-PRZ*I_ESP_Gx3y_S_M1_vrr;
        Double I_ESP_Hx2y2z_S_vrr = PAX*I_ESP_G2y2z_S_vrr-PRX*I_ESP_G2y2z_S_M1_vrr;
        Double I_ESP_Hxy3z_S_vrr = PAY*I_ESP_Gx3z_S_vrr-PRY*I_ESP_Gx3z_S_M1_vrr;
        Double I_ESP_Hx4z_S_vrr = PAX*I_ESP_G4z_S_vrr-PRX*I_ESP_G4z_S_M1_vrr;
        Double I_ESP_H5y_S_vrr = PAY*I_ESP_G4y_S_vrr-PRY*I_ESP_G4y_S_M1_vrr+4*oned2z*I_ESP_F3y_S_vrr-4*oned2z*I_ESP_F3y_S_M1_vrr;
        Double I_ESP_H4yz_S_vrr = PAZ*I_ESP_G4y_S_vrr-PRZ*I_ESP_G4y_S_M1_vrr;
        Double I_ESP_H3y2z_S_vrr = PAZ*I_ESP_G3yz_S_vrr-PRZ*I_ESP_G3yz_S_M1_vrr+oned2z*I_ESP_F3y_S_vrr-oned2z*I_ESP_F3y_S_M1_vrr;
        Double I_ESP_H2y3z_S_vrr = PAY*I_ESP_Gy3z_S_vrr-PRY*I_ESP_Gy3z_S_M1_vrr+oned2z*I_ESP_F3z_S_vrr-oned2z*I_ESP_F3z_S_M1_vrr;
        Double I_ESP_Hy4z_S_vrr = PAY*I_ESP_G4z_S_vrr-PRY*I_ESP_G4z_S_M1_vrr;
        Double I_ESP_H5z_S_vrr = PAZ*I_ESP_G4z_S_vrr-PRZ*I_ESP_G4z_S_M1_vrr+4*oned2z*I_ESP_F3z_S_vrr-4*oned2z*I_ESP_F3z_S_M1_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_H_S_aa
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_H_S_aa_coefs = alpha*alpha;
        I_ESP_H5x_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H5x_S_vrr;
        I_ESP_H4xy_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H4xy_S_vrr;
        I_ESP_H4xz_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H4xz_S_vrr;
        I_ESP_H3x2y_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H3x2y_S_vrr;
        I_ESP_H3xyz_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H3xyz_S_vrr;
        I_ESP_H3x2z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H3x2z_S_vrr;
        I_ESP_H2x3y_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H2x3y_S_vrr;
        I_ESP_H2x2yz_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H2x2yz_S_vrr;
        I_ESP_H2xy2z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H2xy2z_S_vrr;
        I_ESP_H2x3z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H2x3z_S_vrr;
        I_ESP_Hx4y_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_Hx4y_S_vrr;
        I_ESP_Hx3yz_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_Hx3yz_S_vrr;
        I_ESP_Hx2y2z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_Hx2y2z_S_vrr;
        I_ESP_Hxy3z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_Hxy3z_S_vrr;
        I_ESP_Hx4z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_Hx4z_S_vrr;
        I_ESP_H5y_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H5y_S_vrr;
        I_ESP_H4yz_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H4yz_S_vrr;
        I_ESP_H3y2z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H3y2z_S_vrr;
        I_ESP_H2y3z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H2y3z_S_vrr;
        I_ESP_Hy4z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_Hy4z_S_vrr;
        I_ESP_H5z_S_aa += SQ_ESP_H_S_aa_coefs*I_ESP_H5z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_F_S_a
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        Double SQ_ESP_F_S_a_coefs = alpha;
        I_ESP_F3x_S_a += SQ_ESP_F_S_a_coefs*I_ESP_F3x_S_vrr;
        I_ESP_F2xy_S_a += SQ_ESP_F_S_a_coefs*I_ESP_F2xy_S_vrr;
        I_ESP_F2xz_S_a += SQ_ESP_F_S_a_coefs*I_ESP_F2xz_S_vrr;
        I_ESP_Fx2y_S_a += SQ_ESP_F_S_a_coefs*I_ESP_Fx2y_S_vrr;
        I_ESP_Fxyz_S_a += SQ_ESP_F_S_a_coefs*I_ESP_Fxyz_S_vrr;
        I_ESP_Fx2z_S_a += SQ_ESP_F_S_a_coefs*I_ESP_Fx2z_S_vrr;
        I_ESP_F3y_S_a += SQ_ESP_F_S_a_coefs*I_ESP_F3y_S_vrr;
        I_ESP_F2yz_S_a += SQ_ESP_F_S_a_coefs*I_ESP_F2yz_S_vrr;
        I_ESP_Fy2z_S_a += SQ_ESP_F_S_a_coefs*I_ESP_Fy2z_S_vrr;
        I_ESP_F3z_S_a += SQ_ESP_F_S_a_coefs*I_ESP_F3z_S_vrr;

        /************************************************************
         * shell quartet name: SQ_ESP_P_S
         * doing contraction work for VRR part 
         * totally 0 integrals are omitted 
         ************************************************************/
        I_ESP_Px_S += I_ESP_Px_S_vrr;
        I_ESP_Py_S += I_ESP_Py_S_vrr;
        I_ESP_Pz_S += I_ESP_Pz_S_vrr;
    }

    /************************************************************
     * shell quartet name: SQ_ESP_F_S_dax_dax
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_S_aa
     * RHS shell quartet name: SQ_ESP_F_S_a
     * RHS shell quartet name: SQ_ESP_F_S_a
     * RHS shell quartet name: SQ_ESP_P_S
     ************************************************************/
    abcd[iGrid*60+0] = 4.0E0*I_ESP_H5x_S_aa-2.0E0*3*I_ESP_F3x_S_a-2.0E0*4*I_ESP_F3x_S_a+3*2*I_ESP_Px_S;
    abcd[iGrid*60+1] = 4.0E0*I_ESP_H4xy_S_aa-2.0E0*2*I_ESP_F2xy_S_a-2.0E0*3*I_ESP_F2xy_S_a+2*1*I_ESP_Py_S;
    abcd[iGrid*60+2] = 4.0E0*I_ESP_H4xz_S_aa-2.0E0*2*I_ESP_F2xz_S_a-2.0E0*3*I_ESP_F2xz_S_a+2*1*I_ESP_Pz_S;
    abcd[iGrid*60+3] = 4.0E0*I_ESP_H3x2y_S_aa-2.0E0*1*I_ESP_Fx2y_S_a-2.0E0*2*I_ESP_Fx2y_S_a;
    abcd[iGrid*60+4] = 4.0E0*I_ESP_H3xyz_S_aa-2.0E0*1*I_ESP_Fxyz_S_a-2.0E0*2*I_ESP_Fxyz_S_a;
    abcd[iGrid*60+5] = 4.0E0*I_ESP_H3x2z_S_aa-2.0E0*1*I_ESP_Fx2z_S_a-2.0E0*2*I_ESP_Fx2z_S_a;
    abcd[iGrid*60+6] = 4.0E0*I_ESP_H2x3y_S_aa-2.0E0*1*I_ESP_F3y_S_a;
    abcd[iGrid*60+7] = 4.0E0*I_ESP_H2x2yz_S_aa-2.0E0*1*I_ESP_F2yz_S_a;
    abcd[iGrid*60+8] = 4.0E0*I_ESP_H2xy2z_S_aa-2.0E0*1*I_ESP_Fy2z_S_a;
    abcd[iGrid*60+9] = 4.0E0*I_ESP_H2x3z_S_aa-2.0E0*1*I_ESP_F3z_S_a;

    /************************************************************
     * shell quartet name: SQ_ESP_F_S_dax_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_S_aa
     * RHS shell quartet name: SQ_ESP_F_S_a
     * RHS shell quartet name: SQ_ESP_F_S_a
     * RHS shell quartet name: SQ_ESP_P_S
     ************************************************************/
    abcd[iGrid*60+10] = 4.0E0*I_ESP_H4xy_S_aa-2.0E0*3*I_ESP_F2xy_S_a;
    abcd[iGrid*60+11] = 4.0E0*I_ESP_H3x2y_S_aa-2.0E0*1*I_ESP_F3x_S_a-2.0E0*2*I_ESP_Fx2y_S_a+2*1*I_ESP_Px_S;
    abcd[iGrid*60+12] = 4.0E0*I_ESP_H3xyz_S_aa-2.0E0*2*I_ESP_Fxyz_S_a;
    abcd[iGrid*60+13] = 4.0E0*I_ESP_H2x3y_S_aa-2.0E0*2*I_ESP_F2xy_S_a-2.0E0*1*I_ESP_F3y_S_a+2*I_ESP_Py_S;
    abcd[iGrid*60+14] = 4.0E0*I_ESP_H2x2yz_S_aa-2.0E0*1*I_ESP_F2xz_S_a-2.0E0*1*I_ESP_F2yz_S_a+1*I_ESP_Pz_S;
    abcd[iGrid*60+15] = 4.0E0*I_ESP_H2xy2z_S_aa-2.0E0*1*I_ESP_Fy2z_S_a;
    abcd[iGrid*60+16] = 4.0E0*I_ESP_Hx4y_S_aa-2.0E0*3*I_ESP_Fx2y_S_a;
    abcd[iGrid*60+17] = 4.0E0*I_ESP_Hx3yz_S_aa-2.0E0*2*I_ESP_Fxyz_S_a;
    abcd[iGrid*60+18] = 4.0E0*I_ESP_Hx2y2z_S_aa-2.0E0*1*I_ESP_Fx2z_S_a;
    abcd[iGrid*60+19] = 4.0E0*I_ESP_Hxy3z_S_aa;

    /************************************************************
     * shell quartet name: SQ_ESP_F_S_dax_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_S_aa
     * RHS shell quartet name: SQ_ESP_F_S_a
     * RHS shell quartet name: SQ_ESP_F_S_a
     * RHS shell quartet name: SQ_ESP_P_S
     ************************************************************/
    abcd[iGrid*60+20] = 4.0E0*I_ESP_H4xz_S_aa-2.0E0*3*I_ESP_F2xz_S_a;
    abcd[iGrid*60+21] = 4.0E0*I_ESP_H3xyz_S_aa-2.0E0*2*I_ESP_Fxyz_S_a;
    abcd[iGrid*60+22] = 4.0E0*I_ESP_H3x2z_S_aa-2.0E0*1*I_ESP_F3x_S_a-2.0E0*2*I_ESP_Fx2z_S_a+2*1*I_ESP_Px_S;
    abcd[iGrid*60+23] = 4.0E0*I_ESP_H2x2yz_S_aa-2.0E0*1*I_ESP_F2yz_S_a;
    abcd[iGrid*60+24] = 4.0E0*I_ESP_H2xy2z_S_aa-2.0E0*1*I_ESP_F2xy_S_a-2.0E0*1*I_ESP_Fy2z_S_a+1*I_ESP_Py_S;
    abcd[iGrid*60+25] = 4.0E0*I_ESP_H2x3z_S_aa-2.0E0*2*I_ESP_F2xz_S_a-2.0E0*1*I_ESP_F3z_S_a+2*I_ESP_Pz_S;
    abcd[iGrid*60+26] = 4.0E0*I_ESP_Hx3yz_S_aa;
    abcd[iGrid*60+27] = 4.0E0*I_ESP_Hx2y2z_S_aa-2.0E0*1*I_ESP_Fx2y_S_a;
    abcd[iGrid*60+28] = 4.0E0*I_ESP_Hxy3z_S_aa-2.0E0*2*I_ESP_Fxyz_S_a;
    abcd[iGrid*60+29] = 4.0E0*I_ESP_Hx4z_S_aa-2.0E0*3*I_ESP_Fx2z_S_a;

    /************************************************************
     * shell quartet name: SQ_ESP_F_S_day_day
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_S_aa
     * RHS shell quartet name: SQ_ESP_F_S_a
     * RHS shell quartet name: SQ_ESP_F_S_a
     * RHS shell quartet name: SQ_ESP_P_S
     ************************************************************/
    abcd[iGrid*60+30] = 4.0E0*I_ESP_H3x2y_S_aa-2.0E0*1*I_ESP_F3x_S_a;
    abcd[iGrid*60+31] = 4.0E0*I_ESP_H2x3y_S_aa-2.0E0*1*I_ESP_F2xy_S_a-2.0E0*2*I_ESP_F2xy_S_a;
    abcd[iGrid*60+32] = 4.0E0*I_ESP_H2x2yz_S_aa-2.0E0*1*I_ESP_F2xz_S_a;
    abcd[iGrid*60+33] = 4.0E0*I_ESP_Hx4y_S_aa-2.0E0*2*I_ESP_Fx2y_S_a-2.0E0*3*I_ESP_Fx2y_S_a+2*1*I_ESP_Px_S;
    abcd[iGrid*60+34] = 4.0E0*I_ESP_Hx3yz_S_aa-2.0E0*1*I_ESP_Fxyz_S_a-2.0E0*2*I_ESP_Fxyz_S_a;
    abcd[iGrid*60+35] = 4.0E0*I_ESP_Hx2y2z_S_aa-2.0E0*1*I_ESP_Fx2z_S_a;
    abcd[iGrid*60+36] = 4.0E0*I_ESP_H5y_S_aa-2.0E0*3*I_ESP_F3y_S_a-2.0E0*4*I_ESP_F3y_S_a+3*2*I_ESP_Py_S;
    abcd[iGrid*60+37] = 4.0E0*I_ESP_H4yz_S_aa-2.0E0*2*I_ESP_F2yz_S_a-2.0E0*3*I_ESP_F2yz_S_a+2*1*I_ESP_Pz_S;
    abcd[iGrid*60+38] = 4.0E0*I_ESP_H3y2z_S_aa-2.0E0*1*I_ESP_Fy2z_S_a-2.0E0*2*I_ESP_Fy2z_S_a;
    abcd[iGrid*60+39] = 4.0E0*I_ESP_H2y3z_S_aa-2.0E0*1*I_ESP_F3z_S_a;

    /************************************************************
     * shell quartet name: SQ_ESP_F_S_day_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_S_aa
     * RHS shell quartet name: SQ_ESP_F_S_a
     * RHS shell quartet name: SQ_ESP_F_S_a
     * RHS shell quartet name: SQ_ESP_P_S
     ************************************************************/
    abcd[iGrid*60+40] = 4.0E0*I_ESP_H3xyz_S_aa;
    abcd[iGrid*60+41] = 4.0E0*I_ESP_H2x2yz_S_aa-2.0E0*1*I_ESP_F2xz_S_a;
    abcd[iGrid*60+42] = 4.0E0*I_ESP_H2xy2z_S_aa-2.0E0*1*I_ESP_F2xy_S_a;
    abcd[iGrid*60+43] = 4.0E0*I_ESP_Hx3yz_S_aa-2.0E0*2*I_ESP_Fxyz_S_a;
    abcd[iGrid*60+44] = 4.0E0*I_ESP_Hx2y2z_S_aa-2.0E0*1*I_ESP_Fx2y_S_a-2.0E0*1*I_ESP_Fx2z_S_a+1*I_ESP_Px_S;
    abcd[iGrid*60+45] = 4.0E0*I_ESP_Hxy3z_S_aa-2.0E0*2*I_ESP_Fxyz_S_a;
    abcd[iGrid*60+46] = 4.0E0*I_ESP_H4yz_S_aa-2.0E0*3*I_ESP_F2yz_S_a;
    abcd[iGrid*60+47] = 4.0E0*I_ESP_H3y2z_S_aa-2.0E0*1*I_ESP_F3y_S_a-2.0E0*2*I_ESP_Fy2z_S_a+2*1*I_ESP_Py_S;
    abcd[iGrid*60+48] = 4.0E0*I_ESP_H2y3z_S_aa-2.0E0*2*I_ESP_F2yz_S_a-2.0E0*1*I_ESP_F3z_S_a+2*I_ESP_Pz_S;
    abcd[iGrid*60+49] = 4.0E0*I_ESP_Hy4z_S_aa-2.0E0*3*I_ESP_Fy2z_S_a;

    /************************************************************
     * shell quartet name: SQ_ESP_F_S_daz_daz
     * code section is: DERIV
     * totally 0 integrals are omitted 
     * RHS shell quartet name: SQ_ESP_H_S_aa
     * RHS shell quartet name: SQ_ESP_F_S_a
     * RHS shell quartet name: SQ_ESP_F_S_a
     * RHS shell quartet name: SQ_ESP_P_S
     ************************************************************/
    abcd[iGrid*60+50] = 4.0E0*I_ESP_H3x2z_S_aa-2.0E0*1*I_ESP_F3x_S_a;
    abcd[iGrid*60+51] = 4.0E0*I_ESP_H2xy2z_S_aa-2.0E0*1*I_ESP_F2xy_S_a;
    abcd[iGrid*60+52] = 4.0E0*I_ESP_H2x3z_S_aa-2.0E0*1*I_ESP_F2xz_S_a-2.0E0*2*I_ESP_F2xz_S_a;
    abcd[iGrid*60+53] = 4.0E0*I_ESP_Hx2y2z_S_aa-2.0E0*1*I_ESP_Fx2y_S_a;
    abcd[iGrid*60+54] = 4.0E0*I_ESP_Hxy3z_S_aa-2.0E0*1*I_ESP_Fxyz_S_a-2.0E0*2*I_ESP_Fxyz_S_a;
    abcd[iGrid*60+55] = 4.0E0*I_ESP_Hx4z_S_aa-2.0E0*2*I_ESP_Fx2z_S_a-2.0E0*3*I_ESP_Fx2z_S_a+2*1*I_ESP_Px_S;
    abcd[iGrid*60+56] = 4.0E0*I_ESP_H3y2z_S_aa-2.0E0*1*I_ESP_F3y_S_a;
    abcd[iGrid*60+57] = 4.0E0*I_ESP_H2y3z_S_aa-2.0E0*1*I_ESP_F2yz_S_a-2.0E0*2*I_ESP_F2yz_S_a;
    abcd[iGrid*60+58] = 4.0E0*I_ESP_Hy4z_S_aa-2.0E0*2*I_ESP_Fy2z_S_a-2.0E0*3*I_ESP_Fy2z_S_a+2*1*I_ESP_Py_S;
    abcd[iGrid*60+59] = 4.0E0*I_ESP_H5z_S_aa-2.0E0*3*I_ESP_F3z_S_a-2.0E0*4*I_ESP_F3z_S_a+3*2*I_ESP_Pz_S;
  }
}
