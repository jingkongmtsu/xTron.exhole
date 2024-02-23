#include<iostream>
#include "xcfunccalc.h"
#include "functionallist.h"
#include "batchvar.h"
#include "xcfunc.h"
using namespace xcfunc;
using namespace batchvar;
using namespace xcfunccalc;


void XCFuncCalc::calcB05NDPar(Double* F, Double* D1F, UInt nF, UInt nD1F) const
{
	for (UInt i = 0; i < ng; i++)
	{
		b05NDPar_1P(nF, nD1F, i).addToBatch(F, 0, D1F, i);
	}
}



void XCFuncCalc::calcB05NDPar1(Double* F, Double* D1F, UInt nF, UInt nD1F) const
{
	for (UInt i = 0; i < ng; i++)
	{
		b05NDPar1_1P(nF, nD1F, i).addToBatch(F, 0, D1F, i);
	}
}

XCFuncCalcResult1P XCFuncCalc::b05NDPar_1P(UInt nF, UInt nD1F, UInt i) const
{
	//!!The following copying is rather inefficient.
	double pB05 = static_cast<double>(xcfunc->getBecke05PVal());
	double qB05 = static_cast<double>(xcfunc->getBecke05QVal());
	if ( xcfunc->getB05NDPar1Spin() == 0 )
	{
		XCFuncCalcResult1P ndpar(*this, nF, 0, nD1F);
		b05ndpar1p_(&NA,&NB,&thresh,&pB05,&qB05,&rhoA[i],&rhoB[i],&GAA[i],&GBB[i],
	            &TA[i],&TB[i],&LA[i],&LB[i],&EXRA[i],&EXRB[i],&(*hirWts)[i],
	            ndpar.F0, ndpar.DRhoA, ndpar.DGAA, ndpar.DTA, ndpar.DLA, ndpar.DUA,
	            ndpar.DRhoB, ndpar.DGBB, ndpar.DTB, ndpar.DLB, ndpar.DUB);
		return ndpar;
	}
	else
	{
		//This version of ndpar is for calling the version of b05ndpar1p for one
		//spin only (xx). For the beta spin, it switches the alpha and beta parts
		//just like in the formula.  I think it is correct but the numbers change
		//slightly at the 3rd scf iterations for bn/b05 job.  Quite weird.
		XCFuncCalcResult1P ndparA(*this, nF, 0, nD1F);
		b05ndpar1pxx_(&NA,&NB,&thresh,&pB05,&qB05,&rhoA[i],&rhoB[i],&GAA[i],&GBB[i],
		              &TA[i],&TB[i],&LA[i],&LB[i],&EXRA[i],&EXRB[i],&(*hirWts)[i],
  	              ndparA.F0, 
		              ndparA.DRhoA, ndparA.DGAA, ndparA.DTA, ndparA.DLA, ndparA.DUA,
  	              ndparA.DRhoB, ndparA.DGBB, ndparA.DTB, ndparA.DLB, ndparA.DUB);
		XCFuncCalcResult1P ndparB(*this, nF, 0, nD1F);
		//Now switch alpha and beta.
  	b05ndpar1pxx_(&NB,&NA,&thresh,&pB05,&qB05,&rhoB[i],&rhoA[i],&GBB[i],&GAA[i],
  	              &TB[i],&TA[i],&LB[i],&LA[i],&EXRB[i],&EXRA[i],&(*hirWts)[i],
  	              ndparB.F0, 
  	              ndparB.DRhoB, ndparB.DGBB, ndparB.DTB, ndparB.DLB, ndparB.DUB,
		              ndparB.DRhoA, ndparB.DGAA, ndparB.DTA, ndparB.DLA, ndparB.DUA);
		ndparA += ndparB;
		return ndparA;
	}
}

XCFuncCalcResult1P XCFuncCalc::b05NDPar1_1P(UInt nF, UInt nD1F, UInt i) const
{
	//!!The following copying is rather inefficient.
	Double cap = xcfunc->getKP14CNDPARCapVal();
	double pB05 = static_cast<double>(xcfunc->getBecke05PVal());

	XCFuncCalcResult1P ndpar = b05NDPar_1P(nF, nD1F, i); 

	XCFuncCalcResult1P ndopp(*this, nF, 0, nD1F);
	b05ndopp1p_(&thresh,&pB05,&rhoA[i],&rhoB[i],&GAA[i],&GBB[i],&TA[i],&TB[i],
	            &LA[i],&LB[i],&EXRA[i],&EXRB[i],&(*hirWts)[i],
	            ndopp.F0, ndopp.DRhoA, ndopp.DGAA, ndopp.DTA, ndopp.DLA, ndopp.DUA,
	            ndopp.DRhoB, ndopp.DGBB, ndopp.DTB, ndopp.DLB, ndopp.DUB);

	ndopp *= cap;

	return minSmooth2(ndpar, ndopp);
}

//minSmooth2 does dimensionless min(r1, r2) with z=(f2-f1)/sqrt(f1^2+f2^2).
//Note, it does not do anything to the energy profiling components. Those
//are presumably recalculated after the smoothening.
XCFuncCalcResult1P XCFuncCalc::minSmooth2(const XCFuncCalcResult1P& r1, 
                   const XCFuncCalcResult1P& r2) const
{
	Double P = 500.0E0;  //!!Parameter for smoothening.
	XCFuncCalcResult1P rv(*r1.xcfc, r1.f0.size(), r1.f1.size(), r1.d1F.size());  //return value.

	Double F1 = r1.f0[0];
	Double F2 = r2.f0[0];
	if (F2 > ZERO) F2=ZERO;    //This line is not smooth!

	Double T = F1*F1+F2*F2;
	Double Z(0), DZ1(0), DZ2(0);
	if (T > thresh)   //Nothing to do. IF (T > THRESH)
	{
		//Get the drv of Z wrt F1 and F2.
		Double T1 = sqrt(T); //sqrt is to keep z dimensionless.
		Z  = (F2-F1)/T1;
		Double T2 = pow(T1,3);
		DZ1= -F2*(F1+F2)/T2;
		DZ2= F1*(F1+F2)/T2;
	}
	//Drv of Z wrt functional variables.
	XVector<Double> DZ = DZ1*r1.d1F + DZ2*r2.d1F; //DZDRA=DZ1*DF1DRA + DZ2*DF2DRA..

	//Get the HPZ and DHZDZ (drv of HPZ wrt Z) value
	Double HPZ(ZERO), DHPZDZ(ZERO);
	{
		//one thing to note, that PZ below could be very large then
		//it may break the expression of HPZ therefore we use log to test it
		Double LIMIT_Z = -log(thresh);
		Double PZ    = P*Z;
		if (abs(PZ) < LIMIT_Z) 
   	{
			Double TP   = exp(P*Z);
			HPZ = ONE/(ONE+TP);
			DHPZDZ = -P*HPZ*HPZ*TP;
		}
		else if (PZ < ZERO) { HPZ = ONE;  }
		else if (PZ > ZERO) { HPZ = ZERO; }
	}

	//Collecting results. Note, only f0[0] is being calculated here. (Not sure what
  //other elements in f0 are for.)
	rv.f0[0] = (F1-F2)*HPZ + F2;
	rv.d1F = HPZ*(r1.d1F-r2.d1F) + (F1-F2)*DHPZDZ*DZ + r2.d1F;
	return rv;
}
