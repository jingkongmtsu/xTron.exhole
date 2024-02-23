#include<iostream>
#include "xcfunccalc.h"
#include "functionallist.h"
#include "batchvar.h"
#include "xcfunc.h"
using namespace xcfunc;
using namespace batchvar;
using namespace xcfunccalc;

void XCFuncCalc::calcKP14C(Double* F, Double* F1, Double* D1F, UInt nF, 
	                         UInt nF1, UInt nD1F) const
{
	//!!Referencing would be the way to go. -jk.
	double pB05 = static_cast<double>(xcfunc->getBecke05PVal());
	double qB05 = static_cast<double>(xcfunc->getBecke05QVal());
	double b = static_cast<double>(xcfunc->getKP14BVal());
	double alpha = static_cast<double>(xcfunc->getKP14AlphaVal());
	double c_ndpar = static_cast<double>(xcfunc->getKP14CNDPARVal());
	Int ndparMethod = xcfunc->getB05NDParMethod();
	int nden_ = static_cast<int>(nden);

	for (UInt i = 0; i < ng; i++)
	{
		//XCFuncCalcResult1P ndpar;
		XCFuncCalcResult1P ndpar(*this, nF, 0, nD1F);
		//!!decision at this deep level is very innefficient! Too many repetitions.
		// the function calls
		//should be set outside of batch loop, at the xcfunc level, perhaps.

		if ( ndparMethod == 0 )
		{
			ndpar = b05NDPar_1P(nF, nD1F, i);
		}
		else if ( ndparMethod == 1 ) 
		{
		  ndpar = b05NDPar1_1P(nF, nD1F, i);
		}
		else
		{
			printf("wrong choice of method for NDpar\n");
			exit(0);
		}

		//!!This is ia purely silly thing due the ad hoc calls of functions with ngrid points
		//for 1 grid point, and the fact that DRA,DRB have ngrid as the leading dim. -jk

		DoubleVec TMP_DRA(3), TMP_DRB(3);
		for (UInt crd = 0; crd < 3; crd++)
		{
			TMP_DRA[crd] = DRA[i+crd*ng];
			TMP_DRB[crd] = DRB[i+crd*ng];
		}
		XCFuncCalcResult1P ndkp(*this, nF, nF1, nD1F);
		kp14ec1p_(infor,&nden_,&NA,&NB,&b,&alpha,&c_ndpar,&thresh,&pB05,&qB05,
		          &rhoA[i],&rhoB[i],&GAA[i],&GBB[i],
		            &TA[i],&TB[i],&LA[i],&LB[i],&EXRA[i],&EXRB[i],&(*hirWts)[i],
		            &TMP_DRA[0],&TMP_DRB[0],ndpar.F0, 
		            ndpar.DRhoA, ndpar.DGAA, ndpar.DTA, ndpar.DLA, ndpar.DUA,
		            ndpar.DRhoB, ndpar.DGBB, ndpar.DTB, ndpar.DLB, ndpar.DUB,
		            &ndkp.f0[0],&ndkp.f1[0],&ndkp.f1[1],
		            ndkp.DRhoA, ndkp.DGAA, ndkp.DTA, ndkp.DLA, ndkp.DUA,
		            ndkp.DRhoB, ndkp.DGBB, ndkp.DTB, ndkp.DLB, ndkp.DUB);
		ndkp.addToBatch(F, F1, D1F, i);
	}
}
