#include<iostream>
#include "xcfunccalc.h"
#include "functionallist.h"
#include "batchvar.h"
#include "xcfunc.h"
using namespace xcfunc;
using namespace batchvar;
using namespace xcfunccalc;


void XCFuncCalc::calcBR94CorrOppX3(Double* F, Double* D1F, UInt nF, UInt nD1F) const
{
	//!!Referencing would be the way to go. -jk.
	int nden_ = static_cast<int>(nden);
	double pB05 = static_cast<double>(xcfunc->getBecke05PVal());
	const Double CAB2(-0.8E0);

	for (UInt i = 0; i < ng; i++)
	{
		//Skip if there is no opp spin electron around.
		if (rhoA[i] < thresh || rhoB[i] < thresh) continue;

		XCFuncCalcResult1P b05_f(*this, nF, 0, nD1F);
		if ( i == (i/194)*194 ) 
		{
		XCFuncCalcResult1P ubrx(*this, nF, 0, nD1F);
		double ubrxA, ubrxB;
		br89hole(&thresh,&rhoA[i], &GAA[i], &TA[i], &LA[i], &ubrxA,
              ubrx.DRhoA, ubrx.DGAA, ubrx.DTA, ubrx.DLA);
		br89hole(&thresh,&rhoB[i], &GBB[i], &TB[i], &LB[i], &ubrxB,
              ubrx.DRhoB, ubrx.DGBB, ubrx.DTB, ubrx.DLB);
		becke05_f(&thresh, &pB05, &(*hirWts)[i],
		          &rhoA[i], &GAA[i], &TA[i], &LA[i], &EXRA[i],
		          &rhoB[i], &GBB[i], &TB[i], &LB[i], &EXRB[i], b05_f.F0, 
		          b05_f.DRhoA, b05_f.DGAA, b05_f.DTA, b05_f.DLA, b05_f.DUA,
		          b05_f.DRhoB, b05_f.DGBB, b05_f.DTB, b05_f.DLB, b05_f.DUB);
		        printf("ubrxA%7.3f  EXRA%7.3f  ubrxB%7.3f  EXRB%7.3f  f%7.3f\n", 
		                ubrxA, EXRA[i]/rhoA[i], ubrxB, EXRB[i]/rhoB[i], b05_f.F0[0]);
		}
		continue;
		XCFuncCalcResult1P br(*this, nF, 0, nD1F);
		{
			XCFuncCalcResult1P ua(*this, nF, 0, nD1F);
			XCFuncCalcResult1P ub(*this, nF, 0, nD1F);
			Double roAm1 = 1E0/rhoA[i];
			Double roAm2 = roAm1*roAm1;
			Double roBm1 = 1E0/rhoB[i];
			Double roBm2 = roBm1*roBm1;
			ua.F0[0] = EXRA[i]*roAm1;
			ub.F0[0] = EXRB[i]*roBm1;

			Double fz, fz_ua, fz_ub;
			br94corroppfzab(&fz, &fz_ua, &fz_ub, &nden_, &thresh, ua.F0, ub.F0);

			*ua.DRhoA = -EXRA[i]*roAm2;
			*ua.DUA += roAm1;

			*ub.DRhoB += -EXRB[i]*roBm2;
			*ub.DUB += roBm1;

			br.F0[0] = CAB2*rhoA[i]*rhoB[i]*fz;
			br.d1F = CAB2*rhoA[i]*rhoB[i]*(fz_ua*ua.d1F + fz_ub*ub.d1F);
			*br.DRhoA += CAB2*rhoB[i]*fz;
			*br.DRhoB += CAB2*rhoA[i]*fz;
		}
		br.addToBatch(F, 0, D1F, i);
		XCFuncCalcResult1P rv(*this, nF, 0, nD1F);
		rv.f0[0] = (ONE - b05_f.f0[0])*br.f0[0];
		rv.d1F = br.d1F - br.f0[0]*b05_f.d1F - b05_f.f0[0]*br.d1F;
	}
}

void XCFuncCalc::calcBR94CorrOppX2(Double* F, Double* D1F, UInt nF, UInt nD1F) const
{
	//!!Referencing would be the way to go. -jk.
	int nden_ = static_cast<int>(nden);
	double pB05 = static_cast<double>(xcfunc->getBecke05PVal());
	const Double CAB2(-0.8E0);

	for (UInt i = 0; i < ng; i++)
	{
		//Skip if there is no opp spin electron around.
		if (rhoA[i] < thresh || rhoB[i] < thresh) continue;
		XCFuncCalcResult1P b05_f(*this, nF, 0, nD1F);
		becke05_f(&thresh, &pB05, &(*hirWts)[i],
		          &rhoA[i], &GAA[i], &TA[i], &LA[i], &EXRA[i],
		          &rhoB[i], &GBB[i], &TB[i], &LB[i], &EXRB[i], b05_f.F0, 
		          b05_f.DRhoA, b05_f.DGAA, b05_f.DTA, b05_f.DLA, b05_f.DUA,
		          b05_f.DRhoB, b05_f.DGBB, b05_f.DTB, b05_f.DLB, b05_f.DUB);
		XCFuncCalcResult1P br(*this, nF, 0, nD1F);
		{
			XCFuncCalcResult1P ua(*this, nF, 0, nD1F);
			XCFuncCalcResult1P ub(*this, nF, 0, nD1F);
			Double roAm1 = 1E0/rhoA[i];
			Double roAm2 = roAm1*roAm1;
			Double roBm1 = 1E0/rhoB[i];
			Double roBm2 = roBm1*roBm1;
			ua.F0[0] = EXRA[i]*roAm1;
			ub.F0[0] = EXRB[i]*roBm1;

			Double fz, fz_ua, fz_ub;
			br94corroppfzab(&fz, &fz_ua, &fz_ub, &nden_, &thresh, ua.F0, ub.F0);

			*ua.DRhoA = -EXRA[i]*roAm2;
			*ua.DUA += roAm1;

			*ub.DRhoB += -EXRB[i]*roBm2;
			*ub.DUB += roBm1;

			br.F0[0] = CAB2*rhoA[i]*rhoB[i]*fz;
			br.d1F = CAB2*rhoA[i]*rhoB[i]*(fz_ua*ua.d1F + fz_ub*ub.d1F);
			*br.DRhoA += CAB2*rhoB[i]*fz;
			*br.DRhoB += CAB2*rhoA[i]*fz;
		}
		br.addToBatch(F, 0, D1F, i);
		XCFuncCalcResult1P rv(*this, nF, 0, nD1F);
		rv.f0[0] = (ONE - b05_f.f0[0])*br.f0[0];
		rv.d1F = br.d1F - br.f0[0]*b05_f.d1F - b05_f.f0[0]*br.d1F;
	}
}

void XCFuncCalc::calcBR94CorrOppX1(Double* F, Double* D1F, UInt nF, UInt nD1F) const
{
	//!!Referencing would be the way to go. -jk.
	int nden_ = static_cast<int>(nden);
	double pB05 = static_cast<double>(xcfunc->getBecke05PVal());
	const Double CAB2(-0.8E0);

	for (UInt i = 0; i < ng; i++)
	{
		//Skip if there is no opp spin electron around.
		if (rhoA[i] < thresh || rhoB[i] < thresh) continue;
		XCFuncCalcResult1P b05_f(*this, nF, 0, nD1F);
		becke05_f(&thresh, &pB05, &(*hirWts)[i],
		          &rhoA[i], &GAA[i], &TA[i], &LA[i], &EXRA[i],
		          &rhoB[i], &GBB[i], &TB[i], &LB[i], &EXRB[i], b05_f.F0, 
		          b05_f.DRhoA, b05_f.DGAA, b05_f.DTA, b05_f.DLA, b05_f.DUA,
		          b05_f.DRhoB, b05_f.DGBB, b05_f.DTB, b05_f.DLB, b05_f.DUB);
		XCFuncCalcResult1P br(*this, nF, 0, nD1F);
		{
			XCFuncCalcResult1P ua(*this, nF, 0, nD1F);
			XCFuncCalcResult1P ub(*this, nF, 0, nD1F);
			Double roAm1 = 1E0/rhoA[i];
			Double roAm2 = roAm1*roAm1;
			Double roBm1 = 1E0/rhoB[i];
			Double roBm2 = roBm1*roBm1;
			ua.F0[0] = EXRA[i]*roAm1 + b05_f.F0[0]*EXRB[i]*roBm1;
			ub.F0[0] = EXRB[i]*roBm1 + b05_f.F0[0]*EXRA[i]*roAm1;

			Double fz, fz_ua, fz_ub;
			br94corroppfzab(&fz, &fz_ua, &fz_ub, &nden_, &thresh, ua.F0, ub.F0);

			ua.d1F = EXRB[i]*roBm1*b05_f.d1F;
			*ua.DRhoA += -EXRA[i]*roAm2;
			*ua.DUA += roAm1;
			*ua.DRhoB += -b05_f.F0[0]*EXRB[i]*roBm2;
			*ua.DUB += b05_f.f0[0]*roBm1;

			ub.d1F = EXRA[i]*roAm1*b05_f.d1F;
			*ub.DRhoB += -EXRB[i]*roBm2;
			*ub.DUB += roBm1;
			*ub.DRhoA += -b05_f.F0[0]*EXRB[i]*roAm2;
			*ub.DUA += b05_f.f0[0]*roAm1;

			br.F0[0] = CAB2*rhoA[i]*rhoB[i]*fz;
			br.d1F = CAB2*rhoA[i]*rhoB[i]*(fz_ua*ua.d1F + fz_ub*ub.d1F);
			*br.DRhoA += CAB2*rhoB[i]*fz;
			*br.DRhoB += CAB2*rhoA[i]*fz;
		}
		XCFuncCalcResult1P rv(*this, nF, 0, nD1F);
		rv.f0[0] = (ONE - b05_f.f0[0])*br.f0[0];
		rv.d1F = br.d1F - br.f0[0]*b05_f.d1F - b05_f.f0[0]*br.d1F;
		rv.addToBatch(F, 0, D1F, i);
	}
}

void XCFuncCalc::calcBR94CorrOppX0(Double* F, Double* D1F, UInt nF, UInt nD1F) const
{
	//!!Referencing would be the way to go. -jk.
	int nden_ = static_cast<int>(nden);
	double pB05 = static_cast<double>(xcfunc->getBecke05PVal());

	for (UInt i = 0; i < ng; i++)
	{
		XCFuncCalcResult1P br94COpp(*this, nF, 0, nD1F);
    br94corropp1p(&nden_,&NA,&NB,&thresh,
		          &rhoA[i],&rhoB[i],&GAA[i],&GBB[i],
		            &TA[i],&TB[i],&LA[i],&LB[i],
                  br94COpp.F0,
                  br94COpp.DRhoA, br94COpp.DGAA, br94COpp.DTA, br94COpp.DLA,
                  br94COpp.DRhoB, br94COpp.DGBB, br94COpp.DTB, br94COpp.DLB);

		XCFuncCalcResult1P b05_f(*this, nF, 0, nD1F);
		becke05_f(&thresh, &pB05, &(*hirWts)[i],
		          &rhoA[i], &GAA[i], &TA[i], &LA[i], &EXRA[i],
		          &rhoB[i], &GBB[i], &TB[i], &LB[i], &EXRB[i],
		          b05_f.F0, b05_f.DRhoA, b05_f.DGAA, b05_f.DTA, b05_f.DLA, b05_f.DUA,
		          b05_f.DRhoB, b05_f.DGBB, b05_f.DTB, b05_f.DLB, b05_f.DUB);

		XCFuncCalcResult1P rv(*this, nF, 0, nD1F);

		rv.f0[0] = (ONE - b05_f.f0[0])*br94COpp.f0[0];
		rv.d1F = br94COpp.d1F - br94COpp.f0[0]*b05_f.d1F - b05_f.f0[0]*br94COpp.d1F;
		rv.addToBatch(F, 0, D1F, i);
	}
}

void XCFuncCalc::calcBR94CorrOpp(Double* F, Double* D1F, UInt nF, UInt nD1F) const
{
	//!!Referencing would be the way to go. -jk.
	int nden_ = static_cast<int>(nden);
	for (UInt i = 0; i < ng; i++)
	{
		XCFuncCalcResult1P br94COpp(*this, nF, 0, nD1F);
    br94corropp1p(&nden_,&NA,&NB,&thresh,
		          &rhoA[i],&rhoB[i],&GAA[i],&GBB[i],
		            &TA[i],&TB[i],&LA[i],&LB[i],
                  br94COpp.F0,
                  br94COpp.DRhoA, br94COpp.DGAA, br94COpp.DTA, br94COpp.DLA,
                  br94COpp.DRhoB, br94COpp.DGBB, br94COpp.DTB, br94COpp.DLB);
		br94COpp.addToBatch(F, 0, D1F, i);
	}
}
