#include<iostream>
#include "xcfunccalc.h"
#include "functionallist.h"
#include "batchvar.h"
#include "xcfunc.h"

using namespace xcfunc;
using namespace batchvar;
using namespace xcfunccalc;

XCFuncCalc::XCFuncCalc(const Int* inf, UInt ngg, UInt nd, Int na, Int nb,
                       const BatchVar& bvar, Double thr,
                       const XCVar& xcvar, const XCFunc& xcf,
                       const DoubleVec& hwts)
: xcfunc(&xcf), hirWts(&hwts), infor(inf), d1Vars(N_FUNC_DERIV_1), ng(ngg), 
  nden(nd), GAA(ngg,0E0), GBB(ngg,0E0), GAB(ngg,0E0)
{
	bvar.varForFort(rhoA, rhoB, DRA, DRB, TA, TB, LA, LB, EXRA, EXRB, xcvar);
	init_func_deriv_1_(infor,&d1Vars.front());
	NA = static_cast<int>(na);
	NB = static_cast<int>(nb);
	thresh = static_cast<double>(thr);
	for (UInt i = 0; i < ng; i++)
	{
		for (UInt crd = 0; crd < 3; crd++)
		{
			GAA[i] += pow((DRA + crd*ng)[i],2);
			GBB[i] += pow((DRB + crd*ng)[i],2);
			GAB[i] += (DRA + crd*ng)[i]*(DRB + crd*ng)[i];
		}
	}
}

XCFuncCalcResult1P::XCFuncCalcResult1P(const XCFuncCalc& xf, UInt nF, UInt nF1, UInt nD1F)
: xcfc(&xf), f0(nF,0E0), f1(nF1,0E0),  d1F(nD1F,0E0)
{
	assignPointers();
}

XCFuncCalcResult1P::XCFuncCalcResult1P()
: xcfc(0), f0(0,0E0), f1(0,0E0),  d1F(0,0E0) {}

XCFuncCalcResult1P::XCFuncCalcResult1P(const XCFuncCalcResult1P& r1p)
: xcfc(r1p.xcfc), f0(r1p.f0), f1(r1p.f1), d1F(r1p.d1F)
{
  assignPointers();
}

XCFuncCalcResult1P& XCFuncCalcResult1P::operator=(const XCFuncCalcResult1P& r1p)
{
	xcfc = r1p.xcfc;
	f0 = r1p.f0;
	f1 = r1p.f1;
	d1F = r1p.d1F;
  assignPointers();
	return *this;
}

void XCFuncCalcResult1P::assignPointers()
{
	const vector<Int>& d1Vars = xcfc->d1Vars;
	F0    = &f0[0];
	F1		= &f1[0];
	DRhoA = &d1F[d1Vars[ID_RA]-1];
	DGAA  = &d1F[d1Vars[ID_GAA]-1];
	DTA   = &d1F[d1Vars[ID_TA]-1];
	DLA   = &d1F[d1Vars[ID_LA]-1];
	DUA   = &d1F[d1Vars[ID_EXA]-1];
	DRhoB = &d1F[d1Vars[ID_RB]-1];
	DGBB  = &d1F[d1Vars[ID_GBB]-1];
	DTB   = &d1F[d1Vars[ID_TB]-1];
	DLB   = &d1F[d1Vars[ID_LB]-1];
	DUB   = &d1F[d1Vars[ID_EXB]-1];
	DGAB  = &d1F[d1Vars[ID_GAB]-1];
}

XCFuncCalcResult1P& XCFuncCalcResult1P::operator *= (const Double& f)
{
	f0 *= f; f1 *=f; d1F *= f;
	return *this;
}

XCFuncCalcResult1P& XCFuncCalcResult1P::operator += (const XCFuncCalcResult1P& r)
{
	f0 += r.f0; f1 +=r.f1; d1F += r.d1F;
	return *this;
}

void XCFuncCalcResult1P::addToBatch(double* F, double* F1, double* D1F, UInt ig)
{
	UInt ng = xcfc->ng;
	for (UInt icmp = 0; icmp < f0.size(); icmp++) { F[ig+icmp*ng] = f0[icmp]; }
	for (UInt icmp = 0; icmp < f1.size(); icmp++) { F1[ig+icmp*ng] = f1[icmp]; }
	for (UInt icmp = 0; icmp < d1F.size(); icmp++) { D1F[ig+icmp*ng] = d1F[icmp]; }
}

void XCFuncCalc::calcB05NDOpp(double* F, double* D1F, UInt nF, UInt nD1F) const
{
	//The reason not to make this a data member is that this will
	//show where the data come from.
	double pB05 = static_cast<double>(xcfunc->getBecke05PVal());
	for (UInt i = 0; i < ng; i++)
	{
		XCFuncCalcResult1P r1p(*this, nF, 0, nD1F);
		b05ndopp1p_(&thresh,&pB05,&rhoA[i],&rhoB[i],&GAA[i],&GBB[i],&TA[i],&TB[i], 
		            &LA[i],&LB[i],&EXRA[i],&EXRB[i],&(*hirWts)[i],
		            r1p.F0, r1p.DRhoA, r1p.DGAA, r1p.DTA, r1p.DLA, r1p.DUA, 
		            r1p.DRhoB, r1p.DGBB, r1p.DTB, r1p.DLB, r1p.DUB);
		r1p.addToBatch(F, 0, D1F, i);
	}
}


