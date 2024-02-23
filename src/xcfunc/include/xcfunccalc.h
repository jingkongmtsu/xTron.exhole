#ifndef XCFUNCCALC_H
#define XCFUNCCALC_H

#include<stdlib.h>
#include<vector>
#include "libgen.h"
#include "scalablevec.h"
#include "xvector.h"

namespace batchvar { class BatchVar; }
namespace xcvar { class XCVar; }
namespace xcfunc { class XCFunc; }

namespace xcfunccalc {
	//These indexes are 1 less than the one in fderiv1.inc.
	const Int ID_RA  = 0;
	const Int ID_RB  = 1;
	const Int ID_GAA = 2;
	const Int ID_GAB = 3;
	const Int ID_GBB = 4;
	const Int ID_TA  = 5;
	const Int ID_TB  = 6;
	const Int ID_LA  = 7;
	const Int ID_LB  = 8;
	const Int ID_EXA = 9;
	const Int ID_EXB = 10;
	const Int N_FUNC_DERIV_1 = 11;

	using namespace std;
	using namespace batchvar;
	using namespace xcvar;
	using namespace xcfunc;

	class XCFuncCalc;

	class XCFuncCalcResult1P
	{
		friend class XCFuncCalc;
		private:
			const XCFuncCalc* xcfc;
			XVector<double> f0;   //For functional values.
			XVector<double> f1;   //For energy profiling.
			XVector<double> d1F;

			//The following pointers are alias to d1F for fortran calls.
  		double* F0;
  		double* F1;
			double* DRhoA;
			double* DGAA;
			double* DTA;
			double* DLA;
			double* DUA;
			double* DRhoB;
			double* DGBB;
			double* DTB;
			double* DLB;
			double* DUB;
			double* DGAB;

			void assignPointers();  //For initialization and copy constructor, etc.

		public:
			XCFuncCalcResult1P();
			XCFuncCalcResult1P(const XCFuncCalc&, UInt nF, UInt nF1, UInt nD1F);
			XCFuncCalcResult1P(const XCFuncCalcResult1P& r1p);
			void addToBatch(double* F, double* F1, double* D1F, UInt ig);
			XCFuncCalcResult1P& operator *= (const Double& f);
			XCFuncCalcResult1P& operator += (const XCFuncCalcResult1P&);
			XCFuncCalcResult1P& operator=(const XCFuncCalcResult1P& r1p);
	};

	class XCFuncCalc {
		friend class XCFuncCalcResult1P;
		private:
			//Input variables.
			const XCFunc* xcfunc;
			const DoubleVec* hirWts;
			const Int* infor;           //For fortran calls such as kp14ecp1.
			vector<Int> d1Vars;         //Still use the fortran indexing for consistency.
			UInt ng;
			UInt nden;
			Int NA;
			Int NB;
			double thresh;
			const Double *rhoA;
			const Double *rhoB; 
			const Double *DRA;
			const Double *DRB;
			const Double *TA;
			const Double *TB;
			const Double *LA;
			const Double *LB;
			const Double *EXRA;
			const Double *EXRB;

			//The following is a bit waste of memory but makes coding a bit easier.
			XVector<double> GAA;
			XVector<double> GBB;
			XVector<double> GAB;

		public:
			XCFuncCalc(const Int* infor, UInt ngg, UInt nd, Int na, Int nb,
                 const BatchVar& bvar, Double thr,
                 const XCVar& xcvar, const XCFunc& xcf,
                 const DoubleVec& hwts);
			void calcB05NDOpp(Double* F, Double* D1F, UInt nF, UInt nD1F) const;
			void calcB05NDPar(Double* F, Double* D1F, UInt nF, UInt nD1F) const;
			void calcB05NDPar1(Double* F, Double* D1F, UInt nF, UInt nD1F) const;
			void calcKP14C(Double* F, Double* F1, Double* D1F, UInt nF, UInt nF1, 
			               UInt nD1F) const;
			XCFuncCalcResult1P b05NDPar1_1P(UInt nF, UInt nD1F, UInt i) const;
			XCFuncCalcResult1P b05NDPar_1P (UInt nF, UInt nD1F, UInt i) const;
			void calcB05NDParXX(Double* F, Double* D1F, UInt nF, UInt nD1F) const;
			void calcBR94CorrOpp(Double* F, Double* D1F, UInt nF, UInt nD1F) const;
			void calcBR94CorrOppX0(Double* F, Double* D1F, UInt nF, UInt nD1F) const;
			void calcBR94CorrOppX1(Double* F, Double* D1F, UInt nF, UInt nD1F) const;
			void calcBR94CorrOppX2(Double* F, Double* D1F, UInt nF, UInt nD1F) const;
			void calcBR94CorrOppX3(Double* F, Double* D1F, UInt nF, UInt nD1F) const;
			XCFuncCalcResult1P minSmooth2(const XCFuncCalcResult1P& r1,
                   const XCFuncCalcResult1P& r2) const;
	};
}

#endif
