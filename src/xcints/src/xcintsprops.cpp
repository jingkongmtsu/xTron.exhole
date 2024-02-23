/**
 * CPP files corresponding to the xcenergyinfor.h
 * \author  Fenglai liu and Jing kong
 */
#include<iostream>
#include<cstdio>
#include<cmath>
#include "blas1.h"
#include "excep.h"
#include "xcfunc.h"
#include "xcintsinfor.h"
#include "batchgrid.h"
#include "batchvar.h"
#include "xcints.h"
#include "xcintsprops.h"
#include "denmtrx.h"
#include "functionallist.h"

using namespace excep;
using namespace blas;
using namespace xcfunc;
using namespace xcintsinfor;
using namespace batchgrid;
using namespace batchvar;
using namespace xcints;
using namespace xcintsprops;
using namespace denmtrx;
//using namespace xcenergyinfor;
using namespace std;

//The output of this function is to a file, for now.
BatchProp* XCHoleAvg::makeBatchProp(const XCInput& xcInp,
           const XCIterInput& xcItInp, const BatchIterBasicData& bibd) const
{
  //cout << "makebatchprop 1: thread::id = " << get_id() << endl;
  //Reform: this must be made a unique_ptr.
  BatchHoleAvg* bHoleAvg = new BatchHoleAvg(*this, xcInp, xcItInp, bibd);
  makeBatchHole(*bHoleAvg, xcInp, xcItInp, bibd);
  return bHoleAvg;
}

void XCHoleAvg::accumulate(const BatchProp* bProp)
{
  const BatchHoleAvg* bHoleAvg = (const BatchHoleAvg*)bProp;
	if (xcItInp->infor.param("XCHoleAvg.dumpBatchHoles") == "true" ) 
	{ bHoleAvg->dump(); }
  for (UInt ispin = 0; ispin < holeAvg.size(); ispin++)
  {
    for (UInt isval = 0; isval < holeAvg[ispin].size(); isval++)
    {
      holeAvg[ispin][isval] += bHoleAvg->holeAvg[ispin][isval];
    }
  }
}

BatchHoleAvg::BatchHoleAvg(const XCHoleAvg& hAvg, const XCInput& xi, 
	const XCIterInput& xii, const BatchIterBasicData& bibd) 
	: BatchProp(xi,xii,bibd), xcHoleAvg(hAvg), 
	  holes(hAvg.sValues.size(), vector<DoubleVec>(xi.denMtrx.getNSpin(), 
    DoubleVec(bibd.bg.getNGrids()))), 
  	holeAvg(xi.denMtrx.getNSpin(), DoubleVec(hAvg.holeAvg.size())) {}

XCHoleAvg::XCHoleAvg(UInt nSpin, UInt nSValue, const string& fo, 
	const XCIterInput& xcii) : XCProp(xcii), holeAvg(nSpin, DoubleVec(nSValue)), 
	sValues(nSValue), fout(fo)
{
	//Clean up the batch hole file since it is appended.
	FILE* myfile = fopen(fout.c_str(),"w"); fclose(myfile); //Clear file.
  Double radian = Double(1.50);
  for(UInt irad=0; irad<nSValue; irad++) {
    Double t = irad+1;
    Double u = nSValue-irad;
    Double tdu = t/u;
    sValues[irad] = radian*pow(tdu,TWO);
    //sValueWeightVec[irad] = TWO*pow(radian,THREE)*(nSValue+1)*
		//pow(tdu,FIVE)/pow(u,TWO);
  }
}

void BatchHoleAvg::dump() const
{
	FILE* myfile = fopen(xcHoleAvg.fout.c_str(),"a");
	UInt nGrid = bibd.bg.getNGrids();
	UInt nSValue = holes.size();
	UInt posRA = xcItInp.xcvar.getVarPos(xcvarinfor::RA);
	UInt posDRA = xcItInp.xcvar.getVarPos(xcvarinfor::DAX);
	UInt posTA = xcItInp.xcvar.getVarPos(xcvarinfor::TA);
	UInt posLA = xcItInp.xcvar.getVarPos(xcvarinfor::LA);
	//for (UInt ispin = 0; ispin < nSpin; ispin++)
  // Looping over ispin leads to problem, maybe mismatch of variable dimensions
	UInt ispin = 0;
	for (UInt ig = 0; ig < nGrid; ig++)
	{
		const Double* crd = bibd.bg.getGridCoord() + 3*ig;
		fprintf(myfile,"Coord: %12.6f  %12.6f  %12.6f\n", crd[0], crd[1], crd[2]);
		Double rA = bibd.bvar.getVar(posRA)[ig];
		const Double* dRA = bibd.bvar.getVar(posDRA);
		Double GAA = 0;
    for (UInt ic = 0; ic < 3; ic++) GAA += pow((dRA + ic*nGrid)[ig],2);
		Double tA = bibd.bvar.getVar(posTA)[ig];
		Double lA = bibd.bvar.getVar(posLA)[ig];
  	fprintf(myfile,"Rho  GAA  TA  LA\n");
		fprintf(myfile,"%20.14f %20.14f %20.14f %20.14f\n", rA, GAA, tA, lA);
  	for (UInt isval = 0; isval < nSValue; isval++)
		{
    	fprintf(myfile,"%20.14f  %20.14f\n",
				      xcHoleAvg.sValues[isval], holes[isval][ispin][ig]);
		}
	}
  fclose(myfile);
}

void XCHoleAvg::dump() const
{
	string foutAvg(fout + ".Avg"); 
	FILE* myfile = fopen(foutAvg.c_str(),"w");
  //for (UInt ispin = 0; ispin < holeAvg.size(); ispin++)
	UInt ispin = 0;
  {
  	for (UInt isval = 0; isval < holeAvg[ispin].size(); isval++)
    {
    	//fprintf(myfile,"xholeVal[%u] =  %12.30lf\n",isval, holeAvg[ispin][isval]);
    	//fprintf(myfile,"sValue = %12.30lf\n",sValues[isval]);
    	fprintf(myfile,"%20.14f  %20.14f\n",sValues[isval],holeAvg[ispin][isval]);
    }
  }
  fclose(myfile);
}

void BatchHoleAvg::integrate()
{
	UInt nSpin = xcInp.denMtrx.getNSpin();
	UInt nGrid = bibd.bg.getNGrids();
	UInt nSValue = holes.size();
	vector<UInt> posDen(nSpin);
	posDen[0] = xcItInp.xcvar.getVarPos(RA);
	posDen[nSpin-1] = xcItInp.xcvar.getVarPos(RB);
	for (UInt isval = 0; isval < nSValue; isval++)
	{
		//for (UInt ispin = 0; ispin < nSpin; ispin++)
  // Looping over ispin leads to problem, maybe mismatch of variable dimensions
		UInt ispin = 0;
		{
			holeAvg[ispin][isval] = 0;
			const Double* wts = bibd.bg.getGridWts();
			const Double* den = bibd.bvar.getVar(posDen[ispin]);
			for (UInt ig = 0; ig < nGrid; ig++)
			{
				holeAvg[ispin][isval] += wts[ig]*den[ig]*holes[isval][ispin][ig];
			}
		}
	}
}

MetaGGAHole::MetaGGAHole(UInt nSpin, UInt nSValue, const string& fo, 
	                       const XCIterInput& xcii)
  : XCHoleAvg(nSpin, nSValue, fo, xcii)
{
	//Change here for different MetaGGA holes.
	cpMetaGGAHole = &br89xhole_s;
}

HyperGGAHole::HyperGGAHole(UInt nSpin, UInt nSValue, const string& fo,
	                         const XCIterInput& xcii)
  : XCHoleAvg(nSpin, nSValue, fo, xcii)
{
	//Change here for different HyperGGA holes.
	cpHyperGGAHole = &brrelaxhole_s;
}

void MetaGGAHole::makeBatchHole(BatchHoleAvg& bHoleAvg, const XCInput& xcInp,
     const XCIterInput& xcItInp, const BatchIterBasicData& bibd) const
{
  const Double *rhoA(0), *rhoB(0), *DRA(0), *DRB(0), *TA(0), *TB(0), 
	             *LA(0), *LB(0), *EXRA(0), *EXRB(0);
  bibd.bvar.varForFort(rhoA, rhoB, DRA, DRB, TA, TB, LA, LB, EXRA, EXRB, 
	                     xcItInp.xcvar);
  Int ng_     = static_cast<Int>(bibd.bg.getNGrids());
  Int nDen_   = static_cast<Int>(xcItInp.infor.getNDen());

	for ( UInt is = 0; is < sValues.size(); is++ )
  {
		Double sValue = sValues[is];
		vector<DoubleVec>& xhole = bHoleAvg.holes[is];
		(*cpMetaGGAHole)(&xhole[0].front(), &sValue, rhoA, rhoB, DRA, DRB, LA, LB,
      TA, TB, &ng_, &nDen_);
	}
}

void HyperGGAHole::makeBatchHole(BatchHoleAvg& bHoleAvg, const XCInput& xcInp,
     const XCIterInput& xcItInp, const BatchIterBasicData& bibd) const
{
  const Double *rhoA(0), *rhoB(0), *DRA(0), *DRB(0), *TA(0), *TB(0), 
	             *LA(0), *LB(0), *EXRA(0), *EXRB(0);
  bibd.bvar.varForFort(rhoA, rhoB, DRA, DRB, TA, TB, LA, LB, EXRA, EXRB, 
	                     xcItInp.xcvar);
  Int ng_     = static_cast<Int>(bibd.bg.getNGrids());
  Int nDen_   = static_cast<Int>(xcItInp.infor.getNDen());
  double thresh_   = static_cast<Double>(xcItInp.infor.getFuncTol());
	for ( UInt is = 0; is < sValues.size(); is++ )
  {
		Double sValue = sValues[is];
		vector<DoubleVec>& xhole = bHoleAvg.holes[is];
		(*cpHyperGGAHole)(&xhole[0].front(), &sValue, &thresh_, rhoA, rhoB, 
		                  DRA, DRB, LA, LB, TA, TB, EXRA, EXRB, &ng_, 
		                  &nDen_);
	}
}
