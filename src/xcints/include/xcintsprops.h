/**
 * \file    xcintsprops.h
 * \brief   this is used to analyze the XC energy component
 * \author  Jing Kong
 */
#ifndef XCINTSPROPS_H
#define XCINTSPROPS_H
#include <string>
#include <iostream>
#include "libgen.h"
#include "scalablevec.h"

namespace xcintsinfor { class XCIntsinfor; }
namespace xcints { 
	class TBB_XCIntsMatrix; 
	class BatchIterBasicData;
}
namespace xcfunc { class XCFunc; }
namespace denmtrx { class DenMtrx; }
namespace sigshellpairinfor { class SigMolShellPairInfor; };
namespace atomdenmtrx { class AtomDenMtrx; }
namespace denmtrx { class DenMtrx; }
namespace xcvar { class XCVar; }
namespace molegrids { class MoleGrids;}
namespace batchgrid { class BatchGrid; };
namespace batchbasis { class BatchBasis; };
namespace sigatombasis { class SigAtomBasis; };
namespace batchgrid { class BatchGrid; };
namespace batchbasis { class BatchBasis; };
namespace sigatombasis { class SigAtomBasis; };
namespace batchvar { class BatchVar; };

namespace xcintsprops {
	using namespace std;
	using namespace batchgrid;
	using namespace xcintsinfor;
	using namespace xcfunc;
	using namespace xcints;
	using namespace xcvar;
	using namespace atomdenmtrx;
	using namespace denmtrx;
	using namespace xcvar;
	using namespace molegrids;
	using namespace denmtrx;
	using namespace sigshellpairinfor;
	using namespace batchbasis;
	using namespace sigatombasis;
  using namespace batchgrid;
  using namespace batchbasis;
  using namespace sigatombasis;
  using namespace batchvar;

  //Reform. This class will turn into part of BatchIter class as pointers to be refreshed
  //at each iteration.
  class BatchIterBasicData {
    public:
      const BatchGrid& bg;
      const SigAtomBasis& sigList;
      const BatchBasis& bbasis;
      const BatchVar& bvar;

      BatchIterBasicData(const BatchGrid& bg0, const SigAtomBasis& sigList0,
        const BatchBasis& bbasis0, const BatchVar& bvar0) : bg(bg0), sigList(sigList0),
        bbasis(bbasis0), bvar(bvar0) {};
      ~BatchIterBasicData() {};
  };

  //Collection of input data needed for initializing the xc job.
  class XCInput {
    public:
      const XCFunc& xcfunc;                  ///< xc functional information
      const MolShell& ms;                    ///< input shell
      const DenMtrx& denMtrx;                ///< density matrix

      XCInput(const XCFunc& xcfunc0, const MolShell& ms0, const DenMtrx& den)
        : xcfunc(xcfunc0), ms(ms0), denMtrx(den) {}
  };

  //Collection of input data needed for looping with BatchIter loop, corresponding
	//to the data members of XCInts that are not copied from ctor args.
  class XCIterInput {
    public:
      const XCIntJobInfor& infor;            ///< information center for real calculation
      const SigMolShellPairInfor& spData;    ///< shell pair information data
      const MoleGrids& molGrids;             ///< molecular grids
      const MolShellSize& molShellSize;      ///< molecular shell size data
      const XCVar& xcvar;                    ///< xc variable information
      const DoubleVec& cartDen;              ///< Cartesian form density matrix in
                                             ///< significant shell pair form
      const AtomDenMtrx& atomDenMtrx;        ///< density matrix data for free atom


      XCIterInput(const SigMolShellPairInfor& spData0,
          const MoleGrids& molGrids0, const MolShellSize& molShellSize0, 
          const XCVar& xcvar0, const XCIntJobInfor& infor0, 
          const DoubleVec& cartDen0, const AtomDenMtrx& atomDenMtrx0)
        : spData(spData0), molGrids(molGrids0),molShellSize(molShellSize0),
          xcvar(xcvar0),infor(infor0),cartDen(cartDen0),
          atomDenMtrx(atomDenMtrx0) {};
  };

	//A container class for grid point wise data.
	class BatchProp {
    protected:
    	const XCInput& xcInp;
    	const XCIterInput& xcItInp;
    	const BatchIterBasicData& bibd;
		public:
			BatchProp(const XCInput& xi, const XCIterInput& xii,
			          const BatchIterBasicData& bi) 
			  : xcInp(xi), xcItInp(xii), bibd(bi) {}
			virtual void integrate() = 0;
			virtual ~BatchProp() {};
	};

	class XCProp {
		protected:
    	const XCIterInput* xcItInp;
		public:
			XCProp() : xcItInp(0) {}
			XCProp(const XCIterInput& xcii) : xcItInp(&xcii) {}
			virtual BatchProp* makeBatchProp(const XCInput&, const XCIterInput&,
				const BatchIterBasicData& bibd) const = 0;
			//virtual void deleteBatchProp() = 0;
			virtual void accumulate(const BatchProp*) = 0;
			//virtual void batchDump(const BatchHoleAvg&, const XCIterInput&, 
			//	             const BatchIterBasicData&) const = 0;
			virtual void dump() const = 0;
			virtual ~XCProp() {};
	};

	class XCHoleAvg;  //Forward declaration.
	class HFXHole;  //Forward declaration.
	class MetaGGAHole;  //Forward declaration.
	class HyperGGAHole;  //Forward declaration.

	class BatchHoleAvg: public BatchProp {
		friend class XCHoleAvg;
		friend class HFXHole;
		friend class MetaGGAHole; 
		friend class HyperGGAHole; 
		private:
			vector<vector<DoubleVec> > holes;   //holes[nsvalue][nspin][ngrid]
			vector<DoubleVec> holeAvg;					//holes[nspin][nsvalue]
			const XCHoleAvg& xcHoleAvg;         //for printing svalue only.
		public:
			BatchHoleAvg(const XCHoleAvg& hAvg, const XCInput& in, const XCIterInput& xii,
				const BatchIterBasicData& bibd);
			void integrate();
			void dump() const;
			~BatchHoleAvg() {}
	};

	class XCHoleAvg : public XCProp {
		friend class BatchHoleAvg;
		friend class TBB_XCIntsMatrix;
		protected:
			vector<DoubleVec> holeAvg;
			DoubleVec sValues;
			string fout;  //This name is for batch holes. fout.avg is for the avg hole.
		public:
			XCHoleAvg() : XCProp(),  holeAvg(0,DoubleVec(0)), sValues(0), fout(0) {}
			XCHoleAvg(UInt nSpin, UInt nSValue, const string& fo, const XCIterInput& xcii); 
			BatchProp* makeBatchProp(const XCInput&, const XCIterInput&,
				const BatchIterBasicData& bibd) const;
			virtual void makeBatchHole(BatchHoleAvg& bHoleAvg, const XCInput& xcInp,
                        const XCIterInput& xcItInp,
                        const BatchIterBasicData& bibd) const = 0;
			//void deleteBatchProp();
			void accumulate(const BatchProp*);
			//void batchDump(const BatchHoleAvg&, const XCIterInput&, 
			//	             const BatchIterBasicData&) const;
			void dump() const;
			virtual ~XCHoleAvg() {};
	};

	class HFXHole : public XCHoleAvg {
		typedef void (HFXHole::*CPtrHFXHole)(BatchHoleAvg& bHoleAvg,
                const XCInput& xcInp, const XCIterInput& xcItInp,
                const BatchIterBasicData& bibd) const;
		private:
			CPtrHFXHole cpHFXHole;
			void hfxHole(BatchHoleAvg& bHoleAvg, const XCInput& xcInp,
                   const XCIterInput& xcItInp,
                   const BatchIterBasicData& bibd) const;
		public:
			HFXHole() : XCHoleAvg(), cpHFXHole(0) {}
			HFXHole(UInt nSpin, UInt nSValue, const string& fo, const XCIterInput& xcii) 
			  : XCHoleAvg(nSpin, nSValue, fo, xcii), cpHFXHole(&HFXHole::hfxHole) {} 
			void makeBatchHole(BatchHoleAvg& bHoleAvg, const XCInput& xcInp,
                        const XCIterInput& xcItInp,
                        const BatchIterBasicData& bibd) const;
			//void deleteBatchProp();
			virtual ~HFXHole() {};
	};

	typedef void (*CPtrMetaGGAHole)(double* hx, const double* sValue, 
		               const double* RhoA, const double* RhoB, const double* DRA, 
		               const double* DRB, const double* LapA, const double* LapB, 
		               const double* TauA, const double* TauB, const Int* ND, 
		               const Int* NDEN);
	class MetaGGAHole : public XCHoleAvg {
		private:
			CPtrMetaGGAHole cpMetaGGAHole;
		public:
			MetaGGAHole() : XCHoleAvg(), cpMetaGGAHole(0) {}
			MetaGGAHole(UInt nSpin, UInt nSValue, const string& fo, const XCIterInput& xcii);
			void makeBatchHole(BatchHoleAvg& bHoleAvg, const XCInput&, const XCIterInput&,
				                 const BatchIterBasicData& bibd) const;
			virtual ~MetaGGAHole() {};
	};

	typedef void (*CPtrHyperGGAHole)(double* hx, const double* sValue, 
	                 const double* thresh,
		               const double* RhoA, const double* RhoB, const double* DRA, 
		               const double* DRB, const double* LapA, const double* LapB, 
		               const double* TauA, const double* TauB, const double* EXA,
	                 const double* EXB, const Int* ND, const Int* NDEN);
	class HyperGGAHole : public XCHoleAvg {
		private:
			CPtrHyperGGAHole cpHyperGGAHole;
		public:
			HyperGGAHole() : XCHoleAvg(), cpHyperGGAHole(0) {}
			HyperGGAHole(UInt nSpin, UInt nSValue, const string& fo, const XCIterInput& xcii);
			void makeBatchHole(BatchHoleAvg& bHoleAvg, const XCInput&, const XCIterInput&,
				                 const BatchIterBasicData& bibd) const;
			virtual ~HyperGGAHole() {};
	};

	//Reform.  Temporary. Future XCInts class.  My own tbb class for now. 
	// Should be put into xcints class in future. In the future, change the 
	// const& to real thing and intialize at the constructor for xcinput, and 
	// use unique_ptr for xciterinp to make_unique in the do call, perhaps.
	class XCIntsOper {
		private:
			const XCInput& xcInp;
			const XCIterInput& xcItInp;
			XCProp& xcProp;
		public:
			XCIntsOper(const XCInput& xcin, const XCIterInput& xcitin, XCProp& xcp)
			  : xcInp(xcin), xcItInp(xcitin), xcProp(xcp) {}
			void operator() (const blocked_range<UInt>& r) const;
	};
}
#endif
