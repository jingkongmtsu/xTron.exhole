//
// This is the template to set up fractional spin calculation
//
#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<fstream>
#include<string> 
#include<vector>
#include<cmath>
#include<boost/lexical_cast.hpp>
#include <boost/timer/timer.hpp>
#include "libgen.h"
#include "globalinfor.h"
#include "excep.h"
#include "textread.h"
#include "molecule.h"
#include "denmtrx.h"
#include "fracspininfor.h"
#include "shell.h"
#include "scf.h"
#include "parameterparsing.h"
using namespace molecule;
using namespace excep;
using namespace textread;
using namespace globalinfor;
using namespace shell;
using namespace denmtrx;
using namespace fracspininfor; 
using namespace scf;
using namespace std;
using namespace parameterparsing;

void frac1Job(const string& inf, const GlobalInfor& globalInfor)
{
	// let's see how many molecules in the input file
	// right now we only support two cases:
	// 1  only one molecule;
	// 2  two molecules, one is for generating guess 
	//    and the second is for real calculation
	UInt nMol = 0;
	// open the file
	ifstream in;
	in.open(inf.c_str(),ios::in);
	if (!in) {
		string infor = "failed to open the input file";
		Excep excep("frac","frac",EXCEPTION_FILE_MISSING,infor);
		handleExcep(excep);
	}

	// counting that how many sections we have
	string key = "molecule";
	nMol = 0;
	string line;
	WordConvert w;
	in.seekg(0,ios::beg);
	while(getline(in,line)) {
		LineParse l(line);
		if (l.isSec() && w.compare(l.getSecName(), key)) {
			nMol++;
		}
	}
	in.close();

	if (nMol <= 0) {
		string infor = "No %molecule section";
		Excep excep("frac","frac",EXCEPTION_FUNCTION_NOT_AVAILABLE,infor);
		handleExcep(excep);
	}

	// basic printtings
	cout << "----------------------------------------------" << endl;
	cout << "          FIRST JOB                           " << endl;
	cout << "----------------------------------------------" << endl;

	// this is real job
	Molecule molecule0(inf,1);

	// Get frac info, and a modified molecule. Note, the molecule
	// will be forced to be unrestricted.
	FracSpinInfor fracSpinInfor(globalInfor,molecule0);
	if (!fracSpinInfor.doFracSpinJob()) { cout << endl << "ERROR: No fracspin info." << endl; abort(); }
	const Molecule& newMol = fracSpinInfor.getMol();
	newMol.print();

	//should use newMol here but it probably does not make a difference.
	MolShell s0(inf,molecule0);
	cout << "number of basis sets " << s0.getNBas() << endl;

	SCF scf(globalInfor,newMol,s0);
	DenMtrx den(globalInfor,newMol,s0,s0,newMol.getNSpin()); 
	//Forming initial guess. 
	{
		MO mo(globalInfor,s0,newMol,newMol.getNSpin());
		mo.recover(scf.getSCFParam().getGuessMOPath());
		mo.formFracSpinMO(fracSpinInfor);
		den.formDenMtrx(mo);
	}
	scf.doSCF(globalInfor,newMol,s0,den);
}
