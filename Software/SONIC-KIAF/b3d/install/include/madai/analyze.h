#ifndef __ANALYZE_H__
#define __ANALYZE_H__

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <sys/stat.h>
#include "coral.h"
#include "coralutils.h"
#include "H5Cpp.h"
#include "b3d.h"
#include "qualifier.h"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

using namespace std;

class CAnalyze{
public:
	CB3D *b3d;
	string qualifier,run_name;
	int npartsmax,neventsmax,nsample;
	parameterMap parmap;
	CompType *ptype;
	//CompType ptype(sizeof(CPartH5));
	CAnalyze(string &run_name_set);
	~CAnalyze();
	void SetQualifier(string qualname);
	string input_dataroot,output_dataroot;
	string h5_infilename,output_filename,vizfilename,detailsdirname;
	H5File *h5file,*vizfile;
	double phenix_pion_yield;
	double phenix_pion_yield_min,phenix_pion_yield_max;
	FILE *output;
	CPartH5 *partH5;
	CPart *part;
	bool CALCGARRAYS;
	bool STAR_ACCEPTANCE;
	int ReadDataH5(int ievent); // returns nparts for given event
	int ReadVizData(double tau);
	void CalcSpectra_STAR();
	double CalcSpectra_PHENIX();
	void CalcSpectra_PHENIXppbar();
	void Calc3DSpectra();
	void CalcV2();
	void CalcFlucVn();
	void CalcFlucVnGamma();
	void CalcFlucQn();
	void CalcFlucQnEbE();
	void CalcHBT(int num_threads);
	void GetHBTpars(CPartH5 *pH5,double &tau,double &rout,double &rside,double &rlong);
	void CalcGamma();
	void CalcGamma_BothSigns();
	void CalcBalance();
	double legendre(int ell,double x);
	void Consolidate(string run_name);
	CResList *reslist;
	
	bool bPrintSpectra, bPrintV2;
};

class CDistill{
public:
	void Distill(int nruns);
};



#endif
