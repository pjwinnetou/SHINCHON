#include "runglauber_v3.2.C"

void Run_glauber(){

	gSystem->Load("libMathMore.so");
	gROOT->ProcessLine(".L runglauber_v3.2.C+g");

	bool bHIST = true;

	//runAndOutputLemonTree(10, 0.5, "Pb", "Pb", 68.0, 0.4, 0.0, 2.0, 1, "MCGlauber-PbPb-5.02TeV-b0-2fm.root");
	runAndOutputLemonTree(50, 0.4, "He3", "Au", 42.0, 0.4, 0.0, 2.0, bHIST, "MCGlauber-He3Au-200GeV-b0-2fm.root");


}
