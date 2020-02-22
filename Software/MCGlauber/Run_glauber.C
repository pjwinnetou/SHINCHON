//#include "runglauber_v3.2.C"

void Run_glauber(){

	gSystem->Load("libMathMore.so");
	gROOT->ProcessLine(".L runglauber_v3.2.C+g");

	bool bHIST = true;

	//runAndOutputLemonTree(100, 0.5, "Pb", "Pb", 68.0, 0.4, 4.0, 6.0, 1, "MCGlauber-PbPb-5.02TeV-b4-6fm.root");
	//runAndOutputLemonTree(100, 0.5, "Pb", "Pb", 68.0, 0.4, 6.0, 8.0, 1, "MCGlauber-PbPb-5.02TeV-b6-8fm.root");
	runAndOutputLemonTree(100, 0.5, "p", "Pb", 68.0, 0.4, 0.0, 2.0, 1, "MCGlauber-pPb-5.02TeV-b0-2fm.root");
	//runAndOutputLemonTree(50, 0.4, "He3", "Au", 42.0, 0.4, 0.0, 2.0, bHIST, "MCGlauber-He3Au-200GeV-b0-2fm.root");


}
