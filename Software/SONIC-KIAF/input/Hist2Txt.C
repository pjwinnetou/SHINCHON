const int nbins = 100; //TODO: change me
const double max_x = 10; // TODO: change me
const bool saveQuarkDists = false; // TODO: whether to save wounded nucleon distributions in individual root files.

// Output file from runglauber should be called "outFile_${system}_${firstEvent}_${lastEvent}.root".
const string name_system = "PbPb"; //TODO: e.g. "He3Au" or "pPb". For accessing and saving the right files/directories.
const int energy = 5020;
const int bmin = 0;
const int bmax = 1;
const int firstEvent = 0;
const int lastEvent = 100;
const int scaler = 5;

// TODO: select an energy scaling value
//const double e0 = (scaler*0.1+0.1)*0.00150022 * TMath::Power(140.*max_x/(5.*nbins), 4); // 200 GeV (RHIC) 
//const double e0 = (scaler*0.1+0.3)*0.00150022 * TMath::Power(140.*max_x/(5.*nbins), 4); // >5 TeV (LHC) 
//const double e0 = (scaler*0.1+0.6)*0.00150022 * TMath::Power(140.*max_x/(5.*nbins), 4); // >5 TeV (LHC) //b 12-13 fm 
const double e0 = (scaler*0.1+1.2)*0.00150022 * TMath::Power(140.*max_x/(5.*nbins), 4); // >5 TeV (LHC) //PbPb b 0-1 fm 

void Hist2Txt () {

	gSystem->Exec(Form("mkdir -p initedFiles_%s_%dGeV_b%d_%dfm_scaler%03d",name_system.c_str(),energy,bmin,bmax,scaler));

	// inFile contains the histogram with the wounded nucleon distribution for all events.
	TFile* inFile = new TFile(Form("MCGlauber-%s-%dGeV-b%d-%dfm.root",name_system.c_str(),energy,bmin,bmax), "READ");
	cout << "OPEN: " << inFile->GetName() << endl;
	TFile* outRootFile = NULL;

  double max = 0; // for finding the maximum 
	for (int eventNum = firstEvent; eventNum < lastEvent; eventNum++) {
		TH2D* initedHist = (TH2D*)inFile->Get(Form("inited_event%i",eventNum));
		//initedHist->Scale (e0);

		ofstream outFile;
		char outFileName[500];
		sprintf(outFileName, "./initedFiles_%s_%dGeV_b%d_%dfm_scaler%03d/event%d.dat",name_system.c_str(),energy,bmin,bmax,scaler,eventNum);
		outFile.open(outFileName);

    for (int xbin = 1; xbin <= nbins; xbin++) {
      for (int ybin = 1; ybin <= nbins; ybin++) {
				double content = initedHist->GetBinContent(xbin, ybin); 

				//if ( content<0.001 ){
				if ( 0 ){
					initedHist->SetBinContent(xbin, ybin, 0.001*e0);
				}else{
					initedHist->SetBinContent(xbin, ybin, content*e0);
				}

				content = initedHist->GetBinContent(xbin, ybin); 

        outFile << content << "\t";
        if (content > max) max = content;
      }
      outFile << "\n";
    }
    outFile.close();

    if (saveQuarkDists) {
     // outRootFile copies the histogram with the distribution of wounded nucleons for a single event.
     // This can be used by diffusion.C in generating the initial positions of heavy qqbar pairs.
			TFile* outRootFile = new TFile(Form("./initedFiles_%s_%dGeV_b%d_%dfm_scaler%03d/event%d.root",name_system.c_str(),energy,bmin,bmax,scaler,eventNum), "RECREATE");
     outRootFile->cd();
     initedHist->Write();
    }
  }
  cout << "Max val: " << max << " GeV" << endl; // prints out the maximum initial energy density found

  inFile->Close();
  delete inFile;
}
