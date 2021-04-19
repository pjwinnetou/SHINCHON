#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

using namespace std;

void ScanV0(const char *set="pO8160GeV"){

	int npart;
	int part_pid[10000];
	float part_pt[10000], part_phi[10000], part_eta[10000];

	ifstream flist;

	char fname[300];
	sprintf(fname,"file_%s.lst",set);
	//sprintf(fname,"file.lst");

	flist.open(fname);

	//TH1D *hdNchdeta1 = new TH1D("hdNchdeta1","",200,0,200);
	//TH1D *hdNchdeta2 = new TH1D("hdNchdeta2","",200,0,200);

	TH1D *hdNchdeta1 = new TH1D("hdNchdeta1","",5000,0,5000);
	TH1D *hdNchdeta2 = new TH1D("hdNchdeta2","",5000,0,5000);

	while ( flist >> fname ){

		TFile *infile = new TFile(fname,"read");

		TTree *T = (TTree*)infile->Get("T");

		if ( !T ){
			infile->Close();
			delete infile;
			continue;
		}

		cout << "OPEN: " << fname << endl;

		T->SetBranchAddress("npart",&npart);
		T->SetBranchAddress("part_pid",part_pid);
		T->SetBranchAddress("part_pt",part_pt);
		T->SetBranchAddress("part_phi",part_phi);
		T->SetBranchAddress("part_eta",part_eta);

		int nentries = T->GetEntries();

		if ( nentries<1000 ){
			cout << "SKIP, # entries: " << nentries << endl;
			infile->Close();
			delete infile;
			continue;
		}

		double Nch_all = 0.0;

		for (int ien=0; ien<nentries; ien++){

			T->GetEntry(ien);

			double Nch = 0.0;

			for (int ipart=0; ipart<npart; ipart++){
				if ( abs(part_pid[ipart])==211 || abs(part_pid[ipart])==321 || abs(part_pid[ipart])==2212 ){
					Nch++;
				}
			}//ipart

			hdNchdeta1->Fill(Nch/2.0);

			Nch_all += Nch;

		}//ien

		hdNchdeta2->Fill(Nch_all/nentries/2.0);


		infile->Close();
		delete infile;

	}//

	flist.close();

	TFile *outfile = new TFile(Form("pre_outfile_SONIC_%s.root",set),"recreate");
	//TFile *outfile = new TFile("pre_outfile_SONIC.root","recreate");

	hdNchdeta1->Write();
	hdNchdeta2->Write();

	outfile->Close();

}
