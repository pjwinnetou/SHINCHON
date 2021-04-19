#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TProfile.h>

using namespace std;

void ScanV1(const char *set="pAu200GeV"){

	int npart;
	int part_pid[10000];
	float part_pt[10000], part_phi[10000], part_eta[10000];

	TFile *infile_pre = new TFile(Form("pre_outfile_SONIC_%s.root",set),"read");
	TH1D *hdNchdeta = (TH1D*)infile_pre->Get("hdNchdeta2");

	TProfile *hv2[5];
	TProfile *hv3[5];
	TProfile *hv4[5];

	for (int icent=0; icent<5; icent++){
		hv2[icent] = new TProfile(Form("hv2_%d",icent),"",15,0,3);
		hv3[icent] = new TProfile(Form("hv3_%d",icent),"",15,0,3);
		hv4[icent] = new TProfile(Form("hv4_%d",icent),"",15,0,3);
	}

	ifstream flist;

	char fname[300];
	sprintf(fname,"file_%s.lst",set);
	//sprintf(fname,"file.lst");

	flist.open(fname);

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

			Nch_all += Nch;

		}//ien

		Nch_all /= (nentries*2.0);

		int frac = hdNchdeta->Integral(hdNchdeta->FindBin(Nch_all), hdNchdeta->GetNbinsX());
		int cent = int(100*frac/hdNchdeta->Integral()) + 1;

		int cent_index = -1;

		if ( cent<=10 ){
			cent_index = 0;
		}else if ( cent<=20 ){
			cent_index = 1;
		}else if ( cent<=30 ){
			cent_index = 2;
		}else if ( cent<=40 ){
			cent_index = 3;
		}else if ( cent<=50 ){
			cent_index = 4;
		}

		if ( cent_index<0 ){
			infile->Close();
			delete infile;
			continue;
		}

		cout << Nch_all << " " << cent << endl;

		double cos_phi[3][15] = {0.}, sin_phi[3][15] = {0.};
		double Nch_pt[15] = {0.};

		for (int ien=0; ien<nentries; ien++){

			T->GetEntry(ien);

			for (int ipart=0; ipart<npart; ipart++){
				if ( part_pt[ipart]>3.0 ) continue;

				int ipt = int(part_pt[ipart]/0.2);

				if ( !(ipt>=0 && ipt<15) ) continue;

				if ( abs(part_pid[ipart])==211 || abs(part_pid[ipart])==321 || abs(part_pid[ipart])==2212 ){

					for (int ii=0; ii<3; ii++){
						cos_phi[ii][ipt] += cos((ii+2)*part_phi[ipart]);
						sin_phi[ii][ipt] += sin((ii+2)*part_phi[ipart]);
					}

					Nch_pt[ipt]++;
				}//chg particle

			}//ipart

		}//ien

		for (int ipt=0; ipt<15; ipt++){

			if ( Nch_pt[ipt]<1 ) continue;

			double v2 = sqrt(cos_phi[0][ipt]*cos_phi[0][ipt] + sin_phi[0][ipt]*sin_phi[0][ipt])/Nch_pt[ipt];
			double v3 = sqrt(cos_phi[1][ipt]*cos_phi[1][ipt] + sin_phi[1][ipt]*sin_phi[1][ipt])/Nch_pt[ipt];
			double v4 = sqrt(cos_phi[2][ipt]*cos_phi[2][ipt] + sin_phi[2][ipt]*sin_phi[2][ipt])/Nch_pt[ipt];

			hv2[cent_index]->Fill(hv2[cent_index]->GetBinCenter(ipt+1),v2);
			hv3[cent_index]->Fill(hv3[cent_index]->GetBinCenter(ipt+1),v3);
			hv4[cent_index]->Fill(hv4[cent_index]->GetBinCenter(ipt+1),v4);
		}

		infile->Close();
		delete infile;

	}//

	flist.close();

	TFile *outfile = new TFile(Form("outfile_vn_SONIC_%s.root",set),"recreate");
	//TFile *outfile = new TFile("outfile_vn_SONIC.root","recreate");

	for (int icent=0; icent<5; icent++){
		hv2[icent]->Write();
		hv3[icent]->Write();
		hv4[icent]->Write();
	}

	outfile->Close();

}
