#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>

#include <fstream>
#include <iostream>
#include <set>
#include <string>

using namespace std;

void Make_hist(const float length=20., const int ngrid=100, const int event=0, const int nA=197, const int nB=197, const int nQA=0, const int nQB=0){

	gStyle->SetOptStat(0);

	const bool bSAVE = true;
	const bool bDRAW = false;

	float xmax = length/2.;
	float ymax = length/2.;
	//int ngrid = 4096;

	float width = 2*xmax / (ngrid);

	TH2D *h2_time[1000];
	TH2D *h2_time_bin[1000];
	TH1D *h1_time = new TH1D("h1_time","",1000,0,1000);

	char buf[300];
	char fname[300];

	ifstream fdata;
	ifstream flist;
	flist.open("file.lst");

	float eta, x, y, u; 
	float dum;

	float max[1000] = {0.0};

	int count = 0;

	while ( flist >> fname ){

		cout << "open: " << fname << endl;

		string buf_fname = fname;
		string buf_sub = buf_fname.substr(buf_fname.find("r_")+2,5);

		float time = atof(buf_sub.c_str());
		cout << buf_sub << " " << time << endl;

		//continue;

		h1_time->Fill(count, time);

		h2_time[count] = new TH2D(Form("h2_evt%05d_t%05d",event,count),"",ngrid,-xmax,xmax,ngrid,-ymax,ymax);
		h2_time_bin[count] = new TH2D(Form("h2_bin_evt%05d_t%05d",event,count),"",ngrid,-xmax,xmax,ngrid,-ymax,ymax);

		fdata.open(fname);
		fdata.getline(buf,300);

		cout << buf << endl;

		while ( fdata >> x >> y >> u
				){

			//int xbin = h2_time[count]->GetXaxis()->FindBin(x+0.5*width);
			//int ybin = h2_time[count]->GetYaxis()->FindBin(y+0.5*width);
			int xbin = h2_time[count]->GetXaxis()->FindBin(x+0.1*width);
			int ybin = h2_time[count]->GetYaxis()->FindBin(y+0.1*width);

			h2_time[count]->SetBinContent(xbin, ybin, u);
			h2_time_bin[count]->Fill(x+0.1*width, y+0.1*width);

			if ( u>max[count] ) max[count] = u;
		}

		count++;
		fdata.close();

		if ( count>=1000 ) break;

	}//flist

	cout << "grid: " << sqrt(h2_time_bin[0]->Integral()) << ", " << h2_time[0]->Integral() << endl; 


	if ( bSAVE ){
		TFile *outfile = new TFile("outfile.root","recreate");

		for (int ii=0; ii<count; ii++){
			h2_time[ii]->Write();
		}
		h1_time->Write();
		outfile->Close();
	}


}
