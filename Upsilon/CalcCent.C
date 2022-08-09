#include "Style.h"

#define pif pair<int,float>

bool cmp(const pif& a, const pif& b) {
	if (a.second == b.second) return a.first < b.first;
	return a.second > b.second;
}

void CalcCent(){

	gStyle->SetOptStat(0);
	gStyle->SetPalette(55);

	const int nset = 3;
	const bool bDRAW = true;
	const bool bWRITE = false;

	string system[nset] = {"pPb", "pO", "OO"};

	TFile *infileMult[nset];
	TH1D *hMult[3];
	TH1D *hCent[3];
	TH2D *hCentMult[3];

	for (int iset=0; iset<nset; iset++){

		infileMult[iset] = new TFile(Form("totalout_%s_fixed_v1.root",system[iset].c_str()),"read");

		hMult[iset] = new TH1D(Form("hMult_%s",system[iset].c_str()),Form("hMult_%s",system[iset].c_str()),1000,0,1000);
		hCent[iset] = new TH1D(Form("hCent_%s",system[iset].c_str()),Form("hCent_%s",system[iset].c_str()),1000,0,1000);

		hCentMult[iset] = new TH2D(Form("hCentMult_%s",system[iset].c_str()),Form("hCentMult_%s",system[iset].c_str()),100,0,100,100,0,100);

		map<int, float> map_mult;

		for (int ii=0; ii<1000; ii++){

			TH1D *_hMult = (TH1D*)infileMult[iset]->Get(Form("hMultDist_%d_0",ii));

			hMult[iset]->SetBinContent(ii+1, _hMult->GetMean());
			map_mult.insert(make_pair(ii, _hMult->GetMean()));

		}//ii

		vector<pif> vec_map(map_mult.begin(), map_mult.end());
		sort(vec_map.begin(), vec_map.end(), cmp);

		int count = 0;
		for (auto tmp_map: vec_map){

			int cent = int(count*0.1) + 1;

			hCent[iset]->SetBinContent(tmp_map.first+1, cent);

			//cout << cent << " " << tmp_map.second << " " << tmp_map.first << endl;
			count++;
		}

		//infileMult[iset]->Close();

		for (int ii=0; ii<1000; ii++){

			hCentMult[iset]->Fill(hCent[iset]->GetBinContent(ii+1), hMult[iset]->GetBinContent(ii+1));

		}

	}//iset

	if ( bDRAW ){

		TCanvas *c1 = new TCanvas("c1","c1",1.1*3*500,500);
		c1->Divide(3,1);

		for (int iset=0; iset<nset; iset++){
			c1->cd(iset+1);
			SetPadStyle(1);

			htmp = (TH1D*)hCentMult[iset];
			SetHistoStyle("Centrality","Multiplicity","",22,20);
			htmp->GetXaxis()->SetTitleOffset(1.1);
			htmp->GetYaxis()->SetTitleOffset(1.1);

			hCentMult[iset]->Draw("colz same");

		}//

	}//bDRAW


	if ( bWRITE ){
		TFile *outfile = new TFile("outfile_Centrality_pPb_pO_OO_8TeV.root","recreate");

		for (int iset=0; iset<nset; iset++){
			hMult[iset]->Write();
			hCent[iset]->Write();
			hCentMult[iset]->Write();
		}

		outfile->Close();
	}


}
