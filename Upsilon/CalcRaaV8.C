#include "Style.h"

#include <TRandom.h>
#include <TRandom3.h>
#include <TH1.h>
#include <TFile.h>
#include <TProfile.h>
#include <TH2.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TProfile2D.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TH3.h>
#include <TVector2.h>

#include <iostream>
#include <fstream>

using namespace std;

Double_t fTsallis1S_v2(Double_t *x, Double_t *fpar);

void CalcRaaV8(const string system = "pPb"){

	gStyle->SetOptStat(0);
	gStyle->SetPalette(55);

	gRandom = new TRandom3(0);

	const bool bDRAW_Config = false; 
	const bool bDRAW_Result = false; 
	const bool bSAVE = true; 
	const bool b2D = false;

	const bool bPreQGP = false; //should be off
	const bool bPreRES = false;
	const bool bPTONLY = true;

	const int run_i = 0;
	const int run_f = 1000;
	const int nrun = 1000;
	const float const_hbarc = 197.5; //MeV fm
	const float const_mY1S = 9.46; //GeV
	const float const_mY2S = 10.02; //GeV
	const float const_mY3S = 10.36; //GeV

	//arXiv:1605.03561
	//arXiv:2007.03939
	//const int const_Opt = 0;
	//const float const_TmaxY1S = 600.0; //MeV
	//const float const_TmaxY2S = 230.0; //MeV
	//const float const_TmaxY3S = 170.0; //MeV
	//const float const_tau0Y1S = 0.2; //fm/c
	//const float const_tau0Y2S = 0.4; //fm/c
	//const float const_tau0Y3S = 0.6; //fm/c
	//const float const_Tf = 160.0; //arXiv:2007.03939
	//const float const_Tf = 192.0; //arXiv:1605.03561

	//arXiv:1706.08670
	const int const_Opt = 1;
	const float const_tau0Y1S = 0.5; //fm/c
	const float const_tau0Y2S = 1.0; //fm/c
	const float const_tau0Y3S = 1.5; //fm/c
	const float const_TmaxY1S = 600.0; //MeV
	const float const_TmaxY2S = 240.0; //MeV
	const float const_TmaxY3S = 190.0; //MeV
	const float const_Tf = 170.0;
	
	//nSAMP*nColl Ups per event
	const int nSAMP = 10; //times Ncoll
	//const int nSAMP = 0; //times Ncoll

	const float npartmax = 50;
	const float nmultmax = 100;

	ifstream fdata;

	char buf[500];
	vector<float> T[30];
	vector<float> Gdiss[30];
	vector<float> Cregen[30];

	float f_tmp[20];

	//thermal width Y(1S)
	fdata.open("Gdiss0.dat");
	fdata.getline(buf,500);

	while ( fdata 
			>> f_tmp[0] >> f_tmp[1] >> f_tmp[2] >> f_tmp[3]
			>> f_tmp[4] >> f_tmp[5] >> f_tmp[6] >> f_tmp[7]
			>> f_tmp[8] >> f_tmp[9] >> f_tmp[10] >> f_tmp[11]
			){
		for (int ipt=0; ipt<11; ipt++){
			T[ipt].push_back(f_tmp[0]);
			Gdiss[ipt].push_back(f_tmp[ipt+1]);
		}
	}

	fdata.close();

	fdata.open("Gdiss1.dat");
	fdata.getline(buf,500);

	while ( fdata 
			>> f_tmp[0] >> f_tmp[1] >> f_tmp[2] >> f_tmp[3]
			>> f_tmp[4] >> f_tmp[5] >> f_tmp[6] >> f_tmp[7]
			>> f_tmp[8] >> f_tmp[9] >> f_tmp[10]
			){
		for (int ipt=0; ipt<10; ipt++){
			T[ipt+11].push_back(f_tmp[0]);
			Gdiss[ipt+11].push_back(f_tmp[ipt+1]);
		}
	}

	fdata.close();

	TGraphErrors *gGdiss[21];
	TF1 *fGdiss[21];

	for (int ipt=0; ipt<21; ipt++){
		gGdiss[ipt] = new TGraphErrors;
		for (unsigned int iT=0; iT<T[ipt].size(); iT++){
			gGdiss[ipt]->SetPoint(iT, T[ipt][iT], Gdiss[ipt][iT]);
		}

		fGdiss[ipt] = new TF1(Form("fGdiss_%d",ipt),"pol5",0,520);
		gGdiss[ipt]->Fit(fGdiss[ipt],"R0Q");
	}

	//thermal width Y(2S)
	for (int ii=0; ii<30; ii++){
		T[ii].clear();
		Gdiss[ii].clear();
	}
	fdata.open("diss_2s.dat");
	fdata.getline(buf,500);

	while ( fdata 
			>> f_tmp[0] >> f_tmp[1] >> f_tmp[2] >> f_tmp[3]
			>> f_tmp[4] >> f_tmp[5] >> f_tmp[6] >> f_tmp[7]
			>> f_tmp[8] >> f_tmp[9] >> f_tmp[10] >> f_tmp[11]
			){
		for (int ipt=0; ipt<11; ipt++){
			if ( f_tmp[ipt+1]>0 ){
				T[ipt].push_back(f_tmp[0]);
				Gdiss[ipt].push_back(f_tmp[ipt+1]);
			}//
		}//ipt
	}

	TGraphErrors *gGdiss2S[11];
	TF1 *fGdiss2S[11];

	for (int ipt=0; ipt<11; ipt++){
		gGdiss2S[ipt] = new TGraphErrors;
		for (unsigned int iT=0; iT<T[ipt].size(); iT++){
			gGdiss2S[ipt]->SetPoint(iT, T[ipt][iT], Gdiss[ipt][iT]);
		}

		fGdiss2S[ipt] = new TF1(Form("fGdiss2S_%d",ipt),"pol2",0,230);
		gGdiss2S[ipt]->Fit(fGdiss2S[ipt],"R0Q");
	}

	fdata.close();

	//thermal width Y(3S)
	for (int ii=0; ii<30; ii++){
		T[ii].clear();
		Gdiss[ii].clear();
	}
	fdata.open("diss_3s.dat");
	fdata.getline(buf,500);

	while ( fdata 
			>> f_tmp[0] >> f_tmp[1] >> f_tmp[2]
			){
		for (int ipt=0; ipt<2; ipt++){
			if ( f_tmp[ipt+1]>0 ){
				T[ipt].push_back(f_tmp[0]);
				Gdiss[ipt].push_back(f_tmp[ipt+1]);
			}//
		}//ipt
	}

	TGraphErrors *gGdiss3S[11];
	TF1 *fGdiss3S[11];

	for (int ipt=0; ipt<2; ipt++){
		gGdiss3S[ipt] = new TGraphErrors;
		for (unsigned int iT=0; iT<T[ipt].size(); iT++){
			gGdiss3S[ipt]->SetPoint(iT, T[ipt][iT], Gdiss[ipt][iT]);
		}

		fGdiss3S[ipt] = new TF1(Form("fGdiss3S_%d",ipt),"pol2",0,180);
		gGdiss3S[ipt]->Fit(fGdiss3S[ipt],"R0Q");
	}

	fdata.close();

	//regeneration term
	for (int ii=0; ii<30; ii++){
		T[ii].clear();
	}
	fdata.open("regen0.dat");
	fdata.getline(buf,500);

	while ( fdata 
			>> f_tmp[0] >> f_tmp[1] >> f_tmp[2] >> f_tmp[3]
			>> f_tmp[4] >> f_tmp[5] >> f_tmp[6] >> f_tmp[7]
			>> f_tmp[8] >> f_tmp[9]
			){
		for (int ipt=0; ipt<9; ipt++){
			T[ipt].push_back(f_tmp[0]);
			Cregen[ipt].push_back(f_tmp[ipt+1]);
		}
	}

	fdata.close();

	fdata.open("regen1.dat");
	fdata.getline(buf,500);

	while ( fdata 
			>> f_tmp[0] >> f_tmp[1] >> f_tmp[2] >> f_tmp[3]
			>> f_tmp[4] >> f_tmp[5] >> f_tmp[6] >> f_tmp[7]
			>> f_tmp[8] >> f_tmp[9]
			){
		for (int ipt=0; ipt<9; ipt++){
			T[9+ipt].push_back(f_tmp[0]);
			Cregen[9+ipt].push_back(f_tmp[ipt+1]);
		}
	}

	fdata.close();

	fdata.open("regen2.dat");
	fdata.getline(buf,500);

	while ( fdata 
			>> f_tmp[0] >> f_tmp[1] >> f_tmp[2] >> f_tmp[3]
			>> f_tmp[4] >> f_tmp[5] >> f_tmp[6] >> f_tmp[7]
			>> f_tmp[8]
			){
		for (int ipt=0; ipt<8; ipt++){
			T[18+ipt].push_back(f_tmp[0]);
			Cregen[18+ipt].push_back(f_tmp[ipt+1]);
		}
	}

	fdata.close();

	TGraphErrors *gCregen[21];
	TF1 *fCregen[21];

	for (int ipt=0; ipt<21; ipt++){
		gCregen[ipt] = new TGraphErrors;
		for (unsigned int iT=0; iT<T[ipt].size(); iT++){
			gCregen[ipt]->SetPoint(iT, T[ipt][iT], Cregen[ipt][iT]);
		}

		fCregen[ipt] = new TF1(Form("fCregen_%d",ipt),"pol5",0,520);
		gCregen[ipt]->Fit(fCregen[ipt],"R0Q");
	}

	//dN/dt vs pT for various T
	TH1D *hRregen[45];
	
	for (int iT=0; iT<45; iT++){
		float Temp = 170 + (10*iT);
		hRregen[iT] = new TH1D(Form("hRregen_T%03dMeV",int(Temp)),"",21,-0.5,20.5);

		for (int ipt=0; ipt<21; ipt++){
			hRregen[iT]->SetBinContent(ipt+1, fCregen[ipt]->Eval(Temp)/const_hbarc);
		}//ipt
	}//iT

	TCanvas *c0;
	TCanvas *c0_2S;
	TCanvas *c0_3S;
	TCanvas *c0_;

	if ( bDRAW_Config ){
		c0 = new TCanvas("c0","c0",1.2*500,500);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(160,0,600,450);
		SetHistoStyle("T [MeV]","#Gamma_{diss} [MeV]");

		for (int ipt=0; ipt<21; ipt++){
			gGdiss[ipt]->SetLineWidth(3);
			gGdiss[ipt]->SetLineColorAlpha(kGray,0.3);
			gGdiss[ipt]->SetMarkerColor(kGray);
			gGdiss[ipt]->SetMarkerStyle(24);
			if ( ipt==0 ){
				gGdiss[ipt]->SetLineColor(2);
				gGdiss[ipt]->SetMarkerColor(2);
			}else if ( ipt==5 ){
				gGdiss[ipt]->SetLineColor(kGreen+2);
				gGdiss[ipt]->SetMarkerColor(kGreen+2);
			}else if ( ipt==10 ){
				gGdiss[ipt]->SetLineColor(kBlue);
				gGdiss[ipt]->SetMarkerColor(kBlue);
			}else if ( ipt==15 ){
				gGdiss[ipt]->SetLineColor(kMagenta);
				gGdiss[ipt]->SetMarkerColor(kMagenta);
			}else if ( ipt==20 ){
				gGdiss[ipt]->SetLineColor(1);
				gGdiss[ipt]->SetMarkerColor(1);
			}
			gGdiss[ipt]->Draw("P");

			fGdiss[ipt]->SetRange(160,600);
			fGdiss[ipt]->SetLineColor(1);
			fGdiss[ipt]->SetLineWidth(1);
			fGdiss[ipt]->SetLineStyle(7);
			fGdiss[ipt]->Draw("same");
		}

		TLegend *leg = new TLegend(0.15,0.65,0.4,0.95);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(43);
		leg->SetTextSize(18);
		leg->AddEntry("","#Upsilon(1S)","h");
		leg->AddEntry(gGdiss[0],"p_{T}=0 GeV/c","p");
		leg->AddEntry(gGdiss[5],"p_{T}=5 GeV/c","p");
		leg->AddEntry(gGdiss[10],"p_{T}=10 GeV/c","p");
		leg->AddEntry(gGdiss[15],"p_{T}=15 GeV/c","p");
		leg->AddEntry(gGdiss[20],"p_{T}=20 GeV/c","p");
		leg->Draw();

	}

	if ( bDRAW_Config ){
		c0_2S = new TCanvas("c0_2S","c0_2S",1.2*500,500);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(160,0,300,450);
		SetHistoStyle("T [MeV]","#Gamma_{diss} [MeV]");

		for (int ipt=0; ipt<11; ipt++){
			gGdiss2S[ipt]->SetLineWidth(2);
			gGdiss2S[ipt]->SetLineColorAlpha(kGray,0.3);
			gGdiss2S[ipt]->SetMarkerColor(kGray);
			gGdiss2S[ipt]->SetMarkerStyle(24);
			if ( ipt==0 ){
				gGdiss2S[ipt]->SetLineColor(2);
				gGdiss2S[ipt]->SetMarkerColor(2);
			}else if ( ipt==2 ){
				gGdiss2S[ipt]->SetLineColor(kGreen+2);
				gGdiss2S[ipt]->SetMarkerColor(kGreen+2);
			}else if ( ipt==4 ){
				gGdiss2S[ipt]->SetLineColor(kBlue);
				gGdiss2S[ipt]->SetMarkerColor(kBlue);
			}else if ( ipt==6 ){
				gGdiss2S[ipt]->SetLineColor(kMagenta);
				gGdiss2S[ipt]->SetMarkerColor(kMagenta);
			}else if ( ipt==8 ){
				gGdiss2S[ipt]->SetLineColor(1);
				gGdiss2S[ipt]->SetMarkerColor(1);
			}
			gGdiss2S[ipt]->Draw("P");

			fGdiss2S[ipt]->SetRange(160,300);
			fGdiss2S[ipt]->SetLineColor(1);
			fGdiss2S[ipt]->SetLineWidth(1);
			fGdiss2S[ipt]->SetLineStyle(7);
			fGdiss2S[ipt]->Draw("same");
		}

		TLegend *leg = new TLegend(0.15,0.65,0.4,0.95);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(43);
		leg->SetTextSize(18);
		leg->AddEntry("","#Upsilon(2S)","h");
		leg->AddEntry(gGdiss2S[0],"p_{T}=0 GeV/c","p");
		leg->AddEntry(gGdiss2S[2],"p_{T}=4 GeV/c","p");
		leg->AddEntry(gGdiss2S[4],"p_{T}=8 GeV/c","p");
		leg->AddEntry(gGdiss2S[6],"p_{T}=12 GeV/c","p");
		leg->AddEntry(gGdiss2S[8],"p_{T}=16 GeV/c","p");
		leg->Draw();
	}

	if ( bDRAW_Config ){
		c0_3S = new TCanvas("c0_3S","c0_3S",1.2*500,500);
		SetPadStyle();

		htmp = (TH1D*)gPad->DrawFrame(160,0,185,450);
		SetHistoStyle("T [MeV]","#Gamma_{diss} [MeV]");

		for (int ipt=0; ipt<2; ipt++){
			gGdiss3S[ipt]->SetLineWidth(2);
			gGdiss3S[ipt]->SetLineColorAlpha(kGray,0.3);
			gGdiss3S[ipt]->SetMarkerColor(kGray);
			gGdiss3S[ipt]->SetMarkerStyle(24);
			if ( ipt==0 ){
				gGdiss3S[ipt]->SetLineColor(2);
				gGdiss3S[ipt]->SetMarkerColor(2);
			}else if ( ipt==1 ){
				gGdiss3S[ipt]->SetLineColor(kGreen+2);
				gGdiss3S[ipt]->SetMarkerColor(kGreen+2);
			}
			gGdiss3S[ipt]->Draw("P");

			fGdiss3S[ipt]->SetRange(160,180);
			fGdiss3S[ipt]->SetLineColor(1);
			fGdiss3S[ipt]->SetLineWidth(1);
			fGdiss3S[ipt]->SetLineStyle(7);
			fGdiss3S[ipt]->Draw("same");
		}

		TLegend *leg = new TLegend(0.15,0.65,0.4,0.95);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(43);
		leg->SetTextSize(18);
		leg->AddEntry("","#Upsilon(3S)","h");
		leg->AddEntry(gGdiss2S[0],"p_{T}=0 GeV/c","p");
		leg->AddEntry(gGdiss2S[1],"p_{T}=2 GeV/c","p");
		leg->AddEntry("","","");
		leg->AddEntry("","","");
		leg->AddEntry("","","");
		leg->Draw();
	}

	if ( bDRAW_Config ){
		c0_ = new TCanvas("c0_","c0_",1.2*2*500,500);
		c0_->Divide(2,1);

		c0_->cd(1);
		SetPadStyle();
		gPad->SetTopMargin(0.05);

		htmp = (TH1D*)gPad->DrawFrame(170,0,600,1.5e-8);
		SetHistoStyle("T [MeV]","C_{reg} [MeV]");

		for (int ipt=0; ipt<21; ipt++){
			gCregen[ipt]->SetLineWidth(3);
			gCregen[ipt]->SetLineColorAlpha(kGray,0.5);
			gCregen[ipt]->SetMarkerColor(kGray);
			gCregen[ipt]->SetMarkerStyle(24);
			if ( ipt==0 ){
				gCregen[ipt]->SetLineColor(2);
				gCregen[ipt]->SetMarkerColor(2);
			}else if ( ipt==5 ){
				gCregen[ipt]->SetLineColor(kGreen+2);
				gCregen[ipt]->SetMarkerColor(kGreen+2);
			}else if ( ipt==10 ){
				gCregen[ipt]->SetLineColor(kBlue);
				gCregen[ipt]->SetMarkerColor(kBlue);
			}else if ( ipt==15 ){
				gCregen[ipt]->SetLineColor(kMagenta);
				gCregen[ipt]->SetMarkerColor(kMagenta);
			}else if ( ipt==20 ){
				gCregen[ipt]->SetLineColor(1);
				gCregen[ipt]->SetMarkerColor(1);
			}
			gCregen[ipt]->Draw("P");

			fCregen[ipt]->SetRange(160,600);
			fCregen[ipt]->SetLineColor(1);
			fCregen[ipt]->SetLineWidth(1);
			fCregen[ipt]->SetLineStyle(7);
			fCregen[ipt]->Draw("same");
		}

		{
			TLegend *leg = new TLegend(0.15,0.65,0.4,0.95);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextFont(43);
			leg->SetTextSize(18);
			leg->AddEntry("","#Upsilon(1S)","h");
			leg->AddEntry(gCregen[0],"p_{T}=0 GeV/c","p");
			leg->AddEntry(gCregen[5],"p_{T}=5 GeV/c","p");
			leg->AddEntry(gCregen[10],"p_{T}=10 GeV/c","p");
			leg->AddEntry(gCregen[15],"p_{T}=15 GeV/c","p");
			leg->AddEntry(gCregen[20],"p_{T}=20 GeV/c","p");
			leg->Draw();
		}

		c0_->cd(2);
		SetPadStyle();
		gPad->SetTopMargin(0.05);

		htmp = (TH1D*)gPad->DrawFrame(0,0,20,2.0e-10);
		SetHistoStyle("p_{T} (GeV/c)","dN_{reg}/dt (fm/c)^{-1}");

		for (int ipt=0; ipt<44; ipt++){
			hRregen[ipt]->SetLineWidth(1);
			hRregen[ipt]->SetLineColorAlpha(kGray,0.3);
			if ( ipt==3 ){
				hRregen[ipt]->SetLineColor(2);
				hRregen[ipt]->SetMarkerColor(2);
			}else if ( ipt==13 ){
				hRregen[ipt]->SetLineColor(kGreen+2);
				hRregen[ipt]->SetMarkerColor(kGreen+2);
			}else if ( ipt==23 ){
				hRregen[ipt]->SetLineColor(kBlue);
				hRregen[ipt]->SetMarkerColor(kBlue);
			}else if ( ipt==33 ){
				hRregen[ipt]->SetLineColor(kMagenta);
				hRregen[ipt]->SetMarkerColor(kMagenta);
			}else if ( ipt==43 ){
				hRregen[ipt]->SetLineColor(1);
				hRregen[ipt]->SetMarkerColor(1);
			}
			hRregen[ipt]->Draw("same");
		}
		hRregen[ 3]->Draw("same");
		hRregen[13]->Draw("same");
		hRregen[23]->Draw("same");
		hRregen[33]->Draw("same");
		hRregen[43]->Draw("same");

		{
			TLegend *leg = new TLegend(0.15,0.65,0.4,0.95);
			leg->SetFillStyle(0);
			leg->SetBorderSize(0);
			leg->SetTextFont(43);
			leg->SetTextSize(18);
			leg->AddEntry("","#Upsilon(1S)","h");
			leg->AddEntry(hRregen[ 3],"T=200 MeV","L");
			leg->AddEntry(hRregen[13],"T=300 MeV","L");
			leg->AddEntry(hRregen[23],"T=400 MeV","L");
			leg->AddEntry(hRregen[33],"T=500 MeV","L");
			leg->AddEntry(hRregen[43],"T=600 MeV","L");
			leg->Draw();
		}
	}

	//return;

	//Upsilon beta vs. p 
	TF1 *fP1S = new TF1("fP1S","[0]*x/sqrt(1-x*x)",0,1);
	fP1S->SetParameter(0, const_mY1S);

	TF1 *fP2S = new TF1("fP2S","[0]*x/sqrt(1-x*x)",0,1);
	fP2S->SetParameter(0, const_mY2S);

	TF1 *fP3S = new TF1("fP3S","[0]*x/sqrt(1-x*x)",0,1);
	fP3S->SetParameter(0, const_mY3S);

	TCanvas *c1_;
	if ( bDRAW_Config ){
		c1_ = new TCanvas("c1_","c1_",1.2*500,500);
		SetPadStyle();
		htmp = (TH1D*)gPad->DrawFrame(0,0,1,25);
		SetHistoStyle("#beta","p [GeV]");

		fP1S->SetLineColor(1);
		fP1S->SetLineWidth(3);
		fP1S->Draw("same");

		fP2S->SetLineColor(2);
		fP2S->SetLineWidth(3);
		fP2S->Draw("same");

		fP3S->SetLineColor(4);
		fP3S->SetLineWidth(3);
		fP3S->Draw("same");
	}

	//return;
	
	//Glauber
	TFile *infileGlauber;
	if ( system=="pPb" ){
		infileGlauber = new TFile("/alice/data/junleekim/SHINCHON/MCGlauber/MCGlauber-pPb-8160GeV-b0-10fm.root","read");
	}else if ( system=="pO" ){
		infileGlauber = new TFile("/alice/data/junleekim/SHINCHON/MCGlauber/MCGlauber-pO-8160GeV-b0-10fm.root","read");
	}else if ( system=="OO" ){
		infileGlauber = new TFile("/alice/data/junleekim/SHINCHON/MCGlauber/MCGlauber-OO-8160GeV-b0-10fm.root","read");
	}
	TTree *TGlauber = (TTree*)infileGlauber->Get("lemon");
	int Gnpart, Gncoll;
	TGlauber->SetBranchAddress("npart",&Gnpart);
	TGlauber->SetBranchAddress("ncoll",&Gncoll);

	/*
	TFile *infileMult;
	if ( system=="pPb" ){
		infileMult = new TFile("totalout_pPb_fixed_v1.root","read");
	}else if ( system=="pO" ){
		infileMult = new TFile("totalout_pO_fixed_v1.root","read");
	}else if ( system=="OO" ){
		infileMult = new TFile("totalout_OO_fixed_v1.root","read");
	}
	*/

	TFile *infileCent = new TFile("outfile_Centrality_pPb_pO_OO_8TeV.root","read");
	TH1D *hMult = (TH1D*)infileCent->Get(Form("hMult_%s",system.c_str()));
	TH1D *hCent = (TH1D*)infileCent->Get(Form("hCent_%s",system.c_str()));


	//QGP T profile
	TH2D *hTHydro[300];
	TGraphErrors *gTHydro[nrun];
	TF1 *fT2[nrun];
	double timeHydro[nrun][300] = {0.0}; 
	float errTHydro[nrun][300] = {0.0};
	float freezeT[nrun];
	float Npart[nrun];

	TProfile2D *hprofRAA_xy[nrun];
	TProfile2D *hprofRAA_xy_rl[nrun];

	TProfile *hprofRAA1S_Npart = new TProfile("hprofRAA1S_Npart","",50,0,npartmax);
	TProfile *hprofRAA2S_Npart = new TProfile("hprofRAA2S_Npart","",50,0,npartmax);
	TProfile *hprofRAA3S_Npart = new TProfile("hprofRAA3S_Npart","",50,0,npartmax);

	TProfile *hprofRAA1S_pT = new TProfile("hprofRAA1S_pT","",20,0,20);
	TProfile *hprofRAA2S_pT = new TProfile("hprofRAA2S_pT","",20,0,20);
	TProfile *hprofRAA3S_pT = new TProfile("hprofRAA3S_pT","",20,0,20);

	TProfile *hprofRAA1S_pT_rl = new TProfile("hprofRAA1S_pT_rl","",20,0,20);

	TProfile *hprofRAA1S_Npart_rl = new TProfile("hprofRAA1S_Npart_rl","",50,0,npartmax);
	TProfile *hprofRAA2S_Npart_rl = new TProfile("hprofRAA2S_Npart_rl","",50,0,npartmax);
	TProfile *hprofRAA3S_Npart_rl = new TProfile("hprofRAA3S_Npart_rl","",50,0,npartmax);

	TH1D *htau_form_Y1S = new TH1D("htau_form_Y1S","",100,0,2.5);
	TH1D *htau_form_Y2S = new TH1D("htau_form_Y2S","",100,0,2.5);
	TH1D *htau_form_Y3S = new TH1D("htau_form_Y3S","",100,0,2.5);

	TH1D *hpT1S_i = new TH1D("hpT_i","",100,0,20);
	TH1D *hpT1S_f = new TH1D("hpT_f","",100,0,20);
	TH1D *hpT1S_r = new TH1D("hpT_r","",100,0,20);
	hpT1S_r->Sumw2();

	TProfile *hprofRAA1S_Npart_allpT = new TProfile("hprofRAA1S_Npart_allpT","",50,0,npartmax);
	TProfile *hprofRAA2S_Npart_allpT = new TProfile("hprofRAA2S_Npart_allpT","",50,0,npartmax);
	TProfile *hprofRAA3S_Npart_allpT = new TProfile("hprofRAA3S_Npart_allpT","",50,0,npartmax);

	TProfile *hprofRAA1S_Mult_allpT = new TProfile("hprofRAA1S_Mult_allpT","",50,0,nmultmax);
	TProfile *hprofRAA2S_Mult_allpT = new TProfile("hprofRAA2S_Mult_allpT","",50,0,nmultmax);
	TProfile *hprofRAA3S_Mult_allpT = new TProfile("hprofRAA3S_Mult_allpT","",50,0,nmultmax);

	TProfile *hprofRAA1S_Multst_allpT = new TProfile("hprofRAA1S_Multst_allpT","",50,0,20);
	TProfile *hprofRAA2S_Multst_allpT = new TProfile("hprofRAA2S_Multst_allpT","",50,0,20);
	TProfile *hprofRAA3S_Multst_allpT = new TProfile("hprofRAA3S_Multst_allpT","",50,0,20);

	TProfile *hprofsT_Npart = new TProfile("hprofsT_Npart","",50,0,npartmax);
	TProfile *hprofsT_Mult = new TProfile("hprofsT_Mult","",50,0,nmultmax);

//****
	const int nPtBin = 10;
	double PtBin[nPtBin] = {
		0, 2, 4, 6, 8,
		12, 25, 50, 100, 1e4 };

	const int nNpartBin = 10;
	double NpartBin[nNpartBin] = {
		20, 50, 80, 100, 150,
		200, 250, 300, 350, 1e3 };
	const int nPhiBin = 33;
	double PhiBin[nPhiBin];
	for(int phi=0;phi<nPhiBin;phi++) PhiBin[phi] = (TMath::Pi()*2.0 / nPhiBin) * phi;

	//      TProfile3D* hFinalState = new TProfile3D("hFinalState","hFinalState",
	TH3D* hFinalState = new TH3D("hFinalState","hFinalState",
			nPtBin-1, PtBin, nNpartBin-1, NpartBin, nPhiBin-1, PhiBin);
	//**** temporary bin definition.

	//Init
	TF1 *fInitialUpsilon = new TF1("fInitialUpsilon",fTsallis1S_v2,0,30,3);
	fInitialUpsilon -> SetParameters(  1.06450e+00 ,  7.97649e-01 , 100);

	TF1 *fInitY = new TF1("fInitY","[2]*x/TMath::Power(((x/[1])*(x/[1])+1.),[0])",0,30);
	fInitY->SetParameters(2.44, 6.05, 1);

	TF1 *fInitb = new TF1("fInitb","[2]*x/TMath::Power(((x/[1])*(x/[1])+1.),[0])",0,30);
	fInitb->SetParameters(2.85, 6.07, 1);

	if ( bDRAW_Config ){
		TCanvas *ctmp = new TCanvas("ctmp","ctmp");
		fInitb->SetLineStyle(2);
		fInitb->Draw("");
		fInitY->Draw("same");
	}

	cout << "before: " << fInitb->Integral(0,30) << " " << fInitY->Integral(0,30) << endl;

	fInitb->SetParameter(2, 1./fInitb->Integral(0,30));
	fInitY->SetParameter(2, 1./fInitY->Integral(0,30));

	cout << "after: " << fInitb->Integral(0,30) << " " << fInitY->Integral(0,30) << endl;

	double nY_tot = 0.0;

	//return;

	for (int irun=run_i; irun<run_f; irun++){

		if ( irun==918 ) continue;

		//TH1D *hMult = (TH1D*)infileMult->Get(Form("hMultDist_%d_0",irun));

		float Mult = hMult->GetBinContent(irun+1);
		float Cent = hCent->GetBinContent(irun+1);

		cout << "Scan event #" << irun << ", Mult: " << Mult << ", Cent: " << Cent << endl;

		if ( Cent>90 ) continue;

		//Glauber info
		TGlauber->GetEntry(irun);
		Npart[irun] = Gnpart;

		TH2D *hGlauber = (TH2D*)infileGlauber->Get(Form("inited_event%d",irun));

		float GmeanXw = 0, GmeanYw = 0;
		float GmeanX2w = 0, GmeanY2w = 0;
		float GmeanXYw = 0;
		float Gsumw = 0;

		float Gnbinsx = hGlauber->GetNbinsX();
		float Gnbinsy = hGlauber->GetNbinsY();

		for (int ix=0; ix<Gnbinsx; ix++){
			for (int iy=0; iy<Gnbinsy; iy++){
				float xx = hGlauber->GetXaxis()->GetBinCenter(ix+1);
				float yy = hGlauber->GetYaxis()->GetBinCenter(iy+1);

				float ww = hGlauber->GetBinContent(ix+1, iy+1);

				GmeanXw += xx*ww;
				GmeanYw += yy*ww;
				GmeanX2w += xx*xx*ww;
				GmeanY2w += yy*yy*ww;
				GmeanXYw += xx*yy*ww;

				Gsumw += ww; 
			}   
		}   

		GmeanXw /= Gsumw;
		GmeanYw /= Gsumw;
		GmeanX2w /= Gsumw;
		GmeanY2w /= Gsumw;
		GmeanXYw /= Gsumw;

		float GvarXw = GmeanX2w - GmeanXw*GmeanXw;
		float GvarYw = GmeanY2w - GmeanYw*GmeanYw;
		float GvarXYw = GmeanXYw - GmeanXw*GmeanYw;
		float GsTwsq = GvarXw*GvarYw - GvarXYw*GvarXYw;

		if ( GsTwsq<=0 ){
			cout << "Invalid sTsq: " << GsTwsq << ", " << system << ", run: " << irun << endl; 
			continue;
		}

		float GsTw = TMath::Pi()*sqrt(GsTwsq); 

		hprofsT_Npart->Fill(Npart[irun], GsTw);
		hprofsT_Mult->Fill(Mult, GsTw);

		if ( Mult/GsTw>2.4 && Mult/GsTw<2.8 ){
			cout << "######## " << Mult/GsTw << endl;
		}

		//continue;

		if ( b2D ){
			hprofRAA_xy[irun] = new TProfile2D(Form("hprofRAA_xy_run%05d",irun),"",100,-15,15,100,-15,15);
			hprofRAA_xy_rl[irun] = new TProfile2D(Form("hprofRAA_xy_rl_run%05d",irun),"",100,-15,15,100,-15,15);
		}

		//TFile *infileHydro = new TFile(Form("/alice/data/junleekim/SHINCHON/superSONIC_profile_pPb_8160GeV/superSONIC_profile_pPb_8160GeV_event%05d.root",irun),"read");
		//TFile *infileHydro = new TFile(Form("/alice/data/junleekim/SHINCHON/superSONIC_profile_OO_8160GeV/superSONIC_profile_OO_8160GeV_event%05d.root",irun),"read");
		//TFile *infileHydro = new TFile(Form("/alice/data/shlim/SONIC/SONIC_profile_PbPb5TeV_0_18fm_t0_0_3_v2/SONIC_profile_PbPb5TeV_0_18fm_event%05d.root",irun),"read");
		TFile *infileHydro;
		if ( system=="pPb" ){
			infileHydro = new TFile(Form("/alice/data/junleekim/SHINCHON/superSONIC_profile_pPb_8160GeV_fine/superSONIC_profile_pPb_8160GeV_fine_event%05d.root",irun),"read");
		}else if ( system=="pO" ){
			infileHydro = new TFile(Form("/alice/data/junleekim/SHINCHON/superSONIC_profile_pO_8160GeV_fine/superSONIC_profile_pO_8160GeV_fine_event%05d.root",irun),"read");
		}else if ( system=="OO" ){
			infileHydro = new TFile(Form("/alice/data/junleekim/SHINCHON/superSONIC_profile_OO_8160GeV_fine/superSONIC_profile_OO_8160GeV_fine_event%05d.root",irun),"read");
		}
		TH1D *htimeHydro = (TH1D*)infileHydro->Get("Time");
		int ntimeHydro = (int)htimeHydro->GetEntries();

		//Load histograms
		for (int it=0; it<ntimeHydro; it++){
			hTHydro[it] = (TH2D*)infileHydro->Get(Form("T_%d",it));
			timeHydro[irun][it] = htimeHydro->GetBinContent(it+1);
		}

		const int nY = nSAMP * Gncoll;

		nY_tot += nY;

		//Upsilon
		for (int iY=0; iY<nY; iY++){

			//Momentum
			//double pT = 20.0*gRandom->Rndm();
			//double pT = 3.0;
			//double pT = fInitialUpsilon->GetRandom();
			double pT = fInitY->GetRandom();
			double phi = (gRandom->Rndm()-0.5)*TMath::TwoPi(); 
			double px = pT*cos(phi);
			double py = pT*sin(phi);

			double rap = 2.5*gRandom->Rndm();
			double mT = sqrt(pT*pT + const_mY1S*const_mY1S);
			double pz = mT*sinh(rap);

			int pTbin = int(pT);
			if ( pTbin>=20 ) continue;

			//Beta
			double bx = fP1S->GetX(fabs(px)); 
			double by = fP1S->GetX(fabs(py)); 

			if ( px<0 ) bx *= -1;
			if ( py<0 ) by *= -1;

			TLorentzVector lvec;
			if ( bPTONLY ){
				lvec.SetPxPyPzE(px, py, 0, mT);
			}else{
				lvec.SetPxPyPzE(px, py, pz, sqrt(const_mY1S*const_mY1S + pT*pT + pz*pz));
			}

			//formation time
			double tau_formY1S = const_tau0Y1S*lvec.Gamma();
			htau_form_Y1S->Fill(tau_formY1S);

			//Position
			double vx, vy;
			hGlauber->GetRandom2(vx, vy);

			double vx0 = vx;
			double vy0 = vy;

			double modF = 1.0;
			int nFO = 0;

			double det_flg = 1.0;

			//Pre-QGP
			if ( bPreQGP )
			{
				//approximation of 20% higher T in pre-hydro
				float TPre = 1.2*hTHydro[0]->GetBinContent(hTHydro[0]->FindBin(vx, vy))*1000.;
				if ( TPre>const_Tf && tau_formY1S<0.3 ){
					float GdissPre0 = fGdiss[pTbin]->Eval(TPre);
					float GdissPre1 = fGdiss[pTbin+1]->Eval(TPre);
					float GdissPre = GdissPre0 + (GdissPre1 - GdissPre0)*(pT - pTbin);

					float dt = 0.3 - tau_formY1S;
					modF = exp(-(dt)*GdissPre/const_hbarc);
					if( exp(-(dt)*GdissPre/const_hbarc) < gRandom->Rndm() ) det_flg = 0.0;

					//cout << tau_form << " " << modF << endl;
				}
			}

			vx += bx*0.3;
			vy += by*0.3;

			//Time evolution
			for (int it=0; it<ntimeHydro-1; it++){

				if ( nFO>=10 ) break;

				float dt = timeHydro[irun][it+1] - timeHydro[irun][it];
				float dx = bx*dt;
				float dy = by*dt;

				float THydro0 = hTHydro[it]->GetBinContent(hTHydro[it]->FindBin(vx, vy))*1000.;
				float THydro1 = hTHydro[it+1]->GetBinContent(hTHydro[it+1]->FindBin(vx+dx, vy+dy))*1000.;

				float Gdiss0_1S = fGdiss[pTbin]->Eval((THydro0+THydro1)/2);
				float Gdiss1_1S = fGdiss[pTbin+1]->Eval((THydro0+THydro1)/2);
				float Gdiss_1S = Gdiss0_1S + (Gdiss1_1S - Gdiss0_1S)*(pT - pTbin); 

				if ( (THydro0+THydro1)/2<const_Tf ){
					vx += dx;
					vy += dy;
					nFO++;
					continue;
				}

				if ( timeHydro[irun][it]<tau_formY1S ){
					//before formation time
					if ( bPreRES ){
						Gdiss_1S *= timeHydro[irun][it]/tau_formY1S;
					}else{
						vx += dx;
						vy += dy;
						continue;
					}

				}else{
					//after formation time
					if ( (THydro0+THydro1)/2>const_TmaxY1S ){
						//cout << "High enough for Y(1S), T:" << (THydro0+THydro1)/2 << endl;
						modF = 0.0;
						det_flg = 0.0;
						break;
					}

				}

				modF *= exp(-(dt)*Gdiss_1S/const_hbarc);
				if( exp(-(dt)*Gdiss_1S/const_hbarc) < gRandom->Rndm() ) det_flg = 0.0;

				vx += dx;
				vy += dy;

			}//it

			if ( pT>2.0 && pT<4.0 ){
				hprofRAA1S_Npart->Fill(Npart[irun], modF);
				hprofRAA1S_Npart_rl->Fill(Npart[irun], det_flg );
				if ( b2D ){
					hprofRAA_xy[irun]->Fill(vx0, vy0, modF);
					hprofRAA_xy_rl[irun]->Fill(vx0, vy0, det_flg );
				}
			}
			hFinalState->Fill( pT, Npart[irun], TVector2(vx,vy).Phi(), det_flg );
			
			hprofRAA1S_Npart_allpT->Fill(Npart[irun], modF);
			hprofRAA1S_Mult_allpT->Fill(Mult, modF);
			hprofRAA1S_Multst_allpT->Fill(Mult/GsTw, modF);
			hprofRAA1S_pT->Fill(pT, modF);
			hprofRAA1S_pT_rl->Fill(pT, det_flg );

			hpT1S_i->Fill(pT);
			if ( det_flg ){
				hpT1S_f->Fill(pT);
			}
		}//iY

		//Upsilon 2S
		for (int iY=0; iY<nY; iY++){

			//double pT = fInitialUpsilon->GetRandom();
			double pT = fInitY->GetRandom();
			double phi = (gRandom->Rndm()-0.5)*TMath::TwoPi(); 
			double px = pT*cos(phi);
			double py = pT*sin(phi);

			double rap = 2.5*gRandom->Rndm();
			double mT = sqrt(pT*pT + const_mY2S*const_mY2S);
			double pz = mT*sinh(rap);

			if ( pT>=20.0 ) continue;
			int pTbin = int(pT/2);

			//Beta
			double bx = fP2S->GetX(fabs(px)); 
			double by = fP2S->GetX(fabs(py)); 

			if ( px<0 ) bx *= -1;
			if ( py<0 ) by *= -1;

			TLorentzVector lvec;
			if ( bPTONLY ){
				lvec.SetPxPyPzE(px, py, 0, mT);
			}else{
				lvec.SetPxPyPzE(px, py, pz, sqrt(const_mY2S*const_mY2S + pT*pT + pz*pz));
			}

			//formation time
			double tau_formY1S = const_tau0Y1S*lvec.Gamma();
			double tau_formY2S = const_tau0Y2S*lvec.Gamma();
			htau_form_Y2S->Fill(tau_formY2S);

			//Position
			double vx, vy;
			hGlauber->GetRandom2(vx, vy);

			double vx0 = vx;
			double vy0 = vy;

			double modF = 1.0;
			int nFO = 0;

			double det_flg = 1.0;

			vx += bx*0.3;
			vy += by*0.3;

			//Time evolution
			for (int it=0; it<ntimeHydro-1; it++){

				if ( nFO>=10 ) break;

				float dt = timeHydro[irun][it+1] - timeHydro[irun][it];
				float dx = bx*dt;
				float dy = by*dt;

				float THydro0 = hTHydro[it]->GetBinContent(hTHydro[it]->FindBin(vx, vy))*1000.;
				float THydro1 = hTHydro[it+1]->GetBinContent(hTHydro[it+1]->FindBin(vx+dx, vy+dy))*1000.;

				float Gdiss0_1S = fGdiss[pTbin]->Eval((THydro0+THydro1)/2);
				float Gdiss1_1S = fGdiss[pTbin+1]->Eval((THydro0+THydro1)/2);
				float Gdiss_1S = Gdiss0_1S + (Gdiss1_1S - Gdiss0_1S)*(pT - pTbin*2)/2.0; 

				float Gdiss0_2S = fGdiss2S[pTbin]->Eval((THydro0+THydro1)/2);
				float Gdiss1_2S = fGdiss2S[pTbin+1]->Eval((THydro0+THydro1)/2);
				float Gdiss_2S = Gdiss0_2S + (Gdiss1_2S - Gdiss0_2S)*(pT - pTbin*2)/2.0; 

				if ( (THydro0+THydro1)/2<const_Tf ){
					vx += dx;
					vy += dy;
					nFO++;
					continue;
				}

				if ( timeHydro[irun][it]<tau_formY2S ){
					//before formation time
					if ( bPreRES ){
						Gdiss_1S *= timeHydro[irun][it]/tau_formY1S;
						Gdiss_2S *= timeHydro[irun][it]/tau_formY2S;
					}else{
						vx += dx;
						vy += dy;
						continue;
					}

				}else{
					//after formation time
					if ( (THydro0+THydro1)/2>const_TmaxY2S ){
						modF = 0.0;
						det_flg = 0.0;
						break;
					}

				}

				/*
				if ( timeHydro[irun][it]<tau_formY1S ){
					modF *= exp(-(dt)*Gdiss_1S/const_hbarc);
					if( exp(-(dt)*Gdiss_1S/const_hbarc) < gRandom->Rndm() ) det_flg = 0.0;
				}else{
					modF *= exp(-(dt)*Gdiss_2S/const_hbarc);
					if( exp(-(dt)*Gdiss_2S/const_hbarc) < gRandom->Rndm() ) det_flg = 0.0;
				}
				*/
				modF *= exp(-(dt)*Gdiss_2S/const_hbarc);
				if( exp(-(dt)*Gdiss_2S/const_hbarc) < gRandom->Rndm() ) det_flg = 0.0;

				vx += dx;
				vy += dy;

			}//it

			if ( pT>2.0 && pT<4.0 ){
				hprofRAA2S_Npart->Fill(Npart[irun], modF);
				hprofRAA2S_Npart_rl->Fill(Npart[irun], det_flg );
			}

			hprofRAA2S_Npart_allpT->Fill(Npart[irun], modF);
			hprofRAA2S_Mult_allpT->Fill(Mult, modF);
			hprofRAA2S_Multst_allpT->Fill(Mult/GsTw, modF);
			hprofRAA2S_pT->Fill(pT, modF);
		}//iY

		//Upsilon 3S
		for (int iY=0; iY<nY; iY++){

			//double pT = fInitialUpsilon->GetRandom();
			double pT = fInitY->GetRandom();
			double phi = (gRandom->Rndm()-0.5)*TMath::TwoPi(); 
			double px = pT*cos(phi);
			double py = pT*sin(phi);

			double rap = 2.5*gRandom->Rndm();
			double mT = sqrt(pT*pT + const_mY3S*const_mY3S);
			double pz = mT*sinh(rap);

			if ( pT>=20.0 ) continue;
			int pTbin = int(pT/2);

			//Beta
			double bx = fP3S->GetX(fabs(px)); 
			double by = fP3S->GetX(fabs(py)); 

			if ( px<0 ) bx *= -1;
			if ( py<0 ) by *= -1;

			TLorentzVector lvec;
			if ( bPTONLY ){
				lvec.SetPxPyPzE(px, py, 0, mT);
			}else{
				lvec.SetPxPyPzE(px, py, pz, sqrt(const_mY3S*const_mY3S + pT*pT + pz*pz));
			}

			//formation time
			double tau_formY1S = const_tau0Y1S*lvec.Gamma();
			double tau_formY2S = const_tau0Y2S*lvec.Gamma();
			double tau_formY3S = const_tau0Y3S*lvec.Gamma();
			htau_form_Y3S->Fill(tau_formY3S);

			//Position
			double vx, vy;
			hGlauber->GetRandom2(vx, vy);

			double vx0 = vx;
			double vy0 = vy;

			double modF = 1.0;
			int nFO = 0;

			double det_flg = 1.0;

			vx += bx*0.3;
			vy += by*0.3;

			//Time evolution
			for (int it=0; it<ntimeHydro-1; it++){

				if ( nFO>=10 ) break;

				float dt = timeHydro[irun][it+1] - timeHydro[irun][it];
				float dx = bx*dt;
				float dy = by*dt;

				float THydro0 = hTHydro[it]->GetBinContent(hTHydro[it]->FindBin(vx, vy))*1000.;
				float THydro1 = hTHydro[it+1]->GetBinContent(hTHydro[it+1]->FindBin(vx+dx, vy+dy))*1000.;

				float Gdiss0_1S = fGdiss[pTbin]->Eval((THydro0+THydro1)/2);
				float Gdiss1_1S = fGdiss[pTbin+1]->Eval((THydro0+THydro1)/2);
				float Gdiss_1S = Gdiss0_1S + (Gdiss1_1S - Gdiss0_1S)*(pT - pTbin*2)/2.0; 

				float Gdiss0_2S = fGdiss2S[pTbin]->Eval((THydro0+THydro1)/2);
				float Gdiss1_2S = fGdiss2S[pTbin+1]->Eval((THydro0+THydro1)/2);
				float Gdiss_2S = Gdiss0_2S + (Gdiss1_2S - Gdiss0_2S)*(pT - pTbin*2)/2.0; 

				float GdissRef = fGdiss3S[0]->Eval((THydro0+THydro1)/2);
				float Gdiss0_3S = GdissRef*(fGdiss2S[pTbin]->Eval((THydro0+THydro1)/2))/(fGdiss2S[0]->Eval((THydro0+THydro1)/2));
				float Gdiss1_3S = GdissRef*(fGdiss2S[pTbin+1]->Eval((THydro0+THydro1)/2))/(fGdiss2S[0]->Eval((THydro0+THydro1)/2));
				float Gdiss_3S = Gdiss0_3S + (Gdiss1_3S - Gdiss0_3S)*(pT - pTbin*2)/2.0; 

				if ( (THydro0+THydro1)/2<const_Tf ){
					vx += dx;
					vy += dy;
					nFO++;
					continue;
				}

				if ( timeHydro[irun][it]<tau_formY3S ){
					//before formation time
					if ( bPreRES ){
						Gdiss_1S *= timeHydro[irun][it]/tau_formY1S;
						Gdiss_2S *= timeHydro[irun][it]/tau_formY2S;
						Gdiss_3S *= timeHydro[irun][it]/tau_formY3S;
					}else{
						vx += dx;
						vy += dy;
						continue;
					}

				}else{
					//after formation time
					if ( (THydro0+THydro1)/2>const_TmaxY3S ){
						modF = 0.0;
						det_flg = 0.0;
						break;
					}

				}
				
				/*
				if ( timeHydro[irun][it]<tau_formY1S ){
					modF *= exp(-(dt)*Gdiss_1S/const_hbarc);
					if( exp(-(dt)*Gdiss_1S/const_hbarc) < gRandom->Rndm() ) det_flg = 0.0;
				}else if ( timeHydro[irun][it]<tau_formY2S ){
					modF *= exp(-(dt)*Gdiss_2S/const_hbarc);
					if( exp(-(dt)*Gdiss_2S/const_hbarc) < gRandom->Rndm() ) det_flg = 0.0;
				}else{
					modF *= exp(-(dt)*Gdiss_3S/const_hbarc);
					if( exp(-(dt)*Gdiss_3S/const_hbarc) < gRandom->Rndm() ) det_flg = 0.0;
				}
				*/
				modF *= exp(-(dt)*Gdiss_3S/const_hbarc);
				if( exp(-(dt)*Gdiss_3S/const_hbarc) < gRandom->Rndm() ) det_flg = 0.0;

				vx += dx;
				vy += dy;

			}//it

			if ( pT>2.0 && pT<4.0 ){
				hprofRAA3S_Npart->Fill(Npart[irun], modF);
				hprofRAA3S_Npart_rl->Fill(Npart[irun], det_flg );
			}

			hprofRAA3S_Npart_allpT->Fill(Npart[irun], modF);
			hprofRAA3S_Mult_allpT->Fill(Mult, modF);
			hprofRAA3S_Multst_allpT->Fill(Mult/GsTw, modF);
			hprofRAA3S_pT->Fill(pT, modF);
		}//iY

		/*
		//Upsilon 1S (regeneration)
		for (int ib=0; ib<nY; ib++){

			//double pT = fInitialUpsilon->GetRandom();
			double pT = fInitb->GetRandom();
			double phi = (gRandom->Rndm()-0.5)*TMath::TwoPi(); 
			double px = pT*cos(phi);
			double py = pT*sin(phi);

			//Beta
			double bx = fP1S->GetX(fabs(px)); 
			double by = fP1S->GetX(fabs(py)); 

			//Position
			double vx, vy;
			hGlauber->GetRandom2(vx, vy);

			double vx0 = vx;
			double vy0 = vy;

			vx += bx*0.3;
			vy += by*0.3;

			//Time evolution
			for (int it=0; it<ntimeHydro-1; it++){

				float dt = timeHydro[irun][it+1] - timeHydro[irun][it];
				float dx = bx*dt;
				float dy = by*dt;

				float THydro0 = hTHydro[it]->GetBinContent(hTHydro[it]->FindBin(vx, vy))*1000.;
				float THydro1 = hTHydro[it+1]->GetBinContent(hTHydro[it+1]->FindBin(vx+dx, vy+dy))*1000.;

				if ( (THydro0+THydro1)/2>const_TmaxY1S || (THydro0+THydro1)/2<Tf ){
					continue;
				}

				int Tbin = int((THydro0+THydro1)/2 - 170.0)/10;

				int pTregen = hRregen[Tbin]->GetRandom(); 
				float Rregen = hRregen[Tbin]->GetBinContent(hRregen[Tbin]->FindBin(pTregen));

				hpT1S_r->Fill(pTregen, Rregen*dt);

				//cout << (THydro0+THydro1)/2 << " " << Tbin << endl;

			}//it
		}//ib
		*/

		infileHydro->Close();
		delete infileHydro;

	}

	if ( bDRAW_Result ){

		TCanvas *c1 = new TCanvas("c1","c1",1.2*1*500,500);
		SetPadStyle();
		gPad->SetTicks();
		htmp = (TH1D*)htau_form_Y1S;
		htmp->GetXaxis()->SetTitle("#tau_{form} [fm/c]");
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetTitle("");
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetYaxis()->SetTitleSize(0.05);

		htau_form_Y1S->SetLineColor(1);
		htau_form_Y1S->SetLineWidth(2);
		htau_form_Y1S->Draw("");

		htau_form_Y2S->SetLineColor(2);
		htau_form_Y2S->SetLineWidth(2);
		htau_form_Y2S->Draw("same");

		htau_form_Y3S->SetLineColor(4);
		htau_form_Y3S->SetLineWidth(2);
		htau_form_Y3S->Draw("same");

		TLegend *leg = new TLegend(0.5,0.75,0.9,0.95);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(43);
		leg->SetTextSize(18);
		leg->AddEntry(htau_form_Y1S,Form("Y(1S), #LT#tau_{form}#GT=%4.2f",htau_form_Y1S->GetMean()),"L");
		leg->AddEntry(htau_form_Y2S,Form("Y(2S), #LT#tau_{form}#GT=%4.2f",htau_form_Y2S->GetMean()),"L");
		leg->AddEntry(htau_form_Y3S,Form("Y(3S), #LT#tau_{form}#GT=%4.2f",htau_form_Y3S->GetMean()),"L");
		leg->Draw();
	}

	if ( bDRAW_Result ){

		TCanvas *c2 = new TCanvas("c2","c2",1.2*2*500,500);
		c2->Divide(2,1);

		c2->cd(1);
		SetPadStyle();
		gPad->SetTicks();
		htmp = (TH1D*)gPad->DrawFrame(0,0,npartmax,1.2);
		//htmp->GetXaxis()->SetTitle("N_{part}");
		htmp->GetXaxis()->SetTitle("Multiplicity");
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetTitle("R_{AA}");
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetYaxis()->SetTitleSize(0.05);

		/*
		hprofRAA1S_Npart_allpT->Rebin();
		hprofRAA2S_Npart_allpT->Rebin();
		hprofRAA3S_Npart_allpT->Rebin();
		hprofRAA1S_Npart_allpT->SetLineColor(1);
		hprofRAA1S_Npart_allpT->SetLineWidth(2);
		hprofRAA1S_Npart_allpT->Draw("same");
		hprofRAA2S_Npart_allpT->SetLineColor(2);
		hprofRAA2S_Npart_allpT->SetLineWidth(2);
		hprofRAA2S_Npart_allpT->Draw("same");
		hprofRAA3S_Npart_allpT->SetLineColor(4);
		hprofRAA3S_Npart_allpT->SetLineWidth(2);
		hprofRAA3S_Npart_allpT->Draw("same");
		*/

		hprofRAA1S_Mult_allpT->Rebin();
		hprofRAA2S_Mult_allpT->Rebin();
		hprofRAA3S_Mult_allpT->Rebin();
		hprofRAA1S_Mult_allpT->SetLineColor(1);
		hprofRAA1S_Mult_allpT->SetLineWidth(2);
		hprofRAA1S_Mult_allpT->Draw("same");
		hprofRAA2S_Mult_allpT->SetLineColor(2);
		hprofRAA2S_Mult_allpT->SetLineWidth(2);
		hprofRAA2S_Mult_allpT->Draw("same");
		hprofRAA3S_Mult_allpT->SetLineColor(4);
		hprofRAA3S_Mult_allpT->SetLineWidth(2);
		hprofRAA3S_Mult_allpT->Draw("same");


		TGraphErrors *gRAA1S_FD = new TGraphErrors;
		for (int ii=0; ii<hprofRAA1S_Npart_allpT->GetNbinsX(); ii++){
			float xx = hprofRAA1S_Npart_allpT->GetBinCenter(ii+1);
			if ( xx>410 ) continue;

			float yy1S = hprofRAA1S_Npart_allpT->GetBinContent(ii+1);
			float yy2S = hprofRAA2S_Npart_allpT->GetBinContent(ii+1);
			float yy3S = hprofRAA3S_Npart_allpT->GetBinContent(ii+1);

			float yy = yy1S*0.67 + yy2S*0.26 + yy3S*0.07;

			gRAA1S_FD->SetPoint(ii, xx, yy);
		}
		gRAA1S_FD->SetLineColorAlpha(1, 0.5);
		gRAA1S_FD->SetLineWidth(3);
		gRAA1S_FD->SetLineStyle(7);
		//gRAA1S_FD->Draw("C");

		TGraphErrors *gRAA2S_FD = new TGraphErrors;
		for (int ii=0; ii<hprofRAA2S_Npart_allpT->GetNbinsX(); ii++){
			float xx = hprofRAA2S_Npart_allpT->GetBinCenter(ii+1);
			if ( xx>410 ) continue;

			float yy2S = hprofRAA2S_Npart_allpT->GetBinContent(ii+1);
			float yy3S = hprofRAA3S_Npart_allpT->GetBinContent(ii+1);

			float yy = yy2S*0.60 + yy3S*0.40;

			gRAA2S_FD->SetPoint(ii, xx, yy);
		}
		gRAA2S_FD->SetLineColorAlpha(2, 0.5);
		gRAA2S_FD->SetLineWidth(3);
		gRAA2S_FD->SetLineStyle(7);
		//gRAA2S_FD->Draw("C");

		TLegend *leg = new TLegend(0.4,0.68,0.9,0.93);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(43);
		leg->SetTextSize(18);
		leg->AddEntry("","Pb+Pb #sqrt{s_{NN}}=5.02 TeV","h");
		leg->AddEntry(hprofRAA1S_Npart_allpT,Form("Y(1S)"),"L");
		leg->AddEntry(hprofRAA2S_Npart_allpT,Form("Y(2S)"),"L");
		leg->AddEntry(hprofRAA3S_Npart_allpT,Form("Y(3S)"),"L");
		//leg->AddEntry(gRAA1S_FD,Form("Y(1S), Feed down corrected"),"L");
		//leg->AddEntry(gRAA2S_FD,Form("Y(2S), Feed down corrected"),"L");
		leg->AddEntry("","","");
		leg->AddEntry("","","");
		leg->Draw();


		c2->cd(2);
		SetPadStyle();
		gPad->SetTicks();
		htmp = (TH1D*)gPad->DrawFrame(0,0,20,1.2);
		htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetTitle("R_{AA}");
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetYaxis()->SetTitleSize(0.05);

		//hprofRAA1S_pT->Rebin(5);
		//hprofRAA2S_pT->Rebin(5);
		//hprofRAA3S_pT->Rebin(5);

		hprofRAA1S_pT->SetLineColor(1);
		hprofRAA1S_pT->SetLineWidth(2);
		hprofRAA1S_pT->Draw("same");

		hprofRAA2S_pT->SetLineColor(2);
		hprofRAA2S_pT->SetLineWidth(2);
		hprofRAA2S_pT->Draw("same");

		hprofRAA3S_pT->SetLineColor(4);
		hprofRAA3S_pT->SetLineWidth(2);
		hprofRAA3S_pT->Draw("same");

		TGraphErrors *gRAA1S_pT_FD = new TGraphErrors;
		for (int ii=0; ii<hprofRAA1S_pT->GetNbinsX(); ii++){
			float xx = hprofRAA1S_pT->GetBinCenter(ii+1);

			float yy1S = hprofRAA1S_pT->GetBinContent(ii+1);
			float yy2S = hprofRAA2S_pT->GetBinContent(ii+1);
			float yy3S = hprofRAA3S_pT->GetBinContent(ii+1);

			float yy = yy1S*0.67 + yy2S*0.26 + yy3S*0.07;

			gRAA1S_pT_FD->SetPoint(ii, xx, yy);
		}
		gRAA1S_pT_FD->SetLineColorAlpha(1, 0.5);
		gRAA1S_pT_FD->SetLineWidth(3);
		gRAA1S_pT_FD->SetLineStyle(7);
		//gRAA1S_pT_FD->Draw("C");

		TGraphErrors *gRAA2S_pT_FD = new TGraphErrors;
		for (int ii=0; ii<hprofRAA2S_pT->GetNbinsX(); ii++){
			float xx = hprofRAA2S_pT->GetBinCenter(ii+1);

			float yy2S = hprofRAA2S_pT->GetBinContent(ii+1);
			float yy3S = hprofRAA3S_pT->GetBinContent(ii+1);

			float yy = yy2S*0.60 + yy3S*0.40;

			gRAA2S_pT_FD->SetPoint(ii, xx, yy);
		}
		gRAA2S_pT_FD->SetLineColorAlpha(2, 0.5);
		gRAA2S_pT_FD->SetLineWidth(3);
		gRAA2S_pT_FD->SetLineStyle(7);
		//gRAA2S_pT_FD->Draw("C");

		leg->Draw();

		/*

		TH1D *hRAA1S_pT_diss = (TH1D*)hpT1S_f->Clone("hRAA1S_pT_diss");
		hRAA1S_pT_diss->Sumw2();
		hRAA1S_pT_diss->Divide(hpT1S_i);
		hRAA1S_pT_diss->SetMarkerStyle(24);
		hRAA1S_pT_diss->Draw("p same");
		*/

		/*
		TH1D *hRAA_pT_all = (TH1D*)hpT_final->Clone("hRAA_pT_all");
		hRAA_pT_all->Sumw2();
		hpT_regen->Scale(1000.*hpT_init->Integral());
		hRAA_pT_all->Add(hpT_regen);
		hRAA_pT_all->Divide(hpT_init);
		hRAA_pT_all->SetMarkerStyle(25);
		hRAA_pT_all->SetMarkerColor(2);
		hRAA_pT_all->SetLineColor(2);
		hRAA_pT_all->Draw("p same");
		*/

		/*
		hpT1S_i->Rebin(5);
		hpT1S_f->Rebin(5);
		hpT1S_r->Rebin(5);

		TCanvas *c3 = new TCanvas("c3","c3",1.2*500,500);
		SetPadStyle();
		gPad->SetTicks();
		gPad->SetLogy();
		htmp = (TH1D*)gPad->DrawFrame(0,0.5,20,2*hpT1S_i->GetMaximum());
		htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetTitle("N");
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetYaxis()->SetTitleSize(0.05);

		hpT1S_i->SetMarkerStyle(20);
		hpT1S_i->SetMarkerColor(1);
		hpT1S_i->SetLineColor(1);
		hpT1S_i->Draw("p same");

		hpT1S_f->SetMarkerStyle(24);
		hpT1S_f->SetMarkerColor(1);
		hpT1S_f->SetLineColor(1);
		hpT1S_f->Draw("p same");

		hpT1S_r->Scale(9.93*nY_tot/(2e-3));
		hpT1S_r->SetMarkerStyle(21);
		hpT1S_r->SetMarkerSize(0.8);
		hpT1S_r->SetMarkerColor(1);
		hpT1S_r->SetLineColor(1);
		hpT1S_r->Draw("p same");
		
		cout << hpT1S_i->Integral() << endl;
		cout << hpT1S_f->Integral() << endl;
		cout << hpT1S_r->Integral() << endl;
		*/

	}

	if ( bSAVE ){

		TFile *outfile = new TFile(Form("outfile_RaaV8_%s_%04d_%04d.root",system.c_str(),run_i,run_f),"recreate");

		hprofRAA1S_Npart->Write();
		hprofRAA1S_pT->Write();
		hprofRAA2S_pT->Write();
		hprofRAA3S_pT->Write();
		hprofRAA1S_Npart_allpT->Write();
		hprofRAA2S_Npart_allpT->Write();
		hprofRAA3S_Npart_allpT->Write();

		hprofRAA1S_Mult_allpT->Write();
		hprofRAA2S_Mult_allpT->Write();
		hprofRAA3S_Mult_allpT->Write();

		hprofRAA1S_Multst_allpT->Write();
		hprofRAA2S_Multst_allpT->Write();
		hprofRAA3S_Multst_allpT->Write();

		hprofsT_Npart->Write();
		hprofsT_Mult->Write();

		if ( b2D ){
			for (int irun=run_i; irun<run_f; irun++){
				hprofRAA_xy[irun]->Write();
			}
		}
		hFinalState->Write();
	}

	return;



}

Double_t fTsallis1S_v2(Double_t *x, Double_t *fpar)
{
  Float_t xx = x[0];
  Double_t Y1Smass = 9.46;
  Double_t q = fpar[0];
  Double_t T = fpar[1];
  Double_t c = fpar[2];
  Double_t mT = TMath::Sqrt(Y1Smass*Y1Smass+xx*xx);
  Double_t pow = TMath::Power((1+(q-1)*mT/T),(-q/(q-1)));

  Double_t f = c*mT*xx*pow;
  return f;
}
