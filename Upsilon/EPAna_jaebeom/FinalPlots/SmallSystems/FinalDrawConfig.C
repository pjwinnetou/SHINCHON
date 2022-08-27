#include "../../Style_jaebeom_woFrameLegend.h"
#include "../../SHINCHONLegend_raaCent.C"
#include "../../tdrstyle.C"

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

#include <iostream>
#include <fstream>

using namespace std;

Double_t fTsallis1S_v2(Double_t *x, Double_t *fpar);

void FinalDrawConfig(){

	gStyle->SetOptStat(0);
	gStyle->SetPalette(55);
	setTDRStyle();

	gRandom = new TRandom3(0);

	const bool bSAVE = false; 

	const bool bPreQGP = false; //should be off
	const bool bPreRES = false;

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
	const int nSAMP = 1; //times Ncoll

	const float npartmax = 450;

	ifstream fdata;

	char buf[500];
	vector<float> T[30];
	vector<float> Gdiss[30];
	vector<float> Cregen[30];

	float f_tmp[20];

	//thermal width Y(1S)
	fdata.open("../../../Gdiss0.dat");
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

	fdata.open("../../../Gdiss1.dat");
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
	fdata.open("../../../diss_2s.dat");
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
	fdata.open("../../../diss_3s.dat");
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
	fdata.open("../../../regen0.dat");
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

	fdata.open("../../../regen1.dat");
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

	fdata.open("../../../regen2.dat");
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

	TCanvas *c0_;

	{
		TCanvas *c0_1S = new TCanvas("c0","c0",1.2*500,500);
		gPad->SetLeftMargin(0.13);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(160,0,600,450);
		htmp->GetXaxis()->SetTitle("T (MeV)");
		htmp->GetYaxis()->SetTitle("#Gamma_{diss} (MeV)");
    htmp->GetYaxis()->SetTitleOffset(1.0);
    htmp->GetXaxis()->SetTitleOffset(1.05);
    htmp->GetXaxis()->SetTitleSize(0.055);
    htmp->GetXaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetNdivisions(510);

		for (int ipt=0; ipt<21; ipt++){
			gGdiss[ipt]->SetLineWidth(3);
			gGdiss[ipt]->SetLineColorAlpha(kGray,0.3);
			gGdiss[ipt]->SetMarkerColor(kGray);
			gGdiss[ipt]->SetMarkerStyle(24);

			fGdiss[ipt]->SetRange(160,600);
			fGdiss[ipt]->SetLineColorAlpha(1, 0.3);
			fGdiss[ipt]->SetLineWidth(2);

			if ( ipt==0 ){
				fGdiss[ipt]->SetLineColor(2);
				gGdiss[ipt]->SetLineColor(2);
				gGdiss[ipt]->SetMarkerColor(2);
			}else if ( ipt==5 ){
				fGdiss[ipt]->SetLineColor(kGreen+2);
				gGdiss[ipt]->SetLineColor(kGreen+2);
				gGdiss[ipt]->SetMarkerColor(kGreen+2);
			}else if ( ipt==10 ){
				fGdiss[ipt]->SetLineColor(kBlue);
				gGdiss[ipt]->SetLineColor(kBlue);
				gGdiss[ipt]->SetMarkerColor(kBlue);
			}else if ( ipt==15 ){
				fGdiss[ipt]->SetLineColor(kMagenta);
				gGdiss[ipt]->SetLineColor(kMagenta);
				gGdiss[ipt]->SetMarkerColor(kMagenta);
			}else if ( ipt==20 ){
				fGdiss[ipt]->SetLineColor(1);
				gGdiss[ipt]->SetLineColor(1);
				gGdiss[ipt]->SetMarkerColor(1);
			}
			//gGdiss[ipt]->Draw("P");

			fGdiss[ipt]->Draw("same");
		}

		TLegend *leg = new TLegend(0.15,0.93-0.055*6,0.5,0.93);
		SetLegendStyle(leg);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(43);
		leg->SetTextSize(20);
		leg->AddEntry("","#varUpsilon(1S)","h");
		leg->AddEntry(fGdiss[0],"p_{T}=0 GeV/c","L");
		leg->AddEntry(fGdiss[5],"p_{T}=5 GeV/c","L");
		leg->AddEntry(fGdiss[10],"p_{T}=10 GeV/c","L");
		leg->AddEntry(fGdiss[15],"p_{T}=15 GeV/c","L");
		leg->AddEntry(fGdiss[20],"p_{T}=20 GeV/c","L");
		leg->Draw();

		if ( bSAVE ){
			c0_1S->cd();
			c0_1S->SaveAs("ThermalWidth_Y1S.pdf");
		}

	}

	{
		TCanvas *c0_2S = new TCanvas("c0_2S","c0_2S",1.2*500,500);
		gPad->SetLeftMargin(0.13);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(160,0,250,450);
		htmp->GetXaxis()->SetTitle("T (MeV)");
		htmp->GetYaxis()->SetTitle("#Gamma_{diss} (MeV)");
    htmp->GetYaxis()->SetTitleOffset(1.0);
    htmp->GetXaxis()->SetTitleOffset(1.05);
    htmp->GetXaxis()->SetTitleSize(0.055);
    htmp->GetXaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetNdivisions(510);

		for (int ipt=0; ipt<11; ipt++){
			gGdiss2S[ipt]->SetLineWidth(2);
			gGdiss2S[ipt]->SetLineColorAlpha(kGray,0.3);
			gGdiss2S[ipt]->SetMarkerColor(kGray);
			gGdiss2S[ipt]->SetMarkerStyle(24);

			fGdiss2S[ipt]->SetRange(160,300);
			fGdiss2S[ipt]->SetLineColorAlpha(1, 0.3);
			fGdiss2S[ipt]->SetLineWidth(2);
			//fGdiss2S[ipt]->SetLineStyle(7);

			if ( ipt==0 ){
				fGdiss2S[ipt]->SetLineColor(2);
				gGdiss2S[ipt]->SetLineColor(2);
				gGdiss2S[ipt]->SetMarkerColor(2);
			}else if ( ipt==2 ){
				fGdiss2S[ipt]->SetLineColor(kGreen+2);
				gGdiss2S[ipt]->SetLineColor(kGreen+2);
				gGdiss2S[ipt]->SetMarkerColor(kGreen+2);
			}else if ( ipt==4 ){
				fGdiss2S[ipt]->SetLineColor(kBlue);
				gGdiss2S[ipt]->SetLineColor(kBlue);
				gGdiss2S[ipt]->SetMarkerColor(kBlue);
			}else if ( ipt==6 ){
				fGdiss2S[ipt]->SetLineColor(kMagenta);
				gGdiss2S[ipt]->SetLineColor(kMagenta);
				gGdiss2S[ipt]->SetMarkerColor(kMagenta);
			}else if ( ipt==8 ){
				fGdiss2S[ipt]->SetLineColor(1);
				gGdiss2S[ipt]->SetLineColor(1);
				gGdiss2S[ipt]->SetMarkerColor(1);
			}
			//gGdiss2S[ipt]->Draw("P");

			fGdiss2S[ipt]->Draw("same");
		}

		TLegend *leg = new TLegend(0.15,0.93-0.055*6,0.5,0.93);
		SetLegendStyle(leg);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(43);
		leg->SetTextSize(20);
		leg->AddEntry("","#varUpsilon(2S)","h");
		leg->AddEntry(fGdiss2S[0],"p_{T}=0 GeV/c","l");
		leg->AddEntry(fGdiss2S[2],"p_{T}=4 GeV/c","l");
		leg->AddEntry(fGdiss2S[4],"p_{T}=8 GeV/c","l");
		leg->AddEntry(fGdiss2S[6],"p_{T}=12 GeV/c","l");
		leg->AddEntry(fGdiss2S[8],"p_{T}=16 GeV/c","l");
		leg->Draw();

		if ( bSAVE ){
			c0_2S->cd();
			c0_2S->SaveAs("ThermalWidth_Y2S.pdf");
		}
	}

	{
		TCanvas *c0_3S = new TCanvas("c0_3S","c0_3S",1.2*500,500);
		gPad->SetLeftMargin(0.13);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(160,0,190,450);
		htmp->GetXaxis()->SetTitle("T (MeV)");
		htmp->GetYaxis()->SetTitle("#Gamma_{diss} (MeV)");
    htmp->GetYaxis()->SetTitleOffset(1.0);
    htmp->GetXaxis()->SetTitleOffset(1.05);
    htmp->GetXaxis()->SetTitleSize(0.055);
    htmp->GetXaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetNdivisions(510);

		for (int ipt=0; ipt<2; ipt++){
			gGdiss3S[ipt]->SetLineWidth(2);
			gGdiss3S[ipt]->SetLineColorAlpha(kGray,0.3);
			gGdiss3S[ipt]->SetMarkerColor(kGray);
			gGdiss3S[ipt]->SetMarkerStyle(24);

			fGdiss3S[ipt]->SetRange(160,190);
			fGdiss3S[ipt]->SetLineColorAlpha(1, 0.3);
			fGdiss3S[ipt]->SetLineWidth(2);
			//fGdiss3S[ipt]->SetLineStyle(7);

			if ( ipt==0 ){
				fGdiss3S[ipt]->SetLineColor(2);
				gGdiss3S[ipt]->SetLineColor(2);
				gGdiss3S[ipt]->SetMarkerColor(2);
			}else if ( ipt==1 ){
				fGdiss3S[ipt]->SetLineColor(kGreen+2);
				gGdiss3S[ipt]->SetLineColor(kGreen+2);
				gGdiss3S[ipt]->SetMarkerColor(kGreen+2);
			}
			//gGdiss3S[ipt]->Draw("P");

			fGdiss3S[ipt]->Draw("same");
		}

		TLegend *leg = new TLegend(0.15,0.93-0.055*6,0.5,0.93);
		SetLegendStyle(leg);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(43);
		leg->SetTextSize(20);
		leg->AddEntry("","#varUpsilon(3S)","h");
		leg->AddEntry(fGdiss2S[0],"p_{T}=0 GeV/c","l");
		leg->AddEntry(fGdiss2S[1],"p_{T}=2 GeV/c","l");
		leg->AddEntry("","","");
		leg->AddEntry("","","");
		leg->AddEntry("","","");
		leg->Draw();

		if ( bSAVE ){
			c0_3S->cd();
			c0_3S->SaveAs("ThermalWidth_Y3S.pdf");
		}
	}

	/*
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
	*/

	//return;

	//Upsilon beta vs. p 
	TF1 *fP1S = new TF1("fP1S","[0]*x/sqrt(1-x*x)",0,1);
	fP1S->SetParameter(0, const_mY1S);

	TF1 *fP2S = new TF1("fP2S","[0]*x/sqrt(1-x*x)",0,1);
	fP2S->SetParameter(0, const_mY2S);

	TF1 *fP3S = new TF1("fP3S","[0]*x/sqrt(1-x*x)",0,1);
	fP3S->SetParameter(0, const_mY3S);

	{
		TCanvas *c1_ = new TCanvas("c1_","c1_",1.2*500,500);
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,1,25);
		htmp->GetXaxis()->SetTitle("#beta");
		htmp->GetYaxis()->SetTitle("p (GeV)");
    htmp->GetYaxis()->SetTitleOffset(1.0);
    htmp->GetXaxis()->SetTitleOffset(1.05);
    htmp->GetXaxis()->SetTitleSize(0.055);
    htmp->GetXaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetNdivisions(510);

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
