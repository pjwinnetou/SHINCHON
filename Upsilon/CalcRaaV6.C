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

#include <iostream>
#include <fstream>

using namespace std;

Double_t fTsallis1S_v2(Double_t *x, Double_t *fpar);

void CalcRaaV6(){

	gStyle->SetOptStat(0);
	gStyle->SetPalette(55);

	gRandom = new TRandom3(0);

	const bool bDRAW = true; 
	const bool bSAVE = false; 
	const bool b2D = false;

	const int run_i = 0;
	const int run_f = 1;
	const int nrun = 1000;
	const float const_hbarc = 197.5; //MeV fm
	const float const_mY = 9.46; //GeV
	const float const_mY2S = 10.02; //GeV
	const int nSAMP = 1; //times Ncoll

	const float npartmax = 450;
	const float Tf = 192.0;

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

	for (int ipt=0; ipt<21; ipt++){
		gCregen[ipt] = new TGraphErrors;
		for (unsigned int iT=0; iT<T[ipt].size(); iT++){
			gCregen[ipt]->SetPoint(iT, T[ipt][iT], Cregen[ipt][iT]);
		}
	}

	TCanvas *c0;
	TCanvas *c0_2S;
	TCanvas *c0_;
	if ( bDRAW ){
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

		TLegend *leg = new TLegend(0.15,0.70,0.4,0.95);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(43);
		leg->SetTextSize(18);
		leg->AddEntry(gGdiss[0],"p_{T}=0 GeV/c","p");
		leg->AddEntry(gGdiss[5],"p_{T}=5 GeV/c","p");
		leg->AddEntry(gGdiss[10],"p_{T}=10 GeV/c","p");
		leg->AddEntry(gGdiss[15],"p_{T}=15 GeV/c","p");
		leg->AddEntry(gGdiss[20],"p_{T}=20 GeV/c","p");
		leg->Draw();

	}

	if ( bDRAW ){
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

		TLegend *leg = new TLegend(0.15,0.7,0.4,0.95);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(43);
		leg->SetTextSize(18);
		leg->AddEntry(gGdiss2S[0],"p_{T}=0 GeV/c","p");
		leg->AddEntry(gGdiss2S[2],"p_{T}=4 GeV/c","p");
		leg->AddEntry(gGdiss2S[4],"p_{T}=8 GeV/c","p");
		leg->AddEntry(gGdiss2S[6],"p_{T}=12 GeV/c","p");
		leg->AddEntry(gGdiss2S[8],"p_{T}=16 GeV/c","p");
		leg->Draw();
	}

	if ( bDRAW ){
		c0_ = new TCanvas("c0_","c0_",1.2*500,500);
		SetPadStyle();
		gPad->SetTopMargin(0.05);

		htmp = (TH1D*)gPad->DrawFrame(170,0,500,1e-8);
		SetHistoStyle("T [MeV]","C_{reg} [MeV]");

		for (int ipt=0; ipt<21; ipt++){
			gCregen[ipt]->SetLineWidth(3);
			gCregen[ipt]->SetLineColorAlpha(kGray,0.5);
			if ( ipt==0 ){
				gCregen[ipt]->SetLineColor(2);
			}else if ( ipt==5 ){
				gCregen[ipt]->SetLineColor(kGreen+2);
			}else if ( ipt==10 ){
				gCregen[ipt]->SetLineColor(kBlue);
			}else if ( ipt==15 ){
				gCregen[ipt]->SetLineColor(kMagenta);
			}else if ( ipt==20 ){
				gCregen[ipt]->SetLineColor(1);
			}
			gCregen[ipt]->Draw("C");
		}
	}

	//return;

	//Upsilon beta vs. p 
	TF1 *fP = new TF1("fP","[0]*x/sqrt(1-x*x)",0,1);
	fP->SetParameter(0, const_mY);

	TF1 *fP2S = new TF1("fP2S","[0]*x/sqrt(1-x*x)",0,1);
	fP2S->SetParameter(0, const_mY2S);

	TCanvas *c1_;
	if ( bDRAW ){
		c1_ = new TCanvas("c1_","c1_",1.2*500,500);
		SetPadStyle();
		htmp = (TH1D*)gPad->DrawFrame(0,0,1,25);
		SetHistoStyle("#beta","p [GeV]");

		fP->SetLineColor(1);
		fP->SetLineWidth(3);
		fP->Draw("same");

		fP2S->SetLineColor(2);
		fP2S->SetLineWidth(3);
		fP2S->Draw("same");
	}

	//return;
	
	//Glauber
	TFile *infileGlauber = new TFile("/alice/home/shlim/work/SHINCHON/SHINCHON/Software/SONIC-KIAF/input/MCGlauber-PbPb-5020GeV-b0-18fm.root","read");
	//TFile *infileGlauber = new TFile("/alice/data/junleekim/SHINCHON/MCGlauber/MCGlauber-OO-8160GeV-b0-10fm.root","read");
	//TFile *infileGlauber = new TFile("/alice/data/junleekim/SHINCHON/MCGlauber/MCGlauber-pPb-8160GeV-b0-10fm.root","read");
	//TFile *infileGlauber = new TFile("/alice/home/jseo/MCGlauber-PbPb-5020GeV-b0-18fm.root","read");
	TTree *TGlauber = (TTree*)infileGlauber->Get("lemon");
	int Gnpart, Gncoll;
	TGlauber->SetBranchAddress("npart",&Gnpart);
	TGlauber->SetBranchAddress("ncoll",&Gncoll);

	//QGP T profile
	TH2D *hTHydro[300];
	TGraphErrors *gTHydro[nrun];
	TF1 *fT2[nrun];
	double timeHydro[nrun][300] = {0.0}; 
	float errTHydro[nrun][300] = {0.0};
	float freezeT[nrun];
	float Npart[nrun];

	TProfile *hprofRAA_Npart = new TProfile("hprofRAA_Npart","",50,0,npartmax);
	TProfile *hprofRAA_pT = new TProfile("hprofRAA_pT","",100,0,20);
	TProfile2D *hprofRAA_xy[nrun];

	TProfile *hprofRAA2S_Npart = new TProfile("hprofRAA2S_Npart","",50,0,npartmax);
	TProfile *hprofRAA2S_pT = new TProfile("hprofRAA2S_pT","",100,0,20);


	TProfile *hprofRAA_Npart_rl = new TProfile("hprofRAA_Npart_rl","",50,0,npartmax);
	TProfile *hprofRAA_pT_rl = new TProfile("hprofRAA_pT_rl","",100,0,20);
	TProfile2D *hprofRAA_xy_rl[nrun];

	TProfile *hprofRAA2S_Npart_rl = new TProfile("hprofRAA2S_Npart_rl","",50,0,npartmax);


	TH1D *hpT_init = new TH1D("hpT_init","",100,0,20);
	TH1D *hpT_final = new TH1D("hpT_final","",100,0,20);
	TH1D *hpT_regen = new TH1D("hpT_regen","",100,0,20);
	hpT_regen->Sumw2();

	TProfile *hprofRAA_Npart_allpT = new TProfile("hprofRAA_Npart_allpT","",50,0,npartmax);
	TProfile *hprofRAA2S_Npart_allpT = new TProfile("hprofRAA2S_Npart_allpT","",50,0,npartmax);

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
	//Init

	for (int irun=run_i; irun<run_f; irun++){

		cout << "Scan event #" << irun << endl;

		if ( b2D ){
			hprofRAA_xy[irun] = new TProfile2D(Form("hprofRAA_xy_run%05d",irun),"",100,-15,15,100,-15,15);
			hprofRAA_xy_rl[irun] = new TProfile2D(Form("hprofRAA_xy_rl_run%05d",irun),"",100,-15,15,100,-15,15);
		}

		//TFile *infileHydro = new TFile(Form("/alice/data/junleekim/SHINCHON/superSONIC_profile_pPb_8160GeV/superSONIC_profile_pPb_8160GeV_event%05d.root",irun),"read");
		//TFile *infileHydro = new TFile(Form("/alice/data/junleekim/SHINCHON/superSONIC_profile_OO_8160GeV/superSONIC_profile_OO_8160GeV_event%05d.root",irun),"read");
		TFile *infileHydro = new TFile(Form("/alice/data/shlim/SONIC/SONIC_profile_PbPb5TeV_0_18fm_t0_0_3_v2/SONIC_profile_PbPb5TeV_0_18fm_event%05d.root",irun),"read");
		//TFile *infileHydro = new TFile(Form("SONIC_profile_PbPb5TeV_0_18fm_event%05d.root",irun),"read");
		TH1D *htimeHydro = (TH1D*)infileHydro->Get("Time");
		int ntimeHydro = (int)htimeHydro->GetEntries();

		//Glauber info
		TGlauber->GetEntry(irun);
		Npart[irun] = Gnpart;

		TH2D *hGlauber = (TH2D*)infileGlauber->Get(Form("inited_event%d",irun));

		//Load histograms
		for (int it=0; it<ntimeHydro; it++){
			hTHydro[it] = (TH2D*)infileHydro->Get(Form("T_%d",it));
			timeHydro[irun][it] = htimeHydro->GetBinContent(it+1);
		}

		const int nY = nSAMP * Gncoll;

		//Upsilon
		for (int iY=0; iY<nY; iY++){

			//Momentum
			//double pT = 20.0*gRandom->Rndm();
			//double pT = 3.0;
			double pT = fInitialUpsilon->GetRandom();
			double phi = (gRandom->Rndm()-0.5)*TMath::TwoPi(); 
			double px = pT*cos(phi);
			double py = pT*sin(phi);

			int pTbin = int(pT);
			if ( pTbin>=20 ) continue;

			//Beta
			double bx = fP->GetX(fabs(px)); 
			double by = fP->GetX(fabs(py)); 

			if ( px<0 ) bx *= -1;
			if ( py<0 ) by *= -1;

			TLorentzVector lvec;
			lvec.SetPxPyPzE(px, py, 0, sqrt(const_mY*const_mY + pT*pT));

			//formation time
			double tau_form = 0.2*lvec.Gamma();

			//Position
			double vx, vy;
			hGlauber->GetRandom2(vx, vy);

			double vx0 = vx;
			double vy0 = vy;

			double modF = 1.0;
			int nFO = 0;

			double det_flg = 1.0;

			double Yregen = 0.0;

			//Pre-Hydro
			/*
			{
				float TPre = 1.15*hTHydro[0]->GetBinContent(hTHydro[0]->FindBin(vx, vy))*1000.;
				if ( TPre>170.0 && tau_form<0.3 ){
					float GdissPre0 = fGdiss[pTbin]->Eval(TPre);
					float GdissPre1 = fGdiss[pTbin+1]->Eval(TPre);
					float GdissPre = GdissPre0 + (GdissPre1 - GdissPre0)*(pT - pTbin);
					modF = exp(-(0.2)*GdissPre/const_hbarc);
					if( exp(-(0.2)*GdissPre/const_hbarc) < gRandom->Rndm() ) det_flg = 0.0;
				}
			}
			*/

			vx += bx*0.3;
			vy += by*0.3;

			//Time evolution
			for (int it=0; it<ntimeHydro-1; it++){

				if ( nFO>=10 ) break;

				float dt = timeHydro[irun][it+1] - timeHydro[irun][it];
				float dx = bx*dt;
				float dy = by*dt;

				if ( tau_form>timeHydro[irun][it] ){
					vx += dx;
					vy += dy;
					continue;
				}

				float THydro0 = hTHydro[it]->GetBinContent(hTHydro[it]->FindBin(vx, vy))*1000.;
				float THydro1 = hTHydro[it+1]->GetBinContent(hTHydro[it+1]->FindBin(vx+dx, vy+dy))*1000.;

				if ( (THydro0+THydro1)/2>600.0 ){
					cout << "High enough for Y(1S), T:" << (THydro0+THydro1)/2 << endl;
					modF = 0.0;
					det_flg = 0.0;
					break;
				}

				if ( (THydro0+THydro1)/2<Tf ){
					vx += dx;
					vy += dy;
					nFO++;
					continue;
				}

				float Gdiss0 = fGdiss[pTbin]->Eval((THydro0+THydro1)/2);
				float Gdiss1 = fGdiss[pTbin+1]->Eval((THydro0+THydro1)/2);
				float Gdiss = Gdiss0 + (Gdiss1 - Gdiss0)*(pT - pTbin); 
				modF *= exp(-(dt)*Gdiss/const_hbarc);
				if( exp(-(dt)*Gdiss/const_hbarc) < gRandom->Rndm() ) det_flg = 0.0;

				float Cregen0 = gCregen[pTbin]->Eval((THydro0+THydro1)/2);
				float Cregen1 = gCregen[pTbin+1]->Eval((THydro0+THydro1)/2);
				float Cregen = Cregen0 + (Cregen1 - Cregen0)*(pT - pTbin); 

				Yregen += dt*Cregen/const_hbarc;

				vx += dx;
				vy += dy;

			}//it

			//cout << Yregen << endl;

			if ( pT>2.0 && pT<4.0 ){
				hprofRAA_Npart->Fill(Npart[irun], modF);
				hprofRAA_Npart_rl->Fill(Npart[irun], det_flg );
				if ( b2D ){
					hprofRAA_xy[irun]->Fill(vx0, vy0, modF);
					hprofRAA_xy_rl[irun]->Fill(vx0, vy0, det_flg );
				}
			}
			hFinalState->Fill( pT, Npart[irun], TVector2(vx,vy).Phi(), det_flg );
			
			hprofRAA_Npart_allpT->Fill(Npart[irun], modF);

			hprofRAA_pT->Fill(pT, modF);
			hprofRAA_pT_rl->Fill(pT, det_flg );

			hpT_init->Fill(pT);
			if ( det_flg ){
				hpT_final->Fill(pT);
			}
			hpT_regen->Fill(pT, Yregen);
		}//iY

		//Upsilon 2S
		for (int iY=0; iY<nY; iY++){

			double pT = fInitialUpsilon->GetRandom();
			double phi = (gRandom->Rndm()-0.5)*TMath::TwoPi(); 
			double px = pT*cos(phi);
			double py = pT*sin(phi);

			if ( pT>=20.0 ) continue;
			int pTbin = int(pT/2);

			//Beta
			double bx = fP2S->GetX(fabs(px)); 
			double by = fP2S->GetX(fabs(py)); 

			if ( px<0 ) bx *= -1;
			if ( py<0 ) by *= -1;

			TLorentzVector lvec;
			lvec.SetPxPyPzE(px, py, 0, sqrt(const_mY2S*const_mY2S + pT*pT));

			//formation time
			double tau_form = 0.4*lvec.Gamma();

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

				if ( tau_form>timeHydro[irun][it] ){
					vx += dx;
					vy += dy;
					continue;
				}

				float THydro0 = hTHydro[it]->GetBinContent(hTHydro[it]->FindBin(vx, vy))*1000.;
				float THydro1 = hTHydro[it+1]->GetBinContent(hTHydro[it+1]->FindBin(vx+dx, vy+dy))*1000.;

				if ( (THydro0+THydro1)/2>230.0 ){
					modF = 0.0;
					det_flg = 0.0;
					break;
				}

				if ( (THydro0+THydro1)/2<Tf ){
					vx += dx;
					vy += dy;
					nFO++;
					continue;
				}

				float Gdiss0 = fGdiss2S[pTbin]->Eval((THydro0+THydro1)/2);
				float Gdiss1 = fGdiss2S[pTbin+1]->Eval((THydro0+THydro1)/2);
				float Gdiss = Gdiss0 + (Gdiss1 - Gdiss0)*(pT - pTbin*2)/2.0; 
				modF *= exp(-(dt)*Gdiss/const_hbarc);
				if( exp(-(dt)*Gdiss/const_hbarc) < gRandom->Rndm() ) det_flg = 0.0;

				vx += dx;
				vy += dy;

			}//it

			if ( pT>2.0 && pT<4.0 ){
				hprofRAA2S_Npart->Fill(Npart[irun], modF);
				hprofRAA2S_Npart_rl->Fill(Npart[irun], det_flg );
			}

			hprofRAA2S_Npart_allpT->Fill(Npart[irun], modF);

			hprofRAA2S_pT->Fill(pT, modF);


		}//iY

		infileHydro->Close();
		delete infileHydro;

	}

	if ( bDRAW ){

		TCanvas *c1 = new TCanvas("c1","c1",1.2*2*500,500);
		c1->Divide(2,1);

		c1->cd(1);
		SetPadStyle();
		gPad->SetRightMargin(0.13);
		htmp = (TH1D*)gPad->DrawFrame(-10,-10,10,10);
		htmp->GetXaxis()->SetTitle("x [fm]");
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetTitle("y [fm]");
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetYaxis()->SetTitleSize(0.05);

		TH2D *hGlauber = (TH2D*)infileGlauber->Get("inited_event0");
		hGlauber->GetZaxis()->SetTitle("");
		hGlauber->Draw("colz same");

		c1->cd(2);
		SetPadStyle();
		gPad->SetRightMargin(0.13);
		htmp = (TH1D*)gPad->DrawFrame(-10,-10,10,10);
		htmp->GetXaxis()->SetTitle("x [fm]");
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetTitle("y [fm]");
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetYaxis()->SetTitleSize(0.05);

		//hprofRAA_xy[0]->SetMaximum(1);
		//hprofRAA_xy[0]->SetMinimum(0);
		//hprofRAA_xy[0]->Draw("colz same");
	}

	if ( bDRAW ){

		TCanvas *c2 = new TCanvas("c2","c2",1.2*2*500,500);
		c2->Divide(2,1);

		c2->cd(1);
		SetPadStyle();
		gPad->SetRightMargin(0.13);
		htmp = (TH1D*)gPad->DrawFrame(0,0,npartmax,1.2);
		htmp->GetXaxis()->SetTitle("N_{part}");
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetTitle("R_{AA}");
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetYaxis()->SetTitleSize(0.05);

		hprofRAA_Npart_allpT->Rebin();
		hprofRAA2S_Npart_allpT->Rebin();

		hprofRAA_Npart_allpT->SetLineColor(1);
		hprofRAA_Npart_allpT->SetLineWidth(2);
		hprofRAA_Npart_allpT->Draw("same");

		hprofRAA2S_Npart_allpT->SetLineColor(2);
		hprofRAA2S_Npart_allpT->SetLineWidth(2);
		hprofRAA2S_Npart_allpT->Draw("same");


		c2->cd(2);
		SetPadStyle();
		gPad->SetRightMargin(0.13);
		htmp = (TH1D*)gPad->DrawFrame(0,0,20,1.2);
		htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetTitle("R_{AA}");
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetYaxis()->SetTitleSize(0.05);

		hprofRAA_pT->Rebin();
		hprofRAA2S_pT->Rebin();

		hprofRAA_pT->SetLineColor(1);
		hprofRAA_pT->SetLineWidth(2);
		hprofRAA_pT->Draw("same");

		hprofRAA2S_pT->SetLineColor(2);
		hprofRAA2S_pT->SetLineWidth(2);
		hprofRAA2S_pT->Draw("same");

		/*
		TH1D *hRAA_pT_diss = (TH1D*)hpT_final->Clone("hRAA_pT_diss");
		hRAA_pT_diss->Sumw2();
		hRAA_pT_diss->Divide(hpT_init);
		hRAA_pT_diss->SetMarkerStyle(24);
		hRAA_pT_diss->Draw("p same");

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

	}

	if ( bSAVE ){

		TFile *outfile = new TFile(Form("outfile_RaaV2_%04d_%04d.root",run_i,run_f),"recreate");

		hprofRAA_Npart->Write();
		hprofRAA_pT->Write();
		hprofRAA_Npart_allpT->Write();

		if ( b2D ){
			for (int irun=run_i; irun<run_f; irun++){
				hprofRAA_xy[irun]->Write();
			}
		}
		hFinalState->Write();
	}


	//hprofRAA_xy[0]->Draw("colz");

	return;

	//Calculate RAA



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
