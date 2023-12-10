#include "../Style.h"
//#include "Style_jaebeom.h"

#include <TClonesArray.h>
#include <TClass.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TProfile.h>
#include <TH2.h>
#include <TH3.h>
#include <TLegend.h>
#include <TTree.h>
#include <TVector2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TProfile2D.h>
#include <TF2.h>
#include <TLorentzVector.h>
#include "Math/Integrator.h"
#include "Math/Functor.h"

#include <iostream>
#include <fstream>
#include <string>
#include <dirent.h>
#include <memory>
#include <cstddef>

using namespace std;

double getMom(double v, double m){
  double mom = m*v / sqrt(1-v*v);
  return mom;
}

double getGamma(double v){
  double gamma = 1 / sqrt(1-v*v);
  return gamma;
}

Double_t fTsallis_v2(Double_t *x, Double_t *fpar){
  Float_t xx = x[0];
  Double_t q = fpar[0];
  Double_t T = fpar[1];
  Double_t c = fpar[2];
  Double_t Ymass = fpar[3];
  Double_t mT = TMath::Sqrt(Ymass*Ymass + xx*xx);
  Double_t pow = TMath::Power((1+(q-1)*mT/T),(-q/(q-1)));

  Double_t f = c*mT*xx*pow;
  return f;
}

void makeUpsTree( string Collision_system, int kInitPos = 1, int runId=0){


  gStyle->SetOptStat(0);
  gStyle->SetPalette(55);

  gRandom = new TRandom3(0);
  ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-6); //to converge

  const bool bDRAW = false; 
  const bool bSAVE = true; 
  const bool b2D = false;
  const bool bPreQGP = true;
  const bool bPreRES = false;

  const float PreHydroTempRatio = 1.20;

  const int run_i = runId;
  const int run_f = runId+10;
  const int nrun = 10;

/*
  const int run_i = 0;
  const int run_f = 1000;
  const int nrun = (run_f-run_i);
*/
  const float const_hbarc = 197.5; //MeV fm
  const int nstates = 3;

  const float const_mY[nstates] = {
	9.46 ,10.02, 10.36 }; //GeV
  const float const_tau0Y[nstates] = {
	0.5, 1.0, 1.5 }; //fm/c
//  0.3, 0.6, 0.9 };
  const float const_TmaxY[nstates] = {
	600.0, 240.0, 190.0}; //MeV

  const int nSAMP = 100; //times Ncoll

  const float npartmax = 450;
  const float Tf = 170.0;
  double sigs = 0.5;

  ifstream fdata;

  char buf[500];
  vector<float> T[21];
  vector<float> Gdiss[21];

  float f_tmp[20];

  fdata.open("/alice/data/junleekim/SHINCHON/DissConst/Gdiss0.dat");
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

  fdata.open("/alice/data/junleekim/SHINCHON/DissConst/Gdiss1.dat");
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

  TGraphErrors *gGdiss1S[21];
  TF1 *fGdiss1S[21];

  for (int ipt=0; ipt<21; ipt++){
    gGdiss1S[ipt] = new TGraphErrors;
    for (unsigned int iT=0; iT<T[ipt].size(); iT++){
      gGdiss1S[ipt]->SetPoint(iT, T[ipt][iT], Gdiss[ipt][iT]);
    }
    fGdiss1S[ipt] = new TF1(Form("fGdiss1S_%d",ipt),"pol5",0,520);
    gGdiss1S[ipt]->Fit(fGdiss1S[ipt],"R0Q");
  }

  for (int ii=0; ii<21; ii++){
    T[ii].clear();
    Gdiss[ii].clear();
  }


  fdata.open("/alice/data/junleekim/SHINCHON/DissConst/diss_2s.dat");
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
      }
    }
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


  for (int ii=0; ii<21; ii++){
    T[ii].clear();
    Gdiss[ii].clear();
  }

  fdata.open("/alice/data/junleekim/SHINCHON/DissConst/diss_3s.dat");
  fdata.getline(buf,500);

  while ( fdata
	>> f_tmp[0] >> f_tmp[1] >> f_tmp[2]
  ){
    for (int ipt=0; ipt<2; ipt++){
      if ( f_tmp[ipt+1]>0 ){
	T[ipt].push_back(f_tmp[0]);
	Gdiss[ipt].push_back(f_tmp[ipt+1]);
      }
    }
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


/*
  TCanvas *c0;
  if ( bDRAW ){
    c0 = new TCanvas("c0","c0",1.2*500,500);
    SetPadStyle();

    htmp = (TH1D*)gPad->DrawFrame(170,0,500,200);
    SetHistoStyle("T [MeV]","#Gamma_{diss} [MeV]");

    for (int ipt=0; ipt<21; ipt++){
      gGdiss1S[ipt]->SetLineWidth(3);
      gGdiss1S[ipt]->SetLineColorAlpha(kGray,0.5);
      if ( ipt==0 ){
        gGdiss1S[ipt]->SetLineColor(2);
      }else if ( ipt==5 ){
        gGdiss1S[ipt]->SetLineColor(kGreen+2);
      }else if ( ipt==10 ){
        gGdiss1S[ipt]->SetLineColor(kBlue);
      }else if ( ipt==15 ){
        gGdiss1S[ipt]->SetLineColor(kMagenta);
      }else if ( ipt==20 ){
        gGdiss1S[ipt]->SetLineColor(1);
      }
      gGdiss1S[ipt]->Draw("C");
    }

    TLegend *leg = new TLegend(0.2,0.6,0.6,0.9);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(43);
    leg->SetTextSize(20);
    leg->AddEntry(gGdiss1S[0],"p_{T}=0 GeV/c","L");
    leg->AddEntry(gGdiss1S[5],"p_{T}=5 GeV/c","L");
    leg->AddEntry(gGdiss1S[10],"p_{T}=10 GeV/c","L");
    leg->AddEntry(gGdiss1S[15],"p_{T}=15 GeV/c","L");
    leg->AddEntry(gGdiss1S[20],"p_{T}=20 GeV/c","L");
    leg->Draw();
  }
*/

  //Upsilon beta vs. p 
  TF1 *fP1S = new TF1("fP1S","[0]*x/sqrt(1-x*x)",0,1);
  fP1S->SetParameter(0, const_mY[0]);

  TF1 *fP2S = new TF1("fP2S","[0]*x/sqrt(1-x*x)",0,1);
  fP2S->SetParameter(0, const_mY[1]);

  TF1 *fP3S = new TF1("fP3S","[0]*x/sqrt(1-x*x)",0,1);
  fP3S->SetParameter(0, const_mY[2]);

/*
  TCanvas *c1_;
  if ( bDRAW ){
    c1_ = new TCanvas("c1_","c1_",1.2*500,500);
    SetPadStyle();
    htmp = (TH1D*)gPad->DrawFrame(0,0,1,25);
    SetHistoStyle("#beta","p [GeV]");

    fP1S->SetLineColor(1);
    fP1S->SetLineWidth(3);
    fP1S->Draw("same");
  }
*/

  //return;
  double phiN_ = 2;

  //Glauber
  const int Necc = 10;
  const int nNucl = 500;
  
  string GlbFileName ("MCGlauber-PbPb-5020GeV-b0-18fm-bin100-v3.root");
  int nNucl_proj = -1;
  int nNucl_targ = -1;

 
  if(GlbFileName.find("PbPb") != string::npos){nNucl_proj=208; nNucl_targ=208;}
  else if(GlbFileName.find("pPb") != string::npos){nNucl_proj=1; nNucl_targ=208;}
  else if(GlbFileName.find("AuAu") != string::npos){nNucl_proj=197; nNucl_targ=197;}
  else if(GlbFileName.find("dAu") != string::npos){nNucl_proj=2; nNucl_targ=197;}
  else if(GlbFileName.find("OO") != string::npos){nNucl_proj=16; nNucl_targ=16;}
  else if(GlbFileName.find("pO") != string::npos){nNucl_proj=1; nNucl_targ=16;}
  else if(GlbFileName.find("pp") != string::npos){nNucl_proj=1; nNucl_targ=1;}

  TFile* infileGlauber;
  if( ( Collision_system.find("pPb") != string::npos || 
	Collision_system.find("pO")  != string::npos ||
	Collision_system.find("OO")  != string::npos ||
	Collision_system.find("pp")  != string::npos ) && 
    Collision_system.find("trento") == string::npos){
    infileGlauber = new TFile(Form("/alice/data/junleekim/SHINCHON/MCGlauber/MCGlauber-%s-8160GeV-b0-10fm.root",Collision_system.c_str()),"read"); 
  } else if( Collision_system.find("trento") != string::npos ){
    string collname = Collision_system.substr( 0, Collision_system.find("_"));
    infileGlauber = new TFile(Form("/alice/data/junleekim/SHINCHON/MCGlauber/MCGlauber-%s-8160GeV-b0-10fm.root",collname.c_str()),"read");
  } else if( Collision_system.find("AuAu")  != string::npos ){ 
    infileGlauber = new TFile(Form("/alice/data/junleekim/SHINCHON/MCGlauber/MCGlauber-%s-200GeV-b0-18fm.root",Collision_system.c_str()),"read");
  } else{
    infileGlauber = new TFile(Form("../%s",GlbFileName.c_str()),"read");
  }

  TTree *TGlauber = (TTree*)infileGlauber->Get("lemon");

  int Gnpart, Gncoll;
  float b, eccgaus[Necc], eccpoint[Necc];
  float xproj[nNucl], yproj[nNucl], xtarg[nNucl], ytarg[nNucl];
  Bool_t wproj[nNucl], wtarg[nNucl]; 
  TGlauber->SetBranchAddress("npart",&Gnpart);
  TGlauber->SetBranchAddress("ncoll",&Gncoll);
  TGlauber->SetBranchAddress("b",&b);
  TGlauber->SetBranchAddress("eccgaus",eccgaus);
  TGlauber->SetBranchAddress("eccpoint",eccpoint);
  TGlauber->SetBranchAddress("xproj",xproj);
  TGlauber->SetBranchAddress("yproj",yproj);
  TGlauber->SetBranchAddress("xtarg",xtarg);
  TGlauber->SetBranchAddress("ytarg",ytarg);
  TGlauber->SetBranchAddress("wproj",wproj);
  TGlauber->SetBranchAddress("wtarg",wtarg);

  //QGP T profile
  TH2D *hTHydro[1000];
  TGraphErrors *gTHydro[nrun];
  TF1 *fT2[nrun];
  double timeHydro[nrun][1000] = {0.0}; 
  float errTHydro[nrun][1000] = {0.0};
  float freezeT[nrun];

  TProfile *hprofRAA_Npart = new TProfile("hprofRAA_Npart","",40,0,400);
  TProfile *hprofRAA_pT = new TProfile("hprofRAA_pT","",100,0,20);

  TProfile *hprofRAA_Npart_rl[nstates];
  TProfile *hprofRAA_mult_rl[nstates];
  TProfile *hprofRAA_multST_rl[nstates];
  TProfile *hprofRAA_pT_rl[nstates];

  TGraphErrors* gGdiss[nstates][100];
  for(int i=0;i<21;i++){
    gGdiss[0][i] = (TGraphErrors*)gGdiss1S[i]->Clone();
  }
  for(int i=0;i<11;i++){
    if( gGdiss2S[i]->GetN() ) gGdiss[1][i] = (TGraphErrors*)gGdiss2S[i]->Clone();
  }
  for(int i=0;i<2;i++){
    if( gGdiss3S[i]->GetN() ) gGdiss[2][i] = (TGraphErrors*)gGdiss3S[i]->Clone();
  }

  TF1 *fP[nstates];
  fP[0] = (TF1*)fP1S->Clone();
  fP[1] = (TF1*)fP2S->Clone();
  fP[2] = (TF1*)fP3S->Clone();

  for(int s=0;s<nstates;s++){
    hprofRAA_Npart_rl[s] 	= new TProfile(Form("hprofRAA_Npart_rl_%dS",s+1),"",40,0,400);
    hprofRAA_mult_rl[s] 	= new TProfile(Form("hprofRAA_mult_rl_%dS",s+1),"",1500,0,150);
    hprofRAA_multST_rl[s] 	= new TProfile(Form("hprofRAA_multST_rl_%dS",s+1),"",40,0,40);
    hprofRAA_pT_rl[s] 		= new TProfile(Form("hprofRAA_pT_rl_%dS",s+1),"",100,0,20);
  }

  //multiplicity
  TFile* fMultOut = new TFile(Form("/alice/data/junleekim/SHINCHON/mult/totalout_%s_fixed_v1.root",Collision_system.c_str()),"read");
  TH1D* hMult;  //to be changed

  //**** temporary bin definition.
  const int nPtBin = 10;
  double PtBin[nPtBin] = {
    0,   2,  4,   6,   8,
    12, 25, 50, 100, 1e4 };

  const int nNpartBin = 10;
  double PartBin[nNpartBin] = {
    20,   50,  80, 100, 150,
    200, 250, 300, 350, 1000 };

  TH2D* hTempEnergy[nrun];

  // cross check histogram 
  TH1D* hphi_initMom[nstates]; 
  TH1D* hphi_finalMom[nstates];  
  for(int s=0;s<nstates;s++){
    hphi_initMom[s]  = new TH1D(Form("hphi_initMom_%dS",s+1),"",100,-TMath::Pi(),TMath::Pi());
    hphi_finalMom[s] = new TH1D(Form("hphi_finalMom_%dS",s+1),"",100,-TMath::Pi(),TMath::Pi());
  }  

  TH2D* hPosFinal[nrun];
  TH2D* hPosInit[nrun]; 
  TProfile2D* hprofRAA_xy[nrun];
  TProfile2D* hprofRAA_xy_rl[nrun];

  //***** Tree & Output Definition
  string fInitPos;
  if(kInitPos==0) fInitPos = "InitPosZero";
  else if(kInitPos==1) fInitPos = "InitPosGlauber";
  else if(kInitPos==2) fInitPos = "InitPosMean";
  TFile *outfile = new TFile(Form("./outfile_UpsSkim_PhiAng_%s_%.f_%s_%04d_%04d.root",Collision_system.c_str(),phiN_,fInitPos.c_str(),run_i,run_f),"recreate");

  static const long MAXTREESIZE = 1000000000;
  static const long MaxUpsSize = 200000;

  int runid;
  int nUps;
  TClonesArray* Ups4momRaw = new TClonesArray("TLorentzVector", MaxUpsSize);
  TClonesArray* Ups4momEnProfCor = new TClonesArray("TLorentzVector", MaxUpsSize);
  TClonesArray* Ups4momGlauberCor = new TClonesArray("TLorentzVector", MaxUpsSize);
  double IsUpsSurv_modf[MaxUpsSize];
  double IsUpsSurv_prob[MaxUpsSize];
  double vx0_[MaxUpsSize];
  double vy0_[MaxUpsSize];
  double EPangEnProf, EPangGlauber, EPangGlauberh;
  int Npart_;
  int Ncoll_;
  float b_;
  float eccgaus_[Necc];
  float eccpoint_[Necc];


  TTree *tree = new TTree("tree","Upsilon Tree");
  tree->SetMaxTreeSize(MAXTREESIZE);
  tree->Branch("run",&runid,"run/I");
  tree->Branch("nUps",&nUps,"nUps/I");
  tree->Branch("Ups4momRaw", "TClonesArray", &Ups4momRaw, 32000, 0);
  tree->Branch("Ups4momEnProfCor", "TClonesArray", &Ups4momEnProfCor, 32000, 0);
  tree->Branch("Ups4momGlauberCor", "TClonesArray", &Ups4momGlauberCor, 32000, 0);
  tree->Branch("EPangEnProf",&EPangEnProf,"EPangEnProf/D");
  tree->Branch("EPangGlauber",&EPangGlauberh,"EPangGlauber/D");
  tree->Branch("IsUpsSurv_modf",IsUpsSurv_modf,"IsUpsSurv_modf[nUps]/D");
  tree->Branch("IsUpsSurv_prob",IsUpsSurv_prob,"IsUpsSurv_prob[nUps]/D");
  tree->Branch("vx0",vx0_,"vx0[nUps]/D");
  tree->Branch("vy0",vy0_,"vy0[nUps]/D");
  tree->Branch("Npart",&Npart_,"Npart/I");
  tree->Branch("Ncoll",&Ncoll_,"Ncoll/I");
  tree->Branch("b",&b_,"Ncoll/F");
  tree->Branch("eccgaus",eccgaus_,Form("eccgaus[%d]/F",Necc));
  tree->Branch("eccpoint",eccpoint_,Form("eccpoint[%d]/F",Necc));

  //Initial pT distribution from JJS
  TF1 *fInitialUpsilon = new TF1("fInitialUpsilon",fTsallis_v2,0,30,4);
  fInitialUpsilon -> SetParameters(  1.06450e+00 ,  7.97649e-01 , 100, const_mY[0] );
 
  //Glauber smearing function  
  TF2* smearing_function = new TF2("smear_tf2", "TMath::Exp(-(x*x+y*y)/(2.*[0]*[0]))/(2*TMath::Pi()*[0]*[0])", -100*sigs, 100*sigs, -100*sigs, 100*sigs);
  smearing_function->FixParameter(0,sigs);

  const Double_t max_x = 10.0;
  const int nSliceGlb = 300;
  TH2D* hist_glauber_fine[nrun];

  float meanX = 0, meanY = 0;
  float meanX2 = 0, meanY2 = 0;
  float meanXY = 0;
  float sum = 0;

  float varXw, varYw, varXYw;
  double sTsq, sT;

  for (int irun=run_i; irun<run_f; irun++){
    gRandom->SetSeed( irun );

    hMult = (TH1D*)fMultOut->Get(Form("hMultDist_%d_0",irun));

    cout << "Scan event #" << irun << endl;
    hTempEnergy[irun-run_i]= new TH2D(Form("histTempEnergy_%d",irun),";x;y",100,-10,10,100,-10,10);
    
    runid = irun;
    Ups4momRaw->Clear();
    Ups4momEnProfCor->Clear();
    Ups4momGlauberCor->Clear();
    nUps=0;

    if ( b2D ){
      hprofRAA_xy[irun-run_i] = new TProfile2D(Form("hprofRAA_xy_run%05d",irun),"",100,-15,15,100,-15,15);
      hprofRAA_xy_rl[irun-run_i] = new TProfile2D(Form("hprofRAA_xy_rl_run%05d",irun),"",100,-15,15,100,-15,15);
      hPosFinal[irun-run_i] = new TH2D(Form("hPosFinal_run%d",irun),"",100,-15,15,100,-15,15);
      hPosInit[irun-run_i] = new TH2D(Form("hPosInit_run%d",irun),"",100,-15,15,100,-15,15);
    }
    TFile* infileHydro;
    if( Collision_system.find("pPb") != string::npos || 
	Collision_system.find("pO")  != string::npos ||
	Collision_system.find("OO")  != string::npos ||
	Collision_system.find("pp")  != string::npos ){

      infileHydro = new TFile(Form("/alice/data/junleekim/SHINCHON/superSONIC_profile_%s_8160GeV_fine/superSONIC_profile_%s_8160GeV_fine_event%05d.root",Collision_system.c_str(),Collision_system.c_str(),irun),"read");
    } else if( Collision_system.find("AuAu") != string::npos ){
	  infileHydro = new TFile(Form("/alice/data/junleekim/SHINCHON/superSONIC_profile_%s_200GeV_fine/superSONIC_profile_%s_200GeV_fine_event%05d.root",Collision_system.c_str(),Collision_system.c_str(),irun),"read");
    } else{

	}
    TH1D *htimeHydro = (TH1D*)infileHydro->Get("Time");
    int ntimeHydro = (int)htimeHydro->GetEntries();

    //Glauber info
    TGlauber->GetEntry(irun);
    Npart_ = Gnpart;
    Ncoll_ = Gncoll;
    b_ = b;
    for(int iecc=0;iecc<Necc;iecc++){eccgaus_[iecc] = eccgaus[iecc]; eccpoint_[iecc] = eccpoint[iecc];}
    TH2D* hGlauber = (TH2D*)infileGlauber->Get(Form("inited_event%d",irun));

    ifstream fGlauber;
    if( Collision_system.find("trento") == string::npos ){
      if( Collision_system.find("pPb") != string::npos ){
	fGlauber.open(Form("/alice/data/junleekim/SHINCHON/initedFiles_pPb_8160GeV_b0_10fm/event%d_9999.dat",irun));
      } else if( Collision_system.find("pO") != string::npos ){
	fGlauber.open(Form("/alice/data/junleekim/SHINCHON/initedFiles_pO_8160GeV_b0_10fm/event%d_9998.dat",irun));
      } else if( Collision_system.find("OO") != string::npos ){
	fGlauber.open(Form("/alice/data/junleekim/SHINCHON/initedFiles_OO_8160GeV_b0_10fm/event%d_9997.dat",irun));
      } else if( Collision_system.find("pp") != string::npos ){
	fGlauber.open(Form("/alice/data/junleekim/SHINCHON/initedFiles_pp_8160GeV_b0_5fm/event%d_9996.dat",irun));
      } else if( Collision_system.find("AuAu") != string::npos ){
	fGlauber.open(Form("/alice/data/junleekim/SHINCHON/initedFiles_AuAu_200GeV_b0_18fm/event%d_9995.dat",irun));
	  }
    } else if( Collision_system.find("trento") != string::npos ){
      if( Collision_system.find("pPb") != string::npos ){
	if( Collision_system.find("p0") != string::npos ){
	  fGlauber.open(Form("/alice/data/shlim/SONIC/Trento_pPb_8160GeV_p0/event%d.dat",irun));
	} else if( Collision_system.find("p1") != string::npos ){
	  fGlauber.open(Form("/alice/data/shlim/SONIC/Trento_pPb_8160GeV_p1/event%d.dat",irun));
	}
      }
    }

    double ed;
    int index_x=0;
    int index_y=0;
    while( fGlauber >> ed ){
      hGlauber->SetBinContent( index_x+1, index_y+1, ed );
      index_y++;
      if( index_y == 100 ){
        index_y = 0;
        index_x++;
      }
    }
    fGlauber.close();
    
	// sT calculations
    meanX = 0.0;
    meanY = 0.0;
    meanX2 = 0.0;
    meanY2 = 0.0;
    meanXY = 0.0;
    sum = 0.0;

    for(int x=0;x<hGlauber->GetNbinsX();x++){
      for(int y=0;y<hGlauber->GetNbinsY();y++){
		float xx = hGlauber->GetXaxis()->GetBinCenter( x+1 );
        float yy = hGlauber->GetYaxis()->GetBinCenter( y+1 );
        float ww = hGlauber->GetBinContent( x+1, y+1 );

        meanX += xx*ww;
        meanY += yy*ww;
        meanX2 += xx*xx*ww;
        meanY2 += yy*yy*ww;
        meanXY += xx*yy*ww;
        sum += ww;
      }
    }

    meanX /= sum;
    meanY /= sum;
    meanX2 /= sum;
    meanY2 /= sum;
    meanXY /= sum;

    varXw = meanX2 - meanX*meanX;
    varYw = meanY2 - meanY*meanY;
    varXYw = meanXY - meanX*meanY;

    sTsq = varXw*varYw - varXYw*varXYw;
    sT = sTsq>=0 ? TMath::Pi()*sqrt( sTsq ) : -1;

    //Load histograms
    for (int it=0; it<ntimeHydro; it++){
      hTHydro[it] = (TH2D*)infileHydro->Get(Form("T_%d",it));
      timeHydro[irun-run_i][it] = htimeHydro->GetBinContent(it+1);
    }
    
    //EP angle from energy density profile
    int nBinsX = hTHydro[ntimeHydro-1]->GetNbinsX();
    int nBinsY = hTHydro[ntimeHydro-1]->GetNbinsY();
    double sumSinEnProf =0; double sumCosEnProf =0;
    for(int ibinx = 1; ibinx <= nBinsX; ibinx++){
      double bcx = hTHydro[ntimeHydro-1]->GetXaxis()->GetBinCenter(ibinx);
      for(int ibiny = 1; ibiny<=nBinsY; ibiny++){
        double bcy = hTHydro[ntimeHydro-1]->GetXaxis()->GetBinCenter(ibiny);
        double phiAng;
        if(bcx==0) phiAng = 0;
        else if(bcx>0 && bcy>=0) phiAng = TMath::ATan(abs(bcy/bcx));
        else if(bcx<0 && bcy>=0) phiAng = TMath::Pi() - TMath::ATan(abs(bcy/bcx));
        else if(bcx<0 && bcy<=0) phiAng = -TMath::Pi() + TMath::ATan(abs(bcy/bcx));
        else if(bcx>0 && bcy<=0) phiAng = -TMath::ATan(abs(bcy/bcx));
        else{cout << "ERROR : no phi angle calculation " << endl; return;}
        double bcw = TMath::Power(hTHydro[ntimeHydro-1]->GetBinContent(ibinx,ibiny),0.25);
        sumSinEnProf += bcw * TMath::Sin(phiN_*phiAng);
        sumCosEnProf += bcw * TMath::Cos(phiN_*phiAng);
        hTempEnergy[irun-run_i]->Fill(bcx,bcy,bcw);
      }
    }
    EPangEnProf = 1./phiN_ * TMath::ATan(sumSinEnProf/sumCosEnProf);

    //EP angle from Glauber 
    hist_glauber_fine[irun-run_i] = new TH2D(Form("hist_glauber_fineBins_irun%d",irun),";x;y",nSliceGlb, -max_x, max_x, nSliceGlb, -max_x, max_x);
    double sumSinGlb = 0; double sumCosGlb=0;

    if(!(Gncoll>0)){cout << "No Ncoll event!! " << endl; continue;}
    for(int isly=0; isly<nSliceGlb;isly++){
      for(int islx = 0; islx<nSliceGlb;islx++){
        double xval_ = -max_x + max_x/nSliceGlb + 2*max_x/nSliceGlb*islx;
        double yval_ = -max_x + max_x/nSliceGlb + 2*max_x/nSliceGlb*isly;
        long double content=0;

        for(int in=0; in<nNucl; in++){
          if(wproj[in]==0) continue;
          if(in >= nNucl_proj) continue;
          content += smearing_function->Eval(xproj[in] - xval_, yproj[in] - yval_);
        }

        for(int in=0; in<nNucl; in++){
          if(wtarg[in]==0) continue;
          if(in >= nNucl_targ) continue;
          content += smearing_function->Eval(xtarg[in] - xval_, ytarg[in] - yval_);
        }

        double phiAng;
        if(xval_==0) phiAng=0;
        else if(xval_>0 && yval_>=0) phiAng = TMath::ATan(abs(yval_/xval_));
        else if(xval_<0 && yval_>=0) phiAng = TMath::Pi() - TMath::ATan(abs(yval_/xval_));
        else if(xval_<0 && yval_<=0) phiAng = -TMath::Pi() + TMath::ATan(abs(yval_/xval_));
        else if(xval_>0 && yval_<=0) phiAng = -TMath::ATan(abs(yval_/xval_));
        else{cout << "ERROR : no phi angle calculation " << endl; return;}

        sumCosGlb += content*TMath::Cos(phiN_*phiAng);
        sumSinGlb += content*TMath::Sin(phiN_*phiAng);
      }
    }
    EPangGlauber  = 1./phiN_ * TMath::ATan(sumSinGlb/sumCosGlb);
   

    //EP angle from Glauber produced hist 
    int nBinsXglbh = hGlauber->GetNbinsX();
    int nBinsYglbh = hGlauber->GetNbinsY();
    double sumSinGlbh =0; double sumCosGlbh =0;
    for(int ibinx = 1; ibinx <= nBinsXglbh; ibinx++){
      double bcx = hGlauber->GetXaxis()->GetBinCenter(ibinx);
      for(int ibiny = 1; ibiny<=nBinsYglbh; ibiny++){
        double bcy = hGlauber->GetXaxis()->GetBinCenter(ibiny);
        double phiAng;
        if(bcx==0) phiAng = 0;
        else if(bcx>0 && bcy>=0) phiAng = TMath::ATan(abs(bcy/bcx));
        else if(bcx<0 && bcy>=0) phiAng = TMath::Pi() - TMath::ATan(abs(bcy/bcx));
        else if(bcx<0 && bcy<=0) phiAng = -TMath::Pi() + TMath::ATan(abs(bcy/bcx));
        else if(bcx>0 && bcy<=0) phiAng = -TMath::ATan(abs(bcy/bcx));
        else{cout << "ERROR : no phi angle calculation " << endl; return;}
        double bcw = hGlauber->GetBinContent(ibinx,ibiny);
        sumSinGlbh += bcw * TMath::Sin(phiN_*phiAng);
        sumCosGlbh += bcw * TMath::Cos(phiN_*phiAng);
      }
    }
    EPangGlauberh = 1./phiN_ * TMath::ATan(sumSinGlbh/sumCosGlbh);

/*
    cout << "EP Angle (Energy density profile) : " << EPangEnProf << endl;
    cout << "EP Angle (Glauber fine bins) : " << EPangGlauber << endl;
    cout << "EP Angle (Glauber hist) : " << EPangGlauberh << endl;
*/
    
    //Upsilon
    const int nY = nSAMP * Gncoll;
    cout << "nY : " << nY << endl;
    cout << endl;

    for(int s=0;s<nstates;s++){
      for (int iY=0; iY<nY; iY++){
        //Momentum
		fInitialUpsilon->SetParameters( 1.06450e+00 ,  7.97649e-01 , 100, const_mY[s] );
        double pT = fInitialUpsilon->GetRandom();
        double phi = (gRandom->Rndm()-0.5)*TMath::TwoPi(); 
        double px = pT*cos(phi);
        double py = pT*sin(phi);
        double mT = sqrt(pT*pT + const_mY[s]*const_mY[s]);

		int pTbin = int(pT);
		if ( pT>=20 ) continue;
		if( s>0 ) pTbin = int(pT/2);

      	//Beta
      	double bx = fP[s]->GetX(fabs(px)); 
      	double by = fP[s]->GetX(fabs(py)); 

      	if ( px<0 ) bx *= -1;
      	if ( py<0 ) by *= -1;

      	//Position
      	double vx, vy;
      	if(kInitPos==0){vx=0;vy=0;};
      	if(kInitPos==1) hGlauber->GetRandom2(vx, vy);
      	if(kInitPos==2){
          TH1D* hGlauberX = hGlauber->ProjectionX("hGlauberX",1,hGlauber->GetNbinsX());
          TH1D* hGlauberY = hGlauber->ProjectionY("hGlauberY",1,hGlauber->GetNbinsY());
          vx = hGlauberX->GetMean();
          vy = hGlauberY->GetMean();
      	}

      	TLorentzVector lvec;
      	lvec.SetPxPyPzE(px, py, 0, mT);

      	double tau_form = const_tau0Y[s]*lvec.Gamma();

      	double vx0 = vx;
      	double vy0 = vy;

      	double modF = 1.0;
      	int nFO = 0;
      	double det_flg = 1.0;

      	//Pre-Hydro 
      	if( bPreQGP && s==0 ){
          float TPre = PreHydroTempRatio*hTHydro[0]->GetBinContent(hTHydro[0]->FindBin(vx, vy))*1000.;
          if( TPre>Tf && tau_form<0.3 ){
            float GdissPre0 = gGdiss[s][pTbin]->Eval(TPre);
            float GdissPre1 = gGdiss[s][pTbin+1]->Eval(TPre);
            float GdissPre = GdissPre0 + (GdissPre1 - GdissPre0)*(pT - pTbin);
			float dt = 0.3 - tau_form;
            modF = exp(-(dt)*GdissPre/const_hbarc);
            if( exp(-(dt)*GdissPre/const_hbarc) < gRandom->Rndm() ) det_flg = 0.0;
          } else{
            modF = 1.0;
            nFO++;
          }
      	}

      	vx += bx*0.3;
      	vy += by*0.3;


      	//Time evolution
      	double totTime =0;
      	for (int it=0; it<ntimeHydro-1; it++){
          if ( nFO>=5 ) break;

          float dt = timeHydro[irun-run_i][it+1] - timeHydro[irun-run_i][it];
          float dx = bx*dt;
          float dy = by*dt;

		  if ( !bPreRES && timeHydro[irun-run_i][it]<tau_form ){
			vx += dx;
			vy += dy;
			continue;
		  }

          float THydro0 = hTHydro[it]->GetBinContent(hTHydro[it]->FindBin(vx, vy))*1000.;
          float THydro1 = hTHydro[it+1]->GetBinContent(hTHydro[it+1]->FindBin(vx+dx, vy+dy))*1000.;

		  if ( (THydro0+THydro1)/2>const_TmaxY[s] ){
			modF = 0.0;
			det_flg = 0.0;
			break;
		  }

          if ( (THydro0+THydro1)/2<Tf ){
            vx += dx;
            vy += dy;
            totTime += dt;
            nFO++;
            continue;
          }

		  float Gdiss0, Gdiss1, Gdiss, GdissRef;

		  if( s!=2 ){
            Gdiss0 = gGdiss[s][pTbin]->Eval((THydro0+THydro1)/2);
            Gdiss1 = gGdiss[s][pTbin+1]->Eval((THydro0+THydro1)/2);
            Gdiss = Gdiss0 + (Gdiss1 - Gdiss0)*(pT - pTbin*(s+1) );  // for 1s and 2s
		  } else if( s==2 ){
			GdissRef = fGdiss3S[0]->Eval((THydro0+THydro1)/2);
            Gdiss0 = GdissRef*(fGdiss2S[pTbin]->Eval((THydro0+THydro1)/2))/(fGdiss2S[0]->Eval((THydro0+THydro1)/2));
            Gdiss1 = GdissRef*(fGdiss2S[pTbin+1]->Eval((THydro0+THydro1)/2))/(fGdiss2S[0]->Eval((THydro0+THydro1)/2));
            Gdiss = Gdiss0 + (Gdiss1 - Gdiss0)*(pT - pTbin*2)/2.0;
		  }

		  if( bPreRES && timeHydro[irun-run_i][it]<tau_form ){
			Gdiss *= timeHydro[irun-run_i][it]/tau_form;
		  }

          modF *= exp(-(dt)*Gdiss/const_hbarc);
          if( exp(-(dt)*Gdiss/const_hbarc) < gRandom->Rndm() ) det_flg = 0.0;

          vx += dx;
          vy += dy;
          totTime += dt;

        }//it
/*
      	double velfx = (vx-vx0)/totTime;
      	double velfy = (vy-vy0)/totTime;
      	double pxf = getMom(velfx,const_mY[s]);
      	double pyf = getMom(velfy,const_mY[s]);

      	double gammax = getGamma(velfx); 
      	double gammay = getGamma(velfy); 

      	double vxf_ = (vx-vx0)*gammax;
      	double vyf_ = (vy-vy0)*gammay;

      	double vxf = vxf_ + vx0;
      	double vyf = vyf_ + vy0;
*/
/*
      	if( abs(px-pxf) > 1e-5 || abs(py-pyf) > 1e-5){
          cout << "Inconsistent momentum diff > 1e-5 GeV/c" << endl;
          cout << "initial pT, px, py : " << pT << "," << px << "," << py << endl;
          cout << "vx0, vy0 : " << vx0 << "," << vy0 << endl;
          cout << "vx, vy : " << vx << "," << vy << endl;
          cout << "total time : " << totTime << endl;
          cout << "final momentum px,py : " << pxf << "," << pyf << endl;
          cout << endl;
          return;
        }
*/

      	hprofRAA_Npart_rl[s]->Fill( Npart_, det_flg );
      	hprofRAA_mult_rl[s]->Fill( hMult->GetMean(), det_flg );
      	hprofRAA_multST_rl[s]->Fill( hMult->GetMean()/sT, det_flg );
      	hprofRAA_pT_rl[s]->Fill( pT, det_flg );

      	if ( b2D ){
          hprofRAA_xy[irun-run_i]->Fill(vx0, vy0, modF);
          hprofRAA_xy_rl[irun-run_i]->Fill(vx0, vy0, det_flg );
//	  hPosFinal[irun-run_i] -> Fill(vxf,vyf);
//	  hPosInit[irun-run_i] -> Fill(vx0,vy0);
      	}
      
      	//Fill Upsilon position
/*
      	double fphiAng;
      	if(pxf==0) fphiAng = 0;
      	else if(pxf>0 && pyf>=0) fphiAng = TMath::ATan(abs(pyf/pxf));
      	else if(pxf<0 && pyf>=0) fphiAng = TMath::Pi() - TMath::ATan(abs(pyf/pxf));
      	else if(pxf<0 && pyf<=0) fphiAng = -TMath::Pi() + TMath::ATan(abs(pyf/pxf));
      	else if(pxf>0 && pyf<=0) fphiAng = -TMath::ATan(abs(pyf/pxf));
      
      	hphi_initMom[s]->Fill(phi);
      	hphi_finalMom[s]->Fill(fphiAng);

      
      if( abs((fphiAng-fphiAng_)/fphiAng*100)>15){
      cout << "-----COMP START----" << endl;
      cout << "Upsilon : " << iY << " (px,py), (vx,vy) : " << pxf << ", " << pyf << " -- " << vx << ", " << vy << endl;
      cout << "fphiAng(px,py), fphiAng(vx,vy) : " << fphiAng << ", " << fphiAng_ << endl;
      cout << "py/px vs vy/vx : " << pyf/pxf << " - " << vy/vx << endl;
      cout << "--------------" << endl;
      cout << "initial pT, px, py : " << pT << "," << px << "," << py << endl;
      cout << "vx0, vy0 : " << vx0 << "," << vy0 << endl;
      cout << "vx, vy : " << vx << "," << vy << endl;
      cout << "total time : " << totTime << endl;
      cout << "final momentum px,py : " << pxf << "," << pyf << endl;
      cout << "initial phi : " << phi << endl;
      cout << endl;
      cout << " correction " << endl;
      cout << "lorentz factor x : " << 1/sqrt(1-vx/totTime*vx/totTime) << endl;
      cout << "lorentz factor y : " << 1/sqrt(1-vy/totTime*vy/totTime) << endl;
      double gammax = 1/sqrt(1-vx/totTime*vx/totTime);
      double gammay = 1/sqrt(1-vy/totTime*vy/totTime);
      cout << "corrected vx : " << vx*gammax << endl;
      cout << "corrected vy : " << vy*gammay << endl;
      cout << "corr vy/vx : " << (vy*gammay)/(vx*gammax) << endl;
      cout << endl;
      }
      */

/*
      	TLorentzVector* Ups4VRaw = new TLorentzVector;
      	TLorentzVector* Ups4VEnProfCor = new TLorentzVector;
      	TLorentzVector* Ups4VGlauberCor = new TLorentzVector;
      	Ups4VRaw->SetPtEtaPhiM(pT,0,phi,const_mY[s]);
      	Ups4VEnProfCor->SetPtEtaPhiM(Ups4VRaw->Pt(), Ups4VRaw->Eta(), Ups4VRaw->Phi()-EPangEnProf, Ups4VRaw->M());
      	Ups4VGlauberCor->SetPtEtaPhiM(Ups4VRaw->Pt(), Ups4VRaw->Eta(), Ups4VRaw->Phi()-EPangGlauber, Ups4VRaw->M());
      	new((*Ups4momRaw)[nUps])TLorentzVector(*Ups4VRaw);
      	new((*Ups4momEnProfCor)[nUps])TLorentzVector(*Ups4VEnProfCor);
      	new((*Ups4momGlauberCor)[nUps])TLorentzVector(*Ups4VGlauberCor);

      	vx0_[nUps] = vx0;
      	vy0_[nUps] = vy0;
      	IsUpsSurv_modf[nUps] = modF;
      	IsUpsSurv_prob[nUps] = det_flg;
      	nUps++;
*/
      }//iY
    }
    
    infileHydro->Close();
    delete infileHydro;

    tree->Fill();

    outfile->cd();
    if ( b2D ){
      hTempEnergy[irun-run_i]->Write();
    }
  }


/*
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
  }

  if ( bDRAW ){

    TCanvas *c2 = new TCanvas("c2","c2",1.2*2*500,500);
    c2->Divide(2,1);

    c2->cd(1);
    SetPadStyle();
    gPad->SetRightMargin(0.13);
    htmp = (TH1D*)gPad->DrawFrame(0,0,400,1.2);
    htmp->GetXaxis()->SetTitle("N_{part}");
    htmp->GetXaxis()->SetLabelSize(0.04);
    htmp->GetXaxis()->SetTitleSize(0.05);
    htmp->GetYaxis()->SetTitle("R_{AA} at p_{T}=3 GeV/c");
    htmp->GetYaxis()->SetLabelSize(0.04);
    htmp->GetYaxis()->SetTitleSize(0.05);

    hprofRAA_Npart->SetLineWidth(2);
    hprofRAA_Npart->Draw("same");


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

    hprofRAA_pT->SetLineWidth(2);
    hprofRAA_pT->Draw("same");

  }
*/

  if ( bSAVE ){
    outfile->cd();

    for(int s=0;s<nstates;s++){
      hprofRAA_Npart_rl[s]->Write();
      hprofRAA_mult_rl[s]->Write();
      hprofRAA_multST_rl[s]->Write();
      hprofRAA_pT_rl[s]->Write();

      hphi_initMom[s]->Write();
      hphi_finalMom[s]->Write();
    }

    if ( b2D ){
      for (int irun=run_i; irun<run_f; irun++){
        hprofRAA_xy[irun-run_i]->Write();
        hPosFinal[irun-run_i]->Write();
        hPosInit[irun-run_i]->Write();
      }
    }
    tree->Write();
  }
  return;
}
