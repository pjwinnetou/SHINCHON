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

#include <iostream>
#include <fstream>

using namespace std;

Double_t fTsallis1S_v2(Double_t *x, Double_t *fpar);

double getMom(double v, double m)
{
  double mom = m*v / sqrt(1-v*v);
  return mom;
}

double getGamma(double v)
{
  double gamma = 1 / sqrt(1-v*v);
  return gamma;
}

void makeUpsTree(const int kInitPos = 1){

  gStyle->SetOptStat(0);
  gStyle->SetPalette(55);

  gRandom = new TRandom3(0);

  const bool bDRAW = true; 
  const bool bSAVE = true; 
  const bool b2D = true;

  const int run_i = 0;
  const int run_f = 10;
  const int nrun = (run_f-run_i);
  const float const_hbarc = 197.5; //MeV fm
  const float const_mY = 9.46; //GeV
  const int nSAMP = 100; //times Ncoll
  double sigs = 0.5;

  ifstream fdata;

  char buf[500];
  vector<float> T[21];
  vector<float> Gdiss[21];

  float f_tmp[20];

  fdata.open("../Gdiss0.dat");
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

  fdata.open("../Gdiss1.dat");
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

  for (int ipt=0; ipt<21; ipt++){
    gGdiss[ipt] = new TGraphErrors;
    for (unsigned int iT=0; iT<T[ipt].size(); iT++){
      gGdiss[ipt]->SetPoint(iT, T[ipt][iT], Gdiss[ipt][iT]);
    }
  }

  TCanvas *c0;
  if ( bDRAW ){
    c0 = new TCanvas("c0","c0",1.2*500,500);
    SetPadStyle();

    htmp = (TH1D*)gPad->DrawFrame(170,0,500,200);
    SetHistoStyle("T [MeV]","#Gamma_{diss} [MeV]");

    for (int ipt=0; ipt<21; ipt++){
      gGdiss[ipt]->SetLineWidth(3);
      gGdiss[ipt]->SetLineColorAlpha(kGray,0.5);
      if ( ipt==0 ){
        gGdiss[ipt]->SetLineColor(2);
      }else if ( ipt==5 ){
        gGdiss[ipt]->SetLineColor(kGreen+2);
      }else if ( ipt==10 ){
        gGdiss[ipt]->SetLineColor(kBlue);
      }else if ( ipt==15 ){
        gGdiss[ipt]->SetLineColor(kMagenta);
      }else if ( ipt==20 ){
        gGdiss[ipt]->SetLineColor(1);
      }
      gGdiss[ipt]->Draw("C");
    }

    TLegend *leg = new TLegend(0.2,0.6,0.6,0.9);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(43);
    leg->SetTextSize(20);
    leg->AddEntry(gGdiss[0],"p_{T}=0 GeV/c","L");
    leg->AddEntry(gGdiss[5],"p_{T}=5 GeV/c","L");
    leg->AddEntry(gGdiss[10],"p_{T}=10 GeV/c","L");
    leg->AddEntry(gGdiss[15],"p_{T}=15 GeV/c","L");
    leg->AddEntry(gGdiss[20],"p_{T}=20 GeV/c","L");
    leg->Draw();
  }


  //Upsilon beta vs. p 
  TF1 *fP = new TF1("fP","[0]*x/sqrt(1-x*x)",0,1);
  fP->SetParameter(0, const_mY);

  TCanvas *c1_;
  if ( bDRAW ){
    c1_ = new TCanvas("c1_","c1_",1.2*500,500);
    SetPadStyle();
    htmp = (TH1D*)gPad->DrawFrame(0,0,1,25);
    SetHistoStyle("#beta","p [GeV]");

    fP->SetLineColor(1);
    fP->SetLineWidth(3);
    fP->Draw("same");
  }

  //return;

  //Glauber
  const int Necc = 10;
  const int nNucl = 500;
  TFile *infileGlauber = new TFile("../MCGlauber-PbPb-5020GeV-b0-18fm.root","read");
  TTree *TGlauber = (TTree*)infileGlauber->Get("lemon");
  int Gnpart, Gncoll;
  float b, eccgaus[Necc], eccpoint[Necc];
  float xproj[nNucl], yproj[nNucl], xtarg[nNucl], ytarg[nNucl];
  TGlauber->SetBranchAddress("npart",&Gnpart);
  TGlauber->SetBranchAddress("ncoll",&Gncoll);
  TGlauber->SetBranchAddress("b",&b);
  TGlauber->SetBranchAddress("eccgaus",eccgaus);
  TGlauber->SetBranchAddress("eccpoint",eccpoint);
  TGlauber->SetBranchAddress("xproj",xproj);
  TGlauber->SetBranchAddress("yproj",yproj);
  TGlauber->SetBranchAddress("xtarg",xtarg);
  TGlauber->SetBranchAddress("ytarg",ytarg);

  //QGP T profile
  TH2D *hTHydro[300];
  TGraphErrors *gTHydro[nrun];
  TF1 *fT2[nrun];
  double timeHydro[nrun][300] = {0.0}; 
  float errTHydro[nrun][300] = {0.0};
  float freezeT[nrun];

  TProfile *hprofRAA_Npart = new TProfile("hprofRAA_Npart","",40,0,400);
  TProfile *hprofRAA_pT = new TProfile("hprofRAA_pT","",100,0,20);
  TProfile2D *hprofRAA_xy[nrun];


  TProfile *hprofRAA_Npart_rl = new TProfile("hprofRAA_Npart_rl","",40,0,400);
  TProfile *hprofRAA_pT_rl = new TProfile("hprofRAA_pT_rl","",100,0,20);
  TProfile2D *hprofRAA_xy_rl[nrun];

  //**** temporary bin definition.
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

  TH2D* hTempEnergy[nrun];

  // cross check histogram 
  TH1D* hphi_initMom   = new TH1D("hphi_initMom","",100,-TMath::Pi(),TMath::Pi());
  TH1D* hphi_finalMom  = new TH1D("hphi_finalMom","",100,-TMath::Pi(),TMath::Pi());
  
  TH2D* hPosFinal[nrun];
  TH2D* hPosInit[nrun]; 


  //***** Tree & Output Definition
  string fInitPos;
  if(kInitPos==0) fInitPos = "InitPosZero";
  else if(kInitPos==1) fInitPos = "InitPosGlauber";
  else if(kInitPos==2) fInitPos = "InitPosMean";
  TFile *outfile = new TFile(Form("./outfile_UpsSkim_%s_%04d_%04d.root",fInitPos.c_str(),run_i,run_f),"recreate");

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
  TF1 *fInitialUpsilon = new TF1("fInitialUpsilon",fTsallis1S_v2,0,30,3);
  fInitialUpsilon -> SetParameters(  1.06450e+00 ,  7.97649e-01 , 100);
  
  //Glauber smearing function  
  TF2* smearing_function = new TF2("smear_tf2", "TMath::Exp(-(x*x+y*y)/(2.*[0]*[0]))/(2*TMath::Pi()*[0]*[0])", 0, 10*sigs, 0, 10*sigs);
  smearing_function->FixParameter(0,sigs);


  for (int irun=run_i; irun<run_f; irun++){

    cout << "Scan event #" << irun << endl;
    hTempEnergy[irun]= new TH2D(Form("histTempEnergy_%d",irun),";x;y",100,-15,15,100,-15,15);
    
    runid = irun;
    Ups4momRaw->Clear();
    Ups4momEnProfCor->Clear();
    Ups4momGlauberCor->Clear();
    nUps=0;

    if ( b2D ){
      hprofRAA_xy[irun] = new TProfile2D(Form("hprofRAA_xy_run%05d",irun),"",100,-15,15,100,-15,15);
      hprofRAA_xy_rl[irun] = new TProfile2D(Form("hprofRAA_xy_rl_run%05d",irun),"",100,-15,15,100,-15,15);
      hPosFinal[irun] = new TH2D(Form("hPosFinal_run%d",irun),"",100,-15,15,100,-15,15);
      hPosInit[irun] = new TH2D(Form("hPosInit_run%d",irun),"",100,-15,15,100,-15,15);
    }

    TFile *infileHydro = new TFile(Form("../SONIC_profile_PbPb5TeV_0_18fm_event%05d.root",irun),"read");
    TH1D *htimeHydro = (TH1D*)infileHydro->Get("Time");
    int ntimeHydro = (int)htimeHydro->GetEntries();

    //Glauber info
    TGlauber->GetEntry(irun);
    Npart_ = Gnpart;
    Ncoll_ = Gncoll;
    b_ = b;
    for(int iecc=0;iecc<Necc;iecc++){eccgaus_[iecc] = eccgaus[iecc]; eccpoint_[iecc] = eccpoint[iecc];}

    TH2D *hGlauber = (TH2D*)infileGlauber->Get(Form("inited_event%d",irun));

    //Load histograms
    for (int it=0; it<ntimeHydro; it++){
      hTHydro[it] = (TH2D*)infileHydro->Get(Form("T_%d",it));
      timeHydro[irun][it] = htimeHydro->GetBinContent(it+1);
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
        sumSinEnProf += bcw * TMath::Sin(2*phiAng);
        sumCosEnProf += bcw * TMath::Cos(2*phiAng);
        hTempEnergy[irun]->Fill(bcx,bcy,bcw);
      }
    }
    EPangEnProf = 0.5 * TMath::ATan(sumSinEnProf/sumCosEnProf);



    //EP angle from Glauber 
    double lowcut = 10e-5; double highcut = 10000;
    const Double_t max_x = 15.0;
    const int nSliceGlb = 300;
    double sumSinGlb = 0; double sumCosGlb=0;

    if(!(Gncoll>0)){cout << "No Ncoll event!! " << endl; continue;}
    for(int islx = 0; islx<nSliceGlb;islx++){
      double xval_ = -max_x + max_x/nSliceGlb + 2*max_x/nSliceGlb*islx;
      for(int isly=0; isly<nSliceGlb;isly++){
        double yval_ = -max_x + max_x/nSliceGlb + 2*max_x/nSliceGlb*isly;
        long double content=0;

        for(int in=0; in<nNucl; in++){
          if(xproj[in]==0 || yproj[in]==0 || abs(xproj[in])<lowcut || abs(xproj[in])>highcut || abs(yproj[in])<lowcut || abs(yproj[in]) > highcut) continue;
          content += smearing_function->Eval(xproj[in] - xval_, yproj[in] - yval_);
        }

        for(int in=0; in<nNucl; in++){
          if(xtarg[in]==0 || ytarg[in]==0 || abs(xtarg[in])<lowcut || abs(xtarg[in])>highcut || abs(ytarg[in])<lowcut || abs(ytarg[in]) > highcut) continue;
          content += smearing_function->Eval(xtarg[in] - xval_, ytarg[in] - yval_);
        }


        double phiAng;
        if(xval_==0) phiAng=0;
        else if(xval_>0 && yval_>=0) phiAng = TMath::ATan(abs(yval_/xval_));
        else if(xval_<0 && yval_>=0) phiAng = TMath::Pi() - TMath::ATan(abs(yval_/xval_));
        else if(xval_<0 && yval_<=0) phiAng = -TMath::Pi() + TMath::ATan(abs(yval_/xval_));
        else if(xval_>0 && yval_<=0) phiAng = -TMath::ATan(abs(yval_/xval_));
        else{cout << "ERROR : no phi angle calculation " << endl; return;}

        sumCosGlb += content*TMath::Cos(2*phiAng);
        sumSinGlb += content*TMath::Sin(2*phiAng);

      }
    }
    EPangGlauber  = 0.5 * TMath::ATan(sumSinGlb/sumCosGlb);
   


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
        sumSinGlbh += bcw * TMath::Sin(2*phiAng);
        sumCosGlbh += bcw * TMath::Cos(2*phiAng);
      }
    }
    EPangGlauberh = 0.5 * TMath::ATan(sumSinGlbh/sumCosGlbh);

    cout << "EP Angle (Energy density profile) : " << EPangEnProf << endl;
    cout << "EP Angle (Glauber fine bins) : " << EPangGlauber << endl;
    cout << "EP Angle (Glauber hist) : " << EPangGlauberh << endl;

		const int nY = nSAMP * Gncoll;
    cout << "nY : " << nY << endl;
    cout << endl;

    //Upsilon
    for (int iY=0; iY<nY; iY++){

      //Momentum
      //double pT = 20.0*gRandom->Rndm();
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

      //Position
      double vx, vy;
      if(kInitPos==0){vx=0;vy=0;};
      if(kInitPos==1) hGlauber->GetRandom2(vx, vy);
      else if(kInitPos==2){
        TH1D* hGlauberX = hGlauber->ProjectionX("hGlauberX",1,hGlauber->GetNbinsX());
        TH1D* hGlauberY = hGlauber->ProjectionY("hGlauberY",1,hGlauber->GetNbinsY());
        vx = hGlauberX->GetMean();
        vy = hGlauberY->GetMean();
      }

      double vx0 = vx;
      double vy0 = vy;

      double modF = 1.0;
      int nFO = 0;

      double det_flg = 1.0;

      //Pre-Hydro 
      {
        float TPre = 1.20*hTHydro[0]->GetBinContent(hTHydro[0]->FindBin(vx, vy))*1000.;
        if ( TPre>170.0 ){
          float GdissPre0 = gGdiss[pTbin]->Eval(TPre);
          float GdissPre1 = gGdiss[pTbin+1]->Eval(TPre);
          float GdissPre = GdissPre0 + (GdissPre1 - GdissPre0)*(pT - pTbin);
          modF = exp(-(0.2)*GdissPre/const_hbarc);
          if( exp(-(0.2)*GdissPre/const_hbarc) < gRandom->Rndm() ) det_flg = 0.0;
        }else{
          modF = 1.0;
          nFO++;
        }
      }

      //Time evolution
      double totTime =0;
      for (int it=0; it<ntimeHydro-1; it++){

        if ( nFO>=10 ) break;

        float dt = timeHydro[irun][it+1] - timeHydro[irun][it];
        float dx = bx*dt;
        float dy = by*dt;

        float THydro0 = hTHydro[it]->GetBinContent(hTHydro[it]->FindBin(vx, vy))*1000.;
        float THydro1 = hTHydro[it+1]->GetBinContent(hTHydro[it+1]->FindBin(vx+dx, vy+dy))*1000.;

        if ( (THydro0+THydro1)/2<170.0 ){
          vx += dx;
          vy += dy;
          totTime += dt;
          nFO++;
          continue;
        }

        float Gdiss0 = gGdiss[pTbin]->Eval((THydro0+THydro1)/2);
        float Gdiss1 = gGdiss[pTbin+1]->Eval((THydro0+THydro1)/2);
        float Gdiss = Gdiss0 + (Gdiss1 - Gdiss0)*(pT - pTbin); 
        modF *= exp(-(dt)*Gdiss/const_hbarc);
        if( exp(-(dt)*Gdiss/const_hbarc) < gRandom->Rndm() ) det_flg = 0.0;

        vx += dx;
        vy += dy;
        totTime += dt;

      }//it

      double velfx = (vx-vx0)/totTime;
      double velfy = (vy-vy0)/totTime;
      double pxf = getMom(velfx,const_mY);
      double pyf = getMom(velfy,const_mY);

      double gammax = getGamma(velfx); 
      double gammay = getGamma(velfy); 

      double vxf_ = (vx-vx0)*gammax;
      double vyf_ = (vy-vy0)*gammay;

      double vxf = vxf_ + vx0;
      double vyf = vyf_ + vy0;


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

      hprofRAA_Npart->Fill(Npart_, modF);
      hprofRAA_Npart_rl->Fill(Npart_, det_flg );
      if ( b2D ){
        hprofRAA_xy[irun]->Fill(vx0, vy0, modF);
        hprofRAA_xy_rl[irun]->Fill(vx0, vy0, det_flg );
      }

      hprofRAA_pT->Fill(pT, modF);
      hprofRAA_pT_rl->Fill(pT, det_flg );
      
      //Fill Upsilon position
      hPosFinal[irun] -> Fill(vxf,vyf);   
      hPosInit[irun] -> Fill(vx0,vy0);   
  

      double fphiAng;
      if(pxf==0) fphiAng = 0;
      else if(pxf>0 && pyf>=0) fphiAng = TMath::ATan(abs(pyf/pxf));
      else if(pxf<0 && pyf>=0) fphiAng = TMath::Pi() - TMath::ATan(abs(pyf/pxf));
      else if(pxf<0 && pyf<=0) fphiAng = -TMath::Pi() + TMath::ATan(abs(pyf/pxf));
      else if(pxf>0 && pyf<=0) fphiAng = -TMath::ATan(abs(pyf/pxf));
      
      hphi_initMom->Fill(phi);
      hphi_finalMom->Fill(fphiAng);

      /*
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

      TLorentzVector* Ups4VRaw = new TLorentzVector;
      TLorentzVector* Ups4VEnProfCor = new TLorentzVector;
      TLorentzVector* Ups4VGlauberCor = new TLorentzVector;
      Ups4VRaw->SetPtEtaPhiM(pT,0,phi,const_mY);
      Ups4VEnProfCor->SetPtEtaPhiM(Ups4VRaw->Pt(), Ups4VRaw->Eta(), Ups4VRaw->Phi()-EPangEnProf, Ups4VRaw->M());
      Ups4VGlauberCor->SetPtEtaPhiM(Ups4VRaw->Pt(), Ups4VRaw->Eta(), Ups4VRaw->Phi()-EPangGlauberh + TMath::Pi()/2, Ups4VRaw->M());
      new((*Ups4momRaw)[nUps])TLorentzVector(*Ups4VRaw);
      new((*Ups4momEnProfCor)[nUps])TLorentzVector(*Ups4VEnProfCor);
      new((*Ups4momGlauberCor)[nUps])TLorentzVector(*Ups4VGlauberCor);

      vx0_[nUps] = vx0;
      vy0_[nUps] = vy0;
      IsUpsSurv_modf[nUps] = modF;
      IsUpsSurv_prob[nUps] = det_flg;
      nUps++;
    }//iY

    infileHydro->Close();
    delete infileHydro;

    tree->Fill();

    outfile->cd();
    hTempEnergy[irun]->Write();
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

  if ( bSAVE ){

    outfile->cd();
    hprofRAA_Npart->Write();
    hprofRAA_pT->Write();

    if ( b2D ){
      for (int irun=run_i; irun<run_f; irun++){
        hprofRAA_xy[irun]->Write();
        hPosFinal[irun]->Write();
        hPosInit[irun]->Write();
      }
    }
    tree->Write();
    hphi_initMom->Write();
    hphi_finalMom->Write();
  }

  return;
}

Double_t fTsallis1S_v2(Double_t *x, Double_t *fpar)
{
  Float_t xx = x[0];
  Double_t Y1Smass = 10.023;
  Double_t q = fpar[0];
  Double_t T = fpar[1];
  Double_t c = fpar[2];
  Double_t mT = TMath::Sqrt(Y1Smass*Y1Smass+xx*xx);
  Double_t pow = TMath::Power((1+(q-1)*mT/T),(-q/(q-1)));

  Double_t f = c*mT*xx*pow;
  return f;
}


