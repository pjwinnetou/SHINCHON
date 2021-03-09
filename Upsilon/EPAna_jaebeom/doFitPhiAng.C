#include <iostream>
#include "TFile.h"
#include "TROOT.h"
#include "TF1.h"
#include "Style_jaebeom.h"

using namespace std;

void AssignFunction(Double_t CosOrd(Double_t *x, Double_t *p), int ir);

Double_t CosOrd2(Double_t *x, Double_t *p)
{ 
  Double_t xx = x[0];
  Double_t A  = p[0];
  Double_t v1 = p[1];
  Double_t v2 = p[2];
  Double_t out = A*(1+2*(v1*TMath::Cos(xx)+v2*TMath::Cos(2*xx)));
  return out;
} 

Double_t CosOrd3(Double_t *x, Double_t *p)
{ 
  Double_t xx = x[0];
  Double_t A  = p[0];
  Double_t v1 = p[1];
  Double_t v2 = p[2];
  Double_t v3 = p[3];
  Double_t out = A*(1+2*(v1*TMath::Cos(xx)+v2*TMath::Cos(2*xx)+v3*TMath::Cos(3*xx)));
  return out;
} 

Double_t CosOrd4(Double_t *x, Double_t *p)
{ 
  Double_t xx = x[0];
  Double_t A  = p[0];
  Double_t v1 = p[1];
  Double_t v2 = p[2];
  Double_t v3 = p[3];
  Double_t v4 = p[4];
  Double_t out = A*(1+2*(v1*TMath::Cos(xx)+v2*TMath::Cos(2*xx)+v3*TMath::Cos(3*xx)+v4*TMath::Cos(4*xx)));
  return out;
} 

Double_t CompOrd(Double_t *x, Double_t *p)
{ 
  Double_t xx = x[0];
  Double_t A  = p[0];
  Double_t v1 = p[1];
  Double_t out = A*(1+2*(v1*TMath::Cos(xx)));
  return out;
} 

  

void doFitPhiAng(int nrun = 10, int InitPos=1, const int fitorder = 3, int kDataSel = 0, int phiN=3)
{

  if(phiN==3 && fitorder<3){cout << "less fit order than phi" << endl; return;}

  const char* initPosStr;
  if(InitPos==0) initPosStr = "InitPosZero";
  else if(InitPos==1) initPosStr = "InitPosGlauber";
  else if(InitPos==2) initPosStr = "InitPosMean";

  const int nData = 3;
  enum fData {CMS_502, CMS_276, STAR_200,};
  const char* fDataStr[nData] = {"CMS502","CMS276","STAR200"};

  TFile* rf = new TFile(Form("res/SavedPlots_%s_PhiAng%d_%s_nRun%d.root",fDataStr[kDataSel],phiN,initPosStr,nrun),"read");

  TFile* wf = new TFile(Form("res/FitResults_v%d_fitorder%d_%s_%s_nRun%d.root",phiN,fitorder,fDataStr[kDataSel],initPosStr,nrun),"recreate");

  const int nPt = 4;
  const int nFitOrder = 3;
  TH1D* hPhiPt[nPt];
  
  TF1* f1[nPt];
  TF1* f1_v_comp[nFitOrder][nPt];

  TCanvas* c1[nPt];

  Double_t CosOrd(Double_t *x, Double_t *p);

  int npixel = 1300;
  const int nParComp = 2;
  double xmax = 3.2;
  double ymaxExt = 1.2;
  int fColor[nFitOrder+1] = {kBlue+1, kGreen+1, kMagenta+2, kOrange+2};
  int fStyle[nFitOrder+1] = {kDashed, kDotted, 10, 6};


  for(int ipt=0; ipt<nPt; ipt++){
    c1[ipt] = new TCanvas(Form("c_%d",ipt),"",600,600);
    hPhiPt[ipt] = (TH1D*) rf->Get(Form("hPhi_allEvt_EPGlauber_pt%d",ipt));
    SetHistStyle(hPhiPt[ipt],0,0);
    hPhiPt[ipt]->SetMinimum(0);
    hPhiPt[ipt]->SetMaximum(hPhiPt[ipt]->GetMaximum()*ymaxExt);
   
    //AssignFunction(CosOrd,fitorder-2);
    if(fitorder==2) f1[ipt] = new TF1(Form("fit_pt%d",ipt),CosOrd2,-xmax,xmax,fitorder+1);
    else if(fitorder==3) f1[ipt] = new TF1(Form("fit_pt%d",ipt),CosOrd3,-xmax,xmax,fitorder+1);
    else if(fitorder==4) f1[ipt] = new TF1(Form("fit_pt%d",ipt),CosOrd4,-xmax,xmax,fitorder+1);
    f1[ipt]->SetLineColor(kRed-4);
    f1[ipt]->SetLineWidth(3);
    f1[ipt]->SetNpx(npixel);
    
    for(int ifit=0;ifit<fitorder;ifit++){
      f1_v_comp[ifit][ipt] = new TF1(Form("fit_v%dcompt_pt%d",ifit,ipt),CompOrd,-xmax,xmax,nParComp);
    }
    
    hPhiPt[ipt]->Fit(Form("fit_pt%d",ipt),"R");

    c1[ipt]->cd();
    wf->cd();
    hPhiPt[ipt]->Draw("pe");
    f1[ipt]->Draw("same");

    hPhiPt[ipt]->Write();
    f1[ipt]->Write();
    for(int ifit=0;ifit<fitorder;ifit++){
      f1_v_comp[ifit][ipt] -> FixParameter(0, f1[ipt]->GetParameter(0));
      f1_v_comp[ifit][ipt] -> FixParameter(1, f1[ipt]->GetParameter(ifit+1));
      f1_v_comp[ifit][ipt] -> SetLineColor(fColor[ifit]);
      f1_v_comp[ifit][ipt] -> SetLineStyle(fStyle[ifit]);
      f1_v_comp[ifit][ipt] -> SetLineWidth(2);
      f1_v_comp[ifit][ipt]->SetNpx(npixel);
      f1_v_comp[ifit][ipt]->Draw("same");
      f1_v_comp[ifit][ipt]->Write();
    }
  }
    
  //integrated pt
  TCanvas* c1_int = new TCanvas("c_int","",600,600);
  TH1D* hPhiPt_int = (TH1D*) rf->Get("hPhi_allEvt_EPGlauber");
  SetHistStyle(hPhiPt_int,0,0);
  hPhiPt_int->SetMinimum(0);
  hPhiPt_int->SetMaximum(hPhiPt_int->GetMaximum()*ymaxExt);

  TF1* f1_int; 
  //AssignFunction(CosOrd,fitorder-2);
  if(fitorder==2) f1_int = new TF1("fit_ptInt",CosOrd2,-xmax,xmax,fitorder+1);
  else if(fitorder==3) f1_int = new TF1("fit_ptInt",CosOrd3,-xmax,xmax,fitorder+1);
  else if(fitorder==4) f1_int = new TF1("fit_ptInt",CosOrd4,-xmax,xmax,fitorder+1);
  f1_int->SetLineColor(kRed-4);
  f1_int->SetLineWidth(3);
  f1_int->SetNpx(npixel);

  TF1* f1_v_comp_int[nFitOrder];
  for(int ifit=0;ifit<nFitOrder;ifit++){
    f1_v_comp_int[ifit] = new TF1(Form("fit_v%dcompt_ptInt",ifit),CompOrd,-xmax,xmax,nParComp);
  }

  hPhiPt_int->Fit("fit_ptInt","R");

  c1_int->cd();
  hPhiPt_int->Draw("pe");
  f1_int->Draw("same");

  hPhiPt_int->Write();
  f1_int->Write();
  for(int ifit=0;ifit<fitorder;ifit++){
    f1_v_comp_int[ifit] -> FixParameter(0, f1_int->GetParameter(0));
    f1_v_comp_int[ifit] -> FixParameter(1, f1_int->GetParameter(ifit+1));
    f1_v_comp_int[ifit] -> SetLineColor(fColor[ifit]);
    f1_v_comp_int[ifit] -> SetLineStyle(fStyle[ifit]);
    f1_v_comp_int[ifit] -> SetLineWidth(2);
    f1_v_comp_int[ifit]->SetNpx(npixel);
    f1_v_comp_int[ifit]->Draw("same");
    f1_v_comp_int[ifit]->Write();
  }
}
  
void AssignFunction(Double_t CosOrd(Double_t *x, Double_t *p), int ir){
  if(ir==0) CosOrd = CosOrd2;
  else if(ir==1) CosOrd = CosOrd3;
  else if(ir==2) CosOrd = CosOrd4;
};
