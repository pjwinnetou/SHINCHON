#include <iostream>
#include "TFile.h"
#include "TROOT.h"
#include "TF1.h"
#include "Style_jaebeom.h"

using namespace std;

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

Double_t CompOrd1(Double_t *x, Double_t *p)
{ 
  Double_t xx = x[0];
  Double_t A  = p[0];
  Double_t v1 = p[1];
  Double_t out = A*(1+2*(v1*TMath::Cos(xx)));
  return out;
} 

Double_t CompOrd2(Double_t *x, Double_t *p)
{ 
  Double_t xx = x[0];
  Double_t A  = p[0];
  Double_t v1 = p[1];
  Double_t out = A*(1+2*(v1*TMath::Cos(2*xx)));
  return out;
} 

Double_t CompOrd3(Double_t *x, Double_t *p)
{ 
  Double_t xx = x[0];
  Double_t A  = p[0];
  Double_t v1 = p[1];
  Double_t out = A*(1+2*(v1*TMath::Cos(3*xx)));
  return out;
} 


void doFitV2(int nrun = 10, int InitPos=1, const int fitorder = 3)
{
  const char* initPosStr;
  if(InitPos==0) initPosStr = "InitPosZero";
  else if(InitPos==1) initPosStr = "InitPosGlauber";
  else if(InitPos==2) initPosStr = "InitPosMean";

  TFile* rf = new TFile(Form("res/SavedPlots_%s_nRun%d.root",initPosStr,nrun),"read");

  TFile* wf = new TFile(Form("res/FitResults_v2_fitorder%d_%s_nRun%d.root",fitorder,initPosStr,nrun),"recreate");

  const int nPt = 4;
  TH1D* hPhiPt[nPt];
  
  TF1* f1[nPt];
  TF1* f1_v1comp[nPt];
  TF1* f1_v2comp[nPt];
  TF1* f1_v3comp[nPt];

  TCanvas* c1[nPt];

  int npixel = 1300;

  for(int ipt=0; ipt<nPt; ipt++){
    c1[ipt] = new TCanvas(Form("c_%d",ipt),"",600,600);
    hPhiPt[ipt] = (TH1D*) rf->Get(Form("hPhi_allEvt_EPGlauber_pt%d",ipt));
    SetHistStyle(hPhiPt[ipt],0,0);
    hPhiPt[ipt]->SetMinimum(0);
    hPhiPt[ipt]->SetMaximum(hPhiPt[ipt]->GetMaximum()*1.2);
    
    if(fitorder==2){ 
      f1[ipt] = new TF1(Form("fit_pt%d",ipt),CosOrd2,-3.2,3.2,fitorder+1);
      f1_v1comp[ipt] = new TF1(Form("fit_v1compt_pt%d",ipt),CompOrd1,-3.2,3.2,2);
      f1_v2comp[ipt] = new TF1(Form("fit_v2compt_pt%d",ipt),CompOrd2,-3.2,3.2,2);
    }
    else if(fitorder==3){
      f1[ipt] = new TF1(Form("fit_pt%d",ipt),CosOrd3,-3.2,3.2,fitorder+1);
      f1_v1comp[ipt] = new TF1(Form("fit_v1compt_pt%d",ipt),CompOrd1,-3.2,3.2,2);
      f1_v2comp[ipt] = new TF1(Form("fit_v2compt_pt%d",ipt),CompOrd2,-3.2,3.2,2);
      f1_v3comp[ipt] = new TF1(Form("fit_v3compt_pt%d",ipt),CompOrd3,-3.2,3.2,2);
    }
    f1[ipt]->SetLineColor(kRed-4);
    f1[ipt]->SetLineWidth(3);
    f1[ipt]->SetNpx(npixel);

    hPhiPt[ipt]->Fit(Form("fit_pt%d",ipt),"R");

    f1_v1comp[ipt]->FixParameter(0, f1[ipt]->GetParameter(0));
    f1_v1comp[ipt]->FixParameter(1, f1[ipt]->GetParameter(1));
    f1_v2comp[ipt]->FixParameter(0, f1[ipt]->GetParameter(0));
    f1_v2comp[ipt]->FixParameter(1, f1[ipt]->GetParameter(2));
    if(fitorder==3){
      f1_v3comp[ipt]->FixParameter(0,f1[ipt]->GetParameter(0));
      f1_v3comp[ipt]->FixParameter(1,f1[ipt]->GetParameter(3));
      f1_v3comp[ipt]->SetLineColor(kMagenta+2);
      f1_v3comp[ipt]->SetLineWidth(2);
      f1_v3comp[ipt]->SetLineStyle(10);
      f1_v3comp[ipt]->SetNpx(npixel);
    }

    f1_v1comp[ipt]->SetLineColor(kBlue+1);
    f1_v1comp[ipt]->SetLineWidth(2);
    f1_v1comp[ipt]->SetLineStyle(kDashed);
    f1_v1comp[ipt]->SetNpx(npixel);

    f1_v2comp[ipt]->SetLineColor(kGreen+1);
    f1_v2comp[ipt]->SetLineWidth(2);
    f1_v2comp[ipt]->SetLineStyle(kDotted);
    f1_v2comp[ipt]->SetNpx(npixel);


    c1[ipt]->cd();
    wf->cd();
    hPhiPt[ipt]->Draw("pe");
    f1[ipt]->Draw("same");
    f1_v1comp[ipt]->Draw("same");
    f1_v2comp[ipt]->Draw("same");

    hPhiPt[ipt]->Write();
    f1[ipt]->Write();
    f1_v1comp[ipt]->Write();
    f1_v2comp[ipt]->Write();
    if(fitorder==3){
      f1_v3comp[ipt]->Draw("same");
      f1_v3comp[ipt]->Write();
    }
  }
    
  //integrated pt
  TCanvas* c1_int = new TCanvas("c_int","",600,600);
  TH1D* hPhiPt_int = (TH1D*) rf->Get("hPhi_allEvt_EPGlauber");
  SetHistStyle(hPhiPt_int,0,0);
  hPhiPt_int->SetMinimum(0);
  hPhiPt_int->SetMaximum(hPhiPt_int->GetMaximum()*1.2);

  TF1* f1_int; TF1* f1_v1comp_int; TF1* f1_v2comp_int; TF1* f1_v3comp_int;
  if(fitorder==2){ 
    f1_int = new TF1("fit_ptInt",CosOrd2,-3.2,3.2,fitorder+1);
    f1_v1comp_int = new TF1("fit_v1compt_ptInt",CompOrd1,-3.2,3.2,2);
    f1_v2comp_int = new TF1("fit_v2compt_ptInt",CompOrd2,-3.2,3.2,2);
  }
  else if(fitorder==3){
    f1_int = new TF1("fit_ptInt",CosOrd3,-3.2,3.2,fitorder+1);
    f1_v1comp_int = new TF1("fit_v1compt_ptInt",CompOrd1,-3.2,3.2,2);
    f1_v2comp_int = new TF1("fit_v2compt_ptInt",CompOrd2,-3.2,3.2,2);
    f1_v3comp_int = new TF1("fit_v3compt_ptInt",CompOrd3,-3.2,3.2,2);
  }
  f1_int->SetLineColor(kRed-4);
  f1_int->SetLineWidth(3);
  f1_int->SetNpx(npixel);

  hPhiPt_int->Fit("fit_ptInt","R");

  f1_v1comp_int->FixParameter(0, f1_int->GetParameter(0));
  f1_v1comp_int->FixParameter(1, f1_int->GetParameter(1));
  f1_v2comp_int->FixParameter(0, f1_int->GetParameter(0));
  f1_v2comp_int->FixParameter(1, f1_int->GetParameter(2));
  if(fitorder==3){
    f1_v3comp_int->FixParameter(0,f1_int->GetParameter(0));
    f1_v3comp_int->FixParameter(1,f1_int->GetParameter(3));
    f1_v3comp_int->SetLineColor(kMagenta+2);
    f1_v3comp_int->SetLineWidth(2);
    f1_v3comp_int->SetLineStyle(10);
    f1_v3comp_int->SetNpx(npixel);
  }

  f1_v1comp_int->SetLineColor(kBlue+1);
  f1_v1comp_int->SetLineWidth(2);
  f1_v1comp_int->SetLineStyle(kDashed);
  f1_v1comp_int->SetNpx(npixel);

  f1_v2comp_int->SetLineColor(kGreen+1);
  f1_v2comp_int->SetLineWidth(2);
  f1_v2comp_int->SetLineStyle(kDotted);
  f1_v2comp_int->SetNpx(npixel);


  c1_int->cd();
  hPhiPt_int->Draw("pe");
  f1_int->Draw("same");
  f1_v1comp_int->Draw("same");
  f1_v2comp_int->Draw("same");

  hPhiPt_int->Write();
  f1_int->Write();
  f1_v1comp_int->Write();
  f1_v2comp_int->Write();

  if(fitorder==3){
    f1_v3comp_int->Draw("same");
    f1_v3comp_int->Write();
  }

}
