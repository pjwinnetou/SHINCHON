#include <iostream>
#include "../Style_jaebeom.h"

using namespace std;

void draw_vn_pt_compare(int fitor = 3, int nrun = 10, int kInitPos=1, int kDataSel=0,bool isLine = false, int phiN=3){

  setTDRStyle();
  writeExtraText= false;
  int iPeriod = -1;
  int iPos = 33;
  bool drawInner = true;

  gStyle->SetEndErrorSize(0);
  gStyle->SetOptFit(kTRUE);
 
  const char* fPosStr;
  if(kInitPos==0) fPosStr = "InitPosZero";
  else if(kInitPos==1) fPosStr = "InitPosGlauber";
  else if(kInitPos==2) fPosStr = "InitPosMean";
  
  const int nData = 3;
  enum fData {CMS_502, CMS_276, STAR_200,};
  const char* fDataStr[nData] = {"CMS502","CMS276","STAR200"};

  const int nPt = 4;
  double ptBin[nPt+1] = {0, 3, 6, 12, 20};
  double exsys[nPt] =  {1.5, 1.5, 3, 4};
  double ptBin_m[nPt] = {1.5, 4.5, 9, 16};

  double ymin = -0.005; double ymax = 0.023; double xmin = 0; double xmax = 20;

  //// read input file : value & stat.
  // Fit value
  TFile *fvnres = new TFile(Form("../res/FitResults_v%d_fitorder%d_%s_%s_nRun%d.root",phiN,fitor,fDataStr[kDataSel],fPosStr,nrun),"read");
  TFile *fsv = new TFile(Form("../res/SavedPlots_%s_PhiAng%d_%s_nRun%d.root",fDataStr[kDataSel],phiN,fPosStr,nrun),"read");

  TF1* ftot[nPt];

  TGraphErrors* gvn_fit = new TGraphErrors();
  TGraphErrors* gvn_fit_err = new TGraphErrors();

  for(int ipt=0; ipt<nPt; ipt++){
    ftot[ipt] = (TF1*) fvnres -> Get(Form("fit_pt%d",ipt));
    gvn_fit->SetPoint(ipt, ptBin_m[ipt],ftot[ipt]->GetParameter(phiN));
    gvn_fit->SetPointError(ipt, exsys[ipt], ftot[ipt]->GetParError(phiN));
    if(isLine){
      gvn_fit_err->SetPoint(ipt, ptBin_m[ipt],ftot[ipt]->GetParameter(phiN));
      gvn_fit_err->SetPointError(ipt, 0, ftot[ipt]->GetParError(phiN));
      gvn_fit->SetPointError(ipt, 0, 0);
    }
  }
  

  //Cos NPi
  TH1D* hvnre_cosglb = (TH1D*) fsv->Get(Form("hist_allEvt_pt_cos%dPi_EPGlauber",phiN));
  TGraphErrors* gvn_fit_glb = new TGraphErrors();
  TGraphErrors* gvn_fit_glb_err = new TGraphErrors();
  for(int ipt=0; ipt<nPt; ipt++){
    gvn_fit_glb->SetPoint(ipt, ptBin_m[ipt], hvnre_cosglb->GetBinContent(ipt+1));
    gvn_fit_glb->SetPointError(ipt, exsys[ipt], hvnre_cosglb->GetBinError(ipt+1));
    if(isLine){
      gvn_fit_glb_err->SetPoint(ipt, ptBin_m[ipt], hvnre_cosglb->GetBinContent(ipt+1));
      gvn_fit_glb_err->SetPointError(ipt, 0, hvnre_cosglb->GetBinError(ipt+1));
      gvn_fit_glb->SetPointError(ipt, 0,0);
    }
  }
  
  TH1D* hvnre_cosenp = (TH1D*) fsv->Get(Form("hist_allEvt_pt_cos%dPi_EPEnProf",phiN));
  TGraphErrors* gvn_fit_enp = new TGraphErrors();
  TGraphErrors* gvn_fit_enp_err = new TGraphErrors();
  for(int ipt=0; ipt<nPt; ipt++){
    gvn_fit_enp->SetPoint(ipt, ptBin_m[ipt], hvnre_cosenp->GetBinContent(ipt+1));
    gvn_fit_enp->SetPointError(ipt, exsys[ipt], hvnre_cosenp->GetBinError(ipt+1));
    if(isLine){
      gvn_fit_enp_err->SetPoint(ipt, ptBin_m[ipt], hvnre_cosenp->GetBinContent(ipt+1));
      gvn_fit_enp_err->SetPointError(ipt, 0, hvnre_cosenp->GetBinError(ipt+1));
      gvn_fit_enp->SetPointError(ipt, 0,0);
    }
  }

  SetGraphStyle(gvn_fit,1,0);
  SetGraphStyle(gvn_fit_glb,0,1);
  SetGraphStyle(gvn_fit_enp,2,2);
  
  if(isLine){ 
    SetGraphStyleSys(gvn_fit_err,1);
    SetGraphStyleSys(gvn_fit_glb_err,0);
    SetGraphStyleSys(gvn_fit_enp_err,2);
    gvn_fit->SetLineWidth(4);
    gvn_fit->SetMarkerSize(0);
    gvn_fit_err->SetFillColorAlpha(kBlue-3,0.4);
    gvn_fit_err->SetMarkerSize(0);
    gvn_fit_glb->SetLineWidth(4);
    gvn_fit_glb->SetMarkerSize(0);
    gvn_fit_glb_err->SetFillColorAlpha(kRed-4,0.4);
    gvn_fit_glb_err->SetMarkerSize(0);
    gvn_fit_enp->SetLineWidth(4);
    gvn_fit_enp->SetMarkerSize(0);
    gvn_fit_enp_err->SetFillColorAlpha(kGreen-4,0.4);
    gvn_fit_enp_err->SetMarkerSize(0);
  }

  SetGraphAxis(gvn_fit,"p_{T}^{#varUpsilon} (GeV/c)",Form("#it{v_{%d}}",phiN));

  //Draw
  TCanvas *c1= new TCanvas("c1","c1",600,600);
  SetCanvasSquare2(c1);

  c1->cd();
  gvn_fit->GetYaxis()->SetTitleOffset(1.4);
  gvn_fit->GetXaxis()->SetTitleOffset(1.1);
  gvn_fit->GetXaxis()->SetTitleSize(0.055);
  gvn_fit->GetXaxis()->SetLabelSize(0.043);
  gvn_fit->GetYaxis()->SetLabelSize(0.043);
  gvn_fit->GetYaxis()->SetNdivisions(510);
  gvn_fit->GetXaxis()->SetLimits(xmin,xmax);
  gvn_fit->GetXaxis()->SetRangeUser(xmin,xmax);
  gvn_fit->SetMinimum(ymin);
  gvn_fit->SetMaximum(ymax);

  if(!isLine) gvn_fit->Draw("AP");
  else if(isLine){
    gvn_fit->Draw("AL");
    gvn_fit_err->Draw("LE3 same");
    gvn_fit_glb_err->Draw("LE3 same");
    gvn_fit_enp_err->Draw("LE3 same");
  }
  gvn_fit_glb->Draw("P");
  gvn_fit_enp->Draw("P");
  dashedLine(xmin,0.,xmax,0.,1,1);

  double legposx1 = 0.275;
  double legposx2 = 0.480;
  double legposy1 = 0.65;
  double legposy2 = 0.82;
  double labtextsize=0.038;

  TLegend* leg= new TLegend(legposx1,legposy1,legposx2,legposy2);
  SetLegendStyle(leg);
  leg->SetTextSize(labtextsize);
  leg->SetHeader("#varUpsilon(1S) (#eta#approx0)");
  if(!isLine){
    leg->AddEntry(gvn_fit,"fit","pe");
    leg->AddEntry(gvn_fit_glb,Form("<cos(%d#phi)> (EP Glauber)",phiN),"pe");
    leg->AddEntry(gvn_fit_enp,Form("<cos(%d#phi)> (EP Energy Density)",phiN),"pe");
  }
  else if(isLine){
    leg->AddEntry(gvn_fit,"fit","l");
    leg->AddEntry(gvn_fit_glb,Form("<cos(%d#phi)> (EP Glauber)",phiN),"l");
    leg->AddEntry(gvn_fit_enp,Form("<cos(%d#phi)> (EP Energy Density)",phiN),"l");
  }
  leg->Draw("same");

  double lab_posx = 0.275; double lab_posy = 0.87; double lab_pos_diff = 0.269;
  drawGlobText("PbPb, #sqrt{s_{NN}} = 5.02 TeV", lab_posx, lab_posy, 1, labtextsize);
  SHINCHONLegend(c1,iPeriod,iPos,drawInner);
  
  string savedir = Form("plots_%s",fPosStr);
  void * dirf = gSystem->OpenDirectory(savedir.c_str());
  if(dirf) gSystem->FreeDirectory(dirf);
  else {gSystem->mkdir(savedir.c_str(), kTRUE);}

  c1->SaveAs(Form("%s/compare_v%d_vs_pt_%s_isLine%d.pdf",savedir.c_str(),phiN,fDataStr[kDataSel],isLine));

	return;
} 
