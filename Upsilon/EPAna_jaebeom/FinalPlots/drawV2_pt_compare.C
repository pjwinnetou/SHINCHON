#include <iostream>
#include "../Style_jaebeom.h"

using namespace std;

void drawV2_pt_compare(int fitor = 3, int nrun = 10, int kInitPos=1, bool isLine = true){

  setTDRStyle();
  writeExtraText= true;
  int iPeriod = 21;
  int iPos = 33;
  bool drawInner = true;

  gStyle->SetEndErrorSize(0);
  gStyle->SetOptFit(kTRUE);
 
  const char* fPosStr;
  if(kInitPos==0) fPosStr = "InitPosZero";
  else if(kInitPos==1) fPosStr = "InitPosGlauber";
  else if(kInitPos==2) fPosStr = "InitPosMean";

  const int nPt = 4;
  double ptBin[nPt+1] = {0, 3, 6, 12, 20};
  double exsys[nPt] =  {1.5, 1.5, 3, 4};
  double ptBin_m[nPt] = {1.5, 4.5, 9, 16};

  double ymin = -0.020; double ymax = 0.021; double xmin = 0; double xmax = 20;

  //// read input file : value & stat.
  // Fit value
  TFile *fV2res = new TFile(Form("../res/FitResults_v2_fitorder%d_%s_nRun%d.root",fitor,fPosStr,nrun),"read");
  TFile *fsv = new TFile(Form("../res/SavedPlots_%s_nRun%d.root",fPosStr,nrun),"read");

  TF1* ftot[nPt];

  TGraphErrors* gv2_fit = new TGraphErrors();
  TGraphErrors* gv2_fit_err = new TGraphErrors();

  for(int ipt=0; ipt<nPt; ipt++){
    ftot[ipt] = (TF1*) fV2res -> Get(Form("fit_pt%d",ipt));
    gv2_fit->SetPoint(ipt, ptBin_m[ipt],ftot[ipt]->GetParameter(2));
    gv2_fit->SetPointError(ipt, exsys[ipt], ftot[ipt]->GetParError(2));
    if(isLine){
      gv2_fit_err->SetPoint(ipt, ptBin_m[ipt],ftot[ipt]->GetParameter(2));
      gv2_fit_err->SetPointError(ipt, 0, ftot[ipt]->GetParError(2));
      gv2_fit->SetPointError(ipt, 0, 0);
    }
  }
  

  //Cos 2Pi
  TH1D* hV2re_cosglb = (TH1D*) fsv->Get("hist_allEvt_pt_cos2Pi_EPGlauber");
  TGraphErrors* gv2_fit_glb = new TGraphErrors();
  TGraphErrors* gv2_fit_glb_err = new TGraphErrors();
  for(int ipt=0; ipt<nPt; ipt++){
    gv2_fit_glb->SetPoint(ipt, ptBin_m[ipt], hV2re_cosglb->GetBinContent(ipt+1));
    gv2_fit_glb->SetPointError(ipt, exsys[ipt], hV2re_cosglb->GetBinError(ipt+1));
    if(isLine){
      gv2_fit_glb_err->SetPoint(ipt, ptBin_m[ipt], hV2re_cosglb->GetBinContent(ipt+1));
      gv2_fit_glb_err->SetPointError(ipt, 0, hV2re_cosglb->GetBinError(ipt+1));
      gv2_fit_glb->SetPointError(ipt, 0,0);
    }

  }
  TH1D* hV2re_cosenp = (TH1D*) fsv->Get("hist_allEvt_pt_cos2Pi_EPEnProf");
  TGraphErrors* gv2_fit_enp = new TGraphErrors();
  TGraphErrors* gv2_fit_enp_err = new TGraphErrors();
  for(int ipt=0; ipt<nPt; ipt++){
    gv2_fit_enp->SetPoint(ipt, ptBin_m[ipt], hV2re_cosenp->GetBinContent(ipt+1));
    gv2_fit_enp->SetPointError(ipt, exsys[ipt], hV2re_cosenp->GetBinError(ipt+1));
    if(isLine){
      gv2_fit_enp_err->SetPoint(ipt, ptBin_m[ipt], hV2re_cosenp->GetBinContent(ipt+1));
      gv2_fit_enp_err->SetPointError(ipt, 0, hV2re_cosenp->GetBinError(ipt+1));
      gv2_fit_enp->SetPointError(ipt, 0,0);
    }
  }

  SetGraphStyle(gv2_fit,1,0);
  SetGraphStyle(gv2_fit_glb,0,1);
  SetGraphStyle(gv2_fit_enp,2,2);
  
  
  if(isLine){ 
    SetGraphStyleSys(gv2_fit_err,1);
    SetGraphStyleSys(gv2_fit_glb_err,0);
    SetGraphStyleSys(gv2_fit_enp_err,2);
    gv2_fit->SetLineWidth(4);
    gv2_fit->SetMarkerSize(0);
    gv2_fit_err->SetFillColorAlpha(kBlue-3,0.4);
    gv2_fit_err->SetMarkerSize(0);
    gv2_fit_glb->SetLineWidth(4);
    gv2_fit_glb->SetMarkerSize(0);
    gv2_fit_glb_err->SetFillColorAlpha(kRed-4,0.4);
    gv2_fit_glb_err->SetMarkerSize(0);
    gv2_fit_enp->SetLineWidth(4);
    gv2_fit_enp->SetMarkerSize(0);
    gv2_fit_enp_err->SetFillColorAlpha(kGreen-4,0.4);
    gv2_fit_enp_err->SetMarkerSize(0);
  }

  SetGraphAxis(gv2_fit,"p_{T}^{#varUpsilon} (GeV/c)","#it{v_{2}}");

  //Draw
  TCanvas *c1= new TCanvas("c1","c1",600,600);
  SetCanvasSquare2(c1);

  c1->cd();
  gv2_fit->GetYaxis()->SetTitleOffset(1.4);
  gv2_fit->GetXaxis()->SetTitleOffset(1.1);
  gv2_fit->GetXaxis()->SetTitleSize(0.055);
  gv2_fit->GetXaxis()->SetLabelSize(0.043);
  gv2_fit->GetYaxis()->SetLabelSize(0.043);
  gv2_fit->GetYaxis()->SetNdivisions(510);
  gv2_fit->GetXaxis()->SetLimits(xmin,xmax);
  gv2_fit->GetXaxis()->SetRangeUser(xmin,xmax);
  gv2_fit->SetMinimum(ymin);
  gv2_fit->SetMaximum(ymax);

  if(!isLine) gv2_fit->Draw("AP");
  else if(isLine){
    gv2_fit->Draw("AL");
    gv2_fit_err->Draw("LE3 same");
    gv2_fit_glb_err->Draw("LE3 same");
    gv2_fit_enp_err->Draw("LE3 same");
  }
  gv2_fit_glb->Draw("P");
  gv2_fit_enp->Draw("P");
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
    leg->AddEntry(gv2_fit,"fit","pe");
    leg->AddEntry(gv2_fit_glb,"<cos(2#phi)> (EP Glauber)","pe");
    leg->AddEntry(gv2_fit_enp,"<cos(2#phi)> (EP Energy Density)","pe");
  }
  else if(isLine){
    leg->AddEntry(gv2_fit,"fit","l");
    leg->AddEntry(gv2_fit_glb,"<cos(2#phi)> (EP Glauber)","l");
    leg->AddEntry(gv2_fit_enp,"<cos(2#phi)> (EP Energy Density)","l");
  }
  leg->Draw("same");

  double lab_posx = 0.275; double lab_posy = 0.87; double lab_pos_diff = 0.269;
  drawGlobText("PbPb, #sqrt{s_{NN}} = 5.02 TeV", lab_posx, lab_posy, 1, labtextsize);
  SHINCHONLegend(c1,iPeriod,iPos,drawInner);
  
  string savedir = Form("plots_%s",fPosStr);
  void * dirf = gSystem->OpenDirectory(savedir.c_str());
  if(dirf) gSystem->FreeDirectory(dirf);
  else {gSystem->mkdir(savedir.c_str(), kTRUE);}

  c1->SaveAs(Form("%s/compare_v2_vs_pt_cos2Pi_isLine%d.pdf",savedir.c_str(),isLine));

	return;
} 
