#include <iostream>
#include "../Style_jaebeom.h"

using namespace std;

void draw_v2v3_pt_compare(int fitorv2 = 3, int fitorv3 = 3, int nrun = 10, int kInitPos=1, int kDataSel=0,bool isLine = true){
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
  TFile *fv2res = new TFile(Form("../res/FitResults_v2_fitorder%d_%s_%s_nRun%d.root",fitorv2,fDataStr[kDataSel],fPosStr,nrun),"read");
  TFile *fv3res = new TFile(Form("../res/FitResults_v3_fitorder%d_%s_%s_nRun%d.root",fitorv3,fDataStr[kDataSel],fPosStr,nrun),"read");
  TFile *fsv2 = new TFile(Form("../res/SavedPlots_%s_PhiAng2_%s_nRun%d.root",fDataStr[kDataSel],fPosStr,nrun),"read");
  TFile *fsv3 = new TFile(Form("../res/SavedPlots_%s_PhiAng3_%s_nRun%d.root",fDataStr[kDataSel],fPosStr,nrun),"read");

  TF1* ftotv2[nPt];
  TF1* ftotv3[nPt];

  TGraphErrors* gv2_fit = new TGraphErrors();
  TGraphErrors* gv2_fit_err = new TGraphErrors();

  TGraphErrors* gv3_fit = new TGraphErrors();
  TGraphErrors* gv3_fit_err = new TGraphErrors();


  for(int ipt=0; ipt<nPt; ipt++){
    ftotv2[ipt] = (TF1*) fv2res -> Get(Form("fit_pt%d",ipt));
    ftotv3[ipt] = (TF1*) fv3res -> Get(Form("fit_pt%d",ipt));
    gv2_fit->SetPoint(ipt, ptBin_m[ipt],ftotv2[ipt]->GetParameter(2));
    gv2_fit->SetPointError(ipt, exsys[ipt], ftotv2[ipt]->GetParError(2));
    gv3_fit->SetPoint(ipt, ptBin_m[ipt],ftotv3[ipt]->GetParameter(3));
    gv3_fit->SetPointError(ipt, exsys[ipt], ftotv3[ipt]->GetParError(3));
    if(isLine){
      gv2_fit_err->SetPoint(ipt, ptBin_m[ipt],ftotv2[ipt]->GetParameter(2));
      gv2_fit_err->SetPointError(ipt, 0, ftotv2[ipt]->GetParError(2));
      gv2_fit->SetPointError(ipt, 0, 0);
      gv3_fit_err->SetPoint(ipt, ptBin_m[ipt],ftotv3[ipt]->GetParameter(3));
      gv3_fit_err->SetPointError(ipt, 0, ftotv3[ipt]->GetParError(3));
      gv3_fit->SetPointError(ipt, 0, 0);
    }
  }
  
  SetGraphStyle(gv2_fit,1,0);
  SetGraphStyle(gv3_fit,0,1);
  
  if(isLine){ 
    SetGraphStyleSys(gv2_fit_err,1);
    SetGraphStyleSys(gv3_fit_err,0);
    gv2_fit->SetLineWidth(4);
    gv2_fit->SetMarkerSize(0);
    gv2_fit_err->SetFillColorAlpha(kBlue-3,0.4);
    gv2_fit_err->SetMarkerSize(0);
    gv3_fit->SetLineWidth(4);
    gv3_fit->SetMarkerSize(0);
    gv3_fit_err->SetFillColorAlpha(kRed-4,0.4);
    gv3_fit_err->SetMarkerSize(0);
  }

  SetGraphAxis(gv2_fit,"p_{T}^{#varUpsilon} (GeV/c)","#it{v_{2}}");
  SetGraphAxis(gv3_fit,"p_{T}^{#varUpsilon} (GeV/c)","#it{v_{3}}");

  //Draw
  TCanvas *c1= new TCanvas("c1","c1",600,600);
  SetCanvasSquare2(c1);

  c1->cd();
  gv2_fit->GetYaxis()->SetTitleOffset(1.4);
  gv2_fit->GetXaxis()->SetTitleOffset(1.1);
  gv2_fit->GetXaxis()->SetTitleSize(0.055);
  gv2_fit->GetXaxis()->SetLabelSize(0.043);
  gv2_fit->GetYaxis()->SetLabelSize(0.043);
  gv2_fit->GetYaxis()->SetTitle("#it{v_{n}}");
  gv2_fit->GetYaxis()->SetNdivisions(510);
  gv2_fit->GetXaxis()->SetLimits(xmin,xmax);
  gv2_fit->GetXaxis()->SetRangeUser(xmin,xmax);
  gv2_fit->SetMinimum(ymin);
  gv2_fit->SetMaximum(ymax);

  if(!isLine) gv2_fit->Draw("AP");
  else if(isLine){
    gv2_fit->Draw("AL");
    gv2_fit_err->Draw("LE3 same");
    gv3_fit_err->Draw("LE3 same");
  }
  gv3_fit->Draw("P");
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
    leg->AddEntry(gv2_fit,"v_{2}","pe");
    leg->AddEntry(gv3_fit,"v_{3}","pe");
  }
  else if(isLine){
    leg->AddEntry(gv2_fit,"v_{2}","l");
    leg->AddEntry(gv3_fit,"v_{3}","l");
  }
  leg->Draw("same");

  double lab_posx = 0.275; double lab_posy = 0.87; double lab_pos_diff = 0.269;
  drawGlobText("PbPb, #sqrt{s_{NN}} = 5.02 TeV", lab_posx, lab_posy, 1, labtextsize);
  SHINCHONLegend(c1,iPeriod,iPos,drawInner);
  
  string savedir = Form("plots_%s",fPosStr);
  void * dirf = gSystem->OpenDirectory(savedir.c_str());
  if(dirf) gSystem->FreeDirectory(dirf);
  else {gSystem->mkdir(savedir.c_str(), kTRUE);}

  c1->SaveAs(Form("%s/compare_v2v3_vs_pt_%s_isLine%d.pdf",savedir.c_str(),fDataStr[kDataSel],isLine));

	return;
} 
