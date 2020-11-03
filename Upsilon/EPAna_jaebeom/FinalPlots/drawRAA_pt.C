#include <iostream>
#include "../Style_jaebeom.h"

using namespace std;

void drawRAA_pt(int nrun = 10, int kInitPos=1, int kDataSel=0, bool isLine = true){

  setTDRStyle();
  writeExtraText= false;
  int iPeriod = -1;
  int iPos = 33;
  bool drawInner = true;

  gStyle->SetEndErrorSize(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
 
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

  double ymin = 0; double ymax = 1.0; double xmin = 0; double xmax = 20;

  //// read input file : value & stat.
  TFile* rf = new TFile(Form("../res/SavedPlots_%s_%s_nRun%d.root",fDataStr[kDataSel],fPosStr,nrun),"read");
  TH1D* hPtNum = (TH1D*) rf->Get("hPt_allEvt_EPGlauber");
  TH1D* hPtDen = (TH1D*) rf->Get("hPt_allEvt_noRAA_EPGlauber");
  TH1D* hRAA_pt = (TH1D*) hPtNum->Clone("hRAA_pt");
  hRAA_pt->Divide(hPtDen);

  TGraphErrors* gRAA = new TGraphErrors();
  TGraphErrors* gRAA_err = new TGraphErrors();

  int iptg=0;
  for(int ipt=0; ipt<hRAA_pt->GetNbinsX(); ipt++){
    if(hRAA_pt->GetBinContent(ipt+1) ==0 ) continue;
    gRAA->SetPoint(iptg, hRAA_pt->GetBinCenter(ipt+1),hRAA_pt->GetBinContent(ipt+1));
    gRAA->SetPointError(iptg, hRAA_pt->GetBinWidth(ipt+1)/2, hRAA_pt->GetBinError(ipt+1));
    if(isLine){
      gRAA_err->SetPoint(iptg, hRAA_pt->GetBinCenter(ipt+1),hRAA_pt->GetBinContent(ipt+1));
      gRAA_err->SetPointError(iptg, hRAA_pt->GetBinWidth(ipt+1)/2, hRAA_pt->GetBinError(ipt+1));
      gRAA->SetPointError(iptg, 0, 0); 
    }
    iptg++;
  }
  
  SetGraphStyle(gRAA,1,0);
  if(isLine){ 
    SetGraphStyleSys(gRAA_err,1);
    gRAA->SetLineWidth(3);
    gRAA->SetMarkerSize(0);
    gRAA_err->SetMarkerSize(0);
    gRAA_err->SetFillColorAlpha(kBlue-3,0.4);
  }
  SetGraphAxis(gRAA,"p_{T}^{#varUpsilon} (GeV/c)","R_{AA}");

  TCanvas *c1= new TCanvas("c1","c1",600,600);
  SetCanvasSquare2(c1);

  c1->cd();
  gRAA->GetYaxis()->SetTitleOffset(1.4);
  gRAA->GetXaxis()->SetTitleOffset(1.1);
  gRAA->GetXaxis()->SetTitleSize(0.055);
  gRAA->GetXaxis()->SetLabelSize(0.043);
  gRAA->GetYaxis()->SetLabelSize(0.043);
  gRAA->GetYaxis()->SetNdivisions(510);
  gRAA->GetXaxis()->SetLimits(xmin,xmax);
  gRAA->GetXaxis()->SetRangeUser(xmin,xmax);
  gRAA->SetMinimum(ymin);
  gRAA->SetMaximum(ymax);

  if(!isLine) gRAA->Draw("AP");
  else if(isLine){
    gRAA->Draw("AL");
    gRAA_err->Draw("LE3 same");
  }
  dashedLine(xmin,0.,xmax,0.,1,1);

  double legposx1 = 0.275;
  double legposx2 = 0.480;
  double legposy1 = 0.71;
  double legposy2 = 0.81;
  double labtextsize=0.038;

  TLegend* leg= new TLegend(legposx1,legposy1,legposx2,legposy2);
  SetLegendStyle(leg);
  leg->SetTextSize(labtextsize);
  if(!isLine) leg->AddEntry(gRAA,"#varUpsilon(1S) (#eta#approx0)","pe");
  else if(isLine) leg->AddEntry(gRAA,"#varUpsilon(1S) (#eta#approx0)","l");
  leg->Draw("same");

  double lab_posx = 0.275; double lab_posy = 0.83; double lab_pos_diff = 0.269;
  drawGlobText("PbPb, #sqrt{s_{NN}} = 5.02 TeV", lab_posx, lab_posy, 1, labtextsize);
  SHINCHONLegend(c1,iPeriod,iPos,drawInner);
  
  string savedir = Form("plots_%s",fPosStr);
  void * dirf = gSystem->OpenDirectory(savedir.c_str());
  if(dirf) gSystem->FreeDirectory(dirf);
  else {gSystem->mkdir(savedir.c_str(), kTRUE);}

  c1->SaveAs(Form("%s/RAA_vs_pt_%s_isLine%d.pdf",savedir.c_str(),fDataStr[kDataSel],isLine));

  return;

} 
