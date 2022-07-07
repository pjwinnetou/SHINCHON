#ifndef header_H
#define header_H

#include <TF1.h>
#include <TCut.h>
#include <TChain.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <iostream>
#include <TLine.h>
#include <TMath.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>

double br_1P_0to1Sgamma = 0.0194;
double br_1P_1to1Sgamma = 0.352;
double br_1P_2to1Sgamma = 0.18;

double br_2P_0to1Sgamma = 3.8e-05;
double br_2P_1to1Sgamma = 0.099;
double br_2P_2to1Sgamma = 0.066;

double br_2P_0to2Sgamma = 0.0138;
double br_2P_1to2Sgamma = 0.181;
double br_2P_2to2Sgamma = 0.089;

double br_2Sto1S = 17.85 + 8.6 + 2.9e-04;// + 6.9*br_1P_1to1Sgamma + 7.15*br_1P_2to1Sgamma + 3.8*br_1P_0to1Sgamma;
double br_3Sto1S = 4.37 + 2.20;// + 10.6*br_2Sto1S + 13.1*(br_2P_2to1Sgamma + br_2P_2to2Sgamma*br_2Sto1S) + 12.6*(br_2P_1to1Sgamma + br_2P_1to2Sgamma*br_2Sto1S) + 5.9*(br_2P_0to1Sgamma + br_2P_0to2Sgamma*br_2Sto1S);
double br_3Sto2S = 10.6;

struct branchingRatioToDimuons { double Y1S, Y2S, Y3S;};
branchingRatioToDimuons brToDM = {2.48, 1.93, 2.18};

struct branchingRatio { double Y2Sto1S, Y3Sto1S, Y3Sto2S;};
branchingRatio br = {br_2Sto1S, br_3Sto1S, br_3Sto2S};

int fExpCMSnS =0;
int fExpLHCbnS =1;
TString fDir[2] = {"Table 6","Table 19"};
TString fHistName[2] = {"Hist1D_y","Hist1D_y"};
TString fn[2] = {"YnS_CMS_7TeV_SR_corr_vs_pT_absy1p2.root","YnS_LHCb_8TeV_SR_corr_vs_pT_y2to4p5.root"};
TString dir = "/home/CMS/DataFiles/QuarkoniaFeedDown/";

double ptLow = 0; double ptHigh=50;
double yrangeLow = 0; double yrangeHigh=50;

double getErrorPropaDivide(double a, double da, double b, double db){
  double err = a/b*TMath::Sqrt((da/a)*(da/a) + (db/b)*(db/b));
  return err;
};

void GraphColorSetting(int state, TGraph *g1=0, TGraph *g2=0, TGraph *g3=0, TGraph *g4=0, TGraph *g5=0, TGraph *g6=0, TGraph *g7=0) 
{
  int col;
  if(state==1){
    col = TColor::GetColor("#0229be");
    g1->SetLineColor(col);
    g1->SetMarkerColor(col);
    g1->SetMarkerStyle(kFullTriangleUp);
    g1->SetMarkerSize(1.2);

    col = TColor::GetColor("#50b7ef");
    g2->SetLineColor(col);
    g2->SetMarkerColor(col);
    g2->SetMarkerStyle(kFullTriangleDown);
    g2->SetMarkerSize(1.2);

    col = TColor::GetColor("#de8900");
    g3->SetLineColor(col);
    g3->SetMarkerColor(col);
    g3->SetMarkerStyle(kFullTriangleUp);
    g3->SetMarkerSize(1.2);

    col = TColor::GetColor("#fad855");
    g4->SetLineColor(col);
    g4->SetMarkerColor(col);
    g4->SetMarkerStyle(kFullTriangleDown);
    g4->SetMarkerSize(1.2);

    col = TColor::GetColor("#e93300");
    g5->SetLineColor(col);
    g5->SetMarkerColor(col);
    g5->SetMarkerStyle(kFullCircle);
    g5->SetMarkerSize(1.2);

    col = TColor::GetColor("#a52ff7");
    g6->SetLineColor(col);
    g6->SetMarkerColor(col);
    g6->SetMarkerStyle(kFullCircle);
    g6->SetMarkerSize(1.2);

    col = TColor::GetColor("#3cd600");
    g7->SetLineColor(col);
    g7->SetMarkerColor(col);
    g7->SetMarkerStyle(kFullCircle);
    g7->SetMarkerSize(1.2);
  }
  else if(state==2){
    col = TColor::GetColor("#e73ca6");
    g1->SetLineColor(col);
    g1->SetMarkerColor(col);
    g1->SetMarkerStyle(kFullTriangleUp);
    g1->SetMarkerSize(1.2);

    col = TColor::GetColor("#ec1818");
    g2->SetLineColor(col);
    g2->SetMarkerColor(col);
    g2->SetMarkerStyle(kFullTriangleDown);
    g2->SetMarkerSize(1.2);

    col = TColor::GetColor("#4199dc");
    g3->SetLineColor(col);
    g3->SetMarkerColor(col);
    g3->SetMarkerStyle(kFullCircle);
    g3->SetMarkerSize(1.2);

    col = TColor::GetColor("#a3d928");
    g4->SetLineColor(col);
    g4->SetMarkerColor(col);
    g4->SetMarkerStyle(kFullCircle);
    g4->SetMarkerSize(1.2);
  }
  else if(state==3){
    col = TColor::GetColor("#9b43ed");
    g1->SetLineColor(col);
    g1->SetMarkerColor(col);
    g1->SetMarkerStyle(kFullCircle);
    g1->SetMarkerSize(1.2);
  }
};

void DefaultGraphPoint(TGraphErrors* g1=0, TGraphErrors* g2=0, int n=0)
{
  for(int ip=0; ip<g2->GetN(); ip++){
    double px = g2->GetPointX(ip);
    double py = g2->GetPointY(ip);
    double ex = g2->GetErrorX(ip);
    double ey = g2->GetErrorY(ip);
    g1->SetPoint(n+ip,px,py);
    g1->SetPointError(n+ip,ex,ey);
  }
};

void DefaultGraphPoint(TGraphErrors* g1=0, TGraphAsymmErrors* g2=0, int n=0)
{
  for(int ip=0; ip<g2->GetN(); ip++){
    double px = g2->GetPointX(ip);
    double py = g2->GetPointY(ip);
    double ex = g2->GetErrorX(ip);
    double ey = g2->GetErrorY(ip);
    g1->SetPoint(n+ip,px,py);
    g1->SetPointError(n+ip,ex,ey);
  }
};

void SetGraphAxisRange(TGraph* g1=0, double ptLow=0, double ptHigh=50, double yrangeLow=0, double yrangeHigh=50)
{
  g1->GetXaxis()->SetLimits(ptLow,ptHigh);
  g1->GetXaxis()->SetRangeUser(ptLow,ptHigh);
  g1->GetYaxis()->SetLimits(yrangeLow,yrangeHigh);
  g1->GetYaxis()->SetRangeUser(yrangeLow,yrangeHigh);
};

void SetGraphAxisSetting(TGraph* g1=0, const char* xtitle = "", const char* ytitle="")
{
  g1->GetXaxis()->SetTitle("p_{T}^{#varUpsilon(1S)} (GeV/c)");
  g1->GetYaxis()->SetTitle("Feed-down fraction (%)");
  g1->GetXaxis()->CenterTitle();
  g1->GetYaxis()->CenterTitle();
  g1->GetYaxis()->SetTitleOffset(1.2);
  g1->GetXaxis()->SetTitleOffset(1.2);
  g1->GetYaxis()->SetTitleSize(0.047);
  g1->GetXaxis()->SetTitleSize(0.047);
  g1->GetYaxis()->SetLabelSize(0.035);
  g1->GetXaxis()->SetLabelSize(0.035);
};

void SetCanvasMarginSetting(TCanvas* c1)
{
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.04);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.04);
  c1->SetTicks(1,1);
};

void SetLegendStyle(TLegend* l) {
  l->SetFillColor(0);
  l->SetFillStyle(4000);
  l->SetBorderSize(0);
  l->SetMargin(0.2);
  l->SetTextSize(0.029);
  l->SetTextFont(42);
};


#endif
