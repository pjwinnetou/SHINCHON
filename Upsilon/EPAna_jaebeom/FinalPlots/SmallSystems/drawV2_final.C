#include <iostream>
#include "../../Style_jaebeom.h"

using namespace std;

void drawV2_final(){
    setTDRStyle();
    writeExtraText= false;
    int iPeriod = -1;
    int iPos = 33;
	int Ystate = 1;
    bool drawInner = true;

    TFile *fV2respPb = new TFile(Form("v2input/v2_vs_pt_pPb_isLine1_%is_FDall.root",Ystate),"read");
    TFile *fV2respO = new TFile(Form("v2input/v2_vs_pt_pO_isLine1_%is_FDall.root",Ystate),"read");
    TFile *fV2resOO = new TFile(Form("v2input/v2_vs_pt_OO_isLine1_%is_FDall.root",Ystate),"read");

    double ymin = -0.005; double ymax = 0.015; double xmin = 0; double xmax = 20;

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Exclusion graphs");
    TMultiGraph *mgerror = new TMultiGraph();
    mgerror->SetTitle("Exclusion graphs");

    TGraphErrors* gv2_pPb = (TGraphErrors*) fV2respPb->Get("gv2");
    TGraphErrors* gv2_pPb_error = (TGraphErrors*) fV2respPb->Get("gv2_error");
    TGraphErrors* gv2_pO = (TGraphErrors*) fV2respO->Get("gv2");
    TGraphErrors* gv2_pO_error = (TGraphErrors*) fV2respO->Get("gv2_error");
    TGraphErrors* gv2_OO = (TGraphErrors*) fV2resOO->Get("gv2");
    TGraphErrors* gv2_OO_error = (TGraphErrors*) fV2resOO->Get("gv2_error");

    TCanvas *c1= new TCanvas("c1","c1",600,600);
  SetCanvasSquare2(c1);

    SetGraphStyleSys(gv2_pPb,1); SetGraphStyleSys(gv2_pPb_error,1); gv2_pPb_error->SetFillColorAlpha(kBlue-3,0.4);
    SetGraphStyleSys(gv2_pO,4);  SetGraphStyleSys(gv2_pO_error,4); gv2_pO_error->SetFillColorAlpha(kGray+2,0.4);
    SetGraphStyleSys(gv2_OO,0); SetGraphStyleSys(gv2_OO_error,0); gv2_OO_error->SetFillColorAlpha(kRed-4,0.4);

    mg->Add(gv2_pPb); mgerror->Add(gv2_pPb_error);
    mg->Add(gv2_pO); mgerror->Add(gv2_pO_error);
    mg->Add(gv2_OO); mgerror->Add(gv2_OO_error);


  c1->cd();
  mgerror->GetYaxis()->SetTitleOffset(1.4);
  mgerror->GetXaxis()->SetTitleOffset(1.1);
  mgerror->GetXaxis()->SetTitleSize(0.055);
  mgerror->GetXaxis()->SetLabelSize(0.043);
  mgerror->GetYaxis()->SetLabelSize(0.043);
  mgerror->GetYaxis()->SetNdivisions(510);
  mgerror->GetXaxis()->SetLimits(xmin,xmax);
  mgerror->GetXaxis()->SetRangeUser(xmin,xmax);
  mgerror->SetMinimum(ymin);
  mgerror->SetMaximum(ymax);
  mgerror->GetXaxis()->SetTitle("p_{T}^{#varUpsilon} (GeV/c)");
  mgerror->GetYaxis()->SetTitle("#it{v_{2}}");
  mg->Draw("AL"); mgerror->Draw("LE3 SAME"); 

  //SetGraphAxis(mg,"p_{T}^{#varUpsilon} (GeV/c)","#it{v_{2}}");
  //gv2_pPb->Draw("SAME");

   

  dashedLine(xmin,0.,xmax,0.,1,1);

  double legposx1 = 0.275;
  double legposx2 = 0.55;
  double legposy1 = 0.51;
  double legposy2 = 0.71;
  double labtextsize=0.04;
  TLegend* leg= new TLegend(legposx1,legposy1,legposx2,legposy2);
  SetLegendStyle(leg);
  leg->SetTextSize(labtextsize);
  leg->AddEntry(gv2_pO,"p+O","f");
  leg->AddEntry(gv2_OO,"O+O","f");
  leg->AddEntry(gv2_pPb,"p+Pb","f");
  leg->Draw("same");

  double lab_posx = 0.275; double lab_posy = 0.73; double lab_pos_diff = 0.269;
  drawGlobText(Form("Y(%iS), #sqrt{s_{NN}} = 8.16 TeV",Ystate), lab_posx, lab_posy, 1, labtextsize);
  SHINCHONLegend(c1,iPeriod,iPos,drawInner);

    
return;

}
