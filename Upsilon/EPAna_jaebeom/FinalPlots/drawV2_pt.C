#include <iostream>
#include "../Style_jaebeom.h"

using namespace std;

void drawV2_pt(int fitor = 3, int nrun = 10, int kInitPos=1, int kDataSel=0, bool isLine = false){

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

  double ymin = -0.005; double ymax = 0.015; double xmin = 0; double xmax = 20;

  //// read input file : value & stat.
  TFile *fV2res = new TFile(Form("../res/FitResults_v2_fitorder%d_%s_%s_nRun%d.root",fitor,fDataStr[kDataSel],fPosStr,nrun),"read");

  TH1D* hV2res[nPt];
  TF1* ftot[nPt];
  TF1* f_v1[nPt];
  TF1* f_v2[nPt];
  TF1* f_v3[nPt];
  if(fitor==2) delete *f_v3;

  TGraphErrors* gv2_fit = new TGraphErrors();
  TGraphErrors* gv2_fit_err = new TGraphErrors();

  for(int ipt=0; ipt<nPt; ipt++){
    hV2res[ipt] = (TH1D*) fV2res->Get(Form("hPhi_allEvt_EPGlauber_pt%d",ipt));
    SetHistStyle(hV2res[ipt],0,0);
    SetHistAxis(hV2res[ipt],"#phi (rad)","#frac{dN}{d#phi}");
    ftot[ipt] = (TF1*) fV2res -> Get(Form("fit_pt%d",ipt));
    f_v1[ipt] = (TF1*) fV2res -> Get(Form("fit_v1compt_pt%d",ipt));
    f_v2[ipt] = (TF1*) fV2res -> Get(Form("fit_v2compt_pt%d",ipt));
    if(fitor==3) f_v3[ipt] = (TF1*) fV2res -> Get(Form("fit_v3compt_pt%d",ipt));
    gv2_fit->SetPoint(ipt, ptBin_m[ipt],ftot[ipt]->GetParameter(2));
    gv2_fit->SetPointError(ipt, exsys[ipt], ftot[ipt]->GetParError(2));
    if(isLine){
    gv2_fit_err->SetPoint(ipt, ptBin_m[ipt],ftot[ipt]->GetParameter(2));
    gv2_fit_err->SetPointError(ipt, exsys[ipt], ftot[ipt]->GetParError(2));
    gv2_fit->SetPointError(ipt, 0,0);
    }
    cout << gv2_fit->GetPointX(ipt) << ", " << gv2_fit->GetPointY(ipt) << endl;

  }
  
  TH1D* hV2re_int = (TH1D*) fV2res -> Get("hPhi_allEvt_EPGlauber");
  TF1* ftot_int = (TF1*) fV2res->Get("fit_ptInt");
  TF1* f_v1_int = (TF1*) fV2res->Get("fit_v1compt_ptInt");
  TF1* f_v2_int = (TF1*) fV2res->Get("fit_v2compt_ptInt");
  TF1* f_v3_int; if(fitor==3) f_v3_int = (TF1*) fV2res->Get("fit_v3compt_ptInt");

  SetGraphStyle(gv2_fit,1,0);
  if(isLine){ 
    SetGraphStyleSys(gv2_fit_err,1);
    gv2_fit->SetLineWidth(3);
    gv2_fit->SetMarkerSize(0);
    gv2_fit_err->SetMarkerSize(0);
    gv2_fit_err->SetFillColorAlpha(kBlue-3,0.4);
  }
  SetGraphAxis(gv2_fit,"p_{T}^{#varUpsilon} (GeV/c)","#it{v_{2}}");

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
  if(!isLine) leg->AddEntry(gv2_fit,"#varUpsilon(1S) (#eta#approx0)","pe");
  else if(isLine) leg->AddEntry(gv2_fit,"#varUpsilon(1S) (#eta#approx0)","l");
  leg->Draw("same");

  double lab_posx = 0.275; double lab_posy = 0.83; double lab_pos_diff = 0.269;
  drawGlobText("PbPb, #sqrt{s_{NN}} = 5.02 TeV", lab_posx, lab_posy, 1, labtextsize);
  SHINCHONLegend(c1,iPeriod,iPos,drawInner);
  
  string savedir = Form("plots_%s",fPosStr);
  void * dirf = gSystem->OpenDirectory(savedir.c_str());
  if(dirf) gSystem->FreeDirectory(dirf);
  else {gSystem->mkdir(savedir.c_str(), kTRUE);}

  c1->SaveAs(Form("%s/v2_vs_pt_%s_isLine%d.pdf",savedir.c_str(),fDataStr[kDataSel],isLine));

  //Draw Each Fit 
  legposy1 = 0.32; legposy2 = 0.57;
  lab_posy = 0.60;
  labtextsize = 0.035;
  TCanvas* cf[nPt];
  TLegend* legf[nPt];
  for(int ipt=0; ipt<nPt; ipt++){
    cf[ipt] = new TCanvas(Form("cf_pt%d",ipt),"",600,600);
    SetCanvasSquare3(cf[ipt]);
    hV2res[ipt]->SetMinimum(hV2res[ipt]->GetMinimum()*0.8);
    hV2res[ipt]->SetMaximum(hV2res[ipt]->GetMaximum()*1.1);
    hV2res[ipt]->Draw("pe");
    ftot[ipt]->Draw("same");
    f_v1[ipt]->Draw("same");
    f_v2[ipt]->Draw("same");
    if(fitor==3) f_v3[ipt]->Draw("same");

    legf[ipt]= new TLegend(legposx1,legposy1,legposx2,legposy2);
    SetLegendStyle(legf[ipt]);
    legf[ipt]->SetTextSize(labtextsize);
    legf[ipt]->AddEntry(hV2res[ipt],"#varUpsilon(1S) (#eta#approx0)","pe");
    legf[ipt]->AddEntry(ftot[ipt],"Total fit","l");
    legf[ipt]->AddEntry(f_v1[ipt],"#it{v_{1}} component","l");
    legf[ipt]->AddEntry(f_v2[ipt],"#it{v_{2}} component","l");
    if(fitor==3) legf[ipt]->AddEntry(f_v3[ipt],"#it{v_{3}} component","l");
    legf[ipt]->Draw("same");
    drawGlobText("PbPb, #sqrt{s_{NN}} = 5.02 TeV", lab_posx, lab_posy, 1, labtextsize);
    SHINCHONLegend(cf[ipt],iPeriod,iPos,drawInner);
    cf[ipt]->SaveAs(Form("%s/v2_fit_%s_pt%d.pdf",savedir.c_str(),fDataStr[kDataSel],ipt));
  }





	return;
} 
