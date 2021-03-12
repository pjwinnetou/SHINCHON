#include <iostream>
#include "../Style_jaebeom.h"

using namespace std;

void draw_vn_pt(int fitor = 3, int nrun = 10, int kInitPos=1, int kDataSel=0, bool isLine = false, int phiN=3){

  setTDRStyle();
  writeExtraText= false;
  int iPeriod = -1;
  int iPos = 33;
  bool drawInner = true;

  if(fitor<phiN){cout << "fitorder < vn order " << endl; return;}

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
  TFile *fvnres = new TFile(Form("../res/FitResults_v%d_fitorder%d_%s_%s_nRun%d.root",phiN,fitor,fDataStr[kDataSel],fPosStr,nrun),"read");

  TH1D* hvnres[nPt];
  TF1* ftot[nPt];
  TF1* f_v[fitor][nPt];

  TGraphErrors* gvn_fit = new TGraphErrors();
  TGraphErrors* gvn_fit_err = new TGraphErrors();

  for(int ipt=0; ipt<nPt; ipt++){
    hvnres[ipt] = (TH1D*) fvnres->Get(Form("hPhi_allEvt_EPGlauber_pt%d",ipt));
    SetHistStyle(hvnres[ipt],0,0);
    SetHistAxis(hvnres[ipt],"#phi (rad)","#frac{dN}{d#phi}");
    ftot[ipt] = (TF1*) fvnres -> Get(Form("fit_pt%d",ipt));
    gvn_fit->SetPoint(ipt, ptBin_m[ipt],ftot[ipt]->GetParameter(phiN));
    gvn_fit->SetPointError(ipt, exsys[ipt], ftot[ipt]->GetParError(phiN));
    if(isLine){
    gvn_fit_err->SetPoint(ipt, ptBin_m[ipt],ftot[ipt]->GetParameter(phiN));
    gvn_fit_err->SetPointError(ipt, exsys[ipt], ftot[ipt]->GetParError(phiN));
    gvn_fit->SetPointError(ipt, 0,0);
    }
    for(int ifit=0;ifit<fitor;ifit++){
      f_v[ifit][ipt] = (TF1*) fvnres -> Get(Form("fit_v%dcompt_pt%d",ifit,ipt));
    }
    cout << gvn_fit->GetPointX(ipt) << ", " << gvn_fit->GetPointY(ipt) << endl;

  }
  
  TH1D* hVnre_int = (TH1D*) fvnres -> Get("hPhi_allEvt_EPGlauber");
  TF1* ftot_int = (TF1*) fvnres->Get("fit_ptInt");
  TF1* f_v1_int = (TF1*) fvnres->Get("fit_v1compt_ptInt");
  TF1* f_v2_int = (TF1*) fvnres->Get("fit_v2compt_ptInt");
  TF1* f_v3_int; TF1* f_v4_int; 
  if(fitor>2){
    f_v3_int = (TF1*) fvnres->Get("fit_v3compt_ptInt");
    if(phiN==3 && fitor==4) f_v4_int = (TF1*) fvnres->Get("fit_v4compt_ptInt"); 
  }

  SetGraphStyle(gvn_fit,1,0);
  if(isLine){ 
    SetGraphStyleSys(gvn_fit_err,1);
    gvn_fit->SetLineWidth(3);
    gvn_fit->SetMarkerSize(0);
    gvn_fit_err->SetMarkerSize(0);
    gvn_fit_err->SetFillColorAlpha(kBlue-3,0.4);
  }
  SetGraphAxis(gvn_fit,"p_{T}^{#varUpsilon} (GeV/c)",Form("#it{v_{%d}}",phiN));

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
  if(!isLine) leg->AddEntry(gvn_fit,"#varUpsilon(1S) (#eta#approx0)","pe");
  else if(isLine) leg->AddEntry(gvn_fit,"#varUpsilon(1S) (#eta#approx0)","l");
  leg->Draw("same");

  double lab_posx = 0.275; double lab_posy = 0.83; double lab_pos_diff = 0.269;
  drawGlobText("PbPb, #sqrt{s_{NN}} = 5.02 TeV", lab_posx, lab_posy, 1, labtextsize);
  SHINCHONLegend(c1,iPeriod,iPos,drawInner);
  
  string savedir = Form("plots_%s",fPosStr);
  void * dirf = gSystem->OpenDirectory(savedir.c_str());
  if(dirf) gSystem->FreeDirectory(dirf);
  else {gSystem->mkdir(savedir.c_str(), kTRUE);}

  c1->SaveAs(Form("%s/v%d_vs_pt_%s_isLine%d.pdf",savedir.c_str(),phiN,fDataStr[kDataSel],isLine));

  //Draw Each Fit 
  legposy1 = 0.32; legposy2 = 0.57;
  lab_posy = 0.60;
  labtextsize = 0.035;
  TCanvas* cf[nPt];
  TLegend* legf[nPt];
  for(int ipt=0; ipt<nPt; ipt++){
    cf[ipt] = new TCanvas(Form("cf_pt%d",ipt),"",600,600);
    SetCanvasSquare3(cf[ipt]);
    hvnres[ipt]->SetMinimum(hvnres[ipt]->GetMinimum()*0.8);
    hvnres[ipt]->SetMaximum(hvnres[ipt]->GetMaximum()*1.1);
    hvnres[ipt]->Draw("pe");
    ftot[ipt]->Draw("same");
    for(int ifit=0;ifit<fitor;ifit++){
      f_v[ifit][ipt]->Draw("same");
    }

    legf[ipt]= new TLegend(legposx1,legposy1,legposx2,legposy2);
    SetLegendStyle(legf[ipt]);
    legf[ipt]->SetTextSize(labtextsize);
    legf[ipt]->AddEntry(hvnres[ipt],"#varUpsilon(1S) (#eta#approx0)","pe");
    legf[ipt]->AddEntry(ftot[ipt],"Total fit","l");
    for(int ifit=0;ifit<fitor;ifit++){
      f_v[ifit][ipt]->Draw("same");
      legf[ipt]->AddEntry(f_v[ifit][ipt],Form("#it{v_{%d}} component",ifit+1),"l");
    }
    legf[ipt]->Draw("same");
    drawGlobText("PbPb, #sqrt{s_{NN}} = 5.02 TeV", lab_posx, lab_posy, 1, labtextsize);
    SHINCHONLegend(cf[ipt],iPeriod,iPos,drawInner);
    cf[ipt]->SaveAs(Form("%s/v%d_fit_%s_pt%d.pdf",savedir.c_str(),phiN,fDataStr[kDataSel],ipt));
  }





	return;
} 
