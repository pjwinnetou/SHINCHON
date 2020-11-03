#include <iostream>
#include "../Style_jaebeom.h"

using namespace std;

void drawRAA_pt_compare_data(int nrun = 10, int kInitPos=1, bool isLine = true){

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

  string savedir = Form("plots_%s",fPosStr);
  void * dirf = gSystem->OpenDirectory(savedir.c_str());
  if(dirf) gSystem->FreeDirectory(dirf);
  else {gSystem->mkdir(savedir.c_str(), kTRUE);}

  const int nData = 3;
  enum fData {CMS_502, CMS_276, STAR_200,};
  const char* fDataStr[nData] = {"CMS502","CMS276","STAR200"};

  const int nPt = 4;
  double ptBin[nPt+1] = {0, 3, 6, 12, 20};
  double exsys[nPt] =  {1.5, 1.5, 3, 4};
  double ptBin_m[nPt] = {1.5, 4.5, 9, 16};

  double ymin = 0; double ymax = 1.1; double xmin = 0; double xmax = 30;

  //// read input file : value & stat.
  TFile* rf[nData];
  TH1D* hPtNum[nData];
  TH1D* hPtDen[nData]; 
  TH1D* hRAA_pt[nData]; 

  TGraphErrors* gRAA[nData];
  TGraphErrors* gRAA_err[nData];

  for(int id=0;id<nData;id++){

    rf[id]= new TFile(Form("../res/SavedPlots_%s_%s_nRun%d.root",fDataStr[id],fPosStr,nrun),"read");
    hPtNum[id]= (TH1D*) rf[id]->Get("hPt_allEvt_EPGlauber");
    hPtDen[id]= (TH1D*) rf[id]->Get("hPt_allEvt_noRAA_EPGlauber");
    hRAA_pt[id]= (TH1D*) hPtNum[id]->Clone("hRAA_pt");
    hRAA_pt[id]->Divide(hPtDen[id]);
    gRAA[id]= new TGraphErrors();
    gRAA_err[id]= new TGraphErrors();

    int iptg=0;
    for(int ipt=0; ipt<hRAA_pt[id]->GetNbinsX(); ipt++){
      if(hRAA_pt[id]->GetBinContent(ipt+1) ==0 ) continue;
      gRAA[id]->SetPoint(iptg, hRAA_pt[id]->GetBinCenter(ipt+1),hRAA_pt[id]->GetBinContent(ipt+1));
      gRAA[id]->SetPointError(iptg, hRAA_pt[id]->GetBinWidth(ipt+1)/2, hRAA_pt[id]->GetBinError(ipt+1));
      if(isLine){
        gRAA_err[id]->SetPoint(iptg, hRAA_pt[id]->GetBinCenter(ipt+1),hRAA_pt[id]->GetBinContent(ipt+1));
        gRAA_err[id]->SetPointError(iptg, hRAA_pt[id]->GetBinWidth(ipt+1)/2, hRAA_pt[id]->GetBinError(ipt+1));
        gRAA[id]->SetPointError(iptg, 0, 0); 
      }
      iptg++;
    }

    SetGraphStyle(gRAA[id],1,0);
    if(isLine){ 
      SetGraphStyleSys(gRAA_err[id],1);
      gRAA[id]->SetLineWidth(3);
      gRAA[id]->SetMarkerSize(0);
      gRAA_err[id]->SetMarkerSize(0);
      gRAA_err[id]->SetFillColorAlpha(kBlue-3,0.4);
    }
    SetGraphAxis(gRAA[id],"p_{T}^{#varUpsilon} (GeV/c)","R_{AA}");
    gRAA[id]->GetYaxis()->SetTitleOffset(1.4);
    gRAA[id]->GetXaxis()->SetTitleOffset(1.1);
    gRAA[id]->GetXaxis()->SetTitleSize(0.055);
    gRAA[id]->GetXaxis()->SetLabelSize(0.043);
    gRAA[id]->GetYaxis()->SetLabelSize(0.043);
    gRAA[id]->GetYaxis()->SetNdivisions(510);
    gRAA[id]->GetXaxis()->SetLimits(xmin,xmax);
    gRAA[id]->GetXaxis()->SetRangeUser(xmin,xmax);
    gRAA[id]->SetMinimum(ymin);
    gRAA[id]->SetMaximum(ymax);
  }

  //Read Data
  //CMS 5.02 TeV
  TFile* fCMS_502 = new TFile("../DataFiles/CMS_RAA_pt.root","read");
  TH1D* hCMS_RAA_502 = (TH1D*) fCMS_502->Get("Table 13/Hist1D_y1");
  TH1D* hCMS_RAA_stat_502 = (TH1D*) fCMS_502->Get("Table 13/Hist1D_y1_e1");
  TH1D* hCMS_RAA_sys_502 = (TH1D*) fCMS_502->Get("Table 13/Hist1D_y1_e2");

  TGraphErrors* gRAA_cms_502 = new TGraphErrors();
  TGraphAsymmErrors* gRAA_cms_sys_502 = new TGraphAsymmErrors();
  for(int ipt=0; ipt<hCMS_RAA_502->GetNbinsX(); ipt++){
    gRAA_cms_502->SetPoint(ipt, hCMS_RAA_502->GetBinCenter(ipt+1), hCMS_RAA_502->GetBinContent(ipt+1));
    gRAA_cms_502->SetPointError(ipt, 0, hCMS_RAA_stat_502->GetBinContent(ipt+1));
    gRAA_cms_sys_502->SetPoint(ipt, hCMS_RAA_502->GetBinCenter(ipt+1), hCMS_RAA_502->GetBinContent(ipt+1));
    gRAA_cms_sys_502->SetPointError(ipt, hCMS_RAA_502->GetBinWidth(ipt+1)/2,hCMS_RAA_502->GetBinWidth(ipt+1)/2, hCMS_RAA_sys_502->GetBinContent(ipt+1), hCMS_RAA_sys_502->GetBinContent(ipt+1));
  }
  SetGraphStyleOpen(gRAA_cms_502,4,0,0);
  SetGraphStyleSys2(gRAA_cms_sys_502,0);

  //CMS 2.76 TeV
  TFile* fCMS_276 = new TFile("../DataFiles/CMS_RAA_Y1S_pT_2p76TeV.root","read");
  TH1D* hCMS_RAA_276 = (TH1D*) fCMS_276->Get("Table 11/Hist1D_y1");
  TH1D* hCMS_RAA_stat_276 = (TH1D*) fCMS_276->Get("Table 11/Hist1D_y1_e1");
  TH1D* hCMS_RAA_sys_276 = (TH1D*) fCMS_276->Get("Table 11/Hist1D_y1_e2");

  TGraphErrors* gRAA_cms_276 = new TGraphErrors();
  TGraphAsymmErrors* gRAA_cms_sys_276 = new TGraphAsymmErrors();
  for(int ipt=0; ipt<hCMS_RAA_276->GetNbinsX(); ipt++){
    gRAA_cms_276->SetPoint(ipt, hCMS_RAA_276->GetBinCenter(ipt+1), hCMS_RAA_276->GetBinContent(ipt+1));
    gRAA_cms_276->SetPointError(ipt, 0, hCMS_RAA_stat_276->GetBinContent(ipt+1));
    gRAA_cms_sys_276->SetPoint(ipt, hCMS_RAA_276->GetBinCenter(ipt+1), hCMS_RAA_276->GetBinContent(ipt+1));
    gRAA_cms_sys_276->SetPointError(ipt, hCMS_RAA_276->GetBinWidth(ipt+1)/2,hCMS_RAA_276->GetBinWidth(ipt+1)/2, hCMS_RAA_sys_276->GetBinContent(ipt+1), hCMS_RAA_sys_276->GetBinContent(ipt+1));
  }
  SetGraphStyleOpen(gRAA_cms_276,4,0,0);
  SetGraphStyleSys2(gRAA_cms_sys_276,0);
  
  //200 GeV
  TFile* fSTAR_200 = new TFile("../DataFiles/STAR_RAA_pt_Y1S.root");
  TGraphErrors* grSTAR_200 = (TGraphErrors*) fSTAR_200 ->Get("g_RAA");
  TGraphErrors* grSTARsys_200 = (TGraphErrors*) fSTAR_200 ->Get("g_RAA_sys");
  
  SetGraphStyleOpen(grSTAR_200,4,0,0);
  SetGraphStyleSys2(grSTARsys_200,0);


  //DRAW
  TCanvas *c1[nData];
  double legposx1 = 0.275;
  double legposx2 = 0.480;
  double legposy1 = 0.60;
  double legposy2 = 0.79;
  double labtextsize=0.038;
  double lab_posx = 0.275; double lab_posy = 0.83; double lab_pos_diff = 0.269;

  TLegend* leg= new TLegend(legposx1,legposy1,legposx2,legposy2);

  //Draw CMS 5.02 TeV
  c1[CMS_502] = new TCanvas("c_CMS_502","c_CMS_502",600,600);
  SetCanvasSquare2(c1[CMS_502]);

  c1[CMS_502]->cd();

  if(!isLine) gRAA[CMS_502]->Draw("AP");
  else if(isLine){
    gRAA[CMS_502]->Draw("AL");
    gRAA_err[CMS_502]->Draw("LE3 same");
  }
  gRAA_cms_sys_502->Draw("5");
  gRAA_cms_502->Draw("P");
  dashedLine(xmin,0.,xmax,0.,1,1);

  SetLegendStyle(leg);
  leg->SetTextSize(labtextsize);
  leg->SetHeader("#varUpsilon(1S)");
  if(!isLine) leg->AddEntry(gRAA[CMS_502],"SHINCHON (#eta#approx0)","pe");
  else if(isLine) leg->AddEntry(gRAA[CMS_502],"SHINCHON (#eta#approx0)","l");
  leg->AddEntry(gRAA_cms_502,"CMS (|y|<2.4), PLB790(2019)10","pe");
  leg->Draw("same");

  drawGlobText("PbPb, #sqrt{s_{NN}} = 5.02 TeV", lab_posx, lab_posy, 1, labtextsize);
  SHINCHONLegend(c1[CMS_502],iPeriod,iPos,drawInner);
  
  c1[CMS_502]->SaveAs(Form("%s/compare_RAA_vs_pt_data_%s_isLine%d.pdf",savedir.c_str(),fDataStr[CMS_502],isLine));


  //Draw CMS 2.76 TeV
  c1[CMS_276] = new TCanvas("c_CMS_276","c_CMS_276",600,600);
  SetCanvasSquare2(c1[CMS_276]);

  c1[CMS_276]->cd();

  if(!isLine) gRAA[CMS_276]->Draw("AP");
  else if(isLine){
    gRAA[CMS_276]->Draw("AL");
    gRAA_err[CMS_276]->Draw("LE3 same");
  }
  gRAA_cms_sys_276->Draw("5");
  gRAA_cms_276->Draw("P");
  dashedLine(xmin,0.,xmax,0.,1,1);

  leg->Clear();
  SetLegendStyle(leg);
  leg->SetTextSize(labtextsize);
  leg->SetHeader("#varUpsilon(1S)");
  if(!isLine) leg->AddEntry(gRAA[CMS_276],"SHINCHON (#eta#approx0)","pe");
  else if(isLine) leg->AddEntry(gRAA[CMS_276],"SHINCHON (#eta#approx0)","l");
  leg->AddEntry(gRAA_cms_276,"CMS (|y|<2.4), PLB770(2017)357","pe");
  leg->Draw("same");

  drawGlobText("PbPb, #sqrt{s_{NN}} = 2.76 TeV", lab_posx, lab_posy, 1, labtextsize);
  SHINCHONLegend(c1[CMS_276],iPeriod,iPos,drawInner);
  
  c1[CMS_276]->SaveAs(Form("%s/compare_RAA_vs_pt_data_%s_isLine%d.pdf",savedir.c_str(),fDataStr[CMS_276],isLine));


  //Draw STAR 200 GeV
  c1[STAR_200] = new TCanvas("c_STAR_200","c_STAR_200",600,600);
  SetCanvasSquare2(c1[STAR_200]);

  c1[STAR_200]->cd();

  xmin = 0; xmax = 12;
  gRAA[STAR_200]->GetXaxis()->SetLimits(xmin,xmax);
  gRAA[STAR_200]->GetXaxis()->SetRangeUser(xmin,xmax);
  ymin = 0; ymax = 1.5;
  gRAA[STAR_200]->SetMinimum(ymin);
  gRAA[STAR_200]->SetMaximum(ymax);
  
  if(!isLine) gRAA[STAR_200]->Draw("AP");
  else if(isLine){
    gRAA[STAR_200]->Draw("AL");
    gRAA_err[STAR_200]->Draw("LE3 same");
  }
  grSTARsys_200->Draw("5");
  grSTAR_200->Draw("P");
  dashedLine(xmin,0.,xmax,0.,1,1);

  leg->Clear();
  SetLegendStyle(leg);
  leg->SetTextSize(labtextsize);
  leg->SetHeader("#varUpsilon(1S)");
  if(!isLine) leg->AddEntry(gRAA[STAR_200],"SHINCHON (#eta#approx0)","pe");
  else if(isLine) leg->AddEntry(gRAA[STAR_200],"SHINCHON (#eta#approx0)","l");
  leg->AddEntry(grSTAR_200,"STAR (|y|<0.5), Preliminary","pe");
  leg->Draw("same");

  drawGlobText("AuAu, #sqrt{s_{NN}} = 200 GeV", lab_posx, lab_posy, 1, labtextsize);
  SHINCHONLegend(c1[STAR_200],iPeriod,iPos,drawInner);
  
  c1[STAR_200]->SaveAs(Form("%s/compare_RAA_vs_pt_data_%s_isLine%d.pdf",savedir.c_str(),fDataStr[STAR_200],isLine));

  return;

} 
