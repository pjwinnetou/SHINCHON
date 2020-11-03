#include <iostream>
#include "../Style_jaebeom_woFrameLegend.h"
#include "../SHINCHONLegend_raaCent.C"
#include "../tdrstyle.C"

using namespace std;

void drawRAA_nPart_compare_data(int nrun = 10, int kInitPos=1, bool isLine=true)
{
  setTDRStyle();
  writeExtraText= false;
  int iPeriod = -1;
  int iPos = 33;
  bool drawInner = true;
  
  double xmax = 420.0;
  double xmin = 0.0;
  double xmin_int = 0.5;
  double xmax_int = 1.5;
  double ymin =0;
  double ymax = 1.3;
  double boxw = 6.5; // for syst. box (vs cent)
  double boxw_int = 0.09;
  double xlonger = 120;
  
  //Legend
  double legposx1 = 0.40;
  double legposx2 = 0.59;
  double legposy1 = 0.57;
  double legposy2 = 0.75;
  double labtextsize=0.0385;
  
  //text draw
  double lab_posx = 0.275; double lab_posy = 0.87; double lab_pos_diff = 0.269;
  double sz_allign = 0.034;
  double sz_init = 0.574; double sz_step = 0.0558;

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

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile *rf[nData];
  TH1D* hRAA_npart[nData];
  TH1D* hRAA_int[nData];
  TGraphErrors* gRAA_npart[nData];
  TGraphErrors* gRAA_npart_err[nData];
  TGraphErrors* gRAA_int[nData];
  for(int idata=0; idata<nData; idata++){
    rf[idata] = new TFile(Form("../res/SavedPlots_%s_%s_nRun%d.root",fDataStr[idata],fPosStr,nrun),"read");
    hRAA_npart[idata] = (TH1D*) rf[idata]->Get("hNpart_allEvt_RAA");
    hRAA_int[idata] = (TH1D*) rf[idata]->Get("hRAA_int");
    gRAA_npart[idata] = new TGraphErrors();
    gRAA_npart_err[idata] = new TGraphErrors();
    int ipgr=0;
    for(int ip=0; ip<hRAA_npart[idata]->GetNbinsX();ip++){
      if(hRAA_npart[idata]->GetBinContent(ip+1)==0) continue;
      gRAA_npart[idata]->SetPoint(ipgr,hRAA_npart[idata]->GetBinCenter(ip+1),hRAA_npart[idata]->GetBinContent(ip+1));
      gRAA_npart[idata]->SetPointError(ipgr,hRAA_npart[idata]->GetBinWidth(ip+1)/2,hRAA_npart[idata]->GetBinError(ip+1));
      if(isLine){
        gRAA_npart_err[idata]->SetPoint(ipgr,hRAA_npart[idata]->GetBinCenter(ip+1),hRAA_npart[idata]->GetBinContent(ip+1));
        gRAA_npart_err[idata]->SetPointError(ipgr,hRAA_npart[idata]->GetBinWidth(ip+1)/2,hRAA_npart[idata]->GetBinError(ip+1));
        gRAA_npart[idata]->SetPointError(ipgr,0,0);
      }
      ipgr++;
    }
    gRAA_int[idata]= new TGraphErrors();
    gRAA_int[idata]->SetPoint(0,1,hRAA_int[idata]->GetBinContent(1));
    gRAA_int[idata]->SetPointError(0,boxw_int,hRAA_int[idata]->GetBinError(1));

    SetGraphStyle(gRAA_npart[idata],1,0);
    SetGraphStyle(gRAA_int[idata],1,0);
    if(isLine){ 
      SetGraphStyleSys(gRAA_npart[idata],1);
      gRAA_npart[idata]->SetLineWidth(3);
      gRAA_npart[idata]->SetMarkerSize(0);
      gRAA_int[idata]->SetLineWidth(3);
      gRAA_int[idata]->SetMarkerSize(0);
      gRAA_npart_err[idata]->SetMarkerSize(0);
      gRAA_npart_err[idata]->SetFillColorAlpha(kBlue-3,0.4);
    }
    SetGraphAxis(gRAA_npart[idata],"N_{part}","R_{AA}");
    gRAA_npart[idata]->GetYaxis()->SetTitleOffset(1.4);
    gRAA_npart[idata]->GetXaxis()->SetTitleOffset(1.1);
    gRAA_npart[idata]->GetXaxis()->SetTitleSize(0.055);
    gRAA_npart[idata]->GetXaxis()->SetLabelSize(0.043);
    gRAA_npart[idata]->GetYaxis()->SetLabelSize(0.043);
    gRAA_npart[idata]->GetYaxis()->SetNdivisions(510);
    gRAA_npart[idata]->GetXaxis()->SetLimits(xmin,xmax);
    gRAA_npart[idata]->GetXaxis()->SetRangeUser(xmin,xmax);
    gRAA_npart[idata]->SetMinimum(ymin);
    gRAA_npart[idata]->SetMaximum(ymax);
  
    gRAA_int[idata]->GetXaxis()->SetLimits(xmin_int,xmax_int);
    gRAA_int[idata]->SetMinimum(ymin);
    gRAA_int[idata]->SetMaximum(ymax);
    gRAA_int[idata]->GetXaxis()->SetNdivisions(101);
    gRAA_int[idata]->GetXaxis()->SetLabelSize(0);
    gRAA_int[idata]->GetYaxis()->SetTickLength(0.03*600/xlonger);
    gRAA_int[idata]->GetYaxis()->SetLabelSize(0);
  }

  //data

  // 5.02 TeV
  TFile* fCMS_502 = new TFile("../DataFiles/CMS_RAA_npart.root");
  TH1D* hCMS_RAA_502 = (TH1D*) fCMS_502->Get("Table 19/Hist1D_y1"); 
  TH1D* hCMS_RAA_est_502 = (TH1D*) fCMS_502->Get("Table 19/Hist1D_y1_e1"); 
  TH1D* hCMS_RAA_essl_502 = (TH1D*) fCMS_502->Get("Table 19/Hist1D_y1_e2minus"); 
  TH1D* hCMS_RAA_essh_502 = (TH1D*) fCMS_502->Get("Table 19/Hist1D_y1_e2plus"); 
  TH1D* hCMS_Npart_502 = (TH1D*) fCMS_502->Get("Table 19/Hist1D_y2");

  TGraphErrors* grCMS_502 = new TGraphErrors();
  TGraphAsymmErrors* grCMSsys_502 = new TGraphAsymmErrors();
  for(int ip=0; ip<hCMS_RAA_502->GetNbinsX(); ip++){
    grCMS_502->SetPoint(ip, hCMS_Npart_502->GetBinContent(hCMS_RAA_502->GetNbinsX()-ip), hCMS_RAA_502->GetBinContent(hCMS_RAA_502->GetNbinsX()-ip));
    grCMS_502->SetPointError(ip, 0 , hCMS_RAA_est_502->GetBinError(hCMS_RAA_502->GetNbinsX()-ip));
    grCMSsys_502->SetPoint(ip, hCMS_Npart_502->GetBinContent(hCMS_RAA_502->GetNbinsX()-ip), hCMS_RAA_502->GetBinContent(hCMS_RAA_502->GetNbinsX()-ip));
    grCMSsys_502->SetPointError(ip, boxw, boxw , hCMS_RAA_essl_502->GetBinError(hCMS_RAA_502->GetNbinsX()-ip), hCMS_RAA_essh_502->GetBinError(hCMS_RAA_502->GetNbinsX()-ip));
  }
  
  TGraphErrors* grCMSint_502 = new TGraphErrors();
  TGraphAsymmErrors* grCMSintsys_502 = new TGraphAsymmErrors();
  double CMSRAAInt_502 = 0.376; double CMSRAAIntErr_502 = 0.013; double CMSRAAIntErrSysUp_502 = 0.034; double CMSRAAIntErrSysDo_502 = 0.035;
  grCMSint_502->SetPoint(0, 1, CMSRAAInt_502);
  grCMSint_502->SetPointError(0, 0, CMSRAAIntErr_502);
  grCMSintsys_502->SetPoint(0, 1, CMSRAAInt_502);
  grCMSintsys_502->SetPointError(0, boxw_int, boxw_int, CMSRAAIntErrSysUp_502, CMSRAAIntErrSysDo_502);

  SetGraphStyleOpen(grCMS_502,4,0,0);
  SetGraphStyleSys2(grCMSsys_502,0);
  SetGraphStyleOpen(grCMSint_502,4,0,0);
  SetGraphStyleSys2(grCMSintsys_502,0);

  //2.76 TeV
  TFile* fCMS_276 = new TFile("../DataFiles/CMS_RAA_npart_Y1S_2p76TeV.root");
  TGraphErrors* grCMS_276 = (TGraphErrors*) fCMS_276->Get("g_RAA");
  TGraphErrors* grCMSsys_276 = (TGraphErrors*) fCMS_276->Get("g_RAA_sys");
  for(int ip=0; ip<grCMS_276->GetN(); ip++){grCMSsys_276->SetPointError(ip,boxw, grCMSsys_276->GetErrorY(ip));}
  
  TGraphErrors* grCMSint_276 = new TGraphErrors();
  TGraphErrors* grCMSintsys_276 = new TGraphErrors();
  double CMSRAAInt_276 = 0.453; double CMSRAAIntErr_276 = 0.014; double CMSRAAIntErrSys_276 = 0.046;
  grCMSint_276->SetPoint(0, 1, CMSRAAInt_276);
  grCMSint_276->SetPointError(0, 0, CMSRAAIntErr_276);
  grCMSintsys_276->SetPoint(0, 1, CMSRAAInt_276);
  grCMSintsys_276->SetPointError(0, boxw_int, CMSRAAIntErrSys_276);
  
  SetGraphStyleOpen(grCMS_276,4,0,0);
  SetGraphStyleSys2(grCMSsys_276,0);
  SetGraphStyleOpen(grCMSint_276,4,0,0);
  SetGraphStyleSys2(grCMSintsys_276,0);

  //200 GeV
  TFile* fSTAR_200 = new TFile("../DataFiles/STAR_RAA_npart_Y1S.root");
  TGraphErrors* grSTAR_200 = (TGraphErrors*) fSTAR_200 ->Get("g_RAA");
  TGraphErrors* grSTARsys_200 = (TGraphErrors*) fSTAR_200 ->Get("g_RAA_sys");
  TGraphAsymmErrors* grSTARsys200_ncollUnc = (TGraphAsymmErrors*) fSTAR_200->Get("g_RAA_nCollUnc");
  for(int ip=0; ip<grSTAR_200->GetN(); ip++){grSTARsys_200->SetPointError(ip,boxw, grSTARsys_200->GetErrorY(ip)); grSTARsys200_ncollUnc->SetPointError(ip,boxw*0.6,boxw*0.6,grSTARsys200_ncollUnc->GetErrorYlow(ip),grSTARsys200_ncollUnc->GetErrorYhigh(ip));}
  
  SetGraphStyleOpen(grSTAR_200,4,0,0);
  SetGraphStyleSys2(grSTARsys_200,0);
  SetGraphStyleSys2(grSTARsys200_ncollUnc,1);


  
  //// draw  
  // several dataset draw
  TCanvas* c1[nData];
  
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.0384);
  globtex->SetTextAlign(22); //center-center
  globtex->SetTextSize(0.038*600./xlonger);

  //==============================================
  //==============================================
  //===============CMS 5.02 TeV===================
  //==============================================
  //==============================================
  c1[CMS_502] = new TCanvas("c_CMS_502","c_CMS_502",600+xlonger,600);
  
  TPad* pad_diff = new TPad("pad_diff", "",0, 0, 600/(600.+xlonger), 1.0); // vs centrality
  pad_diff->SetLeftMargin(0.2);
  pad_diff->SetRightMargin(0);
  pad_diff->SetBottomMargin(0.14);
  pad_diff->SetTopMargin(0.067);
  TPad* pad_int = new TPad("pad_int", "",600/(600.+xlonger), 0, 1.0, 1.0); // centrality-integrated
  pad_int->SetBottomMargin(0.14);
  pad_int->SetTopMargin(0.067);
  pad_int->SetLeftMargin(0);
  pad_int->SetRightMargin(0.032*600/xlonger);
  TLegend *leg= new TLegend(legposx1,legposy1,legposx2,legposy2);

  //// --- 1st pad!!!   
  c1[CMS_502]->cd();
  pad_diff->Draw(); 
  pad_diff->cd(); 
  if(!isLine) gRAA_npart[CMS_502]->Draw("AP");
  else if(isLine){
    gRAA_npart[CMS_502]->Draw("AL");
    gRAA_npart_err[CMS_502]->Draw("LE3 same");
  }
  grCMSsys_502->Draw("5");
  grCMS_502->Draw("P");
  dashedLine(xmin,1.,xmax,1.,1,1);
  
  //// legend
  SetLegendStyle(leg);
  leg->SetTextSize(labtextsize);
  leg->SetHeader("#varUpsilon(1S)");
  if(!isLine) leg->AddEntry(gRAA_npart[CMS_502],"SHINCHON (#eta#approx0)","pe");
  else if(isLine) leg->AddEntry(gRAA_npart[CMS_502],"SHINCHON (#eta#approx0)","l");
  leg->AddEntry(grCMS_502,"CMS (|y|<2.4), PLB790(2019)10","pe");
  leg->Draw("same");

  drawGlobText("PbPb, #sqrt{s_{NN}} = 5.02 TeV", lab_posx, lab_posy, 1, labtextsize);

  pad_diff->Update();
  SHINCHONLegend_raaCent(pad_diff, iPeriod, iPos);
  pad_diff->Update();

  //// --- 2nd pad!!!   
  c1[CMS_502]->cd();
  pad_int->Draw(); 
  pad_int->cd(); 
  
  //// for int
 
  if(!isLine) gRAA_int[CMS_502]->Draw("AP"); 
  else if(isLine){
    gRAA_int[CMS_502]->SetFillColorAlpha(kBlue-3,0.4);
    gRAA_int[CMS_502]->Draw("ALE");
  }
  grCMSintsys_502->Draw("5");
  grCMSint_502->Draw("P");
  dashedLine(0.,1.,xmax_int,1.,1,1);
  
  pad_int->Update();
  //// draw text
  globtex->DrawLatex(0.5*(1-0.032*600/xlonger), sz_init-sz_step-sz_allign, "<Npart>");

	c1[CMS_502]->Update();
  c1[CMS_502]->SaveAs(Form("%s/compare_RAA_vs_npart_data_%s_isLine%d.pdf",savedir.c_str(),fDataStr[CMS_502],isLine));

  //==============================================
  //==============================================
  //===============CMS 2.76 TeV===================
  //==============================================
  //==============================================
  c1[CMS_276] = new TCanvas("c_CMS_276","c_CMS_276",600+xlonger,600);

  ymin = 0; ymax = 1.7;
  gRAA_npart[CMS_276]->SetMinimum(ymin);
  gRAA_npart[CMS_276]->SetMaximum(ymax);
  gRAA_int[CMS_276]->SetMinimum(ymin);
  gRAA_int[CMS_276]->SetMaximum(ymax);
  //// --- 1st pad!!!   
  c1[CMS_276]->cd();
  pad_diff->Draw(); 
  pad_diff->cd(); 
  if(!isLine) gRAA_npart[CMS_276]->Draw("AP");
  else if(isLine){
    gRAA_npart[CMS_276]->Draw("AL");
    gRAA_npart_err[CMS_276]->Draw("LE3 same");
  }
  grCMSsys_276->Draw("5");
  grCMS_276->Draw("P");
  dashedLine(xmin,1.,xmax,1.,1,1);
  
  //// legend
  leg->Clear();
  legposy1=0.62; legposy2=0.80;
  leg->SetY1NDC(legposy1); leg->SetY2NDC(legposy2);
  SetLegendStyle(leg);
  leg->SetTextSize(labtextsize);
  leg->SetHeader("#varUpsilon(1S)");
  if(!isLine) leg->AddEntry(gRAA_npart[CMS_276],"SHINCHON (#eta#approx0)","pe");
  else if(isLine) leg->AddEntry(gRAA_npart[CMS_276],"SHINCHON (#eta#approx0)","l");
  leg->AddEntry(grCMS_276,"CMS (|y|<2.4), PLB770(2017)357","pe");
  leg->Draw("same");

  drawGlobText("PbPb, #sqrt{s_{NN}} = 2.76 TeV", lab_posx, lab_posy, 1, labtextsize);

  pad_diff->Update();
  SHINCHONLegend_raaCent(pad_diff, iPeriod, iPos);
  pad_diff->Update();

  //// --- 2nd pad!!!   
  c1[CMS_276]->cd();
  pad_int->Draw(); 
  pad_int->cd(); 
  
  //// for int
  if(!isLine) gRAA_int[CMS_276]->Draw("AP"); 
  else if(isLine){
    gRAA_int[CMS_276]->SetFillColorAlpha(kBlue-3,0.4);
    gRAA_int[CMS_276]->Draw("ALE");
  }
  grCMSintsys_276->Draw("5");
  grCMSint_276->Draw("P");
  dashedLine(0.,1.,xmax_int,1.,1,1);
  
  pad_int->Update();
  //// draw text
  globtex->DrawLatex(0.5*(1-0.032*600/xlonger), sz_init-sz_step-sz_allign, "<Npart>");

	c1[CMS_276]->Update();
  c1[CMS_276]->SaveAs(Form("%s/compare_RAA_vs_npart_data_%s_isLine%d.pdf",savedir.c_str(),fDataStr[CMS_276],isLine));

  //==============================================
  //==============================================
  //===============STAR 200 GeV===================
  //==============================================
  //==============================================
  c1[STAR_200] = new TCanvas("c_STAR_200","c_STAR_200",600+xlonger,600);
  
  ymin = 0; ymax = 1.2;
  gRAA_npart[STAR_200]->SetMinimum(ymin);
  gRAA_npart[STAR_200]->SetMaximum(ymax);
  gRAA_int[STAR_200]->SetMinimum(ymin);
  gRAA_int[STAR_200]->SetMaximum(ymax);

  //// --- 1st pad!!!   
  c1[STAR_200]->cd();
  pad_diff->Draw(); 
  pad_diff->cd(); 
  if(!isLine) gRAA_npart[STAR_200]->Draw("AP");
  else if(isLine){
    gRAA_npart[STAR_200]->Draw("AL");
    gRAA_npart_err[STAR_200]->Draw("LE3 same");
  }
  grSTARsys200_ncollUnc->Draw("5");
  grSTARsys_200->Draw("5");
  grSTAR_200->Draw("P");
  dashedLine(xmin,1.,xmax,1.,1,1);
  
  //// legend
  leg->Clear();
  legposy1=0.61; legposy2=0.79;
  leg->SetY1NDC(legposy1); leg->SetY2NDC(legposy2);
  SetLegendStyle(leg);
  leg->SetTextSize(labtextsize);
  leg->SetHeader("#varUpsilon(1S)");
  if(!isLine) leg->AddEntry(gRAA_npart[STAR_200],"SHINCHON (#eta#approx0)","pe");
  else if(isLine) leg->AddEntry(gRAA_npart[STAR_200],"SHINCHON (#eta#approx0)","l");
  leg->AddEntry(grSTAR_200,"STAR (|y|<0.5), Preliminary","pe");
  leg->Draw("same");

  drawGlobText("AuAu, #sqrt{s_{NN}} = 200 GeV", lab_posx, lab_posy, 1, labtextsize);

  pad_diff->Update();
  SHINCHONLegend_raaCent(pad_diff, iPeriod, iPos);
  pad_diff->Update();
  
  //// --- 2nd pad!!!   
  c1[STAR_200]->cd();
  pad_int->Draw(); 
  pad_int->cd(); 
  
  //// for int
  if(!isLine) gRAA_int[STAR_200]->Draw("AP"); 
  else if(isLine){
    gRAA_int[STAR_200]->SetFillColorAlpha(kBlue-3,0.4);
    gRAA_int[STAR_200]->Draw("ALE");
  }
  dashedLine(0.,1.,xmax_int,1.,1,1);
  //// draw text
  globtex->DrawLatex(0.5*(1-0.032*600/xlonger), sz_init-sz_step-sz_allign, "<Npart>");
	c1[STAR_200]->Update();
  c1[STAR_200]->SaveAs(Form("%s/compare_RAA_vs_npart_data_%s_isLine%d.pdf",savedir.c_str(),fDataStr[STAR_200],isLine));


	return;

} // end of main func.

