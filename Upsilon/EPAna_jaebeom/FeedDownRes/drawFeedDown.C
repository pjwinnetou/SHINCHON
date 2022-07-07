#include <iostream>
#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "header.h"

using namespace std;

void drawFeedDown()
{
  const int ns = 3;
  TFile *f[ns];
  for(int i=0; i<ns; i++){f[i] =new TFile(Form("chibY%dS_LHCb_8TeV_FD.root",i+1),"read");}
  TGraphAsymmErrors *g1P1S =  (TGraphAsymmErrors*) f[0]->Get("Table 4/Graph1D_y1");
  TGraphAsymmErrors *g2P1S =  (TGraphAsymmErrors*) f[0]->Get("Table 4/Graph1D_y2");
  TGraphAsymmErrors *g3P1S =  (TGraphAsymmErrors*) f[0]->Get("Table 4/Graph1D_y3");
  TGraphAsymmErrors *g2P2S =  (TGraphAsymmErrors*) f[1]->Get("Table 5/Graph1D_y2");
  TGraphAsymmErrors *g3P2S =  (TGraphAsymmErrors*) f[1]->Get("Table 5/Graph1D_y3");
  TGraphAsymmErrors *g3P3S =  (TGraphAsymmErrors*) f[2]->Get("Table 6/Graph1D_y3");
  
  TFile* fCMS = new TFile("Plot_YnS_CMS_7TeV_SR_corr_vs_pT_absy1p2.root","read");
  TGraphErrors* gCMS12 = (TGraphErrors*) fCMS->Get("graph2sTo1s"); 
  TGraphErrors* gCMS13 = (TGraphErrors*) fCMS->Get("graph3sTo1s"); 
  TGraphErrors* gCMS23 = (TGraphErrors*) fCMS->Get("graph3sTo2s"); 
  
  TFile* fLHCb = new TFile("Plot_YnS_LHCb_8TeV_SR_corr_vs_pT_y2to4p5.root","read");
  TGraphErrors* gLHCb12 = (TGraphErrors*) fLHCb->Get("graph2sTo1s"); 
  TGraphErrors* gLHCb13 = (TGraphErrors*) fLHCb->Get("graph3sTo1s"); 
  TGraphErrors* gLHCb23 = (TGraphErrors*) fLHCb->Get("graph3sTo2s"); 
  
  //Graph for 2S->1S
  TGraphErrors* g12f = new TGraphErrors();
  g12f->SetName("g12f");
  DefaultGraphPoint(g12f, gCMS12, 0);
  DefaultGraphPoint(g12f, gLHCb12, g12f->GetN());
  
  TCanvas* c2sto1s = new TCanvas("c2sto1s","",700,700);
  c2sto1s->cd();
  SetGraphAxisRange(g12f,ptLow, ptHigh,yrangeLow,yrangeHigh);
  g12f->Draw("AP");
  TF1* f12s = new TF1("frac2sTo1s","pol4",ptLow,ptHigh);
  f12s->SetParameters(6.9,0.1,0.02,-0.05,0.003);
  f12s->SetParLimits(0,2,10);
  g12f->Fit("frac2sTo1s","R");

  //Graph for 3S->1S
  TGraphErrors* g13f = new TGraphErrors();
  g13f->SetName("g13f");
  DefaultGraphPoint(g13f, gCMS13, 0);
  DefaultGraphPoint(g13f, gLHCb13, g13f->GetN());
  
  TCanvas* c3sto1s = new TCanvas("c3sto1s","",700,700);
  c3sto1s->cd();
  SetGraphAxisRange(g13f,ptLow, ptHigh,yrangeLow,yrangeHigh);
  g13f->Draw("AP");
  TF1* f13s = new TF1("frac3sTo1s","pol4",ptLow,ptHigh);
  f13s->SetParameters(0.6,-0.01,0.02,-0.005,0.0003);
  f13s->SetParLimits(0,0,1);
  f13s->SetParLimits(1,-1,1);
  f13s->SetParLimits(2,-1,1);
  f13s->SetParLimits(3,-1,1);
  f13s->SetParLimits(4,-1,1);
  g13f->Fit("frac3sTo1s","R");

  //Graph for 3S->2S
  TGraphErrors* g23f = new TGraphErrors();
  g23f->SetName("g23f");
  DefaultGraphPoint(g23f, gCMS23, 0);
  DefaultGraphPoint(g23f, gLHCb23, g23f->GetN());
  
  TCanvas* c3sto2s = new TCanvas("c3sto2s","",700,700);
  c3sto2s->cd();
  SetGraphAxisRange(g23f,ptLow, ptHigh,yrangeLow,yrangeHigh);
  g23f->Draw("AP");
  TF1* f23s = new TF1("frac3sTo2s","pol4",ptLow,ptHigh);
  f23s->SetParameters(0.04,-0.02,0.03,-0.005,0.0003);
  f23s->SetParLimits(0,0,10);
  f23s->SetParLimits(1,-1,1);
  f23s->SetParLimits(2,-1,1);
  f23s->SetParLimits(3,-1,1);
  f23s->SetParLimits(4,-1,1);
  g23f->Fit("frac3sTo2s","R");

  //Graph for nP->1S
  TGraphErrors* g1p1sf = new TGraphErrors();
  g1p1sf->SetName("g1p1sf");
  DefaultGraphPoint(g1p1sf, g1P1S, 0);
  
  TCanvas* c1pto1s = new TCanvas("c1pto1s","",700,700);
  c1pto1s->cd();
  SetGraphAxisRange(g1p1sf,ptLow, ptHigh,yrangeLow,yrangeHigh);
  g1p1sf->Draw("AP");
  TF1* f1p1s = new TF1("frac1pTo1s","[0]*(erf((x+[1])/[2])+1)",ptLow,ptHigh);
  f1p1s->SetParameters(10,5,0.2);
  f1p1s->SetParLimits(0,0,100);
  f1p1s->SetParLimits(1,-20,50);
  f1p1s->SetParLimits(2,-30,50);
  g1p1sf->Fit("frac1pTo1s","REM");

  TGraphErrors* g2p1sf = new TGraphErrors();
  g2p1sf->SetName("g2p1sf");
  DefaultGraphPoint(g2p1sf, g2P1S, 0);

  TCanvas* c2pto1s = new TCanvas("c2pto1s","",700,700);
  c2pto1s->cd();
  SetGraphAxisRange(g2p1sf,ptLow, ptHigh,yrangeLow,yrangeHigh);
  g2p1sf->Draw("AP");
  TF1* f2p1s = new TF1("frac2pTo1s","[0]*(erf((x+[1])/[2])+1)",ptLow,ptHigh);
  f2p1s->SetParameters(5,-4,0.5);
  f2p1s->SetParLimits(0,0,150);
  f2p1s->SetParLimits(1,-85,50);
  f2p1s->SetParLimits(2,-30,100);
  g2p1sf->Fit("frac2pTo1s","REM");

  TGraphErrors* g3p1sf = new TGraphErrors();
  g3p1sf->SetName("g3p1sf");
  DefaultGraphPoint(g3p1sf, g3P1S, 0);

  TCanvas* c3pto1s = new TCanvas("c3pto1s","",700,700);
  c3pto1s->cd();
  SetGraphAxisRange(g3p1sf,ptLow, ptHigh,yrangeLow,yrangeHigh);
  g3p1sf->Draw("AP");
  TF1* f3p1s = new TF1("frac3pTo1s","pol1",ptLow,ptHigh);
  g3p1sf->Fit("frac3pTo1s","REM");

  //2,3P -> 2S
  TGraphErrors* g2p2sf = new TGraphErrors();
  g2p2sf->SetName("g2p2sf");
  DefaultGraphPoint(g2p2sf, g2P2S, 0);

  TCanvas* c2pto2s = new TCanvas("c2pto2s","",700,700);
  c2pto2s->cd();
  SetGraphAxisRange(g2p2sf,ptLow, ptHigh,yrangeLow,yrangeHigh);
  g2p2sf->Draw("AP");
  TF1* f2p2s = new TF1("frac2pTo2s","[0]*(erf((x+[1])/[2])+1)",ptLow,ptHigh);
  f2p2s->FixParameter(1,f1p1s->GetParameter(1));
  f2p2s->FixParameter(2,f1p1s->GetParameter(2));
  g2p2sf->Fit("frac2pTo2s","REM");

  TGraphErrors* g3p2sf = new TGraphErrors();
  g3p2sf->SetName("g3p2sf");
  DefaultGraphPoint(g3p2sf, g3P2S, 0);
  
  TCanvas* c3pto2s = new TCanvas("c3pto2s","",700,700);
  c3pto2s->cd();
  SetGraphAxisRange(g3p2sf,ptLow, ptHigh,yrangeLow,yrangeHigh);
  g3p2sf->Draw("AP");
  TF1* f3p2s = new TF1("frac3pTo2s","pol0",ptLow,ptHigh);
  g3p2sf->Fit("frac3pTo2s","REM");

  
  // 3P -> 3S
  TGraphErrors* g3p3sf = new TGraphErrors();
  g3p3sf->SetName("g3p3sf");
  DefaultGraphPoint(g3p3sf, g3P3S, 0);
  
  TCanvas* c3pto3s = new TCanvas("c3pto3s","",700,700);
  c3pto3s->cd();
  SetGraphAxisRange(g3p3sf,ptLow, ptHigh,yrangeLow,yrangeHigh);
  g3p3sf->Draw("AP");
  TF1* f3p3s = new TF1("frac3pTo3s","[0]*(erf((x+[1])/[2])+1)",ptLow,ptHigh);
  f3p3s->FixParameter(1,f1p1s->GetParameter(1));
  f3p3s->FixParameter(2,f1p1s->GetParameter(2));
  g3p3sf->Fit("frac3pTo3s","REM");

  //Set Color
  GraphColorSetting(1, gCMS12, gLHCb12, gCMS13, gLHCb13, g1P1S, g2P1S, g3P1S);
  GraphColorSetting(2, gCMS23, gLHCb23, g2P2S, g3P2S);
  GraphColorSetting(3, g3P3S);

  //draw 1S
  SetGraphAxisRange(gCMS12,ptLow, ptHigh,yrangeLow,85);
  SetGraphAxisSetting(gCMS12,"p_{T}^{#varUpsilon(1S)} (GeV/c)","Feed-down fraction (%)");

  TLegend* leg = new TLegend(0.20,0.72,0.40,0.92);
  SetLegendStyle(leg);
  leg->AddEntry(gCMS12,"CMS 7 TeV #varUpsilon(2S) #rightarrow #varUpsilon(1S)","pe");
  leg->AddEntry(gLHCb12,"LHCb 8 TeV #varUpsilon(2S) #rightarrow #varUpsilon(1S)","pe");
  leg->AddEntry(gCMS13,"CMS 7 TeV #varUpsilon(3S) #rightarrow #varUpsilon(1S)","pe");
  leg->AddEntry(gLHCb13,"LHCb 8 TeV #varUpsilon(3S) #rightarrow #varUpsilon(1S)","pe");
  leg->SetTextSize(0.025);

  TLegend* leg2 = new TLegend(0.58,0.77,0.78,0.92);
  SetLegendStyle(leg2);
  leg2->AddEntry(g1P1S,"LHCb 8 TeV #chi_{b}(1P) #rightarrow #varUpsilon(1S)","pe");
  leg2->AddEntry(g2P1S,"LHCb 8 TeV #chi_{b}(2P) #rightarrow #varUpsilon(1S)","pe");
  leg2->AddEntry(g3P1S,"LHCb 8 TeV #chi_{b}(3P) #rightarrow #varUpsilon(1S)","pe");
  leg2->SetTextSize(0.025);

  gStyle->SetEndErrorSize(0);
  TCanvas *c1S = new TCanvas("c1S","",700,700);
  SetCanvasMarginSetting(c1S);
  c1S->cd();
  gCMS12->Draw("AP");
  gLHCb12->Draw("P");
  g1P1S->Draw("P");
  g2P1S->Draw("P");
  g3P1S->Draw("P");
  gLHCb13->Draw("P");
  gCMS13->Draw("P");
  leg->Draw("same");
  leg2->Draw("same");

  int nfpxl=1500;
  f12s->SetNpx(nfpxl);
  f13s->SetNpx(nfpxl);
  f1p1s->SetNpx(nfpxl);
  f2p1s->SetNpx(nfpxl);
  f3p1s->SetNpx(nfpxl);
  f12s->SetLineColor(kBlue-2);
  f12s->Draw("same");
  f13s->SetLineColor(kOrange);
  f13s->Draw("same");
  f1p1s->SetLineColor(kPink-2);
  f1p1s->Draw("same");
  f2p1s->SetLineColor(kMagenta+2);
  f2p1s->Draw("same");
  f3p1s->SetLineColor(kGreen+2);
  f3p1s->Draw("same");
  TF1 *f1stot = new TF1("f1stot",[&](double*x, double *p){ return f12s->Eval(x[0]) + f13s->Eval(x[0]) + f1p1s->Eval(x[0]) + f2p1s->Eval(x[0]) + f3p1s->Eval(x[0]); },ptLow,ptHigh,0);
  f1stot->SetLineColor(kBlack);
  f1stot->SetNpx(nfpxl);
  f1stot->Draw("same");
  
  TLegend* legf = new TLegend(0.48,0.52,0.78,0.59);
  SetLegendStyle(legf);
  legf->AddEntry(f1stot,"Total feed-down fraction to #varUpsilon(1S)","l");
  legf->SetTextSize(0.025);
  legf->Draw("same");
  c1S->Update();
  c1S->Modified();
  c1S->SaveAs("FeedDownTo1S.pdf");
  c1S->SaveAs("macro_FeedDownTo1S.C");

  

  //draw 2S
  SetGraphAxisRange(gCMS23,ptLow, ptHigh,yrangeLow,65);
  SetGraphAxisSetting(gCMS23,"p_{T}^{#varUpsilon(2S)} (GeV/c)","Feed-down fraction (%)");

  TLegend* leg2s = new TLegend(0.20,0.77,0.40,0.92);
  SetLegendStyle(leg2s);
  leg2s->AddEntry(gCMS23,"CMS 7 TeV #varUpsilon(3S) #rightarrow #varUpsilon(2S)","pe");
  leg2s->AddEntry(gLHCb23,"LHCb 8 TeV #varUpsilon(3S) #rightarrow #varUpsilon(2S)","pe");
  leg2s->SetTextSize(0.025);

  TLegend* leg2_2s = new TLegend(0.58,0.77,0.78,0.92);
  SetLegendStyle(leg2_2s);
  leg2_2s->AddEntry(g2P2S,"LHCb 8 TeV #chi_{b}(2P) #rightarrow #varUpsilon(2S)","pe");
  leg2_2s->AddEntry(g3P2S,"LHCb 8 TeV #chi_{b}(3P) #rightarrow #varUpsilon(2S)","pe");
  leg2_2s->SetTextSize(0.025);

  gStyle->SetEndErrorSize(0);
  TCanvas *c2S = new TCanvas("c2S","",700,700);
  SetCanvasMarginSetting(c2S);
  c2S->cd();
  gCMS23->Draw("AP");
  gLHCb23->Draw("P");
  g2P2S->Draw("P");
  g3P2S->Draw("P");
  leg2s->Draw("same");
  leg2_2s->Draw("same");

  f23s->SetNpx(nfpxl);
  f2p2s->SetNpx(nfpxl);
  f3p2s->SetNpx(nfpxl);
  f23s->SetLineColor(kRed-2);
  f23s->Draw("same");
  f2p2s->SetLineColor(kBlue+2);
  f2p2s->Draw("same");
  f3p2s->SetLineColor(kGreen+2);
  f3p2s->Draw("same");
  
  TF1 *f2stot = new TF1("f2stot",[&](double*x, double *p){ return f23s->Eval(x[0]) + f2p2s->Eval(x[0]) + f3p2s->Eval(x[0]); },ptLow,ptHigh,0);
  f2stot->SetLineColor(kBlack);
  f2stot->SetNpx(nfpxl);
  f2stot->Draw("same");
  
  TF1 *f3stot = new TF1("f3stot",[&](double*x, double *p){ return f3p3s->Eval(x[0]); },ptLow,ptHigh,0);
  f3stot->SetLineColor(kBlack);
  f3stot->SetNpx(nfpxl);
  
  TLegend* legf2s = new TLegend(0.48,0.68,0.78,0.75);
  SetLegendStyle(legf2s);
  legf2s->AddEntry(f2stot,"Total feed-down fraction to #varUpsilon(2S)","l");
  legf2s->SetTextSize(0.025);
  legf2s->Draw("same");
  c2S->SaveAs("FeedDownTo2S.pdf");
  c1S->SaveAs("macro_FeedDownTo2S.C");
  

  //draw 3S
  SetGraphAxisRange(g3P3S,ptLow, ptHigh,yrangeLow,65);
  SetGraphAxisSetting(g3P3S,"p_{T}^{#varUpsilon(3S)} (GeV/c)","Feed-down fraction (%)");
  g3P3S->SetTitle(" ");

  TLegend* leg2_3s = new TLegend(0.58,0.77,0.78,0.92);
  SetLegendStyle(leg2_3s);
  leg2_3s->AddEntry(g3P3S,"LHCb 8 TeV #chi_{b}(3P) #rightarrow #varUpsilon(3S)","pe");
  leg2_3s->SetTextSize(0.025);

  gStyle->SetEndErrorSize(0);
  TCanvas *c3S = new TCanvas("c3S","",700,700);
  SetCanvasMarginSetting(c3S);
  c3S->cd();
  g3P3S->Draw("AP");
  leg2_3s->Draw("same");

  f3p3s->SetLineColor(kMagenta+2);
  f3p3s->Draw("same");
 /* 
  TF1 *f2stot = new TF1("f2stot",[&](double*x, double *p){ return f23s->Eval(x[0]) + f2p2s->Eval(x[0]) + f3p2s->Eval(x[0]); },ptLow,ptHigh,0);
  f2stot->SetLineColor(kBlack);
  f2stot->SetNpx(nfpxl);
  f2stot->Draw("same");
  */
  c3S->SaveAs("FeedDownTo3S.pdf");
  c3S->SaveAs("macro_FeedDownTo3S.C");

  TFile* wf = new TFile("Results_FD_Bottomonium.root","recreate");
  wf->cd();
  f12s->Write();
  f13s->Write();
  f23s->Write();
  f1p1s->Write();
  f2p1s->Write();
  f3p1s->Write();
  f2p2s->Write();
  f3p2s->Write();
  f3p3s->Write();
  f1stot->Write();
  f2stot->Write();
  f3stot->Write();

  
}

