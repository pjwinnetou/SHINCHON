#include <iostream>
#include "TFile.h"
#include "TROOT.h"
#include "TF1.h"
#include "Style_jaebeom.h"

using namespace std;

Double_t CosOrd2(Double_t *x, Double_t *p)
{ 
  Double_t xx = x[0];
  Double_t A  = p[0];
  Double_t v1 = p[1];
  Double_t v2 = p[2];
  Double_t out = A*(1+2*(v1*TMath::Cos(xx)+v2*TMath::Cos(2*xx)));
  return out;
} 

Double_t CosOrd3(Double_t *x, Double_t *p)
{ 
  Double_t xx = x[0];
  Double_t A  = p[0];
  Double_t v1 = p[1];
  Double_t v2 = p[2];
  Double_t v3 = p[3];
  Double_t out = A*(1+2*(v1*TMath::Cos(xx)+v2*TMath::Cos(2*xx)+v3*TMath::Cos(3*xx)));
  return out;
} 

Double_t CompOrd1(Double_t *x, Double_t *p)
{ 
  Double_t xx = x[0];
  Double_t A  = p[0];
  Double_t v1 = p[1];
  Double_t out = A*(1+2*(v1*TMath::Cos(xx)));
  return out;
} 

Double_t CompOrd2(Double_t *x, Double_t *p)
{ 
  Double_t xx = x[0];
  Double_t A  = p[0];
  Double_t v1 = p[1];
  Double_t out = A*(1+2*(v1*TMath::Cos(2*xx)));
  return out;
} 

Double_t CompOrd3(Double_t *x, Double_t *p)
{ 
  Double_t xx = x[0];
  Double_t A  = p[0];
  Double_t v1 = p[1];
  Double_t out = A*(1+2*(v1*TMath::Cos(3*xx)));
  return out;
} 


void doFitV2(int nrun = 1000, int InitPos=1, const int fitorder = 2, int kDataSel = 1, int phiN=2)
{
  const char* initPosStr;
  if(InitPos==0) initPosStr = "InitPosZero";
  else if(InitPos==1) initPosStr = "InitPosGlauber";
  else if(InitPos==2) initPosStr = "InitPosMean";

  vector <TString> nameOfFDContribution = {"2S", "2S3S", "2S3S1P", "2S3S1P2P", "all"};

  TFile* fdFraction = new TFile("FeedDownRes/Results_FD_Bottomonium.root","read");
  vector <TString> nameOfFDFraction = {"frac2sTo1s", "frac1pTo1s", "frac3sTo1s", "frac2pTo1s", "frac3pTo1s",
                                      "frac3sTo2s", "frac2pTo2s", "frac3pTo2s",
                                      "frac3pTo3s",
                                      "f1stot", "f2stot", "f3stot"};
    
  const int nfd = nameOfFDFraction.size();
  TF1* ffracFD[nfd];
  TH1D* hfracFD[nfd];
  const int nPt = 4;
  double ptBin[nPt+1]  = {0, 3, 6, 10, 20};
  
  for(int i=0; i<nfd; i++){
    ffracFD[i] = (TF1*) fdFraction->Get(Form("%s",nameOfFDFraction[i].Data()));
    hfracFD[i] = new TH1D(Form("%s",nameOfFDFraction[i].Data()),"",nPt,ptBin);
    hfracFD[i]->Eval(ffracFD[i]);
  }

  const int nStates = 3;
  const int nData = 3;
  enum fData {CMS_502, CMS_276, STAR_200,};
  //const char* fDataStr[nData] = {"CMS502","CMS276","STAR200"};
  const char* fDataStr[nData] = {"pp","pPb","OO"};

  for(int k=0; k<nameOfFDContribution.size(); k++){
    TH1D* hPhiPt[nStates][nPt];
    
    TF1* f1[nPt];
    TF1* f1_v1comp[nPt];
    TF1* f1_v2comp[nPt];
    TF1* f1_v3comp[nPt];

    TCanvas* c1[nStates][nPt];

    int npixel = 1300;

    double fdScaler1S = 0; 
    double fdScaler2S = 0;
    double fdScaler3S = 0;

    for(int ipt=0; ipt<nPt; ipt++){
      for (int istate = 0; istate<nStates; istate++){
        TFile* rf = new TFile(Form("res/SavedPlots_%s_PhiAng%d_%s_nRun%d_%ds.root",fDataStr[kDataSel],phiN,initPosStr,nrun,istate+1),"read");
      
        hPhiPt[istate][ipt] = (TH1D*) rf->Get(Form("hPhi_allEvt_EPGlauber_pt%d",ipt));
        SetHistStyle(hPhiPt[istate][ipt],0,0);
        hPhiPt[istate][ipt]->SetMinimum(0);
        hPhiPt[istate][ipt]->SetMaximum(hPhiPt[istate][ipt]->GetMaximum()*1.2);
      }

      if(strcmp(nameOfFDContribution[k].Data(),"2S")==0) {
        fdScaler2S = hfracFD[0]->GetBinContent(ipt+1)/100;
        fdScaler1S = fdScaler2S + fdScaler3S;
      }
      if(strcmp(nameOfFDContribution[k].Data(),"2S3S")==0) {
        fdScaler2S = hfracFD[0]->GetBinContent(ipt+1)/100;
        fdScaler3S = hfracFD[2]->GetBinContent(ipt+1)/100;
        fdScaler1S = fdScaler2S + fdScaler3S;
      }
      if(strcmp(nameOfFDContribution[k].Data(),"2S3S1P")==0) {
        fdScaler2S = hfracFD[0]->GetBinContent(ipt+1)/100 + hfracFD[1]->GetBinContent(ipt+1)/100;
        fdScaler3S = hfracFD[2]->GetBinContent(ipt+1)/100;
        fdScaler1S = fdScaler2S + fdScaler3S;
      }
      if(strcmp(nameOfFDContribution[k].Data(),"2S3S1P2P")==0) {
        fdScaler2S = hfracFD[0]->GetBinContent(ipt+1)/100 + hfracFD[1]->GetBinContent(ipt+1)/100;
        fdScaler3S = hfracFD[2]->GetBinContent(ipt+1)/100 + hfracFD[3]->GetBinContent(ipt+1)/100;
        fdScaler1S = fdScaler2S + fdScaler3S;
      }
      if(strcmp(nameOfFDContribution[k].Data(),"all")==0) {
        fdScaler2S = hfracFD[0]->GetBinContent(ipt+1)/100 + hfracFD[1]->GetBinContent(ipt+1)/100;
        fdScaler3S = hfracFD[2]->GetBinContent(ipt+1)/100 + hfracFD[3]->GetBinContent(ipt+1)/100 + hfracFD[4]->GetBinContent(ipt+1)/100;
        fdScaler1S = fdScaler2S + fdScaler3S; // = hfracFD[9]->GetBinContent(ipt+1)/100
      }
  
      hPhiPt[0][ipt]->Scale(1-fdScaler1S);
      hPhiPt[1][ipt]->Scale(fdScaler2S);
      hPhiPt[2][ipt]->Scale(fdScaler3S);    
    }

    for (int istate = 0; istate<nStates; istate++){
      TFile* rf = new TFile(Form("res/SavedPlots_%s_PhiAng%d_%s_nRun%d_%ds.root",fDataStr[kDataSel],phiN,initPosStr,nrun,istate+1),"read");
      TFile* wf = new TFile(Form("res/FitResults_v2_fitorder%d_%s_%s_nRun%d_%ds_FD%s.root",fitorder,fDataStr[kDataSel],initPosStr,nrun,istate+1,nameOfFDContribution[k].Data()),"recreate");

      for(int ipt=0; ipt<nPt; ipt++){
        c1[istate][ipt] = new TCanvas(Form("c_%d_%ds",ipt,istate),"",600,600);
        
        for(int i = istate+1; i<nStates; i++)  hPhiPt[istate][ipt]->Add(hPhiPt[i][ipt]);

            if(fitorder==2){ 
            f1[ipt] = new TF1(Form("fit_pt%d",ipt),CosOrd2,-3.2,3.2,fitorder+1);
            f1_v1comp[ipt] = new TF1(Form("fit_v1compt_pt%d",ipt),CompOrd1,-3.2,3.2,2);
            f1_v2comp[ipt] = new TF1(Form("fit_v2compt_pt%d",ipt),CompOrd2,-3.2,3.2,2);
          }
          else if(fitorder==3){
            f1[ipt] = new TF1(Form("fit_pt%d",ipt),CosOrd3,-3.2,3.2,fitorder+1);
            f1_v1comp[ipt] = new TF1(Form("fit_v1compt_pt%d",ipt),CompOrd1,-3.2,3.2,2);
            f1_v2comp[ipt] = new TF1(Form("fit_v2compt_pt%d",ipt),CompOrd2,-3.2,3.2,2);
            f1_v3comp[ipt] = new TF1(Form("fit_v3compt_pt%d",ipt),CompOrd3,-3.2,3.2,2);
          }
          f1[ipt]->SetLineColor(kRed-4);
          f1[ipt]->SetLineWidth(3);
          f1[ipt]->SetNpx(npixel);

          hPhiPt[istate][ipt]->Fit(Form("fit_pt%d",ipt),"R");

          f1_v1comp[ipt]->FixParameter(0, f1[ipt]->GetParameter(0));
          f1_v1comp[ipt]->FixParameter(1, f1[ipt]->GetParameter(1));
          f1_v2comp[ipt]->FixParameter(0, f1[ipt]->GetParameter(0));
          f1_v2comp[ipt]->FixParameter(1, f1[ipt]->GetParameter(2));
          if(fitorder==3){
            f1_v3comp[ipt]->FixParameter(0,f1[ipt]->GetParameter(0));
            f1_v3comp[ipt]->FixParameter(1,f1[ipt]->GetParameter(3));
            f1_v3comp[ipt]->SetLineColor(kMagenta+2);
            f1_v3comp[ipt]->SetLineWidth(2);
            f1_v3comp[ipt]->SetLineStyle(10);
            f1_v3comp[ipt]->SetNpx(npixel);
          }

          f1_v1comp[ipt]->SetLineColor(kBlue+1);
          f1_v1comp[ipt]->SetLineWidth(2);
          f1_v1comp[ipt]->SetLineStyle(kDashed);
          f1_v1comp[ipt]->SetNpx(npixel);

          f1_v2comp[ipt]->SetLineColor(kGreen+1);
          f1_v2comp[ipt]->SetLineWidth(2);
          f1_v2comp[ipt]->SetLineStyle(kDotted);
          f1_v2comp[ipt]->SetNpx(npixel);


          c1[istate][ipt]->cd();
          wf->cd();
          hPhiPt[istate][ipt]->Draw("pe");
          f1[ipt]->Draw("same");
          f1_v1comp[ipt]->Draw("same");
          f1_v2comp[ipt]->Draw("same");

          hPhiPt[istate][ipt]->Write();
          f1[ipt]->Write();
          f1_v1comp[ipt]->Write();
          f1_v2comp[ipt]->Write();
          if(fitorder==3){
            f1_v3comp[ipt]->Draw("same");
            f1_v3comp[ipt]->Write();
          }
        }
      
      //integrated pt
      TCanvas* c1_int = new TCanvas("c_int","",600,600);
      TH1D* hPhiPt_int = (TH1D*) rf->Get("hPhi_allEvt_EPGlauber");
      SetHistStyle(hPhiPt_int,0,0);
      hPhiPt_int->SetMinimum(0);
      hPhiPt_int->SetMaximum(hPhiPt_int->GetMaximum()*1.2);

      TF1* f1_int; TF1* f1_v1comp_int; TF1* f1_v2comp_int; TF1* f1_v3comp_int;
      if(fitorder==2){ 
        f1_int = new TF1("fit_ptInt",CosOrd2,-3.2,3.2,fitorder+1);
        f1_v1comp_int = new TF1("fit_v1compt_ptInt",CompOrd1,-3.2,3.2,2);
        f1_v2comp_int = new TF1("fit_v2compt_ptInt",CompOrd2,-3.2,3.2,2);
      }
      else if(fitorder==3){
        f1_int = new TF1("fit_ptInt",CosOrd3,-3.2,3.2,fitorder+1);
        f1_v1comp_int = new TF1("fit_v1compt_ptInt",CompOrd1,-3.2,3.2,2);
        f1_v2comp_int = new TF1("fit_v2compt_ptInt",CompOrd2,-3.2,3.2,2);
        f1_v3comp_int = new TF1("fit_v3compt_ptInt",CompOrd3,-3.2,3.2,2);
      }
      f1_int->SetLineColor(kRed-4);
      f1_int->SetLineWidth(3);
      f1_int->SetNpx(npixel);

      hPhiPt_int->Fit("fit_ptInt","R");

      f1_v1comp_int->FixParameter(0, f1_int->GetParameter(0));
      f1_v1comp_int->FixParameter(1, f1_int->GetParameter(1));
      f1_v2comp_int->FixParameter(0, f1_int->GetParameter(0));
      f1_v2comp_int->FixParameter(1, f1_int->GetParameter(2));
      if(fitorder==3){
        f1_v3comp_int->FixParameter(0,f1_int->GetParameter(0));
        f1_v3comp_int->FixParameter(1,f1_int->GetParameter(3));
        f1_v3comp_int->SetLineColor(kMagenta+2);
        f1_v3comp_int->SetLineWidth(2);
        f1_v3comp_int->SetLineStyle(10);
        f1_v3comp_int->SetNpx(npixel);
      }

      f1_v1comp_int->SetLineColor(kBlue+1);
      f1_v1comp_int->SetLineWidth(2);
      f1_v1comp_int->SetLineStyle(kDashed);
      f1_v1comp_int->SetNpx(npixel);

      f1_v2comp_int->SetLineColor(kGreen+1);
      f1_v2comp_int->SetLineWidth(2);
      f1_v2comp_int->SetLineStyle(kDotted);
      f1_v2comp_int->SetNpx(npixel);


      c1_int->cd();
      wf->cd();
      hPhiPt_int->Draw("pe");
      f1_int->Draw("same");
      f1_v1comp_int->Draw("same");
      f1_v2comp_int->Draw("same");

      hPhiPt_int->Write();
      f1_int->Write();
      f1_v1comp_int->Write();
      f1_v2comp_int->Write();

      if(fitorder==3){
        f1_v3comp_int->Draw("same");
        f1_v3comp_int->Write();
      }
    }
  }

}
