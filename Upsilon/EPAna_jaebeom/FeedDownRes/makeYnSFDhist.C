#include <iostream>
#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "header.h"

using namespace std;

void makeYnSFDhist(int iExp = fExpCMSnS)
{
  std::unique_ptr<TFile> myFile( TFile::Open(dir+fn[iExp]));
  TDirectory *tdir = (TDirectory*) myFile->Get(fDir[iExp].Data());
  
  TGraphErrors* g12 = new TGraphErrors();
  g12->SetName("graph2sTo1s");
  TGraphErrors* g13 = new TGraphErrors();
  g13->SetName("graph3sTo1s");
  TGraphErrors* g23 = new TGraphErrors();
  g23->SetName("graph3sTo2s");
  
  TIter next(tdir->GetListOfKeys()); TObject *obj;
  while ((obj = next())) {
    obj = ((TKey*) obj) -> ReadObj();
    if(iExp == fExpCMSnS){
      if(TString(obj->GetName()).Contains(fHistName[iExp].Data())) {
        auto h21 = (TH1F*) tdir->Get(fHistName[iExp]+"1");
        auto h21_err_stat = (TH1F*) tdir->Get(fHistName[iExp] + "1_e1");
        auto h21_err_sys = (TH1F*) tdir->Get(fHistName[iExp] + "1_e2");
        auto h31 = (TH1F*) tdir->Get(fHistName[iExp]+"2");
        auto h31_err_stat = (TH1F*) tdir->Get(fHistName[iExp] + "2_e1");
        auto h31_err_sys = (TH1F*) tdir->Get(fHistName[iExp] + "2_e2");

        for(int ip=0; ip<h21->GetNbinsX(); ip++){
          double c = h21->GetBinContent(ip+1);
          double px = h21->GetBinCenter(ip+1);
          double ex = h21->GetBinWidth(ip+1);
          double ey = h21_err_stat->GetBinContent(ip+1);
          double eys = h21_err_sys->GetBinContent(ip+1);

          double corr = brToDM.Y1S/brToDM.Y2S * br.Y2Sto1S;
          c = c*corr;
          double eyt = sqrt(ey*ey + eys*eys)*corr;
          g12->SetPoint(ip, px, c);
          g12->SetPointError(ip, ex/2, eyt);
        }
        for(int ip=0; ip<h31->GetNbinsX(); ip++){
          double c = h31->GetBinContent(ip+1);
          double px = h31->GetBinCenter(ip+1);
          double ex = h31->GetBinWidth(ip+1);
          double ey = h31_err_stat->GetBinContent(ip+1);
          double eys = h31_err_sys->GetBinContent(ip+1);

          double corr = brToDM.Y1S/brToDM.Y3S * br.Y3Sto1S;
          c = c*corr;
          double eyt = sqrt(ey*ey + eys*eys)*corr;
          g13->SetPoint(ip, px, c);
          g13->SetPointError(ip, ex/2, eyt);
        }
        for(int ip=0; ip<g13->GetN(); ip++){
          double px = g12->GetPointX(ip);
          double ex = g12->GetErrorX(ip);
          double c1 = g12->GetPointY(ip);
          double e1 = g12->GetErrorY(ip);
          double c2 = g13->GetPointY(ip);
          double e2 = g13->GetErrorY(ip);

          double c = c2/c1; 
          double err = getErrorPropaDivide(c2,e2,c1,e1);
          double corr = br.Y3Sto2S * (br.Y2Sto1S/br.Y3Sto1S);
          c = c*corr;
          err = err*corr;
          g23->SetPoint(ip, px, c);
          g23->SetPointError(ip, ex/2, err);
        }
      }
    }
    else if(iExp == fExpLHCbnS){
      if(TString(obj->GetName()).Contains(fHistName[iExp].Data())) {
        auto h21 = (TH1F*) tdir->Get(fHistName[iExp]+"1");
        auto h21_err = (TH1F*) tdir->Get(fHistName[iExp] + "1_e1");
        auto h31 = (TH1F*) tdir->Get(fHistName[iExp]+"2");
        auto h31_err = (TH1F*) tdir->Get(fHistName[iExp] + "2_e1");
        auto h32 = (TH1F*) tdir->Get(fHistName[iExp]+"3");
        auto h32_err = (TH1F*) tdir->Get(fHistName[iExp] + "3_e1");

        for(int ip=0; ip<h21->GetNbinsX(); ip++){
          double c = h21->GetBinContent(ip+1);
          double px = h21->GetBinCenter(ip+1);
          double ex = h21->GetBinWidth(ip+1);
          double ey = h21_err->GetBinContent(ip+1);

          double corr = brToDM.Y1S/brToDM.Y2S * br.Y2Sto1S;
          c = c*corr;
          ey = ey*corr;
          g12->SetPoint(ip, px, c);
          g12->SetPointError(ip, ex/2, ey);
        }
        for(int ip=0; ip<h31->GetNbinsX(); ip++){
          double c = h31->GetBinContent(ip+1);
          double px = h31->GetBinCenter(ip+1);
          double ex = h31->GetBinWidth(ip+1);
          double ey = h31_err->GetBinContent(ip+1);

          double corr = brToDM.Y1S/brToDM.Y3S * br.Y3Sto1S;
          c = c*corr;
          ey = ey*corr;
          g13->SetPoint(ip, px, c);
          g13->SetPointError(ip, ex/2, ey);
        }
        for(int ip=0; ip<h32->GetNbinsX(); ip++){
          double c = h32->GetBinContent(ip+1);
          double px = h32->GetBinCenter(ip+1);
          double ex = h32->GetBinWidth(ip+1);
          double ey = h32_err->GetBinContent(ip+1);

          double corr = brToDM.Y2S/brToDM.Y3S * br.Y3Sto2S;
          c = c*corr;
          ey = ey*corr;
          g23->SetPoint(ip, px, c);
          g23->SetPointError(ip, ex/2, ey);
        }
      }
    }
  }

  TFile* rf = new TFile(Form("Plot_%s",fn[iExp].Data()),"recreate");
  rf->cd();
  g12->Write();
  g13->Write();
  g23->Write();
}
