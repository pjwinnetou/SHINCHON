#include <iostream>
#include "TROOT.h"
#include "TF1.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TreeSetting.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "Style_jaebeom.h"

using namespace std;

void makePlots(int nrun= 10, int kInitPos = 1)
{
  const char* fPosStr;
  if(kInitPos==0) fPosStr = "InitPosZero";
  else if(kInitPos==1) fPosStr = "InitPosGlauber";
  else if(kInitPos==2) fPosStr = "InitPosMean";
  TFile* rf = new TFile(Form("outfile_UpsSkim_%s_0000_%04d.root",fPosStr,nrun),"read");
  TTree* tree = (TTree*) rf->Get("tree");
  if(nrun!= tree->GetEntries()){cout << "ERROR!! :::: Number of entries and runs inconsistent!!" << endl;return;}

  string savedir = "res";
  void * dirf = gSystem->OpenDirectory(savedir.c_str());
  if(dirf) gSystem->FreeDirectory(dirf);
  else {gSystem->mkdir(savedir.c_str(), kTRUE);}
  TFile* wf = new TFile(Form("%s/SavedPlots_%s_nRun%d.root",savedir.c_str(),fPosStr,nrun),"recreate");

  SetTree settree_;
  settree_.TreeSetting(tree);

  const int nPt = 4;
  double ptBin[nPt+1]  = {0, 3, 6, 10, 20};
  double varVect[nHist][nVar];

  map<TString, map<TString,TH1D*>> hist;
  map<TString, map<TString,TH1D*>> hist_noRAA;
  map<TString, map<TString,TH1D*>> hist_ptBin;
  map<TString, map<TString,TH1D*>> hist_ptBin_noRAA;
  
  Sets BinHistSets;

  int it=0;
  for(const auto& t : histName){
    int ivar=0;
    for(auto& var : varName){
      BinHistSets.BinAndHistSetting(hist[it][ivar],Form("allEvt_EP%s",t),Form("%s",var),-1);
      BinHistSets.BinAndHistSetting(hist_noRAA[it][ivar],Form("allEvt_noRAA_EP%s",t),Form("%s",var),-1);
      ivar++;
    }
    for(int ipt=0;ipt<nPt; ipt++){
      BinHistSets.BinAndHistSetting(hist_ptBin[it][ipt],Form("allEvt_EP%s_pt%d",t,ipt),"phi",-1);
      BinHistSets.BinAndHistSetting(hist_ptBin_noRAA[it][ipt],Form("allEvt_noRAA_EP%s_pt%d",t,ipt),"phi",-1);
    }
    it++;
  }

  TH1D* hPhi_allEvt_pt_cos2Pi_EnProf = new TH1D("hist_allEvt_pt_cos2Pi_EPEnProf",";Bin Num;<cos2Pi>",nPt+1,0,5);
  TH1D* hPhi_allEvt_pt_cos2Pi_Glauber = new TH1D("hist_allEvt_pt_cos2Pi_EPGlauber",";Bin Num;<cos2Pi>",nPt+1,0,5);

  TH1D* hRAA_npart;
  BinHistSets.BinAndHistSetting(hRAA_npart,"allEvt_RAA","Npart",-1);
  TH1D* hRAA_npart_ = (TH1D*) hRAA_npart->Clone("hNpart_noRAA");

  TH1D* hRAA_int = new TH1D("hRAA_int",";;",1,0,1);

  const static int CountMax = 100000000;

  vector<double> AvgCos2Pi_en(nPt+1,0);
  vector<double> AvgCos2PiErr_en(nPt+1,0);
  vector<double> weight_s_en(nPt+1,0);
  vector<vector<double>> Cos2Pi_raw_en(nPt+1,vector<double> (CountMax,0));
  vector<vector<double>> weight_en(nPt+1,vector<double> (CountMax,0));
  vector<unsigned int> count_en(nPt+1,0);

  vector<double> AvgCos2Pi_gl(nPt+1,0);
  vector<double> AvgCos2PiErr_gl(nPt+1,0);
  vector<double> weight_s_gl(nPt+1,0);
  vector<vector<double>> Cos2Pi_raw_gl(nPt+1,vector<double> (CountMax,0));
  vector<vector<double>> weight_gl(nPt+1,vector<double> (CountMax,0));
  vector<unsigned int> count_gl(nPt+1,0);

  TLorentzVector* UpsEnProfCor = new TLorentzVector;
  TLorentzVector* UpsGlauberCor = new TLorentzVector;
  TLorentzVector* UpsRaw = new TLorentzVector;
  
  long int RAAint_den = 0;
  long int RAAint_num = 0;
  for(int i=0; i<tree->GetEntries(); i++)
  {
    tree->GetEntry(i);
    for(int iups=0; iups<nUps; iups++)
    {
      UpsEnProfCor = (TLorentzVector*) Ups4momEnProfCor->At(iups);
      UpsGlauberCor = (TLorentzVector*) Ups4momGlauberCor->At(iups);
      UpsRaw = (TLorentzVector*) Ups4momRaw->At(iups);
      
      //Cos 2Pi for EnProf
      for(int ipt=0;ipt<nPt;ipt++)
      {
        if(UpsEnProfCor->Pt()>ptBin[ipt] && UpsEnProfCor->Pt()<ptBin[ipt+1]){
          AvgCos2Pi_en[ipt] += TMath::Cos(2*UpsEnProfCor->Phi())*IsUpsSurv_prob[iups];
          Cos2Pi_raw_en[ipt][count_en[ipt]] = TMath::Cos(2*UpsEnProfCor->Phi())*IsUpsSurv_prob[iups];
          weight_en[ipt][count_en[ipt]] = IsUpsSurv_prob[iups];
          weight_s_en[ipt] += IsUpsSurv_prob[iups];
          count_en[ipt]++;
        }
      }
      
      AvgCos2Pi_en[nPt] += TMath::Cos(2*UpsEnProfCor->Phi())*IsUpsSurv_prob[iups];
      Cos2Pi_raw_en[nPt][count_en[nPt]] = TMath::Cos(2*UpsEnProfCor->Phi())*IsUpsSurv_prob[iups];
      weight_en[nPt][count_en[nPt]] = IsUpsSurv_prob[iups];
      weight_s_en[nPt]++;
      count_en[nPt]++;
      
      //Cos 2Pi for Glauber
      for(int ipt=0;ipt<nPt;ipt++)
      {
        if(UpsGlauberCor->Pt()>ptBin[ipt] && UpsGlauberCor->Pt()<ptBin[ipt+1]){
          AvgCos2Pi_gl[ipt] += TMath::Cos(2*UpsGlauberCor->Phi())*IsUpsSurv_prob[iups];
          Cos2Pi_raw_gl[ipt][count_gl[ipt]] = TMath::Cos(2*UpsGlauberCor->Phi())*IsUpsSurv_prob[iups];
          weight_gl[ipt][count_gl[ipt]] = IsUpsSurv_prob[iups];
          weight_s_gl[ipt] += IsUpsSurv_prob[iups];
          count_gl[ipt]++;
        }
      }
      
      AvgCos2Pi_gl[nPt] += TMath::Cos(2*UpsGlauberCor->Phi())*IsUpsSurv_prob[iups];
      Cos2Pi_raw_gl[nPt][count_gl[nPt]] = TMath::Cos(2*UpsGlauberCor->Phi())*IsUpsSurv_prob[iups];
      weight_gl[nPt][count_gl[nPt]] = IsUpsSurv_prob[iups];
      weight_s_gl[nPt]++;
      count_gl[nPt]++;
      
      //Fill histogram
      settree_.VarVectSet(varVect, iups);
      it=0;
      for(const auto& t : histName){
        int ivar=0;
        for(auto& var : varName){
          hist[it][ivar]->Fill(varVect[it][ivar],IsUpsSurv_prob[iups]);
          hist_noRAA[it][ivar]->Fill(varVect[it][ivar]);
          ivar++;
        }
        for(int ipt=0;ipt<nPt; ipt++){
          if(varVect[it][0] > ptBin[ipt] && varVect[it][0] < ptBin[ipt+1]){
            hist_ptBin[it][ipt]->Fill(varVect[it][1],IsUpsSurv_prob[iups]);
            hist_ptBin_noRAA[it][ipt]->Fill(varVect[it][1]);
          }
        }
        it++;
      }
      
      hRAA_npart->Fill(Npart,IsUpsSurv_prob[iups]);
      hRAA_npart_->Fill(Npart);

      //cal RAA int 
      RAAint_den++;
      if(IsUpsSurv_prob[iups]>0) RAAint_num++; 

    }
  }

  //Calculate cos2Pi error
  for(int ipt=0; ipt<=nPt;ipt++){
    AvgCos2Pi_en[ipt] = AvgCos2Pi_en[ipt]/weight_s_en[ipt];
    for(int ic=0; ic<count_en[ipt]; ic++){AvgCos2PiErr_en[ipt] += pow((Cos2Pi_raw_en[ipt][ic] - AvgCos2Pi_en[ipt])*weight_en[ipt][ic],2);}
    AvgCos2PiErr_en[ipt] = sqrt(AvgCos2PiErr_en[ipt])/weight_s_en[ipt];
    hPhi_allEvt_pt_cos2Pi_EnProf->SetBinContent(ipt+1, AvgCos2Pi_en[ipt]);
    hPhi_allEvt_pt_cos2Pi_EnProf->SetBinError(ipt+1,AvgCos2PiErr_en[ipt]);
  }

  for(int ipt=0; ipt<=nPt;ipt++){
    AvgCos2Pi_gl[ipt] = AvgCos2Pi_gl[ipt]/weight_s_gl[ipt];
    for(int ic=0; ic<count_gl[ipt]; ic++){AvgCos2PiErr_gl[ipt] += pow((Cos2Pi_raw_gl[ipt][ic] - AvgCos2Pi_gl[ipt])*weight_gl[ipt][ic],2);}
    AvgCos2PiErr_gl[ipt] = sqrt(AvgCos2PiErr_gl[ipt])/weight_s_gl[ipt];
    hPhi_allEvt_pt_cos2Pi_Glauber->SetBinContent(ipt+1, AvgCos2Pi_gl[ipt]);
    hPhi_allEvt_pt_cos2Pi_Glauber->SetBinError(ipt+1,AvgCos2PiErr_gl[ipt]);
  }


  //Npart RAA
  hRAA_npart->Divide(hRAA_npart_);

  //RAA int
  double RAAint_den_err = sqrt(RAAint_den);
  double RAAint_num_err = sqrt(RAAint_num);
  double RAAInt;
  double RAAIntErr;
 
  DivideValue(RAAint_num, RAAint_num_err, RAAint_den, RAAint_den_err, RAAInt, RAAIntErr);
  hRAA_int->SetBinContent(1, RAAInt);
  hRAA_int->SetBinError(1, RAAIntErr);


  wf->cd();
  
  hPhi_allEvt_pt_cos2Pi_EnProf->Write();
  hPhi_allEvt_pt_cos2Pi_Glauber->Write();
  hRAA_npart->Write();
  hRAA_int->Write();
      
  it=0;
  for(const auto& t : histName){
    int ivar=0;
    for(auto& var : varName){
      hist[it][ivar]->Write();
      hist_noRAA[it][ivar]->Write();
      ivar++;
    }
    for(int ipt=0;ipt<nPt; ipt++){
        hist_ptBin[it][ipt]->Write();
        hist_ptBin_noRAA[it][ipt]->Write();
    }
    it++;
  }

}
