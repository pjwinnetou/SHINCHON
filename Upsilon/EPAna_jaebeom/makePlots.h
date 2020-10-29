#ifndef makePlots_h
#define makePlots_h

#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TChain.h"

using namespace std;
  
double ptBins[] = {0, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 9, 10, 11, 12, 13, 14, 15, 18, 20};
int nPtbins = sizeof(ptBins)/sizeof(double)-1;

double NpartBins[] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420};
int nNpartBins = sizeof(NpartBins)/sizeof(double)-1;

const int nPhibins=20;
double phiBins[nPhibins+1];// = {-3.2,-2.56,-1.92,-1.28,-0.64,0,0.64,1.28,1.92,2.56,3.2};

const int nHist = 2;
const int nVar = 2;
const char *histName[nHist] = {"EnProf","Glauber"};
const char *varName[nVar] = {"pt","phi"};


class Sets
{
  public:
    Sets(){};

    float getDeltaR(float eta1,float phi1,float eta2,float phi2)
    {
      float deltaPhi = TMath::Abs(phi1-phi2);
      float deltaEta = eta1-eta2;
      if(deltaPhi > TMath::Pi())
        deltaPhi = TMath::TwoPi() - deltaPhi;
      return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
    }

    virtual ~Sets();
    virtual void BinAndHistSetting(TH1D*& h1, const char* title, const char* var, int irun);
};

Sets::~Sets()
{
}


void Sets::BinAndHistSetting(TH1D*& h1, const char* title, const char* var, int irun)
{
  for(int ib=0;ib<=nPhibins; ib++){phiBins[ib]= -TMath::Pi() + 2*TMath::Pi()/nPhibins*ib;}

  if(strcmp(var, "pt") == 0){
    if(irun==-1) h1 = new TH1D(Form("hPt_%s",title),";p_{T} (GeV/c);",nPtbins,ptBins);
    else h1 = new TH1D(Form("hPt_%s_run%d",title,irun),";p_{T} (GeV/c);",nPtbins,ptBins);
    cout << "Set TH1D :: " << h1->GetName() << " -- address : " << h1 << endl;
  }
  else if(strcmp(var, "phi") == 0){
    if(irun==-1)  h1 = new TH1D(Form("hPhi_%s",title),";#phi (rad);",nPhibins,phiBins);
    else  h1 = new TH1D(Form("hPhi_%s_run%d",title,irun),";#phi (rad);",nPhibins,phiBins);
    cout << "Set TH1D :: " << h1->GetName() << " -- address : " << h1 << endl;
  }
  else if(strcmp(var, "Npart") == 0){
    if(irun==-1)  h1 = new TH1D(Form("hNpart_%s",title),";<N_{part}>;",nNpartBins,NpartBins);
    else  h1 = new TH1D(Form("hNpart_%s_run%d",title,irun),";<N_{part}>;",nNpartBins,NpartBins);
    cout << "Set TH1D :: " << h1->GetName() << " -- address : " << h1 << endl;
  }
  else{cout << "No matching variable!!" << endl;return;}
};

#endif 
