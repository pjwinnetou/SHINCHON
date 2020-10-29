#ifndef TreeSetting_h
#define TreeSetting_h

#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TChain.h"
#include "makePlots.h"

using namespace std;
  
  static const long MaxUpsSize = 200000;
  Int_t nUps = 0;
  TClonesArray *Ups4momRaw = 0;
  TClonesArray *Ups4momEnProfCor = 0;
  TClonesArray *Ups4momGlauberCor = 0;
  double EPangEnProf, EPangGlauber;
  double IsUpsSurv_modf[MaxUpsSize];
  double IsUpsSurv_prob[MaxUpsSize];
  double vx0[MaxUpsSize];
  double vy0[MaxUpsSize];
  int Npart;
  int Ncoll;
  float b;
  float eccgaus[10];
  float eccpoint[10];
  
class SetTree
{
  public:
    SetTree(){};

    float getDeltaR(float eta1,float phi1,float eta2,float phi2)
    {
      float deltaPhi = TMath::Abs(phi1-phi2);
      float deltaEta = eta1-eta2;
      if(deltaPhi > TMath::Pi())
        deltaPhi = TMath::TwoPi() - deltaPhi;
      return TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
    }

    virtual ~SetTree();
    virtual void TreeSetting(TTree* tree);

    void VarVectSet(double varVect[nHist][nVar], int iups);
};

SetTree::~SetTree()
{
}


void SetTree::TreeSetting(TTree* tree)
{
  //set the branches you will use further
  tree->SetBranchAddress("nUps",&nUps);
  tree->SetBranchAddress("Ups4momRaw",&Ups4momRaw);
  tree->SetBranchAddress("Ups4momEnProfCor",&Ups4momEnProfCor);
  tree->SetBranchAddress("Ups4momGlauberCor",&Ups4momGlauberCor);
  tree->SetBranchAddress("EPangEnProf",&EPangEnProf);
  tree->SetBranchAddress("EPangGlauber",&EPangGlauber);
  tree->SetBranchAddress("IsUpsSurv_modf",IsUpsSurv_modf);
  tree->SetBranchAddress("IsUpsSurv_prob",IsUpsSurv_prob);
  tree->SetBranchAddress("vx0",vx0);
  tree->SetBranchAddress("vy0",vy0);
  tree->SetBranchAddress("Npart",&Npart);
  tree->SetBranchAddress("Ncoll",&Ncoll);
  tree->SetBranchAddress("b",&b);
  tree->SetBranchAddress("eccgaus",eccgaus);
  tree->SetBranchAddress("eccpoint",eccpoint);
};

void SetTree::VarVectSet(double varVect[nHist][nVar], int iups)
{
  TLorentzVector* lorentzen = new TLorentzVector;
  TLorentzVector* lorentzgl = new TLorentzVector;

  lorentzen = (TLorentzVector*) Ups4momEnProfCor->At(iups); 
  lorentzgl = (TLorentzVector*) Ups4momGlauberCor->At(iups); 
   
  varVect[0][0] = lorentzen->Pt(); 
  varVect[0][1] = lorentzen->Phi(); 
  varVect[1][0] = lorentzgl->Pt(); 
  varVect[1][1] = lorentzgl->Phi(); 
  lorentzen->Clear();
  lorentzgl->Clear();
}

#endif 
