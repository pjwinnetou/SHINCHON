#include <iostream>
#include "TROOT.h"
#include "TSystem.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include <fstream>

using namespace std;

int main(int argc, char **argv)
{
  const int nPoint = 8 ;
  double npart, RAA, RAA_stat, RAA_sys;
 
  fstream openFile("CMS_RAA_npart_Y1S_2p76TeV.dat");
  
  TFile *wfile = new TFile("CMS_RAA_npart_Y1S_2p76TeV.root","recreate");
  TTree *tree = new TTree("tree","RAA tree");
         tree -> Branch("npart",&npart,"npart/D");
         tree -> Branch("RAA",&RAA,"RAA/D");
         tree -> Branch("RAA_stat",&RAA_stat,"RAA_stat/D");
         tree -> Branch("RAA_sys",&RAA_sys,"RAA_sys/D");

  for(int i=0; i<nPoint; i++)
  {
    openFile >> npart >> RAA >> RAA_stat >> RAA_sys;
    tree->Fill();
  }
  openFile.close();
  
  const int nEntries = tree->GetEntries(); 
  TGraphErrors *g_RAA = new TGraphErrors();
  TGraphErrors *g_RAA_sys = new TGraphErrors();

  for(int i=0; i<nEntries; i++){
    tree->GetEntry(i);
    g_RAA->SetPoint(i, npart, RAA);
    g_RAA->SetPointError(i, 0, RAA_stat);

    g_RAA_sys->SetPoint(i, npart, RAA);
    g_RAA_sys->SetPointError(i, 0, RAA_sys);
  }

  g_RAA->SetName("g_RAA");
  g_RAA_sys->SetName("g_RAA_sys");

  wfile ->  cd();
  tree  ->  Write();
  g_RAA ->  Write();
  g_RAA_sys ->  Write();
  wfile ->  Close();

  return 0;
}
   
