#include <iostream>
#include "TROOT.h"
#include "TSystem.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include <fstream>

using namespace std;

int main(int argc, char **argv)
{
  const int nPoint = 3 ;
  double pT, RAA, RAA_stat, RAA_sys;
  double pT_width;
  fstream openFile("STAR_RAA_pt_Y1S.dat");
  
  TFile *wfile = new TFile("STAR_RAA_pt_Y1S.root","recreate");
  TTree *tree = new TTree("tree","RAA tree");
         tree -> Branch("pT",&pT,"pT/D");
         tree -> Branch("pT_width",&pT_width,"pT_width/D");
         tree -> Branch("RAA",&RAA,"RAA/D");
         tree -> Branch("RAA_stat",&RAA_stat,"RAA_stat/D");
         tree -> Branch("RAA_sys",&RAA_sys,"RAA_sys/D");

  for(int i=0; i<nPoint; i++)
  {
    openFile >> pT >> pT_width >> RAA >> RAA_stat >> RAA_sys;
    tree->Fill();
  }
  openFile.close();
  
  const int nEntries = tree->GetEntries(); 
  TGraphErrors *g_RAA = new TGraphErrors();
  TGraphErrors *g_RAA_sys = new TGraphErrors();

  for(int i=0; i<nEntries; i++){
    tree->GetEntry(i);
    g_RAA->SetPoint(i, pT, RAA);
    g_RAA->SetPointError(i, 0, RAA_stat);

    g_RAA_sys->SetPoint(i, pT, RAA);
    g_RAA_sys->SetPointError(i, pT_width, RAA_sys);
    
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
   
