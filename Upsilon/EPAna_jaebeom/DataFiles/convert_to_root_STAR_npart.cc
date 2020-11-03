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
  const int nPoint = 3 ;
  double npart, RAA, RAA_stat, RAA_sys, RAA_ncoll_uncup,RAA_ncoll_uncdo;
 
  fstream openFile("STAR_RAA_npart_Y1S.dat");
  
  TFile *wfile = new TFile("STAR_RAA_npart_Y1S.root","recreate");
  TTree *tree = new TTree("tree","RAA tree");
         tree -> Branch("npart",&npart,"npart/D");
         tree -> Branch("RAA",&RAA,"RAA/D");
         tree -> Branch("RAA_stat",&RAA_stat,"RAA_stat/D");
         tree -> Branch("RAA_sys",&RAA_sys,"RAA_sys/D");
         tree -> Branch("RAA_ncoll_uncup",&RAA_ncoll_uncup,"RAA_ncoll_uncup/D");
         tree -> Branch("RAA_ncoll_uncdo",&RAA_ncoll_uncdo,"RAA_ncoll_uncdo/D");

  for(int i=0; i<nPoint; i++)
  {
    openFile >> npart >> RAA >> RAA_stat >> RAA_sys >> RAA_ncoll_uncup >> RAA_ncoll_uncdo;
    tree->Fill();
  }
  openFile.close();
  
  const int nEntries = tree->GetEntries(); 
  TGraphErrors *g_RAA = new TGraphErrors();
  TGraphErrors *g_RAA_sys = new TGraphErrors();
  TGraphAsymmErrors *g_RAA_ncollUnc = new TGraphAsymmErrors();

  for(int i=0; i<nEntries; i++){
    tree->GetEntry(i);
    g_RAA->SetPoint(i, npart, RAA);
    g_RAA->SetPointError(i, 0, RAA_stat);

    g_RAA_sys->SetPoint(i, npart, RAA);
    g_RAA_sys->SetPointError(i, 0, RAA_sys);
    
    g_RAA_ncollUnc->SetPoint(i, npart, RAA);
    g_RAA_ncollUnc->SetPointError(i, 0, 0, RAA_ncoll_uncdo, RAA_ncoll_uncup);
  }

  g_RAA->SetName("g_RAA");
  g_RAA_sys->SetName("g_RAA_sys");
  g_RAA_ncollUnc->SetName("g_RAA_nCollUnc");

  wfile ->  cd();
  tree  ->  Write();
  g_RAA ->  Write();
  g_RAA_sys ->  Write();
  g_RAA_ncollUnc ->  Write();
  wfile ->  Close();

  return 0;
}
   
