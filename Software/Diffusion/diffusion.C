// j.nagle 
//
// 12/21/2012 
// re-edited with a clean version and modified to take input function of temperature for D (not just a/T)
// 01/15/2014
// adding tracking of change in delta-phi versus time

#include "TROOT.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TBox.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TNtuple.h"
#include "TArrow.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TSystem.h"

#include <iostream>
#include <fstream>

using namespace std;

void SetGraphProps(TGraph* g,
		   Int_t linecolor,
		   Int_t markercolor,
		   Int_t markerstyle,
		   Double_t markersize);

void diffusion(
         const string    inFileName = "hydro.root", // input nagle-hydro.root file
         const bool      plotThings = true, // whether to create canvases for final results. Non-essential for large batch jobs
         const bool      yesbeauty = false, // only beauty, otherwise only charm
	       const int       NQuarks = 1000,
	       int             N_timesteps = 5000, // just give it larger than the 508 steps
	       const TString   output_filename = "outfile.root",
         const double    etaOverS = 1,

         const int       eventNum = 1, // parameterA is 3/2pi * etaOverS
	       //double    parameterA = 3.0/(2.*pi), //  / TMath::TwoPi(), // A = DT
	       // for reference, 3/2pi is equivalent to eta/s = 1/4pi
	       // this parameter is overloaded for negative values to use the Dlookup
	       // set parameterA to -1, -2, -3 for temperature dependent cases

	       const double    timeStepScale = 1.0,   // If > 1, skip timesteps, coasting in between. 
	       const bool      onepanel = false,      // if you just want the hydro+quark animation
	       const TString   save = "none",         // Canvases: "all", "none", or "last"
	       const double    ScaleDdown = 1.0,      // ignore this 
	       const double    NonInteractTime = 0.0  // how long before you allow the quarks to interact
	       ) 
{

  gROOT->Reset();
  if (gSystem->Getenv("TMPDIR"))
    gSystem->SetBuildDir(gSystem->Getenv("TMPDIR"));
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //==============================================================================================
  TFile* inFile = new TFile(inFileName.c_str(), "read");
  //==============================================================================================
  // j.nagle - this is how it was from the b6.5 file ??? arggg...
  // 4/18/2018 Jeff Ouellette: static_cast doesn't work right in root 5.. for root 6 this works:
  // .... TH2D *quarkxy = static_cast <TH2D *> (fin2->Get("ncoll_rp"));
  //TH2D* quarkxy= static_cast <TH2D*> (fin2->Get(inNcollDistName.c_str()));
  //
  // for root 5 use this instead:
  TH2D* quarkxy = (TH2D*)inFile->Get("ED_0");
  if (! quarkxy) cout << "Did not find quark x,y distribution TH2D" << endl;

  TH1D* htime = (TH1D*)inFile->Get("Time");
  if (! htime) cout << "Did not find time distribution TH1D" << endl;
  //==============================================================================================


  const double pi = TMath::Pi();
  const double ptMin = 0.0;
  const double ptMax = 30.0;
  const int numPtBins = 300;
  const int numPhiBins = 120;

  const double radii[5] = {0, 0.5, 1.0, 1.7, 2.5};

  const double QuarkMass                   = (yesbeauty)? 4.20 : 1.27;
	const int    NGridX                      = quarkxy->GetNbinsX();   // # of x and y grid cells
	const int    NGridY                      = quarkxy->GetNbinsY(); 
  double time_step_length = htime->GetBinContent(2) - htime->GetBinContent(1);
  double time_step_length_inverseGeV = time_step_length * (1.0/0.1975);
  const double temperature_cutoff          = 0.1;
  const double parameterA                  = etaOverS * 1.0/(4*pi);

  const bool sampleRapidity                = false; // samples a near-central rapidity
  const bool MomentumDependent             = true; // scales D with total momentum

  cout << "DIFFUSION CODE RUNNING with parameter A = " << etaOverS << " x 1/4pi" << endl;

  time_step_length            *= timeStepScale;
  time_step_length_inverseGeV *= timeStepScale;

  TF1 *fkt = new TF1("fkt","TMath::Exp(-x*x/(2.*[0]*[0]))",-10.0,10.0);
  fkt->SetParameter(0,1.0);

  //==============================================================================================
  TH1D *Dlookup; // temperature dependence of D parameter (default as in Teaney and Moore is a/T)
  if (parameterA > 0) {
    Dlookup = new TH1D("Dlookup","Dlookup",100,0.0,0.500); // up to 500 MeV/c
    for (int i=1;i<=100;i++) 
      Dlookup->SetBinContent(i, parameterA / Dlookup->GetBinCenter(i));  // just D=a/T
  } else {
    // read in Dlookup from an external file !!!!
    TFile *inFileDlookup = new TFile("DlookupModel.root");
    if (parameterA == -1) Dlookup  = (TH1D*) inFileDlookup->Get("DlookupA");
    if (parameterA == -2) Dlookup  = (TH1D*) inFileDlookup->Get("DlookupB");
    if (parameterA == -3) Dlookup  = (TH1D*) inFileDlookup->Get("DlookupC");
    if (parameterA == -4) Dlookup  = (TH1D*) inFileDlookup->Get("DlookupD");
    if (parameterA == -5) Dlookup  = (TH1D*) inFileDlookup->Get("DlookupPerfect");
  }
  //==============================================================================================

  //==============================================================================================
  // input for charm and beauty initial pT shapes
  // Alternative to using TH1s from heavy_quark_pt.root
  TF1* ptc = new TF1("ptc", "[0]*x*TMath::Power(x*x + pow([1],2), [2])", 0., 20.);
  ptc->SetParameters(3.3, 2.1, -3.9);    // parameters from arXiv:1205.2396v1 (Cao,Qin, Bass)
  ptc->SetLineColor(kRed);

  TF1* ptb = new TF1("ptb", "[0]*x*TMath::Power(x*x + pow([1],2), [2])", 0., 20.);
  ptb->SetParameters(900, 7.5, -4.9);
  ptb->SetLineColor(kBlue);
  //==============================================================================================

  //==============================================================================================
  // Hydro histos for each timestep - to be read from inFile
  const int numhists = N_timesteps+1;
  TH2D *htemperature[numhists];
  TH2D* hbetax[numhists];
  TH2D* hbetay[numhists];

  TH1D* hheavypt[numhists][5];
  TH1D* hheavyraa[numhists][5];
  TH2D* hheavyptphi[numhists][5];
  TH2D* hheavyraa_particle[5];
  TH2D* hheavydpt_particle[5];

  TH1D* hdeltaphiorig = new TH1D("hdeltaphiorig","hdeltaphiorig",30,0.0,pi);
  TH1D* hdeltaphiorig1 = new TH1D("hdeltaphiorig1","hdeltaphiorig1",30,0.0,pi);
  TH1D* hdeltaphiorig2 = new TH1D("hdeltaphiorig2","hdeltaphiorig2",30,0.0,pi);
  TH1D* hdeltaphiorig3 = new TH1D("hdeltaphiorig3","hdeltaphiorig3",30,0.0,pi);
  TH1D* hdeltaphi = new TH1D("hdeltaphi","hdeltaphi",30,0.0,pi);
  hdeltaphi->SetXTitle("charm-anticharm #Delta #phi (rad)");
  TH1D* hdeltaphi1 = new TH1D("hdeltaphi1","hdeltaphi1",30,0.0,pi);
  TH1D* hdeltaphi2 = new TH1D("hdeltaphi2","hdeltaphi2",30,0.0,pi);
  TH1D* hdeltaphi3 = new TH1D("hdeltaphi3","hdeltaphi3",30,0.0,pi);

  TLorentzVector* tquark = new TLorentzVector();
  TLorentzVector* tquark1 = new TLorentzVector();
  TLorentzVector* tquark2 = new TLorentzVector();
  TLatex ltx, ltxn;  ltxn.SetNDC();
  TObjArray* boosts[numhists]; // TArrow* arrays
  TRandom3 ran(0);                  // 0 means new seed every run
  
	TCanvas* cc = new TCanvas("cc", "cc", 10,10,1200,400);
	if (onepanel) cc->SetWindowSize(500,500);  
	if (!onepanel) cc->Divide(3,1);

  // Read in hydro histograms: T, E density, beta, and pt.
  for (int itime=0; itime<N_timesteps; itime++) {
    
		cout << Form("Reading hydro inputs for timestep %d/%d", itime, N_timesteps) << endl;
    
    // check if this information exists in hydro file or not
    if (! inFile->Get(Form("ED_%d",    itime))) {
			N_timesteps = itime - 1;
      cout << "Last hydro timestep = " << N_timesteps << endl;
      break;
    }

    htemperature[itime]   = (TH2D*)inFile->Get(Form("T_%d",    itime));
    htemperature[itime]->SetXTitle("x coordinate [fm]");
    htemperature[itime]->SetYTitle("y coordinate [fm]");
    hbetax[itime]         = (TH2D*)inFile->Get(Form("Vx_%d",   itime));
    hbetay[itime]         = (TH2D*)inFile->Get(Form("Vy_%d",   itime));
    
    // For graphics: boost vectors for every 10th x,y cell
    boosts[itime] = new TObjArray();
    for(int ybin=1; ybin<=NGridY; ybin++) {
      for(int xbin=1; xbin<=NGridX; xbin++) {  	  
	      if(xbin%10==0 && ybin%10==0) {
	        double x = htemperature[itime]->GetXaxis()->GetBinCenter(xbin);
	        double y = htemperature[itime]->GetYaxis()->GetBinCenter(ybin);
	        double beta_x = hbetax[itime]->GetBinContent(xbin,ybin);
	        double beta_y = hbetay[itime]->GetBinContent(xbin,ybin);
	        double length = TMath::Sqrt(beta_x*beta_x + beta_y*beta_y);
	        TArrow *boost_vector = 
	          new TArrow(x, y, x+beta_x, y+beta_y, 0.005*length, ">");
	        boosts[itime]->Add(boost_vector);
	      }
      }
    } // end loop over x,y, grid
  } // end loop over timesteps

	for (int ir = 0; ir < 5; ir++) {
		for (int itime = 0; itime <= N_timesteps; itime++) { // inclusive with N_timesteps to get endpoints of interval
			if (plotThings || itime==0 || itime == N_timesteps) {
				hheavypt[itime][ir] = new TH1D (Form ("hheavypt_itime%i_radius%i_%i", itime, ir, ir+1), "", numPtBins, ptMin, ptMax);
				hheavyraa[itime][ir] = new TH1D (Form ("hheavyraa_itime%i_radius%i_%i", itime, ir, ir+1), "", numPtBins, ptMin, ptMax);
				hheavyptphi[itime][ir] = new TH2D (Form ("hheavyptphi_itime%i_radius%i_%i", itime, ir, ir+1), "", numPtBins, ptMin, ptMax, numPhiBins, -pi, pi);
			}
		}
		hheavyraa_particle[ir] = new TH2D (Form ("hheavyraa_particle_radius%i_%i", ir, ir+1), "", numPtBins, ptMin, ptMax, 200, 0.0, 10.0);
		hheavydpt_particle[ir] = new TH2D (Form ("hheavydpt_particle_radius%i_%i", ir, ir+1), "", numPtBins, ptMin, ptMax, 200, 0.0, 10.0);
	}

	//return;

  //--------------------------------------------------------------------------------------
  // NOW STEP QUARKS THROUGH THE TIME EVOLUTION
  //--------------------------------------------------------------------------------------

  // for drawing some example quarks moving around - put in double array of TGraphs
  const int npoints = N_timesteps;
  TGraph *trackheavy[100][npoints];  // map the first 100 for examples
  int maxdraw_trackheavy = 20; // this really draws maxdraw/2 pairs
  for (int it1=0;it1<100;it1++) {
    for (int it2=0;it2<N_timesteps;it2++) {
      trackheavy[it1][it2] = new TGraph();
      trackheavy[it1][it2]->SetMarkerStyle(20);
      trackheavy[it1][it2]->SetMarkerSize(0.4);
      trackheavy[it1][it2]->SetLineColor(it1+1);
      trackheavy[it1][it2]->SetMarkerColor(it1+1);
      trackheavy[it1][it2]->SetLineColor(1);
      trackheavy[it1][it2]->SetMarkerColor(it1%2? kGray : kBlack);
    }
  }

  double qxstore = 0.0;
  double qystore = 0.0;
  double qptstore = 0.0;
  double qphistore = 0.0;
  double qphifinalstore = 0.0;
  double qptfinalstore = 0.0;


  cout << "Staring loop over quarks ..." << endl;
  for(int iquark=1; iquark<=NQuarks; iquark++) {
    
    cout << Form("Quark # %d/%d\r", iquark, NQuarks) << endl;
    
    // start with quark initial information
    // in odd numbered events, calculate both momentum for charm and anti-charm
    // and then use the anti-charm for the even numbered event.
    double qpt, qphi, qx, qy, qrapidity, qpz;
    if (iquark%2!=0) {

      qx = 0.0; 
      qy = 0.0;
      // Get initial x,y coordinate from the Ncoll spatial distribution
      quarkxy->GetRandom2(qx,qy);

      qpt = 0;
      qphi = ran.Uniform(-TMath::Pi(), TMath::Pi());
			if (yesbeauty) qpt = ptb->GetRandom();
			else qpt = ptc->GetRandom();

      if (sampleRapidity) qrapidity = ran.Uniform(-0.25, 0.25);
      else qrapidity = 0;

      qpz = sqrt(QuarkMass*QuarkMass + qpt*qpt) * sinh (qrapidity);

      tquark1->SetPxPyPzE( cos(qphi)*qpt, sin(qphi)*qpt , qpz, sqrt(QuarkMass*QuarkMass + qpt*qpt + qpz*qpz));
      tquark2->SetPxPyPzE( cos(qphi+TMath::Pi())*qpt, sin(qphi+TMath::Pi())*qpt , qpz, sqrt(QuarkMass*QuarkMass + qpt*qpt + qpz*qpz));

      // CONSIDER OPTION OF A KT KICK (SO NOT BACK TO BACK)
      // 1.  RANDOMLY SELECT A KT KICK FROM A GAUSSIAN with sigma = 1 GeV (for example)
      // 2.  in the rest frame of the c/cbar the total mass = 2 x charm quark mass
      // 3.  convert momentum kick in this rest frame into beta boost
      // 4.  then boost each c and cbar using TLorenzVectors

      double kt = fkt->GetRandom();
			//kt = 0.0;
      
      double DiQuarkMass = 2.0 * QuarkMass;
      double ktBetaBoost = TMath::Abs( kt / sqrt(pow(kt,2)+pow(DiQuarkMass,2)) );
      
      double BoostPhi = ran.Uniform(-TMath::Pi(), TMath::Pi());
      tquark1->Boost(ktBetaBoost*sin(BoostPhi),ktBetaBoost*cos(BoostPhi),0.0); 
      tquark2->Boost(ktBetaBoost*sin(BoostPhi),ktBetaBoost*cos(BoostPhi),0.0); 
      // then extract px, py again...
      qpt = tquark1->Pt();
      qphi = tquark1->Phi();// double check....
	
      // store anti-charm information to be grabbed in the next event
      qxstore = qx; // spatial coordinates
      qystore = qy;
      qptstore = tquark2->Pt(); // momentum and phi
      qphistore = tquark2->Phi();

      double deltaphi = tquark1->Phi() - tquark2->Phi();
      if (deltaphi < 0) deltaphi *= -1.0;
      if (deltaphi > pi) deltaphi = 2.0*pi - deltaphi;
      hdeltaphiorig->Fill(deltaphi);
      if (tquark1->Pt()>1.0 && tquark2->Pt()>1.0) hdeltaphiorig1->Fill(deltaphi);
      if (tquark1->Pt()>2.0 && tquark2->Pt()>2.0) hdeltaphiorig2->Fill(deltaphi);
      if (tquark1->Pt()>3.0 && tquark2->Pt()>3.0) hdeltaphiorig3->Fill(deltaphi);

    } else {
      // every other event make the same but opposite direction
      qx = qxstore;
      qy = qystore;
      qpt = qptstore;
      qphi = qphistore;
    }

    // all of the below code uses qx, qy, qpt, qphi... (so need to grab from alternating events above)

    // put into TLorentzVector
    tquark->SetPxPyPzE( cos(qphi)*qpt, sin(qphi)*qpt , 0.0, 
			sqrt(QuarkMass*QuarkMass + qpt*qpt));

    const double radius = sqrt (pow(qx,2) + pow(qy,2));
    int ir = 0;
    while (ir < 5 && radii[ir] < radius) ir++;

    const double init_pt = tquark->Pt(); // pt, phi at the beginning of the time evolution
    const double init_phi = tquark->Phi(); // returns between -pi and pi

		//continue;

    // now loop over time steps for diffusion
    for (int istep=0; istep < N_timesteps; istep++) {

      if (istep % (int) timeStepScale) continue;

      if (iquark <= maxdraw_trackheavy) {
	      for (int it1=istep;it1<N_timesteps;it1++) 
	        trackheavy[iquark-1][it1]->SetPoint(istep, qx, qy);
      }
      
      double local_temp  = 0;
      double local_betax = 0;
      double local_betay = 0;

      // from qx,qy coordinate, determine what cell the quark is in... 
      int ibinx = htemperature[istep]->GetXaxis()->FindBin(qx);
      int ibiny = htemperature[istep]->GetYaxis()->FindBin(qy);
      
      // make sure we are inside the grid bounaries
      if (ibinx>=1 && ibinx <=NGridX && ibiny >=1 && ibiny<=NGridY) {
	      local_temp  = htemperature[istep]->GetBinContent(ibinx,ibiny);
	      local_betax = hbetax[istep]->GetBinContent(ibinx,ibiny);
	      local_betay = hbetay[istep]->GetBinContent(ibinx,ibiny);
      }

      if (local_temp >= temperature_cutoff) {
	
	      // boost into rest frame of cell
	      tquark->Boost(-local_betax,-local_betay,0.0);
	      
	      double Dparameter = Dlookup->GetBinContent(Dlookup->FindBin(local_temp));

	      // if we wanted a momentum dependence to D - take arXiv 1005.0769v1 Figure 20 Riek,Rapp
	      // posit a 0.5 drop in coupling strength in going from 0 GeV to 5 GeV (simple linear for now)
	      if (MomentumDependent) {
	        // D parameter should be increasing with larger momentum (i.e. weaker coupling)
	        double scaleValue = 1.0 + 2.0*(tquark->P()/10.0); //  at P=0, sV=1; at P=5, sV = 2.0
	        Dparameter = Dparameter * scaleValue;
	      }

	      double px_kick = ran.Gaus(0, sqrt(2.0 * pow(local_temp,2) / (Dparameter * time_step_length_inverseGeV))) * time_step_length_inverseGeV;
	      double py_kick = ran.Gaus(0, sqrt(2.0 * pow(local_temp,2) / (Dparameter * time_step_length_inverseGeV))) * time_step_length_inverseGeV;
	      double pz_kick = ran.Gaus(0, sqrt(2.0 * pow(local_temp,2) / (Dparameter * time_step_length_inverseGeV))) * time_step_length_inverseGeV;

	      // could rotate and have a different kappa in and out of direction of motion...
	      
	      double quarkbetax = tquark->Px() / tquark->E();
	      double quarkbetay = tquark->Py() / tquark->E();
	      double quarkbetaz = tquark->Pz() / tquark->E();

	      // drag calculation
	      double drag_loss_x = time_step_length_inverseGeV * pow(local_temp,1) * quarkbetax / Dparameter;
	      double drag_loss_y = time_step_length_inverseGeV * pow(local_temp,1) * quarkbetay / Dparameter;
	      double drag_loss_z = time_step_length_inverseGeV * pow(local_temp,1) * quarkbetaz / Dparameter;

	      // put is all back together 
	      double newpx = tquark->Px() - drag_loss_x + px_kick;
	      double newpy = tquark->Py() - drag_loss_y + py_kick;
	      double newpz = tquark->Pz() - drag_loss_z + pz_kick;

	      // can set a non-interacting time
	      if (((double) istep)*time_step_length <  NonInteractTime) {
	        newpx = tquark->Px();  newpy = tquark->Py();  newpz = tquark->Pz();
	      }

	      tquark->SetPxPyPzE( newpx, newpy , newpz, sqrt(pow(QuarkMass,2)+pow(newpx,2) + pow(newpy,2) + pow(newpz,2)) );
	        
	      // boost back to normal frame (out of fluid rest frame)
	      tquark->Boost(+local_betax, +local_betay, 0.0);
	  
      } // even if below temperature cutoff still fill and just linearly propagate

      if (istep == 0) {
        hheavyptphi[0][0]->Fill (init_pt, init_phi);
        hheavypt[0][0]->Fill (init_pt);

        if (0 <= ir && ir < 4) {
          hheavyptphi[0][ir]->Fill (init_pt, init_phi);
          hheavypt[0][ir]->Fill (init_pt);
        }
      }

      if (plotThings || istep == N_timesteps-1) {
        hheavyptphi[istep+1][0]->Fill (tquark->Pt(), tquark->Phi()); // fill kinematics at the end of this time step
        hheavypt[istep+1][0]->Fill (tquark->Pt());

        if (0 <= ir && ir < 4) {
          hheavyptphi[istep+1][ir]->Fill (tquark->Pt(), tquark->Phi());
          hheavypt[istep+1][ir]->Fill (tquark->Pt());
        }
      }

      if (istep == N_timesteps-1) { // only fill these plots at the last timestep
        hheavyraa_particle[0]->Fill(init_pt, tquark->Pt()/init_pt);
        hheavydpt_particle[0]->Fill(init_pt, tquark->Pt()-init_pt);

        if (0 <= ir && ir < 4) {
          hheavyraa_particle[ir]->Fill(init_pt, tquark->Pt()/init_pt);
          hheavydpt_particle[ir]->Fill(init_pt, tquark->Pt()-init_pt);
        }
      }

      // update position of quark
      if (tquark->P() > 0) {
	      qx = qx + (tquark->Px()/tquark->E())* time_step_length;
	      qy = qy + (tquark->Py()/tquark->E())* time_step_length;
      } 
      
    } // goto next step

    if (iquark%2!=0) {
      qphifinalstore = tquark->Phi();
      qptfinalstore = tquark->Pt();
    } else {
      double deltaphi = qphifinalstore - tquark->Phi();
      if (deltaphi < 0) deltaphi *= -1.0;
      if (deltaphi > pi) deltaphi = 2.0*pi - deltaphi;
      hdeltaphi->Fill(deltaphi);
      if (tquark->Pt() > 1.0 && qptfinalstore > 1.0) hdeltaphi1->Fill(deltaphi);
      if (tquark->Pt() > 2.0 && qptfinalstore > 2.0) hdeltaphi2->Fill(deltaphi);
      if (tquark->Pt() > 3.0 && qptfinalstore > 3.0) hdeltaphi3->Fill(deltaphi);
    }

  } // end loop over quarks

  // output file
  //==============================================================
  TFile* outFile = new TFile(output_filename,"recreate");

  // plotting 
  //==============================================================
  if (plotThings) {
    cout << "\nDrawing..." << endl;
    
    TH1D *hratio = new TH1D("hratio","hratio",50,0.0,5.0);

    // NOW DRAW THINGS INTO THE CANVAS
    for (int itime=0; itime<N_timesteps; itime++) {

      if (itime % (int) timeStepScale) continue;
      
      // Panel 1: Flow field and quark trajectories
      cc->cd(1);

      Printf("step %d", itime);

      htemperature[itime]->Draw("colz");
      boosts[itime]->Draw();

      for(int iquark=0; iquark<maxdraw_trackheavy; iquark++) {
        if (!trackheavy[iquark][itime])
	        Error("", "no trackheavy[%d][%d]", iquark, itime);
        else
	        trackheavy[iquark][itime]->Draw("p,same");
      }
      ltx.DrawLatex(-9.0, 10.2, "Temperature");
      
      // Panel 2: pt distribution, original and modified
      if (!onepanel) cc->cd(2);
      TH1D* hheavypt_orig = hheavypt[0][0];
      hheavypt_orig->SetLineWidth(3);
      hheavypt_orig->SetLineColor(1);
      hheavypt_orig->SetXTitle("Transverse Momentum (GeV/c)");
      hheavypt_orig->SetYTitle("dN/dp_{T} (Initial=black,Current=red)");
      
      // determine the y-axis scale for this figure....
      double maxyrange = hheavypt_orig->GetBinContent(hheavypt_orig->GetMaximumBin());
      if (maxyrange < hheavypt[itime][0]->GetBinContent(hheavypt[itime][0]->GetMaximumBin())) 
        maxyrange = hheavypt[itime][0]->GetBinContent(hheavypt[itime][0]->GetMaximumBin());
      maxyrange = 1.1 * maxyrange;
      hheavypt_orig->SetMaximum(maxyrange);
      
      if (!onepanel) hheavypt_orig->DrawCopy();
      hheavypt[itime][0]->SetLineWidth(3);
      hheavypt[itime][0]->SetLineColor(2);
      if (!onepanel) hheavypt[itime][0]->DrawCopy("same");
      
      TLine *t1 = new TLine(hheavypt_orig->GetMean(),0.0,hheavypt_orig->GetMean(),10000.0);
      if (!onepanel) t1->Draw("same");
      TLine *t2 = new TLine(hheavypt[itime][0]->GetMean(),0.0,hheavypt[itime][0]->GetMean(),10000.0);
      t2->SetLineColor(2);
      if (!onepanel) t2->Draw("same");
      
      if (!onepanel) ltxn.DrawLatex(0.2, 0.9, Form("Time (fm/c): %s", htemperature[itime]->GetTitle()));

      // Panel 3: pt ratio plot ("R_AA")
      if (!onepanel) cc->cd(3);
      for (int ibin=1;ibin<=50;ibin++) {
        if (hheavypt_orig->GetBinContent(ibin) > 0) {
          hratio->SetBinContent(ibin, hheavypt[itime][0]->GetBinContent(ibin) / hheavypt_orig->GetBinContent(ibin));
	        hheavyraa[itime][0]->SetBinContent(ibin, hheavypt[itime][0]->GetBinContent(ibin) / hheavypt_orig->GetBinContent(ibin));
        }
      }

      hratio->SetLineWidth(3);
      hratio->SetLineColor(4);
      hratio->SetMinimum(0.0);
      hratio->SetMaximum(2.0);
      hratio->SetXTitle("Transverse Momentum (GeV/c)");
      hratio->SetYTitle("R_{AA} for this time step");
      maxyrange = hratio->GetBinContent(hratio->GetMaximumBin());
      if (maxyrange < 1.5)
        maxyrange = 1.5;
      maxyrange = 1.1 * maxyrange;
      hratio->SetMaximum(2.0); // try just leaving fixed

      if (!onepanel) hratio->DrawCopy();

      cc->Modified();
      cc->Update();
      cc->Draw();

      // options for saving out information
      if (save.Contains("last")) { 
        if (itime==N_timesteps-(int)timeStepScale) {
	        cc->Print(Form("VisualOut/a%.3g.gif",parameterA));
	        cc->SaveAs(Form("VisualOut/a%.3g.C",parameterA));
        }
        else if (itime==N_timesteps) {
	        cc->Print(Form("VisualOut/a%.3g.gif",parameterA));
	        cc->SaveAs(Form("VisualOut/a%.3g.C",parameterA));
	      }
      }
      else if (save.Contains("all")) { 
        cc->Print(Form("VisualOut/canvas_%03d.gif",itime));
      }
      
    } // timesteps

    TCanvas *c2 = new TCanvas();
    c2->cd();
    hdeltaphi1->SetLineWidth(4);
    hdeltaphi1->SetLineColor(1);
    hdeltaphi2->SetLineWidth(4);
    hdeltaphi2->SetLineColor(2);
    hdeltaphi3->SetLineWidth(4);
    hdeltaphi3->SetLineColor(3);
    hdeltaphi->SetMinimum(0.0);
    hdeltaphi->SetLineColor(4);
    hdeltaphi->SetFillColor(5);
    hdeltaphi->Scale(1.0/hdeltaphi->Integral());
    hdeltaphi->SetMinimum(0.0);
    hdeltaphi->DrawCopy("l,f");

    hdeltaphi1->Scale(1.0/hdeltaphi1->Integral());
    hdeltaphi1->DrawCopy("l,same");
    hdeltaphi2->Scale(1.0/hdeltaphi2->Integral());
    hdeltaphi2->DrawCopy("l,same");
    hdeltaphi3->Scale(1.0/hdeltaphi3->Integral());
    hdeltaphi3->DrawCopy("l,same");

    cc->Write("maincanvas");
    c2->Write("subcanvas");
  }

	//return;
  
  // write out pT distribution in time steps
  // write out raa in time steps
  // write out canvas of final display ?????
	//
  
  hdeltaphiorig->Write();
  hdeltaphiorig1->Write();
  hdeltaphiorig2->Write();
  hdeltaphiorig3->Write();
  hdeltaphi->Write();
  hdeltaphi1->Write();
  hdeltaphi2->Write();
  hdeltaphi3->Write();

  for (int ir = 0; ir < 5; ir++) {
    for (int itime = 0; itime <= N_timesteps; itime++) {
      if (plotThings || itime == 0 || itime == N_timesteps) {
        hheavypt[itime][ir]->Write();
        hheavyptphi[itime][ir]->Write();
        hheavyraa[itime][ir]->Write();
      }
    }
    hheavyraa_particle[ir]->Write();
    hheavydpt_particle[ir]->Write();
  }

  outFile->Close();

}  

void SetGraphProps(TGraph* g,
		   Int_t linecolor,
		   Int_t markercolor,
		   Int_t markerstyle,
		   Double_t markersize) 
{
  g->SetLineColor(linecolor);
  g->SetMarkerColor(markercolor);
  g->SetMarkerStyle(markerstyle);
  g->SetMarkerSize(markersize);
  g->SetLineWidth(2);
}
