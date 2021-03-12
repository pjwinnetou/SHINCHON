#include <iostream>
#include "Style_jaebeom.h"
#include "TreeSetting.h"

using namespace std;

void drawBasicPlots(int nrun=10, int kInitPos=1, int phiN = 3)
{

  setTDRStyle();
  writeExtraText= false;
  int iPeriod = 999;
  int iPos = 33;

  gStyle->SetEndErrorSize(0);
  gStyle->SetOptFit(kTRUE);
  
  const char* fPosStr;
  if(kInitPos==0) fPosStr = "InitPosZero";
  else if(kInitPos==1) fPosStr = "InitPosGlauber";
  else if(kInitPos==2) fPosStr = "InitPosMean";

  TFile* rfout = new TFile(Form("outfile_UpsSkim_PhiAng%d_%s_0000_%04d.root",phiN,fPosStr,nrun),"read");
  TTree* tree = (TTree*) rfout->Get("tree");
  if(nrun!= tree->GetEntries()){cout << "ERROR!! :::: Number of entries and runs inconsistent!!" << endl;return;}
  SetTree settree_;
  settree_.TreeSetting(tree);
  
  //Glauber
  TFile* glf = new TFile("../MCGlauber-PbPb-5020GeV-b0-18fm-bin100-v3.root","read");

  string savedir = Form("res/plots/%s_nRun%d",fPosStr,nrun);
  void * dirf = gSystem->OpenDirectory(savedir.c_str());
  if(dirf) gSystem->FreeDirectory(dirf);
  else {gSystem->mkdir(savedir.c_str(), kTRUE);}

  double lab_posx_ = 0.322; double lab_posy_ = 0.664; double lab_pos_diff_ = 0.076;
  TLatex* globtex_int = new TLatex();
  globtex_int->SetNDC();
  globtex_int->SetTextAlign(12); //left-center
  globtex_int->SetTextFont(42);
  globtex_int->SetTextSize(0.036);

  TCanvas* c_tempEnergy[nrun];
  TCanvas* c_PosInit[nrun];
  TCanvas* c_PosFinal[nrun];
  TH1D* h_tempEnergy[nrun];
  TH1D* h_PosInit[nrun];
  TH1D* h_PosFinal[nrun];

  TH1D* hphi_initMom = (TH1D*) rfout->Get("hphi_initMom");
  TH1D* hphi_finalMom = (TH1D*) rfout->Get("hphi_finalMom");

  SetHistAxis(hphi_initMom,"#phi (rad)","#frac{dN}{d#phi}");
  SetHistAxis(hphi_finalMom,"#phi (rad)","#frac{dN}{d#phi}");

  hphi_initMom->Sumw2();
  hphi_finalMom->Sumw2();
  hphi_initMom->GetYaxis()->SetMaxDigits(3);

  TH1D* hphi_rat = (TH1D*) hphi_initMom->Clone("hphi_rat");
  hphi_rat->Divide(hphi_finalMom);
  SetHistAxis(hphi_rat,"#phi (rad)","ratio");
  for(int ibin=1;ibin<=hphi_initMom->GetNbinsX();ibin++){
    double val1 = hphi_initMom->GetBinContent(ibin);
    double val2 = hphi_finalMom->GetBinContent(ibin);
    double valErr1 = hphi_initMom->GetBinError(ibin);
    double valErr2 = hphi_finalMom->GetBinError(ibin);
    double err = getErrorPropaDivide(val1, valErr1, val2, valErr2);
    hphi_rat->SetBinError(ibin,err);
  }

  //Make Canvas
  double ymax = 1.05; double ymin = 0.95; bool setLogopt = false;
  TCanvas* cr = makeHistRatioCanvas(hphi_rat, hphi_initMom, hphi_finalMom, ymax, ymin, setLogopt);
  TPad *pad1 = (TPad*) cr->GetPad(0); 
  pad1->Draw();
  pad1->cd();
  
  double posx1 = 0.4; double posy1 = 0.78; double posx2 = 0.6; double posy2 = 0.9;
  TLegend* leg = new TLegend(posx1, posy1, posx2, posy2);
  SetLegendStyle(leg);
  leg->AddEntry(hphi_initMom,"Initial momentum #phi","pe");
  leg->AddEntry(hphi_finalMom,"Final momentum #phi","pe");
  leg->Draw("same"); 
  
  drawGlobText(fPosStr,0.21,0.872,1,0.027);
  pad1->Update();
  cr->Update();
  
  TPad *pad2 = (TPad*) cr->GetPad(1); 
  pad2->Draw();
  pad2->cd();
  TLegend* legr = new TLegend(posx1, posy1+0.05, posx2, posy2+0.02);
  SetLegendStyle(legr);
  legr->SetTextSize(0.07);
  legr->AddEntry(hphi_rat,"Ratio","pe");
  legr->Draw("same"); 
  
  SHINCHONLegend( pad1, 2, iPos,1);
  cr->Modified();
  cr->Update();
  cr->SaveAs(Form("%s/hist_UpsPhiDist.pdf",savedir.c_str()));


  for(int ih=0; ih<nrun; ih++){
    c_tempEnergy[ih] = new TCanvas(Form("c_tempEnergy_%d",ih),"",600,600);
    c_PosInit[ih]    = new TCanvas(Form("c_PosInit_%d",ih),"",600,600);
    c_PosFinal[ih]   = new TCanvas(Form("c_PosFinal_%d",ih),"",600,600);
    SetCanvasSquare(c_tempEnergy[ih]);
    SetCanvasSquare(c_PosInit[ih]);
    SetCanvasSquare(c_PosFinal[ih]);

    tree->GetEntry(ih);
    h_tempEnergy[ih] = (TH1D*) rfout->Get(Form("histTempEnergy_%d",ih));
    h_PosInit[ih] = (TH1D*) rfout->Get(Form("hPosInit_run%d",ih));
    h_PosFinal[ih] = (TH1D*) rfout->Get(Form("hPosFinal_run%d",ih));
    
    SetHistAxis(h_tempEnergy[ih],"x (fm)","y (fm)");
    SetHistAxis(h_PosInit[ih],"x (fm)","y (fm)");
    SetHistAxis(h_PosFinal[ih],"x (fm)","y (fm)");

    c_tempEnergy[ih]->cd();
    h_tempEnergy[ih]->Draw("colz");
    globtex_int->DrawLatex(lab_posx_,lab_posy_,"#bf{Energy density profile}"); 
    globtex_int->DrawLatex(lab_posx_,lab_posy_-lab_pos_diff_,Form("Npart = %d",Npart)); 
    globtex_int->DrawLatex(lab_posx_,lab_posy_-2*lab_pos_diff_,Form("b = %.2f fm",b)); 
    globtex_int->DrawLatex(lab_posx_,lab_posy_-3*lab_pos_diff_,Form("#Psi_{%d,EP}(En profile) = %.2f #circ",phiN,EPangEnProf*180/TMath::Pi())); 
    globtex_int->DrawLatex(lab_posx_,lab_posy_-4*lab_pos_diff_,Form("#Psi_{%d,EP}(Glauber) = %.2f #circ",phiN,EPangGlauber*180/TMath::Pi())); 
    SHINCHONLegend( c_tempEnergy[ih], iPeriod, iPos, 1 );
  
    c_PosInit[ih]->cd();
    h_PosInit[ih]->Draw("colz");
    globtex_int->DrawLatex(lab_posx_,lab_posy_,"#bf{Initial #varUpsilon(1S) position}"); 
    globtex_int->DrawLatex(lab_posx_,lab_posy_-lab_pos_diff_,Form("Npart = %d",Npart)); 
    globtex_int->DrawLatex(lab_posx_,lab_posy_-2*lab_pos_diff_,Form("b = %.2f fm",b)); 
    globtex_int->DrawLatex(lab_posx_,lab_posy_-3*lab_pos_diff_,Form("#Psi_{%d,EP}(En profile) = %.2f #circ",phiN,EPangEnProf*180/TMath::Pi())); 
    globtex_int->DrawLatex(lab_posx_,lab_posy_-4*lab_pos_diff_,Form("#Psi_{%d,EP}(Glauber) = %.2f #circ",phiN,EPangGlauber*180/TMath::Pi())); 
    SHINCHONLegend( c_PosInit[ih], iPeriod, iPos,0 );

    c_PosFinal[ih]->cd();
    h_PosFinal[ih]->Draw("colz");
    globtex_int->DrawLatex(lab_posx_,lab_posy_,"#bf{Final #varUpsilon(1S) position}"); 
    globtex_int->DrawLatex(lab_posx_,lab_posy_-lab_pos_diff_,Form("Npart = %d",Npart)); 
    globtex_int->DrawLatex(lab_posx_,lab_posy_-2*lab_pos_diff_,Form("b = %.2f fm",b)); 
    globtex_int->DrawLatex(lab_posx_,lab_posy_-3*lab_pos_diff_,Form("#Psi_{2,EP}(En profile) = %.2f #circ",EPangEnProf*180/TMath::Pi())); 
    globtex_int->DrawLatex(lab_posx_,lab_posy_-4*lab_pos_diff_,Form("#Psi_{2,EP}(Glauber) = %.2f #circ",EPangGlauber*180/TMath::Pi())); 
    SHINCHONLegend( c_PosFinal[ih], iPeriod, iPos, 1 );

    c_tempEnergy[ih]->SaveAs(Form("%s/hist_TempEn_phiAng%d_run%d.pdf",savedir.c_str(),phiN,ih));
    c_PosInit[ih]->SaveAs(Form("%s/hist_UpsInitPos_phiAng%d_run%d.pdf",savedir.c_str(),phiN,ih));
    c_PosFinal[ih]->SaveAs(Form("%s/hist_UpsFinalPos_phiAng%d_run%d.pdf",savedir.c_str(),phiN,ih));

    TH2D* hinited_event = (TH2D*) glf->Get(Form("inited_event%d",ih));
    SetHistAxis(hinited_event,"x (fm)","y (fm)");
    TCanvas* c_glb = new TCanvas("c_glb","",600,600);
    SetCanvasSquare(c_glb);
    c_glb->cd();
    hinited_event->Draw("colz");
    globtex_int->DrawLatex(lab_posx_,lab_posy_,Form("Npart = %d, Ncoll = %d",Npart,Ncoll)); 
    globtex_int->DrawLatex(lab_posx_,lab_posy_-lab_pos_diff_,Form("b = %.2f fm",b)); 
    globtex_int->DrawLatex(lab_posx_,lab_posy_-2*lab_pos_diff_,Form("#epsilon_{%d} = %.3f",phiN,eccgaus[2])); 
    globtex_int->DrawLatex(lab_posx_,lab_posy_-3*lab_pos_diff_,Form("#Psi_{%d,EP}(Glauber) = %.2f #circ",phiN,EPangGlauber*180/TMath::Pi())); 
    SHINCHONLegend( c_glb, iPeriod, iPos, 1); 
    c_glb->SaveAs(Form("%s/hist_GlauberEvt_phiAng%d_run%d.pdf",savedir.c_str(),phiN,ih));

    delete hinited_event;
    delete c_glb;
  }

}
