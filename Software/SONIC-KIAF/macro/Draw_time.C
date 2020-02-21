#include "/phenix/u/shlim/Style.h"

void Draw_time(const int grp=49, const int run=0){

	gStyle->SetPalette(55);
	gStyle->SetOptStat(0);

	const bool bLOG = true;

	//const int ntime = 87;
	const int ntime = 172;

	const float Qs = 1.0;
	const float m = 0.0;

	//TFile *infile = new TFile(Form("running_grp%d/IPGlasma_flat_useNucleus0_grid1024_g2mu0.025_m0.15_run%05d.root",grp,run),"read");
	//TFile *infile = new TFile(Form("running_grp%d/IPGlasma_flat_useNucleus0_grid1024_g2mu0.10_m0.15_run%05d.root",grp,run),"read");
	//TFile *infile = new TFile(Form("running_grp%d/IPGlasma_dAu_useNucleus1_grid256_g2mu1.60_m0.0_run%05d.root",grp,run),"read");
	TFile *infile = new TFile(Form("running_grp%d/IPGlasma_dAu_useNucleus1_grid512_g2mu0.40_m0.0_run%05d.root",grp,run),"read");

	TH2D *h2d[ntime];

	for (int it=0; it<ntime; it++){
		h2d[it] = (TH2D*)infile->Get(Form("h2_evt%05d_t%05d",run,it));
	}

	TH1D *h1_time = (TH1D*)infile->Get("h1_time");

	float max = h2d[0]->GetMaximum();
	float min = h2d[ntime-1]->GetMinimum();
	float integral = h2d[0]->Integral();

	TCanvas *c1 = new TCanvas("c1","c1",1.2*1*500,500);

	//c1->cd(1);
	SetPadStyle();
	gPad->SetRightMargin(0.135);
	gPad->SetLeftMargin(0.13);
	gPad->SetTopMargin(0.08);
	if ( bLOG ) gPad->SetLogz();

	/*
	c1->cd(2);
	SetPadStyle();
	gPad->SetRightMargin(0.135);
	gPad->SetLeftMargin(0.13);
	gPad->SetTopMargin(0.08);
	if ( bLOG ) gPad->SetLogz();
	*/


	gSystem->Exec(Form("mkdir -p running_grp%d/gif_run%05d",grp,run));

	for (int it=0; it<ntime; it++){
	//for (int it=0; it<1; it++){
		//c1->cd(1);
		c1->Clear();

		//c1->cd(2);
		//c1->Clear();

		if ( h2d[it]->Integral()<0.01*integral ) continue;

		//htmp = (TH1D*)gPad->DrawFrame(-15,-15,15,15);
		htmp = (TH1D*)gPad->DrawFrame(-7.5,-7.5,7.5,7.5);
		SetHistoStyle("x [fm]", "y [fm]");
		htmp->GetXaxis()->CenterTitle();
		htmp->GetYaxis()->CenterTitle();
		htmp->GetYaxis()->SetTitleOffset(1.0);

		if ( bLOG ){
			h2d[it]->SetMinimum(min);
			h2d[it]->SetMaximum(max);
		}else{
			h2d[it]->SetMaximum(1.02*max);
			h2d[it]->SetMinimum(0);
		}
		htmp = (TH1D*)h2d[it];
		SetHistoStyle("","","",0.04,0.04);

		float time = h1_time->GetBinContent(it);
		TLatex *tex = new TLatex(0,6.5,Form("Q_{s}=%g GeV, m=%g GeV, t=%7.5f fm/c",Qs,m,time));
		tex->SetTextFont(63);
		tex->SetTextSize(20);
		tex->SetTextAlign(21);
		tex->Draw();

		h2d[it]->Draw("colz same");
		if ( bLOG ){
			c1->SaveAs(Form("running_grp%d/gif_run%05d/time_log_%05d.gif",grp,run,it));
		}else{
			c1->SaveAs(Form("running_grp%d/gif_run%05d/time_linear_%05d.gif",grp,run,it));
		}
	}



}
