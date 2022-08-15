#include <iostream>
#include "../../Style_jaebeom_woFrameLegend.h"
#include "../../SHINCHONLegend_raaCent.C"
#include "../../tdrstyle.C"

using namespace std;

void FinalDrawRAA_PbPb(){

	const bool bSAVE = true;

	gStyle->SetOptStat(0);
	gStyle->SetPalette(55);
	TColor::InvertPalette();
	setTDRStyle();

	//TFile *infile = new TFile("outfile_RaaV7_0000_1000.root","read");
	TFile *infile = new TFile("outfile_RaaV7_0000_1000_corrected.root","read");

	TProfile *hprofRAA1S_Npart = (TProfile*)infile->Get("hprofRAA1S_Npart_allpT");
	TProfile *hprofRAA2S_Npart = (TProfile*)infile->Get("hprofRAA2S_Npart_allpT");
	TProfile *hprofRAA3S_Npart = (TProfile*)infile->Get("hprofRAA3S_Npart_allpT");

	TH1D *hprofRAA1S_Npart_c = (TH1D*)infile->Get("hRpA_pPb_mult_corr_1S");
	TH1D *hprofRAA2S_Npart_c = (TH1D*)infile->Get("hRpA_pPb_mult_corr_2S");
	TH1D *hprofRAA3S_Npart_c = (TH1D*)infile->Get("hRpA_pPb_mult_corr_3S");

	//hprofRAA1S_Npart->Rebin(10);
	//hprofRAA2S_Npart->Rebin(10);
	//hprofRAA3S_Npart->Rebin(10);

	TGraphErrors *gRAA1S_Npart = new TGraphErrors;
	TGraphErrors *gRAA2S_Npart = new TGraphErrors;
	TGraphErrors *gRAA3S_Npart = new TGraphErrors;

	TGraphErrors *gRAA1S_Npart_c = new TGraphErrors;
	TGraphErrors *gRAA2S_Npart_c = new TGraphErrors;
	TGraphErrors *gRAA3S_Npart_c = new TGraphErrors;

	for (int ix=0; ix<hprofRAA1S_Npart->GetNbinsX(); ix++){
		float xx = hprofRAA1S_Npart->GetBinCenter(ix+1);
		float yy = hprofRAA1S_Npart->GetBinContent(ix+1);
		float yyerr = hprofRAA1S_Npart->GetBinError(ix+1);

		gRAA1S_Npart->SetPoint(ix, xx, yy);
		gRAA1S_Npart->SetPointError(ix, 0, yyerr);

		xx = hprofRAA2S_Npart->GetBinCenter(ix+1);
		yy = hprofRAA2S_Npart->GetBinContent(ix+1);
		yyerr = hprofRAA2S_Npart->GetBinError(ix+1);

		gRAA2S_Npart->SetPoint(ix, xx, yy);
		gRAA2S_Npart->SetPointError(ix, 0, yyerr);

		xx = hprofRAA3S_Npart->GetBinCenter(ix+1);
		yy = hprofRAA3S_Npart->GetBinContent(ix+1);
		yyerr = hprofRAA3S_Npart->GetBinError(ix+1);

		gRAA3S_Npart->SetPoint(ix, xx, yy);
		gRAA3S_Npart->SetPointError(ix, 0, yyerr);

		xx = hprofRAA1S_Npart_c->GetBinCenter(ix+1);
		yy = hprofRAA1S_Npart_c->GetBinContent(ix+1);
		yyerr = hprofRAA1S_Npart_c->GetBinError(ix+1);

		gRAA1S_Npart_c->SetPoint(ix, xx, yy);
		gRAA1S_Npart_c->SetPointError(ix, 0, yyerr);

		xx = hprofRAA2S_Npart_c->GetBinCenter(ix+1);
		yy = hprofRAA2S_Npart_c->GetBinContent(ix+1);
		yyerr = hprofRAA2S_Npart_c->GetBinError(ix+1);

		gRAA2S_Npart_c->SetPoint(ix, xx, yy);
		gRAA2S_Npart_c->SetPointError(ix, 0, yyerr);

		xx = hprofRAA3S_Npart_c->GetBinCenter(ix+1);
		yy = hprofRAA3S_Npart_c->GetBinContent(ix+1);
		yyerr = hprofRAA3S_Npart_c->GetBinError(ix+1);

		gRAA3S_Npart_c->SetPoint(ix, xx, yy);
		gRAA3S_Npart_c->SetPointError(ix, 0, yyerr);
	}

	TFile *infile1 = new TFile("../../CMSPbPbPublicResults/RAA_Y1S_Npart_PbPb5TeV_HEPData-ins1674529-v2-Table_19.root","read");
	TDirectoryFile *tdf = (TDirectoryFile*)infile1->Get("Table 19");

	TH1D *Hist1D_y1 = (TH1D*)tdf->Get("Hist1D_y1");
	TH1D *Hist1D_y1_e1 = (TH1D*)tdf->Get("Hist1D_y1_e1");
	TH1D *Hist1D_y1_e2p = (TH1D*)tdf->Get("Hist1D_y1_e2plus");
	TH1D *Hist1D_y1_e2m = (TH1D*)tdf->Get("Hist1D_y1_e2minus");
	TH1D *Hist1D_y2 = (TH1D*)tdf->Get("Hist1D_y2");

	TGraphErrors *gCMS = new TGraphErrors;
	TGraphAsymmErrors *gCMSsys = new TGraphAsymmErrors;

	for (int ix=0; ix<Hist1D_y1->GetNbinsX(); ix++){
		float xx = Hist1D_y2->GetBinContent(ix+1);
		float yy = Hist1D_y1->GetBinContent(ix+1);
		float yy_err = Hist1D_y1_e1->GetBinContent(ix+1);
		float yy_syserr0 = Hist1D_y1_e2m->GetBinContent(ix+1);
		float yy_syserr1 = Hist1D_y1_e2p->GetBinContent(ix+1);

		gCMS->SetPoint(ix, xx, yy);
		gCMS->SetPointError(ix, 0, yy_err);

		gCMSsys->SetPoint(ix, xx, yy);
		gCMSsys->SetPointError(ix, 5, 5, fabs(yy_syserr0), yy_syserr1);
	}

	TFile *infile2 = new TFile("../../CMSPbPbPublicResults/RAA_HIN21007_prelim.root","read");
	TGraphErrors *gCMS2 = (TGraphErrors*)infile2->Get("RAA_2S");
	TGraphErrors *gCMS3 = (TGraphErrors*)infile2->Get("RAA_3S");
	TGraphErrors *gCMS2sys = (TGraphErrors*)infile2->Get("RAA_2S_sys");
	TGraphErrors *gCMS3sys = (TGraphErrors*)infile2->Get("RAA_3S_sys");

	//return;
	{
		//w/o FD correction
		TCanvas *c1 = new TCanvas("c1","c1",1.2*500,500);
		gPad->SetLeftMargin(0.12);

		//SetPadStyle();

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,400,1.2);
		SetHistAxis(htmp,"N_{part}","R_{AA}");
    htmp->GetYaxis()->SetTitleOffset(1.0);
    htmp->GetXaxis()->SetTitleOffset(1.1);
    htmp->GetXaxis()->SetTitleSize(0.055);
    htmp->GetXaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetNdivisions(510);

		gCMSsys->SetFillColorAlpha(1,0.3);
		gCMSsys->SetLineWidth(0);
		gCMSsys->Draw("2");

		gCMS->SetMarkerStyle(20);
		gCMS->SetLineWidth(2);
		gCMS->Draw("P");

		gCMS2sys->SetFillColorAlpha(4,0.3);
		gCMS2sys->SetLineWidth(0);
		gCMS2sys->Draw("2");

		gCMS2->SetMarkerStyle(20);
		gCMS2->SetMarkerColor(4);
		gCMS2->SetLineWidth(2);
		gCMS2->SetLineColor(4);
		gCMS2->Draw("P");

		gCMS3sys->SetFillColorAlpha(8,0.3);
		gCMS3sys->SetLineWidth(0);
		gCMS3sys->Draw("2");

		gCMS3->SetMarkerStyle(20);
		gCMS3->SetMarkerColor(kGreen+2);
		gCMS3->SetLineWidth(2);
		gCMS3->SetLineColor(kGreen+2);
		gCMS3->Draw("P");

		gRAA1S_Npart->SetLineWidth(4);
		gRAA1S_Npart->SetLineStyle(2);
		gRAA1S_Npart->SetLineColorAlpha(1, 0.3);
		gRAA1S_Npart->SetFillColorAlpha(1, 0.3);
		gRAA1S_Npart->Draw("L3");

		gRAA2S_Npart->SetLineWidth(4);
		gRAA2S_Npart->SetLineStyle(2);
		gRAA2S_Npart->SetLineColorAlpha(4, 0.3);
		gRAA2S_Npart->SetFillColorAlpha(4, 0.3);
		gRAA2S_Npart->Draw("L3");

		gRAA3S_Npart->SetLineWidth(4);
		gRAA3S_Npart->SetLineStyle(2);
		gRAA3S_Npart->SetLineColorAlpha(kGreen+2, 0.5);
		gRAA3S_Npart->SetFillColorAlpha(kGreen+2, 0.5);
		gRAA3S_Npart->Draw("L3");

		TLegend *leg = new TLegend(0.5,0.6,0.8,0.9);
		SetLegendStyle(leg);
		//leg->SetFillStyle(0);
		//leg->SetBorderSize(0);
		leg->SetTextFont(43);
		leg->SetTextSize(20);
		leg->SetMargin(0.3);
		leg->AddEntry("","Pb-Pb #sqrt{s_{NN}}=5.02 TeV","h");
		leg->SetNColumns(2);
		leg->AddEntry("","SHINCHON","");
		leg->AddEntry("","CMS","");
		leg->AddEntry(gRAA1S_Npart,"#varUpsilon(1S)","LF");
		leg->AddEntry(gCMS,"#varUpsilon(1S)","P");
		leg->AddEntry(gRAA2S_Npart,"#varUpsilon(1S)","LF");
		leg->AddEntry(gCMS2,"#varUpsilon(2S)","P");
		leg->AddEntry(gRAA3S_Npart,"#varUpsilon(1S)","LF");
		leg->AddEntry(gCMS3,"#varUpsilon(3S)","P");
		leg->Draw();

	}

	{
		//w/o FD correction
		TCanvas *c2 = new TCanvas("c2","c2",1.2*500,500);
		gPad->SetLeftMargin(0.13);

		//SetPadStyle();

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,400,1.2);
		SetHistAxis(htmp,"N_{part}","R_{AA}");
    htmp->GetYaxis()->SetTitleOffset(1.0);
    htmp->GetXaxis()->SetTitleOffset(1.1);
    htmp->GetXaxis()->SetTitleSize(0.055);
    htmp->GetXaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetNdivisions(510);

		gCMSsys->Draw("2");
		gCMS->Draw("P");

		gCMS2sys->Draw("2");
		gCMS2->Draw("P");

		gCMS3sys->Draw("2");
		gCMS3->Draw("P");

		gRAA1S_Npart_c->SetLineWidth(4);
		gRAA1S_Npart_c->SetLineStyle(2);
		gRAA1S_Npart_c->SetLineColorAlpha(1, 0.3);
		gRAA1S_Npart_c->SetFillColorAlpha(1, 0.3);
		gRAA1S_Npart_c->Draw("L3");

		gRAA2S_Npart_c->SetLineWidth(4);
		gRAA2S_Npart_c->SetLineStyle(2);
		gRAA2S_Npart_c->SetLineColorAlpha(4, 0.3);
		gRAA2S_Npart_c->SetFillColorAlpha(4, 0.3);
		gRAA2S_Npart_c->Draw("L3");

		gRAA3S_Npart_c->SetLineWidth(4);
		gRAA3S_Npart_c->SetLineStyle(2);
		gRAA3S_Npart_c->SetLineColorAlpha(kGreen+2, 0.5);
		gRAA3S_Npart_c->SetFillColorAlpha(kGreen+2, 0.5);
		gRAA3S_Npart_c->Draw("L3");

		TLegend *leg1 = new TLegend(0.5,0.6,0.8,0.9);
		SetLegendStyle(leg1);
		//leg1->SetFillStyle(0);
		//leg1->SetBorderSize(0);
		leg1->SetTextFont(43);
		leg1->SetTextSize(20);
		leg1->SetMargin(0.3);
		leg1->AddEntry("","Pb-Pb #sqrt{s_{NN}}=5.02 TeV","h");
		leg1->SetNColumns(2);
		leg1->AddEntry("","SHINCHON","");
		leg1->AddEntry("","CMS","");
		leg1->AddEntry(gRAA1S_Npart_c,"#varUpsilon(1S)","LF");
		leg1->AddEntry(gCMS,"#varUpsilon(1S)","P");
		leg1->AddEntry(gRAA2S_Npart_c,"#varUpsilon(1S)","LF");
		leg1->AddEntry(gCMS2,"#varUpsilon(2S)","P");
		leg1->AddEntry(gRAA3S_Npart_c,"#varUpsilon(1S)","LF");
		leg1->AddEntry(gCMS3,"#varUpsilon(3S)","P");
		leg1->Draw();

		if ( bSAVE ){
			c2->SaveAs("RAA_Npart_PbPb_5TeV.pdf");
		}

	}


}
