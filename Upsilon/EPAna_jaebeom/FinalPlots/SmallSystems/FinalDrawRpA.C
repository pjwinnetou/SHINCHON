#include <iostream>
#include "../../Style_jaebeom_woFrameLegend.h"
#include "../../SHINCHONLegend_raaCent.C"
#include "../../tdrstyle.C"

using namespace std;

void FinalDrawRpA(){

	const bool bSAVE = true;

	gStyle->SetOptStat(0);
	gStyle->SetPalette(55);
	TColor::InvertPalette();
	setTDRStyle();

	const int nsyst = 3;
	const int nstat = 3; 

	const char fname[nsyst][100] = {
		"pPb",
		"pO",
		"OO" };

	const char dname[nsyst][100] = {
		"p#font[122]{-}Pb",
		"p#font[122]{-}O",
		"O#font[122]{-}O" };

	const float multmax[nsyst] = {50, 25, 75};

	int nColor[nstat] = {1, 4, kGreen+2}; 

	TFile* fin = new TFile("NMFdata/FDCOR_out.root","read");

	TGraphErrors* gRpPbMult[nsyst][nstat];
	TGraphErrors* gRpPbMultST[nsyst][nstat];
	TGraphErrors* gRpPbPt[nsyst][nstat];

	for(int i=0;i<nsyst;i++){
		for(int j=0;j<nstat;j++){
			gRpPbMult[i][j] = (TGraphErrors*)fin->Get(Form("g_R%s_multdep_fdcor_%dS",fname[i],j+1));
			gRpPbPt[i][j] = (TGraphErrors*)fin->Get(Form("g_R%s_ptdep_fdcor_%dS",fname[i],j+1));
			gRpPbMultST[i][j] = (TGraphErrors*)fin->Get(Form("g_R%s_multSTdep_fdcor_%dS",fname[i],j+1));		
		}
	}

	TCanvas *c1[nsyst];

	for(int i=0;i<nsyst;i++){
		c1[i] = new TCanvas(Form("c1_%d",i),Form("c1_%d",i),1.2*500,500);
		gPad->SetLeftMargin(0.12);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,multmax[i],1.5);
		htmp->GetXaxis()->SetTitle("dN_{ch}/d#eta");
		htmp->GetYaxis()->SetTitle("Nuclear modification factor");
    htmp->GetYaxis()->SetTitleOffset(1.0);
    htmp->GetXaxis()->SetTitleOffset(1.05);
    htmp->GetXaxis()->SetTitleSize(0.055);
    htmp->GetXaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetNdivisions(510);

		TLegend *leg = new TLegend(0.45,0.93-0.06*4,0.85,0.93);
		SetLegendStyle(leg);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(43);
		leg->SetTextSize(20);
		leg->SetMargin(0.3);
		leg->AddEntry( (TObject*)0, Form("SHINCHON, %s #sqrt{#it{s}_{NN}} = 8 TeV",dname[i]), "h");

		for(int j=0;j<nstat;j++){
			gRpPbMult[i][j]->SetFillColorAlpha(nColor[j], 0.5);
			gRpPbMult[i][j]->SetLineColorAlpha(nColor[j], 0.5);
			gRpPbMult[i][j]->SetLineWidth(4);
			gRpPbMult[i][j]->SetLineStyle(1);

			gRpPbMult[i][j]->Draw("L3");
			leg->AddEntry( gRpPbMult[i][j], Form("#varUpsilon(%dS)",j+1), "LF");
		}

		leg->Draw();

		if ( bSAVE ){
			c1[i]->cd();
			c1[i]->SaveAs(Form("NMFplots/RpPb_mult_%s.pdf",fname[i]));
		}
	}

	TCanvas *c2[nstat];

	for(int i=0;i<nstat;i++){
		c2[i] = new TCanvas(Form("c2_%d",i),Form("c2_%d",i),1.2*500,500);
		gPad->SetLeftMargin(0.12);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,multmax[i],1.5);
		htmp->GetXaxis()->SetTitle("dN_{ch}/d#eta");
		htmp->GetYaxis()->SetTitle("Nuclear modification factor");
    htmp->GetYaxis()->SetTitleOffset(1.0);
    htmp->GetXaxis()->SetTitleOffset(1.05);
    htmp->GetXaxis()->SetTitleSize(0.055);
    htmp->GetXaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetNdivisions(510);

		TLegend *leg = new TLegend(0.45,0.93-0.06*4,0.85,0.93);
		SetLegendStyle(leg);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(43);
		leg->SetTextSize(20);
		leg->SetMargin(0.3);
		leg->AddEntry( (TObject*)0, Form("SHINCHON, #sqrt{#it{s}_{NN}} = 8 TeV, #varUpsilon(%dS)",i+1), "h");

		for(int j=0;j<nsyst;j++){

			TGraphErrors *gcopy = new TGraphErrors(*gRpPbMult[j][i]);
			gcopy->SetFillColorAlpha(nColor[j], 0.5);
			gcopy->SetLineColorAlpha(nColor[j], 0.5);
			gcopy->SetLineWidth(4);
			gcopy->Draw("L3");
			leg->AddEntry(gcopy, Form("%s",dname[j]), "LF");
		}

		leg->Draw();

		if ( bSAVE ){
			c2[i]->cd();
			c2[i]->SaveAs(Form("NMFplots/RpPb_mult_%dS.pdf",i+1));
		}
	}

	TCanvas *c3[nsyst];

	for(int i=0;i<nsyst;i++){
		c3[i] = new TCanvas(Form("c3_%d",i),Form("c3_%d",i),1.2*500,500);
		gPad->SetLeftMargin(0.12);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,20,1.5);
		htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		htmp->GetYaxis()->SetTitle("Nuclear modification factor");
    htmp->GetYaxis()->SetTitleOffset(1.0);
    htmp->GetXaxis()->SetTitleOffset(1.05);
    htmp->GetXaxis()->SetTitleSize(0.055);
    htmp->GetXaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetNdivisions(510);

		TLegend *leg = new TLegend(0.2,0.93-0.06*4,0.6,0.93);
		SetLegendStyle(leg);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(43);
		leg->SetTextSize(20);
		leg->SetMargin(0.3);
		leg->AddEntry( (TObject*)0, Form("SHINCHON, %s #sqrt{#it{s}_{NN}} = 8 TeV",dname[i]), "h");

		for(int j=0;j<nstat;j++){
			gRpPbPt[i][j]->SetFillColorAlpha(nColor[j], 0.5);
			gRpPbPt[i][j]->SetLineColorAlpha(nColor[j], 0.5);
			gRpPbPt[i][j]->SetLineWidth(4);
			gRpPbPt[i][j]->SetLineStyle(1);

			gRpPbPt[i][j]->Draw("L3");
			leg->AddEntry( gRpPbMult[i][j], Form("#varUpsilon(%dS)",j+1), "LF");
		}

		leg->Draw();

		if ( bSAVE ){
			c3[i]->cd();
			c3[i]->SaveAs(Form("NMFplots/RpPb_Pt_%s.pdf",fname[i]));
		}
	}

	TFile* fCMS = new TFile("NMFdata//HEPData-ins2037640-v1-root.root","read");
	TDirectoryFile* dir;
	TH1F* hRpPb_CMS_cntl[nstat];
	TH1F* hRpPb_CMS_stat[nstat];
	TH1F* hRpPb_CMS_syst[nstat];

	for(int i=0;i<nstat;i++){
		dir = (TDirectoryFile*)fCMS->GetDirectory(Form("Table %d",13+i));
		hRpPb_CMS_cntl[i] = (TH1F*)dir->Get("Hist1D_y1");
		hRpPb_CMS_stat[i] = (TH1F*)dir->Get("Hist1D_y1_e1");
		hRpPb_CMS_syst[i] = (TH1F*)dir->Get("Hist1D_y1_e2");

		for(int j=0;j<hRpPb_CMS_cntl[i]->GetNbinsX();j++){
			hRpPb_CMS_stat[i]->SetBinContent( j+1, hRpPb_CMS_cntl[i]->GetBinContent(j+1) );
			hRpPb_CMS_syst[i]->SetBinContent( j+1, hRpPb_CMS_cntl[i]->GetBinContent(j+1) );
		}

		hRpPb_CMS_stat[i]->SetLineColor(nColor[i]);
		hRpPb_CMS_syst[i]->SetLineColor(nColor[i]);
		hRpPb_CMS_stat[i]->SetMarkerColor(nColor[i]);
		hRpPb_CMS_syst[i]->SetMarkerColor(nColor[i]);
		hRpPb_CMS_stat[i]->SetMarkerStyle(20);
		hRpPb_CMS_syst[i]->SetFillStyle(0);
	}

	{
		TCanvas *c5 = new TCanvas("c5","c5",1.2*500,500);
		gPad->SetLeftMargin(0.12);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,30,1.8);
		htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		htmp->GetYaxis()->SetTitle("Nuclear modification factor");
		htmp->GetYaxis()->SetTitleOffset(1.0);
		htmp->GetXaxis()->SetTitleOffset(1.05);
		htmp->GetXaxis()->SetTitleSize(0.055);
		htmp->GetXaxis()->SetLabelSize(0.043);
		htmp->GetYaxis()->SetLabelSize(0.043);
		htmp->GetYaxis()->SetNdivisions(510);

		TLegend *leg1 = new TLegend(0.5,0.63,0.9,0.93);
		SetLegendStyle(leg1);
		leg1->SetFillStyle(0);
		leg1->SetBorderSize(0);
		leg1->SetTextFont(43);
		leg1->SetTextSize(20);
		leg1->SetMargin(0.3);
		leg1->AddEntry("","p#font[122]{-}Pb #sqrt{#it{s}_{NN}} = 8 TeV","h");
		leg1->SetNColumns(2);
		leg1->AddEntry("","SHINCHON","");
		leg1->AddEntry("","CMS","");

		for(int i=0;i<nstat;i++){
			hRpPb_CMS_stat[i]->Draw("same");
			hRpPb_CMS_syst[i]->Draw("same,e2");

			gRpPbPt[0][i]->Draw("l3");

			leg1->AddEntry( gRpPbPt[0][i], Form("#varUpsilon(%dS)",i+1), "LF");
			leg1->AddEntry( hRpPb_CMS_stat[i], Form("#varUpsilon(%dS)",i+1), "P");
		}

		leg1->Draw();

		if ( bSAVE ){
			c5->cd();
			c5->SaveAs("NMFplots/RpPb_Pt_pPb_CMS.pdf");
		}

	}

 return;

 TCanvas *c;

 TLegend* legMultST = new TLegend(0.332, 0.749, 0.849, 0.949);
 legMultST->SetFillColorAlpha(0,0);
 legMultST->SetLineWidth(0.0);

 for(int i=0;i<nsyst;i++){
        legMultST->Clear();
        legMultST->AddEntry( (TObject*)0, Form("SHINCHON, %s #sqrt{#it{s}_{NN}} = 8.16 TeV",dname[i]), "");
        for(int j=0;j<nstat;j++){
                gRpPbMultST[i][j]->SetFillColor(j+1);
                gRpPbMultST[i][j]->SetMinimum(0.0);
                if( j==0 ) gRpPbMultST[i][j]->Draw("ae3");
                gRpPbMultST[i][j]->Draw("e3");
                legMultST->AddEntry( gRpPbMultST[i][j], Form("#varUpsilon(%dS)",j+1), "f");
        }
        legMultST->Draw();
	c->SaveAs(Form("NMFplots/R%s_multST.pdf",fname[i]));
 }

 TLegend* legMultST1 = new TLegend(0.311, 0.732, 0.834, 0.947);
 legMultST1->SetFillColorAlpha(0,0);
 legMultST1->SetLineWidth(0.0);

 for(int j=0;j<nstat;j++){
        legMultST1->Clear();
        legMultST1->AddEntry( (TObject*)0, Form("SHINCHON, #sqrt{#it{s}_{NN}} = 8.16 TeV, #varUpsilon(%dS)",j+1), "");
        for(int i=0;i<nsyst;i++){
                gRpPbMultST[i][j]->SetFillColor(i+1);
                if( i==0 ) gRpPbMultST[i][j]->Draw("ae3");
                gRpPbMultST[i][j]->Draw("e3");
                legMultST1->AddEntry( gRpPbMultST[i][j], Form("%s",dname[i]), "f");
        }
        legMultST1->Draw();
	c->SaveAs(Form("NMFplots/R_multST_%dS.pdf",j+1));
 }



 TLegend* legPt = new TLegend(0.156, 0.721, 0.680, 0.937);
 legPt->SetFillColorAlpha(0,0);
 legPt->SetLineWidth(0.0);

 TLegend* legPt1 = new TLegend(0.156, 0.721, 0.680, 0.937);
 legPt1->SetFillColorAlpha(0,0);
 legPt1->SetLineWidth(0.0);

 for(int j=0;j<nstat;j++){
        legPt1->Clear();
        legPt1->AddEntry( (TObject*)0, Form("SHINCHON, #sqrt{#it{s}_{NN}} = 8.16 TeV, #varUpsilon(%dS)",j+1), "");
        for(int i=0;i<nsyst;i++){
                gRpPbPt[i][j]->SetFillColor(i+1);
                if( i==0 ) gRpPbPt[i][j]->Draw("ae3");
                gRpPbPt[i][j]->Draw("e3");
                legPt1->AddEntry( gRpPbPt[i][j], Form("%s",dname[i]), "f");
        }
        legPt1->Draw();
	c->SaveAs(Form("NMFplots/R_Pt_%dS.pdf",j+1));
 }


 c->SaveAs("NMFplots/RpPb_Pt_compCMS.pdf");



}
