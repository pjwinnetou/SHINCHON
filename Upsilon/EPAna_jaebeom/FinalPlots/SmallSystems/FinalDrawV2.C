#include <iostream>
//#include "../../Style_jaebeom.h"
#include "../../Style_jaebeom_woFrameLegend.h"
#include "../../SHINCHONLegend_raaCent.C"
#include "../../tdrstyle.C"

using namespace std;

void FinalDrawV2(){

	const bool bSAVE = true;

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

	int nColor[nstat] = {1, 4, kGreen+2}; 

	TFile *infile[nsyst][nstat];
	TGraphErrors *gv2[nsyst][nstat];
	TGraphErrors *gv2_err[nsyst][nstat];
	
	for (int isys=0; isys<nsyst; isys++){
		for (int ista=0; ista<nstat; ista++){

			infile[isys][ista] = new TFile(Form("v2input/v2_vs_pt_%s_isLine1_%ds_FDall.root",fname[isys],ista+1),"read");

			gv2[isys][ista] = (TGraphErrors*) infile[isys][ista]->Get("gv2");
			gv2_err[isys][ista] = (TGraphErrors*) infile[isys][ista]->Get("gv2_error");

			gv2_err[isys][ista]->SetLineColorAlpha(nColor[isys], 0.3);
			gv2_err[isys][ista]->SetFillColorAlpha(nColor[isys], 0.3);
			gv2_err[isys][ista]->SetLineWidth(4);
			gv2_err[isys][ista]->SetLineStyle(1);

		}//ista
	}//isys

	double ymin = -0.01; double ymax = 0.03; double xmin = 0; double xmax = 20;

	TCanvas *c2[nstat];

	for(int i=0;i<nstat;i++){
		c2[i] = new TCanvas(Form("c2_%d",i),Form("c2_%d",i),1.2*500,500);
		gPad->SetLeftMargin(0.14);

		TH1D *htmp = (TH1D*)gPad->DrawFrame(xmin, ymin, xmax, ymax);
		htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		htmp->GetYaxis()->SetTitle("v_{2}");
    htmp->GetYaxis()->SetTitleOffset(1.1);
    htmp->GetXaxis()->SetTitleOffset(1.05);
    htmp->GetXaxis()->SetTitleSize(0.055);
    htmp->GetXaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetNdivisions(510);

		dashedLine(xmin,0.,xmax,0.,1,1);

		gv2_err[0][i]->Draw("L3");
		gv2_err[1][i]->Draw("L3");
		gv2_err[2][i]->Draw("L3");

		TLegend *leg = new TLegend(0.45,0.93-0.06*4,0.85,0.93);
		SetLegendStyle(leg);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(43);
		leg->SetTextSize(20);
		leg->SetMargin(0.3);
		leg->AddEntry( (TObject*)0, Form("SHINCHON, #sqrt{#it{s}_{NN}} = 8 TeV, #varUpsilon(%dS)",i+1), "h");
		leg->AddEntry(gv2_err[0][i], Form("%s",dname[0]), "LF");
		leg->AddEntry(gv2_err[1][i], Form("%s",dname[1]), "LF");
		leg->AddEntry(gv2_err[2][i], Form("%s",dname[2]), "LF");
		leg->Draw();

		if ( bSAVE ){
			c2[i]->cd();
			c2[i]->SaveAs(Form("NMFplots/V2_Pt_Y%dS.pdf",i+1));
		}

	}//

	return;

	/*
	int iPeriod = -1;
	int iPos = 33;
	int Ystate = 1;
	bool drawInner = true;

	TMultiGraph *mg = new TMultiGraph();
	mg->SetTitle("Exclusion graphs");
	TMultiGraph *mgerror = new TMultiGraph();
	mgerror->SetTitle("Exclusion graphs");

	TGraphErrors* gv2_pO = (TGraphErrors*) fV2respO->Get("gv2");
	TGraphErrors* gv2_pO_error = (TGraphErrors*) fV2respO->Get("gv2_error");
	TGraphErrors* gv2_OO = (TGraphErrors*) fV2resOO->Get("gv2");
	TGraphErrors* gv2_OO_error = (TGraphErrors*) fV2resOO->Get("gv2_error");

	TCanvas *c1= new TCanvas("c1","c1",600,600);
	SetCanvasSquare2(c1);

	SetGraphStyleSys(gv2_pPb,1); SetGraphStyleSys(gv2_pPb_error,1); gv2_pPb_error->SetFillColorAlpha(kBlue-3,0.4);
	SetGraphStyleSys(gv2_pO,4);  SetGraphStyleSys(gv2_pO_error,4); gv2_pO_error->SetFillColorAlpha(kGray+2,0.4);
	SetGraphStyleSys(gv2_OO,0); SetGraphStyleSys(gv2_OO_error,0); gv2_OO_error->SetFillColorAlpha(kRed-4,0.4);

	mg->Add(gv2_pPb); mgerror->Add(gv2_pPb_error);
	mg->Add(gv2_pO); mgerror->Add(gv2_pO_error);
	mg->Add(gv2_OO); mgerror->Add(gv2_OO_error);


	c1->cd();
	mgerror->GetYaxis()->SetTitleOffset(1.4);
	mgerror->GetXaxis()->SetTitleOffset(1.1);
	mgerror->GetXaxis()->SetTitleSize(0.055);
	mgerror->GetXaxis()->SetLabelSize(0.043);
	mgerror->GetYaxis()->SetLabelSize(0.043);
	mgerror->GetYaxis()->SetNdivisions(510);
	mgerror->GetXaxis()->SetLimits(xmin,xmax);
	mgerror->GetXaxis()->SetRangeUser(xmin,xmax);
	mgerror->SetMinimum(ymin);
	mgerror->SetMaximum(ymax);
	mgerror->GetXaxis()->SetTitle("p_{T}^{#varUpsilon} (GeV/c)");
	mgerror->GetYaxis()->SetTitle("#it{v_{2}}");
	mg->Draw("AL"); mgerror->Draw("LE3 SAME"); 

	//SetGraphAxis(mg,"p_{T}^{#varUpsilon} (GeV/c)","#it{v_{2}}");
	//gv2_pPb->Draw("SAME");



	double legposx1 = 0.275;
	double legposx2 = 0.55;
	double legposy1 = 0.51;
	double legposy2 = 0.71;
	double labtextsize=0.04;
	TLegend* leg= new TLegend(legposx1,legposy1,legposx2,legposy2);
	SetLegendStyle(leg);
	leg->SetTextSize(labtextsize);
	leg->AddEntry(gv2_pO,"p+O","f");
	leg->AddEntry(gv2_OO,"O+O","f");
	leg->AddEntry(gv2_pPb,"p+Pb","f");
	leg->Draw("same");

	double lab_posx = 0.275; double lab_posy = 0.73; double lab_pos_diff = 0.269;
	drawGlobText(Form("Y(%iS), #sqrt{s_{NN}} = 8.16 TeV",Ystate), lab_posx, lab_posy, 1, labtextsize);
	SHINCHONLegend(c1,iPeriod,iPos,drawInner);
	*/


	return;

}
