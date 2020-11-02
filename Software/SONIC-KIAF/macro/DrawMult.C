void DrawMult(){


	float tmp_xx, tmp_yy, tmp_yyerr;

	TGraphErrors *gALICE_5020GeV = new TGraphErrors;
	TGraphErrors *gALICE_2760GeV = new TGraphErrors;
	TGraphErrors *gPHENIX_200GeV = new TGraphErrors;

	ifstream fdata;

	fdata.open("ALICE_PbPb_5020GeV.dat");
	int ndata = 0;
	while ( fdata >> tmp_xx >> tmp_yy >> tmp_yyerr ){
		gALICE_5020GeV->SetPoint(ndata, tmp_xx, tmp_yy);
		gALICE_5020GeV->SetPointError(ndata, 0, tmp_yyerr);
		ndata++;
	}
	fdata.close();

	fdata.open("ALICE_PbPb_2760GeV.dat");
	ndata = 0;
	while ( fdata >> tmp_xx >> tmp_yy >> tmp_yyerr ){
		gALICE_2760GeV->SetPoint(ndata, tmp_xx, tmp_yy);
		gALICE_2760GeV->SetPointError(ndata, 0, tmp_yyerr);
		ndata++;
	}
	fdata.close();

	float PHENIX_npart[12] = {
		350.9, 297.9, 251.0, 211.0,
		176.3, 146.8, 120.9, 98.3,
		78.7, 61.9, 47.6, 35.6};
	float PHENIX_dNchdeta[12] = {
		3.92, 3.77, 3.64, 3.52,
		3.43, 3.35, 3.26, 3.18, 
		3.14, 3.08, 2.98, 2.93};
	float PHENIX_dNchdeta_err[12] = {
		0.22, 0.21, 0.21, 0.21,
		0.22, 0.25, 0.28, 0.31,
		0.34, 0.38, 0.44, 0.56};

	for (int ii=0; ii<12; ii++){
		gPHENIX_200GeV->SetPoint(ii, PHENIX_npart[ii], PHENIX_dNchdeta[ii]);
		gPHENIX_200GeV->SetPointError(ii, 0, PHENIX_dNchdeta_err[ii]);
	}

	TFile *infile0 = new TFile("SONIC-ana-PbPb5TeV.root","read");
	TProfile *hprof0 = (TProfile*)infile0->Get("hprof_tuned");
	hprof0->Rebin();

	TGraphErrors *gSONIC0 = new TGraphErrors;
	for (int ii=0; ii<hprof0->GetNbinsX(); ii++){
		float xx = hprof0->GetBinCenter(ii+1);
		float yy = hprof0->GetBinContent(ii+1);
		float yyerr = 0.1;
		if ( yy<0.1 ) continue;
		gSONIC0->SetPoint(ii, xx, yy);
		gSONIC0->SetPointError(ii, 0, yyerr);
	}

	TFile *infile1 = new TFile("SONIC-ana-PbPb2760GeV.root","read");
	TProfile *hprof1 = (TProfile*)infile1->Get("hprof_tuned");
	hprof1->Rebin();

	TGraphErrors *gSONIC1 = new TGraphErrors;
	for (int ii=0; ii<hprof1->GetNbinsX(); ii++){
		float xx = hprof1->GetBinCenter(ii+1);
		float yy = hprof1->GetBinContent(ii+1);
		float yyerr = 0.1;
		if ( yy<0.1 ) continue;
		gSONIC1->SetPoint(ii, xx, yy);
		gSONIC1->SetPointError(ii, 0, yyerr);
	}

	TFile *infile2 = new TFile("SONIC-ana-AuAu200GeV.root","read");
	TProfile *hprof2 = (TProfile*)infile2->Get("hprof_tuned");
	hprof2->Rebin();

	TGraphErrors *gSONIC2 = new TGraphErrors;
	for (int ii=0; ii<hprof2->GetNbinsX(); ii++){
		float xx = hprof2->GetBinCenter(ii+1);
		float yy = hprof2->GetBinContent(ii+1);
		float yyerr = 0.1;
		if ( yy<0.1 ) continue;
		gSONIC2->SetPoint(ii, xx, yy);
		gSONIC2->SetPointError(ii, 0, yyerr);
	}


	{
		TCanvas *c1 = new TCanvas("c1","c1",1.2*500,500);
		gPad->SetMargin(0.14,0.05,0.12,0.03);
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,420,15); 
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitle("#LTN_{part}#GT");
		htmp->GetYaxis()->SetTitle("#frac{2}{#LTN_{part}#GT} #LTdN_{ch}/d#eta#GT");

		gSONIC0->SetLineColorAlpha(4,0.3);
		gSONIC0->SetFillColorAlpha(4,0.3);
		gSONIC0->Draw("3");

		gSONIC1->SetLineColorAlpha(2,0.3);
		gSONIC1->SetFillColorAlpha(2,0.3);
		gSONIC1->Draw("3");

		gSONIC2->SetLineColorAlpha(1,0.3);
		gSONIC2->SetFillColorAlpha(1,0.3);
		gSONIC2->Draw("3");

		gALICE_5020GeV->SetMarkerStyle(20);
		gALICE_5020GeV->SetMarkerColor(4);
		gALICE_5020GeV->SetLineColor(4);
		gALICE_5020GeV->SetLineWidth(2);
		gALICE_5020GeV->Draw("p");

		gALICE_2760GeV->SetMarkerStyle(21);
		gALICE_2760GeV->SetMarkerColor(2);
		gALICE_2760GeV->SetLineColor(2);
		gALICE_2760GeV->SetLineWidth(2);
		gALICE_2760GeV->Draw("p");

		gPHENIX_200GeV->SetMarkerStyle(24);
		gPHENIX_200GeV->SetMarkerColor(1);
		gPHENIX_200GeV->SetLineColor(1);
		gPHENIX_200GeV->SetLineWidth(2);
		gPHENIX_200GeV->Draw("p");


		TLegend *leg = new TLegend(0.15,0.75,0.7,0.9);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.04);
		leg->AddEntry(gALICE_5020GeV,"ALICE, Pb+Pb 5.02 TeV","P");
		leg->AddEntry(gALICE_2760GeV,"ALICE, Pb+Pb 2.76 TeV","P");
		leg->AddEntry(gPHENIX_200GeV,"PHENIX, Au+Au 200 GeV","P");
		leg->Draw();
	}


}
