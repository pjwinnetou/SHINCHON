TGraphErrors* convertHist( TH1* h){
 TGraphErrors* g = new TGraphErrors();
 int np=0;
 for(int i=0;i<h->GetNbinsX();i++){
        if( h->GetBinContent(i+1) == 0 ) continue;

        g->SetPoint( np, h->GetBinCenter(i+1), h->GetBinContent(i+1) );
        g->SetPointError( np, h->GetBinWidth(i+1)/2., h->GetBinError(i+1) );
        np++;
 }

 return g;
}


void RAA_default(){
// TFile* fpb = new TFile("data/FDCOR_out_PbPb.root","read");
 TFile* fpb = new TFile("data/FDCOR_out_PbPb_trento_p1.root","read");
 TFile* fau = new TFile("data/FDCOR_out_AuAu_default.root","read");


 TGraphErrors* gAuNpart[3];
 gAuNpart[0] = (TGraphErrors*)fau->Get("g_RpPb_multdep_fdcor_1S_default");
 gAuNpart[1] = (TGraphErrors*)fau->Get("g_RpPb_multdep_fdcor_2S_default");
 gAuNpart[2] = (TGraphErrors*)fau->Get("g_RpPb_multdep_fdcor_3S_default");


 TProfile* hPbNpart[3];
// hPbNpart[0] = (TProfile*)fpb->Get("hRpA_pPb_mult_corr_1S");
// hPbNpart[1] = (TProfile*)fpb->Get("hRpA_pPb_mult_corr_2S");
// hPbNpart[2] = (TProfile*)fpb->Get("hRpA_pPb_mult_corr_3S");

 TGraphErrors* gPbNpart[3];
// gPbNpart[0] = (TGraphErrors*)convertHist(hPbNpart[0]);
// gPbNpart[1] = (TGraphErrors*)convertHist(hPbNpart[1]);
// gPbNpart[2] = (TGraphErrors*)convertHist(hPbNpart[2]);
 gPbNpart[0] = (TGraphErrors*)fpb->Get("g_RAA_multdep_fdcor_1S_trento_p1");
 gPbNpart[1] = (TGraphErrors*)fpb->Get("g_RAA_multdep_fdcor_2S_trento_p1");
 gPbNpart[2] = (TGraphErrors*)fpb->Get("g_RAA_multdep_fdcor_3S_trento_p1");
 TCanvas* c = new TCanvas("c","c",800,700);
 gPad->SetLeftMargin(0.15);
 gPad->SetBottomMargin(0.15);
 gPad->SetRightMargin(0.03);
 gPad->SetTopMargin(0.03);
 gPad->SetTicks();

 gPbNpart[0]->SetLineStyle(5);
 gPbNpart[0]->SetLineWidth(5);

 gPbNpart[0]->GetXaxis()->SetTitleFont(43);
 gPbNpart[0]->GetXaxis()->SetLabelFont(43);
 gPbNpart[0]->GetYaxis()->SetTitleFont(43);
 gPbNpart[0]->GetYaxis()->SetLabelFont(43);

 gPbNpart[0]->GetXaxis()->SetTitleSize(32);
 gPbNpart[0]->GetXaxis()->SetLabelSize(28);
 gPbNpart[0]->GetYaxis()->SetTitleSize(32);
 gPbNpart[0]->GetYaxis()->SetLabelSize(28);

 gPbNpart[0]->SetTitle(";N_{part};Nuclear modification factor");
 gPbNpart[0]->SetMaximum(1.3);
 gPbNpart[0]->SetMinimum(0.00);
 gPbNpart[0]->GetXaxis()->SetNdivisions(505);
 gPbNpart[0]->GetYaxis()->SetNdivisions(505);

 gPbNpart[0]->Draw("ACex0");

 gAuNpart[0]->SetLineWidth(5);
 gAuNpart[0]->Draw("Cex0");

 int lc[3] = { 1, 30, 46};

 for(int i=0;i<3;i++){
	gPbNpart[i]->SetLineWidth(5);
	gAuNpart[i]->SetLineWidth(5);

	gPbNpart[i]->SetLineColor( lc[i] );
	gAuNpart[i]->SetLineColor( lc[i] );

	gPbNpart[i]->SetLineStyle(5);

	gPbNpart[i]->Draw("Cex0");
	gAuNpart[i]->Draw("Cex0");
 } 

//
 TFile *infile1 = new TFile("data/RAA_Y1S_Npart_PbPb5TeV_HEPData-ins1674529-v2-Table_19.root","read");
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

 TFile *infile2 = new TFile("data/RAA_HIN21007_prelim.root","read");
 TGraphErrors *gCMS2 = (TGraphErrors*)infile2->Get("RAA_2S");
 TGraphErrors *gCMS3 = (TGraphErrors*)infile2->Get("RAA_3S");
 TGraphErrors *gCMS2sys = (TGraphErrors*)infile2->Get("RAA_2S_sys");
 TGraphErrors *gCMS3sys = (TGraphErrors*)infile2->Get("RAA_3S_sys");

 gCMSsys->SetFillColorAlpha(lc[0],0.3);
 gCMSsys->SetLineWidth(0);
 gCMS->SetMarkerStyle(20);
 gCMS->SetMarkerColor(lc[0]);
 gCMS->SetLineColor(lc[0]);
 gCMS->SetMarkerSize(2);

 gCMS2sys->SetFillColorAlpha(lc[1],0.3);
 gCMS2sys->SetLineWidth(0);
 gCMS2->SetMarkerStyle(20);
 gCMS2->SetMarkerColor(lc[1]);
 gCMS2->SetLineColor(lc[1]);
 gCMS2->SetMarkerSize(2);

 gCMS3sys->SetFillColorAlpha(lc[2],0.3);
 gCMS3sys->SetLineWidth(0);
 gCMS3->SetMarkerStyle(20);
 gCMS3->SetMarkerColor(lc[2]);
 gCMS3->SetLineColor(lc[2]);        
 gCMS3->SetMarkerSize(2);

 gCMSsys->Draw("2");
 gCMS->Draw("P");

 gCMS2sys->Draw("2");
 gCMS2->Draw("P");

 gCMS3sys->Draw("2");
 gCMS3->Draw("P");


 TGraphErrors* gSTAR_1S_stat = new TGraphErrors();
 TGraphErrors* gSTAR_1S_syst = new TGraphErrors();
 TGraphErrors* gSTAR_2S_stat = new TGraphErrors();
 TGraphErrors* gSTAR_2S_syst = new TGraphErrors();
 TGraphErrors* gSTAR_3S_stat = new TGraphErrors();
 TGraphErrors* gSTAR_3S_syst = new TGraphErrors();

 gSTAR_1S_stat->SetPoint( 0, 81, 0.56);
 gSTAR_1S_stat->SetPoint( 1, 203, 0.36);
 gSTAR_1S_stat->SetPoint( 2, 325, 0.40);
 gSTAR_1S_stat->SetPointError( 0, 0, 0.07);
 gSTAR_1S_stat->SetPointError( 1, 0, 0.05);
 gSTAR_1S_stat->SetPointError( 2, 0, 0.06);

 gSTAR_1S_syst->SetPoint( 0, 81, 0.56);
 gSTAR_1S_syst->SetPoint( 1, 203, 0.36);
 gSTAR_1S_syst->SetPoint( 2, 325, 0.40);
 gSTAR_1S_syst->SetPointError( 0, 5, 0.10);
 gSTAR_1S_syst->SetPointError( 1, 5, 0.03);
 gSTAR_1S_syst->SetPointError( 2, 5, 0.05);

 gSTAR_2S_stat->SetPoint( 0, 81, 0.48);
 gSTAR_2S_stat->SetPoint( 1, 203, 0.26);
 gSTAR_2S_stat->SetPoint( 2, 325, 0.09);
 gSTAR_2S_stat->SetPointError( 0, 0, 0.16);
 gSTAR_2S_stat->SetPointError( 1, 0, 0.11);
 gSTAR_2S_stat->SetPointError( 2, 0, 0.14);

 gSTAR_2S_syst->SetPoint( 0, 81, 0.48);
 gSTAR_2S_syst->SetPoint( 1, 203, 0.26);
 gSTAR_2S_syst->SetPoint( 2, 325, 0.09);
 gSTAR_2S_syst->SetPointError( 0, 5, 0.09);
 gSTAR_2S_syst->SetPointError( 1, 5, 0.03);
 gSTAR_2S_syst->SetPointError( 2, 5, 0.02);

 
 gSTAR_1S_syst->SetFillColorAlpha(lc[0],0.3);
 gSTAR_1S_syst->SetLineWidth(0);
 gSTAR_1S_stat->SetMarkerStyle(29);
 gSTAR_1S_stat->SetMarkerColor(lc[0]);
 gSTAR_1S_stat->SetLineColor(lc[0]);
 gSTAR_1S_stat->SetMarkerSize(2);

 gSTAR_2S_syst->SetFillColorAlpha(lc[1],0.3);
 gSTAR_2S_syst->SetLineWidth(0);
 gSTAR_2S_stat->SetMarkerStyle(29);
 gSTAR_2S_stat->SetMarkerColor(lc[1]);
 gSTAR_2S_stat->SetLineColor(lc[1]);
 gSTAR_2S_stat->SetMarkerSize(2);

 gSTAR_1S_syst->Draw("2");
 gSTAR_1S_stat->Draw("P");
 gSTAR_2S_syst->Draw("2");
 gSTAR_2S_stat->Draw("P");

// data

 for(int i=0;i<3;i++){
    gPbNpart[i]->Draw("Cex0");
    gAuNpart[i]->Draw("Cex0");
 }

 TLegend* leg = new TLegend(0.245, 0.703, 0.943, 0.943);
 leg->SetTextFont(43);
 leg->SetTextSize(24);
 leg->SetFillStyle(0);
 leg->SetLineWidth(0);
 leg->SetNColumns(4);
 leg->AddEntry( (TObject*)0, "Au#font[122]{-}Au, 200 GeV", "");
 leg->AddEntry( (TObject*)0, "", "");
 leg->AddEntry( (TObject*)0, "Pb#font[122]{-}Pb, 5.02 TeV", "");
 leg->AddEntry( (TObject*)0, "", "");

 leg->AddEntry( (TObject*)0, "SHINCHON", "");
 leg->AddEntry( (TObject*)0, "STAR", "");
 leg->AddEntry( (TObject*)0, "SHINCHON", "");
 leg->AddEntry( (TObject*)0, "CMS", "");

 leg->AddEntry( gAuNpart[0], "#varUpsilon(1S)", "l");
 leg->AddEntry( gSTAR_1S_stat, "#varUpsilon(1S)", "p");
 leg->AddEntry( gPbNpart[0], "#varUpsilon(1S)", "l");
 leg->AddEntry( gCMS, "#varUpsilon(1S)", "p");

 leg->AddEntry( gAuNpart[1], "#varUpsilon(2S)", "l");
 leg->AddEntry( gSTAR_2S_stat, "#varUpsilon(2S)", "p");
 leg->AddEntry( gPbNpart[1], "#varUpsilon(2S)", "l");
 leg->AddEntry( gCMS2, "#varUpsilon(2S)", "p");

 leg->AddEntry( gAuNpart[2], "#varUpsilon(3S)", "l");
 leg->AddEntry( (TObject*)0, "", "");
 leg->AddEntry( gPbNpart[2], "#varUpsilon(3S)", "l");
 leg->AddEntry( gCMS3, "#varUpsilon(3S)", "p");

 leg->Draw();

 c->SaveAs("figs/syst_comp.pdf");
}
