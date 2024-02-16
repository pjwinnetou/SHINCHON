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


void compareNMF_InitGeo_PbPb(){

 TCanvas* c = new TCanvas("c","c",800,700);
 gPad->SetLeftMargin(0.15);
 gPad->SetBottomMargin(0.15);
 gPad->SetRightMargin(0.03);
 gPad->SetTopMargin(0.03);
 gPad->SetTicks();

 TLegend* leg = new TLegend(0.245, 0.703, 0.943, 0.943);
 leg->SetTextFont(43);
 leg->SetTextSize(24);
 leg->SetFillStyle(0);
 leg->SetLineWidth(0);
 leg->SetNColumns(3);

 TFile* fdata_np = new TFile("data/ExpNMF_npart_Out.root","read");
 TGraphErrors* gCMS_1S_stat = (TGraphErrors*)fdata_np->Get("gCMS1stat");
 TGraphAsymmErrors* gCMS_1S_syst = (TGraphAsymmErrors*)fdata_np->Get("gCMS1syst");
 TGraphErrors* gCMS_2S_stat = (TGraphErrors*)fdata_np->Get("gCMS2stat");
 TGraphErrors* gCMS_2S_syst = (TGraphErrors*)fdata_np->Get("gCMS2syst");
 TGraphErrors* gCMS_3S_stat = (TGraphErrors*)fdata_np->Get("gCMS3stat");
 TGraphErrors* gCMS_3S_syst = (TGraphErrors*)fdata_np->Get("gCMS3syst");

 TFile* fdata_pt = new TFile("data/ExpNMF_pt_CMS_Out.root","read");
 TGraphErrors* gCMS_1S_pt_stat = (TGraphErrors*)fdata_pt->Get("gCMS_1S_pt_stat");
 TGraphErrors* gCMS_1S_pt_syst = (TGraphErrors*)fdata_pt->Get("gCMS_1S_pt_syst");
 TGraphErrors* gCMS_2S_pt_stat = (TGraphErrors*)fdata_pt->Get("gCMS_2S_pt_stat");
 TGraphErrors* gCMS_2S_pt_syst = (TGraphErrors*)fdata_pt->Get("gCMS_2S_pt_syst");
 TGraphErrors* gCMS_3S_pt_stat = (TGraphErrors*)fdata_pt->Get("gCMS_3S_pt_stat");
 TGraphErrors* gCMS_3S_pt_syst = (TGraphErrors*)fdata_pt->Get("gCMS_3S_pt_syst");

 TFile* ftp0 = new TFile("data/FDCOR_out_PbPb_trento_p0.root","read");
 TFile* ftp1 = new TFile("data/FDCOR_out_PbPb_trento_p1.root","read");

 const int nstat=3;
 int lc[nstat] = { 1, 30, 46};

 TGraphErrors* g_glb_pt[nstat];
 TGraphErrors* g_glb_npart[nstat];
 TGraphErrors* g_tp0_pt[nstat];
 TGraphErrors* g_tp0_npart[nstat];
 TGraphErrors* g_tp1_pt[nstat];
 TGraphErrors* g_tp1_npart[nstat];

 for(int i=0;i<nstat;i++){
	g_tp0_pt[i] = (TGraphErrors*)ftp0->Get(Form("g_RAA_ptdep_fdcor_%dS_trento_p0",i+1));
	g_tp0_npart[i] = (TGraphErrors*)ftp0->Get(Form("g_RAA_multdep_fdcor_%dS_trento_p0",i+1));

	g_tp1_pt[i] = (TGraphErrors*)ftp1->Get(Form("g_RAA_ptdep_fdcor_%dS_trento_p1",i+1));
	g_tp1_npart[i] = (TGraphErrors*)ftp1->Get(Form("g_RAA_multdep_fdcor_%dS_trento_p1",i+1));

	g_tp0_pt[i]->SetLineColor( lc[i] );
    g_tp0_npart[i]->SetLineColor( lc[i] );
	g_tp1_pt[i]->SetLineColor( lc[i] );
    g_tp1_npart[i]->SetLineColor( lc[i] );

    g_tp0_pt[i]->SetLineWidth(5);
    g_tp0_npart[i]->SetLineWidth(5);
    g_tp1_pt[i]->SetLineWidth(5);
    g_tp1_npart[i]->SetLineWidth(5);

	g_tp0_pt[i]->SetLineStyle(5);
	g_tp0_npart[i]->SetLineStyle(5);
 }

 g_tp1_npart[0]->GetXaxis()->SetTitleFont(43);
 g_tp1_npart[0]->GetXaxis()->SetLabelFont(43);
 g_tp1_npart[0]->GetYaxis()->SetTitleFont(43);
 g_tp1_npart[0]->GetYaxis()->SetLabelFont(43);

 g_tp1_npart[0]->GetXaxis()->SetTitleSize(32);
 g_tp1_npart[0]->GetXaxis()->SetLabelSize(28);
 g_tp1_npart[0]->GetYaxis()->SetTitleSize(32);
 g_tp1_npart[0]->GetYaxis()->SetLabelSize(28);

 g_tp1_npart[0]->SetTitle(";N_{part};Nuclear modification factor");
 g_tp1_npart[0]->GetXaxis()->SetNdivisions(505);
 g_tp1_npart[0]->GetYaxis()->SetNdivisions(505);

 g_tp1_npart[0]->SetMaximum(1.3);
 g_tp1_npart[0]->SetMinimum(0);
 g_tp1_npart[0]->Draw("ACEX0");
 for(int i=0;i<nstat;i++){
    g_tp1_npart[i]->Draw("CEX0");
	g_tp0_npart[i]->Draw("CEX0");
 }

 gCMS_1S_syst->SetFillColorAlpha(lc[0],0.3);
 gCMS_1S_syst->SetLineWidth(0);
 gCMS_1S_stat->SetMarkerStyle(29);
 gCMS_1S_stat->SetMarkerColor(lc[0]);
 gCMS_1S_stat->SetLineColor(lc[0]);
 gCMS_1S_stat->SetMarkerSize(2);

 gCMS_2S_syst->SetFillColorAlpha(lc[1],0.3);
 gCMS_2S_syst->SetLineWidth(0);
 gCMS_2S_stat->SetMarkerStyle(29);
 gCMS_2S_stat->SetMarkerColor(lc[1]);
 gCMS_2S_stat->SetLineColor(lc[1]);
 gCMS_2S_stat->SetMarkerSize(2);

 gCMS_1S_syst->Draw("2");
 gCMS_1S_stat->Draw("P");
 gCMS_2S_syst->Draw("2");
 gCMS_2S_stat->Draw("P");
 gCMS_3S_syst->Draw("2");
 gCMS_3S_stat->Draw("P");


 leg->SetHeader("Pb#font[122]{-}Pb, 5020 GeV");
 leg->AddEntry( g_tp1_npart[0], "#varUpsilon(1S)", "l");
 leg->AddEntry( g_tp1_npart[1], "#varUpsilon(2S)", "l");
 leg->AddEntry( g_tp1_npart[2], "#varUpsilon(3S)", "l");
 leg->AddEntry( g_tp0_npart[0], "Trento #it{p} = 0", "l");
 leg->AddEntry( g_tp1_npart[0], "Trento #it{p} = 1", "l");
 leg->AddEntry( gCMS_1S_stat, "CMS", "p");
 leg->Draw();

 c->SaveAs("figs/Comp_NMF_npart_trento_PbPb.pdf");


 g_tp1_pt[0]->GetXaxis()->SetTitleFont(43);
 g_tp1_pt[0]->GetXaxis()->SetLabelFont(43);
 g_tp1_pt[0]->GetYaxis()->SetTitleFont(43);
 g_tp1_pt[0]->GetYaxis()->SetLabelFont(43);

 g_tp1_pt[0]->GetXaxis()->SetTitleSize(32);
 g_tp1_pt[0]->GetXaxis()->SetLabelSize(28);
 g_tp1_pt[0]->GetYaxis()->SetTitleSize(32);
 g_tp1_pt[0]->GetYaxis()->SetLabelSize(28);

 g_tp1_pt[0]->SetTitle(";#it{p}_{T} (GeV/#it{c});Nuclear modification factor");
 g_tp1_pt[0]->GetXaxis()->SetNdivisions(505);
 g_tp1_pt[0]->GetYaxis()->SetNdivisions(505);

 g_tp1_pt[0]->SetMaximum(1.3);
 g_tp1_pt[0]->SetMinimum(0);
 g_tp1_pt[0]->Draw("ACEX0");
 for(int i=0;i<nstat;i++){
    g_tp0_pt[i]->Draw("CEX0");
    g_tp1_pt[i]->Draw("CEX0");
 }

 gCMS_1S_pt_syst->SetFillColorAlpha(lc[0],0.3);
 gCMS_1S_pt_syst->SetLineWidth(0);
 gCMS_1S_pt_stat->SetMarkerStyle(29);
 gCMS_1S_pt_stat->SetMarkerColor(lc[0]);
 gCMS_1S_pt_stat->SetLineColor(lc[0]);
 gCMS_1S_pt_stat->SetMarkerSize(2);

 gCMS_2S_pt_syst->SetFillColorAlpha(lc[1],0.3);
 gCMS_2S_pt_syst->SetLineWidth(0);
 gCMS_2S_pt_stat->SetMarkerStyle(29);
 gCMS_2S_pt_stat->SetMarkerColor(lc[1]);
 gCMS_2S_pt_stat->SetLineColor(lc[1]);
 gCMS_2S_pt_stat->SetMarkerSize(2);

 gCMS_3S_pt_syst->SetFillColorAlpha(lc[2],0.3);
 gCMS_3S_pt_syst->SetLineWidth(0);
 gCMS_3S_pt_stat->SetMarkerStyle(29);
 gCMS_3S_pt_stat->SetMarkerColor(lc[2]);
 gCMS_3S_pt_stat->SetLineColor(lc[2]);
 gCMS_3S_pt_stat->SetMarkerSize(2);

 gCMS_1S_pt_syst->Draw("2");
 gCMS_1S_pt_stat->Draw("P");
 gCMS_2S_pt_syst->Draw("2");
 gCMS_2S_pt_stat->Draw("P");
 gCMS_3S_pt_syst->Draw("2");
 gCMS_3S_pt_stat->Draw("P");

 leg->Clear();
 leg->SetHeader("Pb#font[122]{-}Pb, 5020 GeV");
 leg->AddEntry( g_tp1_npart[0], "#varUpsilon(1S)", "l");
 leg->AddEntry( g_tp1_npart[1], "#varUpsilon(2S)", "l");
 leg->AddEntry( g_tp1_npart[2], "#varUpsilon(3S)", "l");
 leg->AddEntry( g_tp0_npart[0], "Trento #it{p} = 0", "l");
 leg->AddEntry( g_tp1_npart[0], "Trento #it{p} = 1", "l");
 leg->AddEntry( gCMS_1S_pt_stat, "CMS", "p");
 leg->Draw();



 c->SaveAs("figs/Comp_NMF_npart_trento_PbPb_pt.pdf");

}
