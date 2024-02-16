{
 TCanvas* c = new TCanvas("c","c",800,700);
 gPad->SetLeftMargin(0.15);
 gPad->SetBottomMargin(0.15);
 gPad->SetRightMargin(0.03);
 gPad->SetTopMargin(0.03);
 gPad->SetTicks();
 gStyle->SetOptStat(0);

 TLegend* leg = new TLegend(0.245, 0.703, 0.943, 0.943);
 leg->SetTextFont(43);
 leg->SetTextSize(24);
 leg->SetFillStyle(0);
 leg->SetLineWidth(0);
 leg->SetNColumns(3);

 TFile* fdata_pt = new TFile("data/ExpNMF_pt_Out.root","read");
 TGraphErrors* gSTAR_1S_pt_stat = (TGraphErrors*)fdata_pt->Get("gSTAR_1S_pt_stat");
 TGraphErrors* gSTAR_1S_pt_syst = (TGraphErrors*)fdata_pt->Get("gSTAR_1S_pt_syst");
 TGraphErrors* gSTAR_2S_pt_stat = (TGraphErrors*)fdata_pt->Get("gSTAR_2S_pt_stat");
 TGraphErrors* gSTAR_2S_pt_syst = (TGraphErrors*)fdata_pt->Get("gSTAR_2S_pt_syst");

 TFile* fglb = new TFile("data/FDCOR_out_AuAu_default.root","read");
 TFile* ftp0 = new TFile("data/FDCOR_out_AuAu_trento_p0.root","read");
 TFile* ftp1 = new TFile("data/FDCOR_out_AuAu_trento_p1.root","read");

 const int nstat=3;
 int lc[nstat] = { 1, 30, 46};

 TGraphErrors* g_glb_pt[nstat];
 TGraphErrors* g_tp0_pt[nstat];
 TGraphErrors* g_tp1_pt[nstat];

 TProfile* h_glb_pt[nstat];
 TProfile* h_tp0_pt[nstat];
 TProfile* h_tp1_pt[nstat];

 for(int i=0;i<nstat;i++){
    g_glb_pt[i] = (TGraphErrors*)fglb->Get(Form("g_RpPb_ptdep_fdcor_%dS_default",i+1));
    g_tp0_pt[i] = (TGraphErrors*)ftp0->Get(Form("g_RAA_ptdep_fdcor_%dS_trento_p0",i+1));
    g_tp1_pt[i] = (TGraphErrors*)ftp1->Get(Form("g_RAA_ptdep_fdcor_%dS_trento_p1",i+1));

	h_glb_pt[i] = new TProfile(Form("h_glb_pt_%d",i),"",10,0,20);
	h_tp0_pt[i] = new TProfile(Form("h_tp0_pt_%d",i),"",10,0,20);
	h_tp1_pt[i] = new TProfile(Form("h_tp1_pt_%d",i),"",10,0,20);

	for(int j=0;j<g_glb_pt[i]->GetN();j++){
		h_glb_pt[i]->Fill( g_glb_pt[i]->GetX()[j], g_glb_pt[i]->GetY()[j]);
	}
	for(int j=0;j<g_tp0_pt[i]->GetN();j++){
		h_tp0_pt[i]->Fill( g_tp0_pt[i]->GetX()[j], g_tp0_pt[i]->GetY()[j]);
	}
	for(int j=0;j<g_tp1_pt[i]->GetN();j++){
        h_tp1_pt[i]->Fill( g_tp1_pt[i]->GetX()[j], g_tp1_pt[i]->GetY()[j]);
    }

    h_glb_pt[i]->SetLineColor( lc[i] );
    h_tp0_pt[i]->SetLineColor( lc[i] );
    h_tp1_pt[i]->SetLineColor( lc[i] );

    h_glb_pt[i]->SetLineWidth(5);
    h_tp0_pt[i]->SetLineWidth(5);
    h_tp1_pt[i]->SetLineWidth(5);

    h_tp0_pt[i]->SetLineStyle(5);
    h_tp1_pt[i]->SetLineStyle(7);
 }

 h_glb_pt[0]->GetXaxis()->SetTitleFont(43);
 h_glb_pt[0]->GetXaxis()->SetLabelFont(43);
 h_glb_pt[0]->GetYaxis()->SetTitleFont(43);
 h_glb_pt[0]->GetYaxis()->SetLabelFont(43);

 h_glb_pt[0]->GetXaxis()->SetTitleSize(32);
 h_glb_pt[0]->GetXaxis()->SetLabelSize(28);
 h_glb_pt[0]->GetYaxis()->SetTitleSize(32);
 h_glb_pt[0]->GetYaxis()->SetLabelSize(28);

 h_glb_pt[0]->SetTitle(";#it{p}_{T} (GeV/#it{c});Nuclear modification factor");
 h_glb_pt[0]->GetXaxis()->SetNdivisions(505);
 h_glb_pt[0]->GetYaxis()->SetNdivisions(505);

 h_glb_pt[0]->SetMaximum(1.3);
 h_glb_pt[0]->SetMinimum(0);
 h_glb_pt[0]->Draw("C");
 for(int i=0;i<nstat;i++){
    h_glb_pt[i]->Draw("C,same");
    h_tp0_pt[i]->Draw("C,same");
    h_tp1_pt[i]->Draw("C,same");
 }

 gSTAR_1S_pt_syst->SetFillColorAlpha(lc[0],0.3);
 gSTAR_1S_pt_syst->SetLineWidth(0);
 gSTAR_1S_pt_stat->SetMarkerStyle(29);
 gSTAR_1S_pt_stat->SetMarkerColor(lc[0]);
 gSTAR_1S_pt_stat->SetLineColor(lc[0]);
 gSTAR_1S_pt_stat->SetMarkerSize(2);

 gSTAR_2S_pt_syst->SetFillColorAlpha(lc[1],0.3);
 gSTAR_2S_pt_syst->SetLineWidth(0);
 gSTAR_2S_pt_stat->SetMarkerStyle(29);
 gSTAR_2S_pt_stat->SetMarkerColor(lc[1]);
 gSTAR_2S_pt_stat->SetLineColor(lc[1]);
 gSTAR_2S_pt_stat->SetMarkerSize(2);

 gSTAR_1S_pt_syst->Draw("2");
 gSTAR_1S_pt_stat->Draw("P");
 gSTAR_2S_pt_syst->Draw("2");
 gSTAR_2S_pt_stat->Draw("P");

 leg->SetHeader("Au#font[122]{-}Au, 200 GeV");
 leg->AddEntry( h_glb_pt[0], "#varUpsilon(1S)", "l");
 leg->AddEntry( h_glb_pt[1], "#varUpsilon(2S)", "l");
 leg->AddEntry( h_glb_pt[2], "#varUpsilon(3S)", "l");
 leg->AddEntry( h_glb_pt[0], "Glauber", "l");
 leg->AddEntry( h_tp0_pt[0], "Trento #it{p} = 0", "l");
 leg->AddEntry( h_tp1_pt[0], "Trento #it{p} = 1", "l");
 leg->AddEntry( gSTAR_1S_pt_stat, "STAR", "p");


 leg->Draw();
 c->SaveAs("figs/Comp_NMF_pt_trento.pdf");

}
