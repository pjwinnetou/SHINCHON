{
 TCanvas* c = new TCanvas("c","c",800,700);
 gPad->SetLeftMargin(0.15);
 gPad->SetBottomMargin(0.15);
 gPad->SetRightMargin(0.03);
 gPad->SetTopMargin(0.03);
 gPad->SetTicks();

 TFile* fp0_AuAu = new TFile("data/FDCOR_out_AuAu_trento_p0.root","read");
 TFile* fp1_AuAu = new TFile("data/FDCOR_out_AuAu_trento_p1.root","read");
 TFile* fp0_PbPb = new TFile("data/FDCOR_out_PbPb_trento_p0.root","read");
 TFile* fp1_PbPb = new TFile("data/FDCOR_out_PbPb_trento_p1.root","read");

 const int nstat = 3;
 TGraphErrors* gp0_AuAu_np[nstat];
 TGraphErrors* gp1_AuAu_np[nstat];
 TGraphErrors* gp0_PbPb_np[nstat];
 TGraphErrors* gp1_PbPb_np[nstat];

 TGraphErrors* gRatio_AuAu[nstat];
 TGraphErrors* gRatio_PbPb[nstat];
 
 TLegend* leg = new TLegend(0.484, 0.64, 0.984, 0.940);
 leg->SetTextFont(43);
 leg->SetTextSize(32);
 leg->SetLineWidth(0);
 leg->SetFillStyle(0);


 for(int i=0;i<nstat;i++){
	gp0_AuAu_np[i] = (TGraphErrors*)fp0_AuAu->Get(Form("g_RAA_multdep_fdcor_%dS_trento_p0",i+1));
	gp1_AuAu_np[i] = (TGraphErrors*)fp1_AuAu->Get(Form("g_RAA_multdep_fdcor_%dS_trento_p1",i+1));
	gp0_AuAu_np[i]->SetName(Form("g_RAA_multdep_fdcor_%dS_AuAu_trento_p0",i+1));
	gp1_AuAu_np[i]->SetName(Form("g_RAA_multdep_fdcor_%dS_AuAu_trento_p1",i+1));

	gp0_PbPb_np[i] = (TGraphErrors*)fp0_PbPb->Get(Form("g_RAA_multdep_fdcor_%dS_trento_p0",i+1));
	gp1_PbPb_np[i] = (TGraphErrors*)fp1_PbPb->Get(Form("g_RAA_multdep_fdcor_%dS_trento_p1",i+1));
	gp0_PbPb_np[i]->SetName(Form("g_RAA_multdep_fdcor_%dS_PbPb_trento_p0",i+1));
	gp1_PbPb_np[i]->SetName(Form("g_RAA_multdep_fdcor_%dS_PbPb_trento_p1",i+1));

	gRatio_AuAu[i] = new TGraphErrors();
	gRatio_PbPb[i] = new TGraphErrors();
	for(int j=0;j<gp0_AuAu_np[i]->GetN();j++){
		gRatio_AuAu[i]->SetPoint( j, gp0_AuAu_np[i]->GetX()[j], gp0_AuAu_np[i]->GetY()[j] / gp1_AuAu_np[i]->GetY()[j] );
	}
	for(int j=0;j<gp0_PbPb_np[i]->GetN();j++){
        gRatio_PbPb[i]->SetPoint( j, gp0_PbPb_np[i]->GetX()[j], gp0_PbPb_np[i]->GetY()[j] / gp1_PbPb_np[i]->GetY()[j] );
    }

	gRatio_AuAu[i]->GetXaxis()->SetTitleFont(43);
	gRatio_AuAu[i]->GetYaxis()->SetTitleFont(43);
	gRatio_AuAu[i]->GetXaxis()->SetLabelFont(43);
	gRatio_AuAu[i]->GetYaxis()->SetLabelFont(43);

	gRatio_AuAu[i]->GetXaxis()->SetTitleSize(32);
	gRatio_AuAu[i]->GetYaxis()->SetTitleSize(32);
	gRatio_AuAu[i]->GetXaxis()->SetLabelSize(28);
	gRatio_AuAu[i]->GetYaxis()->SetLabelSize(28);

	gRatio_AuAu[i]->SetTitle(";#it{N}_{part};R_{AA}^{p=0}/R_{AA}^{p=1}");
	gRatio_AuAu[i]->SetMarkerStyle(20);
	gRatio_AuAu[i]->SetMinimum(0.0);
	gRatio_AuAu[i]->SetMaximum(2.0);

	gRatio_PbPb[i]->SetMarkerStyle(22);
	gRatio_PbPb[i]->SetMarkerColor(2);
	gRatio_PbPb[i]->SetLineColor(2);
	
	gRatio_AuAu[i]->Draw("AP");
	gRatio_PbPb[i]->Draw("P");

	leg->Clear();
	leg->AddEntry( (TObject*)0, Form("SHINCHON, #varUpsilon(%dS)",i+1), "");	
	leg->AddEntry( gRatio_AuAu[i], "Au#font[122]{-}Au 200 GeV", "p");
	leg->AddEntry( gRatio_PbPb[i], "Pb#font[122]{-}Pb 5020 GeV", "p");
	leg->Draw();

	c->SaveAs(Form("figs/NDFRatio_%dstat.pdf",i));
 }


}
