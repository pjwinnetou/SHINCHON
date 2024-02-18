{
 TCanvas* c = new TCanvas("c","c",800,700);
 gPad->SetLeftMargin(0.15);
 gPad->SetBottomMargin(0.15);
 gPad->SetRightMargin(0.03);
 gPad->SetTopMargin(0.03);
 gPad->SetTicks();


 TFile* fSTAR = new TFile("data/NchSTAROut.root","read");
 TGraphErrors* gSTAR = (TGraphErrors*)fSTAR->Get("gSTAR200");
 TFile* fALICE = new TFile("data/NchALICEOut.root","read");
 TGraphErrors* gALICE = (TGraphErrors*)fALICE->Get("gALICE5020");

 gSTAR->SetMaximum(14);

 gSTAR->GetXaxis()->SetTitleFont(43);
 gSTAR->GetYaxis()->SetTitleFont(43);
 gSTAR->GetXaxis()->SetLabelFont(43);
 gSTAR->GetYaxis()->SetLabelFont(43);

 gSTAR->GetXaxis()->SetTitleSize(32);
 gSTAR->GetYaxis()->SetTitleSize(32);
 gSTAR->GetXaxis()->SetLabelSize(28);
 gSTAR->GetYaxis()->SetLabelSize(28);

 gSTAR->SetTitle(";#it{N}_{part};(d#it{N}_{ch}/d#eta)/(#it{N}_{part}/2)");
 gSTAR->GetXaxis()->SetRangeUser(0,450);
 gSTAR->Draw("AP");

 gALICE->SetMarkerSize(1.5);
 gALICE->GetFunction("f1")->Delete();
 gALICE->Draw("P");

 TFile* fAuAu_glb = new TFile("data/Nch_glauber_AuAu_Out.root","read");
 TFile* fAuAu_tp0 = new TFile("data/Nch_Trento_p0_AuAu_Out.root","read");
 TFile* fAuAu_tp1 = new TFile("data/Nch_Trento_p1_AuAu_Out.root","read");

 TFile* fPbPb_tp0 = new TFile("data/Nch_Trento_p0_PbPb_Out.root","read");
 TFile* fPbPb_tp1 = new TFile("data/Nch_Trento_p1_PbPb_Out.root","read");

 TGraphErrors* gAuAu_glb = (TGraphErrors*)fAuAu_glb->Get("gMultScaled_Glb");
 TGraphErrors* gAuAu_tp0 = (TGraphErrors*)fAuAu_tp0->Get("gMultScaled_tp0");
 TGraphErrors* gAuAu_tp1 = (TGraphErrors*)fAuAu_tp1->Get("gMultScaled_tp1");

 TGraphErrors* gPbPb_tp0 = (TGraphErrors*)fPbPb_tp0->Get("gMultScaled_tp0");
 TGraphErrors* gPbPb_tp1 = (TGraphErrors*)fPbPb_tp1->Get("gMultScaled_tp1");

 const int nbins_npart = 9;
 double binnings[nbins_npart] = {
	0, 28, 67,
	82, 124, 180,
	255, 331, 450 };

 TProfile* hAuAu_glb = new TProfile("hAuAu_glb","",nbins_npart-1,binnings);
 TProfile* hAuAu_tp0 = new TProfile("hAuAu_tp0","",nbins_npart-1,binnings);
 TProfile* hAuAu_tp1 = new TProfile("hAuAu_tp1","",nbins_npart-1,binnings);

 TProfile* hPbPb_tp0 = new TProfile("hPbPb_tp0","",nbins_npart-1,binnings);
 TProfile* hPbPb_tp1 = new TProfile("hPbPb_tp1","",nbins_npart-1,binnings);

 TGraphErrors* gAddress[5];
 gAddress[0] = gAuAu_glb;
 gAddress[1] = gAuAu_tp0;
 gAddress[2] = gAuAu_tp1;
 gAddress[3] = gPbPb_tp0;
 gAddress[4] = gPbPb_tp1;

 TProfile* hAddress[5];
 hAddress[0] = hAuAu_glb;
 hAddress[1] = hAuAu_tp0;
 hAddress[2] = hAuAu_tp1;
 hAddress[3] = hPbPb_tp0;
 hAddress[4] = hPbPb_tp1;

 for(int i=0;i<5;i++){
	for(int j=0;j<gAddress[i]->GetN();j++){
		hAddress[i]->Fill( gAddress[i]->GetX()[j], gAddress[i]->GetY()[j] );
	}
 }


 hAuAu_glb->SetFillColorAlpha(37,0.7);
 hAuAu_glb->Draw("same,e3");
 hAuAu_tp0->SetFillColorAlpha(30,0.7);
 hAuAu_tp0->Draw("same,e3");
 hAuAu_tp1->SetFillColorAlpha(41,0.7);
 hAuAu_tp1->Draw("same,e3");

 hPbPb_tp0->SetFillColorAlpha(30,0.7);
 hPbPb_tp0->Draw("same,e3");
 hPbPb_tp1->SetFillColorAlpha(41,0.7);
 hPbPb_tp1->Draw("same,e3");

 hAuAu_glb->SetLineColor(37); hAuAu_glb->SetLineWidth(2);
 hAuAu_tp0->SetLineColor(30); hAuAu_tp0->SetLineWidth(2);
 hAuAu_tp1->SetLineColor(41); hAuAu_tp1->SetLineWidth(2);
 hPbPb_tp0->SetLineColor(30); hPbPb_tp0->SetLineWidth(2);
 hPbPb_tp1->SetLineColor(41); hPbPb_tp1->SetLineWidth(2);

 TLegend* leg = new TLegend(0.201, 0.642, 0.701, 0.943);
 leg->SetTextFont(43);
 leg->SetTextSize(28);
 leg->SetLineWidth(0);
 leg->SetFillStyle(0);

 leg->AddEntry( gSTAR, "STAR, Au#font[122]{-}Au, #sqrt{#it{s}_{NN}} = 200 GeV", "pl");
 leg->AddEntry( gALICE, "ALICE, Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5020 GeV", "pl");
 leg->AddEntry( hAuAu_tp0, "Trento #it{p} = 0", "l");
 leg->AddEntry( hAuAu_tp1, "Trento #it{p} = 1", "l");

 leg->Draw();


 c->SaveAs("figs/Nch.pdf");


}
