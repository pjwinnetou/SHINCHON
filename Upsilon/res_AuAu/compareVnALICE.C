{
 TCanvas* c = new TCanvas("c","c",800,700);
 gPad->SetLeftMargin(0.15);
 gPad->SetBottomMargin(0.15);
 gPad->SetRightMargin(0.03);
 gPad->SetTopMargin(0.03);
 gPad->SetTicks();

 TLegend* leg = new TLegend(0.184, 0.64, 0.684, 0.940);
 leg->SetTextFont(43);
 leg->SetTextSize(32);
 leg->SetLineWidth(0);
 leg->SetFillStyle(0);

 const int nmult = 8;
 char cname[nmult][100] = {
    "0#font[122]{-}5%",
    "5#font[122]{-}10%",
    "10#font[122]{-}20%",
    "20#font[122]{-}30%",
    "30#font[122]{-}40%",
    "40#font[122]{-}50%",
    "50#font[122]{-}60%",
    "60#font[122]{-}70%" };

 TFile* fdata = new TFile("data/VnALICE5020GeV.root","read");
 TGraphErrors* gdata[nmult];

 TFile* ftp0 = new TFile("data/VnOut_PbPb_tp0.root","read");
 TFile* ftp1 = new TFile("data/VnOut_PbPb_tp1.root","read");

 TGraphErrors* gtp0[nmult];
 TGraphErrors* gtp1[nmult];

 for(int i=0;i<nmult;i++){
	gdata[i] = (TGraphErrors*)fdata->Get(Form("gVnALICE_%d",i));
	gtp0[i] = (TGraphErrors*)ftp0->Get(Form("gVn_tp0_%d",i));
	gtp1[i] = (TGraphErrors*)ftp1->Get(Form("gVn_tp1_%d",i));

	gtp0[i]->SetLineColorAlpha(36,0.7);
	gtp1[i]->SetLineColorAlpha(46,0.7);
 }

 for(int i=0;i<nmult;i++){
	gdata[i]->Draw("AP");
	gtp0[i]->Draw("c"); 
	gtp1[i]->Draw("c"); 	

	leg->Clear();
	leg->AddEntry( (TObject*)0, Form("Pb#font[122]{-}Pb, 5020 GeV, %s",cname[i]), "");
    leg->AddEntry( gdata[i], "ALICE, v_{2}{2}, JHEP 07 (2018) 103", "p");
	leg->AddEntry( gtp0[i], "Trento #it{p} = 0", "l");
	leg->AddEntry( gtp1[i], "Trento #it{p} = 1", "l");

    leg->Draw();

    c->SaveAs(Form("figs/vn/vn_ALICE_%d.pdf",i));
 }


}
