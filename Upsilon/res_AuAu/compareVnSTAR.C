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

 const int nmult = 6;
 char cname[nmult][100] = {
    "0#font[122]{-}10%",
    "10#font[122]{-}20%",
    "20#font[122]{-}30%",
    "30#font[122]{-}40%",
    "40#font[122]{-}50%",
    "50#font[122]{-}60%" };

 TFile* fdata = new TFile("data/VnSTAR200GeV.root","read");
 TGraphErrors* gdata[nmult];

 TFile* fglb = new TFile("data/VnOut_AuAu_glb.root","read");
 TFile* ftp0 = new TFile("data/VnOut_AuAu_tp0.root","read");
 TFile* ftp1 = new TFile("data/VnOut_AuAu_tp1.root","read");

 TGraphErrors* gglb[nmult];
 TGraphErrors* gtp0[nmult];
 TGraphErrors* gtp1[nmult];

 for(int i=0;i<nmult;i++){
	gdata[i] = (TGraphErrors*)fdata->Get(Form("gVnSTAR_%d",i));
	gglb[i] = (TGraphErrors*)fglb->Get(Form("gVn_glb_%d",i));
	gtp0[i] = (TGraphErrors*)ftp0->Get(Form("gVn_tp0_%d",i));
	gtp1[i] = (TGraphErrors*)ftp1->Get(Form("gVn_tp1_%d",i));

	gglb[i]->SetLineColorAlpha(30,0.7);
	gtp0[i]->SetLineColorAlpha(36,0.7);
	gtp1[i]->SetLineColorAlpha(46,0.7);
 }

 for(int i=0;i<nmult;i++){
	gdata[i]->Draw("AP");
	gglb[i]->Draw("c");
	gtp0[i]->Draw("c"); 
	gtp1[i]->Draw("c"); 	

	leg->Clear();
	leg->AddEntry( (TObject*)0, Form("Au#font[122]{-}Au, 200 GeV, %s",cname[i]), "");
    leg->AddEntry( gdata[i], "PHENIX, v_{2}{2}, PRC 80 (2009) 024909", "p");
    leg->AddEntry( gglb[i], "Glauber", "l");
	leg->AddEntry( gtp0[i], "Trento #it{p} = 0", "l");
	leg->AddEntry( gtp1[i], "Trento #it{p} = 1", "l");

    leg->Draw();

    c->SaveAs(Form("figs/vn/vn_STAR_%d.pdf",i));
 }


}
