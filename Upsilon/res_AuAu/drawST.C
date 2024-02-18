{
 TCanvas* c = new TCanvas("c","c",800,700);
 gPad->SetLeftMargin(0.15);
 gPad->SetBottomMargin(0.15);
 gPad->SetRightMargin(0.03);
 gPad->SetTopMargin(0.03);
 gPad->SetTicks();

 TLegend* leg = new TLegend(0.2,0.6,0.9,0.9);
 leg->SetTextFont(43);
 leg->SetTextSize(32);
 leg->SetLineWidth(0);
 leg->SetFillStyle(0);

 TFile* fglb = new TFile("data/OutST_glauber.root","read");
 TProfile* hglb = (TProfile*)fglb->Get("hprofST");

 TFile* ftp0 = new TFile("data/UpsSkim_AuAu_trento_p0_Hydro.root","read");
 TProfile* htp0 = (TProfile*)ftp0->Get("hprofST");
 TProfile* htp0T = (TProfile*)ftp0->Get("hprofHydroT");

 TFile* ftp1 = new TFile("data/UpsSkim_AuAu_trento_p1_Hydro.root","read");
 TProfile* htp1 = (TProfile*)ftp1->Get("hprofST");
 TProfile* htp1T = (TProfile*)ftp1->Get("hprofHydroT");

 TFile* ftp0_PbPb = new TFile("data/UpsSkim_PbPb_trento_p0_Hydro.root","read");
 TProfile* htp0_PbPb = (TProfile*)ftp0_PbPb->Get("hprofST");
 TProfile* htp0T_PbPb = (TProfile*)ftp0_PbPb->Get("hprofHydroT");

 TFile* ftp1_PbPb = new TFile("data/UpsSkim_PbPb_trento_p1_Hydro.root","read");
 TProfile* htp1_PbPb = (TProfile*)ftp1_PbPb->Get("hprofST");
 TProfile* htp1T_PbPb = (TProfile*)ftp1_PbPb->Get("hprofHydroT");


 htp0->SetMarkerStyle(20);
 htp0->SetMarkerColorAlpha(36,0.6);
 htp0->SetLineColor(36);

 htp1->SetMarkerStyle(20);
 htp1->SetMarkerColorAlpha(46,0.6);
 htp1->SetLineColor(46);


 htp0_PbPb->SetMarkerStyle(21);
 htp0_PbPb->SetMarkerColorAlpha(30,0.6);
 htp0_PbPb->SetLineColor(30);

 htp1_PbPb->SetMarkerStyle(21);
 htp1_PbPb->SetMarkerColorAlpha(41,0.6);
 htp1_PbPb->SetLineColor(41);


 hglb->SetMarkerStyle(20);
 hglb->SetMarkerColorAlpha(1,0.6);
 hglb->SetLineColor(1);

 htp0->SetMaximum(40);

 htp0->GetXaxis()->SetTitleFont(43);
 htp0->GetYaxis()->SetTitleFont(43);
 htp0->GetXaxis()->SetLabelFont(43);
 htp0->GetYaxis()->SetLabelFont(43);

 htp0->GetXaxis()->SetTitleSize(32);
 htp0->GetYaxis()->SetTitleSize(32);
 htp0->GetXaxis()->SetLabelSize(28);
 htp0->GetYaxis()->SetLabelSize(28);

 htp0->GetXaxis()->SetNdivisions(505);
 htp0->GetYaxis()->SetNdivisions(505);

 htp0->SetTitle(";#it{N}_{part};#it{S}_{T}");
 htp0->Draw();
 htp1->Draw("same");
 hglb->Draw("same");

 leg->AddEntry( (TObject*)0, Form("SHINCHON, Au#font[122]{-}Au, 200 GeV"), "");
 leg->AddEntry( htp0, "Trento #it{p} = 0", "pl");
 leg->AddEntry( htp1, "Trento #it{p} = 1", "pl");
 leg->AddEntry( hglb, "Glauber", "pl");


 leg->Draw();

 c->SaveAs("figs/comp_ST.pdf");

 htp0_PbPb->Draw("same");
 htp1_PbPb->Draw("same"); 

 TLegend* leg1 = new TLegend(0.334, 0.161, 0.734, 0.355);
 leg1->SetTextFont(43);
 leg1->SetTextSize(32);
 leg1->SetLineWidth(0);
 leg1->SetFillStyle(0);

 leg1->AddEntry( (TObject*)0, Form("SHINCHON, Pb#font[122]{-}Pb, 5020 GeV"), "");
 leg1->AddEntry( htp0_PbPb, "Trento #it{p} = 0", "pl");
 leg1->AddEntry( htp1_PbPb, "Trento #it{p} = 1", "pl");
 leg1->Draw("same");

 c->SaveAs("figs/comp_ST_all.pdf");



 htp0T->GetXaxis()->SetTitleFont(43);
 htp0T->GetYaxis()->SetTitleFont(43);
 htp0T->GetXaxis()->SetLabelFont(43);
 htp0T->GetYaxis()->SetLabelFont(43);

 htp0T->GetXaxis()->SetTitleSize(32);
 htp0T->GetYaxis()->SetTitleSize(32);
 htp0T->GetXaxis()->SetLabelSize(28);
 htp0T->GetYaxis()->SetLabelSize(28);

 htp0T->GetXaxis()->SetNdivisions(505);
 htp0T->GetYaxis()->SetNdivisions(505);

 htp0T->SetTitle(";#it{N}_{part};T_{max,hydro} (fm/#it{c})");
 htp0T->Rebin(20);
 htp1T->Rebin(20);
 htp0T_PbPb->Rebin(20);
 htp1T_PbPb->Rebin(20);

 htp0T->SetMaximum(17);


 htp0T->SetMarkerStyle(20);
 htp0T->SetMarkerColorAlpha(36,0.6);
 htp0T->SetLineColor(36);

 htp1T->SetMarkerStyle(20);
 htp1T->SetMarkerColorAlpha(46,0.6);
 htp1T->SetLineColor(46);


 htp0T_PbPb->SetMarkerStyle(21);
 htp0T_PbPb->SetMarkerColorAlpha(30,0.6);
 htp0T_PbPb->SetLineColor(30);

 htp1T_PbPb->SetMarkerStyle(21);
 htp1T_PbPb->SetMarkerColorAlpha(41,0.6);
 htp1T_PbPb->SetLineColor(41);


 htp0T->Draw();
 htp1T->Draw("same");
 htp0T_PbPb->Draw("same");
 htp1T_PbPb->Draw("same");

 TLegend* leg2 = new TLegend(0.2,0.7,0.9,0.9);
 leg2->SetTextFont(43);
 leg2->SetTextSize(28);
 leg2->SetLineWidth(0);
 leg2->SetFillStyle(0);

 leg2->AddEntry( htp0, "Au#font[122]{-}Au Trento #it{p} = 0", "pl");
 leg2->AddEntry( htp1, "Au#font[122]{-}Au Trento #it{p} = 1", "pl");
 leg2->AddEntry( htp0_PbPb, "Pb#font[122]{-}Pb Trento #it{p} = 0", "pl");
 leg2->AddEntry( htp1_PbPb, "Pb#font[122]{-}Pb Trento #it{p} = 1", "pl");
 leg2->Draw();

 c->SaveAs("figs/hydroT.pdf");
}
