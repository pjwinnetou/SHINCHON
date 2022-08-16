void draw_RpA(){


 TFile* fin = new TFile("NMFdata/FDCOR_out.root","read");

 TCanvas* c = new TCanvas("c","c",800,600);
 gPad->SetLeftMargin(0.13);
 gPad->SetBottomMargin(0.13);
 gPad->SetTopMargin(0.03);
 gPad->SetRightMargin(0.03);
 gPad->SetTicks();
 gStyle->SetOptStat(0);

 const int nsyst = 3;
 const int nstat = 3; 

 const char fname[nsyst][100] = {
	"pPb",
	"pO",
	"OO" };

 const char dname[nsyst][100] = {
        "p#font[122]{-}Pb",
        "pO",
        "OO" };



 TGraphErrors* gRpPbMult[nsyst][nstat];
 TGraphErrors* gRpPbMultST[nsyst][nstat];
 TGraphErrors* gRpPbPt[nsyst][nstat];


 for(int i=0;i<nsyst;i++){
	for(int j=0;j<nstat;j++){
		gRpPbMult[i][j] = (TGraphErrors*)fin->Get(Form("g_R%s_multdep_fdcor_%dS",fname[i],j+1));
		gRpPbPt[i][j] = (TGraphErrors*)fin->Get(Form("g_R%s_ptdep_fdcor_%dS",fname[i],j+1));
		gRpPbMultST[i][j] = (TGraphErrors*)fin->Get(Form("g_R%s_multSTdep_fdcor_%dS",fname[i],j+1));		

		gRpPbMult[i][j]->GetYaxis()->SetTitle("Nuclear modification factor");
		gRpPbPt[i][j]->GetYaxis()->SetTitle("Nuclear modification factor");
		gRpPbMultST[i][j]->GetYaxis()->SetTitle("Nuclear modification factor");

		gRpPbMult[i][j]->GetXaxis()->SetTitle("N_{ch}");
		gRpPbPt[i][j]->GetXaxis()->SetTitle("#font[12]{p}_{T} (GeV/#font[12]{c})");
		gRpPbMultST[i][j]->GetXaxis()->SetTitle("N_{ch}/#LT#font[12]{s}_{T}#GT");

		gRpPbMult[i][j]->GetXaxis()->SetTitleSize(0.055);
		gRpPbMultST[i][j]->GetXaxis()->SetTitleSize(0.055);
		gRpPbPt[i][j]->GetXaxis()->SetTitleSize(0.055);

                gRpPbMult[i][j]->GetYaxis()->SetTitleSize(0.055);
                gRpPbMultST[i][j]->GetYaxis()->SetTitleSize(0.055);
                gRpPbPt[i][j]->GetYaxis()->SetTitleSize(0.055);

                gRpPbMult[i][j]->GetXaxis()->SetLabelSize(0.05);
                gRpPbMultST[i][j]->GetXaxis()->SetLabelSize(0.05);
                gRpPbPt[i][j]->GetXaxis()->SetLabelSize(0.05);

                gRpPbMult[i][j]->GetYaxis()->SetLabelSize(0.05);
                gRpPbMultST[i][j]->GetYaxis()->SetLabelSize(0.05);
                gRpPbPt[i][j]->GetYaxis()->SetLabelSize(0.05);

		gRpPbMult[i][j] ->SetMaximum(1.5);
		gRpPbMultST[i][j] ->SetMaximum(1.5);
		gRpPbPt[i][j] ->SetMaximum(1.5);

	}
 }


 TLegend* legMult = new TLegend(0.332, 0.749, 0.849, 0.949);
 legMult->SetFillColorAlpha(0,0);
 legMult->SetLineWidth(0.0);

 for(int i=0;i<nsyst;i++){
	legMult->Clear();
	legMult->AddEntry( (TObject*)0, Form("SHINCHON, %s #sqrt{#it{s}_{NN}} = 8.16 TeV",dname[i]), "");
	for(int j=0;j<nstat;j++){
		gRpPbMult[i][j]->SetFillColor(j+1);
		gRpPbMult[i][j]->SetMinimum(0.0);
		if( j==0 ) gRpPbMult[i][j]->Draw("ae3");
		gRpPbMult[i][j]->Draw("e3");
		legMult->AddEntry( gRpPbMult[i][j], Form("#varUpsilon(%dS)",j+1), "f");
	}
	legMult->Draw();
	c->SaveAs(Form("NMFplots/R%s_mult.pdf",fname[i]));
 }

 TLegend* legMult1 = new TLegend(0.311, 0.732, 0.834, 0.917);
 legMult1->SetFillColorAlpha(0,0);
 legMult1->SetLineWidth(0.0);

 for(int j=0;j<nstat;j++){
	legMult1->Clear();
	legMult1->AddEntry( (TObject*)0, Form("SHINCHON, #sqrt{#it{s}_{NN}} = 8.16 TeV, #varUpsilon(%dS)",j+1), "");
	for(int i=0;i<nsyst;i++){
		gRpPbMult[i][j]->SetFillColor(i+1);
		if( i==0 ) gRpPbMult[i][j]->Draw("ae3");
		gRpPbMult[i][j]->Draw("e3");
		legMult1->AddEntry( gRpPbMult[i][j], Form("%s",dname[i]), "f");
	}
	legMult1->Draw();
	c->SaveAs(Form("NMFplots/R_mult_%dS.pdf",j+1));
 }



 TLegend* legMultST = new TLegend(0.332, 0.749, 0.849, 0.949);
 legMultST->SetFillColorAlpha(0,0);
 legMultST->SetLineWidth(0.0);

 for(int i=0;i<nsyst;i++){
        legMultST->Clear();
        legMultST->AddEntry( (TObject*)0, Form("SHINCHON, %s #sqrt{#it{s}_{NN}} = 8.16 TeV",dname[i]), "");
        for(int j=0;j<nstat;j++){
                gRpPbMultST[i][j]->SetFillColor(j+1);
                gRpPbMultST[i][j]->SetMinimum(0.0);
                if( j==0 ) gRpPbMultST[i][j]->Draw("ae3");
                gRpPbMultST[i][j]->Draw("e3");
                legMultST->AddEntry( gRpPbMultST[i][j], Form("#varUpsilon(%dS)",j+1), "f");
        }
        legMultST->Draw();
	c->SaveAs(Form("NMFplots/R%s_multST.pdf",fname[i]));
 }

 TLegend* legMultST1 = new TLegend(0.311, 0.732, 0.834, 0.947);
 legMultST1->SetFillColorAlpha(0,0);
 legMultST1->SetLineWidth(0.0);

 for(int j=0;j<nstat;j++){
        legMultST1->Clear();
        legMultST1->AddEntry( (TObject*)0, Form("SHINCHON, #sqrt{#it{s}_{NN}} = 8.16 TeV, #varUpsilon(%dS)",j+1), "");
        for(int i=0;i<nsyst;i++){
                gRpPbMultST[i][j]->SetFillColor(i+1);
                if( i==0 ) gRpPbMultST[i][j]->Draw("ae3");
                gRpPbMultST[i][j]->Draw("e3");
                legMultST1->AddEntry( gRpPbMultST[i][j], Form("%s",dname[i]), "f");
        }
        legMultST1->Draw();
	c->SaveAs(Form("NMFplots/R_multST_%dS.pdf",j+1));
 }



 TLegend* legPt = new TLegend(0.156, 0.721, 0.680, 0.937);
 legPt->SetFillColorAlpha(0,0);
 legPt->SetLineWidth(0.0);

 for(int i=0;i<nsyst;i++){
        legPt->Clear();
        legPt->AddEntry( (TObject*)0, Form("SHINCHON, %s #sqrt{#it{s}_{NN}} = 8.16 TeV",dname[i]), "");
        for(int j=0;j<nstat;j++){
                gRpPbPt[i][j]->SetFillColor(j+1);
                gRpPbPt[i][j]->SetMinimum(0.0);
                if( j==0 ) gRpPbPt[i][j]->Draw("ae3");
                gRpPbPt[i][j]->Draw("e3");
                legPt->AddEntry( gRpPbPt[i][j], Form("#varUpsilon(%dS)",j+1), "f");
        }
        legPt->Draw();
	c->SaveAs(Form("NMFplots/R%s_Pt.pdf",fname[i]));
 }

 TLegend* legPt1 = new TLegend(0.156, 0.721, 0.680, 0.937);
 legPt1->SetFillColorAlpha(0,0);
 legPt1->SetLineWidth(0.0);

 for(int j=0;j<nstat;j++){
        legPt1->Clear();
        legPt1->AddEntry( (TObject*)0, Form("SHINCHON, #sqrt{#it{s}_{NN}} = 8.16 TeV, #varUpsilon(%dS)",j+1), "");
        for(int i=0;i<nsyst;i++){
                gRpPbPt[i][j]->SetFillColor(i+1);
                if( i==0 ) gRpPbPt[i][j]->Draw("ae3");
                gRpPbPt[i][j]->Draw("e3");
                legPt1->AddEntry( gRpPbPt[i][j], Form("%s",dname[i]), "f");
        }
        legPt1->Draw();
	c->SaveAs(Form("NMFplots/R_Pt_%dS.pdf",j+1));
 }

 TFile* fCMS = new TFile("NMFdata//HEPData-ins2037640-v1-root.root","read");
 TDirectoryFile* dir;
 TH1F* hRpPb_CMS_cntl[nstat];
 TH1F* hRpPb_CMS_stat[nstat];
 TH1F* hRpPb_CMS_syst[nstat];

 for(int i=0;i<nstat;i++){
	dir = (TDirectoryFile*)fCMS->GetDirectory(Form("Table %d",13+i));
	hRpPb_CMS_cntl[i] = (TH1F*)dir->Get("Hist1D_y1");
	hRpPb_CMS_stat[i] = (TH1F*)dir->Get("Hist1D_y1_e1");
	hRpPb_CMS_syst[i] = (TH1F*)dir->Get("Hist1D_y1_e2");

	for(int j=0;j<hRpPb_CMS_cntl[i]->GetNbinsX();j++){
		hRpPb_CMS_stat[i]->SetBinContent( j+1, hRpPb_CMS_cntl[i]->GetBinContent(j+1) );
		hRpPb_CMS_syst[i]->SetBinContent( j+1, hRpPb_CMS_cntl[i]->GetBinContent(j+1) );
	}

	hRpPb_CMS_stat[i]->SetLineColor(i+1);
	hRpPb_CMS_syst[i]->SetLineColor(i+1);
	hRpPb_CMS_stat[i]->SetMarkerColor(i+1);
	hRpPb_CMS_syst[i]->SetMarkerColor(i+1);
	hRpPb_CMS_stat[i]->SetMarkerStyle(20+i);
	hRpPb_CMS_syst[i]->SetFillStyle(0);
 }

 legPt->Clear();
 legPt->SetNColumns(2);
 legPt->SetHeader( Form("%s #sqrt{#it{s}_{NN}} = 8.16 TeV",dname[0]) );
 legPt->AddEntry( (TObject*)0, "SHINCHON", "");
 legPt->AddEntry( (TObject*)0, "CMS", "");

 gRpPbPt[0][0]->Draw("ae3");
 for(int i=0;i<nstat;i++){
	hRpPb_CMS_stat[i]->Draw("same");
	hRpPb_CMS_syst[i]->Draw("same,e2");
 }
 for(int i=0;i<nstat;i++){
	gRpPbPt[0][i]->SetFillColor(i+1);
	gRpPbPt[0][i]->SetLineColor(i+1);
	gRpPbPt[0][i]->Draw("e3");
 }
 for(int i=0;i<nstat;i++){
	legPt->AddEntry( gRpPbPt[0][i], Form("#varUpsilon(%dS)",i+1), "f");
	legPt->AddEntry( hRpPb_CMS_stat[i], Form("#varUpsilon(%dS)",i+1), "pl");
 }

 legPt->Draw();

 c->SaveAs("NMFplots/RpPb_Pt_compCMS.pdf");



}
