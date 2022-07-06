TGraphErrors* convertHist( TH1D* h){
 TGraphErrors* g = new TGraphErrors();
 int np=0;
 for(int i=0;i<h->GetNbinsX();i++){
        if( h->GetBinContent(i+1) == 0 ) continue;

        g->SetPoint( np, h->GetBinCenter(i+1), h->GetBinContent(i+1) );
        g->SetPointError( np, h->GetBinWidth(i+1), h->GetBinError(i+1) );
        np++;
 }

 return g;
}


Double_t fTsallis_v2(Double_t *x, Double_t *fpar){
  Float_t xx = x[0];
  Double_t q = fpar[0];
  Double_t T = fpar[1];
  Double_t c = fpar[2];
  Double_t Ymass = fpar[3];
  Double_t mT = TMath::Sqrt(Ymass*Ymass + xx*xx);
  Double_t pow = TMath::Power((1+(q-1)*mT/T),(-q/(q-1)));

  Double_t f = c*mT*xx*pow;
  return f;
}

void Feeddonw_Correction_v1(){

 TFile* fin = new TFile("../data/outfile_UpsSkim_PhiAng_pPb_2_InitPosGlauber_0000_1000.root","read"); //CHANGE DIR
 TFile* ffd_corr = new TFile("../data/Results_FD_Bottomonium.root","read"); //CHANGE DIR

 TF1* frac2sTo1s = (TF1*)ffd_corr->Get("frac2sTo1s");
 TF1* frac3sTo1s = (TF1*)ffd_corr->Get("frac3sTo1s");
 TF1* frac1pTo1s = (TF1*)ffd_corr->Get("frac1pTo1s");
 TF1* frac2pTo1s = (TF1*)ffd_corr->Get("frac2pTo1s");
 TF1* frac3pTo1s = (TF1*)ffd_corr->Get("frac3pTo1s");

 TF1* frac3sTo2s = (TF1*)ffd_corr->Get("frac3sTo2s");
 TF1* frac2pTo2s = (TF1*)ffd_corr->Get("frac2pTo2s");
 TF1* frac3pTo2s = (TF1*)ffd_corr->Get("frac3pTo2s");

 TF1* frac3pTo3s = (TF1*)ffd_corr->Get("frac3pTo3s");


 const int ntot_fdc = 9;
 TF1* fFeeddonw[ntot_fdc];
 fFeeddonw[0] = (TF1*)frac2sTo1s->Clone();
 fFeeddonw[1] = (TF1*)frac3sTo1s->Clone();
 fFeeddonw[2] = (TF1*)frac1pTo1s->Clone();
 fFeeddonw[3] = (TF1*)frac2pTo1s->Clone();
 fFeeddonw[4] = (TF1*)frac3pTo1s->Clone();

 fFeeddonw[5] = (TF1*)frac3sTo2s->Clone();
 fFeeddonw[6] = (TF1*)frac2pTo2s->Clone();
 fFeeddonw[7] = (TF1*)frac3pTo2s->Clone();

 fFeeddonw[8] = (TF1*)frac3pTo3s->Clone();

 int iparStat[ntot_fdc] = {
	2, 3, 2, 3, 3,
	3, 3, 3,
	3 };

 int fparStat[ntot_fdc] = {
	1, 1, 1, 1, 1,
	2, 2, 2,
	3 };


 TProfile* hRpA_pPb_pt_1S = (TProfile*)fin->Get("hprofRAA_pT_rl_1S");
 TProfile* hRpA_pPb_pt_2S = (TProfile*)fin->Get("hprofRAA_pT_rl_2S");
 TProfile* hRpA_pPb_pt_3S = (TProfile*)fin->Get("hprofRAA_pT_rl_3S");
 TProfile* hRpA_pPb_pt[3];

 hRpA_pPb_pt[0] = (TProfile*)hRpA_pPb_pt_1S->Clone(Form("hRpA_pPb_pT_1S"));
 hRpA_pPb_pt[1] = (TProfile*)hRpA_pPb_pt_2S->Clone(Form("hRpA_pPb_pT_2S"));
 hRpA_pPb_pt[2] = (TProfile*)hRpA_pPb_pt_3S->Clone(Form("hRpA_pPb_pT_3S"));


 const int nstat = 3;
 TH1D* hRpA_pPb_pt_corr[nstat];

 double owncont;
 double fdcont[nstat];

 for(int s=0;s<nstat;s++){
	hRpA_pPb_pt_corr[s] = new TH1D(Form("hRpA_pPb_pt_corr_%dS",s+1),"",100,0,20);

	for(int p=0;p<hRpA_pPb_pt_corr[s]->GetNbinsX();p++){
		owncont = 100.;
		for(int ss=0;ss<nstat;ss++){
			fdcont[ss] = 0.;
		}

		for(int i=0;i<ntot_fdc;i++){
			if( fparStat[i] != s+1 ) continue;
			owncont -= fFeeddonw[i]->Integral( hRpA_pPb_pt_corr[s]->GetBinCenter(p+1) - hRpA_pPb_pt_corr[s]->GetBinWidth(p+1)/2.0,
							   hRpA_pPb_pt_corr[s]->GetBinCenter(p+1) + hRpA_pPb_pt_corr[s]->GetBinWidth(p+1)/2.0 ) / hRpA_pPb_pt_corr[s]->GetBinWidth(p+1);
			if( iparStat[i] < fparStat[i] ) continue;
			fdcont[ iparStat[i]-1 ] +=
				fFeeddonw[i]->Integral( hRpA_pPb_pt_corr[s]->GetBinCenter(p+1) - hRpA_pPb_pt_corr[s]->GetBinWidth(p+1)/2.0,
                                                        hRpA_pPb_pt_corr[s]->GetBinCenter(p+1) + hRpA_pPb_pt_corr[s]->GetBinWidth(p+1)/2.0 ) / hRpA_pPb_pt_corr[s]->GetBinWidth(p+1);
		}
		hRpA_pPb_pt_corr[s]->SetBinContent( p+1, owncont*hRpA_pPb_pt[s]->GetBinContent( p+1 ) );
		hRpA_pPb_pt_corr[s]->SetBinError( p+1, owncont*hRpA_pPb_pt[s]->GetBinError( p+1 ) );
		for(int ss=0;ss<nstat;ss++){
			hRpA_pPb_pt_corr[s]->AddBinContent( p+1, fdcont[ss]*hRpA_pPb_pt[ss]->GetBinContent( p+1 ) );
			hRpA_pPb_pt_corr[s]->SetBinError( p+1, sqrt( pow( hRpA_pPb_pt_corr[s]->GetBinError(p+1),2) + pow( fdcont[ss]*hRpA_pPb_pt[ss]->GetBinError(p+1),2) ) );
		}
	}
	hRpA_pPb_pt_corr[s]->Scale(1./100. );
 }

// pt dep;

 TH1D* hRpA_pPb_pt_corr_2St1S = new TH1D("hRpA_pPb_pt_corr_2St1S","",100,0,20); //dedicated study
 double fdcont_2Sto1S;

 for(int p=0;p<hRpA_pPb_pt_corr_2St1S->GetNbinsX();p++){
	owncont = 100.;
	fdcont_2Sto1S = 0.;
	for(int i=0;i<ntot_fdc;i++){
		if( fparStat[i] == 1 && iparStat[i] == 2 ){
			owncont -= fFeeddonw[i]->Integral( hRpA_pPb_pt_corr_2St1S->GetBinCenter(p+1) - hRpA_pPb_pt_corr_2St1S->GetBinWidth(p+1)/2.0,
							   hRpA_pPb_pt_corr_2St1S->GetBinCenter(p+1) + hRpA_pPb_pt_corr_2St1S->GetBinWidth(p+1)/2.0 ) / hRpA_pPb_pt_corr_2St1S->GetBinWidth(p+1);
			fdcont_2Sto1S += fFeeddonw[i]->Integral( hRpA_pPb_pt_corr_2St1S->GetBinCenter(p+1) - hRpA_pPb_pt_corr_2St1S->GetBinWidth(p+1)/2.0,
								 hRpA_pPb_pt_corr_2St1S->GetBinCenter(p+1) + hRpA_pPb_pt_corr_2St1S->GetBinWidth(p+1)/2.0 ) / hRpA_pPb_pt_corr_2St1S->GetBinWidth(p+1); 
		}
	}
	hRpA_pPb_pt_corr_2St1S->SetBinContent( p+1, owncont*hRpA_pPb_pt[0]->GetBinContent( p+1 ) );
	hRpA_pPb_pt_corr_2St1S->SetBinError( p+1, owncont*hRpA_pPb_pt[0]->GetBinError( p+1 ) );

	hRpA_pPb_pt_corr_2St1S->AddBinContent( p+1, fdcont_2Sto1S*hRpA_pPb_pt[1]->GetBinContent( p+1 ) );
	hRpA_pPb_pt_corr_2St1S->SetBinError( p+1, sqrt( pow( hRpA_pPb_pt_corr_2St1S->GetBinError(p+1),2) +
		pow( fdcont_2Sto1S*hRpA_pPb_pt[1]->GetBinError(p+1),2) ) );
 }
 hRpA_pPb_pt_corr_2St1S->Scale(1./100.);


// dedicated pt dep study

 double invMass[ntot_fdc] = {
        10.02, 10.36, 9.892, 10.26, 10.51, //2s 3s 1p 2p 3p
        10.36, 10.26, 10.51,
        10.51 };

 TF1 *fInitialUpsilon = new TF1("fInitialUpsilon",fTsallis_v2,0,30,4);
 double Feeddonw_corr_meanpt[ntot_fdc];
 double ptmin, ptmax;

 int LARGE=10000;
 for(int i=0;i<ntot_fdc;i++){
	fInitialUpsilon->SetParameters(  1.06450e+00 ,  7.97649e-01 , 100, invMass[i] );	

	Feeddonw_corr_meanpt[i] = 0.;
	for(int j=0;j<LARGE;j++){
		ptmin = (double)j*20./LARGE;
		ptmax = (double)j*20./LARGE + 20./LARGE;

		Feeddonw_corr_meanpt[i] += fInitialUpsilon->Integral( ptmin, ptmax ) *
			fFeeddonw[i]->Integral( ptmin, ptmax ) / (ptmax - ptmin);
	}
	Feeddonw_corr_meanpt[i] /= fInitialUpsilon->Integral( 0, 20 );
	Feeddonw_corr_meanpt[i] /= 100.;
	cout << Feeddonw_corr_meanpt[i] << endl;
 }

 
 TProfile* hRpA_pPb_mult_1S = (TProfile*)fin->Get("hprofRAA_mult_rl_1S");
 TProfile* hRpA_pPb_mult_2S = (TProfile*)fin->Get("hprofRAA_mult_rl_2S");
 TProfile* hRpA_pPb_mult_3S = (TProfile*)fin->Get("hprofRAA_mult_rl_3S");
 TProfile* hRpA_pPb_mult[3];

 hRpA_pPb_mult[0] = (TProfile*)hRpA_pPb_mult_1S->Clone(Form("hRpA_pPb_mult_1S"));
 hRpA_pPb_mult[1] = (TProfile*)hRpA_pPb_mult_2S->Clone(Form("hRpA_pPb_mult_2S"));
 hRpA_pPb_mult[2] = (TProfile*)hRpA_pPb_mult_3S->Clone(Form("hRpA_pPb_mult_3S"));


 TH1D* hRpA_pPb_mult_corr[nstat];
 for(int s=0;s<nstat;s++){
	hRpA_pPb_mult_corr[s] = new TH1D(Form("hRpA_pPb_mult_corr_%dS",s+1),"",100,0,100);
	for(int p=0;p<hRpA_pPb_mult_corr[s]->GetNbinsX();p++){
		hRpA_pPb_mult_corr[s]->SetBinContent( p+1, hRpA_pPb_mult[s]->GetBinContent(p+1) );
		hRpA_pPb_mult_corr[s]->SetBinError( p+1, hRpA_pPb_mult[s]->GetBinError(p+1) );
		for(int i=0;i<ntot_fdc;i++){
			if( fparStat[i] != s+1 ) continue;
			hRpA_pPb_mult_corr[s]->AddBinContent( p+1, -hRpA_pPb_mult[s]->GetBinContent(p+1)*Feeddonw_corr_meanpt[i] );
			hRpA_pPb_mult_corr[s]->SetBinError( p+1, sqrt( pow( hRpA_pPb_mult_corr[s]->GetBinError(p+1),2) +
				pow( hRpA_pPb_mult[s]->GetBinError(p+1)*Feeddonw_corr_meanpt[i],2) ) );

			if( iparStat[i] < fparStat[i] ) continue;
			hRpA_pPb_mult_corr[s]->AddBinContent( p+1, hRpA_pPb_mult[iparStat[i]-1]->GetBinContent(p+1)*Feeddonw_corr_meanpt[i] );
			hRpA_pPb_mult_corr[s]->SetBinError( p+1, sqrt( pow( hRpA_pPb_mult_corr[s]->GetBinError(p+1),2) +
				pow( hRpA_pPb_mult[iparStat[i]-1]->GetBinError(p+1)*Feeddonw_corr_meanpt[i],2) ) );
		}
	}
 }
// mult dep



 TCanvas* c = new TCanvas("c","c",800,600);
 gPad->SetLeftMargin(0.13);
 gPad->SetBottomMargin(0.13);
 gPad->SetTopMargin(0.03);
 gPad->SetRightMargin(0.03);
 gPad->SetTicks();
 gStyle->SetOptStat(0);


 TGraphErrors* gpt[nstat];
 TGraphErrors* gpt_fdc[nstat];
 TGraphErrors* gmult[nstat];
 TGraphErrors* gmult_fdc[nstat];

 TLegend* leg = new TLegend(0.463, 0.754, 0.993, 0.944);
 leg->SetFillColorAlpha(0,0);
 leg->SetLineWidth(0.0);


 for(int s=0;s<nstat;s++){
	gpt[s] = new TGraphErrors();
	gpt_fdc[s] = new TGraphErrors();
	gmult[s] = new TGraphErrors();
	gmult_fdc[s] = new TGraphErrors();

	gpt[s] = convertHist( hRpA_pPb_pt[s] );
	gpt_fdc[s] = convertHist( hRpA_pPb_pt_corr[s] );
	gmult[s] = convertHist( hRpA_pPb_mult[s] );
	gmult_fdc[s] = convertHist( hRpA_pPb_mult_corr[s] );

	gpt[s]->SetFillColorAlpha(kBlack,0.5);
	gpt[s]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
	gpt[s]->GetYaxis()->SetTitle("R_{pPb}");
	gpt[s]->SetMinimum(0.1);
	gpt[s]->SetMaximum(1.3);

	gpt_fdc[s]->SetFillColorAlpha(kRed,0.5);


        gmult[s]->SetFillColorAlpha(kBlack,0.5);
        gmult[s]->GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#eta#GT");
        gmult[s]->GetYaxis()->SetTitle("R_{pPb}");
        gmult[s]->SetMinimum(0.1);
        gmult[s]->SetMaximum(1.3);

        gmult_fdc[s]->SetFillColorAlpha(kRed,0.5);

	gpt[s]->Draw("ae3");
	gpt_fdc[s]->Draw("e3");

	leg->Clear();
	leg->AddEntry( (TObject*)0, "SHINCHON", "");
	leg->AddEntry( (TObject*)0, Form("#varUpsilon(%dS), p#font[122]{-}Pb #sqrt{#it{s}_{NN}} = 8.16 TeV",s+1), "");
	leg->AddEntry( gpt[s], "Before feed-down correction","f");
	leg->AddEntry( gpt_fdc[s], "After feed-down correction","f");
	leg->Draw();

	c->SaveAs(Form("figs/gpt_%dS.pdf",s+1));	



        gmult[s]->Draw("ae3");
        gmult_fdc[s]->Draw("e3");
	leg->Draw();

        c->SaveAs(Form("figs/gmult_%dS.pdf",s+1));
 }

 TFile* fout = new TFile("FDCOR_out.root","recreate");
 for(int s=0;s<nstat;s++){
        gpt[s]->SetName(Form("g_RpPb_ptdep_%dS",s+1));
        gpt_fdc[s]->SetName(Form("g_RpPb_ptdep_fdcor_%dS",s+1));
        gmult[s]->SetName(Form("g_RpPb_multdep_%dS",s+1));
        gmult_fdc[s]->SetName(Form("g_RpPb_multdep_fdcor_%dS",s+1));
 
        gpt[s]->Write();
        gpt_fdc[s]->Write();
        gmult[s]->Write();
        gmult_fdc[s]->Write();
 }


 TGraphErrors* g_pt_corr_2St1S = (TGraphErrors*)convertHist( hRpA_pPb_pt_corr_2St1S );
 g_pt_corr_2St1S->SetFillColorAlpha(kRed, 0.5 );
 gpt[1]->SetFillColorAlpha(kBlue, 0.5 );
 gpt[0]->Draw("ae3");
 gpt[1]->Draw("e3");
 g_pt_corr_2St1S->Draw("e3");
 gpt_fdc[0]->SetFillColorAlpha(kMagenta, 0.5 );
 gpt_fdc[0]->Draw("e3");


 TLegend* leg1 = new TLegend(0.387, 0.187, 0.917, 0.476);
 leg1->SetFillColorAlpha(0,0);
 leg1->SetLineWidth(0.0);
 leg1->AddEntry( (TObject*)0, "SHINCHON", "");
 leg1->AddEntry( (TObject*)0, Form("p#font[122]{-}Pb #sqrt{#it{s}_{NN}} = 8.16 TeV"), "");
 leg1->AddEntry( gpt[0], "#varUpsilon(1S), uncorrected", "f");
 leg1->AddEntry( gpt[1], "#varUpsilon(2S), uncorrected", "f");
 leg1->AddEntry( g_pt_corr_2St1S, "#varUpsilon(1S), corrected with #varUpsilon(2S)", "f");
 leg1->AddEntry( gpt_fdc[0], "#varUpsilon(1S), fully corrected", "f");
 leg1->Draw();

 c->SaveAs(Form("figs/g_pt_corr_2Sto1S.pdf"));


 fout->cd();
 g_pt_corr_2St1S->SetName("g_RpPb_ptdep_1S_corw2S");
 g_pt_corr_2St1S->Write();


}
