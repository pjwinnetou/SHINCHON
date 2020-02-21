void Draw_time(){

	gStyle->SetPalette(55);
	gStyle->SetOptStat(0);

	const bool bLOG = false;

	//const int ntime = 87;
	//const int ntime = 172;

	TFile *infile = new TFile("outfile.root","read");

	TH1D *h1_time = (TH1D*)infile->Get("h1_time");

	const int ntime = int(h1_time->GetEntries());

	TH2D *h2d[ntime];

	for (int it=0; it<ntime; it++){
		h2d[it] = (TH2D*)infile->Get(Form("h2_evt%05d_t%05d",0,it));
	}

	float max = h2d[0]->GetMaximum();
	float min = h2d[ntime-1]->GetMinimum();
	float integral = h2d[0]->Integral();

	TCanvas *c1 = new TCanvas("c1","c1",1.2*1*500,500);

	//c1->cd(1);
	//SetPadStyle();
	gPad->SetRightMargin(0.13);
	gPad->SetLeftMargin(0.09);
	gPad->SetTopMargin(0.08);
	gPad->SetBottomMargin(0.1);
	if ( bLOG ) gPad->SetLogz();

	/*
	c1->cd(2);
	SetPadStyle();
	gPad->SetRightMargin(0.135);
	gPad->SetLeftMargin(0.13);
	gPad->SetTopMargin(0.08);
	if ( bLOG ) gPad->SetLogz();
	*/


	gSystem->Exec("mkdir -p gif/gif_run00");

	for (int it=0; it<ntime; it++){
	//for (int it=0; it<1; it++){
		//c1->cd(1);
		c1->Clear();

		//c1->cd(2);
		//c1->Clear();

		//if ( h2d[it]->Integral()<0.01*integral ) continue;

		//htmp = (TH1D*)gPad->DrawFrame(-15,-15,15,15);
		TH1D *htmp = (TH1D*)gPad->DrawFrame(-10.,-10.,10.,10.);
		htmp->GetXaxis()->SetTitle("x [fm]");
		htmp->GetXaxis()->SetLabelFont(43);
		htmp->GetXaxis()->SetLabelSize(20);
		htmp->GetXaxis()->SetTitleFont(43);
		htmp->GetXaxis()->SetTitleSize(20);

		htmp->GetYaxis()->SetTitle("y [fm]");
		htmp->GetYaxis()->SetLabelFont(43);
		htmp->GetYaxis()->SetLabelSize(20);
		htmp->GetYaxis()->SetTitleFont(43);
		htmp->GetYaxis()->SetTitleSize(20);

		htmp->GetXaxis()->CenterTitle();
		htmp->GetYaxis()->CenterTitle();
		htmp->GetYaxis()->SetTitleOffset(1.0);

		htmp->GetZaxis()->SetLabelFont(43);
		htmp->GetZaxis()->SetLabelSize(20);
		htmp->GetZaxis()->SetTitleFont(43);
		htmp->GetZaxis()->SetTitleSize(20);

		if ( bLOG ){
			h2d[it]->SetMinimum(min);
			h2d[it]->SetMaximum(max);
		}else{
			h2d[it]->SetMaximum(1.05*max);
			h2d[it]->SetMinimum(0);
		}
		htmp = (TH1D*)h2d[it];
		//SetHistoStyle("","","",0.04,0.04);

		TLatex *tex = new TLatex(10.,10.5,"T [GeV]");
		tex->SetTextFont(43);
		tex->SetTextSize(20);
		//tex->SetTextAlign(21);
		tex->Draw();

		/*
		float time = h1_time->GetBinContent(it);
		TLatex *tex = new TLatex(0,6.5,Form("Q_{s}=%g GeV, m=%g GeV, t=%7.5f fm/c",Qs,m,time));
		tex->SetTextFont(63);
		tex->SetTextSize(20);
		tex->SetTextAlign(21);
		tex->Draw();
		*/

		h2d[it]->Draw("colz same");
		if ( bLOG ){
			c1->SaveAs(Form("gif/gif_run00/time_log_%05d.gif",it));
		}else{
			c1->SaveAs(Form("gif/gif_run00/time_linear_%05d.gif",it));
		}
	}



}
