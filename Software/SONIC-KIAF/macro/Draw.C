void Draw(){

	ifstream fdata;
	fdata.open("unid_vn_fluc.dat");

	float tmp_pt, tmp_v0, tmp_v1, tmp_v2, tmp_v3, tmp_v4, tmp_v5;

	TGraphErrors *gv2 = new TGraphErrors;
	TGraphErrors *gv3 = new TGraphErrors;

	int ii = 0;

	while ( fdata >> tmp_pt >> tmp_v0 >> tmp_v1 >> tmp_v2 >> tmp_v3 >> tmp_v4 >> tmp_v5 ){

		gv2->SetPoint(ii, tmp_pt, tmp_v2);
		gv3->SetPoint(ii, tmp_pt, tmp_v3);

		ii++;
	}

	TFile *infile = new TFile("/Users/shlim/Work/Software/SHINCHON/Software/SONIC-Ubuntu/input/MCGlauber-He3Au-200GeV-b0-2fm.root","read");
	TH2D *h2d = (TH2D*)infile->Get("inited_event27");

	TCanvas *c1 = new TCanvas ("c1","c1",800,400);
	c1->Divide(2,1);

	c1->cd(1);
	h2d->Draw("colz");

	c1->cd(2);
	TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,3,0.2);
	htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	htmp->GetYaxis()->SetTitle("v_{n}");
	gv2->SetLineWidth(5);
	gv2->SetLineColor(1);
	gv2->Draw("C");
	
	gv3->SetLineWidth(5);
	gv3->SetLineColor(4);
	gv3->Draw("C");




}
