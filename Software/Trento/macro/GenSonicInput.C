void GenSonicInput(){

	const int nset = 1;
	//const float SF[nset] = {0.763906*1.15, 0.285441};
	//const float SF[nset] = {0.0500457*1.15, 0.0211591};
	const float SF[nset] = {0.0500457*1.15};

	TFile *infile[nset];

	//infile[0] = new TFile("outfile_pPb_p0.root","read");
	//infile[1] = new TFile("outfile_pPb_p1.root","read");
	infile[0] = new TFile("outfile_pAu_p0.root","read");
	//infile[1] = new TFile("outfile_pAu_p1.root","read");

	ofstream fdata;

	char fname[500];

	for (int iset=0; iset<1; iset++){

		for(int ii=0; ii<5000; ii++){

			sprintf(fname,"SONIC_Trento_pAu200GeV_p%d/event%d.dat",iset,ii);
			fdata.open(fname);

			TH2D *h2 = (TH2D*)infile[iset]->Get(Form("h2d_event%d",ii));
			h2->Scale(SF[iset]);

			for (int ix=0; ix<100; ix++){
				for (int iy=0; iy<100; iy++){

					float val = h2->GetBinContent(ix+1, iy+1);
					if ( val<1e-6 ){
						val = 1e-6;
					}
					fdata << val << " ";
				}

				fdata << endl;
			}

			fdata.close();

		}//ii
	}//iset


}
