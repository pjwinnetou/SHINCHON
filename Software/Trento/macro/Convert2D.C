void Convert2D(){

	const int nevent = 5000;

	const float xmax = 7.5;
	const int nbin = 100;

	char fname[300];
	char buf[300];

	TH2D *h2d[nevent];

	int npart;
	float b;

	TFile *outfile = new TFile("outfile_pAu_p0.root","recreate");

	TTree *T = new TTree("T","T");
	T->Branch("npart",&npart,"npart/I");
	T->Branch("b",&b,"b/F");

	ifstream fdata;

	for (int ievt=0; ievt<nevent; ievt++){

		h2d[ievt] = new TH2D(Form("h2d_event%d",ievt),"",nbin,-xmax,xmax,nbin,-xmax,xmax);

		sprintf(fname, "../config/pAu_events_p0/%04d.dat", ievt);
		cout << "OPEN: " << fname << endl;

		fdata.open(fname);

		for (int ii=0; ii<8; ii++){
			fdata.getline(buf, 300);
			if ( ii==1 || ii==2 ){
				strtok(buf,"=");
				char *tok = strtok(NULL,"=");

				if ( ii==1 ){
					b = atof(tok);
				}else if ( ii==2 ){
					npart = atoi(tok);
				}
				//cout << tok << endl;
			}
			//cout << buf << endl;
		}

		cout << npart << " " << b << endl;

		T->Fill();

		float val;

		for (int ix=0; ix<nbin; ix++){
			for (int iy=0; iy<nbin; iy++){
				fdata >> val; 

				h2d[ievt]->SetBinContent(ix+1, iy+1, val);
			}//iy
		}//ix

		fdata.close();

	}//ievt

	for (int ievt=0; ievt<nevent; ievt++){
		h2d[ievt]->Write();
	}

	T->Write();

	outfile->Close();

	return;


}
