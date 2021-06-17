void CalOverlap(){

	int nparta, npartb;
	float b;
	float xproj[500], yproj[500];
	float xtarg[500], ytarg[500];
	bool wproj[500], wtarg[500];

	float eccgaus[10];

	const int nset = 1;

	//const int nmaxA[nset] = {1, 1, 16};
	//const int nmaxB[nset] = {208, 16, 16};
	const int nmaxA[nset] = {1};
	const int nmaxB[nset] = {197};

	TFile *infile[nset];

	//infile[0] = new TFile("/alice/data/junleekim/SHINCHON/MCGlauber/MCGlauber-pPb-8160GeV-b0-10fm.root","read");
	//infile[1] = new TFile("/alice/data/junleekim/SHINCHON/MCGlauber/MCGlauber-pO-8160GeV-b0-10fm.root","read");
	//infile[2] = new TFile("/alice/data/junleekim/SHINCHON/MCGlauber/MCGlauber-OO-8160GeV-b0-10fm.root","read");

	infile[0] = new TFile("MCGlauber-pAu-200GeV-b0-15fm-bin100-v2.root","read");


	TH1D *hsT[nset];
	TProfile *hsT_b[nset];

	TH1D *hsTw[nset];
	TProfile *hsTw_b[nset];
	TProfile *hsTw_npart[nset];

	TProfile *he2w_b[nset];
	TProfile *he2w_npart[nset];

	TProfile *he3w_b[nset];
	TProfile *he3w_npart[nset];


	for (int iset=0; iset<nset; iset++){

		hsT[iset] = new TH1D(Form("hsT_set%d",iset),"",100,0,10);
		hsT_b[iset] = new TProfile(Form("hsT_b_set%d",iset),"",50,0,10);

		hsTw[iset] = new TH1D(Form("hsTw_set%d",iset),"",100,0,10);
		hsTw_b[iset] = new TProfile(Form("hsTw_b_set%d",iset),"",50,0,10);
		hsTw_npart[iset] = new TProfile(Form("hsTw_npart_set%d",iset),"",50,0,50);

		he2w_b[iset] = new TProfile(Form("he2w_b_set%d",iset),"",50,0,10);
		he2w_npart[iset] = new TProfile(Form("he2w_npart_set%d",iset),"",50,0,50);

		he3w_b[iset] = new TProfile(Form("he3w_b_set%d",iset),"",50,0,10);
		he3w_npart[iset] = new TProfile(Form("he3w_npart_set%d",iset),"",50,0,50);

		TTree *T = (TTree*)infile[iset]->Get("lemon");
		T->SetBranchAddress("nparta",&nparta);
		T->SetBranchAddress("npartb",&npartb);
		T->SetBranchAddress("b",&b);
		T->SetBranchAddress("xproj",xproj);
		T->SetBranchAddress("yproj",yproj);
		T->SetBranchAddress("wproj",wproj);
		T->SetBranchAddress("xtarg",xtarg);
		T->SetBranchAddress("ytarg",ytarg);
		T->SetBranchAddress("wtarg",wtarg);
		T->SetBranchAddress("eccgaus",eccgaus);

		int nentries = T->GetEntries();

		for (int ien=0; ien<nentries; ien++){

			float meanX = 0, meanY = 0; 
			float meanX2 = 0, meanY2 = 0; 
			float meanXY = 0;
			float sumA = 0, sumB = 0;

			T->GetEntry(ien);

			for (int iA=0; iA<nmaxA[iset]; iA++){
				if ( !wproj[iA] ) continue;

				meanX += xproj[iA];
				meanY += yproj[iA];
				meanX2 += xproj[iA]*xproj[iA];
				meanY2 += yproj[iA]*yproj[iA];
				meanXY += xproj[iA]*yproj[iA];
				sumA++;
			}

			if ( int(sumA)!=nparta ){
				cout << "inconsistent parta, " << sumA << " " << nparta << endl;
			}

			for (int iB=0; iB<nmaxB[iset]; iB++){
				if ( !wtarg[iB] ) continue;

				meanX += xtarg[iB];
				meanY += ytarg[iB];
				meanX2 += xtarg[iB]*xtarg[iB];
				meanY2 += ytarg[iB]*ytarg[iB];
				meanXY += xtarg[iB]*ytarg[iB];
				sumB++;
			}

			if ( int(sumB)!=npartb ){
				cout << "inconsistent partb, " << sumB << " " << npartb << endl;
			}

			meanX /= (sumA + sumB);
			meanY /= (sumA + sumB);
			meanX2 /= (sumA + sumB);
			meanY2 /= (sumA + sumB);
			meanXY /= (sumA + sumB);

			float varX = meanX2 - meanX*meanX;
			float varY = meanY2 - meanY*meanY;
			float varXY = meanXY - meanX*meanY;

			float sTsq = varX*varY - varXY*varXY;
			float sT = sTsq>=0 ? TMath::Pi()*sqrt(sTsq) : -1;

			if ( sT>=0 ){
				hsT[iset]->Fill(sT);
				hsT_b[iset]->Fill(b, sT);
			}

			float meanXw = 0, meanYw = 0;
			float meanX2w = 0, meanY2w = 0;
			float meanXYw = 0;
			float sumw = 0;

			TH2D *h2 = (TH2D*)infile[iset]->Get(Form("inited_event%d",ien));

			float nbinsx = h2->GetNbinsX();
			float nbinsy = h2->GetNbinsY();

			for (int ix=0; ix<nbinsx; ix++){
				for (int iy=0; iy<nbinsy; iy++){
					float xx = h2->GetXaxis()->GetBinCenter(ix+1);
					float yy = h2->GetYaxis()->GetBinCenter(iy+1);

					float ww = h2->GetBinContent(ix+1, iy+1);

					meanXw += xx*ww;
					meanYw += yy*ww;
					meanX2w += xx*xx*ww;
					meanY2w += yy*yy*ww;
					meanXYw += xx*yy*ww;

					sumw += ww;
				}
			}

			meanXw /= sumw;
			meanYw /= sumw;
			meanX2w /= sumw;
			meanY2w /= sumw;
			meanXYw /= sumw;

			float varXw = meanX2w - meanXw*meanXw;
			float varYw = meanY2w - meanYw*meanYw;
			float varXYw = meanXYw - meanXw*meanYw;

			float sTwsq = varXw*varYw - varXYw*varXYw;
			float sTw = sTwsq>=0 ? TMath::Pi()*sqrt(sTwsq) : -1;

			if ( sTw>=0 ){
				hsTw[iset]->Fill(sTw);
				hsTw_b[iset]->Fill(b, sTw);
				hsTw_npart[iset]->Fill(nparta+npartb, sTw);
			}

			he2w_b[iset]->Fill(b, eccgaus[2]);
			he2w_npart[iset]->Fill(nparta+npartb, eccgaus[2]);

			he3w_b[iset]->Fill(b, eccgaus[3]);
			he3w_npart[iset]->Fill(nparta+npartb, eccgaus[3]);


		}//ien

	}//iset

	TFile *outfile = new TFile("outfile_CalOverlap_pAu.root","recreate");

	for (int iset=0; iset<nset; iset++){
		hsT[iset]->Write();
		hsT_b[iset]->Write();

		hsTw[iset]->Write();
		hsTw_b[iset]->Write();
		hsTw_npart[iset]->Write();

		he2w_b[iset]->Write();
		he2w_npart[iset]->Write();

		he3w_b[iset]->Write();
		he3w_npart[iset]->Write();

	}

}
