#include "Style.h"

void CalcRaaV1(){

	gStyle->SetOptStat(0);
	gStyle->SetPalette(55);

	const bool bDRAW = true; 
	const bool bSAVE = bDRAW && true; 

	const int nrun = 1000;
	const float const_hbarc = 197.5; //MeV fm
	const float const_mY = 9.46; //GeV
	//const float const_mY = 4.18; //GeV

	//const int ntime = 7;
	//float f_time[ntime] = {0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0};
	const int ntime = 2;
	float f_time[ntime] = {0.5, 10.0};

	ifstream fdata;

	char buf[500];
	vector<float> T[21];
	vector<float> Gdiss[21];

	float f_tmp[20];

	fdata.open("Gdiss0.dat");
	fdata.getline(buf,500);

	while ( fdata 
			>> f_tmp[0] >> f_tmp[1] >> f_tmp[2] >> f_tmp[3]
			>> f_tmp[4] >> f_tmp[5] >> f_tmp[6] >> f_tmp[7]
			>> f_tmp[8] >> f_tmp[9] >> f_tmp[10] >> f_tmp[11]
			){
		for (int ipt=0; ipt<11; ipt++){
			T[ipt].push_back(f_tmp[0]);
			Gdiss[ipt].push_back(f_tmp[ipt+1]);
		}
	}

	fdata.close();

	fdata.open("Gdiss1.dat");
	fdata.getline(buf,500);

	while ( fdata 
			>> f_tmp[0] >> f_tmp[1] >> f_tmp[2] >> f_tmp[3]
			>> f_tmp[4] >> f_tmp[5] >> f_tmp[6] >> f_tmp[7]
			>> f_tmp[8] >> f_tmp[9] >> f_tmp[10]
			){
		for (int ipt=0; ipt<10; ipt++){
			T[ipt+11].push_back(f_tmp[0]);
			Gdiss[ipt+11].push_back(f_tmp[ipt+1]);
		}
	}

	fdata.close();

	TGraphErrors *gGdiss[21];

	for (int ipt=0; ipt<21; ipt++){
		gGdiss[ipt] = new TGraphErrors;
		for (unsigned int iT=0; iT<T[ipt].size(); iT++){
			gGdiss[ipt]->SetPoint(iT, T[ipt][iT], Gdiss[ipt][iT]);
		}
	}

	TCanvas *c0;
	if ( bDRAW ){
		c0 = new TCanvas("c0","c0",1.2*500,500);
		htmp = (TH1D*)gPad->DrawFrame(170,0,500,200);
		htmp->GetXaxis()->SetTitle("T [MeV]");
		htmp->GetYaxis()->SetTitle("#Gamma_{diss} [MeV]");

		for (int ipt=0; ipt<21; ipt++){
			gGdiss[ipt]->SetLineWidth(3);
			gGdiss[ipt]->SetLineColorAlpha(kGray,0.5);
			if ( ipt==0 ){
				gGdiss[ipt]->SetLineColor(2);
			}else if ( ipt==5 ){
				gGdiss[ipt]->SetLineColor(kGreen+2);
			}else if ( ipt==10 ){
				gGdiss[ipt]->SetLineColor(kBlue);
			}else if ( ipt==15 ){
				gGdiss[ipt]->SetLineColor(kMagenta);
			}
			gGdiss[ipt]->Draw("C");
		}
	}

	//return;

	//Upsilon 
	TF1 *fP = new TF1("fP","[0]*x/sqrt(1-x*x)",0,1);
	fP->SetParameter(0, const_mY);

	TCanvas *c1_;
	if ( bDRAW ){
		c1_ = new TCanvas("c1_","c1_",1.2*500,500);
		SetPadStyle();
		htmp = (TH1D*)gPad->DrawFrame(0,0,1,25);
		SetHistoStyle("#beta","p [GeV]");
		//htmp->GetXaxis()->SetTitle("#beta");
		//htmp->GetYaxis()->SetTitle("p [GeV]");

		fP->SetLineColor(1);
		fP->SetLineWidth(3);
		fP->Draw("same");
	}

	//return;
	
	//Glauber
	TFile *infileGlauber = new TFile("/alice/home/shlim/work/SHINCHON/SHINCHON/Software/SONIC-KIAF/input/MCGlauber-PbPb-5020GeV-b0-18fm.root","read");
	TTree *TGlauber = (TTree*)infileGlauber->Get("lemon");
	int Gnpart, Gncoll;
	TGlauber->SetBranchAddress("npart",&Gnpart);
	TGlauber->SetBranchAddress("ncoll",&Gncoll);

	//QGP T profile
	//TFile *infileHydro = new TFile("SONIC_profile_PbPb5TeV_event00008.root","read");
	TH2D *hTHydro[300];
	TGraphErrors *gTHydro[nrun];
	TF1 *fT2[nrun];
	float timeHydro[nrun][300] = {0.0}; 
	float meanTHydro[nrun][300] = {0.0};
	float errTHydro[nrun][300] = {0.0};
	float freezeT[nrun];
	float Npart[nrun];

	float modFactor[nrun];

	TGraphErrors *gInitTHydro = new TGraphErrors;
	TGraphErrors *gMeanTHydro = new TGraphErrors;
	TGraphErrors *gFreezeTHydro = new TGraphErrors;
	TGraphErrors *gRAAHydro = new TGraphErrors;
	TGraphErrors *gRAA[nrun][ntime];

	for (int irun=0; irun<nrun; irun++){
		TFile *infileHydro = new TFile(Form("/alice/home/shlim/work/SHINCHON/SHINCHON/Software/Diffusion/macro/SONIC_profile_PbPb5TeV_0_18fm_event%05d.root",irun),"read");
		TH1D *htimeHydro = (TH1D*)infileHydro->Get("Time");
		int ntimeHydro = (int)htimeHydro->GetEntries();

		for (int it=0; it<ntimeHydro; it++){
			hTHydro[it] = (TH2D*)infileHydro->Get(Form("T_%d",it));

			timeHydro[irun][it] = htimeHydro->GetBinContent(it+1);

			int tmpCount = 0;
			float tmpT = 0.0;
			for (int ix=0; ix<100; ix++){
				for (int iy=0; iy<100; iy++){
					float dx = hTHydro[it]->GetXaxis()->GetBinCenter(ix+1);
					float dy = hTHydro[it]->GetYaxis()->GetBinCenter(iy+1);

					if ( sqrt(dx*dx+dy*dy)<timeHydro[irun][it] ){
						tmpT += hTHydro[it]->GetBinContent(ix+1, iy+1);
						tmpCount++;
					}

				}//iy
			}//ix

			tmpT /= tmpCount;
			meanTHydro[irun][it] = tmpT*1000;
			errTHydro[irun][it] = 0.02*meanTHydro[irun][it];

			if ( meanTHydro[irun][it]>170.0 ){
				freezeT[irun] = htimeHydro->GetBinContent(it+1);
			}

			//cout << timeHydro[it] << " " << meanT << " " << ncount << endl;

		}//it

		infileHydro->Close();
		delete infileHydro;

		gTHydro[irun] = new TGraphErrors(ntimeHydro, &timeHydro[irun][0], &meanTHydro[irun][0], 0, &errTHydro[irun][0]); 

		fT2[irun] = new TF1("fT2","[0]*pow(0.3/x,[1])",0.5,1.2);
		fT2[irun]->SetParameters(450, 0.333);
		gTHydro[irun]->Fit(fT2[irun],"R0Q");

		TGlauber->GetEntry(irun);
		Npart[irun] = Gnpart;

		cout << irun << " " << Npart[irun] << " " << meanTHydro[irun][0] << " " << freezeT[irun] << endl; 

		gInitTHydro->SetPoint(irun, Npart[irun], meanTHydro[irun][0]);
		gInitTHydro->SetPointError(irun, 0, 0.02*meanTHydro[irun][0]);
		gFreezeTHydro->SetPoint(irun, Npart[irun], freezeT[irun]);
		gFreezeTHydro->SetPointError(irun, 0, 0.02*freezeT[irun]);

		float meanTPre = fT2[irun]->Integral(0.3,0.5)/0.2;

		float meanGdissPre = gGdiss[3]->Eval(meanTPre);
		modFactor[irun] = exp(-(0.2)*meanGdissPre/const_hbarc);

		for (int it=0; it<ntimeHydro; it++){
			if ( meanTHydro[irun][it]<170.0 ) break;

			float meanGdiss = gGdiss[3]->Eval(meanTHydro[irun][it]);
			float dt = timeHydro[irun][it+1] - timeHydro[irun][it];
			modFactor[irun] *= exp(-(dt)*meanGdiss/const_hbarc);
		}

		gRAAHydro->SetPoint(irun, Npart[irun], modFactor[irun]);
		gRAAHydro->SetPointError(irun, 0, 0.02*modFactor[irun]);

		continue;

		/*
		for (int it=0; it<ntime; it++){

			gRAA[irun][it] = new TGraphErrors;
			gRAA[irun][it]->SetLineWidth(3);

			int ncount = 0;
			float intTimePost = 0.0;
			float meanTPre = fT2[irun]->Integral(0.3,0.5)/0.2;
			float meanTPost = 0.0;

			for (int jt=0; jt<ntimeHydro; jt++){
				if ( timeHydro[irun][jt]>f_time[it] ){
					break;
				}
				if ( meanTHydro[irun][jt]>170.0 ){
					ncount++;
					intTimePost += (timeHydro[irun][jt] - timeHydro[irun][jt-1]);
					meanTPost += meanTHydro[irun][jt];
				}
			}//jt

		}//it

		if ( ncount==0 ){
			meanTPost = meanTHydro[irun][0];
		}else{
			meanTPost /= ncount;
		}

		//float meanT = (meanTPre*0.2 + meanTPost*(f_time[it]-0.5))/(f_time[it]-0.3);
		float meanT = (meanTPre*0.2 + meanTPost*(intTimePost-0.5))/(f_time[it]-0.3);

		for (int ipt=0; ipt<21; ipt++){
			float meanGdiss = gGdiss[ipt]->Eval(meanT);
			//float meanRAA = exp(-(f_time[it]-0.3)*meanGdiss/const_hbarc);
			float meanRAA = exp(-(intTimePost+0.5-0.3)*meanGdiss/const_hbarc);

			if ( ipt==0 ){
				//cout << meanGdiss << " " << meanRAA << endl;
			}

			gRAA[irun][it]->SetPoint(ipt, ipt, meanRAA); 
		}
		*/


	}

	TF1 *fT = new TF1("fT","[0]*pow([1]/x,0.333)",0.3,15);
	fT->SetParameters(550, 0.3);

	if ( bDRAW ){
		TCanvas *c1 = new TCanvas("c1","c1",1.2*500,500);
		SetPadStyle();
		htmp = (TH1D*)gPad->DrawFrame(0,100,12,600);
		SetHistoStyle("t [fm/c]","T [MeV]");

		fT->SetLineColor(1);
		fT->SetLineWidth(3);
		fT->Draw("same");

		gTHydro[0]->SetLineWidth(1);
		gTHydro[0]->SetLineColor(2);
		gTHydro[0]->SetLineStyle(7);
		gTHydro[0]->Draw("C");

		fT2[0]->SetRange(0.3,1.2);
		fT2[0]->Draw("same");

		TCanvas *c2 = new TCanvas("c2","c2",1.2*500,500);
		SetPadStyle();
		htmp = (TH1D*)gPad->DrawFrame(0,0,400,600);
		SetHistoStyle("#LTN_{part}#GT","T_{init} [MeV]");

		gInitTHydro->SetLineStyle(7);
		gInitTHydro->SetLineColor(1);
		gInitTHydro->Draw("p");

		TCanvas *c3 = new TCanvas("c3","c3",1.2*500,500);
		SetPadStyle();
		htmp = (TH1D*)gPad->DrawFrame(0,0,400,10);
		SetHistoStyle("#LTN_{part}#GT","#tau_{FO} [fm]");

		gFreezeTHydro->SetLineStyle(7);
		gFreezeTHydro->SetLineColor(1);
		gFreezeTHydro->Draw("p");

		TCanvas *c4 = new TCanvas("c4","c4",1.2*500,500);
		SetPadStyle();
		htmp = (TH1D*)gPad->DrawFrame(0,0,400,1.1);
		SetHistoStyle("#LTN_{part}#GT","R_{AA} p_{T}=3 GeV");

		gRAAHydro->SetLineStyle(7);
		gRAAHydro->SetLineColor(1);
		gRAAHydro->Draw("p");
	}

	if ( bSAVE ){

		TFile *outfile = new TFile("outfile.root","recreate");

		gInitTHydro->Write("gInitTHydro");
		gFreezeTHydro->Write("gFreezeTHydro");
		gRAAHydro->Write("gRAAHydro");

	}

	//Calculate RAA



	return;


	/*

	//return;


		cout << it << " " << ncount << " " << intTimePost << " " << meanTPre << " " << meanTPost << " " << meanT << endl; 


	}//it

	//return;

	TCanvas *c2 = new TCanvas("c2","c2",1.2*500,500);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,20,1);
	SetHistoStyle("p_{T} [GeV]", "R_{AA}");

	gRAA[0]->SetLineColor(kMagenta);
	gRAA[1]->SetLineColor(4);
	gRAA[2]->SetLineColor(kGreen+2);
	gRAA[3]->SetLineColor(2);
	gRAA[4]->SetLineColor(1);
	gRAA[5]->SetLineColor(kOrange+2);
	gRAA[6]->SetLineColor(kCyan+2);

	for (int it=0; it<ntime; it++){
		gRAA[it]->Draw("C");
	}

	if ( bSAVE ){

		c1->cd();
		c1->SaveAs("plots/Draw1_c1.png");

		c1_->cd();
		c1_->SaveAs("plots/Draw1_c1_.png");

		c11->cd();
		c11->SaveAs("plots/Draw1_c11.png");

		c2->cd();
		c2->SaveAs("plots/Draw1_c2.png");

	}
	*/

	//return;



}
