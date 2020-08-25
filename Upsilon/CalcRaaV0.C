#include "Style.h"

void CalcRaaV0(){

	gStyle->SetOptStat(0);
	gStyle->SetPalette(55);

	const bool bSAVE = false; 

	const float const_hbarc = 197.5; //MeV fm
	const float const_mY = 9.46; //GeV
	//const float const_mY = 4.18; //GeV

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

	TCanvas *c0 = new TCanvas("c0","c0",1.2*500,500);
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

	//return;

	//Upsilon 
	TF1 *fP = new TF1("fP","[0]*x/sqrt(1-x*x)",0,1);
	fP->SetParameter(0, const_mY);

	TCanvas *c1_ = new TCanvas("c1_","c1_",1.2*500,500);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,1,25);
	SetHistoStyle("#beta","p [GeV]");
	//htmp->GetXaxis()->SetTitle("#beta");
	//htmp->GetYaxis()->SetTitle("p [GeV]");

	fP->SetLineColor(1);
	fP->SetLineWidth(3);
	fP->Draw("same");

	//return;

	//QGP T profile
	//TFile *infileHydro = new TFile("SONIC_profile_PbPb5TeV_event00008.root","read");
	TFile *infileHydro = new TFile("SONIC_profile_PbPb5TeV_0_18fm_event00000.root","read");
	TH1D *htimeHydro = (TH1D*)infileHydro->Get("Time");
	TH2D *hTHydro[500];
	int ntimeHydro = (int)htimeHydro->GetEntries();
	float timeHydro[500];   
	float meanTHydro[500];
	float errTHydro[500];

	for (int it=0; it<ntimeHydro; it++){
		hTHydro[it] = (TH2D*)infileHydro->Get(Form("T_%d",it));

		timeHydro[it] = htimeHydro->GetBinContent(it+1);

		int ncount = 0;
		float meanT = 0.0;
		for (int ix=0; ix<100; ix++){
			for (int iy=0; iy<100; iy++){
				float dx = hTHydro[it]->GetXaxis()->GetBinCenter(ix+1);
				float dy = hTHydro[it]->GetYaxis()->GetBinCenter(iy+1);

				if ( sqrt(dx*dx+dy*dy)<timeHydro[it] ){
					meanT += hTHydro[it]->GetBinContent(ix+1, iy+1);
					ncount++;
				}

			}//iy
		}//ix

		meanT /= ncount;
		meanTHydro[it] = meanT*1000;
		errTHydro[it] = 0.02*meanTHydro[it];

		//cout << timeHydro[it] << " " << meanT << endl;

	}//it

	TF1 *fT = new TF1("fT","[0]*pow([1]/x,0.333)",0.3,15);
	fT->SetParameters(550, 0.3);

	TCanvas *c1 = new TCanvas("c1","c1",1.2*500,500);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,100,12,600);
	SetHistoStyle("t [fm/c]","T [MeV]");

	fT->SetLineColor(1);
	fT->SetLineWidth(3);
	fT->Draw("same");

	TGraphErrors *gTHydro = new TGraphErrors(ntimeHydro, &timeHydro[0], &meanTHydro[0], 0, &errTHydro[0]); 
	gTHydro->SetLineWidth(1);
	gTHydro->SetLineColor(2);
	gTHydro->SetLineStyle(7);
	gTHydro->Draw("C");

	//return;

	TF1 *fT2 = new TF1("fT2","[0]*pow(0.3/x,[1])",0.5,1.2);
	fT2->SetParameters(450, 0.333);
	gTHydro->Fit(fT2,"R0");
	fT2->SetRange(0.3,1.2);
	fT2->Draw("same");

	TCanvas *c11 = new TCanvas("c11","c11",1.2*500,500);
	SetPadStyle();
	gPad->SetRightMargin(0.12);
	htmp = (TH1D*)gPad->DrawFrame(-10,-10,10,10);
	SetHistoStyle("x [fm]","y [fm]","T [GeV]");
	htmp->GetXaxis()->CenterTitle();
	htmp->GetYaxis()->CenterTitle();
	hTHydro[0]->Draw("colz same");

	//return;

	//Calculate RAA
	const int ntime = 5;
	float f_time[ntime] = {0.5, 1.0, 2.0, 3.0, 5.0};

	TGraphErrors *gRAA[ntime];

	for (int it=0; it<ntime; it++){

		gRAA[it] = new TGraphErrors;
		gRAA[it]->SetLineWidth(3);

		int ncount = 0;
		float meanTPre = fT2->Integral(0.3,0.5)/0.2;
		float meanTPost = 0.0;

		for (int jt=0; jt<ntimeHydro; jt++){
			if ( timeHydro[jt]>f_time[it] ){
				break;
			}
			ncount++;
			meanTPost += meanTHydro[jt];
		}//

		if ( ncount==0 ){
			meanTPost = meanTHydro[0];
		}else{
			meanTPost /= ncount;
		}

		float meanT = (meanTPre*0.2 + meanTPost*(f_time[it]-0.5))/(f_time[it]-0.3);

		cout << it << " " << meanTPre << " " << meanTPost << " " << meanT << endl; 

		for (int ipt=0; ipt<21; ipt++){
			float meanGdiss = gGdiss[ipt]->Eval(meanT);
			float meanRAA = exp(-(f_time[it]-0.3)*meanGdiss/const_hbarc);

			if ( ipt==0 ){
				cout << meanGdiss << " " << meanRAA << endl;
			}

			gRAA[it]->SetPoint(ipt, ipt, meanRAA); 
		}

	}//it

	TCanvas *c2 = new TCanvas("c2","c2",1.2*500,500);
	SetPadStyle();
	htmp = (TH1D*)gPad->DrawFrame(0,0,20,1);
	SetHistoStyle("p_{T} [GeV]", "R_{AA}");

	gRAA[0]->SetLineColor(kMagenta);
	gRAA[1]->SetLineColor(4);
	gRAA[2]->SetLineColor(kGreen+2);
	gRAA[3]->SetLineColor(2);
	gRAA[4]->SetLineColor(1);

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

	//return;



}
