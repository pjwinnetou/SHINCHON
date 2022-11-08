#include <iostream>
#include "../../Style_jaebeom_woFrameLegend.h"
#include "../../SHINCHONLegend_raaCent.C"
#include "../../tdrstyle.C"

using namespace std;

void FinalDrawSONIC(){

	const bool bSAVE = false;

	gStyle->SetOptStat(0);
	setTDRStyle();

	const int nset = 3;
	const int nColor[nset] = {1, 2, 4};

	const char *setname[nset] = {"pO","OO","pPb"};

	TFile *infile_pre[nset];
	TFile *infile_vn[nset];

	TH1D *hdNchdeta[nset];
	TH1D *hvn[nset][3];

	TGraphErrors *gvn[nset][3];

	for (int iset=0; iset<nset; iset++){

		infile_pre[iset] = new TFile(Form("pre_outfile_SONIC_%s8160GeV.root",setname[iset]),"read");
		infile_vn[iset] = new TFile(Form("outfile_vn_SONIC_%s8160GeV.root",setname[iset]),"read");

		hdNchdeta[iset] = (TH1D*)infile_pre[iset]->Get("hdNchdeta1");
		hdNchdeta[iset]->SetLineColor(nColor[iset]);
		hdNchdeta[iset]->SetLineWidth(2);

		for (int ii=0; ii<3; ii++){
			hvn[iset][ii] = (TH1D*)infile_vn[iset]->Get(Form("hv%d",ii+2));
			hvn[iset][ii]->SetLineColor(nColor[iset]);
			hvn[iset][ii]->SetMarkerColor(nColor[iset]);
			hvn[iset][ii]->SetMarkerStyle(24);

			gvn[iset][ii] = new TGraphErrors;
			gvn[iset][ii]->SetLineColorAlpha(nColor[iset],0.3);
			gvn[iset][ii]->SetFillColorAlpha(nColor[iset],0.3);

			for (int ip=0; ip<hvn[iset][ii]->GetNbinsX()-1; ip++){
				float xx = hvn[iset][ii]->GetBinCenter(ip+1);
				float yy = hvn[iset][ii]->GetBinContent(ip+1);
				float xx_err = 0.5*hvn[iset][ii]->GetBinWidth(ip+1);
				float yy_err = hvn[iset][ii]->GetBinError(ip+1);

				gvn[iset][ii]->SetPoint(ip, xx, yy);
				gvn[iset][ii]->SetPointError(ip, xx_err, yy_err);
			}
		}



	}

	TCanvas *c1 = new TCanvas("c1","c1",1.1*2*400,2*400);
	c1->Divide(2,2);

	{
		c1->cd(1);
		gPad->SetMargin(0.14,0.03,0.12,0.03);
		gPad->SetLogy();
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,1,150,2*hdNchdeta[0]->GetMaximum());
		htmp->GetYaxis()->SetTitle("N_{event}");
		htmp->GetXaxis()->SetTitle("dN_{ch}/d#eta |_{#eta=0}");
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		hdNchdeta[0]->Draw("same");
		hdNchdeta[1]->Draw("same");
		hdNchdeta[2]->Draw("same");

		TLegend *leg = new TLegend(0.55,0.75,0.9,0.95);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.05);
		leg->AddEntry("","SONIC 8.16 TeV","");
		leg->AddEntry(hdNchdeta[0],"p+O","L");
		leg->AddEntry(hdNchdeta[1],"O+O","L");
		leg->AddEntry(hdNchdeta[2],"p+Pb","L");
		leg->Draw();

	}

	for (int ii=0; ii<3; ii++){

		c1->cd(ii+2);

		gPad->SetMargin(0.14,0.03,0.12,0.03);
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,3,0.15);
		htmp->GetYaxis()->SetTitle(Form("v_{%d}",ii+2));
		htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		htmp->GetYaxis()->SetTitleSize(0.05);
		htmp->GetYaxis()->SetLabelSize(0.04);
		htmp->GetXaxis()->SetTitleSize(0.05);
		htmp->GetXaxis()->SetLabelSize(0.04);

		hvn[0][ii]->Draw("p same");
		hvn[1][ii]->Draw("p same");
		hvn[2][ii]->Draw("p same");

	}

	TCanvas *c2 = new TCanvas("c2","c2",1.2*500,500);

	{
		gPad->SetMargin(0.14,0.03,0.12,0.03);
		gPad->SetLogy();
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,1,150,8*hdNchdeta[0]->GetMaximum());
		htmp->GetYaxis()->SetTitle("N_{event}");
		htmp->GetXaxis()->SetTitle("dN_{ch}/d#eta");
    htmp->GetYaxis()->SetTitleOffset(1.0);
    htmp->GetXaxis()->SetTitleOffset(1.05);
    htmp->GetXaxis()->SetTitleSize(0.055);
    htmp->GetXaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetNdivisions(510);

		hdNchdeta[0]->Draw("same");
		hdNchdeta[1]->Draw("same");
		hdNchdeta[2]->Draw("same");

		TLegend *leg = new TLegend(0.6,0.95-0.055*5,0.9,0.95);
		SetLegendStyle(leg);
		leg->SetTextFont(43);
		leg->SetTextSize(20);
		leg->AddEntry("","MC Glauber+SONIC","h");
		leg->AddEntry("","#sqrt{s_{NN}}=8 TeV","h");
		leg->AddEntry(hdNchdeta[0],Form("p+O, #LTdN_{ch}/d#eta#GT=%4.1f",hdNchdeta[0]->GetMean()),"L");
		leg->AddEntry(hdNchdeta[1],Form("O+O, #LTdN_{ch}/d#eta#GT=%4.1f",hdNchdeta[1]->GetMean()),"L");
		leg->AddEntry(hdNchdeta[2],Form("p+Pb, #LTdN_{ch}/d#eta#GT=%4.1f",hdNchdeta[2]->GetMean()),"L");
		leg->Draw();
	}

	if ( bSAVE ){
		c2->cd();
		c2->SaveAs("SONIC_mult_small.pdf");
	}

	//return;

	TCanvas *c3 = new TCanvas("c3","c3",1.2*500,500);

	{
		gPad->SetMargin(0.14,0.03,0.12,0.03);
		TH1D *htmp = (TH1D*)gPad->DrawFrame(0,0,3,0.15);
		htmp->GetYaxis()->SetTitle(Form("v_{%d}",2));
		htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    htmp->GetYaxis()->SetTitleOffset(1.0);
    htmp->GetXaxis()->SetTitleOffset(1.05);
    htmp->GetXaxis()->SetTitleSize(0.055);
    htmp->GetXaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetLabelSize(0.043);
    htmp->GetYaxis()->SetNdivisions(510);

		//hvn[0][ii]->Draw("p same");
		//hvn[1][ii]->Draw("p same");
		//hvn[2][ii]->Draw("p same");

		gvn[0][0]->Draw("e3");
		gvn[1][0]->Draw("e3");
		gvn[2][0]->Draw("e3");

		TLegend *leg = new TLegend(0.2,0.95-0.055*5,0.6,0.95);
		SetLegendStyle(leg);
		leg->SetTextFont(43);
		leg->SetTextSize(20);
		leg->AddEntry("","MC Glauber+SONIC","h");
		leg->AddEntry("","#sqrt{s_{NN}}=8 TeV","h");
		leg->AddEntry(gvn[0][0],"0-5% p+O","F");
		leg->AddEntry(gvn[1][0],"0-5% O+O","F");
		leg->AddEntry(gvn[2][0],"0-5% p+Pb","F");
		leg->Draw();
	}

	if ( bSAVE ){
		c3->cd();
		c3->SaveAs("SONIC_v2_small.pdf");
	}


}
