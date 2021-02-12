#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include <fstream>
#include <iostream>

using namespace std;

void convert_tree(const char *ifname="output/default/cent0to5/oscar.dat", const char *ofname="outfile.root"){

	int npart = 0;
	int part_pid[2000];
	float part_eta[2000], part_phi[2000], part_pt[2000];

	TTree *T = new TTree("T","T");
	T->Branch("npart",&npart,"npart/I");
	T->Branch("part_pid",part_pid,"part_pid[npart]/I");
	T->Branch("part_eta",part_eta,"part_eta[npart]/F");
	T->Branch("part_phi",part_phi,"part_phi[npart]/F");
	T->Branch("part_pt",part_pt,"part_pt[npart]/F");

	int    evtnumber;
	int    testnum;
	int    nlist;
	double impactpar;
	int    npartproj;
	int    nparttarg;
	int    npartprojelas;
	int    npartprojinelas;
	int    nparttargelas;
	int    nparttarginelas;
	double junk;

	ifstream fdata;
	fdata.open(ifname);

	char buf[300];
	fdata.getline(buf,300);
	fdata.getline(buf,300);
	fdata.getline(buf,300);

	while ( fdata >> evtnumber >> nlist >> junk >> junk ){

		//cout << evtnumber << " " << impactpar << endl;

		npart = 0;
		for (int ii=0; ii<2000; ii++){
			part_pid[ii] = 0;
			part_eta[ii] = part_phi[ii] = part_pt[ii] = 0.0;
		}

		for (int ipart=0; ipart<nlist; ipart++){

			int		 partindex;
			int    partid;
			float  pv[4];
			float  mass;
			double space[4];

			fdata >> partindex >> partid >> pv[0] >> pv[1] >> pv[2] >> pv[3] >> mass >> space[0] >> space[1] >> space[2] >> space[3];

			//MeV to GeV
			TVector3 vec(pv[0]/1000., pv[1]/1000., pv[2]/1000.);

			//if ( (abs(partid)==211 || abs(partid)==321 || abs(partid)==2212) && vec.Pt()>0 ){
			if ( vec.Pt()>0 ){
				part_pid[npart] = partid;
				part_eta[npart] = vec.Eta();
				part_phi[npart] = vec.Phi();
				part_pt[npart] = vec.Pt();
				npart++;
			}
		}//ipart

		T->Fill();

	}//while

	fdata.close();


	TFile *outfile = new TFile(ofname,"recreate");
	T->Write();
	outfile->Close();

}
