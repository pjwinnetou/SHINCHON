#include <fstream>
#include <iostream>

#include <TString.h>
#include <TH2D.h>
#include <TList.h>
#include <TSystemDirectory.h>
#include <TFile.h>

using namespace std;

void make_hist_profile(const int runnum=1, const int count=0){

	//TFile *outfile = new TFile(Form("superSONIC_profile_pp13TeV_event%05d.root",runnum),"recreate");
	TFile *outfile = new TFile(Form("superSONIC_profile_pp13TeV_event%05d.root",count),"recreate");

	bool bfirst = true;
	bool bfirst_ed = true;
	float tmax = 1.0;
	float edmax = 1.0;

	const int ngrid = 100;
	const int xmax = 7.5;

	const int tmax_ind = 500;

	TString path = Form("/gpfs/mnt/gpfs02/phenix/hhj/hhj1/shlim/work/98.mc_study/00.SONIC/job/pp_13TeV_3_grp004/wrk_%04d/snapshot",runnum);
	TSystemDirectory dir("snapshot", path);
	TList *files = dir.GetListOfFiles();

	TString fname;
	TSystemFile *file;

	char buf[300];

	TH2D *h2d_t = NULL;
	TH2D *h2d_ed = NULL;
	TH2D *h2d_ux = NULL;
	TH2D *h2d_uy = NULL;
	TH1D *h1d_time = new TH1D("Time",";index;time [fm/c]",tmax_ind+1,-0.5,tmax_ind+0.5);

	ifstream fdata;

	int index_t = 0, index_u = 0, index_ed = 0;

	if ( files ){

		files->Sort();

		TIter next(files);

		while ( (file = (TSystemFile*)next()) ){

			fname = file->GetName();
			if ( file->IsDirectory() ) continue; 
			if ( !fname.EndsWith("dat") ) continue;

			size_t ind_begin = fname.Index("_");
			size_t ind_end = fname.Index(".dat");

			//cout << ind_begin << " " << ind_end << endl;
			TString str_time = fname(ind_begin+1, ind_end-ind_begin-1);

			float time = str_time.Atof();
			str_time.ReplaceAll(".","_");

			cout << Form("%05d",runnum) << " " << fname << " " << Form("%4.3f",time) << endl;
			//cout << fname.(ind_begin, ind_end-ind_begin) << endl;

			if ( fname.BeginsWith("Tcontour") ){

				h2d_t = new TH2D(Form("T_%d",index_t),";x [fm];y [fm]",ngrid,-xmax,xmax,ngrid,-xmax,xmax);
				float tmp_xx, tmp_yy, tmp_tt;

				TString fname_fullpath = path + "/" + fname; 
				fdata.open(fname_fullpath.Data());

				fdata.getline(buf,100);
				//cout << buf << endl;
				for (int iy=0; iy<ngrid; iy++){
					for (int ix=0; ix<ngrid; ix++){
						fdata >> tmp_xx >> tmp_yy >> tmp_tt;
						h2d_t->SetBinContent(ix+1, iy+1, tmp_tt);
					}
				}//
				fdata.close();

				if ( bfirst ){
					tmax = h2d_t->GetMaximum();
					bfirst = false;
				}

				h2d_t->SetMinimum(0);
				h2d_t->SetMaximum(1.05*tmax);

				outfile->cd();
				h2d_t->Write();

				h1d_time->Fill(index_t, time);

				index_t++;
			}else if ( fname.BeginsWith("Vcontour") ){

				h2d_ux = new TH2D(Form("Vx_%d",index_u),";x [fm];y [fm]",ngrid,-xmax,xmax,ngrid,-xmax,xmax);
				h2d_uy = new TH2D(Form("Vy_%d",index_u),";x [fm];y [fm]",ngrid,-xmax,xmax,ngrid,-xmax,xmax);
				float tmp_xx, tmp_yy, tmp_ux, tmp_uy, tmp_gg;

				TString fname_fullpath = path + "/" + fname; 
				fdata.open(fname_fullpath.Data());

				fdata.getline(buf,100);
				for (int iy=0; iy<ngrid; iy++){
					for (int ix=0; ix<ngrid; ix++){
						fdata >> tmp_xx >> tmp_yy >> tmp_ux >> tmp_uy >> tmp_gg;
						h2d_ux->SetBinContent(ix+1, iy+1, tmp_ux);
						h2d_uy->SetBinContent(ix+1, iy+1, tmp_uy);
					}
				}//
				fdata.close();

				outfile->cd();
				h2d_ux->Write();
				h2d_uy->Write();

				index_u++;
			}else if ( fname.BeginsWith("EDcontour") ){

				h2d_ed = new TH2D(Form("ED_%d",index_ed),";x [fm];y [fm]",ngrid,-xmax,xmax,ngrid,-xmax,xmax);
				float tmp_xx, tmp_yy, tmp_ed;

				TString fname_fullpath = path + "/" + fname; 
				fdata.open(fname_fullpath.Data());

				fdata.getline(buf,100);
				//cout << buf << endl;
				for (int iy=0; iy<ngrid; iy++){
					for (int ix=0; ix<ngrid; ix++){
						fdata >> tmp_xx >> tmp_yy >> tmp_ed;
						h2d_ed->SetBinContent(ix+1, iy+1, tmp_ed);
					}
				}//
				fdata.close();

				if ( bfirst_ed ){
					edmax = h2d_ed->GetMaximum();
					bfirst_ed = false;
				}

				h2d_ed->SetMinimum(0);
				h2d_ed->SetMaximum(1.05*edmax);

				outfile->cd();
				h2d_ed->Write();

				index_ed++;
			}

			if ( index_t>tmax_ind || index_u>tmax_ind || index_ed>tmax_ind ){
				break;
			}

		}//while

	}//files

	outfile->cd();
	h1d_time->Write();
	outfile->Close();

}
