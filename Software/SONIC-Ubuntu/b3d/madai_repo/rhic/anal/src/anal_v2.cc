#ifndef __ANALYZE_V2_CC__
#define __ANALYZE_V2_CC__
#include "analyze.h"

void CAnalyze::CalcV2(){
	bool WRITEDETAILS=parameter::getB(parmap,"WRITEDETAILS",false);
	double error;
	bool WRITE_V2=true;
	FILE* detailsfile;
	string detailsfilename;
		//if(phenix_pion_yield>phenix_pion_yield_min && phenix_pion_yield<phenix_pion_yield_max){
	const int NSPECIES=3;
	int NPT;
	const int STAR_PION_NPTBINS=25;
	double STAR_PION_DELPT=80.0, STAR_PION_PTMIN=0.0,STAR_PION_PTMAX=STAR_PION_NPTBINS*STAR_PION_DELPT;
	const int STAR_KAON_NPTBINS=25;
	double STAR_KAON_DELPT=80.0, STAR_KAON_PTMIN=0.0,STAR_KAON_PTMAX=STAR_KAON_DELPT*STAR_KAON_NPTBINS;
	const int STAR_PROTON_NPTBINS=25;
	double STAR_PROTON_DELPT=80.0, STAR_PROTON_PTMIN=0.0,STAR_PROTON_PTMAX=STAR_PROTON_DELPT*STAR_PROTON_NPTBINS;
	double starv2_pion_ptweight=0.0,starv2_pion_ptweightnorm=0.0;
	double starv2_kaon_ptweight=0.0,starv2_kaon_ptweightnorm=0.0;
	double starv2_proton_ptweight=0.0,starv2_proton_ptweightnorm=0.0;
	parameter::set(parmap,string("B3D_RESONANCES_DECAYS_FILE"),string("progdata/resinfo/decays_pdg_weak.dat"));
	parameter::set(parmap,string("B3D_RESONANCES_INFO_FILE"),string("progdata/resinfo/resonances_pdg_weak.dat"));
	reslist->ReadResInfo();
	
	int ievent=1,ipart,nparts,nevents,ispecies,ID,imoment;
	int ipt;
	NPT=STAR_PION_NPTBINS;
	if(STAR_KAON_NPTBINS>NPT) NPT=STAR_KAON_NPTBINS;
	if(STAR_PROTON_NPTBINS>NPT) NPT=STAR_PROTON_NPTBINS;
	double *starv2[NSPECIES],*starv2norm[NSPECIES];
	for(ispecies=0;ispecies<NSPECIES;ispecies++){
		starv2[ispecies]=new double[NPT];
		starv2norm[ispecies]=new double[NPT];
		for(ipt=0;ipt<NPT;ipt++){
			starv2[ispecies][ipt]=0.0;
			starv2norm[ispecies][ipt]=0.0;
		}
	}
	double pt,cos2phi,sigma;
	string sname;
	CPartH5 *ph5;
	printf("CAnalyze::CalcV2, opening %s\n",h5_infilename.c_str());
	h5file = new H5File(h5_infilename,H5F_ACC_RDONLY);
	nevents=int(h5file->getNumObjs());
	printf("nevents=%d\n",nevents);
	if(nevents>neventsmax) nevents=neventsmax;
	for(ievent=1;ievent<=nevents;ievent++){
		nparts=ReadDataH5(ievent);
		for(ipart=0;ipart<nparts;ipart++){
			ph5=&partH5[ipart];
			pt=sqrt(ph5->px*ph5->px+ph5->py*ph5->py);
			ID=ph5->ID;
			ispecies=-1;
			if(ID==211 || ID==-211 || ID==111) ispecies=0;
			else if(ID==310 || ID==130 || ID==321 || ID==-321) ispecies=1;
			else if(ID==2112 || ID==2212 || ID==-2112 || ID==-2212) ispecies=2;
				//else if(ID==3334 || ID==-3334) ispecies=3;
			if(ispecies==0){
					// if(pt>STAR_PION_PTMIN && pt<STAR_PION_PTMAX){
				ipt=int(lrint(floor(pt/STAR_PION_DELPT)));
				if(ipt<NPT){
					cos2phi=ph5->px*ph5->px-ph5->py*ph5->py;
					cos2phi=cos2phi/(pt*pt);
					starv2norm[ispecies][ipt]+=1.0;
					starv2[ispecies][ipt]+=cos2phi;
				}
			}
			if(ispecies==1){
					// if(pt>STAR_KAON_PTMIN && pt<STAR_KAON_PTMAX){
				ipt=int(lrint(floor(pt/STAR_KAON_DELPT)));
				if(ipt<NPT){
					cos2phi=ph5->px*ph5->px-ph5->py*ph5->py;
					cos2phi=cos2phi/(pt*pt);
					starv2norm[ispecies][ipt]+=1.0;
					starv2[ispecies][ipt]+=cos2phi;
				}
			}
			if(ispecies==2){
					// if(pt>STAR_PROTON_PTMIN && pt<STAR_PROTON_PTMAX){
				ipt=int(lrint(floor(pt/STAR_PROTON_DELPT)));
				if(ipt<NPT){
					cos2phi=ph5->px*ph5->px-ph5->py*ph5->py;
					cos2phi=cos2phi/(pt*pt);
					starv2norm[ispecies][ipt]+=1.0;
					starv2[ispecies][ipt]+=cos2phi;
				}
			}
		}
	}
	delete h5file;
	
	if (bPrintV2){
		string mFileName = output_dataroot + string("/fStPiV2.dat");
		FILE *fStPiV2 = fopen(mFileName.c_str(),"w");
		for (ipt=0;ipt<STAR_PION_NPTBINS;ipt++) {
			pt=(0.5+ipt)*STAR_PION_DELPT;
				//			if (pt>STAR_PION_PTMIN && pt<STAR_PION_PTMAX)
			fprintf(fStPiV2,"%g %g\n",pt/1000.,starv2[0][ipt]/starv2norm[0][ipt]);
		}
		fclose(fStPiV2);
		
		mFileName = output_dataroot + string("/fStKV2.dat");
		FILE *fStKV2 = fopen(mFileName.c_str(),"w");
		for (ipt=0;ipt<STAR_KAON_NPTBINS;ipt++) {
			pt=(0.5+ipt)*STAR_KAON_DELPT;
			if (pt>STAR_KAON_PTMIN && pt<STAR_KAON_PTMAX)
				fprintf(fStKV2,"%g %g\n",pt/1000.,starv2[1][ipt]/starv2norm[1][ipt]);
		}
		fclose(fStKV2);
		
		mFileName = output_dataroot + string("/fStPrV2.dat");
		FILE *fStPrV2 = fopen(mFileName.c_str(),"w");
		for (ipt=0;ipt<STAR_PROTON_NPTBINS;ipt++) {
			pt=(0.5+ipt)*STAR_PROTON_DELPT;
			if (pt>STAR_PROTON_PTMIN && pt<STAR_PROTON_PTMAX)
				fprintf(fStPrV2,"%g %g\n",pt/1000.,starv2[2][ipt]/starv2norm[2][ipt]);
		}
		fclose(fStPrV2);
	}
	for(ispecies=0;ispecies<NSPECIES;ispecies++){
		if(ispecies==0){
			for(ipt=0;ipt<STAR_PION_NPTBINS;ipt++){
				starv2[ispecies][ipt]=starv2[ispecies][ipt]/starv2norm[ispecies][ipt];
				pt=(0.5+ipt)*STAR_PION_DELPT;
				printf("%5.2f %g\n",pt,starv2[0][ipt]);
				if(pt>STAR_PION_PTMIN && pt<STAR_PION_PTMAX){
					starv2_pion_ptweight+=pt*starv2[ispecies][ipt];
					starv2_pion_ptweightnorm+=pt;
				}
			}
			starv2_pion_ptweight=starv2_pion_ptweight/starv2_pion_ptweightnorm;
			sigma=sqrt(0.0001+pow(0.075*starv2_pion_ptweight,2));
			fprintf(output,"double STAR_V2_PION_PTWEIGHT %g %g\n",starv2_pion_ptweight,sigma);
			//printf("double STAR_V2_PION_PTWEIGHT %g %g\n",starv2_pion_ptweight,sigma);
		}
		if(ispecies==1){
			for(ipt=0;ipt<STAR_KAON_NPTBINS;ipt++){
				starv2[ispecies][ipt]=starv2[ispecies][ipt]/starv2norm[ispecies][ipt];
				pt=(0.5+ipt)*STAR_KAON_DELPT;
				printf("%5.2f %g\n",pt,starv2[1][ipt]);
				if(pt>STAR_KAON_PTMIN && pt<STAR_KAON_PTMAX){
					starv2_kaon_ptweight+=pt*starv2[ispecies][ipt];
					starv2_kaon_ptweightnorm+=pt;
				}
			}
			starv2_kaon_ptweight=starv2_kaon_ptweight/starv2_kaon_ptweightnorm;
			sigma=sqrt(0.00001+pow(0.075*starv2_kaon_ptweight,2));
			fprintf(output,"double STAR_V2_KAON_PTWEIGHT %g %g\n",starv2_kaon_ptweight,sigma);
		}
		if(ispecies==2){
			for(ipt=0;ipt<STAR_PROTON_NPTBINS;ipt++){
				starv2[ispecies][ipt]=starv2[ispecies][ipt]/starv2norm[ispecies][ipt];
				pt=(0.5+ipt)*STAR_PROTON_DELPT;
				printf("%5.2f %g\n",pt,starv2[2][ipt]);
				if(pt>STAR_PROTON_PTMIN && pt<STAR_PROTON_PTMAX){
					starv2_proton_ptweight+=pt*starv2[ispecies][ipt];
					starv2_proton_ptweightnorm+=pt;
				}
			}
			starv2_proton_ptweight=starv2_proton_ptweight/starv2_proton_ptweightnorm;
			sigma=sqrt(0.0001+pow(0.075*starv2_proton_ptweight,2));
			fprintf(output,"double STAR_V2_PROTON_PTWEIGHT %g %g\n",starv2_proton_ptweight,sigma);
		}
		if(WRITE_V2){
			if(ispecies==0){
				for(ipt=0;ipt<STAR_PION_NPTBINS;ipt++){
					pt=(0.5+ipt)*STAR_PROTON_DELPT;
					if(pt>STAR_PION_PTMIN && pt<STAR_PION_PTMAX){
						error=100.0/sqrt(2.0*starv2norm[0][ipt]);
						fprintf(output,"double STAR_PION_V2_pt%d %g %g\n",int(lrint(pt)),100.0*starv2[0][ipt],error);
					}
				}
			}
			if(ispecies==1){
				for(ipt=0;ipt<STAR_KAON_NPTBINS;ipt++){
					pt=(0.5+ipt)*STAR_KAON_DELPT;
					if(pt>STAR_KAON_PTMIN && pt<STAR_KAON_PTMAX){
						error=100.0/sqrt(2.0*starv2norm[0][ipt]);
						fprintf(output,"double STAR_KAON_V2_pt%d %g %g\n",int(lrint(pt)),100.0*starv2[1][ipt],error);
					}
				}
			}
			if(ispecies==2){
				for(ipt=0;ipt<STAR_PROTON_NPTBINS;ipt++){
					pt=(0.5+ipt)*STAR_PROTON_DELPT;
					if(pt>STAR_PROTON_PTMIN && pt<STAR_PROTON_PTMAX){
						error=100.0/sqrt(2.0*starv2norm[0][ipt]);
						fprintf(output,"double STAR_PROTON_V2_pt%d %g %g\n",int(lrint(pt)),100.0*starv2[2][ipt],error);
					}
				}
			}
		}
		if(WRITEDETAILS){
			if(ispecies==0) detailsfilename=detailsdirname+"/star_v22_pion.dat";
			if(ispecies==1) detailsfilename=detailsdirname+"/star_v22_kaon.dat";
			if(ispecies==2) detailsfilename=detailsdirname+"/star_v22_ppbar.dat";
			detailsfile=fopen(detailsfilename.c_str(),"w");
			if(ispecies==0){
				fprintf(detailsfile,"# STAR_PION_V2\n");
				for(ipt=0;ipt<STAR_PION_NPTBINS;ipt++){
					pt=(0.5+ipt)*STAR_PROTON_DELPT;
					if(pt>STAR_PION_PTMIN && pt<STAR_PION_PTMAX){
						fprintf(detailsfile,"%6.4f %g\n",pt/1000.0,100.0*starv2[0][ipt]);
					}
				}
			}
			if(ispecies==1){
				fprintf(detailsfile,"# STAR_KAON_V2\n");
				for(ipt=0;ipt<STAR_KAON_NPTBINS;ipt++){
					pt=(0.5+ipt)*STAR_PROTON_DELPT;
					if(pt>STAR_KAON_PTMIN && pt<STAR_KAON_PTMAX){
						fprintf(detailsfile,"%6.4f %g\n",pt/1000.0,100.0*starv2[1][ipt]);
					}
				}
			}
			if(ispecies==2){
				fprintf(detailsfile,"# STAR_PROTON_V2\n");
				for(ipt=0;ipt<STAR_PROTON_NPTBINS;ipt++){
					pt=(0.5+ipt)*STAR_PROTON_DELPT;
					if(pt>STAR_PROTON_PTMIN && pt<STAR_PROTON_PTMAX){
						fprintf(detailsfile,"%6.4f %g\n",pt/1000.0,100.0*starv2[2][ipt]);
					}
				}
			}
			fflush(detailsfile);
			fclose(detailsfile);
		}
	}
}

int mysign(double x)
{
  int res;
  if (x>0)
    res=1;
  if (x<0)
    res=-1;
  if (x==0)
    res=0;
					       
  return res;
}


void CAnalyze::CalcFlucVn(){
  const int nmax=5;
  printf("Calculating vn for lumpy output up to n=%i\n",nmax);
  
  bool WRITEDETAILS=parameter::getB(parmap,"WRITEDETAILS",false);
  double error;
  bool WRITE_V2=true;
  FILE* detailsfile;
  string detailsfilename;
  const int NSPECIES=3;
  int NPT;
  const int STAR_PION_NPTBINS=25;
  double STAR_PION_DELPT=160.0, STAR_PION_PTMIN=0.0,STAR_PION_PTMAX=STAR_PION_NPTBINS*STAR_PION_DELPT;
  const int STAR_KAON_NPTBINS=25;
  double STAR_KAON_DELPT=160.0, STAR_KAON_PTMIN=0.0,STAR_KAON_PTMAX=STAR_KAON_DELPT*STAR_KAON_NPTBINS;
  const int STAR_PROTON_NPTBINS=25;
  double STAR_PROTON_DELPT=160.0, STAR_PROTON_PTMIN=0.0,STAR_PROTON_PTMAX=STAR_PROTON_DELPT*STAR_PROTON_NPTBINS;
  double starv2_pion_ptweight=0.0,starv2_pion_ptweightnorm=0.0;
  double starv2_kaon_ptweight=0.0,starv2_kaon_ptweightnorm=0.0;
  double starv2_proton_ptweight=0.0,starv2_proton_ptweightnorm=0.0;
  parameter::set(parmap,string("B3D_RESONANCES_DECAYS_FILE"),string("progdata/resinfo/decays_pdg_weak.dat"));
  parameter::set(parmap,string("B3D_RESONANCES_INFO_FILE"),string("progdata/resinfo/resonances_pdg_weak.dat"));
  reslist->ReadResInfo();
  
  double DELRAPIDITY;
  DELRAPIDITY=2.0*parameter::getD(parmap,"B3D_ETAMAX",1.0);
  int ievent=1,ipart,nparts,nevents,ispecies,ID,imoment;
  int ipt;
  NPT=STAR_PION_NPTBINS;
  if(STAR_KAON_NPTBINS>NPT) NPT=STAR_KAON_NPTBINS;
  if(STAR_PROTON_NPTBINS>NPT) NPT=STAR_PROTON_NPTBINS;
  double **starvncos[NSPECIES],**starvnsin[NSPECIES],*starvnnorm[NSPECIES];
  double *starvnint[NSPECIES];
  for(ispecies=0;ispecies<NSPECIES;ispecies++){
    starvncos[ispecies]=new double*[NPT];
    starvnsin[ispecies]=new double*[NPT];
    starvnnorm[ispecies]=new double[NPT];
    starvnint[ispecies]=new double[nmax+1];
    for (int n=0;n<=nmax;n++)
      starvnint[ispecies][n]=0.0;
    for(ipt=0;ipt<NPT;ipt++){
      starvnnorm[ispecies][ipt]=0.0;
      starvncos[ispecies][ipt]=new double[nmax+1];
      starvnsin[ispecies][ipt]=new double[nmax+1];
      for (int n=1;n<=nmax;n++)
	{
	  starvncos[ispecies][ipt][n]=0.0;
	  starvnsin[ispecies][ipt][n]=0.0;
	}
    }
  }
  
  double pt,phi,sigma;
  string sname;
  CPartH5 *ph5;
  printf("CAnalyze::CalcFlucVn, opening %s\n",h5_infilename.c_str());
  h5file = new H5File(h5_infilename,H5F_ACC_RDONLY);
  nevents=int(h5file->getNumObjs());
  printf("nevents=%d\n",nevents);
  if(nevents>neventsmax) nevents=neventsmax;
  for(ievent=1;ievent<=nevents;ievent++){
    nparts=ReadDataH5(ievent);
    for(ipart=0;ipart<nparts;ipart++){
      ph5=&partH5[ipart];
      pt=sqrt(ph5->px*ph5->px+ph5->py*ph5->py);
      ID=ph5->ID;
      ispecies=-1;
      if(ID==211 || ID==-211) ispecies=0; //pi^+ or pi^-, but not pi^0
      else if(ID==321 || ID==-321) ispecies=1; //K^+ or K^-
      else if(ID==2212 || ID==-2212) ispecies=2;//p or pbar    
      if(ispecies==0){
	ipt=int(lrint(floor(pt/STAR_PION_DELPT)));
	if(ipt<NPT){
	  phi=atan2(ph5->py,ph5->px);
	  //cos2phi=cos2phi/(pt*pt);
	  starvnnorm[ispecies][ipt]+=1.0;
	  for (int n=1;n<=nmax;n++)
	    {
	      starvncos[ispecies][ipt][n]+=cos(n*phi);
	      starvnsin[ispecies][ipt][n]+=sin(n*phi);
	    }
	}
      }
      if(ispecies==1){
	ipt=int(lrint(floor(pt/STAR_KAON_DELPT)));
	if(ipt<NPT){
	  phi=atan2(ph5->py,ph5->px);
	  //cos2phi=cos2phi/(pt*pt);
	  starvnnorm[ispecies][ipt]+=1.0;
	  for (int n=1;n<=nmax;n++)
	    {
	      starvncos[ispecies][ipt][n]+=cos(n*phi);
	      starvnsin[ispecies][ipt][n]+=sin(n*phi);
	    }
	}
      }
      if(ispecies==2){
	// if(pt>STAR_PROTON_PTMIN && pt<STAR_PROTON_PTMAX){
	ipt=int(lrint(floor(pt/STAR_PROTON_DELPT)));
	if(ipt<NPT){
	  phi=atan2(ph5->py,ph5->px);
	  //cos2phi=cos2phi/(pt*pt);
	  starvnnorm[ispecies][ipt]+=1.0;
	  for (int n=1;n<=nmax;n++)
	    {
	      starvncos[ispecies][ipt][n]+=cos(n*phi);
	      starvnsin[ispecies][ipt][n]+=sin(n*phi);
	    }
	}
      }
    }    
  }
  delete h5file;
  printf("CAnalyze::CalcFlucVn, done with reading file\n");
  
  double psi;
  int sign;
  for(ispecies=0;ispecies<NSPECIES;ispecies++){
    printf("attempting species %i of %i\n",ispecies,NSPECIES);
    if(WRITE_V2){
      if(ispecies==0){
	for(ipt=0;ipt<STAR_PION_NPTBINS;ipt++){
	  pt=(0.5+ipt)*STAR_PION_DELPT;
	  if(pt>STAR_PION_PTMIN && pt<STAR_PION_PTMAX){
	    
	    starvnint[ispecies][0]+=starvnnorm[ispecies][ipt];
	    for (int n=1;n<=nmax;n++)
	      {		
		psi=atan2(starvnsin[ispecies][ipt][n],starvncos[ispecies][ipt][n]);
		
		starvncos[ispecies][ipt][n]=sqrt(starvnsin[ispecies][ipt][n]*starvnsin[ispecies][ipt][n]+starvncos[ispecies][ipt][n]*starvncos[ispecies][ipt][n])/starvnnorm[ispecies][ipt];
		
                
		starvnint[ispecies][n]+=starvncos[ispecies][ipt][n]*starvnnorm[ispecies][ipt];
		error=100.0/sqrt(2.0*starvnnorm[0][ipt]);
		//fprintf(output,"double STAR_PION_V%i_pt%d %g %g psi=%g\n",n,int(lrint(pt)),100.0*starvncos[0][ipt][n],error,psi);
	      }
	  }
	}             
	for (int n=1;n<=nmax;n++)
	  starvnint[ispecies][n]/=starvnint[ispecies][0];
	//printf("finished calculating pions, v0=%f v2=%f\n",starvnint[ispecies][0]*STAR_PION_DELPT/1.e6/nevents,starvnint[ispecies][2]*starvnint[ispecies][0]*STAR_PION_DELPT/1.e6/nevents);
      }
      
      if(ispecies==1){
	for(ipt=0;ipt<STAR_KAON_NPTBINS;ipt++){
	  pt=(0.5+ipt)*STAR_KAON_DELPT;
	  if(pt>STAR_KAON_PTMIN && pt<STAR_KAON_PTMAX){
	    starvnint[ispecies][0]+=starvnnorm[ispecies][ipt];
	    for (int n=1;n<=nmax;n++)
	      {		
		psi=atan2(starvnsin[ispecies][ipt][n],starvncos[ispecies][ipt][n]);
		starvncos[ispecies][ipt][n]=sqrt(starvnsin[ispecies][ipt][n]*starvnsin[ispecies][ipt][n]+starvncos[ispecies][ipt][n]*starvncos[ispecies][ipt][n])/starvnnorm[ispecies][ipt];
		starvnint[ispecies][n]+=starvncos[ispecies][ipt][n]*starvnnorm[ispecies][ipt];
		error=100.0/sqrt(2.0*starvnnorm[0][ipt]);
		//fprintf(output,"double STAR_KAON_V%i_pt%d %g %g psi=%g\n",n,int(lrint(pt)),100.0*starvncos[0][ipt][n],error,psi);
	      }
	  }
	}
	for (int n=1;n<=nmax;n++)
	  {
	    starvnint[ispecies][n]/=starvnint[ispecies][0];
	    //printf("species=%i n=%i intvn=%f\n",ispecies,n,starvnint[ispecies][n]);
	  }
	//printf("finished calculating kaons\n");
      }
      
      if(ispecies==2){
	for(ipt=0;ipt<STAR_PROTON_NPTBINS;ipt++){
	  pt=(0.5+ipt)*STAR_PROTON_DELPT;
	  if(pt>STAR_PROTON_PTMIN && pt<STAR_PROTON_PTMAX){
	    starvnint[ispecies][0]+=starvnnorm[ispecies][ipt];
	    for (int n=1;n<=nmax;n++)
	      {		
		psi=atan2(starvnsin[ispecies][ipt][n],starvncos[ispecies][ipt][n]);
		starvncos[ispecies][ipt][n]=sqrt(starvnsin[ispecies][ipt][n]*starvnsin[ispecies][ipt][n]+starvncos[ispecies][ipt][n]*starvncos[ispecies][ipt][n])/starvnnorm[ispecies][ipt];
		starvnint[ispecies][n]+=starvncos[ispecies][ipt][n]*starvnnorm[ispecies][ipt];
		error=100.0/sqrt(2.0*starvnnorm[0][ipt]);
		//fprintf(output,"double STAR_PROTON_V%i_pt%d %g %g psi=%g\n",n,int(lrint(pt)),100.0*starvncos[0][ipt][n],error,psi);
	      }
	  }
	}
	for (int n=1;n<=nmax;n++)
	  starvnint[ispecies][n]/=starvnint[ispecies][0];
      }	    
    }
    if(WRITEDETAILS){
      if(ispecies==0) detailsfilename=detailsdirname+"/pion_vn_fluc.dat";
      if(ispecies==1) detailsfilename=detailsdirname+"/kaon_vn_fluc.dat";
      if(ispecies==2) detailsfilename=detailsdirname+"/ppbar_vn_fluc.dat";
      detailsfile=fopen(detailsfilename.c_str(),"w");
      double myfac=1.e6/4./M_PI/DELRAPIDITY/STAR_PION_DELPT;
      if(ispecies==0){
	//printf("writing details pion_vn, myfac=%g bins=%g\n",myfac,STAR_PION_NPTBINS);
	fprintf(detailsfile,"# PION_Vn\n");
	for(ipt=0;ipt<STAR_PION_NPTBINS;ipt++){
	  pt=(0.5+ipt)*STAR_PION_DELPT;
	  if(pt>STAR_PION_PTMIN && pt<STAR_PION_PTMAX){
	    if (nmax>4)
	      fprintf(detailsfile,"%6.4f %g\t%g\t%g\t%g\t%g\t%g\n",pt/1000.0,starvnnorm[ispecies][ipt]/nevents/pt*myfac,starvncos[ispecies][ipt][1],starvncos[ispecies][ipt][2],starvncos[ispecies][ipt][3],starvncos[ispecies][ipt][4],starvncos[ispecies][ipt][5]);
	  }
	}
      }
      if(ispecies==1){
	printf("writing details kaon_vn\n");
	fprintf(detailsfile,"# KAON_Vn\n");
	for(ipt=0;ipt<STAR_KAON_NPTBINS;ipt++){
	  pt=(0.5+ipt)*STAR_PROTON_DELPT;
	  if(pt>STAR_KAON_PTMIN && pt<STAR_KAON_PTMAX){
	    if (nmax>4)
	      fprintf(detailsfile,"%6.4f %g\t%g\t%g\t%g\t%g\t%g\n",pt/1000.0,starvnnorm[ispecies][ipt]/nevents/pt*myfac,starvncos[ispecies][ipt][1],starvncos[ispecies][ipt][2],starvncos[ispecies][ipt][3],starvncos[ispecies][ipt][4],starvncos[ispecies][ipt][5]);
	  }
	}
      }
      if(ispecies==2){
	printf("writing details proton_vn\n");
	fprintf(detailsfile,"# PROTON_Vn\n");
	for(ipt=0;ipt<STAR_PROTON_NPTBINS;ipt++){
	  pt=(0.5+ipt)*STAR_PROTON_DELPT;
	  if(pt>STAR_PROTON_PTMIN && pt<STAR_PROTON_PTMAX){
	    if (nmax>4)
	      fprintf(detailsfile,"%6.4f %g\t%g\t%g\t%g\t%g\t%g\n",pt/1000.0,starvnnorm[ispecies][ipt]/nevents/pt*myfac,starvncos[ispecies][ipt][1],starvncos[ispecies][ipt][2],starvncos[ispecies][ipt][3],starvncos[ispecies][ipt][4],starvncos[ispecies][ipt][5]);
	  }
	}
      }
      fflush(detailsfile);
      fclose(detailsfile);
    }
  }

  detailsfilename=detailsdirname+"/star_vint_fluc.dat";
  detailsfile=fopen(detailsfilename.c_str(),"w");
  fprintf(detailsfile,"# pion, kaon and ppbar integrated vn\n");
  for(ispecies=0;ispecies<NSPECIES;ispecies++)
    {
      for (int n=1;n<=nmax;n++)
	{
	  fprintf(detailsfile,"%g\t",starvnint[ispecies][n]);
	}
      fprintf(detailsfile,"\n");
    }
  fflush(detailsfile);
  fclose(detailsfile);

}

void CAnalyze::CalcFlucVnGamma(){
  const int nmax=5;
  printf("Calculating vn for lumpy output up to n=%i\n",nmax);
  
  bool WRITEDETAILS=parameter::getB(parmap,"WRITEDETAILS",false);
  double error;
  bool WRITE_V2=true;
  FILE* detailsfile;
  string detailsfilename;
  const int NSPECIES=5;
  int NPT;
  const int STAR_PION_NPTBINS=25;
  double STAR_PION_DELPT=160.0, STAR_PION_PTMIN=0.0,STAR_PION_PTMAX=STAR_PION_NPTBINS*STAR_PION_DELPT;
  const int STAR_KAON_NPTBINS=25;
  double STAR_KAON_DELPT=160.0, STAR_KAON_PTMIN=0.0,STAR_KAON_PTMAX=STAR_KAON_DELPT*STAR_KAON_NPTBINS;
  const int STAR_PROTON_NPTBINS=25;
  double STAR_PROTON_DELPT=160.0, STAR_PROTON_PTMIN=0.0,STAR_PROTON_PTMAX=STAR_PROTON_DELPT*STAR_PROTON_NPTBINS;
  double starv2_pion_ptweight=0.0,starv2_pion_ptweightnorm=0.0;
  double starv2_kaon_ptweight=0.0,starv2_kaon_ptweightnorm=0.0;
  double starv2_proton_ptweight=0.0,starv2_proton_ptweightnorm=0.0;
  parameter::set(parmap,string("B3D_RESONANCES_DECAYS_FILE"),string("progdata/resinfo/decays_pdg_weak.dat"));
  parameter::set(parmap,string("B3D_RESONANCES_INFO_FILE"),string("progdata/resinfo/resonances_pdg_weak.dat"));
  reslist->ReadResInfo();
  
  double DELRAPIDITY;
  DELRAPIDITY=2.0*parameter::getD(parmap,"B3D_ETAMAX",1.0);
  int ievent=1,ipart,nparts,nevents,ispecies,ID,imoment;
  int ipt;
  NPT=STAR_PION_NPTBINS;
  if(STAR_KAON_NPTBINS>NPT) NPT=STAR_KAON_NPTBINS;
  if(STAR_PROTON_NPTBINS>NPT) NPT=STAR_PROTON_NPTBINS;
  double **starvncos[NSPECIES],**starvnsin[NSPECIES],*starvnnorm[NSPECIES];
  double *starvnint[NSPECIES];
  for(ispecies=0;ispecies<NSPECIES;ispecies++){
    starvncos[ispecies]=new double*[NPT];
    starvnsin[ispecies]=new double*[NPT];
    starvnnorm[ispecies]=new double[NPT];
    starvnint[ispecies]=new double[nmax+1];
    for (int n=0;n<=nmax;n++)
      starvnint[ispecies][n]=0.0;
    for(ipt=0;ipt<NPT;ipt++){
      starvnnorm[ispecies][ipt]=0.0;
      starvncos[ispecies][ipt]=new double[nmax+1];
      starvnsin[ispecies][ipt]=new double[nmax+1];
      for (int n=1;n<=nmax;n++)
	{
	  starvncos[ispecies][ipt][n]=0.0;
	  starvnsin[ispecies][ipt][n]=0.0;
	}
    }
  }
  
  double pt,phi,sigma;
  string sname;
  CPartH5 *ph5;
  printf("CAnalyze::CalcFlucVn, opening %s\n",h5_infilename.c_str());
  h5file = new H5File(h5_infilename,H5F_ACC_RDONLY);
  nevents=int(h5file->getNumObjs());
  printf("nevents=%d\n",nevents);
  if(nevents>neventsmax) nevents=neventsmax;
  for(ievent=1;ievent<=nevents;ievent++){
    nparts=ReadDataH5(ievent);
    for(ipart=0;ipart<nparts;ipart++){
      ph5=&partH5[ipart];
      pt=sqrt(ph5->px*ph5->px+ph5->py*ph5->py);
      ID=ph5->ID;
      ispecies=-1;
      if(ID==211 || ID==-211) ispecies=0; //pi^+ or pi^-, but not pi^0
      else if(ID==321 || ID==-321) ispecies=1; //K^+ or K^-
      else if(ID==2212 || ID==-2212) ispecies=2;//p or pbar 
      else if(ID==22) ispecies=3; //photons
      if(ispecies==0){
	ipt=int(lrint(floor(pt/STAR_PION_DELPT)));
	if(ipt<NPT){
	  phi=atan2(ph5->py,ph5->px);
	  //cos2phi=cos2phi/(pt*pt);
	  starvnnorm[ispecies][ipt]+=1.0;
	  for (int n=1;n<=nmax;n++)
	    {
	      starvncos[ispecies][ipt][n]+=cos(n*phi);
	      starvnsin[ispecies][ipt][n]+=sin(n*phi);
	    }
	}
      }
      if(ispecies==1){
	ipt=int(lrint(floor(pt/STAR_KAON_DELPT)));
	if(ipt<NPT){
	  phi=atan2(ph5->py,ph5->px);
	  //cos2phi=cos2phi/(pt*pt);
	  starvnnorm[ispecies][ipt]+=1.0;
	  for (int n=1;n<=nmax;n++)
	    {
	      starvncos[ispecies][ipt][n]+=cos(n*phi);
	      starvnsin[ispecies][ipt][n]+=sin(n*phi);
	    }
	}
      }
      if(ispecies==2){
	// if(pt>STAR_PROTON_PTMIN && pt<STAR_PROTON_PTMAX){
	ipt=int(lrint(floor(pt/STAR_PROTON_DELPT)));
	if(ipt<NPT){
	  phi=atan2(ph5->py,ph5->px);
	  //cos2phi=cos2phi/(pt*pt);
	  starvnnorm[ispecies][ipt]+=1.0;
	  for (int n=1;n<=nmax;n++)
	    {
	      starvncos[ispecies][ipt][n]+=cos(n*phi);
	      starvnsin[ispecies][ipt][n]+=sin(n*phi);
	    }
	}
      }
      if(ispecies==3){
	ipt=int(lrint(floor(pt/STAR_PION_DELPT)));
	if(ipt<NPT){
	  phi=atan2(ph5->py,ph5->px);
	  starvnnorm[ispecies][ipt]+=1.0;
	  for (int n=1;n<=nmax;n++)
	    {
	      starvncos[ispecies][ipt][n]+=cos(n*phi);
	      starvnsin[ispecies][ipt][n]+=sin(n*phi);
	    }
	}
      }
      //species 4: unidentified (pion, kaon, proton)
      if(ID==211 || ID==-211 || ID==321 || ID==-321 || ID==2212 || ID==-2212 )
	{
	  ispecies=4;
	  ipt=int(lrint(floor(pt/STAR_PION_DELPT)));
	  if(ipt<NPT){
	    phi=atan2(ph5->py,ph5->px);
	    //cos2phi=cos2phi/(pt*pt);
	    starvnnorm[ispecies][ipt]+=1.0;
	    for (int n=1;n<=nmax;n++)
	      {
		starvncos[ispecies][ipt][n]+=cos(n*phi);
		starvnsin[ispecies][ipt][n]+=sin(n*phi);
	      }
	  }
	}
    }    
  }
  delete h5file;
  printf("CAnalyze::CalcFlucVn, done with reading file\n");
  
  double psi;
  int sign;
  for(ispecies=0;ispecies<NSPECIES;ispecies++){
    printf("attempting species %i of %i\n",ispecies,NSPECIES);
    if(WRITE_V2){
      if(ispecies==0){
	for(ipt=0;ipt<STAR_PION_NPTBINS;ipt++){
	  pt=(0.5+ipt)*STAR_PION_DELPT;
	  if(pt>STAR_PION_PTMIN && pt<STAR_PION_PTMAX){
	    
	    starvnint[ispecies][0]+=starvnnorm[ispecies][ipt];
	    for (int n=1;n<=nmax;n++)
	      {		
		psi=atan2(starvnsin[ispecies][ipt][n],starvncos[ispecies][ipt][n]);
		
		starvncos[ispecies][ipt][n]=sqrt(starvnsin[ispecies][ipt][n]*starvnsin[ispecies][ipt][n]+starvncos[ispecies][ipt][n]*starvncos[ispecies][ipt][n])/starvnnorm[ispecies][ipt];
		
                
		starvnint[ispecies][n]+=starvncos[ispecies][ipt][n]*starvnnorm[ispecies][ipt];
		error=100.0/sqrt(2.0*starvnnorm[0][ipt]);
		//fprintf(output,"double STAR_PION_V%i_pt%d %g %g psi=%g\n",n,int(lrint(pt)),100.0*starvncos[0][ipt][n],error,psi);
	      }
	  }
	}             
	for (int n=1;n<=nmax;n++)
	  starvnint[ispecies][n]/=starvnint[ispecies][0];
	//printf("finished calculating pions, v0=%f v2=%f\n",starvnint[ispecies][0]*STAR_PION_DELPT/1.e6/nevents,starvnint[ispecies][2]*starvnint[ispecies][0]*STAR_PION_DELPT/1.e6/nevents);
      }
      
      if(ispecies==1){
	for(ipt=0;ipt<STAR_KAON_NPTBINS;ipt++){
	  pt=(0.5+ipt)*STAR_KAON_DELPT;
	  if(pt>STAR_KAON_PTMIN && pt<STAR_KAON_PTMAX){
	    starvnint[ispecies][0]+=starvnnorm[ispecies][ipt];
	    for (int n=1;n<=nmax;n++)
	      {		
		psi=atan2(starvnsin[ispecies][ipt][n],starvncos[ispecies][ipt][n]);
		starvncos[ispecies][ipt][n]=sqrt(starvnsin[ispecies][ipt][n]*starvnsin[ispecies][ipt][n]+starvncos[ispecies][ipt][n]*starvncos[ispecies][ipt][n])/starvnnorm[ispecies][ipt];
		starvnint[ispecies][n]+=starvncos[ispecies][ipt][n]*starvnnorm[ispecies][ipt];
		error=100.0/sqrt(2.0*starvnnorm[0][ipt]);
		//fprintf(output,"double STAR_KAON_V%i_pt%d %g %g psi=%g\n",n,int(lrint(pt)),100.0*starvncos[0][ipt][n],error,psi);
	      }
	  }
	}
	for (int n=1;n<=nmax;n++)
	  {
	    starvnint[ispecies][n]/=starvnint[ispecies][0];
	    //printf("species=%i n=%i intvn=%f\n",ispecies,n,starvnint[ispecies][n]);
	  }
	//printf("finished calculating kaons\n");
      }
      
      if(ispecies==2){
	for(ipt=0;ipt<STAR_PROTON_NPTBINS;ipt++){
	  pt=(0.5+ipt)*STAR_PROTON_DELPT;
	  if(pt>STAR_PROTON_PTMIN && pt<STAR_PROTON_PTMAX){
	    starvnint[ispecies][0]+=starvnnorm[ispecies][ipt];
	    for (int n=1;n<=nmax;n++)
	      {		
		psi=atan2(starvnsin[ispecies][ipt][n],starvncos[ispecies][ipt][n]);
		starvncos[ispecies][ipt][n]=sqrt(starvnsin[ispecies][ipt][n]*starvnsin[ispecies][ipt][n]+starvncos[ispecies][ipt][n]*starvncos[ispecies][ipt][n])/starvnnorm[ispecies][ipt];
		starvnint[ispecies][n]+=starvncos[ispecies][ipt][n]*starvnnorm[ispecies][ipt];
		error=100.0/sqrt(2.0*starvnnorm[0][ipt]);
		//fprintf(output,"double STAR_PROTON_V%i_pt%d %g %g psi=%g\n",n,int(lrint(pt)),100.0*starvncos[0][ipt][n],error,psi);
	      }
	  }
	}
	for (int n=1;n<=nmax;n++)
	  starvnint[ispecies][n]/=starvnint[ispecies][0];
      }

      if(ispecies==3){
	for(ipt=0;ipt<STAR_PION_NPTBINS;ipt++){
	  pt=(0.5+ipt)*STAR_PION_DELPT;
	  if(pt>STAR_PION_PTMIN && pt<STAR_PION_PTMAX){
	    starvnint[ispecies][0]+=starvnnorm[ispecies][ipt];
	    for (int n=1;n<=nmax;n++)
	      {		
		psi=atan2(starvnsin[ispecies][ipt][n],starvncos[ispecies][ipt][n]);
		starvncos[ispecies][ipt][n]=sqrt(starvnsin[ispecies][ipt][n]*starvnsin[ispecies][ipt][n]+starvncos[ispecies][ipt][n]*starvncos[ispecies][ipt][n])/starvnnorm[ispecies][ipt];
		starvnint[ispecies][n]+=starvncos[ispecies][ipt][n]*starvnnorm[ispecies][ipt];
		error=100.0/sqrt(2.0*starvnnorm[0][ipt]);
	      }
	  }
	}
	for (int n=1;n<=nmax;n++)
	  starvnint[ispecies][n]/=starvnint[ispecies][0];
      }

      
      if(ispecies==4){
	for(ipt=0;ipt<STAR_PION_NPTBINS;ipt++){
	  pt=(0.5+ipt)*STAR_PION_DELPT;
	  if(pt>STAR_PION_PTMIN && pt<STAR_PION_PTMAX){
	    
	    starvnint[ispecies][0]+=starvnnorm[ispecies][ipt];
	    for (int n=1;n<=nmax;n++)
	      {		
		psi=atan2(starvnsin[ispecies][ipt][n],starvncos[ispecies][ipt][n]);
		
		starvncos[ispecies][ipt][n]=sqrt(starvnsin[ispecies][ipt][n]*starvnsin[ispecies][ipt][n]+starvncos[ispecies][ipt][n]*starvncos[ispecies][ipt][n])/starvnnorm[ispecies][ipt];
		
                
		starvnint[ispecies][n]+=starvncos[ispecies][ipt][n]*starvnnorm[ispecies][ipt];
		error=100.0/sqrt(2.0*starvnnorm[0][ipt]);
		//fprintf(output,"double STAR_PION_V%i_pt%d %g %g psi=%g\n",n,int(lrint(pt)),100.0*starvncos[0][ipt][n],error,psi);
	      }
	  }
	}             
	for (int n=1;n<=nmax;n++)
	  starvnint[ispecies][n]/=starvnint[ispecies][0];
      }
      
    }
    if(WRITEDETAILS){
      if(ispecies==0) detailsfilename=detailsdirname+"/pion_vn_fluc.dat";
      if(ispecies==1) detailsfilename=detailsdirname+"/kaon_vn_fluc.dat";
      if(ispecies==2) detailsfilename=detailsdirname+"/ppbar_vn_fluc.dat";
      if(ispecies==3) detailsfilename=detailsdirname+"/gamma_vn_fluc.dat";
      if(ispecies==4) detailsfilename=detailsdirname+"/unid_vn_fluc.dat";
      detailsfile=fopen(detailsfilename.c_str(),"w");
      double myfac=1.e6/4./M_PI/DELRAPIDITY/STAR_PION_DELPT;
      if(ispecies==0){
	//printf("writing details pion_vn, myfac=%g bins=%g\n",myfac,STAR_PION_NPTBINS);
	fprintf(detailsfile,"# PION_Vn\n");
	for(ipt=0;ipt<STAR_PION_NPTBINS;ipt++){
	  pt=(0.5+ipt)*STAR_PION_DELPT;
	  if(pt>STAR_PION_PTMIN && pt<STAR_PION_PTMAX){
	    if (nmax>4)
	      fprintf(detailsfile,"%6.4f %g\t%g\t%g\t%g\t%g\t%g\n",pt/1000.0,starvnnorm[ispecies][ipt]/nevents/pt*myfac,starvncos[ispecies][ipt][1],starvncos[ispecies][ipt][2],starvncos[ispecies][ipt][3],starvncos[ispecies][ipt][4],starvncos[ispecies][ipt][5]);
	  }
	}
      }
      if(ispecies==1){
	printf("writing details kaon_vn\n");
	fprintf(detailsfile,"# KAON_Vn\n");
	for(ipt=0;ipt<STAR_KAON_NPTBINS;ipt++){
	  pt=(0.5+ipt)*STAR_PROTON_DELPT;
	  if(pt>STAR_KAON_PTMIN && pt<STAR_KAON_PTMAX){
	    if (nmax>4)
	      fprintf(detailsfile,"%6.4f %g\t%g\t%g\t%g\t%g\t%g\n",pt/1000.0,starvnnorm[ispecies][ipt]/nevents/pt*myfac,starvncos[ispecies][ipt][1],starvncos[ispecies][ipt][2],starvncos[ispecies][ipt][3],starvncos[ispecies][ipt][4],starvncos[ispecies][ipt][5]);
	  }
	}
      }
      if(ispecies==2){
	printf("writing details proton_vn\n");
	fprintf(detailsfile,"# PROTON_Vn\n");
	for(ipt=0;ipt<STAR_PROTON_NPTBINS;ipt++){
	  pt=(0.5+ipt)*STAR_PROTON_DELPT;
	  if(pt>STAR_PROTON_PTMIN && pt<STAR_PROTON_PTMAX){
	    if (nmax>4)
	      fprintf(detailsfile,"%6.4f %g\t%g\t%g\t%g\t%g\t%g\n",pt/1000.0,starvnnorm[ispecies][ipt]/nevents/pt*myfac,starvncos[ispecies][ipt][1],starvncos[ispecies][ipt][2],starvncos[ispecies][ipt][3],starvncos[ispecies][ipt][4],starvncos[ispecies][ipt][5]);
	  }
	}
      }
      if(ispecies==3){
	printf("writing details gamma_vn\n");
	fprintf(detailsfile,"# GAMMA_Vn\n");
	for(ipt=0;ipt<STAR_PION_NPTBINS;ipt++){
	  pt=(0.5+ipt)*STAR_PION_DELPT;
	  if(pt>STAR_PION_PTMIN && pt<STAR_PION_PTMAX){
	    if (nmax>4)
	      fprintf(detailsfile,"%6.4f %g\t%g\t%g\t%g\t%g\t%g\n",pt/1000.0,starvnnorm[ispecies][ipt]/nevents/pt*myfac,starvncos[ispecies][ipt][1],starvncos[ispecies][ipt][2],starvncos[ispecies][ipt][3],starvncos[ispecies][ipt][4],starvncos[ispecies][ipt][5]);
	  }
	}
      }
      if(ispecies==4){
	fprintf(detailsfile,"# UNIDENTIFIED CHARGED HADRONS Vn\n");
	for(ipt=0;ipt<STAR_PION_NPTBINS;ipt++){
	  pt=(0.5+ipt)*STAR_PION_DELPT;
	  if(pt>STAR_PION_PTMIN && pt<STAR_PION_PTMAX){
	    if (nmax>4)
	      fprintf(detailsfile,"%6.4f %g\t%g\t%g\t%g\t%g\t%g\n",pt/1000.0,starvnnorm[ispecies][ipt]/nevents/pt*myfac,starvncos[ispecies][ipt][1],starvncos[ispecies][ipt][2],starvncos[ispecies][ipt][3],starvncos[ispecies][ipt][4],starvncos[ispecies][ipt][5]);
	  }
	}
      }
      fflush(detailsfile);
      fclose(detailsfile);
    }
  }

  detailsfilename=detailsdirname+"/star_vint_fluc.dat";
  detailsfile=fopen(detailsfilename.c_str(),"w");
  fprintf(detailsfile,"# pion, kaon and ppbar integrated vn\n");
  for(ispecies=0;ispecies<NSPECIES;ispecies++)
    {
      for (int n=1;n<=nmax;n++)
	{
	  fprintf(detailsfile,"%g\t",starvnint[ispecies][n]);
	}
      fprintf(detailsfile,"\n");
    }
  fflush(detailsfile);
  fclose(detailsfile);

}

#endif
