#ifndef __ANALYZE_SPECTRA_CC__
#define __ANALYZE_SPECTRA_CC__
#include "analyze.h"
//
void CAnalyze::CalcSpectra(){
	const int NSPECIES=3;
	
	double PHENIX_PTMAX[NSPECIES]={1000.0,1300.0,1600.0},PHENIX_PTMIN[NSPECIES]={200.0,375.0,600.0};
	const int PHENIX_NPTBINS=16,NPHI=12;
	double PHENIX_pcaval[NSPECIES]={0.0},PHENIX_pcapt[NSPECIES]={240.0,375.0,575.0};
	double PHENIX_DELPT=100.0;
	double binnedPHENIX_spectra[NSPECIES][PHENIX_NPTBINS]={0.0};
	double PHENIX_Rx[PHENIX_NPTBINS][NPHI]={0.0},PHENIX_Ry[PHENIX_NPTBINS][NPHI]={0.0},PHENIX_PHENIX_Rz2[PHENIX_NPTBINS][NPHI]={0.0};
	double PHENIX_Rx2[PHENIX_NPTBINS][NPHI]={0.0},PHENIX_Ry2[PHENIX_NPTBINS][NPHI]={0.0},PHENIX_Rz2[PHENIX_NPTBINS][NPHI]={0.0};
	double wR_tot_PHENIX[PHENIX_NPTBINS][NPHI]={0.0},wR_PHENIX[PHENIX_NPTBINS][NPHI]={0.0};
	double Rout_PHENIX[PHENIX_NPTBINS][NPHI],Rlong_PHENIX[PHENIX_NPTBINS][NPHI],Rside_PHENIX[PHENIX_NPTBINS][NPHI];
	
	
	double STAR_PTMAX[NSPECIES]={800.0,700.0,1000.0},STAR_PTMIN[NSPECIES]={200.0,200.0,350.0};
	const int STAR_NPTBINS=20;
	double STAR_pcaval[NSPECIES]={0.0},STAR_pcapt[NSPECIES]={240.0,375.0,575.0};
	double STAR_DELPT=50.0;
	double binnedSTAR_spectra[NSPECIES][STAR_NPTBINS]={0.0};
	double STAR_Rx[STAR_NPTBINS][NPHI]={0.0},STAR_Ry[STAR_NPTBINS][NPHI]={0.0},STAR_STAR_STAR_Rz2[STAR_NPTBINS][NPHI]={0.0};
	double STAR_Rx2[STAR_NPTBINS][NPHI]={0.0},STAR_Ry2[STAR_NPTBINS][NPHI]={0.0},STAR_Rz2[STAR_NPTBINS][NPHI]={0.0};
	double wR_tot_STAR[STAR_NPTBINS][NPHI]={0.0},wR_STAR[STAR_NPTBINS][NPHI]={0.0};
	double Rout_STAR[STAR_NPTBINS][NPHI],Rlong_STAR[STAR_NPTBINS][NPHI],Rside_STAR[STAR_NPTBINS][NPHI];

	
	int ievent=1,ipart,nparts,nevents,ispecies,ID,ipt,iphi,n;
	double yield,pt,meanpt,spread,sigma,sigmapt,sigmaspread,etot=0.0,mass,pz,rapidity,et,phi;
	double sigmspread,DELRAPIDITY=2.0,pionspectra;
	double xval,spectranorm,xsum,degen;
	
	double rout,rlong,rside,tau;
	double w,wtot;
	double alpha,boseweight;
	nsample=parameter::getI(parmap,"B3D_NSAMPLE",1);
	CPartH5 *ph5;
	printf("CAnalyze::CalcSpectra, opening %s, ",h5_infilename.c_str());
	h5file=new H5File(h5_infilename,H5F_ACC_RDONLY);
	
	nevents=int(h5file->getNumObjs());
	if(nevents>neventsmax) nevents=neventsmax;
	printf("nevents=%d\n",nevents);
	double meanpt_pion=0.0,meanpt_kaon=0.0,meanpt_proton=0.0,meanpt_omega=0.0;
	long long int npions=0,nkaons=0,nprotons=0,nomegas=0;
	for(ievent=1;ievent<=nevents;ievent++){
		nparts=ReadDataH5(ievent);
		for(ipart=0;ipart<nparts;ipart++){
			ph5=&partH5[ipart];
			pt=sqrt(ph5->px*ph5->px+ph5->py*ph5->py);
			et=sqrt(ph5->mass*ph5->mass+pt*pt);
			pz=et*sinh(ph5->rapidity-ph5->eta);
			phi=acos(fabs(ph5->px/pt));
			etot+=sqrt(et*et+pz*pz);
			//printf("ID=%d, px=%g py=%g, rapidity=%g, tau=%g, x=%g, y=%g, eta=%g\n",ph5->ID,ph5->px,ph5->py,ph5->rapidity,ph5->tau,ph5->x,ph5->y,ph5->eta);
			ID=ph5->ID;
			ispecies=-1;
			if(ID==211 || ID==-211 || ID==111) ispecies=0;
			else if(ID==310 || ID==130 || ID==321 || ID==-321) ispecies=1;
			else if(ID==2112 || ID==2212 || ID==-2112 || ID==-2212) ispecies=2;
			else if(abs(ID)!=311 && abs(ID)!=3122 && abs(ID)!=22 && abs(ID)!=3222 && abs(ID)!=3312 && abs(ID)!=3112 && abs(ID)!=3322 && abs(ID)!=3334) printf("ID=%d\n",ID);
			//else if(ID==3334 || ID==-3334) ispecies=3;

			if(ispecies==0){
				npions+=1;
				meanpt_pion+=pt;
			}
			if(ispecies==1) {
				nkaons+=1;
				meanpt_kaon+=pt;
			}
			if(ispecies==2){
				nprotons+=1;
				meanpt_proton+=pt;	
			}
			if(ispecies==3){
				nomegas+=1;
				meanpt_omega+=pt;
			}
			/** FOR PHENIX */
			if(ispecies>=0){
				ipt=lrint(floor(pt/PHENIX_DELPT));
				iphi=lrint(floor(2.0*NPHI*phi/PI));
				if(iphi>=NPHI || iphi<0){
					printf("iphi=%d???\n",iphi);
					exit(1);
				}
				if(ipt<PHENIX_NPTBINS){
					binnedPHENIX_spectra[ispecies][ipt]+=1.0;
				}
				if(ispecies==0 && ipt<PHENIX_NPTBINS){
					GetHBTpars(ph5,tau,rout,rside,rlong);
					//printf("rout=%g, rside=%g, rlong=%g, tau=%g\n",rout,rside,rlong,tau);
					tau=ph5->tau;
					wR_tot_PHENIX[ipt][iphi]+=1.0;
					w=exp(-pow(tau-20.0,2)/800.0);
					wR_PHENIX[ipt][iphi]+=w;
					PHENIX_Rx2[ipt][iphi]+=w*rout*rout;
					PHENIX_Ry2[ipt][iphi]+=w*rside*rside;
					PHENIX_Rz2[ipt][iphi]+=w*rlong*rlong;
					PHENIX_Rx[ipt][iphi]+=w*rout;
				}
			}
			/** NOW FOR STAR */
			if(ispecies>=0){
				ipt=lrint(floor(pt/STAR_DELPT));
				iphi=lrint(floor(2.0*NPHI*phi/PI));
				if(iphi>=NPHI || iphi<0){
					printf("iphi=%d???\n",iphi);
					exit(1);
				}
				if(ipt<STAR_NPTBINS){
					binnedSTAR_spectra[ispecies][ipt]+=1.0;
				}
				if(ispecies==0 && ipt<STAR_NPTBINS){
					GetHBTpars(ph5,tau,rout,rside,rlong);
					//printf("rout=%g, rside=%g, rlong=%g, tau=%g\n",rout,rside,rlong,tau);
					tau=ph5->tau;
					wR_tot_STAR[ipt][iphi]+=1.0;
					w=exp(-pow(tau-20.0,2)/800.0);
					wR_STAR[ipt][iphi]+=w;
					STAR_Rx2[ipt][iphi]+=w*rout*rout;
					STAR_Ry2[ipt][iphi]+=w*rside*rside;
					STAR_Rz2[ipt][iphi]+=w*rlong*rlong;
					STAR_Rx[ipt][iphi]+=w*rout;
				}
			}

		}
	}
	delete h5file;
	printf("SPECTRA yields: Npions=%g, Nkaons=%g, Nprotons=%g\n",double(npions)/double(nsample*nevents),double(nkaons)/double(nsample*nevents),double(nprotons)/double(nsample*nevents)); 
	printf("meanpt: pion %g, kaon %g, proton %g\n",meanpt_pion/double(npions),meanpt_kaon/double(nkaons),meanpt_proton/double(nprotons));
	
/** CALC binned spectra, FIRST -- FOR STAR */
	
	// Now Calculating HBT Radii
	for(ipt=0;ipt<PHENIX_NPTBINS;ipt++){
		rout=rlong=rside=w=wtot=0.0;
		for(iphi=0;iphi<NPHI;iphi++){
			PHENIX_Rx2[ipt][iphi]=PHENIX_Rx2[ipt][iphi]/wR_PHENIX[ipt][iphi];
			PHENIX_Ry2[ipt][iphi]=PHENIX_Ry2[ipt][iphi]/wR_PHENIX[ipt][iphi];
			PHENIX_Rz2[ipt][iphi]=PHENIX_Rz2[ipt][iphi]/wR_PHENIX[ipt][iphi];
			PHENIX_Rx[ipt][iphi]=PHENIX_Rx[ipt][iphi]/wR_PHENIX[ipt][iphi];
			Rout_PHENIX[ipt][iphi]=sqrt(PHENIX_Rx2[ipt][iphi]-PHENIX_Rx[ipt][iphi]*PHENIX_Rx[ipt][iphi]);
			Rside_PHENIX[ipt][iphi]=sqrt(PHENIX_Ry2[ipt][iphi]);
			Rlong_PHENIX[ipt][iphi]=sqrt(PHENIX_Rz2[ipt][iphi]);
			rout+=Rout_PHENIX[ipt][iphi]*wR_PHENIX[ipt][iphi];
			rside+=Rside_PHENIX[ipt][iphi]*wR_PHENIX[ipt][iphi];
			rlong+=Rlong_PHENIX[ipt][iphi]*wR_PHENIX[ipt][iphi];
			w+=wR_PHENIX[ipt][iphi];
			wtot+=wR_tot_PHENIX[ipt][iphi];
		}
		rout=rout/w; rside=rside/w; rlong=rlong/w;
	}
	//
	
	degen=3.0;
	for(ipt=0;ipt<PHENIX_NPTBINS;ipt++){
		pt=(0.5+ipt)*PHENIX_DELPT;
		if(ispecies==0){
			binnedPHENIX_spectra[0][ipt]=0.0;
			pionspectra=0.0;
			for(iphi=0;iphi<NPHI;iphi++){
				alpha=pow(2.0*PI,1.5)*HBARC*HBARC*HBARC/(Rout_PHENIX[ipt][iphi]*Rlong_PHENIX[ipt][iphi]*Rside_PHENIX[ipt][iphi]);
				alpha*=wR_PHENIX[ipt][iphi]/(sqrt(pt*pt+139.57*139.57));
				alpha=double(NPHI)*alpha/(2.0*PI*PHENIX_DELPT*pt*DELRAPIDITY*degen*double(nevents*nsample));
				boseweight=0.0;
				for(n=1;n<10;n+=1){
					boseweight+=pow(double(n),-1.5)*pow(alpha,n-1);
				}
				boseweight=(1.0-wR_PHENIX[ipt][iphi]/wR_tot_PHENIX[ipt][iphi])+wR_PHENIX[ipt][iphi]*boseweight/wR_tot_PHENIX[ipt][iphi];
				//printf("ipt=%d, iphi=%d, boseweight=%g\n",ipt,iphi,boseweight);
				//boseweight=1.0;
				pionspectra=(wR_PHENIX[ipt][iphi]*boseweight+(wR_tot_PHENIX[ipt][iphi]-wR_PHENIX[ipt][iphi]));
				binnedPHENIX_spectra[0][ipt]+=pionspectra;
			}
			//printf("pt=%g, alpha=%g, spectra0=%g, spectra=%g, boseweight=%g\n",pt,alpha,binnedPHENIX_spectra[0][ipt],pionspectra,pionspectra/binnedPHENIX_spectra[0][ipt]);
		}
	}
	
	for(ispecies=0;ispecies<NSPECIES;ispecies++){
		if(ispecies==0) degen=3.0;
		if(ispecies==1) degen=4.0;
		if(ispecies==2) degen=4.0;
		if(ispecies==3) degen=2.0;
		spectranorm=0.0;
		meanpt=0.0;
		for(ipt=0;ipt<PHENIX_NPTBINS;ipt++){
			pt=PHENIX_DELPT*(0.5+ipt);
			if(pt>PHENIX_PTMIN[ispecies] && pt<PHENIX_PTMAX[ispecies]){
				spectranorm+=binnedPHENIX_spectra[ispecies][ipt]/double(nevents);
				meanpt+=binnedPHENIX_spectra[ispecies][ipt]*pt/double(nevents);
			}
		}
		sigma=0.05*spectranorm/degen;
		meanpt=meanpt/spectranorm;
		if(ispecies==0){
			fprintf(output,"double PHENIX_SPECTRA_PION_YIELD %g %g\n",spectranorm/degen,sigma);	
			fprintf(output,"double PHENIX_SPECTRA_PION_MEANPT %g %g\n",meanpt,25.0);
		}
		if(ispecies==1){
			fprintf(output,"double PHENIX_SPECTRA_KAON_YIELD %g %g\n",spectranorm/degen,sigma);
			fprintf(output,"double PHENIX_SPECTRA_KAON_MEANPT %g %g\n",meanpt,30.0);
		}
		if(ispecies==2){
			fprintf(output,"double PHENIX_SPECTRA_PROTON_YIELD %g %g\n",spectranorm/degen,sigma);
			fprintf(output,"double PHENIX_SPECTRA_PROTON_MEANPT %g %g\n",meanpt,40.0);
		}
		
		xsum=0.0;
		for(ipt=0;ipt<PHENIX_NPTBINS;ipt++){
			pt=PHENIX_DELPT*(0.5+ipt);
			if(pt>PHENIX_PTMIN[ispecies] && pt<PHENIX_PTMAX[ispecies]){
				xval=binnedPHENIX_spectra[ispecies][ipt]/spectranorm;
				xsum+=xval;
				PHENIX_pcaval[ispecies]+=xval*pt*sinh(pt/PHENIX_pcapt[ispecies]);
				sigma=0.05*xval;
				/*
				if(ispecies==0) fprintf(output,"double PHENIX_SPECTRA_PION_pt%d %g %g\n",int(lrint(pt)),xval,sigma);
				if(ispecies==1) fprintf(output,"double PHENIX_SPECTRA_KAON_pt%d %g %g\n",int(lrint(pt)),xval,sigma);
				if(ispecies==2) fprintf(output,"double PHENIX_SPECTRA_PROTON_pt%d %g %g\n",int(lrint(pt)),xval,sigma);*/
			}
		}
		if(ispecies==0){
			fprintf(output,"double PHENIX_SPECTRA_PION_PCAVAL %g %g\n",PHENIX_pcaval[ispecies],20.0);
		}
		if(ispecies==1) fprintf(output,"double PHENIX_SPECTRA_KAON_PCAVAL %g %g\n",PHENIX_pcaval[ispecies],30.0);
		if(ispecies==2) fprintf(output,"double PHENIX_SPECTRA_PROTON_PCAVAL %g %g\n",PHENIX_pcaval[ispecies],0.2*PHENIX_pcaval[ispecies]);
		//printf("ispecies=%d, xsum=%g\n",ispecies,xsum);
	}
	
	/** NOW, CALC FOR STAR */
	
	// Now Calculating HBT Radii
	for(ipt=0;ipt<STAR_NPTBINS;ipt++){
		rout=rlong=rside=w=wtot=0.0;
		for(iphi=0;iphi<NPHI;iphi++){
			STAR_Rx2[ipt][iphi]=STAR_Rx2[ipt][iphi]/wR_STAR[ipt][iphi];
			STAR_Ry2[ipt][iphi]=STAR_Ry2[ipt][iphi]/wR_STAR[ipt][iphi];
			STAR_Rz2[ipt][iphi]=STAR_Rz2[ipt][iphi]/wR_STAR[ipt][iphi];
			STAR_Rx[ipt][iphi]=STAR_Rx[ipt][iphi]/wR_STAR[ipt][iphi];
			Rout_STAR[ipt][iphi]=sqrt(STAR_Rx2[ipt][iphi]-STAR_Rx[ipt][iphi]*STAR_Rx[ipt][iphi]);
			Rside_STAR[ipt][iphi]=sqrt(STAR_Ry2[ipt][iphi]);
			Rlong_STAR[ipt][iphi]=sqrt(STAR_Rz2[ipt][iphi]);
			rout+=Rout_STAR[ipt][iphi]*wR_STAR[ipt][iphi];
			rside+=Rside_STAR[ipt][iphi]*wR_STAR[ipt][iphi];
			rlong+=Rlong_STAR[ipt][iphi]*wR_STAR[ipt][iphi];
			w+=wR_STAR[ipt][iphi];
			wtot+=wR_tot_STAR[ipt][iphi];
		}
		rout=rout/w; rside=rside/w; rlong=rlong/w;
	}
	//
	
	degen=3.0;
	for(ipt=0;ipt<STAR_NPTBINS;ipt++){
		pt=(0.5+ipt)*STAR_DELPT;
		if(ispecies==0){
			binnedSTAR_spectra[0][ipt]=0.0;
			pionspectra=0.0;
			for(iphi=0;iphi<NPHI;iphi++){
				alpha=pow(2.0*PI,1.5)*HBARC*HBARC*HBARC/(Rout_STAR[ipt][iphi]*Rlong_STAR[ipt][iphi]*Rside_STAR[ipt][iphi]);
				alpha*=wR_STAR[ipt][iphi]/(sqrt(pt*pt+139.57*139.57));
				alpha=double(NPHI)*alpha/(2.0*PI*STAR_DELPT*pt*DELRAPIDITY*degen*double(nevents*nsample));
				boseweight=0.0;
				for(n=1;n<10;n+=1){
					boseweight+=pow(double(n),-1.5)*pow(alpha,n-1);
				}
				boseweight=(1.0-wR_STAR[ipt][iphi]/wR_tot_STAR[ipt][iphi])+wR_STAR[ipt][iphi]*boseweight/wR_tot_STAR[ipt][iphi];
				//printf("ipt=%d, iphi=%d, boseweight=%g\n",ipt,iphi,boseweight);
				//boseweight=1.0;
				pionspectra=(wR_STAR[ipt][iphi]*boseweight+(wR_tot_STAR[ipt][iphi]-wR_STAR[ipt][iphi]));
				binnedSTAR_spectra[0][ipt]+=pionspectra;
			}
			//printf("pt=%g, alpha=%g, spectra0=%g, spectra=%g, boseweight=%g\n",pt,alpha,binnedSTAR_spectra[0][ipt],pionspectra,pionspectra/binnedSTAR_spectra[0][ipt]);
		}
	}
	
	for(ispecies=0;ispecies<NSPECIES;ispecies++){
		if(ispecies==0) degen=3.0;
		if(ispecies==1) degen=4.0;
		if(ispecies==2) degen=4.0;
		if(ispecies==3) degen=2.0;
		spectranorm=0.0;
		meanpt=0.0;
		for(ipt=0;ipt<STAR_NPTBINS;ipt++){
			pt=STAR_DELPT*(0.5+ipt);
			if(pt>STAR_PTMIN[ispecies] && pt<STAR_PTMAX[ispecies]){
				spectranorm+=binnedSTAR_spectra[ispecies][ipt]/double(nevents);
				meanpt+=binnedSTAR_spectra[ispecies][ipt]*pt/double(nevents);
			}
		}
		sigma=0.05*spectranorm/degen;
		meanpt=meanpt/spectranorm;
		if(ispecies==0){
			fprintf(output,"double STAR_SPECTRA_PION_YIELD %g %g\n",spectranorm/degen,sigma);	
			fprintf(output,"double STAR_SPECTRA_PION_MEANPT %g %g\n",meanpt,25.0);
		}
		if(ispecies==1){
			fprintf(output,"double STAR_SPECTRA_KAON_YIELD %g %g\n",spectranorm/degen,sigma);
			fprintf(output,"double STAR_SPECTRA_KAON_MEANPT %g %g\n",meanpt,30.0);
		}
		if(ispecies==2){
			fprintf(output,"double STAR_SPECTRA_PROTON_YIELD %g %g\n",spectranorm/degen,sigma);
			fprintf(output,"double STAR_SPECTRA_PROTON_MEANPT %g %g\n",meanpt,40.0);
		}
		
		xsum=0.0;
		for(ipt=0;ipt<STAR_NPTBINS;ipt++){
			pt=STAR_DELPT*(0.5+ipt);
			if(pt>STAR_PTMIN[ispecies] && pt<STAR_PTMAX[ispecies]){
				xval=binnedSTAR_spectra[ispecies][ipt]/spectranorm;
				xsum+=xval;
				STAR_pcaval[ispecies]+=xval*pt*sinh(pt/STAR_pcapt[ispecies]);
				sigma=0.05*xval;
				/*
				 if(ispecies==0) fprintf(output,"double STAR_SPECTRA_PION_pt%d %g %g\n",int(lrint(pt)),xval,sigma);
				 if(ispecies==1) fprintf(output,"double STAR_SPECTRA_KAON_pt%d %g %g\n",int(lrint(pt)),xval,sigma);
				 if(ispecies==2) fprintf(output,"double STAR_SPECTRA_PROTON_pt%d %g %g\n",int(lrint(pt)),xval,sigma);*/
			}
		}
		if(ispecies==0){
			fprintf(output,"double STAR_SPECTRA_PION_PCAVAL %g %g\n",STAR_pcaval[ispecies],20.0);
		}
		if(ispecies==1) fprintf(output,"double STAR_SPECTRA_KAON_PCAVAL %g %g\n",STAR_pcaval[ispecies],30.0);
		if(ispecies==2) fprintf(output,"double STAR_SPECTRA_PROTON_PCAVAL %g %g\n",STAR_pcaval[ispecies],0.15*STAR_pcaval[ispecies]);
		//printf("ispecies=%d, xsum=%g\n",ispecies,xsum);
	}
	
}


#endif
