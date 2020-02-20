#ifndef __ANALYZE_SPECTRA_CC__
#define __ANALYZE_SPECTRA_CC__
#include "analyze.h"
	//
double CAnalyze::CalcSpectra_PHENIX(){
	string detailsfilename;
	bool WRITEDETAILS=parameter::getB(parmap,"WRITEDETAILS",false);
	bool WRITE_SPECTRA=true;
	FILE *detailsfile;
	const int NSPECIES=2,NPHI=12;
	double PHENIX_PTMAX[NSPECIES]={4000,4000},PHENIX_PTMIN[NSPECIES]={0,0};
	const int PHENIX_NPTBINS=80;
	double pcaval,PHENIX_pcapt[NSPECIES]={240.0,375.0};
	double PHENIX_DELPT=50.0;
	double binnedPHENIX_spectra[NSPECIES][PHENIX_NPTBINS]={0.0};
	double PHENIX_Rx[PHENIX_NPTBINS][NPHI]={0.0},PHENIX_Ry[PHENIX_NPTBINS][NPHI]={0.0},PHENIX_PHENIX_Rz2[PHENIX_NPTBINS][NPHI]={0.0};
	double PHENIX_Rx2[PHENIX_NPTBINS][NPHI]={0.0},PHENIX_Ry2[PHENIX_NPTBINS][NPHI]={0.0},PHENIX_Rz2[PHENIX_NPTBINS][NPHI]={0.0};
	double wR_tot_PHENIX[PHENIX_NPTBINS][NPHI]={0.0},wR_PHENIX[PHENIX_NPTBINS][NPHI]={0.0};
	double Rout_PHENIX[PHENIX_NPTBINS][NPHI],Rlong_PHENIX[PHENIX_NPTBINS][NPHI],Rside_PHENIX[PHENIX_NPTBINS][NPHI];
	
	int ievent=1,ipart,nparts,nevents,ispecies,ID,ipt,iphi,n;
	double yield,pt,meanpt,spread,sigma,sigmapt,sigmaspread,etot=0.0,mass,pz,rapidity,et,phi;
	double DELRAPIDITY,pionspectra;
	double dmult,spectranorm,xsum,degen,phenixpionyield;
	
	double rout,rlong,rside,tau;
	double w,wtot;
	double alpha,boseweight;
	
	parameter::set(parmap,string("B3D_RESONANCES_DECAYS_FILE"),string("progdata/resinfo/decays_pdg_weak.dat"));
	parameter::set(parmap,string("B3D_RESONANCES_INFO_FILE"),string("progdata/resinfo/resonances_pdg_weak.dat"));
	
		//parameter::set(reslist->b3d->parmap,string("B3D_RESONANCES_INFO_FILE"),string("progdata/resinfo/resonances_pdg_weak.dat"));
		//parameter::set(reslist->b3d->parmap,string("B3D_RESONANCES_INFO_FILE"),string("progdata/resinfo/resonances_pdg_weak.dat"));
	
	reslist->ReadResInfo();
	
	nsample=parameter::getI(parmap,"B3D_NSAMPLE",1);
	DELRAPIDITY=2.0*parameter::getD(parmap,"B3D_ETAMAX",1.0);
	CPartH5 *ph5;
	printf("CAnalyze::CalcSpectra, opening %s\n",h5_infilename.c_str());
	h5file=new H5File(h5_infilename,H5F_ACC_RDONLY);
	
	nevents=int(h5file->getNumObjs());
		//printf("nevents=%d\n",nevents);
	if(nevents>neventsmax) nevents=neventsmax;
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
				//printf("ID=%d, px=%g py=%g, rapidity=%g, tau=%g, x=%g, y=%g, eta=%g\n",ph5->ID,
				//ph5->px,ph5->py,ph5->rapidity,ph5->tau,ph5->x,ph5->y,ph5->eta);
			ID=ph5->ID;
			ispecies=-1;
			if(ID==211 || ID==-211) ispecies=0;
			else if(ID==321 || ID==-321) ispecies=1;
			else if(ID==2112 || ID==2212 || ID==-2112 || ID==-2212) ispecies=2;
			else if(ID==3334 || ID==-3334) ispecies=3;
				//			else if(ID!=22 && ID!=-311){
				//				printf("why is particle here? ID=%d\n",ID);
				//			}
				//else if(abs(ID)!=311 &ma& abs(ID)!=3122 && abs(ID)!=22 && abs(ID)!=3222 
				// && abs(ID)!=3312 && abs(ID)!=3112 && abs(ID)!=3322 && abs(ID)!=3334) 
				//printf("ID=%d\n",ID);
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
			
			if(ispecies==0 || ispecies==1){
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
						//					printf("rout=%g, rside=%g, rlong=%g, tau=%g\n",rout,rside,rlong,tau);
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
		}
	}
	delete h5file;
		//printf("SPECTRA yields: Npions=%g, Nkaons=%g, Nprotons=%g\n",double(npions)/double(nsample*nevents),
		//double(nkaons)/double(nsample*nevents),double(nprotons)/double(nsample*nevents)); 
		//printf("meanpt: pion %g, kaon %g, proton %g\n",meanpt_pion/double(npions),
		//meanpt_kaon/double(nkaons),meanpt_proton/double(nprotons));
	
	FILE *fPhPiSp, *fPhKSp, *fPhPrSp, *fSPiSp, *fSKSp, *fSPrSp;
	if (bPrintSpectra){
		string mFileName = output_dataroot + string("/fPhPiSp.dat");
		fPhPiSp = fopen(mFileName.c_str(),"w");
		
		mFileName = output_dataroot + string("/fPhKSp.dat");
		fPhKSp = fopen(mFileName.c_str(),"w");
		
			//		mFileName = output_dataroot + string("/fPhPrSp.dat");
			//		fPhPrSp = fopen(mFileName.c_str(),"w");
		
		mFileName = output_dataroot + string("/fSPiSp.dat");
		fSPiSp = fopen(mFileName.c_str(),"w");
		
		mFileName = output_dataroot + string("/fSKSp.dat");
		fSKSp = fopen(mFileName.c_str(),"w");
		
		mFileName = output_dataroot + string("/fSPrSp.dat");
		fSPrSp = fopen(mFileName.c_str(),"w");
	}
	
	degen=2.0;
	
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
				boseweight=1.0;
				pionspectra=(wR_PHENIX[ipt][iphi]*boseweight+(wR_tot_PHENIX[ipt][iphi]-wR_PHENIX[ipt][iphi]));
				binnedPHENIX_spectra[0][ipt]+=pionspectra;
			}
				//printf("pt=%g, alpha=%g, spectra0=%g, spectra=%g, boseweight=%g\n",pt,alpha,
				//binnedPHENIX_spectra[0][ipt],pionspectra,pionspectra/binnedPHENIX_spectra[0][ipt]);
		}
	}
	
	for(ispecies=0;ispecies<NSPECIES;ispecies++){
		if(ispecies==0) degen=2.0;
		if(ispecies==1) degen=2.0;
		if(ispecies==2) degen=4.0;
		if(ispecies==3) degen=2.0;
		yield=pcaval=meanpt=0.0;
		if(WRITEDETAILS){
			if(ispecies==0) detailsfilename=detailsdirname+"/phenix_spectra_pion.dat";
			if(ispecies==1) detailsfilename=detailsdirname+"/phenix_spectra_kaon.dat";
			detailsfile=fopen(detailsfilename.c_str(),"w");
			if(ispecies==0) fprintf(detailsfile,"# PION SPECTRA \n");
			if(ispecies==1) fprintf(detailsfile,"# KAON SPECTRA \n");
		}
		for(ipt=0;ipt<PHENIX_NPTBINS;ipt++){
			pt=PHENIX_DELPT*(0.5+ipt);
			if(pt>PHENIX_PTMIN[ispecies] && pt<PHENIX_PTMAX[ispecies]){
				dmult=binnedPHENIX_spectra[ispecies][ipt]/double(nevents);
				yield+=dmult;
				meanpt+=dmult*pt;
				pcaval+=dmult*pt*sinh(pt/PHENIX_pcapt[ispecies]);
				sigma=0.1*dmult*pt/1000.0;
				if(WRITE_SPECTRA){
				  if(ispecies==0) fprintf(output,"double PHENIX_SPECTRA_PION_pt%d %g %g and bins=%i\n",int(lrint(pt)),dmult,sigma,PHENIX_NPTBINS);
					if(ispecies==1) fprintf(output,"double PHENIX_SPECTRA_KAON_pt%d %g %g\n",int(lrint(pt)),dmult,sigma);
				}
					//printf("pt=%g, dmult=%g, pcaweight=%g\n",pt,0.5*binnedPHENIX_spectra[ispecies][ipt]/double(nevents),sinh(pt/PHENIX_pcapt[ispecies]));
			}
			if(WRITEDETAILS && (ispecies==0 || ispecies==1)){
				if(pt>PHENIX_PTMIN[ispecies] && pt<PHENIX_PTMAX[ispecies]){
					fprintf(detailsfile,"%g %g\n",pt/1000.,
						binnedPHENIX_spectra[ispecies][ipt]/(degen*PHENIX_DELPT*DELRAPIDITY*2*PI*pt/(1.E6)*double(nevents)));
				}
			}
			
			if (bPrintSpectra) {
				if (ispecies==0 && pt>PHENIX_PTMIN[0])
					fprintf(fPhPiSp,"%g %g\n",pt/1000.,
						binnedPHENIX_spectra[0][ipt]/(degen*PHENIX_DELPT*DELRAPIDITY*2*PI*pt/(1.E6)*double(nevents)));
				if (ispecies==1 && pt>PHENIX_PTMIN[1])
					fprintf(fPhKSp,"%g %g\n",pt/1000.,
						binnedPHENIX_spectra[1][ipt]/(degen*PHENIX_DELPT*DELRAPIDITY*2*PI*pt/(1.E6)*double(nevents)));
					//				if (ispecies==2 && pt>PHENIX_PTMIN[2])
					//					fprintf(fPhPrSp,"%g %g\n",pt/1000.,
					//							binnedPHENIX_spectra[2][ipt]/(degen*PHENIX_DELPT*DELRAPIDITY*2*PI*pt/(1.E6)*double(nevents)));
			}			
		}
		meanpt=meanpt/yield;
		pcaval=pcaval/yield;
		yield=yield/(DELRAPIDITY);
			//printf("yield=%g, meanpt=%g, pcaval=%g\n",yield,meanpt,pcaval);
		

		if(ispecies==0){
			if(qualifier=="cent0to5")	phenix_pion_yield=yield;
			phenixpionyield=yield;
			fprintf(output,"double PHENIX_SPECTRA_PION_YIELD %g %g\n",yield,0.06*yield);	
			fprintf(output,"double PHENIX_SPECTRA_PION_MEANPT %g %g\n",meanpt,0.06*meanpt);
			fprintf(output,"double PHENIX_SPECTRA_PION_PCAVAL %g %g\n",pcaval,0.06*pcaval);
		}
		if(ispecies==1){
			fprintf(output,"double PHENIX_SPECTRA_KAON_YIELD %g %g\n",yield,0.06*yield);
			fprintf(output,"double PHENIX_SPECTRA_KAON_MEANPT %g %g\n",meanpt,0.06*meanpt);
			fprintf(output,"double PHENIX_SPECTRA_KAON_PCAVAL %g %g\n",pcaval,0.06*pcaval);
		}

		if(WRITEDETAILS){
			fflush(detailsfile);
			fclose(detailsfile);
		}
	}


	if (bPrintSpectra) {
		fclose(fPhPiSp);
		fclose(fPhKSp);
			//		fclose(fPhPrSp);
		fclose(fSPiSp);
		fclose(fSKSp);
		fclose(fSPrSp);
	}
	return phenixpionyield;
}

void CAnalyze::CalcSpectra_PHENIXppbar(){
	string detailsfilename;
	bool WRITEDETAILS=parameter::getB(parmap,"WRITEDETAILS",false);
	bool WRITE_SPECTRA=true;
	FILE* detailsfile;
	if(WRITEDETAILS) detailsfile=fopen(detailsfilename.c_str(),"w");
	const int NSPECIES=1;
	double PHENIX_PTMAX=4000.0,PHENIX_PTMIN=0.0;
	const int PHENIX_NPTBINS=80,NPHI=12;
	double PHENIX_pcaval=0.0,PHENIX_pcapt=575.0;
	double PHENIX_DELPT=50.0;
	double binnedPHENIX_spectra[PHENIX_NPTBINS]={0.0};

	int ievent=1,ipart,nparts,nevents,ID,ipt,iphi,n;
	double yield,pt,meanpt,spread,sigma,sigmapt,sigmaspread,etot=0.0,mass,pz,rapidity,et,phi;
	double DELRAPIDITY,pionspectra;
	double dmult,spectranorm,xsum,degen;
	DELRAPIDITY=2.0*parameter::getD(parmap,"B3D_ETAMAX",1.0);

	parameter::set(parmap,string("B3D_RESONANCES_DECAYS_FILE"),string("progdata/resinfo/decays_pdg_weak.dat"));
	parameter::set(parmap,string("B3D_RESONANCES_INFO_FILE"),string("progdata/resinfo/resonances_pdg_weak.dat"));
	reslist->ReadResInfo();
	CResInfo *resinfoptr;
	reslist->GetResInfoptr(3122,resinfoptr);
	resinfoptr->decay=0;
	reslist->GetResInfoptr(-3122,resinfoptr);
	resinfoptr->decay=0;


	nsample=parameter::getI(parmap,"B3D_NSAMPLE",1);
	CPartH5 *ph5;
	printf("CAnalyze::CalcSpectra, opening %s\n",h5_infilename.c_str());
	h5file=new H5File(h5_infilename,H5F_ACC_RDONLY);

	nevents=int(h5file->getNumObjs());
	printf("nevents=%d\n",nevents);
	if(nevents>neventsmax) nevents=neventsmax;
	for(ievent=1;ievent<=nevents;ievent++){
		nparts=ReadDataH5(ievent);
		for(ipart=0;ipart<nparts;ipart++){
			ph5=&partH5[ipart];
			ID=ph5->ID;
			if(abs(ID)==2212 || abs(ID)==2112){
				pt=sqrt(ph5->px*ph5->px+ph5->py*ph5->py);
				ipt=lrint(floor(pt/PHENIX_DELPT));
				if(ipt<PHENIX_NPTBINS){
					binnedPHENIX_spectra[ipt]+=1.0;
				}
			}
		}
	}
	delete h5file;

	degen=4.0;
	spectranorm=0.0;
	yield=0.0;
	meanpt=0.0;
	for(ipt=0;ipt<PHENIX_NPTBINS;ipt++){
		pt=PHENIX_DELPT*(0.5+ipt);
		if(pt>PHENIX_PTMIN && pt<PHENIX_PTMAX){
			dmult=binnedPHENIX_spectra[ipt]/double(nevents);
			yield+=dmult;
			meanpt+=dmult*pt;
			PHENIX_pcaval+=dmult*pt*sinh(pt/PHENIX_pcapt);
			sigma=0.1*dmult*pt/1000.0;
			if(WRITE_SPECTRA)
				fprintf(output,"double PHENIX_SPECTRA_PPBAR_pt%d %g %g\n",int(lrint(pt)),dmult,sigma);
		}
	}
	meanpt=meanpt/yield;
	PHENIX_pcaval+=PHENIX_pcaval/yield;
	yield=0.5*yield/DELRAPIDITY ; // dN/dy of p+pbar in acceptance, factor of 0.5 for neutrons
	fprintf(output,"double PHENIX_SPECTRA_PPBAR_YIELD %g %g\n",yield,0.06*yield);
	fprintf(output,"double PHENIX_SPECTRA_PPBAR_MEANPT %g %g\n",meanpt,0.06*meanpt);
	fprintf(output,"double PHENIX_SPECTRA_PPBAR_PCAVAL %g %g\n",PHENIX_pcaval,0.1*PHENIX_pcaval);
	
	
	if (bPrintSpectra){
		FILE *fPhPrSp;
		string mFileName = output_dataroot + string("/fPhPrSp.dat");
		fPhPrSp = fopen(mFileName.c_str(),"w");
		
		for(ipt=0;ipt<PHENIX_NPTBINS;ipt++){
			pt=PHENIX_DELPT*(0.5+ipt);
			if(pt>PHENIX_PTMIN && pt<PHENIX_PTMAX){			
				fprintf(fPhPrSp,"%g %g\n",pt/1000.,
					binnedPHENIX_spectra[ipt]/(degen*PHENIX_DELPT*DELRAPIDITY*2*PI*pt/(1.E6)*double(nevents)));
				
			}
		}
		fclose(fPhPrSp);
	}
	if(WRITEDETAILS){
		detailsfilename=detailsdirname+"/phenix_spectra_ppbar.dat";
		detailsfile=fopen(detailsfilename.c_str(),"w");
		fprintf(detailsfile,"# PPBAR SPECTRA \n");
		for(ipt=0;ipt<PHENIX_NPTBINS;ipt++){
			pt=PHENIX_DELPT*(0.5+ipt);
			if(pt>PHENIX_PTMIN && pt<PHENIX_PTMAX){
				fprintf(detailsfile,"%g %g\n",pt/1000.0,
					binnedPHENIX_spectra[ipt]/(degen*PHENIX_DELPT*DELRAPIDITY*2*PI*pt/(1.E6)*double(nevents)));
			}
		}
		fflush(detailsfile);
		fclose(detailsfile);
	}
}

#endif
