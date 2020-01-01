#ifndef __ANALYZE_SPECTRA_CC__
#define __ANALYZE_SPECTRA_CC__
#include "analyze.h"
	//
void CAnalyze::Calc3DSpectra(){
	string detailsfilename;
	bool WRITEDETAILS=parameter::getB(parmap,"WRITEDETAILS",false);
	FILE *detailsfile;
	const int NSPECIES=3,NPHI=12;
	const int NPTBINS=100;
	double DELPT=25.0,DELY=1.0;
	double binnedspectra[NSPECIES][NPTBINS]={0.0};
	double Rx[NPTBINS][NPHI]={0.0},Ry[NPTBINS][NPHI]={0.0},Rz2[NPTBINS][NPHI]={0.0};
	double Rx2[NPTBINS][NPHI]={0.0},Ry2[NPTBINS][NPHI]={0.0},Rz2[NPTBINS][NPHI]={0.0};
	double wR_tot[NPTBINS][NPHI]={0.0},wR[NPTBINS][NPHI]={0.0};
	double Rout[NPTBINS][NPHI],Rlong[NPTBINS][NPHI],Rside[NPTBINS][NPHI];
	
	int ievent=1,ipart,nparts,nevents,ispecies,ID,ipt,iphi,n;
	double yield,pt,meanpt,spread,sigma,sigmapt,sigmaspread,etot=0.0,mass,pz,rapidity,et,phi;
	double sigmspread,DELRAPIDITY,pionspectra;
	double dmult,spectranorm,xsum,degen,y;
	
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
			y=asinh(pz/et);
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
			
			if(ispecies==0 || ispecies==1 || ispecies==2){
				ipt=lrint(floor(pt/DELPT));
				iphi=lrint(floor(2.0*NPHI*phi/PI));
				if(iphi>=NPHI || iphi<0){
					printf("iphi=%d???\n",iphi);
					exit(1);
				}
				if(ipt<NPTBINS){
					binnedspectra[ispecies][ipt]+=1.0;
				}
				if(ispecies==0 && ipt<NPTBINS){
					GetHBTpars(ph5,tau,rout,rside,rlong);
						//					printf("rout=%g, rside=%g, rlong=%g, tau=%g\n",rout,rside,rlong,tau);
					tau=ph5->tau;
					wR_tot[ipt][iphi]+=1.0;
					w=exp(-pow(tau-20.0,2)/800.0);
					wR[ipt][iphi]+=w;
					Rx2[ipt][iphi]+=w*rout*rout;
					Ry2[ipt][iphi]+=w*rside*rside;
					Rz2[ipt][iphi]+=w*rlong*rlong;
					Rx[ipt][iphi]+=w*rout;
				}
			}
		}
	}
	delete h5file;
		//printf("SPECTRA yields: Npions=%g, Nkaons=%g, Nprotons=%g\n",double(npions)/double(nsample*nevents),
		//double(nkaons)/double(nsample*nevents),double(nprotons)/double(nsample*nevents)); 
		//printf("meanpt: pion %g, kaon %g, proton %g\n",meanpt_pion/double(npions),
		//meanpt_kaon/double(nkaons),meanpt_proton/double(nprotons));
	
	degen=2.0;
	
	for(ipt=0;ipt<NPTBINS;ipt++){
		pt=(0.5+ipt)*DELPT;
		if(ispecies==0){
			binnedspectra[0][ipt]=0.0;
			pionspectra=0.0;
			for(iphi=0;iphi<NPHI;iphi++){
				alpha=pow(2.0*PI,1.5)*HBARC*HBARC*HBARC/(Rout[ipt][iphi]*Rlong[ipt][iphi]*Rside[ipt][iphi]);
				alpha*=wR[ipt][iphi]/(sqrt(pt*pt+139.57*139.57));
				alpha=double(NPHI)*alpha/(2.0*PI*DELPT*pt*DELRAPIDITY*degen*double(nevents*nsample));
				boseweight=0.0;
				for(n=1;n<10;n+=1){
					boseweight+=pow(double(n),-1.5)*pow(alpha,n-1);
				}
				boseweight=(1.0-wR[ipt][iphi]/wR_tot[ipt][iphi])+wR[ipt][iphi]*boseweight/wR_tot[ipt][iphi];
					//printf("ipt=%d, iphi=%d, boseweight=%g\n",ipt,iphi,boseweight);
				boseweight=1.0;
				pionspectra=(wR[ipt][iphi]*boseweight+(wR_tot[ipt][iphi]-wR[ipt][iphi]));
				binnedspectra[0][ipt]+=pionspectra;
			}
				//printf("pt=%g, alpha=%g, spectra0=%g, spectra=%g, boseweight=%g\n",pt,alpha,
				//binnedspectra[0][ipt],pionspectra,pionspectra/binnedspectra[0][ipt]);
		}
	}
	
	for(ispecies=0;ispecies<NSPECIES;ispecies++){
		if(ispecies==0) degen=2.0;
		if(ispecies==1) degen=2.0;
		if(ispecies==2) degen=4.0;
		if(ispecies==3) degen=2.0;
		yield=meanpt=0.0;
		if(WRITEDETAILS){
			if(ispecies==0) detailsfilename=detailsdirname+"/phenix_spectra_pion.dat";
			if(ispecies==1) detailsfilename=detailsdirname+"/phenix_spectra_kaon.dat";
			detailsfile=fopen(detailsfilename.c_str(),"w");
			if(ispecies==0) fprintf(detailsfile,"# PION SPECTRA \n");
			if(ispecies==1) fprintf(detailsfile,"# KAON SPECTRA \n");
		}
		for(ipt=0;ipt<NPTBINS;ipt++){
			pt=DELPT*(0.5+ipt);
			if(pt>PTMIN[ispecies] && pt<PTMAX[ispecies]){
				dmult=binnedspectra[ispecies][ipt]/double(nevents);
				yield+=dmult;
				meanpt+=dmult*pt;
					//printf("pt=%g, dmult=%g, pcaweight=%g\n",pt,0.5*binnedspectra[ispecies][ipt]/double(nevents),sinh(pt/pcapt[ispecies]));
			}
			if(WRITEDETAILS && (ispecies==0 || ispecies==1)){
				if(pt>PTMIN[ispecies] && pt<PTMAX[ispecies]){
					fprintf(detailsfile,"%g %g\n",pt/1000.,
									binnedspectra[0][ipt]/(degen*DELPT*DELRAPIDITY*2*PI*pt/(1.E6)*double(nevents)));
				}
			}
			
			if (bPrintSpectra) {
				if (ispecies==0 && pt>PTMIN[0])
					fprintf(fPhPiSp,"%g %g\n",pt/1000.,
									binnedspectra[0][ipt]/(degen*DELPT*DELRAPIDITY*2*PI*pt/(1.E6)*double(nevents)));
				if (ispecies==1 && pt>PTMIN[1])
					fprintf(fPhKSp,"%g %g\n",pt/1000.,
									binnedspectra[1][ipt]/(degen*DELPT*DELRAPIDITY*2*PI*pt/(1.E6)*double(nevents)));
					//				if (ispecies==2 && pt>PTMIN[2])
					//					fprintf(fPhPrSp,"%g %g\n",pt/1000.,
					//							binnedspectra[2][ipt]/(degen*DELPT*DELRAPIDITY*2*PI*pt/(1.E6)*double(nevents)));
			}			
		}
		meanpt=meanpt/yield;
		yield=yield/(DELRAPIDITY);
			//printf("yield=%g, meanpt=%g\n",yield,meanpt);
		
		if(WRITEDETAILS) fclose(detailsfile);	
	}
	
	
}

#endif
