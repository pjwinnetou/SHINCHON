#ifndef __ANALYZE_QN_CC__
#define __ANALYZE_QN_CC__
#include "analyze.h"


void CAnalyze::CalcFlucQnEbE(){
  const int nmax=5;
  printf("Using EbE Q-cumulants for Qn's up to n=%i\n",nmax);
  
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
  
  int ievent=1,ipart,nparts,nevents,ispecies,ID,imoment;
  int ipt;
  NPT=STAR_PION_NPTBINS;
  if(STAR_KAON_NPTBINS>NPT) NPT=STAR_KAON_NPTBINS;
  if(STAR_PROTON_NPTBINS>NPT) NPT=STAR_PROTON_NPTBINS;
  double Qnreal[nmax+1],Qnimag[nmax+1],*qnreal[nmax+1],*qnimag[nmax+1];
  double *twopart[nmax+1],normi[nmax+1];
  //so far only for pions
  for (int n=0;n<nmax+1;n++)
    {
      qnreal[n]=new double[NPT];
      qnimag[n]=new double[NPT];
      twopart[n]=new double[NPT];
      Qnreal[n]=0.0;
      Qnimag[n]=0.0;
      normi[n]=0.0;
      for (int ipt=0;ipt<NPT;ipt++)
	{
	  qnreal[n][ipt]=0.0;
	  qnimag[n][ipt]=0.0;
	  twopart[n][ipt]=0.0;
	}
      
    }

  
  
  double pt,phi,sigma;
  string sname;
  CPartH5 *ph5;
  printf("CAnalyze::CalcFlucQn, opening %s\n",h5_infilename.c_str());
  h5file = new H5File(h5_infilename,H5F_ACC_RDONLY);
  nevents=int(h5file->getNumObjs());
  printf("nevents=%d\n",nevents);
  if(nevents>neventsmax) nevents=neventsmax;
  for(ievent=1;ievent<=nevents;ievent++){
    nparts=ReadDataH5(ievent);
    //single event
    for (int n=0;n<nmax;n++)
      {
	Qnreal[n]=0.;
	Qnimag[n]=0.;
	for (int j=0;j<NPT;j++)
	  {
	    qnreal[n][j]=0.;
	    qnimag[n][j]=0.;
	  }
      }


    for(ipart=0;ipart<nparts;ipart++){
      ph5=&partH5[ipart];
      pt=sqrt(ph5->px*ph5->px+ph5->py*ph5->py);
      ID=ph5->ID;
      ipt=int(lrint(floor(pt/STAR_PION_DELPT)));
      phi=atan2(ph5->py,ph5->px);
      ispecies=-1;
      if(ID==211 || ID==-211) ispecies=0; //pi^+ or pi^-, but not pi^0
      else if(ID==321 || ID==-321) ispecies=1; //K^+ or K^-
      else if(ID==2212 || ID==-2212) ispecies=2;//p or pbar    

      
       for (int n=0;n<nmax;n++)
	 {
	   Qnreal[n]+=cos(n*phi);
	   Qnimag[n]+=sin(n*phi);
	   if((ipt<NPT)&&(ispecies==0)) //pions in a certain pt bin
	     {
	       qnreal[n][ipt]+=cos(n*phi);
	       qnimag[n][ipt]+=sin(n*phi);
	     }
	 }

       /*
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

       */
    }
    double r2;
    for (int n=0;n<nmax;n++)
      {
	normi[n]+=(Qnreal[n]*Qnreal[n]+Qnimag[n]*Qnimag[n]-Qnreal[0]);
	for (int j=0;j<NPT;j++)
	  {
	    r2=(qnreal[n][j]*Qnreal[n]+qnimag[n][j]*Qnimag[n]-qnreal[0][j]);
	    twopart[n][j]+=r2;	   
	  }
      }
  }
  delete h5file;
  printf("CAnalyze::CalcFlucQn, done with reading file\n");
  
  double psi;
  int sign;
  for(ispecies=0;ispecies<NSPECIES;ispecies++){
    if(WRITEDETAILS){
      if(ispecies==0) detailsfilename=detailsdirname+"/pion_vn-q-ebe.dat";
      if(ispecies==1) detailsfilename=detailsdirname+"/kaon_vn-q-ebe.dat";
      if(ispecies==2) detailsfilename=detailsdirname+"/ppbar_vn-q-ebe.dat";
      detailsfile=fopen(detailsfilename.c_str(),"w");


      printf("attempting species %i of %i\n",ispecies,NSPECIES);
      if(WRITE_V2){
	if(ispecies==0){
	  fprintf(detailsfile,"# PION_Vn\n");
	  for(ipt=0;ipt<STAR_PION_NPTBINS;ipt++){
	    pt=(0.5+ipt)*STAR_PION_DELPT;
	    if(pt>STAR_PION_PTMIN && pt<STAR_PION_PTMAX){
	      fprintf(detailsfile,"%6.4f\t",pt/1000.0);
	      for (int n=1;n<=nmax;n++)
		{		
		  fprintf(detailsfile,"%g\t",twopart[n][ipt]/twopart[0][ipt]/sqrt(normi[n]/normi[0]));
		}
	      fprintf(detailsfile,"\n");
	    }
	  }           	
	}
      }
      fflush(detailsfile);
      fclose(detailsfile);
    }
  }  
}

void CAnalyze::CalcFlucQn(){
  const int nmax=5;
  printf("Using Q-cumulants for Qn's up to n=%i\n",nmax);
  
  bool WRITEDETAILS=parameter::getB(parmap,"WRITEDETAILS",false);
  double error;
  bool WRITE_V2=true;
  FILE* detailsfile;
  FILE* intvnfile;
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
  
  int ievent=1,ipart,nparts,nevents,ispecies,ID,imoment;
  int ipt;
  NPT=STAR_PION_NPTBINS;
  if(STAR_KAON_NPTBINS>NPT) NPT=STAR_KAON_NPTBINS;
  if(STAR_PROTON_NPTBINS>NPT) NPT=STAR_PROTON_NPTBINS;
  double Qnreal[nmax+1],Qnimag[nmax+1],**qnreal[NSPECIES],**qnimag[NSPECIES];
  
  for (int n=0;n<nmax+1;n++)
    {
      Qnreal[n]=0.0;
      Qnimag[n]=0.0;
    }

  for(ispecies=0;ispecies<NSPECIES;ispecies++){
    qnreal[ispecies]=new double*[nmax+1];
    qnimag[ispecies]=new double*[nmax+1];
    for (int n=0;n<nmax+1;n++)
      {
	qnreal[ispecies][n]=new double[NPT];
	qnimag[ispecies][n]=new double[NPT];
	for (int ipt=0;ipt<NPT;ipt++)
	  {
	    qnreal[ispecies][n][ipt]=0.0;
	    qnimag[ispecies][n][ipt]=0.0;
	  }      
      }
  }
  
  double pt,phi,sigma;
  string sname;
  CPartH5 *ph5;
  printf("CAnalyze::CalcFlucQn, opening %s\n",h5_infilename.c_str());
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
      ipt=int(lrint(floor(pt/STAR_PION_DELPT)));
      phi=atan2(ph5->py,ph5->px);
      ispecies=-1;
      if(ID==211 || ID==-211) ispecies=0; //pi^+ or pi^-, but not pi^0
      else if(ID==321 || ID==-321) ispecies=1; //K^+ or K^-
      else if(ID==2212 || ID==-2212) ispecies=2;//p or pbar    

      
      for (int n=0;n<nmax;n++)
	{
	  Qnreal[n]+=cos(n*phi);
	  Qnimag[n]+=sin(n*phi);
	  if((ipt<NPT)&&(ispecies<NSPECIES)&&(ispecies>=0))
	    {
	      qnreal[ispecies][n][ipt]+=cos(n*phi);
	      qnimag[ispecies][n][ipt]+=sin(n*phi);
	    }
	}

    }    
  }
  delete h5file;
  printf("CAnalyze::CalcFlucQn, done with reading file\n");
  
  double psi;
  int sign;
  if(WRITEDETAILS){
    detailsfilename=detailsdirname+"/int_vn.dat";
    intvnfile=fopen(detailsfilename.c_str(),"w");
    fprintf(intvnfile,"#all vn\n");
    //fprintf(intvnfile,"%i\t %g\n",0,Qnreal[0]/nevents);
    for (int n=1;n<=nmax;n++)
      {
	fprintf(intvnfile,"%i\t",n);
	fprintf(intvnfile,"%g\t",sqrt((Qnreal[n]*Qnreal[n]+Qnimag[n]*Qnimag[n]-Qnreal[0])/(Qnreal[0]*(Qnreal[0]-1))));	  
	fprintf(intvnfile,"\n");
      }
    fflush(intvnfile);
    fclose(intvnfile);      
  }
  for(ispecies=0;ispecies<NSPECIES;ispecies++){
    if(WRITEDETAILS){
      if(ispecies==0) detailsfilename=detailsdirname+"/pion_vn-q.dat";
      if(ispecies==1) detailsfilename=detailsdirname+"/kaon_vn-q.dat";
      if(ispecies==2) detailsfilename=detailsdirname+"/ppbar_vn-q.dat";
      detailsfile=fopen(detailsfilename.c_str(),"w");
      
      printf("attempting species %i of %i\n",ispecies,NSPECIES);
      if(WRITE_V2){

	if(ispecies==0){
	  fprintf(detailsfile,"# PION_Vn\n");
	  for(ipt=0;ipt<STAR_PION_NPTBINS;ipt++){
	    pt=(0.5+ipt)*STAR_PION_DELPT;
	    if(pt>STAR_PION_PTMIN && pt<STAR_PION_PTMAX){
	      fprintf(detailsfile,"%6.4f\t",pt/1000.0);
	      for (int n=1;n<=nmax;n++)
		{		
		  double r2,i2,normi;
		  r2=(qnreal[ispecies][n][ipt]*Qnreal[n]+qnimag[ispecies][n][ipt]*Qnimag[n]-qnreal[ispecies][0][ipt])/(qnreal[ispecies][0][ipt]*(Qnreal[0]-1));
		  //i2=(-qnreal[n][ipt]*Qnimag[n]+qnimag[n][ipt]*Qnreal[n])/(qnreal[0][ipt]*(Qnreal[0]-1));
		  normi=(Qnreal[n]*Qnreal[n]+Qnimag[n]*Qnimag[n]-Qnreal[0])/(Qnreal[0]*(Qnreal[0]-1));
		  fprintf(detailsfile,"%g\t",r2/sqrt(normi));
		}
	      fprintf(detailsfile,"\n");
	    }
	  }           	
	}

	if(ispecies==1){
	  fprintf(detailsfile,"# KAON_Vn\n");
	  for(ipt=0;ipt<STAR_KAON_NPTBINS;ipt++){
	    pt=(0.5+ipt)*STAR_KAON_DELPT;
	    if(pt>STAR_KAON_PTMIN && pt<STAR_KAON_PTMAX){
	      fprintf(detailsfile,"%6.4f\t",pt/1000.0);
	      for (int n=1;n<=nmax;n++)
		{		
		  double r2,i2,normi;
		  r2=(qnreal[ispecies][n][ipt]*Qnreal[n]+qnimag[ispecies][n][ipt]*Qnimag[n]-qnreal[ispecies][0][ipt])/(qnreal[ispecies][0][ipt]*(Qnreal[0]-1));
		  normi=(Qnreal[n]*Qnreal[n]+Qnimag[n]*Qnimag[n]-Qnreal[0])/(Qnreal[0]*(Qnreal[0]-1));
		  fprintf(detailsfile,"%g\t",r2/sqrt(normi));
		}
	      fprintf(detailsfile,"\n");
	    }
	  }           	
	}

	if(ispecies==2){
	  fprintf(detailsfile,"# PROTON_Vn\n");
	  for(ipt=0;ipt<STAR_PROTON_NPTBINS;ipt++){
	    pt=(0.5+ipt)*STAR_PROTON_DELPT;
	    if(pt>STAR_PROTON_PTMIN && pt<STAR_PROTON_PTMAX){
	      fprintf(detailsfile,"%6.4f\t",pt/1000.0);
	      for (int n=1;n<=nmax;n++)
		{		
		  double r2,i2,normi;
		  r2=(qnreal[ispecies][n][ipt]*Qnreal[n]+qnimag[ispecies][n][ipt]*Qnimag[n]-qnreal[ispecies][0][ipt])/(qnreal[ispecies][0][ipt]*(Qnreal[0]-1));
		  normi=(Qnreal[n]*Qnreal[n]+Qnimag[n]*Qnimag[n]-Qnreal[0])/(Qnreal[0]*(Qnreal[0]-1));
		  fprintf(detailsfile,"%g\t",r2/sqrt(normi));
		}
	      fprintf(detailsfile,"\n");
	    }
	  }           	
	}


      }
      fflush(detailsfile);
      fclose(detailsfile);
    }
  } 
}


#endif
