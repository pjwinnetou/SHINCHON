#ifndef __ANALYZE_HBT_CC__
#define __ANALYZE_HBT_CC__
#include "analyze.h"
#include "coral.h"
#include <omp.h>

int get_istarbin(int ikt); // Defined in this file.

void CAnalyze::CalcHBT(int num_threads){
	if (num_threads <= 0) {
		printf("ERROR: Argument num_threads to CalcHBT is less than or equal to 0. Exiting.\n");
	}
  double Rx=6.0,Ry=6.8,Rz=7.6,lambda=0.775;
  bool WRITEDETAILS=parameter::getB(parmap,"WRITEDETAILS",false);
  //string HBTPAIRS="KK";
  string HBTPAIRS="PIPI";
  string detailsfilename=detailsdirname+"/star_hbt_pion.dat";
  FILE* detailsfile;
  if(WRITEDETAILS){
    detailsfile=fopen(detailsfilename.c_str(),"w");
    fprintf(detailsfile,"# STAR HBT RADII\n");
  }
  if(qualifier=="cent0to5"){
    Rx=10.0; Ry=6.8; Rz=7.6;
  }
  else if(qualifier=="cent5to10"){
    Rx=9.7; Ry=6.6; Rz=7.4;
  }
  else if(qualifier=="cent0to10"){
    Rx=10.2; Ry=6.5; Rz=7.5;
  }
  else if(qualifier=="cent10to15"){
    Rx=8.6; Ry=5.9; Rz=6.5;
  }
  else if(qualifier=="cent10to20"){
    Rx=7.4; Ry=5.3; Rz=5.6;
  }
  else if(qualifier=="cent15to20"){
    Rx=8.2; Ry=5.7; Rz=6.6;
  }
  else if(qualifier=="cent20to30"){
    Rx=7.8; Ry=5.3; Rz=6.2;
  }
  else if(qualifier=="cent30to40"){
    Rx=6.9; Ry=4.6; Rz=5.3;
  }
  else if(qualifier=="cent40to50"){
    Rx=5.8; Ry=4.0; Rz=4.8;
  }
  else{
    printf("Don't recognize qualifier, can't initialize radii for hbt radius search\n");
    exit(1);
  }
  //printf("qualifier=%s, R=(%g,%g,%g)\n",qualifier.c_str(),Rx,Ry,Rz);


  //Rx*=0.7; Ry*=0.7; Rz*=0.7;



  int *idlista,*idlistb,nida,nidb;
  if(HBTPAIRS=="PIPI"){
    idlista=new int[3];
    idlistb=new int[3];
    idlista[0]=idlistb[0]=211;
    idlista[1]=idlistb[1]=-211;
    idlista[2]=idlistb[2]=111;
    nida=nidb=3;
  }
  if(HBTPAIRS=="KK"){
    idlista=new int [4];
    idlistb=new int [4];
    idlista[0]=idlistb[0]=311;
    idlista[1]=idlistb[1]=-311;
    idlista[2]=idlistb[2]=321;
    idlista[3]=idlistb[3]=-321;
    nida=nidb=4;
  }
  int nmc,istarbin;
  char sdirname[120],cfdirname[120];
  double kt,spectra,w,m=139.57,sigma,ktbar;
  double wtot[8];
  CMCPRList ***list=NULL;
  // Create Array for storing source info
  //CCHArray *cf=new CCHArray("coralparameters/apars_hdf5_cf.dat");
  //CWaveFunction_pipluspiplus_sqwell *wf=new CWaveFunction_pipluspiplus_sqwell("coralparameters/wfparameters.dat");
  //CWaveFunction_generic *wf=new CWaveFunction_generic("coralparameters/wfparameters.dat",1,493.7,493.7,1.0);
  //CKernel *kernel=new CKernel("coralparameters/wfparameters.dat");
  //kernel->Calc(wf);
  //kernel->WriteData("/Users/scottepratt/data/kdata/pipluspiplus");
  //kernel->ReadData("/Users/scottepratt/data/kdata/pipluspiplus");
	
  
  // Initialize Source Calc Object
  CSourceCalc_HDF5_MultiBin *scalc=new CSourceCalc_HDF5_MultiBin();
  parameter::set(scalc->spars,"DELPT",25.0);
  parameter::set(scalc->spars,"NPTBINS",18);
  parameter::set(scalc->spars,"NPHIBINS",6);
  parameter::set(scalc->spars,"PTMIN",150.0);
  parameter::set(scalc->spars,"PHIMIN_DEG",0.0);
  parameter::set(scalc->spars,"PHIMAX_DEG",360.0);
  parameter::set(scalc->spars,"YMIN",-1.0);
  parameter::set(scalc->spars,"YMAX",1.0);
  parameter::set(scalc->spars,"NPARTSMAX",20000);
  parameter::set(scalc->spars,"HDF5filename",h5_infilename);
  parameter::set(scalc->spars,"NEVENTSMAX",neventsmax);
  parameter::set(scalc->spars,"ETA_GAUSS",1.2);
  scalc->NPHIBINS=parameter::getI(scalc->spars,"NPHIBINS",6);
  scalc->NPTBINS=parameter::getI(scalc->spars,"NPTBINS",18);
  scalc->DELPT=parameter::getD(scalc->spars,"DELPT",25.0);
  scalc->PTMIN=parameter::getD(scalc->spars,"PTMIN",150.0);
  scalc->PTMAX=scalc->PTMIN+scalc->NPTBINS*scalc->DELPT;
  scalc->DELPHI=90.0/double(scalc->NPHIBINS);
  parameter::PrintPars(scalc->spars);
	
  scalc->SetIDs(idlista,nida,idlistb,nidb);
  int nptbins=scalc->NPTBINS;
  int nphibins=scalc->NPHIBINS;
  int nptbins_fp = nptbins; 
  double **Rout=new double *[nptbins];
  double **Rside=new double *[nptbins];
  double **Rlong=new double *[nptbins];
  double *Routbar=new double[nptbins];
  double *Rsidebar=new double[nptbins];
  double *Rlongbar=new double[nptbins];
  for(int ikt=0;ikt<nptbins;ikt++){
    Rout[ikt]=new double[nphibins];
    Rside[ikt]=new double[nphibins];
    Rlong[ikt]=new double[nphibins];
  }
  for(istarbin=0;istarbin<8;istarbin++){
    wtot[istarbin] = 0.0;
    Routbar[istarbin]=0.0;
    Rsidebar[istarbin]=0.0;
    Rlongbar[istarbin]=0.0;
  }
  double xoff,yoff,zoff,gamma; // Where are xoff, yoff, and zoff used??
	
  scalc->CalcS(list,list);
  istarbin=0;
	

  // Starting parallelisation -- RDW,MH
  printf("Beginning parallel region.\nNumber of threads used: %d.\n", num_threads);
	omp_set_num_threads(num_threads); // Set number of threads to equal the number of cores available.
#pragma omp parallel for                                  				\
  default(shared)                                                 \
  private(kt, w, nmc, spectra, gamma, ktbar, sigma, istarbin)  		\
  firstprivate(m, nptbins_fp, nphibins, Rx, Ry, Rz, lambda)													
  for (int ikt=0; ikt<nptbins; ikt++) {
    // Loop over momentum bins.
    kt=scalc->PTMIN+scalc->DELPT*(0.5+ikt); // Calculate momentum of ikt'th bin.
    printf("XXXXXXXXXXXX kt=%g XXXXXXXXXXXXXXX\n",kt);
    gamma=sqrt(kt*kt+m*m)/m;
		istarbin = get_istarbin(ikt); // Figure out which istarbin (0 through 7) we're in.
    for (int iphi=0;iphi<scalc->NPHIBINS;iphi++) {
      // Loop over angle phi bins.
      nmc=list[ikt][iphi]->GetNMC();
      spectra=double(nmc)/(kt*sqrt(kt*kt+m*m));
      w=double(nmc)*spectra; // Calculate the weight for this particular ikt and iphi.
#pragma omp critical
      {
        wtot[istarbin] += w; // Add this weight to the sum of all the weights for a particular istarbin.
      }
      printf("---- kt=%g, ikt=%d, istarbin=%d, iphi=%d, nmc=%d -----\n",kt,ikt,istarbin,iphi,nmc);
      scalc->CalcEffGaussParsPureBose(list[ikt][iphi],lambda,Rx,Ry,Rz);
      Rout[ikt][iphi]=Rx;
      Rside[ikt][iphi]=Ry;
      Rlong[ikt][iphi]=Rz;
#pragma omp critical
			{
      	Routbar[istarbin]+=Rout[ikt][iphi]*w/gamma;
      	Rsidebar[istarbin]+=Rside[ikt][iphi]*w;
      	Rlongbar[istarbin]+=Rlong[ikt][iphi]*w;
			}
      sigma=0.5;
    }
  }
	sigma = 0.5;
	int istarbin_max = get_istarbin(nptbins - 1);
	for (istarbin = 0; istarbin <= istarbin_max; istarbin++) {
		// Divide by the sum of all weights for a particular istarbin.
		Routbar[istarbin]=Routbar[istarbin]/wtot[istarbin]; 
    Rsidebar[istarbin]=Rsidebar[istarbin]/wtot[istarbin];
    Rlongbar[istarbin]=Rlongbar[istarbin]/wtot[istarbin];
    ktbar=scalc->PTMIN+(0.5+istarbin)*100.0;
		printf("----------- istarbin=%d, ktbar=%g ---------------\n",istarbin,ktbar);	
    printf("ROUT_PION_STAR = %g +/- %g\n",Routbar[istarbin],sigma);
    printf("RSIDE_PION_STAR = %g +/- %g\n",Rsidebar[istarbin],sigma);
    printf("RLONG_PION_STAR = %g +/- %g\n",Rlongbar[istarbin],sigma);
  	if(WRITEDETAILS) fprintf(detailsfile,"%6.2f %6.4f %6.4f %6.4f\n",ktbar,Routbar[istarbin],Rsidebar[istarbin],Rlongbar[istarbin]);
	}
  sigma=0.5;
  Rx=(Routbar[0]+Routbar[1]+Routbar[2]+Routbar[3])/4.0;
  Ry=(Rsidebar[0]+Rsidebar[1]+Rsidebar[2]+Rsidebar[3])/4.0;
  Rz=(Rlongbar[0]+Rlongbar[1]+Rlongbar[2]+Rlongbar[3])/4.0;
  fprintf(output,"double STAR_ROUT_PION %g %g\n",Rx,sigma);
  fprintf(output,"double STAR_RSIDE_PION %g %g\n",Ry,sigma);
  fprintf(output,"double STAR_RLONG_PION %g %g\n",Rz,sigma);

  for(int ikt=0;ikt<nptbins;ikt++){
    for(int iphi=0;iphi<scalc->NPHIBINS;iphi++){
      delete list[ikt][iphi];
    }
    delete [] list[ikt];
  }
  delete [] list;
  //Asourcebarbar->ScaleArray(1.0/wtottot);
  if(WRITEDETAILS){
    fflush(detailsfile);
    fclose(detailsfile);
  }
}

void CAnalyze::GetHBTpars(CPartH5 *ph5,double &tau,double &rout,double &rside,double &rlong){
  double tcompare=15.0;
  double x=ph5->x,y=ph5->y,px=ph5->px,py=ph5->py,m=ph5->mass,pt,et,t;
  double eta=ph5->eta,rapidity=ph5->rapidity;
  tau=ph5->tau;
  rlong=tau*sinh(eta-rapidity);
  t=tau*cosh(eta-rapidity);
  pt=sqrt(px*px+py*py);
  et=sqrt(pt*pt+m*m);
  rout=(px*x+py*y)/pt;
  //rout=rout-(pt/et)*(t-tcompare);
  rside=(px*y-py*x)/pt;
}

int get_istarbin(int ikt) {
  // Figures out which istarbin (0 through 8) a given ikt is in.
  // ikt within range: [0,3], [4,7], [8,11], [12, 17], [18, 23], [24, 29], [30, 35], [36, 41]
  // istarbin:         0,     1,     2,      3,        4,        5,        6,        7
  int istarbin;
  if (ikt < 4) {
    istarbin = 0;
  }
  else if (ikt < 8) {
    istarbin = 1;
  }
  else if (ikt < 12) {
    istarbin = 2;
  }
  else if (ikt < 18) {
    istarbin = 3;
  }
  else if (ikt < 24) {
    istarbin = 4;
  }
  else if (ikt < 30) {
    istarbin = 5;
  }
  else if (ikt < 36) {
    istarbin = 6;
  }
  else if (ikt < 42) {
    istarbin = 7;
  }
  else {
    printf("ERROR: istarbin is not 0 through 7.\n");
    istarbin = -1;
  }
  return istarbin;
}

#endif
