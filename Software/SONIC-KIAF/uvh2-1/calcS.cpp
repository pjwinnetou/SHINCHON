#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>
#include <cmath>

using namespace std;

const double fmtoGeV=5.0677;
const double pi=M_PI;

// these global vars are initialized from parameters file
// defaults set here are overridden by that file

int NUMT=8;
long int STEPS=4000,UPDATE=100,SNAPUPDATE=1000;
double B=0.0,AT=0.05,EPS=0.001;
double ETAOS=0.3;
double TSTART=0.5,TF=0.1,TINIT=1.0;
double L1COEF;
double L2COEF;
double IC=0;
int PTASIZE,PHIPASIZE;
int FREEZE=0;
char EOSNAME[255];
double PTMAX, TRIEPS = 0.0, TRIANGLE = 0.0;
double QUADEPS = 0.0, QUADANGLE = 0.0;
double QUINTEPS = 0.0, QUINTANGLE = 0.0;
double SEXEPS = 0.0, SEXANGLE = 0.0;
double SEPTEPS = 0.0, SEPTANGLE = 0.0;
double BIEPS = 0.0, BIANGLE = 0.0;
double MONOEPS = 0.0, MONOANGLE = 0.0;

double SCAL=1.0;

double SMOOTH;
int SMOOTHING;
char ETANAME[255];
char ZETANAME[255];
char BETANAME[255];
char LAMBDANAME[255];
long int B3DEVENTS;

char B3D_PPBAR_ANNIHILATION[255];

int preeqflow;

int FULL = 0;
int PCE = 0, IFLOW = 0, NS = 0;

double RNUC, ANUC, SIGMANN, TANORM;

//controls value of tau_Pi
double COEFF=3.0;

//flags
int wflag=0;
int reachedTf=0;
int bflag=0;

int dof=37;

//also for convenience
long int globali;
double globalx;
//defined in loadeos
long int length;

//radius of nucleus in fm
//double Rnuc = 6.6; //Lead
double Rnuc=6.37; //Gold

//wood-saxon parameter in fm;
//double anuc = 0.55; //Lead
double anuc=0.54;  //Gold

//inelastic nucleon-nucleon cross-section in mb
//double sigmaNN = 60; //LHC (5.5 A*TeV)
double sigmaNN = 40; //RHIC (200 A*GeV)

//normalization of thickness function
// double TAnorm = 2*208/1286.8; //Lead
double TAnorm = 2*197./1175.22; //Gold



//for the equation of state
double *eoT4,*cs2i,*poT4,*Ti;

// these hold the current values
double ***u,**e,**pixy,**pixx,**piyy;

double **cyminp;

double atuomas=Rnuc/125;
double g2mua=0.43328;
double ttuomas;

//overall time
double t = 0;

//center of lattice 
int Middle;




//for convenience 
//global definition of ut(sx,sy) -- in order not to have it calculated
//a gazillion times
double globut;

double thf[4];
double dtpixx[4],dtpixy[4],dtpiyy[4];
double mypixt,mypiyt,mypitt,mypiee;

// output files
fstream cym,inited,initux,inituy,initpixx,initpixy,initpiyy,itime;

//splines -- for fancy freeze-out
gsl_interp_accel *wac;
gsl_spline *workspline;
//splines -- for equation of state
gsl_interp_accel *pacc,*Tacc,*cs2acc;
gsl_spline *pspline,*Tspline,*cs2spline;

//to know where to stop interpolation
double lowestE;

// initialize global arrays
void allocateMemory() {

cout << "==> Allocating memory\n";


 u = new double**[2];

 for (int i=0;i<2;i++) 
   u[i] = new double*[NUMT+2];

 for (int i=0;i<2;i++)
   for (int j=0;j<NUMT+2;j++) 
     u[i][j] = new double[NUMT+2];

	 
 e = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
    e[i] = new double[NUMT+2];
 
 pixy = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   pixy[i] = new double[NUMT+2];


 pixx = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   pixx[i] = new double[NUMT+2];

 piyy = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   piyy[i] = new double[NUMT+2];

 //allocating memory
  cyminp = new double*[500];
  for (int i=0;i<500;i++) 
    cyminp[i] = new double[500];


}


void loadeos()
{
  fstream eosf;
  char eosfile[255];

  extern char EOSNAME[255];

  sprintf(eosfile,"input/%s.dat",EOSNAME);

  printf("===> Info: Loading EOS file info from %s\n",eosfile);

  //ideal equation of state
  //eosf.open("idEOS.dat", ios::in);

  //qcd equation of state from Mikko Laine and York Schroeder, 
  //hep-ph/0603048
  //eosf.open("qcdEOS.dat",ios::in);
  //interpolated equation of state -- numerically more stable
  //eosf.open("qcdIEOS.dat",ios::in);
  eosf.open(eosfile,ios::in);

  if (eosf.is_open())
    {
      eosf >> length;
      cout << "Read length of " << length << endl;

      

      eoT4 = new double[length];
      cs2i = new double[length];
      poT4 = new double[length];
      Ti = new double[length];

      for (int i=1;i<=length;i++)
	{
	  eosf >> Ti[i-1];
	  eosf >> eoT4[i-1];
	  eosf >> poT4[i-1];
	  eosf >> cs2i[i-1];
	}

      eosf.close();
    }
  else
    cout << "Could not open EOS file" << endl; 

  //interpolate
  pacc=gsl_interp_accel_alloc (); 
  Tacc=gsl_interp_accel_alloc (); 
  cs2acc=gsl_interp_accel_alloc (); 
  pspline=gsl_spline_alloc (gsl_interp_cspline, length);
  Tspline=gsl_spline_alloc (gsl_interp_cspline, length);
  cs2spline=gsl_spline_alloc (gsl_interp_cspline, length);
  
  
  double warrx[length];
  double warry[length];

  for (int i=0;i<length;i++)
    {
      warrx[i]=eoT4[i]*Ti[i]*Ti[i]*Ti[i]*Ti[i];
      warry[i]=poT4[i]*Ti[i]*Ti[i]*Ti[i]*Ti[i];
    }

  lowestE=warrx[0];

  gsl_spline_init (pspline,warrx,warry,length);
  
  for (int i=0;i<length;i++)
    {
      warry[i]=Ti[i];
    }  
  gsl_spline_init (Tspline,warrx,warry,length);
  
  for (int i=0;i<length;i++)
    {
      warry[i]=cs2i[i];
    }  
  gsl_spline_init (cs2spline,warrx,warry,length);

}



double Stot2()
{
  double stot=0;
  double T,p;
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	if (e[sx][sy]/pow(AT,4) >= lowestE)
	{
	  T = AT*gsl_spline_eval(Tspline,e[sx][sy]/pow(AT,4),Tacc);
	  //sqrtl(sqrtl(e[sx][sy]/eoT4[0]));
	  //p = AT*AT*AT*AT*e[sx][sy]*gsl_spline_eval(cs2spline,lowestE,cs2acc);
	  p = AT*AT*AT*AT*gsl_spline_eval(pspline,e[sx][sy]/pow(AT,4),pacc);
	  stot+=(e[sx][sy]+p)/T;
	}
      }
  return stot*TINIT*fmtoGeV/AT;
}

double ecc(int n)
{
  double temp[2]={0};
  //holds c.o.m. position w.r.t. the lab frame  
  double r_com[2]={0};
  //total energy density summing over the entire lattice
  double e_tot=0;
  for(int sx=1;sx<NUMT+1;sx++)
    for(int sy=1;sy<NUMT+1;sy++)
      {
	r_com[0]+= double(sx*e[sx][sy]);
	r_com[1]+= double(sy*e[sx][sy]);	
	e_tot+= e[sx][sy];
      }
  //mean x position/x c.o.m. position w.r.t. lab frame in lattice units
  r_com[0]/= e_tot;
  //mean y position/x c.o.m. position w.r.t. lab frame in lattice units
  r_com[1]/= e_tot;
  //These are correct
  //cout << r_com[0] << "\t" << r_com[1] << "\t" << e_tot << endl;

  // calculating the r2
  double r2_tot=0.0;
  for(int sx=1;sx<NUMT+1;sx++)
    for(int sy=1;sy<NUMT+1;sy++)
      {
	double r2=0,sc[2]={0};
	//position w.r.t. c.o.m. frame
	sc[0] = sx - r_com[0];
	sc[1] = sy - r_com[1];
	//squared distance to c.o.m.
	r2 = sc[0]*sc[0] + sc[1]*sc[1];
	//average squared distance weighted by the energy density
	r2_tot+= e[sx][sy] * r2;
	//angular expressions to calculate the moments; this is the actual eccentricity calculation
	temp[0]+= e[sx][sy] * r2 * cos( n * atan2( sc[1] , sc[0] ) );
	temp[1]+= e[sx][sy] * r2 * sin( n * atan2( sc[1] , sc[0] ) );
      }
  for(int i=0;i<2;i++)
    {
      temp[i]/= e_tot;
      temp[i]*= temp[i];
      //cout << temp[i] << endl;
    }
  r2_tot/= e_tot;
  //this returns the n-th order eccentricity; note that this does not coincide with the geometrical eccentriciy \sqrt{1-b^2/a^2}
  return sqrt(temp[0]+temp[1])/r2_tot;
}


void outputEcc(int n)
{
  ofstream exporting ("eccentricities.txt" , ofstream::app);

  if (exporting.is_open())
      {
	exporting << ecc(n) << endl;
      	exporting.close();
      }
  else cout << "Unable to open file " << "eccentricities.txt" << endl;

}

int main() 
{
  
  extern void readParameters(const char*);
  
  readParameters("data/params.txt");  

  allocateMemory();
  
  loadeos();

  
  
  fstream inited,itime;
  inited.open("input/inited.dat", ios::in);
  itime.open("input/time.dat",ios::in);
  itime >> TINIT;

  printf("Initial Time=%f\n",TINIT);
  
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	inited >> e[sx][sy];
	//This scales the energy density
	e[sx][sy]*=SCAL;
      }

   inited.close();
   itime.close();

   double ss=Stot2();
   outputEcc(2);
   cout << "Total Entropy " << ss << endl;
   cout << "Total dNch/dY estimator " << ss/6.5 << endl;
   cout << "Eccentricity eps_2 " << ecc(2) << endl;

 
  return 0;
}
