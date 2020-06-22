/*
Causal Viscous Hydro Code for Non-Central Heavy Ion Collisions

by

Ulrike Romatschke and Paul Romatschke

version 0.2 February 2009
...and Matt Luzum
version 0.3 April 2010
version 0.4 September 2010
version 0.5 April 2013
version 1.0 May 2013
version 1.1 June 2013
version 1.3 November 2013
version 1.4 July 2014
version 1.5 July 2014
version 1.7 August 2014

The causal viscous hydro equations have been
derived in 

R. Baier, P. Romatschke, D.T. Son, A. Starinets, \
M. Stephanov, arXiv:0712.2451,
(JHEP 0804:100,2008).

Results from this code were used in 

   P.~Romatschke and U.~Romatschke, arXiv:0706.1522,
   (Phys. Rev. Lett.99, 172301,2007).

The setup and tests are documented in

   M.~Luzum and P.~Romatschke, arXiv:0804.4015
   (Phys.Rev.C78, 034915, 2008)

Some parts of the code, in particular the paramreader module, 
were salvaged and modified from unrelated work by one of the 
authors (PR) with Michael Strickland (MS) and in fact were 
originally written by MS.

If you use this code or some part of it, be sure to reference these
articles. Also, if you use the code package as it is, be sure to
refer to the following articles:

* Equation of State: 
M.~Laine and Y.~Schroder, Phys.\ Rev.\  D {\bf 73} (2006) 085009 [arXiv:hep-ph/0603048].

* Resonance Decay Routines: 
J.~Sollfrank, P.~Koch, U.W.~Heinz,  Z.\ Phys.\  C {\bf 52} (1991) 593.
J.~Sollfrank, P.~Koch and U.~W.~Heinz, hys.\ Lett.\  B {\bf 252} (1990) 256.

* If you use the Color Glass Condensate initial conditions, refer to 

H.~J.~Drescher, A.~Dumitru, A.~Hayashigaki and Y.~Nara, Phys.\ Rev.\  C {\bf 74} (2006) 044905 [arXiv:nucl-th/0605012].

Permission to copy this code is granted provided you keep this disclaimer.

*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>
#include <myspline.h>

//custom defined output
//Betz-Gyulassy
//#define BG 1
//Adare-McCumber-Nagle-Romatschke
#define AMNR 1

using namespace std;

const double fmtoGeV=5.0677;

//smoothing to be able to run lumpy initial conditions
int SMOOTHING=0;
double SMOOTH=0.01;

double SCAL=1.0;

// these global vars are initialized from parameters file
// defaults set here are overridden by that file

int NUMT=8;
long int STEPS=4000,UPDATE=100,SNAPUPDATE=1000;
double AT=0.05,EPS=0.001,B=0.0;
double ETAOS=0.3;
double TSTART=0.5,TF=0.1,TINIT=1.0;
double IC;
int PTASIZE,PHIPASIZE;

//controls value of tau_Pi
double COEFF=3.0;

//controls value of lambda_1
double L1COEF=2.0;
double L2COEF=1.0;

//create freeze-out surface
int FREEZE=1;

char EOSNAME[255];
char ETANAME[255];
char ZETANAME[255];
char BETANAME[255];
char LAMBDANAME[255];

long int B3DEVENTS;
char B3D_PPBAR_ANNIHILATION[255];

double PTMAX, TRIEPS = 0.0, TRIANGLE = 0.0;
double QUADEPS = 0.0, QUADANGLE = 0.0;
double QUINTEPS = 0.0, QUINTANGLE = 0.0;
double SEXEPS = 0.0, SEXANGLE = 0.0;
double SEPTEPS = 0.0, SEPTANGLE = 0.0;
double BIEPS = 0.0, BIANGLE = 0.0;
double MONOEPS = 0.0, MONOANGLE = 0.0;
double PCE;
int NS=0, IFLOW=0;

int preeqflow;
int FULL = 0;

double RNUC, ANUC, SIGMANN, TANORM;

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
double Rnuc=6.4;
//wood-saxon parameter in fm;
double anuc=0.54;


//for the equation of state
double *eoT4,*cs2i,*poT4,*Ti;

// these hold the current values
double ***u,**e,**pixy,**pixx,**piyy,**pib;

// these hold the updated values
double ***U,**E,**Pixy,**Pixx,**Piyy,**Pib;

// these hold the values from the last UPDATE time step
double ***ulast,**elast,**pixylast,**pixxlast,**piyylast,**pilast;

// these hold the past values
double ***upast;

//overall time
double t = 0;

//these are global for convenience; used in doInc
double ****dtmat;

double ****vec;/*bei Pauli rhs*/

//center of lattice 
int Middle;




//for convenience 
//global definition of ut(sx,sy) -- in order not to have it calculated
//a gazillion times
double **globut;

double ***thf;
double ***dtpixx,***dtpixy,***dtpiyy,***dtpi;
double **mypixt,**mypiyt,**mypitt,**mypiee;

// output files
fstream freeze_out;
fstream meta;
fstream ecces;
fstream Tzetaos;

//custom output files

#ifdef BG
fstream bgout;
#endif


//splines -- for fancy freeze-out
gsl_interp_accel *wac;
gsl_spline *workspline;
//splines -- for equation of state
gsl_spline *pspline,*Tspline,*cs2spline;
//splines -- for eta/s
gsl_spline *etaspline,*betaspline,*lambdaspline,*zetaspline;

myspline mpspline,mTspline;

//to know where to stop interpolation
double lowestE,loweta,higheta,lowbeta,highbeta,lowlambda,highlambda,lowzeta,highzeta;

double eos(double mye,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc);
double T(int sx,int sy,gsl_interp_accel *Tacc);
double T(int sx,int sy);
double Tlast(int sx,int sy,gsl_interp_accel *Tacc);


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

 pib = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   pib[i] = new double[NUMT+2];

//////////////////////
 ulast = new double**[2];

 for (int i=0;i<2;i++) 
   ulast[i] = new double*[NUMT+2];

 for (int i=0;i<2;i++)
   for (int j=0;j<NUMT+2;j++) 
     ulast[i][j] = new double[NUMT+2];

	 
 elast = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
    elast[i] = new double[NUMT+2];
 
 pixylast = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   pixylast[i] = new double[NUMT+2];


 pixxlast = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   pixxlast[i] = new double[NUMT+2];

 piyylast = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   piyylast[i] = new double[NUMT+2];

 pilast = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   pilast[i] = new double[NUMT+2];
//////////////////////

 U = new double**[2];

 for (int i=0;i<2;i++) 
   U[i] = new double*[NUMT+2];

 for (int i=0;i<2;i++)
   for (int j=0;j<NUMT+2;j++) 
     U[i][j] = new double[NUMT+2];

	 
 E = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   E[i] = new double[NUMT+2];


 Pixy = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   Pixy[i] = new double[NUMT+2];


 Pixx = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   Pixx[i] = new double[NUMT+2];


 Piyy = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   Piyy[i] = new double[NUMT+2];

 Pib = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++) 
   Pib[i] = new double[NUMT+2];

 upast = new double**[2];

 for (int i=0;i<2;i++) 
   upast[i] = new double*[NUMT+2];

 for (int i=0;i<2;i++)
   for (int j=0;j<NUMT+2;j++) 
     upast[i][j] = new double[NUMT+2];

 dtmat= new double***[NUMT+2];
 vec=new double***[NUMT+2];

 for (int i=0;i<NUMT+2;i++)
   { 
     dtmat[i]=new double**[NUMT+2];
     vec[i]=new double**[NUMT+2];
     for (int j=0;j<NUMT+2;j++) 
       {
	 dtmat[i][j]=new double*[3];
	 vec[i][j]=new double*[3];
	 
	 for (int k=0;k<3;k++) 
	   {
	     dtmat[i][j][k] = new double[3];
	     vec[i][j][k] = new double[3];
	   }
	 //dtmat = new double*[3];
	 //for (int i=0;i<3;i++) dtmat[i] = new double[3];
	 // vec = new double*[3];
	 //for (int i=0;i<3;i++) vec[i] = new double[1];
       }}
  
 globut = new double*[NUMT+2];
 thf = new double**[NUMT+2];
 dtpixx = new double**[NUMT+2];
 dtpixy = new double**[NUMT+2];
 dtpiyy = new double**[NUMT+2];
 dtpi   = new double**[NUMT+2];
 mypixt = new double*[NUMT+2];
 mypiyt = new double*[NUMT+2];
 mypitt = new double*[NUMT+2];
 mypiee = new double*[NUMT+2];

 for (int i=0;i<NUMT+2;i++)
   { 
     globut[i] = new double[NUMT+2];
     mypixt[i] = new double[NUMT+2];
     mypiyt[i] = new double[NUMT+2];
     mypitt[i] = new double[NUMT+2];
     mypiee[i] = new double[NUMT+2];
     thf[i] = new double*[NUMT+2];
     dtpixx[i] = new double*[NUMT+2];
     dtpixy[i] = new double*[NUMT+2];
     dtpi[i]   = new double*[NUMT+2];
     dtpiyy[i] = new double*[NUMT+2];
     for (int j=0;j<NUMT+2;j++) 
       {
	 thf[i][j] = new double[4];
	 dtpixx[i][j] = new double[4];
	 dtpixy[i][j] = new double[4];
	 dtpi[i][j] = new double[4];
	 dtpiyy[i][j] = new double[4];
       }
   }

}

//enforce periodic BoundaryConditions

void enforcePBCs()
{
  
  for(int sy=1;sy<=NUMT;sy++)
    {
      u[0][0][sy]=u[0][NUMT][sy];
      u[0][NUMT+1][sy]=u[0][1][sy];     
    }
  for(int sx=0;sx<=NUMT+1;sx++)
    {
      u[0][sx][0]=u[0][sx][NUMT];
      u[0][sx][NUMT+1]=u[0][sx][1];
    }
  for(int sy=1;sy<=NUMT;sy++)
    {
      u[1][0][sy]=u[1][NUMT][sy];
      u[1][NUMT+1][sy]=u[1][1][sy];     
    }
  for(int sx=0;sx<=NUMT+1;sx++)
    {
      u[1][sx][0]=u[1][sx][NUMT];
      u[1][sx][NUMT+1]=u[1][sx][1];
   }
  for(int sy=1;sy<=NUMT;sy++)
    {
      e[0][sy]=e[NUMT][sy];
      e[NUMT+1][sy]=e[1][sy]; 
    }
  for(int sx=0;sx<=NUMT+1;sx++)
    {
      e[sx][0]=e[sx][NUMT];
      e[sx][NUMT+1]=e[sx][1];
    }

 for(int sy=1;sy<=NUMT;sy++)
    {
      pixx[0][sy]=pixx[NUMT][sy];
      pixx[NUMT+1][sy]=pixx[1][sy];     
    }

  for(int sx=0;sx<=NUMT+1;sx++)
    {
      pixx[sx][0]=pixx[sx][NUMT];
      pixx[sx][NUMT+1]=pixx[sx][1];
    }

 for(int sy=1;sy<=NUMT;sy++)
    {
      pixy[0][sy]=pixy[NUMT][sy];
      pixy[NUMT+1][sy]=pixy[1][sy];     
    }
  for(int sx=0;sx<=NUMT+1;sx++)
    {
      pixy[sx][0]=pixy[sx][NUMT];
      pixy[sx][NUMT+1]=pixy[sx][1];
    }

 for(int sy=1;sy<=NUMT;sy++)
    {
      piyy[0][sy]=piyy[NUMT][sy];
      piyy[NUMT+1][sy]=piyy[1][sy];     
    }
  for(int sx=0;sx<=NUMT+1;sx++)
    {
      piyy[sx][0]=piyy[sx][NUMT];
      piyy[sx][NUMT+1]=piyy[sx][1];
    }
  for(int sy=1;sy<=NUMT;sy++)
    {
      pib[0][sy]=pib[NUMT][sy];
      pib[NUMT+1][sy]=pib[1][sy];     
    }
}


//enforce Neumann BoundaryConditions

void enforceNBCs2()
{
  
  for(int sy=1;sy<=NUMT;sy++)
    {
      u[0][0][sy]=u[0][1][sy];
      u[0][NUMT+1][sy]=u[0][NUMT][sy];     
    }
  for(int sx=0;sx<=NUMT+1;sx++)
    {
      u[0][sx][0]=u[0][sx][1];
      u[0][sx][NUMT+1]=u[0][sx][NUMT];
    }
  for(int sy=1;sy<=NUMT;sy++)
    {
      u[1][0][sy]=u[1][NUMT][sy];
      u[1][NUMT+1][sy]=u[1][1][sy];     
    }
  for(int sx=0;sx<=NUMT+1;sx++)
    {
      u[1][sx][0]=u[1][sx][1];
      u[1][sx][NUMT+1]=u[1][sx][NUMT];
    }
  for(int sy=1;sy<=NUMT;sy++)
    {
      e[0][sy]=e[1][sy];
      e[NUMT+1][sy]=e[NUMT][sy]; 
    }
  for(int sx=0;sx<=NUMT+1;sx++)
    {
      e[sx][0]=e[sx][1];
      e[sx][NUMT+1]=e[sx][NUMT];
    }

  for(int sy=1;sy<=NUMT;sy++)
    {
      pixx[0][sy]=pixx[1][sy];
      pixx[NUMT+1][sy]=pixx[NUMT][sy];     
    }

  for(int sx=0;sx<=NUMT+1;sx++)
    {
      pixx[sx][0]=pixx[sx][1];
      pixx[sx][NUMT+1]=pixx[sx][NUMT];
    }

  for(int sy=1;sy<=NUMT;sy++)
    {
      pixy[0][sy]=pixy[1][sy];
      pixy[NUMT+1][sy]=pixy[NUMT][sy];     
    }
  for(int sx=0;sx<=NUMT+1;sx++)
    {
      pixy[sx][0]=pixy[sx][1];
      pixy[sx][NUMT+1]=pixy[sx][NUMT];
    }

  for(int sy=1;sy<=NUMT;sy++)
    {
      piyy[0][sy]=piyy[1][sy];
      piyy[NUMT+1][sy]=piyy[NUMT][sy];     
    }
  for(int sx=0;sx<=NUMT+1;sx++)
    {
      piyy[sx][0]=piyy[sx][1];
      piyy[sx][NUMT+1]=piyy[sx][NUMT];
    }
}

void enforceNBCs()
{
  
  for(int i=1;i<=NUMT;i++)
    {
      for(int j=0;j<=1;j++)
	{
	  u[j][0][i]=u[0][1][i];
	  u[j][NUMT+1][i]=u[0][NUMT][i];  
	  u[j][i][0]=u[0][i][1];
	  u[j][i][NUMT+1]=u[0][i][NUMT];  
	}
      e[0][i]=e[1][i];
      e[NUMT+1][i]=e[NUMT][i]; 
      e[i][0]=e[i][1];
      e[i][NUMT+1]=e[i][NUMT]; 
      
      pixx[0][i]=pixx[1][i];
      pixx[NUMT+1][i]=pixx[NUMT][i]; 
      pixx[i][0]=pixx[i][1];
      pixx[i][NUMT+1]=pixx[i][NUMT]; 
      
      pixy[0][i]=pixy[1][i];
      pixy[NUMT+1][i]=pixy[NUMT][i]; 
      pixy[i][0]=pixy[i][1];
      pixy[i][NUMT+1]=pixy[i][NUMT]; 
      
      piyy[0][i]=piyy[1][i];
      piyy[NUMT+1][i]=piyy[NUMT][i]; 
      piyy[i][0]=piyy[i][1];
      piyy[i][NUMT+1]=piyy[i][NUMT];
    }
}

//Copy fields at each UPDATE time step
void copyUPDATE()
{
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	for (int i=0;i<2;i++) 
	  ulast[i][sx][sy]=U[i][sx][sy];
	elast[sx][sy]=E[sx][sy];
	pixylast[sx][sy]=Pixy[sx][sy];
	pixxlast[sx][sy]=Pixx[sx][sy];
	piyylast[sx][sy]=Piyy[sx][sy];
	pilast[sx][sy]=Pib[sx][sy];
      }
}

void smeare(double **pe)
{
  // smooth small densities -MPM
   //======================================================
      // j.nagle - 11/01/2013 - copied in McCumber's smoothing
      // for lumpy condition hydro to run 
      // Set to -1.0 to turn off ....
      // p.romatschke - 10.9.2014 modified to run larger systems
      //double SMOOTH = 0.001;  
      // turn smooth on and off via the params file. -MPM

	for (int sx=1;sx<=NUMT;sx++)
		for (int sy=1;sy<=NUMT;sy++)
		{
			if(pe[sx][sy]/pow(AT,4)<SMOOTH)
			{
				/*
				int sx_up = sx+1;
				int sx_dn = sx-1;
				int sy_up = sy+1;
				int sy_dn = sy-1;

				if (sx_up > NUMT) sx_up = NUMT;
				if (sx_dn < 1)   sx_dn = 1;
				if (sy_up > NUMT) sy_up = NUMT;
				if (sy_dn < 1)   sy_dn = 1;

				pe[sx][sy] = (pe[sx_up][sy_up]+pe[sx_up][sy]+pe[sx_up][sy_dn]+pe[sx][sy_up]+pe[sx][sy]+pe[sx][sy_dn]+pe[sx_dn][sy_up]+pe[sx_dn][sy]+pe[sx_dn][sy_dn])/9.0;
				*/
				double tmp_e = 0.0, tmp_norm = 0.0;
				for (int ii=0; ii<5; ii++){
					for (int jj=0; jj<5; jj++){
						int sxx = sx - 2 + ii;
						int syy = sy - 2 + jj;

						if (sxx>NUMT || sxx<1) continue;
						if (syy>NUMT || syy<1) continue;

						tmp_e += pe[sxx][syy];
						tmp_norm += 1.0;
					}//ii
				}//jj

				pe[sx][sy] = tmp_e/tmp_norm;

			}
		}
}

void smearu(double ***pu, double **pe)
{
	for (int sx=1;sx<=NUMT;sx++)
		for (int sy=1;sy<=NUMT;sy++)
		{
			if(pe[sx][sy]/pow(AT,4)<SMOOTH)
			{
				/*
				int sx_up = sx+1;
				int sx_dn = sx-1;
				int sy_up = sy+1;
				int sy_dn = sy-1;
				if (sx_up > NUMT) sx_up = NUMT;
				if (sx_dn < 1)   sx_dn = 1;
				if (sy_up > NUMT) sy_up = NUMT;
				if (sy_dn < 1)   sy_dn = 1;

				pu[0][sx][sy] = (pu[0][sx_up][sy_up]+pu[0][sx_up][sy]+pu[0][sx_up][sy_dn]+pu[0][sx][sy_up]+pu[0][sx][sy]+pu[0][sx][sy_dn]+pu[0][sx_dn][sy_up]+pu[0][sx_dn][sy]+pu[0][sx_dn][sy_dn])/9.0;
				pu[1][sx][sy] = (pu[1][sx_up][sy_up]+pu[0][sx_up][sy]+pu[1][sx_up][sy_dn]+pu[1][sx][sy_up]+pu[1][sx][sy]+pu[1][sx][sy_dn]+pu[1][sx_dn][sy_up]+pu[1][sx_dn][sy]+pu[1][sx_dn][sy_dn])/9.0;
				*/

				double tmp_u0 = 0.0, tmp_u1 = 0.0, tmp_norm = 0.0;
				for (int ii=0; ii<5; ii++){
					for (int jj=0; jj<5; jj++){
						int sxx = sx - 2 + ii;
						int syy = sy - 2 + jj;

						if (sxx>NUMT || sxx<1) continue;
						if (syy>NUMT || syy<1) continue;

						tmp_u0 += pu[0][sxx][syy];
						tmp_u1 += pu[1][sxx][syy];
						tmp_norm += 1.0;
					}//ii
				}//jj

				pu[0][sx][sy] = tmp_u0 / tmp_norm;
				pu[1][sx][sy] = tmp_u1 / tmp_norm;
			}
		}
}

//prepares next time-step
void copyDown() 
{
  if (SMOOTHING)
    {
      smeare(E);
      //smearu(U,E);	
    }//end if smooth======================================================
      
  for (int sx=1;sx<=NUMT;sx++)
    for (int sy=1;sy<=NUMT;sy++)
      {
	for (int i=0;i<2;i++) 
	  {
	    upast[i][sx][sy]=u[i][sx][sy];
	    u[i][sx][sy]=U[i][sx][sy];
	  }
	e[sx][sy]=E[sx][sy];
	pixy[sx][sy]=Pixy[sx][sy];
	pixx[sx][sy]=Pixx[sx][sy];
	piyy[sx][sy]=Piyy[sx][sy];
	pib[sx][sy]=Pib[sx][sy];
      }


    
  enforcePBCs();
  // enforceNBCs();
}


//load eta/s t-dependence
void loadeta()
{
  fstream eosf;

  char eosfile[255];
  
  //extern char ETANAME[255];

  sprintf(eosfile,"input/%s.dat",ETANAME);
  double *dummyx;
  double *dummyy;

  eosf.open(eosfile,ios::in);
  if (eosf.is_open())
    {
      eosf >> length;
      //cout << "Read length of " << length << endl;
      

      dummyx = new double[length];
      dummyy = new double[length];

      for (int i=1;i<=length;i++)
	{
	  eosf >> dummyx[i-1];
	  eosf >> dummyy[i-1];
	  dummyx[i-1]*=AT; //temperature in lattice units
	  if (i==1)
	    loweta=dummyx[i-1];
	  if (i==length)
	    higheta=dummyx[i-1];
	  //printf("x=%f y=%f\n",dummyx[i-1],dummyy[i-1]);
	}

      eosf.close();

      //interpolate
      
      etaspline=gsl_spline_alloc (gsl_interp_cspline, length);
      
      gsl_spline_init (etaspline,dummyx,dummyy,length);  

      //printf("Successfully interpolated between T=%f and %f\n",loweta/AT,higheta/AT);
      printf("===> Info: Successfully loaded eta/s file of length %li\n",length);
      /*
      gsl_interp_accel *macc;
      macc=gsl_interp_accel_alloc ();
      printf("just checking eta/s(%f)=%f=%f\n",dummyx[100],gsl_spline_eval(etaspline,dummyx[100],macc),dummyy[100]);
      gsl_interp_accel_free (macc);
      */
    }
  else
    printf("Could not open eta/s file %s\n",ETANAME);

  


}

//load zeta/s t-dependence
void loadzeta()
{
  fstream eosf;

  char eosfile[255];
  
  //extern char ETANAME[255];

  sprintf(eosfile,"input/%s.dat",ZETANAME);
  double *dummyx;
  double *dummyy;
 

  eosf.open(eosfile,ios::in);
  if (eosf.is_open())
    {
      eosf >> length;
      //cout << "Read length of " << length << endl;
      

      dummyx = new double[length];
      dummyy = new double[length];

      for (int i=1;i<=length;i++)
	{
	  eosf >> dummyx[i-1];
	  eosf >> dummyy[i-1];
	  dummyx[i-1]*=AT; //temperature in lattice units
	  if (i==1)
	    lowzeta=dummyx[i-1];
	  if (i==length)
	    highzeta=dummyx[i-1];
	  //printf("x=%f y=%f\n",dummyx[i-1],dummyy[i-1]);
	}

      eosf.close();

      //interpolate
      
      zetaspline=gsl_spline_alloc (gsl_interp_cspline, length);
      
      gsl_spline_init (zetaspline,dummyx,dummyy,length);  

      //printf("Successfully interpolated between T=%f and %f\n",loweta/AT,higheta/AT);
      printf("===> Info: Successfully loaded zeta/s file of length %li\n",length);
      /*
      gsl_interp_accel *macc;
      macc=gsl_interp_accel_alloc ();
      printf("just checking eta/s(%f)=%f=%f\n",dummyx[100],gsl_spline_eval(etaspline,dummyx[100],macc),dummyy[100]);
      gsl_interp_accel_free (macc);
      */
    }
  else
    printf("Could not open zeta/s file %s\n",ZETANAME);

  


}


//load taupi/eta/2 t-dependence
void loadbeta2()
{
  fstream eosf;

  char eosfile[255];
  
  sprintf(eosfile,"input/%s.dat",BETANAME);
  double *dummyx;
  double *dummyy;

  eosf.open(eosfile,ios::in);
  if (eosf.is_open())
    {
      eosf >> length;
      //cout << "Read length of " << length << endl;
      

      dummyx = new double[length];
      dummyy = new double[length];

      for (int i=1;i<=length;i++)
	{
	  eosf >> dummyx[i-1];
	  eosf >> dummyy[i-1];
	  dummyx[i-1]*=AT; //temperature in lattice units
	  if (i==1)
	    lowbeta=dummyx[i-1];
	  if (i==length)
	    highbeta=dummyx[i-1];
	  //printf("x=%f y=%f\n",dummyx[i-1],dummyy[i-1]);
	}

      eosf.close();

      //interpolate
      
      betaspline=gsl_spline_alloc (gsl_interp_cspline, length);
      
      gsl_spline_init (betaspline,dummyx,dummyy,length);  

      //printf("Successfully interpolated between T=%f and %f\n",lowbeta/AT,highbeta/AT);
      printf("===> Info: Successfully loaded beta2 file of length %li\n",length);
    }
  else
    cout << "Could not open beta2 file" << endl; 
}


//load taupi/eta/2 t-dependence
void loadlambda1()
{
  fstream eosf;

  char eosfile[255];
  
  sprintf(eosfile,"input/%s.dat",LAMBDANAME);
  double *dummyx;
  double *dummyy;

  eosf.open(eosfile,ios::in);
  if (eosf.is_open())
    {
      eosf >> length;
      //cout << "Read length of " << length << endl;
      

      dummyx = new double[length];
      dummyy = new double[length];

      for (int i=1;i<=length;i++)
	{
	  eosf >> dummyx[i-1];
	  eosf >> dummyy[i-1];
	  dummyx[i-1]*=AT; //temperature in lattice units
	  if (i==1)
	    lowlambda=dummyx[i-1];
	  if (i==length)
	    highlambda=dummyx[i-1];
	  //printf("x=%f y=%f\n",dummyx[i-1],dummyy[i-1]);
	}

      eosf.close();

      //interpolate
      
      lambdaspline=gsl_spline_alloc (gsl_interp_cspline, length);
      
      gsl_spline_init (lambdaspline,dummyx,dummyy,length);  
      printf("===> Info: Successfully loaded lambda1 file of length %li\n",length);
      //printf("Successfully interpolated between T=%f and %f\n",lowbeta/AT,highbeta/AT);
    }
  else
    cout << "Could not open lambda1 file" << endl; 

}


//load equation of state
void loadeos()
{
  fstream eosf;

  char eosfile[255];
  

  sprintf(eosfile,"input/%s.dat",EOSNAME);
  eosf.open(eosfile,ios::in);

  if (eosf.is_open())
    {
      eosf >> length;
      //cout << "Read length of " << length << endl;

      

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

      printf("===> Info: Successfully loaded EoS file of length %li\n",length);
    }
  else
    cout << "Could not open EOS file" << endl; 

  //interpolate
  
  

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

  mpspline.alloc(warrx,warry,length);
 


  gsl_spline_init (pspline,warrx,warry,length);
  
  for (int i=0;i<length;i++)
    {
      warry[i]=Ti[i];
    }  

  mTspline.alloc(warrx,warry,length);
  gsl_spline_init (Tspline,warrx,warry,length);
  
  for (int i=0;i<length;i++)
    {
      warry[i]=cs2i[i];
    }  
  gsl_spline_init (cs2spline,warrx,warry,length);

}

//return interpolated temperature in lattice units
double getintT(int i,double x)
{
  return (Ti[i]+x*(Ti[i+1]-Ti[i]))*AT;
}

//Wood-Saxon routine, in physical units
double WS(double x, double y, double z)
{
  
  //cout <<"R= " << R << endl;

  double temp;
  temp=x*x+y*y+z*z;
  temp=sqrt(temp);
  temp-=Rnuc;
  temp/=anuc;
  return 1/(1+exp(temp));
}

//transverse density T_A, in physical units
double TA(double x, double y)
{
  double temp=0;
  for (int i=0;i<1200;i++)
    temp+=WS(x,y,(i+0.5)*Rnuc/400.)*Rnuc/400.;

  //return result normalized to a gold nucleus
  return 2*temp*197./1175.22;
}


//number density of wounded nucleons
//from Kolb et. al, hep-ph/0103234
double getwnuc(double xx,double yy,double b)
{
  
  double mTAp=TA(xx+b/2.,yy);
  double mTAm=TA(xx-b/2.,yy);

  //return mTA*2*(1.-pow(1-mTA/197.*4.,197.));
  double temp=0;
  temp+=mTAp*(1.-exp(-mTAm*4.));
  temp+=mTAm*(1.-exp(-mTAp*4.));
  return temp;
}

// for use in generating Gubser initial conditions
void RungeKutta4(double *y, double x0, double h, int N,
	void (*dydx)(double x, const double *y, double *f)  )
{
    int i;

    double * f0 = NULL;
    double * f1 = NULL;
    double * f2 = NULL;
    double * f3 = NULL;

    double * yTemp = NULL;

    f0 = (double *) malloc( N * sizeof(double) );
    f1 = (double *) malloc( N * sizeof(double) );
    f2 = (double *) malloc( N * sizeof(double) );
    f3 = (double *) malloc( N * sizeof(double) );

    yTemp = (double *) malloc( N * sizeof(double) );

    dydx( x0, y, f0);

    for (i=0;i<N;++i) yTemp[i] = y[i] + h*f0[i]/2.;
    dydx( x0+h/2., yTemp, f1);

    for (i=0;i<N;++i) yTemp[i] = y[i] + h*f1[i]/2.;
    dydx( x0+h/2., yTemp, f2);

    for (i=0;i<N;++i) yTemp[i] = y[i] + h*f2[i];
    dydx( x0+h, yTemp, f3);

    for (i=0;i<N;++i) y[i] = y[i] + h*(f0[i]+2.*f1[i]+2.*f2[i]+f3[i])/6.;

/*
    for (i=0;i<N;++i) yTemp[i] = y[i] + h*f0[i]/2.;
    dydx( x0+h/2., yTemp, f1);

    for (i=0;i<N;++i) y[i] = y[i] + h*f1[i];
*/

    free(f0);
    free(f1);
    free(f2);
    free(f3);
    free(yTemp);
    return;
}

void dydx(double x, const double *y, double *f)
{

double eta_s, taupi, b;


b=5.;

eta_s = 0.2;
taupi = b*eta_s/y[0];

f[0] = (-2.*tanh(x)/3. + y[1]*tanh(x)/3.)*y[0] ;
f[1] = -y[1]/taupi + 4./3.*eta_s*tanh(x)/(taupi*y[0]) - 4./3.*y[1]*y[1]*tanh(x);

}

// generate initial conditions for Gubser flow
void Get_Gubser(double& eden,double& pi_vector,double rho){

int i,imax;
double Temp,Temp0;
double pi,pi0;
double t,t0, dt, psi[2];

//  condicoes iniciais  //

t0=0.;
Temp0 = 1.2;
pi0 = 0.;

t = t0;
Temp = Temp0;
pi = pi0;

//  definicoes necessarias  //

dt = 0.0001;
if(rho<0){
dt = -0.0001;
}

imax = int( fabs(rho)/dt );

//  Main Loop  //

//for(i=0; i<imax; i++){

while( fabs(t)<fabs(rho) ){  
    psi[0] = Temp;
    psi[1] = pi;

    RungeKutta4( psi, t, dt, 2, dydx);

    t = t + dt;

    Temp = psi[0];
    pi = psi[1];

    }

eden = 3.*(16.+7./2.*3.*2.5)*(3.14159265*3.14159265/90.)*pow(Temp,4.);
pi_vector = pi;

//if((eden != eden) | (pi != pi) | (eden <= 0)) cout << "eden = " << eden << ", pi = " << pi << endl;

// printf("Done! -- %lf    %lf    %lf\n", eden, pi_vector, rho);
}

//reads in initial energy density profile
//from inited.dat (generate by either initE.cpp or your
//favorite routine) as well as flow velocities and shear tensor
void setInitialConditions()
{

  double sig;
  //extern double randGauss(double); // generates a gaussian random number with mean 0

  //Middle
  Middle=(NUMT-1)/2+1;
  //Middle=1;

  printf("mi %i\n",Middle);

  //freeze-out temperature
  TF=TF*AT;

  cout << "TF=" << TF/AT << endl;

  cout << "TSTART=" << TSTART << endl;

  //load equation of state

  loadeos();
  loadeta();
  loadzeta();
  loadbeta2();
  loadlambda1();

  if(int(IC) != -9 ) //all runs except Gubser
    {
      fstream inited,initux,inituy,initpixx,initpixy,initpiyy,itime,initpi;
      itime.open("input/time.dat",ios::in);
      inited.open("input/inited.dat", ios::in);
      initux.open("input/initux.dat", ios::in);
      inituy.open("input/inituy.dat", ios::in);  
      initpixx.open("input/initpixx.dat", ios::in);
      initpixy.open("input/initpixy.dat", ios::in);
      initpiyy.open("input/initpiyy.dat", ios::in);
      initpi.open("input/initpi.dat", ios::in);
  
 
      for (int sx=1;sx<=NUMT;sx++)
	for (int sy=1;sy<=NUMT;sy++)
	  {
	    inited >> e[sx][sy];
	    e[sx][sy]*=SCAL;
	    initux >> u[0][sx][sy];
	    inituy >> u[1][sx][sy];
	    initpixx >> pixx[sx][sy];
	    initpixy >> pixy[sx][sy];
	    initpiyy >> piyy[sx][sy];
	    initpi >> pib[sx][sy];
	  }
 
      itime >> TINIT;

      //convert fm/c to lattice units
      t=TINIT*fmtoGeV/AT;

  

      wac=gsl_interp_accel_alloc (); 
      workspline=gsl_spline_alloc (gsl_interp_cspline, Middle);
  
      inited.close();
      initux.close();
      inituy.close();
      initpixx.close();
      initpixy.close();
      initpiyy.close();
      initpi.close();
    }
  enforcePBCs();
  //  enforceNBCs();

  //if this option is chosen, pre-flow is being generated
  //from ed distribution
  if (preeqflow==2)
    {
      printf("===> Info: preeqflow=2 option chosen, ignoring velocity input files\n");
      


      /*   for (int sx=1;sx<=NUMT;sx++)
	for (int sy=1;sy<=NUMT;sy++)
	  {
	    u[0][sx][sy]=0.;
	    u[1][sx][sy]=0.;
	    }*/

      
      
      if (SMOOTHING)
	for (int a=0;a<40;a++)
	smeare(e);
      
      for (int sx=2;sx<=NUMT-1;sx++)
	for (int sy=2;sy<=NUMT-1;sy++)
	  {
	    u[0][sx][sy]=-0.5*(e[sx+1][sy]-e[sx-1][sy])/e[sx][sy]/3.*t;
	    u[1][sx][sy]=-0.5*(e[sx][sy+1]-e[sx][sy-1])/e[sx][sy]/3.*t;
	  }

      /*
      for (int sx=1;sx<=NUMT;sx++)
	{
	  int sy=1;
	  u[0][sx][sy]=-0.5*(e[sx+1][sy]-e[sx-1][sy])/e[sx][sy]/3.*t;
	  u[1][sx][sy]=-(e[sx][2]-e[sx][1])/e[sx][sy]/3.*t;
	}
      for (int sx=1;sx<=NUMT;sx++)
	{
	  int sy=NUMT;
	  u[0][sx][sy]=-0.5*(e[sx+1][sy]-e[sx-1][sy])/e[sx][sy]/3.*t;
	  u[1][sx][sy]=-(e[sx][NUMT]-e[sx][NUMT-1])/e[sx][sy]/3.*t;
	}
      for (int sy=1;sy<=NUMT;sy++)
	{
	  int sx=1;
	  u[0][sx][sy]=-(e[sx+1][sy]-e[sx][sy])/e[sx][sy]/3.*t;
	  u[1][sx][sy]=-0.5*(e[sx][sy+1]-e[sx][sy-1])/e[sx][sy]/3.*t;
	}
      for (int sy=1;sy<=NUMT;sy++)
	{
	  int sx=NUMT;
	  u[0][sx][sy]=-(e[sx][sy]-e[sx-1][sy])/e[sx][sy]/3.*t;
	  u[1][sx][sy]=-0.5*(e[sx][sy+1]-e[sx][sy-1])/e[sx][sy]/3.*t;
	  }*/

      for (int sx=1;sx<=NUMT;sx++)
	{
	  int sy=1;
	  u[0][sx][sy]=u[0][sx][sy+1];
	  u[1][sx][sy]=u[1][sx][sy+1];
	}
      for (int sx=1;sx<=NUMT;sx++)
	{
	  int sy=NUMT;
	  u[0][sx][sy]=u[0][sx][sy-1];
	  u[1][sx][sy]=u[1][sx][sy-1];
	}
      for (int sy=1;sy<=NUMT;sy++)
	{
	  int sx=1;
	  u[0][sx][sy]=u[0][sx+1][sy];
	  u[1][sx][sy]=u[1][sx+1][sy];
	}
      for (int sy=1;sy<=NUMT;sy++)
	{
	  int sx=NUMT;
	  u[0][sx][sy]=u[0][sx-1][sy];
	  u[1][sx][sy]=u[1][sx-1][sy];
	}

      
      double stopval=1.;
      for (int sx=1;sx<=NUMT;sx++)
	for (int sy=1;sy<=NUMT;sy++)
	  {
	    if (fabs(u[0][sx][sy])>stopval)
	      u[0][sx][sy]*=stopval/fabs(u[0][sx][sy]);
	    if (fabs(u[1][sx][sy])>stopval)
	      u[1][sx][sy]*=stopval/fabs(u[1][sx][sy]);	   
	      } 
      
      if (SMOOTHING)
	for (int a=0;a<10;a++)
	//for (int a=0;a<20;a++)
	smearu(u,e);
    }
}

//gets index associated with energy density -- internal use only
long int geti(double mye)
{
  long int i;
  mye/=AT*AT*AT*AT;
  for (i=0;i<length;i++)
    {
      if (eoT4[i]*Ti[i]*Ti[i]*Ti[i]*Ti[i]>mye)
	break;
    }
  return (i-1);
}

//get fraction associated with index i-- internal use only
double getx(long int i,double mye)
{
  mye/=AT*AT*AT*AT;
  
  double temp;

  temp=(eoT4[i+1]-eoT4[i])/(Ti[i+1]-Ti[i])*(Ti[i+1]+Ti[i])/2;
  temp+=2*(eoT4[i+1]+eoT4[i]);
  temp*=(Ti[i+1]+Ti[i])/2;
  temp*=(Ti[i+1]+Ti[i])/2;
  temp*=(Ti[i+1]+Ti[i])/2;

  return (mye-eoT4[i]*Ti[i]*Ti[i]*Ti[i]*Ti[i])/temp/(Ti[i+1]-Ti[i]);
}


//-----------------------------------------------------------------------------
//declaration of functions, shells for functions (for faster running)
//and spatial derivatives

double ut(int sx,int sy)
{
  double temp=1.;
  temp+=u[0][sx][sy]*u[0][sx][sy];
  temp+=u[1][sx][sy]*u[1][sx][sy];
  return sqrt(temp);
}

double umu(int mu, int sx,int sy)
{
  if (mu<2)
    return u[mu][sx][sy];
  if (mu==2)
    return ut(sx,sy);
  else
    return 0;
}



//returns correct pi
double pi (int delta,int beta,int sx,int sy)
{

  //unconventional notation: 0=x, 1=y, 2=tau, 3=eta

  int phi=0;
  int deltap,betap;

  if (beta==2)
    betap=0;
  if (beta==3)
    betap=3;
  if (beta==1)
    betap=2;
  if (beta==0)
    betap=1;
  if (delta==2)
    deltap=0;
  if (delta==3)
    deltap=3;
  if (delta==1)
    deltap=2;
  if (delta==0)
    deltap=1;
  

  if (deltap>betap)
    {
      phi=betap;
      betap=deltap;
      deltap=phi;
    }
  if (deltap<betap)
    {
      if (betap==1 && deltap==0)
	{
	  return u[0][sx][sy]/ut(sx,sy)*pixx[sx][sy]+u[1][sx][sy]/ut(sx,sy)*pixy[sx][sy];
	}
      if (betap==2)
	{
	  if (deltap==0)
	    return u[0][sx][sy]/ut(sx,sy)*pixy[sx][sy]+u[1][sx][sy]/ut(sx,sy)*piyy[sx][sy];
	  if (deltap==1)
	    return pixy[sx][sy];
	}
      else return 0;
    }

  if (deltap==betap)
    {
      if (deltap==0)
	return u[0][sx][sy]/ut(sx,sy)*u[0][sx][sy]/ut(sx,sy)*pixx[sx][sy]+2*u[0][sx][sy]/ut(sx,sy)*u[1][sx][sy]/ut(sx,sy)*pixy[sx][sy]+u[1][sx][sy]/ut(sx,sy)*u[1][sx][sy]/ut(sx,sy)*piyy[sx][sy];
      if (deltap==1)
	return pixx[sx][sy];
      if (deltap==2)
	return piyy[sx][sy];
      if (deltap==3)
	return (-1)/(t*t)*((-1)*u[0][sx][sy]/ut(sx,sy)*u[0][sx][sy]/ut(sx,sy)*pixx[sx][sy]-2*u[0][sx][sy]/ut(sx,sy)*u[1][sx][sy]/ut(sx,sy)*pixy[sx][sy]-u[1][sx][sy]/ut(sx,sy)*u[1][sx][sy]/ut(sx,sy)*piyy[sx][sy]+pixx[sx][sy]+piyy[sx][sy]);
    }
}


//to make things faster
double pishell (int delta,int beta,int sx,int sy)
{
  int phi=0;
  if (delta>beta)
    {
      phi=beta;
      beta=delta;
      delta=phi;
    }
  if (delta==0&&beta==0)
    return pixx[sx][sy];
  if (delta==0&&beta==1)
    return pixy[sx][sy];
  if (delta==0&&beta==2)
    return mypixt[sx][sy];
  if (delta==0&&beta==3)
    return 0;
  if (delta==1&&beta==1)
    return piyy[sx][sy];
  if (delta==1&&beta==2)
    return mypiyt[sx][sy];
  if (delta==1&&beta==3)
    return 0;
  if (delta==2&&beta==2)
    return mypitt[sx][sy];
  if (delta==2&&beta==3)
    return 0;
  if (delta==3&&beta==3)
    return mypiee[sx][sy];
  
}


//this provides dx u[i]
double dxu(int i,int sx,int sy)
{
  double temp=0;
  if(sx!=1)
    {
      if(sx==NUMT)
	{
	  temp=(umu(i,sx,sy)-umu(i,sx-1,sy));
	}
      else
	temp=(umu(i,sx+1,sy)-umu(i,sx-1,sy))/2;
    }
  else
    temp=(umu(i,sx+1,sy)-umu(i,sx,sy));
  return temp;
}

//this provides dy u[i]
double dyu(int i,int sx,int sy)
{
  double temp=0;
  if(sy!=1)
    {
      if(sy==NUMT)
	{
	  temp=(umu(i,sx,sy)-umu(i,sx,sy-1));
	}
      else
	temp=(umu(i,sx,sy+1)-umu(i,sx,sy-1))/2;
    }
  else
    temp=(umu(i,sx,sy+1)-umu(i,sx,sy));
  return temp;
}

//provides \partial_i u^\mu
double diumu(int i,int mu,int sx,int sy)
{
  double temp=0;
  if (mu!=2)
    {
      if (i==0)
	temp=dxu(mu,sx,sy);
      if (i==1)
	temp=dyu(mu,sx,sy);
    }
  if (mu==2)
    {
      if (i==0)
	temp=umu(0,sx,sy)*dxu(0,sx,sy)+umu(1,sx,sy)*dxu(1,sx,sy);
      if (i==1)
	temp=umu(0,sx,sy)*dyu(0,sx,sy)+umu(1,sx,sy)*dyu(1,sx,sy);
      temp/=umu(2,sx,sy);
    }
  return temp;
}

//this provides dx e
double dxe(int sx,int sy)
{
  double temp=0;
  if(sx!=1)
    {
      if(sx==NUMT)
	{
	  temp=(e[sx][sy]-e[sx-1][sy]);
	}
      else
	temp=(e[sx+1][sy]-e[sx-1][sy])/2;
    }
  else
    temp=(e[sx+1][sy]-e[sx][sy]);
  return temp;
}


//this provides dy e
double dye(int sx,int sy)
{
  double temp=0;
  if(sy!=1)
    {
      if(sy==NUMT)
	{
	  temp=(e[sx][sy]-e[sx][sy-1]);
	}
      else
	temp=(e[sx][sy+1]-e[sx][sy-1])/2;
    }
  else
    temp=(e[sx][sy+1]-e[sx][sy]);
  return temp;
}

//this provides dx pi_bulk
double dxpi(int sx,int sy)
{
  double temp=0;
  if(sx!=1)
    {
      if(sx==NUMT)
	{
	  temp=(pib[sx][sy]-pib[sx-1][sy]);
	}
      else
	temp=(pib[sx+1][sy]-pib[sx-1][sy])/2;
    }
  else
    temp=(pib[sx+1][sy]-pib[sx][sy]);
  return temp;
}

//this provides dy pi_bulk
double dypi(int sx,int sy)
{
  double temp=0;
  if(sy!=1)
    {
      if(sy==NUMT)
	{
	  temp=(pib[sx][sy]-pib[sx][sy-1]);
	}
      else
	temp=(pib[sx][sy+1]-pib[sx][sy-1])/2;
    }
  else
    temp=(pib[sx][sy+1]-pib[sx][sy]);
  return temp;
}

//term f
//di pi^mu,alpha
double djpi(int j,int mu,int alpha,int sx,int sy)
{
  double temp=0;
  if (j==0)
    {
      if(sx!=1)
	{
	  if(sx==NUMT)
	    {
	      temp=(pi(mu,alpha,sx,sy)-pi(mu,alpha,sx-1,sy));
	    }
	  else
	    temp=(pi(mu,alpha,sx+1,sy)-pi(mu,alpha,sx-1,sy))/2;
	}
      else
	temp=(pi(mu,alpha,sx+1,sy)-pi(mu,alpha,sx,sy));
    }
  if (j==1)
    {
      if(sy!=1)
	{
	  if(sy==NUMT)
	    {
	      temp=(pi(mu,alpha,sx,sy)-pi(mu,alpha,sx,sy-1));
	    }
	  else
	    temp=(pi(mu,alpha,sx,sy+1)-pi(mu,alpha,sx,sy-1))/2;
	}
      else
	temp=(pi(mu,alpha,sx,sy+1)-pi(mu,alpha,sx,sy));
    }
  return temp;
}

//returns all non zero gamma
double gamma (int alpha, int beta, int delta)
{

  //unconventional notation: 0=x, 1=y, 2=tau, 3=eta
  if (beta==3)
    {
    if (alpha==2 &&delta==3)
      return t;
    else
      if (alpha==3 && delta==2)
	return 1/t;
      else 
	return 0;
    }
  if (beta==2 && alpha==3 && delta==3)
    {
      return 1/t;
    }
  else
    return 0;
}


//metric 2 upper indices
double g(int alpha,int beta)
{
  if (alpha!=beta)
    return 0.;
  else
    {
      if (alpha==0 || alpha==1)
	return -1.;
      if (alpha==2)
	return 1.;
      if (alpha==3)
	return (-1)/(t*t);
    }
}

//metric 2 lower indices
double gdown(int alpha,int beta)
{
  if (alpha!=beta)
    return 0.;
  else
    {
      if ((alpha==0)||(alpha==1))
	return -1.;
      if (alpha==2)
	return 1.;
      if (alpha==3)
	return -(t*t);
    }
}

//Delta^{mu kappa}
double Delta(int mu,int kappa,int sx,int sy)
{
   return g(mu,kappa)- umu(mu,sx,sy)*umu(kappa,sx,sy);
}

//----------------------------------------------------------------------
//thermodynamic functions

double eos(double mye,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc)
{
  double temp=0;

  double phys=mye/AT/AT/AT/AT;

  if (phys>lowestE)
    {
      temp=gsl_spline_eval(pspline,phys,pacc);
      temp*=AT*AT;
      temp*=AT*AT;
    }
  else
    {
      temp=gsl_spline_eval(cs2spline,lowestE,cs2acc);
      temp*=phys;
      temp*=AT*AT;
      temp*=AT*AT;
    }

  return temp;
}


//Speed of Sound squared = dp/depsilon
double cs2(int sx,int sy,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc)
{

  double temp;
  double phys=e[sx][sy]/AT/AT/AT/AT;

  if (phys>lowestE)
    {
      temp=gsl_spline_eval(cs2spline,phys,cs2acc);
    }
  else
    {
      temp=gsl_spline_eval(cs2spline,lowestE,cs2acc);
    }
  return temp;
}

//provides Temperature*lattice spacing
double T(int sx,int sy,gsl_interp_accel *Tacc)
{
  if (e[sx][sy]<0)
    {
      printf("Negative e at sx=%i sy=%i\n",sx,sy);
    } 

  double temp;
  double phys=e[sx][sy]/AT/AT/AT/AT;

  if (phys>lowestE)
    {
      temp=gsl_spline_eval(Tspline,phys,Tacc);
      temp*=AT;
    }
  else
    {
      temp=sqrtl(sqrtl(e[sx][sy]/eoT4[0]));
    }
  return temp;

}

//provides Temperature*lattice spacing
double T(int sx,int sy)
{
  gsl_interp_accel *macc;
  double res=0;
  res=T(sx,sy,macc);
  gsl_interp_accel_free (macc);
  return res;

}

double Tlast(int sx,int sy,gsl_interp_accel *Tacc)
{

  double temp;
  double phys=elast[sx][sy]/AT/AT/AT/AT;

  if (phys>lowestE)
    {
      temp=gsl_spline_eval(Tspline,phys,Tacc);
      temp*=AT;
    }
  else
    {
      temp=sqrtl(sqrtl(elast[sx][sy]/eoT4[0]));
    }
  return temp;

}

double etaos(double TT)
{
  double result;
  //printf("this is etaos, high=%f low=%f cur=%f\n",loweta,higheta,TT);

  if ((TT>loweta)&&(TT<higheta))
    {
      gsl_interp_accel *macc;
      result=gsl_spline_eval(etaspline,TT,macc);
      gsl_interp_accel_free (macc);
    }
  else
    result=ETAOS;

  return result;
}

double zetaoeta(double TT)
{
  double result;
  //printf("this is zeta over eta, high=%f low=%f cur=%f\n",lowzeta,highzeta,TT);

  if ((TT>lowzeta)&&(TT<highzeta))
    {
      gsl_interp_accel *macc;
      result=gsl_spline_eval(zetaspline,TT,macc)/etaos(TT);
      gsl_interp_accel_free (macc);
    }
  else
    result=0.0;

  return result;
}

double coeff(double TT)
{
  double result;

  if ((TT>lowbeta)&&(TT<highbeta))
    {
      gsl_interp_accel *macc;
      result=gsl_spline_eval(betaspline,TT,macc);
      gsl_interp_accel_free (macc);
    }
  else
    result=COEFF;

  return result;
}

double l1coeff(double TT)
{
  double result;

  if ((TT>lowlambda)&&(TT<highlambda))
    {
      gsl_interp_accel *macc;
      result=gsl_spline_eval(lambdaspline,TT,macc);
      gsl_interp_accel_free (macc);
    }
  else
    result=L1COEF;

  return result;
}

//Provides tau_Pi/lattice spacing;
double taupi(int sx,int sy,gsl_interp_accel *Tacc)
{
  double temp=coeff(T(sx,sy,Tacc))*2.0*etaos(T(sx,sy,Tacc))/T(sx,sy,Tacc);
  //printf("E %f\n",e[sx][sy]);
  if (isnan(temp)!=0)
    {
      cout << "Error in taupi\n";
    }
  return temp;
}

//eta/taupi
double etataupi(int sx,int sy,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc)
{
  double temp=0.5/coeff(T(sx,sy))*(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc));
  return temp;
}




//-----------------------------------------------------------------------
//functions that are used to fill in the update matrix Eq.(9,10) of 
//nucl-th/0610108


//convention:
//coefficients denote
//a[0]-> coeff multiplying partial_t u[0]
//a[1]-> coeff multiplying partial_t u[1]
//a[2]-> coeff multiplying partial_t e
//a[3]-> remainder


//this provides D u^mu
void Dumu(double *a,int mu,int sx,int sy)
{
  //unconventional notation: 0=x, 1=y, 2=tau, 3=eta
 if (mu==0)
   {
     a[0]=globut[sx][sy];
     a[1]=0;
     a[2]=0;
     a[3]=u[0][sx][sy]*dxu(0,sx,sy);
     a[3]+=u[1][sx][sy]*dyu(0,sx,sy);
}

 if (mu==1)
   {
     a[0]=0;
     a[1]=globut[sx][sy];
     a[2]=0;
     a[3]=u[0][sx][sy]*dxu(1,sx,sy);
     a[3]+=u[1][sx][sy]*dyu(1,sx,sy);
   }
 
 if (mu==2)
   {
     a[0]=u[0][sx][sy];
     a[1]=u[1][sx][sy];
     a[2]=0;
     a[3]=u[0][sx][sy]*dxu(2,sx,sy);
     a[3]+=u[1][sx][sy]*dyu(2,sx,sy);
   }
 if (mu==3)
   {
     a[0]=0;
     a[1]=0;
     a[2]=0;
     a[3]=0;
   }
}



//Second term Grad p:


//this provides \nabla^\mu p
void Nablap(double *a,int mu,int sx,int sy,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc)
{
  double cc2=cs2(sx,sy,pacc,cs2acc);

  //unconventional notation: 0=x, 1=y, 2=tau
  if (mu==0)
    {
      a[0]=0;
      a[1]=0;
      a[2]=(-1)*cc2*u[0][sx][sy]*globut[sx][sy];
      a[3]=(-1)*cc2*u[0][sx][sy]*u[0][sx][sy]*dxe(sx,sy);
      a[3]-=cc2*u[0][sx][sy]*u[1][sx][sy]*dye(sx,sy);
      a[3]-=cc2*dxe(sx,sy);
    }
  if (mu==1)
    {
      a[0]=0;
      a[1]=0;
      a[2]=(-1)*cc2*u[1][sx][sy]*globut[sx][sy];
      a[3]=(-1)*cc2*u[1][sx][sy]*u[0][sx][sy]*dxe(sx,sy);
      a[3]-=cc2*u[1][sx][sy]*u[1][sx][sy]*dye(sx,sy);
      a[3]-=cc2*dye(sx,sy);
    }
 
}


//Third term:


//Second part of third term with gammas:


//g*Delta*gammas

double gDeltagamma(int mu,int sx,int sy)
{
  double temp=0;
  double h1=0;
  double h2=0;
  for (int kappa=0;kappa<=3;kappa++)
    for (int alpha=0;alpha<=3;alpha++)
      for (int beta=0;beta<=3;beta++)
	for (int delta=0;delta<=3;delta++)
	  {
	    h1=gamma(alpha,beta,delta);
	    h2=gamma(beta,beta,delta);
	    if (h1!=0)
	      temp+=Delta(mu,kappa,sx,sy)*gdown(alpha,kappa)*(gamma(alpha,beta,delta)*pishell(delta,beta,sx,sy));
	    if (h2!=0)
	      temp+=Delta(mu,kappa,sx,sy)*gdown(alpha,kappa)*(gamma(beta,beta,delta)*pishell(alpha,delta,sx,sy));
	  }
  return temp;
}

//First part of third term:

//From eq. (13)

//d_t pi^i,alpha
//term a
//double vorticity(.....)
//{
//}
//term b
//<Grad u> upper indices!

double theta3(int sx,int sy)
{
  //unconventional notation: 0=x, 1=y, 2=tau, 3=eta
  double temp=0;
  for (int kappa=0;kappa<=3;kappa++)
    {
      {
	for (int alpha=0;alpha<=3;alpha++)
	  for (int i=0;i<=1;i++)
	    temp+=gdown(alpha,kappa)*Delta(i,kappa,sx,sy)*diumu(i,alpha,sx,sy);
      }
      temp+=gdown(3,kappa)*Delta(3,kappa,sx,sy)*gamma(3,2,3)*globut[sx][sy];
    }
  /*
    if ((sx==Middle+8)&&(sy==Middle))
    {
    double t2=0;
    for (int kappa=0;kappa<=3;kappa++)
    for (int alpha=0;alpha<=3;alpha++)
    for (int i=0;i<=1;i++)
    t2+=gdown(alpha,kappa)*Delta(i,kappa,sx,sy)*diumu(i,alpha,sx,sy);
    
    cout << "theta " << temp;
    cout << "\t where" << t2;
    cout << endl; 
      
    }
  */
  return temp;
}

double theta0(int sx,int sy)
{
  double theta0=Delta(2,2,sx,sy)*u[0][sx][sy]/globut[sx][sy]-Delta(2,0,sx,sy);
  return theta0;
}

double theta1(int sx ,int sy)
{
  double theta1=-Delta(2,1,sx,sy)+Delta(2,2,sx,sy)*u[1][sx][sy]/globut[sx][sy];;
  return theta1;
}


//<nabla_x u_x>
void gradxux(double *a,int sx,int sy)
{
  
  //unconventional notation: 0=x, 1=y, 2=tau, 3=eta
  a[0]=2*Delta(2,0,sx,sy)-2/3.*Delta(0,0,sx,sy)*thf[sx][sy][0];
  a[1]=(-1)*2/3.*Delta(0,0,sx,sy)*thf[sx][sy][1];
  a[2]=0;
  a[3]=2*(Delta(0,0,sx,sy)*dxu(0,sx,sy)+Delta(0,1,sx,sy)*dyu(0,sx,sy));
  a[3]-=2/3.*Delta(0,0,sx,sy)*thf[sx][sy][3];


}

void gradxuy(double *a,int sx,int sy)
{
  //unconventional notation: 0=x, 1=y, 2=tau, 3=eta
  a[0]=Delta(2,1,sx,sy)-2/3.*Delta(0,1,sx,sy)*thf[sx][sy][0];
  a[1]=Delta(2,0,sx,sy)-2/3.*Delta(0,1,sx,sy)*thf[sx][sy][1];
  a[2]=0;
  a[3]=Delta(0,0,sx,sy)*dxu(1,sx,sy);
  a[3]+=Delta(0,1,sx,sy)*dyu(1,sx,sy);
  a[3]+=Delta(0,1,sx,sy)*dxu(0,sx,sy);
  a[3]+=Delta(1,1,sx,sy)*dyu(0,sx,sy);
  a[3]-=2/3.*Delta(0,1,sx,sy)*thf[sx][sy][3];
}

//<nabla_x u_t>
void gradxut(double *a,int sx,int sy)
{
 
  
  //oder
  
  double grxux[4];
  gradxux(grxux,sx,sy);
  double grxuy[4];
  gradxuy(grxuy,sx,sy);
  
  a[0]=-(umu(0,sx,sy)/globut[sx][sy]*grxux[0]+umu(1,sx,sy)/globut[sx][sy]*grxuy[0]);
  a[1]=-(umu(0,sx,sy)/globut[sx][sy]*grxux[1]+umu(1,sx,sy)/globut[sx][sy]*grxuy[1]);
  a[2]=0;
  a[3]=-(umu(0,sx,sy)/globut[sx][sy]*grxux[3]+umu(1,sx,sy)/globut[sx][sy]*grxuy[3]);
}

void gradyuy(double *a,int sx,int sy)
{
  //unconventional notation: 0=x, 1=y, 2=tau, 3=eta
  a[0]=(-1)*2/3.*Delta(1,1,sx,sy)*thf[sx][sy][0];
  a[1]=2*Delta(2,1,sx,sy)-2/3.*Delta(1,1,sx,sy)*thf[sx][sy][1];
  a[2]=0;
  a[3]=2*(Delta(0,1,sx,sy)*dxu(1,sx,sy)+Delta(1,1,sx,sy)*dyu(1,sx,sy));
  a[3]-=2/3.*Delta(1,1,sx,sy)*thf[sx][sy][3];
}

//<nabla_y u_t>
void gradyut(double *a,int sx,int sy)
{

  
  //oder
  
  double grxuy[4];
  gradxuy(grxuy,sx,sy);
  double gryuy[4];
  gradyuy(gryuy,sx,sy);
    
  a[0]=-(umu(0,sx,sy)/globut[sx][sy]*grxuy[0]+umu(1,sx,sy)/globut[sx][sy]*gryuy[0]);
  a[1]=-(umu(0,sx,sy)/globut[sx][sy]*grxuy[1]+umu(1,sx,sy)/globut[sx][sy]*gryuy[1]);
  a[2]=0;
  a[3]=-(umu(0,sx,sy)/globut[sx][sy]*grxuy[3]+umu(1,sx,sy)/globut[sx][sy]*gryuy[3]);
} 

void gradtut(double *a,int sx,int sy)
{
 
  
  //oder
  
  double grxux[4];
  gradxux(grxux,sx,sy);
  double grxuy[4];
  gradxuy(grxuy,sx,sy);
  double gryuy[4];
  gradyuy(gryuy,sx,sy);
  
  a[0]=umu(0,sx,sy)/globut[sx][sy]*umu(0,sx,sy)/globut[sx][sy]*grxux[0];
  a[0]+=2*umu(0,sx,sy)/globut[sx][sy]*umu(1,sx,sy)/globut[sx][sy]*grxuy[0];
  a[0]+=umu(1,sx,sy)/globut[sx][sy]*umu(1,sx,sy)/globut[sx][sy]*gryuy[0];
  a[1]=umu(0,sx,sy)/globut[sx][sy]*umu(0,sx,sy)/globut[sx][sy]*grxux[1];
  a[1]+=2*umu(0,sx,sy)/globut[sx][sy]*umu(1,sx,sy)/globut[sx][sy]*grxuy[1];
  a[1]+=umu(1,sx,sy)/globut[sx][sy]*umu(1,sx,sy)/globut[sx][sy]*gryuy[1];
  a[2]=0;
  a[3]=umu(0,sx,sy)/globut[sx][sy]*umu(0,sx,sy)/globut[sx][sy]*grxux[3];
  a[3]+=2*umu(0,sx,sy)/globut[sx][sy]*umu(1,sx,sy)/globut[sx][sy]*grxuy[3];
  a[3]+=umu(1,sx,sy)/globut[sx][sy]*umu(1,sx,sy)/globut[sx][sy]*gryuy[3];
} 

void gradeue(double *a,int sx,int sy)
{
  //unconventional notation: 0=x, 1=y, 2=tau, 3=eta
  a[0]=(-1)*2/3.*t*t*t*t*Delta(3,3,sx,sy)*thf[sx][sy][0];
  a[1]=(-1)*2/3.*t*t*t*t*Delta(3,3,sx,sy)*thf[sx][sy][1];
  a[2]=0;
  a[3]=2*t*t*t*t*Delta(3,3,sx,sy)*gamma(3,2,3)*globut[sx][sy]-2/3.*t*t*t*t*Delta(3,3,sx,sy)*thf[sx][sy][3];
}

//extra terms

//gives D_i u_j
double diuj(int i,int j, int sx,int sy)
{
  double temp;
  if ((i==0)&&(j==0))
    temp=dxu(j,sx,sy);
  if ((i==0)&&(j==1))
    temp=dxu(j,sx,sy);
  if ((i==0)&&(j==3))
    temp=0;
  if ((i==1)&&(j==0))
    temp=dyu(j,sx,sy);
  if ((i==1)&&(j==1))
    temp=dyu(j,sx,sy);
  if ((i==1)&&(j==3))
    temp=0;
  if ((i==3)&&(j==0))
    temp=0;
  if ((i==3)&&(j==1))
    temp=0;
  if ((i==3)&&(j==3))
    temp=-ut(sx,sy)*t;

  return temp;
}


//Omega_{i j}
double prevor(double *a,int i,int j,int sx,int sy)
{
  double bi[4],bj[4];
  Dumu(bi,i,sx,sy);
  Dumu(bj,j,sx,sy);
  
  a[0]=0.5*(umu(i,sx,sy)*bj[0]-umu(j,sx,sy)*bi[0]);
  a[1]=0.5*(umu(i,sx,sy)*bj[1]-umu(j,sx,sy)*bi[1]);
  a[2]=0.5*(umu(i,sx,sy)*bj[2]-umu(j,sx,sy)*bi[2]);
  a[3]=0.5*(umu(i,sx,sy)*bj[3]-umu(j,sx,sy)*bi[3]);
  a[3]-=0.5*(diuj(i,j,sx,sy)-diuj(j,i,sx,sy));

}


//Omega_{xy}
void vorticityxy(double *a,int sx,int sy)
{
  double bx[4],by[4];
  Dumu(bx,0,sx,sy);
  Dumu(by,1,sx,sy);

  a[0]=0.5*(u[1][sx][sy]*bx[0]-u[0][sx][sy]*by[0]);
  a[1]=0.5*(u[1][sx][sy]*bx[1]-u[0][sx][sy]*by[1]);
  a[2]=0.5*(u[1][sx][sy]*bx[2]-u[0][sx][sy]*by[2]);
  a[3]=0.5*(u[1][sx][sy]*bx[3]-u[0][sx][sy]*by[3]);
  a[3]-=0.5*(dxu(1,sx,sy)-dyu(0,sx,sy));

  
}


//gives \Omega_{mu,nu}
double vorticity(double *a,int mu, int nu,int sx,int sy)
{
  if ((mu==2)||(nu==2))
    {
      if ((mu==2)&&(nu==2))
	{
	  a[0]=0;
	  a[1]=0;
	  a[2]=0;
	  a[3]=0;
	}
      if ((mu==2)&&(nu!=2))
	{
	  double bx[4],by[4];
	  prevor(bx,nu,0,sx,sy);
	  prevor(by,nu,1,sx,sy);
	  a[0]=(u[0][sx][sy]*bx[0]+u[1][sx][sy]*by[0])/ut(sx,sy);
	  a[1]=(u[0][sx][sy]*bx[1]+u[1][sx][sy]*by[1])/ut(sx,sy);
	  a[2]=(u[0][sx][sy]*bx[2]+u[1][sx][sy]*by[2])/ut(sx,sy);
	  a[3]=(u[0][sx][sy]*bx[3]+u[1][sx][sy]*by[3])/ut(sx,sy);
	}
      if ((mu!=2)&&(nu==2))
	{
	  double bx[4],by[4];
	  prevor(bx,0,mu,sx,sy);
	  prevor(by,1,mu,sx,sy);
	  a[0]=(u[0][sx][sy]*bx[0]+u[1][sx][sy]*by[0])/ut(sx,sy);
	  a[1]=(u[0][sx][sy]*bx[1]+u[1][sx][sy]*by[1])/ut(sx,sy);
	  a[2]=(u[0][sx][sy]*bx[2]+u[1][sx][sy]*by[2])/ut(sx,sy);
	  a[3]=(u[0][sx][sy]*bx[3]+u[1][sx][sy]*by[3])/ut(sx,sy);
	}
    }
  else 
    {
      prevor(a,mu,nu,sx,sy);
    }

//  if ((sx==15)&&(sy==15))
  //   {
  //  printf("\n here for %i %i",mu,nu);
  //  printf("\t made a[3]=%.12g\n",a[3]);
  //}
}




void terma(double *a,int i,int alpha, int sx,int sy)
{

  //5/2 D ln T
  //a[0]=0;
  //a[1]=0;
  //a[2]=5/8./e[sx][sy]*ut(sx,sy);
  //a[3]=5/8./e[sx][sy]*(u[0][sx][sy]*dxe(sx,sy)+u[1][sx][sy]*dye(sx,sy));

  

  //subtract 1/2 theta
  
  //a[0]-=0.5*theta0(sx,sy);
  //a[1]-=0.5*theta1(sx,sy);
  //a[3]-=0.5*theta3(sx,sy);
  
  //  -4/3 theta

  a[0]=-4./3.*theta0(sx,sy);
  a[1]=-4./3.*theta1(sx,sy);
  a[2]=0;
  a[3]=-4./3.*theta3(sx,sy);


  //multiply by pi^ialpha

  a[0]*=pishell(i,alpha,sx,sy);
  a[1]*=pishell(i,alpha,sx,sy);
  a[2]*=pishell(i,alpha,sx,sy);
  a[3]*=pishell(i,alpha,sx,sy);
  

  
  //subtract vorticity
  

  double b[4];
  //double b2[4];
  
  vorticity(b,0,1,sx,sy);
  

  
  if (i==alpha)
    {
      for (int rho=0;rho<4;rho++)
	{

	  if (alpha!=rho)
	    {
	      vorticity(b,alpha,rho,sx,sy);
  

	      a[0]-=L2COEF*pishell(rho,i,sx,sy)*b[0]*g(alpha,alpha);
	      a[1]-=L2COEF*pishell(rho,i,sx,sy)*b[1]*g(alpha,alpha);
	      a[2]-=L2COEF*pishell(rho,i,sx,sy)*b[2]*g(alpha,alpha);
	      a[3]-=L2COEF*pishell(rho,i,sx,sy)*b[3]*g(alpha,alpha);
	    }
	}
    }
  else
    {
      for (int rho=0;rho<4;rho++)
	{
	  if (alpha!=rho)
	    {
	      vorticity(b,alpha,rho,sx,sy);
	      
	      a[0]-=0.5*L2COEF*pishell(rho,i,sx,sy)*b[0]*g(alpha,alpha);
	      a[1]-=0.5*L2COEF*pishell(rho,i,sx,sy)*b[1]*g(alpha,alpha);
	      a[2]-=0.5*L2COEF*pishell(rho,i,sx,sy)*b[2]*g(alpha,alpha);
	      a[3]-=0.5*L2COEF*pishell(rho,i,sx,sy)*b[3]*g(alpha,alpha);
	    }
	}
  
      for (int rho=0;rho<4;rho++)
	{
	  if (i!=rho)
	    {
	      vorticity(b,i,rho,sx,sy);
  
	      a[0]-=0.5*L2COEF*pishell(rho,alpha,sx,sy)*b[0]*g(i,i);
	      a[1]-=0.5*L2COEF*pishell(rho,alpha,sx,sy)*b[1]*g(i,i);
	      a[2]-=0.5*L2COEF*pishell(rho,alpha,sx,sy)*b[2]*g(i,i);
	      a[3]-=0.5*L2COEF*pishell(rho,alpha,sx,sy)*b[3]*g(i,i);
	    }
	}
    }
  
}



//\eta/taupi/u^tau <\nabla^i u^\alpha>
void termb(double *a,int i,int alpha,int sx,int sy,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc)
{
  a[0]=0;
  a[1]=0;
  a[2]=0;
  a[3]=0;
  
  double grxux[4];
  double grxuy[4];
  double gryuy[4];
  double grxut[4];
  double gryut[4];
 
  int alphap,ip;

  alphap=alpha;
  ip=i;

  gradxux(grxux,sx,sy);
  gradxuy(grxuy,sx,sy);
  gradyuy(gryuy,sx,sy);
  gradxut(grxut,sx,sy);
  gradyut(gryut,sx,sy);


  int phi=0;
  if(alphap<ip)
    {
      phi=alphap;
      alphap=ip;
      ip=phi;
    }

  for(int zaehler=0;zaehler<4;zaehler++)
    {    
      
      if(ip==0 && alphap==0)
	{
	  a[zaehler]+=etataupi(sx,sy,pacc,cs2acc)*grxux[zaehler]/globut[sx][sy];
	  //if (zaehler==0&&(sx==Middle+8)&&sy==Middle)
	  //  cout << "here tb" << a[zaehler]/e[sx][sy] << endl;
	}
      if(ip==0 && alphap==1)
	{
	  a[zaehler]+=etataupi(sx,sy,pacc,cs2acc)*grxuy[zaehler]/globut[sx][sy];
	}
      if(ip==1 && alphap==1)
	{
	  a[zaehler]+=etataupi(sx,sy,pacc,cs2acc)*gryuy[zaehler]/globut[sx][sy];
	}
      if(ip==0 && alphap==2)
	{
	  a[zaehler]-=etataupi(sx,sy,pacc,cs2acc)*grxut[zaehler]/globut[sx][sy];
	}
      if(ip==1 && alphap==2)
	{
	  a[zaehler]-=etataupi(sx,sy,pacc,cs2acc)*gryut[zaehler]/globut[sx][sy];
	}
      //printf("gr %f\n", grxux[3]);    
    }
  
  //  if ((sx==Middle+8)&&(sy==Middle)&&i==0&&alpha==0)
  // cout << "terb " << grxux[3]/globut  << endl;

}

//term c
//1/(tau_Pi u^tau) Pi^{i alpha}
double termc(int i,int alpha,int sx,int sy,gsl_interp_accel *Tacc)
{
  double temp=0;
  temp+=pishell(i,alpha,sx,sy)/(taupi(sx,sy,Tacc)*globut[sx][sy]);
  //printf("tau %f pi %f\n",taupi(sx,sy),pi(i,alpha,sx,sy));
return temp;
}


//term d
//1/u^tau (u^i Pi^alpha_kappa+u^alpha Pi^i_kappa) D u^kappa
void  termd(double *a,int i,int alpha,int sx,int sy)
{
  a[0]=0;
  a[1]=0;
  a[2]=0;
  a[3]=0;
  double Du[4];
	
  for (int kappa=0;kappa<=3;kappa++)
    {
      Dumu(Du,kappa,sx,sy);
      for (int beta=0;beta<=3;beta++)
	for (int zaehler=0;zaehler<4;zaehler++)
	  a[zaehler]+=1/globut[sx][sy]*(umu(i,sx,sy)*gdown(kappa,beta)*pishell(alpha,beta,sx,sy)+umu(alpha,sx,sy)*gdown(kappa,beta)*pishell(i,beta,sx,sy))*Du[zaehler];
    }
}

//term e
//vanishes
/*double terme(int i,int alpha,int sx,int sy)
{
  double temp=0;
  for (int kappa=0;kappa<=3;kappa++)
    for (int beta=0;beta<=3;beta++)
      temp+=umu(kappa,sx,sy)/globut*(gamma(i,kappa,beta)*pi(beta,alpha,sx,sy)+gamma(alpha,kappa,beta)*pi(i,beta,sx,sy));
  return temp;
  }*/



//u^j/u^tau \partial_j Pi^{i alpha}
double termf(int i,int alpha,int sx,int sy)
{
  double temp=0;
  for(int j=0;j<=1;j++)
    temp+=umu(j,sx,sy)/globut[sx][sy]*djpi(j,i,alpha,sx,sy);
  return temp;
}

//shear-shear self coupling
double termg(int i,int alpha, int sx,int sy,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc,gsl_interp_accel *Tacc)
{
  double temp=0;

  for(int j=0;j<=3;j++)
    for (int k=0;k<=3;k++)
      temp+=pishell(k,j,sx,sy)*pishell(k,j,sx,sy)*gdown(j,j)*gdown(k,k);
  
  temp*=Delta(i,alpha,sx,sy);
  temp/=-3.;

  for(int j=0;j<=3;j++)
    temp+=pishell(i,j,sx,sy)*pishell(alpha,j,sx,sy)*gdown(j,j);
    
  temp*=l1coeff(T(sx,sy,Tacc));
  temp/=globut[sx][sy]*taupi(sx,sy,Tacc)*(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc));
  //temp/=4*M_PI*ETAOS*taupi(sx,sy)*(e[sx][sy]+eos(e[sx][sy])); //wrong
				   
  return temp;
}

//d_t pi^{i,alpha} summary
void dtpiialpha(double *a,int i,int alpha,int sx,int sy,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc,gsl_interp_accel *Tacc)
{
  double terd[4];
  termd(terd,i,alpha,sx,sy);
  double terb[4];
  termb(terb,i,alpha,sx,sy,pacc,cs2acc);
  
  double tera[4];
  terma(tera,i,alpha,sx,sy);

  //  tera[0]=0.0;
  // tera[1]=0.0;
  //tera[2]=0.0;
  //tera[3]=0.0;
  

  a[0]=terb[0]-terd[0]+tera[0];
  a[1]=terb[1]-terd[1]+tera[1];
  a[2]=terb[2]-terd[2]+tera[2];
  a[3]=terb[3]-termf(i,alpha,sx,sy)-terd[3]-termc(i,alpha,sx,sy,Tacc)+tera[3];
  a[3]-=termg(i,alpha,sx,sy,pacc,cs2acc,Tacc);

   if (isnan(a[3])!=0)
    {
      cout << "Problem in dtpiialpha" << endl;
      cout << "More specific: " << terb[3] << "\t" << terd[3] << "\t";
      cout << termf(i,alpha,sx,sy) << "\t" << termc(i,alpha,sx,sy,Tacc) << endl;
      printf("sx=%i sy=%i a=%p i=%i alpha=%i globut=%f\n",sx,sy,a,i,alpha,globut[sx][sy]);
    }


}

//\partial_\tau \Pi^{i \alpha}
void dtpishell(double *a,int i,int alpha,int sx,int sy,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc,gsl_interp_accel *Tacc)
{
  int phi;
  if (alpha<i)
    {
      phi=alpha;
      alpha=i;
      i=phi;
    }
  if ((i==0)&&(alpha==0))
    {
      a[0]=dtpixx[sx][sy][0];
      a[1]=dtpixx[sx][sy][1];
      a[2]=dtpixx[sx][sy][2];
      a[3]=dtpixx[sx][sy][3];
    }
  if ((i==0)&&(alpha==1))
    {
      a[0]=dtpixy[sx][sy][0];
      a[1]=dtpixy[sx][sy][1];
      a[2]=dtpixy[sx][sy][2];
      a[3]=dtpixy[sx][sy][3];
    }
  if ((i==0)&&(alpha==2))
    dtpiialpha(a,i,alpha,sx,sy,pacc,cs2acc,Tacc);
  if ((i==1)&&(alpha==1))
    {
      a[0]=dtpiyy[sx][sy][0];
      a[1]=dtpiyy[sx][sy][1];
      a[2]=dtpiyy[sx][sy][2];
      a[3]=dtpiyy[sx][sy][3];
    }
  if ((i==1)&&(alpha==2))
    dtpiialpha(a,i,alpha,sx,sy,pacc,cs2acc,Tacc);

}

//d_t pi^alpha,t 
double help(int alpha,int sx,int sy)
{
  double temp=0;
  for(int i=0;i<=1;i++)
    temp+=pishell(i,alpha,sx,sy)*umu(i,sx,sy)/globut[sx][sy];
  return temp;
}

//u^j/u^tau \partial_\tau \Pi^{j \alpha}
void help2(double *a,int alpha,int sx,int sy,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc,gsl_interp_accel *Tacc)
{
  double dtpiialp[4];
  
  
  a[0]=0;
  a[1]=0;
  a[2]=0;
  a[3]=0;
  
  for(int j=0;j<2;j++)
    {
      dtpishell(dtpiialp,j,alpha,sx,sy,pacc,cs2acc,Tacc);     
      for(int zaehler=0;zaehler<4;zaehler++)
	a[zaehler]+=umu(j,sx,sy)/globut[sx][sy]*dtpiialp[zaehler];
	  //printf("dt %f\n",dtpiialp[3]);
    }

      // }
  
}

//\partial_tau \Pi^{\alpha tau} 
void dtpialphat(double *a,int alpha,int sx,int sy,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc,gsl_interp_accel *Tacc)
{
  double hel2[4];
  help2(hel2,alpha,sx,sy,pacc,cs2acc,Tacc);

  //qcout << "a[1]" << a[1] << endl;

  a[0]=hel2[0]+pishell(0,alpha,sx,sy)/globut[sx][sy]-help(alpha,sx,sy)*umu(0,sx,sy)/(globut[sx][sy]*globut[sx][sy]);
  a[1]=hel2[1]+pishell(1,alpha,sx,sy)/globut[sx][sy]-help(alpha,sx,sy)*umu(1,sx,sy)/(globut[sx][sy]*globut[sx][sy]);
  a[2]=0;
  a[3]=hel2[3];
  

}

//First part of third term

//g_{alpha kappa} \Delta^{\mu \kappa} \partial_\beta \Pi^{\alpha \beta}

void gDeltapi(double *a,int mu,int sx,int sy,gsl_interp_accel *pacc,gsl_interp_accel *cs2acc,gsl_interp_accel *Tacc)
{
  double temp0=0;
  double temp1=0;
  double temp3=0;
  double dtpialpt[4];
  for(int alpha=0;alpha<=3;alpha++)
    {     
      dtpialphat(dtpialpt,alpha,sx,sy,pacc,cs2acc,Tacc);
      for(int kappa=0;kappa<=3;kappa++) 
	{
	  temp0+=gdown(alpha,kappa)*Delta(mu,kappa,sx,sy)*dtpialpt[0];
	  temp1+=gdown(alpha,kappa)*Delta(mu,kappa,sx,sy)*dtpialpt[1];
	  temp3+=gdown(alpha,kappa)*Delta(mu,kappa,sx,sy)*(dtpialpt[3]+djpi(0,alpha,0,sx,sy)+djpi(1,alpha,1,sx,sy));
	}
    }

  a[0]=temp0;
  a[1]=temp1;
  a[2]=0;
  a[3]=temp3;
  
}



//Eq. (2)

//fourth term De

void De(double *a,int sx,int sy)
{
  a[0]=0;
  a[1]=0;
  a[2]=globut[sx][sy];
  a[3]=umu(0,sx,sy)*dxe(sx,sy)+umu(1,sx,sy)*dye(sx,sy);
}


//fifth term (e+p)Nabla_mu u^Mu

//Nabla_mu u^mu
void nablamuumu(double *a,int sx,int sy)
{
  a[0]=thf[sx][sy][0];//-Delta(0,2,sx,sy)+Delta(2,2,sx,sy)*umu(0,sx,sy)/globut;
  a[1]=thf[sx][sy][1];//-Delta(1,2,sx,sy)+Delta(2,2,sx,sy)*umu(1,sx,sy)/globut;
  a[2]=0;
  a[3]=thf[sx][sy][3];//-Delta(0,0,sx,sy)*dxu(0,sx,sy)-Delta(0,1,sx,sy)*dyu(0,sx,sy)-Delta(1,0,sx,sy)*dxu(1,sx,sy)-Delta(1,1,sx,sy)*dyu(1,sx,sy)+Delta(2,0,sx,sy)*dxu(2,sx,sy)+Delta(2,1,sx,sy)*dyu(2,sx,sy);
}

//sixth term 1/2 pi^mu,nu <grad_mu u_nu>
void pigradu(double *a,int sx,int sy)
{
  a[0]=0;
  a[1]=0;
  a[2]=0;
  a[3]=0;
  
  double grxux[4];
  double grxuy[4];
  double gryuy[4];
  double grxut[4];
  double gryut[4];
  double grtut[4];
  double greue[4];

  gradxux(grxux,sx,sy);
  gradxuy(grxuy,sx,sy);
  gradyuy(gryuy,sx,sy);
  gradxut(grxut,sx,sy);
  gradyut(gryut,sx,sy);
  gradtut(grtut,sx,sy);
  gradeue(greue,sx,sy);

  int mup,nup;
  
  for(int zaehler=0;zaehler<=3;zaehler++)
    for(int mu=0;mu<=3;mu++)
        for(int nu=0;nu<=3;nu++)
	  { 
	  mup=mu;
	  nup=nu;
	  
	  int phi=0;
	  if(mu>nu)
	    {
	      phi=mup;
	      mup=nup;
	      nup=phi;
	    }
	  
	  if(mup==0 && nup==0)
	    {
	      a[zaehler]+=1/2.*pishell(mup,nup,sx,sy)*grxux[zaehler];
	    }
	  
	  if(mup==0 && nup==1)
	    {
	      a[zaehler]+=1/2.*pishell(mup,nup,sx,sy)*grxuy[zaehler];
	    }
	  if(mup==0 && nup==2)
	    {
	      a[zaehler]+=1/2.*pishell(mup,nup,sx,sy)*grxut[zaehler];
	    }
	  if(mup==1 && nup==1)
	    {
	      a[zaehler]+=1/2.*pishell(mup,nup,sx,sy)*gryuy[zaehler];
	    }
	  if(mup==1 && nup==2)
	    {
	      a[zaehler]+=1/2.*pishell(mup,nup,sx,sy)*gryut[zaehler];
	    }
	  if(mup==2 && nup==2)
	    {
	      a[zaehler]+=1/2.*pishell(mup,nup,sx,sy)*grtut[zaehler];
	    }
	  if(mup==3 && nup==3)
	    {
	      a[zaehler]+=1/2.*pishell(mup,nup,sx,sy)*greue[zaehler];
	    }
	  }
}


//gauss-jordan elimination procedure
//provides 
//void gaussj(double **a, int n,double **b, int m)
#include "GJE.cpp"
#include "diags.cpp"
#include "bulkvisc.cpp"



//main update routine
inline void doInc(double eps) 
{
  //double Dumx[4],Dumy[4],Nabpx[4],Nabpy[4],gDeltpix[4],gDeltpiy[4],Des[4],nablau[4],pigru[4];

  //cout << "am here " << endl;

	//time_t time1,time2;

	int debug=0;

	int sx,sy;
	long int position;
	int nthreads, tid;
	int chunk=1;
	gsl_interp_accel *pacc,*Tacc,*cs2acc;


#pragma omp parallel shared(nthreads,chunk,u,e,pixy,pib,pixx,piyy,U,E,Pixy,Pixx,Piyy,Pib,t,globut,thf,dtpixx,mypixt,mypiyt,mypitt,mypiee,position,Middle,lowestE) private(sx,sy,tid,pacc,Tacc,cs2acc)
	{

		double Dumx[4],Dumy[4],Nabpx[4],Nabpy[4],gDeltpix[4],gDeltpiy[4],Des[4],nablau[4],pigru[4];

		pacc=gsl_interp_accel_alloc (); 
		Tacc=gsl_interp_accel_alloc (); 
		cs2acc=gsl_interp_accel_alloc (); 

		tid = omp_get_thread_num();
		if (tid == 0)
		{
			nthreads = omp_get_num_threads();
			//printf("===> Info: Number of threads = %d\n", nthreads);
		}
		//printf("===> Info: Thread %d starting...\n",tid);

#pragma omp for schedule(dynamic,chunk)
		for(position=0;position<NUMT*NUMT;position++)
		{
			sx=position%NUMT+1;
			sy=position/NUMT+1;

			//printf("position =%i sx=%i sy=%i thread=%i\n",position,sx,sy,tid);
			//for (sx=1;sx<=NUMT;sx++)
			//for (sy=1;sy<=NUMT;sy++)
			//{
			//these are here to fix the internal indices once
			//instead of having to recalc every time -- faster!
			//globali=geti(e[sx][sy]);
			//globalx=getx(globali,e[sx][sy]);

			globut[sx][sy]=ut(sx,sy);

			//printf("sx=%i sy=%i globut=%f\n",sx,sy,globut[sx][sy]);

			thf[sx][sy][0]=theta0(sx,sy);
			thf[sx][sy][1]=theta1(sx,sy);
			thf[sx][sy][3]=theta3(sx,sy);

			mypixt[sx][sy]=pi(0,2,sx,sy);
			mypiyt[sx][sy]=pi(1,2,sx,sy);
			mypitt[sx][sy]=pi(2,2,sx,sy);
			mypiee[sx][sy]=pi(3,3,sx,sy);

			dtpiialpha(dtpixx[sx][sy],0,0,sx,sy,pacc,cs2acc,Tacc);
			dtpiialpha(dtpixy[sx][sy],0,1,sx,sy,pacc,cs2acc,Tacc);
			dtpiialpha(dtpiyy[sx][sy],1,1,sx,sy,pacc,cs2acc,Tacc);



			Dumu(Dumx,0,sx,sy); 
			Dumu(Dumy,1,sx,sy);
			Nablap(Nabpx,0,sx,sy,pacc,cs2acc);
			Nablap(Nabpy,1,sx,sy,pacc,cs2acc);
			gDeltapi(gDeltpix,0,sx,sy,pacc,cs2acc,Tacc);
			gDeltapi(gDeltpiy,1,sx,sy,pacc,cs2acc,Tacc);
			De(Des,sx,sy);	
			nablamuumu(nablau,sx,sy);
			pigradu(pigru,sx,sy);




			//Matrix



			dtmat[sx][sy][0][0]=(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)-pib[sx][sy])*Dumx[0]-Nabpx[0]
				//+gDeltagamma(0,sx,sy)
				+gDeltpix[0]-u[0][sx][sy]*nablau[0]*zetaoeta(T(sx,sy,Tacc))*etataupi(sx,sy,pacc,cs2acc);
			dtmat[sx][sy][0][1]=(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)-pib[sx][sy])*Dumx[1]-Nabpx[1]
				//+gDeltagamma(0,sx,sy)
				+gDeltpix[1]-u[0][sx][sy]*nablau[1]*zetaoeta(T(sx,sy,Tacc))*etataupi(sx,sy,pacc,cs2acc);
			dtmat[sx][sy][0][2]=(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)-pib[sx][sy])*Dumx[2]-Nabpx[2]//+gDeltagamma(0,sx,sy)
				+gDeltpix[2];
			dtmat[sx][sy][1][0]=(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)-pib[sx][sy])*Dumy[0]-Nabpy[0]//+gDeltagamma(1,sx,sy)
				+gDeltpiy[0]-u[1][sx][sy]*nablau[0]*zetaoeta(T(sx,sy,Tacc))*etataupi(sx,sy,pacc,cs2acc);
			dtmat[sx][sy][1][1]=(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)-pib[sx][sy])*Dumy[1]-Nabpy[1]//+gDeltagamma(1,sx,sy)
				+gDeltpiy[1]-u[1][sx][sy]*nablau[1]*zetaoeta(T(sx,sy,Tacc))*etataupi(sx,sy,pacc,cs2acc);
			dtmat[sx][sy][1][2]=(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)-pib[sx][sy])*Dumy[2]-Nabpy[2]//+gDeltagamma(1,sx,sy)
				+gDeltpiy[2];
			dtmat[sx][sy][2][0]=Des[0]+(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)-pib[sx][sy])*nablau[0]-pigru[0];
			dtmat[sx][sy][2][1]=Des[1]+(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)-pib[sx][sy])*nablau[1]-pigru[1];
			dtmat[sx][sy][2][2]=Des[2]+(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)-pib[sx][sy])*nablau[2]-pigru[2];


			//Vector (written like a 3*0 Matrix)
			vec[sx][sy][0][0]=Nabpx[3]-(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)-pib[sx][sy])*Dumx[3]-gDeltagamma(0,sx,sy)-gDeltpix[3]+u[0][sx][sy]*nablau[3]*zetaoeta(T(sx,sy,Tacc))*etataupi(sx,sy,pacc,cs2acc)-u[0][sx][sy]*pib[sx][sy]/taupi(sx,sy,Tacc)+dxpi(sx,sy);

			if (isnan(vec[sx][sy][0][0])!=0)
			{
				cout << "Problem in v00" << endl;
				cout << "More specific 1 : " << Nabpx[3] << "\t";
				cout << "2 : " << (e[sx][sy]+eos(e[sx][sy],pacc,cs2acc))*Dumx[3] << "\t";
				cout << "3 : " << gDeltagamma(0,sx,sy) << "\t";
				cout << "4 : " << gDeltpix[3] << endl;
			}
			vec[sx][sy][1][0]=Nabpy[3]-(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)-pib[sx][sy])*Dumy[3]-gDeltagamma(1,sx,sy)-gDeltpiy[3]+u[1][sx][sy]*nablau[3]*zetaoeta(T(sx,sy,Tacc))*etataupi(sx,sy,pacc,cs2acc)-u[1][sx][sy]*pib[sx][sy]/taupi(sx,sy,Tacc)+dypi(sx,sy);
			vec[sx][sy][2][0]=(Des[3]+(e[sx][sy]+eos(e[sx][sy],pacc,cs2acc)-pib[sx][sy])*nablau[3])*(-1.0)+pigru[3];

			int check;
			check=gaussj(dtmat[sx][sy],3,vec[sx][sy],1);
			//printf("cs %f\n",cs2(sx,sy));  

			if (check==0)
			{
				for (int i=0;i<2;i++)
					U[i][sx][sy]=u[i][sx][sy]+eps*vec[sx][sy][i][0];
				//printf("vecx %f   vecy %f\n",vec[0][0],vec[1][0]);


				E[sx][sy]=e[sx][sy]+eps*vec[sx][sy][2][0];
				//if smoothing, prevent negative energy densities
				if (SMOOTHING)
				{
					if (E[sx][sy]<0)
						E[sx][sy]=mTspline.low*pow(AT,4);
					if (E[sx][sy]>mTspline.high*pow(AT,4))
						E[sx][sy]=mTspline.high*pow(AT,4);
					//E[sx][sy]=mTspline.low*pow(AT,4);
				}
				Pixx[sx][sy]=pixx[sx][sy]+eps*(dtpixx[sx][sy][0]*vec[sx][sy][0][0]+dtpixx[sx][sy][1]*vec[sx][sy][1][0]+dtpixx[sx][sy][2]*vec[sx][sy][2][0]+dtpixx[sx][sy][3]);


				Pixy[sx][sy]=pixy[sx][sy]+eps*(dtpixy[sx][sy][0]*vec[sx][sy][0][0]+dtpixy[sx][sy][1]*vec[sx][sy][1][0]+dtpixy[sx][sy][2]*vec[sx][sy][2][0]+dtpixy[sx][sy][3]);

				Piyy[sx][sy]=piyy[sx][sy]+eps*(dtpiyy[sx][sy][0]*vec[sx][sy][0][0]+dtpiyy[sx][sy][1]*vec[sx][sy][1][0]+dtpiyy[sx][sy][2]*vec[sx][sy][2][0]+dtpiyy[sx][sy][3]);

				dtpi[sx][sy][0]=nablau[0]*zetaoeta(T(sx,sy,Tacc))*etataupi(sx,sy,pacc,cs2acc);
				dtpi[sx][sy][1]=nablau[1]*zetaoeta(T(sx,sy,Tacc))*etataupi(sx,sy,pacc,cs2acc);
				dtpi[sx][sy][3]=nablau[3]*zetaoeta(T(sx,sy,Tacc))*etataupi(sx,sy,pacc,cs2acc);
				dtpi[sx][sy][3]-=u[0][sx][sy]*dxpi(sx,sy)+u[1][sx][sy]*dypi(sx,sy);
				dtpi[sx][sy][3]-=pib[sx][sy]/taupi(sx,sy,Tacc);



				Pib[sx][sy]=pib[sx][sy]+eps*(dtpi[sx][sy][0]*vec[sx][sy][0][0]+dtpi[sx][sy][1]*vec[sx][sy][1][0]+dtpi[sx][sy][3])/globut[sx][sy];

				//Pixx[sx][sy]=0.0;
				//Pixy[sx][sy]=0.0;
				//Piyy[sx][sy]=0.0;
			}
			else
			{
				cout << "Error! at " << sx  << "  " << sy << endl;
				snapshot(t,pacc,cs2acc,Tacc);
				bflag=1;
			}
		}

		gsl_interp_accel_free (pacc);
		gsl_interp_accel_free (Tacc);
		gsl_interp_accel_free (cs2acc);
		}//end of parallel section

  
}


//main driver routine
void Evolve() 
{

   printf("     Time[fm/c]   T_cent [GeV] \t Ecc(spatial)\t Ecc(momentum)\n");
  //setting step sizes to maximum step size
  double eps = EPS;
  long int i=0;
  //for (long int i=1;i<=STEPS;i++) 
  //{
  while((reachedTf==0)&&(wflag==0)&&(t/fmtoGeV*AT<10.0)) 
    {
      i++;
      // evolve fields eps forward in time storing updated fields in captial vars
	  
	  doInc(eps); 
      
	  gsl_interp_accel *pacc,*cs2acc,*Tacc;
	  pacc=gsl_interp_accel_alloc (); 
	  cs2acc=gsl_interp_accel_alloc (); 
	  Tacc=gsl_interp_accel_alloc (); 
	  

	  // measurements and data dump
	  //if ( (i>1 && (i-1)%UPDATE==0)) 
	  //{
	  if ((i-1)%UPDATE==0) 
	    {
	      outputMeasurements(t,pacc,cs2acc,Tacc);
	      
	      if (FREEZE > 1) 
		{
		  //Copy fields to memory to compare with next UPDATE time step
		  copyUPDATE();
		}
	      
	    }
	  if ((i-1)%SNAPUPDATE==0) 
	    {      
	      snapshot(t,pacc,cs2acc,Tacc); 
	      //snapshotBulkvisc(t,pacc,cs2acc,Tacc); //mh
#ifdef BG
	      betzgyulassy(t,pacc,cs2acc,Tacc);	      
#endif
#ifdef AMNR
	      amnr(t,pacc,cs2acc,Tacc);
#endif
	    }
	  
	  //}
	  
	  gsl_interp_accel_free (pacc);
	  gsl_interp_accel_free (Tacc);
	  gsl_interp_accel_free (cs2acc);

	  //copy fields from capital vars to lowercase vars
	  copyDown();
	  
	  // increment time
	  t += eps;
					       
	  if (bflag==1)
	    break;

	  }
}




void printDivider() 
{
  int dwidth = 13;
  for (int i=0;i<8*dwidth;i++) cout << "-"; cout << endl;
  return;
}

void cleanup()
{
  gsl_spline_free (workspline);
  gsl_interp_accel_free (wac);
  gsl_spline_free (pspline);
  gsl_spline_free (etaspline);
  gsl_spline_free (betaspline);
  gsl_spline_free (lambdaspline);
  gsl_spline_free (Tspline);
  gsl_spline_free (cs2spline);
}

void generatehadronparameters()
{
  FILE * pFile;
  

  pFile = fopen ("parameters/default/fixed.param","w");
  
  if (pFile!=NULL)
    {
      fprintf(pFile,"###################\n");
      fprintf(pFile,"bool B3D_PRHYDRO true\n");
      fprintf(pFile,"double B3D_PR_OMEGA_XYSPACING %f\n",AT/fmtoGeV);
      fprintf(pFile,"double B3D_PR_OMEGA_TAUSPACING %f\n",AT*EPS*UPDATE/fmtoGeV);
      fprintf(pFile,"int B3D_PR_NPRCELLSMAX 10000000\n");
      fprintf(pFile,"int B3D_NEVENTSMAX %li\n",B3DEVENTS);
      fprintf(pFile,"bool B3D_ANNIHILATION_CHECK %s\n",B3D_PPBAR_ANNIHILATION);
      fprintf(pFile,"int B3D_NPARTSMAX 20000\n");
      fprintf(pFile,"double HYDRO_FOTEMP %f\n",TF/AT*1000);
      fprintf(pFile,"int B3D_NACTIONSMAX 100000\n");
      fprintf(pFile,"###################\n");
      
      fprintf(pFile,"bool B3D_OSUHYDRO false\n");
      fprintf(pFile,"double B3D_TAUCOLLMAX 50.0\n");
      fprintf(pFile,"double B3D_ETAMAX 1.0\n");
      fprintf(pFile,"double B3D_XYMAX 20.0\n");
      //fprintf(pFile,"double B3D_XYMAX 40.0\n");
      fprintf(pFile,"int B3D_NETA 4\n");
      //fprintf(pFile,"int B3D_NXY 20\n");
      fprintf(pFile,"int B3D_NXY 10\n");
      fprintf(pFile,"int B3D_NSAMPLE 1\n");
      fprintf(pFile,"string B3D_RESONANCES_DECAYS_FILE progdata/resinfo/decays_pdg.dat\n");
      fprintf(pFile,"string B3D_RESONANCES_INFO_FILE progdata/resinfo/resonances_pdg.dat\n");
      fprintf(pFile,"string B3D_INPUT_DATAROOT output\n");
      fprintf(pFile,"string B3D_OUTPUT_DATAROOT output\n");
      fprintf(pFile,"bool B3D_BJORKEN true\n");
      fprintf(pFile,"bool B3D_ERROR_PRINT false\n");
      
      fprintf(pFile,"double B3D_SIGMAMAX 30.0\n");
      fprintf(pFile,"double B3D_SIGMADEFAULT 1.0\n");
      fprintf(pFile,"bool B3D_BJMAKER_GAUSSIAN false\n");
      fprintf(pFile,"bool B3D_BJMAKER_BALANCE false\n");
      fprintf(pFile,"bool B3D_STAR_ACCEPTANCE false\n");
      fprintf(pFile,"bool B3D_VIZWRITE false\n");
      fprintf(pFile,"double B3D_SIGMAINELASTIC 0.0\n");
      fprintf(pFile,"double B3D_INELASTIC_Q0 200\n");
      fprintf(pFile,"bool B3D_INELASTIC false\n");
      fprintf(pFile,"string B3D_INELASTIC_INFO_FILE inelastic.tmpdouble GLAUBER_PP_ENERGY_FRAC 1.0\n");
      fprintf(pFile,"#######################\n");
      
      
      fclose (pFile);
    }
  else
    printf("==> Info: opening b3d params file failed!\n");


  //dparams.close();
}

int main() 
{
 
  
  extern void readParameters(const char*);
  
  printDivider();

  printf("This is VH2-2.1\n");
  
  readParameters("data/params.txt");
  

  printDivider();

  allocateMemory();
 
  setInitialConditions();

  printf("Programmstart\n");  

  //generate paramter file for hadronic afterburner
  generatehadronparameters();
  
  freeze_out.open("data/freezeout_bulk.dat", ios::out);
  Tzetaos.open("data/Tzetaos_bulk.dat", ios::out);


#ifdef BG
  printf("Opening custom file for Betz/Gyulassy\n");
  bgout.open("data/betzgyulassy.dat",ios::out);
  bgout << "#t [fm/c] \t x [fm] \t y [fm] \t T [GeV] \t u^x \t u^y \t e [GeV^4]\n";
#endif



  meta.open("data/meta.dat", ios::out);
  meta << "# tau [fm]\t" << "T [GeV]\t" << "epsilon [GeV4]\t" << "Phi[GeV4]\t" << "\n";

  ecces.open("data/ecc.dat", ios::out);
  ecces << "#1-tau\t2-e_x\t3-e_p\t4-totalptx\t5-totalpty\t"
	<<  "6-mome1c\t7-mome1s\t8-mome2c\t9-mome2s\t10-mome3c\t11-mome3s\t"
	<<  "12-eps1c\t13-eps1s\t14-eps2c\t15-eps2s\t16-eps3c\t17-eps3s\t"
	<<  "18-eps4c\t19-eps4s\t20-eps5c\t21-eps5s\t22-eps6c\t23-eps5s\n";

  if(SMOOTHING)
    printf("Using SMOOTHING to be able to run lumpy conditions\n");
  
  /*
  gsl_interp_accel *Tacc;
  Tacc=gsl_interp_accel_alloc (); 
  printf("sanity test: 0.08=%f higher=%f\n",etaos(0.2,Tacc),etaos(0.16,Tacc));
  gsl_interp_accel_free (Tacc);*/

  Evolve();
 
  freeze_out.close();
  Tzetaos.close();
 
#ifdef BG 
  bgout.close();
#endif
 

  meta.close();
  ecces.close();

  

  
  
  

  cleanup();

  if (bflag==0)
    {  
      //UVH2+1 succeeded.
      //This returns '0' to the system; status in bash: 'echo "$?"'
      cout << "UVH2+1 finished successfully." << endl;
      return 0;
    }
  else
    {
      //UVH2+1 failed.
      //This returns '3' to the system; status in bash: 'echo "$?"'
      cout << "Aborted because encountered nan." << endl;
      return 3;
    }
  //if necessary, show system vh2 is done
  //int stat=system("cp data/params.txt logdir/vh2-is-done.log");

  //return stat;
}
