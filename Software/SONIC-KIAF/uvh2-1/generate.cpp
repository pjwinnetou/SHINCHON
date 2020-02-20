#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <math.h>

using namespace std;

double eta200,eta300,etainf;
double zeta;

int ABORT=0;

void loadparams()
{
  char dummychar[255];
  int dummyint;
  double dummydouble;

  fstream parameters;
  parameters.open("input/transport_params.dat", ios::in);
  if(parameters.is_open())
    {
      parameters >> dummychar;
      parameters >> dummydouble;
      eta200=dummydouble;
      parameters >> dummychar;
      parameters >> dummydouble;
      eta300=dummydouble;
      parameters >> dummychar;
      parameters >> dummydouble;
      etainf=dummydouble;
      //parameters >> dummychar;
      //parameters >> dummydouble;
      //etalimT=dummydouble;
      parameters >> dummychar;
      parameters >> dummydouble;
      zeta=dummydouble;
    }
  else
    {
      printf("\nERROR: params file params.txt could not be opened\n");
      ABORT=1;
    }
  parameters.close();
}



void generate_eta()
{
  double T;
  double Tlow=0.1;
  double Thigh=1;
  double Tstep=0.002;
  double etaos;
  //double Tc=0.17;
  fstream etadyn;
  etadyn.open("input/etaos-dyn.dat", ios::out);
  if(etadyn.is_open())
    {      
      T=Tlow;
      int len=int((Thigh-Tlow)/Tstep)+1;
      etadyn << len << "\n";
      //fprintf(etadyn,"%i\n",len);    
      while(T<Thigh+Tstep)
	{
	  if (T<0.2)
	    etaos=eta200;
	  if ((T>0.2)&&(T<0.3))
	    etaos=eta300;
	  if (T>0.4)
	    etaos=etainf;
	  etadyn << T << "\t" << etaos << "\n";
	  T+=Tstep;
	}
    }
  else
    {
      printf("\nERROR: could not generate input/etaos-dyn.dat\n");
    }
  etadyn.close();
}

void generate_zeta()
{
  double T;
  double Tlow=0.1;
  double Thigh=1;
  double Tstep=0.002;
  double zetaos;
  double Tc=0.17;
  fstream zetadyn;
  zetadyn.open("input/zetaos-dyn.dat", ios::out);
  if(zetadyn.is_open())
    {      
      T=Tlow;
      int len=int((Thigh-Tlow)/Tstep)+1;
      zetadyn << len << "\n";
      while(T<Thigh+Tstep)
	{
	  zetaos=zeta;
	  zetadyn << T << "\t" << zetaos << "\n";
	  T+=Tstep;
	}
    }
  else
    {
      printf("\nERROR: could not generate input/zetaos-dyn.dat\n");
    }
  zetadyn.close();
}

int main(int argc, char* argv[])
{
  printf("Generating input files for vh2 from input/transport_params.dat\n");
  loadparams();
  if (!ABORT)
    {
      printf("Successfully loaded parameters\n");
      generate_eta();
      generate_zeta();
    }
}
