#include "b3d.h"
#include "analyze.h"
#include "qualifier.h"
using namespace std;

int main(int argc, char *argv[]){
	if (argc != 2) {
		printf("Usage: analyze run_name\n");
		exit(-1);
	}
	int iqual;
	string run_name=argv[1];
	CAnalyze *anal=new CAnalyze(run_name);
	parameter::set(anal->parmap,"WRITEDETAILS",string("true"));
	CQualifiers qualifiers;
	double ppytemp,phenixpionyield;
	qualifiers.Read("qualifiers.dat");
	printf("nqualifiers=%d\n",qualifiers.nqualifiers);
	phenixpionyield=403;
	for(iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		printf("---------- %s --------------\n",(qualifiers.qualifier[iqual]).c_str());
		anal->SetQualifier(qualifiers.qualifier[iqual]);
		//ppytemp=anal->CalcSpectra_PHENIX();
		//if(qualifiers.qualifier[iqual]=="cent0to5"){
		//	phenixpionyield=ppytemp;
		//	printf("For centrality 0-5%%, PHENIX_PION_YIELD=%g\n",phenixpionyield);
		//	}
		//anal->CalcSpectra_PHENIXppbar();
		//anal->CalcSpectra_STAR();
		//anal->CalcSpectra_PHENIX();
		//anal->CalcV2();
		anal->CalcFlucVnGamma();
		//anal->CalcFlucQn();
		//anal->CalcFlucQnEbE();
		anal->CalcHBT();
 		//anal->CalcGamma();//not functional
	}
	delete anal;
	return 0;
}
