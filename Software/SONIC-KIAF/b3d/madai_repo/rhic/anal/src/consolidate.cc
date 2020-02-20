#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <cstdlib>
#include <cstring>
#include "analyze.h"

using namespace std;

void CAnalyze::Consolidate(string run_name){
	int iqual;
	FILE *input,*output;
	string s2,filename;
	double value,error;
	string command="mkdir -p model_results/"+run_name;
	system(command.c_str());
	filename="model_results/"+run_name+"/results.dat";
	output=fopen(filename.c_str(),"w");
	char c1[20],c2[80],c3[40],c4[40];
	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.dat");
	printf("nqualifiers=%d\n",qualifiers.nqualifiers);
	
	for(iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		filename="analysis/"+run_name+"/"+qualifiers.qualifier[iqual]+"/results.dat";
		printf("opening up %s\n",filename.c_str());
		input=fopen(filename.c_str(),"r");
		while(!feof(input)){
			fscanf(input,"%s %s %s %s",c1,c2,c3,c4);
			if(!feof(input)){
				s2=string(c2);
				value=atof(c3);
				error=atof(c4);
				//error=0.06*value;
				if(s2=="STAR_RSIDE_PION" || s2=="STAR_ROUT_PION" || s2=="STAR_RLONG_PION") error=sqrt(0.04+0.06*0.06*value*value);
				if(s2=="STAR_V2_PION_PTWEIGHT") error=sqrt(0.0013*0.0013+0.06*value*0.06*value);
				if(s2=="STAR_V2_KAON_PTWEIGHT" || s2=="STAR_V2_PROTON_PTWEIGHT") error=sqrt(0.004*0.004+0.12*value*0.12*value);
				if(s2=="PHENIX_SPECTRA_PION_MEANPT" || s2=="PHENIX_SPECTRA_PROTON_MEANPT" || s2=="PHENIX_SPECTRA_KAON_MEANPT") error=0.03*value;
				
				s2=qualifiers.qualifier[iqual]+"_"+string(c2);
				fprintf(output,"%s %s %g %g\n",c1,s2.c_str(),value,error);
			}
		}
		fclose(input);
	}
	fclose(output);
}



