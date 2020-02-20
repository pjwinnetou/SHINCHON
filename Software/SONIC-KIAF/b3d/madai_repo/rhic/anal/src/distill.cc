#ifndef __distill_CC__
#define __distill_CC__

#include "analyze.h"

void CDistill::Distill(int nruns){
	int irun,iqual;
	string qualifier[9]={"cent0to5","cent0to10","cent5to10","cent10to15","cent10to20","cent15to20","cent20to30","cent30to40","cent40to50"};
	string infilename,outfilename;
	string type,varname;
	double value,error;
	char dummy[100],dummy1[100],dummy2[100],dummy3[100],*stopstring,idummy[10];
	string irunstring;
	
	double rout=0.0,rside=0.0,rlong=0.0;
	double phenixpionmeanpt=0.0,phenixkaonmeanpt=0.0,phenixprotonmeanpt=0.0;
	double pionv2=0.0,kaonv2=0.0,protonv2=0.0;
	double phenixpionyield=0.0,phenixkaonyield=0.0,phenixprotonyield=0.0;
	double phenixpionpcaval=0.0,phenixkaonpcaval=0.0,phenixprotonpcaval=0.0;
	int nb=4;
	double wb=1.0/double(nb);
	FILE *infile,*outfile;
	for(irun=1;irun<=nruns;irun++){
		sprintf(idummy,"%d",irun);
		irunstring=idummy;
		
		rout=rside=rlong=phenixpionmeanpt=phenixkaonmeanpt=phenixprotonmeanpt=0.0;
		pionv2=kaonv2=protonv2=phenixpionyield=phenixkaonyield=phenixprotonyield=0.0;
		phenixpionpcaval=phenixkaonpcaval=phenixprotonpcaval=0.0;
		
		for(iqual=0;iqual<9;iqual++){
			infilename="analysis/run"+irunstring+"/"+qualifier[iqual]+"/results.dat";
			printf("opening %s\n",infilename.c_str());
			infile=fopen(infilename.c_str(),"r");
			do{
				fscanf(infile,"%s",dummy);
				if(dummy[0]!='#'){
					if(!feof(infile)){
						fscanf(infile,"%s %s %s",dummy1,dummy2,dummy3);
						type=dummy;
						varname=dummy1;
						value=strtof(dummy2,&stopstring);
						error=strtof(dummy3,&stopstring);
						if(irun==1 && iqual==0) printf("%s = %g +/- %g\n",varname.c_str(),value,error);
						/** Integrated HBT Radii */
						if(varname=="STAR_ROUT_PION_kt200" && qualifier[iqual]=="cent0to10") rout+=0.25*wb*value;
						if(varname=="STAR_ROUT_PION_kt200" && qualifier[iqual]=="cent10to20") rout+=0.25*wb*value;
						if(varname=="STAR_ROUT_PION_kt200" && qualifier[iqual]=="cent20to30") rout+=0.25*wb*value;
						if(varname=="STAR_ROUT_PION_kt200" && qualifier[iqual]=="cent30to40") rout+=0.25*wb*value;
						
						if(varname=="STAR_ROUT_PION_kt300" && qualifier[iqual]=="cent0to10") rout+=0.25*wb*value;
						if(varname=="STAR_ROUT_PION_kt300" && qualifier[iqual]=="cent10to20") rout+=0.25*wb*value;
						if(varname=="STAR_ROUT_PION_kt300" && qualifier[iqual]=="cent20to30") rout+=0.25*wb*value;
						if(varname=="STAR_ROUT_PION_kt300" && qualifier[iqual]=="cent30to40") rout+=0.25*wb*value;
						
						if(varname=="STAR_ROUT_PION_kt400" && qualifier[iqual]=="cent0to10") rout+=0.25*wb*value;
						if(varname=="STAR_ROUT_PION_kt400" && qualifier[iqual]=="cent10to20") rout+=0.25*wb*value;
						if(varname=="STAR_ROUT_PION_kt400" && qualifier[iqual]=="cent20to30") rout+=0.25*wb*value;
						if(varname=="STAR_ROUT_PION_kt400" && qualifier[iqual]=="cent30to40") rout+=0.25*wb*value;
						
						if(varname=="STAR_ROUT_PION_kt500" && qualifier[iqual]=="cent0to10") rout+=0.25*wb*value;
						if(varname=="STAR_ROUT_PION_kt500" && qualifier[iqual]=="cent10to20") rout+=0.25*wb*value;
						if(varname=="STAR_ROUT_PION_kt500" && qualifier[iqual]=="cent20to30") rout+=0.25*wb*value;
						if(varname=="STAR_ROUT_PION_kt500" && qualifier[iqual]=="cent30to40") rout+=0.25*wb*value;
						
						if(varname=="STAR_RSIDE_PION_kt200" && qualifier[iqual]=="cent0to10") rside+=0.25*wb*value;
						if(varname=="STAR_RSIDE_PION_kt200" && qualifier[iqual]=="cent10to20") rside+=0.25*wb*value;
						if(varname=="STAR_RSIDE_PION_kt200" && qualifier[iqual]=="cent20to30") rside+=0.25*wb*value;
						if(varname=="STAR_RSIDE_PION_kt200" && qualifier[iqual]=="cent30to40") rside+=0.25*wb*value;

						if(varname=="STAR_RSIDE_PION_kt300" && qualifier[iqual]=="cent0to10") rside+=0.25*wb*value;
						if(varname=="STAR_RSIDE_PION_kt300" && qualifier[iqual]=="cent10to20") rside+=0.25*wb*value;
						if(varname=="STAR_RSIDE_PION_kt300" && qualifier[iqual]=="cent20to30") rside+=0.25*wb*value;
						if(varname=="STAR_RSIDE_PION_kt300" && qualifier[iqual]=="cent30to40") rside+=0.25*wb*value;

						if(varname=="STAR_RSIDE_PION_kt400" && qualifier[iqual]=="cent0to10") rside+=0.25*wb*value;
						if(varname=="STAR_RSIDE_PION_kt400" && qualifier[iqual]=="cent10to20") rside+=0.25*wb*value;
						if(varname=="STAR_RSIDE_PION_kt400" && qualifier[iqual]=="cent20to30") rside+=0.25*wb*value;
						if(varname=="STAR_RSIDE_PION_kt400" && qualifier[iqual]=="cent30to40") rside+=0.25*wb*value;

						if(varname=="STAR_RSIDE_PION_kt500" && qualifier[iqual]=="cent0to10") rside+=0.25*wb*value;
						if(varname=="STAR_RSIDE_PION_kt500" && qualifier[iqual]=="cent10to20") rside+=0.25*wb*value;
						if(varname=="STAR_RSIDE_PION_kt500" && qualifier[iqual]=="cent20to30") rside+=0.25*wb*value;
						if(varname=="STAR_RSIDE_PION_kt500" && qualifier[iqual]=="cent30to40") rside+=0.25*wb*value;

						if(varname=="STAR_RLONG_PION_kt200" && qualifier[iqual]=="cent0to10") rlong+=0.25*wb*value;
						if(varname=="STAR_RLONG_PION_kt200" && qualifier[iqual]=="cent10to20") rlong+=0.25*wb*value;
						if(varname=="STAR_RLONG_PION_kt200" && qualifier[iqual]=="cent20to30") rlong+=0.25*wb*value;
						if(varname=="STAR_RLONG_PION_kt200" && qualifier[iqual]=="cent30to40") rlong+=0.25*wb*value;
						
						if(varname=="STAR_RLONG_PION_kt300" && qualifier[iqual]=="cent0to10") rlong+=0.25*wb*value;
						if(varname=="STAR_RLONG_PION_kt300" && qualifier[iqual]=="cent10to20") rlong+=0.25*wb*value;
						if(varname=="STAR_RLONG_PION_kt300" && qualifier[iqual]=="cent20to30") rlong+=0.25*wb*value;
						if(varname=="STAR_RLONG_PION_kt300" && qualifier[iqual]=="cent30to40") rlong+=0.25*wb*value;

						if(varname=="STAR_RLONG_PION_kt400" && qualifier[iqual]=="cent0to10") rlong+=0.25*wb*value;
						if(varname=="STAR_RLONG_PION_kt400" && qualifier[iqual]=="cent10to20") rlong+=0.25*wb*value;
						if(varname=="STAR_RLONG_PION_kt400" && qualifier[iqual]=="cent20to30") rlong+=0.25*wb*value;
						if(varname=="STAR_RLONG_PION_kt400" && qualifier[iqual]=="cent30to40") rlong+=0.25*wb*value;

						if(varname=="STAR_RLONG_PION_kt500" && qualifier[iqual]=="cent0to10") rlong+=0.25*wb*value;
						if(varname=="STAR_RLONG_PION_kt500" && qualifier[iqual]=="cent10to20") rlong+=0.25*wb*value;
						if(varname=="STAR_RLONG_PION_kt500" && qualifier[iqual]=="cent20to30") rlong+=0.25*wb*value;
						if(varname=="STAR_RLONG_PION_kt500" && qualifier[iqual]=="cent30to40") rlong+=0.25*wb*value;

						/** Spectra */
						if(varname=="PHENIX_SPECTRA_PION_MEANPT" && qualifier[iqual]=="cent0to10") phenixpionmeanpt+=wb*value;
						if(varname=="PHENIX_SPECTRA_PION_MEANPT" && qualifier[iqual]=="cent10to20") phenixpionmeanpt+=wb*value;
						if(varname=="PHENIX_SPECTRA_PION_MEANPT" && qualifier[iqual]=="cent20to30") phenixpionmeanpt+=wb*value;
						if(varname=="PHENIX_SPECTRA_PION_MEANPT" && qualifier[iqual]=="cent30to40") phenixpionmeanpt+=wb*value;
						
						if(varname=="PHENIX_SPECTRA_KAON_MEANPT" && qualifier[iqual]=="cent0to10") phenixkaonmeanpt+=wb*value;
						if(varname=="PHENIX_SPECTRA_KAON_MEANPT" && qualifier[iqual]=="cent10to20") phenixkaonmeanpt+=wb*value;
						if(varname=="PHENIX_SPECTRA_KAON_MEANPT" && qualifier[iqual]=="cent20to30") phenixkaonmeanpt+=wb*value;
						if(varname=="PHENIX_SPECTRA_KAON_MEANPT" && qualifier[iqual]=="cent30to40") phenixkaonmeanpt+=wb*value;

						if(varname=="PHENIX_SPECTRA_PROTON_MEANPT" && qualifier[iqual]=="cent0to10") phenixprotonmeanpt+=wb*value;
						if(varname=="PHENIX_SPECTRA_PROTON_MEANPT" && qualifier[iqual]=="cent10to20") phenixprotonmeanpt+=wb*value;
						if(varname=="PHENIX_SPECTRA_PROTON_MEANPT" && qualifier[iqual]=="cent20to30") phenixprotonmeanpt+=wb*value;
						if(varname=="PHENIX_SPECTRA_PROTON_MEANPT" && qualifier[iqual]=="cent30to40") phenixprotonmeanpt+=wb*value;
						
						if(varname=="PHENIX_SPECTRA_PION_PCAVAL" && qualifier[iqual]=="cent0to10") phenixpionpcaval+=wb*value;
						if(varname=="PHENIX_SPECTRA_PION_PCAVAL" && qualifier[iqual]=="cent10to20") phenixpionpcaval+=wb*value;
						if(varname=="PHENIX_SPECTRA_PION_PCAVAL" && qualifier[iqual]=="cent20to30") phenixpionpcaval+=wb*value;
						if(varname=="PHENIX_SPECTRA_PION_PCAVAL" && qualifier[iqual]=="cent30to40") phenixpionpcaval+=wb*value;
						
						if(varname=="PHENIX_SPECTRA_KAON_PCAVAL" && qualifier[iqual]=="cent0to10") phenixkaonpcaval+=wb*value;
						if(varname=="PHENIX_SPECTRA_KAON_PCAVAL" && qualifier[iqual]=="cent10to20") phenixkaonpcaval+=wb*value;
						if(varname=="PHENIX_SPECTRA_KAON_PCAVAL" && qualifier[iqual]=="cent20to30") phenixkaonpcaval+=wb*value;
						if(varname=="PHENIX_SPECTRA_KAON_PCAVAL" && qualifier[iqual]=="cent30to40") phenixkaonpcaval+=wb*value;
						
						if(varname=="PHENIX_SPECTRA_PROTON_PCAVAL" && qualifier[iqual]=="cent0to10") phenixprotonpcaval+=wb*value;
						if(varname=="PHENIX_SPECTRA_PROTON_PCAVAL" && qualifier[iqual]=="cent10to20") phenixprotonpcaval+=wb*value;
						if(varname=="PHENIX_SPECTRA_PROTON_PCAVAL" && qualifier[iqual]=="cent20to30") phenixprotonpcaval+=wb*value;
						if(varname=="PHENIX_SPECTRA_PROTON_PCAVAL" && qualifier[iqual]=="cent30to40") phenixprotonpcaval+=wb*value;

						/** YIELDS */
						if(varname=="PHENIX_SPECTRA_PION_YIELD" && qualifier[iqual]=="cent0to10") phenixpionyield+=wb*value;
						if(varname=="PHENIX_SPECTRA_PION_YIELD" && qualifier[iqual]=="cent10to20") phenixpionyield+=wb*value;
						if(varname=="PHENIX_SPECTRA_PION_YIELD" && qualifier[iqual]=="cent20to30") phenixpionyield+=wb*value;
						if(varname=="PHENIX_SPECTRA_PION_YIELD" && qualifier[iqual]=="cent30to40") phenixpionyield+=wb*value;
						
						if(varname=="PHENIX_SPECTRA_KAON_YIELD" && qualifier[iqual]=="cent0to10") phenixkaonyield+=wb*value;
						if(varname=="PHENIX_SPECTRA_KAON_YIELD" && qualifier[iqual]=="cent10to20") phenixkaonyield+=wb*value;
						if(varname=="PHENIX_SPECTRA_KAON_YIELD" && qualifier[iqual]=="cent20to30") phenixkaonyield+=wb*value;
						if(varname=="PHENIX_SPECTRA_KAON_YIELD" && qualifier[iqual]=="cent30to40") phenixkaonyield+=wb*value;
						
						if(varname=="PHENIX_SPECTRA_PROTON_YIELD" && qualifier[iqual]=="cent0to10") phenixprotonyield+=wb*value;
						if(varname=="PHENIX_SPECTRA_PROTON_YIELD" && qualifier[iqual]=="cent10to20") phenixprotonyield+=wb*value;
						if(varname=="PHENIX_SPECTRA_PROTON_YIELD" && qualifier[iqual]=="cent20to30") phenixprotonyield+=wb*value;
						if(varname=="PHENIX_SPECTRA_PROTON_YIELD" && qualifier[iqual]=="cent30to40") phenixprotonyield+=wb*value;
						
						/** V2 */
						if(varname=="STAR_V2_PION_PTWEIGHT" && qualifier[iqual]=="cent0to10") pionv2+=wb*value;
						if(varname=="STAR_V2_PION_PTWEIGHT" && qualifier[iqual]=="cent10to20") pionv2+=wb*value;
						if(varname=="STAR_V2_PION_PTWEIGHT" && qualifier[iqual]=="cent20to30") pionv2+=wb*value;
						if(varname=="STAR_V2_PION_PTWEIGHT" && qualifier[iqual]=="cent30to40") pionv2+=wb*value;
						
						if(varname=="STAR_V2_KAON_PTWEIGHT" && qualifier[iqual]=="cent0to10") kaonv2+=wb*value;
						if(varname=="STAR_V2_KAON_PTWEIGHT" && qualifier[iqual]=="cent10to20") kaonv2+=wb*value;
						if(varname=="STAR_V2_KAON_PTWEIGHT" && qualifier[iqual]=="cent20to30") kaonv2+=wb*value;
						if(varname=="STAR_V2_KAON_PTWEIGHT" && qualifier[iqual]=="cent30to40") kaonv2+=wb*value;
						
						if(varname=="STAR_V2_PROTON_PTWEIGHT" && qualifier[iqual]=="cent0to10") protonv2+=wb*value;
						if(varname=="STAR_V2_PROTON_PTWEIGHT" && qualifier[iqual]=="cent10to20") protonv2+=wb*value;
						if(varname=="STAR_V2_PROTON_PTWEIGHT" && qualifier[iqual]=="cent20to30") protonv2+=wb*value;
						if(varname=="STAR_V2_PROTON_PTWEIGHT" && qualifier[iqual]=="cent30to40") protonv2+=wb*value;

					}
					
				}
				else{
					fgets(dummy,100,infile);
				}
			} while(!feof(infile));
			fclose(infile);
		}
		string command="mkdir -p  analysis/run"+irunstring+"/allb";
		system(command.c_str());
		outfilename="analysis/run"+irunstring+"/allb/results.dat";
		printf("writing to %s\n",outfilename.c_str());
		outfile=fopen(outfilename.c_str(),"w");
		fprintf(outfile,"double PIONV2 %g 0.005\n",pionv2);
		fprintf(outfile,"double KAONV2 %g 0.005\n",kaonv2);
		fprintf(outfile,"double PROTONV2 %g 0.005\n",protonv2);
		fprintf(outfile,"double PHENIX_PION_YIELD %g %g\n",phenixpionyield,0.1*phenixpionyield);
		fprintf(outfile,"double PHENIX_KAON_YIELD %g %g\n",phenixkaonyield,0.1*phenixkaonyield);
		fprintf(outfile,"double PHENIX_PROTON_YIELD %g %g\n",phenixprotonyield,0.1*phenixprotonyield);
		fprintf(outfile,"double PHENIX_PION_SPECTRA_PCAVAL %g %g\n",phenixpionpcaval,0.15*phenixpionpcaval);
		fprintf(outfile,"double PHENIX_KAON_SPECTRA_PCAVAL %g %g\n",phenixkaonpcaval,0.15*phenixkaonpcaval);
		fprintf(outfile,"double PHENIX_PROTON_SPECTRA_PCAVAL %g %g\n",phenixprotonpcaval,0.15*phenixprotonpcaval);
		fprintf(outfile,"double ROUT %g %g\n",rout,0.07*rout);
		fprintf(outfile,"double RSIDE %g %g\n",rside,0.07*rside);
		fprintf(outfile,"double RLONG %g %g\n",rlong,0.07*rlong);
		fclose(outfile);
	}
}


#endif
