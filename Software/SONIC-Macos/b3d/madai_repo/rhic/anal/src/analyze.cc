#ifndef __ANALYZE_CC__
#define __ANALYZE_CC__

#include "analyze.h"

CAnalyze::CAnalyze(string &run_name_set){
	run_name=run_name_set;
	string parsfilename="parameters/"+run_name+"/fixed.param";
	parameter::ReadParsFromFile(parmap,parsfilename);
	parsfilename="parameters/"+run_name+"/stats.param";
	parameter::ReadParsFromFile(parmap,parsfilename);

	//parameter::ReadParsFromFile(parmap,"parameters/analysis.param");

	ptype=new CompType(sizeof(CPartH5));
	
	neventsmax=parameter::getI(parmap,"B3D_NEVENTSMAX",40000);
	npartsmax=parameter::getI(parmap,"B3D_NPARTSMAX",3000);
	nsample=parameter::getI(parmap,"B3D_NSAMPLE",1);
	npartsmax*=nsample;
	
	bPrintSpectra = parameter::getB(parmap,"ANAL_PRINTSPEC",false);
	bPrintV2 = parameter::getB(parmap,"ANAL_PRINTV2",false);
	
	ptype->insertMember("listid", HOFFSET(CPartH5,listid), PredType::NATIVE_INT);
	ptype->insertMember("ID", HOFFSET(CPartH5,ID), PredType::NATIVE_INT);
	ptype->insertMember("x", HOFFSET(CPartH5,x), PredType::NATIVE_DOUBLE);
	ptype->insertMember("y", HOFFSET(CPartH5,y), PredType::NATIVE_DOUBLE);
	ptype->insertMember("eta", HOFFSET(CPartH5,eta), PredType::NATIVE_DOUBLE);
	ptype->insertMember("tau", HOFFSET(CPartH5,tau), PredType::NATIVE_DOUBLE);
	ptype->insertMember("px", HOFFSET(CPartH5,px), PredType::NATIVE_DOUBLE);
	ptype->insertMember("py", HOFFSET(CPartH5,py), PredType::NATIVE_DOUBLE);
	ptype->insertMember("rapidity", HOFFSET(CPartH5,rapidity), PredType::NATIVE_DOUBLE);
	ptype->insertMember("mass", HOFFSET(CPartH5,mass), PredType::NATIVE_DOUBLE);
	
	partH5=new CPartH5[npartsmax];
	part=new CPart[npartsmax];
	STAR_ACCEPTANCE=parameter::getB(parmap,"B3D_STAR_ACCEPTANCE","false");
	CALCGARRAYS=false;
	output=NULL;
	
	phenix_pion_yield_min=325;
	phenix_pion_yield_max=500;
	
	reslist=new CResList();
	b3d=new CB3D();
	reslist->parmap=&parmap;
	b3d->reslist=reslist;
	b3d->BJORKEN=false;
	b3d->ETAMAX=1.0E20;
	reslist->b3d=b3d;
}

void CAnalyze::SetQualifier(string qualifiername){
	if(output!=NULL) fclose(output);
	input_dataroot="output/"+run_name+"/"+qualifiername;
	output_dataroot="analysis/"+run_name+"/"+qualifiername;
	string command="mkdir -p "+output_dataroot;
	system(command.c_str());
	h5_infilename=input_dataroot+"/b3d.h5";
	output_filename=output_dataroot+"/results.dat";
	printf("h5_infilename=%s\n",h5_infilename.c_str());
	printf("output_filename=%s\n",output_filename.c_str());
	output=fopen(output_filename.c_str(),"w");
	qualifier=qualifiername;
	detailsdirname="analysis/"+run_name+"/"+qualifier+"/details";
	bool WRITEDETAILS=parameter::getB(parmap,string("WRITEDETAILS"),false);
	if(WRITEDETAILS){
		string command="mkdir -p "+detailsdirname;
		system(command.c_str());
	}
	//if(h5file!=NULL) delete h5file;
	//h5_infilename=parameter::getS(parmap,"B3D_H5_INFILENAME","b3d.h5");	
	//vizfilename=parameter::getS(parmap,"B3D_VIZ_INFILENAME","b3dviz.h5");
	//vizfile = new H5File(infilename,H5F_ACC_RDONLY);
}

CAnalyze::~CAnalyze(){
	if(output!=NULL) fclose(output);
}

int CAnalyze::ReadDataH5(int ievent){
	int nparts,ipart=0;
	char eventno[20];
	H5D_space_status_t status;
	sprintf(eventno,"event%d",ievent);
	hsize_t dim[1];
	DataSet *dataset = new DataSet (h5file->openDataSet(eventno));
	//dataset->getSpaceStatus(status);
	//hsize_t datasetsize=dataset->getStorageSize();
	//printf("ievent=%d, status=%d, size=%d\n",ievent,int(status),int(datasetsize));
	DataSpace filespace = dataset->getSpace();
	int rank=filespace.getSimpleExtentDims(dim);
	nparts=dim[0];
	part->b3d=b3d;
	if(nparts>npartsmax){
		printf("Increase B3D_NPARTSMAX, nparts=%d, npartsmax=%d\n",nparts,npartsmax);
		exit(1);
	}
	dataset->read(partH5,*ptype);
	delete dataset;
	CPart *dptr;
	CPartH5 *ph5;
	
	CPart *mother,**daughter=new CPart *[5];
	CResInfo **daughterresinfo=new CResInfo *[5];
	
	int nbodies,ibody,ntry,nparts0=nparts;
	double mtot,mothermass;
	while(ipart<nparts){
		mother=&part[ipart];
		if(ipart<nparts0) part[ipart].Copy(&partH5[ipart]);
		if(fabs(mother->mass-mother->resinfo->mass)>1.0E-3){
			printf("in ReadDataH5, mother mass out of whack, %g!=%g\n",mother->mass,mother->resinfo->mass);
			mother->Print();
			exit(1);
		}
		mothermass=mother->GetMass();
		nbodies=0;
		if(mother->resinfo->decay){
				//printf("this mother won't last, ID=%d\n",mother->resinfo->code);
			ntry=0;
			do{
				mtot=0.0;
				if(ntry<25) mother->resinfo->DecayGetResInfoptr(nbodies,daughterresinfo);
				for(ibody=0;ibody<nbodies;ibody++){
					mtot+=daughterresinfo[ibody]->mass;
				}
				if(ntry>25){
					printf("FATAL: action_perform_decay, ntry too big, mothermass=%g\n",mother->GetMass());
					mother->Print();
					exit(1);
				}
				ntry++;
			}while(mtot>mothermass);
			for(ibody=0;ibody<nbodies;ibody++){
				daughter[ibody]=&part[nparts+ibody];
				daughter[ibody]->resinfo=daughterresinfo[ibody];
			}
			b3d->Decay(mother,nbodies,daughter);
			mother->Copy(daughter[nbodies-1]);
			for(ibody=0;ibody<nbodies;ibody++){
				dptr=daughter[ibody];
				if(fabs(dptr->mass-dptr->resinfo->mass)>1.0E-3){
					printf("in ReadDataH5, masses out of whack, %g!=%g\n",dptr->mass,dptr->resinfo->mass);
					dptr->Print();
					exit(1);
				}
				if(ibody<nbodies-1)ph5=&partH5[nparts+ibody];
				else ph5=&partH5[ipart];
				ph5->px=dptr->p[1];
				ph5->py=dptr->p[2];
				ph5->rapidity=dptr->y;
				ph5->eta=dptr->eta;
				ph5->mass=dptr->mass;
				ph5->x=dptr->r[1];
				ph5->y=dptr->r[2];
				ph5->tau=dptr->tau0;
				ph5->ID=dptr->resinfo->code;
			}
			nparts+=nbodies-1;
		}
		else ipart+=1;
	}
		//	printf("nparts=%d, nparts0=%d\n",nparts,nparts0);

	delete [] daughterresinfo;
	delete [] daughter;
	
	
	//printf("READ IN %d PARTS\n",nparts);
	return nparts;
}
/**
int CAnalyze::ReadVizData(double tau){
	int nparts,ipart,rank;
	string infilename=input_dataroot+"/"+qualifier+"/"+vizfilename;
	printf("will read %s for viz data\n",infilename.c_str());
	vizfile = new H5File(infilename,H5F_ACC_RDONLY);

	char setname[40];
	H5D_space_status_t status;
	sprintf(setname,"px_tau%g",tau);
	hsize_t pdim[1];
	DataSet *dataset = new DataSet (vizfile->openDataSet(setname));
	//dataset->getSpaceStatus(status);
	//hsize_t datasetsize=dataset->getStorageSize();
	//printf("ievent=%d, status=%d, size=%d\n",ievent,int(status),int(datasetsize));
	DataSpace filespace = dataset->getSpace();
	rank=filespace.getSimpleExtentDims(pdim);
	nparts=pdim[0];
	printf("For px: nparts=%d\n",nparts);
	double *px=new double[nparts];
	if(nparts>npartsmax){
		printf("Increase NPARTSMAX, nparts=%d\n",nparts);
		exit(1);
	}
	dataset->read(px,PredType::NATIVE_DOUBLE);
	delete dataset;
	
	
	printf("------------------\n");
	sprintf(setname,"xyz_tau%g",tau);
	DataSet *xyzdataset = new DataSet (vizfile->openDataSet(setname));
	//dataset->getSpaceStatus(status);
	//hsize_t datasetsize=dataset->getStorageSize();
	//printf("ievent=%d, status=%d, size=%d\n",ievent,int(status),int(datasetsize));
	hsize_t xyzdim[]={nparts,3};
	DataSpace xyzfilespace = xyzdataset->getSpace();
	rank=xyzfilespace.getSimpleExtentDims(xyzdim);
	printf("xyz rank=%d\n",rank);
	nparts=xyzdim[0];
	printf("For Reading set %s, nparts=%d, dimension=%d\n",setname,nparts,int(xyzdim[1]));
	double (*xyz)[3]=new double[nparts][3];
	if(nparts>npartsmax){
		printf("Increase NPARTSMAX, nparts=%d\n",nparts);
		exit(1);
	}
	xyzdataset->read(xyz,PredType::NATIVE_DOUBLE);
	delete xyzdataset;
	
	
	printf("READ IN %d PARTS\n",nparts);
	for(ipart=0;ipart<nparts;ipart++){
		printf("ipart=%d: px=%g, xyz=(%g,%g,%g)\n",ipart,px[ipart],xyz[ipart][0],xyz[ipart][1],xyz[ipart][2]);
	}
	
	delete vizfile;
	delete [] px;
	delete [] xyz;
	return nparts;
}
*/

double CAnalyze::legendre(int ell,double x){
	if(ell==0) return 1.0;
	else if(ell==1) return x;
	else if(ell==2) return 0.5*(3.0*x*x-1.0);
	else if(ell==3) return 0.5*(5.0*x*x*x-3.0*x);
	else if(ell==4) return 0.125*(35.0*x*x*x*x-30.0*x*x+3.0);
	else if(ell==5) return 0.125*(63.0*pow(x,5)-70.0*x*x*x+15.0*x);
	else if(ell==6) return 0.0625*(231.0*pow(x,6)-315.0*pow(x,4)+105.0*x*x-5.0);
	else{
		printf("calculating CAnalyze::legendre(int ell,double x) for ell=%d (only goes to 6)\n",ell);
		return 0.0;
	}
}

#endif
