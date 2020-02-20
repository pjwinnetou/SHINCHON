#include <stdlib.h>
#include <omp.h>
#include "b3d.h"
#include "qualifier.h"
#include <sstream>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]){
	if (argc != 2) 
	  {
	    printf("Usage: b3d run_name\n");
	    exit(-1);
	  }
	
	//This variable is set by the PBS file's environment variable OMP_NUM_THREADS                                                
	int num_threads = omp_get_max_threads();
	cout << num_threads << endl;

	double dnchdy=0.0;
	int nevents_written = 0; // The number of events written to an output file.
	int nparts;
	int iqual,neventsmax;
	string run_name=argv[1];
	CB3D *b3d[num_threads];
	CQualifiers qualifiers[num_threads];
	// CQualifiers qualifiers;
	// qualifiers.Read("qualifiers.dat");
	//printf("checkpt 0\n");
	for (int ithread = 0; ithread < num_threads; ithread++) {
		qualifiers[ithread].Read("qualifiers.dat");
		b3d[ithread] = new CB3D(run_name, ithread);
	}
	neventsmax=parameter::getI(b3d[0]->parmap,"B3D_NEVENTSMAX",10);
	//b3d->randy->reset(-time(NULL));	
	//printf("checkpt 1\n");

	for(int iqual=0;iqual<qualifiers[0].nqualifiers;iqual++) {
		dnchdy=0;
		nevents_written = 0;
		b3d[0]->SetQualifier(qualifiers[0].qualifier[iqual]);
		qualifiers[0].SetPars(&(b3d[0]->parmap), iqual);
		b3d[0]->hydrotob3d->ReadInputPR();
		for (int ithread = 1; ithread < num_threads; ithread++) {
 			//printf("iqual %d, ithread %d, checkpt 0\n", iqual, ithread);
			/*
			stringstream ss;
			ss.fill('0');
			ss.width(10);
			ss.setf(std::ios::right, std::ios::adjustfield);
			ss << ithread;
			string ithread_str = ss.str();
			cout << "ithread = " << ithread_str << ".\n";
			qualifiers[ithread].qualifier[iqual] += "_" + ithread_str;
			*/
			//b3d[ithread]->SetQualifier(qualifiers[ithread].qualifier[iqual]);
			b3d[ithread]->qualifier = qualifiers[ithread].qualifier[iqual];
			b3d[ithread]->h5outfile = b3d[0]->h5outfile;
			b3d[ithread]->h5infile = b3d[0]->h5infile;
 			//printf("iqual %d, ithread %d, checkpt 1\n", iqual, ithread);
			qualifiers[ithread].SetPars(&(b3d[ithread]->parmap),iqual);
 			//printf("iqual %d, ithread %d, checkpt 2\n", iqual, ithread);
			b3d[ithread]->hydrotob3d->ReadInputPR();
 			//printf("iqual %d, ithread %d, checkpt 3\n", iqual, ithread);
		}
		omp_set_num_threads(num_threads);
#pragma omp parallel for \
	default(shared) \
	private(nparts) \
	firstprivate(num_threads, neventsmax, iqual)
		for (int ithread = 0; ithread < num_threads; ithread++) {
			int neventsthread = neventsmax / num_threads; // Assign to each thread a number of events equal to neventsmax/num_threads.
			int startevent = 0; // Stores the event the thread should start on.
			
			// Find the starting event for this thread.
			for (int i = 0; i < ithread; i++) {
				startevent += neventsthread;
				if (i < (neventsmax % num_threads)) {
					startevent++;
				}
			}

			if (ithread < (neventsmax % num_threads)) { // If there are remaining events to be assigned, give 1 of the remaining events to each thread, starting at thread 0, until all events are assigned.
				neventsthread++;
			}
			int endevent = startevent + neventsthread; // This is the event the thread should end on.
			// Now, loop through all events assigned to a given thread.
			for(int ievent=startevent;ievent<endevent;ievent++){
 				//printf("iqual %d, ithread %d, ievent %d, checkpt 0\n", iqual, ithread, ievent);
				nparts=b3d[ithread]->hydrotob3d->MakeEventPR();
 				//printf("iqual %d, ithread %d, ievent %d, checkpt 1\n", iqual, ithread, ievent);
				b3d[ithread]->PerformAllActions();
 				//printf("iqual %d, ithread %d, ievent %d, checkpt 2\n", iqual, ithread, ievent);
#pragma omp critical
				{
					// b3d[ithread]->ievent_write = nevents_written;
					b3d[ithread]->ievent_write = ievent;
      		dnchdy+=b3d[ithread]->WriteDataH5(); // Write to the output file.
					nevents_written++; // Update total number of events written.
 					//printf("iqual %d, ithread %d, ievent %d, checkpt 3\n", iqual, ithread, ievent);
      		printf("##### finished event %d ##### dNch/dy=%g #####\n",b3d[ithread]->ievent_write,double(dnchdy)/double(nevents_written));
					printf("nscatter=%g, ninelastic=%g, nmerges=%g, ncellexits=%g\n",double(b3d[ithread]->nscatter),double(b3d[ithread]->ninelastic),double(b3d[ithread]->nmerge),double(b3d[ithread]->nexit));
      		cout << "Elastic Scattering: " << b3d[ithread]->nscatter << " Inelastic Scatter: " << b3d[ithread]->ninelastic << " Merges: " << b3d[ithread]->nmerge << endl;
				}
				//cout << "Total collisions: " << b3d->ncollision << endl;
				// printf("nactivate=%lld, nexit/N=%g, nscatter/N=%g, nmerge/N=%g, ndecay/N=%g,\n", b3d->nactivate,double(b3d->nexit)/N,double(b3d->nscatter)/N,double(b3d->nmerge)/N,double(b3d->ndecay)/N);
      	// printf("ncheck=%lld, nscatter=%lld, npass=%lld\n",b3d->ncheck,b3d->nscatter,b3d->npass);
			}
		}
	}
  
	return 0;
}
