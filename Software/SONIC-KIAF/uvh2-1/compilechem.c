/* functions.c            functions for reso.c main program calculations 
*   Josef Sollfrank                  Nov. .98             */

// Commented and adjusted by Evan Frodermann, Aug 2005

#include <iostream>
#include <fstream>
using namespace std;

 #include	"readin.h" 
  double stablecp[500][319];
  
  
//*********************************************************************************


double getchem(int t, int i) 
{
  if (particle[i].baryon == -1) i--;
  if (particle[i].stable)
  {
    int stableindex = stablelookup[i];
    return stablecp[t][stableindex];
  }
  else
  {
    double chem = 0;
    for (int channel=0; channel< particle[i].decays; channel++)
    { 
	for (int daughter=0; daughter < abs(decay[i][channel].numpart); daughter++)
	{

	  int daughtval = decay[i][channel].part[daughter];

	  int daughterindex = partid[MHALF + daughtval];
	  double branch = decay[i][channel].branch/abs(decay[i][channel].numpart);
	  chem += branch*getchem(t,daughterindex);
// 	  chem += 1;
	}
    }
    return chem;
  }
     
}


/**************************************************************************
*									  *
*   readin() 								  *
*									  *
*   reads in the particle data file and stores the datas in the arrays	  *
*   particle.* and decay.* and fills up the rest of data (antibaryons)	  *
***********************************************************************   *
**************************************************************************/

// void   readin(filename, particlemax, decaymax)
// void readin() {
int main() {


//      char filename[FILEDIM];
     char filename[] = { 'p', 'd', 'g', '0', '5', '.', 'd', 'a', 't','\0' };
     int ntemp = 500; //number of temperature steps in chemical potential filw
//      string filename = "pdg_weak.dat";
//      int particlemax, decaymax;

//        {
	 int i=0, j=0, k, h;
	 FILE	*dat;
	 double dummy1;
	 
	 for(k=0;k<MAXINTV;k++) 
	   partid[k] = -1; 

	 dat = fopen(filename,"r");
	 if(dat == NULL){
	   printf(" NO file: %s  available ! \n", filename);
	   printf(" GOOD BYE AND HAVE A NICE DAY! \n");
	   exit(0);
	}
	 // Read in the particle data from the specified resonance table 
	 // Save the data is the structure particle[pn]

	while(fscanf(dat,"%li%s%lf%lf%i%i%i%i%i%i%i%i", &particle[i].monval,
	       particle[i].name, &particle[i].mass, &particle[i].width,
	       &particle[i].gspin, &particle[i].baryon,
	       &particle[i].strange, &particle[i].charm,
	       &particle[i].bottom, &particle[i].gisospin,
	       &particle[i].charge, &particle[i].decays)==12)
	  { 
// cout << i << endl;
	    //in resoweak
	    // Read in the unused portion of the data file 
	    //fscanf(dat,"%lf%lf%lf", &dummy1, &dummy1, &dummy1);

	    //dummy1=0.0;

	    //in Pasi's file: no such numbers...
	    partid[MHALF + particle[i].monval] = i;
	    particle[i].stable = 0;


// 	    printf("%i   %li %s %lf %lf %i %i %i %i %i %i %i %i\n", i, particle[i].monval,
// 	    	      particle[i].name, particle[i].mass, particle[i].width,
// 	  	      particle[i].gspin, particle[i].baryon,
// 	    	      particle[i].strange, particle[i].charm,
// 	          particle[i].bottom, particle[i].gisospin,
// 	    	   particle[i].charge, particle[i].decays);


	    /* read in the decays */
	    // These decays are saved in a seperate data set, decay[i].
	   for(k=0;k<particle[i].decays;k++) 
	     {
	       h=fscanf(dat,"%i%i%lf%i%i%i%i%i",
			&decay[i][k].reso, &decay[i][k].numpart, &decay[i][k].branch, 
			&decay[i][k].part[0], &decay[i][k].part[1], &decay[i][k].part[2],
			&decay[i][k].part[3], &decay[i][k].part[4]);
	       //  printf("      %i %i %lf %i %i %i %i %i\n", decay[j].reso,
	       //	      decay[j].numpart, decay[j].branch, decay[j].part[0],
	       //	      decay[j].part[1], decay[j].part[2],decay[j].part[3],
	       //      decay[j].part[4]);
	      
	        if (h != 8) {
	          printf("Error in scanf decay \n");
                  printf(" GOOD BYE AND HAVE A NICE DAY! \n");
	          exit(0);
	        }
		if (decay[i][k].numpart == 1) {
		  particle[i].stable = 1;
		  stablelookuprev[j] = i;
		  stablelookup[i] = j;
// 		  cout << "Stable particle: " << particle[i].name << endl;
		  j++;
		}
// 		j++; // Add one to the decay counting variable "j"
	    }
	   //printf("\n");
	   
	   

	   /* setting of additional parameters */
	   
	    if (particle[i].baryon == 1)
	      {
		i++;// If the particle is a baryon, add a particle for the anti-baryon
		    // Add one to the counting variable "i" for the number of particles for the antibaryon
		particle[i].monval = -particle[i-1].monval;
		strcpy(particle[i].name, "  Anti-");
		strncat(particle[i].name, particle[i-1].name, 18);
		particle[i].mass     =  particle[i-1].mass;
		particle[i].width    =  particle[i-1].width;
		particle[i].gspin    =  particle[i-1].gspin;
		particle[i].baryon   = -particle[i-1].baryon;
		particle[i].strange  = -particle[i-1].strange;
		particle[i].charm    = -particle[i-1].charm;
		particle[i].bottom   = -particle[i-1].bottom;
		particle[i].gisospin =  particle[i-1].gisospin;
		particle[i].charge   = -particle[i-1].charge;
		particle[i].decays   = particle[i-1].decays;
	 	partid[MHALF + particle[i].monval] = i;
	        particle[i].stable =  particle[i-1].stable;
		if (particle[i].stable == 1) {
		  stablelookuprev[j] = i;
		  stablelookup[i] = j;
// 		  cout << "Stable particle: " << particle[i].name << endl;
		  j++;
		}
	    }
	    
	    i++; // Add one to the counting variable "i" for the meson/baryon
	  }
	fclose(dat);

	//printf("last particle: %i %s \n",particle[i-1].monval,
        //                  particle[i-1].name);
	//printf("# particle %5i; # decays %5i; \n\n",i,j);

// 	particlemax = i;   // Set the maxparticle variable to be the value of "i"
// 	decaymax  = j;     // Set the maxdecays variable to be the value of "j"
int particlemax = i;
int stableparticlemax = j;
cout << j << " stable particles\n";
cout << particlemax << " total resonances\n" << endl;

        if((i) > NUMPARTICLE){
          printf("Array for particles to small!!\n");
          printf(" GOOD BYE AND HAVE A NICE DAY! \n");
	  exit(0);
	} 

  fstream stablecpfile;
  stablecpfile.open("pasichem.dat", ios::in);
  for (int i=0;i<ntemp;i++)
    for (int j=0;j<stableparticlemax;j++)
      {
	stablecpfile >> stablecp[i][j];
      }
      cout << stablecp[0][0];
  stablecpfile.close();
  double allchemtable[particlemax][ntemp];
  fstream chemout;
  chemout.open("allchem.dat", ios::out);

    
  for (j=0; j<ntemp; j++)
//   for (j=1; j<2; j++)
  {
    for (i=0; i<particlemax; i++)
//     for (i=0; i<1; i++) 
    {

// 	}
      double chem = getchem(j, i);
      chemout << chem << "\t";
//       cout << "chem = " << chem << endl;
    }
    chemout << endl;
  }
  chemout.close();
  return 0;
	
      }

// int main() {
//   fstream stablecpfile;
//   readin();
//   stablecpfile.open("pasichem.dat", ios::in);
//   double stablecp[500][319];
//   for (int i=0;i<500;i++)
//     for (int j=0;j<319;j++)
//       {
// 	
// 	stablecpfile >> stablecp[i][j];
// 	
//       }
//       cout << stablecp[0][1] << endl;
//   return 0;
// }