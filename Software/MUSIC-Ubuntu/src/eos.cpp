// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

// Copyright 2011 @ Bjoern Schenke, Sangyong Jeon, and Charles Gale
#include<fstream>
#include<iomanip>

#include "util.h"
#include "eos.h"
#include "data.h"
using namespace std;

#define cs2 (1.0/3.0)

EOS::EOS(InitData *para_in) {
    parameters_ptr = para_in;
    util = new Util;
    initialize_eos();
    eps_max = 1e5;  // [1/fm^4]
}

// destructor
EOS::~EOS() {
    if (parameters_ptr->whichEOS == 1) {
        util->mtx_free(pressure1, NBNP1 + 1, NEPP1 + 1);
        util->mtx_free(pressure2, NBNP2 + 1, NEPP2 + 1);
        util->mtx_free(temperature1, NBNP1 + 1, NEPP1 + 1);
        util->mtx_free(temperature2, NBNP2 + 1, NEPP2 + 1);
        util->mtx_free(mu1, NBNP1 + 1, NEPP1 + 1);
        util->mtx_free(mu2, NBNP2 + 1, NEPP2 + 1);
        util->mtx_free(cs2_1, NBNP1 + 1, NEPP1 + 1);
        util->mtx_free(cs2_2, NBNP2 + 1, NEPP2 + 1);
    } else if (parameters_ptr->whichEOS >= 2
                    && parameters_ptr->whichEOS <= 7) {
        util->mtx_free(pressure1, NBNP1 + 1, NEPP1 + 1);
        util->mtx_free(pressure2, NBNP2 + 1, NEPP2 + 1);
        util->mtx_free(pressure3, NBNP3 + 1, NEPP3 + 1);
        util->mtx_free(pressure4, NBNP4 + 1, NEPP4 + 1);
        util->mtx_free(pressure5, NBNP5 + 1, NEPP5 + 1);
        util->mtx_free(pressure6, NBNP6 + 1, NEPP6 + 1);
        util->mtx_free(pressure7, NBNP7 + 1, NEPP7 + 1);
        util->mtx_free(entropyDensity1, NBNP1 + 1, NEPP1 + 1);
        util->mtx_free(entropyDensity2, NBNP2 + 1, NEPP2 + 1);
        util->mtx_free(entropyDensity3, NBNP3 + 1, NEPP3 + 1);
        util->mtx_free(entropyDensity4, NBNP4 + 1, NEPP4 + 1);
        util->mtx_free(entropyDensity5, NBNP5 + 1, NEPP5 + 1);
        util->mtx_free(entropyDensity6, NBNP6 + 1, NEPP6 + 1);
        util->mtx_free(entropyDensity7, NBNP7 + 1, NEPP7 + 1);
        util->mtx_free(QGPfraction1, NBNP1 + 1, NEPP1 + 1);
        util->mtx_free(QGPfraction2, NBNP2 + 1, NEPP2 + 1);
        util->mtx_free(QGPfraction3, NBNP3 + 1, NEPP3 + 1);
        util->mtx_free(QGPfraction4, NBNP4 + 1, NEPP4 + 1);
        util->mtx_free(QGPfraction5, NBNP5 + 1, NEPP5 + 1);
        util->mtx_free(QGPfraction6, NBNP6 + 1, NEPP6 + 1);
        util->mtx_free(QGPfraction7, NBNP7 + 1, NEPP7 + 1);
        util->mtx_free(temperature1, NBNP1 + 1, NEPP1 + 1);
        util->mtx_free(temperature2, NBNP2 + 1, NEPP2 + 1);
        util->mtx_free(temperature3, NBNP3 + 1, NEPP3 + 1);
        util->mtx_free(temperature4, NBNP4 + 1, NEPP4 + 1);
        util->mtx_free(temperature5, NBNP5 + 1, NEPP5 + 1);
        util->mtx_free(temperature6, NBNP6 + 1, NEPP6 + 1);
        util->mtx_free(temperature7, NBNP7 + 1, NEPP7 + 1);
        util->mtx_free(cs2_1, NBNP1 + 1, NEPP1 + 1);
        util->mtx_free(cs2_2, NBNP2 + 1, NEPP2 + 1);
        util->mtx_free(cs2_3, NBNP3 + 1, NEPP3 + 1);
        util->mtx_free(cs2_4, NBNP4 + 1, NEPP4 + 1);
        util->mtx_free(cs2_5, NBNP5 + 1, NEPP5 + 1);
        util->mtx_free(cs2_6, NBNP6 + 1, NEPP6 + 1);
        util->mtx_free(cs2_7, NBNP7 + 1, NEPP7 + 1);
    } else if (parameters_ptr->whichEOS != 0) {
        music_message << "No EOS for whichEOS = " << parameters_ptr->whichEOS
             << ". Use EOS_to_use = 0 (ideal gas) 1 (AZHYDRO EOS-Q), "
             << "2 (s95p-v1), 3 (s95p-PCE150-v1), 4 (s95p-PCE155-v1), "
             << "5 (s95p-PCE160-v1), 6 (s95p-PCE165-v1),"
             << "7 (s95p-v1.2)";
        music_message.flush("error");
        exit(1);
    }
    delete util;
    if (parameters_ptr->whichEOS >= 2 && parameters_ptr->whichEOS <= 7) {
        gsl_interp_free(interp_s2e);
        gsl_interp_accel_free(accel_s2e);
    }
}

void EOS::checkForReadError(FILE *file, const char* name) {
    if (!(file)) {
        music_message << "file " << name << "not found. Exiting...";
        music_message.flush("error");
        exit(0);
    }
}

void EOS::initialize_eos() {
    if (parameters_ptr->whichEOS == 0) {
        music_message.info("Using the ideal gas EOS");
        init_eos0();
    } else if (parameters_ptr->whichEOS == 1) {
        music_message.info("Using EOS-Q from AZHYDRO");
        init_eos();
    } else if (parameters_ptr->whichEOS == 2) {
        music_message.info("Using lattice EOS from Huovinen/Petreczky");
        init_eos2();
    } else if (parameters_ptr->whichEOS == 3) {
        music_message << "Using lattice EOS from Huovinen/Petreczky with "
                      << "partial chemical equilibrium (PCE) "
                      << "chem. f.o. at 150 MeV";
        music_message.flush("info");
        init_eos3(1);
    } else if (parameters_ptr->whichEOS == 4) {
        music_message << "Using lattice EOS from Huovinen/Petreczky with "
                      << "partial chemical equilibrium (PCE) "
                      << "chem. f.o. at 155 MeV";
        music_message.flush("info");
        init_eos3(2);
    } else if (parameters_ptr->whichEOS == 5) {
        music_message << "Using lattice EOS from Huovinen/Petreczky with "
                      << "partial chemical equilibrium (PCE) "
                      << "chem. f.o. at 160 MeV";
        music_message.flush("info");
        init_eos3(3);
    } else if (parameters_ptr->whichEOS == 6) {
        music_message << "Using lattice EOS from Huovinen/Petreczky with "
                      << "partial chemical equilibrium (PCE) "
                      << "chem. f.o. at 165 MeV";
        music_message.flush("info");
        init_eos3(4);
    } else if (parameters_ptr->whichEOS == 7) {
        music_message << "Using lattice EOS from Huovinen/Petreczky s95p-v1.2"
                      << "(for UrQMD)";
        music_message.flush("info");
        init_eos7();
    } else {
        music_message << "No EOS for whichEOS = " << parameters_ptr->whichEOS
             << ". Use EOS_to_use = 0 (ideal gas) 1 (AZHYDRO EOS-Q), "
             << "2 (s95p-v1), 3 (s95p-PCE150-v1), 4 (s95p-PCE155-v1), "
             << "5 (s95p-PCE160-v1), 6 (s95p-PCE165-v1),"
             << "7 (s95p-v1.2)";
        music_message.flush("error");
        exit(1);
    }
    
    if (whichEOS >= 2 && whichEOS < 10) {
        eps_max = (EPP7 + deltaEPP7*(NEPP7-1))/hbarc;  // [1/fm^4]
    }
}

void EOS::init_eos0() {
    whichEOS = 0;
}

void EOS::init_eos() {
// read the azhydro pressure, temperature, and 
// baryon chemical potential from file
    whichEOS = 1;
    music_message.info("reading EOS...");
    int i, j;
    FILE *eos_p1, *eos_p2;
    FILE *eos_T1, *eos_T2;
    FILE *eos_mu1, *eos_mu2;
    const char* EOSPATH = "HYDROPROGRAMPATH";
    char * pre_envPath= getenv(EOSPATH);
    std::string envPath;
    std::string eos_p1_name;
    std::string eos_p2_name;
    std::string eos_T1_name;
    std::string eos_T2_name;
    std::string eos_mu1_name;
    std::string eos_mu2_name;
      
    if (pre_envPath == 0) {
	    envPath=".";
    }
    else {
	    envPath=pre_envPath;
    }
    
    music_message << "from path " << envPath.c_str() << "/EOS";
    music_message.flush("info");
    eos_p1_name=envPath+"/EOS/EOS-Q/aa1_p.dat";
    eos_p2_name=envPath+"/EOS/EOS-Q/aa2_p.dat";
    eos_T1_name=envPath+"/EOS/EOS-Q/aa1_t.dat";
    eos_T2_name=envPath+"/EOS/EOS-Q/aa2_t.dat";
    eos_mu1_name=envPath+"/EOS/EOS-Q/aa1_mb.dat";
    eos_mu2_name=envPath+"/EOS/EOS-Q/aa2_mb.dat";
  
  eos_p1 = fopen(eos_p1_name.c_str(), "r");
  eos_p2 = fopen(eos_p2_name.c_str(), "r");
  eos_T1 = fopen(eos_T1_name.c_str(), "r");
  eos_T2 = fopen(eos_T2_name.c_str(), "r");
  eos_mu1 = fopen(eos_mu1_name.c_str(), "r");
  eos_mu2 = fopen(eos_mu2_name.c_str(), "r");
  
  checkForReadError(eos_p1,eos_p1_name.c_str());
  checkForReadError(eos_p2,eos_p2_name.c_str());
  checkForReadError(eos_T1,eos_T1_name.c_str());
  checkForReadError(eos_T2,eos_T2_name.c_str());
  checkForReadError(eos_mu1,eos_mu1_name.c_str());
  checkForReadError(eos_mu2,eos_mu2_name.c_str());
 
  int fscan_status = 0;
  // read the first two lines:
  // first value of rhob, first value of epsilon
  // deltaRhob, deltaEpsilon, number of rhob steps, number of epsilon steps
  fscan_status = fscanf(eos_p1, "%lf %lf", &BNP1, &EPP1);
  fscan_status = fscanf(eos_p1, "%lf %lf %d %d", &deltaBNP1, &deltaEPP1,
                        &NBNP1, &NEPP1);
  fscan_status = fscanf(eos_p2, "%lf %lf", &BNP2, &EPP2);
  fscan_status = fscanf(eos_p2, "%lf %lf %d %d", &deltaBNP2, &deltaEPP2,
                        &NBNP2, &NEPP2);
  fscan_status = fscanf(eos_T1, "%lf %lf", &BNP1, &EPP1);
  fscan_status = fscanf(eos_T1, "%lf %lf %d %d", &deltaBNP1, &deltaEPP1,
                        &NBNP1, &NEPP1);
  fscan_status = fscanf(eos_T2, "%lf %lf", &BNP2, &EPP2);
  fscan_status = fscanf(eos_T2, "%lf %lf %d %d", &deltaBNP2, &deltaEPP2,
                        &NBNP2, &NEPP2);
  fscan_status = fscanf(eos_mu1, "%lf %lf", &BNP1, &EPP1);
  fscan_status = fscanf(eos_mu1, "%lf %lf %d %d", &deltaBNP1, &deltaEPP1,
                        &NBNP1, &NEPP1);
  fscan_status = fscanf(eos_mu2, "%lf %lf", &BNP2, &EPP2);
  fscan_status = fscanf(eos_mu2, "%lf %lf %d %d", &deltaBNP2, &deltaEPP2,
                        &NBNP2, &NEPP2);
  if (fscan_status == 0 && parameters_ptr->echo_level > 9) {
      music_message << "fscanf_status = " << fscan_status;
      music_message.flush("warning");
  }
 
  // allocate memory for pressure arrays
  pressure1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  pressure2=util->mtx_malloc(NBNP2+1,NEPP2+1);

  // allocate memory for temperature arrays
  // we assume here that all files have the same structure (i.e., same EPP1, deltaEPP1, etc....)
  temperature1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  temperature2=util->mtx_malloc(NBNP2+1,NEPP2+1);

  // allocate memory for baryon chemical potential arrays
  mu1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  mu2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  
  // allocate memory for baryon chemical potential arrays
  cs2_1 = util->mtx_malloc(NBNP1+1, NEPP1+1);
  cs2_2 = util->mtx_malloc(NBNP2+1, NEPP2+1);

  // read pressure, temperature and chemical potential values
  for(j=0;j<=NEPP1;j++)
    for(i=0;i<=NBNP1;i++)
      {
          fscan_status = fscanf(eos_p1,"%lf",&pressure1[i][j]);
          fscan_status = fscanf(eos_T1,"%lf",&temperature1[i][j]);
          fscan_status = fscanf(eos_mu1,"%lf",&mu1[i][j]);
      }

  for(j=0;j<=NEPP2;j++)
    for(i=0;i<=NBNP2;i++)
      {
          fscan_status = fscanf(eos_p2,"%lf",&pressure2[i][j]);
          fscan_status = fscanf(eos_T2,"%lf",&temperature2[i][j]);
          fscan_status = fscanf(eos_mu2,"%lf",&mu2[i][j]);
      }

  build_velocity_of_sound_sq_matrix();
  fclose(eos_p1);
  fclose(eos_p2);
  fclose(eos_T1);
  fclose(eos_T2);
  fclose(eos_mu1);
  fclose(eos_mu2);
  music_message.info("done reading");
}

void EOS::init_eos2()
{
  // read the lattice EOS pressure, temperature, and 
  // baryon chemical potential from file
  music_message.info("reading EOS...");
  whichEOS = 2; 
//   int bytes_read;
  int i, j;
  FILE *eos_d1, *eos_d2, *eos_d3, *eos_d4, *eos_d5, *eos_d6, *eos_d7;
  FILE *eos_T1, *eos_T2, *eos_T3, *eos_T4, *eos_T5, *eos_T6, *eos_T7;
  const char* EOSPATH = "HYDROPROGRAMPATH";
  char * pre_envPath= getenv(EOSPATH);
  std::string envPath;
  std::string eos_d1_name;
  std::string eos_d2_name;
  std::string eos_d3_name;
  std::string eos_d4_name;
  std::string eos_d5_name;
  std::string eos_d6_name;
  std::string eos_d7_name;
  std::string eos_T1_name;
  std::string eos_T2_name;
  std::string eos_T3_name;
  std::string eos_T4_name;
  std::string eos_T5_name;
  std::string eos_T6_name;
  std::string eos_T7_name;
  double eps, baryonDensity; //dummies for now

  if (pre_envPath == 0) {
          envPath=".";
  }
  else {
	envPath=pre_envPath;
  }

  music_message << "from path " << envPath.c_str() << "/EOS/s95p-v1";
  music_message.flush("info");

  eos_d1_name=envPath+"/EOS/s95p-v1/s95p-v1_dens1.dat";
  eos_d2_name=envPath+"/EOS/s95p-v1/s95p-v1_dens2.dat";
  eos_d3_name=envPath+"/EOS/s95p-v1/s95p-v1_dens3.dat";
  eos_d4_name=envPath+"/EOS/s95p-v1/s95p-v1_dens4.dat";
  eos_d5_name=envPath+"/EOS/s95p-v1/s95p-v1_dens5.dat";
  eos_d6_name=envPath+"/EOS/s95p-v1/s95p-v1_dens6.dat";
  eos_d7_name=envPath+"/EOS/s95p-v1/s95p-v1_dens7.dat";

  eos_T1_name=envPath+"/EOS/s95p-v1/s95p-v1_par1.dat";
  eos_T2_name=envPath+"/EOS/s95p-v1/s95p-v1_par2.dat";
  eos_T3_name=envPath+"/EOS/s95p-v1/s95p-v1_par3.dat";
  eos_T4_name=envPath+"/EOS/s95p-v1/s95p-v1_par4.dat";
  eos_T5_name=envPath+"/EOS/s95p-v1/s95p-v1_par5.dat";
  eos_T6_name=envPath+"/EOS/s95p-v1/s95p-v1_par6.dat";
  eos_T7_name=envPath+"/EOS/s95p-v1/s95p-v1_par7.dat";
      
  eos_d1 = fopen(eos_d1_name.c_str(), "r");
  eos_d2 = fopen(eos_d2_name.c_str(), "r");
  eos_d3 = fopen(eos_d3_name.c_str(), "r");
  eos_d4 = fopen(eos_d4_name.c_str(), "r");
  eos_d5 = fopen(eos_d5_name.c_str(), "r");
  eos_d6 = fopen(eos_d6_name.c_str(), "r");
  eos_d7 = fopen(eos_d7_name.c_str(), "r");
  eos_T1 = fopen(eos_T1_name.c_str(), "r");
  eos_T2 = fopen(eos_T2_name.c_str(), "r");
  eos_T3 = fopen(eos_T3_name.c_str(), "r");
  eos_T4 = fopen(eos_T4_name.c_str(), "r");
  eos_T5 = fopen(eos_T5_name.c_str(), "r");
  eos_T6 = fopen(eos_T6_name.c_str(), "r");
  eos_T7 = fopen(eos_T7_name.c_str(), "r");
  
  checkForReadError(eos_d1,eos_d1_name.c_str());
  checkForReadError(eos_d2,eos_d2_name.c_str());
  checkForReadError(eos_d3,eos_d3_name.c_str());
  checkForReadError(eos_d4,eos_d4_name.c_str());
  checkForReadError(eos_d5,eos_d5_name.c_str());
  checkForReadError(eos_d6,eos_d6_name.c_str());
  checkForReadError(eos_d7,eos_d7_name.c_str());
  checkForReadError(eos_T1,eos_T1_name.c_str());
  checkForReadError(eos_T2,eos_T2_name.c_str());
  checkForReadError(eos_T3,eos_T3_name.c_str());
  checkForReadError(eos_T4,eos_T4_name.c_str());
  checkForReadError(eos_T5,eos_T5_name.c_str());
  checkForReadError(eos_T6,eos_T6_name.c_str());
  checkForReadError(eos_T7,eos_T7_name.c_str());
 
  int fscan_status = 0;
  //read the first two lines with general info:
  // lowest value of epsilon
  // deltaEpsilon, number of epsilon steps (i.e. # of lines)
  fscan_status = fscanf(eos_T1,"%lf",&EPP1);
  fscan_status = fscanf(eos_T1,"%lf %d",&deltaEPP1,&NEPP1);
  fscan_status = fscanf(eos_T2,"%lf",&EPP2);
  fscan_status = fscanf(eos_T2,"%lf %d",&deltaEPP2,&NEPP2);
  fscan_status = fscanf(eos_T3,"%lf",&EPP3);
  fscan_status = fscanf(eos_T3,"%lf %d",&deltaEPP3,&NEPP3);
  fscan_status = fscanf(eos_T4,"%lf",&EPP4);
  fscan_status = fscanf(eos_T4,"%lf %d",&deltaEPP4,&NEPP4);
  fscan_status = fscanf(eos_T5,"%lf",&EPP5);
  fscan_status = fscanf(eos_T5,"%lf %d",&deltaEPP5,&NEPP5);
  fscan_status = fscanf(eos_T6,"%lf",&EPP6);
  fscan_status = fscanf(eos_T6,"%lf %d",&deltaEPP6,&NEPP6);
  fscan_status = fscanf(eos_T7,"%lf",&EPP7);
  fscan_status = fscanf(eos_T7,"%lf %d",&deltaEPP7,&NEPP7);
  fscan_status = fscanf(eos_d1,"%lf",&EPP1);
  fscan_status = fscanf(eos_d1,"%lf %d",&deltaEPP1,&NEPP1);
  fscan_status = fscanf(eos_d2,"%lf",&EPP2);
  fscan_status = fscanf(eos_d2,"%lf %d",&deltaEPP2,&NEPP2);
  fscan_status = fscanf(eos_d3,"%lf",&EPP3);
  fscan_status = fscanf(eos_d3,"%lf %d",&deltaEPP3,&NEPP3);
  fscan_status = fscanf(eos_d4,"%lf",&EPP4);
  fscan_status = fscanf(eos_d4,"%lf %d",&deltaEPP4,&NEPP4);
  fscan_status = fscanf(eos_d5,"%lf",&EPP5);
  fscan_status = fscanf(eos_d5,"%lf %d",&deltaEPP5,&NEPP5);
  fscan_status = fscanf(eos_d6,"%lf",&EPP6);
  fscan_status = fscanf(eos_d6,"%lf %d",&deltaEPP6,&NEPP6);
  fscan_status = fscanf(eos_d7,"%lf",&EPP7);
  fscan_status = fscanf(eos_d7,"%lf %d",&deltaEPP7,&NEPP7);
  if (fscan_status == 0 && parameters_ptr->echo_level > 9) {
      music_message << "fscanf_status = " << fscan_status;
      music_message.flush("warning");
  }
  
  //NEPP1 -= 1;
  //NEPP2 -= 1;
  //NEPP3 -= 1;
  //NEPP4 -= 1;
  //NEPP5 -= 1;
  //NEPP6 -= 1;
  //NEPP7 -= 1;
 
  // no rho_b dependence at the moment
  NBNP1=0; 
  NBNP2=0;
  NBNP3=0;
  NBNP4=0;
  NBNP5=0;
  NBNP6=0;
  NBNP7=0;

  // allocate memory for pressure arrays
  pressure1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  pressure2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  pressure3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  pressure4=util->mtx_malloc(NBNP4+1,NEPP4+1);
  pressure5=util->mtx_malloc(NBNP5+1,NEPP5+1);
  pressure6=util->mtx_malloc(NBNP6+1,NEPP6+1);
  pressure7=util->mtx_malloc(NBNP7+1,NEPP7+1);

  // allocate memory for entropy density arrays
  entropyDensity1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  entropyDensity2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  entropyDensity3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  entropyDensity4=util->mtx_malloc(NBNP4+1,NEPP4+1);
  entropyDensity5=util->mtx_malloc(NBNP5+1,NEPP5+1);
  entropyDensity6=util->mtx_malloc(NBNP6+1,NEPP6+1);
  entropyDensity7=util->mtx_malloc(NBNP7+1,NEPP7+1);

  // allocate memory for QGP fraction arrays
  QGPfraction1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  QGPfraction2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  QGPfraction3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  QGPfraction4=util->mtx_malloc(NBNP4+1,NEPP4+1);
  QGPfraction5=util->mtx_malloc(NBNP5+1,NEPP5+1);
  QGPfraction6=util->mtx_malloc(NBNP6+1,NEPP6+1);
  QGPfraction7=util->mtx_malloc(NBNP7+1,NEPP7+1);

  // allocate memory for temperature arrays
  temperature1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  temperature2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  temperature3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  temperature4=util->mtx_malloc(NBNP4+1,NEPP4+1);
  temperature5=util->mtx_malloc(NBNP5+1,NEPP5+1);
  temperature6=util->mtx_malloc(NBNP6+1,NEPP6+1);
  temperature7=util->mtx_malloc(NBNP7+1,NEPP7+1);
  
  // allocate memory for velocity of sound squared arrays
  cs2_1 = util->mtx_malloc(NBNP1+1, NEPP1+1);
  cs2_2 = util->mtx_malloc(NBNP2+1, NEPP2+1);
  cs2_3 = util->mtx_malloc(NBNP3+1, NEPP3+1);
  cs2_4 = util->mtx_malloc(NBNP4+1, NEPP4+1);
  cs2_5 = util->mtx_malloc(NBNP5+1, NEPP5+1);
  cs2_6 = util->mtx_malloc(NBNP6+1, NEPP6+1);
  cs2_7 = util->mtx_malloc(NBNP7+1, NEPP7+1);

  // allocate memory for baryon chemical potential arrays
  // currently always zero
  /*   mu1=mtx_malloc(NBNP1+1,NEPP1+1); */
  /*   mu2=mtx_malloc(NBNP2+1,NEPP2+1); */
  /*   mu3=mtx_malloc(NBNP3+1,NEPP3+1); */
  /*   mu4=mtx_malloc(NBNP4+1,NEPP4+1); */
  
  // read pressure, temperature and chemical potential values
  // files have it backwards, so I start with maximum j and count down
  i=0;
  for(j=NEPP1-1; j>=0; j--)
    {
      fscan_status = fscanf(eos_d1,"%lf",&eps);
      fscan_status = fscanf(eos_d1,"%lf",&pressure1[i][j]);
      fscan_status = fscanf(eos_d1,"%lf",&entropyDensity1[i][j]);
      fscan_status = fscanf(eos_d1,"%lf",&baryonDensity);
      fscan_status = fscanf(eos_d1,"%lf",&QGPfraction1[i][j]);
      fscan_status = fscanf(eos_T1,"%lf",&temperature1[i][j]);
      fscan_status = fscanf(eos_T1,"%lf",&eps); //dummy
      fscan_status = fscanf(eos_T1,"%lf",&eps); //dummy
    }

  for(j=NEPP2-1; j>=0; j--)
    {
      fscan_status = fscanf(eos_d2,"%lf",&eps);
      fscan_status = fscanf(eos_d2,"%lf",&pressure2[i][j]);
      fscan_status = fscanf(eos_d2,"%lf",&entropyDensity2[i][j]);
      fscan_status = fscanf(eos_d2,"%lf",&baryonDensity);
      fscan_status = fscanf(eos_d2,"%lf",&QGPfraction2[i][j]);
      fscan_status = fscanf(eos_T2,"%lf",&temperature2[i][j]);
      fscan_status = fscanf(eos_T2,"%lf",&eps); //dummy
      fscan_status = fscanf(eos_T2,"%lf",&eps); //dummy
    }

  for(j=NEPP3-1; j>=0; j--)
    {
      fscan_status = fscanf(eos_d3,"%lf",&eps);
      fscan_status = fscanf(eos_d3,"%lf",&pressure3[i][j]);
      fscan_status = fscanf(eos_d3,"%lf",&entropyDensity3[i][j]);
      fscan_status = fscanf(eos_d3,"%lf",&baryonDensity);
      fscan_status = fscanf(eos_d3,"%lf",&QGPfraction3[i][j]);
      fscan_status = fscanf(eos_T3,"%lf",&temperature3[i][j]);
      fscan_status = fscanf(eos_T3,"%lf",&eps); //dummy
      fscan_status = fscanf(eos_T3,"%lf",&eps); //dummy
    }

  for(j=NEPP4-1; j>=0; j--)
    {
      fscan_status = fscanf(eos_d4,"%lf",&eps);
      fscan_status = fscanf(eos_d4,"%lf",&pressure4[i][j]);
      fscan_status = fscanf(eos_d4,"%lf",&entropyDensity4[i][j]);
      fscan_status = fscanf(eos_d4,"%lf",&baryonDensity);
      fscan_status = fscanf(eos_d4,"%lf",&QGPfraction4[i][j]);
      fscan_status = fscanf(eos_T4,"%lf",&temperature4[i][j]);
      fscan_status = fscanf(eos_T4,"%lf",&eps); //dummy
      fscan_status = fscanf(eos_T4,"%lf",&eps); //dummy
    }

  for(j=NEPP5-1; j>=0; j--)
    {
      fscan_status = fscanf(eos_d5,"%lf",&eps);
      fscan_status = fscanf(eos_d5,"%lf",&pressure5[i][j]);
      fscan_status = fscanf(eos_d5,"%lf",&entropyDensity5[i][j]);
      fscan_status = fscanf(eos_d5,"%lf",&baryonDensity);
      fscan_status = fscanf(eos_d5,"%lf",&QGPfraction5[i][j]);
      fscan_status = fscanf(eos_T5,"%lf",&temperature5[i][j]);
      fscan_status = fscanf(eos_T5,"%lf",&eps); //dummy
      fscan_status = fscanf(eos_T5,"%lf",&eps); //dummy
    } 

  for(j=NEPP6-1; j>=0; j--)
    {
      fscan_status = fscanf(eos_d6,"%lf",&eps);
      fscan_status = fscanf(eos_d6,"%lf",&pressure6[i][j]);
      fscan_status = fscanf(eos_d6,"%lf",&entropyDensity6[i][j]);
      fscan_status = fscanf(eos_d6,"%lf",&baryonDensity);
      fscan_status = fscanf(eos_d6,"%lf",&QGPfraction6[i][j]);
      fscan_status = fscanf(eos_T6,"%lf",&temperature6[i][j]);
      fscan_status = fscanf(eos_T6,"%lf",&eps); //dummy
      fscan_status = fscanf(eos_T6,"%lf",&eps); //dummy
    } 

  for(j=NEPP7-1; j>=0; j--)
    {
      fscan_status = fscanf(eos_d7,"%lf",&eps);
      fscan_status = fscanf(eos_d7,"%lf",&pressure7[i][j]);
      fscan_status = fscanf(eos_d7,"%lf",&entropyDensity7[i][j]);
      fscan_status = fscanf(eos_d7,"%lf",&baryonDensity);
      fscan_status = fscanf(eos_d7,"%lf",&QGPfraction7[i][j]);
      fscan_status = fscanf(eos_T7,"%lf",&temperature7[i][j]);
      fscan_status = fscanf(eos_T7,"%lf",&eps); //dummy
      fscan_status = fscanf(eos_T7,"%lf",&eps); //dummy
    } 

//test if the reading worked:
/*   for(j=0; j<NEPP1; j++) */
/*     { */
/*       fprintf(stderr,"%lf %lf %lf %lf\n",pressure1[i][j],entropyDensity1[i][j],QGPfraction1[i][j],temperature1[i][j]); */
/*     } */
/*   for(j=0; j<NEPP2; j++) */
/*     { */
/*       fprintf(stderr,"%lf %lf %lf %lf\n",pressure2[i][j],entropyDensity2[i][j],QGPfraction2[i][j],temperature2[i][j]); */
/*     } */
/*   for(j=0; j<NEPP3; j++) */
/*     { */
/*       fprintf(stderr,"%lf %lf %lf %lf\n",pressure3[i][j],entropyDensity3[i][j],QGPfraction3[i][j],temperature3[i][j]); */
/*     } */
/*   for(j=0; j<NEPP4; j++) */
/*     { */
/*       fprintf(stderr,"%lf %lf %lf %lf\n",pressure4[i][j],entropyDensity4[i][j],QGPfraction4[i][j],temperature4[i][j]); */
/*     } */

  fclose(eos_d1);
  fclose(eos_d2);
  fclose(eos_d3);
  fclose(eos_d4);
  fclose(eos_d5);
  fclose(eos_d6);
  fclose(eos_d7);
  fclose(eos_T1);
  fclose(eos_T2);
  fclose(eos_T3);
  fclose(eos_T4);
  fclose(eos_T5);
  fclose(eos_T6);
  fclose(eos_T7);

  music_message.info("Done reading EOS.");

  build_velocity_of_sound_sq_matrix();

  //Allocate and fill arrays for get_s2e()
  //First, create new arrays that will permit to extract information we need from each EOS table in a unified way
  const int number_of_EOS_tables = 7;
  double lowest_eps_list[number_of_EOS_tables] = {EPP1,EPP2,EPP3,EPP4,EPP5,EPP6,EPP7};
  double delta_eps_list[number_of_EOS_tables] = {deltaEPP1,deltaEPP2,deltaEPP3,deltaEPP4,deltaEPP5,deltaEPP6,deltaEPP7};
  int nb_elements_list[number_of_EOS_tables] = {NEPP1,NEPP2,NEPP3,NEPP4,NEPP5,NEPP6,NEPP7};
  double ** s_list[number_of_EOS_tables] = {entropyDensity1,entropyDensity2,entropyDensity3,entropyDensity4,entropyDensity5,entropyDensity6,entropyDensity7};
  //Compute the maximum number of entropyDensity elements (the actual number of entropyDensity elements we will use 
  //will be smaller since we skip duplicate entropyDensities)
  int s_list_rho0_maximal_length=0;
  for(int kk=0;kk<number_of_EOS_tables;kk++) s_list_rho0_maximal_length+=nb_elements_list[kk];
  double last=0.0, tmp_s;
  s_list_rho0_length=0;
  //Allocate memory
  eps_list_rho0= new double [s_list_rho0_maximal_length];
  s_list_rho0= new double [s_list_rho0_maximal_length];

  for(int table_id=0; table_id<number_of_EOS_tables;table_id++) {
	  for(i=0;i<nb_elements_list[table_id];i++) {
		tmp_s=s_list[table_id][0][i];
		if (tmp_s > last) {
			eps_list_rho0[s_list_rho0_length]=(lowest_eps_list[table_id]+i*delta_eps_list[table_id])/hbarc;
			s_list_rho0[s_list_rho0_length]=tmp_s;
			s_list_rho0_length++;
		}
		last=tmp_s;
	  }
  }
  //Initialise the gsl interpolator
  interp_s2e = gsl_interp_alloc(gsl_interp_cspline, s_list_rho0_length);
  gsl_interp_init(interp_s2e, s_list_rho0, eps_list_rho0, s_list_rho0_length);
  accel_s2e = gsl_interp_accel_alloc();
}

void EOS::init_eos3(int selector)
{
  // read the lattice EOS pressure, temperature, and 
  // baryon chemical potential from file
  music_message.info("reading EOS...");
  whichEOS = 3; 
  int i, j;
  std::string envPath;
  if (getenv("HYDROPROGRAMPATH") != 0) {
    envPath=getenv("HYDROPROGRAMPATH");
  }
  else {
	envPath="";
  }
  
  stringstream spath;
  stringstream slocalpath;
  stringstream streos_d1_name;
  stringstream streos_d2_name;
  stringstream streos_d3_name;
  stringstream streos_d4_name;
  stringstream streos_d5_name;
  stringstream streos_d6_name;
  stringstream streos_d7_name;
  stringstream streos_T1_name;
  stringstream streos_T2_name;
  stringstream streos_T3_name;
  stringstream streos_T4_name;
  stringstream streos_T5_name;
  stringstream streos_T6_name;
  stringstream streos_T7_name;
  string eos_d1_name;
  string eos_d2_name;
  string eos_d3_name;
  string eos_d4_name;
  string eos_d5_name;
  string eos_d6_name;
  string eos_d7_name;
  string eos_T1_name;
  string eos_T2_name;
  string eos_T3_name;
  string eos_T4_name;
  string eos_T5_name;
  string eos_T6_name;
  string eos_T7_name;
  double eps, baryonDensity; //dummies for now

  spath << envPath;
  if(selector == 1)
    {
      spath << "/EOS/s95p-PCE-v1/";
      slocalpath << "./EOS/s95p-PCE-v1/";
    }
  else if(selector ==2)
    {
      spath << "/EOS/s95p-PCE155/";
      slocalpath << "./EOS/s95p-PCE155/";
     }  
  else if(selector ==3)
    {
      spath << "/EOS/s95p-PCE160/";
      slocalpath << "./EOS/s95p-PCE160/";
    }
  else if(selector ==4)
    {
      spath << "/EOS/s95p-PCE165-v0/";
      slocalpath << "./EOS/s95p-PCE165-v0/";
    }

  string path = spath.str();
  string localpath = slocalpath.str();


  if (envPath != "") // if path is set in the environment
    {
      streos_d1_name << path << "dens1.dat";
      eos_d1_name = streos_d1_name.str();
      streos_d2_name << path << "dens2.dat";
      eos_d2_name = streos_d2_name.str();
      streos_d3_name << path << "dens3.dat";
      eos_d3_name = streos_d3_name.str();
      streos_d4_name << path << "dens4.dat";
      eos_d4_name = streos_d4_name.str();
      streos_d5_name << path << "dens5.dat";
      eos_d5_name = streos_d5_name.str();
      streos_d6_name << path << "dens6.dat";
      eos_d6_name = streos_d6_name.str();
      streos_d7_name << path << "dens7.dat";
      eos_d7_name = streos_d7_name.str();

      streos_T1_name << path << "par1.dat";
      eos_T1_name = streos_T1_name.str();
      streos_T2_name << path << "par2.dat";
      eos_T2_name = streos_T2_name.str();
      streos_T3_name << path << "par3.dat";
      eos_T3_name = streos_T3_name.str();
      streos_T4_name << path << "par4.dat";
      eos_T4_name = streos_T4_name.str();
      streos_T5_name << path << "par5.dat";
      eos_T5_name = streos_T5_name.str();
      streos_T6_name << path << "par6.dat";
      eos_T6_name = streos_T6_name.str();
      streos_T7_name << path << "par7.dat";
      eos_T7_name = streos_T7_name.str();
      
      music_message << "from path " << path;
      music_message.flush("info");
    }
  else //if path is not set in the environment use current folder and then /EOS/s95p-PCE160 subfolder
    {
      streos_d1_name << localpath << "dens1.dat";
      eos_d1_name = streos_d1_name.str();
      streos_d2_name << localpath << "dens2.dat";
      eos_d2_name = streos_d2_name.str();
      streos_d3_name << localpath << "dens3.dat";
      eos_d3_name = streos_d3_name.str();
      streos_d4_name << localpath << "dens4.dat";
      eos_d4_name = streos_d4_name.str();
      streos_d5_name << localpath << "dens5.dat";
      eos_d5_name = streos_d5_name.str();
      streos_d6_name << localpath << "dens6.dat";
      eos_d6_name = streos_d6_name.str();
      streos_d7_name << localpath << "dens7.dat";
      eos_d7_name = streos_d7_name.str();

      streos_T1_name << localpath << "par1.dat";
      eos_T1_name = streos_T1_name.str();
      streos_T2_name << localpath << "par2.dat";
      eos_T2_name = streos_T2_name.str();
      streos_T3_name << localpath << "par3.dat";
      eos_T3_name = streos_T3_name.str();
      streos_T4_name << localpath << "par4.dat";
      eos_T4_name = streos_T4_name.str();
      streos_T5_name << localpath << "par5.dat";
      eos_T5_name = streos_T5_name.str();
      streos_T6_name << localpath << "par6.dat";
      eos_T6_name = streos_T6_name.str();
      streos_T7_name << localpath << "par7.dat";
      eos_T7_name = streos_T7_name.str();
       
      music_message << "from path " << localpath;
      music_message.flush("info");
    
    }
 

  ifstream eos_d1(eos_d1_name.c_str());
  ifstream eos_d2(eos_d2_name.c_str());
  ifstream eos_d3(eos_d3_name.c_str());
  ifstream eos_d4(eos_d4_name.c_str());
  ifstream eos_d5(eos_d5_name.c_str());
  ifstream eos_d6(eos_d6_name.c_str());
  ifstream eos_d7(eos_d7_name.c_str());
  ifstream eos_T1(eos_T1_name.c_str());
  ifstream eos_T2(eos_T2_name.c_str());
  ifstream eos_T3(eos_T3_name.c_str());
  ifstream eos_T4(eos_T4_name.c_str());
  ifstream eos_T5(eos_T5_name.c_str());
  ifstream eos_T6(eos_T6_name.c_str());
  ifstream eos_T7(eos_T7_name.c_str());

  //checkForReadError(eos_d1,eos_d1_name);
  //checkForReadError(eos_d2,eos_d2_name);
  //checkForReadError(eos_d3,eos_d3_name);
  //checkForReadError(eos_d4,eos_d4_name);
  //checkForReadError(eos_T1,eos_T1_name);
  //checkForReadError(eos_T2,eos_T2_name);
  //checkForReadError(eos_T3,eos_T3_name);
  //checkForReadError(eos_T4,eos_T4_name);
  
  //read the first two lines with general info:
  // lowest value of epsilon
  // deltaEpsilon, number of epsilon steps (i.e. # of lines)
  eos_T1 >> EPP1;
  eos_T1 >> deltaEPP1;
  eos_T1 >> NEPP1;
  eos_T2 >> EPP2;
  eos_T2 >> deltaEPP2;
  eos_T2 >> NEPP2;
  eos_T3 >> EPP3;
  eos_T3 >> deltaEPP3;
  eos_T3 >> NEPP3;
  eos_T4 >> EPP4;
  eos_T4 >> deltaEPP4;
  eos_T4 >> NEPP4;
  eos_T5 >> EPP5;
  eos_T5 >> deltaEPP5;
  eos_T5 >> NEPP5;
  eos_T6 >> EPP6;
  eos_T6 >> deltaEPP6;
  eos_T6 >> NEPP6;
  eos_T7 >> EPP7;
  eos_T7 >> deltaEPP7;
  eos_T7 >> NEPP7;
  eos_d1 >> EPP1;
  eos_d1 >> deltaEPP1;
  eos_d1 >> NEPP1;
  eos_d2 >> EPP2;
  eos_d2 >> deltaEPP2;
  eos_d2 >> NEPP2;
  eos_d3 >> EPP3;
  eos_d3 >> deltaEPP3;
  eos_d3 >> NEPP3;
  eos_d4 >> EPP4;
  eos_d4 >> deltaEPP4;
  eos_d4 >> NEPP4;
  eos_d5 >> EPP5;
  eos_d5 >> deltaEPP5;
  eos_d5 >> NEPP5;
  eos_d6 >> EPP6;
  eos_d6 >> deltaEPP6;
  eos_d6 >> NEPP6;
  eos_d7 >> EPP7;
  eos_d7 >> deltaEPP7;
  eos_d7 >> NEPP7;
  
  // fscanf(eos_d1,"%lf",&EPP1);
  // fscanf(eos_d1,"%lf %d",&deltaEPP1,&NEPP1);
  // cout << "check 4" << endl;

  // no rho_b dependence at the moment
  NBNP1=0; 
  NBNP2=0;
  NBNP3=0;
  NBNP4=0;
  NBNP5=0;
  NBNP6=0;
  NBNP7=0;

  // allocate memory for pressure arrays
  pressure1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  pressure2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  pressure3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  pressure4=util->mtx_malloc(NBNP4+1,NEPP4+1);
  pressure5=util->mtx_malloc(NBNP5+1,NEPP5+1);
  pressure6=util->mtx_malloc(NBNP6+1,NEPP6+1);
  pressure7=util->mtx_malloc(NBNP7+1,NEPP7+1);
 
  // allocate memory for entropy density arrays
  entropyDensity1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  entropyDensity2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  entropyDensity3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  entropyDensity4=util->mtx_malloc(NBNP4+1,NEPP4+1);
  entropyDensity5=util->mtx_malloc(NBNP5+1,NEPP5+1);
  entropyDensity6=util->mtx_malloc(NBNP6+1,NEPP6+1);
  entropyDensity7=util->mtx_malloc(NBNP7+1,NEPP7+1);

  // allocate memory for QGP fraction arrays
  QGPfraction1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  QGPfraction2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  QGPfraction3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  QGPfraction4=util->mtx_malloc(NBNP4+1,NEPP4+1);
  QGPfraction5=util->mtx_malloc(NBNP5+1,NEPP5+1);
  QGPfraction6=util->mtx_malloc(NBNP6+1,NEPP6+1);
  QGPfraction7=util->mtx_malloc(NBNP7+1,NEPP7+1);

  // allocate memory for temperature arrays
  temperature1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  temperature2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  temperature3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  temperature4=util->mtx_malloc(NBNP4+1,NEPP4+1);
  temperature5=util->mtx_malloc(NBNP5+1,NEPP5+1);
  temperature6=util->mtx_malloc(NBNP6+1,NEPP6+1);
  temperature7=util->mtx_malloc(NBNP7+1,NEPP7+1);
  
  // allocate memory for velocity of sound squared arrays
  cs2_1 = util->mtx_malloc(NBNP1+1, NEPP1+1);
  cs2_2 = util->mtx_malloc(NBNP2+1, NEPP2+1);
  cs2_3 = util->mtx_malloc(NBNP3+1, NEPP3+1);
  cs2_4 = util->mtx_malloc(NBNP4+1, NEPP4+1);
  cs2_5 = util->mtx_malloc(NBNP5+1, NEPP5+1);
  cs2_6 = util->mtx_malloc(NBNP6+1, NEPP6+1);
  cs2_7 = util->mtx_malloc(NBNP7+1, NEPP7+1);

  // allocate memory for baryon chemical potential arrays
  // currently always zero
  /*   mu1=mtx_malloc(NBNP1+1,NEPP1+1); */
  /*   mu2=mtx_malloc(NBNP2+1,NEPP2+1); */
  /*   mu3=mtx_malloc(NBNP3+1,NEPP3+1); */
  /*   mu4=mtx_malloc(NBNP4+1,NEPP4+1); */
  
  // read pressure, temperature and chemical potential values
  // files have it backwards, so I start with maximum j and count down


  i=0;
  for(j=NEPP1-1; j>=0; j--)
    {
      eos_d1 >> eps;
      eos_d1 >> pressure1[i][j];
      eos_d1 >> entropyDensity1[i][j];
      eos_d1 >> baryonDensity;
      eos_d1 >> QGPfraction1[i][j];
      //      fscanf(eos_d1,"%lf",&eps);
      //fscanf(eos_d1,"%lf",&pressure1[i][j]);
      //fscanf(eos_d1,"%lf",&entropyDensity1[i][j]);
      //fscanf(eos_d1,"%lf",&baryonDensity);
      //fscanf(eos_d1,"%lf",&QGPfraction1[i][j]);
      eos_T1 >> temperature1[i][j];
      eos_T1 >> eps;
      eos_T1 >> eps;
      //fscanf(eos_T1,"%lf",&temperature1[i][j]);
      //fscanf(eos_T1,"%lf",&eps); //dummy
      //fscanf(eos_T1,"%lf",&eps); //dummy
    }

  for(j=NEPP2-1; j>=0; j--)
    {
      eos_d2 >> eps;
      eos_d2 >> pressure2[i][j];
      eos_d2 >> entropyDensity2[i][j];
      eos_d2 >> baryonDensity;
      eos_d2 >> QGPfraction2[i][j];
      eos_T2 >> temperature2[i][j];
      eos_T2 >> eps;
      eos_T2 >> eps;
    }
  for(j=NEPP3-1; j>=0; j--)
    {
      eos_d3 >> eps;
      eos_d3 >> pressure3[i][j];
      eos_d3 >> entropyDensity3[i][j];
      eos_d3 >> baryonDensity;
      eos_d3 >> QGPfraction3[i][j];
      eos_T3 >> temperature3[i][j];
      eos_T3 >> eps;
      eos_T3 >> eps;
    }
  for(j=NEPP4-1; j>=0; j--)
    {
      eos_d4 >> eps;
      eos_d4 >> pressure4[i][j];
      eos_d4 >> entropyDensity4[i][j];
      eos_d4 >> baryonDensity;
      eos_d4 >> QGPfraction4[i][j];
      eos_T4 >> temperature4[i][j];
      eos_T4 >> eps;
      eos_T4 >> eps;
    }
  for(j=NEPP5-1; j>=0; j--)
    {
      eos_d5 >> eps;
      eos_d5 >> pressure5[i][j];
      eos_d5 >> entropyDensity5[i][j];
      eos_d5 >> baryonDensity;
      eos_d5 >> QGPfraction5[i][j];
      eos_T5 >> temperature5[i][j];
      eos_T5 >> eps;
      eos_T5 >> eps;
    }
  for(j=NEPP6-1; j>=0; j--)
    {
      eos_d6 >> eps;
      eos_d6 >> pressure6[i][j];
      eos_d6 >> entropyDensity6[i][j];
      eos_d6 >> baryonDensity;
      eos_d6 >> QGPfraction6[i][j];
      eos_T6 >> temperature6[i][j];
      eos_T6 >> eps;
      eos_T6 >> eps;
    }
  for(j=NEPP7-1; j>=0; j--)
    {
      eos_d7 >> eps;
      eos_d7 >> pressure7[i][j];
      eos_d7 >> entropyDensity7[i][j];
      eos_d7 >> baryonDensity;
      eos_d7 >> QGPfraction7[i][j];
      eos_T7 >> temperature7[i][j];
      eos_T7 >> eps;
      eos_T7 >> eps;
    }

 //  //test if the reading worked:
//   for(j=0; j<NEPP1; j++) 
//     { 
//       fprintf(stderr,"%lf %lf %lf %lf\n",pressure1[i][j],entropyDensity1[i][j],QGPfraction1[i][j],temperature1[i][j]); 
//     } 
//   for(j=0; j<NEPP2; j++) 
//     { 
//       fprintf(stderr,"%lf %lf %lf %lf\n",pressure2[i][j],entropyDensity2[i][j],QGPfraction2[i][j],temperature2[i][j]); 
//     } 
//   for(j=0; j<NEPP3; j++) 
//     { 
//       fprintf(stderr,"%lf %lf %lf %lf\n",pressure3[i][j],entropyDensity3[i][j],QGPfraction3[i][j],temperature3[i][j]); 
//     } 
//   cout << " check before 4" << endl;
//   for(j=0; j<NEPP4; j++) 
//     { 
//       fprintf(stderr,"%lf %lf %lf %lf\n",pressure4[i][j],entropyDensity4[i][j],QGPfraction4[i][j],temperature4[i][j]); 
//     } 
//   cout << " check before 5" << endl;
//   for(j=0; j<NEPP5; j++) 
//     { 
//       fprintf(stderr,"%lf %lf %lf %lf\n",pressure5[i][j],entropyDensity5[i][j],QGPfraction5[i][j],temperature5[i][j]); 
//     } 
//   cout << " check before 6" << endl;
//   for(j=0; j<NEPP6; j++) 
//     { 
//       fprintf(stderr,"%lf %lf %lf %lf\n",pressure6[i][j],entropyDensity6[i][j],QGPfraction6[i][j],temperature6[i][j]); 
//     } 
//   cout << " check before 7" << endl;
//   for(j=0; j<NEPP7; j++) 
//     { 
//       fprintf(stderr,"%lf %lf %lf %lf\n",pressure7[i][j],entropyDensity7[i][j],QGPfraction7[i][j],temperature7[i][j]); 
//     } 

  
  eos_d1.close();
  eos_d2.close();
  eos_d3.close();
  eos_d4.close();
  eos_d5.close();
  eos_d6.close();
  eos_d7.close();
  eos_T1.close();
  eos_T2.close();
  eos_T3.close();
  eos_T4.close();
  eos_T5.close();
  eos_T6.close();
  eos_T7.close();

  music_message.info("Done reading EOS.");

  build_velocity_of_sound_sq_matrix();

  //Allocate and fill arrays for get_s2e()
  //First, create new arrays that will permit to extract information we need from each EOS table in a unified way
  const int number_of_EOS_tables = 7;
  double lowest_eps_list[number_of_EOS_tables] = {EPP1,EPP2,EPP3,EPP4,EPP5,EPP6,EPP7};
  double delta_eps_list[number_of_EOS_tables] = {deltaEPP1,deltaEPP2,deltaEPP3,deltaEPP4,deltaEPP5,deltaEPP6,deltaEPP7};
  int nb_elements_list[number_of_EOS_tables] = {NEPP1,NEPP2,NEPP3,NEPP4,NEPP5,NEPP6,NEPP7};
  double ** s_list[number_of_EOS_tables] = {entropyDensity1,entropyDensity2,entropyDensity3,entropyDensity4,entropyDensity5,entropyDensity6,entropyDensity7};
  //Compute the maximum number of entropyDensity elements (the actual number of entropyDensity elements we will use 
  //will be smaller since we skip duplicate entropyDensities)
  int s_list_rho0_maximal_length=0;
  for(int kk=0;kk<number_of_EOS_tables;kk++) s_list_rho0_maximal_length+=nb_elements_list[kk];
  double last=0.0, tmp_s;
  s_list_rho0_length=0;
  //Allocate memory
  eps_list_rho0= new double [s_list_rho0_maximal_length];
  s_list_rho0= new double [s_list_rho0_maximal_length];

  for(int table_id=0; table_id<number_of_EOS_tables;table_id++) {
	  for(i=0;i<nb_elements_list[table_id];i++) {
		tmp_s=s_list[table_id][0][i];
		if (tmp_s > last) {
			eps_list_rho0[s_list_rho0_length]=(lowest_eps_list[table_id]+i*delta_eps_list[table_id])/hbarc;
			s_list_rho0[s_list_rho0_length]=tmp_s;
			s_list_rho0_length++;
		}
		last=tmp_s;
	  }
  }
  //Initialise the gsl interpolator
  interp_s2e = gsl_interp_alloc(gsl_interp_cspline, s_list_rho0_length);
  gsl_interp_init(interp_s2e, s_list_rho0, eps_list_rho0, s_list_rho0_length);
  accel_s2e = gsl_interp_accel_alloc();
}

void EOS::init_eos7()
{
  // read the lattice EOS pressure, temperature, and 
  // baryon chemical potential from file
  music_message.info("reading EOS s95p-v1.2 (for UrQMD) ...");
  whichEOS = 7; 
//   int bytes_read;
  int i, j;
  FILE *eos_d1, *eos_d2, *eos_d3, *eos_d4, *eos_d5, *eos_d6, *eos_d7;
  FILE *eos_T1, *eos_T2, *eos_T3, *eos_T4, *eos_T5, *eos_T6, *eos_T7;
  const char* EOSPATH = "HYDROPROGRAMPATH";
  char * pre_envPath= getenv(EOSPATH);
  std::string envPath;
  std::string eos_d1_name;
  std::string eos_d2_name;
  std::string eos_d3_name;
  std::string eos_d4_name;
  std::string eos_d5_name;
  std::string eos_d6_name;
  std::string eos_d7_name;
  std::string eos_T1_name;
  std::string eos_T2_name;
  std::string eos_T3_name;
  std::string eos_T4_name;
  std::string eos_T5_name;
  std::string eos_T6_name;
  std::string eos_T7_name;
  double eps, baryonDensity; //dummies for now

  if (pre_envPath == 0) {
	  envPath=".";
  }
  else {
	  envPath=pre_envPath;
  }

  music_message << "from path " << envPath.c_str() << "/EOS/s95p-v1.2";
  music_message.flush("info");

  eos_d1_name=envPath+"/EOS/s95p-v1.2/s95p-v1.2_dens1.dat";
  eos_d2_name=envPath+"/EOS/s95p-v1.2/s95p-v1.2_dens2.dat";
  eos_d3_name=envPath+"/EOS/s95p-v1.2/s95p-v1.2_dens3.dat";
  eos_d4_name=envPath+"/EOS/s95p-v1.2/s95p-v1.2_dens4.dat";
  eos_d5_name=envPath+"/EOS/s95p-v1.2/s95p-v1.2_dens5.dat";
  eos_d6_name=envPath+"/EOS/s95p-v1.2/s95p-v1.2_dens6.dat";
  eos_d7_name=envPath+"/EOS/s95p-v1.2/s95p-v1.2_dens7.dat";


  eos_T1_name=envPath+"/EOS/s95p-v1.2/s95p-v1.2_par1.dat";
  eos_T2_name=envPath+"/EOS/s95p-v1.2/s95p-v1.2_par2.dat";
  eos_T3_name=envPath+"/EOS/s95p-v1.2/s95p-v1.2_par3.dat";
  eos_T4_name=envPath+"/EOS/s95p-v1.2/s95p-v1.2_par4.dat";
  eos_T5_name=envPath+"/EOS/s95p-v1.2/s95p-v1.2_par5.dat";
  eos_T6_name=envPath+"/EOS/s95p-v1.2/s95p-v1.2_par6.dat";
  eos_T7_name=envPath+"/EOS/s95p-v1.2/s95p-v1.2_par7.dat";
  
  eos_d1 = fopen(eos_d1_name.c_str(), "r");
  eos_d2 = fopen(eos_d2_name.c_str(), "r");
  eos_d3 = fopen(eos_d3_name.c_str(), "r");
  eos_d4 = fopen(eos_d4_name.c_str(), "r");
  eos_d5 = fopen(eos_d5_name.c_str(), "r");
  eos_d6 = fopen(eos_d6_name.c_str(), "r");
  eos_d7 = fopen(eos_d7_name.c_str(), "r");
  eos_T1 = fopen(eos_T1_name.c_str(), "r");
  eos_T2 = fopen(eos_T2_name.c_str(), "r");
  eos_T3 = fopen(eos_T3_name.c_str(), "r");
  eos_T4 = fopen(eos_T4_name.c_str(), "r");
  eos_T5 = fopen(eos_T5_name.c_str(), "r");
  eos_T6 = fopen(eos_T6_name.c_str(), "r");
  eos_T7 = fopen(eos_T7_name.c_str(), "r");
  
  checkForReadError(eos_d1,eos_d1_name.c_str());
  checkForReadError(eos_d2,eos_d2_name.c_str());
  checkForReadError(eos_d3,eos_d3_name.c_str());
  checkForReadError(eos_d4,eos_d4_name.c_str());
  checkForReadError(eos_d5,eos_d5_name.c_str());
  checkForReadError(eos_d6,eos_d6_name.c_str());
  checkForReadError(eos_d7,eos_d7_name.c_str());
  checkForReadError(eos_T1,eos_T1_name.c_str());
  checkForReadError(eos_T2,eos_T2_name.c_str());
  checkForReadError(eos_T3,eos_T3_name.c_str());
  checkForReadError(eos_T4,eos_T4_name.c_str());
  checkForReadError(eos_T5,eos_T5_name.c_str());
  checkForReadError(eos_T6,eos_T6_name.c_str());
  checkForReadError(eos_T7,eos_T7_name.c_str());

  int fscan_status = 0;
  // read the first two lines with general info:
  // lowest value of epsilon
  // deltaEpsilon, number of epsilon steps (i.e. # of lines)
  fscan_status = fscanf(eos_T1,"%lf",&EPP1);
  fscan_status = fscanf(eos_T1,"%lf %d",&deltaEPP1,&NEPP1);
  fscan_status = fscanf(eos_T2,"%lf",&EPP2);
  fscan_status = fscanf(eos_T2,"%lf %d",&deltaEPP2,&NEPP2);
  fscan_status = fscanf(eos_T3,"%lf",&EPP3);
  fscan_status = fscanf(eos_T3,"%lf %d",&deltaEPP3,&NEPP3);
  fscan_status = fscanf(eos_T4,"%lf",&EPP4);
  fscan_status = fscanf(eos_T4,"%lf %d",&deltaEPP4,&NEPP4);
  fscan_status = fscanf(eos_T5,"%lf",&EPP5);
  fscan_status = fscanf(eos_T5,"%lf %d",&deltaEPP5,&NEPP5);
  fscan_status = fscanf(eos_T6,"%lf",&EPP6);
  fscan_status = fscanf(eos_T6,"%lf %d",&deltaEPP6,&NEPP6);
  fscan_status = fscanf(eos_T7,"%lf",&EPP7);
  fscan_status = fscanf(eos_T7,"%lf %d",&deltaEPP7,&NEPP7);
  fscan_status = fscanf(eos_d1,"%lf",&EPP1);
  fscan_status = fscanf(eos_d1,"%lf %d",&deltaEPP1,&NEPP1);
  fscan_status = fscanf(eos_d2,"%lf",&EPP2);
  fscan_status = fscanf(eos_d2,"%lf %d",&deltaEPP2,&NEPP2);
  fscan_status = fscanf(eos_d3,"%lf",&EPP3);
  fscan_status = fscanf(eos_d3,"%lf %d",&deltaEPP3,&NEPP3);
  fscan_status = fscanf(eos_d4,"%lf",&EPP4);
  fscan_status = fscanf(eos_d4,"%lf %d",&deltaEPP4,&NEPP4);
  fscan_status = fscanf(eos_d5,"%lf",&EPP5);
  fscan_status = fscanf(eos_d5,"%lf %d",&deltaEPP5,&NEPP5);
  fscan_status = fscanf(eos_d6,"%lf",&EPP6);
  fscan_status = fscanf(eos_d6,"%lf %d",&deltaEPP6,&NEPP6);
  fscan_status = fscanf(eos_d7,"%lf",&EPP7);
  fscan_status = fscanf(eos_d7,"%lf %d",&deltaEPP7,&NEPP7);
  if (fscan_status == 0 && parameters_ptr->echo_level > 9) {
      music_message << "fscanf_status = " << fscan_status;
      music_message.flush("warning");
  }
  
  //NEPP1 -= 1;
  //NEPP2 -= 1;
  //NEPP3 -= 1;
  //NEPP4 -= 1;
  //NEPP5 -= 1;
  //NEPP6 -= 1;
  //NEPP7 -= 1;
 
  // no rho_b dependence at the moment
  NBNP1=0; 
  NBNP2=0;
  NBNP3=0;
  NBNP4=0;
  NBNP5=0;
  NBNP6=0;
  NBNP7=0;

  // allocate memory for pressure arrays
  pressure1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  pressure2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  pressure3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  pressure4=util->mtx_malloc(NBNP4+1,NEPP4+1);
  pressure5=util->mtx_malloc(NBNP5+1,NEPP5+1);
  pressure6=util->mtx_malloc(NBNP6+1,NEPP6+1);
  pressure7=util->mtx_malloc(NBNP7+1,NEPP7+1);

  // allocate memory for entropy density arrays
  entropyDensity1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  entropyDensity2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  entropyDensity3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  entropyDensity4=util->mtx_malloc(NBNP4+1,NEPP4+1);
  entropyDensity5=util->mtx_malloc(NBNP5+1,NEPP5+1);
  entropyDensity6=util->mtx_malloc(NBNP6+1,NEPP6+1);
  entropyDensity7=util->mtx_malloc(NBNP7+1,NEPP7+1);

  // allocate memory for QGP fraction arrays
  QGPfraction1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  QGPfraction2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  QGPfraction3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  QGPfraction4=util->mtx_malloc(NBNP4+1,NEPP4+1);
  QGPfraction5=util->mtx_malloc(NBNP5+1,NEPP5+1);
  QGPfraction6=util->mtx_malloc(NBNP6+1,NEPP6+1);
  QGPfraction7=util->mtx_malloc(NBNP7+1,NEPP7+1);

  // allocate memory for temperature arrays
  temperature1=util->mtx_malloc(NBNP1+1,NEPP1+1);
  temperature2=util->mtx_malloc(NBNP2+1,NEPP2+1);
  temperature3=util->mtx_malloc(NBNP3+1,NEPP3+1);
  temperature4=util->mtx_malloc(NBNP4+1,NEPP4+1);
  temperature5=util->mtx_malloc(NBNP5+1,NEPP5+1);
  temperature6=util->mtx_malloc(NBNP6+1,NEPP6+1);
  temperature7=util->mtx_malloc(NBNP7+1,NEPP7+1);
  
  // allocate memory for velocity of sound squared arrays
  cs2_1 = util->mtx_malloc(NBNP1+1, NEPP1+1);
  cs2_2 = util->mtx_malloc(NBNP2+1, NEPP2+1);
  cs2_3 = util->mtx_malloc(NBNP3+1, NEPP3+1);
  cs2_4 = util->mtx_malloc(NBNP4+1, NEPP4+1);
  cs2_5 = util->mtx_malloc(NBNP5+1, NEPP5+1);
  cs2_6 = util->mtx_malloc(NBNP6+1, NEPP6+1);
  cs2_7 = util->mtx_malloc(NBNP7+1, NEPP7+1);

  // allocate memory for baryon chemical potential arrays
  // currently always zero
  /*   mu1=mtx_malloc(NBNP1+1,NEPP1+1); */
  /*   mu2=mtx_malloc(NBNP2+1,NEPP2+1); */
  /*   mu3=mtx_malloc(NBNP3+1,NEPP3+1); */
  /*   mu4=mtx_malloc(NBNP4+1,NEPP4+1); */
  
  // read pressure, temperature and chemical potential values
  // files have it backwards, so I start with maximum j and count down
  i=0;
  for(j=NEPP1-1; j>=0; j--)
    {
      fscan_status = fscanf(eos_d1,"%lf",&eps);
      fscan_status = fscanf(eos_d1,"%lf",&pressure1[i][j]);
      fscan_status = fscanf(eos_d1,"%lf",&entropyDensity1[i][j]);
      fscan_status = fscanf(eos_d1,"%lf",&baryonDensity);
      fscan_status = fscanf(eos_d1,"%lf",&QGPfraction1[i][j]);
      fscan_status = fscanf(eos_T1,"%lf",&temperature1[i][j]);
      fscan_status = fscanf(eos_T1,"%lf",&eps); //dummy
      fscan_status = fscanf(eos_T1,"%lf",&eps); //dummy
    }

  for(j=NEPP2-1; j>=0; j--)
    {
      fscan_status = fscanf(eos_d2,"%lf",&eps);
      fscan_status = fscanf(eos_d2,"%lf",&pressure2[i][j]);
      fscan_status = fscanf(eos_d2,"%lf",&entropyDensity2[i][j]);
      fscan_status = fscanf(eos_d2,"%lf",&baryonDensity);
      fscan_status = fscanf(eos_d2,"%lf",&QGPfraction2[i][j]);
      fscan_status = fscanf(eos_T2,"%lf",&temperature2[i][j]);
      fscan_status = fscanf(eos_T2,"%lf",&eps); //dummy
      fscan_status = fscanf(eos_T2,"%lf",&eps); //dummy
    }

  for(j=NEPP3-1; j>=0; j--)
    {
      fscan_status = fscanf(eos_d3,"%lf",&eps);
      fscan_status = fscanf(eos_d3,"%lf",&pressure3[i][j]);
      fscan_status = fscanf(eos_d3,"%lf",&entropyDensity3[i][j]);
      fscan_status = fscanf(eos_d3,"%lf",&baryonDensity);
      fscan_status = fscanf(eos_d3,"%lf",&QGPfraction3[i][j]);
      fscan_status = fscanf(eos_T3,"%lf",&temperature3[i][j]);
      fscan_status = fscanf(eos_T3,"%lf",&eps); //dummy
      fscan_status = fscanf(eos_T3,"%lf",&eps); //dummy
    }

  for(j=NEPP4-1; j>=0; j--)
    {
      fscan_status = fscanf(eos_d4,"%lf",&eps);
      fscan_status = fscanf(eos_d4,"%lf",&pressure4[i][j]);
      fscan_status = fscanf(eos_d4,"%lf",&entropyDensity4[i][j]);
      fscan_status = fscanf(eos_d4,"%lf",&baryonDensity);
      fscan_status = fscanf(eos_d4,"%lf",&QGPfraction4[i][j]);
      fscan_status = fscanf(eos_T4,"%lf",&temperature4[i][j]);
      fscan_status = fscanf(eos_T4,"%lf",&eps); //dummy
      fscan_status = fscanf(eos_T4,"%lf",&eps); //dummy
    }

  for(j=NEPP5-1; j>=0; j--)
    {
      fscan_status = fscanf(eos_d5,"%lf",&eps);
      fscan_status = fscanf(eos_d5,"%lf",&pressure5[i][j]);
      fscan_status = fscanf(eos_d5,"%lf",&entropyDensity5[i][j]);
      fscan_status = fscanf(eos_d5,"%lf",&baryonDensity);
      fscan_status = fscanf(eos_d5,"%lf",&QGPfraction5[i][j]);
      fscan_status = fscanf(eos_T5,"%lf",&temperature5[i][j]);
      fscan_status = fscanf(eos_T5,"%lf",&eps); //dummy
      fscan_status = fscanf(eos_T5,"%lf",&eps); //dummy
    } 

  for(j=NEPP6-1; j>=0; j--)
    {
      fscan_status = fscanf(eos_d6,"%lf",&eps);
      fscan_status = fscanf(eos_d6,"%lf",&pressure6[i][j]);
      fscan_status = fscanf(eos_d6,"%lf",&entropyDensity6[i][j]);
      fscan_status = fscanf(eos_d6,"%lf",&baryonDensity);
      fscan_status = fscanf(eos_d6,"%lf",&QGPfraction6[i][j]);
      fscan_status = fscanf(eos_T6,"%lf",&temperature6[i][j]);
      fscan_status = fscanf(eos_T6,"%lf",&eps); //dummy
      fscan_status = fscanf(eos_T6,"%lf",&eps); //dummy
    } 

  for(j=NEPP7-1; j>=0; j--)
    {
      fscan_status = fscanf(eos_d7,"%lf",&eps);
      fscan_status = fscanf(eos_d7,"%lf",&pressure7[i][j]);
      fscan_status = fscanf(eos_d7,"%lf",&entropyDensity7[i][j]);
      fscan_status = fscanf(eos_d7,"%lf",&baryonDensity);
      fscan_status = fscanf(eos_d7,"%lf",&QGPfraction7[i][j]);
      fscan_status = fscanf(eos_T7,"%lf",&temperature7[i][j]);
      fscan_status = fscanf(eos_T7,"%lf",&eps); //dummy
      fscan_status = fscanf(eos_T7,"%lf",&eps); //dummy
    } 

//test if the reading worked:
/*   for(j=0; j<NEPP1; j++) */
/*     { */
/*       fprintf(stderr,"%lf %lf %lf %lf\n",pressure1[i][j],entropyDensity1[i][j],QGPfraction1[i][j],temperature1[i][j]); */
/*     } */
/*   for(j=0; j<NEPP2; j++) */
/*     { */
/*       fprintf(stderr,"%lf %lf %lf %lf\n",pressure2[i][j],entropyDensity2[i][j],QGPfraction2[i][j],temperature2[i][j]); */
/*     } */
/*   for(j=0; j<NEPP3; j++) */
/*     { */
/*       fprintf(stderr,"%lf %lf %lf %lf\n",pressure3[i][j],entropyDensity3[i][j],QGPfraction3[i][j],temperature3[i][j]); */
/*     } */
/*   for(j=0; j<NEPP4; j++) */
/*     { */
/*       fprintf(stderr,"%lf %lf %lf %lf\n",pressure4[i][j],entropyDensity4[i][j],QGPfraction4[i][j],temperature4[i][j]); */
/*     } */

  fclose(eos_d1);
  fclose(eos_d2);
  fclose(eos_d3);
  fclose(eos_d4);
  fclose(eos_d5);
  fclose(eos_d6);
  fclose(eos_d7);
  fclose(eos_T1);
  fclose(eos_T2);
  fclose(eos_T3);
  fclose(eos_T4);
  fclose(eos_T5);
  fclose(eos_T6);
  fclose(eos_T7);

  music_message.info("Done reading EOS.");

  build_velocity_of_sound_sq_matrix();

  //Allocate and fill arrays for get_s2e()
  //First, create new arrays that will permit to extract information we need from each EOS table in a unified way
  const int number_of_EOS_tables = 7;
  double lowest_eps_list[number_of_EOS_tables] = {EPP1,EPP2,EPP3,EPP4,EPP5,EPP6,EPP7};
  double delta_eps_list[number_of_EOS_tables] = {deltaEPP1,deltaEPP2,deltaEPP3,deltaEPP4,deltaEPP5,deltaEPP6,deltaEPP7};
  int nb_elements_list[number_of_EOS_tables] = {NEPP1,NEPP2,NEPP3,NEPP4,NEPP5,NEPP6,NEPP7};
  double ** s_list[number_of_EOS_tables] = {entropyDensity1,entropyDensity2,entropyDensity3,entropyDensity4,entropyDensity5,entropyDensity6,entropyDensity7};
  //Compute the maximum number of entropyDensity elements (the actual number of entropyDensity elements we will use 
  //will be smaller since we skip duplicate entropyDensities)
  int s_list_rho0_maximal_length=0;
  for(int kk=0;kk<number_of_EOS_tables;kk++) s_list_rho0_maximal_length+=nb_elements_list[kk];
  double last=0.0, tmp_s;
  s_list_rho0_length=0;
  //Allocate memory
  eps_list_rho0= new double [s_list_rho0_maximal_length];
  s_list_rho0= new double [s_list_rho0_maximal_length];

  for(int table_id=0; table_id<number_of_EOS_tables;table_id++) {
	  for(i=0;i<nb_elements_list[table_id];i++) {
		tmp_s=s_list[table_id][0][i];
		if (tmp_s > last) {
			eps_list_rho0[s_list_rho0_length]=(lowest_eps_list[table_id]+i*delta_eps_list[table_id])/hbarc;
			s_list_rho0[s_list_rho0_length]=tmp_s;
			s_list_rho0_length++;
		}
		last=tmp_s;
	  }
  }
  //Initialise the gsl interpolator
  interp_s2e = gsl_interp_alloc(gsl_interp_cspline, s_list_rho0_length);
  gsl_interp_init(interp_s2e, s_list_rho0, eps_list_rho0, s_list_rho0_length);
  accel_s2e = gsl_interp_accel_alloc();
}

double EOS::interpolate_pressure(double e, double rhob)
{
  double p, pa, pb, pa1, pa2, pb1, pb2;
  //use linear interpolation
  int ie1, ie2, inb1, inb2;
  double frace, fracnb;
  
  //if(rhob>0.00001) fprintf(stderr,"rhob=%lf\n", rhob);
  
  e*=hbarc; // in the files epsilon is in GeV/fm^3

  if(e>24.) return (e-4.*.3642)/3./hbarc;

  if(e<EPP2) //use first file for small epsilon values
    {
      if(e<EPP1) 
	{
	  ie1 = 0;
	  ie2 = 1;
	  inb1 = floor((rhob-BNP1)/deltaBNP1);
	  inb2 = floor((rhob-BNP1)/deltaBNP1+1);
	  frace = e/(EPP1);
	}
      else
	{
	  ie1 = floor((e-EPP1)/deltaEPP1);
	  ie2 = floor((e-EPP1)/deltaEPP1+1);
	  inb1 = floor((rhob-BNP1)/deltaBNP1);
	  inb2 = floor((rhob-BNP1)/deltaBNP1+1);
	  frace = (e-(ie1*deltaEPP1+EPP1))/deltaEPP1; 
	}
      
      if(ie1>NEPP1 || inb1>NBNP1 || ie1<0 || inb1<0 )
	{
	  music_message.error("ERROR in inperpolate_pressure. out of range.");
	  fprintf(stderr,"ie1=%d,NEPP1=%d; inb1=%d, NBNP1=%d, e=%f, rhob=%f \n", ie1, NEPP1, inb1, NBNP1,e,rhob);
	  fprintf(stderr,"ie1=%d,BNP1=%lf; inb1=%d, NBNP1=%d, e=%f, rhob=%f \n", ie1, BNP1, inb1, NBNP1,e,rhob);
	 
	  exit(0);
	}
      if(ie2>NEPP1 || inb2>NBNP1 || ie2<0 || inb2<0 )
	{
	  music_message.error("ERROR in inperpolate_pressure. out of range.");
	  fprintf(stderr,"ie2=%d,NEPP1=%d; inb2=%d, NBNP1=%d, e=%f, rhob=%f \n", ie2, NEPP1, inb2, NBNP1,e,rhob);
	  exit(0);
	}

      pa1 = pressure1[inb1][ie1];
      pa2 = pressure1[inb2][ie1];

      fracnb = (rhob-(inb1*deltaBNP1+BNP1))/deltaBNP1; 
      
      pa = pa1*(1-fracnb) + pa2*fracnb;
      
      pb1 = pressure1[inb1][ie2];
      pb2 = pressure1[inb2][ie2];
      
      pb = pb1*(1-fracnb) + pb2*fracnb;
      
     
      if(e<EPP1) 
	{
	  p = pa*(frace);
	  //if (p<0) fprintf(stderr,"pa=%lf\n", pa);
	  //if (p<0) fprintf(stderr,"p=%lf\n", p);
	}
      else
	{
	  p = pa*(1-frace) + pb*frace;
	}
    }
  else // use second file for larger epsilon values
    {
      ie1 = floor((e-EPP2)/deltaEPP2);
      ie2 = floor((e-EPP2)/deltaEPP2+1);
      
      inb1 = floor((rhob-BNP2)/deltaBNP2);
      inb2 = floor((rhob-BNP2)/deltaBNP2+1);
      
      if(ie1>NEPP2 || inb1>NBNP2 || ie1<0 || inb1<0 )
	{
	  fprintf(stderr,"ERROR in inperpolate_pressure. out of range.\n");
	  fprintf(stderr,"ie1=%d,NEPP2=%d; inb1=%d, NBNP2=%d\n", ie1, NEPP2, inb1, NBNP2);
	  fprintf(stderr,"e=%f,EPP2=%f; maxe=%f, deltaEPP2=%f\n", e, EPP2, NEPP2*deltaEPP2+EPP2, deltaEPP2);
	  fprintf(stderr,"ie1=%d,BNP1=%lf; inb1=%d, NBNP1=%d, e=%f, rhob=%f \n", ie1, BNP1, inb1, NBNP1,e,rhob);
	  exit(0);
	}
      if(ie2>NEPP2 || inb2>NBNP2 || ie2<0 || inb2<0 )
	{
	  fprintf(stderr,"ERROR in inperpolate_pressure. out of range.\n");
	  fprintf(stderr,"ie2=%d,NEPP2=%d; inb2=%d, NBNP2=%d\n", ie2, NEPP2, inb2, NBNP2);
	  exit(0);
	}
/*       fprintf(stderr,"ie1=%d\n", ie1); */
/*       fprintf(stderr,"ie2=%d\n", ie2); */
/*       fprintf(stderr,"e=%lf\n", e); */
      
/*       fprintf(stderr,"inb1=%d\n", inb1); */
/*       fprintf(stderr,"inb2=%d\n", inb2); */
/*       fprintf(stderr,"nb=%lf\n", rhob); */
      
      pa1 = pressure2[inb1][ie1];
      pa2 = pressure2[inb2][ie1];
      
      //     fprintf(stderr,"e1=%lf\n", ie1*deltaEPP2+EPP2);
      //fprintf(stderr,"e2=%lf\n", (ie1+1)*deltaEPP2+EPP2);
      
      
      frace = (e-(ie1*deltaEPP2+EPP2))/deltaEPP2; 
      //fprintf(stderr,"frace=%lf\n", frace);
      fracnb = (rhob-(inb1*deltaBNP2+BNP2))/deltaBNP2; 
      
      pa = pa1*(1-fracnb) + pa2*fracnb;
      
      pb1 = pressure2[inb1][ie2];
      pb2 = pressure2[inb2][ie2];
      
      pb = pb1*(1-fracnb) + pb2*fracnb;
      
   /*    fprintf(stderr,"fracnb=%lf\n", fracnb); */
/*       fprintf(stderr,"inb1=%d\n", inb1); */
/*       fprintf(stderr,"inb2=%d\n", inb2); */
/*       fprintf(stderr,"nb1=%lf\n", (inb1*deltaBNP2+BNP2)); */
/*       fprintf(stderr,"nb2=%lf\n", (inb2*deltaBNP2+BNP2)); */
/*       fprintf(stderr,"nb=%lf\n", rhob); */
      
      p = pa*(1-frace) + pb*frace;
      
      //fprintf(stderr,"pa(ie1)=%lf\n", pa1);
      //fprintf(stderr,"pa(ie2)=%lf\n", pa2);
      
      //fprintf(stderr,"pa(e)=%lf\n", pa);  
      //fprintf(stderr,"pb(e)=%lf\n", pb);
      //fprintf(stderr,"p(e)=%lf\n", p);
  
    }
  return p/hbarc;
}

double EOS::interpolate2(double e, double rhob, int selector)
{
    double p, pa, pb;

    //use linear interpolation
    int ie1, ie2, NEps;
    double frace;
    double eps0, deltaEps;
    double **array; 
    //selector = 0 : pressure
    //selector = 1 : temperature
    //selector = 2 : entropy density
    //selector = 3 : QGP fraction
    //selector = 4 : velocity of sound squared
    
    //if(rhob>0.00001) fprintf(stderr,"rhob=%lf\n", rhob);
    //fprintf(stderr,"e=%f\n",e);
    
    e*=hbarc; // in the files epsilon is in GeV/fm^3
     
    //fprintf(stderr,"e=%f\n",e);

    if(e<EPP2) //use first file for small epsilon values
    {
      if(e<EPP1) 
	{
	  ie1 = 0;
	  ie2 = 1;
	  frace = e/(EPP1);
	}
      else
	{
	  ie1 = floor((e-EPP1)/deltaEPP1);
	  ie2 = floor((e-EPP1)/deltaEPP1+1);
	  frace = (e-(ie1*deltaEPP1+EPP1))/deltaEPP1; 
	}
      
      if(ie1>NEPP1)
	{
	  music_message.error("ERROR interpolate2: someting's wrong.");
	  fprintf(stderr,"ie1=%d,NEPP1=%d\n", ie1, NEPP1);
	  exit(0);
	}
      if(ie2>NEPP1)
	{
	  music_message.error("ERROR interpolate2: someting's wrong.");
	  fprintf(stderr,"ie2=%d,NEPP1=%d\n", ie2, NEPP1);
	  exit(0);
	}

      switch (selector) 
	{
	    case 0: array = pressure1; break;
	    case 1: array = temperature1; break;
	    case 2: array = entropyDensity1; break;
	    case 3: array = QGPfraction1; break;
	    case 4: array = cs2_1; break;
        default: music_message.error("ERROR in interpolate2 - selector must be 0,1,2,3, or 4"); exit(1);
	}

      pa = array[0][ie1];
      pb = array[0][ie2];
      
      if(e<EPP1) 
	{
	  p = pa*(frace);
	  //if (p<0) fprintf(stderr,"pa=%lf\n", pa);
	  //if (p<0) fprintf(stderr,"p=%lf\n", p);
	}
      else
	{
	  p = pa*(1-frace) + pb*frace;
	}
    }
    else // use other files for larger epsilon values
    {
      if (e<EPP3)
	{
	    eps0 = EPP2;
	    NEps = NEPP2;
	    deltaEps = deltaEPP2;
	    switch (selector) 
	    {
	        case 0: array = pressure2; break;
	        case 1: array = temperature2; break;
	        case 2: array = entropyDensity2; break;
	        case 3: array = QGPfraction2; break;
	        case 4: array = cs2_2; break;
            default: music_message.error("ERROR in interpolate2 - selector must be 0,1,2,3, or 4"); exit(1);
	    }
	}
      else if (e<EPP4)
	{
	    eps0 = EPP3;
	    NEps = NEPP3;
	    deltaEps = deltaEPP3;
	    switch (selector) 
	    {
	        case 0: array = pressure3; break;
	        case 1: array = temperature3; break;
	        case 2: array = entropyDensity3; break;
	        case 3: array = QGPfraction3; break;
	        case 4: array = cs2_3; break;
            default: music_message.error("ERROR in interpolate2 - selector must be 0,1,2,3, or 4"); exit(1);
	    }
	}
      else if (e<EPP5)
	{
	    eps0 = EPP4;
	    NEps = NEPP4;
	    deltaEps = deltaEPP4;
	    switch (selector) 
	    {
	        case 0: array = pressure4; break;
	        case 1: array = temperature4; break;
	        case 2: array = entropyDensity4; break;
	        case 3: array = QGPfraction4; break;
	        case 4: array = cs2_4; break;
            default: music_message.error("ERROR in interpolate2 - selector must be 0,1,2,3, or 4"); exit(1);
	    }
	}
      else if (e<EPP6)
	{
	    eps0 = EPP5;
	    NEps = NEPP5;
	    deltaEps = deltaEPP5;
	    switch (selector) 
	    {
	        case 0: array = pressure5; break;
	        case 1: array = temperature5; break;
	        case 2: array = entropyDensity5; break;
	        case 3: array = QGPfraction5; break;
	        case 4: array = cs2_5; break;
            default: music_message.error("ERROR in interpolate2 - selector must be 0,1,2,3, or 4"); exit(1);
	    }
	}
      else if (e<EPP7)
	{
	    eps0 = EPP6;
	    NEps = NEPP6;
	    deltaEps = deltaEPP6;
	    switch (selector) 
	    {
	        case 0: array = pressure6; break;
	        case 1: array = temperature6; break;
	        case 2: array = entropyDensity6; break;
	        case 3: array = QGPfraction6; break;
	        case 4: array = cs2_6; break;
            default: music_message.error("ERROR in interpolate2 - selector must be 0,1,2,3, or 4"); exit(1);
	    }
	}
      else
	{
	    eps0 = EPP7;
	    NEps = NEPP7;
	    deltaEps = deltaEPP7;
	    switch (selector) 
	    {
	        case 0: array = pressure7; break;
	        case 1: array = temperature7; break;
	        case 2: array = entropyDensity7; break;
	        case 3: array = QGPfraction7; break;
	        case 4: array = cs2_7; break;
            default: music_message.error("ERROR in interpolate2 - selector must be 0,1,2,3, or 4"); exit(1);
	    }
	}

      ie1 = floor((e-eps0)/deltaEps);
      ie2 = floor((e-eps0)/deltaEps+1);

      if(ie1>NEps)
	{
	  music_message.error("ERROR in inperpolate2. out of range.");
	  fprintf(stderr,"ie1=%d,NEPP2=%d\n", ie1, NEps);
	  fprintf(stderr,"e=%f,eps0=%f; maxe=%f, deltaEps=%f\n", e, eps0, NEps*deltaEps+eps0, deltaEps);
	  exit(0);
	}
      if(ie2>NEps)
	{
	  music_message.error("ERROR in inperpolate2. out of range.");
	  fprintf(stderr,"ie2=%d,NEPP2=%d\n", ie2, NEps);
	  exit(0);
	}

      //if (ie1 ==NEps)
	//cout << " last e=" << e << " T=" << T1*pow(e,T2)/hbarc << endl;

      pa = array[0][ie1];
           
      frace = (e-(ie1*deltaEps+eps0))/deltaEps; 
      
      pb = array[0][ie2];
      
      p = pa*(1-frace) + pb*frace;
      //fprintf(stderr,"p=%f\n",p);
    }
    switch (selector) 
    {
        case 0: p/=hbarc; break;
        case 1: p/=hbarc; break;
        case 2: break;
        case 3: break;
        case 4: break;
        default: music_message.error("ERROR in interpolate2 - selector must be 0,1,2,3, or 4"); exit(1);
    }

    return p;
}


double EOS::get_dpOverde(double e, double rhob)
{
  double dp, pL, pR;
  //use linear interpolation
  double eLeft, eRight;
  
  e*=hbarc;

  //fprintf(stderr,"EPP2=%lf\n", EPP2);
  if(e<EPP2) //use first file for small epsilon values
    {
      eLeft=e-deltaEPP1/2.;
      eRight=e+deltaEPP1/2.;
      if(eLeft<EPP1)
	{
	  eLeft=EPP1;
	  eRight=EPP1+deltaEPP1;
	}
      
      if(eRight>EPP1+NEPP1*deltaEPP1)
	{
	  //eLeft=EPP1+(NEPP1-1)*deltaEPP1;
	  eRight=EPP1+NEPP1*deltaEPP1;
	  //fprintf(stderr,"e=%lf, EPP2=%lf\n", e, EPP2); 
	  //fprintf(stderr,"eLeft=%lf\n", eLeft); 
	  //fprintf(stderr,"eRight=%lf\n", eRight); 
	  //if(eRight-eLeft<deltaEPP1) fprintf(stderr,"delta=%lf, EPP1=%lf\n", eRight-eLeft, deltaEPP1); 
	}

      pL = get_pressure(eLeft/hbarc, rhob);
      pR = get_pressure(eRight/hbarc,  rhob);
      
      dp = pR-pL;
      
       if(dp<0) 
       { 
 	  fprintf(stderr,"1dp/de=%lf\n", dp/((eRight-eLeft)/hbarc)); 
 	  fprintf(stderr,"1pL=%lf\n", pL); 
 	  fprintf(stderr,"1pR=%lf\n", pR); 
 	  fprintf(stderr,"1e=%lf\n", e);  
 	  fprintf(stderr,"1eLeft=%lf\n", eLeft/hbarc);  
 	  fprintf(stderr,"1eRight=%lf\n", eRight/hbarc);  
       } 
      //if(dp<0) dp=0.;
	
      return dp/((eRight-eLeft)/hbarc);
    } 
  else
    {
      eLeft=e-deltaEPP2;
      eRight=e+deltaEPP2;
      if(eLeft<EPP2)
	{
	  eLeft=EPP2;
	  eRight=EPP2+deltaEPP2;
	}
      
      if(eRight>EPP2+NEPP2*deltaEPP2)
	{
	  eLeft=EPP2+(NEPP2-1)*deltaEPP2;
	  eRight=EPP2+NEPP2*deltaEPP2;
	}
      
      pL = get_pressure(eLeft/hbarc, rhob);
      pR = get_pressure(eRight/hbarc, rhob);
      
      dp = pR-pL;
      
      if(dp<0 && pL*hbarc>0.1 && eLeft<24.)
	{
 	  fprintf(stderr,"*2dp=%lf\n", dp/hbarc);
 	  fprintf(stderr,"2dp/de=%lf\n", dp/((eRight-eLeft)/hbarc));
 	  fprintf(stderr,"2pL=%lf\n", pL*hbarc);
 	  fprintf(stderr,"2pR=%lf\n", pR*hbarc);
 	  fprintf(stderr,"2rhob=%lf\n", rhob);
 	  fprintf(stderr,"2e=%lf\n", e);
 	  fprintf(stderr,"EPP2=%lf\n", EPP2);
 	  fprintf(stderr,"2eLeft=%lf\n", eLeft);
 	  fprintf(stderr,"2eRight=%lf\n", eRight);
 	}  
      //fprintf(stderr,"dp/de=%lf\n", dp/((eRight-eLeft)/hbarc));  
      if(dp<0) dp=0.;

      return dp/((eRight-eLeft)/hbarc);
    }
}

double EOS::get_dpOverdrhob(double e, double rhob)
{
  double dp, pL, pR;
  //use linear interpolation
  double rhoLeft, rhoRight;
  
  e*=hbarc;

  //fprintf(stderr,"EPP2=%lf\n", EPP2);
  if(e<EPP2) //use first file for small epsilon values
    {
      rhoLeft=rhob-deltaBNP1/2.;
      rhoRight=rhob+deltaBNP1/2.;
      
      if(rhoLeft<BNP1)
	{
	  rhoLeft=BNP1;
	  rhoRight=BNP1+deltaBNP1;
	}
      
      if(rhoRight>BNP1+NBNP1*deltaBNP1)
	{
	  rhoRight=BNP1+NBNP1*deltaBNP1;
	}
/*       if(rhob>0) */
/* 	{ */
/* 	  fprintf(stderr,"rhoLeft=%lf\n", rhoLeft); */
/* 	  fprintf(stderr,"rhoRight=%lf\n", rhoRight); */
/* 	  fprintf(stderr,"rho=%lf\n", rhob); */
/* 	  fprintf(stderr,"BNP1=%lf\n", BNP1); */
/* 	} */
      pL = get_pressure(e/hbarc, rhoLeft);
      pR = get_pressure(e/hbarc, rhoRight);

      dp = pR-pL;
      return dp/((rhoRight-rhoLeft));
    } 
  else
    {
      rhoLeft=rhob-deltaBNP2/2.;
      rhoRight=rhob+deltaBNP2/2.;
      if(rhoLeft<BNP2)
	{
	  rhoLeft=BNP2;
	}
      
      if(rhoRight>BNP2+NBNP2*deltaBNP2)
	{
	  rhoRight=BNP2+NBNP2*deltaBNP2;
	}

      pL = get_pressure(e/hbarc, rhoLeft); 
      pR = get_pressure(e/hbarc, rhoRight);

/*       fprintf(stderr,"rhoLeft=%lf\n", rhoLeft); */
/*       fprintf(stderr,"rhoRight=%lf\n", rhoRight); */
/*       fprintf(stderr,"rho=%lf\n", rhob); */
/*       fprintf(stderr,"BNP2=%lf\n", BNP2); */
/*       fprintf(stderr,"pL=%lf\n", pL*hbarc); */
/*       fprintf(stderr,"pR=%lf\n", pR*hbarc); */
      
      dp = pR-pL;
      //      fprintf(stderr,"dp=%lf\n", dp);
      //if(dp/(rhoRight-rhoLeft)<1.) fprintf(stderr,"dp/drho2=%lf\n", dp/((rhoRight-rhoLeft)));
      
      return dp/((rhoRight-rhoLeft));
    }
}

double EOS::get_dpOverde2(double e, double rhob)
{

   //The energy has to be in GeV in what follows
   e=e*hbarc;

   double eps0, deltaEps;
   int Neps;

   if(e < EPP1)
   {
       eps0 = 0.0;
       deltaEps = EPP1;
       Neps = 2;
   }
   else if(e<EPP2)
   {
       eps0 = EPP1;
       deltaEps = deltaEPP1;
       Neps = NEPP1;
   }
   else if (e<EPP3)
   {
       eps0 = EPP2;
       deltaEps = deltaEPP2;
       Neps = NEPP2;
   }
   else if (e<EPP4)
   {
       eps0 = EPP3;
       deltaEps = deltaEPP3;
       Neps = NEPP3;
   }
   else if (e<EPP5)
   {
       eps0 = EPP4;
       deltaEps = deltaEPP4;
       Neps = NEPP4;
   }
   else if (e<EPP6)
   {
       eps0 = EPP5;
       deltaEps = deltaEPP5;
       Neps = NEPP5;
   }
   else if (e<EPP7)
   {
       eps0 = EPP6;
       deltaEps = deltaEPP6;
       Neps = NEPP6;
   }
   else 
   {
       eps0 = EPP7;
       deltaEps = deltaEPP7;
       Neps = NEPP7;
   }
   double eps_end = eps0 + (Neps - 1)*deltaEps;

   double eLeft = e - deltaEps*0.5;    // GeV/fm^3
   double eRight = e + deltaEps*0.5;   // GeV/fm^3
  
   // deal with boundary, avoid to exceed the table
   if(eLeft < (eps0 + 1e-6))
       eLeft = eps0;
   if(eRight > (eps_end - 1e-6))
       eRight = eps_end;

   double pL = get_pressure(eLeft/hbarc, rhob);  // 1/fm^4
   double pR = get_pressure(eRight/hbarc, rhob); // 1/fm^4
      
   double dpde = (pR - pL)*hbarc/(eRight - eLeft);
      
   //if(dpde > 1./3.) 
   //{ 
   //    fprintf(stderr, "dp/de=%lf\n", dpde); 
   //    fprintf(stderr, "pL=%lf\n", pL); 
   //    fprintf(stderr, "pR=%lf\n", pR); 
   //    fprintf(stderr, "e=%lf\n", e/hbarc);  
   //    fprintf(stderr, "eLeft=%lf\n", eLeft/hbarc);  
   //    fprintf(stderr, "eRight=%lf\n", eRight/hbarc);  
   //} 
   //if(dpde < 0) 
   //{ 
   //    fprintf(stderr, "dp/de=%lf\n", dpde); 
   //    fprintf(stderr, "pL=%lf\n", pL); 
   //    fprintf(stderr, "pR=%lf\n", pR); 
   //    fprintf(stderr, "e=%lf\n", e);  
   //    fprintf(stderr, "eLeft=%lf\n", eLeft/hbarc);  
   //    fprintf(stderr, "eRight=%lf\n", eRight/hbarc);  
   //} 
   return dpde;
}

double EOS::get_dpOverde3(double e, double rhob)
{
   double eLeft = 0.9*e;
   double eRight = 1.1*e;

   double pL = get_pressure(eLeft, rhob);  // 1/fm^4
   double pR = get_pressure(eRight, rhob); // 1/fm^4
      
   double dpde = (pR - pL)/(eRight - eLeft);
      
   //if(dpde > 1./3.) 
   //{ 
   //    fprintf(stderr, "dp/de=%lf\n", dpde); 
   //    fprintf(stderr, "pL=%lf\n", pL); 
   //    fprintf(stderr, "pR=%lf\n", pR); 
   //    fprintf(stderr, "e=%lf\n", e/hbarc);  
   //    fprintf(stderr, "eLeft=%lf\n", eLeft/hbarc);  
   //    fprintf(stderr, "eRight=%lf\n", eRight/hbarc);  
   //} 
   //if(dpde < 0) 
   //{ 
   //    fprintf(stderr, "dp/de=%lf\n", dpde); 
   //    fprintf(stderr, "pL=%lf\n", pL); 
   //    fprintf(stderr, "pR=%lf\n", pR); 
   //    fprintf(stderr, "e=%lf\n", e);  
   //    fprintf(stderr, "eLeft=%lf\n", eLeft/hbarc);  
   //    fprintf(stderr, "eRight=%lf\n", eRight/hbarc);  
   //} 
   return dpde;
}

double EOS::get_dpOverdrhob2(double e, double rhob)
{
    double local_ed = e*hbarc;      // GeV/fm^3
    double local_rhob = rhob;       // 1/fm^3
    
    double deltaRhob;               // in 1/fm^3
    double rhob_end;
    if (local_ed < EPP2)
    {
        deltaRhob = deltaBNP1;
        rhob_end = BNP1 + (NBNP1 - 1)*deltaRhob;
    }
    else if (local_ed < EPP3)
    {
        deltaRhob = deltaBNP2;
        rhob_end = BNP2 + (NBNP2 - 1)*deltaRhob;
    }
    else if (local_ed < EPP4)
    {
        deltaRhob = deltaBNP3;
        rhob_end = BNP3 + (NBNP3 - 1)*deltaRhob;
    }
    else if (local_ed < EPP5)
    {
        deltaRhob = deltaBNP4;
        rhob_end = BNP4 + (NBNP4 - 1)*deltaRhob;
    }
    else if (local_ed < EPP6)
    {
        deltaRhob = deltaBNP5;
        rhob_end = BNP5 + (NBNP5 - 1)*deltaRhob;
    }
    else if (local_ed < EPP7)
    {
        deltaRhob = deltaBNP6;
        rhob_end = BNP6 + (NBNP6 - 1)*deltaRhob;
    }
    else 
    {
        deltaRhob = deltaBNP7;
        rhob_end = BNP7 + (NBNP7 - 1)*deltaRhob;
    }
    
    double rhobLeft = local_rhob - deltaRhob*0.1;
    double rhobRight = local_rhob + deltaRhob*0.1;

    // avoid exceeding the table
    if(rhobRight > rhob_end)
        rhobRight = rhob_end;
   
    double pL = get_pressure(e, rhobLeft);  // 1/fm^4
    double pR = get_pressure(e, rhobRight); // 1/fm^4
      
    double dpdrho = (pR - pL)/(rhobRight - rhobLeft);  // 1/fm
   
    return dpdrho;   // in 1/fm
}

double EOS::get_cs2(double e, double rhob)
{
    double f;
    if (whichEOS==0)
        f = cs2;
    else if (whichEOS==1)
        f = interpolate(e, rhob, 2);
    else if (whichEOS>=2 && whichEOS < 10)
        f = interpolate2(e, rhob, 4);
    else
    {
        music_message << __func__ <<  "whichEOS = " << whichEOS
                      << "is out of range!";
        music_message.flush("error");
        exit(0);
    }
    return f;
}

void EOS::build_velocity_of_sound_sq_matrix()
{
    // whiichEOS == 0 : ideal gas EoS, cs^2 = 1/3 always
    if (whichEOS==1)  // 2 tables
    {
        fill_cs2_matrix(EPP1, deltaEPP1, NEPP1, BNP1, deltaBNP1, NBNP1, cs2_1);
        fill_cs2_matrix(EPP2, deltaEPP2, NEPP2, BNP2, deltaBNP2, NBNP2, cs2_2);
    }
    else if (whichEOS>=2 && whichEOS < 10)  // 7 tables
    {
        fill_cs2_matrix(EPP1, deltaEPP1, NEPP1, 0.0, 0.0, 1, cs2_1);
        fill_cs2_matrix(EPP2, deltaEPP2, NEPP2, 0.0, 0.0, 1, cs2_2);
        fill_cs2_matrix(EPP3, deltaEPP3, NEPP3, 0.0, 0.0, 1, cs2_3);
        fill_cs2_matrix(EPP4, deltaEPP4, NEPP4, 0.0, 0.0, 1, cs2_4);
        fill_cs2_matrix(EPP5, deltaEPP5, NEPP5, 0.0, 0.0, 1, cs2_5);
        fill_cs2_matrix(EPP6, deltaEPP6, NEPP6, 0.0, 0.0, 1, cs2_6);
        fill_cs2_matrix(EPP7, deltaEPP7, NEPP7, 0.0, 0.0, 1, cs2_7);
    }
    else
    {
        music_message << __func__ <<  "whichEOS = " << whichEOS
                      << "is out of range!";
        music_message.flush("error");
        exit(0);
    }
}

void EOS::fill_cs2_matrix(double e0, double de, int ne, double rhob0, double drhob, int nrhob, double** cs2_ptr)
{
    for(int i = 0; i < nrhob; i++)
    {
        double rhob_local = rhob0 + i*drhob;
        for(int j = 0; j < ne; j++)
        {
            double epsilon_local = (e0 + j*de)/hbarc;
            double cs2_local = calculate_velocity_of_sound_sq(epsilon_local, rhob_local);
            cs2_ptr[i][j] = cs2_local;
        }
    }
}

double EOS::calculate_velocity_of_sound_sq(double e, double rhob)
{
    double v_min = 0.01;
    double v_max = 1./3;
    double dpde = p_e_func(e, rhob);
    double dpdrho = p_rho_func(e, rhob);
    double pressure = get_pressure(e, rhob);
    double v_sound = dpde + rhob/(e + pressure + 1e-15)*dpdrho;

    if(v_sound < v_min)
        v_sound = v_min;
    if(v_sound > v_max)
    {
        //fprintf(stderr, "EOS::get_velocity_of_sound:velocity of sound is larger than 1/3!\n");
        //fprintf(stderr, "v_sound = %lf \n", v_sound);
        //fprintf(stderr, "e = %lf \n", e*hbarc);
        //fprintf(stderr, "rhob = %lf \n", rhob);
        //fprintf(stderr, "pressure = %lf \n", pressure);
        //fprintf(stderr, "dP/de = %lf \n", dpde);
        //fprintf(stderr, "dP/drho = %lf \n", dpdrho);
        //exit(0);
        v_sound = v_max;
    }
    /*
    if(v_sound < 0.)
    {
        fprintf(stderr, "EOS::get_velocity_of_sound:velocity of sound is negative!\n");
        fprintf(stderr, "v_sound = %lf \n", v_sound);
        fprintf(stderr, "e = %lf \n", e*hbarc);
        fprintf(stderr, "rhob = %lf \n", rhob);
        fprintf(stderr, "pressure = %lf \n", pressure);
        fprintf(stderr, "dP/de = %lf \n", dpde);
        fprintf(stderr, "dP/drho = %lf \n", dpdrho);
    }
    */
    return(v_sound);
}

double EOS::get_pressure(double e, double rhob)
// return pressure in [1/fm^4]
{
    double f;
    if (whichEOS==0)
        f = cs2*e;
    else if (whichEOS==1)
        f = interpolate_pressure(e,rhob);
    else if (whichEOS>=2 && whichEOS < 10)
        f = interpolate2(e,rhob,0);    //selector 0 means get pressure 
    else
    {
        music_message << __func__ <<  "whichEOS = " << whichEOS
                      << "is out of range!";
        music_message.flush("error");
        exit(0);
    }
    return f;
}/* get_pressure */


double EOS::p_rho_func(double e, double rhob)
// return dP/drho_b (in 1/fm)
{
    double f;
    if (whichEOS==0) 
        f = 0.0;
    else if (whichEOS==1)
        f = get_dpOverdrhob(e, rhob);
    else if (whichEOS>=2 && whichEOS < 10)
        f = 0.0;
    else
    {
        music_message << __func__ <<  "whichEOS = " << whichEOS
                      << "is out of range!";
        music_message.flush("error");
        exit(0);
    }
    return f;
}/* p_rho_func */

double EOS::p_e_func(double e, double rhob)
// return dP/de
{
    double f;
    if (whichEOS==0)
        f = cs2;
    else if (whichEOS==1)
        f = get_dpOverde(e, rhob);
    else if (whichEOS>=2 && whichEOS < 10)
        f = get_dpOverde2(e, rhob);
    else
    {
        music_message << __func__ <<  "whichEOS = " << whichEOS
                      << "is out of range!";
        music_message.flush("error");
        exit(0);
    }
    return f;
}/* p_e_func */


double EOS::interpolate(double e, double rhob, int selector)
{
  // selector=0 : return temperature
  // selector=1 : return baryon chemical potential
  // selector=2 : return cs^2
  double T, pa, pb, pa1, pa2, pb1, pb2;
  //use linear interpolation
  int ie1, ie2, inb1, inb2;
  double frace, fracnb;
  
  e*=hbarc; // in the files epsilon is in GeV/fm^3

  if(e>24.)
    {
      if (selector==0)
	{
	  return pow(((e-.3642)*120./(pow(PI,2.))/169.*(pow(hbarc,3.))),0.25)/hbarc;
	}    
      else if (selector==1)
	{
	  return 0.;
	}
    }
  if(e<EPP2) //use first file for small epsilon values
    {
      if(e<EPP1) 
	{
	  ie1 = 0;
	  ie2 = 1;
	  inb1 = floor((rhob-BNP1)/deltaBNP1);
	  inb2 = floor((rhob-BNP1)/deltaBNP1+1);
	  frace = e/(EPP1);
  	}
      else
	{
	  ie1 = floor((e-EPP1)/deltaEPP1);
	  ie2 = floor((e-EPP1)/deltaEPP1+1);
	  inb1 = floor((rhob-BNP1)/deltaBNP1);
	  inb2 = floor((rhob-BNP1)/deltaBNP1+1);
	  frace = (e-(ie1*deltaEPP1+EPP1))/deltaEPP1; 
  	}
      
      if(ie1>NEPP1 || inb1>NBNP1)
	{
	  music_message.error("ERROR in inperpolate_temperature. out of range.");
	  fprintf(stderr,"ie1=%d,NEPP1=%d; inb1=%d, NBNP1=%d\n", ie1, NEPP1, inb1, NBNP1);
	  exit(0);
	}
      if(ie2>NEPP1 || inb2>NBNP1)
	{
	  music_message.error("ERROR in inperpolate_temperature. out of range.");
	  fprintf(stderr,"ie2=%d,NEPP1=%d; inb2=%d, NBNP1=%d\n", ie2, NEPP1, inb2, NBNP1);
	  exit(0);
	}
      
      if (selector==0)
	{
	  pa1 = temperature1[inb1][ie1];
	  pa2 = temperature1[inb2][ie1];
	}
      else if (selector==1)
	{
	  pa1 = mu1[inb1][ie1];
	  pa2 = mu1[inb2][ie1];
	}
      else if (selector == 2)
      {
	  pa1 = cs2_1[inb1][ie1];
	  pa2 = cs2_1[inb2][ie1];
      }
   	else
      {
        music_message << "selector = " << selector << " is out of range.";
        music_message.flush("error");
        exit(0);
      }

      if (pa1<0) pa1=0.;
      if (pa2<0) pa2=0.;
    
      fracnb = (rhob-(inb1*deltaBNP1+BNP1))/deltaBNP1; 
      
      pa = pa1*(1-fracnb) + pa2*fracnb;
      
      if (selector==0)
	{
	  pb1 = temperature1[inb1][ie2];
	  pb2 = temperature1[inb2][ie2];
	}
      else if (selector==1)
	{
	  pb1 = mu1[inb1][ie2];
	  pb2 = mu1[inb2][ie2];
	}
      else if (selector == 2)
      {
	  pb1 = cs2_1[inb1][ie2];
	  pb2 = cs2_1[inb2][ie2];
      }
   	else
      {
        music_message << "selector = " << selector << " is out of range.";
        music_message.flush("error");
        exit(0);
      }

      if (pb1<0) pb1=0.;
      if (pb2<0) pb2=0.;
      
      pb = pb1*(1-fracnb) + pb2*fracnb;
     
      //fprintf(stderr,"fracnb=%lf\n", fracnb);

      if(e<EPP1) 
	{
	  T = pa*(frace);
	}
      else 
	{
	  T = pa*(1-frace) + pb*frace;
	  
	}
    }
  else // use second file for larger epsilon values
    {
      ie1 = floor((e-EPP2)/deltaEPP2);
      ie2 = floor((e-EPP2)/deltaEPP2+1);
      
      inb1 = floor((rhob-BNP2)/deltaBNP2);
      inb2 = floor((rhob-BNP2)/deltaBNP2+1);
      
      if(ie1>NEPP2 || inb1>NBNP2)
	{
	  music_message.error("ERROR in inperpolate_temperature. out of range.");
	  fprintf(stderr,"ie1=%d,NEPP2=%d; inb1=%d, NBNP2=%d\n", ie1, NEPP2, inb1, NBNP2);
	  fprintf(stderr,"e=%f,EPP2=%f; maxe=%f, deltaEPP2=%f\n", e, EPP2, NEPP2*deltaEPP2+EPP2, deltaEPP2);
	  exit(0);
	}
      if(ie2>NEPP2 || inb2>NBNP2)
	{
	  music_message.error("ERROR in inperpolate_temperature. out of range.");
	  fprintf(stderr,"ie2=%d,NEPP2=%d; inb2=%d, NBNP2=%d\n", ie2, NEPP2, inb2, NBNP2);
	  exit(0);
	}

      if (selector==0)
	{
	  pa1 = temperature2[inb1][ie1];
	  pa2 = temperature2[inb2][ie1];
	}
      else if (selector==1)
	{
	  pa1 = mu2[inb1][ie1];
	  pa2 = mu2[inb2][ie1];
	}
      else if (selector == 2)
      {
	  pa1 = cs2_1[inb1][ie1];
	  pa2 = cs2_1[inb2][ie1];
      }
   	else
      {
        music_message << "selector = " << selector << " is out of range.";
        music_message.flush("error");
        exit(0);
      }
      
      if (pa1<0) pa1=0.;
      if (pa2<0) pa2=0.;
      
      frace = (e-(ie1*deltaEPP2+EPP2))/deltaEPP2; 
      fracnb = (rhob-(inb1*deltaBNP2+BNP2))/deltaBNP2; 
      
      pa = pa1*(1-fracnb) + pa2*fracnb;
      
      if (selector==0)
	{
	  pb1 = temperature2[inb1][ie2];
	  pb2 = temperature2[inb2][ie2];
	}     
      else if (selector==1)
	{
	  pb1 = mu2[inb1][ie2];
	  pb2 = mu2[inb2][ie2];
	}
      else if (selector == 2)
      {
	  pb1 = cs2_1[inb1][ie2];
	  pb2 = cs2_1[inb2][ie2];
      }
   	else
      {
        music_message << "selector = " << selector << " is out of range.";
        music_message.flush("error");
        exit(0);
      }

      if (pb1<0) pb1=0.;
      if (pb2<0) pb2=0.;
      
      pb = pb1*(1-fracnb) + pb2*fracnb;

      T = pa*(1-frace) + pb*frace;
    }
  return T/hbarc;
}

double EOS::T_from_eps_ideal_gas(double eps)
{

	//Define number of colours and of flavours
 	const double Nc=3, Nf=2.5;
	
	return pow(90.0/M_PI/M_PI*(eps/3.0)/(2*(Nc*Nc-1)+7./2*Nc*Nf),.25);

}
double EOS::s2e_ideal_gas(double s)
{

	//Define number of colours and of flavours
 	const double Nc=3, Nf=2.5;

	//e=T*T*T*T*(M_PI*M_PI*3.0*(2*(Nc*Nc-1)+7./2*Nc*Nf)/90.0);
	//s = 4 e / (3 T)
	//s =4/3 T*T*T*(M_PI*M_PI*3.0*(2*(Nc*Nc-1)+7./2*Nc*Nf)/90.0);
	//T = pow(3. * s / 4. / (M_PI*M_PI*3.0*(2*(Nc*Nc-1)+7./2*Nc*Nf)/90.0), 1./3.);
	return 3. / 4. * s * pow(3. * s / 4. / (M_PI*M_PI*3.0*(2*(Nc*Nc-1)+7./2*Nc*Nf)/90.0), 1./3.); //in 1/fm^4

}

double EOS::get_entropy(double epsilon, double rhob)
// return entropy density in [1/fm^3]
{
    double f;
    double P, T, mu;
    
    if (whichEOS >=0 && whichEOS <=1)
    {
        P = get_pressure(epsilon, rhob);
        T = get_temperature(epsilon,rhob);
        mu = get_mu(epsilon, rhob);
   
        if (T!=0)
            f = (epsilon + P - mu*rhob)/(T + 1e-15);
        else
            f = 0.;
    }
    else if (whichEOS >= 2 && whichEOS < 10)
        f = interpolate2(epsilon,rhob,2);
    else
    {
        music_message << __func__ << "whichEOS = " << whichEOS 
                      << " is out of range!";
        music_message.flush("error");
        exit(0);
    }
    return f;
}/* get_entropy */

double EOS::ssolve(double e, double rhob, double s)
{
  // takes e in GeV/fm^3 and passes it on in 1/fm^4 ...
  double P, T, mu;
  P = get_pressure(e/hbarc, rhob);
  T = get_temperature(e/hbarc, rhob);
  mu = get_mu(e/hbarc, rhob);
      
//   fprintf(stderr,"T=%f\n",T*hbarc);
//   fprintf(stderr,"P=%f\n",P*hbarc);
//   fprintf(stderr,"e=%f\n",e);
//   fprintf(stderr,"s=%f\n",s);
//   fprintf(stderr,"ssolve=%f\n",T*s-e/hbarc-P+mu*rhob);

  return T*s-e/hbarc-P+mu*rhob;
}

double EOS::Tsolve(double e, double rhob, double T)
{
  // takes e in GeV/fm^3 and passes it on in 1/fm^4 ...
  return T-get_temperature(e/hbarc,rhob);
}

double EOS::findRoot(double (EOS::*func)(double, double, double), double rhob, double s, double e1, double e2, double eacc)
{
  int j, jmax;
  jmax=40;
  double emid, de, value, f, fmid;
  fmid = (this->*func)(e2, rhob, s);
  f = (this->*func)(e1, rhob, s);
 
  //  fprintf(stderr,"fmid=%f\n",fmid);
  //fprintf(stderr,"fabs(f)=%f\n",fabs(f));
  //fprintf(stderr,"eacc=%f\n",eacc);
  
  if(f*fmid>=0)
    {
      if( fabs(f) < eacc )
	{
	 return 0.;
	}
      fprintf(stderr,"root must be bracketed in findRoot\n");
      fprintf(stderr,"f=%f\n",f);
      fprintf(stderr,"fmid=%f\n",fmid);
      fprintf(stderr,"fabs(f)=%f\n",fabs(f));
      fprintf(stderr,"eacc=%f\n",eacc);
    }
     
  if (f<0)
    {
      value=e1;
      de=e2-e1;
    }
  else
    {
      value=e2;
      de=e1-e2;
    }
  for(j=1; j<=jmax; j++)
    {
      de*=0.5;
      emid = value+de;
      fmid = (this->*func)(emid, rhob, s);
      //fprintf(stderr,"fmid(emid)=%f\n",fmid);
      //fprintf(stderr,"emid=%f\n",emid);
      //fprintf(stderr,"value=%f\n",value);
      if (fmid<=0.) value = emid;
      if (fabs(de)<eacc || fmid==0.) return value/hbarc;
    }
  fprintf(stderr,"too many bisections in findRoot\n");
  return 0.;
}


double EOS::get_temperature(double eps, double rhob)
// return temperature in [1/fm]
{
    double T;
    if (whichEOS==0)
        T = T_from_eps_ideal_gas(eps);
    else if (whichEOS==1)
        T = interpolate(eps, rhob, 0);
    else if (whichEOS>=2 && whichEOS < 10)
        T = interpolate2(eps, rhob, 1);
    else
    {
        music_message << __func__ << "whichEOS = " << whichEOS 
                      << " is out of range!";
        music_message.flush("error");
        exit(0);
    }
    return max(T, 1e-15);
}


double EOS::get_mu(double eps, double rhob)
// return mu_B in [1/fm]
{
    double mu;
    if (whichEOS==0)
        mu=0.0;
    else if (whichEOS==1)
        mu = 0.0;
    else if (whichEOS>=2 && whichEOS < 10)
        mu = 0.0;
    else
    {
        music_message << __func__ << "whichEOS = " << whichEOS 
                      << " is out of range!";
        music_message.flush("error");
        exit(0);
    }
    return mu;
}


double EOS::get_muS(double eps, double rhob)
// return mu_B in [1/fm]
{
    double mu;
    if (whichEOS==0)
        mu=0.0;
    else if (whichEOS==1)
        mu = 0.0;
    else if (whichEOS>=2 && whichEOS < 10)
        mu = 0.0;
    else
    {
        music_message << __func__ << "whichEOS = " << whichEOS 
                      << " is out of range!";
        music_message.flush("error");
        exit(0);
    }
    return mu;
}

double EOS::get_qgp_frac(double eps, double rhob) {

 double frac;

 if (whichEOS==0)
   {
     frac=-1.0;
   }
 else if (whichEOS==1)
   {
     frac=(eps*hbarc-0.45)/(1.6-0.45); // e=(1-QGPfrac)*e_H + QGPfrac*e_QGP in the mixed phase
     if (frac>1.) frac = 1;
     else if (frac<0.) frac=0.;
   }
 else if (whichEOS>=2)
   {
     frac = interpolate2(eps, rhob, 3);
   }
   	else {
        music_message << __func__ << "whichEOS = " << whichEOS 
                      << " is out of range!";
        music_message.flush("error");
        exit(0);
    }

  return frac;

}


double EOS::get_s2e(double s, double rhob)
{ //s - entropy density in 1/fm^3
   double e; //epsilon - energy density
   int status;
   
   if (whichEOS==0)
      e = s2e_ideal_gas(s);
   else if (whichEOS>=2 && whichEOS<=6)
   {
      if(s < entropyDensity1[0][0])
      {
         double frac = s/entropyDensity1[0][0];
         e = frac*EPP1;
      }
      else
      {
         status = gsl_interp_eval_e(interp_s2e, s_list_rho0, eps_list_rho0, s, accel_s2e, &e);
         if (status == GSL_EDOM)
         {
            music_message << "Error: can't get energy from entropy, entropy s="
                          << s
                          <<" fm^-3 is outside the current tabulation of the EOS...\n";
            music_message.flush("error");
            exit(1);
         }
      }
   }
   else
   {
      music_message << __func__ << "whichEOS = " << whichEOS
                    << "is out of range";
      music_message.flush("error");
      exit(0);
   }
   return e; //in 1/fm^4
}

