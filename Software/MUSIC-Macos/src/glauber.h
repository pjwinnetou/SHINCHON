// Glauber routines
// Copyright (C) 2017 Sangyong Jeon, Joseph Kapusta, Chun Shen & Björn Schenke

#ifndef glauber_h //avoids multiple inclusions of the header file
#define glauber_h

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "data.h"
#include "./pretty_ostream.h"
#include <sstream>
#include <fstream>
#include <gsl/gsl_rng.h>



#define TOL (1.0e-6)
#define tiny (1.0e-10)
#define limit 10000


//! This class generates MC-Glauber initial conditions for hydrodynamics

class Glauber{

 private:
    Util *util;
    pretty_ostream music_message;
    init_data *parameter_list;
  
 public:
    typedef double (*ptr_func)(double);
  
    typedef struct nucleus {
        string name;
        double A;
        double Z;
        int AnumFunc;
        int AnumFuncIntegrand;
        int DensityFunc;
        double w_WS;
        double a_WS;
        double R_WS;
        double rho_WS;
    } Nucleus;
  
    typedef struct data {
        double SigmaNN; 
        Nucleus Target;
        Nucleus Projectile;
        double SCutOff;
        int InterMax; 
    } Data;
  
    double AnumR, NuInS_S;
    Nucleus *Nuc_WS; 
    Data LexusData;  
    ptr_func tempFunc;
    double b; // impact parameter
    double currentTAB;
    double currentA;

 public:
    Glauber(init_data *parameter_list_in);
    ~Glauber();
    double nucleusA() const {return(currentA);}
    int IsFile(char *file_name);
    int setNucleusParameters(Nucleus *nucleus, string nucleus_name);
    void PrintLexusData();
    void PrintNucleusData(Nucleus *nucleus);
    int LinearFindXorg(double x, double *Vx, int ymax);
    double FourPtInterpolate(double x, double *Vx, double *Vy, double h, int x_org, int ymax);
    void MakeCoeff(double *a, double *b, double *c, double *d, 
      	 double *Vy, double *Vx, double h, int x_org);
    double VInterpolate(double x, double *Vx, double *Vy, int ymax);
    int FindXorg(double x, double *Vx, int ymax);
    double *MakeVx(double down, double up, int maxi_num);
    double *MakeVy(string st, double *vx, int maxi_num);
    double *ReadInVx(char *, int maxi_num, int quiet);
    double *ReadInVy(char *, int maxi_num, int quiet);
    
    double InterNuPInSP(double s);
    double InterNuTInST(double s);
    void CalcRho(Nucleus *nucleus);
    double NuInS(double s);
    
    double Anum3Fermi(double R_WS);
    double Anum3FermiInt(double xi);
    double NuInt3Fermi(double xi);
    double Anum3Gauss(double R_WS);
    double Anum3GaussInt(double xi);
    double NuInt3Gauss(double xi);
    double Anum2HO(double R_WS);
    double Anum2HOInt(double xi);
    double NuInt2HO(double xi);
    
    double integral (int id, double down, double up, double tol, int *count);
    double qnc7(int id, double tol, double down, double dx, double *f_of, 
            double pre_sum, double area, int *count);
    double OLSIntegrand(double s);
    double TAB();
    double PAB(double x, double y);
    void initGlauber(double SigmaNN, string Target, string Projectile, double b,
                     int size, int rank);
    double areaTA(double x, double A);
    MCGlauberReturnValue SampleTARejection(gsl_rng *random);
};
#endif

