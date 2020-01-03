// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

#ifndef INIT_H
#define INIT_H

#include <stdio.h>
#include <iostream>
#include "data.h"
#include "grid.h"
#include "glauber.h"
#include "./pretty_ostream.h"
#include <vector>
#include <time.h>
#include <cmath>
#include <gsl/gsl_rng.h>

//! This class initializes hydrodynamics on a grid

class Init
{
    private:
        gsl_rng * random; 
        EOS *eos;
        Util * util;
        Glauber *glauber;
        pretty_ostream music_message;

        // list of x and y coordinates of nucleons in nucleus A      
        vector<MCGlauberReturnValue> nucleusA;  
        // list of x and y coordinates of nucleons in nucleus B 
        vector<MCGlauberReturnValue> nucleusB;  

        // support for JETSCAPE
        std::vector<double> initial_entropy_density;
        std::vector<double> initial_energy_density;
        std::vector<double> initial_u_tau;
        std::vector<double> initial_u_x;
        std::vector<double> initial_u_y;
        std::vector<double> initial_u_eta;
        std::vector<double> initial_pi_00;
        std::vector<double> initial_pi_01;
        std::vector<double> initial_pi_02;
        std::vector<double> initial_pi_03;
        std::vector<double> initial_pi_11;
        std::vector<double> initial_pi_12;
        std::vector<double> initial_pi_13;
        std::vector<double> initial_pi_22;
        std::vector<double> initial_pi_23;
        std::vector<double> initial_pi_33;
        std::vector<double> initial_bulk_pi;

    public:
        Init(EOS *eos, Glauber* glauber);  //constructor
        ~Init();  //destructor

        void sampleTA();
        void get_entropy_density_vector(vector<double> entropy_density_in);
        void get_pre_equilibrium_vectors(std::vector<double> e_in,
                                         std::vector<double> u_tau_in,
                                         std::vector<double> u_x_in,
                                         std::vector<double> u_y_in,
                                         std::vector<double> u_eta_in,
                                         std::vector<double> pi_00_in,
                                         std::vector<double> pi_01_in,
                                         std::vector<double> pi_02_in,
                                         std::vector<double> pi_03_in,
                                         std::vector<double> pi_11_in,
                                         std::vector<double> pi_12_in,
                                         std::vector<double> pi_13_in,
                                         std::vector<double> pi_22_in,
                                         std::vector<double> pi_23_in,
                                         std::vector<double> pi_33_in,
                                         std::vector<double> Bulk_pi_in);

        void InitArena(InitData *DATA, Grid ****arena, Grid ****Lneighbor, 
                       Grid ****Rneighbor, int size, int rank);
        void LinkNeighbors(InitData *DATA, Grid ****arena, int size, int rank);
        int InitTJb(InitData *DATA, Grid ****arena, Grid ****Lneighbor, 
                    Grid ****Rneighbor, int size, int rank);

        // The following two functions return T_A and T_B. 
        // Normalization: \int r T_A(r) dr dphi = A
        double TATarget(InitData *DATA, double r);
        double TAProjectile(InitData *DATA, double r);
    
        double eta_profile_normalisation(InitData *DATA, double eta);
};

#endif
