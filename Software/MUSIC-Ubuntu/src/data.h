// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

#ifndef SRC_DATA_H
#define SRC_DATA_H

#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

struct MCGlauberReturnValue {
    double x;
    double y;
    int collided;
};


//! this data structure contains all the parameters used in the simulation general parameters

typedef struct init_data {
    double **gmunu; /* metric */

    int mode;               //!< 1: do everything (mode 2, 3 and 4); 
                            //!< 2: do hydro evolution only; 
                            //!< 3: do calculation of thermal spectra only;
                            //!< 4: do resonance decays only
    int seed;
    int echo_level;

    //! MPI control
    int size;
    int rank;

    //! parameters for initial condition
    string Target;
    string Projectile;
    string initName;     // filename of the initial density profile

    double epsilon0;
    double rhoB0;
    double eta_fall_off;
    double eta_flat;
    double SigmaNN;
    double b;

    int Initial_profile;
    int initializeEntropy;  //!< 1: MC-Glauber outputs the entropy density
                            //!< 0: MC-Glauber outputs the energy denisty
    double hard;            //!< fraction of binary collisons scaling in Mc-Glauber
                            //!< initial distribution
    double bmin;            //!< minimal and maximal impact parameter
    double bmax;            //!< when using MC-Glauber (Initial_Distribution 3), for sampling
    double sigma0;
    double sFactor;
    int initial_eta_profile;

    // grid information
    int nx;
    int ny;
    int neta;
    int nt;

    double x_size;       /*!< in fermi -x_size/2 < x < x_size/2 */
    double y_size;       /*!< in fermi, ditto */
    double eta_size;     /*!< ditto */
    double tau_size;     /*!< tau_0 < tau < tau0+tau_size  */
    
    double delta_x;
    double delta_y;
    double delta_eta;
    double delta_tau;
    
    // hydro parameters
    bool boost_invariant;
    double tau0;
    int whichEOS;
    int store_hydro_info_in_memory;
    int outputEvolutionData;    //!< whether to output the evolution data
    int outputBinaryEvolution;  //!< whether to output evolution data in binary
                                //!< format (1) or in text format (0)
    int output_evolution_every_N_timesteps;
    int output_evolution_every_N_x;
    int output_evolution_every_N_y;
    int output_evolution_every_N_eta;
    bool output_hydro_params_header;

    int turn_on_shear;
    int turn_on_bulk;
    double shear_to_s;
    int viscosity_flag;
    int T_dependent_shear_to_s;
    int include_second_order_terms;

    // parameters for numerical solver
    int reconst_type;
    int rk_order;
    double minmod_theta;
  

    // parameters for Cooper-Fyre
    int doFreezeOut;  
    int pseudofreeze;
    double epsilonFreeze;   //!< freeze-out energy density in GeV/fm^3
    int freezeOutMethod;    //!< 2: 4d-triangulation

    int include_deltaf;
    int include_deltaf_qmu;
    int include_deltaf_bulk;

    double TFO;            //!< freeze-out temperature. Used if useEpsFO=0
    int useEpsFO;          //!< if 1, use energy density value to define freeze
                           //!< out condition, if 0 use temperature in TFO

    // skip step for freeze out surface finder to save CPU-time/memory/diskspace
    int facTau;
    int fac_x;
    int fac_eta;

    // for calculation of spectra
    int NumberOfParticlesToInclude;  //!< # of resonances to include.
                                     //!< maximum=319 (all up to 2 GeV) 
    int particleSpectrumNumber;      //!< number of particle for which
                                     //!< the spectrum is to be computed.
                                     //!< 0: all particles
    double max_pseudorapidity;
    int pseudo_steps;
    int phi_steps;
    double min_pt;
    double max_pt;
    int pt_steps;
  
    // parameters for mode 13 and 14
    double dNdy_y_min;    //!< rapidity range for dN/dy as a function of y
    double dNdy_y_max;
    //! pseudo-rapidity range for dN/dy as a function of y
    double dNdy_eta_min;
    double dNdy_eta_max;
    int dNdy_nrap;
    //! the integrated rapidity range for dN/dypTdpT
    double dNdyptdpt_y_min;
    double dNdyptdpt_y_max;
    double dNdyptdpt_eta_min;
    double dNdyptdpt_eta_max;

} InitData;


#endif
