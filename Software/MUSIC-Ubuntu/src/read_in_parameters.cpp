// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen


#include <iostream>
#include <cstring>
#include "./read_in_parameters.h"

using namespace std;

ReadInParameters::ReadInParameters() {
    util = new Util();
    emoji_face = new emoji();
}

ReadInParameters::~ReadInParameters() {
    delete util;
    delete emoji_face;
}


void ReadInParameters::read_in_parameters(InitData *parameter_list,
                                          string input_file_in) {
    input_file = input_file_in;
    // this function reads in parameters
    int m, n;
    string tempinput;

    // Impact_parameter: in fm, >=0
    double tempb = 0;
    tempinput = util->StringFind4(input_file, "Impact_parameter");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempb;
    parameter_list->b = tempb;

    // echo_level controls the mount of
    // warning message output during the evolution
    double temp_echo_level = 9;
    tempinput = util->StringFind4(input_file, "echo_level");
    if(tempinput != "empty") istringstream ( tempinput ) >> temp_echo_level;
    parameter_list->echo_level = temp_echo_level;

    // Target, Projectile
    string tempTarget = "Pb";
    tempinput = util->StringFind4(input_file, "Target");
    if (tempinput != "empty")
        tempTarget.assign(tempinput);
    parameter_list->Target.assign(tempTarget);
    string tempProjectile = "Pb";
    tempinput = util->StringFind4(input_file, "Projectile");
    if (tempinput != "empty")
        tempProjectile.assign(tempinput);
    parameter_list->Projectile.assign(tempProjectile);

    // SigmaNN: nucleon-nucleon cross section in mb
    double tempSigmaNN = 60.;
    tempinput = util->StringFind4(input_file, "SigmaNN");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempSigmaNN;
    parameter_list->SigmaNN = tempSigmaNN;
    
    // Initial_proinput_file: 
    int tempInitial_profile = 1;
    tempinput = util->StringFind4(input_file, "Initial_profile");
    if(tempinput != "empty") istringstream(tempinput) >> tempInitial_profile;
    parameter_list->Initial_profile = tempInitial_profile;

    // boost-invariant
    int temp_boost_invariant = 0;
    tempinput = util->StringFind4(input_file, "boost_invariant");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_boost_invariant;
    if (temp_boost_invariant == 0) {
        parameter_list->boost_invariant = false;
    } else {
        parameter_list->boost_invariant = true;
    }

    //Select the profile to use in eta for the energy/entropy initialisation
    //1 for Hirano's central plateau + Gaussian decay
    //2 for a Woods-Saxon profile
    int tempinitial_eta_profile = 1;
    tempinput = util->StringFind4(input_file, "initial_eta_profile");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempinitial_eta_profile;
    parameter_list->initial_eta_profile = tempinitial_eta_profile;
    
    
    //initialize_with_entropy: only for Glauber initialization
    //0: scale with energy density
    //1: scale with entropy density
    int tempinitializeEntropy = 0;
    tempinput = util->StringFind4(input_file, "initialize_with_entropy");
    if(tempinput != "empty") istringstream ( tempinput ) >> tempinitializeEntropy;
    parameter_list->initializeEntropy = tempinitializeEntropy;
    
    // T_freeze: freeze out temperature
    // only used with use_eps_for_freeze_out = 0
    double tempTFO = 0.12;
    int tfoset = 0;
    tempinput = util->StringFind4(input_file, "T_freeze");
    if (tempinput != "empty") {
        istringstream(tempinput) >> tempTFO;
        // if only freeze out temperature is set, freeze out by temperature
        parameter_list->useEpsFO = 0;
        tfoset = 1;
    }
    parameter_list->TFO = tempTFO;
    
    // epsilon_freeze: freeze-out energy density in GeV/fm^3
    // only used with use_eps_for_freeze_out = 1
    double tempepsilonFreeze = 0.12;
    tempinput = util->StringFind4(input_file, "epsilon_freeze");
    if (tempinput != "empty") {
        istringstream(tempinput) >> tempepsilonFreeze;
        // if epsilon_freeze is set, freeze out by epsilon
        parameter_list->useEpsFO = 1;
    }
    parameter_list->epsilonFreeze = tempepsilonFreeze;

    //use_eps_for_freeze_out: 
    // 0: freeze out at constant temperature T_freeze
    // 1: freeze out at constant energy density epsilon_freeze
    // if set in input input_file, overide above defaults
    int tempuseEpsFO = parameter_list->useEpsFO;
    tempinput = util->StringFind4(input_file, "use_eps_for_freeze_out");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempuseEpsFO;
    parameter_list->useEpsFO = tempuseEpsFO;
    if (tfoset == 1 && parameter_list->useEpsFO == 1) {
        music_message << "T_freeze set but overridden"
                      << " -- freezing out by energy density at " 
                      << parameter_list->epsilonFreeze << " GeV/fm^3";
        music_message.flush("warning");
    }
    
    //particle_spectrum_to_compute:
    // 0: Do all up to number_of_particles_to_include
    // any natural number: Do the particle with this (internal) ID
    int tempparticleSpectrumNumber = 0;
    tempinput = util->StringFind4(input_file, "particle_spectrum_to_compute");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempparticleSpectrumNumber;
    parameter_list->particleSpectrumNumber = tempparticleSpectrumNumber;
    
    // mode: 
    // 1: Does everything. Evolution. Computation of thermal spectra.
    //    Resonance decays. Observables.
    // 2: Evolution only.
    // 3: Compute all thermal spectra only.
    // 4: Resonance decays only.
    // 13: Compute observables from previously-computed thermal spectra
    // 14: Compute observables from post-decay spectra
    int tempmode = 1;
    tempinput = util->StringFind4(input_file, "mode");
    if (tempinput != "empty") {
        istringstream(tempinput) >> tempmode;
    } else {
        music_message.error("Must specify mode. Exiting.");
        exit(1);
    }
    parameter_list->mode = tempmode;
    
    //EOS_to_use:
    // 0: ideal gas
    // 1: EOS-Q from azhydro
    // 2: lattice EOS from Huovinen and Petreczky
    // 3: lattice EOS from Huovinen and Petreczky
    //    with partial chemical equilibrium (PCE) at 150 MeV
    //    (see https://wiki.bnl.gov/TECHQM/index.php/QCD_Equation_of_State)
    // 4: PCE EOS with chemical freeze out at 155 MeV
    // 5: PCE EOS at 160 MeV
    // 6: PCE EOS at 165 MeV
    int tempwhichEOS = 2;
    tempinput = util->StringFind4(input_file, "EOS_to_use");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempwhichEOS;
    parameter_list->whichEOS = tempwhichEOS;
    
    // number_of_particles_to_include:
    // This determines up to which particle in the list spectra
    // should be computed (mode=3) or resonances should be included (mode=4)
    // current maximum = 319
    int tempNumberOfParticlesToInclude = 2;
    tempinput = util->StringFind4(input_file,
                                  "number_of_particles_to_include");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempNumberOfParticlesToInclude;
    parameter_list->NumberOfParticlesToInclude = tempNumberOfParticlesToInclude;
    
    // freeze_out_method:
    // 2: Schenke's more complex method
    int tempfreezeOutMethod = 2;
    tempinput = util->StringFind4(input_file, "freeze_out_method");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempfreezeOutMethod;
    parameter_list->freezeOutMethod = tempfreezeOutMethod;

    // rho_b_max: maximum baryon density for zero impact parameter.
    // The shape of the ρB distribution is the same as that
    // for the energy/entropy density distribution
    // any number > 0.
    // (don’t use 0, except if the EOS does not support finite baryon
    // chemical potential, in which case it does not matter what you put)
    double temprhoB0 = 0.000000001;
    tempinput = util->StringFind4(input_file, "rho_b_max");
    if(tempinput != "empty") istringstream ( tempinput ) >> temprhoB0;
    parameter_list->rhoB0 = temprhoB0;
    
    // binary_collision_scaling_fraction: entropy or energy density is proportional to [(1-f)N_{wounded} + f N{binary colls}}
    // for Glauber initialization
    double temphard = 0.05;
    tempinput = util->StringFind4(input_file,
                                  "binary_collision_scaling_fraction");
    if (tempinput != "empty")
        istringstream(tempinput) >> temphard;
    parameter_list->hard = temphard;
    
    // average_surface_over_this_many_time_steps:
    // Only save every N timesteps for finding freeze out surface
    int tempfacTau = 1;
    tempinput = util->StringFind4(input_file,
                                  "average_surface_over_this_many_time_steps");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempfacTau;
    parameter_list->facTau = tempfacTau;
    
    int tempfac_x = 1;
    tempinput = util->StringFind4(input_file, "freeze_Ncell_x_step");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempfac_x;
    parameter_list->fac_x = tempfac_x;
    
    int tempfac_eta = 1;
    tempinput = util->StringFind4(input_file, "freeze_Ncell_eta_step");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempfac_eta;
    parameter_list->fac_eta = tempfac_eta;
    
    // Grid_size_in_*
    // number of cells in x,y direction
    int tempnx = 10;
    tempinput = util->StringFind4(input_file, "Grid_size_in_x");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempnx;
    parameter_list->nx = tempnx;
    int tempny = 10;
    tempinput = util->StringFind4(input_file, "Grid_size_in_y");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempny;
    parameter_list->ny = tempny;
    
    // Grid_size_in_eta
    // number of cells in eta direction.
    // One cell is positioned at eta=0, 
    // half the cells are at negative eta,
    // the rest (one fewer) are at positive eta
    int tempneta = 1;
    tempinput = util->StringFind4(input_file, "Grid_size_in_eta");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempneta;
    parameter_list->neta = tempneta;
    
    // grid_size_in_fm:
    // total length of box in x,y direction in fm (minus delta_*)
    double tempx_size = 25.;
    tempinput = util->StringFind4(input_file, "X_grid_size_in_fm");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempx_size;
    parameter_list->x_size = tempx_size;
    double tempy_size = 25.;
    tempinput = util->StringFind4(input_file, "Y_grid_size_in_fm");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempy_size;
    parameter_list->y_size = tempy_size;
    
    
    // Eta_grid_size:  total length of box in eta direction (minus delta_eta)
    // e.g., neta=8 and eta_size=8 has 8 cells that run from eta=-4 to eta=3
    double tempeta_size = 8.;
    tempinput = util->StringFind4(input_file, "Eta_grid_size");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempeta_size;
    parameter_list->eta_size = tempeta_size;
    
    // Total_evolution_time_tau
    // Maximum evolution time in [fm].
    // Evolution may stop earlier if a freeze-out hypersurface is computed [ parameter "Do_FreezeOut_Yes_1_No_0" ]
    // and all cells are frozen out.
    double temptau_size = 50.;
    tempinput = util->StringFind4(input_file, "Total_evolution_time_tau");
    if (tempinput != "empty")
        istringstream(tempinput) >> temptau_size;
    parameter_list->tau_size = temptau_size;
    
    // Initial_time_tau_0:  in fm
    double temptau0 = 0.4;
    tempinput = util->StringFind4(input_file, "Initial_time_tau_0");
    if (tempinput != "empty")
        istringstream(tempinput) >> temptau0;
    parameter_list->tau0 = temptau0;
    
    if (parameter_list->Initial_profile == 0) {
        /* x-grid, for instance, runs from 0 to nx */
        parameter_list->delta_x = (parameter_list->x_size
                                /static_cast<double>(parameter_list->nx - 1));
        parameter_list->delta_y = (parameter_list->y_size
                                /static_cast<double>(parameter_list->ny - 1));
        parameter_list->delta_eta = (parameter_list->eta_size
                                /static_cast<double>(parameter_list->neta));
    } else {
        /* x-grid, for instance, runs from 0 to nx */
        parameter_list->delta_x = (parameter_list->x_size
                                /static_cast<double>(parameter_list->nx));
        parameter_list->delta_y = (parameter_list->y_size
                                /static_cast<double>(parameter_list->ny));
        parameter_list->delta_eta = (parameter_list->eta_size
                                /static_cast<double>(parameter_list->neta));
    }

    
    music_message << "DeltaX = " << parameter_list->delta_x << " fm";
    music_message.flush("info");
    music_message << "DeltaY = " << parameter_list->delta_y << " fm";
    music_message.flush("info");
    music_message << "DeltaETA = " << parameter_list->delta_eta;
    music_message.flush("info");
    
    // Delta_Tau: 
    // time step to use in [fm].
    double tempdelta_tau = 0.02;
    tempinput = util->StringFind4(input_file, "Delta_Tau");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempdelta_tau;
    parameter_list->delta_tau = tempdelta_tau;
    music_message << "DeltaTau = " << parameter_list->delta_tau << " fm";
    music_message.flush("info");
    
    // output_evolution_data:  
    // 1: output energy, flow velocities, shear stress tensor and bulk pressure information
    // Output at every grid point at every time step by default, but can be controlled
    // with parameters "output_evolution_every_N_timesteps", ... below
    int tempoutputEvolutionData = 0;
    tempinput = util->StringFind4(input_file, "output_evolution_data");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempoutputEvolutionData;
    parameter_list->outputEvolutionData = tempoutputEvolutionData;
    
    int temp_store_hydro_info_in_memory = 0;
    tempinput = util->StringFind4(input_file, "store_hydro_info_in_memory");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_store_hydro_info_in_memory;
    parameter_list->store_hydro_info_in_memory =
                                            temp_store_hydro_info_in_memory;

    // The evolution is outputted every
    // "output_evolution_every_N_timesteps" timesteps
    int temp_evo_N_tau = 1;
    tempinput = util->StringFind4(input_file,
                                  "output_evolution_every_N_timesteps");
    if (tempinput != "empty")
        istringstream(tempinput) >> temp_evo_N_tau;
    parameter_list->output_evolution_every_N_timesteps = temp_evo_N_tau;
    
    int temp_evo_N_x = 1;
    tempinput = util->StringFind4(input_file, "output_evolution_every_N_x");
    if(tempinput != "empty") istringstream ( tempinput ) >> temp_evo_N_x;
    parameter_list->output_evolution_every_N_x = temp_evo_N_x;
    
    int temp_evo_N_y = 1;
    tempinput = util->StringFind4(input_file, "output_evolution_every_N_y");
    if(tempinput != "empty") istringstream ( tempinput ) >> temp_evo_N_y;
    parameter_list->output_evolution_every_N_y = temp_evo_N_y;

    int temp_evo_N_eta = 1;
    tempinput = util->StringFind4(input_file, "output_evolution_every_N_eta");
    if(tempinput != "empty") istringstream ( tempinput ) >> temp_evo_N_eta;
    parameter_list->output_evolution_every_N_eta = temp_evo_N_eta;
    
    // outputBinaryEvolution: output the above information in binary format (float, not double) to reduce diskspace usage
    int tempoutputBinaryEvolution = 0;
    tempinput = util->StringFind4(input_file, "outputBinaryEvolution");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempoutputBinaryEvolution;
    parameter_list->outputBinaryEvolution = tempoutputBinaryEvolution;
   
    
    parameter_list->nt = static_cast<int>(
            parameter_list->tau_size/(parameter_list->delta_tau) + 0.5);
    music_message << __func__ << ": Time step size = "
                  << parameter_list->delta_tau;
    music_message.flush("info");
    music_message << __func__ << ": Number of time steps required = "
                  << parameter_list->nt;
    music_message.flush("info");
  
    // Maximum_energy_density:
    // determines the maximum energy density at zero impact parameter
    // given in [GeV/fm3] (for initialize_with_entropy = 0)
    // or the maximum entropy density at zero impact parameter given in [1/fm3]
    // (for initialize_with_entropy = 1)
    // for Glauber initial conditions
    double tempepsilon0 = 21.67;
    tempinput = util->StringFind4(input_file, "Maximum_energy_density");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempepsilon0;
    parameter_list->epsilon0 = tempepsilon0;
    
    // Eta_fall_off:
    // width of half-Gaussian on each side of a central plateau in eta
    double tempeta_fall_off  = 0.4;
    tempinput = util->StringFind4(input_file, "Eta_fall_off");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempeta_fall_off ;
    parameter_list->eta_fall_off  = tempeta_fall_off;
    
    // Eta_plateau_size:
    // width of the flat region symmetrical around eta=0
    double tempeta_flat = 20.0;
    tempinput = util->StringFind4(input_file, "Eta_plateau_size");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempeta_flat;
    parameter_list->eta_flat = tempeta_flat;

    // s_factor:  normalisation of energy or entropy density at initialisation (see init.cpp for details)
    //   parameter_list->sFactor = util->DFind(input_file, "s_factor");
    double tempsFactor   = 20.;
    tempinput = util->StringFind4(input_file, "s_factor");
    if (tempinput != "empty") istringstream ( tempinput ) >> tempsFactor  ;
    parameter_list->sFactor   = tempsFactor;
    
    // for calculation of spectra:
    // max_pseudorapidity:
    // spectra calculated from pseudorapidity -eta to eta
    double tempmax_pseudorapidity = 5.0;
    tempinput = util->StringFind4(input_file, "max_pseudorapidity");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempmax_pseudorapidity;
    parameter_list->max_pseudorapidity = tempmax_pseudorapidity;
    
    // pseudo_steps:
    // steps in pseudorapidity in calculation of spectra
    int temppseudo_steps = 51;
    tempinput = util->StringFind4(input_file, "pseudo_steps");
    if (tempinput != "empty")
        istringstream(tempinput) >> temppseudo_steps;
    parameter_list->pseudo_steps = temppseudo_steps; 
    
    // phi_steps
    // steps in azimuthal angle in calculation of spectra
    int tempphi_steps = 48;
    tempinput = util->StringFind4(input_file, "phi_steps");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempphi_steps  ;
    parameter_list->phi_steps = tempphi_steps; 
    
    // min_pt: 
    // spectra calculated from this to max_pt transverse momentum in GeV
    double tempmin_pt   = 0.0;
    tempinput = util->StringFind4(input_file, "min_pt");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempmin_pt  ;
    parameter_list->min_pt = tempmin_pt;
      
    // max_pt: 
    // spectra calculated from min_pt to this transverse momentum in GeV
    double tempmax_pt   = 3.0;
    tempinput = util->StringFind4(input_file, "max_pt");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempmax_pt;
    parameter_list->max_pt = tempmax_pt;
    
    // pt_steps:
    // steps in transverse momentum in calculation of spectra
    int temppt_steps   = 60;
    tempinput = util->StringFind4(input_file, "pt_steps");
    if (tempinput != "empty")
        istringstream(tempinput) >> temppt_steps  ;
    parameter_list->pt_steps = temppt_steps;   
    
    // pseudofreeze
    // Calculate spectra at fixed,
    // equally-spaced grid in pseudorapidity, pt, and phi
    int temppseudofreeze = 1;
    // Hard-coded in the public version
//    tempinput = util->StringFind4(input_file, "pseudofreeze");
//    if (tempinput != "empty")
//        istringstream(tempinput) >> temppseudofreeze;
    parameter_list->pseudofreeze = temppseudofreeze;
    
    // Runge_Kutta_order:  must be 1 or 2
    int temprk_order = 2;
    // Hard-coded in the public version
//    tempinput = util->StringFind4(input_file, "Runge_Kutta_order");
//    if (tempinput != "empty")
//        istringstream(tempinput) >> temprk_order;
    parameter_list->rk_order = temprk_order;
    
    // reconstruction type
    // details...
    int tempreconst_type = 1;
    // Hard-coded in the public version
//    tempinput = util->StringFind4(input_file, "reconst_type");
//    if (tempinput != "empty")
//        istringstream (tempinput) >> tempreconst_type;
    parameter_list->reconst_type = tempreconst_type;
    
    // in case of using Initial_Distribution 3 (MC-Glauber initial conditions),
    // these are the limits between which to sample the impact parameter
    double tempbmin   = 0.;
    tempinput = util->StringFind4(input_file, "bmin");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempbmin  ;
    parameter_list->bmin   = tempbmin;
    if (tempinput != "empty" && parameter_list->Initial_profile != 3) {
        music_message.warning("bmin unused when Initial_profile != 3");
    }

    double tempbmax   = 6.73;
    tempinput = util->StringFind4(input_file, "bmax");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempbmax  ;
    parameter_list->bmax = tempbmax;
    if (tempinput != "empty" && parameter_list->Initial_profile != 3) {
        music_message.warning("bmax unused when Initial_profile != 3");
    }

    // Minmod_Theta: theta parameter in the min-mod like limiter
    double tempminmod_theta   = 1.8;
    // Hard-coded in the public version
    //tempinput = util->StringFind4(input_file, "Minmod_Theta");
    //if (tempinput != "empty")
    //    istringstream(tempinput) >> tempminmod_theta  ;
    parameter_list->minmod_theta = tempminmod_theta;
    
    // Viscosity_Flag_Yes_1_No_0:   set to 0 for ideal hydro
    int tempviscosity_flag = 1;
    tempinput = util->StringFind4(input_file, "Viscosity_Flag_Yes_1_No_0");
    if(tempinput != "empty") istringstream ( tempinput ) >> tempviscosity_flag;
    parameter_list->viscosity_flag = tempviscosity_flag;
    
    // Include_Shear_Visc_Yes_1_No_0
    int tempturn_on_shear = 0;
    tempinput = util->StringFind4(input_file, "Include_Shear_Visc_Yes_1_No_0");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempturn_on_shear;
    parameter_list->turn_on_shear = tempturn_on_shear;

    // T_dependent_Shear_to_S_ratio:
    // if 1, ignore constant eta/s
    // and use hard-coded T-dependent shear viscosity (see dissipative.cpp)
    int tempT_dependent_shear_to_s = 0;
    tempinput = util->StringFind4(input_file, "T_dependent_Shear_to_S_ratio");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempT_dependent_shear_to_s;
    parameter_list->T_dependent_shear_to_s = tempT_dependent_shear_to_s;

    //Shear_to_S_ratio:  constant eta/s
    double tempshear_to_s = 0.08;
    tempinput = util->StringFind4(input_file, "Shear_to_S_ratio");
    if (tempinput != "empty") {
        istringstream(tempinput) >> tempshear_to_s;
    } else if (parameter_list->turn_on_shear == 1
                && parameter_list->T_dependent_shear_to_s == 0) {
        music_message.error("please define Shear_to_S_ratio!");
        exit(1);
    }
    parameter_list->shear_to_s = tempshear_to_s;

    // Include_Bulk_Visc_Yes_1_No_0
    int tempturn_on_bulk = 0;
    tempinput = util->StringFind4(input_file, "Include_Bulk_Visc_Yes_1_No_0");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempturn_on_bulk;
    parameter_list->turn_on_bulk = tempturn_on_bulk;
    
    // Include secord order terms
    int tempturn_on_second_order = 0;
    tempinput = util->StringFind4(input_file, "Include_second_order_terms");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempturn_on_second_order;
    parameter_list->include_second_order_terms = tempturn_on_second_order;
    
    // Include_deltaf:  
    // Looks like 0 sets delta_f=0, 1 uses standard quadratic ansatz,
    // and 2 is supposed to use p^(2-alpha)
    int tempinclude_deltaf = 1;
    tempinput = util->StringFind4(input_file, "Include_deltaf");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempinclude_deltaf;
    parameter_list->include_deltaf = tempinclude_deltaf;
    
    int tempinclude_deltaf_bulk = 0;
    tempinput = util->StringFind4(input_file, "Include_deltaf_bulk");
    if(tempinput != "empty") istringstream ( tempinput ) >> tempinclude_deltaf_bulk;
    parameter_list->include_deltaf_bulk = tempinclude_deltaf_bulk;
    
    // Do_FreezeOut_Yes_1_No_0
    // set to 0 to bypass freeze out surface finder
    int tempdoFreezeOut = 1;
    tempinput = util->StringFind4(input_file, "Do_FreezeOut_Yes_1_No_0");
    if(tempinput != "empty") istringstream ( tempinput ) >> tempdoFreezeOut;
    parameter_list->doFreezeOut = tempdoFreezeOut;
    
    // sigma_0:  width of MC-Glauber Gaussian energy deposition
    double tempsigma0   = 0.4;
    tempinput = util->StringFind4(input_file, "sigma_0");
    if(tempinput != "empty") istringstream ( tempinput ) >> tempsigma0  ;
    parameter_list->sigma0   = tempsigma0;
    
    // Initial_Distribution_input_filename
    string tempinitName = "initial/initial_ed.dat";
    tempinput = util->StringFind4(input_file, "Initial_Distribution_input_filename");
    std::string tempinput2 = util->StringFind4(input_file, "Initial_Distribution_Filename");
    if(tempinput != "empty") tempinitName.assign(tempinput);
    //for backward compatibility
    else if (tempinput2 != "empty") tempinitName.assign(tempinput2); 
    parameter_list->initName.assign(tempinitName);
    
    /* initialize the metric, mostly plus */
    parameter_list->gmunu = util->mtx_malloc(4, 4);
    for (m=0; m<4; m++) {
        for (n=0; n<4; n++) {
            if (m == n)
                (parameter_list->gmunu)[m][n] = 1.0;
            else
                (parameter_list->gmunu)[m][n] = 0.0;
            if (m==0 && n==0)
                (parameter_list->gmunu)[m][n] *= -1.0;
        }
    }  /* m */

    //  Make MUSIC output additionnal hydro information
    //  0 for false (do not output), 1 for true
    //bool tempoutput_hydro_debug_info = false;
    //tempinput = util->StringFind4(input_file, "output_hydro_debug_info");
    //if (tempinput != "empty")
    //    istringstream(tempinput) >> tempoutput_hydro_debug_info;
    //parameter_list->output_hydro_debug_info = tempoutput_hydro_debug_info;
    
    //Make MUSIC output a C header input_file containing
    //informations about the hydro parameters used
    //0 for false (do not output), 1 for true
    bool tempoutput_hydro_params_header = false;
    tempinput = util->StringFind4(input_file, "output_hydro_params_header");
    if (tempinput != "empty")
        istringstream(tempinput) >> tempoutput_hydro_params_header;
    parameter_list->output_hydro_params_header = tempoutput_hydro_params_header;

    // initial parameters for modes 13/14
    double temp_dNdy_y_min = -0.5;
    tempinput = util->StringFind4(input_file, "dNdy_y_min");
    if(tempinput != "empty") istringstream ( tempinput ) >> temp_dNdy_y_min;
    parameter_list->dNdy_y_min = temp_dNdy_y_min;
    
    double temp_dNdy_y_max = 0.5;
    tempinput = util->StringFind4(input_file, "dNdy_y_max");
    if(tempinput != "empty") istringstream ( tempinput ) >> temp_dNdy_y_max;
    parameter_list->dNdy_y_max = temp_dNdy_y_max;
    
    double temp_dNdy_eta_min = -2.0;
    tempinput = util->StringFind4(input_file, "dNdy_eta_min");
    if(tempinput != "empty") istringstream ( tempinput ) >> temp_dNdy_eta_min;
    parameter_list->dNdy_eta_min = temp_dNdy_eta_min;
    
    double temp_dNdy_eta_max = 2.0;
    tempinput = util->StringFind4(input_file, "dNdy_eta_max");
    if(tempinput != "empty") istringstream ( tempinput ) >> temp_dNdy_eta_max;
    parameter_list->dNdy_eta_max = temp_dNdy_eta_max;
    
    int temp_dNdy_nrap = 30;
    tempinput = util->StringFind4(input_file, "dNdy_nrap");
    if(tempinput != "empty") istringstream ( tempinput ) >> temp_dNdy_nrap;
    parameter_list->dNdy_nrap = temp_dNdy_nrap;
    
    double temp_dNdyptdpt_y_min = -0.5;
    tempinput = util->StringFind4(input_file, "dNdyptdpt_y_min");
    if(tempinput != "empty") istringstream ( tempinput ) >> temp_dNdyptdpt_y_min;
    parameter_list->dNdyptdpt_y_min = temp_dNdyptdpt_y_min;
    
    double temp_dNdyptdpt_y_max = 0.5;
    tempinput = util->StringFind4(input_file, "dNdyptdpt_y_max");
    if(tempinput != "empty") istringstream ( tempinput ) >> temp_dNdyptdpt_y_max;
    parameter_list->dNdyptdpt_y_max = temp_dNdyptdpt_y_max;
    
    double temp_dNdyptdpt_eta_min = -0.5;
    tempinput = util->StringFind4(input_file, "dNdyptdpt_eta_min");
    if(tempinput != "empty") istringstream ( tempinput ) >> temp_dNdyptdpt_eta_min;
    parameter_list->dNdyptdpt_eta_min = temp_dNdyptdpt_eta_min;
    
    double temp_dNdyptdpt_eta_max = 0.5;
    tempinput = util->StringFind4(input_file, "dNdyptdpt_eta_max");
    if(tempinput != "empty") istringstream ( tempinput ) >> temp_dNdyptdpt_eta_max;
    parameter_list->dNdyptdpt_eta_max = temp_dNdyptdpt_eta_max;

    music_message.info("Done read_in_parameters.");
    check_parameters(parameter_list);
}

void ReadInParameters::check_parameters(InitData *parameter_list) {
    music_message.info("Checking input parameter list ... ");
    if (parameter_list->b < 0) {
        music_message.error("Impact parameter must be greater than zero");
        exit(1);
    }

    if (parameter_list->SigmaNN<0) {
        music_message.error(
                    "NN cross section must be greater than zero (in mb)");
        exit(1);
    }

    if (parameter_list->Initial_profile < 0) {
        music_message << "Initial profile" << parameter_list->Initial_profile
                      << "not defined";
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list->initial_eta_profile > 2
            || parameter_list->initial_eta_profile < 0) {
        music_message << "Initial eta profile"
                      << parameter_list->Initial_profile
                      << "not defined";
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list->initializeEntropy > 1
            || parameter_list->initializeEntropy < 0) {
        music_message.error("Must initialize with entropy or energy");
        exit(1);
    }

    if (parameter_list->TFO < 0.0 || parameter_list->TFO > 0.2) {
        music_message << "T_freeze = " << parameter_list->TFO
                      << " is not physical!";
        music_message.flush("error");
        exit(1);
    }
    
    if (parameter_list->epsilonFreeze <= 0) {
        music_message.error(
                "Freeze out energy density must be greater than zero");
        exit(1);
    }
    
    if (parameter_list->useEpsFO > 1 || parameter_list->useEpsFO < 0) {
        music_message << "Error: did not set either freeze out energy density "
                      << "or temperature, or invalid option for "
                      << "use_eps_for_freeze_out:"
                      << parameter_list->useEpsFO;
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list->whichEOS > 7 || parameter_list->whichEOS < 0) {
        music_message << "EOS_to_use unspecified or invalid option: "
                      << parameter_list->whichEOS;
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list->whichEOS > 1 && parameter_list->whichEOS < 7
            && parameter_list->NumberOfParticlesToInclude > 320) {
        music_message << "Invalid option for number_of_particles_to_include:"
                      << parameter_list->NumberOfParticlesToInclude;
        music_message.flush("error");
        exit(1);
    }
    
    if (parameter_list->whichEOS == 7
            && parameter_list->NumberOfParticlesToInclude > 327) {
        music_message << "Invalid option for number_of_particles_to_include:"
                      << parameter_list->NumberOfParticlesToInclude;
        music_message.flush("error");
        exit(1);
    }
    
    if (parameter_list->freezeOutMethod != 2) {
        music_message << "Invalid option for freeze_out_method: "
                      << parameter_list->freezeOutMethod;
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list->hard < 0. || parameter_list->hard > 1.) {
        music_message << "Glauber WN/BC mixing ratio is unphysical: alpha = "
                      << parameter_list->hard;
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list->facTau <= 0) {
        music_message << "average_surface_over_this_many_time_steps <= 0: "
                      << parameter_list->facTau;
        music_message.flush("error");
        exit(1);
    }

    double freeze_dtau = parameter_list->facTau*parameter_list->delta_tau;
    if (freeze_dtau > 1.) {
        music_message << "freeze-out time setp is too large! "
                      << "freeze_dtau = " << freeze_dtau
                      << ", hydro_dtau = " << parameter_list->delta_tau
                      << ", average_surface_over_this_many_time_steps = "
                      << parameter_list->facTau;
        music_message.flush("warning");
    }
    
    if (parameter_list->fac_x <= 0) {
        music_message << "freeze out every x step <= 0: "
                      << parameter_list->fac_x;
        music_message.flush("error");
        exit(1);
    }
    
    if (parameter_list->fac_eta <= 0) {
        music_message << "freeze out every eta step <= 0: "
                      << parameter_list->fac_eta;
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list->nx != parameter_list->ny) {
        music_message.warning(
                "Grid size in x is not equal to grid size in y!");
    }

    if (parameter_list->neta < 32 && !parameter_list->boost_invariant) {
        music_message << "Grid size in eta = " << parameter_list->neta 
                      << "is too small for a (3+1)-d run! "
                      << "Please increase Grid_size_in_eta to be "
                      << "larger than 32 at least!";
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list->boost_invariant && parameter_list->neta > 1) {
        music_message << "Grid size in eta is set to "
                      << parameter_list->neta << " for a (2+1)-d simulation! "
                      << "This is redundant! Reset neta to 1!";
        music_message.flush("warning");
        parameter_list->neta = 1;
    }
    
    if (parameter_list->delta_tau > 0.1) {
        music_message << "Delta_Tau = " << parameter_list->delta_tau
                      << " maybe too large! "
                      << "Please choose a dtau < 0.1 fm.";
        music_message.flush("warning");

        bool reset_dtau_use_CFL_condition = true;
        int temp_CFL_condition = 1;
        string tempinput = util->StringFind4(input_file,
                                      "reset_dtau_use_CFL_condition");
        if (tempinput != "empty")
            istringstream(tempinput) >> temp_CFL_condition;
        if (temp_CFL_condition == 0)
            reset_dtau_use_CFL_condition = false;

        if (reset_dtau_use_CFL_condition) {
            music_message.info("reset dtau using CFL condition.");
            double dtau_CFL = mini(
                    parameter_list->delta_x/10.0,
                    parameter_list->tau0*parameter_list->delta_eta/10.0);
            parameter_list->delta_tau = dtau_CFL;
            parameter_list->nt = static_cast<int>(
                parameter_list->tau_size/(parameter_list->delta_tau) + 0.5);
            music_message << __PRETTY_FUNCTION__ << "Time step size = "
                          << parameter_list->delta_tau;
            music_message.flush("info");
            music_message << __PRETTY_FUNCTION__
                          << "Number of time steps required = "
                          << parameter_list->nt;
            music_message.flush("info");
        } else {
            exit(1);
        }
    }

    if (parameter_list->min_pt > parameter_list->max_pt) {
        music_message << "min_pt = " << parameter_list->min_pt << " > "
                      << "max_pt = " << parameter_list->max_pt;
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list->phi_steps < 20) {
        music_message << "phi_steps = " << parameter_list->phi_steps
                      << " is too small for computing vn, "
                      << "please increase to at least 40!";
        music_message.flush("error");
        exit(1);
    }
    
    if (parameter_list->rk_order > 2 || parameter_list->rk_order < 0) {
        music_message << "Invalid option for Runge_Kutta_order: "
                      << parameter_list->rk_order;
        music_message.flush("error");
        exit(1);
    }
    if (parameter_list->rk_order != 2) {
        music_message << "Runge-Kutta order = " << parameter_list->rk_order;
        music_message.flush("info");
    }
    
    if (parameter_list->reconst_type != 1
            && parameter_list->reconst_type != 0) {
        music_message << "unrecognized reconst_type: "
                      << parameter_list->reconst_type;
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list->bmin > parameter_list->bmax) {
        music_message << "bmin = " << parameter_list->bmin << " > "
                      << "bmax = " << parameter_list->bmax;
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list->minmod_theta < 1.
            || parameter_list->minmod_theta > 2.) {
        music_message << "minmod = " << parameter_list->minmod_theta
                      << " is out of allowed range [1., 2]";
        music_message.flush("error");
        exit(1);
    }

    if (parameter_list->turn_on_shear == 0 && parameter_list->shear_to_s > 0) {
        music_message << "non-zero eta/s = " << parameter_list->shear_to_s
                      << " is set with "
                      << "Include_Shear_Visc = "
                      << parameter_list->turn_on_shear
                      << ". Please check you want to run ideal hydro!";
        music_message.flush("warning");
    }

    if (parameter_list->turn_on_shear == 0
            && parameter_list->include_deltaf == 1) {
        music_message << "hydro with zero shear viscosity does not need "
                      << "shear delta f in Cooper-Frye! ";
        music_message << "input Include_deltaf = "
                      << parameter_list->include_deltaf
                      << ". Now rewrite it to 0!";
        music_message.flush("warning");
        parameter_list->include_deltaf = 0;
    }

    if (parameter_list->turn_on_bulk == 0
            && parameter_list->include_deltaf_bulk == 1) {
        music_message << "hydro with zero bulk viscosity does not need "
                      << "bulk delta f in Cooper-Frye!";
        music_message << "input Include_deltaf_bulk = "
                      << parameter_list->include_deltaf_bulk
                      << ". Now rewrite it to 0!";
        music_message.flush("warning");
        parameter_list->include_deltaf_bulk = 0;
    }

    if (parameter_list->turn_on_bulk == 1
            && parameter_list->include_deltaf_bulk == 0 && ((parameter_list->mode == 1)||(parameter_list->mode==3)) ) {
        music_message << "hydro with bulk viscosity should also include "
                      << "a bulk delta f in Cooper-Frye!";
        music_message << "input Include_deltaf_bulk = "
                      << parameter_list->include_deltaf_bulk
                      << ". Unless you're doing this on purpose, go back to your input file and set \"Include_deltaf_bulk  1\"!";
        music_message.flush("warning");
    }

    if (parameter_list->output_evolution_every_N_timesteps <= 0) {
        music_message.error("output_evolution_every_N_timesteps < 0!");
        exit(1);
    }
    
    if (parameter_list->output_evolution_every_N_x <= 0) {
        music_message.error("output_evolution_every_N_x < 0!");
        exit(1);
    }
    
    if (parameter_list->output_evolution_every_N_y <= 0) {
        music_message.error("output_evolution_every_N_y < 0!");
        exit(1);
    }
    
    if (parameter_list->output_evolution_every_N_eta <= 0) {
        music_message.error("output_evolution_every_N_eta < 0!");
        exit(1);
    }

    if (parameter_list->dNdy_y_min > parameter_list->dNdy_y_max) {
        music_message << "dNdy_y_min = " << parameter_list->dNdy_y_min << " < " 
                      << "dNdy_y_max = " << parameter_list->dNdy_y_max << "!";
        music_message.flush("error");
        exit(1);
    }
    
    if (parameter_list->dNdy_eta_min > parameter_list->dNdy_eta_max) {
        music_message << "dNdy_eta_min = " << parameter_list->dNdy_eta_min
                      << " < " 
                      << "dNdy_eta_max = " << parameter_list->dNdy_eta_max
                      << "!";
        music_message.flush("error");
        exit(1);
    }
    
    if (parameter_list->dNdyptdpt_y_min > parameter_list->dNdyptdpt_y_max) {
        music_message << "dNdyptdpt_y_min = "
                      << parameter_list->dNdyptdpt_y_min << " < " 
                      << "dNdyptdpt_y_max = "
                      << parameter_list->dNdyptdpt_y_max << "!";
        music_message.flush("error");
        exit(1);
    }
    
    if (parameter_list->dNdyptdpt_eta_min
            > parameter_list->dNdyptdpt_eta_max) {
        music_message << "dNdyptdpt_eta_min = "
                      << parameter_list->dNdyptdpt_eta_min << " < " 
                      << "dNdyptdpt_eta_max = "
                      << parameter_list->dNdyptdpt_eta_max << "!";
        music_message.flush("error");
        exit(1);
    }

    music_message << "Finished checking input parameter list. "
                  << "Everything looks reasonable so far "
                  << emoji_face->success();
    music_message.flush("info");
}
