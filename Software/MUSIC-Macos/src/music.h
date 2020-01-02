// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

// Copyright

#ifndef SRC_MUSIC_H_
#define SRC_MUSIC_H_

#include <cstring>
#include <vector>

#include "./util.h"
#include "./grid.h"
#include "./data.h"
#include "./init.h"
#include "./eos.h"
#include "./freeze.h"
#include "./evolve.h"
#include "./read_in_parameters.h"
#include "./pretty_ostream.h"
#include "./HydroinfoMUSIC.h"

#include "mpi.h"
#include <stdio.h>
#include <sys/stat.h>  // for mkdir


using namespace std;

//! This is a wrapper class for MUSIC hydro

class MUSIC {
 private:
    int mode;            //!< records running mode

    int rank;
    int size;

    Grid ***arena;
    Grid ***Lneighbor;
    Grid ***Rneighbor;
    
    string input_file;
    InitData DATA;
    
    ReadInParameters reader;

    Init *init;
    Glauber *glauber;
    EOS *eos;
    Evolve *evolve;
    Freeze *freeze;
    HydroinfoMUSIC *hydro_info_ptr;

    pretty_ostream music_message;

 public:
    MUSIC(int argc, char *argv[]);
    ~MUSIC();

    int get_running_mode () {return(mode);}

    int initialize_hydro();  //!< Initialize hydro
    int initialize_hydro_from_vector(std::vector<double> entropy_density,
                                     double dx);
    int initialize_hydro_from_pre_equilibrium_vectors(
        const double dx, const double dz, const double z_max, const int nz,
        std::vector<double> e_in,
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

    int run_hydro();         //!< run hydrodynamic simulations

    int run_Cooper_Frye(int running_mode);  //!< perform Cooper-Frye freeze-out

    //! This function prints out welcome message
    void welcome_message();

    //! This function displays the program logo
    void display_logo(int selector);

    //! This function prints out code desciprtion and copyright information
    void display_code_description_and_copyright();

    //! This function returns hydrodynamic fields informaiton at a given
    //! space-time point (t, x, y, z)
    void get_hydro_info(double x, double y, double z, double t,
                        fluidCell* info);

    //! This function cleans up the memory which stores hydro evolution history
    void clear_hydro_info_from_memory();
};

#endif  // SRC_MUSIC_H_
