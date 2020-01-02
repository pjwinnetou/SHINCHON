// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

// Copyright @ 2011 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#ifndef SRC_EVOLVE_H_
#define SRC_EVOLVE_H_

#include <mpi.h>
#include <time.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>

#include "./util.h"
#include "./data.h"
#include "./grid.h"
#include "./grid_info.h"
#include "./eos.h"
#include "./advance.h"
#include "./u_derivative.h"
#include "./pretty_ostream.h"
#include "./HydroinfoMUSIC.h"

//! this is a control class for the hydrodynamic evolution

class Evolve {
 private:
    EOS *eos; // declare EOS object
    Grid *grid; // declare Grid object
    Grid_info *grid_info;
    Util *util;
    Advance *advance;
    U_derivative *u_derivative;
    pretty_ostream music_message;

    InitData *DATA_ptr;
    bool boost_invariant;

    // simulation information
    int rk_order;
    int grid_nx, grid_ny, grid_neta;
    
    double SUM, SUM2;
    int warnings;
    int cells;
    int weirdCases;
    int facTau;
  
 public:
    Evolve(EOS *eos, InitData *DATA_in);//constructor
    ~Evolve();//destructor
    int EvolveIt(InitData *DATA, Grid ***arena,
                 Grid ***Lneighbor, Grid ***Rneighbor, int size, int rank,
                 HydroinfoMUSIC *hydro_info_ptr);
    
    int AdvanceRK(double tau, InitData *DATA, Grid ***arena, 
                  Grid ***Lneighbor, Grid ***Rneighbor, int size, int rank);
    
    int UpdateArena(double tau, Grid ***arena);
    
    int FindFreezeOutSurface2(double tau, InitData *DATA, Grid ***arena, 
                               int size, int rank);

    void initial_prev_variables(Grid ***arena);
    void storePreviousEpsilon(Grid ***arena);
    void storePreviousW(Grid ***arena);
};

#endif  // SRC_EVOLVE_H_
