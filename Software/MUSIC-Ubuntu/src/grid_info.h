// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Maxime Dion, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-Francois Paquet, Bjoern Schenke, Chun Shen, Clint Young

#ifndef GRID_INFO_H
#define GRID_INFO_H
#include "data.h"
#include "eos.h"
#include "grid.h"
#include "HydroinfoMUSIC.h"
#include <mpi.h>
#include <iostream>
#include <iomanip>

//! This class contains many helper functions for outputing useful information
//! about hydrodynamic evolution
class Grid_info
{
 private:
    InitData* DATA_ptr;
        
 public:
    Grid_info(InitData* DATA_in);
    ~Grid_info();
    
    void OutputEvolutionDataXYEta(Grid ***arena, InitData *DATA, EOS *eos,
                                  double tau, int size, int rank,
                                  HydroinfoMUSIC *hydro_info_ptr);
    void Output_hydro_information_header(InitData *DATA, EOS *eos);

    //! This function outputs quantities for Gubser flow check
    void Gubser_flow_check_file(Grid ***arena, EOS *eos_ptr, double tau);

    //! This function outputs the evolution of hydrodynamic variables at a
    //! give fluid cell
    void monitor_fluid_cell(Grid ***arena, int ix, int iy, int ieta,
                            double tau);
};
#endif
