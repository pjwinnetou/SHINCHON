// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

#ifndef U_DERIVATIVE_H
#define U_DERIVATIVE_H

#include "data.h"
#include "grid.h"
#include <string.h>
#include <iostream>

//! This class contains functions to compute derivatives of flow velocity
class U_derivative{

 private:
  Minmod *minmod;
  
 public:
  // Sangyong Nov 18 2014: added EOS *eos in the argument
  U_derivative(InitData* DATA_in);//constructor
  ~U_derivative();
  int UpdateDSpatial(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
		     Grid *Lneighbor2, Grid *Rneighbor2, int rk_flag, int size, int rank);
  int MakedU(double tau, InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int rk_flag, int size, int rank);
  int MakeDSpatial(double tau, InitData *DATA, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
		   Grid *Lneighbor2, Grid *Rneighbor2, int rk_flag, int size, int rank);
  int MakeDTau(double tau, InitData *DATA, Grid *grid_pt, int rk_flag);
  int UpdateDTau(double tau, InitData *DATA, Grid *grid_pt, int rk_flag);
  int UpdatedU(double tau, InitData *DATA, Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, int rk_flag, int size, int rank);
};
#endif
