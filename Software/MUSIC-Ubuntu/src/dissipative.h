// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

#ifndef DISSIPATIVE_H
#define DISSIPATIVE_H

#include "util.h"
#include "grid.h"
#include "data.h"
#include "minmod.h"
#include "./pretty_ostream.h"
#include <iostream>

//! This class computes the all viscous source terms for T^\mu\nu evolution and for viscous stress tensor

class Diss{
 private:
  EOS *eos;
  Minmod *minmod;
  pretty_ostream music_message;

 public:
  Diss(EOS *eosIn, InitData* DATA_in);
  ~Diss();
  
  double MakeWSource(double tau, int alpha, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
		     Grid *Lneighbor2, Grid *Rneighbor2, InitData *DATA, int rk_flag, int size, int rank);
  
  int Make_uWRHS
    (double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
     Grid *Lneighbor2, Grid *Rneighbor2, double **w_rhs, InitData *DATA, int rk_flag, int size, int rank);
  
  void Get_uWmns(double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
			 Grid *Lneighbor2, Grid *Rneighbor2, 
		 int mu, int nu, int direc,
		 double *g, double *f, double *gp1, double *fp1,
		 double *gp2, double *fp2, double *gm1, double *fm1, 
		 double *gm2, double *fm2, InitData *DATA, int rk_flag, int size, int rank);
  
  double Make_uWSource
    (double tau, Grid *grid_pt, int mu, int nu, InitData *DATA, int rk_flag); 
  
  int Make_uPRHS
    (double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
     Grid *Lneighbor2, Grid *Rneighbor2, 
     double *p_rhs, InitData *DATA, int rk_flag, int size, int rank);
  void Get_uPis(double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
			 Grid *Lneighbor2, Grid *Rneighbor2, int direc,
		double *g, double *f, double *gp1, double *fp1,
		double *gp2, double *fp2, double *gm1, double *fm1, 
		double *gm2, double *fm2, InitData *DATA, int rk_flag, int size, int rank); 
  
  double Make_uPiSource
    (double tau, Grid *grid_pt, InitData *DATA, int rk_flag);

};
#endif
