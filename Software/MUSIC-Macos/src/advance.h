// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

// Copyright

#ifndef SRC_ADVANCE_H_
#define SRC_ADVANCE_H_

#include "data.h"
#include "grid.h"
#include "reconst.h"
#include "dissipative.h"
#include "minmod.h"
#include "./pretty_ostream.h"
#include <iostream>

//! This class performs propagation for hydrodynamic variables using the KT algorithm.

//! advance routines separate for 
//! T^{0 nu} \del T^{i\nu} (T)
//! W
//! T^{0 nu} with W source (TS)
//! W with source (WS)

class Advance{
 private:
  
    Util *util;
    Reconst *reconst; // declare Reconst object
    Diss *diss; // dissipative object
    EOS *eos;
    Grid *grid;
    Minmod *minmod;
    pretty_ostream music_message;

    int grid_nx, grid_ny, grid_neta;
    int rk_order;
  
    typedef struct bdry_cells {
        Grid *grid_p_h_L;
        Grid *grid_p_h_R;
        Grid *grid_m_h_L;
        Grid *grid_m_h_R;
    
        double **qiphL;
        double **qiphR;
        double **qimhL;
        double **qimhR;
    
    } BdryCells;
  
  
    typedef struct nbrs {
        double **qip1;
        double **qip2;
        double **qim1;
        double **qim2;
    } NbrQs;

 public:
    Advance(EOS *eosIn, Grid *grid, InitData* DATA_in);
    ~Advance();

    int AdvanceIt(double tau_init, InitData *DATA, Grid ***arena,
                  Grid ***Lneighbor, Grid ***Rneighbor, int rk_flag,
                  int size, int rank);

   void MPISendReceiveW(Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, 
                        int size, int rank, int rk_flag);
   void MPISendReceiveT(Grid ***arena, Grid ***Lneighbor, Grid ***Rneighbor, 
                        int size, int rank, int rk_flag);

   int AdvanceLocalT(double tau_init, InitData *DATA, Grid *grid_pt, 
                     Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2, 
		         int rk_flag, int size, int rank);
   int AdvanceLocalW(double tau_init, InitData *DATA, Grid *grid_pt, 
                     Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2, 
		         int rk_flag, int size, int rank);
 
   int FirstRKStepT(double tau_it, InitData *DATA, Grid *grid_pt, 
                    Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2, 
                    int rk_flag, double *qi, double *rhs,
                    double **qirk, Grid *grid_rk, int size, int rank);
  
   int FirstRKStepW(double tau_it, InitData *DATA, Grid *grid_pt,
                    Grid *Lneighbor, Grid *Rneighbor,
                    Grid *Lneighbor2, Grid *Rneighbor2,
                    int rk_flag, int size, int rank);

   void UpdateTJbRK(Grid *grid_rk, Grid *grid_pt, int rk_flag);
   int QuestRevert(double tau, Grid *grid_pt, int rk_flag, InitData *DATA, int size, int rank);
   void TestW(double tau, Grid *grid_pt, int rk_flag);
   void ProjectSpin2W(double tau, Grid *grid_pt, int rk_flag, InitData *DATA, int size, int rank);
   void ProjectSpin2WS(double tau, Grid *grid_pt, int rk_flag, InitData *DATA, int size, int rank);
  
   void MakeDeltaQI(double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, 
		   Grid *Lneighbor2, Grid *Rneighbor2, double *qi, double *rhs, 
		   InitData *DATA, int rk_flag, int size, int rank);
   void GetQIs(double tau, Grid *grid_pt, Grid *Lneighbor, Grid *Rneighbor, Grid *Lneighbor2, Grid *Rneighbor2,
	      double *qi, NbrQs *NbrCells, int rk_flag, InitData *DATA, int size, int rank);
   int MakeQIHalfs(double *qi, NbrQs *NbrCells, BdryCells *HalfwayCells, 
		  Grid *grid_pt, InitData *DATA);
   int ConstHalfwayCells(double tau, BdryCells *HalfwayCells,
                         double *qi, Grid *grid_pt, int rk_flag);
   void MakeKTCurrents(double tau, double **DFmmp, Grid *grid_pt, 
		      BdryCells *HalfwayCells, int rk_flag);
   void MakeMaxSpeedAs(double tau, BdryCells *HalfwayCells, double aiph[], double aimh[], int rk_flag);
   double MaxSpeed (double tau, int direc, Grid *grid_p, int rk_flag);
   void InitNbrQs(NbrQs *NbrCells);
   void InitTempGrids(BdryCells *HalfwayCells, int rk_order);

};  
#endif  // SRC_ADVANCE_H_
