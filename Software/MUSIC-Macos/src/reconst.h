// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

#ifndef RECONST_H
#define RECONST_H

#include "./data.h"
#include "./grid.h"
#include "./eos.h"
#include "./pretty_ostream.h"
#include <iostream>
#include <cstring>

//! This class handles reconstruct hydrodynamical variables from T^\mu\nu tensor
class Reconst {
 private:
    EOS *eos;
    InitData *DATA_ptr;
    pretty_ostream music_message;

    int echo_level;
    int max_iter;
    double abs_err;
    double rel_err;
    double eos_eps_max;

 public:
    Reconst(EOS *eos, InitData *DATA_in);//constructor
    ~Reconst();//destructor
    
    double GuessEps(double T00, double K00, double cs2);

    //! This is a shell function for reconstruct T^\mu\nu from T^\tau\mu
    void ReconstIt_shell(Grid *grid_p, int direc, double tau, double **uq,
                         Grid *grid_pt, int rk_flag);
    
    //! reconstruct TJb from q[0] - q[4]
    //! use iteration to solve  eps and rhob
    int ReconstIt(Grid *grid_p, int direc, double tau,
                  double **uq, Grid *grid_pt, int rk_flag);

    //! reconstruct TJb from q[0] - q[4]
    //! reconstruct velocity first for finite mu_B case
    //! use iteration to solve v and u0
    int ReconstIt_velocity_iteration(Grid *grid_p, int direc, double tau,
                                     double **uq, Grid *grid_pt,
                                     int rk_flag);

    //! this function returns f(v) = M/(M0 + P)
    double reconst_velocity_f(double v, double T00, double M, double J0);

    //! this function returns f(u0) = (M0+P)/sqrt((M0+P)^2 - M^2)
    double reconst_u0_f(double u0, double T00, double K00, double M,
                        double J0);

    //! This function reverts the grid information back its values
    //! at the previous time step
    void revert_grid(Grid *grid_current, Grid *grid_prev, int rk_flag);

    //! This function regulate the grid information
    void regulate_grid(Grid *grid_cell, double elocal, int rk_flag);
};
#endif
