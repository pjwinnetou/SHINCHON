// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

#ifndef GRID_H
#define GRID_H
#include "data.h"
#include "eos.h"
#include <mpi.h>
#include <iostream>
#include <iomanip>

//! This class contains all hydrodynamic variables on one grid point

class Grid
{
 public:
    double epsilon;
    double p;
    double rhob;
    //! stress energy tensor plus baryon current
    /*! TJb[flag][alpha][mu] */
    /*! flag = 0 is the actual values. flag != 0 are the intermediate values
      ! for the Runge-Kutta step */
    /*! alpha = 4 means Jb */
    double ***TJb; 
    
    //! temporary values for the final RK update
    double epsilon_t;
    //! temporary values for the final RK update
    double p_t;
    //! temporary values for the final RK update
    double rhob_t;
    
    //! store the epsilon and rhob at previous time step
    double prev_epsilon;
    //! store the epsilon and rhob at previous time step
    double prev_rhob;

    //! u[flag][mu]: flag=0 is the actual values. flag !=0 are for RK steps
    double **u;
    
    /* to include shear viscosity */
    /* we need to calculate partial_tau u[mu] */
    //! u[mu] from the previous time step including the rk flag
    double **prev_u;
    
    Grid **nbr_p_1; 
    Grid **nbr_m_1; 
    Grid **nbr_p_2; 
    Grid **nbr_m_2; 
    
    //! This is u^mu partial_mu u^nu
    double **a;
    
    //! this is the expansion rate partial_mu u^mu
    double *theta_u;

    //! the velocity shear tensor
    double ***sigma;
    
    /* we need to calculate partial_tau u[mu] */
    //! dU[flag][m][n] = u^{m,n} = partial^n u^m with the rk flag
    // note that they are superscripted. So partial^t = -partial_t
    double ***dUsup; 
    //! shear part of the TJb with the rk_flag
    double ***Wmunu;
    //! shear part of the TJb at previous time step
    double ***prevWmunu; 
    
    //! bulk part of the TJb with the rk_flag
    double ***Pimunu;
    //! bulk part of the TJb at previous time step
    double ***prevPimunu; 
    
    //! bulk pressure */
    double *pi_b;

    double T;
    double mu;
    
    int position[4];

    //! the following variables are for hyper-surface finder to determine freeze-out surface
    //! they are only updated every freeze-out step not every evolution time step
    double epsilon_prev;
    double rhob_prev;
    double u_prev[4];
    double pi_b_prev;
    //! the one for the freeze-out surface finder for interpolation
    double **W_prev;
    
    Grid();
    ~Grid() {};
    Grid *grid_v_malloc(int );
    Grid **grid_m_malloc(int , int );
    Grid ***grid_c_malloc(int , int , int );
    void print_all_information();
};
#endif
