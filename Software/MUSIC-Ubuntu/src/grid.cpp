// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

#include "util.h"
#include "grid.h"
#include <iostream>

using namespace std;

Grid::Grid() {
    prev_epsilon = 0.0;
    epsilon_t = 0.0;
    epsilon = 0.0;
    epsilon_prev = 0.0;

    prev_rhob = 0.0;
    rhob_t = 0.0;
    rhob = 0.0;
    rhob_prev = 0.0;

    p = 0.0;
    p_t = 0.0;

    pi_b_prev = 0.0;

    T = 0.0;
    mu = 0.0;

    TJb = NULL;
    u = NULL;
    prev_u = NULL;

    nbr_p_1 = NULL;
    nbr_m_1 = NULL;
    nbr_p_2 = NULL;
    nbr_m_2 = NULL;
    
    a = NULL;
    theta_u = NULL;
    sigma = NULL;
    dUsup = NULL;

    Wmunu = NULL;
    prevWmunu = NULL;
    W_prev = NULL;

    Pimunu = NULL;
    prevPimunu = NULL;
    pi_b = NULL;


    for (int i = 0; i < 4; i++) {
        position[i] = 0;
        u_prev[i] = 0.0;
    }
}

Grid *Grid::grid_v_malloc(int n1) {
    Grid *d1_ptr;
  
    /* pointer to the n1 array */
    d1_ptr = new Grid[n1];
    
    return d1_ptr;
}/* grid_v_malloc */


Grid **Grid::grid_m_malloc(int n1, int n2) {
    Grid **d1_ptr, *tmp_ptr;

    tmp_ptr = (Grid *)malloc(sizeof(Grid)*n1*n2);
    d1_ptr = (Grid **) malloc (sizeof(Grid *)*n1);

    for(int i = 0; i < n1; i++) {
       d1_ptr[i] = &(tmp_ptr[i*n2]);
    }
     
    return d1_ptr;
}/* grid_m_malloc */

Grid ***Grid::grid_c_malloc(int n1, int n2, int n3) {
   Grid ***d1_ptr;

   d1_ptr = new Grid **[n1];

   for(int i = 0; i < n1; i++) {
      d1_ptr[i] = new Grid *[n2];
      for (int j = 0; j < n2; j++) {
         d1_ptr[i][j] = new Grid[n3];
      }
   } 
   
   return d1_ptr;
}/* grid_c_malloc */

void Grid::print_all_information() {
    cout << "Grid cell position: " << position[1] << ", " << position[2]
         << ", " << position[3] << endl;
    cout << "e = " << epsilon << ", p = " << p << ", rhob = " << rhob << endl;
    cout << "e_t = " << epsilon_t << ", p_t = " << p_t
         << ", rhob_t = " << rhob_t << endl;
    cout << "prev_e = " << prev_epsilon << ", prev_rhob = " << prev_rhob
         << endl;
    cout << "T = " << T << ", mu = " << mu << endl;
    cout << "Quantities for freeze-out: " << endl;
    cout << "e_prev = " << epsilon_prev << ", rhob_prev = " << rhob_prev
         << endl;
    cout << "u_prev^mu = " << u_prev[0] << ", " << u_prev[1] << ", "
         << u_prev[2] << ", " << u_prev[3] << endl;
    cout << "W_prev[" << 0 << "] = " << endl;
    for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
            cout << W_prev[j][k] << ", ";
        }
        cout << endl;
    }
         
    int i = 0;
    cout << "prev_u^mu[" << i << "]= " << prev_u[i][0] << ", "
         << prev_u[i][1] << ", " << prev_u[i][2] << ", " << prev_u[i][3]
         << endl;
    cout << "prevWmunu[" << i << "] = " << endl;
    for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
            cout << prevWmunu[i][j][k] << ", ";
        }
        cout << endl;
    }
    cout << "prevPimunu[" << i << "] = " << endl;
    for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
            cout << prevPimunu[i][j][k] << ", ";
        }
        cout << endl;
    }

    for (i = 0; i < 3; i++) {
        cout << "rk_flag = " << i << endl;
        cout << "u^mu[" << i << "] = " << u[i][0] << ", " << u[i][1] << ", "
             << u[i][2] << ", " << u[i][3] << endl;
        cout << "a^mu[" << i << "] = " << a[i][0] << ", " << a[i][1] << ", "
             << a[i][2] << ", " << a[i][3] << endl;
        cout << "theta[" << i << "] = " << theta_u[i] << endl;
        cout << "sigma[" << i << "] = " << endl;
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                cout << sigma[i][j][k] << ", ";
            }
            cout << endl;
        }
        cout << "dUsup[" << i << "] = " << endl;
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                cout << dUsup[i][j][k] << ", ";
            }
            cout << endl;
        }
        cout << "TJb[" << i << "] = " << endl;
        for (int j = 0; j < 5; j++) {
            for (int k = 0; k < 4; k++) {
                cout << TJb[i][j][k] << ", ";
            }
            cout << endl;
        }
        cout << "Wmunu[" << i << "] = " << endl;
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                cout << Wmunu[i][j][k] << ", ";
            }
            cout << endl;
        }
        cout << "Pimunu[" << i << "] = " << endl;
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                cout << Pimunu[i][j][k] << ", ";
            }
            cout << endl;
        }
    }
}
