// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

// Copyright @ 2011 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#include "./util.h"
#include "./grid.h"
#include "./init.h"
#include "./eos.h"
#include <vector>


using namespace std;

Init::Init(EOS *eosIn, Glauber *glauberIn) {
    eos = eosIn;
    util = new Util;
    glauber = glauberIn;
    random = gsl_rng_alloc(gsl_rng_ranlxs2);
}

// destructor
Init::~Init() {
    gsl_rng_free(random);
    delete util;
}

void Init::InitArena(InitData *DATA, Grid ****arena, Grid ****Lneighbor, 
                     Grid ****Rneighbor, int size, int rank) {
    Grid *helperGrid;
    helperGrid = new Grid;
    music_message.info("initArena");
    
    if (DATA->Initial_profile == 0) {
        music_message << "Using Initial_profile=" << DATA->Initial_profile;
        music_message.flush("info");
        DATA->nx = DATA->nx - 1;
        DATA->ny = DATA->ny - 1;
        DATA->delta_x = 0.05;
        DATA->delta_y = 0.05;
        music_message << "nx=" << DATA->nx+1 << ", ny=" << DATA->ny+1;
        music_message.flush("info");
        music_message << "dx=" << DATA->delta_x << ", dy=" << DATA->delta_y;
        music_message.flush("info");
    } else if (DATA->Initial_profile == 8) {
        music_message.info(DATA->initName);
        ifstream profile(DATA->initName.c_str());
        string dummy;
        int nx, ny, neta;
        double deta, dx, dy, dummy2;
        
        // read the first line with general info
        profile >> dummy >> dummy >> dummy2 
                >> dummy >> neta >> dummy >> nx >> dummy >> ny 
                >> dummy >> deta >> dummy >> dx >> dummy >> dy ;
        profile.close();
        music_message << "Using Initial_profile=" << DATA->Initial_profile 
                      << ". Overwriting lattice dimensions:";
        music_message.flush("info");

        DATA->nx = nx-1;
        DATA->ny = ny-1;
        DATA->delta_x = dx;
        DATA->delta_y = dy;

        music_message << "neta=" << neta << ", nx=" << nx << ", ny=" << ny;
        music_message.flush("info");
        music_message << "deta=" << DATA->delta_eta << ", dx=" << DATA->delta_x 
                      << ", dy=" << DATA->delta_y;
        music_message.flush("info");
    } else if (DATA->Initial_profile == 41 || DATA->Initial_profile == 42) {
        // initial condition from the JETSCAPE framework
        music_message << "Using Initial_profile=" << DATA->Initial_profile 
                      << ". Overwriting lattice dimensions:";
        music_message.flush("info");

        int nx = static_cast<int>(sqrt(initial_entropy_density.size()));
        if (DATA->Initial_profile == 42) {
            nx = static_cast<int>(sqrt(initial_energy_density.size()/DATA->neta));
        }
        int ny = nx;
        DATA->nx = nx - 1;
        DATA->ny = ny - 1;
        DATA->x_size = DATA->delta_x*(nx - 1);
        DATA->y_size = DATA->delta_y*(ny - 1);

        music_message << "neta = " << DATA->neta
                      << ", nx = " << nx << ", ny = " << ny;
        music_message.flush("info");
        music_message << "deta=" << DATA->delta_eta
                      << ", dx=" << DATA->delta_x 
                      << ", dy=" << DATA->delta_y;
        music_message.flush("info");
        music_message << "x_size = " << DATA->x_size
                      << ", y_size = " << DATA->y_size
                      << ", eta_size = " << DATA->eta_size;
        music_message.flush("info");
    }

    // initialize arena
    *arena = helperGrid->grid_c_malloc(DATA->nx+1, DATA->ny+1, DATA->neta);
    *Lneighbor = helperGrid->grid_c_malloc(DATA->nx+1, DATA->ny+1, 2);
    *Rneighbor = helperGrid->grid_c_malloc(DATA->nx+1, DATA->ny+1, 2);
    
    music_message.info("Grid allocated.");
    InitTJb(DATA, arena, Lneighbor, Rneighbor, size, rank);
    
    music_message << "rank " << rank << " ok.";
    music_message.flush("info");
    LinkNeighbors(DATA, arena, size, rank);
    delete helperGrid;
}


void Init::get_entropy_density_vector(vector<double> entropy_density_in) {
    initial_entropy_density = entropy_density_in;
}


void Init::get_pre_equilibrium_vectors(std::vector<double> e_in,
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
                                       std::vector<double> Bulk_pi_in) {
    initial_energy_density = e_in;
    initial_u_tau          = u_tau_in;
    initial_u_x            = u_x_in;
    initial_u_y            = u_y_in;
    initial_u_eta          = u_eta_in;
    initial_pi_00          = pi_00_in;
    initial_pi_01          = pi_01_in;
    initial_pi_02          = pi_02_in;
    initial_pi_03          = pi_03_in;
    initial_pi_11          = pi_11_in;
    initial_pi_12          = pi_12_in;
    initial_pi_13          = pi_13_in;
    initial_pi_22          = pi_22_in;
    initial_pi_23          = pi_23_in;
    initial_pi_33          = pi_33_in;
    initial_bulk_pi        = Bulk_pi_in;
}


void Init::LinkNeighbors(InitData *DATA, Grid ****arena, int size, int rank) {
    int nx = DATA->nx;
    int ny = DATA->ny;
    int neta = DATA->neta;

    /* allocate memory */
    for (int ix = 0; ix <= nx; ix++) {
        for (int iy = 0; iy <= ny; iy++) {
            for (int ieta = 0; ieta < neta; ieta++) {
                (*arena)[ix][iy][ieta].nbr_p_1 = new Grid *[4];
                (*arena)[ix][iy][ieta].nbr_m_1 = new Grid *[4];
                (*arena)[ix][iy][ieta].nbr_p_2 = new Grid *[4];
                (*arena)[ix][iy][ieta].nbr_m_2 = new Grid *[4];
    
                (*arena)[ix][iy][ieta].position[1] = ix;
                (*arena)[ix][iy][ieta].position[2] = iy;
                (*arena)[ix][iy][ieta].position[3] = ieta;
            }  /* ieta */
        }  /* iy */
    }  /* ix */

    for (int ix = 0; ix <= nx; ix++) {
        for (int iy = 0; iy <= ny; iy++) {
            for (int ieta = 0; ieta < neta; ieta++) {
                if (ix != nx)
                    (*arena)[ix][iy][ieta].nbr_p_1[1] = 
                                                    &(*arena)[ix+1][iy][ieta]; 
                else 
                    (*arena)[ix][iy][ieta].nbr_p_1[1] =
                                                    &(*arena)[nx][iy][ieta]; 

                if(ix < nx - 1)
                    (*arena)[ix][iy][ieta].nbr_p_2[1] = 
                                                    &(*arena)[ix+2][iy][ieta]; 
                else 
                    (*arena)[ix][iy][ieta].nbr_p_2[1] =
                                                    &(*arena)[nx][iy][ieta]; 
         
                if(ix != 0)
                    (*arena)[ix][iy][ieta].nbr_m_1[1] = 
                                                    &(*arena)[ix-1][iy][ieta]; 
                else
                    (*arena)[ix][iy][ieta].nbr_m_1[1] =
                                                    &(*arena)[0][iy][ieta]; 

                if(ix > 1)
                    (*arena)[ix][iy][ieta].nbr_m_2[1] = 
                                                    &(*arena)[ix-2][iy][ieta]; 
                else
                    (*arena)[ix][iy][ieta].nbr_m_2[1] =
                                                    &(*arena)[0][iy][ieta]; 
         
                if(iy != ny)
                    (*arena)[ix][iy][ieta].nbr_p_1[2] = 
                                                    &(*arena)[ix][iy+1][ieta]; 
                else
                    (*arena)[ix][iy][ieta].nbr_p_1[2] =
                                                    &(*arena)[ix][ny][ieta]; 
        
                if(iy < ny - 1)
                    (*arena)[ix][iy][ieta].nbr_p_2[2] = 
                                                    &(*arena)[ix][iy+2][ieta]; 
                else
                    (*arena)[ix][iy][ieta].nbr_p_2[2] =
                                                    &(*arena)[ix][ny][ieta]; 
         
                if(iy != 0)
                    (*arena)[ix][iy][ieta].nbr_m_1[2] = 
                                                    &(*arena)[ix][iy-1][ieta]; 
                else
                    (*arena)[ix][iy][ieta].nbr_m_1[2] =
                                                    &(*arena)[ix][0][ieta]; 
        
                if(iy > 1)
                    (*arena)[ix][iy][ieta].nbr_m_2[2] = 
                                                    &(*arena)[ix][iy-2][ieta]; 
                else
                    (*arena)[ix][iy][ieta].nbr_m_2[2] =
                                                    &(*arena)[ix][0][ieta]; 

                // do not care which rank it is 
                // - that is dealt with in evolve.cpp
                if(ieta != neta-1)
                    (*arena)[ix][iy][ieta].nbr_p_1[3] = 
                                                    &(*arena)[ix][iy][ieta+1]; 
                else 
                    (*arena)[ix][iy][ieta].nbr_p_1[3] =
                                                    &(*arena)[ix][iy][neta-1]; 
        
                if(ieta < neta-2)
                    (*arena)[ix][iy][ieta].nbr_p_2[3] = 
                                                    &(*arena)[ix][iy][ieta+2]; 
                else 
                    (*arena)[ix][iy][ieta].nbr_p_2[3] =
                                                    &(*arena)[ix][iy][neta-1]; 
         
                if(ieta != 0)
                    (*arena)[ix][iy][ieta].nbr_m_1[3] = 
                                                    &(*arena)[ix][iy][ieta-1]; 
                else
                    (*arena)[ix][iy][ieta].nbr_m_1[3] =
                                                    &(*arena)[ix][iy][0]; 
        
                if(ieta > 1)
                    (*arena)[ix][iy][ieta].nbr_m_2[3] = 
                                                    &(*arena)[ix][iy][ieta-2]; 
                else
                    (*arena)[ix][iy][ieta].nbr_m_2[3] =
                                                    &(*arena)[ix][iy][0]; 
     
            }/* ieta */
        }/* iy */
    }/* ix */
}/* LinkNeighbors */


void Init::sampleTA()
{
    int A = static_cast<int>(glauber->nucleusA());

    MCGlauberReturnValue rv, rv2;

    for (int i = 0; i < A; i++) // get all nucleon coordinates
    {
        rv = glauber->SampleTARejection(random);
        rv2 = glauber->SampleTARejection(random);
        nucleusA.push_back(rv);
        nucleusB.push_back(rv2);
    }
}  

int Init::InitTJb(InitData *DATA, Grid ****arena, Grid ****Lneighbor, 
                  Grid ****Rneighbor, int size, int rank) {
    double epsilon = 0.0;
    double rhob = 0.0;
    double p, u[4], x, y, eta;
    int ix, iy, ieta, mu, nu, rk_order;
    int initializeEntropy = DATA->initializeEntropy;

    rk_order = DATA->rk_order;
    music_message << "rk_order=" << rk_order;
    music_message.flush("info");
 
    if (DATA->Initial_profile == 0) {
        // Gubser flow test
        music_message.info(" Perform Gubser flow test ... ");
        music_message.info(" ----- information on initial distribution -----");

        for (ix = 0; ix <= DATA->nx; ix++) {
            for (iy = 0; iy <= DATA->ny; iy++) {
                (*Lneighbor)[ix][iy][0].TJb =
                                    util->cube_malloc(rk_order+1, 5, 4);
                (*Rneighbor)[ix][iy][0].TJb =
                                    util->cube_malloc(rk_order+1, 5, 4);
                (*Lneighbor)[ix][iy][1].TJb =
                                    util->cube_malloc(rk_order+1, 5, 4);
                (*Rneighbor)[ix][iy][1].TJb =
                                    util->cube_malloc(rk_order+1, 5, 4);
                (*Lneighbor)[ix][iy][0].Wmunu =
                                    util->cube_malloc(rk_order+1, 4, 4);
                (*Rneighbor)[ix][iy][0].Wmunu =
                                    util->cube_malloc(rk_order+1, 4, 4);
                (*Lneighbor)[ix][iy][1].Wmunu =
                                    util->cube_malloc(rk_order+1, 4, 4);
                (*Rneighbor)[ix][iy][1].Wmunu =
                                    util->cube_malloc(rk_order+1, 4, 4);
                (*Lneighbor)[ix][iy][0].Pimunu =
                                    util->cube_malloc(rk_order+1, 4, 4);
                (*Rneighbor)[ix][iy][0].Pimunu =
                                    util->cube_malloc(rk_order+1, 4, 4);
                (*Lneighbor)[ix][iy][1].Pimunu =
                                    util->cube_malloc(rk_order+1, 4, 4);
                (*Rneighbor)[ix][iy][1].Pimunu =
                                    util->cube_malloc(rk_order+1, 4, 4);
                (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
                (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
                (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
                (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
                (*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
                (*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
                (*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
                (*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
            }
        }

        string input_filename;
        string input_filename_prev;
        if (DATA->turn_on_shear == 1) {
            input_filename = "tests/Gubser_flow/Initial_Profile.dat";
        } else {
            input_filename = "tests/Gubser_flow/y=0_tau=1.00_ideal.dat";
            input_filename_prev = "tests/Gubser_flow/y=0_tau=0.98_ideal.dat";
        }
        ifstream profile(input_filename.c_str());
        if (!profile.good()) {
            music_message << __PRETTY_FUNCTION__
                          << ": Can not open the initial file: "
                          << input_filename;
            music_message.flush("error");
            exit(1);
        }
        ifstream profile_prev;
        if (DATA->turn_on_shear == 0) {
            profile_prev.open(input_filename_prev.c_str());
            if (!profile_prev.good()) {
                music_message << __PRETTY_FUNCTION__
                              << ": Can not open the initial file: "
                              << input_filename_prev;
                music_message.flush("error");
                exit(1);
            }
        }
        int nx = DATA->nx + 1;
        int ny = DATA->ny + 1;
        double** temp_profile_ed = new double* [nx];
        double** temp_profile_ux = new double* [nx];
        double** temp_profile_uy = new double* [nx];
        double **temp_profile_ed_prev = NULL;
        double **temp_profile_rhob = NULL;
        double **temp_profile_rhob_prev = NULL;
        double **temp_profile_ux_prev = NULL;
        double **temp_profile_uy_prev = NULL;
        double **temp_profile_pixx = NULL;
        double **temp_profile_piyy = NULL;
        double **temp_profile_pixy = NULL;
        double **temp_profile_pi00 = NULL;
        double **temp_profile_pi0x = NULL;
        double **temp_profile_pi0y = NULL;
        double **temp_profile_pi33 = NULL;
        if (DATA->turn_on_shear == 1) {
            temp_profile_pixx = new double* [nx];
            temp_profile_piyy = new double* [nx];
            temp_profile_pixy = new double* [nx];
            temp_profile_pi00 = new double* [nx];
            temp_profile_pi0x = new double* [nx];
            temp_profile_pi0y = new double* [nx];
            temp_profile_pi33 = new double* [nx];
        } else {
            temp_profile_ed_prev = new double* [nx];
            temp_profile_rhob = new double* [nx];
            temp_profile_rhob_prev = new double* [nx];
            temp_profile_ux_prev = new double* [nx];
            temp_profile_uy_prev = new double* [nx];
        }
        for (int i = 0; i < nx; i++) {
            temp_profile_ed[i] = new double[ny];
            temp_profile_ux[i] = new double[ny];
            temp_profile_uy[i] = new double[ny];
            if (DATA->turn_on_shear == 1) {
                temp_profile_pixx[i] = new double[ny];
                temp_profile_pixy[i] = new double[ny];
                temp_profile_piyy[i] = new double[ny];
                temp_profile_pi00[i] = new double[ny];
                temp_profile_pi0x[i] = new double[ny];
                temp_profile_pi0y[i] = new double[ny];
                temp_profile_pi33[i] = new double[ny];
            } else {
                temp_profile_ed_prev[i] = new double[ny];
                temp_profile_rhob[i] = new double[ny];
                temp_profile_rhob_prev[i] = new double[ny];
                temp_profile_ux_prev[i] = new double[ny];
                temp_profile_uy_prev[i] = new double[ny];
            }
        }

        double dummy;
        double u[4];
        int rk_order = DATA->rk_order;
        for (ix = 0; ix < nx; ix++) {
            for (iy = 0; iy < ny; iy++) {
                if (DATA->turn_on_shear == 1) {
                    profile >> dummy >> dummy >> temp_profile_ed[ix][iy]
                            >> temp_profile_ux[ix][iy]
                            >> temp_profile_uy[ix][iy];
                    profile >> temp_profile_pixx[ix][iy]
                            >> temp_profile_piyy[ix][iy]
                            >> temp_profile_pixy[ix][iy]
                            >> temp_profile_pi00[ix][iy]
                            >> temp_profile_pi0x[ix][iy]
                            >> temp_profile_pi0y[ix][iy]
                            >> temp_profile_pi33[ix][iy];
                } else {
                    profile >> dummy >> dummy >> temp_profile_ed[ix][iy]
                            >> temp_profile_rhob[ix][iy]
                            >> temp_profile_ux[ix][iy]
                            >> temp_profile_uy[ix][iy];
                    profile_prev >> dummy >> dummy
                                 >> temp_profile_ed_prev[ix][iy]
                                 >> temp_profile_rhob_prev[ix][iy]
                                 >> temp_profile_ux_prev[ix][iy]
                                 >> temp_profile_uy_prev[ix][iy];
                }
            }
        }
        profile.close();
        if (DATA->turn_on_shear == 0) {
            profile_prev.close();
        }

        for (ieta = 0; ieta < DATA->neta; ieta++) {
            for (ix = 0; ix <= DATA->nx; ix++) {
                for (iy = 0; iy<= DATA->ny; iy++) {
                    double rhob = 0.0;
                    double epsilon = temp_profile_ed[ix][iy];
                    // initial pressure distribution
                    double p = eos->get_pressure(epsilon, rhob);
                    (*arena)[ix][iy][ieta].epsilon = epsilon;
                    (*arena)[ix][iy][ieta].epsilon_t = epsilon;
                    (*arena)[ix][iy][ieta].prev_epsilon = epsilon;
                    (*arena)[ix][iy][ieta].rhob = rhob;
                    (*arena)[ix][iy][ieta].rhob_t = rhob;
                    (*arena)[ix][iy][ieta].prev_rhob = rhob;
                    (*arena)[ix][iy][ieta].p = p;
                    (*arena)[ix][iy][ieta].p_t = p;
                    (*arena)[ix][iy][ieta].T =
                                        eos->get_temperature(epsilon, rhob);
                    (*arena)[ix][iy][ieta].mu = eos->get_mu(epsilon, rhob);
                    (*arena)[ix][iy][ieta].TJb =
                                        util->cube_malloc(rk_order+1, 5, 4);
                    (*arena)[ix][iy][ieta].dUsup =
                                        util->cube_malloc(rk_order+1, 4, 4);
                    (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
                    (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 4);
                    (*arena)[ix][iy][ieta].theta_u =
                                        util->vector_malloc(rk_order+1);
                    (*arena)[ix][iy][ieta].sigma =
                                        util->cube_malloc(rk_order+1, 4, 4);
                    (*arena)[ix][iy][ieta].pi_b =
                                        util->vector_malloc(rk_order+1);
                    (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
                    (*arena)[ix][iy][ieta].Wmunu =
                                        util->cube_malloc(rk_order+1, 4, 4);
                    (*arena)[ix][iy][ieta].prevWmunu =
                                        util->cube_malloc(1, 4, 4);
                    (*arena)[ix][iy][ieta].Pimunu =
                                        util->cube_malloc(rk_order+1, 4, 4);
                    (*arena)[ix][iy][ieta].prevPimunu =
                                        util->cube_malloc(1, 4, 4);
                    (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(4, 4);
                    double utau_local = sqrt(1.
                            + temp_profile_ux[ix][iy]*temp_profile_ux[ix][iy]
                            + temp_profile_uy[ix][iy]*temp_profile_uy[ix][iy]);
                    (*arena)[ix][iy][ieta].u[0][0] = utau_local;
                    (*arena)[ix][iy][ieta].u[0][1] = temp_profile_ux[ix][iy];
                    (*arena)[ix][iy][ieta].u[0][2] = temp_profile_uy[ix][iy];
                    (*arena)[ix][iy][ieta].u[0][3] = 0.0;
                    
                    u[0] = utau_local;
                    u[1] = temp_profile_ux[ix][iy];
                    u[2] = temp_profile_uy[ix][iy];
                    u[3] = 0.0;
                    
                    if (DATA->turn_on_shear == 0) {
                        double utau_prev = sqrt(1.
                                + temp_profile_ux_prev[ix][iy]
                                  *temp_profile_ux_prev[ix][iy]
                                + temp_profile_uy_prev[ix][iy]
                                  *temp_profile_uy_prev[ix][iy]);
                        (*arena)[ix][iy][ieta].prev_u[0][0] = utau_prev;
                        (*arena)[ix][iy][ieta].prev_u[0][1] =
                                                temp_profile_ux_prev[ix][iy];
                        (*arena)[ix][iy][ieta].prev_u[0][2] =
                                                temp_profile_uy_prev[ix][iy];
                        (*arena)[ix][iy][ieta].prev_u[0][3] = 0.0;
                    } else {
                        (*arena)[ix][iy][ieta].prev_u[0][0] = utau_local;
                        (*arena)[ix][iy][ieta].prev_u[0][1] =
                                                    temp_profile_ux[ix][iy];
                        (*arena)[ix][iy][ieta].prev_u[0][2] =
                                                    temp_profile_uy[ix][iy];
                        (*arena)[ix][iy][ieta].prev_u[0][3] = 0.0;
                    }

                    (*arena)[ix][iy][ieta].pi_b[0] = 0.0;
                    
                    if (DATA->turn_on_shear == 1) {
                        (*arena)[ix][iy][ieta].Wmunu[0][0][0] =
                                                    temp_profile_pi00[ix][iy];
                        (*arena)[ix][iy][ieta].Wmunu[0][0][1] =
                                                    temp_profile_pi0x[ix][iy];
                        (*arena)[ix][iy][ieta].Wmunu[0][0][2] =
                                                    temp_profile_pi0y[ix][iy];
                        (*arena)[ix][iy][ieta].Wmunu[0][0][3] = 0.0;
                        (*arena)[ix][iy][ieta].Wmunu[0][1][0] =
                                                    temp_profile_pi0x[ix][iy];
                        (*arena)[ix][iy][ieta].Wmunu[0][1][1] =
                                                    temp_profile_pixx[ix][iy];
                        (*arena)[ix][iy][ieta].Wmunu[0][1][2] =
                                                    temp_profile_pixy[ix][iy];
                        (*arena)[ix][iy][ieta].Wmunu[0][1][3] = 0.0;
                        (*arena)[ix][iy][ieta].Wmunu[0][1][3] = 0.0;
                        (*arena)[ix][iy][ieta].Wmunu[0][2][0] =
                                                    temp_profile_pi0y[ix][iy];
                        (*arena)[ix][iy][ieta].Wmunu[0][2][1] =
                                                    temp_profile_pixy[ix][iy];
                        (*arena)[ix][iy][ieta].Wmunu[0][2][2] =
                                                    temp_profile_piyy[ix][iy];
                        (*arena)[ix][iy][ieta].Wmunu[0][2][3] = 0.0;
                        (*arena)[ix][iy][ieta].Wmunu[0][3][0] = 0.0;
                        (*arena)[ix][iy][ieta].Wmunu[0][3][1] = 0.0;
                        (*arena)[ix][iy][ieta].Wmunu[0][3][2] = 0.0;
                        (*arena)[ix][iy][ieta].Wmunu[0][3][3] =
                                                    temp_profile_pi33[ix][iy];
                        for(int mu = 0; mu < 4; mu++) {
                            for (nu = 0; nu < 4; nu++) {
                                (*arena)[ix][iy][ieta].prevWmunu[0][nu][mu] =
                                    (*arena)[ix][iy][ieta].Wmunu[0][nu][mu];
                            }
                        }
                    } else {
                        for(int mu = 0; mu < 4; mu++) {
                            for (nu = 0; nu < 4; nu++) {
                                (*arena)[ix][iy][ieta].Wmunu[0][nu][mu] = 0.0;
                                (*arena)[ix][iy][ieta].prevWmunu[0][nu][mu] = 0.0;
                            }
                        }
                    }
                    for (int mu = 0; mu < 4; mu++) {
                        /* baryon density */
                        (*arena)[ix][iy][ieta].TJb[0][4][mu] = rhob*u[mu];
                        for (nu = 0; nu < 4; nu++) {
                            (*arena)[ix][iy][ieta].TJb[0][nu][mu] = (
                                (epsilon + p)*u[mu]*u[nu] 
                                + p*(DATA->gmunu)[mu][nu]);
                            (*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = 0.0;
                            (*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = 0.0;
                        }
                    }
                }
            }
        }
        // clean up
        for (int i = 0; i < nx; i++) {
            delete[] temp_profile_ed[i];
            delete[] temp_profile_ux[i];
            delete[] temp_profile_uy[i];
            if (DATA->turn_on_shear == 1) {
                delete[] temp_profile_pixx[i];
                delete[] temp_profile_piyy[i];
                delete[] temp_profile_pixy[i];
                delete[] temp_profile_pi00[i];
                delete[] temp_profile_pi0x[i];
                delete[] temp_profile_pi0y[i];
                delete[] temp_profile_pi33[i];
            } else {
                delete[] temp_profile_ed_prev[i];
                delete[] temp_profile_rhob[i];
                delete[] temp_profile_rhob_prev[i];
                delete[] temp_profile_ux_prev[i];
                delete[] temp_profile_uy_prev[i];
            }
        }
        delete[] temp_profile_ed;
        delete[] temp_profile_ux;
        delete[] temp_profile_uy;
        if (DATA->turn_on_shear == 1) {
            delete[] temp_profile_pixx;
            delete[] temp_profile_piyy;
            delete[] temp_profile_pixy;
            delete[] temp_profile_pi00;
            delete[] temp_profile_pi0x;
            delete[] temp_profile_pi0y;
            delete[] temp_profile_pi33;
        } else {
            delete[] temp_profile_ed_prev;
            delete[] temp_profile_rhob;
            delete[] temp_profile_rhob_prev;
            delete[] temp_profile_ux_prev;
            delete[] temp_profile_uy_prev;
        }
    } else if (DATA->Initial_profile == 1) {
        double epsilon0 = DATA->epsilon0;  // this is in GeV/fm^3
        if (initializeEntropy == 0)
            epsilon0 /= hbarc;  // everything is in fm now. eps is in fm^-4
        // full average initial conditions using Glauber
        double s, r1, r2, W, nBinary, nWounded;
        // impact parameter:
        double b=DATA->b;

      double normT = glauber->LexusData.Target.A;
      double normP = glauber->LexusData.Projectile.A;
      double hard = DATA->hard;
      
      // The value of W at x=y=0 for normalization of W in central collisions:
      double W0 = (
           hard*(TATarget(DATA, 0.)*TAProjectile(DATA, 0.)*DATA->SigmaNN/10.)
         + (1. - hard)*(  TATarget(DATA, 0.)
                          *(1. - pow((1. - (DATA->SigmaNN/10.)
                                            *TAProjectile(DATA, 0.)/normP), 
                                     normP)
                           )
                    + TAProjectile(DATA, 0.)
                          *(1. - pow((1. - (DATA->SigmaNN/10.)
                                            *TATarget(DATA, 0.)/normT), normT)
                           )
                       )
                 ); 
     
     // loop over the whole lattice and initialize values:
     for(ix=0; ix<=DATA->nx; ix++)
       {
     x = DATA->delta_x*(ix*2 - DATA->nx)/2.0;
     for(iy=0; iy<=DATA->ny; iy++)
       {
         y = DATA->delta_y*(iy*2 - DATA->ny)/2.0;
         r1 = sqrt((x+b/2.)*(x+b/2.)+(y*y));
         r2 = sqrt((x-b/2.)*(x-b/2.)+(y*y));
         //number of binary collisions:
         nBinary = TAProjectile(DATA, r2)*TATarget(DATA, r1)*DATA->SigmaNN/10.;
         //number of wounded nucleons:
         nWounded = TATarget(DATA, r1)*(1.-pow((1.-(DATA->SigmaNN/10.)*TAProjectile(DATA, r2)/normP),normP))
           + TAProjectile(DATA, r2)*(1.-pow((1.-(DATA->SigmaNN/10.)*TATarget(DATA, r1)/normT),normT));
         //distribution in the transverse plane, normalized so that maximum value is 1:
         W = ((1.-hard)*nWounded + hard*nBinary)/W0;

         (*Lneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
         (*Rneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
         (*Lneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
         (*Rneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
         (*Lneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 4,4);
         (*Rneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 4,4);
         (*Lneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 4,4);
         (*Rneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 4,4);
         (*Lneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 4,4);
         (*Rneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 4,4);
         (*Lneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 4,4);
         (*Rneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 4,4);
         (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
         (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
         (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
         (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
         (*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
         (*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
         (*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
         (*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
         int fakeEnvelope=0;
         if(fakeEnvelope==1)
           {
         if (ix==0 && iy==0) {
             music_message.warning("USING THE FAKE TESTING ENVELOPE!!!!!!");
         }
         W*=(1.6+1.6*tanh((sqrt(0.001*pow(y,4.))))-1.6*tanh(sqrt(0.2*pow(x,4))));
           }
          for(ieta=0; ieta<DATA->neta; ieta++)
           {
         eta = (DATA->delta_eta)*(ieta+DATA->neta*rank) - (DATA->eta_size)/2.0;

         if ( initializeEntropy==0 )
           {
             // distribution of the initial energy density in eta:
             epsilon = epsilon0*eta_profile_normalisation(DATA, eta);
             // and x,y:
             epsilon *= W;

             if (epsilon<0.00000000001)
               epsilon = 0.00000000001;
             
             // distribution of initial baryon density:
             //rhob = DATA->rhoB0/epsilon0*epsilon; /* in fm^-3 */

             rhob = 0.0; /* in fm^-3 */
           }
         else if ( initializeEntropy==1 )
           {
             s = epsilon0*W;
             s *= eta_profile_normalisation(DATA, eta);
             rhob = DATA->rhoB0/epsilon0*s;
             if ( DATA->whichEOS==1 && s>96.288 )
               {
             epsilon=(pow((4./169./pow(PI,2.)),1./3.)/30.*pow(22.5,(4./3.))*
                  pow(s,(4./3.))*0.197+0.0)/hbarc;
               }
             else if (s>0.2)
               {
             epsilon = eos->findRoot(&EOS::ssolve, rhob, s, 1.15*rhob+0.001, 300.,0.001);
               }
             else 
               {
             epsilon=0;
               }
           }
         // intial pressure distribution:
         p = eos->get_pressure(epsilon, rhob);

         // set all values in the grid element:
         (*arena)[ix][iy][ieta].epsilon = epsilon;
         (*arena)[ix][iy][ieta].rhob = rhob;
         (*arena)[ix][iy][ieta].p = p;

         (*arena)[ix][iy][ieta].T = eos->get_temperature(epsilon, rhob); 
         (*arena)[ix][iy][ieta].mu = eos->get_mu(epsilon, rhob);
         
         (*arena)[ix][iy][ieta].TJb = util->cube_malloc(rk_order+1, 5,4);
         (*arena)[ix][iy][ieta].dUsup = util->cube_malloc(rk_order+1, 4,4);
         (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
         (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 4);
         (*arena)[ix][iy][ieta].theta_u = util->vector_malloc(rk_order+1);
         (*arena)[ix][iy][ieta].sigma = util->cube_malloc(rk_order+1, 4, 4);
         (*arena)[ix][iy][ieta].pi_b = util->vector_malloc(rk_order+1);
         (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
         (*arena)[ix][iy][ieta].Wmunu = util->cube_malloc(rk_order+1, 4,4);
         (*arena)[ix][iy][ieta].prevWmunu = util->cube_malloc(1, 4,4);
         (*arena)[ix][iy][ieta].Pimunu = util->cube_malloc(rk_order+1, 4,4);
         (*arena)[ix][iy][ieta].prevPimunu = util->cube_malloc(1, 4,4);

         (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(4,4);
         
         /* for HIC */
         u[0] = (*arena)[ix][iy][ieta].u[0][0] = 1.0;
         u[3] = (*arena)[ix][iy][ieta].u[0][3] = 0.0;
         u[1] = (*arena)[ix][iy][ieta].u[0][1] = 0.0;
         u[2] = (*arena)[ix][iy][ieta].u[0][2] = 0.0;

         (*arena)[ix][iy][ieta].prev_u[0][0] = 1.0;
         (*arena)[ix][iy][ieta].prev_u[0][3] = 0.0;
         (*arena)[ix][iy][ieta].prev_u[0][1] = 0.0;
         (*arena)[ix][iy][ieta].prev_u[0][2] = 0.0;
         
         (*arena)[ix][iy][ieta].pi_b[0] = 0.0;

         for(mu=0; mu<4; mu++)
           {
             /* baryon density */
             (*arena)[ix][iy][ieta].TJb[0][4][mu] = rhob*u[mu];
             
             for(nu=0; nu<4; nu++)
               {
             (*arena)[ix][iy][ieta].TJb[0][nu][mu] 
               = (epsilon + p)*u[mu]*u[nu] + p*(DATA->gmunu)[mu][nu];

             (*arena)[ix][iy][ieta].Wmunu[0][nu][mu] = (double) 0.0;
             (*arena)[ix][iy][ieta].prevWmunu[0][nu][mu] = (double) 0.0;
             
             (*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = (double) 0.0;
             (*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = (double) 0.0;
               }/* nu */
           }/* mu */
         
           }}}/* ix, iy, ieta */
    } else if (DATA->Initial_profile == 3) {
        // sample initial distribution with a Glauber Monte-Carlo
        // (event-by-event initial condition)
        // keep this in mind: from Nara and Hirano: 0904.2966
        // In our MC-Glauber model, we find the defaultWoods-Saxon
        // distribution is reproduced by a larger radius parameter
        // r0 = 6.42 (4.28) fm and a smaller diffuseness parameter
        // d = 0.44 (0.50) fm (i.e., sharper boundary of a nucleus)
        // for a gold (copper) nucleus than the default parameters.
         
        // sample Nu (which gives the number of nucleons at radial position r)
        // for both nuclei for now use the same nucleus
        
        // initialize random seed using the system time
        double rnum = time(0) + DATA->seed;
        music_message << "random seed=" << rnum-DATA->seed << " + "
                      << DATA->seed << " = " << rnum;
        music_message.flush("info");
        gsl_rng_set(random,static_cast <unsigned long int> (rnum));

        // width of Gaussian in fm around the center
        // of the nucleon-nucleon collision
        double width = DATA->sigma0;
        
        // arrays for sending when using MPI
        double bArray[1];
        double AxArray[300], AyArray[300];
        double BxArray[300], ByArray[300];
        MCGlauberReturnValue rv, rv2; //auxiliary nucleons
        double x, y, xm, ym, dx, dy, dij;
        int i;
        double b = 0.;

        // minimal and maximal impact parameter in fm
        double bmin = DATA->bmin;
        double bmax = DATA->bmax;
   

        // set up number of nucleons, free arrays with nucleons
        // in nucleus A and B
        double A = glauber->nucleusA();
        if (A > 300) {
            music_message.error("A > 300, change the code to make it work");
            exit(1);
        }
        nucleusA.clear();
        nucleusB.clear();


        // sample the distribution of nucleons when on rank 0,
        // send the result to the other ranks
        if (rank == 0) {
            // let only one rank do the sampling
            // and then pass the same distribution to all other ranks
            // first sample a b from the distribution
            // b db = 2b/(b_max^2-b_min^2)
            // uniformly distributed random variable
            double xb =  gsl_rng_uniform(random);
            b = sqrt((bmax*bmax - bmin*bmin)*xb + bmin*bmin);
            // that's all, we have b.
            // find a testing routine at the end of this file
            // populate the lists nucleusA and nucleusB with position data
            // of the nucleons
            sampleTA();
            for (i = 0; i < A; i++) {
                bArray[0] = b;
                AxArray[i] = nucleusA.at(i).x;
                AyArray[i] = nucleusA.at(i).y;
                BxArray[i] = nucleusB.at(i).x;
                ByArray[i] = nucleusB.at(i).y;
            }

            for ( i = 1; i < size; i++) {
                // send to all other ranks
                MPI_Send(bArray, 1, MPI_DOUBLE, i, 5, MPI_COMM_WORLD);
                MPI_Send(AxArray, 300, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
                MPI_Send(AyArray, 300, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
                MPI_Send(BxArray, 300, MPI_DOUBLE, i, 3, MPI_COMM_WORLD);
                MPI_Send(ByArray, 300, MPI_DOUBLE, i, 4, MPI_COMM_WORLD);
            }
        }

        if (rank != 0) {
            // let only one rank do the sampling
            // and then pass the same distribution to all other ranks
            MPI_Recv(bArray, 1, MPI_DOUBLE, 0, 5,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(AxArray, 300, MPI_DOUBLE, 0, 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(AyArray, 300, MPI_DOUBLE, 0, 2,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(BxArray, 300, MPI_DOUBLE, 0, 3,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(ByArray, 300, MPI_DOUBLE, 0, 4,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     
            for (int i = 0; i < A; i++) {
                b = bArray[0];
                rv.x= AxArray[i];
                rv.y = AyArray[i];
                rv2.x = BxArray[i];
                rv2.y = ByArray[i];

                nucleusA.push_back(rv);
                nucleusB.push_back(rv2);
            }
        }
        music_message << "rank " << rank << ": Using impact parameter b=" << b;
        music_message.flush("info");

        // define the radius^2 in fm^2:
        double d2 = DATA->SigmaNN/(PI*10.);          // in fm^2
         
        // determine whether the nucleons collide,
        // depending on their position/overlap:
        for (int i = 0; i < A; i++) {
            // shift the nuclei's position by =/- b/2 or +b/2 respectively
            nucleusB.at(i).x = nucleusB.at(i).x + b/2.;
            nucleusA.at(i).x = nucleusA.at(i).x - b/2.;
        }   
     
        int ibinColl = 0;       // counter for binary collisions
        int NbinColl = 0;       // total number of binary collisions
        double xbinColl[5000];  // x coordinate of binary collisions
        double ybinColl[5000];  // y coordinate
         
        for (int i = 0; i < A; i++) {
            for (int j = 0 ; j < A ;j++) {
                dx = nucleusB.at(j).x - nucleusA.at(i).x;
                dy = nucleusB.at(j).y - nucleusA.at(i).y;
                dij = dx*dx + dy*dy;
                if (dij < d2) {
                    nucleusB.at(j).collided = 1;
                    nucleusA.at(i).collided = 1;
                    xbinColl[ibinColl] =
                                    (nucleusA.at(i).x + nucleusB.at(j).x)/2.;
                    ybinColl[ibinColl] =
                                    (nucleusA.at(i).y + nucleusB.at(j).y)/2.;
                    ibinColl+=1;
                }
            }
        }

        NbinColl = ibinColl;
        music_message << "number of binary collisions = " << NbinColl;
        music_message.flush("info");
     
        int participants = 0;
        if (rank == 0) {
            for (int i = 0; i<A; i++) {
                if (nucleusA.at(i).collided == 1) {
                    participants ++;
                }
                if (nucleusB.at(i).collided == 1) {
                    participants ++;
                }
            }
        }
        music_message << " N_part = " << participants;
        music_message.flush("info");

        double s;
        // impact parameter:
        // normalization of TATarget and TAProjectile is given by
        // Target.A and Projectie.A, respectively
        // the fraction of the hard (binary collisions) contribution
        double hard = DATA->hard;

        // loop over the whole lattice and initialize values:
        for (ix = 0; ix <= DATA->nx; ix++) {
            x = DATA->delta_x*(ix*2 - DATA->nx)/2.0;
            for(iy = 0; iy <= DATA->ny; iy++) {
                y = DATA->delta_y*(iy*2 - DATA->ny)/2.0;
             
                double W = 0.;
                double WbinColl = 0.;
                for (int i = 0; i<A; i++) {
                    if (nucleusB.at(i).collided == 1) {
                        xm = nucleusB.at(i).x;
                        ym = nucleusB.at(i).y;
                        W += exp((-(x - xm)*(x - xm) - (y - ym)*(y - ym))
                                 /(2.*width*width));
                    }
                    if (nucleusA.at(i).collided == 1) {
                        xm = nucleusA.at(i).x;
                        ym = nucleusA.at(i).y;
                        W += exp((-(x - xm)*(x - xm) - (y - ym)*(y - ym))
                                 /(2.*width*width));
                    }
                }
            
                for (int i = 0; i < NbinColl; i++) {
                    xm = xbinColl[i];
                    ym = ybinColl[i];
                    WbinColl += exp((-(x - xm)*(x - xm) - (y - ym)*(y - ym))
                                    /(2.*width*width));
                }
             
                double Wfull = hard*WbinColl + (1. - hard)*W;
                W = Wfull;

                (*Lneighbor)[ix][iy][0].TJb =
                                    util->cube_malloc(rk_order+1, 5, 4);
                (*Rneighbor)[ix][iy][0].TJb =
                                    util->cube_malloc(rk_order+1, 5, 4);
                (*Lneighbor)[ix][iy][1].TJb =
                                    util->cube_malloc(rk_order+1, 5, 4);
                (*Rneighbor)[ix][iy][1].TJb =
                                    util->cube_malloc(rk_order+1, 5, 4);
                (*Lneighbor)[ix][iy][0].Wmunu =
                                    util->cube_malloc(rk_order+1, 4, 4);
                (*Rneighbor)[ix][iy][0].Wmunu =
                                    util->cube_malloc(rk_order+1, 4, 4);
                (*Lneighbor)[ix][iy][1].Wmunu =
                                    util->cube_malloc(rk_order+1, 4, 4);
                (*Rneighbor)[ix][iy][1].Wmunu =
                                    util->cube_malloc(rk_order+1, 4, 4);
                (*Lneighbor)[ix][iy][0].Pimunu =
                                    util->cube_malloc(rk_order+1, 4, 4);
                (*Rneighbor)[ix][iy][0].Pimunu =
                                    util->cube_malloc(rk_order+1, 4, 4);
                (*Lneighbor)[ix][iy][1].Pimunu =
                                    util->cube_malloc(rk_order+1, 4, 4);
                (*Rneighbor)[ix][iy][1].Pimunu =
                                    util->cube_malloc(rk_order+1, 4, 4);
                (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
                (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
                (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
                (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
                (*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
                (*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
                (*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
                (*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
             
                double epsilon0 = DATA->sFactor;
                epsilon = epsilon0*W;
                 
                for (ieta = 0; ieta < DATA->neta; ieta++) {
                    eta = ((DATA->delta_eta)*(ieta+DATA->neta*rank)
                            - (DATA->eta_size)/2.0);
                    if (initializeEntropy == 0) {
                        epsilon = epsilon0*W;
                        epsilon *= eta_profile_normalisation(DATA, eta);

                        // distribution of initial baryon density:
                        rhob = DATA->rhoB0/epsilon0*epsilon;  /* in fm^-3 */
                    }
                    else if (initializeEntropy == 1) {
                        s = epsilon0*W;
                        s *= eta_profile_normalisation(DATA, eta);
                        rhob = DATA->rhoB0/epsilon0*s;
                        if (DATA->whichEOS == 1 && s > 96.288) {
                            epsilon = ((pow((4./169./pow(PI,2.)),1./3.)/30.
                                       *pow(22.5,(4./3.))*pow(s,(4./3.))*0.197)
                                      /hbarc);
                        } else if (s > 0.2) {
                            epsilon = eos->findRoot(
                                    &EOS::ssolve, rhob, s,
                                    1.15*rhob+0.001, 300.,0.001);
                        } else {
                            epsilon = 0;
                        }
                    }
                    if (epsilon < 0.00000000001)
                        epsilon = 0.00000000001;

                    // intial pressure distribution:
                    p = eos->get_pressure(epsilon, rhob);

                    // set all values in the grid element:
                    (*arena)[ix][iy][ieta].epsilon = epsilon;
                    (*arena)[ix][iy][ieta].epsilon_t = epsilon;
                    (*arena)[ix][iy][ieta].prev_epsilon = epsilon;
                    (*arena)[ix][iy][ieta].rhob = rhob;
                    (*arena)[ix][iy][ieta].rhob_t = rhob;
                    (*arena)[ix][iy][ieta].prev_rhob = rhob;
                    (*arena)[ix][iy][ieta].p = p;
                    (*arena)[ix][iy][ieta].p_t = p;

                    (*arena)[ix][iy][ieta].T =
                                        eos->get_temperature(epsilon, rhob);
                    (*arena)[ix][iy][ieta].mu = eos->get_mu(epsilon, rhob);
                    
                    (*arena)[ix][iy][ieta].TJb =
                                        util->cube_malloc(rk_order+1, 5, 4);
                    (*arena)[ix][iy][ieta].dUsup =
                                        util->cube_malloc(rk_order+1, 4, 4);
                    (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
                    (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 4);
                    (*arena)[ix][iy][ieta].theta_u =
                                        util->vector_malloc(rk_order+1);
                    (*arena)[ix][iy][ieta].sigma =
                                        util->cube_malloc(rk_order+1, 4, 4);
                    (*arena)[ix][iy][ieta].pi_b =
                                        util->vector_malloc(rk_order+1);
                    (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
                    (*arena)[ix][iy][ieta].Wmunu =
                                        util->cube_malloc(rk_order+1, 4, 4);
                    (*arena)[ix][iy][ieta].prevWmunu =
                                        util->cube_malloc(1, 4, 4);
                    (*arena)[ix][iy][ieta].Pimunu =
                                        util->cube_malloc(rk_order+1, 4, 4);
                    (*arena)[ix][iy][ieta].prevPimunu =
                                        util->cube_malloc(1, 4, 4);
                    (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(4, 4);
                    
                    /* for HIC */
                    u[0] = (*arena)[ix][iy][ieta].u[0][0] = 1.0;
                    u[3] = (*arena)[ix][iy][ieta].u[0][3] = 0.0;
                    u[1] = (*arena)[ix][iy][ieta].u[0][1] = 0.0;
                    u[2] = (*arena)[ix][iy][ieta].u[0][2] = 0.0;

                    (*arena)[ix][iy][ieta].prev_u[0][0] = 1.0;
                    (*arena)[ix][iy][ieta].prev_u[0][3] = 0.0;
                    (*arena)[ix][iy][ieta].prev_u[0][1] = 0.0;
                    (*arena)[ix][iy][ieta].prev_u[0][2] = 0.0;
             
                    (*arena)[ix][iy][ieta].pi_b[0] = 0.0;
             
                    for (mu = 0; mu < 4; mu++) {
                        /* baryon density */
                        (*arena)[ix][iy][ieta].TJb[0][4][mu] = rhob*u[mu];
                        for (nu = 0; nu < 4; nu++) {
                            (*arena)[ix][iy][ieta].TJb[0][nu][mu] =
                                ((epsilon + p)*u[mu]*u[nu]
                                 + p*(DATA->gmunu)[mu][nu]);

                            (*arena)[ix][iy][ieta].Wmunu[0][nu][mu] = 0.0;
                            (*arena)[ix][iy][ieta].prevWmunu[0][nu][mu] = 0.0;
                 
                            (*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = 0.0;
                            (*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = 0.0;
                        }/* nu */
                    }/* mu */
                }
            }
        }/* ix, iy, ieta */
        // for checking
        //if (rank == 0) {
        //    ofstream of_Glauber("initial_ed.dat");
        //    for (ix = 0; ix <= DATA->nx; ix++) {
        //        for (iy = 0; iy <= DATA->ny; iy++) {
        //            of_Glauber << scientific << setw(18) << setprecision(8)
        //                       << (*arena)[ix][iy][0].epsilon << "  ";
        //        }
        //        of_Glauber << endl;
        //    }
        //    of_Glauber.close();
        //}
				if (1) {
					ofstream of_Glauber("initial_ed.dat");
					for (ix = 0; ix <= DATA->nx; ix++) {
						for (iy = 0; iy <= DATA->ny; iy++) {
							of_Glauber << ix << " " << iy << endl;
							for (ieta = 0; ieta < DATA->neta; ieta++){
							of_Glauber << scientific << setw(15) << setprecision(8)
								<< (*arena)[ix][iy][ieta].epsilon;
							}
							of_Glauber << endl;
						}
					}
					of_Glauber.close();
				}
    } else if (DATA->Initial_profile == 8) {
        // read in the profile from file 
        // - IPGlasma initial conditions with initial flow
        size = DATA->size;

        music_message << "size=" << size;
        music_message.flush("info");
        music_message.info(" ----- information on initial distribution -----");
        music_message << "file name used: " << DATA->initName;
        music_message.flush("info");
    
        ifstream profile(DATA->initName.c_str());
    
        string dummy;
        int nx, ny, neta;
        double dx, dy, deta;
        // read the information line
        profile >> dummy >> dummy >> dummy >> dummy >> neta 
                >> dummy >> nx >> dummy >> ny 
                >> dummy >> deta >> dummy >> dx >> dummy >> dy;

        music_message << "neta=" << DATA->neta << ", nx=" << nx
                      << ", ny=" << ny << ", deta=" << DATA->delta_eta
                      << ", dx=" << dx << ", dy=" << dy;
        music_message.flush("info");

        double density, dummy1, dummy2, dummy3;
        double ux, uy, utau;
        
        double** temp_profile_ed = new double* [nx+1];
        double** temp_profile_utau = new double* [nx+1];
        double** temp_profile_ux = new double* [nx+1];
        double** temp_profile_uy = new double* [nx+1];
        for (int i = 0; i < nx+1; i++) {
            temp_profile_ed[i] = new double [ny+1];
            temp_profile_utau[i] = new double [ny+1];
            temp_profile_ux[i] = new double [ny+1];
            temp_profile_uy[i] = new double [ny+1];
        }

        //read the one slice
        for (ix = 0; ix <= DATA->nx; ix++) {
            for (iy = 0; iy <= DATA->ny; iy++) {
                profile >> dummy1 >> dummy2 >> dummy3 
                        >> density >> utau >> ux >> uy
                        >> dummy  >> dummy  >> dummy  >> dummy;
                temp_profile_ed[ix][iy] = density;
                temp_profile_utau[ix][iy] = utau;
                temp_profile_ux[ix][iy] = ux;
                temp_profile_uy[ix][iy] = uy;
                if (ix == 0 && iy == 0) {
                    DATA->x_size = -dummy2*2;
                    DATA->y_size = -dummy3*2;
                    music_message << "eta_size=" << DATA->eta_size 
                                  << ", x_size=" << DATA->x_size 
                                  << ", y_size=" << DATA->y_size;
                    music_message.flush("info");
                }
            }
        }
        profile.close();

        for(ix = 0; ix <= DATA->nx; ix++)
        {
            for(iy = 0; iy <= DATA->ny; iy++)
            {
                (*Lneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
                (*Rneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
                (*Lneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
                (*Rneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
                (*Lneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 4,4);
                (*Rneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 4,4);
                (*Lneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 4,4);
                (*Rneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 4,4);
                (*Lneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 4,4);
                (*Rneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 4,4);
                (*Lneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 4,4);
                (*Rneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 4,4);
                (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
                (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
                (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
                (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
                (*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
                (*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
                (*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
                (*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
            }
        }
       
        int entropy_flag = DATA->initializeEntropy;
        for(int ieta = 0; ieta < DATA->neta; ieta++)
        {
            double eta = (DATA->delta_eta*(ieta + DATA->neta*rank) 
                          - (DATA->eta_size)/2.0);
            double eta_envelop_ed = eta_profile_normalisation(DATA, eta);
            for(ix = 0; ix <= DATA->nx; ix++)
            {
                for(iy = 0; iy<= DATA->ny; iy++)
                {
                    rhob = 0.0;
                    if(entropy_flag == 0)
                        epsilon = (temp_profile_ed[ix][iy]*eta_envelop_ed
                                   *DATA->sFactor/hbarc);  // 1/fm^4
                    else
                    {
                        double local_sd = (temp_profile_ed[ix][iy]
                                           *DATA->sFactor*eta_envelop_ed);
                        epsilon = eos->get_s2e(local_sd, rhob);
                    }
                    if (epsilon<0.00000000001)
                        epsilon = 0.00000000001;

                    // initial pressure distribution
                    p = eos->get_pressure(epsilon, rhob);
             
                    // set all values in the grid element:
                    (*arena)[ix][iy][ieta].epsilon = epsilon;
                    (*arena)[ix][iy][ieta].epsilon_t = epsilon;
                    (*arena)[ix][iy][ieta].prev_epsilon = epsilon;
                    (*arena)[ix][iy][ieta].rhob = rhob;
                    (*arena)[ix][iy][ieta].rhob_t = rhob;
                    (*arena)[ix][iy][ieta].prev_rhob = rhob;
                    (*arena)[ix][iy][ieta].p = p;
                    (*arena)[ix][iy][ieta].p_t = p;
    
                    (*arena)[ix][iy][ieta].T = eos->get_temperature(epsilon, rhob);
                    (*arena)[ix][iy][ieta].mu = eos->get_mu(epsilon, rhob);
                        
                    (*arena)[ix][iy][ieta].TJb = util->cube_malloc(rk_order+1, 5,4);
                    (*arena)[ix][iy][ieta].dUsup = util->cube_malloc(rk_order+1, 4,4);
                    (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
                    (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 4);
                    (*arena)[ix][iy][ieta].theta_u = util->vector_malloc(rk_order+1);
                    (*arena)[ix][iy][ieta].sigma = util->cube_malloc(rk_order+1, 4, 4);
                    (*arena)[ix][iy][ieta].pi_b = util->vector_malloc(rk_order+1);
                    (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
                    (*arena)[ix][iy][ieta].Wmunu = util->cube_malloc(rk_order+1, 4,4);
                    (*arena)[ix][iy][ieta].prevWmunu = util->cube_malloc(1, 4,4);
                    (*arena)[ix][iy][ieta].Pimunu = util->cube_malloc(rk_order+1, 4,4);
                    (*arena)[ix][iy][ieta].prevPimunu = util->cube_malloc(1, 4,4);
                    (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(4,4);

                    /* for HIC */
                    u[0] = (*arena)[ix][iy][ieta].u[0][0] = temp_profile_utau[ix][iy];
                    u[1] = (*arena)[ix][iy][ieta].u[0][1] = temp_profile_ux[ix][iy];
                    u[2] = (*arena)[ix][iy][ieta].u[0][2] = temp_profile_uy[ix][iy];
                    u[3] = (*arena)[ix][iy][ieta].u[0][3] = 0.0;
                   
                    (*arena)[ix][iy][ieta].prev_u[0][0] = temp_profile_utau[ix][iy];
                    (*arena)[ix][iy][ieta].prev_u[0][1] = temp_profile_ux[ix][iy];
                    (*arena)[ix][iy][ieta].prev_u[0][2] = temp_profile_uy[ix][iy];
                    (*arena)[ix][iy][ieta].prev_u[0][3] = 0.0;
                   
                    (*arena)[ix][iy][ieta].pi_b[0] = 0.0;

                    for(int mu=0; mu<4; mu++)
                    {
                        /* baryon density */
                        (*arena)[ix][iy][ieta].TJb[0][4][mu] = rhob*u[mu];

                        for(nu=0; nu<4; nu++)
                        {
                            (*arena)[ix][iy][ieta].TJb[0][nu][mu] = (
                                (epsilon + p)*u[mu]*u[nu] 
                                + p*(DATA->gmunu)[mu][nu]);
                            (*arena)[ix][iy][ieta].Wmunu[0][nu][mu] = 0.0;
                            (*arena)[ix][iy][ieta].prevWmunu[0][nu][mu] = 0.0;
                    
                            (*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = 0.0;
                            (*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = 0.0;
                        }/* nu */
                    }/* mu */
                }
            }
        }/* ix, iy, ieta */

				if (1) {
					ofstream of_Glauber("initial_ed.dat");
					for (ix = 0; ix <= DATA->nx; ix++) {
						for (iy = 0; iy <= DATA->ny; iy++) {
							of_Glauber << ix << " " << iy << endl;
							for (ieta = 0; ieta < DATA->neta; ieta++){
							of_Glauber << scientific << setw(15) << setprecision(8)
								<< (*arena)[ix][iy][ieta].epsilon;
							}
							of_Glauber << endl;
						}
					}
					of_Glauber.close();
				}
        // clean up
        for(int i = 0; i < nx+1; i++)
        {
            delete[] temp_profile_ed[i];
            delete[] temp_profile_utau[i];
            delete[] temp_profile_ux[i];
            delete[] temp_profile_uy[i];
        }
        delete[] temp_profile_ed;
        delete[] temp_profile_utau;
        delete[] temp_profile_ux;
        delete[] temp_profile_uy;
    } else if (DATA->Initial_profile == 41) {
        // initialize hydro with vectors from JETSCAPE
        size = DATA->size;
        music_message << "size=" << size;
        music_message.flush("info");
        music_message.info(" ----- information on initial distribution -----");
        music_message << "initialized with a JETSCAPE initial condition.";
        music_message.flush("info");
    
        for (ix = 0; ix <= DATA->nx; ix++) {
            for (iy = 0; iy <= DATA->ny; iy++) {
                (*Lneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
                (*Rneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
                (*Lneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
                (*Rneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
                (*Lneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 4,4);
                (*Rneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 4,4);
                (*Lneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 4,4);
                (*Rneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 4,4);
                (*Lneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 4,4);
                (*Rneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 4,4);
                (*Lneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 4,4);
                (*Rneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 4,4);
                (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
                (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
                (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
                (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
                (*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
                (*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
                (*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
                (*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
            }
        }
       
        int entropy_flag = 1;
        for (int ieta = 0; ieta < DATA->neta; ieta++) {
            double eta_envelop_ed = 1.0;
            for (ix = 0; ix <= DATA->nx; ix++) {
                for (iy = 0; iy<= DATA->ny; iy++) {
                    int idx = ix + iy*(DATA->nx + 1);
                    rhob = 0.0;
                    if (entropy_flag == 0) {
                        epsilon = (initial_entropy_density[idx]*eta_envelop_ed
                                   *DATA->sFactor/hbarc);  // 1/fm^4
                    } else {
                        double local_sd = (initial_entropy_density[idx]
                                           *DATA->sFactor*eta_envelop_ed);
                        epsilon = eos->get_s2e(local_sd, rhob);
                    }
                    if (epsilon < 0.00000000001) {
                        epsilon = 0.00000000001;
                    }

                    // initial pressure distribution
                    p = eos->get_pressure(epsilon, rhob);
             
                    // set all values in the grid element:
                    (*arena)[ix][iy][ieta].epsilon = epsilon;
                    (*arena)[ix][iy][ieta].epsilon_t = epsilon;
                    (*arena)[ix][iy][ieta].prev_epsilon = epsilon;
                    (*arena)[ix][iy][ieta].rhob = rhob;
                    (*arena)[ix][iy][ieta].rhob_t = rhob;
                    (*arena)[ix][iy][ieta].prev_rhob = rhob;
                    (*arena)[ix][iy][ieta].p = p;
                    (*arena)[ix][iy][ieta].p_t = p;
    
                    (*arena)[ix][iy][ieta].T = eos->get_temperature(epsilon, rhob);
                    (*arena)[ix][iy][ieta].mu = eos->get_mu(epsilon, rhob);
                        
                    (*arena)[ix][iy][ieta].TJb = util->cube_malloc(rk_order+1, 5,4);
                    (*arena)[ix][iy][ieta].dUsup = util->cube_malloc(rk_order+1, 4,4);
                    (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
                    (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 4);
                    (*arena)[ix][iy][ieta].theta_u = util->vector_malloc(rk_order+1);
                    (*arena)[ix][iy][ieta].sigma = util->cube_malloc(rk_order+1, 4, 4);
                    (*arena)[ix][iy][ieta].pi_b = util->vector_malloc(rk_order+1);
                    (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
                    (*arena)[ix][iy][ieta].Wmunu = util->cube_malloc(rk_order+1, 4,4);
                    (*arena)[ix][iy][ieta].prevWmunu = util->cube_malloc(1, 4,4);
                    (*arena)[ix][iy][ieta].Pimunu = util->cube_malloc(rk_order+1, 4,4);
                    (*arena)[ix][iy][ieta].prevPimunu = util->cube_malloc(1, 4,4);
                    (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(4,4);

                    /* for HIC */
                    u[0] = (*arena)[ix][iy][ieta].u[0][0] = 1.0;
                    u[1] = (*arena)[ix][iy][ieta].u[0][1] = 0.0;
                    u[2] = (*arena)[ix][iy][ieta].u[0][2] = 0.0;
                    u[3] = (*arena)[ix][iy][ieta].u[0][3] = 0.0;
                   
                    (*arena)[ix][iy][ieta].prev_u[0][0] = 1.0;
                    (*arena)[ix][iy][ieta].prev_u[0][1] = 0.0;
                    (*arena)[ix][iy][ieta].prev_u[0][2] = 0.0;
                    (*arena)[ix][iy][ieta].prev_u[0][3] = 0.0;
                   
                    (*arena)[ix][iy][ieta].pi_b[0] = 0.0;

                    for(int mu=0; mu<4; mu++)
                    {
                        /* baryon density */
                        (*arena)[ix][iy][ieta].TJb[0][4][mu] = rhob*u[mu];

                        for(nu=0; nu<4; nu++)
                        {
                            (*arena)[ix][iy][ieta].TJb[0][nu][mu] = (
                                (epsilon + p)*u[mu]*u[nu] 
                                + p*(DATA->gmunu)[mu][nu]);
                            (*arena)[ix][iy][ieta].Wmunu[0][nu][mu] = 0.0;
                            (*arena)[ix][iy][ieta].prevWmunu[0][nu][mu] = 0.0;
                    
                            (*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = 0.0;
                            (*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = 0.0;
                        }
                    }
                }
            }
        }
        // clean up
        initial_entropy_density.clear();
    } else if (DATA->Initial_profile == 42) {
        // initialize hydro with vectors from JETSCAPE
        size = DATA->size;
        music_message << "size=" << size;
        music_message.flush("info");
        music_message.info(" ----- information on initial distribution -----");
        music_message << "initialized with a JETSCAPE initial condition.";
        music_message.flush("info");
    
        for (ix = 0; ix <= DATA->nx; ix++) {
            for (iy = 0; iy <= DATA->ny; iy++) {
                (*Lneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
                (*Rneighbor)[ix][iy][0].TJb = util->cube_malloc(rk_order+1, 5,4);
                (*Lneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
                (*Rneighbor)[ix][iy][1].TJb = util->cube_malloc(rk_order+1, 5,4);
                (*Lneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 4,4);
                (*Rneighbor)[ix][iy][0].Wmunu = util->cube_malloc(rk_order+1, 4,4);
                (*Lneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 4,4);
                (*Rneighbor)[ix][iy][1].Wmunu = util->cube_malloc(rk_order+1, 4,4);
                (*Lneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 4,4);
                (*Rneighbor)[ix][iy][0].Pimunu = util->cube_malloc(rk_order+1, 4,4);
                (*Lneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 4,4);
                (*Rneighbor)[ix][iy][1].Pimunu = util->cube_malloc(rk_order+1, 4,4);
                (*Lneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
                (*Rneighbor)[ix][iy][0].u = util->mtx_malloc(rk_order+1, 4);
                (*Lneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
                (*Rneighbor)[ix][iy][1].u = util->mtx_malloc(rk_order+1, 4);
                (*Lneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
                (*Rneighbor)[ix][iy][0].pi_b = util->vector_malloc(rk_order+1);
                (*Lneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
                (*Rneighbor)[ix][iy][1].pi_b = util->vector_malloc(rk_order+1);
            }
        }
       
        for (int ieta = 0; ieta < DATA->neta; ieta++) {
            for (ix = 0; ix <= DATA->nx; ix++) {
                for (iy = 0; iy<= DATA->ny; iy++) {
                    int idx = ix + (iy + ieta*(DATA->nx + 1))*(DATA->nx + 1);
                    rhob = 0.0;
                    epsilon = (initial_energy_density[idx]
                               *DATA->sFactor/hbarc);  // 1/fm^4
                    if (epsilon < 0.00000000001) {
                        epsilon = 0.00000000001;
                    }

                    // initial pressure distribution
                    p = eos->get_pressure(epsilon, rhob);
             
                    // set all values in the grid element:
                    (*arena)[ix][iy][ieta].epsilon = epsilon;
                    (*arena)[ix][iy][ieta].epsilon_t = epsilon;
                    (*arena)[ix][iy][ieta].prev_epsilon = epsilon;
                    (*arena)[ix][iy][ieta].rhob = rhob;
                    (*arena)[ix][iy][ieta].rhob_t = rhob;
                    (*arena)[ix][iy][ieta].prev_rhob = rhob;
                    (*arena)[ix][iy][ieta].p = p;
                    (*arena)[ix][iy][ieta].p_t = p;
    
                    (*arena)[ix][iy][ieta].T = eos->get_temperature(epsilon, rhob);
                    (*arena)[ix][iy][ieta].mu = eos->get_mu(epsilon, rhob);
                        
                    (*arena)[ix][iy][ieta].TJb = util->cube_malloc(rk_order+1, 5,4);
                    (*arena)[ix][iy][ieta].dUsup = util->cube_malloc(rk_order+1, 4,4);
                    (*arena)[ix][iy][ieta].u = util->mtx_malloc(rk_order+1, 4);
                    (*arena)[ix][iy][ieta].a = util->mtx_malloc(rk_order+1, 4);
                    (*arena)[ix][iy][ieta].theta_u = util->vector_malloc(rk_order+1);
                    (*arena)[ix][iy][ieta].sigma = util->cube_malloc(rk_order+1, 4, 4);
                    (*arena)[ix][iy][ieta].pi_b = util->vector_malloc(rk_order+1);
                    (*arena)[ix][iy][ieta].prev_u = util->mtx_malloc(1, 4);
                    (*arena)[ix][iy][ieta].Wmunu = util->cube_malloc(rk_order+1, 4,4);
                    (*arena)[ix][iy][ieta].prevWmunu = util->cube_malloc(1, 4,4);
                    (*arena)[ix][iy][ieta].Pimunu = util->cube_malloc(rk_order+1, 4,4);
                    (*arena)[ix][iy][ieta].prevPimunu = util->cube_malloc(1, 4,4);
                    (*arena)[ix][iy][ieta].W_prev = util->mtx_malloc(4,4);

                    /* for HIC */
                    u[1] = (*arena)[ix][iy][ieta].u[0][1] = initial_u_x[idx];
                    u[2] = (*arena)[ix][iy][ieta].u[0][2] = initial_u_y[idx];
                    u[3] = (*arena)[ix][iy][ieta].u[0][3] = DATA->tau0*initial_u_eta[idx];
                    u[0] = sqrt(1. + u[1]*u[1] + u[2]*u[2] + u[3]*u[3]);
                    (*arena)[ix][iy][ieta].u[0][0] = u[0];
                   
                    (*arena)[ix][iy][ieta].prev_u[0][0] = u[0]; 
                    (*arena)[ix][iy][ieta].prev_u[0][1] = u[1]; 
                    (*arena)[ix][iy][ieta].prev_u[0][2] = u[2]; 
                    (*arena)[ix][iy][ieta].prev_u[0][3] = u[3]; 
                   
                    (*arena)[ix][iy][ieta].pi_b[0] = initial_bulk_pi[idx]/hbarc;

                    (*arena)[ix][iy][ieta].Wmunu[0][0][0] = initial_pi_00[idx]/hbarc;
                    (*arena)[ix][iy][ieta].Wmunu[0][0][1] = initial_pi_01[idx]/hbarc;
                    (*arena)[ix][iy][ieta].Wmunu[0][0][2] = initial_pi_02[idx]/hbarc;
                    (*arena)[ix][iy][ieta].Wmunu[0][0][3] = DATA->tau0*initial_pi_03[idx]/hbarc;
                    (*arena)[ix][iy][ieta].Wmunu[0][1][1] = initial_pi_11[idx]/hbarc;
                    (*arena)[ix][iy][ieta].Wmunu[0][1][2] = initial_pi_12[idx]/hbarc;
                    (*arena)[ix][iy][ieta].Wmunu[0][1][3] = DATA->tau0*initial_pi_13[idx]/hbarc;
                    (*arena)[ix][iy][ieta].Wmunu[0][2][2] = initial_pi_22[idx]/hbarc;
                    (*arena)[ix][iy][ieta].Wmunu[0][2][3] = DATA->tau0*initial_pi_23[idx]/hbarc;
                    (*arena)[ix][iy][ieta].Wmunu[0][3][3] = DATA->tau0*DATA->tau0*initial_pi_33[idx]/hbarc;

                    for (int mu=0; mu<4; mu++) {
                        /* baryon density */
                        (*arena)[ix][iy][ieta].TJb[0][4][mu] = rhob*u[mu];

                        for (nu=0; nu<4; nu++) {
                            (*arena)[ix][iy][ieta].TJb[0][nu][mu] = (
                                (epsilon + p)*u[mu]*u[nu] 
                                + p*(DATA->gmunu)[mu][nu]);
                            
                            (*arena)[ix][iy][ieta].prevWmunu[0][nu][mu] = (
                                    (*arena)[ix][iy][ieta].Wmunu[0][nu][mu]);
                    
                            (*arena)[ix][iy][ieta].Pimunu[0][nu][mu] = (
                                (*arena)[ix][iy][ieta].pi_b[0]*(
                                    (DATA->gmunu)[nu][mu] + u[mu]*u[nu]));

                            (*arena)[ix][iy][ieta].prevPimunu[0][nu][mu] = (
                                    (*arena)[ix][iy][ieta].Pimunu[0][nu][mu]);
                        }
                    }
                }
            }
        }
        // clean up
        initial_energy_density.clear();
        initial_u_tau.clear();
        initial_u_x.clear();
        initial_u_y.clear();
        initial_u_eta.clear();
        initial_pi_00.clear();
        initial_pi_01.clear();
        initial_pi_02.clear();
        initial_pi_03.clear();
        initial_pi_11.clear();
        initial_pi_12.clear();
        initial_pi_13.clear();
        initial_pi_22.clear();
        initial_pi_23.clear();
        initial_pi_33.clear();
        initial_bulk_pi.clear();
    }

    music_message.info("initial distribution done.");
    return 1;
}


double Init::TATarget(InitData *DATA, double r)
{
    return glauber->InterNuTInST(r)/(DATA->SigmaNN/10.);
}

double Init::TAProjectile(InitData *DATA, double r)
{
    return glauber->InterNuPInSP(r)/(DATA->SigmaNN/10.);
}

double Init::eta_profile_normalisation(InitData *DATA, double eta)
{
    double res;

    if (DATA->boost_invariant) {
    	res=1.0;
    }
    else {
	    //Hirano's plateau + Gaussian fall-off
	    if (DATA->initial_eta_profile == 1)
	    {
		    double exparg1 = (fabs(eta) - DATA->eta_flat/2.0)/DATA->eta_fall_off;
		    double exparg = exparg1*exparg1/2.0;
		    res=exp(-exparg*theta(exparg1));
	    }
	    //Woods-Saxon
	    //The radius is set to be half of DATA->eta_flat
	    //The diffusiveness is set to DATA->eta_fall_off
	    else if (DATA->initial_eta_profile == 2)
	    {
		    double ws_R = DATA->eta_flat/2.0;
		    double ws_a = DATA->eta_fall_off;
		    res = (1.0+exp(-ws_R/ws_a))/(1.0+exp((abs(eta)-ws_R)/ws_a));
	    }
	    else
	    {
		    music_message.error("initial_eta_profile out of range.");
		    exit(0);
	    }
    }
    return res;
}
