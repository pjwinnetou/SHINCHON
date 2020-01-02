// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Maxime Dion, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-Francois Paquet, Bjoern Schenke, Chun Shen, Clint Young

#include "util.h"
#include "grid_info.h"

using namespace std;

Grid_info::Grid_info(InitData* DATA_in)
{
   DATA_ptr = DATA_in;
}

Grid_info::~Grid_info()
{
}


void Grid_info::OutputEvolutionDataXYEta(Grid ***arena, InitData *DATA,
                                  EOS *eos, double tau, int size, int rank,
                                  HydroinfoMUSIC *hydro_info_ptr) {
    Util *util = new Util();
  
    int nx = DATA->nx;
    int ny = DATA->ny;
    int neta = DATA->neta;

    // MPI send and receive
    int sizeOfData = (nx+1)*(ny+1)*(neta);
    int to, from;
  
    double *eps;
    double *rhob;
    double *utau;
    double *ux;
    double *uy;
    double *ueta;
    double *p;
    double *Wtautau, *Wtaux, *Wtauy, *Wtaueta, *Wxx, *Wxy, *Wxeta, *Wyy;
    double *Wyeta, *Wetaeta;
    double *bulk_pressure;

    eps = (double *)malloc(sizeof(double)*sizeOfData);
    rhob = (double *)malloc(sizeof(double)*sizeOfData);
    utau = (double *)malloc(sizeof(double)*sizeOfData);
    ux = (double *)malloc(sizeof(double)*sizeOfData);
    uy = (double *)malloc(sizeof(double)*sizeOfData);
    ueta = (double *)malloc(sizeof(double)*sizeOfData);
    p = (double *)malloc(sizeof(double)*sizeOfData);

    // added by Maxime
    Wtautau = (double *)malloc(sizeof(double)*sizeOfData);
    Wtaux = (double *)malloc(sizeof(double)*sizeOfData);
    Wtauy = (double *)malloc(sizeof(double)*sizeOfData);
    Wtaueta = (double *)malloc(sizeof(double)*sizeOfData);
    Wxx = (double *)malloc(sizeof(double)*sizeOfData);
    Wxy = (double *)malloc(sizeof(double)*sizeOfData);
    Wxeta = (double *)malloc(sizeof(double)*sizeOfData);
    Wyy = (double *)malloc(sizeof(double)*sizeOfData); 
    Wyeta = (double *)malloc(sizeof(double)*sizeOfData);
    Wetaeta = (double *)malloc(sizeof(double)*sizeOfData);
    // end of "added by Maxime"

    bulk_pressure = (double *)malloc(sizeof(double)*sizeOfData);
    int position;
    if (rank > 0) {
        // send all to rank 0
        to = 0;
        for (int ix = 0; ix <= nx; ix++) {
	        for (int iy = 0; iy <= ny; iy++) {
	            for (int ieta = 0; ieta < neta; ieta++) {
		            position = ieta + neta*(ix + (nx + 1)*iy);
		            eps[position] = arena[ix][iy][ieta].epsilon;
		            rhob[position] = arena[ix][iy][ieta].rhob;
		            utau[position] = arena[ix][iy][ieta].u[0][0];
		            ux[position] = arena[ix][iy][ieta].u[0][1];
		            uy[position] = arena[ix][iy][ieta].u[0][2];
		            ueta[position] = arena[ix][iy][ieta].u[0][3];
		            p[position] = arena[ix][iy][ieta].p; 

		            // added by Maxime
		            Wtautau[position] = arena[ix][iy][ieta].Wmunu[0][0][0];
		            Wtaux[position] = arena[ix][iy][ieta].Wmunu[0][0][1];
		            Wtauy[position] = arena[ix][iy][ieta].Wmunu[0][0][2];
		            Wtaueta[position] = arena[ix][iy][ieta].Wmunu[0][0][3];
		            Wxx[position] = arena[ix][iy][ieta].Wmunu[0][1][1];
		            Wxy[position] = arena[ix][iy][ieta].Wmunu[0][1][2];
		            Wxeta[position] = arena[ix][iy][ieta].Wmunu[0][1][3];
		            Wyy[position] = arena[ix][iy][ieta].Wmunu[0][2][2];
		            Wyeta[position] = arena[ix][iy][ieta].Wmunu[0][2][3];
		            Wetaeta[position] = arena[ix][iy][ieta].Wmunu[0][3][3];
		            // end of "added by Maxime"

		            bulk_pressure[position] = arena[ix][iy][ieta].pi_b[0];
		        }
	        }
	    }
        MPI_Send(eps,sizeOfData,MPI_DOUBLE,to,1, MPI_COMM_WORLD);
        MPI_Send(rhob,sizeOfData,MPI_DOUBLE,to,2, MPI_COMM_WORLD);
        MPI_Send(utau,sizeOfData,MPI_DOUBLE,to,3, MPI_COMM_WORLD);
        MPI_Send(ux,sizeOfData,MPI_DOUBLE,to,4, MPI_COMM_WORLD);
        MPI_Send(uy,sizeOfData,MPI_DOUBLE,to,5, MPI_COMM_WORLD);
        MPI_Send(ueta,sizeOfData,MPI_DOUBLE,to,6, MPI_COMM_WORLD);
        MPI_Send(p,sizeOfData,MPI_DOUBLE,to,17, MPI_COMM_WORLD);

        // added by Maxime
        MPI_Send(Wtautau,sizeOfData,MPI_DOUBLE,to,7, MPI_COMM_WORLD);
        MPI_Send(Wtaux,sizeOfData,MPI_DOUBLE,to,8, MPI_COMM_WORLD);
        MPI_Send(Wtauy,sizeOfData,MPI_DOUBLE,to,9, MPI_COMM_WORLD);
        MPI_Send(Wtaueta,sizeOfData,MPI_DOUBLE,to,10, MPI_COMM_WORLD);
        MPI_Send(Wxx,sizeOfData,MPI_DOUBLE,to,11, MPI_COMM_WORLD);
        MPI_Send(Wxy,sizeOfData,MPI_DOUBLE,to,12, MPI_COMM_WORLD);
        MPI_Send(Wxeta,sizeOfData,MPI_DOUBLE,to,13, MPI_COMM_WORLD);
        MPI_Send(Wyy,sizeOfData,MPI_DOUBLE,to,14, MPI_COMM_WORLD);
        MPI_Send(Wyeta,sizeOfData,MPI_DOUBLE,to,15, MPI_COMM_WORLD);
        MPI_Send(Wetaeta,sizeOfData,MPI_DOUBLE,to,16, MPI_COMM_WORLD);
        // end of "added by Maxime"

        MPI_Send(bulk_pressure,sizeOfData,MPI_DOUBLE,to,18, MPI_COMM_WORLD);
    }
  
    if (rank == 0) {
        double ***epsFrom;
        double ***rhobFrom;
        double ***utauFrom;
        double ***uxFrom;
        double ***uyFrom;
        double ***uetaFrom;
        double ***WtautauFrom, ***WtauxFrom, ***WtauyFrom, ***WtauetaFrom;
        double ***WxxFrom, ***WxyFrom, ***WxetaFrom;
        double ***WyyFrom, ***WyetaFrom, ***WetaetaFrom; // added by Maxime
        double ***pFrom;
        double ***bulkPressureFrom;
      
        epsFrom = util->cube_malloc(nx+1,ny+1,size*neta);
        rhobFrom = util->cube_malloc(nx+1,ny+1,size*neta);
        utauFrom = util->cube_malloc(nx+1,ny+1,size*neta);
        uxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
        uyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
        uetaFrom = util->cube_malloc(nx+1,ny+1,size*neta);
        pFrom = util->cube_malloc(nx+1,ny+1,size*neta);

        // added by Maxime
        WtautauFrom = util->cube_malloc(nx+1,ny+1,size*neta);
        WtauxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
        WtauyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
        WtauetaFrom = util->cube_malloc(nx+1,ny+1,size*neta);
        WxxFrom = util->cube_malloc(nx+1,ny+1,size*neta);
        WxyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
        WxetaFrom = util->cube_malloc(nx+1,ny+1,size*neta); 
        WyyFrom = util->cube_malloc(nx+1,ny+1,size*neta);
        WyetaFrom = util->cube_malloc(nx+1,ny+1,size*neta);
        WetaetaFrom = util->cube_malloc(nx+1,ny+1,size*neta);  
        // end of "added by Maxime"   

        bulkPressureFrom = util->cube_malloc(nx+1,ny+1,size*neta);  
      
        for (int ix = 0; ix <= nx; ix++) {
	        for (int iy = 0; iy <= ny; iy++) {
	            for (int ieta = 0; ieta < neta; ieta++) {
		            epsFrom[ix][iy][ieta] = arena[ix][iy][ieta].epsilon;
		            rhobFrom[ix][iy][ieta] = arena[ix][iy][ieta].rhob;
		            utauFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][0];
		            uxFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][1];
		            uyFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][2];
		            uetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].u[0][3];
		            pFrom[ix][iy][ieta] = arena[ix][iy][ieta].p;

                    // added by Maxime
                    WtautauFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][0][0]; 
                    WtauxFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][0][1];
                    WtauyFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][0][2];
                    WtauetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][0][3];
                    WxxFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][1][1];
                    WxyFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][1][2];
                    WxetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][1][3];
                    WyyFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][2][2];
                    WyetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][2][3];
                    WetaetaFrom[ix][iy][ieta] = arena[ix][iy][ieta].Wmunu[0][3][3];
	  		        // end of "added by Maxime"

		            bulkPressureFrom[ix][iy][ieta] = arena[ix][iy][ieta].pi_b[0];
		        }
	        }
	    }
        for (int irank=1; irank<size; irank++) {
	        from = irank;
	        MPI_Recv(eps,sizeOfData,MPI_DOUBLE,from,1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	        MPI_Recv(rhob,sizeOfData,MPI_DOUBLE,from,2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	        MPI_Recv(utau,sizeOfData,MPI_DOUBLE,from,3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	        MPI_Recv(ux,sizeOfData,MPI_DOUBLE,from,4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	        MPI_Recv(uy,sizeOfData,MPI_DOUBLE,from,5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	        MPI_Recv(ueta,sizeOfData,MPI_DOUBLE,from,6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	        MPI_Recv(p,sizeOfData,MPI_DOUBLE,from,17, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	        // added by Maxime
	        MPI_Recv(Wtautau,sizeOfData,MPI_DOUBLE,from,7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	        MPI_Recv(Wtaux,sizeOfData,MPI_DOUBLE,from,8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
 	        MPI_Recv(Wtauy,sizeOfData,MPI_DOUBLE,from,9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	        MPI_Recv(Wtaueta,sizeOfData,MPI_DOUBLE,from,10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	        MPI_Recv(Wxx,sizeOfData,MPI_DOUBLE,from,11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	        MPI_Recv(Wxy,sizeOfData,MPI_DOUBLE,from,12, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	        MPI_Recv(Wxeta,sizeOfData,MPI_DOUBLE,from,13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	        MPI_Recv(Wyy,sizeOfData,MPI_DOUBLE,from,14, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	        MPI_Recv(Wyeta,sizeOfData,MPI_DOUBLE,from,15, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   	        MPI_Recv(Wetaeta,sizeOfData,MPI_DOUBLE,from,16, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	        // end of "added by Maxime"

	        MPI_Recv(bulk_pressure,sizeOfData,MPI_DOUBLE,from,18,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  
	        for (int ix = 0; ix <= nx; ix++) {
	            for (int iy = 0; iy <= ny; iy++) {
		            for (int ieta = 0; ieta < neta; ieta++) {
		                position = ieta+(neta*(ix + ((nx+1)*iy)));
		                epsFrom[ix][iy][ieta+irank*neta] = eps[position];
		                rhobFrom[ix][iy][ieta+irank*neta] = rhob[position];
		                utauFrom[ix][iy][ieta+irank*neta] = utau[position];
		                uxFrom[ix][iy][ieta+irank*neta] = ux[position];
		                uyFrom[ix][iy][ieta+irank*neta] = uy[position];
		                uetaFrom[ix][iy][ieta+irank*neta] = ueta[position];
		                pFrom[ix][iy][ieta+irank*neta] = p[position];

		                // added by Maxime
		                WtautauFrom[ix][iy][ieta+irank*neta] = Wtautau[position];
		                WtauxFrom[ix][iy][ieta+irank*neta] = Wtaux[position];
  		                WtauyFrom[ix][iy][ieta+irank*neta] = Wtauy[position];
		                WtauetaFrom[ix][iy][ieta+irank*neta] = Wtaueta[position];
		                WxxFrom[ix][iy][ieta+irank*neta] = Wxx[position];
 		                WxyFrom[ix][iy][ieta+irank*neta] = Wxy[position];
		                WxetaFrom[ix][iy][ieta+irank*neta] = Wxeta[position];
		                WyyFrom[ix][iy][ieta+irank*neta] = Wyy[position];
		                WyetaFrom[ix][iy][ieta+irank*neta] = Wyeta[position];
		                WetaetaFrom[ix][iy][ieta+irank*neta] = Wetaeta[position];
		                // end of added by Maxime

		            }
		        }
	        }
	    }

      
        // set up file name
        const string out_name_xyeta = "evolution_xyeta.dat";
        const string out_name_W_xyeta = "evolution_Wmunu_over_epsilon_plus_P_xyeta.dat";
        const string out_name_bulk_xyeta = "evolution_bulk_pressure_xyeta.dat";
        string out_open_mode;
        FILE *out_file_xyeta;
        FILE *out_file_W_xyeta;
        FILE *out_file_bulk_xyeta;

        // If it's the first timestep, overwrite the previous file
        if (tau == DATA->tau0) {
           out_open_mode = "w";
        }
        else {
           out_open_mode = "a";	
        }
        //If we output in binary, set the mode accordingly
        if (0 == DATA->outputBinaryEvolution) {
          out_open_mode += "b";
        }
        
        if (DATA->outputEvolutionData == 1) {
            out_file_xyeta =
                fopen(out_name_xyeta.c_str(), out_open_mode.c_str());
            out_file_W_xyeta =
                fopen(out_name_W_xyeta.c_str(), out_open_mode.c_str());
            out_file_bulk_xyeta =
                fopen(out_name_bulk_xyeta.c_str(), out_open_mode.c_str());
        }
      
        // Although it is a little confusing,
        // it is easiest to use (tau, x, y, eta) coordinates
        // and save at these points vx, vy, and vz. -CFY 11/16/2010
        float T1, u01, ux1, uy1, uz1, ueta1, utau1, epsilon1, rhob1, QGPfrac1;
        float entropy1;
        double eta;
        float Wtt, Wtx, Wty, Wtz, Wzz, Wxz, Wyz; // added by Maxime
        float Wtautau, Wtaux, Wtauy, Wtaueta, Wxx, Wxy, Wxeta, Wyy;
        float Wyeta, Wetaeta; // added by Maxime
        float div_factor, pressure1;

        int n_eta_step = DATA_ptr->output_evolution_every_N_eta;
        int n_x_step = DATA_ptr->output_evolution_every_N_x;
        int n_y_step = DATA_ptr->output_evolution_every_N_y;

        // No interpolation is necessary here!
        for (int ieta = 0; ieta < (DATA->neta)*size; ieta += n_eta_step) {
            // Recall that the eta grid is shifted by one cell
            // in the negative direction
            // For a boost-invariant system,
            // we want to output this information at eta=0,
            // not eta=-(DATA->eta_size)/2.0
            if (DATA->boost_invariant) {
                eta = 0.0;
            } else {
                eta = ((double)ieta)*(DATA->delta_eta)-(DATA->eta_size)/2.0;
            }

	        for (int iy = 0; iy <= DATA->ny; iy += n_y_step) {
	            for (int ix = 0; ix <= DATA->nx; ix += n_x_step) {
		            epsilon1 = epsFrom[ix][iy][ieta];
		            pressure1 = pFrom[ix][iy][ieta];
		            rhob1 = rhobFrom[ix][iy][ieta];
		            utau1 = utauFrom[ix][iy][ieta];
		            ux1 = uxFrom[ix][iy][ieta];
		            uy1 = uyFrom[ix][iy][ieta];
		            ueta1 = uetaFrom[ix][iy][ieta];
		
		            u01 = ueta1*sinh(eta)+utau1*cosh(eta); // = gamma factor
		            ux1 = ux1/u01;
		            uy1 = uy1/u01;
		            uz1 = ueta1*cosh(eta)+utau1*sinh(eta);
		            uz1 /= u01;
	
		            T1 = eos->get_temperature(epsilon1, rhob1);
                    entropy1 = eos->get_entropy(epsilon1, rhob1);

                    if (DATA_ptr->store_hydro_info_in_memory == 1) {
                        hydro_info_ptr->dump_ideal_info_to_memory(tau,
                                        epsilon1, pressure1, entropy1, T1,
                                        ux1, uy1, ueta1);
                    }

                    QGPfrac1 = 0.0;
		            const float local_bulk=bulkPressureFrom[ix][iy][ieta];
		            const float cs2=eos->p_e_func(epsilon1,rhob1);
                    div_factor=(epsilon1+pressure1);
                    // I need <del_mu u_nu> not Wmunu
                    // so I divide by the viscosity -Maxime
		            Wtautau = WtautauFrom[ix][iy][ieta]/(div_factor); 
		            Wtaux = WtauxFrom[ix][iy][ieta]/(div_factor);
		            Wtauy = WtauyFrom[ix][iy][ieta]/(div_factor);
		            Wtaueta = WtauetaFrom[ix][iy][ieta]/(div_factor);
		            Wxx = WxxFrom[ix][iy][ieta]/(div_factor);
		            Wxy = WxyFrom[ix][iy][ieta]/(div_factor);
		            Wxeta = WxetaFrom[ix][iy][ieta]/(div_factor);
		            Wyy = WyyFrom[ix][iy][ieta]/(div_factor);
		            Wyeta = WyetaFrom[ix][iy][ieta]/(div_factor);
		            Wetaeta = WetaetaFrom[ix][iy][ieta]/(div_factor);

		            Wtt = pow(cosh(eta),2)*Wtautau + pow(tau*sinh(eta),2)*Wetaeta*pow(tau,-2) + 2*tau*cosh(eta)*sinh(eta)*Wtaueta*pow(tau,-1);
		            Wtx = cosh(eta)*Wtaux + tau*sinh(eta)*Wxeta*pow(tau,-1);
		            Wty = cosh(eta)*Wtauy + tau*sinh(eta)*Wyeta*pow(tau,-1);
		            Wtz = cosh(eta)*sinh(eta)*Wtautau + tau*( pow(cosh(eta),2) + pow(sinh(eta),2) )*Wtaueta*pow(tau,-1) + tau*tau*cosh(eta)*sinh(eta)*Wetaeta*pow(tau,-2);
		            Wxz = sinh(eta)*Wtaux + tau*cosh(eta)*Wxeta*pow(tau,-1);
		            Wyz = sinh(eta)*Wtauy + tau*cosh(eta)*Wyeta*pow(tau,-1);
		            Wzz = pow(sinh(eta),2)*Wtautau + pow(tau*cosh(eta),2)*Wetaeta*pow(tau,-2) + 2*tau*cosh(eta)*sinh(eta)*Wtaueta*pow(tau,-1);
		
                    if (DATA->outputEvolutionData == 1) {
                        // exclude the actual coordinates from the output to save space:
                        if (0 == DATA->outputBinaryEvolution) {
		                    fprintf(out_file_xyeta,"%e %e %e %e %e\n", T1*hbarc, QGPfrac1, ux1, uy1, uz1);
		                    if (DATA->viscosity_flag) {
			                    if (1 == DATA->turn_on_shear) {
			                        fprintf(out_file_W_xyeta,"%e %e %e %e %e %e %e %e %e %e\n",Wtt,Wtx,Wty,Wtz,Wxx,Wxy,Wxz,Wyy,Wyz,Wzz); 
			                    }
			                    if (1 == DATA->turn_on_bulk) {
			                        fprintf(out_file_bulk_xyeta,"%e %e %e\n", local_bulk, epsilon1+pressure1, cs2);
			                    }
		                    }
		                } else {
		                    float array[]={T1*hbarc, QGPfrac1, ux1, uy1, uz1};
		                    fwrite(array,sizeof(float),5,out_file_xyeta);
		                    // Write Wmunu/shear only if viscosity is on
                            // --- no need to fill a file with zeros in the ideal case
		                    if (DATA->viscosity_flag) {
			                    if (1 == DATA->turn_on_shear) {
				                    float array2[]={Wtt,Wtx,Wty,Wtz,Wxx,Wxy,Wxz,Wyy,Wyz,Wzz};
				                    fwrite(array2,sizeof(float),10,out_file_W_xyeta);
			                    }
			                    if (1 == DATA->turn_on_bulk) {
				                    float array2[]={local_bulk, epsilon1+pressure1, cs2};
			                        fwrite(array2,sizeof(float),3,out_file_bulk_xyeta);
			                    }
		                    }
		                }
                    }
	            }/* ix */
	        }/* iy */
        }/* ieta */
        if (DATA->outputEvolutionData == 1) {
            fclose(out_file_xyeta);
            fclose(out_file_W_xyeta);
            fclose(out_file_bulk_xyeta);
        }
        
        /*End of hydro output in tau,x,y,eta*/
        util->cube_free(epsFrom,nx+1,ny+1,size*neta);
        util->cube_free(rhobFrom,nx+1,ny+1,size*neta);
        util->cube_free(utauFrom,nx+1,ny+1,size*neta);
        util->cube_free(uxFrom,nx+1,ny+1,size*neta);
        util->cube_free(uyFrom,nx+1,ny+1,size*neta);
        util->cube_free(uetaFrom,nx+1,ny+1,size*neta);
        util->cube_free(pFrom,nx+1,ny+1,size*neta);
        
        // added by Maxime
        util->cube_free(WtautauFrom,nx+1,ny+1,size*neta);
        util->cube_free(WtauxFrom,nx+1,ny+1,size*neta);
        util->cube_free(WtauyFrom,nx+1,ny+1,size*neta);
        util->cube_free(WtauetaFrom,nx+1,ny+1,size*neta);
        util->cube_free(WxxFrom,nx+1,ny+1,size*neta);
        util->cube_free(WxyFrom,nx+1,ny+1,size*neta);
        util->cube_free(WxetaFrom,nx+1,ny+1,size*neta);
        util->cube_free(WyyFrom,nx+1,ny+1,size*neta);
        util->cube_free(WyetaFrom,nx+1,ny+1,size*neta);
        util->cube_free(WetaetaFrom,nx+1,ny+1,size*neta);
        // end of added by Maxime
        util->cube_free(bulkPressureFrom,nx+1,ny+1,size*neta);
    }
    free(eps);
    free(rhob);
    free(utau);
    free(ux);
    free(uy);
    free(ueta);
    free(p);

    // added by Maxime
    free(Wtautau);
    free(Wtaux);
    free(Wtauy);
    free(Wtaueta);
    free(Wxx);
    free(Wxy);
    free(Wxeta);
    free(Wyy);
    free(Wyeta);
    free(Wetaeta);
    // end of "added by Maxime"
  
    free(bulk_pressure);

    delete(util);
}/* OutputEvolutionDataXYEta */


void Grid_info::Output_hydro_information_header(InitData *DATA, EOS *eos)
{
	string fname = "hydro_info_header_h";

	//Open output file
	ofstream outfile;
	outfile.open(fname.c_str());

      int grid_nx = ceil(((double)(DATA->nx+1))/DATA->output_evolution_every_N_x);
      int grid_ny = ceil(((double)(DATA->ny+1))/DATA->output_evolution_every_N_y);
      int mpi_size;
      MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
      int grid_neta = ceil(((double)(DATA->neta*mpi_size))
                           /DATA->output_evolution_every_N_eta);

	outfile << "const int MUSIC_real_nx=" << grid_nx << ";" << endl;
	outfile << "const int MUSIC_real_ny=" << grid_ny << ";" << endl;
	//DATA->neta is _not_ the actual number of cells in eta, it is the number of cells in eta _per_ processor
	outfile << "const int MUSIC_real_neta=" << grid_neta << ";" << endl;

	//double x_size; /* in fermi -x_size/2 < x < x_size/2 */
	//double y_size; /* in fermi, ditto */
	//double eta_size; /* ditto */

	outfile << "const double MUSIC_tau0=" << DATA->tau0 << ";" << endl;

	outfile << "const double MUSIC_dx=" << DATA->delta_x*DATA->output_evolution_every_N_x << ";" << endl;
	outfile << "const double MUSIC_dy=" << DATA->delta_y*DATA->output_evolution_every_N_y << ";" << endl;
	outfile << "const double MUSIC_deta=" << DATA->delta_eta*DATA->output_evolution_every_N_eta << ";" << endl;
	outfile << "const double MUSIC_effective_dtau=" << DATA->output_evolution_every_N_timesteps*DATA->delta_tau << ";" << endl;

	outfile << "const bool MUSIC_with_shear_viscosity=" << ((DATA->viscosity_flag)&&(DATA->turn_on_shear)) << ";\n";
	outfile << "const bool MUSIC_with_bulk_viscosity=" << ((DATA->viscosity_flag)&&(DATA->turn_on_bulk)) << ";\n";

      outfile << "const double MUSIC_kinetic_FO_temperature_in_GeV=";
      if (DATA->useEpsFO) {
      	outfile << eos->get_temperature(DATA->epsilonFreeze/hbarc,0.0)*hbarc;
      }
      else {
      	outfile << DATA->TFO;
      }
      outfile << ";" << endl;
      outfile << "const int MUSIC_use_temperature_FO=1;" << endl;

	outfile << "const bool MUSIC_outputBinaryEvolution=" << DATA->outputBinaryEvolution << ";" << endl;
	outfile.close();
}


//! This function outputs quantities for Gubser flow check
void Grid_info::Gubser_flow_check_file(Grid ***arena, EOS *eos_ptr,
                                       double tau) {
    ostringstream filename;
    filename << "Gubser_flow_check_tau_" << tau << ".dat";
    ofstream output_file(filename.str().c_str());

    double unit_convert = 0.19733;  // hbarC
    double dx = DATA_ptr->delta_x;
    double x_min = -DATA_ptr->x_size/2.;
    double dy = DATA_ptr->delta_y;
    double y_min = -DATA_ptr->y_size/2.;
    for (int ix = 0; ix <= DATA_ptr->nx; ix++) {
        double x_local = x_min + ix*dx;
        for (int iy = 0; iy <= DATA_ptr->ny; iy++) {
            double y_local = y_min + iy*dy;
            double e_local = arena[ix][iy][0].epsilon;
            double rhob_local = arena[ix][iy][0].rhob;
            double T_local = eos_ptr->get_temperature(e_local, 0.0);
            output_file << scientific << setprecision(8) << setw(18)
                        << x_local << "  " << y_local << "  "
                        << e_local*unit_convert << "  " << rhob_local << "  "
                        << T_local*unit_convert << "  "
                        << arena[ix][iy][0].u[0][1] << "  "
                        << arena[ix][iy][0].u[0][2] << "  "
                        << arena[ix][iy][0].Wmunu[0][1][1]*unit_convert << "  "
                        << arena[ix][iy][0].Wmunu[0][2][2]*unit_convert << "  "
                        << arena[ix][iy][0].Wmunu[0][1][2]*unit_convert << "  "
                        << arena[ix][iy][0].Wmunu[0][3][3]*unit_convert << "  "
                        << endl;
        }
    }
    output_file.close();
}


//! This function outputs the evolution of hydrodynamic variables at a
//! give fluid cell
void Grid_info::monitor_fluid_cell(Grid ***arena, int ix, int iy, int ieta,
                                   double tau) {
    ostringstream filename;
    filename << "monitor_fluid_cell_ix_" << ix << "_iy_" << iy
             << "_ieta_" << ieta << ".dat";
    ofstream output_file(filename.str().c_str(),
                         std::ofstream::out | std::ofstream::app);
    output_file << scientific << setprecision(8) << setw(18)
                << tau << "  " << arena[ix][iy][ieta].epsilon << "  "
                << arena[ix][iy][ieta].rhob << "  "
                << arena[ix][iy][ieta].u[0][0] << "  "
                << arena[ix][iy][ieta].u[0][1] << "  "
                << arena[ix][iy][ieta].u[0][2] << "  "
                << arena[ix][iy][ieta].u[0][3] << "  "
                << arena[ix][iy][ieta].Wmunu[0][0][0] << "  "
                << arena[ix][iy][ieta].Wmunu[0][0][1] << "  "
                << arena[ix][iy][ieta].Wmunu[0][0][2] << "  "
                << arena[ix][iy][ieta].Wmunu[0][0][3] << "  "
                << arena[ix][iy][ieta].Wmunu[0][1][1] << "  "
                << arena[ix][iy][ieta].Wmunu[0][1][2] << "  "
                << arena[ix][iy][ieta].Wmunu[0][1][3] << "  "
                << arena[ix][iy][ieta].Wmunu[0][2][2] << "  "
                << arena[ix][iy][ieta].Wmunu[0][2][3] << "  "
                << arena[ix][iy][ieta].Wmunu[0][3][3] << "  "
                << arena[ix][iy][ieta].pi_b[0] << "  "
                << endl;
    output_file.close();
}
