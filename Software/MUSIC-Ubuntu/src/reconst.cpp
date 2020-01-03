// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

#include "reconst.h"

using namespace std;

Reconst::Reconst(EOS *eosIn, InitData *DATA_in) {
    eos = eosIn;
    DATA_ptr = DATA_in;
    echo_level = DATA_ptr->echo_level;
    eos_eps_max = eos->get_eps_max();
    max_iter = 100;
    abs_err = 1e-16;
    rel_err = 1e-8;
}

Reconst::~Reconst() {
}


//! This is a shell function for reconstruct T^\mu\nu from T^\tau\mu
void Reconst::ReconstIt_shell(Grid *grid_p, int direc, double tau, double **uq,
                             Grid *grid_pt, int rk_flag) {
    int flag = 0;
    // reconst v and u0 through iteractions
    flag = ReconstIt_velocity_iteration(grid_p, direc, tau, uq, grid_pt,
                                        rk_flag);
    if (flag < 0) {
        // reconstruct from v and u0 failed
        // try to reconst e and rhob through iterations
        flag = ReconstIt(grid_p, direc, tau, uq, grid_pt, rk_flag);
    }

    if (flag == -1) {
        revert_grid(grid_p, grid_pt, rk_flag);
    } else if (flag == -2) {
        if (uq[0][direc]/tau < abs_err) {
            regulate_grid(grid_p, abs_err, 0);
        } else {
            regulate_grid(grid_p, uq[0][direc]/tau, 0);
        }
    }
}


int Reconst::ReconstIt(Grid *grid_p, int direc, double tau, double **uq,
                       Grid *grid_pt, int rk_flag) {
    /* reconstruct TJb from q[0] - q[4] */
    double K00, T00, J0, u[4], epsilon, p, h, rhob;
    double epsilon_prev, rhob_prev, p_prev, p_guess, temperr;
    double epsilon_next, rhob_next, p_next, err, tempf, temph, cs2;
    double eps_guess, scalef;
    int iter, alpha;
    double q[5];
    const double RECONST_PRECISION = 1e-8;

    /* prepare for the iteration */
    /* uq = qiphL, qiphR, etc 
       qiphL[alpha][direc] means, for instance, TJ[alpha][0] 
       in the cell at x+dx/2 calculated from the left 
       */

    /* uq are the conserved charges. That is, the ones appearing in
       d_tau (Ttautau/tau) + d_eta(Ttaueta/tau) + d_perp(Tperptau) = -Tetaeta
       d_tau (Ttaueta) + d_eta(Tetaeta) + d_v(tau Tveta) = -Ttaueta/tau 
       d_tau (Ttauv) + d_eta(Tetav) + d_w(tau Twv) = 0
       d_tau (Jtau) + d_eta Jeta + d_perp(tau Jperp) = 0
       */

    /* q[0] = Ttautau/tau, q[1] = Ttaux, q[2] = Ttauy, q[3] = Ttaueta
       q[4] = Jtau */
    /* uq = qiphL, qiphR, qimhL, qimhR, qirk */


    for (alpha=0; alpha<5; alpha++) {
        q[alpha] = uq[alpha][direc];
    }

    K00 = q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    K00 /= (tau*tau);
 
    T00 = q[0]/tau;
    J0 = q[4]/tau;
 
    if ((T00 < abs_err) || ((T00 - K00/T00) < 0.0)) {
        // can't make Tmunu with this. restore the previous value
        // remember that uq are eigher halfway cells or the final q_next
        // at this point, the original values in grid_pt->TJb are not touched.
        return(-2);
    } /* if t00-k00/t00 < 0.0 */

    /* Iteration scheme */
 
    double eps_init = grid_pt->epsilon;
    double rhob_init = grid_pt->rhob;
    //cs2 = eos->p_e_func(eps_init, rhob_init);
    cs2 = eos->get_cs2(eps_init, rhob_init);
    eps_guess = GuessEps(T00, K00, cs2);
    epsilon_next = eps_guess;
     
    if (isnan(epsilon_next)) {
        music_message << "problem " << eps_guess << " T00=" << T00
                      << " K00=" << K00 << " cs2=" << cs2
                      << " q[0]=" << q[0] << " uq[0][" << direc << "]="
                      << uq[0][direc] << " q[1]=" << q[1] << " q[2]=" << q[2];
        music_message.flush("warning");
    }
    p_guess = eos->get_pressure(epsilon_next, rhob_init);
    p_next = p_guess;

    /* rhob = J0*sqrt( (eps + p)/(T00 + p) ) */

    if (J0 == 0.0) {
        rhob_next = 0.0;
    } else {
        rhob_next = J0*sqrt((eps_guess + p_guess)/(T00 + p_guess));
        if (rhob_next > 1e20 || rhob_next < 0) {
            rhob_next = 0.0;
        }
    }
    err = 1.0;
    for (iter=0; iter<100; iter++) {
        if (err < (RECONST_PRECISION)*0.01) {
            if (isnan(epsilon_next)) {
                cout << "problem2" << endl;
            }
            p_next = eos->get_pressure(epsilon_next, rhob_next);
            break;
        } else {
            epsilon_prev = epsilon_next;
            rhob_prev = rhob_next;
        }
   
        if (isnan(epsilon_prev)) {
            cout << "problem3" << endl;
        }

        p_prev = eos->get_pressure(epsilon_prev, rhob_prev);
        epsilon_next = T00 - K00/(T00 + p_prev);
        err = 0.0;
   
        rhob_next = J0*sqrt((epsilon_prev + p_prev)/(T00 + p_prev));
        temperr = fabs((rhob_next-rhob_prev)/(rhob_prev+abs_err));
        if (temperr < 1e20)
            err += temperr;
        else if (temperr > 1e20)
            err += 1000.0; /* big enough */
        else if (isnan(temperr))
            err += 1000.0;
   
        temperr = fabs((epsilon_next-epsilon_prev)/(epsilon_prev+abs_err));
        if (temperr < 1e20)
            err += temperr;
        else if (temperr > 1e20)
            err += 1000.0; /* big enough */
        else if (isnan(temperr))
            err += 1000.0;
    }/* iter */ 

    if (iter == 100 && grid_p->epsilon > 1e-6) {
        music_message << "Reconst didn't converge. "
                      << "e = " << grid_p->epsilon << " 1/fm^4.";
        music_message.flush("warning");
        return(-1);
    } /* if iteration is unsuccessful, revert */

    /* update */
 
    epsilon = grid_p->epsilon = epsilon_next;
    p = grid_p->p = p_next;
    rhob = grid_p->rhob = rhob_next;
    h = p+epsilon;

    /* q[0] = Ttautau/tau, q[1] = Ttaux, q[2] = Ttauy, q[3] = Ttaueta,
       q[4] = Jtau */

    u[0] = sqrt((q[0]/tau + p)/h);

    double check_u0_var = (fabs(u[0] - grid_pt->u[rk_flag][0])
                           /(grid_pt->u[rk_flag][0]));
    if (check_u0_var > 100.) {
        if (grid_pt->epsilon > 1e-6 && echo_level > 5) {
            music_message << "u0 varies more than 100 times compared to "
                          << "its value at previous time step";
            music_message.flush("warning");
            music_message << "e = " << grid_pt->epsilon*hbarc
                          << " GeV/fm^3, u[0] = " << u[0]
                          << ", prev_u[0] = " << grid_pt->u[rk_flag][0];
            music_message.flush("warning");
        }
        return(-1);
    }

    if (epsilon > eos_eps_max) {
        if (echo_level > 5) {
            music_message << "Reconst velocity: e = " << epsilon
                          << " > e_max in the EoS table.";
            music_message.flush("warning");
	        music_message << "e_max = " << eos_eps_max << " [1/fm^4]";
            music_message.flush("warning");
	        music_message << "previous epsilon = " << grid_pt->epsilon
                          << " [1/fm^4]";
            music_message.flush("warning");
        }
        return(-1);
    }

    //remove if for speed
    if (u[0] > 1e20) {
        u[0] = 1.0;
        u[1] = 0.0;
        u[2] = 0.0;
        u[3] = 0.0;
    } else {
        u[1] = q[1]/tau/h/u[0]; 
        u[2] = q[2]/tau/h/u[0]; 
        u[3] = q[3]/tau/h/u[0]; 
    }

    double u_max = 242582597.70489514; // cosh(20)
    if (u[0] > u_max) {
        fprintf(stderr, "Reconst: u[0] = %e is too large.\n", u[0]);
        if(grid_pt->epsilon > 0.3) {
	        fprintf(stderr, "Reconst: u[0] = %e is too large.\n", u[0]);
	        fprintf(stderr, "epsilon = %e\n", grid_pt->epsilon);
	    }
        return(-1);
    }/* if u[0] is too large, revert */

    /* Correcting normalization of 4-velocity */
    temph = u[0]*u[0] - u[1]*u[1] - u[2]*u[2] - u[3]*u[3];
    // Correct velocity when unitarity is not satisfied to numerical accuracy
    if (fabs(temph - 1.0) > abs_err) {
        // If the deviation is too large, exit MUSIC
        if (fabs(temph - 1.0) > 0.1) {
            fprintf(stderr, "In Reconst, reconstructed : u2 = %e\n", temph);
            fprintf(stderr, "Can't happen.\n");
            exit(0);
        } else if(fabs(temph - 1.0) > sqrt(abs_err) && echo_level > 5) {
            // Warn only when the deviation from 1 is relatively large
            fprintf(stderr, "In Reconst, reconstructed : u2 = %e\n", temph);
            fprintf(stderr, "with u[0] = %e\n", u[0]);
            fprintf(stderr, "Correcting it...\n");
        }   

        // Rescaling spatial components of velocity so that unitarity 
        // is exactly satisfied (u[0] is not modified)
        scalef = (u[0]-1.0);
        scalef *= (u[0]+1.0);
        scalef /= (u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
        scalef = sqrt(scalef);
        u[1] *= scalef;
        u[2] *= scalef;
        u[3] *= scalef;
    }/* if u^mu u_\mu != 1 */
    /* End: Correcting normalization of 4-velocity */

    for (int mu=0; mu<4; mu++) {
        tempf = grid_p->TJb[0][4][mu] = rhob*u[mu];
        tempf = grid_p->u[0][mu] = u[mu];
   
        for (int nu=0; nu<4; nu++) {
            tempf = grid_p->TJb[0][nu][mu] 
                  = ((epsilon + p)*u[nu]*u[mu] + p*(DATA_ptr->gmunu)[nu][mu]);
            if (tempf > 1e20) {
                fprintf(stderr, "Update: TJb[0][%d][%d] is %e.\n",
                        nu, mu, grid_p->TJb[0][nu][mu]);
                fprintf(stderr, "Update: epsilon is %e.\n", epsilon);
                exit(0);
            }
        }/* nu */
    }/* mu */
    return(1); /* on successful execution */
}/* Reconst */

//! reconstruct TJb from q[0] - q[4]
//! reconstruct velocity first for finite mu_B case
//! use iteration to solve v and u0
int Reconst::ReconstIt_velocity_iteration(
    Grid *grid_p, int direc, double tau, double **uq, Grid *grid_pt,
    int rk_flag) {

    double v_critical = 0.563624;

    /* prepare for the iteration */
    /* uq = qiphL, qiphR, etc 
       qiphL[alpha][direc] means, for instance, TJ[alpha][0] 
       in the cell at x+dx/2 calculated from the left */
    
    /* uq are the conserved charges. That is, the ones appearing in
       d_tau (Ttautau/tau) + d_eta(Ttaueta/tau) + d_perp(Tperptau) = -Tetaeta
       d_tau (Ttaueta) + d_eta(Tetaeta) + d_v(tau Tveta) = -Ttaueta/tau 
       d_tau (Ttauv) + d_eta(Tetav) + d_w(tau Twv) = 0
       d_tau (Jtau) + d_eta Jeta + d_perp(tau Jperp) = 0 */
    
    /* q[0] = Ttautau/tau, q[1] = Ttaux, q[2] = Ttauy, q[3] = Ttaueta
       q[4] = Jtau */
    /* uq = qiphL, qiphR, qimhL, qimhR, qirk */

    double q[5];
    for (int alpha = 0; alpha < 5; alpha++) {
        q[alpha] = uq[alpha][direc]/tau;
    }

    double K00 = q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    double M = sqrt(K00);
    double T00 = q[0];
    double J0 = q[4];

    if ((T00 < abs_err) || (T00 - K00/T00) < 0.0) {
        // can't make Tmunu with this. restore the previous value 
        // remember that uq are eigher halfway cells or the final q_next 
        // at this point, the original values in grid_pt->TJb are not touched. 
        if (echo_level > 9) {
            music_message.warning("Reconst velocity:: can not find solution!");
            music_message << "T00 = " << T00 << ", K00 = " << K00;
            music_message.flush("warning");
        }
        return(-2);
    }/* if t00-k00/t00 < 0.0 */

    double u[4], epsilon, pressure, rhob;

    int v_status = 1;
    int iter = 0;
    double abs_error_v = 10.0;
    double rel_error_v = 10.0;
    double v_next = 1.0;
    double v_prev = 0.0;
    do {
        iter++;
        v_next = reconst_velocity_f(v_prev, T00, M, J0);
        abs_error_v = fabs(v_next - v_prev);
        rel_error_v = 2.*abs_error_v/(v_next + v_prev + 1e-15);
        v_prev = v_next;
        if (iter > max_iter) {
            v_status = 0;
            break;
        }
    } while (abs_error_v > abs_err && rel_error_v > rel_err);

    double v_solution;
    if (v_status == 1) {
        v_solution = v_next;
    } else {
        if (echo_level > 5) {
            music_message.warning("Reconst velocity:: can not find solution!");
            music_message.warning("output the results at the last iteration:");
            music_message.warning("iter  [lower, upper]  root  err(est)");
            music_message << iter << "   [" << v_prev << ",  " << v_next
                          << "]  " << abs_error_v << "  " << rel_error_v;
            music_message.flush("warning");
        }
        return(-1);
    }/* if iteration is unsuccessful, revert */
   
    // for large velocity, solve u0
    double u0_solution = 1.0;
    if (v_solution > v_critical) {
        double u0_prev = 1./sqrt(1. - v_solution*v_solution);
        int u0_status = 1;
        iter = 0;
        double u0_next;
        double abs_error_u0 = 1.0;
        double rel_error_u0 = 1.0;
        do {
            iter++;
            u0_next = reconst_u0_f(u0_prev, T00, K00, M, J0);
            abs_error_u0 = fabs(u0_next - u0_prev);
            rel_error_u0 = 2.*abs_error_u0/(u0_next + u0_prev + 1e-15);
            u0_prev = u0_next;
            if (iter > max_iter) {
                u0_status = 0;
                break;
            }
        } while (abs_error_u0 > abs_err && rel_error_u0 > rel_err);

        if (u0_status == 1) {
            u0_solution = u0_next;
        } else {
            if (echo_level > 5) {
                music_message.warning(
                        "Reconst velocity:: can not find solution!");
                music_message.warning(
                        "output the results at the last iteration:");
                music_message.warning("iter  [lower, upper]  root  err(est)");
                music_message << iter << "   [" << u0_prev << ",  " << u0_next
                              << "]  " << abs_error_u0 << "  " << rel_error_u0;
                music_message.flush("warning");
            }
            return(-1);
        } /* if iteration is unsuccessful, revert */
    }

    // successfully found velocity, now update everything else
    if (v_solution < v_critical) {
        u[0] = 1./(sqrt(1. - v_solution*v_solution) + v_solution*abs_err);
        epsilon = T00 - v_solution*sqrt(K00);
        rhob = J0/u[0];
    } else {
        u[0] = u0_solution;
        epsilon = T00 - sqrt((1. - 1./(u0_solution*u0_solution))*K00);
        rhob = J0/u0_solution;
    }

    double check_u0_var = (fabs(u[0] - grid_pt->u[rk_flag][0])
                           /(grid_pt->u[rk_flag][0]));
    if (check_u0_var > 100.) {
        if (grid_pt->epsilon > 1e-6 || echo_level > 5) {
            music_message << "u0 varies more than 100 times compared to "
                          << "its value at previous time step";
            music_message.flush("warning");
            music_message << "e = " << grid_pt->epsilon
                          << ", u[0] = " << u[0]
                          << ", prev_u[0] = " << grid_pt->u[rk_flag][0];
            music_message.flush("warning");
        }
        return(-1);
    }

    if (epsilon > eos_eps_max) {
        if (echo_level > 5) {
            music_message << "Reconst velocity: e = " << epsilon
                          << " > e_max in the EoS table.";
            music_message.flush("warning");
	        music_message << "e_max = " << eos_eps_max << " [1/fm^4]";
            music_message.flush("warning");
	        music_message << "previous epsilon = " << grid_pt->epsilon
                          << " [1/fm^4]";
            music_message.flush("warning");
        }
        return(-1);
    }

    grid_p->epsilon = epsilon;
    grid_p->rhob = rhob;

    pressure = eos->get_pressure(epsilon, rhob);
    grid_p->p = pressure;

    double u_max = 242582597.70489514; // cosh(20)
    //remove if for speed
    if (u[0] > 1e20) {
        // check whether velocity is too large
        if (echo_level > 5) {
            music_message << "Reconst velocity: u[0] = " << u[0]
                          << " is crazy!"
                          << "epsilon = " << grid_pt->epsilon;
            music_message.flush("warning");
        }
        return(-1);
    } else if(u[0] > u_max) {
        // check whether velocity is too large
        if (echo_level > 5) {
            music_message << "Reconst velocity: u[0] = " << u[0]
                          << " is too large!"
                          << "epsilon = " << grid_pt->epsilon;
            music_message.flush("warning");
        }
        return(-1);
    } else {
        // individual components of velocity
        double velocity_inverse_factor = u[0]/(T00 + pressure);
        u[1] = q[1]*velocity_inverse_factor; 
        u[2] = q[2]*velocity_inverse_factor; 
        u[3] = q[3]*velocity_inverse_factor; 
    }

    // Correcting normalization of 4-velocity
    double temp_usq = u[0]*u[0] - u[1]*u[1] - u[2]*u[2] - u[3]*u[3];
    // Correct velocity when unitarity is not satisfied to numerical accuracy
    if (fabs(temp_usq - 1.0) > abs_err) {
        // If the deviation is too large, exit MUSIC
        if (fabs(temp_usq - 1.0) > 0.1*u[0]) {
            music_message << "In Reconst velocity, reconstructed: u^2 - 1 = "
                          << temp_usq - 1.0;
            music_message.flush("error");
            music_message << "u[0]=" << u[0] << ", u[1]=" << u[1]
                          << ", u[2]=" << u[2] << ", u[3]=" << u[3];
            music_message.flush("error");
            music_message << "e=" << epsilon << ", rhob=" << rhob
                          << ", p=" << pressure;
            music_message.flush("error");
            music_message << "with q1 = " << q[1] << ", q2 = " << q[2]
                          << ", q3 = " << q[3];
            music_message.flush("error");
            exit(0);
        } else if (fabs(temp_usq - 1.0) > sqrt(abs_err)*u[0]
                   && echo_level > 5) {
            // Warn only when the deviation from 1 is relatively large
            music_message << "In Reconst velocity, reconstructed: u^2 - 1 = "
                          << temp_usq - 1.0;
            music_message.flush("warning");
            double f_res;
            if (v_solution < v_critical)
                f_res = fabs(v_solution
                             - reconst_velocity_f(v_solution, T00, M, J0));
            else
                f_res = fabs(u0_solution
                             - reconst_u0_f(u0_solution, T00, K00, M, J0));
            music_message << "with v = " << v_solution << ", u[0] = " << u[0]
                          << ", res = " << f_res;
            music_message.flush("warning");
            music_message << "with u[1] = " << u[1]
                          << "with u[2] = " << u[2]
                          << "with u[3] = " << u[3];
            music_message.flush("warning");
            music_message << "with T00 = " << T00 << ", K = " << K00;
            music_message.flush("warning");
            music_message << "with q1 = " << q[1] << ", q2 = " << q[2]
                          << ", q3 = " << q[3];
            music_message.flush("warning");
            music_message.warning("Correcting it...");
        }
        // Rescaling spatial components of velocity so that unitarity 
        // is exactly satisfied (u[0] is not modified)
        double scalef = sqrt((u[0]*u[0] - 1.0)
                             /(u[1]*u[1] + u[2]*u[2] + u[3]*u[3] + abs_err));
        u[1] *= scalef;
        u[2] *= scalef;
        u[3] *= scalef;
    } // if u^mu u_\mu != 1 
    // End: Correcting normalization of 4-velocity
   
    for (int mu = 0; mu < 4; mu++) {
        grid_p->TJb[0][4][mu] = rhob*u[mu];
        grid_p->u[0][mu] = u[mu];
        double gfac = (mu == 0 ? -1.0 : 1.0);
        double tempf = (epsilon + pressure)*u[mu]*u[mu] + pressure*gfac;
        grid_p->TJb[0][mu][mu] = tempf;
    }
    for (int mu = 0; mu < 4; mu++) {
        double tempf = (epsilon + pressure)*u[mu];
        for (int nu = mu+1; nu < 4; nu++) {
            double Tmunu = tempf*u[nu];
            grid_p->TJb[0][mu][nu] = Tmunu;
            grid_p->TJb[0][nu][mu] = Tmunu;
        }
    }
    return(1);
}


double Reconst::GuessEps(double T00, double K00, double cs2) {
    double f;
 
    if (cs2 < abs_err) {
        f = ((-K00 + T00*T00)
             *(cs2*K00*T00*T00 + pow(T00, 4) 
               + cs2*cs2*K00*(2*K00 - T00*T00)))
            /pow(T00, 5);
    } else {
        f = ((-1.0 + cs2)*T00 
            + sqrt(-4.0*cs2*K00 + T00*T00
                   + 2.0*cs2*T00*T00
                   + cs2*cs2*T00*T00))/(2.0*cs2);
    }
    return f;
}/*  GuessEps */


//! this function returns f(v) = M/(M0 + P)
double Reconst::reconst_velocity_f(double v, double T00, double M,
                                   double J0) {
    double epsilon = T00 - v*M;
    double rho = J0*sqrt(1 - v*v);
   
    double pressure = eos->get_pressure(epsilon, rho);
    double fv = M/(T00 + pressure);
    return(fv);
}


//! this function returns f(u0) = (M0+P)/sqrt((M0+P)^2 - M^2)
double Reconst::reconst_u0_f(double u0, double T00, double K00, double M,
                             double J0) {
    double epsilon = T00 - sqrt(1. - 1./u0/u0)*M;
    double rho = J0/u0;
    
    double pressure = eos->get_pressure(epsilon, rho);
    double fu = (T00 + pressure)/sqrt((T00 + pressure)*(T00 + pressure) - K00);
    return(fu);
}


//! This function reverts the grid information back its values
//! at the previous time step
void Reconst::revert_grid(Grid *grid_current, Grid *grid_prev, int rk_flag) {
    if (echo_level > 5) {
        music_message.warning("revert grid to its previous value ...");
    }
    grid_current->epsilon = grid_prev->epsilon;
    grid_current->rhob = grid_prev->rhob;
    grid_current->p = grid_prev->p;
    for (int mu = 0; mu < 4; mu++) {
        grid_current->TJb[0][4][mu] = grid_prev->TJb[rk_flag][4][mu];
        grid_current->u[0][mu] = grid_prev->u[rk_flag][mu];
 
        for (int nu = 0; nu < 4; nu++) {
            grid_current->TJb[0][nu][mu] = grid_prev->TJb[rk_flag][nu][mu];
        }
    }
}


//! This function regulate the grid information
void Reconst::regulate_grid(Grid *grid_cell, double elocal, int rk_flag) {
    if (echo_level > 9) {
        music_message.warning("regulate grid to ideal form...");
    }
    grid_cell->epsilon = elocal;
    grid_cell->rhob = 0.0;
    double plocal = eos->get_pressure(elocal, 0.0);
    grid_cell->p = plocal;

    grid_cell->u[rk_flag][0] = 1.0;
    grid_cell->u[rk_flag][1] = 0.0;
    grid_cell->u[rk_flag][2] = 0.0;
    grid_cell->u[rk_flag][3] = 0.0;
    
    grid_cell->TJb[rk_flag][0][0] = elocal;
    grid_cell->TJb[rk_flag][1][1] = plocal;
    grid_cell->TJb[rk_flag][2][2] = plocal;
    grid_cell->TJb[rk_flag][3][3] = plocal;

    for (int mu = 0; mu < 4; mu++) {
        grid_cell->TJb[rk_flag][4][mu] = 0.0;
        for (int nu = mu+1; nu < 4; nu++) {
            grid_cell->TJb[rk_flag][mu][nu] = 0.0;
            grid_cell->TJb[rk_flag][nu][mu] = 0.0;
        }
    }
}
