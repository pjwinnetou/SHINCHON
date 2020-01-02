// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

#ifndef MINMOD_H
#define MINMOD_H

#include "data.h"

//! This is a class define the flux limiter
class Minmod
{
    private:
        double theta_flux;

    public:
        Minmod(InitData* DATA);
        ~Minmod();
        double minmod_dx(double up1, double u, double um1);
        double minmod_theta_dx(double up1, double u, double um1, double theta);
};
#endif
