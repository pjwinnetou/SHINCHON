// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

// Copyright @ 2011 Bjoern Schenke, Sangyong Jeon, and Charles Gale
#include "mpi.h"
#include <stdio.h>
#include <sys/stat.h>  // for mkdir

#include "./music.h"

// main program
int main(int argc, char *argv[]) {
    MUSIC *music_hydro = new MUSIC(argc, argv);

    int running_mode = music_hydro->get_running_mode();
    if (running_mode == 1 || running_mode == 2) {
        music_hydro->initialize_hydro();
        music_hydro->run_hydro();
    }
  
    if (running_mode == 1 || running_mode == 3
            || running_mode == 4 || running_mode >= 5) {
        music_hydro->run_Cooper_Frye(running_mode);
    }

    delete music_hydro;
    return(0);
}/* main */
