// MUSIC - a 3+1D viscous relativistic hydrodynamic code for heavy ion collisions
// Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen

#include "./music.h"
#include <cstdlib>
#include <ctime>
#include <vector>

using namespace std;

MUSIC::MUSIC(int argc, char *argv[]) {
    // you have the option to give a second command line option,
    // which is an integer to be added to the random seed from the current time
    // because on the cluster it will happen that two jobs start
    // at exactly the same time, this makes sure that they dont run 
    // with excatly the same seed
    welcome_message();
    string sseed;
    int seed = 0;
    if (argc > 2) {
        sseed = argv[2];
        seed = atoi(sseed.c_str());
    }
    seed *= 10000;
    DATA.seed = seed;
    
    if (argc >1)
        input_file = *(argv+1);
    else
        input_file = "";

    reader.read_in_parameters(&DATA, input_file);

    music_message.info("Initialize MPI ... ");

    MPI_Init(NULL, NULL);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);       // number of current processor
    MPI_Comm_size(MPI_COMM_WORLD, &size);       // total number of processors

    DATA.rank = rank;
    DATA.size = size;

    music_message << "This is processor " << (rank+1) << "/" << size
                  << ": READY.";
    music_message.flush("info");
    
    if (DATA.neta%size != 0 && DATA.mode < 3) {
        music_message << " Number of cells in eta direction " << DATA.neta
                      << " is not a multiple of the number of processors "
                      << size << ". Exiting.";
        if (DATA.boost_invariant == 1 && DATA.neta < size) {
            music_message << "Simulation is set to boost_invariant mode. "
                          << "Please set size to 1";
        }
        music_message.flush("error");
        exit(1);
    }
    
    // reduce the lattice size on each processor
    // (this is slicing it in 'size' pieces)
    DATA.neta = DATA.neta/size;
    
    eos = new EOS(&DATA);
    mode = DATA.mode;
    init = NULL;
    glauber = NULL;
    evolve = NULL;
    freeze = NULL;
    hydro_info_ptr = new HydroinfoMUSIC;
}

MUSIC::~MUSIC() {
    delete eos;
    if (glauber != NULL) {
        delete glauber;
    }
    if (init != NULL) {
        delete init;
    }
    if (evolve != NULL) {
        delete evolve;
    }
    if (freeze != NULL) {
        delete freeze;
    }
    if (hydro_info_ptr != NULL) {
        delete hydro_info_ptr;
    }
    MPI_Finalize();
}


int MUSIC::initialize_hydro() {
    size = DATA.size;
    rank = DATA.rank;
    ostringstream info_message;
    music_message << "size=" << size << ", rank=" << rank;
    music_message.flush("info");

    int status = 0;
    stringstream ss;
    ss << "bash -c 'rm surface.dat surface{0.." << size-1 << "}.dat'";
    status = system(ss.str().c_str());
    if (DATA.Initial_profile == 1 || DATA.Initial_profile == 3) {
        music_message.info("init Glauber");
        if (glauber != NULL) {
            delete glauber;
        }
        glauber = new Glauber(&DATA);
        glauber->initGlauber(DATA.SigmaNN, DATA.Target, DATA.Projectile,
                             DATA.b, size, rank);
    }
    if (init != NULL) {
        delete init;
    }
    init = new Init(eos, glauber);
    init->InitArena(&DATA, &arena, &Lneighbor, &Rneighbor, size, rank);
    
    return(status);
}


int MUSIC::initialize_hydro_from_vector(std::vector<double> entropy_density,
                                        double dx) {
    size = DATA.size;
    rank = DATA.rank;
    DATA.Initial_profile = 41;
    ostringstream info_message;
    music_message << "size=" << size << ", rank=" << rank;
    music_message.flush("info");

    int status = 0;
    stringstream ss;
    ss << "bash -c 'rm surface.dat surface{0.." << size-1 << "}.dat'";
    status = system(ss.str().c_str());
    if (init != NULL) {
        delete init;
    }
    init = new Init(eos, glauber);
    DATA.delta_x = dx;
    DATA.delta_y = dx;
    init->get_entropy_density_vector(entropy_density);
    init->InitArena(&DATA, &arena, &Lneighbor, &Rneighbor, size, rank);
    return(status);
}


int MUSIC::initialize_hydro_from_pre_equilibrium_vectors(
        const double dx, const double dz, const double z_max, const int nz,
        std::vector<double> e_in,
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

    size = DATA.size;
    rank = DATA.rank;
    DATA.Initial_profile = 42;

    if (nz > 1) {
        DATA.boost_invariant = false;
        DATA.eta_size        = z_max;
        DATA.delta_eta       = dz;
        DATA.neta            = nz;
    }

    ostringstream info_message;
    music_message << "size=" << size << ", rank=" << rank;
    music_message.flush("info");

    int status = 0;
    stringstream ss;
    ss << "bash -c 'rm surface.dat surface{0.." << size-1 << "}.dat'";
    status = system(ss.str().c_str());
    if (init != NULL) {
        delete init;
    }
    init = new Init(eos, glauber);
    DATA.delta_x = dx;
    DATA.delta_y = dx;
    init->get_pre_equilibrium_vectors(
        e_in, u_tau_in, u_x_in, u_y_in, u_eta_in,
        pi_00_in, pi_01_in, pi_02_in, pi_03_in, pi_11_in, pi_12_in, pi_13_in,
        pi_22_in, pi_23_in, pi_33_in, Bulk_pi_in);
    init->InitArena(&DATA, &arena, &Lneighbor, &Rneighbor, size, rank);
    return(status);
}


void MUSIC::get_hydro_info(double x, double y, double z, double t,
                           fluidCell* fluid_cell_info) {
    if (DATA.store_hydro_info_in_memory == 0) {
        music_message << "hydro evolution informaiton is not stored "
                      << "in the momeory! Please set the parameter "
                      << "store_hydro_info_in_memory to 1~";
        music_message.flush("error");
        exit(1);
    }
    hydro_info_ptr->getHydroValues(x, y, z, t, fluid_cell_info);
}

void MUSIC::clear_hydro_info_from_memory() {
    if (DATA.store_hydro_info_in_memory == 0) {
        music_message << "The parameter store_hydro_info_in_memory is 0. "
                      << "No need to clean memory~";
        music_message.flush("warning");
    } else {
        hydro_info_ptr->clean_hydro_event();
    }
}

int MUSIC::run_hydro() {
    if (evolve != NULL) {
        delete evolve;
    }

    evolve = new Evolve(eos, &DATA);
    evolve->EvolveIt(&DATA, arena, Lneighbor, Rneighbor, size, rank,
                     hydro_info_ptr);
    MPI_Barrier(MPI_COMM_WORLD);

    int status = 0;
    // combining surface files from different rank
    stringstream act;
	act << "bash -c 'cat surface{0.." << size-1 << "}.dat > surface.dat'";
	status = system(act.str().c_str());
	  
	//remove extraneous files
	stringstream ass;
	ass << "bash -c 'rm surface{0.." << size-1 << "}.dat'";
	status = system(ass.str().c_str());
    return(status);
}

int MUSIC::run_Cooper_Frye(int running_mode) {
    if (freeze != NULL) {
        delete freeze;
    }
    freeze = new Freeze(&DATA);
    freeze->CooperFrye_pseudo(DATA.particleSpectrumNumber, running_mode,
                              &DATA, eos, size, rank);
    int status = 0;
    if (running_mode == 3) {
	    //remove extraneous files
	    stringstream ass;
	    ass << "bash -c 'rm yptphiSpectra{0.." << size-1 << "}.dat'";
	    status = system(ass.str().c_str());
    }
    return(status);
}

void MUSIC::welcome_message() {
    srand (time(NULL));
    display_logo(rand()%4);
    display_code_description_and_copyright();
}

void MUSIC::display_logo(int selector) {
    switch (selector) {
        case 0:  // 3D Diagonal
            cout << "================================================================" << endl;
            cout << "|           ____                                               |" << endl;
            cout << "|         ,'  , `.               .--.--.      ,---,  ,----..   |" << endl;
            cout << "|      ,-+-,.' _ |         ,--, /  /    '. ,`--.' | /   /   \\  |" << endl;
            cout << "|   ,-+-. ;   , ||       ,'_ /||  :  /`. / |   :  :|   :     : |" << endl;
            cout << "|  ,--.'|'   |  ;|  .--. |  | :;  |  |--`  :   |  '.   |  ;. / |" << endl;
            cout << "| |   |  ,', |  ':,'_ /| :  . ||  :  ;_    |   :  |.   ; /--`  |" << endl;
            cout << "| |   | /  | |  |||  ' | |  . . \\  \\    `. '   '  ;;   | ;     |" << endl;
            cout << "| '   | :  | :  |,|  | ' |  | |  `----.   \\|   |  ||   : |     |" << endl;
            cout << "| ;   . |  ; |--' :  | | :  ' ;  __ \\  \\  |'   :  ;.   | '___  |" << endl;
            cout << "| |   : |  | ,    |  ; ' |  | ' /  /`--'  /|   |  ''   ; : .'| |" << endl;
            cout << "| |   : '  |/     :  | : ;  ; |'--'.     / '   :  |'   | '/  : |" << endl;
            cout << "| ;   | |`-'      '  :  `--'   \\ `--'---'  ;   |.' |   :    /  |" << endl;
            cout << "| |   ;/          :  ,      .-./           '---'    \\   \\ .'   |" << endl;
            cout << "| '---'            `--`----'                         `---`     |" << endl;
            cout << "================================================================" << endl;
            break;
        case 1:  // bloody
            cout << "==============================================" << endl;
            cout << "|  ███▄ ▄███▓ █    ██   ██████  ██▓ ▄████▄   |" << endl;
            cout << "| ▓██▒▀█▀ ██▒ ██  ▓██▒▒██    ▒ ▓██▒▒██▀ ▀█   |" << endl;
            cout << "| ▓██    ▓██░▓██  ▒██░░ ▓██▄   ▒██▒▒▓█    ▄  |" << endl;
            cout << "| ▒██    ▒██ ▓▓█  ░██░  ▒   ██▒░██░▒▓▓▄ ▄██▒ |" << endl;
            cout << "| ▒██▒   ░██▒▒▒█████▓ ▒██████▒▒░██░▒ ▓███▀ ░ |" << endl;
            cout << "| ░ ▒░   ░  ░░▒▓▒ ▒ ▒ ▒ ▒▓▒ ▒ ░░▓  ░ ░▒ ▒  ░ |" << endl;
            cout << "| ░  ░      ░░░▒░ ░ ░ ░ ░▒  ░ ░ ▒ ░  ░  ▒    |" << endl;
            cout << "| ░      ░    ░░░ ░ ░ ░  ░  ░   ▒ ░░         |" << endl;
            cout << "|        ░      ░           ░   ░  ░ ░       |" << endl;
            cout << "|                                  ░         |" << endl;
            cout << "==============================================" << endl;
            break;
        case 2:  // Dancing font
            cout << "====================================================" << endl;
            cout << "|   __  __     _   _   ____                   ____  |" << endl;
            cout << "| U|' \\/ '|uU |\"|u| | / __\"| u      ___    U /\"___| |" << endl;
            cout << "| \\| |\\/| |/ \\| |\\| |<\\___ \\/      |_\"_|   \\| | u   |" << endl;
            cout << "|  | |  | |   | |_| | u___) |       | |     | |/__  |" << endl;
            cout << "|  |_|  |_|  <<\\___/  |____/>>    U/| |\\u    \\____| |" << endl;
            cout << "| <<,-,,-.  (__) )(    )(  (__).-,_|___|_,-._// \\\\  |" << endl;
            cout << "|  (./  \\.)     (__)  (__)      \\_)-' '-(_/(__)(__) |" << endl;
            cout << "====================================================" << endl;
            break;
        case 3:  // STAR Wars
            cout << "=====================================================" << endl;
            cout << "| .___  ___.  __    __       _______. __    ______  |" << endl;
            cout << "| |   \\/   | |  |  |  |     /       ||  |  /      | |" << endl;
            cout << "| |  \\  /  | |  |  |  |    |   (----`|  | |  ,----' |" << endl;
            cout << "| |  |\\/|  | |  |  |  |     \\   \\    |  | |  |      |" << endl;
            cout << "| |  |  |  | |  `--'  | .----)   |   |  | |  `----. |" << endl;
            cout << "| |__|  |__|  \\______/  |_______/    |__|  \\______| |" << endl;
            cout << "=====================================================" << endl;
            break;
    }

}


//! This function prints out code desciprtion and copyright information
void MUSIC::display_code_description_and_copyright() {
    cout << "MUSIC - a 3+1D viscous relativistic hydrodynamic code for "
         << "heavy ion collisions" << endl;
    cout << "Copyright (C) 2017  Gabriel Denicol, Charles Gale, Sangyong Jeon, "
         << "Matthew Luzum, Jean-François Paquet, Björn Schenke, Chun Shen"
         << endl;
}
