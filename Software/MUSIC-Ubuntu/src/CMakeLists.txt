set (SOURCES
    main.cpp
    grid.cpp
    eos.cpp
    evolve.cpp
    init.cpp
    reconst.cpp
    freeze.cpp
    freeze_pseudo.cpp
    minmod.cpp
    glauber.cpp
    advance.cpp
    u_derivative.cpp
    dissipative.cpp
    util.cpp
    grid_info.cpp
    read_in_parameters.cpp
    reso_decay.cpp
    music.cpp
    pretty_ostream.cpp
    HydroinfoMUSIC.cpp
    )

include_directories(${MPI_INCLUDE_PATH})
add_library(music_lib SHARED ${SOURCES})
add_executable(mpihydro ${SOURCES})
if (APPLE)
    set_target_properties (mpihydro PROPERTIES COMPILE_FLAGS "-DAPPLE")
    set_target_properties (music_lib PROPERTIES COMPILE_FLAGS "-DAPPLE")
endif (APPLE)
target_link_libraries (music_lib ${MPI_LIBRARIES} ${GSL_LIBRARIES})
target_link_libraries (mpihydro ${MPI_LIBRARIES} ${GSL_LIBRARIES})

install(TARGETS mpihydro DESTINATION ${CMAKE_HOME_DIRECTORY})
