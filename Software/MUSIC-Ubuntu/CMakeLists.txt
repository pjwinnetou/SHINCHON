cmake_minimum_required(VERSION 2.8)
project(mpihydro)

#for MacOSX
cmake_policy(SET CMP0042 NEW)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")

# Require MPI and GSL for this project:
find_package(MPI REQUIRED)
find_package(GSL REQUIRED)

set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")
set(CMAKE_CXX_FLAGS " -O3 -Wall")

include_directories(${MPI_INCLUDE_PATH} ${GSL_INCLUDE_DIR})
add_subdirectory (src)  
