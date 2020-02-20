B3DMAINDIR = $(PWD)
#B3DMAINDIR = /scratch/public/b3d-2.3
MADAI_REPO = $(B3DMAINDIR)/madai_repo
#prefix of .../rhic directory (root of madai installation)
MADAI_GSLPATH = /alice/home/shlim/sources/gsl-2.6-gcc4.9.3-install
#where gnu scientific library is installed
MADAI_X11PATH = /usr
#root of X11 installation
MADAI_HDF5_HOME = /alice/home/shlim/sources/hdf5-1.10.6-gcc4.9.3-install
#root of HDF5 installation
MADAI_CPP = $(shell which g++)
#compiler
MADAI_INSTALLDIR = $(B3DMAINDIR)/install
#location of where you want things installed
MADAI_CFLAGS = -O2 -g -std=c++11
#MADAI_CFLAGS = -O 
#compiler optimization flags, usually -O2 for linux, -fast for OSX with g++
