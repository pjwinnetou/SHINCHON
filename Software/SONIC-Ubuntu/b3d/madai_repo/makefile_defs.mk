#B3DMAINDIR = /scratch/public/b3d-2.3
B3DMAINDIR = $(PWD)
MADAI_REPO = $(B3DMAINDIR)/madai_repo
#prefix of .../rhic directory (root of madai installation)
MADAI_GSLPATH = /usr/local/include/gsl
#where gnu scientific library is installed
MADAI_X11PATH = /usr
#root of X11 installation
MADAI_HDF5_HOME = /usr/local
#root of HDF5 installation
MADAI_CPP = $(shell which g++)
#compiler
MADAI_INSTALLDIR = $(B3DMAINDIR)/install
#location of where you want things installed
MADAI_CFLAGS = -O2
#MADAI_CFLAGS = -O
#compiler optimization flags, usually -O2 for linux, -fast for OSX with g++
