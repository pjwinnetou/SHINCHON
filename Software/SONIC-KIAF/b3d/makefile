include madai_repo/makefile_defs.mk

HPM = -I${MADAI_INSTALLDIR}/include/madai -I${MADAI_GSLPATH}/include -I${MADAI_HDF5_HOME}/include -L${MADAI_GSLPATH}/lib -L${MADAI_INSTALLDIR}/lib/madai -L${MADAI_HDF5_HOME}/lib
PROGDATA=${MADAI_INSTALLDIR}/progdata/madai

all : b3d analyze

alllibs : corallibs rhicutilslibs b3dlibs

anallibs :
	make -C ${MADAI_REPO}/rhic/anal install

b3dlibs :
	make -C ${MADAI_REPO}/rhic/b3d install

corallibs :
	make -C ${MADAI_REPO}/rhic/coral install

rhicutilslibs :
	make -C ${MADAI_REPO}/rhic/rhicutils install

b3d : progdata_link b3dmain.cc b3dlibs corallibs rhicutilslibs
	${MADAI_CPP} -o b3d b3dmain.cc ${MADAI_CFLAGS} ${HPM} -fopenmp -lb3d -lcoralutils -lrhicutils -lgsl -lgslcblas -lhdf5_hl_cpp  -lhdf5_cpp -lhdf5_hl -lhdf5 -lz -lpthread

analyze : progdata_link analmain.cc anallibs corallibs rhicutilslibs
	${MADAI_CPP} -o analyze analmain.cc ${MADAI_CFLAGS} ${HPM} -fopenmp -lanalyze -lb3d -lcoral -lcoralutils -lrhicutils  -lgsl -lgslcblas -lhdf5_hl_cpp  -lhdf5_cpp -lhdf5_hl -lhdf5 -lz -lpthread

progdata_link :
	ln -s -f ${PROGDATA} .

clean :
	rm -f b3d analyze
	rm -f ${MADAI_REPO}/rhic/anal/build/*
	rm -f ${MADAI_REPO}/rhic/b3d/build/*
	rm -f ${MADAI_REPO}/rhic/coral/build/coral/* 
	rm -f ${MADAI_REPO}/rhic/coral/build/coralutils/* 
	rm -f ${MADAI_REPO}/rhic/coral/build/xgraph/* 
	rm -f ${MADAI_REPO}/rhic/rhicutils/build/*
