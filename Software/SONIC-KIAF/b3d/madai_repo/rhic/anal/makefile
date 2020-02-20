include ../../makefile_defs.mk
# enter make all, make b3d, make analyze, make install, make clean or make uninstall
### the following variables are set in the include statement above
CPP=${MADAI_CPP}
OPT=${MADAI_CFLAGS}
#e.g. /Users/scottepratt/git -- 
CORAL_HOME=${MADAI_REPO}/rhic/coral
#########################################################################
#HOPEFULLY, everything below need not be touched

all : analyze

#you must be sure coral and b3d are installed first
#coral :
#	make -C ${CORAL_HOME} -f ${CORAL_HOME}/makefile;
#b3d :
# make -C ${MADAI_REPO}/rhic/b3d ${MADAI_REPO}/rhic/b3d/makefile

analyze: analyze_dirs lib/libanalyze.a ${ANALYZE_HFILES}

install : ${MADAI_INSTALLDIR}/lib/madai/libanalyze.a ${MADAI_INSTALLDIR}/include/madai/analyze.h

clean :
	rm -f lib/libanalyze.a build/*.o include/*.h

INC=-I include -I ${MADAI_GSLPATH}/include -I ${MADAI_INSTALLDIR}/include/madai -I ${MADAI_HDF5_HOME}/include

uninstall :
	rm -f ${MADAI_INSTALLDIR}/lib/madai/libanalyze.a ${MADAI_INSTALLDIR}/include/madai/analyze.h


ANALYZE_HFILES = include/analyze.h

ANALYZE_OBJFILES = build/analyze.o\
build/anal_balance.o\
build/anal_gamma.o\
build/anal_PHENIXspectra.o\
build/anal_STARspectra.o\
build/anal_v2.o\
build/anal_qn.o\
build/anal_hbt.o\
build/distill.o\
build/consolidate.o

#######################

analyze_dirs :
	mkdir -p include;\
	mkdir -p lib;\
	mkdir -p ${MADAI_INSTALLDIR}/include;\
	mkdir -p ${MADAI_INSTALLDIR}/lib;\

#######################

include/analyze.h : src/analyze.h
	cp -f src/analyze.h include/

######################

lib/libanalyze.a : ${ANALYZE_HFILES} ${ANALYZE_OBJFILES}
	rm -f lib/libanalyze.a;\
	ar -ru lib/libanalyze.a ${ANALYZE_OBJFILES}

build/analyze.o : src/analyze.cc ${ANALYZE_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/analyze.o src/analyze.cc

build/anal_balance.o : src/anal_balance.cc ${ANALYZE_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/anal_balance.o src/anal_balance.cc

build/anal_gamma.o : src/anal_gamma.cc ${ANALYZE_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/anal_gamma.o src/anal_gamma.cc

build/anal_PHENIXspectra.o : src/anal_PHENIXspectra.cc ${ANALYZE_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/anal_PHENIXspectra.o src/anal_PHENIXspectra.cc
	
build/anal_STARspectra.o : src/anal_STARspectra.cc ${ANALYZE_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/anal_STARspectra.o src/anal_STARspectra.cc

build/anal_v2.o : src/anal_v2.cc ${ANALYZE_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/anal_v2.o src/anal_v2.cc

build/anal_qn.o : src/anal_qn.cc ${ANALYZE_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/anal_qn.o src/anal_qn.cc

build/anal_hbt.o : src/anal_hbt.cc ${ANALYZE_HFILES}
	${CPP} -c -fopenmp ${OPT} ${INC} -o build/anal_hbt.o src/anal_hbt.cc

build/distill.o : src/distill.cc ${ANALYZE_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/distill.o src/distill.cc

build/consolidate.o : src/consolidate.cc ${ANALYZE_HFILES}
	${CPP} -c ${OPT} ${INC} -o build/consolidate.o src/consolidate.cc

######################

${MADAI_INSTALLDIR}/include/madai/analyze.h : include/analyze.h
	cp -f include/analyze.h ${MADAI_INSTALLDIR}/include/madai/

${MADAI_INSTALLDIR}/lib/madai/libanalyze.a : lib/libanalyze.a
	cp -f lib/libanalyze.a ${MADAI_INSTALLDIR}/lib/madai/

