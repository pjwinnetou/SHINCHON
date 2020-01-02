# ===========================================================================
#  Makefile for MUSIC                                 Chun Shen Dec. 12, 2014
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install	make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

CC := mpicxx-openmpi-gcc5 
CFLAGS= -Wall -O3 $(shell gsl-config --cflags)

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(CFLAGS) $(shell gsl-config --libs)
SYSTEMFILES     =       $(SRCGNU)

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	mpihydro
endif

SRC		=	main.cpp grid.cpp eos.cpp evolve.cpp init.cpp reconst.cpp \
            freeze.cpp freeze_pseudo.cpp minmod.cpp glauber.cpp \
            advance.cpp u_derivative.cpp dissipative.cpp \
            util.cpp grid_info.cpp read_in_parameters.cpp music.cpp \
			reso_decay.cpp pretty_ostream.cpp HydroinfoMUSIC.cpp

INC		= 	grid.h eos.h evolve.h init.h reconst.h freeze.h \
            minmod.h glauber.h advance.h u_derivative.h dissipative.h \
            util.h reconst.h int.h data.h grid_info.h \
			read_in_parameters.h music.h emoji.h pretty_ostream.h \
			HydroinfoMUSIC.h fluidCell.h

# -------------------------------------------------

OBJDIR		=	obj
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN)
INSTPATH	=	../

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

#%.cpp:
#	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install

all:		mkobjdir $(TARGET)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS) 
#		strip $(TARGET)

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm -r obj

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
minmod.cpp : minmod.h data.h
advance.cpp : advance.h util.h dissipative.h minmod.h grid.h data.h eos.h reconst.h u_derivative.h
dissipative.cpp : dissipative.h minmod.h util.h grid.h data.h eos.h
u_derivative.cpp : u_derivative.h util.h data.h grid.h eos.h 
eos.cpp : eos.h util.h data.h 
evolve.cpp : evolve.h data.h eos.h grid.h reconst.h advance.h util.h dissipative.h minmod.h u_derivative.h
grid.cpp : grid.h util.h data.h eos.h
init.cpp : init.h eos.h grid.h util.h data.h glauber.h 
reconst.cpp : reconst.h data.h eos.h grid.h util.h
util.cpp : util.h 
glauber.cpp : glauber.h util.h data.h 
freeze.cpp : freeze.h data.h eos.h grid.h util.h int.h
freeze_pseudo.cpp : freeze.h data.h eos.h grid.h util.h 
main.cpp : music.h
music.cpp : music.h
read_in_parameters.cpp: read_in_parameters.h
