PWD = `pwd`

DEBUG =
OPTIMIZATION = -O3
FLOWTRACE =
PARALL = -fopenmp
CFLAGS = $(DEBUG) $(OPTIMIZATION) $(FLOWTRACE) $(PARALL)
COMPILER = g++-6
LIBS = -L/usr/local/lib -lm -lgsl -lgslcblas 
INCLUDES = -I. -I/usr/local/include 

.SUFFIXES: .o .cpp .h

.cpp.o :
	$(COMPILER) $(CFLAGS) $(INCLUDES) -c $*.cpp
S1 =\
paramreader.cpp\
UVH2+1.cpp
S2 =\
paramreader.cpp\
convert.cpp\
convertnew.cpp
S3 =\
extract.cpp
S3b =\
paramreader.cpp\
extractfull.cpp
S3c =\
paramreader.cpp\
extractfullangle.cpp
S4 =\
paramreader.cpp\
initE.cpp
S4b=\
paramreader.cpp\
calcS.cpp
S5 =\
paramreader.cpp\
givedist.cpp
S6 =\
paramreader.cpp\
prereso.cpp
S7 =\
paramreader.cpp\
preresofull.cpp
S7b =\
paramreader.cpp\
preresofull.cpp
S8 =\
chempot.cpp
S9 =\
generate.cpp
OBJ1 =\
paramreader.o\
UVH2+1.o
OBJ2 =\
paramreader.o\
convert.o
OBJ2b=\
paramreader.o\
convertnew.o
OBJ2c=\
paramreader.o\
convertfull.o
OBJ2d=\
paramreader.o\
convertlinear.o
OBJ2e=\
paramreader.o\
convertfulltime.o
OBJ3 =\
extract.o
OBJ3b =\
paramreader.o\
extractfull.o
OBJ3c =\
paramreader.o\
extractfullangle.o
OBJ4=\
paramreader.o\
initE.o
OBJ4b=\
paramreader.o\
calcS.o
OBJ5=\
paramreader.o\
givedist.o
OBJ6 =\
paramreader.o\
prereso.o
OBJ7 =\
paramreader.o\
preresofull.o
OBJ7b =\
paramreader.o\
preresofullangle.o
OBJ8 =\
chempot.o
OBJ9 =\
generate.o
EXE1 =\
vh2
EXE2 =\
convert
EXE2b =\
convertnew
EXE2c =\
convertfull
EXE2d =\
convertlinear
EXE2e =\
convertfulltime
EXE3 =\
extract
EXE3b =\
extractfull
EXE3c =\
extractfullangle
EXE4 =\
initE
EXE4b =\
calcS
EXE5 =\
givedist
EXE6 =\
prereso
EXE7 =\
preresofull
EXE7b =\
preresofullangle
EXE8 =\
chempot
EXE9 =\
generate

$(EXE1) : $(OBJ1) 
$(EXE2) : $(OBJ2) 
$(EXE2b) : $(OBJ2b) 
$(EXE2c) : $(OBJ2c) 
$(EXE2d) : $(OBJ2d)
$(EXE2e) : $(OBJ2e)
$(EXE3) : $(OBJ3) 
$(EXE3b) : $(OBJ3b)
$(EXE3c) : $(OBJ3c)
$(EXE4) : $(OBJ4)
$(EXE4b) : $(OBJ4b)
$(EXE5) : $(OBJ5)
$(EXE6) : $(OBJ6)
$(EXE7) : $(OBJ7)
$(EXE7b) : $(OBJ7b)
$(EXE8) : $(OBJ8)
$(EXE9) : $(OBJ9)

vh2:	
	$(COMPILER) $(PARALL) $(OBJ1) -o $(EXE1)-2.1  $(LIBS)

convert:
	$(COMPILER) $(PARALL) $(OBJ2) -o $(EXE2)  $(LIBS)

convertnew:
	$(COMPILER) $(PARALL) $(OBJ2b) -o $(EXE2b)  $(LIBS)

convertfull:
	$(COMPILER) $(PARALL) $(OBJ2c) -o $(EXE2c)  $(LIBS)

convertlinear:
	$(COMPILER) $(PARALL) $(OBJ2d) -o $(EXE2d)  $(LIBS)

convertfulltime:
	$(COMPILER) $(PARALL) $(OBJ2e) -o $(EXE2e)  $(LIBS)

extract:
	$(COMPILER) $(PARALL) $(OBJ3) -o $(EXE3)  $(LIBS)

extractfull:
	$(COMPILER) $(PARALL) $(OBJ3b) -o $(EXE3b)  $(LIBS)

extractfullangle:
	$(COMPILER) $(PARALL) $(OBJ3c) -o $(EXE3c)  $(LIBS)

initE:
	$(COMPILER) $(PARALL) $(OBJ4) -o $(EXE4)  $(LIBS)

calcS:
	$(COMPILER) $(PARALL) $(OBJ4b) -o $(EXE4b)  $(LIBS)

givedist:
	$(COMPILER) $(PARALL) $(OBJ5) -o $(EXE5) $(LIBS)

prereso:
	$(COMPILER) $(PARALL) $(OBJ6) -o $(EXE6) $(LIBS)

preresofull:
	$(COMPILER) $(PARALL) $(OBJ7) -o $(EXE7) $(LIBS)

preresofullangle:
	$(COMPILER) $(PARALL) $(OBJ7b) -o $(EXE7b) $(LIBS)

generate:
	$(COMPILER) $(PARALL) $(OBJ9) -o $(EXE9)  $(LIBS)


# clean up misc files
clean :
	rm -f $(EXE1) $(EXE2) $(EXE2b) $(EXE2c) $(EXE2d) $(EXE2e) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE9) *\.o *~ #* core 
