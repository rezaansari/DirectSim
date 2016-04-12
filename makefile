################################################################################
#####  FranZona/DirectSim                                                #######
##  Top level makefile  Alex Abate, Reza Ansari , April 2016                  ##
################################################################################

include $(SOPHYABASE)/include/sophyamake.inc

#### Automatically set up GALSIM variable if it doesn't exist
ifndef FZGALSIM
	FZGALSIM := ${PWD}
	echo "FZGALSIM = ${GALSIM}"
endif

OBJ = ${FZGALSIM}/obj
EXE = ${FZGALSIM}/exe
LIB = ${FZGALSIM}/lib

all: lib progs baoprogs


lib::
	@echo '----- Doing make in classes/ ...'
	cd classes ; make ; cd ..

progs::
	@echo '----- Doing make in progs/ ...'
	cd progs ; make ; cd ..

baoprogs::
	@echo '----- Doing make in baoprogs/ ...'
	cd baoprogs ; make ; cd ..

clean:
	@echo '----- Cleaning : make clean in classes/ progs/ baoprogs/ ...'
	cd classes ; make clean ; cd ..
	cd progs ; make clean ; cd ..
	cd baoprogs ; make clean ; cd ..

##  Create directories 
depend:
	mkdir -p $(FZGALSIM)/obj
	mkdir -p $(FZGALSIM)/lib
	mkdir -p $(FZGALSIM)/exe

