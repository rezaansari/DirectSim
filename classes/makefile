################################################################################
#####  FranZona/DirectSim                                                #######
##  makefile  for classes/  Alex Abate, Reza Ansari , April 2016              ##
################################################################################

include $(SOPHYABASE)/include/sophyamake.inc

#define GEN3D_FLOAT

OBJ = ${FZGALSIM}/obj
EXE = ${FZGALSIM}/exe
LIB = ${FZGALSIM}/lib

MYCL = ../classes

LISTH := $(MYCL)/cosmocalcs.h $(MYCL)/geneutils.h $(MYCL)/gftdist.h \
$(MYCL)/schechter.h $(MYCL)/sinterp.h $(MYCL)/simdata.h $(MYCL)/reddening.h \
$(MYCL)/sedfilter.h $(MYCL)/genefluct3d.h  $(MYCL)/pkspectrum.h \
$(MYCL)/mass2gal.h $(MYCL)/powerspec.h $(MYCL)/matrix.h $(MYCL)/igm.h \
$(MYCL)/hpoly.h $(MYCL)/shapelets.h $(MYCL)/em.h $(MYCL)/cat2grid.h \
$(MYCL)/fitkbaoscale.h $(MYCL)/chisqstats.h
#### Deactivated files 
# $(MYCL)/sedpca.h $(MYCL)/constcosmo.h $(MYCL)/root_plots.h

LISTLIBOBJ := $(OBJ)/cosmocalcs.o $(OBJ)/geneutils.o $(OBJ)/gftdist.o \
$(OBJ)/schechter.o $(OBJ)/sinterp.o $(OBJ)/simdata.o $(OBJ)/reddening.o \
$(OBJ)/sedfilter.o  $(OBJ)/genefluct3d.o  $(OBJ)/pkspectrum.o $(OBJ)/mass2gal.o \
$(OBJ)/powerspec.o $(OBJ)/matrix.o $(OBJ)/igm.o $(OBJ)/hpoly.o $(OBJ)/shapelets.o\
$(OBJ)/em.o $(OBJ)/cat2grid.o $(OBJ)/fitkbaoscale.o $(OBJ)/chisqstats.o

#### Deactivated files 
#  $(OBJ)/sedpca.o  

### Target definition : default=all , build the library 
# all : $(OBJ)/sedfilter.o
# all : $(LISTLIBOBJ)
all : $(LIB)/libFZGalCatSim.a

clean :
	rm $(LISTLIBOBJ) ; rm $(LIB)/libFZGalCatSim.a 

###  Building library
$(LIB)/libFZGalCatSim.a : $(LISTLIBOBJ)
	$(AR) $(ARFLAGS) $(LIB)/libFZGalCatSim.a $(LISTLIBOBJ)

###  Individual file compilation 
$(OBJ)/cat2grid.o : cat2grid.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/cat2grid.o cat2grid.cc
$(OBJ)/chisqstats.o : chisqstats.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/chisqstats.o chisqstats.cc
$(OBJ)/cosmocalcs.o : cosmocalcs.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/cosmocalcs.o cosmocalcs.cc
$(OBJ)/em.o : em.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/em.o em.cc
$(OBJ)/fitkbaoscale.o : fitkbaoscale.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/fitkbaoscale.o fitkbaoscale.cc
$(OBJ)/genefluct3d.o : genefluct3d.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/genefluct3d.o genefluct3d.cc
$(OBJ)/geneutils.o : geneutils.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/geneutils.o geneutils.cc
$(OBJ)/gftdist.o : gftdist.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/gftdist.o gftdist.cc
$(OBJ)/hpoly.o : hpoly.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/hpoly.o hpoly.cc
$(OBJ)/igm.o : igm.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/igm.o igm.cc
$(OBJ)/mass2gal.o : mass2gal.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/mass2gal.o mass2gal.cc
$(OBJ)/massfunc.o : massfunc.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/massfunc.o massfunc.cc
$(OBJ)/matrix.o : matrix.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/matrix.o matrix.cc
$(OBJ)/pkspectrum.o : pkspectrum.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/pkspectrum.o pkspectrum.cc
$(OBJ)/powerspec.o : powerspec.cc
	$(CXXCOMPILE) -o $(OBJ)/powerspec.o powerspec.cc
$(OBJ)/readsim.o : readsim.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/readsim.o readsim.cc
$(OBJ)/reddening.o : reddening.cc
	$(CXXCOMPILE) -o $(OBJ)/reddening.o reddening.cc
$(OBJ)/schechter.o : schechter.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/schechter.o schechter.cc
$(OBJ)/sedfilter.o : sedfilter.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/sedfilter.o sedfilter.cc
$(OBJ)/sedpca.o : sedpca.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/sedpca.o sedpca.cc
$(OBJ)/selectfunc.o : selectfunc.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/selectfunc.o selectfunc.cc
$(OBJ)/shapelets.o : shapelets.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/shapelets.o shapelets.cc
$(OBJ)/simdata.o : simdata.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/simdata.o simdata.cc
$(OBJ)/sinterp.o : sinterp.cc  $(LISTH) 
	$(CXXCOMPILE) -o $(OBJ)/sinterp.o sinterp.cc
