################################################################################
#####  FranZona/DirectSim                                                #######
##  makefile  for baoprogs/  Alex Abate, Reza Ansari , April 2016             ##
################################################################################

include $(SOPHYABASE)/include/sophyamake.inc

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

FZLIBS := $(LIB)/libFZGalCatSim.a

#   Targets 
all : subfromfull computepsfromarray  addGausszerr  fitkbao  grid_data  getpzconvf getsf sim_mcgrids

clean:
	rm $(OBJ)/subfromfull.o $(OBJ)/computepsfromarray.o  $(OBJ)/addGausszerr.o  $(OBJ)/fitkbao.o  $(OBJ)/grid_data.o \
	$(OBJ)/getpzconvf.o $(OBJ)/getsf.o $(OBJ)/sim_mcgrids.o
	rm $(EXE)/subfromfull $(EXE)/computepsfromarray  $(EXE)/addGausszerr  $(EXE)/fitkbao  $(EXE)/grid_data \
	$(EXE)/getpzconvf $(EXE)/getsf $(EXE)/sim_mcgrids

# GET DATA SUB GRID
subfromfull : $(EXE)/subfromfull
	@echo 'makefile : subfromfull made'
$(EXE)/subfromfull : $(OBJ)/subfromfull.o $(FZLIBS)
	$(CXXLINK) -o $(EXE)/subfromfull $(OBJ)/subfromfull.o $(FZLIBS) $(SOPHYAEXTSLBLIST)
$(OBJ)/subfromfull.o : subfromfull.cc $(LISTH)
	$(CXXCOMPILE) -I$(MYCL)  -o $(OBJ)/subfromfull.o subfromfull.cc

# COMPUTE POWER SPECTRA
computepsfromarray : $(EXE)/computepsfromarray
	@echo 'makefile : computepsfromarray made'
$(EXE)/computepsfromarray : $(OBJ)/computepsfromarray.o $(FZLIBS)
	$(CXXLINK) -o $(EXE)/computepsfromarray $(OBJ)/computepsfromarray.o $(FZLIBS) $(SOPHYAEXTSLBLIST) 
$(OBJ)/computepsfromarray.o : computepsfromarray.cc $(LISTH)
	$(CXXCOMPILE) -I$(MYCL) -o $(OBJ)/computepsfromarray.o computepsfromarray.cc 

# ADD GAUSSIAN Z ERROR TO CATALOG
addGausszerr : $(EXE)/addGausszerr
	@echo 'makefile : addGausszerr made'
$(EXE)/addGausszerr : $(OBJ)/addGausszerr.o $(FZLIBS)
	$(CXXLINK) -o $(EXE)/addGausszerr $(OBJ)/addGausszerr.o $(FZLIBS) $(SOPHYAEXTSLBLIST)
$(OBJ)/addGausszerr.o : addGausszerr.cc $(LISTH)
	$(CXXCOMPILE) -I$(MYCL) -o $(OBJ)/addGausszerr.o addGausszerr.cc 

# FIT K BAO
fitkbao : $(EXE)/fitkbao
	@echo 'makefile : fitkbao made'
$(EXE)/fitkbao : $(OBJ)/fitkbao.o $(FZLIBS)
	$(CXXLINK) -o $(EXE)/fitkbao $(OBJ)/fitkbao.o $(FZLIBS) $(SOPHYAEXTSLBLIST) 

$(OBJ)/fitkbao.o : fitkbao.cc $(LISTH)
	$(CXXCOMPILE) -I$(MYCL) -o $(OBJ)/fitkbao.o fitkbao.cc

# GRID GALAXY DATA
grid_data : $(EXE)/grid_data
	@echo 'makefile : grid_data made'
$(EXE)/grid_data : $(OBJ)/grid_data.o $(FZLIBS)
	$(CXXLINK) -o $(EXE)/grid_data $(OBJ)/grid_data.o $(FZLIBS) $(SOPHYAEXTSLBLIST) 

$(OBJ)/grid_data.o : grid_data.cc $(LISTH)
	$(CXXCOMPILE) -I$(MYCL) -o $(OBJ)/grid_data.o grid_data.cc

# CALCULATE PHOTO-Z CONVOLUTION FUNCTION
getpzconvf : $(EXE)/getpzconvf
	@echo 'makefile : getpzconvf made'
$(EXE)/getpzconvf : $(OBJ)/getpzconvf.o $(FZLIBS)
	$(CXXLINK) -o $(EXE)/getpzconvf $(OBJ)/getpzconvf.o $(FZLIBS) $(SOPHYAEXTSLBLIST) 
$(OBJ)/getpzconvf.o : getpzconvf.cc $(LISTH)
	$(CXXCOMPILE) -I$(MYCL) -o $(OBJ)/getpzconvf.o getpzconvf.cc

# CALCULATE SELECTION FUNCTION OF OBSERVED CATALOG
getsf : $(EXE)/getsf
	@echo 'makefile : getsf made'
$(EXE)/getsf : $(OBJ)/getsf.o $(FZLIBS)
	$(CXXLINK) -o $(EXE)/getsf $(OBJ)/getsf.o $(FZLIBS) $(SOPHYAEXTSLBLIST) 
$(OBJ)/getsf.o : getsf.cc $(LISTH)
	$(CXXCOMPILE) -I$(MYCL) -o $(OBJ)/getsf.o getsf.cc 

# SIMULATE RANDOM CATALOGS
sim_mcgrids : $(EXE)/sim_mcgrids
	@echo 'makefile : sim_mcgrids made'
$(EXE)/sim_mcgrids : $(OBJ)/sim_mcgrids.o $(FZLIBS)
	$(CXXLINK) -o $(EXE)/sim_mcgrids $(OBJ)/sim_mcgrids.o $(FZLIBS) $(SOPHYAEXTSLBLIST) 
$(OBJ)/sim_mcgrids.o : sim_mcgrids.cc $(LISTH)
	$(CXXCOMPILE) -I$(MYCL) -o $(OBJ)/sim_mcgrids.o sim_mcgrids.cc

