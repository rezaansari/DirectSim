################################################################################
#####  FranZona/DirectSim                                                #######
##  makefile  for progs/  Alex Abate, Reza Ansari , April 2016                ##
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
ROOTINC :=
MYLIB :=
ROOTLIB := 


#   Targets 
all :  analyzeBPZ  baseSimulation  calculateKcorrections  cfhtColors  colorDistributions convertSEDS \
      fitLSSTspectra lineOfSightLymanAlpha  lymanAlphaToDensity photoZdist \
      rdlss sdssElColors sdssPicklesLibrary simdensity \
      simulateAbsorberLinesOfSight simulateCFHTobs simulateLSSTobs simulateLSSTobsFromTruth
#### Targets removed (need ROOT)
## addIGMToSED  lineOfSightMagnitude  priorFitter  projectTemplates

clean:
	rm $(OBJ)/addIGMToSED.o $(OBJ)/analyzeBPZ.o  $(OBJ)/baseSimulation.o  $(OBJ)/calculateKcorrections.o \
           $(OBJ)/cfhtColors.o  $(OBJ)/colorDistributions.o $(OBJ)/convertSEDS.o $(OBJ)/fitLSSTspectra.o \
           $(OBJ)/lineOfSightLymanAlpha.o $(OBJ)/lineOfSightMagnitude.o $(OBJ)/lymanAlphaToDensity.o $(OBJ)/photoZdist.o \
           $(OBJ)/priorFitter.o $(OBJ)/projectTemplates.o $(OBJ)/rdlss.o $(OBJ)/sdssElColors.o $(OBJ)/sdssPicklesLibrary.o \
           $(OBJ)/simdensity.o $(OBJ)/simulateAbsorberLinesOfSight.o $(OBJ)/simulateCFHTobs.o \
           $(OBJ)/simulateLSSTobs.o $(OBJ)/simulateLSSTobsFromTruth.o

	rm $(EXE)/addIGMToSED $(EXE)/analyzeBPZ  $(EXE)/baseSimulation  $(EXE)/calculateKcorrections  $(EXE)/cfhtColors  \
           $(EXE)/colorDistributions $(EXE)/convertSEDS $(EXE)/fitLSSTspectra $(EXE)/lineOfSightLymanAlpha \
           $(EXE)/lineOfSightMagnitude $(EXE)/lymanAlphaToDensity $(EXE)/photoZdist \
           $(EXE)/priorFitter $(EXE)/projectTemplates $(EXE)/rdlss $(EXE)/sdssElColors $(EXE)/sdssPicklesLibrary \
           $(EXE)/simdensity $(EXE)/simulateAbsorberLinesOfSight $(EXE)/simulateCFHTobs $(EXE)/simulateLSSTobs \
           $(EXE)/simulateLSSTobsFromTruth

# MAIN PROGS

addIGMToSED : $(EXE)/addIGMToSED
	@echo 'makefile : addIGMToSED made'
analyzeBPZ : $(EXE)/analyzeBPZ
	@echo 'makefile : analyzeBPZ made'
baseSimulation : $(EXE)/baseSimulation
	@echo 'makefile : baseSimulation made'
calculateKcorrections : $(EXE)/calculateKcorrections
	@echo 'makefile : calculateKcorrections made'
cfhtColors : $(EXE)/cfhtColors
	@echo 'makefile : cfhtColors made'
colorDistributions	: $(EXE)/colorDistributions
	@echo 'makefile : colorDistributions made'
convertSEDS	: $(EXE)/convertSEDS
	@echo 'makefile : convertSEDS made'
fitLSSTspectra : $(EXE)/fitLSSTspectra
	@echo 'makefile : fitLSSTspectra made'
lineOfSightLymanAlpha : $(EXE)/lineOfSightLymanAlpha
	@echo 'makefile : lineOfSightLymanAlpha made'
lineOfSightMagnitude : $(EXE)/lineOfSightMagnitude
	@echo 'makefile : lineOfSightMagnitude made'
lsstPicklesLibrary : $(EXE)/lsstPicklesLibrary
	@echo 'makefile : lsstPicklesLibrary made'
lymanAlphaToDensity : $(EXE)/lymanAlphaToDensity
	@echo 'makefile : lymanAlphaToDensity made'
pcaTemplates : $(EXE)/pcaTemplates
	@echo 'makefile : pcaTemplates made'
photoZdist : $(EXE)/photoZdist
	@echo 'makefile : photoZdist made'
priorFitter : $(EXE)/priorFitter
	@echo 'makefile : priorFitter made'
projectTemplates : $(EXE)/projectTemplates
	@echo 'makefile : projectTemplates made'
rdlss : $(EXE)/rdlss
	@echo 'makefile : rdlss made'
sdssElColors : $(EXE)/sdssElColors
	@echo 'makefile : sdssElColors made'
sdssPicklesLibrary : $(EXE)/sdssPicklesLibrary
	@echo 'makefile : sdssPicklesLibrary made'
simdensity : $(EXE)/simdensity
	@echo 'makefile : simdensity made'
simulateAbsorberLinesOfSight : $(EXE)/simulateAbsorberLinesOfSight
	@echo 'makefile : simulateAbsorberLinesOfSight made'
simulateCFHTobs : $(EXE)/simulateCFHTobs
	@echo 'makefile : simulateCFHTobs made'
simulateLSSTobs : $(EXE)/simulateLSSTobs
	@echo 'makefile : simulateLSSTobs made'	
simulateLSSTobsFromTruth : $(EXE)/simulateLSSTobsFromTruth
	@echo 'makefile : simulateLSSTobsFromTruth made'
# TESTING PROGS
test2Dinterp : $(EXE)/test2Dinterp 
	@echo 'makefile :test2Dinterp made'
testbasesim : $(EXE)/testbasesim 
	@echo 'makefile :testbasesim made'
test : $(EXE)/test
	@echo 'makefile : test made'
testEMalgorithm : $(EXE)/testEMalgorithm
	@echo 'makefile : testEMalgorithm made'
testErrors : $(EXE)/testErrors
	@echo 'makefile : testErrors made'
testgoodsmagsim : $(EXE)/testgoodsmagsim
	@echo 'makefile : testgoodsmagsim made'
testKcorrColors : $(EXE)/testKcorrColors 
	@echo 'makefile : testKcorrColors made'
testKcorrMethod : $(EXE)/testKcorrMethod 
	@echo 'makefile : testKcorrMethod made'
testLF : $(EXE)/testLF 
	@echo 'makefile :testLF made'
testLymanAlphaAbs : $(EXE)/testLymanAlphaAbs 
	@echo 'makefile :testLymanAlphaAbs made'
testMadau : $(EXE)/testMadau
	@echo 'makefile : testMadau made'
testMeiksin : $(EXE)/testMeiksin
	@echo 'makefile : testMeiksin made'
testSimReadKcorr : $(EXE)/testSimReadKcorr
	@echo 'makefile : testSimReadKcorr made'
testsimulateIGM : $(EXE)/testsimulateIGM
	@echo 'makefile : testsimulateIGM made'
testSimulation : $(EXE)/testSimulation
	@echo 'makefile : testSimulation made'
testTemplateFitting : $(EXE)/testTemplateFitting
	@echo 'makefile : testTemplateFitting made'
## programs below here have not been CHECKED or maybe even finished...
testpsdenscube : $(EXE)/testpsdenscube
	@echo 'makefile :  testpsdenscube made'

testsimdensity : $(EXE)/testsimdensity
	@echo 'makefile :  testsimdensity made'

# ADD LINE OF SIGHT TRANSMISSON TO SEDS IN A LIBRARY
$(EXE)/addIGMToSED : $(OBJ)/addIGMToSED.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/addIGMToSED $(OBJ)/addIGMToSED.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/addIGMToSED.o : addIGMToSED.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/addIGMToSED.o addIGMToSED.cc

# ANALYZE A BPZ CATALOG
$(EXE)/analyzeBPZ : $(OBJ)/analyzeBPZ.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/analyzeBPZ $(OBJ)/analyzeBPZ.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/analyzeBPZ.o : analyzeBPZ.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/analyzeBPZ.o analyzeBPZ.cc 

# BASE SIMULATION
$(EXE)/baseSimulation : $(OBJ)/baseSimulation.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/baseSimulation $(OBJ)/baseSimulation.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/baseSimulation.o : baseSimulation.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/baseSimulation.o baseSimulation.cc

# CALCULATE K CORRECTION TABLES
$(EXE)/calculateKcorrections : $(OBJ)/calculateKcorrections.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/calculateKcorrections $(OBJ)/calculateKcorrections.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/calculateKcorrections.o : calculateKcorrections.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/calculateKcorrections.o calculateKcorrections.cc

# CALCULATE U-G, I-Z CFHT colors w/ and wo/ host galaxy reddening
$(EXE)/cfhtColors : $(OBJ)/cfhtColors.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/cfhtColors $(OBJ)/cfhtColors.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/cfhtColors.o : cfhtColors.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/cfhtColors.o cfhtColors.cc

# COLOR DISTRIBUTIONS
$(EXE)/colorDistributions : $(OBJ)/colorDistributions.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/colorDistributions $(OBJ)/colorDistributions.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/colorDistributions.o : colorDistributions.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/colorDistributions.o colorDistributions.cc

# CONVERT SED UNITS
$(EXE)/convertSEDS : $(OBJ)/convertSEDS.o
	$(CXXLINK) -o $(EXE)/convertSEDS $(OBJ)/convertSEDS.o $(SOPHYAEXTSLBLIST)

$(OBJ)/convertSEDS.o : convertSEDS.cc
	$(CXXCOMPILE) -o $(OBJ)/convertSEDS.o convertSEDS.cc 

# FIT LSST SPECTRA TO CWWK
$(EXE)/fitLSSTspectra : $(OBJ)/fitLSSTspectra.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/fitLSSTspectra $(OBJ)/fitLSSTspectra.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/fitLSSTspectra.o : fitLSSTspectra.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/fitLSSTspectra.o fitLSSTspectra.cc

# SIMULATE LINE OF SIGHT LYMAN ALPHA TRANSMISSION
$(EXE)/lineOfSightLymanAlpha : $(OBJ)/lineOfSightLymanAlpha.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/lineOfSightLymanAlpha $(OBJ)/lineOfSightLymanAlpha.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/lineOfSightLymanAlpha.o : lineOfSightLymanAlpha.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/lineOfSightLymanAlpha.o lineOfSightLymanAlpha.cc

# U,G MAGNITUDES WITH IGM LINE OF SIGHT
$(EXE)/lineOfSightMagnitude : $(OBJ)/lineOfSightMagnitude.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/lineOfSightMagnitude $(OBJ)/lineOfSightMagnitude.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/lineOfSightMagnitude.o : lineOfSightMagnitude.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/lineOfSightMagnitude.o lineOfSightMagnitude.cc

# CALCULATE LSST COLORS FOR PICKLES' LIBRARY OF STARS 
$(EXE)/lsstPicklesLibrary : $(OBJ)/lsstPicklesLibrary.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/lsstPicklesLibrary $(OBJ)/lsstPicklesLibrary.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/lsstPicklesLibrary.o : lsstPicklesLibrary.cc $(LISTH)
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/lsstPicklesLibrary.o lsstPicklesLibrary.cc

# LYMAN ALPHA ALONG LINE OF SIGHT CONVERTED TO DENSITY
$(EXE)/lymanAlphaToDensity : $(OBJ)/lymanAlphaToDensity.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/lymanAlphaToDensity $(OBJ)/lymanAlphaToDensity.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/lymanAlphaToDensity.o : lymanAlphaToDensity.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/lymanAlphaToDensity.o lymanAlphaToDensity.cc

# TEMPLATE PCA
$(EXE)/pcaTemplates : $(OBJ)/pcaTemplates.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/pcaTemplates $(OBJ)/pcaTemplates.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/pcaTemplates.o : pcaTemplates.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/pcaTemplates.o pcaTemplates.cc 

# PHOTO-Z DISTRIBUTION
$(EXE)/photoZdist : $(OBJ)/photoZdist.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/photoZdist $(OBJ)/photoZdist.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/photoZdist.o : photoZdist.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/photoZdist.o photoZdist.cc 

# PRIOR FITTER
$(EXE)/priorFitter : $(OBJ)/priorFitter.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/priorFitter $(OBJ)/priorFitter.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB) $(MINUIT)

$(OBJ)/priorFitter.o : priorFitter.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/priorFitter.o priorFitter.cc

# PROJECT TEMPLATES
$(EXE)/projectTemplates : $(OBJ)/projectTemplates.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/projectTemplates $(OBJ)/projectTemplates.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/projectTemplates.o : projectTemplates.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/projectTemplates.o projectTemplates.cc 

# SIMULATE CATALOG OF BASIC GALAXY PROPERTIES FROM OVER-DENSITY GRID
$(EXE)/rdlss : $(OBJ)/rdlss.o $(FZLIBS)
	$(CXXLINK) -o $(EXE)/rdlss $(OBJ)/rdlss.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/rdlss.o : rdlss.cc $(LISTH)
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/rdlss.o rdlss.cc

# CALCULATE SDSS COLORS OF ELLIPTICAL GALAXY
$(EXE)/sdssElColors : $(OBJ)/sdssElColors.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/sdssElColors $(OBJ)/sdssElColors.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/sdssElColors.o : sdssElColors.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/sdssElColors.o sdssElColors.cc 

# CALCULATE SDSS COLORS FOR PICKLES' LIBRARY OF STARS 
$(EXE)/sdssPicklesLibrary : $(OBJ)/sdssPicklesLibrary.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/sdssPicklesLibrary $(OBJ)/sdssPicklesLibrary.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/sdssPicklesLibrary.o : sdssPicklesLibrary.cc $(LISTH)
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/sdssPicklesLibrary.o sdssPicklesLibrary.cc
# SIMULATE OVERDENSITY GRID
$(EXE)/simdensity : $(OBJ)/simdensity.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/simdensity $(OBJ)/simdensity.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/simdensity.o : simdensity.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/simdensity.o simdensity.cc 

# SIMULATE LINE OF SIGHT ABSORBER DISTRIBUTIONS 
$(EXE)/simulateAbsorberLinesOfSight : $(OBJ)/simulateAbsorberLinesOfSight.o $(FZLIBS)
	$(CXXLINK) -o $(EXE)/simulateAbsorberLinesOfSight $(OBJ)/simulateAbsorberLinesOfSight.o \
	$(FZLIBS) $(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/simulateAbsorberLinesOfSight.o : simulateAbsorberLinesOfSight.cc $(LISTH)
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/simulateAbsorberLinesOfSight.o \
	simulateAbsorberLinesOfSight.cc 

# SIMULATE CFHT OBSERVATIONS FROM BASE SIMULATION
$(EXE)/simulateCFHTobs : $(OBJ)/simulateCFHTobs.o $(FZLIBS)
	$(CXXLINK) -o $(EXE)/simulateCFHTobs $(OBJ)/simulateCFHTobs.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/simulateCFHTobs.o : simulateCFHTobs.cc $(LISTH)
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/simulateCFHTobs.o \
	simulateCFHTobs.cc 

# SIMULATE LSST OBSERVATIONS FROM BASE SIMULATION
$(EXE)/simulateLSSTobs : $(OBJ)/simulateLSSTobs.o $(FZLIBS)
	$(CXXLINK) -o $(EXE)/simulateLSSTobs $(OBJ)/simulateLSSTobs.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/simulateLSSTobs.o : simulateLSSTobs.cc $(LISTH)
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/simulateLSSTobs.o \
	simulateLSSTobs.cc 

# SIMULATE LSST OBSERVATIONS FROM INPUT IMSIM CATALOG OF TRUE MAGS AND REDSHIFTS
$(EXE)/simulateLSSTobsFromTruth : $(OBJ)/simulateLSSTobsFromTruth.o $(FZLIBS)
	$(CXXLINK) -o $(EXE)/simulateLSSTobsFromTruth $(OBJ)/simulateLSSTobsFromTruth.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/simulateLSSTobsFromTruth.o : simulateLSSTobsFromTruth.cc $(LISTH)
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/simulateLSSTobsFromTruth.o \
	simulateLSSTobsFromTruth.cc 

###################### TESTING PROGRAMS ########################################

# TEST 2D INTERPOLATION
$(EXE)/test2Dinterp : $(OBJ)/test2Dinterp.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/test2Dinterp $(OBJ)/test2Dinterp.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/test2Dinterp.o : test2Dinterp.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/test2Dinterp.o test2Dinterp.cc

# TEST BASESIM
$(EXE)/testbasesim : $(OBJ)/testbasesim.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/testbasesim $(OBJ)/testbasesim.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testbasesim.o : testbasesim.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/testbasesim.o testbasesim.cc

# TEST EXPECTATION-MAXIMIZATION ALGORITHM
$(EXE)/testEMalgorithm : $(OBJ)/testEMalgorithm.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/testEMalgorithm $(OBJ)/testEMalgorithm.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testEMalgorithm.o : testEMalgorithm.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/testEMalgorithm.o testEMalgorithm.cc

# TEST LSST ERRORS
$(EXE)/testErrors : $(OBJ)/testErrors.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/testErrors $(OBJ)/testErrors.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testErrors.o : testErrors.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/testErrors.o testErrors.cc

# TEST GOODS SIM
$(EXE)/testgoodsmagsim : $(OBJ)/testgoodsmagsim.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/testgoodsmagsim $(OBJ)/testgoodsmagsim.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testgoodsmagsim.o : testgoodsmagsim.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/testgoodsmagsim.o testgoodsmagsim.cc

# TEST K CORRECTION 
$(EXE)/testKcorrColors : $(OBJ)/testKcorrColors.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/testKcorrColors $(OBJ)/testKcorrColors.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testKcorrColors.o : testKcorrColors.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/testKcorrColors.o testKcorrColors.cc

# TEST K CORRECTION INTERPOLATION
$(EXE)/testKcorrMethod : $(OBJ)/testKcorrMethod.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/testKcorrMethod $(OBJ)/testKcorrMethod.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testKcorrMethod.o : testKcorrMethod.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/testKcorrMethod.o testKcorrMethod.cc

# TEST LF 
$(EXE)/testLF : $(OBJ)/testLF.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/testLF $(OBJ)/testLF.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testLF.o : testLF.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) -o $(OBJ)/testLF.o testLF.cc 

# TEST Lyman-alpha absorption calculation parts 
$(EXE)/testLymanAlphaAbs : $(OBJ)/testLymanAlphaAbs.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/testLymanAlphaAbs $(OBJ)/testLymanAlphaAbs.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testLymanAlphaAbs.o : testLymanAlphaAbs.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) -o $(OBJ)/testLymanAlphaAbs.o testLymanAlphaAbs.cc 

# TEST MADAU
$(EXE)/testMadau : $(OBJ)/testMadau.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/testMadau $(OBJ)/testMadau.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testMadau.o : testMadau.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/testMadau.o testMadau.cc

# TEST MEIKSIN
$(EXE)/testMeiksin : $(OBJ)/testMeiksin.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/testMeiksin $(OBJ)/testMeiksin.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testMeiksin.o : testMeiksin.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/testMeiksin.o testMeiksin.cc

# TEST SIMULATE DATA USING K CORRECTION TABLES
$(EXE)/testSimReadKcorr : $(OBJ)/testSimReadKcorr.o $(FZLIBS)
	$(CXXLINK) -o $(EXE)/testSimReadKcorr $(OBJ)/testSimReadKcorr.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB) 

$(OBJ)/testSimReadKcorr.o : testSimReadKcorr.cc $(LISTH)
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/testSimReadKcorr.o testSimReadKcorr.cc 

# TEST SIMULATE IGM
$(EXE)/testsimulateIGM : $(OBJ)/testsimulateIGM.o $(FZLIBS)
	$(CXXLINK) -o $(EXE)/testsimulateIGM $(OBJ)/testsimulateIGM.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB) 

$(OBJ)/testsimulateIGM.o : testsimulateIGM.cc $(LISTH)
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/testsimulateIGM.o testsimulateIGM.cc

# TEST SIMULATION
$(EXE)/testSimulation : $(OBJ)/testSimulation.o $(FZLIBS)
	$(CXXLINK) -o $(EXE)/testSimulation $(OBJ)/testSimulation.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB) 

$(OBJ)/testSimulation.o : testSimulation.cc $(LISTH)
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/testSimulation.o testSimulation.cc

# TEST TEMPLATE FITTING
$(EXE)/testTemplateFitting : $(OBJ)/testTemplateFitting.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/testTemplateFitting $(OBJ)/testTemplateFitting.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testTemplateFitting.o : testTemplateFitting.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/testTemplateFitting.o testTemplateFitting.cc

# TEMPORARY TEST CODE
$(EXE)/test : $(OBJ)/test.o $(FZLIBS)
	$(CXXLINK) -o $(EXE)/test $(OBJ)/test.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB) 

$(OBJ)/test.o : test.cc $(LISTH)
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/test.o test.cc

## programs below here have not been CHECKED or maybe even finished...

# TEST POWER SPECTRUM FROM OVER-DENSITY CUBE
$(EXE)/testpsdenscube : $(OBJ)/testpsdenscube.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/testpsdenscube $(OBJ)/testpsdenscube.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testpsdenscube.o : testpsdenscube.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/testpsdenscube.o testpsdenscube.cc 

# TEST DENSITY SIMULATION
$(EXE)/testsimdensity : $(OBJ)/testsimdensity.o $(FZLIBS) 
	$(CXXLINK) -o $(EXE)/testsimdensity $(OBJ)/testsimdensity.o $(FZLIBS) \
	$(SOPHYAEXTSLBLIST) $(MYLIB) $(ROOTLIB)

$(OBJ)/testsimdensity.o : testsimdensity.cc $(LISTH)  
	$(CXXCOMPILE) -I$(MYCL) $(ROOTINC) -o $(OBJ)/testsimdensity.o testsimdensity.cc 
