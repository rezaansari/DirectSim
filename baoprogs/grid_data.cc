/**
  * @file  grid_data.cc
  * @brief grid galaxy data and output arrays of gridded data ready for power 
  *        spectrum computation
  *
  * @todo sky area, cosmology and z dim as radial dir should be read from galaxy 
  *       catalog header instead of supplied to the program as arguments or 
  *       hard coded
  *
  * @author Alex Abate
  * Contact: abate@email.arizona.edu
  *
  */


#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <typeinfo>

// sophya
#include "machdefs.h"
#include "sopnamsp.h"
#include "timing.h"
#include "array.h"
#include "hisprof.h"
#include "histerr.h"
#include "histinit.h"
#include "fiosinit.h"
#include "fitsioserver.h"
#include "fftwserver.h"
#include "swfitsdtable.h"
#include "resusage.h"

// DirectSim
#include "mydefrg.h"
#include "geneutils.h"
#include "cat2grid.h"
#include "powerspec.h"
#include "mass2gal.h"
#include "pkspectrum.h"
#include "fitkbaoscale.h"
#include "chisqstats.h"


void usage(void);
void usage(void) {
		cout << endl<<" Usage: grid_data [...options...]          "<<endl<<endl;
		
		cout << "  Read in a galaxy catalog and output the gridded      "<<endl;
		cout << "  galaxy data. Name of catalog to read is specified with"<<endl;
		cout << "  -C option, root name of file to output gridded data  "<<endl;
		cout << "  to is specifed with -O option                        "<<endl;
		cout << endl;
		
		cout << "  The grid is specified by giving the number of pixels "<<endl;
		cout << "  in each dimension, the pixel size and the redshift of"<<endl;
		cout << "  the center pixel via the -P option                   "<<endl;
		cout << endl;
		
		cout << "  The sky area the catalog covers is set with option -a"<<endl;
		cout << "  this refers to the opening angle in radians that is  "<<endl;
		cout << "  covered by the observation cone. This variable is    "<<endl;
		cout << "  used in determining how many pixels of the gridded   "<<endl;
		cout << "  data should in fact contain data.                    "<<endl;
		cout << endl;
		
		cout << "  If the galaxy catalog was simulated such that the z  "<<endl;
		cout << "  dimension is the radial direction then this is       "<<endl;
		cout << "  indicated with option -r                             "<<endl;
		cout << endl;
		
		cout << "  The names of the redshift columns to read in can be  "<<endl;
		cout << "  supplied with the -z, first argument is the observed "<<endl;
		cout << "  redshift column name and the second (optional)       "<<endl;
		cout << "  argument is the spectroscopic redshift column name.  "<<endl;
		cout << "  both arguments must be separated by a comma.         "<<endl;
		cout << endl;
		
		cout << "  If the catalog is subject to selection effects (not   "<<endl;
		cout << "  all galaxies were observed) then a selection function "<<endl;
		cout << "  must either be computed or read in. To read in the    "<<endl;
		cout << "  selection function pass the root name of the          "<<endl;
		cout << "  selection function file to be read with the -s option."<<endl;
		cout << "  Otherwise if it is to be computed, pass the root name "<<endl;
		cout << "  of the selection function file to be written AND the  "<<endl;
		cout << "  file containing a list of all the redshifts with the  "<<endl;
		cout << "  -s option, separated by a comma.                      "<<endl; 
		cout << endl;
		
		cout << "  Four grids are output to the file passed to the      "<<endl;
		cout << "  program via the -O option:                           "<<endl;
		cout << "  - grid of normalised galaxy number per grid cell     "<<endl;
		cout << "  - grid of weighted normalised galaxy number per grid "<<endl;
		cout << "    cell (weighted by selection function)              "<<endl;
		cout << "  - grid of normalised random catalog per grid cell    "<<endl;
		cout << "    (weighted by selection function)                   "<<endl;
		cout << "  - grid of redshifts at the pixel centers             "<<endl;
		cout << endl;
		cout << "  If there is no selection function correction then the"<<endl;
		cout << "  weighting = 1 and the grid of normalised galaxy      "<<endl;
		cout << "  number and grid of normalised galaxy number are      "<<endl;
        cout << "  identical.                                           "<<endl;
        cout << endl;
        cout << "  A random catalog must be generated to account for the"<<endl;
        cout << "  shot noise. The mean density of this catalog is set  "<<endl;
        cout << "  with the -m option                                   "<<endl;
        cout << endl;
						
		cout << "  This code reads the cosmology from the input catalog. "<<endl;
		cout << endl;

		cout << "  EXAMPLE 1: A galaxy catalog is stored in a file      "<<endl;
		cout << "  catalog.fits that covers a circular area of sky with "<<endl;
		cout << "  radius pi/4. You want the grid to be specified by    "<<endl;
		cout << "  Nx,Ny,Nz=500,500,500 with pixels of size 6 Mpc and   "<<endl;
		cout << "  centered at a redshift of 0.5. The selection function"<<endl;
		cout << "  of the catalog is in file sf.txt. The observed       "<<endl;
		cout << "  redshifts are in the column named 'zp' of the catalog"<<endl;
		cout << "  and the true redshifts are in the column named 'zs'. "<<endl;
		cout << "  The gridded galaxy data must be written to file      "<<endl;
		cout << "  grids.fits:"<<endl;
		cout << endl;
		
		cout << "  $ grid data -C catalog.fits -a 0.7854 -P 500,500,500,0.5,6."<<endl;
		cout << "              -s sf.txt -z zp,zs -O grids.fits   "<<endl<<endl;
		
	cout << " -C : input_catalog : FITS filename containing galaxy catalog"<<endl;
	cout << " -O : out_grids_name : Write gridded data to this FITS file"<<endl;
	cout << " -a : SkyArea : Specify sky area (radians)                 "<<endl;
	cout << " -E : Error : Add a Gaussian error on z-axis of            "<<endl;
	cout << "      sigma = E(1+redshift) (properly translated in Mpc)   "<<endl;
	cout << " -e : Error : Add a Gaussian error on redshift of          "<<endl;
	cout << "      sigma = E(1+redshift) (properly translated in Mpc)   "<<endl;
	cout << " -S : Error : random seed, must be set for simulations     "<<endl;
	cout << " -r : isZRadial: z-dimension of catalog IS radial direction"<<endl;
	cout << " -P : Nx,Ny,Nz,zref,Res : Number of pixels, redshift of    "<<endl;
	cout << "      central pixel, pixel size - to specify how to grid   "<<endl;
	cout << " -z : ZOCol,ZSCol: read OBSERVED redshifts from column named"<<endl;
	cout << "      ZOCol, SPECTRO redshifts from column named ZSCol      "<<endl;
	cout << " -m : nc : Mean density of random grid                      "<<endl;
	cout << " -s : sf_file_root,all_z_file: Do selection function        "<<endl;
	cout << "      correction. If both args are given sf is calculated   "<<endl;
	cout << " -d : debug_out : root stem of output filename objects are  "<<endl;
	cout << "      written to if want to debug                           "<<endl;
	cout << endl;
	}



int main(int narg, char* arg[]) {

	cout << " ==== grid_data_main.cc program , output gridded galaxy data"<<endl;
	cout << "      from galaxy catalog fits file                          ==== " <<endl;
	
	// Make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	InitTim();
	
    // defaults
	// FILES TO READ IN/OUT
	string input_catalog, out_grids_name;
	// CATALOG AND REDSHIFT COLUMN PARAMETERS
	string zcols;		      // list of redshift columns to read in
	bool isZColGiven = false; // z columns given via program argument
	string ZOCol = "z"; // by default OBSERVED redshift col labelled "z" read in
	string ZSCol = "z"; // by default SPECTRO redshift col labelled "z" read in
	double SkyArea = 999;	// Catalog covers angle radius SkyArea [999==full sky]
	bool isZRadial = false; // if true, catalog z-dimension IS radial direction
	// SELECTION FUNCTION CORRECTION PARAMETERS
	string sffiles;			  // list of selection function files
	bool doSFCorr = false;	  // if true, apply selection function correction
	string sf_file_root;      // file name root of selection function (to write/to read)
	string all_z_file;		  // file name of catalog of ALL redshifts in sim
	bool doSFCompute = false; // if true, do selection function computation here
	bool isForceZspec =false; // force selection function to be computed using the SPEC-z
	// GRID SPEC PARS
	double R = 8.;		    // Grid cell size in Mpc (exact)
	long Nx=0,Ny=0,Nz=0;	// Grid has Nx,Ny,Nz pixels (approx)
	double zref=0,PZerrAxis=0,PzerrReds=0;	// Grid centered at zref (exact) + (Cecile) error on z (redshit or coordinate, it depends)
	double nc=1;			// Mean density of random grid
	// DEBUGGING
	string debug_out;
	bool DoDebug = false;
	bool RandomSeed = false;
	//bool SaveArr = false;
	
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"hrSC:O:a:E:e:r:P:z:m:s:d:h:")) != -1) {
	  switch (c) {
	  case 'C' :
	    input_catalog = optarg;
	    break;
	  case 'O' :
	    out_grids_name = optarg;
	    break;
	  case 'a' :
	    sscanf(optarg,"%lf",&SkyArea);
	    break;
	  case 'E' :
	    sscanf(optarg,"%lf",&PZerrAxis);
	    break;
	  case 'e' :
	    sscanf(optarg,"%lf",&PzerrReds);
	    break;
	  case 'S' :
	    RandomSeed = true;
	    cout << "Montecarlo mode"<< endl;
	    break;
	  case 'r' :
	    isZRadial = true;
	    break;
	  case 'P' :
	    sscanf(optarg,"%ld,%ld,%ld,%lf,%lf",&Nx,&Ny,&Nz,&zref,&R);
	    break; 
	  case 'z' :
	    zcols = optarg; // list of z column names to read in
	    isZColGiven = true;
	    break;
	  case 'm' :
	    sscanf(optarg,"%lf",&nc);
	    break;
	  case 's' :
	    sffiles = optarg;
	    doSFCorr = true;
	    break;
	  case 'd' :
	    debug_out = optarg; // filename of debug files
	    DoDebug = true;
	    break;
	  case 'h' :
	  default :
	    usage(); return -1;
	  }
	}    
	cout << "    - finished decoding command line arguments "<<endl<<endl;
				
	// split up string to read in the two selection function filenames, and possibly 3rd string
	if (doSFCorr) {
		string delim=",";
		vector<string> results;
		stringSplit(sffiles,delim,results);
		sf_file_root = results[0];
		if (results.size()>1) {
		    all_z_file = results[1];
		    doSFCompute = true;
		    }
		if (results.size()>2)
		    isForceZspec = true;
		/*vector<string>::iterator i;
		i = results.begin();
		sf_file_root=*i;
		i++;
		all_z_file=*i;
		i++;
		if (i!=results.end())// if any string given for FORCESZ, sets isForceZspec to true
			{ isForceZspec=true; }*/
		}
		
  	// get up to two z column names
	if (isZColGiven) {
		string delim=",";
		vector<string> results;
		stringSplit(zcols, delim, results);
		ZOCol = results[0];
		if (results.size()>1)
		    ZSCol = results[1];
		/*vector<string>::iterator i;
		i = results.begin();
		ZOCol=*i;
		i++;
		if (i!=results.end())
			ZSCol=*i;*/
		}
	
	// Command line arguments
	
	// IN FILE TYPE
	cout << "     *CATALOG DETAILS*"<<endl;
	cout << "     Galaxy catalog read from file "<< input_catalog;
	cout << "     Reading redshifts from columns named "<< ZOCol <<" and ";
	cout << ZSCol << endl;
	if (isZRadial)
		cout << "     Z dimension IS the radial direction"<<endl;
	cout << endl;
		
	// SF CORRECTION
	cout << "     *SELECTION FUNCTION DETAILS*"<<endl;
	if (doSFCorr) {
		cout << "     Correcting for selection function"<<endl;
		cout << "     Reading true redshifts from "<< all_z_file <<endl;
		string endoffilename;
		if (isForceZspec)
			endoffilename = "_specz_nofz.txt";
		else
			endoffilename = "_nofz.txt";
		cout << "     Saving selection function to "<< sf_file_root << endoffilename <<endl;
		cout << "     (unless "<< sf_file_root << endoffilename <<" already exists)"<<endl;
		}
	else 
		cout << "     Not correcting for selection function ";
	cout <<endl;
		
	// GRID STUFF
	cout << "     *GRID DETAILS*"<<endl;
	cout << "     Full grid defined by ...."<<endl;
	cout << "     pixels : Nx,Ny,Nz = "<< Nx <<","<< Ny <<","<< Nz;
	cout << ", of size "<< R <<" Mpc, centered at z="<< zref <<endl;
	cout << "     Mean density of random grid = "<< nc <<endl;
	cout << endl;
	
	// OUTPUT FILES
	cout << "     *OUTPUT DETAILS*"<<endl;
	//if(SaveArr)
	cout << "     Saving full arrays to filename "<< out_grids_name <<endl;
	if (DoDebug)
		cout << "     Output root filename for debugging is "<< debug_out <<endl;
	cout << endl;
	//-- end command line arguments
  
	int rc = 1;  
	try {  // exception handling try bloc at top level
	
	ResourceUsage res;
	res.Update();
	cout << " Memory size (KB):" << res.getMemorySize() << endl;
	cout << " Max memory size (KB):" << res.getMaxMemorySize() << endl;
	cout << " Maximum allowed data segment size (KB):"<< res.getMaxDataSize() <<endl;
	cout << " Resource usage info : \n" << res << endl;
	
	
    // Read in galaxy catalog
	cout <<"0/ Read in file "<< input_catalog <<endl;
	FitsInOutFile fin(input_catalog, FitsInOutFile::Fits_RO);
	fin.MoveAbsToHDU(2);
	SwFitsDataTable galaxy_catalog(fin, 512, false);
	cout <<endl;
	
	
	res.Update();
	cout << " Initialised SwFitsDataTable"<<endl;
	cout << " Memory size (KB):" << res.getMemorySize() << endl;
	cout << " Max memory size (KB):" << res.getMaxMemorySize() << endl;
	cout << " Maximum allowed data segment size (KB):"<<res.getMaxDataSize()<<endl;
	cout << " Resource usage info : \n" << res << endl;
	
	
	// Set cosmology (read this from galaxy catalog header)
	cout << "     Initialise cosmology: (read from catalog)"<<endl;
	
	//////////  modif Adeline : read cosmo in Fits_RO header
	string H0_s, OmegaM_s, OmegaL_s, OmegaB_s, OmegaR_s, wDE_s, wDA_s, Sigma8_s, Ns_s;
	double h, OmegaM, OmegaL, OmegaB, OmegaR, wDE, wDA, Sigma8, n_s;
	H0_s = fin.KeyValue("H0");
	OmegaM_s = fin.KeyValue("OMEGAM0");
	OmegaL_s = fin.KeyValue("OMEGADE0");
	OmegaB_s = fin.KeyValue("OMEGAB0");
	OmegaR_s = fin.KeyValue("OMEGAR0");
	wDE_s = fin.KeyValue("DE_W0");
	wDA_s = fin.KeyValue("DE_WA");
	Sigma8_s = fin.KeyValue("SIGMA8");
	Ns_s = fin.KeyValue("N_S");

	h = atof(H0_s.c_str()) / 100.;
	OmegaM = atof(OmegaM_s.c_str());
	OmegaL = atof(OmegaL_s.c_str());
	OmegaB = atof(OmegaB_s.c_str());
	OmegaR = atof(OmegaR_s.c_str());
	wDE = atof(wDE_s.c_str());
	wDA = atof(wDA_s.c_str());
	Sigma8 = atof(Sigma8_s.c_str());
	n_s = atof(Ns_s.c_str());

	SimpleUniverse su(h, OmegaM, OmegaL);
	su.SetOmegaBaryon(OmegaB);
	su.SetOmegaRadiation(OmegaR);
	su.SetSigma8(Sigma8);
	su.SetSpectralIndex(n_s);
	su.SetFlatUniverse_OmegaLambda(); // Cecile modif - no be sure that it is flat by adjusting OmegaLambda
	cout << "     OmegaK="<< su.OmegaCurv() <<", OmegaM="<< su.OmegaMatter();
	cout << ", OmegaL="<< su.OmegaLambda() <<", OmegaB="<< su.OmegaBaryon();
	cout << ", Omega_rad=" << su.OmegaRadiation() << ", Omega_cdm=" << su.OmegaCDM() <<", H0="<< su.H0() << endl;
	cout << "check flatness: OmegaTot=" << su.OmegaTotal() << endl;
	if (wDE != -1 or wDA !=0)  
	  su.SetDarkEnergy(su.OmegaLambda(),wDE,wDA);
	
	cout << " and w0=" << su.wDE() << ", wA=" << su.waDE() << ", sigma8=" << su.Sigma8() << endl;
	cout << "Spectral index=" << su.Ns() << endl;
	cout << endl;
	
	// Initialize grid data class
	RandomGenerator rg; // need this for cat2grid
	FitsInOutFile fos(out_grids_name, FitsInOutFile::Fits_Create);
	if (RandomSeed) {
	  rg.AutoInit(0);
	  cout << "Seed automatically generated" << endl;
	} else {
	  long seed=1;
	  rg.SetSeed(seed);
	}

	Cat2Grid cat(galaxy_catalog, su, rg, fos, ZOCol, ZSCol, isZRadial);
	if (DoDebug)
		cat.SetDebugOutroot(debug_out);
	cout << "    The number of gals in whole simulation is "<< cat.ReturnNgAll() <<endl;
	
	if (PZerrAxis > 0.)  cat.SetGaussErrAxis(PZerrAxis,zref,RandomSeed);
	if (PzerrReds > 0.)  cat.SetGaussErrRedshift(PzerrReds,zref,RandomSeed);

	// Compute min and max coordinates and min max redshift
	// only actually really need to do this if correcting for selection function
	cout <<"1/ Find minimum and maximum galaxy cartesian coordinates"<<endl;
	double maxdL; // luminosity distance of largest z in catalog: spec-z or phot-z according to TypePZ
	if (Nx == 0)  {
	  maxdL = cat.FindMinMaxCoords();
	  cout << "    Maximum luminosity distance = "<< maxdL <<endl;
	} else cout << "Nothing to do : grid provided"  <<endl;
	
	res.Update();
	cout << " Computed FindMinMaxCoords()"<<endl;
	cout << " Memory size (KB):" << res.getMemorySize() << endl;
	cout << " Max memory size (KB):" << res.getMaxMemorySize() << endl;
	cout << " Maximum allowed data segment size (KB):"<< res.getMaxDataSize() <<endl;
	cout << " Resource usage info : \n" << res << endl;


	// Set grid
	cout <<"2/ Lay grid over simulation"<<endl;
	cout <<"    Input cell size is "<< R <<endl;
	cat.SetGrid(Nx,Ny,Nz,R,zref);
	
	
	// Selection function
	cout <<"3/ Selection function .... ";
	if(doSFCorr) {
		cout<<endl;

		// 1) set selection function file name

		// first set to "normal" SF, using "observed" redshifts 
		// these redshifts could be spec-z, phot-z, gauss-z		
		string sffile = sf_file_root + "_nofz.txt";
		
		// if making sure the redshifts are SPEC-Z, change SF filename	
		if (isForceZspec) {
			sffile = sf_file_root + "_specz_nofz.txt";
			cout<<"    Using selection function computed from SPECTRO-z";
			cout<<" (no matter what column "<< ZOCol <<" is)"<<endl;
			}
			
		// 2) do selection function computation now or read it from a file
		
		if (doSFCompute) {
		    // both spec-z and phot-z sf's are computed here
		    cat.SaveSelecFunc(sf_file_root, all_z_file);
		    // WARNING! there will be a problem if the z column in the all_z_file is not labeled "z"
            }
        else {
		    ifstream inp;
		    inp.open(sffile.c_str(), ifstream::in);
		    inp.close();
		    if(inp.fail()) { 
			    // sffile does NOT exist
			    string emsg = "ERROR! Selection function in file " + sffile;
			    emsg += " does not exist";
                throw ParmError(emsg);
			    }
		    else {
		        // sffile DOES exist
			    cout <<"    Selection function file has already been computed";
			    cout <<" and will be read from file " << sffile.c_str() <<endl;
		        }
		    }
		    
		// 3) set selection function in Cat2Grid
		    
		// whether sffile was *just* computed or *already* computed it is
		// read into Cat2Grid here
		ComputedSelFunc* sfp = new ComputedSelFunc(sffile); 
		cat.SetSelectionFunction(*sfp);
		}
	else 
		cout <<"    .... not being correcting for "<<endl;
	cout<<endl<<endl;
	
	
	// If debugging just output something here not sure what or why
	if (DoDebug) {
		cat.OutputEuclidCat(SkyArea);
	    }
	else {	
	
	
	// Project galaxies onto the grid
	cout << "4/ Project galaxies onto grid & write to the file ..."<<endl;
	cat.GalGrid(SkyArea);
	cout <<"    - zero size the arrays to save space"<<endl;
	cat.ZeroGalArrays();
	cout << endl;
	
	
	cout <<"    Return actual grid specification:"<<endl;
	TVector<r_8> gridv = cat.ReturnGridSpec();
	cout <<"    Nx,Ny,Nz,L="<< gridv(0) <<","<< gridv(1) <<","<< gridv(2);
	cout <<","<< gridv(3) <<endl;
	cout << endl;
	
	
	// Make random galaxy grid with mean density 
	cout << "5/ Make random catalog galaxy grid ..."<<endl;
	cout <<"    Mean density of random grid = "<< nc <<endl;
	cat.RandomGrid(nc);//,SaveArr);
	res.Update();
	cout << "    Memory size increase (KB):" << res.getDeltaMemorySize() << endl;
	cout << "    Resource usage info : \n" << res << endl;
	
	
	// Write headers to FITS file that contains array
	// include input galaxy catalog file name in header
	cat.WriteHeader(input_catalog);
	
	
	}// end of if not debugging
	
  }  // End of try bloc 
  
  
catch (PThrowable & exc) {
	// catching SOPHYA exceptions
    cerr << " grid_data.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
	}
catch (std::exception & e) {  
    // catching standard C++ exceptions
    cerr << " grid_data.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
	}
catch (...) {  
    // catching other exceptions
    cerr << " grid_data.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
	}
cout << " ==== End of grid_data.cc program  Rc= " << rc << endl;
return rc;	
}
