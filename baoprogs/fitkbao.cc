/**
  * @file  fitkbao.cc
  * @brief Given an input power spectrum + errors fit the BAO scale to "wiggles only"
  *        power spectrum
  *
  *
  * @author Alex Abate
  * Contact: abate@email.arizona.edu
  *
  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <functional>
#include <numeric>
#include <algorithm>

// sophya
#include "sopnamsp.h"
#include "histinit.h"
#include "hisprof.h"
#include "histerr.h"
#include "histos.h"
#include "datatable.h"
#include "fitshdtable.h"
#include "swfitsdtable.h"
#include "fitsarrhand.h"
#include "fiosinit.h"
#include "tarray.h"

// DirectSim
#include "geneutils.h"
#include "cosmocalcs.h"
#include "pkspectrum.h"
#include "fitkbaoscale.h"



void usage(void);
void usage(void) {

	cout << endl<<" Usage: fitkbao [...options...]              " << endl<<endl;
	
	cout << "  Given an input power spectrum + errors fit the BAO scale."<<endl;
	cout << endl;
 
	cout << "  The input power spectrum has already been corrected for  "<<endl;
	cout << "  shot noise, selection, photo-z etc and is supplied to the"<<endl;
	cout << "  program with option -P                                   "<<endl;
	cout << endl;

	cout << "  Method: divide observed power spectrum by a reference    "<<endl;
	cout << "  power spectrum and fit a sine wave described by some     "<<endl;
	cout << "  amplitude, and a characteristic scale. To compute the    "<<endl;
	cout << "  reference power spectrum the redshift of the observed "<<endl;
	cout << "  power spectrum must be supplied with option -z, and the  "<<endl;
	cout << "  values of sigma_8 and the spectral index parameters must "<<endl;
	cout << "  supplied with option -c "<<endl;
	cout << endl;
	
	cout << "  Results are written to files starting with the root name "<<endl;
	cout << "  supplied with the -O option. The chi-square values,      "<<endl;
	cout << "  reference power spectrum, best-fit sinusoid and fit      "<<endl;
	cout << "  results are written to files "<<endl;
	cout << endl;
	
	cout << " -P : PSFile: power spectrum file to read in               "<<endl;
	cout << "              (3 columns: k (Mpc^-1), P(k) Mpc^3, err)     "<<endl;
	cout << " -U : cosmogical model read from grid file                 "<<endl;
	cout << " -O : outfile_root: file root name to write results to     "<<endl; 
	cout << " -s : if input file is simulation (no shotnoise, no sigma) "<<endl; 
	cout << endl;
}



int main(int narg, char *arg[]) {

	SophyaInit();
	FitsIOServerInit();
  
	// input power spectrum
	string ps_file, cosmo_file, ref_file;
	// output file
	string outfile_root;
	double maxk = 1;
	double Sigma8, n_s;
	bool simu_mode = false;

	//--- decoding command line arguments 
	char c;
	while ((c = getopt(narg,arg,"hsP:U:R:O:d")) != -1) {
	    switch (c) {
		    case 'P' :
			    ps_file = optarg;
			    break;
		    case 'U' :
			    cosmo_file = optarg;
			    break;
		    case 'R' :
			    ref_file = optarg;
			    break;
		    case 'O' :
			    outfile_root = optarg;
			    break;
		    case 's' :
			    simu_mode = true;
			    break;
		    case 'h' :
		        default :
			    usage(); return -1;
		    }
	    }
	
	
	cout << "     Printing command line arguments ... "<<endl<<endl;
	cout << "     Reading in observed power spectrum from: "<< ps_file <<endl;
	cout << "     Reading cosmological parameters from: "<< cosmo_file <<endl;
	if (simu_mode) cout << " Will read shot noise and sigmaP from the input power spectrum file" << endl;
	cout << "     Saving results to files beginning "<< outfile_root <<endl;
	cout <<endl;
	
    try {
	
	
	// Read in power spectrum file
	cout << "     Read in power spectrum file "<< ps_file <<endl;
	ifstream ifs(ps_file.c_str());
	TArray<r_8> power_spectrum;
	sa_size_t nr, nc;
	power_spectrum.ReadASCII(ifs,nr,nc);
	cout << power_spectrum ;

	// Set cosmology 
	cout << "     Initialise cosmology:"<<endl;
	FitsInOutFile fin(cosmo_file, FitsInOutFile::Fits_RO);   
	fin.MoveAbsToHDU(4);// assuming cosmo from cat_grid

	//Modif Adeline : read cosmo parameters in file header
	string H0_s, OmegaM_s, OmegaL_s, OmegaB_s, OmegaR_s, wDE_s, wDA_s, Sigma8_s, Ns_s;
	double h, OmegaM, OmegaL, OmegaB, OmegaR, wDE, wDA; // Sigma8 et n_s defined before as they can be given as parameters
	H0_s = fin.KeyValue("H0");
	OmegaM_s = fin.KeyValue("OMEGAM0");
	OmegaL_s = fin.KeyValue("OMEGADE0");
	OmegaB_s = fin.KeyValue("OMEGAB0");
	OmegaR_s = fin.KeyValue("OMEGAR0");
	wDE_s = fin.KeyValue("DE_W0");
	wDA_s = fin.KeyValue("DE_WA");
	Sigma8_s = fin.KeyValue("SIGMA8");
	Ns_s = fin.KeyValue("N_S");
	
	h = atof(H0_s.c_str()) / 100;
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
	su.SetFlatUniverse_OmegaLambda(); // Cecile modif - to be sure that it is flat by adjusting OmegaLambda
	cout << "     OmegaK="<< su.OmegaCurv() <<", OmegaM="<< su.OmegaMatter();
	cout << ", OmegaL="<< su.OmegaLambda() <<", OmegaB="<< su.OmegaBaryon();
	cout << ", Omega_rad=" << su.OmegaRadiation() << ", Omega_cdm=" << su.OmegaCDM() <<", H0="<< su.H0() << endl;
	cout << "check flatness: OmegaTot=" << su.OmegaTotal() << endl;
	if (wDE != -1 or wDA !=0)  
	  su.SetDarkEnergy(su.OmegaLambda(),wDE,wDA);
	
	cout << " and w0=" << su.wDE() << ", wA=" << su.waDE() << ", sigma8=" << su.Sigma8() << endl;
	cout << "Spectral index=" << su.Ns() << endl;
	cout << "Previous values are inputs of the simulation" << endl;
	cout << "____________________________________________" << endl << endl;
	
	// Initialise FitBAOScale
	cout << "     Compute chisq:"<<endl;
	double zref = atof(fin.KeyValue("ZREF").c_str()); 
	FitBAOScale fitbao(power_spectrum, su, zref, ref_file, simu_mode);
	fitbao.ComputeChisq(maxk);

	
	// Find best fit scale and 1-sig error
	cout << "     Find best-fit scale and 1-sig error:"<<endl;
	double bestfit, siglow, sighigh;
	int nsig = 1;
	fitbao.BestfitStdDev(bestfit, siglow, sighigh, nsig);
	double errup = sighigh - bestfit;
	double errdown = bestfit - siglow;
	cout <<"      ka = "<< bestfit <<"+"<< errup <<"-"<< errdown <<endl;
	cout <<endl;
	

	// print info to a file
	cout << "     Print chisq and results to files"<<endl;
	string outfile;
	outfile = outfile_root + "_chisq.txt";
	cout << "     Write chi^2 to file "<< outfile <<endl;
	fitbao.WriteChisq(outfile);
	outfile = outfile_root + "_ancillary.txt";
	cout << "     Write reference power spectrum AND best-fit sinusoid model";
	cout << " to file "<< outfile <<endl;
	fitbao.WriteAncillaryInfo(outfile);
	outfile = outfile_root + "_result.txt";
	cout << "     Write results to file "<< outfile <<endl;
	fitbao.WriteResults(outfile);
	cout << endl;


    }// end of try
  
  
catch(PThrowable exc ) {
    cerr << "fitkbao.cc , Catched exception: \n" << exc.what() << endl;
    }
catch(std::exception ex) {
    cerr << "fitkbao.cc , Catched exception ! " << (string)(ex.what()) << endl;
    }
catch(...) {
    cerr << "fitkbao.cc , Catched ... ! " << endl;
    }

cout << "--------------- fitkbao.cc / END --------------------- " << endl;
}// end of main
