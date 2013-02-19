#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include <typeinfo>
#include "timing.h"

#include "array.h"
#include "hisprof.h"
#include "histerr.h"
#include "histinit.h"
#include "fiosinit.h"
#include "fitsioserver.h"
#include "fftwserver.h"
#include "mydefrg.h"
#include "resusage.h"

#include "mass2gal.h"
#include "schechter.h"
#include "cosmocalcs.h"
#include "gftdist.h"
#include "constcosmo.h"

void usage(void);
void usage(void){

	cout << endl<<" Usage: rdlss [...options...]" << endl<<endl;

	cout << "    Program for simulating galaxy catalog from SimLSS output." << endl;
 	cout << "    It simulates galaxy positions, absolute magnitude and galaxy" << endl;
 	cout << "    SED type.  " << endl<<endl;
 
	cout << "    There are 4 possible outputs after simulating the galaxy catalog:" << endl;
	cout << "    1) Output galaxy catalog, either using galaxy clustering or not [DEFAULT]" << endl;
	cout << "    2) Output a HISTO of redshifts [OPTION -H]" << endl;
	cout << "    3) Output catalog of just (true) redshifts [OPTION -Z]" << endl;
	cout << "    4) Simple simulation: ra,dec,z only [OPTION -S]" << endl<<endl;
	
	cout << "    EXAMPLE 1: Simulate a galaxy catalog from an over-density distribution" << endl;
	cout << "    in file odens.fits.  We want the survey to have a circular sky area with "<<endl;
	cout << "    radius = pi/4 radians, and galaxies with randomized positions in the pixels."<<endl;
	cout << "    So not to produce too many galaxies so the absolute magnitude cut is set,"<<endl;
	cout << "    This throws out any galaxies with an unobservable absolute magnitude given"<<endl;
	cout << "    their redshifts. This is simulation number 1 so the id is set to be 1."<<endl;
	cout << "    The output catalog will be saved to out.fits:";
	cout <<endl<<endl;
	
	cout << "    $ rdlss -C odens.fits -O out.fits -a 0.7854 -i 1 -R -M"<<endl<<endl;
	
	cout << "    EXAMPLE 2: Similar to example 1, but write only the true redshifts"<<endl;
	cout << "    to the galaxy catalog and DON'T set the absolute magnitude cut:"<<endl;
	cout <<endl;
	
	cout << "    $ rdlss -C odens.fits -O ztrue.fits -i 1 -R -Z"<<endl<<endl<<endl;

	cout << " -O : OutputFileName: filename that the output is written to"<<endl;
	cout << " -a : SkyArea: Radius of survey sky area in radians (assumed circular),"<<endl;
	cout << "      [default=2PI]"<<endl;
	cout << " -C : MassDistFileName: FITS filename containing SimLSS output"<<endl;
	cout << " -H : Only produce a histogram of the TRUE catalog n(z) [default=NO]"<<endl;
	cout << " -i : Simulation identifier, integer >=1 [default=1]"<<endl;
	cout << " -M : Don't add galaxies with M>AbsMagCut(z) [default=YES]"<<endl; 
	cout << " -N : OutRoot: Output mass cube, ngals cube to files [OutRoot]_ngals.fits"<<endl;
	cout << "      and [OutRoot]_mass.fits [default=NO,OutRoot=tmp]"<<endl;
	cout << " -R : Randomise galaxy positions within cell [default=NO]"<<endl;
	cout << " -S : Only do a simple sim of ra,dec,z [default=NO]"<<endl;
	cout << " -s : ngal: Another type of simple simulation. Simulate ngal galaxies"<<endl;
	cout << "      in each cell [default=NO,ngal=1]"<<endl;
	cout << " -x : xplanes: Number of SimLSS cube planes to remove [default=1]"<<endl;
	cout << " -z : z dimension is 'radial' direction [default=NO]"<<endl;
	cout << " -Z : Only produce a catalog containing the TRUE redshifts (and galaxy IDs)"<<endl;
	cout << "      [default=NO]"<<endl;
	
	cout << endl;
}

int main(int narg, char* arg[])
{
	cout << " ==== rdlss.cc program , reading SimLSS d_rho/rho output  ==== " << endl;
	// Make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	InitTim();
	cout<<endl<<endl;

	cout << " ==== setting defaults ===="<<endl;
	string infile,outfile,outroot="tmp";  
	//double absmc = -13;
	double skyarea=6.3;
	bool CutArea=false;   // throw out gals with phi>skyarea
	bool CutM=false;      // do a cut on absolute magnitude(z) 
	bool OutNgal=false;   // write out ngals, mass cubes
	bool RandPos=false;   // randomize galaxy positions in cells
	bool HistoOnly=false; // only histogram of z
	bool TrueZOnly=false; // only z
	bool SimpleSim=false; // just ra,dec,z
	bool SimpleSim2=false;// no clustering
	bool ZisRadial=false; // z dimension IS radial direction
	int idsim=1;
	int xplanes=1; // SimLSS simulates cube with 1 or 2 extra planes: one dimension is too long
	// xplanes should be the difference between NX and drho.SizeX()
	float ngalcell=1;
	float conv;

  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"hHMRSZza:C:i:N:O:x:s:")) != -1) {
	    switch (c) {
            case 'C' :
                infile = optarg;
                break;
            case 'i' :
                sscanf(optarg,"%d",&idsim);
                break;
            case 'O' :
                outfile	= optarg;
                break;
            case 'M' :
                CutM=true;
                break;
            case 'a' :
                sscanf(optarg,"%lf",&skyarea);
                CutArea = true;
                break;
            case 'N' :
                OutNgal = true;
                outroot = optarg;
                break;
            case 'R' :
                RandPos=true;
                break;
            case 'H' :
                HistoOnly=true;
                break;
            case 'S' :
                SimpleSim=true;
                break;
            case 's' :
                sscanf(optarg,"%f",&ngalcell);
                SimpleSim2=true;
                break;
            case 'x' :
                sscanf(optarg,"%d",&xplanes);
                break;
            case 'z' :
                ZisRadial=true;
                break;
            case 'Z' :
                TrueZOnly=true;
                break;
            case 'h' :
                default :
                usage(); return -1;
		    }
	    }

	if(TrueZOnly&&HistoOnly)
		throw ParmError("ERROR! Cannot have both options to output n(z) histo AND z file");
		
	cout << "     Reading command line arguments ... "<<endl;
	cout << "     SimLSS file is "<<infile<<endl;
	cout << "     Removing "<<xplanes<<" planes from SimLSS cube"<<endl;
	if(!HistoOnly&&!TrueZOnly&&!SimpleSim) {
		cout << "     Output catalog file will be "<<outfile<<endl;
		cout << "     Simulation ID number is "<<idsim<<"00,000,000,000"<<endl; 
		if (CutM)
			cout << "     Removing galaxies with faint absolute magnitudes"<<endl;
		if (ZisRadial)
			cout << "     z-dimension is radial dimension"<<endl;
		if (SimpleSim2)
			cout << "     Catalog simulated with a constant number density of "<<ngalcell<<endl;
		}
	else if (HistoOnly)
		cout << "     Outputting histogram of n(z) ONLY into file "<<outfile<<endl; 
	else if (TrueZOnly)
		cout << "     Outputting true redshifts ONLY into file "<<outfile<<endl; 
	else if (SimpleSim)
		cout << "     Outputting ra, dec, true redshifts ONLY into file "<<outfile<<endl; 
	if (CutArea)
		cout << "     Keeping galaxies with phi<"<<skyarea<<endl; 
	if (RandPos)
		cout << "     Randomising galaxy positions within pixel"<<endl;
  //-- end command line arguments
  
  int rc = 1;  
  try {  // exception handling try bloc at top level
  
	ResourceUsage res;
	res.Update();
	cout << " Memory size (KB):" << res.getMemorySize() << endl;
	cout << " Max memory size (KB):" << res.getMaxMemorySize() << endl;
	cout << " Maximum allowed data segment size (KB):"<<res.getMaxDataSize()<<endl;
	cout << " Resource usage info : \n" << res << endl;
  
	cout << "0/ rdlss.cc: Reading input file= " << infile << endl;  
	FitsInOutFile fin(infile,FitsInOutFile::Fits_RO);
	TArray<r_8> drho;
	fin >> drho;
	cout << drho.Info();
	cout << "    print original drho array size: "<<drho.SizeX()<<"x"<<drho.SizeY()<<"x"<<drho.SizeZ()<<endl<<endl<<endl;
   
	cout << "0.1/ Initialise cosmology: (same as SimLSS)"<<endl;
	double h=0.71, OmegaM=0.267804, OmegaL=0.73;
	SimpleUniverse su(h,OmegaM,OmegaL);
	su.SetFlatUniverse_OmegaMatter();
	cout <<"    OmegaK="<<su.OmegaCurv()<<", OmegaM="<<su.OmegaMatter()<<", OmegaL="<<su.OmegaLambda();
	cout <<", OmegaB="<<su.OmegaBaryon()<<", H0="<<su.H0()<<endl;
	
	/* remove N extra planes */
	cout << "1/ Initialise Mass2Gal: remove planes" << endl;
	RandomGenerator rg;
	Mass2Gal m2g(drho,su,rg,xplanes);
	double mean, sig;
	TArray<r_8> mass;
	m2g.MassArray(mass);
	MeanSigma(mass, mean, sig);
	cout << endl<<"1.1/ RAW DENS CUBE STATS: Mean=" << mean << " Sigma=" << sig << endl<<endl;
	res.Update();
	cout << " Memory size (KB):" << res.getMemorySize() << endl;                                                                 
	cout << " Resource usage info : \n" << res << endl;
	
	/* read in cube and pixel properties from fits header */
	cout << "1.2/ Read in cube properties from fits header" << endl;
	m2g.ReadHeader(fin);
	// xplanes should be the difference between NX and drho.SizeX()
	int NZ=m2g.ReturnNZ();
	int diff = drho.SizeX()-NZ;
	if( xplanes!=abs(diff) )
		throw ParmError("ERROR: removed wrong number of planes from SimLSS cube");
	m2g.SetRandPos(RandPos);

	/* deal with negative mass cells AND drho/rho -> rho/rho^bar */
	cout << " 2/ Clean Negative Mass Cells" << endl;
	sa_size_t nbadc = m2g.CleanNegativeMassCells(); // adds 1 to drho and sets anything <0=0
	m2g.MassArray(mass);
	cout << "    NBadCells=" << nbadc 
	     << " BadPercentage=" << (double)nbadc*100./(double)mass.Size() << "%"<<endl;
	MeanSigma(mass, mean, sig);
	cout << endl<<"2.1/ CLEANED DENS CUBE STATS: Mean=" << mean << " Sigma=" << sig << endl;
	/* double check there are no bad cells */
	sa_size_t nbadc2 = m2g.CheckNegativeMassCells();
	cout <<"    double check there are no bad cells ...."<<nbadc2<<endl;
	cout <<"    check minimum value in mass array is 0"<<endl;
	double min, max;
	mass.MinMax(min,max);
	cout <<"    min of mass = "<<min<<", max of mass = "<<max<<endl<<endl<<endl;
	res.Update();
	cout << " Memory size (KB):" << res.getMemorySize() << endl;
	cout << " Resource usage info : \n" << res << endl;
	
	
	/* convert to mean Ngal by integrating Schechter function */
	cout << "3/ Convert rho/rho^bar To Mean NGal"<<endl;  

	/* Create and setup the galaxy type, magnitude distribution definition object */
	cout << "3.1/ Set up Schechter functions for each galaxy type"<<endl;
	cout <<" ... GOODS B band: Early types, Late types, Starbursts"<<endl;
	cout <<" ... see Table 3 in Dahlen et al 2005"<<endl;
	
	string LFplace;
	char * plf=getenv("SIMBAOLF");
	if (plf==NULL)
		{
		cout <<"ERROR LF LOCATION ENVIRONMENT VARIABLE NOT DEFINED"<<endl;
		return 1;
		}
	else
		{
		LFplace=plf;
		cout <<"    Location of LF file is "<<LFplace<<endl;
		}
	string LFfile=LFplace +	"GOODS_B_LF.txt";// add an option for this

	//******* all this below should be in Schechter class ***********//
	ifstream ifs;
	ifs.open(LFfile.c_str(), ifstream::in);
	if (ifs.fail())
		cout <<"  ERROR: failed to find luminosity function file"<<endl;
	//ifs.open(LFfile);
	TArray<r_4> LFTable;
	sa_size_t nr, nc;
	LFTable.ReadASCII(ifs,nr,nc);
	//cout << LFTable ;
	
	int MstarCol=2, AlphaCol=3, PhiStarCol=4;
	// ALL GALAXIES
	double MstarAz1=LFTable(MstarCol,13),alpAz1=LFTable(AlphaCol,13),phistarAz1=LFTable(PhiStarCol,13)*1e-4;
	double MstarAz2=LFTable(MstarCol,14),alpAz2=LFTable(AlphaCol,14),phistarAz2=LFTable(PhiStarCol,14)*1e-4;
	double MstarAz3=LFTable(MstarCol,15),alpAz3=LFTable(AlphaCol,15),phistarAz3=LFTable(PhiStarCol,15)*1e-4;
	
	// EARLY TYPES
	double MstarEz1=LFTable(MstarCol,1),alpEz1=LFTable(AlphaCol,1),phistarEz1=LFTable(PhiStarCol,1)*1e-4;
	double MstarEz2=LFTable(MstarCol,6),alpEz2=LFTable(AlphaCol,6),phistarEz2=LFTable(PhiStarCol,6)*1e-4;
	double MstarEz3=LFTable(MstarCol,10),alpEz3=LFTable(AlphaCol,10),phistarEz3=LFTable(PhiStarCol,10)*1e-4;
	
	// LATE TYPES
	double MstarLz1=LFTable(MstarCol,3),alpLz1=LFTable(AlphaCol,3),phistarLz1=LFTable(PhiStarCol,3)*1e-4;
	double MstarLz2=LFTable(MstarCol,7),alpLz2=LFTable(AlphaCol,7),phistarLz2=LFTable(PhiStarCol,7)*1e-4;
	double MstarLz3=LFTable(MstarCol,11),alpLz3=LFTable(AlphaCol,11),phistarLz3=LFTable(PhiStarCol,11)*1e-4;
	
	// STARBURST TYPES
	double MstarSz1=LFTable(MstarCol,4),alpSz1=LFTable(AlphaCol,4),phistarSz1=LFTable(PhiStarCol,4)*1e-4;
	double MstarSz2=LFTable(MstarCol,8),alpSz2=LFTable(AlphaCol,8),phistarSz2=LFTable(PhiStarCol,8)*1e-4;
	double MstarSz3=LFTable(MstarCol,12),alpSz3=LFTable(AlphaCol,12),phistarSz3=LFTable(PhiStarCol,12)*1e-4;
	
	string MstarUnits="M-5log10h70";
	string phistarUnits="(Mpc/h70)^-3";
	
	cout <<" ... Schechter parameters"<<endl; 
	cout <<"     z range     Mstar     alpha     phistar     spec type"<<endl;
	cout <<"    0.75<z<1.0  "<<MstarAz3<<"     "<<alpAz3<<"       "<<phistarAz3<<"        All"<<endl;
	cout <<"                "<<MstarEz3<<"     "<<alpEz3<<"       "<<phistarEz3<<"         Early"<<endl;
	cout <<"                "<<MstarLz3<<"     "<<alpLz3<<"       "<<phistarLz3<<"          Late"<<endl;
	cout <<"                "<<MstarSz3<<"     "<<alpSz3<<"        "<<phistarSz3<<"        Starburst"<<endl<<endl;
	//******* all this below should be in Schechter class, or other class ***********//
	
	cout<<"3.2/ Mass to Galaxy number conversion"<<endl;
	Schechter schAz3(phistarAz3,MstarAz3,alpAz3);
	//schechter functions for each type
	Schechter schEz3(phistarEz3,MstarEz3,alpEz3);
	Schechter schLz3(phistarLz3,MstarLz3,alpLz3);
	Schechter schSz3(phistarSz3,MstarSz3,alpSz3);
	
	double schmin=-24, schmax=-13;// units of "M-5log10h70"
	int schnpt=10000;
	cout <<"    ... integrating from Mbright="<<schmin<<" "<<MstarUnits<<" to Mfaint="<<schmax<<" "<<MstarUnits<<", with step="<<schnpt<<endl;
	schAz3.SetInteg(schmin,schmax,schnpt);
	double nz3=schAz3.Integrate();
	cout <<"    ... number density of galaxies: "<<nz3<<" Mpc^-3"<<endl;
	
	double pixVol=m2g.ReturnPixVol();
	cout <<"    pixel volume="<< pixVol<<" Mpc^3"<<endl;
	conv = pixVol*nz3;
	cout <<"    actual gals per pixel vol="<<pixVol*nz3<< endl;
	cout <<"    gals per pixel vol="<<conv<< endl;
	
	if(SimpleSim2)
		conv = ngalcell;
	
	m2g.ConvertToMeanNGal(conv,SimpleSim2); //just multiplies mass_ by conv
	res.Update();
	cout << " Memory size (KB):" << res.getMemorySize() << endl;
	cout << " Resource usage info : \n" << res << endl;
	
	/* write out ngals, mass cubes? */
	if(OutNgal)
		{
		cout <<"    **** Writing out cleaned mass distribution AND ngal distribution ****"<<endl;
		TArray<int_8> ngals;
		m2g.NGalArray(ngals);
		TArray<r_8>   mass2;
		m2g.MassArray(mass2);
		cout <<"    check minimum value in mass array is 0"<<endl;
		mass2.MinMax(min,max);
		cout <<"    min of mass = "<<min<<", max of mass = "<<max<<endl;
		string outngal = outroot +"_ngals.fits";
		FitsInOutFile fos(outngal, FitsInOutFile ::Fits_Create);
		fos << ngals;
		cout <<"    Written ngals array"<<endl;
		
		// have to do below stuff or there is a problem with 
		// writing FT'd array to a FITS file
		// now there shouldn't be a problem because of 
		// doing .PackElements() in constructor
		//TVector<r_8> gridv=m2g.ReturnGridSpec();
		//TArray<r_8> massn;
		//int ndim=3;
		//sa_size_t mydim[ndim];
		//mydim[0]=gridv(0); mydim[1]=gridv(1); mydim[2]=gridv(2);
		//massn.SetSize(ndim, mydim);
		//cout <<"    Size of mass array = "<<
		//for(sa_size_t iz=0; iz<mass2.SizeZ(); iz++) 
		//	for(sa_size_t iy=0; iy<mass2.SizeY(); iy++) 
		//		for(sa_size_t ix=0; ix<mass2.SizeX(); ix++) 
		//			massn(ix, iy, iz) = mass2(ix, iy, iz);
		//cout <<" here2"<<endl;
		string outmass = outroot +"_mass.fits";
		FitsInOutFile fos1(outmass, FitsInOutFile ::Fits_Create);
		fos1 << mass2;
		cout <<"    Written mass array"<<endl;
		cout << endl;
		}
		
	/* create distribution */
	cout <<"4/ Set up Mb-Type 2D distribution"<<endl;
	// set distribution parameters
	int magbin=10000;
	int PrtLevel = 2;
	string IntLFUnits="(Mpc/h70)^-3";
	cout <<"    Renormalise type-specific LFs ..."<<endl;
	//get prob distribution
	// - set schechter distributions
	int type;
	type=1;
	SchechterDist schDistEz3(schAz3,schEz3,schLz3,schSz3,type);
	type=2;
	SchechterDist schDistLz3(schAz3,schEz3,schLz3,schSz3,type);
	type=3;
	SchechterDist schDistSz3(schAz3,schEz3,schLz3,schSz3,type);
	
	cout <<"    Find fraction of each type ..."<<endl;
	// integrate work out type fractions
	double nfE3=schDistEz3.Integrate(schmin,schmax,schnpt);
	double nfL3=schDistLz3.Integrate(schmin,schmax,schnpt);
	double nfS3=schDistSz3.Integrate(schmin,schmax,schnpt);
	double totalnr=nfE3+nfL3+nfS3;
	double fEz3=nfE3/totalnr;
	double fLz3=nfL3/totalnr;
	double fSz3=nfS3/totalnr;
	cout <<"    Type fractions from renormalised LFs are: fE="<<fEz3<<", fL="<<fLz3<<", fS="<<fSz3<<endl;
	GalFlxTypDist gfdz3(rg, PrtLevel); 
	gfdz3.AddGalType(schDistEz3, schmin, schmax, fEz3, magbin, schnpt); 
	gfdz3.AddGalType(schDistLz3, schmin, schmax, fLz3, magbin, schnpt); 
	gfdz3.AddGalType(schDistSz3, schmin, schmax, fSz3, magbin, schnpt);
	
	/* cut magnitude? */
	if (CutM)
		m2g.MaxAbsMag();// this does calculation of what max abs mag is as a function z
	
	/* sim galaxies and write out */
	cout <<"5/ Simulate galaxies"<<endl;
	bool extinct=false;
	
	if(!HistoOnly&&!TrueZOnly&&!SimpleSim)
		m2g.CreateGalCatalog(idsim,outfile,gfdz3,extinct,CutM,skyarea,ZisRadial);
	else if (HistoOnly)
		m2g.CreateNzHisto(outfile,gfdz3,extinct,skyarea);
	else if (TrueZOnly)
		m2g.CreateTrueZFile(idsim,outfile,skyarea);
	else if (SimpleSim)
		m2g.CreateSimpleCatalog(idsim,outfile,skyarea);
		
  }  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " rdlss.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " rdlss.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " rdlss.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ==== End of rdlss.cc program  Rc= " << rc << endl;
  return rc;	
}
