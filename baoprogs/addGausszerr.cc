/**
 * @file  addGausszerr.cc
 * @brief Add Gaussian photo-z errors to distances in a catalog
 *
 *
 * @author Alex Abate
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2009
 * @date 2009
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <fstream>


// SOPHYA
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
#include "ctimer.h"


// DirectSim
#include "geneutils.h"
#include "cat2grid.h"


void usage(void);
void usage(void) {

	cout << endl<<" Usage: addGausszerr [...options...]"           <<endl<<endl;
	
	cout << "  Add a Gaussian redshift error to the galaxies in a       "<<endl;
	cout << "  catalog. The redshift error is added either directly to  "<<endl;
	cout << "  the redshift itself, or to the z-coordinate. If the error"<<endl;
	cout << "  is added to the z-coordinate the redshifts are recomputed"<<endl;
	cout << "  according to the new radial distance to the galaxy. To   "<<endl;
	cout << "  add errors to the z-coordinate instead of the redshift   "<<endl;
	cout << "  coordinate use option -z "<<endl;
	cout << endl;

    cout << "  The size of the error added to the redshift is           "<<endl;
    cout << "  sigma_z*(1+z) which is converted into an equivalent      "<<endl;
    cout << "  comoving distance if it is to be added to the comoving   "<<endl;
    cout << "  distance. The value of sigma_z is supplied to the program"<<endl;
    cout << "  using option -E. The value of 'z' is either the redshift "<<endl;
    cout << "  of the galaxy in question, or the redshift value supplied"<<endl;
    cout << "  to the program using the -Z option                       "<<endl;
	cout << endl;

	cout << "  If the catalog has been simulated with the radial        "<<endl;
	cout << "  dimension parallel to the z-dimension use option -r. This"<<endl;
	cout << "  option has no effect if the error is being added directly"<<endl;
	cout << "  to the redshift.                                         "<<endl;
    cout << endl;

	
	cout << "  EXAMPLE 1: Add a photometric redshift error of size      "<<endl;
	cout << "  0.03*(1+z) to the spectroscopic redshifts in column      "<<endl;
	cout << "  labelled 'zs' in a file called cat.fits and output the   "<<endl;
	cout << "  augmented data to a file called Gausscat.fits            "<<endl;
	cout << endl;
	
	cout << "  $ addGausszerr -C cat.fits -O Gausscat.fits -E 0.03 -c zs"<<endl;
	cout << endl;
	
	cout << "  EXAMPLE 2: As example 1, but add the photometric redshift"<<endl;
	cout << "  error to the z-coordinate (comoving distance in z-dim)   "<<endl;
	cout << "  instead, and make this error constant for all galaxies,  "<<endl;
	cout << "  defined as: sigma_z*(1+zref), where zref=1.              "<<endl;
	cout << endl;
	
	cout << "  $ addGausszerr -C cat.fits -O Gausscat.fits -E 0.03 -c zs"<<endl;
	cout << "                 -z -Z 1.                                  "<<endl;
	cout << endl;
	
	cout << " -C : CatName     FITS filename containing catalog         "<<endl;
	cout << " -O : OutCatName  FITS file containing output catalog with "<<endl;
	cout << "                  Gaussian z errors                        "<<endl;
	cout << " -E : PZerr       Size of photometric redshift error:      "<<endl;
	cout << "                  PZerr*(1+zs) [DEFAULT=0.03]              "<<endl;
	
	cout << " -z : [noarg]     Add error to z-coordinate instead of     "<<endl; 
	cout << "                  redshift [DEFAULT=no]                    "<<endl;
	cout << " -Z : zref        Add constant redshift error with value   "<<endl; 
	cout << "                  PZerr*(1+zref) [DEFAULT=no]              "<<endl;
	cout << " -c : ZSCol       Name of column of spec-z                 "<<endl;
	cout << " -r : [noarg]     z-dimension IS radial direction in       "<<endl;
	cout << "                  catalog [DEFAULT=no]. Has no effect if   "<<endl;
	cout << "                  adding error directly to the redshift    "<<endl;
	cout << endl;
}


int main(int narg, char *arg[]) {

  SophyaInit();
  FitsIOServerInit();
  
  cout << " ==== setting defaults ===="<<endl;
  string InCat, OutCat;
  double PZerr = 0.03;
  bool ADD2Z = false, RadialZ = false;
  double zref = -1; // becomes a value >0 if -Z option used
  string ZCol; 
  bool Zsp = false;
  string ZSCol = "zs"; // by default SPECTRO redshift column labelled "zs" is read in
  
  //--- decoding command line arguments 
  cout << " ==== decoding command line arguments ===="<<endl;
  char c;
  while((c = getopt(narg,arg,"hrzC:O:E:Z:c:")) != -1) {
    switch (c) {
    case 'C' :
      InCat = optarg;
      break;
    case 'O' :
      OutCat = optarg;
      break;
    case 'E' :
      sscanf(optarg,"%lf",&PZerr);
      break;
    case 'r' :
      RadialZ = true;
      break;
    case 'z' :
      ADD2Z = true;
      break;
    case 'Z' :
      //ADD2Z = true;	
      sscanf(optarg,"%lf",&zref);
      break;
    case 'c' :
      ZSCol = optarg; // z column names to read in
      //Zsp = true;
      break;
    case 'h' :
    default :
      usage(); return -1;
    }
  }
  
  //	// get two z column names
  //	if (Zsp)	
  //		{ 
  //		string delim=",";
  //		vector<string> results;
  //		stringSplit(ZCol,delim,results);
  //		vector<string>::iterator i;
  //		i = results.begin();
  //		ZOCol=*i;
  //		i++;
  //		if (i!=results.end())
  //			ZSCol=*i;
  //		}
  
  cout << "     Printing command line arguments ...             "<<endl<<endl;
  
  cout << "     Reading in observed catalogs "<< InCat                 <<endl;
  cout << "     Reading in spec redshifts from column "<< ZSCol        <<endl;
  if (RadialZ) 
    cout << "     Z dimension IS radial direction                   "<<endl;
  cout << "     Catalog with added Gaussian redshift errors will be   "<<endl;
  cout << "     written to "<< OutCat                                  <<endl;
  cout << "     Adding Gaussian photometric redshift error =          "<<endl;
  cout << "     "<< PZerr <<"(1+zs)                                   "<<endl;
  cout << "     Adding error to ";
  if (ADD2Z) {
    cout << "z-coordinate"                                           <<endl;
  }
  else 
    cout << "spectroscopic redshift"<<endl;
  if(zref>0)
    cout << "     - but only equivalent to error at z = "<< zref <<endl;
  cout <<endl;
  
  
  try {
    
    Timer tm("AddGaussZerr");
    
    RandomGenerator rg; // random number generator
    long seed=1;
    rg.SetSeed(seed);
    
    // Read input catalog
    FitsInOutFile fin(InCat,FitsInOutFile::Fits_RO);
    fin.MoveAbsToHDU(2);
    SwFitsDataTable dt(fin,512,false);
    sa_size_t ng = dt.NEntry(); // number of galaxies
    sa_size_t nc = dt.NCols();  // number of columns
    DataTableRow row = dt.EmptyRow();
    cout << "     Number of galaxies in catalog = "<< ng ;
    cout << ", number of columns in catalog = "<< nc <<endl;
    
	     		
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
	      		
    // Get fixed error in comoving z-dimension (if using)
    double ErrFix = 0;
    if (zref>0 && ADD2Z) {
      cout << "     Calculating co-moving distance error at z = "<<zref;
      //su.SetEmissionRedShift(zref);
      ErrFix = ZErr2CoDistErr(su,PZerr,zref);  // change by Reza, 08/07/2014, then Cecile 30/04/15
      // PZerr*(1+zref)*(SpeedOfLight_Cst/su.HZE());
      cout <<", dDc = "<< ErrFix <<endl<<endl;
    }
    
    
    
    // Create z-dcomoving look up table incase needed
    SInterp1D z2dist, dist2z;
    int_8 nz=1000;
    vector<double> zrs, codist;
    double minz=0, maxz=10;
    double dz = (maxz-minz)/(nz-1);
    for(int kk=0; kk<nz; kk++) {
      
      double zs=minz+kk*dz;
      su.SetEmissionRedShift(zs);
      double cod =su.RadialCoordinateMpc(); // radial distance 
      zrs.push_back(zs);
      codist.push_back(cod); 
    }
    z2dist.DefinePoints(zrs,codist,minz,maxz,2*nz);
    dist2z.DefinePoints(codist,zrs,codist[0],codist[codist.size()-1],2*nz);
    
    sa_size_t nobj = 20; // number of objects to print to the screen
    sa_size_t nskip = (sa_size_t)ng/nobj;
    long j=0;
    if (ADD2Z) 
      cout << "     Printing every "<< nskip <<" galaxy to screen "<<endl;

    
    // Create new catalog 
    cout << "     Creating new catalog called "<< OutCat <<endl;
    FitsInOutFile swf(OutCat, FitsInOutFile::Fits_Create);	
    SwFitsDataTable gals(swf, 2048);
    for (int i=0;i<nc;i++) { // add same columns input file has
      
      cout <<"  Col number = "<< i <<", Col name = "<< dt.NomIndex(i);
      cout <<", Col type = "<< dt.GetColumnType(i) <<endl;
      if (dt.GetColumnType(i)<2)
	gals.AddIntegerColumn(dt.NomIndex(i));
      else if ( (dt.GetColumnType(i)>1) && (dt.GetColumnType(i)<3) )
	gals.AddLongColumn(dt.NomIndex(i));
      else if ( (dt.GetColumnType(i)>2) && (dt.GetColumnType(i)<4) )
	gals.AddFloatColumn(dt.NomIndex(i));
      else if ( (dt.GetColumnType(i)>3) && (dt.GetColumnType(i)<5) )
	gals.AddDoubleColumn(dt.NomIndex(i));
      else if ( (dt.GetColumnType(i)>4) && (dt.GetColumnType(i)<6) )
	gals.AddComplexColumn(dt.NomIndex(i));
      else if ( (dt.GetColumnType(i)>5) && (dt.GetColumnType(i)<7) )
	gals.AddDoubleComplexColumn(dt.NomIndex(i));
      else if ( (dt.GetColumnType(i)>6) && (dt.GetColumnType(i)<8) )
	gals.AddStringColumn(dt.NomIndex(i));
			}
    // then add one more that will contain z's w/ Gaussian errors
    gals.AddFloatColumn("zG");
    DataTableRow rowout = gals.EmptyRow();
    
    
    
    // Find required col names (only required if adding to z-dimension)
    sa_size_t Izs = dt.IndexNom(ZSCol);
    sa_size_t Ith,Iph;
    if(!RadialZ) {
      Iph = dt.IndexNom("phi");
      Ith = dt.IndexNom("theta");
    }
    
    // loop over catalog
    cout << "     Starting loop over catalog ... "<<endl;		
    
    for (long i=0;i<ng;i++) {
      
      dt.GetRow(i,row);
      double zs = row[Izs]; // spec-z of galaxy 
      
      
      double zG, th, ph;	
      // if adding error to z-coordinate instead of redshift
      if (ADD2Z) {
	
	// for printing info
	int dij1= i-j;
	int dij = abs(dij1);
	
	// get cartesian galaxy position
	double x,y,z,reds;
	Cat2Grid cat(dt,su,rg,dist2z,z2dist,ZSCol,ZSCol,RadialZ);
	GalRecord grec;
	cat.Row2Record(row,grec);
	cat.Rec2EuclidCoord(grec,x,y,z,reds);
	
	
	if(i<10) {
	  cout <<"    galid="<<i<<": ";
	  grec.Print();
	}
	
	double PZDerr = 0;
	double zold=z;
	
	
	if (zref<0) {
	  // add varying error to z coord
	  //  su.SetEmissionRedShift(reds);
	  //Reza-DEL					PZDerr = su.ZErr2CoDistErr(PZerr);
	  PZDerr = ZErr2CoDistErr(su,PZerr,reds);  // change by Reza, 08/07/2014, then Cecile 30/04/15
	  z = z + PZDerr*rg.Gaussian();
	}
	else // add constant error to z coord
	  z = z + ErrFix*rg.Gaussian();	  
	
	// convert new z-coord back to redshift
	double dc;
	if (RadialZ)
	  dc = z;
	else
	  dc = sqrt(x*x+y*y+z*z);
	zG = dist2z(dc);
	
	// need to also convert th and ph 
	// if not radial direction)
	ph = atan2(y,x); // should be same as old ph
	th = acos(z/dc);

	if (i<10)
	  cout <<", zpG="<< zG <<endl;
	
	if (dij<1) {
	  cout <<"check values: ph="<< ph <<", th="<< th <<endl;
	  cout <<"gal "<< i <<": (x,y,z)=("<< x <<","<< y <<","<< z;
	  cout <<"), zold="<< zold <<", sig_zc="<< PZDerr <<", dc=";
	  cout << dc <<", zs="<< zs <<", zG="<< zG << " diff = " << fabs(zs-zG) << endl;
	  tm.SplitQ();
	  cout <<"Loop "<< i <<", elapsed time ";
	  cout << tm.TotalElapsedTime()<<"s"<<endl;
	  j+=nskip;
	  cout <<"j is now "<< j <<endl;
	}
	
      } // end of case error on z-coordinate

      else { // If adding straight to spec-z
	
	if (zref<0) {
	  double sigz = PZerr*(1+zs); // error at z of galaxy
	  zG = zs + sigz*rg.Gaussian();
	}
	else {
	  double sigz = PZerr*(1+zref); // error at z specified
	  zG = zs + sigz*rg.Gaussian();
	}
	
	// don't think I need to convert the th and ph					
      }
      
      // make sure there are no unphysical redshifts
      if (zG<0)
	zG=0.01;
      
      // fill output row with same values as input row
      for (int n=0; n<row.Size();n++)
	rowout[n]=row[n];
      // add 'redshift' with error to last column
      rowout[row.Size()] = zG;
      
      // recomputed theta and phi if error was added to z-dimension
      if(ADD2Z&&!RadialZ) {
	rowout[Iph] = ph;
	rowout[Ith] = th;
      }
      
      gals.AddRow(rowout);
      
    }// end loop over file
    
    //modified by Adeline : write cosmo parameters in file header
    swf.WriteKey("H0", su.H0()," Cosmo.Param H0");
    swf.WriteKey("OMEGAM0", su.OmegaMatter()," Cosmo.Param OmegaMatter0 ");
    swf.WriteKey("OMEGAB0", su.OmegaBaryon()," Cosmo.Param OmegaBaryon0");
    swf.WriteKey("OMEGAR0", su.OmegaRadiation()," Cosmo.Param OmegaRadiation0");
    swf.WriteKey("OMEGAT0", su.OmegaTotal()," Cosmo.Param OmegaTot0");
    swf.WriteKey("OMEGADE0", su.OmegaLambda(),"  Cosmo.Param OmegaLambda0 (dark energy density)");
    swf.WriteKey("OMEGADK", su.OmegaCurv(),"  Cosmo.Param OmegaK ");
    swf.WriteKey("DE_W0", su.wDE(), " Cosmo.Param w0 (dark energy eq.state)");
    swf.WriteKey("DE_WA",su.waDE() , " Cosmo.Param wA (dark energy eq.state)"); 
    swf.WriteKey("SIGMA8", su.Sigma8(), " Cosmo.Param sigma8_0");
    cout << "Check cosmo parameters : " << endl;
    cout << "  OmegaK="<< su.OmegaCurv() <<", OmegaM="<< su.OmegaMatter();
    cout << ", OmegaL="<< su.OmegaLambda() <<", OmegaB="<< su.OmegaBaryon()  ;
    cout << ", Omega_rad=" << su.OmegaRadiation() << ", H0=" << su.H0() << ", Sig8=" << su.Sigma8() <<endl; 
    cout << ", Omega_curv=" << su.OmegaCurv() << ", DE_W0=" << su.wDE() << ", DE_WA=" << su.waDE() <<endl; 
    cout << endl;
    // end modifications
    			
    }
    catch(PThrowable exc ) {
        cerr << "addGausszerr.cc , Catched exception: \n" << exc.what() << endl;
        }
    catch(std::exception ex) {
        cerr << "addGausszerr.cc , Catched exception ! " << (string)(ex.what()) << endl;
        }
    catch(...) {
        cerr << "addGausszerr.cc , Catched ... ! " << endl;
        }

    cout << "--------------- addGausszerr.cc / END --------------------- " << endl;
}
