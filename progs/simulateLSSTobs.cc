#include <unistd.h>

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
#include "timestamp.h"
#include "ctimer.h"

#include "schechter.h"
#include "cosmocalcs.h"
#include "gftdist.h"
#include "simdata.h"
#define PI 3.141592

void usage(void);
void usage(void)
{
	cout << endl<<" Usage: simulateLSSTobs [...options...]" << endl<<endl;

	cout << " -i : INFILE: input filename of catalog generated by baseSimulation "<<endl;
	cout << " -o : OUTFILE: output filename for observed LSST catalog (will be stored in output/)"<<endl;
	cout << " -s : SEDFILE: reading model galaxy SEDs from file SEDFILE [DEFAULT=CWWK.list]"<<endl;
	cout << " -n : NINTERP: number of interpolations between SEDs [DEFAULT=0]"<<endl;
	cout << " -y : NYEARS: number of years of LSST observations [DEFAULT=10]"<<endl;
	cout << " -t : NELLIPTICAL,NSPIRAL: number of elliptical, spiral SEDs [DEFAULT=1,2]"<<endl;
	cout << " -r : Don't add host galaxy reddening "<<endl;
	cout << " -m : Don't include Madau absorption "<<endl;
	cout << endl;
}
int main(int narg, char* arg[])
{
	cout << " ==== simulateLSSTobs.cc program , to simulate LSST data ==== "<<endl;

	// make sure SOPHYA modules are initialized 
    SophyaInit();  
	FitsIOServerInit();
	InitTim();
	cout<<endl<<endl;

	//--- decoding command line arguments 
	string outfile="output/simulateLSSTobs.fits";
	string infile="output/testbasesim.fits";
	string sedFile = "CWWK.list";
    int nYear = 10;
    int nInterp = 0;
    int nElliptical = 1;
    int nSpiral = 2;
    
    // Number of visits per year (Table 1, Ivezic et al 2008)
    int uVisitsPerYear = 6;
    int gVisitsPerYear = 8;
    int rVisitsPerYear = 18;
    int iVisitsPerYear = 18;
    int zVisitsPerYear = 16;
    int yVisitsPerYear = 16;
    
    bool isAddRedden = true;
    bool isAddMadau = true;
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"hrmi:o:y:n:t:")) != -1)  {
	    switch (c) {
            case 'i' :
		        infile = optarg;
		        break;
	        case 'o' :
		        outfile = optarg;
		        outfile = "output/" + outfile;
		        break;
	        case 's' :
		        sedFile = optarg;
		        break;
	        case 'y' :
		        sscanf(optarg,"%d",&nYear);
		        break;
		    case 'n' :
		        sscanf(optarg,"%d",&nInterp);
		        break;
		    case 't' :
		        sscanf(optarg,"%d,%d",&nElliptical,&nSpiral);
		        break;
		    case 'r' :
		        isAddRedden = false;
		        break;
		    case 'm' :
		        isAddMadau = false;
		        break;
	      case 'h' :
		    default :
		    usage(); return -1;
		    }
	    }
	
	// total number of visits
	int uVisits = uVisitsPerYear*nYear;
    int gVisits = gVisitsPerYear*nYear;
    int rVisits = rVisitsPerYear*nYear;
    int iVisits = iVisitsPerYear*nYear;
    int zVisits = zVisitsPerYear*nYear;
    int yVisits = yVisitsPerYear*nYear;
	
    cout << "     Input base catalog is "<< infile <<endl;
    cout << "     Reading SEDs from file "<< sedFile <<", interpolating "<< nInterp;
    cout << " times between them"<< endl;
    cout << "     Number of ellipticals = "<< nElliptical <<", number of spirals = ";
    cout << nSpiral << endl;
    if (isAddRedden)
        cout <<"     Adding host galaxy reddening "<<endl;
    else
        cout <<"     No host galaxy reddening added "<<endl;
    if (isAddMadau)
        cout <<"     Applying Madau absorption "<<endl;
    else
        cout <<"     No Madau absorption included "<<endl;    
    cout << "     Output catalog is "<< outfile <<endl;
    cout << "     Number of years of LSST observations = "<< nYear <<endl;
    cout << endl;
    //-- end command line arguments
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
	ResourceUsage res;
	InitTim();
	
	// this controls the drawing of random numbers
	RandomGenerator rg;

	// Set cosmology
	double h=0.7, OmegaM=0.3, OmegaL=0.7;
	SimpleUniverse su(h,OmegaM,OmegaL);
	su.SetFlatUniverse_OmegaMatter();
	cout <<"     Set cosmology to: OmegaM="<<OmegaM<<", OmegaL="<<OmegaL;
	cout <<", H0="<<100*h<<endl;
	cout << endl;

	// Read in catalog generated by baseSimulation
	// Contains basic galaxy properties: redshift, absolute magnitude, 
	// broad type (elliptical, spiral, starburst), 
	FitsInOutFile fin(infile,FitsInOutFile::Fits_RO);
	fin.MoveAbsToHDU(2);
	SwFitsDataTable dt(fin,512,false);
	sa_size_t ng=dt.NEntry();
	sa_size_t nc=dt.NCols();
	DataTableRow row=dt.EmptyRow();
	cout <<"     In file "<<infile<<" ... "<<endl;
	cout <<"     Number of columns = "<<nc<<", number of entries = "<<ng;
	cout << endl;

	// Output file
	cout <<"     Creating output file "<<outfile<<endl;
	FitsInOutFile swf(outfile,FitsInOutFile::Fits_Create);
	
	// binary data table
	SwFitsDataTable gals(swf,2048);
	gals.AddFloatColumn("zs");
	gals.AddFloatColumn("am");
	gals.AddFloatColumn("type");
	gals.AddFloatColumn("ext");
	gals.AddFloatColumn("mu");
	gals.AddFloatColumn("mg");
	gals.AddFloatColumn("mr");
	gals.AddFloatColumn("mi");
	gals.AddFloatColumn("mz");
	gals.AddFloatColumn("my");
	gals.AddFloatColumn("muo");
	gals.AddFloatColumn("mgo");
	gals.AddFloatColumn("mro");
	gals.AddFloatColumn("mio");
	gals.AddFloatColumn("mzo");
	gals.AddFloatColumn("myo");
	gals.AddFloatColumn("emu");
	gals.AddFloatColumn("emg");
	gals.AddFloatColumn("emr");
	gals.AddFloatColumn("emi");
	gals.AddFloatColumn("emz");
	gals.AddFloatColumn("emy");
	DataTableRow rowin=gals.EmptyRow();
	cout << endl;

	// Load in filters required
	// wavelength range of the SEDs/filters
	cout <<"     Load in LSST filters"<<endl;
	double lmin=5e-8, lmax=2.5e-6;
	string filterFile = "LSST.filters";
	ReadFilterList readFilterList(filterFile);
	readFilterList.readFilters(lmin,lmax);
	vector<Filter*> LSSTfilters=readFilterList.getFilterArray();
	int nFilter = readFilterList.getNTot();
	cout <<"     Read in "<< nFilter <<" filters"<<endl;
	//Filter restFrameFilter((*GOODSfilters[1]));
	int iU = 0;
	int iG = 1;
	int iR = 2;
	int iI = 3;
	int iZ = 4;
	int iY = 5;
	cout <<endl;
	
	cout <<"     Load in GOODS B filter"<<endl;
	string goodsFilterFile = "GOODSB.filters";
	ReadFilterList readGOODSBfilter(goodsFilterFile);
	readGOODSBfilter.readFilters(lmin,lmax);
	vector<Filter*> goodsBFilter=readGOODSBfilter.getFilterArray();
	
	// Load in SEDs
	ReadSedList readSedList(sedFile);
    readSedList.readSeds(lmin,lmax); // Read out SEDs into array
    int nsed=readSedList.getNSed(); // Get total number of SEDs
    cout <<"     Number of original SEDs = "<<nsed<<endl;
	cout << endl;
	
	
	if (nInterp>0) { // Interpolate SEDs
	    readSedList.interpSeds(nInterp);
        cout <<"     Interpolated SEDs "<<nInterp<<" times "<<endl;
        // Reorder SEDs
        readSedList.reorderSEDs();
        }    
    int ntot = readSedList.getNTot();
    vector<SED*> sedArray=readSedList.getSedArray();
    cout <<"     Final total number of SEDs = "<<ntot<<endl;
	cout << endl;

	// Prepare the class which will calculate the magnitudes
	cout <<"     Initialize class to calculate magnitudes"<<endl;
	SimData simgal(sedArray,LSSTfilters,su,rg,nElliptical,nSpiral);
    cout << endl;
    
    // Add Madau preference
    simgal.setMadau(isAddMadau);

    // Loop over all galaxies in the base catalog
	cout <<"     Start loop over galaxies ..."<<endl;
	Timer tm("TimingGalaxies",false);
	tm.Split();

	for (int i=0; i<ng; i++)
		{
		cout <<"     gal "<<i+1<<" of "<<ng<<endl;
		dt.GetRow(i,row);
		
		// Get galaxy properties: galaxy has redshift zs, absolute magnitude am,
		// SED type type and internal extinction ext.
		double zs=row[0];
		double am=row[1];
		int typ=(int)row[2];
		double type=simgal.SimSED(typ);
		double ext;
		if (isAddRedden)
		    ext = simgal.SimRed(type);
		else
		    ext = 0.;

        // Calculate galaxy magnitude in observed filters ugrizy
        // Galaxy's absolute magnitude is defined in filter (*goodsBFilter[0])
		double uMagTh=simgal.GetMag(zs,type,am,ext,iU,(*goodsBFilter[0]));
		double gMagTh=simgal.GetMag(zs,type,am,ext,iG,(*goodsBFilter[0]));
		double rMagTh=simgal.GetMag(zs,type,am,ext,iR,(*goodsBFilter[0]));
		double iMagTh=simgal.GetMag(zs,type,am,ext,iI,(*goodsBFilter[0]));
		double zMagTh=simgal.GetMag(zs,type,am,ext,iZ,(*goodsBFilter[0]));
		double yMagTh=simgal.GetMag(zs,type,am,ext,iY,(*goodsBFilter[0]));
		
		// The final observations
		// The 1st element is the value of the observed magnitude
		// The 2nd element is the magnitude error
		vector<double> uObservation = simgal.addLSSTuError(uMagTh,uVisits);
		vector<double> gObservation = simgal.addLSSTgError(gMagTh,gVisits);
		vector<double> rObservation = simgal.addLSSTrError(rMagTh,rVisits);
		vector<double> iObservation = simgal.addLSSTiError(iMagTh,iVisits);
		vector<double> zObservation = simgal.addLSSTzError(zMagTh,zVisits);
        vector<double> yObservation = simgal.addLSSTyError(yMagTh,yVisits);

        // Write the data to the FITS file
		rowin[0]=zs;
		rowin[1]=am;
		rowin[2]=type;
		rowin[3]=ext;
		rowin[4]=uMagTh;
		rowin[5]=gMagTh;
		rowin[6]=rMagTh;
		rowin[7]=iMagTh;
		rowin[8]=zMagTh;
		rowin[9]=yMagTh;
		rowin[10]=uObservation[0];
		rowin[11]=gObservation[0];
		rowin[12]=rObservation[0];
		rowin[13]=iObservation[0];
		rowin[14]=zObservation[0];
		rowin[15]=yObservation[0];
		rowin[16]=uObservation[1];
		rowin[17]=gObservation[1];
		rowin[18]=rObservation[1];
		rowin[19]=iObservation[1];
		rowin[20]=zObservation[1];
		rowin[21]=yObservation[1];
		gals.AddRow(rowin);
		}
	cout <<"     End loop"<<endl;
	tm.Split();
	cout <<"     .... done, took "<< tm.PartialElapsedTime()/60 <<" mins";
	
	// Write information on simulation to FITS header
	DVList  dvl;
    dvl("BaseSim") = infile;
	dvl("SedFile") = sedFile;
	dvl("NEllip") = nElliptical;
	dvl("NSpiral") = nSpiral;
	dvl("NInterp") = nInterp;
	dvl("NYearObs") = nYear;
	swf.WriteHeaderRecords(dvl);
	swf.MoveAbsToHDU(2);
	
	}  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " simulateLSSTobs.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " simulateLSSTobs.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " simulateLSSTobs.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of simulateLSSTobs.cc program  Rc= " << rc << endl;
  return rc;	
}
