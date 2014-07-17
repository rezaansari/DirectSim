#include "sedfilter.h"

namespace SOPHYA {

/*******************************************************************************
*                                                                              *
*                                SED CLASSES                                   *
*                                                                              *
*******************************************************************************/


//******* SED ****************************************************************//

// Constructor if reading flux values from a file for interpolation
SED::SED(string& fname, double xmin, double xmax, int npt)
: sed2_(sed2init_)
{
  	string Tplace;
	char * pt=getenv("SEDLOC");
	if (pt==NULL)
		throw ParmError("ERROR SED LOCATION ENVIRONMENT VARIABLE -SEDLOC- NOT DEFINED");
	else
		{
		Tplace=pt;
		cout <<"    Location of template files are "<<Tplace<<endl;
		}
	string filename=Tplace+fname;
	
	// Set interpolation to zero outside of x-range in file 
    bool isSetZero = false; // this needs to be checked
    int nComments = 0;
	sed_.ReadXYFromFile(filename, xmin, xmax, npt, nComments, isSetZero);

    isRead_=true;
	isInterp_ = false;
	isRedden_ = false;
	_test_=2.;
};


// copy constructor
SED::SED(SED const& a)
:  sed_(a.sed_)
{
	Set(a);
	_test_=3.;

};


// Copy SED: this does a DEEP copy
SED& SED::Set(const SED& a)
{
    sed2_=a.sed2_;// missing this line cost me a whole live long day!
	isRead_ = a.isRead_;
    isRedden_=a.isRedden_;
    isInterp_=a.isInterp_;
    a_=a.a_;
    b_=a.b_;
    EBmV_=a.EBmV_;
    RvCard_=a.RvCard_;
    law_=a.law_;
	
};


//******* SED methods ********************************************************//

void SED::readSED(string& filename, double xmin, double xmax, int npt)
{

    // Set interpolation to zero outside of x-range in file 
    bool isSetZero = false; // this needs to be checked
    int nComments = 0;
	sed_.ReadXYFromFile(filename, xmin, xmax, npt, nComments, isSetZero);
	isRead_=true;

};

// Setting up interpolation
// If interpolating between 2 SEDs
// Must take as argument a POINTER to an SED object
void SED::doInterp(SED* sed2,double a,double b)
{
    // copy sed2 into sed2_
    sed2_=sed2;
    
    a_=a; b_=b;

	isInterp_ = true;
	isRedden_ = false;
	
};

// Setting up the reddening
void SED::doRedden(double EBmV, int law, double RvCard)
{
    EBmV_=EBmV;
    law_=law;
    RvCard_=RvCard;

	isRedden_ = true;

};

// Main return function
double SED::returnFlux(double lambda) const
{

    if (isInterp_&&!isRedden_)// if SED is to be interpolated
        return interpSED(lambda);
    else if(isRedden_&&!isInterp_)// if SED is to be reddened
        return addReddening(lambda);
    else if (isRedden_&&isInterp_)// if SED is to be interpolated AND reddened
        return interpAddReddening(lambda);
    else
        return sed_(lambda);

};

double SED::addReddening(double lambda) const
{
    if (!isRedden_)
        throw ParmError("ERROR! Cannot redden spectrum");
        
    Reddening red;
	double k;
	if (law_<1)
		k=red.Cardelli(lambda,RvCard_);
	if (law_>0)
		k=red.Calzetti(lambda);
	
	return  (sed_(lambda)*pow(10,-0.4*k*EBmV_));
};


double SED::interpAddReddening(double lambda) const
{

    if (!isRedden_)
        throw ParmError("ERROR! Cannot redden spectrum");
    if (!isInterp_)
        throw ParmError("ERROR! Cannot interpolate spectrum");
        
    Reddening red;
	double k;
	if (law_<1)
		k=red.Cardelli(lambda,RvCard_);
	if (law_>0)
		k=red.Calzetti(lambda);
	
	double redPart = pow(10,-0.4*k*EBmV_);
	double interpPart = interpSED(lambda);
	return  (interpPart*redPart);

};


double SED::interpSED(double lambda) const
{

    if (!isInterp_)
        throw ParmError("ERROR! Cannot interpolate spectrum");
        
    // sed2_ is a POINTER to a SED object
    double sed1=a_*(sed_(lambda));
    double sed2=b_*(sed2_->returnFlux(lambda));
    return (sed1+sed2);
};

//******* ReadSedList ********************************************************//

ReadSedList::ReadSedList(string sedFile,int prt)
: prt_(prt)
{

    // Initialize these to zero to start
    ntot_=nsed_=0;

	// first get location of SED files
	sedDir_=getSedDirEnviromentVar();
	sedFileFullPath_=sedDir_+sedFile;
		
	// open file, read all lines and count number of SEDs within
	countSeds(sedFileFullPath_);
    cout <<"     There are "<<nsed_<<" SEDs to read in"<<endl;
    cout <<endl;
    
    ntot_=nsed_; // sets both total number of SEDs and number of sed's read in
    
    isInterp_ = isRedden_ = false;

};

//******* ReadSedList methods ************************************************//

string ReadSedList::getSedDirEnviromentVar()
{

    string tPlace;
    char * pt=getenv("SEDLOC");
	if (pt==NULL) {
		string emsg="ERROR SED LOCATION ENVIRONMENT VARIABLE";
		emsg+=" -SEDLOC- NOT DEFINED";
		throw ParmError(emsg);
		}
	else {
		tPlace=pt;
		cout <<"     Location of SED template files is "<<tPlace<<endl;
		}
	cout << endl;
	return tPlace;

};

void ReadSedList::countSeds(string sedFileFullPath)
{
    ifstream ifs;
    ifs.open(sedFileFullPath.c_str(),ifstream::in);
	if (ifs.fail())
		{
		string emsg="error: Unable to open file ";
		emsg+=sedFileFullPath;
		throw ParmError(emsg);
		}
    string line;
    nsed_=0;
   	while ( ifs.good() )
    	{
   		getline(ifs,line);
	    nsed_++;
   		}
   	nsed_-=1;
   	ifs.close();
};

void ReadSedList::readSeds(double lmin,double lmax)
{
    ifstream ifs;
    
    // read in all the SED filenames
	ifs.open(sedFileFullPath_.c_str(),ifstream::in);
	string fileNames[nsed_];
	string line;
	if (prt_ > 0)
	    cout <<"     File contains the following SEDs:"<<endl;
	for (int i=0; i<nsed_; i++)
		{
		getline(ifs,line);
		sedFiles_.push_back(line);
		fileNames[i]=sedDir_+line;
		if (prt_ > 0)
		    cout <<"     "<<fileNames[i]<<endl;
		}
		
	// read in all sed files
	if (prt_ > 0)
	    cout <<"     Reading in all the SED files ... "<<endl;

	for (int i=0; i<nsed_; i++)
		{
		sedArray_.push_back(new SED());// assigned memory for SED pointer
		sedArray_[i]->readSED(fileNames[i],lmin,lmax); 
		}

};

void ReadSedList::interpSeds(int nInterp)
{

    if (isRedden_)
        throw ParmError("ERROR! Make sure interpolation is done first!");

    // If we have nsed's we have nsed-1 spaces inbetween seds
    // Therefore number of SEDs being added is nInterp*(nsed-1)

    if (prt_>0) {
	    cout <<"     Before interpolation, size of SED array is ";
	    cout << sedArray_.size()<<endl;
	    }
	    
	int iSED=ntot_; // starting index of next SED
    for (int i=0; i<(nsed_-1); i++)
		{
		double a=0;

		// loop over number of interps to do
		for (int j=0; j<nInterp; j++)
		    {
		    // need to find values of a and b depending on what nInterp is
		    a+=1./(nInterp+1);
		    double b=1-a;
		
		    // push new SED to end of array
		    sedArray_.push_back(new SED(*(sedArray_[i])));// copies SED in i
		    sedArray_[iSED]->doInterp(sedArray_[i+1],a,b);
		                        
            iSED++;
		    }
		 
		}
		
	int nAdd=nInterp*(nsed_-1);
    ntot_+=nAdd;
    
    if (prt_>0){
	    cout <<"     Size of SED array is now "<<sedArray_.size()<<endl;
	    cout <<"     Total number of SEDs is "<<ntot_<<endl;
	    cout << endl; 
        }
        
     isInterp_ = true; 
};

void ReadSedList::reddenSeds(int nStepRed,double redMax)
{
// This is not properly implemented, the value of redMax and law will change
// depending on the SED involved
// Will need to define another array that encodes these properties e.g. an
// integer array ..
// broadtype_[i] = 0,1,2
// where 0 == elliptical type: redMax cannot be >0.1, law=0 (Cardelli)
//       1 == late type: no limit on redMax, law=0 (Cardelli)
//       1 == starburst type: no limit on redMax, law=1 (Calzetti)
    

    if (prt_>0) {
	    cout <<"     Before reddening, size of SED array is ";
	    cout << sedArray_.size()<<endl;
	    }
	    
	double redStep=redMax/nStepRed;
	int law=0;
	for (int i=0; i<ntot_; i++)
		{
		cout <<"     On SED "<<i+1<<" of "<<ntot_<<endl;
		// loop over number of reddening steps to do
		for (int j=0; j<nStepRed; j++)
		    {
		    
		    double EBmV=(j+1)*redStep;
		    cout <<"     On redden "<<j+1<<" of "<<nStepRed;
		    cout <<", EB-V="<<EBmV<<endl;
		    // push new SED to end of array
		    sedArray_.push_back(new SED(*(sedArray_[i])));
		    int lastIndex=sedArray_.size()-1;
		    sedArray_[lastIndex]->doRedden(EBmV,law);// Not properly implemented!
		    }
		}
		
    ntot_=(ntot_)*nStepRed+ntot_; // reddening all seds nStepRead times

    if (prt_>0){
	    cout <<"     Size of SED array is now "<<sedArray_.size()<<endl;
	    cout <<"     Total number of SEDs is "<<ntot_<<endl; 
        }
        
    isRedden_ = true;

};

void ReadSedList::reorderSEDs()
{
    if (isRedden_)
        throw ParmError("ERROR! Cannot reorder SEDs if reddening has been done!");
    if (!isInterp_)
        throw ParmError("ERROR! No need to reorder if not interpolated!");

    vector<SED*> tmpsedArray;
    
    // If we have nsed's we have nsed-1 spaces inbetween seds
    // Therefore number of SEDs being added is: nInterpd = nBetween*(nsed-1)
    // nBetween is the number of interpolations done between each original SED
    
    
    int nInterpd = ntot_ - nsed_;
    int nBetween = nInterpd/(nsed_ - 1);
    
    //int iBeginInt = nsed_; // index of sedArray_ where the interpolated SEDS begin
    
    for (int i=0; i<(nsed_-1); i++) {
        int ii = i*nBetween;
        tmpsedArray.push_back(sedArray_[i]); 
        for (int j=0; j<nBetween; j++) {
            int ik = nsed_ + ii + j;
            tmpsedArray.push_back(sedArray_[ik]); 
            }
        }
    tmpsedArray.push_back(sedArray_[nsed_-1]); 
    
    int tmp = tmpsedArray.size();
    if ( tmp!=ntot_ )
        throw ParmError("ERROR! Wrong number of SEDs!");
            
    sedArray_.clear();
    sedArray_=tmpsedArray;
};

void ReadSedList::writeSpectra(string outFile,double lmin,double lmax,int nl)
{
    ifstream ifs;
	ofstream outp;
	
	cout <<"     Spectra number = "<<sedArray_.size()<<endl;
	
    double dl=(lmax-lmin)/(nl-1);
	ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int j=0; j<nl; j++)
			{
			double lam=lmin+j*dl;
			//cout <<"     Writing wavelength "<<j+1<<endl;
			outp <<lam<<"  ";
			for (int i=0; i<ntot_; i++)
				{
				//cout <<"     Writing spectrum "<<i+1<<endl;
				double val=sedArray_[i]->returnFlux(lam);
				outp << val <<"  ";
				}
			outp <<endl;
			if ( prt_ & (j<100) ) {
				cout <<lam<<"  ";
				for (int k=0; k<ntot_; k++)
					cout<<sedArray_[k]->returnFlux(lam)<<" ";
				cout <<endl;
				}
			}
		outp.close();
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;
	cout << endl;

};


/*******************************************************************************
*                                                                              *
*                              FILTER CLASSES                                  *
*                                                                              *
*******************************************************************************/


//******* Filter *************************************************************//

Filter::Filter(string& fname, double xmin, double xmax, int nComments, 
                                                      bool zero_outside,int npt)
{
	string Fplace;
	char * pf=getenv("FILTLOC");
	if (pf==NULL) {
	    string emsg = "ERROR FILTER LOCATION ENVIRONMENT VARIABLE -FILTLOC- ";
	    emsg += "NOT DEFINED";
		throw ParmError(emsg);
        }
	else {
		Fplace=pf;
		cout <<"    Location of filter files are "<<Fplace<<endl;
		}

	string filename = Fplace+fname;
	// zero_outside = true: then interpolation sets stuff to zero outside xmin,xmax
	// definitely want this to be true for the filter
    ReadXYFromFile(filename,xmin,xmax,npt,nComments,zero_outside);
};

//******* Filter methods *****************************************************//

void Filter::readFilter(string& fname, double lmin, double lmax, int nComments,
                                            bool zero_outside, int npt)
{
    // zero_outside = true: then interpolation sets stuff to zero outside xmin,xmax
	// definitely want this to be true for the filter
    ReadXYFromFile(fname,lmin,lmax,npt,nComments,zero_outside);

};

//******* ReadFilterList *****************************************************//

ReadFilterList::ReadFilterList(string filterFile,int prt)
{
    // Initialize these to zero to start
    ntot_=0;

	// first get location of filter files
	filterDir_=getFilterDirEnviromentVar();
	filterFileFullPath_=filterDir_+filterFile;
		
	// open file, read all lines and count number of filters within
	countFilters(filterFileFullPath_);
    cout <<"     There are "<<ntot_<<" filters to read in"<<endl;
    cout <<endl;
    
};
    
//******* ReadFilterList methods *********************************************//

string ReadFilterList::getFilterDirEnviromentVar()
{

    string tPlace;
    char * pt=getenv("FILTLOC");
	if (pt==NULL) {
		string emsg="ERROR FILTER LOCATION ENVIRONMENT VARIABLE";
		emsg+=" -FILTLOC- NOT DEFINED";
		throw ParmError(emsg);
		}
	else {
		tPlace=pt;
		cout <<"     Location of Filter template files is "<<tPlace<<endl;
		}
	cout << endl;
	return tPlace;

};
    

void ReadFilterList::countFilters(string filterFileFullPath)
{
    ifstream ifs;
    ifs.open(filterFileFullPath.c_str(),ifstream::in);
	if (ifs.fail())
		{
		string emsg="error: Unable to open file ";
		emsg+=filterFileFullPath;
		throw ParmError(emsg);
		}
    string line;
    ntot_=0;
   	while ( ifs.good() )
    	{
   		getline(ifs,line);
	    ntot_++;
   		}
   	ntot_-=1;
   	ifs.close();

};
    

void ReadFilterList::readFilters(double lmin, double lmax)
{
    ifstream ifs;
    
    // read in all the filter filenames
	ifs.open(filterFileFullPath_.c_str(),ifstream::in);
	string fileNames[ntot_];
	string line;
	if (prt_ > 0)
	    cout <<"     File contains the following filters:"<<endl;
	for (int i=0; i<ntot_; i++) {
		getline(ifs,line);
		fileNames[i]=filterDir_+line;
		if (prt_ > 0)
		    cout <<"     "<<fileNames[i]<<endl;
		}
		
	// read in all sed files
	if (prt_ > 0)
	    cout <<"     Reading in all the filter files ... "<<endl;

	for (int i=0; i<ntot_; i++) {
		filterArray_.push_back(new Filter());// assigned memory for SED pointer
		filterArray_[i]->readFilter(fileNames[i],lmin,lmax); 
		}

};
    

void ReadFilterList::writeFilters(string outFile,double lmin,double lmax,int nl)
{
    ifstream ifs;
	ofstream outp;
	
	cout <<"     Filter number = "<<filterArray_.size()<<endl;
	
    double dl=(lmax-lmin)/(nl-1);
	ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int j=0; j<nl; j++)
			{
			double lam=lmin+j*dl;
			//cout <<"     Writing wavelength "<<j+1<<endl;
			outp <<lam<<"  ";
			for (int i=0; i<ntot_; i++)
				{
				//cout <<"     Writing spectrum "<<i+1<<endl;
				double val=filterArray_[i]->operator()(lam);
				outp << val <<"  ";
				}
			outp <<endl;
			if ( prt_ & (j<100) ) {
				cout <<lam<<"  ";
				for (int k=0; k<ntot_; k++)
					cout<<filterArray_[k]->operator()(lam)<<" ";
				cout <<endl;
				}
			}
		outp.close();
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;
	cout << endl;


};





}// end namespace SOPHYA
