/**
 * @file  gftdist.h
 * @brief Contains a series of classes for generating galaxy distributions;
 *        galaxy redshifts, types, magnitudes. Use CumulDistZ and DrawZ when
 *        simulating a distribution without clustering. Just use CumulDistM,
 *        DrawM, TypeRatio when simulating galaxies with clustering.
 *
 * @todo <CODE>CumulDistZ::Output2File</CODE> <BR>
 *       When the calculated redshift CDF is output to a file the cosmology and
 *       luminosity function information used to calculate it should be included
 *       
 *
 * @todo <CODE>CumulDistZ::SetUp</CODE> <BR>
 *       When the calculated redshift CDF is read from a file into an 
 *       interpolation function, that function should set any values outside of
 *       zmin, zmax to zero.  Making this fix <I>shouldn't</I> affect the 
 *       results anyway though
 *
 * @todo <CODE>DrawZ::Draw</CODE> <BR>
 *       Does it matter that redshifts are drawn from the CDF on the open interval?
 *       This means that the exact values of zmin and zmax will <I>never</I> be
 *       drawn.
 *
 * @todo <CODE>DrawM::SetArray</CODE> <BR>
 *       Find out exactly why some of the interpolated magnitude CDF values have  
 *       the value nan
 *
 * @todo <CODE>DrawM::Draw</CODE> <BR>
 *       Replace with 2D interpolation table
 *
 * @todo Move <CODE>VolElement</CODE> somewhere more appropriate
 *
 * @todo Write wrapper class like <CODE>SimBaseCatalog</CODE> when simulating
 *       galaxies with clustering (ie leave out the DrawZ part)
 *
 * @author Alex Abate, Reza Ansari, ...
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2009
 * @date 2009
 *
 */
#ifndef  GFTDIST_H_SEEN
#define  GFTDIST_H_SEEN


// generic stuff
#include <iostream>
#include <fstream>
#include <stdlib.h>

// sophya stuff
#include "machdefs.h"
#include "sopnamsp.h"
#include "array.h"
// Sophya update v2.3 June 2013 replaces genericfunc with classfunc
//#include "genericfunc.h"
#include "classfunc.h"
#include "pexceptions.h"
#include "randinterf.h" // This define RandomGenerator interface class
#include "ctimer.h"
#include "fiosinit.h"
#include "fitsioserver.h"
#include "swfitsdtable.h"

// DirectSim
#include "mydefrg.h"  // definition of default RandomGenerator
#include "sinterp.h"
#include "geneutils.h"
#include "cosmocalcs.h"
#include "schechter.h"

namespace SOPHYA {


//******* CumulDistZ *********************************************************//
/** @class
  * CumulDistZ class
  *
  * Class to create cumulative distribution to draw redshifts from 
  *
  *
  */
class CumulDistZ {
public:
	/** Default constructor */
	CumulDistZ()
	: lfpars_(lfp_default_), su_(su_default_), zmin_(0), zmax_(0),
	  nptinteg_(0), nptinterp_(0) {    };
	
	/** Constructor when calculating table
	    @param lfpars       class holding LF parameters and their evolution
	    @param su           class that does cosmological calculations
	    @param zmin         minimum redshift of cumulative distribution function
	    @param zmax         maximum redshift of cumulative distribution function
	    @param nptinteg     number of points to use in integration of ?
	    @param nptinterp    number of points to use in interpolation of ?     
	    @param prt          print level                                       */
	CumulDistZ(LFParameters& lfpars, SimpleUniverse& su, double zmin=0.,
			double zmax=6., int nptinteg=1000, int nptinterp=100, int prt=0) 
	: lfpars_(lfpars), su_(su), zmin_(zmin), zmax_(zmax), nptinteg_(nptinteg),
      nptinterp_(nptinterp) {
			SetUp(lfpars_, su_, zmin_, zmax_, nptinteg_, nptinterp_, prt); };
				
	/** Constructor when reading table 
	    @param fname    FITS file to read from
        @param prt      printing level                                        */
	CumulDistZ(string fname, int prt=0)
	: lfpars_(lfp_default_), su_(su_default_), zmin_(0), zmax_(0), nptinteg_(0),
      nptinterp_(0) { SetUp(fname, prt); };

	//METHODS//
	
	/** Return value of cumultive redshift distribution function at redshift z
	    \f$ F_z(z)=\frac{\int_0^z\int_M_1^M_2 \phi(M,z')dV(z')dM}
	                {\int_0^z_{max}\int_M_1^M_2  \phi(M,z')dV(z')dM}          */
	virtual double operator() (double z) const { return schZint_(z); };

    /** Set up for reading distribution from a FITS bintable file 
        @param fname    FITS file to read from
        @param prt      printing level                                        */
	void SetUp(string fname, int prt=0);
	
	/** Set up for doing the actual calculation                
	    @param lfpars       class holding LF parameters and their evolution
	    @param su           class that does cosmological calculations
	    @param zmin         minimum redshift of cumulative distribution function
	    @param zmax         maximum redshift of cumulative distribution function
	    @param nptinteg     number of points to use in integration of \f$\phi(z)\f$
	    @param nptinterp    number of points to use in interpolation of \f$F_z(z)\f$*/
	void SetUp(LFParameters& lfpars, SimpleUniverse& su, double zmin=0.,
			double zmax=6., int nptinteg=1000, int nptinterp=100, int prt=0);
				
	/** Print n values of the distribution between zmin and zmax+0.1 as a check 
	    @param nvals    number of values to print                             */
	void PrintDist(int nvals) {
		cout << "     CumulDistZ::PrintDist Printing "<< nvals <<" values from";
		cout << " the cumulative z dist ... "<<endl;
		//int ntot=zv_.size();
		//int nskip=ntot/nvals;
		double dz = (zmax_+0.1 - zmin_)/(nvals - 1);
		for (int i=0; i<nvals; i++) { 
		    double z = zmin_ + i*dz;
		    double val = schZint_(z);
		    cout << "     "<< z <<"  "<< val <<endl; 
		    }
		    //int j=i*nskip; 
			//cout << "     "<<zv_[j] <<"  "<<scv_[j]<<endl; }
		};

    /** Output 2D array of cumulative redshift distribution to a FITS file
         binary table.  The zmin and zmax of the calculated distribution are
         included in the filename
        @param outfileroot  root name of FITS file to write to                */
	void Output2File(string outfileroot);

	/** Return min and max redshift of the redshift cdf 
	    @param zmin     min redshift
	    @param zmax     max redshift                                          */
	void returnZminZmax(double& zmin, double& zmax)
	    { zmin = zmin_; zmax = zmax_; };

protected:
	LFParameters& lfpars_;      /**< class that stores the LF parameters and evolution */
	SimpleUniverse& su_;        /**< class that calculates cosmological quantities     */
	double zmin_;               /**< min redshift of cdf                      */
	double zmax_;               /**< max redshift of cdf                      */
	int nptinteg_;              /**< number of points in integration of \f$\phi(z)\f$*/
	int nptinterp_;             /**< number of points in interpolation of \f$F_z(z)\f$*/
	vector<double> zv_;         /**< grid of z values for \f$F_z(z)\f$        */
	vector<double> scv_;        /**< grid of F values (\f$F_z(z)\f$)          */
	SInterp1D schZint_;         /**< interpolation of \f$F_z(z)\f$            */
	SimpleUniverse su_default_;
	LFParameters lfp_default_;

};


//******* CumulDistM *********************************************************//
/** @class
  * CumulDistM class
  *
  * Class to create cumulative distribution to draw magnitudes from
  *
  * Returns the value of the cumulative magnitude distribution function given a
  * magnitude M and redshift z, \f$ F_M(M,z) \f$
  *
  */
class CumulDistM {
public:

	/** Default constructor */
  CumulDistM(): lfpars_(lfpars_default_) , cosmoc_(su_default_,0.,10.,1000) , Mmin_(0) , Mmax_(0)
	{    };

	/** Constructor 
	@param lfpars   holds LF parameters and their evolution
	@param su       does cosmological calclations
	@param Mmin     sets lower integration bound
	@param Mmax     sets upper integration bound                              */
	CumulDistM(LFParameters& lfpars, SimpleUniverse& su, double Mmin=-24, double Mmax=-13)
		: lfpars_(lfpars) , cosmoc_(su,0.,10.,1000) , Mmin_(Mmin) , Mmax_(Mmax){    };
		
	// copy constructor
	//CumulDistM(CumulDistM const& a);

    /** Returns value of cumultive magnitude distribution function at magnitude M
        given redshift z
	    \f$ F_M(M,z)=\frac{\int_M_1^M_2 \phi(M',z)dV(z)dM'}
	                {\int_M_1^M_2  \phi(M',z)dV(z)dM'}                        */
	virtual double operator()(double m, double z) {
	  /*   Change par Reza 
		SchechterVol schvol(lfpars_,su_);
		schvol.SetInteg(Mmin_,m);
		double top=schvol.Integrate(z);
		schvol.SetInteg(Mmin_,Mmax_);
		double bot=schvol.Integrate(z);
       	        return top/bot; 
	  */
	  //---- introduit par Reza 
	  double ps,ms,a;
	  lfpars_(z,ps,ms,a);
	  Schechter sch(ps,ms,a);
	  double top=0.;
	  int npt=10000*((m-Mmin_)/(Mmax_-Mmin_));
	  if (npt>1) {
	    sch.SetInteg(Mmin_,m,npt);
	    top=sch.Integrate();
	  }
	  double bot=top;
	  npt=10000*((Mmax_-m)/(Mmax_-Mmin_));
	  if (npt>1) {
	    sch.SetInteg(m,Mmax_,npt);
	    bot+=sch.Integrate();
	  }
	  double retv=0.; 
	  if (bot>1.e-39)  retv=top/bot;
	  return retv;
	  //---- introduit par Reza 
	}
		
	/** Return magnitude integration limits                                   */
    void returnMminMmax(double& Mmin, double& Mmax)
	    { Mmin = Mmin_; Mmax = Mmax_; };

protected:
	LFParameters& lfpars_;      /**< class that stores the LF parameters and evolution */
  // Reza: We use InterpCosmoCalc instead	SimpleUniverse& su_;        /**< class that calculates cosmological quantities     */
        InterpCosmoCalc cosmoc_;        /**< cosmology                                  */
	double Mmin_;               /**< min absolute magnitude of cdf            */
	double Mmax_;               /**< max absolute magnitude of cdf            */
	LFParameters lfpars_default_;
	SimpleUniverse su_default_;
};


//******* DrawZ **************************************************************//
/** @class
  * DrawZ class
  *
  * Class to draw redshifts from cumulative distribution 
  *
  */
class DrawZ {
public:

    /** Constructor
        @param cumz cumulative distribution to draw from
        @param npt  number of points to use in (reverse) interpolation table  */
	DrawZ(CumulDistZ& cumz, RandomGeneratorInterface& rg, int npt=10000)
	: rg_(rg) {
	
	    // Get zmin and zmax
	    cumz.returnZminZmax(zmin_, zmax_);
	    
	    // Reverse interpolation table
	    vector<double> zvals, cumzvals;
		double dz = (zmax_ - zmin_)/(npt - 1);
		for (int i=0; i<npt; i++) { 
			double z=zmin_ + i*dz;
			double cv=cumz(z);	
			zvals.push_back(z);
			cumzvals.push_back(cv);
			//if ( i<10 || i>npt-10)
			//    cout << z <<"  "<< cv <<endl;
			}

        // create interpolation
        // It won't matter if cumzvals is <0 or >1 because redshifts at these 
        // CDF values will never be drawn because a flat random number generator 
        // between 0 and 1 (open interval but shouldn't matter) is used
	    revCumZ_.DefinePoints(cumzvals,zvals,0.,1.,zvals.size()*4);

		};

    /** Draw redshift from redshift cdf                                       */
	double Draw(){
	    // draw a flat random number between 0 and 1 (open interval)
	    double rn=rg_.Flat01();
	    double zs = revCumZ_(rn); // cdf goes from 0 to 1, return z at that value
	    return zs;
	    };

	/** Return min and max z of cdf 
	    @param zmin     min z of cdf
	    @param zmax     max z of cdf                                          */
	void returnZminZmax(double& zmin, double& zmax)
	    { zmin = zmin_; zmax = zmax_; };

protected:
	SInterp1D revCumZ_;             /**< reverse interpolation of \f$F_z(z)\f$ */
	RandomGeneratorInterface& rg_;  /**< random number generation              */
	double zmin_;                   /**< min z of cdf                          */
	double zmax_;                   /**< max z of cdf                          */
};


//******* DrawM **************************************************************//
/** @class
  * DrawM class
  *
  * Class to draw magnitudes from cumulative distribution 
  *
  */
class DrawM {
public:

	/** Default constructor 
	    @param cumm cumulative magnitude distribution \f$ F_M(M,z) \f$ (interp table) */
	DrawM(CumulDistM& cumm)
	: cumm_(cumm) , rg_(rg_default_) , mmin_(1) , mmax_(1), zmin_(-1) , zmax_(-1)
    {  };
	
	/** Constructor when calculating 2D array of cumulative magnitude distribution 
	    @param cumm     cumulative magnitude distribution \f$ F_M(M,z) \f$ (interp table)
	    @param rg       draws random numbers
	    @param mmin     minimum absolute mag value in 2D array
	    @param mmax     maximum absolute mag value in 2D array
	    @param zmin     minimum z value in 2D array
	    @param zmax     maximum z value in 2D array
	    @param nptz     number of z values in 2D array
	    @param nptm     number of absolute mag values in 2D array*/
	DrawM(CumulDistM& cumm, RandomGeneratorInterface& rg, double mmin=-24,
		double mmax=-13, double zmin=0., double zmax=6., int nptz=600, int nptm=600)
	: cumm_(cumm), rg_(rg), mmin_(mmin), mmax_(mmax), zmin_(zmin), zmax_(zmax)
    { SetUp(rg_, mmin_, mmax_, zmin_, zmax_, nptz, nptm); };
				
	/** Constructor when reading 2D array of cumulative magnitude distribution 
	    from a file 
	    @param fname    FITS filename containing pre-calculated 2D array of 
	                    F_M(M,z)                                              */
	DrawM(string fname)
	: cumm_(cumm_default_), rg_(rg_default_), mmin_(1), mmax_(1), zmin_(-1), zmax_(-1)
    { SetUp(fname); };
    
	// Destructor
	//virtual ~DrawM();
	
	/** Calculate 2D array of cumulative magnitude distribution              */
	void SetUp(RandomGeneratorInterface& rg, double mmin=-24, double mmax=-13,
	           double zmin=0., double zmax=6., int nptz=600, int nptm=600);

    /** Read 2D array of cumulative magnitude distribution from a file       */
	void SetUp(string infile);
	
    /** Output 2D array of cumulative magnitude distribution to a FITS file
        The magnitude and redshift range values are written in the filename
        @param outfileroot  root name of FITS file to write to                */
	void Output2File(string outfileroot);

    /** Store cumulative magnitude distribution values in an array           */
  void SetArray();   
    /** Given a redshift z draw a magnitude according the the cumulative distribution
        in the array @cumval_
        THIS SHOULD BE REPLACED WITH A 2D INTERP TABLE                        */
	double Draw(double z);

	/** Return min and max redshift of 2D array used to draw magnitudes from  */
	void returnZminZmax(double& zmin, double& zmax)
	    { zmin = zmin_; zmax = zmax_; };

protected:
	CumulDistM& cumm_;  /**< cumulative magnitude distribution \f$ F_M(M,z) \f$ (interp table) */
	RandomGeneratorInterface& rg_;/**< draws random numbers                       */
	TVector<r_8> mv_;           /**< magnitude grid values                        */
	TVector<r_8> zv_;           /**< redshift grid values                         */
	TArray<r_8> cumval_;        /**< 2D array of cumulative mag distribution      */
	double mmin_;               /**< minimum absolute magnitude value in 2D array */
	double mmax_;               /**< maximum absolute magnitude value in 2D array */
	double dm_;                 /**< step in absolute magnitude in 2D array       */
	double zmin_;               /**< minimum z value in 2D array                  */
	double zmax_;               /**< maximum z value in 2D array                  */
	double dz_;                 /**< step in z in 2D array                        */
	CumulDistM cumm_default_;
	DR48RandGen rg_default_;
};


//******* DrawType ***********************************************************//
/** @class
  * DrawType class
  *
  * Class to draw galaxy type
  *
  */
class DrawType {
public:
    /** Constructor 
        @param tr   ratio of each galaxy type with redshift and magnitude
        @param rg   random number generator                                   */
	DrawType(TypeRatio tr, RandomGeneratorInterface& rg)
	: tr_(tr) , rg_(rg) {  };

	/** Return 
	    @param */
	virtual int operator() (double m, double z) { 
	
	    double rn = rg_.Flat01();
		double Fe, Fs, Fl;
		// return type fractions at the redshift and magnitude
		tr_(m, z, Fe, Fs, Fl);
	    if (rn>=0 && rn<Fe)
		    return 1;
		else if (rn>=Fe && rn<(Fe+Fs))
			return 3;
		else if (rn>=(Fe+Fs))
			return 2;
		else return 0;  // CHECK , is this correct - returning 1  - Reza, April 2016 
		 };
	
protected:
    TypeRatio tr_;                  /**< ratio of each galaxy type with mag, z  */
	RandomGeneratorInterface& rg_;  /**< random number generator                */
};


//******* SimBaseCatalog *****************************************************//
/** @class
  * SimBaseCatalog class
  *
  * Simulates the base catalog of galaxy properties: redshift, magnitude and 
  * type. No clustering
  *
  */
class SimBaseCatalog {
public:

    /** Constructor: simulate galaxy properties without clustering 
        @param drz  redshift distribution
        @param drm  absolute magnitude distribution
        @param drt  galaxy type distribution                                  */
	SimBaseCatalog(DrawZ& drz, DrawM& drm, DrawType& drt)
	: drz_(drz) , drm_(drm) , drt_(drt)
	{  
	
	    double zminRS, zmaxRS, zminMAG, zmaxMAG;
	    drz_.returnZminZmax(zminRS, zmaxRS);
        drm_.returnZminZmax(zminMAG, zmaxMAG);
        
        if ( zminRS<zminMAG || zmaxRS>zmaxMAG )
            throw ParmError("ERROR! Magnitude and redshift CDF's don't have matching redshift ranges");
        else {
            zmin_ = zminRS;
            zmax_ = zmaxRS;
	        }
	};

    /** Simulate ngal galaxies and output to FITS file
        @param ngal         number of galaxies to simulation
        @param outfileroot  FITS file to write galaxies to                    */
	void DoSim(long ngal, string outfileroot);
	
	/** Simulate single galaxy  
	    @param z    simulated galaxy redshift
	    @param am   simulated galaxy absolute magnitude
	    @param typ  simulated galaxy type                                     */
	void DoSim(double& z, double& am, double& typ)
	    {   z=drz_.Draw();
		    am=drm_.Draw(z);
		    typ=(double)drt_(am,z); };

protected:
	DrawZ& drz_;        /**< redshift distribution                            */
	DrawM& drm_;        /**< absolute magnitude distribution                  */
	DrawType& drt_;     /**< galaxy type distribution                         */
	double zmin_;       /**< min redshift of distributions                    */
	double zmax_;       /**< max redshift of distributions                    */

};


//******* NGalZ **************************************************************//
/** @class
  * NGalZ class
  *
  * Class to calculate total number of galaxies between two redshifts and over
  * some solid angle
  *
  */
class NGalZ {
public:

    /** Constructor 
        @param lfpars   LF parameters as a function of z
        @param su       class that does cosmological calculations
        @param npt      number of points to integrate LF with                 */
	NGalZ(LFParameters& lfpars, SimpleUniverse& su, int npt=1000)
       : lfpars_(lfpars) , su_(su) , npt_(npt) {       };
 
    /** Returns the number of galaxies between two redshifts and some solid angle
        @param zmin minimum redshift
        @param zmax maximum redshift
        @param sa   solid angle in steradians                                 */
	virtual long operator() (double zmin, double zmax, double sa) {
	
	    // The below class returns phi(z) = [int phi(M|z) dM]*dV(z)
		SchechterZVol schZ(lfpars_, su_);
		double val=schZ.Integrate();
		long ntot=(long)val*sa;
		return ntot; 
		
		};

    /** Reset the number of points to integrate the LF with */
	void SetInteg(int npt)
		{ npt_ = npt; };

protected:
	LFParameters& lfpars_;  /**< class that holds LF parameters and their evolution */
	SimpleUniverse& su_;    /**< class that does cosmological calculations    */
	int npt_;               /**< number of points to integrate LF with        */
};


//******* VolElement *********************************************************//
/** @class
  * VolElement class
  *
  * Calculate dV/dz as a function of redshift
  *
  * @warning Move this class somewhere more appropriate?
  */
class VolElement : public ClassFunc1D
{
public:
    /** Constructor 
        @param su   does cosmological calculations
        @param sa   solid angle                                               */
    VolElement(SimpleUniverse& su, double sa)
	: su_(su) , sa_(sa) { }

	/** Return dV/dz
	    @param z    redshift                                                  */
    virtual double operator()(double z) const {
	    su_.SetEmissionRedShift(z);
	    double dVdz=su_.HubbleLengthMpc()*pow(su_.RadialCoordinateMpc(),2)*su_.Gz(z)*sa_;
        return (dVdz);
        }
  
protected:
    SimpleUniverse& su_;    /**< does cosmological calculations               */
    double sa_;             /**< solid angle                                  */
};


//******* GalFlxTypDist ******************************************************//
/** @class
  * GalFlxTypDist class
  *
  * This class can be used to define a galaxy distribution dist(magnitude, type) 
  * and then generate (magnitude,type) pair values according to the defined 
  * distribution 
  *
  * @warning This class has now been superceded by the DrawM, DrawType classes
  */
class GalFlxTypDist {
public:
// Constructor, needs a RandomGenerator object
explicit GalFlxTypDist(RandomGeneratorInterface& rg, int prtlev=0);
explicit GalFlxTypDist(int prtlev=0);
virtual ~GalFlxTypDist();

// To be called once magnitude distribution for all 
// types have been defined through AddGalType() call 
// return true if OK, generate an exception (error) if problem
bool Check();

// Adds a galaxy type with the corresponding magnitude distribution
// See the .cc file for more information 
int AddGalType(ClassFunc1D& magf, double magmin, double magmax, double fr=1, 
               int nbinmag=20, int nstep=2000);
// Once all galaxy type-magnitude functions are added using AddGalType, pick a 
// galaxy type and mag at random according to these distributions using:
int GetGalaxy(int& typ, r_8& mag) const ;
// Added by AA Feb 2011:
int GetGalMag(r_8& mag) const ;
void AddSchechter(Schechter& t1,Schechter& t2,Schechter& t3);
void GetGalType(double mag,int& type);

// Return the number of different types defined
size_t NTypes() const { return v_mag_.size(); }
// Return the number of generated galaxies
size_t NGalGen() const { return totngal_; }

// Print some information about the galaxy type list
ostream& Print(ostream& os) const;

// THESE SHOULD BE PUT IN SEPARATE CLASS - they simulate the redshift distribution 
// according to the volume element 
int GetZDist(ClassFunc1D& dVdz, double zmin,double zmax, int nbin=100, int nstep=2000);
int GetZ(double&) const ;

void PrintDist(string);

// .... class variables 
RandomGenerator rgdefault_;
RandomGeneratorInterface& rg_;
int prtlev_ ;						// print/debug level 
vector< TVector<r_8> > v_mag_;		// Magnitude bins for each type
vector< r_8 > v_frac_;				// Fraction of total number of galaxies, for each type
vector< r_8 > v_sumfrac_;			// Sum of fractions for each type Sum [ 0 ... type] 
r_8 tot_frac_ ;						// Sum of type fraction for all types, should be = 1 when all types are defined
mutable size_t totngal_;			// total number of galaxies generated 
TVector<r_8> dVdzbins_;
//Schechter& t1_,t2_,t3_;
bool SchTypes_;
Schechter dt1,dt2,dt3;  // Default Schecter functions
Schechter t1_,t2_,t3_;  // Schechter function objects
};

/*! operator << overloading - Prints the list on the stream \b s */
inline ostream& operator << (ostream& s, GalFlxTypDist const & gfd)
  {  gfd.Print(s);  return(s);  }
  
  

} // End namespace SOPHYA
#endif
