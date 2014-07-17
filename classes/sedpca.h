/**
 * @file  sedpca.h
 * @brief Contains classes that perform PCA (Principal Component Analysis) calculations with SED 
 *
 *
 * @author Alex Abate 
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2013 
 * @date 2013
 *
 */
 
#ifndef SEDPCA_H_SEEN 
#define SEDPCA_H_SEEN 

#include "sedfilter.h"
// check the below is needed
#include "matrix.h"

// ROOT libraries
// removed these root libraries that are not needed
#include "TMinuit.h"
#include "TGraph.h"
#include "TFile.h"
#include "TDecompSVD.h"
#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"
#include "TPrincipal.h"

namespace SOPHYA {


/*******************************************************************************
*                                                                              *
*                              PCA RELATED CLASSES                             *
*                                                                              *
*******************************************************************************/



/** SEDCovMat class
  * 
  * Class to make SED covariance matrix from input point to SED class array
  * Probably defunct now using ROOT TPrincipal for calculating PCA
  *
  */
class SEDCovMat
{
public:
    SEDCovMat(double lmin,double dl,int nl)
    : lmin_(lmin) , dl_(dl) , nl_(nl) {  };

    TArray<double> MakeCovMat(SED **sedarray,int nsed,
                        TVector<double>& meanSED);
    TArray<double> TransposeMult(TArray<double> SEDmatrix);
    
    void SetWavelengths(double lmin,double dl,int nl)
        { lmin_=lmin; dl_=dl; nl_=nl; };
    
protected:
    double lmin_,dl_;
    int nl_;

};



/** TemplatePCA class
  * 
  * Class to calculate principal components of SED templates 
  *
  */
class TemplatePCA {
public:

    /** Constructor for calculating eigenvectors and eigenvalues from the spectra
        in sedArray */
    TemplatePCA(vector<SED*> sedArray, double lmin=5e-8, double lmax=2.5e-6, 
                                                                 int nl=1500);
    
    /** Constructor for projecting the spectra in sedArray onto the eigenvectors 
        read from file eigVectFile */
    TemplatePCA(vector<SED*> sedArray,string eigVectFile, double lmin=5e-8, 
                                                double lmax=2.5e-6,int nl=1500);
                                                                
    /** should add argument to here to reflect different normalization choices*/
    void normalizeSpectra();
    
    /** Find mean of spectra */
    void meanSpectra();
    
    /** add each spectrum to the data matrix */
    void addSpectra();
    
    /** Calculate covariance matrix of the data */
    void calculateCovarianceMatrix();
    
    /** Normalize covariance matrix by its trace */
    void normalizeCovarianceMatrix();
    
    /* Decompose into eigenvalues and eigenvectors */
    void doPCA();
    
    /* Copy the covariance matrix into a ROOT object */
    TMatrixD copyCovMatrixToROOT();
    
    /** Project each spectrum onto eigenvectors 1:nEigKept and find the 
        eigenvalues of the projected spectrum.
        
        Basically this is a loop over reconstructSpectrum method.
        Fills class matrices: reconstructedSpectra_ and eigenvalsProjSpec_.
        
        @param nEigKept number of eigenvalues to project onto */
    void reconstructSpectra(int nEigKept);
    
    /** Project spectrum iSpectrum onto eigenvectors 1:nEigKept 
        @param nEigKept number of eigenvalues to project onto
        @param iSpectrum spectrum to project
        @param reconstructedSpectrum reconstructed spectrum */
    TVector<double> reconstructSpectrum(int nEigKept,int iSpectrum,
                                TMatrix<double>& reconstructedSpectrum);   
                                
    /** Find out how good the fit of the reconstructed to the true spectrum is */
    TVector<double> fitSpectra();                                            
                              
    // These just return stuff as a SOPHYA object
    /** Return meanValues as SOPHYA object */                
    //TVector<double> getMeanValues();
    /** Return covariance matrix as SOPHYA object */
    TMatrix<double> getCovMatrix();
    /** Return eigenvalues as SOPHYA object */
    TVector<double> getEigenValues();
    /** Return eigenvectors as SOPHYA object */
    TMatrix<double> getEigenVectors();
    
    // These just return things
    /** Return reconstructed spectra */
    TMatrix<double> returnRecSpectra(){ return reconstructedSpectra_; };
    /** Return eigenvalues of reconstructed spectra */
    TMatrix<double> returnEigValsProjSpec(){ return eigenvalsProjSpec_; };
    /** Return spectra normalization values */
    TVector<double> returnNormValues(){ return normValues_; };
    
    // Methods to write stuff to a file
    /** Write spectra normalization values to a file */
    void writeNormValues(string outFile);
    /** Write spectra mean values to a file */
    void writeMeanValues(string outFile);
    /** Write covariance matrix to a file */
    void writeCovMatrix(string outFile);
    /** Write data matrix to a file */
    void writeDataMatrix(string outFile);
    /** Write eigenvectors to a file */
    void writeEigenVectors(string outFile);
    /** Write eigenvalues to a file */
    void writeEigenValues(string outFile);
    /** Write eigenvectors and eigenvalues to a file */
    void writeEigenValVecs(string outFile);
    /** Write eigenvalues of projected spectrum to a file */
    void writeEigenValsOfProjSpec(string outFile);
    /** Write reconstructed spectra to a file */
    void writeRecSpec(string outFile);
  
    /** Read eigenvectors from a file */
    TMatrixD readEigenVectors(string inFile);
    
    // Potentially add method(s) to calculate Karnhunen-Loeve angles

protected:
    vector<SED*> sedArray_;             /**< array of SED objects */
    double lmin_;                       /**< min wavelength */
    double dl_;                         /**< step in wavelength */
    int nl_;                            /**< number of wavelengths to in SED data matrix*/
    int nsed_;                          /**< number of SEDs */
    TMatrix<double> dataMatrix_;        /**< matrix of normalized, mean subbed SEDs */
    //TPrincipal* principal_;             /**< class that does PCA calc */
    TVector<double> normValues_;        /**< normalization of each spectrum */
    TVector<double> meanValues_;        /**< mean values of each wl bin (over SEDs)*/
    TMatrix<double> covMatrix_;         /**< data covariance matrix */
    TMatrixD eigenVectors_;             /**< eigenvectors */
    TVectorD eigenValues_;              /**< eigenvalues */
    TMatrix<double> reconstructedSpectra_;/**< rows: wavelength, cols: SEDs */
    TMatrix<double> eigenvalsProjSpec_; /**< doesn't depend on # of eigenvalues kept 
                                             (well length does). Size is (nEigKept,nsed)
                                             Is basically the projection of each SED
                                             onto the first nEigKept eigenvectors */
                                             
};                 

      
}// end namespace sophya

#endif
