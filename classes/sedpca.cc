#include "sedpca.h"

namespace SOPHYA {


/*******************************************************************************
*                                                                              *
*                              PCA RELATED CLASSES                             *
*                                                                              *
*******************************************************************************/

//******* SEDCovMat methods **************************************************//
TArray<double> SEDCovMat::MakeCovMat(SED **sedarray, int nsed,
						TVector<double>& meanSED)
{

	// Initialize size of SED matrix
	TArray<double> SEDmatrix;
	int ndim=2;
	sa_size_t mydim[ndim];
	mydim[0]=nl_; mydim[1]=nsed;
	SEDmatrix.SetSize(ndim, mydim);
	
	// Vector of mean SEDs
	meanSED.SetSize(nsed);
	meanSED=0; // set vector elements to equal zero
	
	// loop over wavelengths and SEDs and fill matrix, and calculate mean
	// SED values
	for (int i=0; i<nl_; i++)
		{
		double lam=lmin_+i*dl_;
		for (int j=0; j<nsed; j++)
			{
			double val=sedarray[j]->returnFlux(lam);
			//cout <<val<<endl;
			SEDmatrix(i,j)=val;
			meanSED(j)+=SEDmatrix(i,j);
			}
		}
	meanSED/=nl_;
	for (int j=0; j<nsed; j++)
		cout <<"mean of SED "<<j<<" is "<<meanSED(j)<<endl;
	
	// Subtract mean
	for (int i=0; i<nl_; i++)
		for (int j=0; j<nsed; j++)
			SEDmatrix(i,j)-=meanSED(j);
			
	/*cout <<"print SEDmatrix part to screen M(0:10,1:6) "<<endl;
	for (int i=0; i<100; i++)
		{
		cout <<" lam = "<<lmin_+i*dl_<<"  ";
		for (int j=0; j<6; j++)
			cout <<SEDmatrix(i,j)<<"  ";
		cout <<endl;
		}*/

	TArray<double> CovMat=TransposeMult(SEDmatrix);
	CovMat/=(nl_-1);
	
	return CovMat;
};


TArray<double> SEDCovMat::TransposeMult(TArray<double> SEDmatrix)
{

	//SVD matinv(SEDmatrix);
	// Do Transpose
	TArray<double> TSEDmatrix=Transpose(SEDmatrix);
	cout <<"     Size of original matrix = "<<SEDmatrix.SizeX()<<"x";
	cout <<SEDmatrix.SizeY()<<endl;
	cout <<"     Size of transposed matrix = "<<TSEDmatrix.SizeX()<<"x";
	cout <<TSEDmatrix.SizeY()<<endl;
	/*for (int i=0; i<6; i++)
		{
		for (int j=100; j<110; j++)
			cout <<TSEDmatrix(i,j)<<"  ";
		cout <<endl;
		}*/
	
	// Do multiplication
	TArray<double> TransMult=Mult(SEDmatrix,TSEDmatrix);

	return TransMult;

};

//******* TemplatePCA ********************************************************//

TemplatePCA::TemplatePCA(vector<SED*> sedArray,double lmin,double lmax,int nl)
: sedArray_(sedArray) , lmin_(lmin) , nl_(nl)
{

    dl_=(lmax-lmin_)/(nl_-1);
    nsed_=sedArray_.size();
    
    dataMatrix_.SetSize(nsed_,nl_);
    cout <<"     Set data matrix size: "<<dataMatrix_.NRows()<<"x";
    cout << dataMatrix_.NCols() << endl;
    cout << endl;

    // Normalize each spectrum such that fmean=f/(sqrt(sum(f^2)))
    // This is not the only option for normalization
    cout <<"     Find normalizations of spectra "<<endl;
    normalizeSpectra();
    cout << endl;
    
    // Find mean value of each spectrum
    cout <<"     Find means of spectra "<<endl;
    meanSpectra();
    cout << endl;
    
    // Add each normalized, mean subtracted spectrum to the data matrix
    cout <<"     Add spectra to data matrix "<<endl;
    addSpectra();
    cout << endl;
    
    // Calculate covariance matrix
    cout <<"     Calculate covariance matrix of data"<<endl;
    calculateCovarianceMatrix();
    cout << endl;

    // Do PCA
    cout <<"     Find eigenvalues and eigenvectors "<<endl;
    doPCA();
    cout << endl;
   
};

TemplatePCA::TemplatePCA(vector<SED*> sedArray, string eigVectFile, double lmin,
                                                        double lmax, int nl)
: sedArray_(sedArray) , lmin_(lmin) , nl_(nl)
{

    dl_=(lmax-lmin_)/(nl_-1);
    nsed_=sedArray_.size();
    
    dataMatrix_.SetSize(nsed_,nl_);
    cout <<"     Set data matrix size "<<dataMatrix_.NRows()<<"x";
    cout << dataMatrix_.NCols() << endl;
    cout << endl;

    // Normalize each spectrum such that fmean=f/(sqrt(sum(f^2)))
    // This is not the only option for normalization
    cout <<"     Find normalization of spectra "<<endl;
    normalizeSpectra();
    cout << endl;
    
    // Find mean value of each spectrum
    cout <<"     Find mean of spectra "<<endl;
    meanSpectra();
    cout << endl;
    
    // Add each normalized, mean subtracted spectrum to the data matrix
    cout <<"     Add spectra to data matrix "<<endl;
    addSpectra();
    cout << endl;
    
    // Don't need to calculate covariance matrix or do PCA
    
    // Read eigenvalues from a file
    cout <<"     Read eigenvectors from file "<<endl;
    eigenVectors_.ResizeTo(nl_,nl_);
    eigenVectors_=readEigenVectors(eigVectFile);
    
    cout <<"    Check eigenvectors ... "<<endl;
    cout <<"    Size of eigenvector matrix = "<< eigenVectors_.GetNrows() <<"x";
    cout << eigenVectors_.GetNcols() <<endl;
    cout <<"    Printing first 10x10 ... "<<endl;
    for (int i=0; i<10; i++) {
        for (int j=0; j<10; j++)
            cout <<eigenVectors_(i,j) <<" ";
        cout << endl;
        }
            
    
    cout << endl;
    
};

//******* TemplatePCA methods ************************************************//

void TemplatePCA::normalizeSpectra()
{
    normValues_.SetSize(nsed_);
    
    for (int i=0; i<nsed_; i++)
        {
        double sumSquared=0;
		for (int j=0; j<nl_; j++)
			{
			double lam=lmin_+j*dl_;
			double val=sedArray_[i]->returnFlux(lam);
			sumSquared+=(val*val);
			}
		normValues_(i)=sqrt(sumSquared);
		}
		
};


void TemplatePCA::meanSpectra()
{
    meanValues_.SetSize(nsed_);

    for (int i=0; i<nsed_; i++)
        {
        double sum=0;
		for (int j=0; j<nl_; j++)
			{
			double lam=lmin_+j*dl_;
			double val=sedArray_[i]->returnFlux(lam);
			sum+=(val/normValues_(i));
			}
		meanValues_(i)=sum/nl_;
		}

};


void TemplatePCA::addSpectra()
{

    for (int i=0; i<nsed_; i++) {
        //cout << "     Adding spectrum "<<i+1<<" of "<<nsed_<<endl;
        for (int j=0; j<nl_; j++) {
            
            //cout <<"     On lambda "<<j+1<<" of "<<nl_<<endl;
			double lam=lmin_+j*dl_;
			double val=sedArray_[i]->returnFlux(lam);
			dataMatrix_(i,j)=(val/normValues_(i))-meanValues_(i);
			
			}
        }

};


void TemplatePCA::calculateCovarianceMatrix()
{

    // Transpose data matrix
    //cout <<"     Transpose data matrix "<<endl;
    TMatrix<double> dataMatrixTransposed;
	dataMatrixTransposed=Transpose(dataMatrix_);
	
	// Make covariance matrix: D^T*D 
	//cout <<"     Make covariance matrix "<<endl;
	covMatrix_=Mult(dataMatrixTransposed,dataMatrix_);
	
	// Normalize covariance matrix by its trace (sum(D_ii))
	//cout <<"     Normalize covariance matrix "<<endl;
    normalizeCovarianceMatrix();
};
	        

void TemplatePCA::normalizeCovarianceMatrix()
{
    double traceCovMatrix=0;
	for (int i=0; i<nl_; i++)
	    traceCovMatrix += covMatrix_(i,i);
	cout <<"     -Trace of the covariance matrix is = "<< traceCovMatrix <<endl;
	cout <<"     -Normalize by this value"<<endl;
	covMatrix_/=traceCovMatrix;

};


void TemplatePCA::doPCA()
{
    // Copy covariance matrix to ROOT object 
    //cout <<"     Copy covariance matrix to ROOT object"<<endl;
    TMatrixD fCovarianceMatrix=copyCovMatrixToROOT();

    // Decompose into eigenvectors and eigenvalues
    //TVectorD fEigenValues(nl_);
    //TMatrixD fEigenVectors(nl_,nl_);
    eigenValues_.ResizeTo(nl_);
    eigenVectors_.ResizeTo(nl_,nl_);
    TMatrixDSym sym; sym.Use(fCovarianceMatrix.GetNrows(),
	                            fCovarianceMatrix.GetMatrixArray());
    TMatrixDSymEigen eigen(sym);
    eigenVectors_ = eigen.GetEigenVectors();
    eigenValues_ = eigen.GetEigenValues();
    
};


TMatrixD TemplatePCA::copyCovMatrixToROOT()
{
       
	TMatrixD fCovarianceMatrix(nl_,nl_);

	// convert to ROOT object
	for (int i=0; i<nl_; i++)
        for (int j=0; j<nl_; j++)
	        fCovarianceMatrix(i,j)=covMatrix_(i,j);
	        
    return fCovarianceMatrix;
    
};


void TemplatePCA::reconstructSpectra(int nEigKept)
{

    // first zero size incase have called method before with a different
    // nEigKept
    reconstructedSpectra_.ZeroSize();
	eigenvalsProjSpec_.ZeroSize();

    reconstructedSpectra_.SetSize(nl_,nsed_);
	eigenvalsProjSpec_.SetSize(nEigKept,nsed_);
    for (int i=0; i<nsed_; i++)
        {
	    int iSpectrum=i;
	    cout <<"     On spectrum "<<iSpectrum+1<<" of "<<nsed_<<endl;
	    TMatrix<double> reconstructedSpectrum;
	    TVector<double> zCut=reconstructSpectrum(nEigKept,iSpectrum,
	                                                     reconstructedSpectrum);
	    cout << "     Size of projected eigenvalues = "<<zCut.Size()<<endl;
	                                                     
	    for (int j=0; j<nl_; j++)
	        reconstructedSpectra_(j,i)=reconstructedSpectrum(j,1);
	        
	    for (int j=0; j<nEigKept; j++) // think this should be nl!!!!
	        eigenvalsProjSpec_(j,i)=zCut(j);
	        
	    cout << endl;
	    }
	cout << endl;

};

TVector<double> TemplatePCA::reconstructSpectrum(int nEigKept, int iSpectrum, 
                                        TMatrix<double>& reconstructedSpectrum)
{

    // Keep only 1st nEigKept
    //cout <<"     Only taking first "<<nEigKept<<" eigenvalues "<<endl;
    
    // Chop the eigenvector matrix at the maximum eigenvalue to keep
	TMatrixT<double> eigenVectorsCut(eigenVectors_.GetNrows(),nEigKept);
	for (int i=0; i<eigenVectors_.GetNrows(); i++)
		for (int j=0; j<nEigKept; j++)
            eigenVectorsCut(i,j)=eigenVectors_(i,j);
            
    /*cout <<"    Printing first nEigKeptxnEigKept of chopped eigenvector matrix "<<endl;
    for (int i=0; i<nEigKept; i++) {
        for (int j=0; j<nEigKept; j++)
            cout <<eigenVectorsCut(i,j) <<" ";
        cout << endl;
        }*/
            
    // Transpose the chopped eigenvector matrix 
	//cout <<"     Transpose cut eigenvector matrix"<<endl;
	cout <<"     Size of eigenvector matrix:";
	int nr=eigenVectorsCut.GetNrows();
	int nc=eigenVectorsCut.GetNcols();
	cout <<nr<<"x"<<nc<<endl;
	TMatrixT<double> eigenVectorsCutTransposed(nc,nr);
	eigenVectorsCutTransposed.Transpose(eigenVectorsCut);
	cout <<"     Size of transposed eigenvector matrix:";
	cout <<eigenVectorsCutTransposed.GetNrows()<<"x";
	cout <<eigenVectorsCutTransposed.GetNcols()<<endl;
	
	/*cout <<"    Printing first nEigKeptxnEigKept of chopped eigenvector matrix";
	cout <<" transposed "<<endl;
    for (int i=0; i<nEigKept; i++) {
        for (int j=0; j<nEigKept; j++)
            cout <<eigenVectorsCutTransposed(i,j) <<" ";
        cout << endl;
        }*/
	
	// check that the size of eigenvector matrix is ok
	if (nr!=nl_) {
	    stringstream ss1,ss2;
	    ss1<<nr; ss2<<nl_;
	    string emsg="ERROR! Length of eigenvectors "+ss1.str()+" is not equal to the";
	    emsg+=" number of wavelengths "+ss2.str();
	    throw ParmError(emsg);
	    }
	
    // Get mean subtracted data 
	TMatrixT<double> spectrumMeanSubtracted(nr,1);// needs to be a TMatrix type
	for (int i=0; i<nl_; i++){
			//double lam=lmin_+i*dl_;
			//double val=sedArray_[iSpectrum]->returnFlux(lam);
			spectrumMeanSubtracted(i,0)=dataMatrix_(iSpectrum,i);
			                            //(val/normValues_(iSpectrum))-
			                            //                      meanValues_(i);
			}

    // Matrix multiply the transposed eigenvectors and the mean subtracted spectrum
	// zCut are the eigenvalues of the spectrum projected onto the nEigKept
	// eigenvectors
	TMatrixT<double> zCut(nEigKept,1);
	zCut.Mult(eigenVectorsCutTransposed,spectrumMeanSubtracted);
	cout <<"     Result size = "<<zCut.GetNrows()<<"x"<<zCut.GetNcols()<<endl;
    
    
	// Rotate data back
	// Multiply each column of the truncated eigenvector matrix by each
	// eigenvalue of the projected spectrum
	TMatrixT<double> rot(nl_,nEigKept);
	for (int i=0; i<nl_; i++)
	    for (int j=0; j<nEigKept; j++)
	        rot(i,j)=eigenVectorsCut(i,j)*zCut(j,0);
	        
    // Sum over eigenvalues (the columns)
	TMatrixT<double> recSpec(nl_,1);
	recSpec=0;
	for (int i=0; i<nl_; i++)
	    for (int j=0; j<nEigKept; j++)
	        recSpec(i,0)+=( rot(i,j) ); //+meanValues_(i)/nEigKept );
	     
	     
	// Add back mean spectrum
	for (int i=0; i<nl_; i++)
	       recSpec(i,0) = recSpec(i,0) + meanValues_(iSpectrum);
	        
	
	// Basically done now!
	
	// Convert to Sophya objects
	TVector<double> eigenvaluesSpectrum(nEigKept);
	for (int i=0; i<nEigKept; i++)
	    eigenvaluesSpectrum(i)=zCut(i,0);
	
	reconstructedSpectrum.SetSize(nl_,2);
	for (int i=0; i<nl_; i++)
	        {
	        double lam=lmin_+i*dl_;
	        reconstructedSpectrum(i,0)=lam; // add in lambda for good measure
	        reconstructedSpectrum(i,1)=recSpec(i,0);
	        }
	        
	return eigenvaluesSpectrum;

};


TVector<double> TemplatePCA::fitSpectra()
{

    TVector<double> fitValues(nsed_);
    for (int i=0; i<nsed_; i++)
        {
        double chisq=0.;
        for (int j=0; j<nl_; j++)
			{
			double lam=lmin_+j*dl_;
			double sp1=sedArray_[i]->returnFlux(lam);
			double sp2=reconstructedSpectra_(j,i)*normValues_(i);
			
            chisq+=((sp1-sp2)*(sp1-sp2));

			}
		fitValues(i)=chisq/(double)nl_;
		}
	return fitValues;

};


TMatrix<double> TemplatePCA::getCovMatrix()
{


    int nr=covMatrix_.NRows();
    int nc=covMatrix_.NCols();
    TMatrix<double> covMatrixSophyaMatrix(nr,nc);
    for (int i=0; i<nr; i++)
			for (int j=0; j<nc; j++)
				covMatrixSophyaMatrix(i,j)=covMatrix_(i,j);
				
	return covMatrixSophyaMatrix;

};

TVector<double> TemplatePCA::getEigenValues()
{

    int nv=eigenValues_.GetNoElements();
    TVector<double> eigenValuesSophyaVector(nv);
    for (int i=0; i<nv; i++)
        eigenValuesSophyaVector(i)=eigenValues_(i);
        
    return eigenValuesSophyaVector;
};

TMatrix<double> TemplatePCA::getEigenVectors()
{

    int nr=eigenVectors_.GetNrows();
    int nc=eigenVectors_.GetNcols();

    TMatrix<double> eigenVectorsSophyaMatrix(nr,nc);
    for (int i=0; i<nr; i++)
			for (int j=0; j<nc; j++)
				eigenVectorsSophyaMatrix(i,j)=eigenVectors_(i,j);
	return eigenVectorsSophyaMatrix;

};

// WRITE TO FILES

// Data normalization values
void TemplatePCA::writeNormValues(string outFile)
{

    ifstream ifs;
	ofstream outp;
    ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int i=0; i<normValues_.Size(); i++)
				outp<< normValues_(i) <<endl;
		outp.close();	
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;   


};


// Data mean values
void TemplatePCA::writeMeanValues(string outFile)
{

    int nMeanValues = meanValues_.Size();
    //cout <<"     Number of mean values = "<<nMeanValues<<endl;

    ifstream ifs;
	ofstream outp;
    ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int i=0; i<nMeanValues; i++)
				outp<< meanValues_(i) <<endl;
		outp.close();	
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;   


};


// Covariance matrix
void TemplatePCA::writeCovMatrix(string outFile)
{

    int nr = covMatrix_.NRows();
    int nc = covMatrix_.NCols();
    //cout <<"     Size of covariance matrix = "<< nr << "x" << nc <<endl;
    
    ifstream ifs;
	ofstream outp;
    ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int i=0; i<nr; i++)
			{
			for (int j=0; j<nc; j++)
				outp<< covMatrix_(i,j)<<"  ";
			outp << endl; 
			}
		outp.close();	
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;

};

// data matrix
void TemplatePCA::writeDataMatrix(string outFile)
{

    
    ifstream ifs;
	ofstream outp;
    ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int i=0; i<nsed_; i++) {
			for (int j=0; j<nl_; j++)
                outp<< dataMatrix_(i,j)<<"  ";
			outp << endl;
			}
		outp.close();	
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;

};

// eigenvectors
void TemplatePCA::writeEigenVectors(string outFile)
{

    ifstream ifs;
	ofstream outp;
    ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int i=0; i<eigenVectors_.GetNrows(); i++)
			{
			for (int j=0; j<eigenVectors_.GetNcols(); j++)
				outp<< eigenVectors_(i,j)<<"  ";
			outp << endl;
			}
		outp.close();	
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;

};

// eigenvalues
void TemplatePCA::writeEigenValues(string outFile)
{

    ifstream ifs;
	ofstream outp;
    ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int i=0; i<eigenValues_.GetNrows(); i++)
				outp<< eigenValues_(i)<<endl;
		outp.close();	
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;



};

// both eigenvectors and eigenvalues
void TemplatePCA::writeEigenValVecs(string outFile)
{

    //const TVectorD* eigenValues=principal_->GetEigenValues();
    
    ifstream ifs;
	ofstream outp;
    ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int i=0; i<eigenVectors_.GetNrows(); i++)
			{
			for (int j=0; j<eigenVectors_.GetNcols(); j++)
				outp<< eigenVectors_(i,j)<<"  ";
			outp <<eigenValues_(i)<< endl;
			}
		outp.close();	
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;

};


void TemplatePCA::writeEigenValsOfProjSpec(string outFile)
{

    ifstream ifs;
	ofstream outp;
    ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int i=0; i<eigenvalsProjSpec_.NRows(); i++)
			{
			for (int j=0; j<nsed_; j++)
			    outp <<eigenvalsProjSpec_(i,j) << "  ";
			outp << endl;
			}
		outp.close();	
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;
	cout << endl;


};

void TemplatePCA::writeRecSpec(string outFile)
{
    ifstream ifs;
	ofstream outp;
    ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int i=0; i<nl_; i++)
			{
			double lam=lmin_+i*dl_;
			outp <<lam<<"  ";
			for (int j=0; j<nsed_; j++)
			    outp << reconstructedSpectra_(i,j)*normValues_(j) << "  ";
			outp << endl;
			}
		outp.close();	
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;
	cout << endl;

};

// Read from file
TMatrixD TemplatePCA::readEigenVectors(string inFile)
{

        ifstream ifs;
        cout <<"     Reading in eigenvectors from file "<<inFile<<endl;
	    ifs.open(inFile.c_str(), ifstream::in);
	    sa_size_t nr, nc;
	    TArray<r_4> spectrum;
        spectrum.ReadASCII(ifs,nr,nc);
        cout <<"    Number of eigenvectors = "<<nc<<endl;
        cout <<"    Length of eigenvectors = "<<nr<<endl;

        //TMatrixT<double> eigenvals(nr,nc);
        TMatrixD eigenvals(nr,nc);
        for (int j=0; j<nr; j++)
            for (int i=0; i<nc; i++)
                eigenvals(j,i)=spectrum(i,j);
                
        return eigenvals;


};

}// end namespace SOPHYA
