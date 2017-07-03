#include <cstdio>
#include <cmath>
#include <iostream>
#include <cassert>
#include <fstream>
#include <tr1/random>
#include "phantom.h"
#include "memory.h"
#include <ctime>

#define Swap(a,b) temp=a; a=b; b=temp;

using namespace std;



phantom::phantom():buffer(NULL),totalScatters(0), bscArray(NULL), bscFreqArray(NULL)
{}


phantom::~phantom()
{
	if (buffer) delete[] buffer;
}
myVector phantom::getPhanSize() {
	return phanSize;
}

/*!This function creates a phantom object whose scatterers are uniformly distributed throughout the volume
 */
void phantom::createUniformPhantom(const myVector& phansize, double density, double soundSpeed, double atten0, double atten1, double atten2, char* fname )
{
	phantom::c0 =soundSpeed;
	phantom::a0 = atten0;
	phantom::a1 = atten1;
	phantom::a2 = atten2;
	totalScatters =  (int)(phansize.x*phansize.y*phansize.z*density * 1e9);
	cout << "The total number of scatterers is: " << totalScatters <<endl;
	phanSize = phansize;
	buffer = new scatterer[totalScatters];
	assert(buffer != NULL);
	
	/* initialize random seed: */
	tr1::mt19937 eng(time(NULL));  // uniform random integer generatot
	tr1::uniform_real<double> dist;
	tr1::variate_generator<tr1::mt19937, tr1::uniform_real<double> > generator(eng, dist);

	for (int i=0; i<totalScatters; i++)
	{
		buffer[i].x = generator()*phansize.x;
		buffer[i].y = generator()*phansize.y;
		buffer[i].z = generator()*phansize.z;


		//x coordinate is lateral, y is elevational, z is axial
	}

	readBscFromFile(fname);

}



/*!  Read in backscatter coefficients from a file.  The file format is binary
 * and all numbers are doubles.  The file contains:
 * frequency spacing (MHz)
 * number of points
 * list of backscatter coefficients 
 * |
 * |
 * |
 * ^
 * */
void phantom::readBscFromFile(char* fname)
{

	
	ifstream bscFile;
	bscFile.open(fname, std::ios::binary);

	if( !bscFile.is_open() ) {	
	cout << "Error opening backscatter file named: " << fname << std::endl;
	exit( EXIT_FAILURE );
	}

	bscFile.read( (char*) &(phantom::freqStep), sizeof(double) );
	cout << "The frequency spacing is " << phantom::freqStep << std::endl;
	double tmp;
	bscFile.read( (char*) &tmp    , sizeof(double) );	
	cout << "The number of backscatter coefficients read will be: " << tmp << std::endl;
	numBsc = static_cast<int>(tmp);

	bscArray = new double[numBsc];
	bscFile.read( (char*) bscArray , sizeof(double)*numBsc );	

	bscFreqArray = new double[numBsc];
	double tempFreq;
	for(int ind = 0; ind < numBsc; ind++)
	{
		tempFreq = phantom::freqStep*ind;
		bscFreqArray[ind] = tempFreq;
	}		
	

}



/*!  Given a frequency in MHz, return the backscatter coefficient value
 */
double phantom::giveBsc(double freq)
{

	int lowInd = static_cast<int>( freq/phantom::freqStep );
    double remainder = freq/phantom::freqStep;
    remainder -= static_cast<int>(remainder);
    double bscDiff = bscArray[lowInd] - bscArray[lowInd + 1];
    double bscCoeff = bscArray[lowInd] + remainder*bscDiff;
	if(bscCoeff < 0 ) bscCoeff = 0;
	return bscCoeff;

} 

/*!  This function saves the phantom data to a binary file.
 */
int phantom::savePhantom(char *filename)
{

	sortScatterer();
	ofstream fpout(filename, ios::binary);


	if(!fpout.is_open() )
    {
       cout << "unable to create phantom file " << filename <<endl;
       return -1;
    }

	fpout.write((char*)&phanSize, sizeof(myVector) );
	fpout.write((char*)&totalScatters, sizeof(int) );
	fpout.write( (char*)&c0, sizeof(double));
	fpout.write( (char*)&a0, sizeof(double));
	fpout.write( (char*)&a1, sizeof(double));
	fpout.write( (char*)&a2, sizeof(double));
	fpout.write((char*)buffer, totalScatters*sizeof(scatterer));
	fpout.write( (char*) &numBsc, sizeof(int) );
	fpout.write( (char*) &freqStep, sizeof(double) );
	fpout.write( (char*) bscArray, numBsc*sizeof( double) );
	fpout.write( (char*) bscFreqArray, numBsc*sizeof( double) );
	fpout.close();
	return(1);
}

/*!  This function reads in a phantom from a file created by savePhantom
 */
int phantom::loadPhantom(char* filename)
{
	ifstream fpin;
	fpin.open(filename, ios::binary);
    
	if(!fpin.is_open() ){
        cout << "Error reading phantom file " << filename << std::endl;
        return 0;
    }

	fpin.read((char*)&phanSize, sizeof(myVector) );
	fpin.read((char*)&totalScatters, sizeof(int) );
	fpin.read((char*)&c0, sizeof(double) );
	fpin.read((char*)&a0, sizeof(double) );
	fpin.read((char*)&a1, sizeof(double) );
	fpin.read((char*)&a2, sizeof(double) );
	
        std::cout<< "The phantom size is: " << phanSize.x*1E3 <<" mm laterally \n" 
                     << phanSize.y*1E3 <<" mm elevationally\n" 
                     <<phanSize.z*1E3 <<" mm axially\n";

        std::cout<< "The number of scatterers in the phantom is: " <<totalScatters << std::endl;
        if(buffer!=NULL) delete[] buffer;
	buffer = new scatterer[totalScatters];
	
	fpin.read((char*)buffer, sizeof(scatterer)*totalScatters);
	fpin.read( (char*) &numBsc, sizeof(int) );
	fpin.read( (char*) &freqStep, sizeof(double) );
        std::cout << "The number of backscatter coefficients stored is equal to: " <<numBsc << std::endl;	
	if(bscArray!=NULL) delete[] bscArray;
        std::cout << "Allocating the backscatter coefficient array \n";	
        bscArray = new double[numBsc];
        
        std::cout << "Reading in the backscatter coefficients" << std::endl;
	fpin.read( (char*) bscArray, numBsc*sizeof( double) );

	std::cout << "Creating a backscatter frequency array" << std::endl;
	if(bscFreqArray!=NULL) delete[] bscFreqArray;
	bscFreqArray = new double[numBsc];

        std::cout << "Reading in the backscatter frequencies" <<std::endl;	
	fpin.read( (char*) bscFreqArray, numBsc*sizeof( double) );

	
	return 1;
}


/*!  This function uses a file containing axial and lateral displacements to move each scatterer contained
 * in a phantom.
 */
void phantom::displaceAnsys(double* u, double* v, int nSize) //nSize is the number of points in the square matrix of displacement estimates
{
	//For all arrays corresponding to images/matrixes the ordering will go like:
	//1   4    7
	//2   5    8
	//3   6    9
	//
	double x,z, localX, localZ, phsize,dispx, dispz;
	//spacing is the coarser spacing between the x and z dimensions
	//same number of points in x and z direction
	(phanSize.x>phanSize.z)? phsize=phanSize.x:phsize=phanSize.z;
	double spacing=(nSize-1)/phsize;
	int i, roundX, roundZ, arrayIdx;

	
    for(i=0; i<totalScatters; i++) {
    	//important to remember coordinates go from 0 to size
	//get x and z for current scatterer in units of points
    	x=buffer[i].x*spacing;
	z= buffer[i].z*spacing;
        
        //for these if statements x and z should run from 0 to nSize - 1 if they are in bounds
        if( x<0 || x>= (nSize-1) || z < 0 || z >=(nSize-1) )
		{
			dispx=0; dispz=0;
		}
		else {
			//linearly interpolating in a 2 by 2 grid of points, ranging from x,z to x+1,z+1
			//u and v are 1-D arrays of displacements, a flattened 2-D array going by row Major indexing
            roundX=(int) x;
            roundZ=(int) z;
            arrayIdx = roundZ + roundX*nSize;        //index of lower left hand coordinate, assuming numbering goes
            localX = x - roundX;                 //1   4   7
            localZ = z - roundZ;                 //2   5   8
                                                 //3   6   9
            //bilinear interpolation using Lagrange polynomials                                                 
		dispx= (.5 - localX)*(.5 - localZ)*u[arrayIdx] + (.5 + localX)*(.5 - localZ)*u[arrayIdx + nSize] + (.5 + localX)*(.5 + localZ)*u[arrayIdx + nSize - 1] + (.5 - localX)*(.5 + localZ)*u[arrayIdx - 1];
		    dispz= (.5 - localX)*(.5 - localZ)*v[arrayIdx] + (.5 + localX)*(.5 - localZ)*v[arrayIdx + nSize] + (.5 + localX)*(.5 + localZ)*v[arrayIdx + nSize - 1] + (.5 - localX)*(.5 + localZ)*v[arrayIdx - 1];
        }

		buffer[i].x+=dispx*phsize;
		buffer[i].z+=dispz*phsize;
	}


}

/*!  This function gets scatterers located between two X coordinates.  It depends on the scatterers
 * contained in a phantom's scatterer array to be sorted by increasing x coordinate
 */
int phantom::getScattersBetween(double start, double end, scatterer*& buf)
{
	assert(start<end);
	int recStart = binSearch(start);
	int recEnd = binSearch(end);
	int recLen = recEnd - recStart;
	buf	= &buffer[recStart];
	return recLen;
}

/*!  This function is a standard binary search.
 */
int phantom::binSearch(double val)
{

	//finds the scatterer with the closest x coordinate to val
	int left = 0;
	int right = totalScatters-1;
	int center;
	if(buffer[0].x>=val)
		return 0;
	if(buffer[totalScatters-1].x<val)
		return totalScatters;
	while(right>=left)
	{
		center = (left + right)/2;
		if (buffer[center].x >= val)
			right = center-1;
		else
		left = center+1;
	}

	return left;
}

/*!  Return the phantom sound speed
 */
double phantom::soundSpeed()
{
	return c0;

}

/*!  Return the phantom's attenuation coefficient at frequency freq
 */
double phantom::attenuation(double freq)
{

	double fMhz = freq/1e6;
	double db = a0+a1*pow(fMhz, a2);
	return db*100.*log(10)/20;
}

/*!  Sort scatterers by increasing x coordinate using a quicksort algorithm
 */
void phantom::sortScatterer() {
	quickSort(buffer, 0, totalScatters-1);
}

/*!  A standard quick sort
 */
void phantom::quickSort(scatterer *A, int F, int L)
{
	int PivotIndex;
	if (F<L)
	{
		if(F==(L-1)) {
			if(A[F].x>A[L].x)
			{
				scatterer temp;
				Swap(A[F], A[L]);
			}
			return;
		}
		partition(A, F, L, PivotIndex);
		quickSort(A, F, PivotIndex-1);
		quickSort(A, PivotIndex+1, L);
	}
}

/*!  A component of the quicksort algorithm
 */
void phantom::partition(scatterer* A, int F, int L, int& PivotIndex)
{
	int pivotind=(int) (F+L)/2;
	scatterer temp;
	Swap(A[F], A[pivotind]);
	scatterer Pivot = A[F];
	int LastS1 = F;
	int FirstUnknown = F+1;

	for(; FirstUnknown<=L;++FirstUnknown)
	{
		if( A[FirstUnknown].x < Pivot.x)
		{
			++LastS1;
			Swap(A[FirstUnknown], A[LastS1]);
		}
	}

	Swap(A[F],A[LastS1]);

	PivotIndex = LastS1;

}
