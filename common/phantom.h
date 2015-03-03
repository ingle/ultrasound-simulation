#ifndef PHANTOM_H
#define PHANTOM_H

/*! \brief A structure that holds the x,y,z position of a scatterer.  More info to be added.
 */
struct scatterer
{
	double x;
	double y;
	double z;

};

/*! \brief An x,y,z vector for working with the phantom class.
 */
struct myVector
{
       double x;
       double y;
       double z;
};
/*! \brief This class encompasses the attenuation, sound speed, and backscatter coefficients of an object to be imaged.  
  */
class phantom
{
	public:
        phantom();
	~phantom();
		
        void createUniformPhantom(const myVector&, double, double, double, double, double, char*);
	
	//saving and loading phantoms for future use
        int savePhantom(char*);
	int loadPhantom(char*);
		
	//displacing the scatterer positions.  Used for compressions and elastography.
        void displaceAnsys(double*, double*, int);
		
        //finding scatterers
        int getScattersBetween(double, double, scatterer*&);
        int binSearch(double);

	//functions for getting info about the phantoms
	myVector getPhanSize();

	//return sound speed or attenuation as a function of frequency
	 double soundSpeed();
	 double attenuation(double freq);

	//sorting scatterers
	void sortScatterer();
	void quickSort(scatterer *A, int F, int L);
	void partition(scatterer* A, int F, int L, int& PivotIndex);

        //obtaining backscatter coefficient
	double giveBsc(double freq);
	 
	private:
	double c0;  //Constant phantom sound speed
	double a0, a1, a2;      //Frequency dependent attenuation
	scatterer* buffer;
	myVector phanSize;
	int totalScatters;
	
	//function for calculating backscatter coefficients eventually, for now just reads in a list
	void readBscFromFile(char*);
	
	//info about backscatter coefficients
	double* bscArray;
	double* bscFreqArray;
	int numBsc;
	double freqStep;
};

#endif
