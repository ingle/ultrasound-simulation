#ifndef __PRESSUREFIELD_H__
#define __PRESSUREFIELD_H__

#define FULL_APERTURE -1
#include <stdio.h> 
#include <assert.h>
#include "phantom.h"

class array;
class fieldBuffer;
class phantom;


//structure used to describe a single rectangular element
struct singleGeom
{
	double width;
	double length;

};



class array
{
	friend class fieldBuffer;


	public:
		array(singleGeom, double, int, double);
		//(geom, spacing between elements, element number, assumed sound speed)


		~array(); //destructor


		void setTransFocus(double, double, double);
		//(focal distance, angle of steering, F number, frequency) 
		//set F<0 to active all the elements

		void setRecFocus(double, double, double);
		//(focal distance, angle of steering, F number, frequency)
		//set F<0 to active all the elements

		double Spacing() { return spacing;};
		//get the spacing
		double trsFnum() { return trnsFnum;};
		double recFnum() { return recvFnum;};
		void setTrsFnum(double);
		void setRecFnum(double);
	private:
		void setFocus(cplx*, double, double, double);
		//common routine to set lateral focus
		//(phase factor[], focal distance, F number, frequency);

		singleGeom geom; //hold the elevational geometry info
		double spacing; //spacing between elements
		int eleCnt; //number of elements
		double trnsFnum, recvFnum;
		cplx *transPhase; //transmit lateral phase factors 
		cplx *recPhase; //receive lateral phase factors
		double assumedSoundSpeed;

};




class fieldBuffer   //incident pressure field associated with a phantom
{
	public:
		fieldBuffer(double, double, const vector&, phantom*, double, array*, double);
		//constructor (transmit focus, beam width, grid step );

		~fieldBuffer(); //destructor

		void calculateBufferField(double);
		//calculate the buffer at frequency (freq)
		
		cplx bufferField(const vector&);		
		//get the pressure field at (location)
		void beamProfile();

		vector giveCenter() {return center; };
		double giveImageDepth() {return size.z;};

		double giveSoundSpeed() {return transducer->assumedSoundSpeed;};
		//get the pressure field contribution due to a single element
		cplx getSingleElementField(const vector&, const singleGeom&, const cplx&);

		
		//given a scatterer, return the coordinates relative to the transducer
		vector phantomCoordinateToPressureCoordinate(const scatterer&, int);
	private:
		double transFocus; //Transmit focus
		vector size; //the field size to be calculated
		vector step; //the grid step
		vector center; //center of the field
		cplx K; //cplx wavenumber
		phantom* target;
		int xLen; //x dimension
		int xLenExtra; //
		int yLen; //y dimension
		int zLen; //z dimension
		int denseFactor; //densefactor is the number of calculated points along each element

		int arrayLineSize; //buffer line size
		int arrayPlaneSize; //buffer plane size


		cplx *singleRowTransField; //single row of the transmit field, constant depth
		cplx *singleRowRecField; //single row of the receive field, constant depth
	
		cplx *arrayField; //resulting buffer field */
		array* transducer;
		double assumedSoundSpeed;
		fresnelInt *fres;
		double phantomGap;
};


#endif
