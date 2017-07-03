#include <iostream>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include "util.h"
#include "pressureField.h"
#include "phantom.h"




using namespace std;
int main(int argc, char* argv[])
{

	//x is lateral, y is elevational, z is axial direction
	if(argc<2)
	{
		cout<<"Error! An input file is needed"<<endl;
		exit(-1);
	}



	singleGeom geom;
	vector step;	
	int count, beamlines;
	double beamspacing,spacing, transfocus, maxfreq, beamWidth, machineSoundSpeed, phantomGap, samplingFrequency;
	char phantomfile[60], outrffile[60], backscatterfile[60];
	fresnelInt *fres;
	array *transducer;
	int success;


	FILE *fpinput;
	cout<<argv[1]<<endl;

	if ((fpinput=fopen(argv[1], "r"))==NULL)
	{
		cout<<"Can't find input file!"<<endl;
		exit(-1);
	}
	
        cout << "Reading Input File" << endl;

	while (fgetc(fpinput)!=':');
	success = fscanf(fpinput, "%lf, %lf", &geom.width, &geom.length);

	while (fgetc(fpinput)!=':');
	success = fscanf(fpinput, "%lf, %d", &spacing, &count);

	while (fgetc(fpinput)!=':');
	success = fscanf(fpinput, "%lf", &transfocus);
	if(transfocus<0)
		cout<<"Focusing disabled."<<endl;
	
	while (fgetc(fpinput)!=':');
	success = fscanf(fpinput, "%lf", &beamWidth);
	
	while (fgetc(fpinput)!=':');
	success = fscanf(fpinput, "%lf, %lf, %lf", &step.x, &step.y, &step.z);

	
	while (fgetc(fpinput)!=':');
	success = fscanf(fpinput, "%d", &beamlines);

	while (fgetc(fpinput)!=':');
	success = fscanf(fpinput, "%lf", &beamspacing);
	
	while (fgetc(fpinput)!=':');
	success = fscanf(fpinput, "%lf", &maxfreq);

	cout << "The maximum simulated frequency is:  " << maxfreq << endl;
	while (fgetc(fpinput)!=':');
	success = fscanf(fpinput, "%s", phantomfile);
	
	
	while (fgetc(fpinput)!=':');
	success = fscanf(fpinput, "%s", outrffile);


	while(fgetc(fpinput)!=':');
	success = fscanf(fpinput, "%lf", &machineSoundSpeed);

        while (fgetc(fpinput)!=':');
	success = fscanf(fpinput, "%lf", &phantomGap);

	fclose(fpinput);


	//Load the phantom and generate the incident pressure field
	phantom target;
	int loaded = target.loadPhantom(phantomfile);
    if(!loaded){
    cout << "Phantom file not loaded" << endl;
    return -1;
    }
	

        //Initialize a class for holding transducer information, pressure field, fresnel integral
	transducer = new array(geom, spacing, count, machineSoundSpeed);
	assert(transducer != NULL);
	fieldBuffer pressure(transfocus, beamWidth, step, &target, machineSoundSpeed, transducer, phantomGap);
	
    
	////////////////////////////////////////////////////
	////////////////////////////////////////////////////
	/////////////Calculate the image FFT////////////////
	////////////////////////////////////////////////////
	////////////////////////////////////////////////////

	//get necessary delta freq and frequency points
	double beamOnTime = (pressure.giveImageDepth()*2)/pressure.giveSoundSpeed();
	//The maximum simulated frequency is the sampling frequency.
	//This is from the relationship deltaF*deltaT = 1/N	
	int freqPoints = static_cast<int>( beamOnTime*maxfreq);
	double freqStep = maxfreq/freqPoints;

	//initialize the coef matrix
	cplx* fftCoef = new cplx[freqPoints*beamlines];
	assert(fftCoef != NULL);


	vector fieldCenter;
	fieldCenter = pressure.giveCenter();

	myVector phanSize;
	phanSize = target.getPhanSize();

	//Need to be sure scatterers are sorted before imaging is performed
	target.sortScatterer();
	//set to 0
	for(int fp=0; fp<freqPoints; fp++)
		for(int il=0; il<beamlines; il++)
				fftCoef[fp +il*freqPoints] = cplxZero;

		time_t t0, t1;
		t0 = time(NULL);

		//loop through freq domain, skip DC frequency as contribution is zero there
		for(int fIndex=1; fIndex<freqPoints; fIndex++)
		{
			double freq = fIndex*freqStep; //Hz
			cout << "The backscatter coefficient at: " << freq/1E6 <<" MHz is: "<< target.giveBsc(freq/1E6) << endl;
			//get the next buffer field
			pressure.calculateBufferField(freq);
			vector loc;
			//loop through image lines
			for(int i=0; i<beamlines; i++)
			{
				//left end, right end, and center of beam
				double leftEnd = i*beamspacing;
				double rightEnd = leftEnd + beamWidth;
				double beamCenter = (leftEnd + rightEnd)/2;
				scatterer *pos;
				int cnt = target.getScattersBetween(leftEnd, rightEnd, pos);

				//loop through each scatterer in beam
				for(int j=0; j<cnt; j++)
				{
				       loc = pressure.phantomCoordinateToPressureCoordinate(pos[j], i);
				       cplx a0 = pressure.bufferField(loc); //get pressure field at location
				       fftCoef[fIndex + i*freqPoints] += a0*target.giveBsc(freq);

					//a0 is pi and ps, incident and scattered pressure multiplied
				}

				//take care of constants
				cplx factor = freq*imUnit;
				fftCoef[fIndex + i*freqPoints] *= factor;
			}

		 //track how long each iteration takes
	    	t1 = time(NULL);
			cout << (fIndex+1) << '/' << freqPoints << " completed: "
			 << freq/1e6 << "MHz, " << t1-t0 << " sec used" << endl;


		}
	
	//////////////////////////////////////
	//////////////////////////////////////
	/////////////SAVE OUTPUT//////////////
	//////////////////////////////////////
	//////////////////////////////////////
	//////////////////////////////////////
	ofstream fp(outrffile, ios::binary );

	if(!fp.is_open() )
	    {
	       cout << "Failure to write RF data file named: " << outrffile << std::endl;
	       exit(EXIT_FAILURE);
	    }

	//First write the freq step, number of points, number of lines as double, int, int
	fp.write( (char*) &freqStep , sizeof(double) );
	fp.write( (char*) &freqPoints, sizeof(int) );
	fp.write( (char*) &beamlines, sizeof(int) );	
	
	double *realSignal=new double[freqPoints*beamlines];
	double *imagSignal=new double[freqPoints*beamlines];


	for (int j=0; j<beamlines; j++)
		{
			int k;
			realSignal[0 + j*freqPoints]=0.0;
			imagSignal[0 + j*freqPoints]=0.0;
			for(k=0; k<freqPoints; k++)
			{
				realSignal[k + j*freqPoints] = fftCoef[k + j*freqPoints].real();
				imagSignal[k + j*freqPoints] = fftCoef[k + j*freqPoints].imag();
			}


		}


	fp.write((char*)realSignal, sizeof(double)*freqPoints*beamlines );
	fp.write((char*)imagSignal, sizeof(double)*freqPoints*beamlines );
	fp.close();
	delete realSignal;
	delete imagSignal;
	delete transducer;

}
