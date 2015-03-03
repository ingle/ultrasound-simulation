#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "phantom.h"

using namespace std;



int main(int argc, char* argv[])
{
	myVector  inSize;
	double density, c0, a0, a1, a2;
	char strPhanFile[60], strBscFile[60];
	int scatterType;  //at the moment an empty feature, need to include backscatter coefficient calculations
	
	FILE *fpinput;


	if(argc<2) {
			cout<<"Error! An input file is needed"<<endl;
			exit(-1);
		}


	if ((fpinput=fopen(argv[1], "r"))==NULL)
		{
			cout<<"Can't find input file!"<<endl;
			exit(-1);
		}


	while (fgetc(fpinput)!=':');
	fscanf(fpinput, "%lf, %lf, %lf", &inSize.x, &inSize.y, &inSize.z);


	while (fgetc(fpinput)!=':');
	fscanf(fpinput, "%lf", &density);

	while (fgetc(fpinput)!=':');
	fscanf(fpinput, "%lf", &c0);

	while (fgetc(fpinput)!=':');
	fscanf(fpinput, "%lf, %lf, %lf", &a0, &a1, &a2);
	
	while (fgetc(fpinput)!=':');
	fscanf(fpinput, "%s", strBscFile);

	while (fgetc(fpinput)!=':');
	fscanf(fpinput, "%s", strPhanFile);

	cout << "The x size is " << inSize.x << endl;
	cout << "The y size is " << inSize.y << endl;
	cout << "The z size is " << inSize.z <<endl;
	cout << "The density is " << density << endl;
	cout << "The filename is " << strPhanFile << endl;
	cout << "The sound Speed is " << c0 << endl;
	cout << "The attenuation is " << a0 <<" " << a1 << "  " << a2 << endl;

	phantom target;
	target.createUniformPhantom(inSize,density, c0, a0, a1, a2, strBscFile );
        target.savePhantom(strPhanFile);

}
