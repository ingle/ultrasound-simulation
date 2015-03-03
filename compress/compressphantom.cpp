#include <iostream>
#include <fstream>
#include "phantom.h"

using namespace std;


int main(int argc, char* argv[])
{

	if(argc<2) {
		cout<<"Error! An input file is needed "<<endl;
		return -1;
	}
	char inphanfile[60], outphanfile[60], displfile[60], line[60];
	int dispSize;

	FILE *fpinput;
	
	if ((fpinput=fopen(argv[1], "r"))==NULL)
		{
			cout<<"Can't find input file!"<<endl;
			return -1;
		}

	while (fgetc(fpinput)!=':');
	fscanf(fpinput, "%s", inphanfile);
	

	while (fgetc(fpinput)!=':');
	fscanf(fpinput, "%s", outphanfile);

	while (fgetc(fpinput)!=':');
	fscanf(fpinput, "%s", displfile);

	while (fgetc(fpinput)!=':');
	fscanf(fpinput, "%d", &dispSize);

	
	phantom target;
	cout<<"Pre Phantom File name is: " <<inphanfile <<endl;
	cout <<"Output Phantom file name will be: " <<outphanfile <<endl;
	cout <<" The displacement matrix has a size of:"  << dispSize <<endl;
	
	target.loadPhantom(inphanfile);


 	  //this bit of the code allocates memory to contain a dispSize by dispSize array of x and y displacements
		double *u, *v;
        u= new double[dispSize*dispSize];
		v= new double[dispSize*dispSize];
	
		ifstream fpdisp;
		fpdisp.open(displfile, ios::binary);
		if(!fpdisp.is_open())
		{
			cout<<"Error! Can't find file "<<displfile<<endl;
			return -1;
		}
		
        //here the displacements are read in from the .dis file
        fpdisp.read((char*)u,dispSize*dispSize*sizeof(double) );
		fpdisp.read((char*)v,dispSize*dispSize*sizeof(double) );
		fpdisp.close();


		target.displaceAnsys(u,v,dispSize);

        target.savePhantom(outphanfile);

		return 0;
	
}


