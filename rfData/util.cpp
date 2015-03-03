#include "util.h"
/*!The logic of Fresnel class is to perform fresnel integral once over wide interval,
 *and then to do 1-D interpolation to get further values of the integral.
 *Fresnel integral is symmetric
 */
fresnelInt::fresnelInt()
{
	
	EPSLON = 6.e-10;
	MAXIT = 500;
	FPMIN = 1.e-50;
	XMIN = 1.5;

        DELTAX = 4e-4;
	LARGELIM = 30.;
	BASESIZE = int(LARGELIM/DELTAX) + 1;
	eps = 3e-15;

	fresBase = new cplx[BASESIZE];
	assert(fresBase != NULL);

	for (int i=0; i<BASESIZE; i++)
	{

		fresBase[i] = fresnel(i*DELTAX);
	}
}

fresnelInt::~fresnelInt()
{
	delete[] fresBase;
}

/*!Code for evaluating fresnel integral straight from numerical recipes 1992 version
 * The Fresnel integrals are defined as follows:
 * \f$ S(x) = /int_0^x{ sin( \frac{\pi}{2} t^2), t } \f$
 * \f$ C(x) = /int_0^x{cos( \frac{\pi}{2} t^2) , t } \f$
 * This code returns a complex number equal to C(x) + i*S(x)
 *
 * By Euler's formula this is equivalent to returning
 * \f$ /int_0^x{ exp(i \frac{ \pi}{2} t^2), t } \f$
 */
cplx fresnelInt::fresnel(double x)
{
	int k, n, odd;
	double a, ax, fact, pix2, sign, sum, sumc, sums, term, test;
	cplx b, cc, d, h, del, cs;
	double s, c;

	ax = fabs(x);
	if (ax < sqrt(FPMIN))  //Input is small enough to call 0
	{
		s = 0.0;
		c = ax;
	}
	else if (ax <= XMIN)  //Use series approximation
	{
		sum = sums = 0;
		sumc = ax;
		sign = 1.0;
		fact = M_PI/2.*ax*ax;
		odd = true;
		term = ax;
		n = 3;
		for (k=1; k<=MAXIT; k++)
		{
			term *= fact/k;
			sum += sign*term/n;
			test = fabs(sum)*EPSLON;
			if (odd)
			{
				sign = -sign;
				sums = sum;
				sum = sumc;
			}
			else
			{
				sumc = sum;
				sum = sums;
			}
			if (term < test) break;
			odd = !odd;
			n += 2;
		}
		assert(k <= MAXIT);
		s = sums;
		c = sumc;
	}
	else      //Evaluate continued fraction using Lentz's method
	{
		pix2 = M_PI*ax*ax;
		b = 1.0 -imUnit*pix2;
		cc = 1.0/FPMIN;
		d = h = cplxOne/b;
		n = -1;
		for (k=2; k<=MAXIT; k++)
		{
			n += 2;
			a = -n*(n+1);
			b += 4;
			d = 1./((a*d)+b);
			cc = b + a/cc;
			del = cc*d;
			h *= del;
			if (fabs(del.real() - 1.0)+fabs(del.imag()) < EPSLON) break;
		}
		if(k >= MAXIT)

		assert(k <= MAXIT);
		h *= (ax-imUnit*ax);
		cs = (.5 + imUnit*.5)*(cplxOne- ( cos(.5*pix2) + imUnit*sin(.5*pix2)) *h);
		c = cs.real();
		s = cs.imag();
	}
	if (x<0)
	{
		c = -c;
		s = -s;
	}
	cplx returnValue;
	returnValue = c + imUnit*s;
	return returnValue;
}




/*!The fast Fresnel function uses linear interpolation of a table to find the
 * fast Fresnel integral
 */
cplx fresnelInt::fastFresnel(double x)
{
	double absx = fabs(x);
	cplx res;

	if (absx > LARGELIM)
	{
		res =	(.5+1./(M_PI*absx)*sin(M_PI/2*absx*absx)) + imUnit*(.5-1./(M_PI*absx)*cos(M_PI/2*absx*absx));
	}
	else //use linear interpolation with pre-calculated look up table
	{
		int index = int(absx/DELTAX);
		res = fresBase[index];
		double ratio = (absx - (index*DELTAX))/DELTAX;
		res = ratio*fresBase[index+1] + (1 - ratio)*fresBase[index];
	}

	if (x < 0)
	{
		res = -1.00*res;
	}
	return res;
}


/*!given two arrays, x and y, representing evenly spaced samples of a function
along with the size of the arrays, return the value of the function at the
intermediate point x0
 */
double fresnelInt::interp1(double *x, double *y, int size, double x0)
{
	int left = 0;
	int right = size-1;
	int center;

	while(right >= left)
	{
		center = (left + right)/2;
		if (x[center] >= x0)
			right = center-1;
		else
			left = center+1;
	}
	if(right<0)
		return(y[0]);
	else if(left>size-1)
		return(y[size-1]);
	else
		return (y[left]+(y[right]-y[left])/(x[right]-x[left])*(x0-x[left]));
}


