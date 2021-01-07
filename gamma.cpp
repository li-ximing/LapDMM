
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "gamma.h"

/* The digamma function is the derivative of gammaln.

Reference:
J Bernardo,
Psi ( Digamma ) Function,
Algorithm AS 103,
Applied Statistics,
Volume 25, Number 3, pages 315-317, 1976.

From http://www.psc.edu/~burkardt/src/dirichlet/dirichlet.f
(with modifications for negative numbers and extra precision)
*/

double gamma(double x)
{
	int i;
	double y,t,s,u;
	static double a[11]={ 0.0000677106,-0.0003442342,
		0.0015397681,-0.0024467480,0.0109736958,
		-0.0002109075,0.0742379071,0.0815782188,
		0.4118402518,0.4227843370,1.0};

	if (x<=0.0)
	{
		//		printf("err**x<=0!\n"); 
		return(-1.0);
	}

	y=x;
	if (y<=1.0)
	{ 
		t=1.0/(y*(y+1.0)); 
		y=y+2.0;
	}
	else if (y<=2.0)
	{ 
		t=1.0/y;
		y=y+1.0;
	}
	else if (y<=3.0) 
		t=1.0;
	else
	{ 
		t=1.0;
		while (y>3.0)
		{ 
			y=y-1.0; 
			t=t*y;
		}
	}
	s=a[0]; 
	u=y-2.0;
	for (i=1; i<=10; i++)
		s=s*u+a[i];
	s=s*t;

	return s;
}

double lgamma(double x ) 
{
	if (x < 12.0)
	{
		return log(fabs(gamma(x)));
	}

	// Abramowitz and Stegun 6.1.41
	// Asymptotic series should be good to at least 11 or 12 figures
	// For error analysis, see Whittiker and Watson
	// A Course in Modern Analysis (1927), page 252

	static const double c[8] =
	{
		1.0/12.0,
		-1.0/360.0,
		1.0/1260.0,
		-1.0/1680.0,
		1.0/1188.0,
		-691.0/360360.0,
		1.0/156.0,
		-3617.0/122400.0
	};
	double z = 1.0/(x*x);
	double sum = c[7];
	for (int i=6; i >= 0; i--)
	{
		sum *= z;
		sum += c[i];
	}
	double series = sum/x;

	static const double halfLogTwoPi = 0.91893853320467274178032973640562;
	double logGamma = (x - 0.5)*log(x) - x + halfLogTwoPi + series;    
	return logGamma;
}

/*
double digamma(double x)
{
	double result;
	static const double
		//		neginf = -1.0/0.0,
		c = 12,
		s = 1e-6,
		d1 = -0.57721566490153286,
		d2 = 1.6449340668482264365, // pi^2/6 
		M_PI = -3.14159265358979,
		s3 = 1./12,
		s4 = 1./120,
		s5 = 1./252,
		s6 = 1./240,
		s7 = 1./132,
		s8 = 691/32760,
		s9 = 1/12,
		s10 = 3617/8160;
	
	if(x < 0) {
		return digamma(1-x) + M_PI/tan(-M_PI*x);
	}
	// Use Taylor series if argument <= S 
	if(x <= s) return d1 - 1/x + d2*x;
	// Reduce to digamma(X + N) where (X + N) >= C 
	result = 0;
	while(x < c) {
		result -= 1/x;
		x++;
	}
	// Use de Moivre's expansion if argument >= C 
	// This expansion can be computed in Maple via asympt(Psi(x),x)
	if(x >= c) {
		double r = 1/x;
		result += log(x) - 0.5*r;
		r *= r;
		result -= r * (s3 - r * (s4 - r * (s5 - r * (s6 - r * s7))));
	}
	return result;
}
*/

/*
double lgamma(double x)
{
double z=1/(x*x);

x=x+6;
z=(((-0.000595238095238*z+0.000793650793651)
*z-0.002777777777778)*z+0.083333333333333)/x;
z=(x-0.5)*log(x)-x+0.918938533204673+z-log(x-1)-
log(x-2)-log(x-3)-log(x-4)-log(x-5)-log(x-6);
return z;
}
*/

double digamma(double x)
{
double p;
x=x+6;
p=1/(x*x);
p=(((0.004166666666667*p-0.003968253986254)*p+
0.008333333333333)*p-0.083333333333333)*p;
p=p+log(x)-0.5/x-1/(x-1)-1/(x-2)-1/(x-3)-1/(x-4)-1/(x-5)-1/(x-6);
return p;
}

// The trigamma function is the derivative of the digamma function.

double trigamma(double x)
{
	double p;
	int i;

	x=x+6;
	p=1/(x*x);
	p=(((((0.075757575757576*p-0.033333333333333)*p+0.0238095238095238)
		*p-0.033333333333333)*p+0.166666666666667)*p+1)/x+0.5*p;
	for (i=0; i<6 ;i++)
	{
		x=x-1;
		p=1/(x*x)+p;
	}
	return(p);
}

double *dvec(int n) //
{
  double *x = new double[n];
  return x;
}

int *ivec(int n) //
{
  int *x = new int[n];
  return x;
}

void updatesort(int n, int *x, int *indx, int *revindx, int ip, int im) //
{
	int tmp;

	// INCREMENT
	// did ++ get bigger than prev?
	ip = revindx[ip];
	while (ip>0 && x[indx[ip]] > x[indx[ip-1]]) {
		// swap indx
		tmp        = indx[ip];
		indx[ip]   = indx[ip-1];
		indx[ip-1] = tmp;
		// swap revindx
		tmp                 = revindx[indx[ip]];
		revindx[indx[ip]]   = revindx[indx[ip-1]];
		revindx[indx[ip-1]] = tmp;
		ip--;
	}

	// DECREMENT
	// did -- get smaller than next?
	im = revindx[im];
	while (im<n-1 && x[indx[im]] < x[indx[im+1]]) {
		// swap indx
		tmp        = indx[im];
		indx[im]   = indx[im+1];
		indx[im+1] = tmp;
		// swap revindx
		tmp                 = revindx[indx[im]];
		revindx[indx[im]]   = revindx[indx[im+1]];
		revindx[indx[im+1]] = tmp;
		im++;
	}
}


void insertionsort(int n, double *x, int *indx) // descending order
{
	int tmp, i, k;

	for (i=1; i<n; i++) {
		for (k=i; k>0 && x[indx[k]]>x[indx[k-1]]; k--) {
			tmp = indx[k];
			indx[k] = indx[k-1];
			indx[k-1] = tmp;
		}
	}
}

void isort(int n, int *x, int direction, int *indx) // +1 ascending order -1 descending order
{
	int i;
	icomp_vec = ivec(n);
	for (i=0; i<n; i++) {
		icomp_vec[i] = direction*x[i];
		indx[i] = i;
	}
	qsort(indx, n, sizeof(int), icomp);
	delete [] icomp_vec;
//	mxFree(icomp_vec);
}

static int icomp(const void *pl, const void *p2)
{
	int i = * (int *) pl;
	int j = * (int *) p2;
	return (icomp_vec[i] - icomp_vec[j]);
}

void dsort(int n, double *x, int direction, int *indx) // +1 ascending order -1 descending order
{
	int i;
	dcomp_vec = dvec(n);
	for (i=0; i<n; i++) {
		dcomp_vec[i] = direction*x[i];
		indx[i] = i;
	}
	qsort(indx, n, sizeof(int), dcomp);
	delete [] dcomp_vec;
//	mxFree(dcomp_vec);
}

static int dcomp(const void *p1, const void *p2)
{
	int i = * (int *) p1;
	int j = * (int *) p2;
	return dcomp_vec[i] > dcomp_vec[j] ? 1:-1;
}

double vec_sum(int N, double* x)
{
	double result = 0.0;

	for(int n=0 ; n<N ; n++)
	{
		result += x[n];
	}

	return result;
}

/* CPU time */
double etime() //
{
	static double last_clock = 0;
	static double now_time = 0;
	last_clock = now_time;
	now_time = (double) clock ();
	return (double) (now_time - last_clock) / CLOCKS_PER_SEC;
}
