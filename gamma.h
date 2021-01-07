
#ifndef __GAMMA_H__
#define __GAMMA_H__

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>

#define isnan(x) _isnan(x)

double gamma(double x);
double lgamma(double x);
double digamma(double x);
double trigamma(double x);

#define psi(x)	digamma(x)
#define ppsi(x)	trigamma(x)

int *ivec(int n);
double *dvec(int n);
void updatesort(int n, int *x, int *indx, int *revindx, int ip, int im);
void insertionsort(int n, double *x, int *indx);
void isort(int n, int *x, int direction, int *indx);
void dsort(int n, double *x, int direction, int *indx);
double vec_sum(int N, double*x);
double etime();

static int icomp(const void *, const void *); /* comparison for isort */
static int dcomp(const void *, const void *); /* comparison for dsort */
static int    *icomp_vec;                     /*  data used for isort */
static double *dcomp_vec;                     /*  data used for dsort */

#endif