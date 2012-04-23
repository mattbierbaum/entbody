#ifndef __UTIL_H__
#define __UTIL_H__

#define pi      3.141592653589
#define BLACK   0
#define RED     1
#define EPSILON DBL_EPSILON

typedef unsigned long long int ullong;

// random number generator functions
void   ran_seed(long j);
double ran_ran2();

// neighbor list functions
void coords_to_index(double *x, int *size, int *index, double L);
int  mod_rvec(int a, int b, int p, int *image);

// random utility functions
double mymod(double a, double b);

#endif
