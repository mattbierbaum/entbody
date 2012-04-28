#ifndef __UTIL_H__
#define __UTIL_H__

#define NMAX    50
#define pi      3.141592653589
#define BLACK   0
#define RED     1
#define EPSILON FLT_EPSILON

typedef unsigned long long int ullong;

// random number generator functions
void   ran_seed(long j);
float ran_ran2();

// neighbor list functions
void coords_to_index(float *x, int *size, int *index, float L);
int  mod_rvec(int a, int b, int p, int *image);

// random utility functions
float mymod(float a, float b);

#endif
