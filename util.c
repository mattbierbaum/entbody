#include <math.h>
#include "util.h"

ullong vseed;
ullong vran;

//=================================================
// random number generator 
//=================================================
void ran_seed(long j){
  vseed = j;  vran = 4101842887655102017LL;
  vran ^= vseed; 
  vran ^= vran >> 21; vran ^= vran << 35; vran ^= vran >> 4;
  vran = vran * 2685821657736338717LL;
}

double ran_ran2(){
    vran ^= vran >> 21; vran ^= vran << 35; vran ^= vran >> 4;
    ullong t = vran * 2685821657736338717LL;
    return 5.42101086242752217e-20*t;
}


//====================================================
// neighbor list helper functions
//====================================================
inline void coords_to_index(double *x, int *size, int *index, double L){   
    index[0] = (int)(x[0]/L  * size[0]);
    index[1] = (int)(x[1]/L  * size[1]);
}

inline int mod_rvec(int a, int b, int p, int *image){
    *image = 1;
    if (b==0) {if (a==0) *image=0; return 0;}
    if (p != 0){
        if (a>b)  return a-b-1;
        if (a<0)  return a+b+1;
    } else {
        if (a>b)  return b;
        if (a<0)  return 0;
    }
    *image = 0;
    return a;
}


//==========================================
// random things
//==========================================
inline double mymod(double a, double b){
  return a - b*(int)(a/b) + b*(a<0);
}
