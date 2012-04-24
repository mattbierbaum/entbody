#include <math.h>
#include "../util.h"

#define INIT_BRAZILNUTS \
do {                                \
    long i;                         \
    double radius  = 1.0;           \
    L = sqrt(5*pi*radius*radius*N); \
                                    \
    for (i=0; i<N; i++){            \
        rad[i] = radius;            \
        x[2*i+0] = L*ran_ran2();    \
        x[2*i+1] = L*ran_ran2();    \
                                    \
        type[i] = BLACK;            \
        if (i %2 ==0)               \
            rad[i] = 2*radius;      \
        v[2*i+0] = 0.0;             \
        v[2*i+1] = 0.0;             \
        if (i==0) {                 \
            type[i] = RED;          \
            rad[i] = 5*radius;      \
        }                           \
     }                              \
} while(0);

#define INIT_RAYLEIGHTAYLOR \
do {                                                \
    long i;                                         \
    double radius  = 1.0;                           \
    L = 2*sqrt(pi*radius*radius*N);                 \
                                                    \
    double f = 0.15;                                \
    for (i=0; i<N; i++){                            \
        rad[i] = radius;                            \
        x[2*i+0] = L - mymod((double)i, L);         \
        x[2*i+1] = L - (double)(i/L);               \
        if (x[2*i+1] < 0)  x[2*i+1] = 0.0;          \
        if (x[2*i+1] >= L) x[2*i+1] = L;            \
                                                    \
        v[2*i+0] = 0.0;                             \
        v[2*i+1] = 0.0;                             \
                                                    \
        type[i] = BLACK;                            \
        if (i < f*N)                                \
            type[i] = RED;                          \
    }                                               \
} while(0);