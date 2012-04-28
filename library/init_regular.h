#include <math.h>
#include "../util.h"

#define INIT_RANDOM \
do {                                \
    long i;                         \
    float radius  = 1.0;            \
    L = INIT_RANDOM_FRACTION*sqrt(pi*radius*radius*N);   \
                                    \
    for (i=0; i<N; i++){            \
        rad[i] = radius;            \
        x[2*i+0] = L*ran_ran2();    \
        x[2*i+1] = L*ran_ran2();    \
                                    \
        type[i] = BLACK;            \
        v[2*i+0] = 0.0;             \
        v[2*i+1] = 0.0;             \
        if (i==0) {                 \
            type[i] = RED;          \
            rad[i] = 1*radius;      \
        }                           \
     }                              \
} while(0);

#define INIT_CUBIC_LATTICE                          \
do {                                                \
    long i;                                         \
    float radius  = 1.0;                            \
    L = sqrt(pi*radius*radius*N);                   \
                                                    \
    float f = INIT_RAYLEIGHTAYLOR_FRACTION;         \
    for (i=0; i<N; i++){                            \
        rad[i] = radius;                            \
        x[2*i+0] = L - mymod((float)2*i, L);        \
        x[2*i+1] = L - (int)(4*i/L);                \
        if (x[2*i+1] < 0)  x[2*i+1] = radius;       \
        if (x[2*i+1] >= L) x[2*i+1] = L-radius;     \
                                                    \
        v[2*i+0] = 0.0;                             \
        v[2*i+1] = 0.0;                             \
                                                    \
        type[i] = BLACK;                            \
    }                                               \
} while(0);
