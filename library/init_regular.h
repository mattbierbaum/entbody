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
    float f = INIT_CUBIC_LATTICE_A;                 \
    L = f*2*radius*((int)sqrt(N));                  \
                                                    \
    float ep = EPSILON;                             \
    float sz = 2*radius*f;                          \
    int side = (int)L/sz;                           \
                                                    \
    for (i=0; i<N; i++){                            \
        rad[i] = radius;                            \
        x[2*i+0] = sz * (int)(i%side) + ep;         \
        x[2*i+1] = sz * (int)(i/side) + ep;         \
        v[2*i+0] = 0.0;                             \
        v[2*i+1] = 0.0;                             \
                                                    \
        type[i] = BLACK;                            \
    }                                               \
} while(0);
