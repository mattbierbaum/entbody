#include <math.h>
#include "../util.h"

#define INIT_RANDOM \
do {                                \
    long i;                         \
    double radius  = 1.0;           \
    L = sqrt(pi*radius*radius*N);   \
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


