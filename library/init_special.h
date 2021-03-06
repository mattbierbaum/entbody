#include <math.h>
#include "../core/util.h"

#define INIT_BRAZILNUTS             \
do {                                \
    long i;                         \
    float radius  = 1.0;            \
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


#define INIT_CENTRAL_CIRCLE                 \
do {                                        \
    long i;                                 \
    float radius = 1.0;                     \
    L = 1.03*sqrt(pi*radius*radius*N);      \
    for (i=0; i<N; i++){                    \
        double tx = L*ran_ran2();           \
        double ty = L*ran_ran2();           \
        double tt = 2*pi*ran_ran2();        \
                                            \
        rad[i] = radius;                    \
        x[2*i+0] = tx;                      \
        x[2*i+1] = ty;                      \
        double dd = sqrt((tx-L/2)*(tx-L/2) +\
            (ty-L/2)*(ty-L/2));             \
        double rad = sqrt(0.15*L*L / pi);   \
        if (dd < rad)                       \
            type[i] = RED;                  \
        v[2*i+0] = 0.0f;                    \
        v[2*i+1] = 0.0f;                    \
    }                                       \
} while(0);                                 \

