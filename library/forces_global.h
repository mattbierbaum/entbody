#include <math.h>
#include "../util.h"

#define FORCE_DAMPING \
do {                                            \
    double speed= FORCE_DAMPING_SPEED;          \
    double damp = FORCE_DAMPING_COEFF;          \
    double vlen = v[2*i+0]*v[2*i+0]             \
                + v[2*i+1]*v[2*i+1];            \
    if (vlen > 1e-6){                           \
        f[2*i+0] -= damp*(vlen-speed)*v[2*i+0]/vlen;    \
        f[2*i+1] -= damp*(vlen-speed)*v[2*i+1]/vlen;    \
    }                                           \
} while (0);


#define FORCE_THERMAL \
do {                                \
    f[2*i+0] += Tglobal*(ran_ran2()-0.5); \
    f[2*i+1] += Tglobal*(ran_ran2()-0.5); \
} while(0);


#define FORCE_KICK \
do {                \
    f[2*i+0] += o[2*i+0];   \
    f[2*i+1] += o[2*i+1];   \
} while(0);


#define FORCE_GRAVITY \
do {                        \
    double g = FORCE_GRAVITY_G;         \
    if (t == RED)           \
        f[2*i+1] += 0.5*g;  \
    else                    \
        f[2*i+1] += g;      \
} while(0);
    
