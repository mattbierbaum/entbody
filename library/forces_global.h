#include <math.h>
#include "../util.h"

#define FORCE_DAMPING \
do {                                            \
    float speed= FORCE_DAMPING_SPEED;           \
    float damp = FORCE_DAMPING_COEFF;           \
    float vlen = vx*vx + vy*vy;                 \
    if (vlen > 1e-6){                           \
        fx -= damp*(vlen-speed)*vx/vlen;        \
        fy -= damp*(vlen-speed)*vy/vlen;        \
    }                                           \
} while (0);


#define FORCE_THERMAL               \
do {                                \
    fx += Tglobal*(ran_ran2()-0.5); \
    fy += Tglobal*(ran_ran2()-0.5); \
} while(0);


#define FORCE_KICK  \
do {                \
    fx += ox;       \
    fy += oy;       \
} while(0);


#define FORCE_GRAVITY           \
do {                            \
    float g = FORCE_GRAVITY_G;  \
    if (type[i] == RED)         \
        fy += 0.1*g;            \
    else                        \
        fy += g;                \
} while(0);
    
