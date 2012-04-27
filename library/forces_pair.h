#include <math.h>
#include "../util.h"

#define FORCE_HERTZ                             \
do {                                            \
    float r0 = trad+rad[tn];                    \
    float l = sqrt(dist);                       \
    float co = FORCE_HERTZ_EPSILON              \
            * (1-l/r0)*(1-l/r0) * (l<r0);       \
    fx += -co * dx[0];                          \
    fy += -co * dx[1];                          \
} while (0); 


#define FORCE_MORSE                         \
do {                                        \
    float e = FORCE_MORSE_EPSILON;          \
    float a = FORCE_MORSE_ALPHA;            \
    float r0 = trad+rad[tn];                \
    float rcut = CONST_CUTOFF_FACTOR*r0;    \
    float fex  = exp(-a*(rcut-r0));         \
    float fcut = 2*e*a * (1-fex)*fex;       \
                                            \
    float l = sqrt(dist);                   \
    float ex = exp(-a*(l-r0));              \
    float co = (2*e*a * (1-ex)*ex-fcut)     \
  *(l<rcut);                                \
                                            \
    fx += co * dx[0]/l;                     \
    fy += co * dx[1]/l;                     \
} while(0);    


#define FORCE_MORSE_2POP                    \
do {                                        \
    float e = FORCE_MORSE_EPSILON;          \
    float a = FORCE_MORSE_ALPHA;            \
    float r0 = trad+rad[tn];                \
    float rcut = CONST_CUTOFF_FACTOR*r0;    \
    float fex  = exp(-a*(rcut-r0));         \
    float fcut = 2*e*a * (1-fex)*fex;       \
                                            \
    float l = sqrt(dist);                   \
    float ex = exp(-a*(l-r0));              \
    float co = (2*e*a * (1-ex)*ex-fcut)     \
  *(l<rcut)*(ttype==type[tn]?               \
        FORCE_MORSE_RELATIVE:1.0);          \
                                            \
    fx += co * dx[0]/l;                     \
    fy += co * dx[1]/l;                     \
} while(0);    


    
