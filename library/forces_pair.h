#include <math.h>
#include "../util.h"

#define FORCE_HERTZ \
do {                                                  \
    double r0 = rad[i]+rad[n];                        \
    double l = sqrt(dist);                            \
    double co = FORCE_HERTZ_EPSILON                   \
            * (1-l/r0)*(1-l/r0) * (l<r0);             \
    f[2*i+0] += -co * dx[0];                          \
    f[2*i+1] += -co * dx[1];                          \
} while (0); 


#define FORCE_MORSE \
do {                                        \
    double e = FORCE_MORSE_EPSILON;         \
    double a = FORCE_MORSE_ALPHA;           \
    double r0 = rad[i]+rad[n];              \
    double rcut = CONST_CUTOFF_FACTOR*r0;   \
    double fex  = exp(-a*(rcut-r0));        \
    double fcut = 2*e*a * (1-fex)*fex;      \
                                            \
    double l = sqrt(dist);                  \
    double ex = exp(-a*(l-r0));             \
    double co = (2*e*a * (1-ex)*ex-fcut)    \
  *(l<rcut);                                \
                                            \
    f[2*i+0] += co * dx[0]/l;               \
    f[2*i+1] += co * dx[1]/l;               \
} while(0);    


#define FORCE_MORSE_2POP \
do {                                        \
    double e = FORCE_MORSE_EPSILON;         \
    double a = FORCE_MORSE_ALPHA;           \
    double r0 = rad[i]+rad[n];              \
    double rcut = CONST_CUTOFF_FACTOR*r0;   \
    double fex  = exp(-a*(rcut-r0));        \
    double fcut = 2*e*a * (1-fex)*fex;      \
                                            \
    double l = sqrt(dist);                  \
    double ex = exp(-a*(l-r0));             \
    double co = (2*e*a * (1-ex)*ex-fcut)    \
  *(l<rcut)*(type[i]==type[n]?FORCE_MORSE_RELATIVE:1.0);     \
                                            \
    f[2*i+0] += co * dx[0]/l;               \
    f[2*i+1] += co * dx[1]/l;               \
} while(0);    


    
