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
    double e = 830.0;                       \
    double a = 0.15;                        \
    double r0 = rad[i]+rad[n];              \
    double rcut = 1.5*r0;                   \
    double fex  = exp(-a*(rcut-r0));        \
    double fcut = 2*e*a * (1-fex)*fex;      \
                                            \
    double l = sqrt(dist);                  \
    double ex = exp(-a*(l-r0));             \
    double co = (2*e*a * (1-ex)*ex-fcut)    \
  *(l<rcut)*(type[i]==type[n]?0.1:1.0);     \
                                            \
    f[2*i+0] += co * dx[0]/l;               \
    f[2*i+1] += co * dx[1]/l;               \
} while(0);    


    
