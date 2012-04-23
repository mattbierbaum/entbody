#include <math.h>
#include "util.h"
#include "forces.h"

inline void force_hertz(double *dx, double dist, double radi, double radj, int typei, int typej, double *f){
    double epsilon = 50.0;
    double r0 = radi+radj;
    double l = sqrt(dist);
    double co = epsilon * (1-l/r0)*(1-l/r0) * (l<r0);
    f[0] += -co * dx[0];
    f[1] += -co * dx[1];
}

inline void force_morse(double *dx, double dist, double radi, double radj, int typei, int typej, double *f){
    double e = 50.0;
    double a = 0.2; 
    double r0 = radi+radj;
    double l = sqrt(dist);
    double ex = exp(-a*(l-r0));
    double co = 2*e*a;
 
    f[0] += co * (1-ex)*ex * dx[0]/l;
    f[1] += co * (1-ex)*ex * dx[1]/l;
}

inline void force_damping(double *v, double *f){
    double damp = 1.0;
    double vlen = v[0]*v[0] + v[1]*v[1];
    if (vlen > 1e-6){
        f[0] -= damp*vlen*v[0]/vlen;
        f[1] -= damp*vlen*v[1]/vlen;
    }
}

inline void force_thermal(double T, double *f){
    f[0] += T*(ran_ran2()-0.5);
    f[1] += T*(ran_ran2()-0.5);
}


inline void force_kick(double *k, double *f){
    f[0] += k[0];
    f[1] += k[1];
}

inline void force_gravity(double *f){
    double g = 1.0;
    f[1] += g;
}
    
