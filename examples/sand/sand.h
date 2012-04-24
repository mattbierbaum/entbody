#define CONST_PBC               {1,0}
#define CONST_PARTICLECOUNT     512*4

#define INIT_RANDOM \
do {                                \
    long i;                         \
    double radius  = 1.0;           \
    L = sqrt(5*pi*radius*radius*N); \
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


#define FORCE_HERTZ_EPSILON     150.0
#define FUNCTION_FORCE_PAIR     FORCE_HERTZ

#define FORCE_DAMPING_COEFF     0.3 
#define FUNCTION_FORCE_GLOBAL   {FORCE_DAMPING FORCE_GRAVITY FORCE_KICK}

