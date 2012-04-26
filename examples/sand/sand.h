#define CONST_PBC               {1,0}
#define CONST_PARTICLECOUNT     512*4
#define CONST_KICKFORCE         10.0

#define ARGV1           FORCE_HERTZ_EPSILON 
#define ARGV1_TYPE      double 
#define ARGV1_DEFAULT   150.0

#define FUNCTION_INIT \
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
            rad[i] = 5*radius;      \
        }                           \
     }                              \
} while(0);


#define FORCE_HERTZ_EPSILON     ARGV1_VAR 
#define FUNCTION_FORCE_PAIR     FORCE_HERTZ

#define FORCE_DAMPING_COEFF     0.1 
#define FORCE_DAMPING_SPEED     0.0
#define FORCE_GRAVITY_G         1.0
#define FUNCTION_FORCE_GLOBAL   {FORCE_DAMPING FORCE_GRAVITY FORCE_KICK}

