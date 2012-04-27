#define CONST_PBC               {1,1}
#define CONST_PARTICLECOUNT     40
#define CONST_COLOR_FACTOR      50000.0
#define CONST_CUTOFF_FACTOR     2.0
//#define CONST_TIME_END          1.0

#define INIT_RANDOM_FRACTION    2.0 
#define FUNCTION_INIT           INIT_RANDOM 

#define FORCE_HERTZ_EPSILON     140.0
#define FUNCTION_FORCE_PAIR     FORCE_HERTZ

#define CONST_KICKFORCE         80.0
#define FORCE_DAMPING_COEFF     0.9
#define FORCE_DAMPING_SPEED     0.0 
#define FUNCTION_FORCE_GLOBAL   {FORCE_DAMPING FORCE_KICK}

