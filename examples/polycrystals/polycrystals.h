#define CONST_PBC               {1,1}
#define CONST_PARTICLECOUNT     512*512
#define CONST_KICKFORCE         10.0f
#define CONST_TIME_END          1.0f

#define INIT_RANDOM_FRACTION    1.0f
#define FUNCTION_INIT           INIT_RANDOM

#define CONST_COLOR_FACTOR      20.0f
#define FORCE_HERTZ_EPSILON     150.0f
#define FUNCTION_FORCE_PAIR     FORCE_HERTZ

#define FORCE_DAMPING_COEFF     1.0f
#define FORCE_DAMPING_SPEED     0.0f
#define FUNCTION_FORCE_GLOBAL   {FORCE_DAMPING FORCE_KICK}

