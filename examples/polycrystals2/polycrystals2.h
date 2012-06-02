#define CONST_PBC               {1,1}
#define CONST_PARTICLECOUNT     512*8
#define CONST_KICK_FORCE        0.1f
#define CONST_FRAME_SKIP        5
#define CONST_COLOR_FACTOR      50.0f

#define INIT_RANDOM_FRACTION    0.97 
#define FUNCTION_INIT           INIT_RANDOM

#define FORCE_HERTZ_EPSILON     80.0f
#define FUNCTION_FORCE_PAIR     FORCE_HERTZ

#define FORCE_DAMPING_COEFF     0.5f
#define FORCE_DAMPING_SPEED     0.0f
#define FUNCTION_FORCE_GLOBAL   {FORCE_DAMPING FORCE_KICK}

