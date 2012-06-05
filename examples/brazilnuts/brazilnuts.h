#define DIM 2
#define CONST_PBC               {1,0}
#define CONST_PARTICLECOUNT     512*4
#define CONST_KICKFORCE         20.0

#define FUNCTION_INIT           INIT_BRAZILNUTS

#define FORCE_HERTZ_EPSILON     150.0
#define FUNCTION_FORCE_PAIR     FORCE_HERTZ

#define FORCE_DAMPING_COEFF     0.1 
#define FORCE_DAMPING_SPEED     0.0
#define FORCE_GRAVITY_G         1.0
#define FUNCTION_FORCE_GLOBAL   {FORCE_DAMPING FORCE_GRAVITY FORCE_KICK}

