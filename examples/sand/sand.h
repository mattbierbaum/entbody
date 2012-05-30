#define CONST_PBC               {1,0}
#define CONST_PARTICLECOUNT     512*1
#define CONST_KICKFORCE         100.0

#define ARGV1           FORCE_HERTZ_EPSILON 
#define ARGV1_TYPE      double 
#define ARGV1_DEFAULT   200.0

#define INIT_RANDOM_FRACTION    2.1
#define FUNCTION_INIT           {INIT_RANDOM rad[0] = 3; type[0] = RED;}

#define FORCE_HERTZ_EPSILON     ARGV1_VAR 
#define FUNCTION_FORCE_PAIR     FORCE_HERTZ

#define FORCE_DAMPING_COEFF     0.1 
#define FORCE_DAMPING_SPEED     0.0
#define FORCE_GRAVITY_G         1.0
#define FUNCTION_FORCE_GLOBAL   {FORCE_DAMPING FORCE_GRAVITY FORCE_KICK}

