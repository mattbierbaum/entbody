#define CONST_PBC               {0,0}
#define CONST_PARTICLECOUNT     10*10 
#define CONST_KICKFORCE         1.0f
#define CONST_FRAME_SKIP        2
#define CONST_COLOR_FACTOR      5.0f

#define INIT_CUBIC_LATTICE_A    1.5
#define FUNCTION_INIT           INIT_CUBIC_LATTICE

#define CONST_CUTOFF_FACTOR     2
#define FORCE_SPRING_CONST      10.0f
#define FUNCTION_FORCE_PAIR     FORCE_SPRING

#define FORCE_DAMPING_COEFF     1.0f
#define FORCE_DAMPING_SPEED     0.0f
#define FUNCTION_FORCE_GLOBAL   {FORCE_DAMPING FORCE_KICK col[i]=0.4;}

