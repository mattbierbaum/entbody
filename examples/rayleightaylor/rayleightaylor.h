#define CONST_PBC               {1,0}
#define CONST_PARTICLECOUNT     512*16
#define CONST_KICKFORCE         1.0
#define CONST_COLOR_FACTOR      1.0

#define INIT_RAYLEIGHTAYLOR_FRACTION 0.3
#define FUNCTION_INIT           {INIT_RAYLEIGHTAYLOR; Tglobal=0.0;}

#define CONST_CUTOFF_FACTOR     1.5
#define FORCE_MORSE_EPSILON     830.0
#define FORCE_MORSE_ALPHA       0.10
#define FORCE_MORSE_RELATIVE    0.1
#define FUNCTION_FORCE_PAIR     FORCE_MORSE_2POP

#define FORCE_DAMPING_COEFF     0.6
#define FORCE_DAMPING_SPEED     0.0
#define FORCE_GRAVITY_G         0.1
#define FUNCTION_FORCE_GLOBAL   {FORCE_DAMPING FORCE_GRAVITY FORCE_KICK FORCE_THERMAL col[i]=0.7;}/* if (Tglobal < 20.0) Tglobal+=0.01/N;}*/
