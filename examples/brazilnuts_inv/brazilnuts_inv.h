#define CONST_PBC               {0,0}
#define CONST_PARTICLECOUNT     512*2
#define CONST_KICKFORCE         10.0f
#define CONST_FRAME_SKIP        10
#define CONST_COLOR_FACTOR      500.0f

#define INIT_RANDOM_FRACTION    4.5f
#define FUNCTION_INIT           {INIT_RANDOM rad[0]=2.0f; for (i=0; i<N; i++) {if (i%2==1) rad[i] = 1.1f;}}

#define CONST_CUTOFF_FACTOR     5.0
#define FORCE_MORSE_EPSILON     2000.0
#define FORCE_MORSE_ALPHA       0.10
#define FUNCTION_FORCE_PAIR     FORCE_MORSE

#define FORCE_DAMPING_COEFF     0.1f
#define FORCE_DAMPING_SPEED     0.0f
#define FORCE_GRAVITY_G         1.0f
#define FUNCTION_FORCE_GLOBAL   {FORCE_DAMPING FORCE_KICK FORCE_GRAVITY}

