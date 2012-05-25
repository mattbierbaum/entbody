#define CONST_PBC               {1,1}
#define CONST_PARTICLECOUNT     512*6
#define CONST_KICKFORCE         1.0
#define CONST_COLOR_FACTOR      1.0
#define CONST_FRAME_SKIP        2

#define INIT_CUBIC_LATTICE_A    2.0 
#define FUNCTION_INIT                           \
        INIT_CUBIC_LATTICE                      \
        for (i=0.3*N; i<0.7*N; i++) type[i] = RED;

#define CONST_CUTOFF_FACTOR     4.5 /*5.5 is really strange*/
#define FORCE_MORSE_EPSILON     2000.0
#define FORCE_MORSE_ALPHA       0.12
#define FORCE_MORSE_RELATIVE    0.7
#define FUNCTION_FORCE_PAIR     FORCE_MORSE_2POP

#define FORCE_DAMPING_COEFF     0.4
#define FORCE_DAMPING_SPEED     0.0
#define FORCE_GRAVITY_G         0.0
#define FUNCTION_FORCE_GLOBAL   {FORCE_DAMPING FORCE_GRAVITY \
                                 FORCE_KICK FORCE_THERMAL } 
