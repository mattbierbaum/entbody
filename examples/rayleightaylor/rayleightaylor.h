#define CONST_PBC               {1,1}
#define CONST_PARTICLECOUNT     512*6
#define CONST_KICKFORCE         1.0
#define CONST_COLOR_FACTOR      1.0
#define CONST_FRAME_SKIP        2

#define INIT_RAYLEIGHTAYLOR_FRACTION 0.3
#define FUNCTION_INIT           {INIT_RAYLEIGHTAYLOR; Tglobal=0.0;}

#define CONST_CUTOFF_FACTOR     4.5
#define FORCE_MORSE_EPSILON     1000.0
#define FORCE_MORSE_ALPHA       0.10
#define FORCE_MORSE_RELATIVE    0.7
#define FUNCTION_FORCE_PAIR     FORCE_MORSE_2POP

#define FORCE_DAMPING_COEFF     0.5
#define FORCE_DAMPING_SPEED     0.0
#define FORCE_GRAVITY_G         0.1
#define FUNCTION_FORCE_GLOBAL               \
    {                                       \
        FORCE_DAMPING FORCE_GRAVITY         \
        FORCE_KICK FORCE_THERMAL            \
        if (ttype == RED && t > 10.0f)      \
            fx -= 0.0;                     \
    }
