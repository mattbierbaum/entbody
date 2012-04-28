#include "../util.h"

#define INPUT_KEYS_QUIT     \
        if (key['q'] == 1)  \
            break;

#define INPUT_KEYS_TEMP         \
        if (key['9'] == 1)      \
            Tglobal -= 0.01;    \
        if (key['0'] == 1)      \
            Tglobal += 0.01;    \
        if (key['8'] == 1){      \
            Tglobal = 0.0;      \
        }

#define KICKFORCE 2.0
#define INPUT_KEYS_WASD                     \
        if (key['w'] == 1){                 \
            if (ttype == RED)               \
                oy = -KICKFORCE;            \
        }                                   \
        if (key['s'] == 1){                 \
            if (ttype == RED)               \
                oy = KICKFORCE;             \
        }                                   \
        if (key['a'] == 1){                 \
            if (ttype == RED)               \
                ox = -KICKFORCE;            \
        }                                   \
        if (key['d'] == 1){                 \
            if (ttype == RED)               \
                ox = KICKFORCE;             \
        }                                   
 
