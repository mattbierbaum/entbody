#include "../core/util.h"

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

//#define KICKFORCE 16.0
#ifndef FORCE_KICK_CONST
#define FORCE_KICK_CONST 16.0f
#endif

#define INPUT_KEYS_WASD                     \
        if (key['w'] == 1){                 \
            if (ttype == RED)               \
                oy = -FORCE_KICK_CONST;     \
        }                                   \
        if (key['s'] == 1){                 \
            if (ttype == RED)               \
                oy = FORCE_KICK_CONST;      \
        }                                   \
        if (key['a'] == 1){                 \
            if (ttype == RED)               \
                ox = -FORCE_KICK_CONST;     \
        }                                   \
        if (key['d'] == 1){                 \
            if (ttype == RED)               \
                ox = FORCE_KICK_CONST;      \
        }                                   
 
