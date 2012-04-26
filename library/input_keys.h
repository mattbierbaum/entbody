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
            printf("temp\n");\
        }


#define INPUT_KEYS_WASD                     \
        if (key['w'] == 1){                 \
            for (i=0; i<N; i++){            \
                if (type[i] == RED)         \
                    o[2*i+1] = -kickforce;  \
            }                               \
        }                                   \
        if (key['s'] == 1){                 \
            for (i=0; i<N; i++){            \
                if (type[i] == RED)         \
                    o[2*i+1] = kickforce;   \
            }                               \
        }                                   \
        if (key['a'] == 1){                 \
            for (i=0; i<N; i++){            \
                if (type[i] == RED)         \
                    o[2*i+0] = -kickforce;  \
            }                               \
        }                                   \
        if (key['d'] == 1){                 \
            for (i=0; i<N; i++){            \
                if (type[i] == RED)         \
                    o[2*i+0] = kickforce;   \
            }                               \
        }                                   
 
