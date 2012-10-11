#define CONST_PBC               {1,1}
#define CONST_PARTICLECOUNT     512*512
#define CONST_NO_NEIGHBORS      

#define INIT_DENSE_BOX                          \
    do {                                        \
        float radius = 1.0;                     \
        L = 100.0f;                             \
        float ep = EPSILON;                     \
        float sz = 1.0f/sqrt(N);                \
        float c  = L/2 - 1.0f/2;                \
        int side = (int)1.0f/sz;                \
                                                \
        for (i=0; i<N; i++){                    \
            rad[i] = 1./5;radius;               \
            x[2*i+0] = sz * (int)(i%side)+c;    \
            x[2*i+1] = sz * (int)(i/side)+c;    \
            v[2*i+0] = v[2*i+1] = 0.0f;         \
            type[i] = BLACK;                    \
        }                                       \
        type[0] = RED;                          \
    } while (0);   

#define INIT_DENSE_RANDOM                       \
    do {                                        \
        float radius = 1.0;                     \
        L = 100.0f;                             \
        float ep = EPSILON;                     \
        float sz = 1.0f/sqrt(N);                \
        float c  = L/2 - 1.0f/2;                \
        int side = (int)1.0f/sz;                \
                                                \
        for (i=0; i<N; i++){                    \
            rad[i] = 1./5;radius;               \
            x[2*i+0] = 2.0*ran_ran2() + L/2;    \
            x[2*i+1] = 2.0*ran_ran2() + L/2;    \
            v[2*i+0] = v[2*i+1] = 0.0f;         \
            type[i] = BLACK;                    \
        }                                       \
        type[0] = RED;                          \
    } while (0);   

#define FUNCTION_INIT           INIT_DENSE_RANDOM

#define FORCE_SOFTENING         0.00001f
#define SPEED_FACTOR            1
#define SWITCHING_TIME          (4.0f/SPEED_FACTOR)
#define FUNCTION_FORCE_GLOBAL                       \
    {                                               \
        int time_seg = (int)(t/SWITCHING_TIME);     \
        float sinkx = time_seg%2?25.0f:75.0f;       \
        float rr = time_seg%2?1.0f:-1.0f;           \
        float ff = SPEED_FACTOR*10.0f;              \
                                                    \
        float sinky = 50.0f;                        \
        float dx = px - sinkx;                      \
        float dy = py - sinky;                      \
        float r = sqrt(dx*dx + dy*dy);              \
        if (i != 0){                                \
            vx = -ff*dx/(r*r + FORCE_SOFTENING);    \
            vy = -ff*dy/(r*r + FORCE_SOFTENING);    \
            vx += - ff*dy/r;                        \
            vy += + ff*dx/r;                        \
        } else {                                    \
            px = sinkx;                             \
            py = sinky;                             \
        }                                           \
    }  
