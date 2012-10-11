#define CONST_PBC               {1,1}
#define CONST_PARTICLECOUNT     512*1
#define CONST_KICKFORCE         20.0
#define CONST_CUTOFF_FACTOR     2 
#define CONST_COLOR_FACTOR      20
#define CONST_FRAME_SKIP        4

#define FUNCTION_INIT           INIT_CENTRAL_CIRCLE

#define FUNCTION_OBJECTS_CREATE                 \
    float wlen=0.0, vlen=0.0, vhappy=0.0, wx=0.0f, wy=0.0f; 
#define FUNCTION_OBJECTS_FREE  

#define FLOCK_STRENGTH          1.0
#define FORCE_FLOCK_ADD                     \
    if (dist > EPSILON &&                   \
        dist < CONST_CUTOFF_FACTOR &&       \
        type[j] == RED && ttype == RED){    \
        wx += vx; wy += vy;                 \
    }

#define FORCE_FLOCK_SUM                     \
    wlen = wx*wx + wy*wy;                   \
    if (ttype == RED && wlen > EPSILON){    \
        fx += FLOCK_STRENGTH * wx / wlen;   \
        fy += FLOCK_STRENGTH * wy / wlen;   \
        wx = 0.0f; wy = 0.0f;               \
    }

#define VHAPPY_RED              1.0
#define VHAPPY_BLACK            0.0
#define DAMP_COEFF              1.0
#define FORCE_PROPULSION                            \
    vlen = vx*vx + vy*vy;                           \
    vhappy = ttype==RED?VHAPPY_RED:VHAPPY_BLACK;    \
    if (vlen > EPSILON){                            \
        fx += DAMP_COEFF*(vhappy - vlen)*vx/vlen;   \
        fy += DAMP_COEFF*(vhappy - vlen)*vy/vlen;   \
    }


#define FORCE_HERTZ_EPSILON     100.0
#define FUNCTION_FORCE_PAIR     {FORCE_HERTZ FORCE_FLOCK_ADD}
#define FUNCTION_FORCE_GLOBAL   {FORCE_KICK FORCE_PROPULSION FORCE_FLOCK_SUM}

