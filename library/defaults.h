#ifndef CONST_PBC
#define CONST_PBC               {1,1}
#endif

#ifndef CONST_PARTICLECOUNT
#define CONST_PARTICLECOUNT     512
#endif

#ifndef CONST_COLOR_FACTOR
#define CONST_COLOR_FACTOR      1.0
#endif

#ifndef CONST_TIME_END         
#define CONST_TIME_END          1e20
#endif

#ifndef CONST_FRAME_SKIP        
#define CONST_FRAME_SKIP        1  
#endif

#ifndef CONST_CUTOFF_FACTOR 
#define CONST_CUTOFF_FACTOR     1.0
#endif

#ifndef FUNCTION_INIT
#define INIT_RANDOM_FRACTION    0.2
#define FUNCTION_INIT           INIT_RANDOM
#endif

#ifndef FUNCTION_FORCE_PAIR 
#define FORCE_HERTZ_EPSILON     150.0
#define FUNCTION_FORCE_PAIR     FORCE_HERTZ
#endif

#ifndef FUNCTION_FORCE_GLOBAL
#define FORCE_DAMPING_COEFF     1.0
#define FUNCTION_FORCE_GLOBAL   {FORCE_DAMPING FORCE_KICKFORCE}
#endif

#ifndef FUNCTION_OBJECTS_CREATE
#define FUNCTION_OBJECTS_CREATE
#endif

#ifndef FUNCTION_OBJECTS_FREE
#define FUNCTION_OBJECTS_FREE
#endif
