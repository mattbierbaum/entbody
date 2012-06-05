#ifndef DIM 
#define DIM                     2
#endif

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
#define INIT_RANDOM_FRACTION    2.0
#define FUNCTION_INIT           INIT_RANDOM
#endif

#ifndef FUNCTION_FORCE_PAIR 
#define FORCE_HERTZ_EPSILON     150.0
#define FUNCTION_FORCE_PAIR     FORCE_HERTZ
#endif

#ifndef FUNCTION_FORCE_GLOBAL
#define FORCE_DAMPING_COEFF     1.0
#define FORCE_DAMPING_SPEED     0.0
#define CONST_KICKFORCE         1.0
#define FUNCTION_FORCE_GLOBAL   {FORCE_DAMPING FORCE_KICK}
#endif

#ifndef FUNCTION_OBJECTS_CREATE
#define FUNCTION_OBJECTS_CREATE
#endif

#ifndef FUNCTION_OBJECTS_FREE
#define FUNCTION_OBJECTS_FREE
#endif

#ifndef NBL_STRUCT
#define NBL_STRUCT nbl_struct_cell
#endif

#ifndef NBL_BUILD
#define NBL_BUILD NBL_STRUCT *nsc = nbl_cell_build(N, L, DIM, R, pbc, x);
#endif

#ifndef NBL_RESET
#define NBL_RESET nbl_cell_reset(N, nsc);
#endif

#ifndef NBL_UPDATE
#define NBL_UPDATE nbl_cell_update(N, L, DIM, R, pbc, x, nsc);
#endif

#ifndef NBL_NEIGHBORS
#define NBL_NEIGHBORS int *neighs; float *rij; int neigh_count = nbl_cell_neighbors(i, &neighs, &rij, DIM, nsc);
#endif

#ifndef NBL_FREE
#define NBL_FREE nbl_cell_free(nsc);
#endif
