#define CONST_PBC               {1,1}
#define CONST_PARTICLECOUNT     512*256
#define CONST_CUTOFF_FACTOR     2
#define CONST_FRAME_SKIP        10

#define INIT_CUBIC_LATTICE_A    1
#define FUNCTION_INIT           {INIT_CUBIC_LATTICE type[34] = RED;}

#define FUNCTION_OBJECTS_CREATE int nz = 0; int na = 0;
#define FUNCTION_FORCE_PAIR     { if (type[j] == RED) nz++; if (type[j] == WHITE) na++; }

#define WHITE 10
#define ZZB 0.8
#define ZZK 0.6
#define FUNCTION_FORCE_GLOBAL   { if (ran_ran2() < (1 - pow(1-ZZB, nz))) type[i] = RED;   \
                                  if (ran_ran2() < (1 - pow(1-ZZK, na))) type[i] = BLACK; \
                                }
