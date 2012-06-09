#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifdef CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <driver_types.h>
#endif

#ifdef HEADER
#include HEADER
#endif
#include "library/defaults.h"

// the order here matters!
#include "core/util.h"
#include "core/neighborlist.h"
#include "library/argv.h"
#include "library/cuda.h"
#include "library/init_regular.h"
#include "library/init_special.h"
#include "library/forces_pair.h"
#include "library/forces_global.h"
#include "library/input_keys.h"

#ifdef PLOT
#include "addons/plot.h"
#endif

#ifdef FPS
#include <time.h>
#endif

void simulate(int s);

#ifdef ARGV1
#define ARGV1_VAR       ARGVNAME(var,ARGV1)
#define ARGV1_CONVERTER ARGV_CONVERTER(ARGV1_TYPE)
ARGV1_TYPE ARGV1_VAR  = ARGV1_DEFAULT;
#endif

#ifdef ARGV2
#define ARGV2_VAR       ARGVNAME(var,ARGV2)
#define ARGV2_CONVERTER ARGV_CONVERTER(ARGV2_TYPE)
ARGV2_TYPE ARGV2_VAR  = ARGV2_DEFAULT;
#endif

#ifdef CUDA
CUDA_DEVICE_FUNCS
#endif
//===================================================
// the main function
//===================================================
int main(int argc, char **argv){
    int seed_in = 0;

    #ifdef CUDA
    CUT_DEVICE_INIT();
    #endif
    
    #ifdef ARGV1
    if (argc > 1)
        ARGV1_VAR = ARGV1_CONVERTER(argv[1]);
    #endif
    #ifdef ARGV2
    if (argc > 2)
        ARGV2_VAR = ARGV2_CONVERTER(argv[2]);
    #endif

    if (argc == 1) 
        simulate(seed_in);
    else if (argc == 2)
        simulate(seed_in);
    else {
        printf("usage:\n");
        printf("\t./entbody [seed]\n");
    }
    return 0;
}

//==================================================================
// the timestep - can be CPU or CUDA!
//==================================================================
CUDA_GLOBAL
void step(float *x, float *v, int *type, float *rad, float *col, int *key, NBL_STRUCT *nsc,
          long N, float L, float R, int *pbc, float dt, float t, float Tglobal, float colfact){

    int i, j;
    float dist;

    float px, py;
    float vx, vy;
    float fx, fy;
    float ox, oy;
 
    int ttype;
    float trad;
    float tcol;
    float R2 = R*R;
    
    FUNCTION_OBJECTS_CREATE
    #ifdef NEIGHBOR_DEBUG
    for (i=0; i<N; i++) col[i] = 10.0;
    #endif
    for (i=0; i<N; i++){
        tcol  = col[i]; ttype = type[i]; trad = rad[i];
        px = x[2*i+0];  py = x[2*i+1];
        vx = v[2*i+0];  vy = v[2*i+1]; 

        fx = 0.0;       fy = 0.0;
        ox = 0.0;       oy = 0.0; 
 
        #ifdef PLOT
        INPUT_KEYS_TEMP
        INPUT_KEYS_WASD
        #endif

        //====================================
        // loop over neighbors
        NBL_NEIGHBORS

        for (j=0; j<neigh_count; j++){
            float *dx = &rij[j*DIM];
            unsigned int tn = neighs[j];
            dist = dx[0]*dx[0] + dx[1]*dx[1];
            FUNCTION_FORCE_PAIR

            #ifdef NEIGHBOR_DEBUG
            if (ttype == RED){
                col[tn] = 0.01;
            }
            #else
            //FIXME - this adds the first force multiple times!
            tcol += fx*fx + fy*fy;
            #endif
        }

        //=====================================
        // global forces    
        FUNCTION_FORCE_GLOBAL

        //=====================================
        // Newton-Stomer-Verlet
        vx += fx * dt;
        vy += fy * dt;

        px += vx * dt;
        py += vy * dt;
        
        //======================================
        // boundary conditions 
        const float restoration = 0.5;
        if (pbc[0] == 1){
            if (px >= L-EPSILON || px < 0)
                px = mymod(px, L);
        }
        else {
            if (px >= L){px = 2*L-px; vx *= -restoration;}
            if (px < 0) {px = -px;    vx *= -restoration;}
            if (px >= L-EPSILON || px < 0){px = mymod(px, L);}
        }

        if (pbc[1] == 1){
            if (py >= L-EPSILON || py < 0)
                py = mymod(py, L);
        }
        else {
            if (py >= L){py = 2*L-py; vy *= -restoration;}
            if (py < 0) {py = -py;    vy *= -restoration;}
            if (py >= L-EPSILON || py < 0){py = mymod(py, L);}
        }
 
        tcol = tcol/colfact; 

        col[i] = tcol;  type[i] = ttype;
        x[2*i+0] = px;  x[2*i+1] = py;
        v[2*i+0] = vx;  v[2*i+1] = vy; 
    }

    FUNCTION_OBJECTS_FREE
}


//==================================================
// simulation
//==================================================
void simulate(int seed){
    ran_seed(seed);

    int    N      = CONST_PARTICLECOUNT;
    int pbc[]     = CONST_PBC; 
    float L       = 0.0;
    float dt      = 0.1;
    float t       = 0.0;
    float Tglobal = 0.0;

    float colfact = CONST_COLOR_FACTOR;

    int i;
    int mem_size_f = sizeof(float)*N;
    int mem_size_i = sizeof(int)*N;
    int mem_size_k = sizeof(int)*256;

    int *type    =   (int*)malloc(mem_size_i);
    float *rad   = (float*)malloc(mem_size_f);
    float *col   = (float*)malloc(mem_size_f);
    for (i=0; i<N; i++){ type[i] = 0; rad[i] = col[i] = 0.0;}

    float *x     = (float*)malloc(2*mem_size_f);
    float *v     = (float*)malloc(2*mem_size_f);
    float *copyx = (float*)malloc(2*mem_size_f);
    for (i=0; i<2*N; i++){x[i] = v[i] = copyx[i] = 0.0;}

    float time_end = CONST_TIME_END;

    #ifdef PLOT 
        int *key;
        plot_init(680);
        plot_clear_screen();
        key = plot_render_particles(x, rad, type, N, L,col);
    #else
        int *key = (int*)malloc(mem_size_k);
        memset(key, 0, mem_size_k);
    #endif

    //==========================================
    // initialize
    FUNCTION_INIT

    // find out what happened in initialization
    float maxr = 0.0;
    for (i=0; i<N; i++)
        if (rad[i] > maxr) maxr = rad[i];
    float R = 2*maxr*CONST_CUTOFF_FACTOR;

    NBL_BUILD

    //==========================================================
    // where the magic happens
    //==========================================================
    #ifdef CUDA
    CUDA_MEMORY_CREATE
    #endif

    int frames = 0;

    #ifdef FPS
    struct timespec start;
    clock_gettime(CLOCK_REALTIME, &start);
    #endif

    const int blocks = 512;
    for (t=0.0; t<time_end; t+=dt){
        #ifndef CUDA
        NBL_RESET
        NBL_UPDATE
        step(x, v, type, rad, col, key, nsc,
             N, L, R, pbc, dt, t, Tglobal, colfact);
        #else
        cudaMemcpy(cu_key, key, mem_size_k, cudaMemcpyHostToDevice);
        step<<<blocks, N/blocks >>>(cu_x, cu_copyx, cu_v, cu_type, cu_rad, cu_col, 
                    cu_cells, cu_count, cu_size, size_total, cu_key,
                    N, L, R, cu_pbc, dt, t, Tglobal, colfact);
        ERROR_CHECK
        #endif

        #ifdef PLOT 
        if (frames % CONST_FRAME_SKIP == 0 && key['h'] != 1){
            plot_clear_screen();
            key = plot_render_particles(x, rad, type, N, L,col);
            INPUT_KEYS_QUIT
        }
        #endif
        frames++;
    }
    // end of the magic, cleanup
    //----------------------------------------------
    #ifdef FPS
    struct timespec end;
    clock_gettime(CLOCK_REALTIME, &end);
    printf("fps = %f\n", frames/((end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec)/1e9));
    #endif

    free(x);
    free(v);
    free(rad);
    free(type);
    free(col);

    NBL_FREE

    #ifdef CUDA
    CUDA_MEMORY_FREE
    #endif

    #ifdef PLOT
    plot_clean(); 
    #else
    free(key);
    #endif
}





