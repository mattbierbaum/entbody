#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

// the order here matters!
#include "util.h"
#include "library/argv.h"
#include "library/cuda.h"
#include "library/init_regular.h"
#include "library/init_special.h"
#include "library/forces_pair.h"
#include "library/forces_global.h"
#include "library/input_keys.h"

#ifdef HEADER
#include HEADER
#endif
#include "library/defaults.h"

#ifdef PLOT
#include "plot.h"
#endif

#ifdef FPS
#include <time.h>
#endif

void simulate(int s);

#ifdef ARGV1
#define ARGV1_VAR       ARGVNAME(var,ARGV1)
#define ARGV1_CONVERTER ARGV_CONVERTER(ARGV1_TYPE)
ARGV1_TYPE ARGV1_VAR = ARGV1_DEFAULT;
#endif

#ifdef ARGV2
#define ARGV2_VAR       ARGVNAME(var,ARGV2)
#define ARGV2_CONVERTER ARGV_CONVERTER(ARGV2_TYPE)
ARGV2_TYPE ARGV2_VAR = ARGV2_DEFAULT;
#endif

//===================================================
// the main function
//===================================================
int main(int argc, char **argv){
    int seed_in = 0;

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
void step(float *x, float *copyx, float *v, int *type, float *rad, float *col, 
          int *cells, int *count, int *size, int size_total, int *key,
          long N, float L, float R, int *pbc, float dt, float Tglobal, float colfact){

    #ifndef CUDA
    int i,j;
    #else 
    int i = blockDim.x*blockIdx.x + threadIdx.x;
    #endif

    //=========================================
    // reset the neighborlists
    int index[2];
    #ifndef CUDA
    for (i=0; i<size_total; i++){
    #endif
        if (i < size_total)
            count[i] = 0;
    #ifndef CUDA
    }
    #else
    __syncthreads();
    #endif

    //=========================================
    // rehash all of the particles into the list
    #ifndef CUDA
    for (i=0; i<N; i++){
    #endif
        index[0] = (int)(x[2*i+0]/L  * size[0]);
        index[1] = (int)(x[2*i+1]/L  * size[1]);
        int t = index[0] + index[1]*size[0];
        cells[NMAX*t + count[t]] = i;
        copyx[2*i+ 0] = x[2*i+0];
        copyx[2*i+ 1] = x[2*i+1];
    
        #ifndef CUDA
        count[t] = count[t]+1;
        #else
        atomicInc(&count[t], 0xffffffff);
        #endif
     #ifndef CUDA
    }   
    #else
    __syncthreads();
    #endif

    //==========================================
    // this is mainly for CUDA optimization
    int tt[2];
    int tix[2];
    int image[2];
    float dx[2];
    int goodcell, ind, tn;
    float dist;

    float px, py;
    float vx, vy;
    float fx, fy;
    float wx, wy;
    float ox, oy;
 
    int ttype;
    float trad;
    float tcol;
    float R2 = R*R;

    //==========================================
    // find forces on all particles
    #ifndef CUDA
    #ifdef OPENMP
    #pragma omp parallel for private(i,dx,index,tt,goodcell,tix,ind,j,tn,image,k,dist)
    #endif
    for (i=0; i<N; i++){
    #endif
        tcol  = col[i]; ttype = type[i]; trad = rad[i];
        px = x[2*i+0];  py = x[2*i+1];
        vx = v[2*i+0];  vy = v[2*i+1]; 

        fx = 0.0;       fy = 0.0;
        wx = 0.0;       wy = 0.0;
        ox = 0.0;       oy = 0.0; 
 
        #ifdef PLOT
        INPUT_KEYS_TEMP
        INPUT_KEYS_WASD
        #endif

        index[0] = (int)(px/L * size[0]);
        index[1] = (int)(py/L * size[1]);

        for (tt[0]=-1; tt[0]<=1; tt[0]++){
        for (tt[1]=-1; tt[1]<=1; tt[1]++){
            goodcell = 1;
            tix[0] = mod_rvec(index[0]+tt[0],size[0]-1,pbc[0],&image[0]);
            tix[1] = mod_rvec(index[1]+tt[1],size[1]-1,pbc[1],&image[1]);
            if (pbc[0] < image[0] || pbc[1] < image[1])  goodcell = 0;

            if (goodcell){
                ind = tix[0] + tix[1]*size[0]; 
                for (j=0; j<count[ind]; j++){
                    tn = cells[NMAX*ind+j];
                    float px2 = copyx[2*tn+0];
                    float py2 = copyx[2*tn+1];

                    dist = 0.0;
                    dx[0] = px2 - px;
                    if (image[0]) dx[0] += L*tt[0];
                    dist += dx[0]*dx[0];
                    
                    dx[1] = py2 - py;
                    if (image[1]) dx[1] += L*tt[1];
                    dist += dx[1]*dx[1];

                    //===============================================
                    // force calculation 
                    if (dist > 1e-6 && dist < R2){
                        FUNCTION_FORCE_PAIR
                        tcol += fx*fx + fy*fy;
                    }
                 }
            }
        } } 

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
    #ifndef CUDA
    }
    #ifdef OPENMP
    #pragma omp barrier
    #endif
    #endif
}


//==================================================
// simulation
//==================================================
void simulate(int seed){
    ran_seed(seed);

    int    N      = CONST_PARTICLECOUNT;
    int pbc[]     = CONST_PBC; 
    float L       = 0.0;
    float dt      = 1e-1;
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
    for (i=0; i<N; i++){ type[i] = rad[i] = col[i] = 0;}

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

    // make boxes for the neighborlist
    int size[2];
    int size_total = 1;
    for (i=0; i<2; i++){
        size[i] = (int)(L / R); 
        size_total *= size[i];
    }

    int *count = (int*)malloc(sizeof(int)*size_total);
    int *cells = (int*)malloc(sizeof(int)*size_total*NMAX);
    for (i=0; i<size_total; i++)
        count[i] = 0;
    for (i=0; i<size_total*NMAX; i++)
        cells[i] = 0;

    //==========================================================
    // where the magic happens
    //==========================================================
    int frames = 0;

    #ifdef FPS
    struct timespec start;
    clock_gettime(CLOCK_REALTIME, &start);
    #endif

    FUNCTION_OBJECTS_CREATE

    for (t=0.0; t<time_end; t+=dt){
        #ifndef CUDA
        memcpy(copyx, x, 2*mem_size_f);
        step(x, copyx, v, type, rad, col, 
            cells, count, size, size_total, key,
            N, L, R, pbc, dt, Tglobal, colfact);
        #else
        cudaMemcpy(copyx, x, 2*mem_size_f);
        cudaMemcpy(key, cu_key, mem_size_k);
        step<<< >>>(cu_x, cu_copyx, cu_v, cu_type, cu_rad, cu_col, 
                    cu_cells, cu_count, cu_size, size_total, cu_key,
                    N, L, R, cu_pbc, dt, Tglobal, colfact);
        #endif

        #ifdef PLOT 
            plot_clear_screen();
            key = plot_render_particles(x, rad, type, N, L,col);
            INPUT_KEYS_QUIT
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

    free(cells);
    free(count);

    free(copyx); 
    free(x);
    free(v);
    free(rad);
    free(type);
    free(col);

    FUNCTION_OBJECTS_FREE

    #ifdef PLOT
    plot_clean(); 
    #else
    free(key);
    #endif
}





