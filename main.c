#include <stdio.h>
#include <stdlib.h>
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
#else
#include "library/defaults.h"
#endif

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
void step(float *x, float *v, int *type, float *rad, float *col, 
          int *cells, int *count, int *size, int size_total, 
          long N, float L, float R, int *pbc, float dt, float colfact){

    #ifndef CUDA
    int i,j;
    #else 
    int i = blockDim.x*blockIdx.x + threadIdx.x;
    #endif

    //=========================================
    // reset the neighborlists
    int index[2];
    #ifndef CUDA
    for (i=0; i<size_total*NMAX; i++){
    #endif
        if (i < size_total)
            count[i] = 0;
        if (i < size_total*NMAX)
            cells[i] = 0;
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
        #ifndef CUDA
        count[t]++;
        unsigned int pos = count[t];
        #else
        unsigned int pos = atomicInc(&count[t], 0xffffffff);
        #endif
        cells[NMAX*t + pos] = i;
    #ifndef CUDA
    }   
    #else
    __syncthreads();
    #endif

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
    float *tempx = (float*)malloc(sizeof(float)*N*2);
    float *tempv = (float*)malloc(sizeof(float)*N*2);
 
    for (i=0; i<N; i++){
    #endif
        tcol  = col[i]; ttype = type[i]; trad = rad[i];
        px = x[2*i+0];  py = x[2*i+1];
        vx = v[2*i+0];  vy = v[2*i+1]; 

        fx = 0.0;       fy = 0.0;
        wx = 0.0;       wy = 0.0;
        ox = 0.0;       oy = 0.0; 
 
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

                    dist = 0.0;
                    dx[0] = x[2*tn+0] - px;
                    if (image[0]) dx[0] += L*tt[0];
                    dist += dx[0]*dx[0];
                    
                    dx[1] = x[2*tn+1] - py;
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

        col[i] = tcol; type[i] = ttype;
        tempx[2*i+0] = px;  tempx[2*i+1] = py;
        tempv[2*i+0] = vx;  tempv[2*i+1] = vy; 
    #ifndef CUDA
    }
    #ifdef OPENMP
    #pragma omp barrier
    #endif
    for (i=0; i<N; i++){
        x[2*i+0] = tempx[2*i+0]; x[2*i+1] = tempx[2*i+1];
        v[2*i+0] = tempv[2*i+0]; v[2*i+1] = tempv[2*i+1];
    }
    free(tempx);
    free(tempv);
    #endif
}


//==================================================
// simulation
//==================================================
void simulate(int seed){
    ran_seed(seed);

    int    N       = CONST_PARTICLECOUNT;
    int pbc[]      = CONST_PBC; 
    float L       = 0.0;
    float dt      = 1e-1;
    float t       = 0.0;
    float Tglobal = 0.0;

    float colfact = 3.0;
    #ifdef CONST_COLOR_FACTOR
    colfact = CONST_COLOR_FACTOR;
    #endif

    int i;
    int *type   = (int*)malloc(sizeof(int)*N);
    float *rad = (float*)malloc(sizeof(float)*N); 
    float *col = (float*)malloc(sizeof(float)*N); 
    for (i=0; i<N; i++){ type[i] = rad[i] = 0;}

    float *x = (float*)malloc(sizeof(float)*2*N);
    float *v = (float*)malloc(sizeof(float)*2*N);
    float *f = (float*)malloc(sizeof(float)*2*N);
    float *w = (float*)malloc(sizeof(float)*2*N);
    float *o = (float*)malloc(sizeof(float)*2*N);
    for (i=0; i<2*N; i++){o[i] = x[i] = v[i] = f[i] = w[i] = 0.0;}

    #ifdef PLOT 
    float time_end = 1e20;
    #else
    float time_end = 1e2;
    #endif

    #ifdef PLOT 
        int *key;
        float kickforce = 1.0;
        #ifdef CONST_KICKFORCE
        kickforce = CONST_KICKFORCE;
        #endif

        //#ifdef POINTS
        //    plot_init((int)(0.92*sqrt(N)));
        //#else
            plot_init(680);
        //#endif 
        plot_clear_screen();
        key = plot_render_particles(x, rad, type, N, L,col);
    #endif

    //==========================================
    // initialize
    FUNCTION_INIT

    // find out what happened in initialization
    float maxr = 0.0;
    for (i=0; i<N; i++)
        if (rad[i] > maxr) maxr = rad[i];
    float R = 2*maxr;
    #ifdef CONST_CUTOFF_FACTOR
    R *= CONST_CUTOFF_FACTOR;
    #endif

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

    #ifdef FUNCTION_OBJECTS_CREATE
    FUNCTION_OBJECTS_CREATE
    #endif

    for (t=0.0; t<time_end; t+=dt){
        step(x, v, type, rad, col, 
            cells, count, size, size_total, 
            N, L, R, pbc, dt, colfact);

        #ifdef PLOT 
            plot_clear_screen();
            key = plot_render_particles(x, rad, type, N, L,col);
        #endif
        frames++;

        #ifdef PLOT
        INPUT_KEYS_QUIT
        INPUT_KEYS_TEMP
        INPUT_KEYS_WASD
        #endif
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
 
    free(x);
    free(v);
    free(f);
    free(w);
    free(o);
    free(rad);
    free(type);
    free(col);

    #ifdef FUNCTION_OBJECTS_FREE
    FUNCTION_OBJECTS_FREE
    #endif

    #ifdef PLOT
    plot_clean(); 
    #endif
}





