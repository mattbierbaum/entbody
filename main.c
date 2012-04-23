#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "util.h"
#include "initialize.h"
#include "forces.h"

#ifdef PLOT
#include "plot.h"
#endif

#ifdef FPS
#include <time.h>
#endif

void simulate(int s);

//===================================================
// the main function
//===================================================
int main(int argc, char **argv){
    int seed_in     = 0;

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


//==================================================
// simulation
//==================================================
void simulate(int seed){
    ran_seed(seed);

    int    NMAX    = 50;
    int    N       = 512*8; 
    double L       = 0.0;

    int pbc[]      = {1,1};
    double dt      = 1e-1;
    double t       = 0.0;
    double Tglobal = 0.0;

    int i, j, k;

    int *type   = (int*)malloc(sizeof(int)*N);
    int *neigh  = (int*)malloc(sizeof(int)*N);
    double *rad = (double*)malloc(sizeof(double)*N); 
    double *col = (double*)malloc(sizeof(double)*N); 
    for (i=0; i<N; i++){ type[i] = neigh[i] = rad[i] = 0;}

    double *x = (double*)malloc(sizeof(double)*2*N);
    double *v = (double*)malloc(sizeof(double)*2*N);
    double *f = (double*)malloc(sizeof(double)*2*N);
    double *w = (double*)malloc(sizeof(double)*2*N);
    double *o = (double*)malloc(sizeof(double)*2*N);
    for (i=0; i<2*N; i++){o[i] = x[i] = v[i] = f[i] = w[i] = 0.0;}

    #ifdef PLOT 
    double time_end = 1e20;
    #else
    double time_end = 1e2;
    #endif

    #ifdef PLOT 
        int *key;
        double kickforce = 10.0;
        #ifdef POINTS
            plot_init((int)(0.92*sqrt(N)));
        #else
            plot_init(680);
        #endif 
        plot_clear_screen();
        key = plot_render_particles(x, rad, type, N, L,col);
    #endif

    //==========================================
    // initialize
    //init_brazilnuts(x, v, rad, type, &L, N);
    init_random(x, v, rad, type, &L, N);

    // find out what happened in initialization
    double maxr = 0.0;
    for (i=0; i<N; i++)
        if (rad[i] > maxr) maxr = rad[i];
    double R = 2*maxr;
    double R2 = R*R;

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

    double colmax = 0.0;
    for (t=0.0; t<time_end; t+=dt){
        int index[2];
        for (i=0; i<size_total; i++)
            count[i] = 0;

        for (i=0; i<N; i++){
            coords_to_index(&x[2*i], size, index, L);
            int t = index[0] + index[1]*size[0];
            cells[NMAX*t + count[t]] = i;
            count[t]++; 
        }

        int tt[2];
        int tix[2];
        int image[2];
        double dx[2];
        int goodcell, ind, n;
        double dist;

        #ifdef OPENMP
        #pragma omp parallel for private(i,dx,index,tt,goodcell,tix,ind,j,n,image,k,dist)
        #endif 
        for (i=0; i<N; i++){
            f[2*i+0] = 0.0;
            f[2*i+1] = 0.0;
            w[2*i+0] = 0.0;
            w[2*i+1] = 0.0;
            neigh[i] = 0;
            
            coords_to_index(&x[2*i], size, index, L);

            for (tt[0]=-1; tt[0]<=1; tt[0]++){
            for (tt[1]=-1; tt[1]<=1; tt[1]++){
                goodcell = 1;    
                for (j=0; j<2; j++){
                    tix[j] = mod_rvec(index[j]+tt[j],size[j]-1,pbc[j],&image[j]);
                    if (pbc[j] < image[j])
                        goodcell=0;
                }

                if (goodcell){
                    ind = tix[0] + tix[1]*size[0]; 

                    for (j=0; j<count[ind]; j++){
                        n = cells[NMAX*ind+j];

                        dist = 0.0;
                        for (k=0; k<2; k++){
                            dx[k] = x[2*n+k] - x[2*i+k];
                    
                            if (image[k])
                                dx[k] += L*tt[k];
                            dist += dx[k]*dx[k];
                        }

                        //===============================================
                        // force calculation 
                        if (dist > 1e-10 && dist < R2){
                            force_hertz(dx, dist, rad[i], rad[n], type[i], type[n], &f[2*i]);
                            col[i] += f[2*i+0]*f[2*i+0] + f[2*i+1]*f[2*i+1]; 
                        }
                    }
                }
            } } 

            //====================================
            // global forces 
            force_damping(&v[2*i], &f[2*i]);
            force_thermal(Tglobal, &f[2*i]);
            force_kick(&o[2*i], &f[2*i]);
            //force_gravity(&f[2*i]);
            if (col[i] > colmax) colmax = col[i];
            o[2*i+0] = 0.0; o[2*i+1] = 0.0;
        }
        #ifdef OPENMP
        #pragma omp barrier
        #endif

        // now integrate the forces since we have found them
        #ifdef OPENMP
        #pragma omp parallel for private(j)
        #endif 
        for (i=0; i<N;i++){
            // Newton-Stomer-Verlet
            #ifdef PLOT
            if (key['h'] != 1){
            #endif
            v[2*i+0] += f[2*i+0] * dt;
            v[2*i+1] += f[2*i+1] * dt;

            x[2*i+0] += v[2*i+0] * dt;
            x[2*i+1] += v[2*i+1] * dt;
            #ifdef PLOT
            }   
            #endif

            // boundary conditions 
            for (j=0; j<2; j++){
                if (pbc[j] == 1){
                    if (x[2*i+j] >= L-EPSILON || x[2*i+j] < 0)
                        x[2*i+j] = mymod(x[2*i+j], L);
                }
                else {
                    const double restoration = 0.5;
                    if (x[2*i+j] >= L){x[2*i+j] = 2*L-x[2*i+j]; v[2*i+j] *= -restoration;}
                    if (x[2*i+j] < 0) {x[2*i+j] = -x[2*i+j];    v[2*i+j] *= -restoration;}
                    if (x[2*i+j] >= L-EPSILON || x[2*i+j] < 0){x[2*i+j] = mymod(x[2*i+j], L);}
                }
            }

            // just check for errors
            if (x[2*i+0] >= L || x[2*i+0] < 0.0 ||
                x[2*i+1] >= L || x[2*i+1] < 0.0)
                printf("out of bounds\n");
            
            col[i] = col[i]/3;//colmax; 
        }
        #ifdef OPENMP
        #pragma omp barrier
        #endif

        #ifdef PLOT 
            plot_clear_screen();
            key = plot_render_particles(x, rad, type, N, L,col);
        #endif
        frames++;

        #ifdef PLOT
        if (key['q'] == 1)
            break;
        if (key['w'] == 1){
            for (i=0; i<N; i++){
                if (type[i] == RED)
                    o[2*i+1] = -kickforce;
            }
        }
        if (key['s'] == 1){
            for (i=0; i<N; i++){
                if (type[i] == RED)
                    o[2*i+1] = kickforce;
            }
        }
        if (key['a'] == 1){
            for (i=0; i<N; i++){
                if (type[i] == RED)
                    o[2*i+0] = -kickforce;
            }
        }
        if (key['d'] == 1){
            for (i=0; i<N; i++){
                if (type[i] == RED)
                    o[2*i+0] = kickforce;
            }
        }
        if (key['9'] == 1)
            Tglobal -= 0.01;
        if (key['0'] == 1)
            Tglobal += 0.01;
        if (key['8'] == 1)
            Tglobal = 0.0;
        if (key['t'] == 1)
            dt *= 1.2;
        if (key['y'] == 1)
            dt /= 1.2;
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
    free(neigh);
    free(rad);
    free(type);
    free(col);

    #ifdef PLOT
    plot_clean(); 
    #endif
}





