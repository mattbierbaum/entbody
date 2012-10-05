#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#ifdef PLOT
#include "plot.h"
#endif

#ifdef FPS
#include <time.h>
#endif

#define pi      3.141592653589
#define BLACK   0
#define RED     1
#define EPSILON DBL_EPSILON

typedef unsigned long long int ullong;

//-----------------------------------------------------------
// some defines and what not 
//------------------------------------------------------------
void   ran_seed(long j);
double ran_ran2();
ullong vseed;
ullong vran;

void simulate(int s);

inline double mymod(double a, double b){
  return a - b*(int)(a/b) + b*(a<0);
}

inline void coords_to_index(double *x, int *size, int *index, double L){   
    index[0] = (int)(x[0]/L  * size[0]);
    index[1] = (int)(x[1]/L  * size[1]);
}

inline int mod_rvec(int a, int b, int p, int *image){
    *image = 1;
    if (b==0) {if (a==0) *image=0; return 0;}
    if (p != 0){
        if (a>b)  return a-b-1;
        if (a<0)  return a+b+1;
    } else {
        if (a>b)  return b;
        if (a<0)  return 0;
    }
    *image = 0;
    return a;
}

//===================================================
// the main function
//===================================================
int main(int argc, char **argv){
    int seed_in     = 0;

    if (argc == 1) 
        simulate(seed_in);
    else if (argc == 4){
        seed_in  = atoi(argv[1]);
        simulate(seed_in);
    }
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

    int    N       = 18*512;
    double radius  = 1.0;
    double L       = sqrt(1.02*pi*radius*radius*N);
    int    Npercell = 50;

    int pbc[] = {1,1};

    double epsilon    = 180.0;
    double damp_coeff = 0.3;
    double kickforce  = 20.0;

    double dt = 1e-1;
    double t  = 0.0;
    double R  = 2*radius; 
    double R2 = R*R;

    int i, j, k;
    int *key;

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
        plot_init(); 
        plot_clear_screen();
        key = plot_render_particles(x, rad, type, N, L,col);
    #endif

    //-------------------------------------------------
    // initialize
    for (i=0; i<N; i++){
        rad[i] = radius;
        x[2*i+0] = L*ran_ran2();
        x[2*i+1] = L*ran_ran2();
    
        v[2*i+0] = 0.0;
        v[2*i+1] = 0.0;
        type[i] = BLACK;

        if (i==0) type[i] = RED;
    }

    //-------------------------------------------------------
    // make boxes for the neighborlist
    int size[2];
    int size_total = 1;
    for (i=0; i<2; i++){
        size[i] = (int)(L / (R)); 
        size_total *= size[i];
    }

    int *count = (int*)malloc(sizeof(int)*2*size_total);
    int *cells = (int*)malloc(sizeof(int)*2*size_total*Npercell);
    for (i=0; i<size_total; i++){
        for (j=0; j<Npercell; j++)
            cells[i*Npercell + j] = 0;
        count[i] = 0;
    }

    //==========================================================
    // where the magic happens
    //==========================================================
    int frames = 0;

    #ifdef FPS
    struct timespec start;
    clock_gettime(CLOCK_REALTIME, &start);
    #endif

    double T = 0.0;
    for (t=0.0; t<time_end; t+=dt){
        double colmax = 0.0;
        int index[2];

        for (i=0; i<size_total; i++)
            count[i] = 0;

        for (i=0; i<N; i++){
            col[i] = 0.0;
            coords_to_index(&x[2*i], size, index, L);
            int t = index[0] + index[1]*size[0];
            cells[t*Npercell+count[t]] = i;
            count[t]++; 
        }

        int tt[2];
        int tix[2];
        int image[2];
        double dx[2];

        for (i=0; i<N; i++){
            f[2*i+0] = 0.0;
            f[2*i+1] = 0.0;
            
            coords_to_index(&x[2*i], size, index, L);

            for (tt[0]=-1; tt[0]<=1; tt[0]++){
            for (tt[1]=-1; tt[1]<=1; tt[1]++){
                int goodcell = 1;    
                for (j=0; j<2; j++){
                    tix[j] = mod_rvec(index[j]+tt[j],size[j]-1,pbc[j],&image[j]);
                    if (pbc[j] < image[j])
                        goodcell=0;
                }

                if (goodcell){
                    int ind = tix[0] + tix[1]*size[0]; 

                    for (j=0; j<count[ind]; j++){
                        int n = cells[ind*Npercell+j];

                        double dist = 0.0;
                        for (k=0; k<2; k++){
                            dx[k] = x[2*n+k] - x[2*i+k];
                    
                            if (image[k])
                                dx[k] += L*tt[k];
                            dist += dx[k]*dx[k];
                        }

                        //===============================================
                        // force calculation - hertz
                        if (dist > 1e-10 && dist < R2){
                            double r0 = R; 
                            double l  = sqrt(dist);
                            double co = epsilon * (1-l/r0)*(1-l/r0) * (l<r0);
                            for (k=0; k<2; k++){
                                f[2*i+k] += - dx[k] * co;
                                col[i] += co*co*dx[k]*dx[k]; 
                            }
                        }
                    }
                }
            } } 

            //====================================
            // damping
            f[2*i+0] -= damp_coeff*v[2*i+0];
            f[2*i+1] -= damp_coeff*v[2*i+1];

            //=====================================
            // noise
            f[2*i+0] += T*(ran_ran2()-0.5);
            f[2*i+1] += T*(ran_ran2()-0.5);

            //=====================================
            // kick force
            f[2*i+0] += o[2*i+0]; o[2*i+0] = 0.0;
            f[2*i+1] += o[2*i+1]; o[2*i+1] = 0.0;
            
            //=======================
            // color norm
            if (col[i] > colmax) colmax = col[i];
        }

        // now integrate the forces since we have found them
        for (i=0; i<N;i++){
            // Newton-Stomer-Verlet
            if (key['h'] != 1){
            v[2*i+0] += f[2*i+0] * dt;
            v[2*i+1] += f[2*i+1] * dt;

            x[2*i+0] += v[2*i+0] * dt;
            x[2*i+1] += v[2*i+1] * dt;
            }
            // boundary conditions 
            for (j=0; j<2; j++){
                if (pbc[j] == 1){
                    if (x[2*i+j] >= L-EPSILON || x[2*i+j] < 0)
                        x[2*i+j] = mymod(x[2*i+j], L);
                }
                else {
                    const double restoration = 1.0;
                    if (x[2*i+j] >= L){x[2*i+j] = 2*L-x[2*i+j]; v[2*i+j] *= -restoration;}
                    if (x[2*i+j] < 0) {x[2*i+j] = -x[2*i+j];    v[2*i+j] *= -restoration;}
                    if (x[2*i+j] >= L-EPSILON || x[2*i+j] < 0){x[2*i+j] = mymod(x[2*i+j], L);}
                }
            }

            // just check for errors
            if (x[2*i+0] >= L || x[2*i+0] < 0.0 ||
                x[2*i+1] >= L || x[2*i+1] < 0.0)
                printf("out of bounds\n");
          
           //col[i] += v[2*i+0]*v[2*i+0] + v[2*i+1]*v[2*i+1]; 
           col[i] = col[i]/8; 
        }

        #ifdef PLOT 
        const int FRAMESKIP = 1;
        if (frames % FRAMESKIP == 0){
            char filename[100];
            sprintf(filename, "screen_%05d.png", frames/FRAMESKIP);
            plot_clear_screen();
            key = plot_render_particles(x, rad, type, N, L,col);
            #ifdef IMAGES
            plot_saveimage(filename);
            #endif
        }
        #endif
        frames++;

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
            T -= 0.01;
        if (key['0'] == 1)
            T += 0.01;
        if (key['8'] == 1)
            T = 0.0;
        if (key['o'] == 1)
            L -= 0.01;
        if (key['p'] == 1)
            L += 0.01;
    }
    // end of the magic, cleanup
    //----------------------------------------------
    #ifdef FPS
    struct timespec end;
    clock_gettime(CLOCK_REALTIME, &end);
    printf("fps = %f\n", frames/(end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec)/1e9);
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

    #ifdef PLOT
    plot_clean(); 
    #endif
}




//=================================================
// extra stuff
//=================================================
void ran_seed(long j){
  vseed = j;  vran = 4101842887655102017LL;
  vran ^= vseed; 
  vran ^= vran >> 21; vran ^= vran << 35; vran ^= vran >> 4;
  vran = vran * 2685821657736338717LL;
}

double ran_ran2(){
    vran ^= vran >> 21; vran ^= vran << 35; vran ^= vran >> 4;
    ullong t = vran * 2685821657736338717LL;
    return 5.42101086242752217e-20*t;
}


