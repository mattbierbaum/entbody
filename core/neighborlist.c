#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "../core/util.h"
#include "neighborlist.h"

nbl_struct_cell *nbl_cell_build(int N, float L, int dim, float R, int *pbc, float *x){
    int i;

    nbl_struct_cell *nsc = (nbl_struct_cell*)malloc(sizeof(nbl_struct_cell));
    nsc->rij     = (float*)malloc(sizeof(float)*N*dim*NMAXTOTAL);
    nsc->nij     = (unsigned int*)malloc(sizeof(unsigned int)*N*NMAXTOTAL);
    nsc->countij = (unsigned int*)malloc(sizeof(unsigned int)*N);
    nsc->size    = (unsigned int*)malloc(sizeof(unsigned int)*dim);

    nsc->size_total = 1;
    for (i=0; i<dim; i++){
        nsc->size[i] = (int)(L / R); 
        nsc->size_total *= nsc->size[i];
    }

    nsc->count   = (unsigned int*)malloc(sizeof(unsigned int)*nsc->size_total);
    nsc->cells   = (unsigned int*)malloc(sizeof(unsigned int)*nsc->size_total*NMAX);
    for (i=0; i<nsc->size_total; i++)
        nsc->count[i] = 0;
    for (i=0; i<nsc->size_total*NMAX; i++)
        nsc->cells[i] = 0;

    return nsc;
}

int nbl_cell_update(int N, float L, int dim, float R, int *pbc, float *x, nbl_struct_cell *nsc){
    int i,j;
    int index[3];

    float dx[3];
    int tt[3], tix[3], image[3];
    int goodcell, ind, tn;

    for (i=0; i<N; i++){
        index[0] = (int)(x[2*i+0]/L  * nsc->size[0]); 
        index[1] = (int)(x[2*i+1]/L  * nsc->size[1]); 
        if (index[0] >= nsc->size[0]) index[0] = nsc->size[0]-1;        
        if (index[1] >= nsc->size[1]) index[1] = nsc->size[1]-1;        
        if (index[0] <  0)            index[0] = 0;        
        if (index[1] <  0)            index[1] = 0;        
         
        int t = index[0] + index[1]*nsc->size[0];
        int tcount = nsc->count[t];
        nsc->count[t] = nsc->count[t]+1;
        nsc->cells[NMAX*t + tcount] = i;
    }   

    for (i=0; i<N; i++){
        index[0] = (int)(x[2*i+0]/L * nsc->size[0]);
        index[1] = (int)(x[2*i+1]/L * nsc->size[1]);

        for (tt[0]=-1; tt[0]<=1; tt[0]++){
        for (tt[1]=-1; tt[1]<=1; tt[1]++){
            goodcell = 1;
            tix[0] = mod_rvec(index[0]+tt[0],nsc->size[0]-1,pbc[0],&image[0]);
            tix[1] = mod_rvec(index[1]+tt[1],nsc->size[1]-1,pbc[1],&image[1]);
            if (pbc[0] < image[0] || pbc[1] < image[1])  goodcell = 0;

            if (goodcell){
                ind = tix[0] + tix[1]*nsc->size[0]; 
                for (j=0; j<nsc->count[ind]; j++){
                    tn = nsc->cells[NMAX*ind+j];

                    dx[0] = x[2*tn+0] - x[2*i+0]; 
                    if (image[0]) dx[0] += L*tt[0];
                    
                    dx[1] = x[2*tn+1] - x[2*i+1];
                    if (image[1]) dx[1] += L*tt[1];

                    double dist = dx[0]*dx[0] + dx[1]*dx[1];
                    if (dist > EPSILON){
                        nsc->rij[NMAXTOTAL*i+nsc->countij[i]+0] = dx[0];
                        nsc->rij[NMAXTOTAL*i+nsc->countij[i]+1] = dx[1];
                        nsc->nij[NMAXTOTAL*i+nsc->countij[i]]   = tn;
                        nsc->countij[i]++;
                    }
                 }
            }
        } } 
    }
    return 0;
}

int nbl_cell_neighbors(int i, int **neighs, float **rij, int dim, nbl_struct_cell *nsc){
    *neighs  = &nsc->nij[NMAXTOTAL*i];
    *rij     = &nsc->rij[dim*NMAXTOTAL*i];
    return nsc->countij[i]; 
}

int nbl_cell_reset(int N, nbl_struct_cell *nsc){
    int i;
    for (i=0; i<nsc->size_total; i++)
        nsc->count[i] = 0;
    for (i=0; i<N; i++)
        nsc->countij[i] = 0;
    return 0;
}

int nbl_cell_free(nbl_struct_cell *nsc){
    free(nsc->size);
    free(nsc->count);
    free(nsc->cells);
    free(nsc->rij);
    free(nsc->nij);
    free(nsc->countij);
    free(nsc);
}

int nbl_cell_mod_rvec(int a, int b, int p, int *image){
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

int nbl_cell_coords_to_index(float *x, int *size, int *index, float L){
    index[0] = (int)(x[0]/L  * size[0]);
    index[1] = (int)(x[1]/L  * size[1]);
}
