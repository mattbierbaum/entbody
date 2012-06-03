#include "neighborlist.h"

nbl_struct_cell *nbl_cell_build(int N, float L, int *pbc, float *x){
    nbl_struct_cell *nsc = (nbl_struct_cell*)malloc(sizeof(nbl_struct_cell));
    nsc->size = (unsigned int*)malloc(sizeof(unsigned int)*DIM);

    nsc->size_total = 1;
    for (i=0; i<DIM; i++){
        nsc->size[i] = (int)(L / R); 
        nsc->size_total *= nsc->size[i];
    }

    nsc->count = (unsigned int*)malloc(sizeof(unsigned int)*nsc->size_total);
    nsc->cells = (unsigned int*)malloc(sizeof(unsigned int)*nsc->size_total*NMAX);
    for (i=0; i<nsc->size_total; i++)
        nsc->count[i] = 0;
    for (i=0; i<nsc->size_total*NMAX; i++)
        nsc->cells[i] = 0;

    return nsc;
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
