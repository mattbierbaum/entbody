#ifndef __NEIGHBORLIST_H__
#define __NEIGHBORLIST_H__

#define NMAX 50

typedef struct {
    unsigned int *size;
    unsigned int size_total;
    unsigned int *count;
    unsigned int *cells;    
} nbl_struct_cell;

nbl_struct_cell *nbl_cell_build(int N, float L, int *pbc, float *x);
int nbl_cell_update(int N, float L, int *pbc, float *x, nbl_struct_cell *nsc);
int nbl_cell_neighbors(int i, int *pbc, int *neighs, int *num, nbl_struct_cell *nsc);
int nbl_cell_mod_rvec(int a, int b, int p, int *image);
int nbl_cell_coords_to_index(float *x, int *size, int *index, float L);
int nbl_cell_free(nbl_struct_cell *nsc);


#endif
