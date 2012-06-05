#ifndef __NEIGHBORLIST_H__
#define __NEIGHBORLIST_H__

#define NMAX        50
#define NMAXTOTAL   (NMAX*9)

typedef struct {
    // counts for single cells
    unsigned int *size;
    unsigned int size_total;
    unsigned int *count;
    unsigned int *cells;    

    // count after cells totalled
    unsigned int *countij; 
    unsigned int *nij;
    float *rij;
} nbl_struct_cell;

nbl_struct_cell *nbl_cell_build(int N, float L, int dim, float R, int *pbc, float *x);
int nbl_cell_update(int N, float L, int dim, float R, int *pbc, float *x, nbl_struct_cell *nsc);
int nbl_cell_neighbors(int i, int **neighs, float **rij, int dim, nbl_struct_cell *nsc);
int nbl_cell_reset(int N, nbl_struct_cell *nsc);
int nbl_cell_free(nbl_struct_cell *nsc);

int nbl_cell_mod_rvec(int a, int b, int p, int *image);
int nbl_cell_coords_to_index(float *x, int *size, int *index, float L);


#endif
