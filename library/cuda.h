
#ifdef CUDA
//=============================================
// these are chosen if there is cuda
//=============================================
#define CUDA_GLOBAL __global__
#define CUDA_DEVICE __device__

#define CUDA_CREATE 
#define CUDA_FREE   
#define CUDA_H2D    CUDA_SAFE_CALL()
#define CUDA_D2H    CUDA_SAFE_CALL()

#else
//=============================================
// blanks are used if we don't need them
//=============================================
#define CUDA_GLOBAL 
#define CUDA_DEVICE 

#define CUDA_CREATE
#define CUDA_FREE
#define CUDA_D2H
#define CUDA_H2D

#endif

#define ERROR_CHECK { cudaError_t err; \
  if ((err = cudaGetLastError()) != cudaSuccess) { \
    printf("CUDA error: %s, line %d\n", cudaGetErrorString(err), __LINE__);}}

//=============================================
// this is too ugly to put in the main func
//=============================================
#define CUDA_MEMORY_CREATE                                      \
    int mem_size2 = sizeof(int)*2;                              \
    int imem_size = sizeof(int)*N;                              \
    int fmem_size = sizeof(float)*N;                            \
    int fmem_siz2 = sizeof(float)*N*2;                          \
    int mem_cell  = sizeof(unsigned int)*size_total;            \
    int mem_cell2 = sizeof(unsigned int)*size_total*NMAX;   \
                                                                \
    unsigned int *cu_count  = NULL;                             \
    unsigned int *cu_cells  = NULL;                             \
    int *cu_size   = NULL;                                      \
    int *cu_type   = NULL;                                      \
    int *cu_key    = NULL;                                      \
    float *cu_rad  = NULL;                                      \
    float *cu_col  = NULL;                                      \
    float *cu_x    = NULL;                                      \
    float *cu_copyx    = NULL;                                      \
    float *cu_v    = NULL;                                      \
    int *cu_pbc    = NULL;                                      \
                                                                \
    cudaMalloc((void**) &cu_pbc,   2*sizeof(int));  \
    cudaMalloc((void**) &cu_count, mem_cell);                   \
    cudaMalloc((void**) &cu_cells, mem_cell2);                  \
    cudaMalloc((void**) &cu_size,  mem_size2);                  \
                                                                \
    cudaMalloc((void**) &cu_key,   mem_size_k);                 \
    cudaMalloc((void**) &cu_type,  imem_size);                  \
    cudaMalloc((void**) &cu_rad,   fmem_size);                  \
    cudaMalloc((void**) &cu_col,   fmem_size);                  \
    cudaMalloc((void**) &cu_x,     fmem_siz2);                  \
    cudaMalloc((void**) &cu_copyx,     fmem_siz2);                  \
    cudaMalloc((void**) &cu_v,     fmem_siz2);                  \
                                                                \
    printf("Copying problem...\n");                             \
    cudaMemcpy(cu_size,  size,  mem_size2, cudaMemcpyHostToDevice);\
    cudaMemcpy(cu_type,  type,  imem_size, cudaMemcpyHostToDevice);\
    cudaMemcpy(cu_rad,   rad,   fmem_size, cudaMemcpyHostToDevice);\
    cudaMemcpy(cu_col,   col,   fmem_size, cudaMemcpyHostToDevice);\
    cudaMemcpy(cu_x,     x,     fmem_siz2, cudaMemcpyHostToDevice);\
    cudaMemcpy(cu_v,     v,     fmem_siz2, cudaMemcpyHostToDevice);\
    cudaMemcpy(cu_pbc, pbc,   2*sizeof(int), cudaMemcpyHostToDevice);\
    cudaMemset(cu_count, 0, mem_cell);                          \
    cudaMemset(cu_cells, 0, mem_cell2);                         \
    ERROR_CHECK                                                 


#define CUDA_MEMORY_FREE \
    cudaFree(cu_count);\
    cudaFree(cu_cells);\
    cudaFree(cu_key);\
    cudaFree(cu_type);\
    cudaFree(cu_rad);\
    cudaFree(cu_col);\
    cudaFree(cu_x);\
    cudaFree(cu_v);\
    cudaFree(cu_copyx);\
    ERROR_CHECK


#define CUDA_DEVICE_FUNCS \
__device__ int mod_rvec(int a, int b, int p, int *image){\
    *image = 1;\
    if (b==0) {if (a==0) *image=0; return 0;}\
    if (p != 0){\
        if (a>b)  return a-b-1;\
        if (a<0)  return a+b+1;\
    } else {\
        if (a>b)  return b;\
        if (a<0)  return 0;\
    }\
    *image = 0;\
    return a;\
}\
__device__ float mymod(float a, float b){\
  return a - b*(int)(a/b) + b*(a<0);\
}



#  define CUDA_SAFE_CALL( call) {                                            \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } }


#  define CUT_DEVICE_INIT() {                                                \
    int deviceCount;                                                         \
    CUDA_SAFE_CALL(cudaGetDeviceCount(&deviceCount));                        \
    if (deviceCount == 0) {                                                  \
        fprintf(stderr, "cutil error: no devices supporting CUDA.\n");       \
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    int dev = 0;                                                             \
    if (dev < 0) dev = 0;                                                    \
    if (dev > deviceCount-1) dev = deviceCount - 1;                          \
    cudaDeviceProp deviceProp;                                               \
    CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, dev));               \
    if (deviceProp.major < 1) {                                              \
        fprintf(stderr, "cutil error: device does not support CUDA.\n");     \
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    CUDA_SAFE_CALL(cudaSetDevice(dev));                                      \
}

