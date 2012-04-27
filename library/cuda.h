
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
