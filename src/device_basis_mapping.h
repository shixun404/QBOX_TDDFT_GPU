#ifndef DEVICE_BASIS_MAPPING_H
#define DEVICE_BASIS_MAPPING_H


#if AUTOTUNER

template<int unroll, int tb, int tb_sm>
__global__ __launch_bounds__(tb,tb_sm) void  zvec_to_vector_kernel(const double *zvec, double *c, const int* ip_, const int ng, const int nvecnp2);
#endif

void cuda_vector_to_zvec (const double *c,
  double *zvec, const int* ip_, const int * im_, const int ng, const int dst_size,cudaStream_t stream, const int batch, const int PLAN);



#if AUTOTUNER
template<int D, int tb, int tb_sm>
__global__ __launch_bounds__(tb,tb_sm) void cu_Z_copy(const int count, const double* x, const int offset_x, const int incx, double * y, const int offset_y, const int incy, const int *offset, const int N_x, const int N_y);
#endif


void cuZcopy(const int count, const double * x, const int offset_x,const int incx, double * y, const int  offset_y, const int incy, const int * offset, const int batchesX, cudaStream_t stream, const int batchY, const int PLAN, const int max_blocks,const int n1, const int n2);


#if AUTOTUNER
template<int D, int tb, int tb_sm>
__global__ __launch_bounds__(tb,tb_sm) void vector_to_zvec_kernel(const double *c, double *zvec, const int* ip_, const int* im_, const int ng, const int n_dest);
#endif


void cuda_zvec_to_vector(const double * zvec, double *c, const int* ip_, const int ng, const int nvecnp2, cudaStream_t stream, const int batch=1);


#if AUTOTUNER
template<int unroll, int tb, int tb_sm>
__global__ __launch_bounds__(tb,tb_sm) void cu_pairwise(const double* src, double* dest, const int N, const int batch);
#endif


void cuPairwise(const double* src, double* dest, const int N, cudaStream_t stream, const int batch=1);

#if AUTOTUNER
template<int unroll, int tb, int tb_sm>
__global__ __launch_bounds__(tb,tb_sm) void cu_dscal(const double scal, double* dest, const int N, const int batch);
#endif


void cuDscal(const double scal, double* dest, const int N, cudaStream_t stream, const int batch=1);

#endif
