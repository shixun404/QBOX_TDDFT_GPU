#ifndef DEVICE_BASIS_MAPPING_H
#define DEVICE_BASIS_MAPPING_H

void cuda_vector_to_zvec (const double *c,
  double *zvec, const int* ip_, const int ng, cudaStream_t stream);

void cuZcopy(const int count, const double * x, const int offset_x,const int incx, double * y, const int  offset_y, const int incy, const int * offset, const int batches, cudaStream_t stream, const int PLAN);


void cuda_zvec_to_vector(const double * zvec, double *c, const int* ip_, const int ng, cudaStream_t stream);

void cuPairwise(const double* src, double* dest, const int N, cudaStream_t stream);

#endif
