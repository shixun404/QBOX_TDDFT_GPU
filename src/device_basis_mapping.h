#ifndef DEVICE_BASIS_MAPPING_H
#define DEVICE_BASIS_MAPPING_H

void cuda_vector_to_zvec (const double *c,
  double *zvec, const int* ip_, const int ng, cudaStream_t stream);


#endif
