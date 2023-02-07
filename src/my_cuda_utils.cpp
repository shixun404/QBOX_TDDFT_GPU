#include <iostream>
#include <unistd.h>
#include "my_cuda_utils.h"



void cuda_error_check (cudaError_t cudaError, const char* const file, int line){

  int         pid;
  size_t      free, total;
  cudaError_t cErr2;

  cErr2 = cudaGetLastError();
  if (cudaError != cudaSuccess || cErr2 != cudaSuccess) {
    pid = getpid();
    fprintf(stderr,"%d CUDA %s Error line %d\n", pid,file, line);
    fprintf(stderr,"%d CUDA  Code Error: %s\n", pid, cudaGetErrorString(cudaError));
    fprintf(stderr,"%d CUDA  Last Error: %s\n", pid, cudaGetErrorString(cErr2));
    cudaMemGetInfo(&free,&total);
    fprintf(stderr,"%d Free: %zu , Total: %zu\n", pid, free, total);
    fflush(stderr);
    exit(-1);
  }


}


void cuda_check_last(const char* const file, const int line)
{
    cudaError_t err{cudaGetLastError()};
    if (err != cudaSuccess)
    {
        std::cerr << "CUDA Runtime Error at: " << file << ":" << line
                  << std::endl;
        std::cerr << cudaGetErrorString(err) << std::endl;
        std::exit(EXIT_FAILURE);
    }
}
