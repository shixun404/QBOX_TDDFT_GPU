#ifndef CUDA_UTILS_H
#define CUDA_UTILS_H
void cuda_error_check (cudaError_t cudaError, const char* const file, int line);
void cuda_check_last (const char* const file, const int line);
#endif
