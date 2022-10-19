#include "device_basis_mapping.h"
#include <stdio.h>



//TODO: Think in another coalescing-friendly version
template <int unroll>
__global__ void vector_to_zvec_kernel(const double *c, double *zvec, const int* ip_, const int ng)
{
        const int global_id=blockDim.x*blockIdx.x*unroll+threadIdx.x*unroll;


        #pragma unroll
        for(int i=0;i<unroll;i++)
        {
                const int id=global_id+i;
                if(id<ng)
                {
                        const int ip=ip_[id];
                        zvec[2*ip]=c[2*id];
                        zvec[2*ip+1]=c[2*id+1];
                }
        }
}


void cuda_vector_to_zvec (const double *c,
  double *zvec, const int* ip_, const int ng, cudaStream_t  stream)
{
        constexpr int THREADS_PER_BLOCK=128;
        constexpr int UNROLL_FACTOR=8;

        //round-up
        const int block_num=(ng+THREADS_PER_BLOCK*UNROLL_FACTOR-1)/(THREADS_PER_BLOCK*UNROLL_FACTOR);


        dim3 threads (THREADS_PER_BLOCK);
        dim3 blocks (block_num);
        vector_to_zvec_kernel<UNROLL_FACTOR><<<blocks,threads/*,0,stream*/>>>(c,zvec,ip_,ng);
	if (cudaGetLastError() != cudaSuccess){
                   fprintf(stderr, "Cuda error device_bm: Failed kernel\n");
                   exit(-1);
                }

	
	return;
}





