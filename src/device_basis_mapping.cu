#include "device_basis_mapping.h"
#include <stdio.h>


//TODO: template it for just doubles and not complex



template <int unroll>
__global__ void zvec_to_vector_kernel(const double *zvec, double *c, const int* ip_, const int ng, const int nvecnp2)
{
	const int global_id=blockDim.x*blockIdx.x*unroll+threadIdx.x*unroll;

	#pragma unroll
        for (int i=0; i<unroll;i++)
	{
		const int id = global_id+i;
		if(id<ng)
		{
			const int ip = ip_[id];
			c[2*id+2*blockIdx.y*ng]=zvec[2*ip+2*blockIdx.y*nvecnp2];
			c[2*id +2*blockIdx.y*ng +1]=zvec[2*ip+2*blockIdx.y*nvecnp2+1];
		}
	}



}
void cuda_zvec_to_vector (const double * zvec, double *c, const int* ip_, const int ng, const int nvecnp2, cudaStream_t stream, const int batch){

	constexpr int THREADS_PER_BLOCK=128;
        constexpr int UNROLL_FACTOR=8;

        //round-up
        const int block_num=(ng+THREADS_PER_BLOCK*UNROLL_FACTOR-1)/(THREADS_PER_BLOCK*UNROLL_FACTOR);

        dim3 threads (THREADS_PER_BLOCK);
        dim3 blocks (block_num, batch);
        zvec_to_vector_kernel<UNROLL_FACTOR><<<blocks,threads,0,stream>>>(zvec,c,ip_,ng,nvecnp2);

}


/*********************************/

template<int D> struct Copy {
        static inline __device__ void 
        copy(double* dst, const double* src, const int *ip, const int *im, const int id, const int n, const int n2);

};

template<> struct Copy<0>{
        static inline __device__ void
        copy(double* dst, const double* src, const int *ip, const int *im, const int id, const int n, const int n2){
                const int index = ip[id];
		dst[2*index+2*blockIdx.y*n2]=src[2*id+2*blockIdx.y*n];
		dst[2*index+1+2*blockIdx.y*n2]=src[2*id+1+2*blockIdx.y*n];
   	}
};



template<> struct Copy<1>{
        static inline __device__ void
        copy(double* dst, const double* src, const int *ip, const int *im, const int id, const int n, const int n2){
                const int indexp = ip[id];
                const int indexm = im[id];

		dst[2*indexp+2*blockIdx.y*n2]=src[2*id+2*blockIdx.y*n];
                dst[2*indexp+1+2*blockIdx.y*n2]=src[2*id+1+2*blockIdx.y*n];

		dst[2*indexm+2*blockIdx.y*n2]=src[2*id+2*blockIdx.y*n];
		dst[2*indexm+1+2*blockIdx.y*n2]=-src[2*id+1+2*blockIdx.y*n];
        }
};


//TODO: Think in another coalescing-friendly version
template <int unroll,int D>
__global__ void vector_to_zvec_kernel(const double *c, double *zvec, const int* ip_, const int* im_, const int ng, const int n_dest)
{
        const int global_id=blockDim.x*blockIdx.x*unroll+threadIdx.x*unroll;

        #pragma unroll
        for(int i=0;i<unroll;i++)
        {
                const int id=global_id+i;
                if(id<ng)
                {
			Copy<D>::copy(zvec,c,ip_,im_,id,ng,n_dest);
                	
		}
        }
}


void cuda_vector_to_zvec (const double *c,
  double *zvec, const int* ip_, const int * im_, const int ng, const int dst_size, cudaStream_t  stream, const int batch, const int PLAN)
{
        constexpr int THREADS_PER_BLOCK=128;
        constexpr int UNROLL_FACTOR=8;

        //round-up
        const int block_num=(ng+THREADS_PER_BLOCK*UNROLL_FACTOR-1)/(THREADS_PER_BLOCK*UNROLL_FACTOR);
	//TODO: check the blocks are ok in the target architecture
	const int block_batch=batch;
        dim3 threads (THREADS_PER_BLOCK);
        dim3 blocks (block_num,block_batch);
	if(PLAN)
        	vector_to_zvec_kernel<UNROLL_FACTOR,1><<<blocks,threads,0,stream>>>(c,zvec,ip_,im_,ng,dst_size);
	else
		vector_to_zvec_kernel<UNROLL_FACTOR,0><<<blocks,threads,0,stream>>>(c,zvec,ip_,im_,ng,dst_size);
	if (cudaGetLastError() != cudaSuccess){
                   fprintf(stderr, "Cuda error device_bm: Failed kernel\n");
                   exit(-1);
                }

	
	return;
}
 

/***********************************/



template<int D> struct Offset {
	static inline __device__ int 
	calculateX(const int offsetX, const int *v, const int id);
	static inline __device__ int
	calculateY(const int offsetY, const int *v, const int id);

};
//Forward case
template<> struct Offset<-1>{
	static inline __device__ int 
	calculateX(const int offsetX, const int *v, const int id){
		return v[id]*2;
	}
	static inline __device__ int	
	calculateY(const int offsetY, const int *v, const int id){
		return id*offsetY*2;
	}

};
//backward case
template<> struct Offset<1>{
        static inline __device__ int
	calculateX(const int offsetX, const int *v, const int id){
                return id*2*offsetX; 
        }   
        static inline __device__ int
        calculateY(const int offsetY, const int *v, const int id){
        	return v[id]*2;
	}
};



//WARNING: It only works for Complex<double> (see stride *2 in code)


template<int D>
__global__ void cu_Z_copy(const int count, const double* x, const int offset_x, const int incx, double * y, const int offset_y, const int incy, const int *offset, const int N_x, const int N_y){

	const int src = Offset<D>::calculateX(offset_x, offset, blockIdx.x);
	const int dest = Offset<D>::calculateY(offset_y,offset, blockIdx.x);
	const int n_iters = ((count+blockDim.x-1)/blockDim.x);


	for(int i=0; i<n_iters; i++){
		const int thx=threadIdx.x + i*blockDim.x;
		if(thx<count){
			y[dest + 2*incy*thx+blockIdx.y*N_y*2] = x[src + 2*incx*thx+blockIdx.y*N_x*2];
			y[dest + 2*incy*thx + 1+blockIdx.y*N_y*2] = x[src + 2*incx*thx + 1+blockIdx.y*N_x*2];
		}
	}	



}

//WARNING: IT does not support negative increments
//https://netlib.org/blas/zcopy.f
void cuZcopy(const int count, const double * x, const int offset_x, const int incx, double * y, const int offset_y, const int incy, const int * offset, const int batchesX, cudaStream_t stream, const int batchY, const int PLAN, const int max_blocks,const int n1, const int n2){

	constexpr int THREADS_PER_BLOCK=128;

	dim3 threads(THREADS_PER_BLOCK);
	//TODO: check batches < max_blocks allowed for the architecture
	if (batchesX>max_blocks)
	{
		 fprintf(stderr, "Cuda error cu_zcopy: More blocks requested than available\n");
		 exit(-1);
	}
	
	dim3 blocks (batchesX,batchY);
 	if (PLAN==-1)	
		cu_Z_copy<-1><<<blocks, threads, 0, stream>>>(count,x,offset_x,incx,y,offset_y,incy,offset,n1,n2);
	if (PLAN==1)
		cu_Z_copy<1><<<blocks, threads, 0, stream>>>(count,x,offset_x,incx,y,offset_y,incy,offset,n1,n2);
	if (cudaGetLastError() != cudaSuccess){
                   fprintf(stderr, "Cuda error cu_zcopy: Failed kernel\n");
                   exit(-1);
        }

	return;

}	

/********************************************/

//TODO: templated for T type

template<int unroll>
__global__ void cu_pairwise(const double* src, double* dest, const int N, const int batch){

	const int thx = blockIdx.x*blockDim.x*unroll + threadIdx.x;

	for(int i=0;i<unroll;i++) 
	{
		const int idx = thx+i*blockDim.x;
		if(idx<N)
		{
			//Complex by scalar
			dest[2*idx+2*blockIdx.y*N] = dest[2*idx+2*blockIdx.y*N]*src[idx];
			dest[2*idx+1+2*blockIdx.y*N]=dest[2*idx+1+2*blockIdx.y*N]*src[idx];

			//Complex by Complex
			/*dest[2*idx+2*blockIdx.y*N] = dest[2*idx+2*blockIdx.y*N]*src[2*idx]-dest[2*idx+1+2*blockIdx.y*N]*src[2*idx+1];
			dest[2*idx+1+2*blockIdx.y*N] = dest[2*idx+2*blockIdx.y*N]*src[2*idx+1] + dest[2*idx+1+2*blockIdx.y*N]*src[2*idx] ;*/
		}
	}
	
}

void cuPairwise(const double* src, double* dest, const int N, cudaStream_t stream, const int batch){
  	constexpr int THREADS_PER_BLOCK=128;
        constexpr int UNROLL_FACTOR=8;

        //round-up
        const int block_num=(N+THREADS_PER_BLOCK*UNROLL_FACTOR-1)/(THREADS_PER_BLOCK*UNROLL_FACTOR);

	//Check total number of threads/blocks meets architecture requirements
        dim3 threads (THREADS_PER_BLOCK);
        dim3 blocks (block_num,batch);

        cu_pairwise<UNROLL_FACTOR><<<blocks,threads,0,stream>>>(src,dest,N,batch);
        if (cudaGetLastError() != cudaSuccess){
                   fprintf(stderr, "Cuda error cu_pairwise: Failed kernel\n");
                   exit(-1);
        }


}




/**************************/



//TODO: templated for T type
template<int unroll>
__global__ void cu_dscal(const double scal, double* dest, const int N, const int batch){

        const int thx = blockIdx.x*blockDim.x*unroll + threadIdx.x;

        for(int i=0;i<unroll;i++)
        {
                const int idx = thx+i*blockDim.x;
                if(idx<N)
                {
                        //Complex by scalar
                        dest[idx+blockIdx.y*N] = dest[idx+blockIdx.y*N]*scal;

                }
        }

}




void cuDscal(const double scal, double * dest, const int N, cudaStream_t stream, const int batch){

	constexpr int THREADS_PER_BLOCK=128;
        constexpr int UNROLL_FACTOR=8;

        //round-up
        const int block_num=(N+THREADS_PER_BLOCK*UNROLL_FACTOR-1)/(THREADS_PER_BLOCK*UNROLL_FACTOR);

        //Check total number of threads/blocks meets architecture requirements
        dim3 threads (THREADS_PER_BLOCK);
        dim3 blocks (block_num,batch);

        cu_dscal<UNROLL_FACTOR><<<blocks,threads,0,stream>>>(scal,dest,N,batch);
        if (cudaGetLastError() != cudaSuccess){
                   fprintf(stderr, "Cuda error cu_pairwise: Failed kernel\n");
                   exit(-1);
        }

}







