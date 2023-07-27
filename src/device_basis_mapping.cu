#include "device_basis_mapping.h"
#include <stdio.h>

//TODO: template it for just doubles and not complex

#if AUTOTUNER

#include "tuning_data_cuda.h"
#include "utilities.hpp"

void log2_x(int32_t m, int32_t & res)
{
        while(m>>=1)
                res++;
}


int getRow (const int32_t unroll_i,const int32_t tb_i , const int32_t tb_sm_i)
{
        int32_t res = 0;

        const int32_t m_j=((MAX_TB-STARTING_TB)/STRIDE_TB)+1;
        const int32_t m_k=((MAX_TB_SM-STARTING_TB_SM)/STRIDE_TB_SM)+1;

        log2_x(unroll_i/STARTING_UNROLL,res);
        const int32_t j = (tb_i-STARTING_TB)/STRIDE_TB;
        const int32_t k = (tb_sm_i-STARTING_TB_SM)/STRIDE_TB_SM;

        return res*m_j*m_k+j*m_k+k;
}





int32_t translate_unroll(const int i){
        return i;
}

int32_t translate_tb(const int32_t i, const int32_t warpSize){
	return i*warpSize;
}

int32_t translate_tb_sm(const int32_t i){
	return i;
}


int getIndex(const int32_t unroll_i, const int32_t tb_i, const int32_t warpSize , const int32_t tb_sm_i){
         return getRow(translate_unroll(unroll_i),translate_tb(tb_i, warpSize),translate_tb_sm(tb_sm_i));

}


	
#endif


#if AUTOTUNER
template<int unroll, int tb, int tb_sm>
__global__ __launch_bounds__(tb,tb_sm) void  zvec_to_vector_kernel(const double *zvec, double *c, const int* ip_, const int ng, const int nvecnp2)
#else
template <int unroll>
__global__ void zvec_to_vector_kernel(const double *zvec, double *c, const int* ip_, const int ng, const int nvecnp2)
#endif
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

#if AUTOTUNER
template <typename T>
int ZvecLauncher(const KernelCfg_zvec<T> &launchSetup, const double * zvec, double *c, const int* ip_, const int ng, const int nvecnp2, cudaStream_t stream, const int batch){

	void (*kerPtr) (const T*, T*,const int *,  const int, const int) = launchSetup.kerPtr;

	if(!kerPtr){
            fprintf(stderr, "ERROR: NO CUDA KERNEL FOUND FOR THE GIVEN VALUES at zVec \n");
            exit(-1); // Invalid configuration, abort execution so autotuner does not measure it
         }

	 const int THREADS_PER_BLOCK=launchSetup.TB;
         const int UNROLL_FACTOR=launchSetup.U;



         //round-up
        const int block_num=(ng+THREADS_PER_BLOCK*UNROLL_FACTOR-1)/(THREADS_PER_BLOCK*UNROLL_FACTOR);

        //Check total number of threads/blocks meets architecture requirements
        dim3 threads (THREADS_PER_BLOCK);
        dim3 blocks (block_num,batch);

        kerPtr<<<blocks,threads,0,stream>>>(zvec,c,ip_,ng,nvecnp2);
        return 0;



}

#endif




void cuda_zvec_to_vector (const double * zvec, double *c, const int* ip_, const int ng, const int nvecnp2, cudaStream_t stream, const int batch){
#if AUTOTUNER
	const int index = getIndex(atoi(DFTuning::getEnv("QBOX_ZVEC_unroll")),atoi(DFTuning::getEnv("QBOX_ZVEC_TB")), 32 ,atoi(DFTuning::getEnv("QBOX_ZVEC_TB_SM")));

	ZvecLauncher(weightTable_zvec[index],zvec,c,ip_,ng,nvecnp2,stream,batch);

#else
	constexpr int THREADS_PER_BLOCK=128;
        constexpr int UNROLL_FACTOR=8;

        //round-up
        const int block_num=(ng+THREADS_PER_BLOCK*UNROLL_FACTOR-1)/(THREADS_PER_BLOCK*UNROLL_FACTOR);

        dim3 threads (THREADS_PER_BLOCK);
        dim3 blocks (block_num, batch);
        zvec_to_vector_kernel<UNROLL_FACTOR><<<blocks,threads,0,stream>>>(zvec,c,ip_,ng,nvecnp2);
#endif

	if (cudaGetLastError() != cudaSuccess){
                   fprintf(stderr, "Cuda error device_bm: Failed kernel\n");
                   exit(-1);
                }
 

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
//TODO: Delete
	//	if((blockIdx.y==1)&&((2*index+2*blockIdx.y*n2)==1024020))
	//		printf("Block %d Th %d (id %d)- mi index %d  a partir de id %d con n2 %d y n %d voy a escribir en %d \n", blockIdx.x, threadIdx.x, id, index, id,n2,n,2*index+2*blockIdx.y*n2);
	
		dst[2*index+(2*blockIdx.y*n2)]=src[2*id+(2*blockIdx.y*n)];
		dst[2*index+1+(2*blockIdx.y*n2)]=src[2*id+1+(2*blockIdx.y*n)];
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

#if AUTOTUNER
template<int D, int tb, int tb_sm>
__global__ __launch_bounds__(tb,tb_sm) void vector_to_zvec_kernel(const double *c, double *zvec, const int* ip_, const int* im_, const int ng, const int n_dest)
#else
template <int D>
__global__ void vector_to_zvec_kernel(const double *c, double *zvec, const int* ip_, const int* im_, const int ng, const int n_dest)
#endif
{
	constexpr int unroll=8;
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

#if AUTOTUNER
template<typename T>
int VecLauncher (const KernelCfg_vec<T> &launchSetup,const double *c, 
  double *zvec, const int* ip_, const int * im_, const int ng, const int dst_size, cudaStream_t  stream, const int batch)
{
	 void (*kerPtr) ( const T*, T*, const int*,const int*, const int, const int) = launchSetup.kerPtr;

	  if(!kerPtr){
            fprintf(stderr, "ERROR: NO CUDA KERNEL FOUND FOR THE GIVEN VALUES at VecToZVec \n");
            exit(-1); // Invalid configuration, abort execution so autotuner does not measure it
         }

         const int THREADS_PER_BLOCK=launchSetup.TB;
	
	constexpr int UNROLL_FACTOR=8;
	const int block_num=(ng+THREADS_PER_BLOCK*UNROLL_FACTOR-1)/(THREADS_PER_BLOCK*UNROLL_FACTOR);

        dim3 threads (THREADS_PER_BLOCK);
        dim3 blocks (block_num,batch);

        kerPtr<<<blocks,threads,0,stream>>>(c,zvec,ip_,im_,ng,dst_size);
        return 0;

}	

#endif

void cuda_vector_to_zvec (const double *c,
  double *zvec, const int* ip_, const int * im_, const int ng, const int dst_size, cudaStream_t  stream, const int batch, const int PLAN)
{

#if AUTOTUNER
	const int index = getIndex(PLAN+1,atoi(DFTuning::getEnv("QBOX_VEC_TB")), 32 ,atoi(DFTuning::getEnv("QBOX_VEC_TB_SM")));
	VecLauncher(weightTable_vec[index], c, zvec, ip_,  im_, ng, dst_size,  stream, batch);

#else
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


#endif

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

#if AUTOTUNER
template<int D, int tb, int tb_sm>
__global__ __launch_bounds__(tb,tb_sm) void cu_Z_copy(const int count, const double* x, const int offset_x, const int incx, double * y, const int offset_y, const int incy, const int *offset, const int N_x, const int N_y){
#else
template<int D>
__global__ void cu_Z_copy(const int count, const double* x, const int offset_x, const int incx, double * y, const int offset_y, const int incy, const int *offset, const int N_x, const int N_y){
#endif
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

#if AUTOTUNER
template<typename T>
int ZcopyLauncher(const KernelCfg_zcopy<T> &launchSetup,const int count, const T* x, const int offset_x, const int incx, T* y, const int offset_y, const int incy, const int * offset, const int batchesX, cudaStream_t stream, const int batchY,const int max_blocks,const int n1, const int n2){

	void (*kerPtr) (int, const T*, const int, const int, T*, const int, const int, const int*, const int, const int) = launchSetup.kerPtr;

	if(!kerPtr){
            fprintf(stderr, "ERROR: NO CUDA KERNEL FOUND FOR THE GIVEN VALUES at cuZcopy \n");
            exit(-1); // Invalid configuration, abort execution so autotuner does not measure it
         }

	 const int THREADS_PER_BLOCK=launchSetup.TB;


	if (batchesX>max_blocks)
        {
                 fprintf(stderr, "Cuda error cu_zcopy: More blocks requested than available\n");
                 exit(-1);
        }

	dim3 threads (THREADS_PER_BLOCK);
	dim3 blocks (batchesX,batchY);

	kerPtr<<<blocks,threads,0,stream>>>(count,x,offset_x,incx,y,offset_y,incy,offset,n1,n2);
	return 0;
}	
#endif


//WARNING: IT does not support negative increments
//https://netlib.org/blas/zcopy.f
void cuZcopy(const int count, const double * x, const int offset_x, const int incx, double * y, const int offset_y, const int incy, const int * offset, const int batchesX, cudaStream_t stream, const int batchY, const int PLAN, const int max_blocks,const int n1, const int n2){




#if AUTOTUNER

        const int index = getIndex(PLAN+1,atoi(DFTuning::getEnv("QBOX_ZCOPY_TB")), 32 ,atoi(DFTuning::getEnv("QBOX_ZCOPY_TB_SM")));
//TODO: Be careful here!! I am reusing the getIndex function with the PLAN parameter, but if the UNROLL_LIMITS are changed, this will likely fail!!!
        ZcopyLauncher(weightTable_zcopy[index], count,x,offset_x,incx,y,offset_y,incy,offset,batchesX,stream,batchY,max_blocks,n1,n2);
#else	
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


#endif

	if (cudaGetLastError() != cudaSuccess){
                   fprintf(stderr, "Cuda error cu_zcopy: Failed kernel\n");
                   exit(-1);
        }

	return;

}	

/********************************************/

//TODO: templated for T type
#if AUTOTUNER
template<int unroll, int tb, int tb_sm>
__global__ __launch_bounds__(tb,tb_sm) void cu_pairwise(const double* src, double* dest, const int N, const int batch){
#else
template<int unroll>
__global__ void cu_pairwise(const double* src, double* dest, const int N, const int batch){
#endif
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


#if AUTOTUNER

template <typename T>
int PairLauncher(const KernelCfg_pair<T> &launchSetup, const T* src, T* dest, const int N, cudaStream_t stream, const int batch){

	 void (*kerPtr) (const T*, T*, const int, const int) = launchSetup.kerPtr;

	  if(!kerPtr){
            fprintf(stderr, "ERROR: NO CUDA KERNEL FOUND FOR THE GIVEN VALUES at cuPairwise \n");
            exit(-1); // Invalid configuration, abort execution so autotuner does not measure it
         }   

	 const int THREADS_PER_BLOCK=launchSetup.TB;
         const int UNROLL_FACTOR=launchSetup.U;



         //round-up
        const int block_num=(N+THREADS_PER_BLOCK*UNROLL_FACTOR-1)/(THREADS_PER_BLOCK*UNROLL_FACTOR);

        //Check total number of threads/blocks meets architecture requirements
        dim3 threads (THREADS_PER_BLOCK);
        dim3 blocks (block_num,batch);

        kerPtr<<<blocks,threads,0,stream>>>(src, dest,N,batch);
	return 0;
}

#endif



void cuPairwise(const double* src, double* dest, const int N, cudaStream_t stream, const int batch){
  
#if AUTOTUNER
	const int index = getIndex(atoi(DFTuning::getEnv("QBOX_PAIR_unroll")),atoi(DFTuning::getEnv("QBOX_PAIR_TB")), 32 ,atoi(DFTuning::getEnv("QBOX_PAIR_TB_SM")));

	PairLauncher(weightTable_pair[index],src,dest,N,stream,batch);
#else
	constexpr int THREADS_PER_BLOCK=128;
        constexpr int UNROLL_FACTOR=8;

        //round-up
        const int block_num=(N+THREADS_PER_BLOCK*UNROLL_FACTOR-1)/(THREADS_PER_BLOCK*UNROLL_FACTOR);

	//Check total number of threads/blocks meets architecture requirements
        dim3 threads (THREADS_PER_BLOCK);
        dim3 blocks (block_num,batch);

        cu_pairwise<UNROLL_FACTOR><<<blocks,threads,0,stream>>>(src,dest,N,batch);
#endif
	
	if (cudaGetLastError() != cudaSuccess){
                   fprintf(stderr, "Cuda error cu_pairwise: Failed kernel\n");
                   exit(-1);
        }


}




/**************************/



//TODO: templated for T type
#if AUTOTUNER
template<int unroll, int tb, int tb_sm>
__global__ __launch_bounds__(tb,tb_sm) void cu_dscal(const double scal, double* dest, const int N, const int batch){
 #else
template<int unroll>
__global__ void cu_dscal(const double scal, double* dest, const int N, const int batch){
#endif

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



#if AUTOTUNER



template <typename T>
int DscalLauncher (const KernelCfg<T> &launchSetup, const T scal, T* dest, const int N,cudaStream_t stream, const int batch){

	void (*kerPtr) (const T, T*, const int, const int) = launchSetup.kerPtr;

	 if(!kerPtr){
            fprintf(stderr, "ERROR: NO CUDA KERNEL FOUND FOR THE GIVEN VALUES at cuDscal \n");
            exit(-1); // Invalid configuration, abort execution so autotuner does not measure it
         }
	
	 /*if(launchSetup.TB > prop.maxThreadsPerBlock){
           fprintf(stderr,"ERROR: CUDA THREADBLOCK SIZE EXCEEDS LIMITS at cuDscal \n"); exit(-1);
         }
	 */
	 const int THREADS_PER_BLOCK=launchSetup.TB;
	 const int UNROLL_FACTOR=launchSetup.U;


	 //round-up
        const int block_num=(N+THREADS_PER_BLOCK*UNROLL_FACTOR-1)/(THREADS_PER_BLOCK*UNROLL_FACTOR);

        //Check total number of threads/blocks meets architecture requirements
        dim3 threads (THREADS_PER_BLOCK);
        dim3 blocks (block_num,batch);
	
	kerPtr<<<blocks,threads,0,stream>>>(scal, dest,N,batch);
		

 	return 0;
}

 #endif

void cuDscal(const double scal, double * dest, const int N, cudaStream_t stream, const int batch){



#if AUTOTUNER

	 const int index = getIndex(atoi(DFTuning::getEnv("QBOX_DSCAL_unroll")),atoi(DFTuning::getEnv("QBOX_DSCAL_TB")), 32 ,atoi(DFTuning::getEnv("QBOX_DSCAL_TB_SM")));
	
	 //const int index = getIndex(DFTuningQB::ExecutorQBox::getDscalU(), DFTuningQB::ExecutorQBox::getDscalTb(), prop.warpSize, DFTuningQB::ExecutorQBox::getDscalTbSm());
	 DscalLauncher(weightTable[index],scal,dest,N,stream,batch );
#else
	constexpr int THREADS_PER_BLOCK=128;
        constexpr int UNROLL_FACTOR=8;
        //round-up
        const int block_num=(N+THREADS_PER_BLOCK*UNROLL_FACTOR-1)/(THREADS_PER_BLOCK*UNROLL_FACTOR);

        //Check total number of threads/blocks meets architecture requirements
        dim3 threads (THREADS_PER_BLOCK);
        dim3 blocks (block_num,batch);

        cu_dscal<UNROLL_FACTOR><<<blocks,threads,0,stream>>>(scal,dest,N,batch);
#endif


	if (cudaGetLastError() != cudaSuccess){
                   fprintf(stderr, "Cuda error cu_pairwise: Failed kernel\n");
                   exit(-1);
        }
}











