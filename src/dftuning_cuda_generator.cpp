#include <iostream>
#include <string>

#include <cuda_runtime.h>
#include <cuda.h>
//const static std::string nameKernel("cu_dscal");


constexpr static int STRIDE_UNROLL=2;
constexpr static int STARTING_UNROLL=1;
constexpr static int MAX_UNROLL=8;


constexpr static int STRIDE_TB=32;
constexpr static int STARTING_TB=32;


constexpr static int STRIDE_TB_SM=1;
constexpr static int STARTING_TB_SM=1;


void print_row(const std::string nameKernel,const int unroll,const int tb , const int tb_sm, const int max_threads_sm)
{
      const bool showNULL=(tb_sm*tb>max_threads_sm);
      if(showNULL)
              std::cout<<"NULL_ROW(),";
      else
              std::cout<<"ROW("<<nameKernel<<","<<unroll<<","<<tb<<","<<tb_sm<<"),";
      std::cout<< "\n";

}








int main(int argc, char* arvg[])
{
        constexpr int dev_id=0;

        int max_tb_sm=1;
        cudaDeviceGetAttribute (&max_tb_sm,cudaDevAttrMaxBlocksPerMultiprocessor,dev_id);
	int max_tb=1024;
	cudaDeviceGetAttribute (&max_tb, cudaDevAttrMaxThreadsPerBlock,dev_id);
/* // Copy this in the generated header file 

template<typename T>
struct KernelCfg{
        void (*kerPtr)(const T, T*, const int, const int ); //Kernel pointer
	
	short U; //Unrolling
        short TB; //Threadblock Size
        short TB_SM;    //Threadblocks per SM
};


template<typename T>
struct KernelCfg_pair{
        void (*kerPtr)(const T*, T*, const int, const int ); //Kernel pointer
        
        short U; //Unrolling
        short TB; //Threadblock Size
        short TB_SM;    //Threadblocks per SM
};



template<typename T>
struct KernelCfg_zcopy{
        void (*kerPtr)(const int, const T*,const int, const int, T*, const int, const int, const int *, const int, const int ); //Kernel pointer
	short P; //Plan        
        short TB; //Threadblock Size
        short TB_SM;    //Threadblocks per SM
};


template<typename T>
struct KernelCfg_vec{
        void (*kerPtr)( const T*, T*, const int*,const int*, const int, const int); //Kernel pointer
        short P; //Plan        
        short TB; //Threadblock Size
        short TB_SM;    //Threadblocks per SM
};


template<typename T>
struct KernelCfg_zvec{
        void (*kerPtr)( const T*, T*, const int*, const int, const int); //Kernel pointer
        short U; //Plan        
        short TB; //Threadblock Size
        short TB_SM;    //Threadblocks per SM
};



*/


	
        int max_threads_sm=1536;
        cudaDeviceGetAttribute(&max_threads_sm,cudaDevAttrMaxThreadsPerMultiProcessor,dev_id);
        std::cout << "#define NULL_ROW() {NULL,0,0,0} \n";
        std::cout << "#define ROW(_Kernel, _U, _TB, _TB_SM) {_Kernel<_U, _TB, _TB_SM>,_U,  _TB, _TB_SM }\n";

        std::cout << "\n\n";
        std::cout << "constexpr static int STRIDE_UNROLL="<<STRIDE_UNROLL <<";\n";
        std::cout << "constexpr static int STARTING_UNROLL="<< STARTING_UNROLL<<";\n";
        std::cout << "constexpr static int MAX_UNROLL="<<MAX_UNROLL <<";\n";
        
	std::cout << "constexpr static int STRIDE_TB="<< STRIDE_TB<<";\n";
        std::cout << "constexpr static int STARTING_TB="<< STARTING_TB<<";\n";
        std::cout << "constexpr static int MAX_TB="<< max_tb<<";\n";
        
	std::cout << "constexpr static int STRIDE_TB_SM="<< STRIDE_TB_SM<<";\n";
        std::cout << "constexpr static int STARTING_TB_SM="<< STARTING_TB_SM <<";\n";
        std::cout << "constexpr static int MAX_TB_SM="<<max_tb_sm<<";\n\n";
/****/
        std::cout << "const static KernelCfg<double> weightTable []={"<<"\n";

        #pragma unroll
        for(int i=STARTING_UNROLL; i<=MAX_UNROLL;i*=STRIDE_UNROLL){
                for(int j=STARTING_TB; j<=max_tb; j+=STRIDE_TB){
                        for (int k=STARTING_TB_SM; k<=max_tb_sm; k+=STRIDE_TB_SM){
                                print_row("cu_dscal",i,j,k,max_threads_sm);
                        }
                }
        }
        std::cout<<"NULL_ROW()\n };\n";
/****/
	std::cout << "const static KernelCfg_pair<double> weightTable_pair []={"<<"\n";

        #pragma unroll
        for(int i=STARTING_UNROLL; i<=MAX_UNROLL;i*=STRIDE_UNROLL){
                for(int j=STARTING_TB; j<=max_tb; j+=STRIDE_TB){
                        for (int k=STARTING_TB_SM; k<=max_tb_sm; k+=STRIDE_TB_SM){
                                print_row("cu_pairwise",i,j,k,max_threads_sm);
                        }
                }
        }
        std::cout<<"NULL_ROW()\n };\n";

/****/
        std::cout << "const static KernelCfg_zcopy<double> weightTable_zcopy []={"<<"\n";

 	#pragma unroll
        for(int i=-1; i<=1;i+=2){
                for(int j=STARTING_TB; j<=max_tb; j+=STRIDE_TB){
                        for (int k=STARTING_TB_SM; k<=max_tb_sm; k+=STRIDE_TB_SM){
                                print_row("cu_Z_copy",i,j,k,max_threads_sm);
                        }
                }
        }
        std::cout<<"NULL_ROW()\n };\n";	
/****/

	std::cout << "const static KernelCfg_zvec<double> weightTable_zvec []={"<<"\n";

        #pragma unroll
        for(int i=STARTING_UNROLL; i<=MAX_UNROLL;i*=STRIDE_UNROLL){
                for(int j=STARTING_TB; j<=max_tb; j+=STRIDE_TB){
                        for (int k=STARTING_TB_SM; k<=max_tb_sm; k+=STRIDE_TB_SM){
                                print_row("zvec_to_vector_kernel",i,j,k,max_threads_sm);
                        }
                }
        }
        std::cout<<"NULL_ROW()\n };\n";
/****/

std::cout << "const static KernelCfg_vec<double> weightTable_vec []={"<<"\n";

        #pragma unroll
        for(int i=0; i<=1;i+=1){
                for(int j=STARTING_TB; j<=max_tb; j+=STRIDE_TB){
                        for (int k=STARTING_TB_SM; k<=max_tb_sm; k+=STRIDE_TB_SM){
                                print_row("vector_to_zvec_kernel",i,j,k,max_threads_sm);
                        }
                }
        }
        std::cout<<"NULL_ROW()\n };\n";

	
}


