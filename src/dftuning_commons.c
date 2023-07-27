#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>




const char* getEnv(const char* env){
        char* value= getenv(env);
        if(value)
                return value;
        printf("NO ENV VARIABLE FOUND FOR %s \n",env);
        exit(-1);
}



//extern "C"
//{
void dftuning_read_configuration (int32_t niters, uint32_t nbatches, int32_t batch_index)
{

   char buffer[1024];
   char *record,*line;
   int32_t i=0,j=0;
   char mat[100][100][100];
   char filename[500];
   // TODO: READ THE INPUT PATH
   snprintf(filename,500,"%s/%d.input",getEnv("DFTUNING_APP_INPUT_TUNING"),batch_index);
   FILE *fstream = fopen(filename,"r");
   
   //General statement for this implementation
   if(nbatches >100)
   {
	   printf("\n DFTUNING: Maximum number of rows per configuration file must be less than 100 configurations");
	   exit(-1);
   }
   if(fstream == NULL)
   {
      printf("\n DFTUNING: configuration file opening failed %s",filename);
      exit(-1);
   }
   while((line=fgets(buffer,sizeof(buffer),fstream))!=NULL)
   {
     record = strtok(line,";");
     while(record != NULL)
     {
     strcpy(&mat[i][j++][0], record);//atoi(record) ;
     record = strtok(NULL,";");
     }
     ++i ;
   }
   fclose(fstream);
   //General check when working with many batches
   if(i!=nbatches)
   {
	 fprintf(stderr,"\n DFTUNING: The number of lines in configuration file does not match with number of batches provided");
	  exit(-1); 
   }

  
   //Specific statement for this implementation. Can be extended in the future.
   if(nbatches>1)
   {	
	fprintf(stderr,"\n In this code, we will use only 1 batch as we cannot fully use the C++ class hierarchy deve    loped in the framework for executing multiple performance configurations in a single app call");
	exit(-1);
   }

   printf("QBOX_DSCAL_unroll %s\n",mat[0][1]);
   setenv("QBOX_DSCAL_unroll",mat[0][1],1);   

   printf("QBOX_DSCAL_TB %s\n",mat[0][2]);
   setenv("QBOX_DSCAL_TB",mat[0][2],1);

   printf("QBOX_DSCAL_TB_SM %s\n",mat[0][3]);
   setenv("QBOX_DSCAL_TB_SM",mat[0][3],1);

}

 void dftuning_set_output( int32_t niters, int32_t nbatches, int32_t batch_index, float time)
{
	char filename[500];
 	snprintf(filename,500,"%s/%d.output",getEnv("DFTUNING_APP_OUTPUT"),batch_index);
 	FILE *fstream = fopen(filename,"w");
	if(fstream==NULL)
	{
		fprintf(stderr,"\n DFTUNING: configuration file opening failed %s",filename);
      		exit(-1);
	}
	fprintf(fstream,"QBox.Time: %lf\n",time/niters);
	fclose(fstream);
}	

//} // end extern C


