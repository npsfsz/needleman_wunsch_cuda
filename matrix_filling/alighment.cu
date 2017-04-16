#include "alighment.cuh"
extern "C" {
    #include "alighment.h"
}

#include "stdlib.h"
#include "stdio.h"
#include <cuda.h>
#include <cuda_profiler_api.h>

//#include "alighment.c"
//#define LAUNCH
//#define DEBUG
//#define HEAD
#define MAX_BLK 32//128
#define MAX_THREAD 2048//2048
#define MAX_THREAD_PER_BLK 1024
void alighment_gpu(char* h_seq_a, char* h_seq_b, int seq_size, double cpu_time)
{

    printf("Kernel warm up\n");
    char* d_seq_a;
    char* d_seq_b;
    char* d_matrix;
    unsigned long int bytes = sizeof(char) * seq_size;
    unsigned long int matrix_bytes = sizeof(int) * (seq_size+1) * (seq_size+1);

	cudaMalloc((void**)&d_seq_a, (size_t)bytes);
    cudaMalloc((void**)&d_seq_b, (size_t)bytes);
    cudaMalloc((void**)&d_matrix, (size_t)matrix_bytes);
    cudaMemcpy((void*)d_seq_a, (const void*)h_seq_a, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy((void*)d_seq_b, (const void*)h_seq_b, bytes, cudaMemcpyHostToDevice);
    //cudaMemset((void*)d_matrix, 0, matrix_bytes);
    int* h_matrix = (int*) malloc(matrix_bytes);
    h_matrix[0] = 0;
        // int mismatch = -1;
        // int match = 1;
    int indel = -1;
    int i, j;
      // work along first row -> [i][0]
    for(i = 1; i < seq_size+1; i++)
    {
      h_matrix[i] = i * indel;
    }

      // work down first column -> [0][j]
    for(j = 1; j < seq_size+1; j++)
    {
      unsigned long index = (unsigned long)j*seq_size+1;
      h_matrix[index] =  j * indel;
    }
    cudaMemcpy((void*)d_matrix, (const void*)h_matrix, matrix_bytes, cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
    cudaFree((void*) d_seq_a);
    cudaFree((void*) d_seq_b);
    cudaFree((void*) d_matrix);
    cudaDeviceSynchronize();
    //cudaProfilerStart();
    //do the workprintf("1\n");
    int block_limit;
    if (seq_size == 4096){
        block_limit = 4;
    }else{
        block_limit = 1;
    }
    for (int max_blocks = block_limit; max_blocks <= MAX_BLK; max_blocks *= 2){//launch with 1, 2, or 4 blocks
		int block_size = seq_size/max_blocks;
		int threads;
		/*
		if (seq_size < 1024){
			threads = seq_size;//assign just enough threads for each block
		}else{
			threads = (MAX_THREAD/max_blocks) > 1024 ? 1024 : (MAX_THREAD/max_blocks);//assign max threads
		}*/
        //printf("%d\n", threads);
		    char* d_seq_a = 0;
            char* d_seq_b = 0;
            char* d_matrix = 0;
            unsigned long int bytes = sizeof(char) * seq_size;
            unsigned long int matrix_bytes = sizeof(int) * (seq_size+1) * (seq_size+1);
  
	        cudaMalloc((void**)&d_seq_a, (size_t)bytes);
            cudaMalloc((void**)&d_seq_b, (size_t)bytes);
            cudaMalloc((void**)&d_matrix, (size_t)matrix_bytes);
            cudaMemcpy((void*)d_seq_a, (const void*)h_seq_a, bytes, cudaMemcpyHostToDevice);
            cudaMemcpy((void*)d_seq_b, (const void*)h_seq_b, bytes, cudaMemcpyHostToDevice);
            //cudaMemset((void*)d_matrix, 0, matrix_bytes);
            int* h_matrix = (int*) malloc(matrix_bytes);
            h_matrix[0] = 0;
                // int mismatch = -1;
                // int match = 1;
            int indel = -1;
            int i, j;
              // work along first row -> [i][0]
            for(i = 1; i < seq_size+1; i++)
            {
              h_matrix[i] = i * indel;
            }

              // work down first column -> [0][j]
            for(j = 1; j < seq_size+1; j++)
            {
              unsigned long index = (unsigned long)j*seq_size+1;
              h_matrix[index] =  j * indel;
            }
            cudaMemcpy((void*)d_matrix, (const void*)h_matrix, matrix_bytes, cudaMemcpyHostToDevice);
            cudaDeviceSynchronize();
 



		int max_block_step = 2 * max_blocks -1;
		int mat_width = seq_size + 1;
		//cudaDeviceSynchronize();
        double avg_gpu_time = 0;
        double total_gpu_time = 0;
		for (int run = 0; run < MAX_RUN; run++){
            //cudaMemset((void*)d_matrix, 0, matrix_bytes);
			double time1 =  CPUtime();
			
		    for (int block_step = 0; block_step < max_block_step; block_step++){


			    int blocks = block_step + 1 < max_blocks ? (block_step + 1) : (max_block_step - block_step);
			    
			    if ((blocks*block_size) > MAX_THREAD ){
			    
			        threads = MAX_THREAD / blocks ;
			        
			        
			    }else{

			        threads = (block_size);

			    }

			    if (threads > MAX_THREAD_PER_BLK){
			        threads = MAX_THREAD_PER_BLK;
			    }
                threads = (threads + 16) & (0xffffffff ^ 31);
                
			
			
			    dim3 grid(blocks, 1 ,1);
			    dim3 block(threads, 1, 1);
                int shared_mem_size = (2*block_size) * sizeof(int); 
			    #ifdef LAUNCH
			    printf("block_step %d, max_blocks %d, threads %d, total threads %d,block_size %d, blocks %d, mem %d\n", block_step,max_blocks, threads, threads*blocks ,  block_size, blocks, shared_mem_size);
			    #endif

			    if (block_step+1 < max_blocks){
			        //printf("1\n");
				    align2<<< grid, block, shared_mem_size >>>(d_seq_a, d_seq_b, d_matrix, block_step, block_size, mat_width);
			    }else if (block_step+1 == max_blocks){
			        //printf("2\n");
				    align4<<< grid, block, shared_mem_size >>>(d_seq_a, d_seq_b, d_matrix, block_step, block_size, mat_width, max_blocks, blocks);
			    }else{
			        //printf("3\n");
				    align3<<< grid, block, shared_mem_size >>>(d_seq_a, d_seq_b, d_matrix, block_step, block_size, mat_width, max_blocks, blocks);
			    }
		    
		    }
cudaDeviceSynchronize();
		    double time2 = CPUtime();
		    total_gpu_time += time2-time1;
		    avg_gpu_time = total_gpu_time / MAX_RUN;

        }

		printf("GPU took %f blocks: %d, speedup %f\n", avg_gpu_time, max_blocks, cpu_time/avg_gpu_time);


			    cudaDeviceSynchronize();
    cudaFree((void*) d_seq_a);
    cudaFree((void*) d_seq_b);
    cudaFree((void*) d_matrix);
		    cudaDeviceSynchronize();
    }
    //cudaDeviceSynchronize();
	//cudaProfilerStop();

    
}


__global__
void align2(char* seq_a, char* seq_b, char* matrix, int block_step, int block_width, int mat_width){
    //this is the alignment function before the longest block anitdiagonal

    /****************
    grid size 1,2,4 largest antidiagonal length 1,2,4
    for each block step(1, 121, 1234321){
	at each block step, let certain block work
	for each thread step(0 to 2*blockwidth -1 ){
	    index = tid
	    while (index + numofthreads < stepsize){
		matrix[index] = max( , , )
		index += numofthreads
	    }
	{
    }
     ****************/
    #ifdef HEAD
    printf("before longest%d\n", threadIdx.x);
    #endif
    int tid = threadIdx.x;
    int num_threads = blockDim.x;
    int antidiagonal_thread_index, antidiagonal_block_index;
    antidiagonal_block_index = blockIdx.x;
    int max_thread_step = 2 * block_width - 1;
    int i, j, mat_index, current_score;
    int thread_step_length;
    extern __shared__ int s[];
    int top_index, left_index;
    //printf("%d\n", max_thread_step);
    

    
    for (int thread_step = 0; thread_step < max_thread_step; thread_step++){

		antidiagonal_thread_index = tid;
		thread_step_length = (block_width - thread_step - 1) >= 0 ? (thread_step + 1) : (max_thread_step - thread_step ) ;
		//this is the actual length for each thread step, it decreases after hitting the longest
//printf("cao ati: %d, tsl %d\n", antidiagonal_thread_index, thread_step_length);
		while (antidiagonal_thread_index < thread_step_length){
			// It's an indexing thing...
			//unsigned long new_index = (unsigned long)j*score_width + i;

			if (thread_step < block_width){
			   	j = (1 + antidiagonal_thread_index) * mat_width //thread part
			       	       	+ (antidiagonal_block_index*block_width) * mat_width;//block part
			   	i = (1 + thread_step - antidiagonal_thread_index) //thread
			       	       	+ (block_step - antidiagonal_block_index) * block_width;//block
			}else{
				j = /*(antidiagonal_thread_index + (thread_step - max_thread_step)+ (block_width + 1)) * mat_width*///thread
				    (antidiagonal_thread_index + thread_step - block_width + 2) * mat_width
				    	+ (antidiagonal_block_index*block_width) * mat_width;//block
				i = (block_width - antidiagonal_thread_index) //thread
					+ (block_step - antidiagonal_block_index) * block_width;//block
			}
			current_score = (seq_a[i-1] == seq_b[j-mat_width]) ? 1 : -1;
			mat_index = j + i;


			#ifdef DEBUG
			printf("thread_step %d, length %d, max %d, thread_id %d, ati %d,i %d, j %d, block_width %d, abi %d\n", thread_step, thread_step_length, max_thread_step, tid,antidiagonal_thread_index, i, j , block_width, antidiagonal_block_index);
			#endif

            /*
			matrix[mat_index] = max(
			                        max(
			                            matrix[mat_index - mat_width], 
			                            matrix[mat_index - 1]
			                        ),
			                        matrix[mat_index - mat_width - 1])
			                        +current_score;
			 */
			s[antidiagonal_thread_index] = j - mat_width + i;
			matrix[mat_index] = max(
			                        max(
			                            s[antidiagonal_thread_index], 
			                            s[antidiagonal_thread_index+1]
			                        ),
			                        matrix[mat_index - mat_width - 1])
			                        +current_score;
			//mat_index = mat_index + num_threads * matrix_width - num_threads;
			antidiagonal_thread_index += num_threads;
		}
    // __syncthreads;
    }



}


__global__
void align3(char* seq_a, char* seq_b, char* matrix, int block_step, int block_width, int mat_width, int max_blocks_num, int block_step_length){
    //this is the alignment function after the longest block anitdiagonal

    /****************
    grid size 1,2,4 largest antidiagonal length 1,2,4
    for each block step(1, 121, 1234321){
	at each block step, let certain block work
	for each thread step(0 to 2*blockwidth -1 ){
	    index = tid
	    while (index + numofthreads < stepsize){
		matrix[index] = max( , , )
		index += numofthreads
	    }
	{
    }
     ****************/
    #ifdef HEAD
    printf("after longest\n");
    #endif
    int tid = threadIdx.x;
    int max_block_step_length = max_blocks_num;
    int num_threads = blockDim.x;
    int antidiagonal_thread_index, antidiagonal_block_index;
    antidiagonal_block_index = blockIdx.x;
    int max_thread_step = 2 * block_width - 1;
    extern __shared__ int s[];

    
    int mat_index, i, j, thread_step_length, current_score;
    
    for (int thread_step = 0; thread_step < max_thread_step; thread_step++){

		antidiagonal_thread_index = tid;
		thread_step_length = (block_width - thread_step - 1) >= 0 ? (thread_step + 1) : (max_thread_step - thread_step ) ;//this is the actual length for each thread step, it decreases after hitting the longest

		while (antidiagonal_thread_index < thread_step_length){
			// It's an indexing thing...
			//unsigned long new_index = (unsigned long)j*score_width + i;

		  
			if (thread_step  < block_width){
			    j = (1 + antidiagonal_thread_index) * mat_width //thread
					    + (antidiagonal_block_index + max_block_step_length - block_step_length)*block_width * mat_width;//block
			    i = (1 + thread_step - antidiagonal_thread_index) //thread
					    + (block_step - antidiagonal_block_index - (max_block_step_length - block_step_length))*block_width;//block
			}else{
			    j = (antidiagonal_thread_index + (thread_step - max_thread_step) + (block_width + 1)) * mat_width//thread
			        			       + (antidiagonal_block_index + max_block_step_length - block_step_length)*block_width * mat_width;//block
			    i = (block_width - antidiagonal_thread_index) //thread
					    + (block_step - antidiagonal_block_index - (max_block_step_length - block_step_length))*block_width;//block
			}
			current_score = (seq_a[i-1] == seq_b[j-mat_index]) ? 1 : -1;
			mat_index = j + i;


			#ifdef DEBUG
			printf("thread_step %d, length %d, max %d, thread_id %d, ati %d,i %d, j %d, block_width %d, abi %d\n", thread_step, thread_step_length, max_thread_step, tid,antidiagonal_thread_index, i, j , block_width, antidiagonal_block_index);
			#endif




            /*
			matrix[mat_index] = max(
			                        max(
			                            matrix[mat_index - mat_width], 
			                            matrix[mat_index - 1]
			                        ),
			                        matrix[mat_index - mat_width - 1])
			                        +current_score;
			 */
			s[antidiagonal_thread_index] = j - mat_width + i;
			matrix[mat_index] = max(
			                        max(
			                            s[antidiagonal_thread_index], 
			                            s[antidiagonal_thread_index+1]
			                        ),
			                        matrix[mat_index - mat_width - 1])
			                        +current_score;
			//mat_index = mat_index + num_threads * matrix_width - num_threads;
			antidiagonal_thread_index += num_threads;
		}
		// __syncthreads;
    }



}


__global__
void align4(char* seq_a, char* seq_b, char* matrix, int block_step, int block_width, int mat_width, int max_blocks_num, int block_step_length){
    //this is the alignment function for the longest block anitdiagonal

    /****************
    grid size 1,2,4 largest antidiagonal length 1,2,4
    for each block step(1, 121, 1234321){
	at each block step, let certain block work
	for each thread step(0 to 2*blockwidth -1 ){
	    index = tid
	    while (index + numofthreads < stepsize){
		matrix[index] = max( , , )
		index += numofthreads
	    }
	{
    }
     ****************/
    #ifdef HEAD
    printf("longest %d %d\n", block_step, threadIdx.x);
    #endif
    int tid = threadIdx.x;
    int max_block_step_length = gridDim.x;
    int num_threads = blockDim.x;
    int antidiagonal_thread_index, antidiagonal_block_index;
    antidiagonal_block_index = blockIdx.x;
    int max_thread_step = 2 * block_width - 1;
    //int block_step_length = gridDim.x;
    int thread_step_length, mat_index, i, j, current_score;
    extern __shared__ int s[];
    for (int thread_step = 0; thread_step < max_thread_step; thread_step++){

		antidiagonal_thread_index = tid;
		thread_step_length = (block_width - thread_step - 1) >= 0 ? (thread_step + 1) : (max_thread_step - thread_step ) ;//this is the actual length for each thread step, it decreases after hitting the longest
		
		while (antidiagonal_thread_index < thread_step_length){

			// It's an indexing thing...
			//unsigned long new_index = (unsigned long)j*score_width + i;

			if (thread_step  < block_width){
			j = (1 + antidiagonal_thread_index) * mat_width + //thread
					(antidiagonal_block_index* block_width) * mat_width;//block
			i = (1 + thread_step - antidiagonal_thread_index) //thread
					+ (block_step - antidiagonal_block_index)*block_width;//block
			}else{
			j = (antidiagonal_thread_index + (thread_step - max_thread_step) + (block_width + 1)) * mat_width//thread
			    			       	+ (antidiagonal_block_index* block_width) * mat_width;//block
			i = (block_width  - antidiagonal_thread_index) //thread
					+ (block_step - antidiagonal_block_index)*block_width;//block
			}
			current_score = (seq_a[i-1] == seq_b[j-mat_width]) ? 1 : -1;
			mat_index = j + i;


			#ifdef DEBUG
			printf("thread_step %d, length %d, max %d, thread_id %d, ati %d,i %d, j %d, block_width %d, abi %d\n", thread_step, thread_step_length, max_thread_step, tid,antidiagonal_thread_index, i, j , block_width, antidiagonal_block_index);
			#endif







            /*
			matrix[mat_index] = max(
			                        max(
			                            matrix[mat_index - mat_width], 
			                            matrix[mat_index - 1]
			                        ),
			                        matrix[mat_index - mat_width - 1])
			                        +current_score;
			 */
			s[antidiagonal_thread_index] = j - mat_width + i;
			matrix[mat_index] = max(
			                        max(
			                            s[antidiagonal_thread_index], 
			                            s[antidiagonal_thread_index+1]
			                        ),
			                        matrix[mat_index - mat_width - 1])
			                        +current_score;
			//mat_index = mat_index + num_threads * matrix_width - num_threads;
			antidiagonal_thread_index += num_threads;
		}
		// __syncthreads;
    }



}


