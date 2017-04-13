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
#define MAX_BLK 4//64
#define MAX_THREAD 2048//2048
void alighment_gpu(char* h_seq_a, char* h_seq_b, int seq_size)
{
    char* d_seq_a;
    char* d_seq_b;
    char* d_matrix;
    unsigned long int bytes = sizeof(int) * seq_size;
    unsigned long int matrix_bytes = sizeof(int) * (seq_size+1) * (seq_size+1);

    // d_seq_a = (char*) malloc(bytes);
    // d_seq_b = (char*) malloc(bytes);
    // d_matrix = (char*) malloc(matrix_bytes);
    cudaMalloc((void**)&d_seq_a, (size_t)bytes);printf("I died before this line\n");
    cudaMalloc((void**)&d_seq_b, (size_t)bytes);
    cudaMalloc((void**)&d_matrix, (size_t)matrix_bytes);
    cudaMemcpy((void*)d_seq_a, (void*)h_seq_a, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy((void*)d_seq_b, (void*)h_seq_b, bytes, cudaMemcpyHostToDevice);
    cudaMemset((void*)d_matrix, 0, matrix_bytes);
    //cudaDeviceSynchronize();
    //cudaProfilerStart();
    //do the work
    for (int max_blocks = 1; max_blocks <= MAX_BLK; max_blocks *= 2){//launch with 1, 2, or 4 blocks
		int block_size = seq_size/max_blocks;
		int threads;
		if (seq_size < block_size){
			threads = seq_size;//assign just enough threads for each block
		}else{
			threads = (MAX_THREAD/max_blocks) > 1024 ? 1024 : (MAX_THREAD/max_blocks);//assign max threads
		}



		int max_block_step = 2 * max_blocks -1;
		int mat_width = seq_size + 1;
		//cudaDeviceSynchronize();

		double time1 =  CPUtime();
		for (int block_step = 0; block_step < max_block_step; block_step++){


			int blocks = block_step + 1 < max_blocks ? (block_step + 1) : (max_block_step - block_step);
			dim3 grid(blocks, 1 ,1);
			dim3 block(threads, 1, 1);

			#ifdef LAUNCH
			printf("max_blocks %d, threads %d, block_step %d, block_size %d, blocks %d\n", max_blocks, threads, block_step, block_size, blocks);
			#endif

			if (block_step+1 < max_blocks){
				align2<<< grid, block >>>(d_seq_a, d_seq_b, d_matrix, block_step, block_size, mat_width);
			}else if (block_step+1 == max_blocks){
				align4<<< grid, block >>>(d_seq_a, d_seq_b, d_matrix, block_step, block_size, mat_width, max_blocks);
			}else{
				align3<<< grid, block >>>(d_seq_a, d_seq_b, d_matrix, block_step, block_size, mat_width, max_blocks);
			}

		}
		cudaDeviceSynchronize();

		double time2 = CPUtime();
		printf("GPU took %f blocks: %d\n", time2-time1, max_blocks);

    }
    //cudaDeviceSynchronize();
	//cudaProfilerStop();
    cudaFree((void*) d_seq_a);
    cudaFree((void*) d_seq_b);
    cudaFree((void*) d_matrix);
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
    #ifdef DEBUG
    printf("before longest%d\n", threadIdx.x);
    #endif
    int tid = threadIdx.x;
    int num_threads = blockDim.x;
    int antidiagonal_thread_index, antidiagonal_block_index;
    antidiagonal_block_index = blockIdx.x;
    int max_thread_step = 2 * block_width - 1;

    for (int thread_step = 0; thread_step < max_thread_step; thread_step++){

		antidiagonal_thread_index = tid;
		int thread_step_length = (block_width - thread_step - 1) >= 0 ? (thread_step + 1) : (max_thread_step - thread_step ) ;
		//this is the actual length for each thread step, it decreases after hitting the longest
		int mat_index;
		while (antidiagonal_thread_index < thread_step_length){
			// It's an indexing thing...
			//unsigned long new_index = (unsigned long)j*score_width + i;
			int i , j ;
			if (thread_step < block_width){
			   	j = (1 + antidiagonal_thread_index) * mat_width //thread part
			       	       	+ (antidiagonal_block_index*block_width) * mat_width;//block part
			   	i = (1 + thread_step - antidiagonal_thread_index) //thread
			       	       	+ (block_step - antidiagonal_block_index) * block_width;//block
			}else{
				j = (antidiagonal_thread_index + (thread_step - max_thread_step)+ (block_width + 1)) * mat_width//thread
				    	+ (antidiagonal_block_index*block_width) * mat_width;//block
				i = ((block_width - thread_step) + thread_step - antidiagonal_thread_index) //thread
					+ (block_step - antidiagonal_block_index) * block_width;//block
			}
			int current_score = (seq_a[i-1] == seq_b[j-1]) ? 1 : -1;
			mat_index = j + i;


			#ifdef DEBUG
			printf("thread_step %d, length %d, max %d, thread_id %d, ati %d,i %d, j %d, block_width %d, abi %d\n", thread_step, thread_step_length, max_thread_step, tid,antidiagonal_thread_index, i, j , block_width, antidiagonal_block_index);
			#endif


			matrix[mat_index] = max_gpu(matrix[mat_index - mat_width - 1], //[i-1][j-1]
						matrix[mat_index - 1], // [i-1][j]
						matrix[mat_index - mat_width])
						+ current_score; //[i][j-1]
			//mat_index = mat_index + num_threads * matrix_width - num_threads;
			antidiagonal_thread_index += num_threads;
		}
    // __syncthreads;
    }



}


__global__
void align3(char* seq_a, char* seq_b, char* matrix, int block_step, int block_width, int mat_width, int max_block_step_length){
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
    #ifdef DEBUG
    printf("after longest\n");
    #endif
    int tid = threadIdx.x;
    int num_threads = blockDim.x;
    int antidiagonal_thread_index, antidiagonal_block_index;
    antidiagonal_block_index = blockIdx.x;
    int max_thread_step = 2 * block_width - 1;
    int block_step_length = gridDim.x;
    for (int thread_step = 0; thread_step < max_thread_step; thread_step++){

		antidiagonal_thread_index = tid;
		int thread_step_length = (block_width - thread_step - 1) >= 0 ? (thread_step + 1) : (max_thread_step - thread_step ) ;//this is the actual length for each thread step, it decreases after hitting the longest
		int mat_index;
		while (antidiagonal_thread_index < thread_step_length){
			// It's an indexing thing...
			//unsigned long new_index = (unsigned long)j*score_width + i;

		   int i , j ;
			if (thread_step + 1 <= block_width){
			j = (1 + antidiagonal_thread_index) * mat_width //thread
					+ (antidiagonal_block_index + (max_block_step_length - block_step_length))*block_width * mat_width;//block
			i = (1 + thread_step - antidiagonal_thread_index) //thread
					+ (block_step - antidiagonal_block_index - (max_block_step_length - block_step_length))*block_width;//block
			}else{
			j = (antidiagonal_thread_index + (thread_step - max_thread_step) //thread
					+ (block_width + 1)) * mat_width + (antidiagonal_block_index + (max_block_step_length - block_step_length))*block_width * mat_width;//block
			i = ((block_width - thread_step) + thread_step - antidiagonal_thread_index) //thread
					+ (block_step - antidiagonal_block_index - (max_block_step_length - block_step_length))*block_width;//block
			}
			int current_score = (seq_a[i-1] == seq_b[j-1]) ? 1 : -1;
			mat_index = j + i;


			#ifdef DEBUG
			printf("thread_step %d, length %d, max %d, thread_id %d, ati %d,i %d, j %d, block_width %d, abi %d\n", thread_step, thread_step_length, max_thread_step, tid,antidiagonal_thread_index, i, j , block_width, antidiagonal_block_index);
			#endif




			matrix[mat_index] = max_gpu(matrix[mat_index - mat_width - 1], //[i-1][j-1]
						matrix[mat_index - 1], // [i-1][j]
						matrix[mat_index - mat_width])
						+ current_score; //[i][j-1]
			//mat_index = mat_index + num_threads * matrix_width - num_threads;
			antidiagonal_thread_index += num_threads;
		}
		// __syncthreads;
    }



}


__global__
void align4(char* seq_a, char* seq_b, char* matrix, int block_step, int block_width, int mat_width, int max_block_step_length){
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
    #ifdef DEBUG
    printf("longest %d %d\n", block_step, threadIdx.x);
    #endif
    int tid = threadIdx.x;
    int num_threads = blockDim.x;
    int antidiagonal_thread_index, antidiagonal_block_index;
    antidiagonal_block_index = blockIdx.x;
    int max_thread_step = 2 * block_width - 1;
    //int block_step_length = gridDim.x;
    for (int thread_step = 0; thread_step < max_thread_step; thread_step++){

		antidiagonal_thread_index = tid;
		int thread_step_length = (block_width - thread_step - 1) >= 0 ? (thread_step + 1) : (max_thread_step - thread_step ) ;//this is the actual length for each thread step, it decreases after hitting the longest
		int mat_index;
		while (antidiagonal_thread_index < thread_step_length){

			// It's an indexing thing...
			//unsigned long new_index = (unsigned long)j*score_width + i;


			int i , j ;
			if (thread_step + 1 <= block_width){
			j = (1 + antidiagonal_thread_index) * mat_width + //thread
					(antidiagonal_block_index* block_width) * mat_width;//block
			i = (1 + thread_step - antidiagonal_thread_index) //thread
					+ (block_step - antidiagonal_block_index)*block_width;//block
			}else{
			j = (antidiagonal_thread_index + (thread_step - max_thread_step) //thread
					+ (block_width + 1)) * mat_width + (antidiagonal_block_index* block_width) * mat_width;//block
			i = ((block_width - thread_step) + thread_step - antidiagonal_thread_index) //thread
					+ (block_step - antidiagonal_block_index)*block_width;//block
			}
			int current_score = (seq_a[i-1] == seq_b[j-1]) ? 1 : -1;
			mat_index = j + i;


			#ifdef DEBUG
			printf("thread_step %d, length %d, max %d, thread_id %d, ati %d,i %d, j %d, block_width %d, abi %d\n", thread_step, thread_step_length, max_thread_step, tid,antidiagonal_thread_index, i, j , block_width, antidiagonal_block_index);
			#endif







			matrix[mat_index] = max_gpu(matrix[mat_index - mat_width - 1], //[i-1][j-1]
						matrix[mat_index - 1], // [i-1][j]
						matrix[mat_index - mat_width])
						+ current_score; //[i][j-1]
			//mat_index = mat_index + num_threads * matrix_width - num_threads;
			antidiagonal_thread_index += num_threads;
		}
		// __syncthreads;
    }



}

__device__
static long max_gpu(long a, long b, long c)
{
  long result = a;

  if(b > result) {
    result = b;
  }
  if(c > result) {
    result = c;
  }

  return result;
}
