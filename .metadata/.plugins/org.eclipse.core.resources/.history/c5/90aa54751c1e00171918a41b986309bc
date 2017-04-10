#include "alighment.h"
__host__
void alignment_gpu(char* h_seq_a, char* h_seq_b, int size)
{
    char* d_seq_a, d_seq_b, d_matrix;
    unsigned long int bytes = sizeof(int) * size;
    unsigned long int matrix_bytes = sizeof(int) * (size+1) * (size+1);
    cudaMalloc((void**)d_seq_a, (size_t)bytes);
    cudaMalloc((void**)d_seq_b, (size_t)bytes);
    cudaMalloc((void**)d_matrix, (size_t)matrix_bytes)
	cudaMemcpy(d_seq_a, h_seq_a, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_seq_b, h_seq_b, bytes, cudaMemcpyHostToDevice);
    cudaMemset(d_matrix, \0, matrix_bytes);
    
    
    //do the work
    for (int blocks = 1; blocks <= 16; blocks *= 4){
        
        if (size < (2048/blocks)){
            threads = size;
        }else{
            threads = 2048/blocks;
        }

        dim3 grid(blocks, 1 ,1);
        dim3 block(threads, 1, 1);
        
        double time1 =  CPUtime();
        align<<< grid, block >>>(d_seq_a, d_seq_b, d_matrix, threads, blocks, size);
        double time2 = CPUtime();
        print("GPU reduction took %lu\n", time2-time1);

    }
    
    cudaFree((void*) d_seq_a);
    cudaFree((void*) d_seq_b);
    cudaFree((void*) d_matrix)
}

__device__
void align(char* seq_a, char* seq_b, char* matrix, int threads, int griddim){

    int tid = threadIdx.x;
    https://github.com/NVlabs/nvbio/blob/master/nvbio/alignment/sw/sw_warp_inl.h
    see above

}








