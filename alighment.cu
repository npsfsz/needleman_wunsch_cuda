#include "alighment.h"
__host__
void alighment_gpu(char* h_seq_a, char* h_seq_b, int size)
{
    char* d_seq_a, d_seq_b;
    unsigned long int bytes = sizeof(int) * size;
    cudaMalloc((void**)d_seq_a, (size_t)bytes);
    cudaMalloc((void**)d_seq_b, (size_t)bytes);
	cudaMemcpy(d_seq_a, h_seq_a, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_seq_b, h_seq_b, bytes, cudaMemcpyHostToDevice);
    
    
    //do the work
    
    
    cudaFree((void*) d_seq_a);
    cudaFree((void*) d_seq_b);
}
