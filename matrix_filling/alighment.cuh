//#include "alighment.h"
//__device__ void align(char* seq_a, char* seq_b, char* matrix, int block_step, int block_width, int mat_width);
__global__ void align2(char* seq_a, char* seq_b, char* matrix, int block_step, int block_width, int mat_width);
__global__ void align3(char* seq_a, char* seq_b, char* matrix, int block_step, int block_width, int mat_width, int max_blocks_num, int block_step_length);
__global__ void align4(char* seq_a, char* seq_b, char* matrix, int block_step, int block_width, int mat_width, int max_blocks_num, int block_step_length);
__device__ static long max_gpu(long a, long b, long c);
