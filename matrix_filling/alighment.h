#define MAX_RUN 100

void alighment_cpu (char* seq_a, char* seq_b);
void generateSequence(char* seq_a, char* seq_b, int size);

void alighment_gpu(char* h_seq_a, char* h_seq_b, int seq_size, double cpu_time);

double CPUtime();

