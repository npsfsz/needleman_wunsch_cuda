#define DEBUG 1
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <limits.h>
#include "alighment.h"
#include <cuda.h>
#define MIN_SEQ_LEN 8//512
#define MAX_SEQ_LEN 16384//16384


static long max(long a, long b, long c)
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

void alighment_cpu (char* seq_a, char* seq_b)
{

  size_t length_a = strlen(seq_a);
  size_t length_b = strlen(seq_b);
  
  // Calculate largest amount of mem needed
//   unsigned int longest_alignment = (unsigned int)length_a +
//                                   (unsigned int)length_b;
  
  unsigned int score_width = length_a+1;
  unsigned int score_height = length_b+1;
  
  // Check sequences aren't too long to align
  if(score_height > ULONG_MAX / score_width)
  {
    fprintf(stderr, "Error: sequences too long to align (%i * %i > %lu)\n",
            score_width, score_height, ULONG_MAX);
    exit(EXIT_FAILURE);
  }
  
  unsigned long arr_size = (unsigned long)score_width * score_height;
  //printf("%d", arr_size);
  // 2d array (length_a x length_b)
  // addressing [a][b]

  // Score having just matched
  int* match_score = (int*) malloc(arr_size * sizeof(int));

  if(match_score == NULL )
  {
    unsigned long num_of_bytes = arr_size * sizeof(int);
    fprintf(stderr, "Couldn't allocate enough memory (%lu bytes required)\n",
            num_of_bytes);
    
    #ifdef DEBUG
    fprintf(stderr, "SeqA length: %i; SeqB length: %i\n", (int)length_a, (int)length_b);
    fprintf(stderr, "arr_size: %lu; int size: %li\n", arr_size, sizeof(int));
    #endif

    exit(EXIT_FAILURE);
  }
  
  #ifdef DEBUG
  printf("Malloc'd score matrices! %lu bytes\n", arr_size*4);
  #endif

  // Fill in traceback matrix

  unsigned int i, j;
  
  // [0][0]
  match_score[0] = 0;
    // int mismatch = -1;
    // int match = 1;
    int indel = -1;
  
  // work along first row -> [i][0]
  for(i = 1; i < score_width; i++)
  {
    match_score[i] = i * indel;
    
  }
  
  // work down first column -> [0][j]
  for(j = 1; j < score_height; j++)
  {
    unsigned long index = (unsigned long)j*score_width;
    match_score[index] =  j * indel;
  }
  
  //
  // Update Dynamic Programming arrays
  //
  for (i = 1; i < score_width; i++)
  {
    for (j = 1; j < score_height; j++)
    {
      // It's an indexing thing...
      unsigned int seq_i = i - 1;
      unsigned int seq_j = j - 1;
      
      // Update match_score[i][j] with position [i-1][j-1]
      // Addressing array must be done with unsigned long
      unsigned long new_index = (unsigned long)j*score_width + i;
      unsigned long old_index1 = (unsigned long)(j-1)*score_width + (i-1);
      // Update gap_a_score[i][j] from position [i][j-1]
      unsigned long old_index2 = (unsigned long)(j-1)*score_width + i;
      // Update gap_b_score[i][j] from position [i-1][j]
      unsigned long old_index3 = (unsigned long)j*score_width + (i-1);
      int current_score = (seq_a[seq_i] == seq_b[seq_j]) ? 1 : -1;
                                     
      //max(1,2,3);
      
      match_score[new_index] = max((long)match_score[old_index1],
                                   (long)match_score[old_index2],
                                   (long)match_score[old_index3])
                               + current_score;
                               
    }
  }
}


void generateSequence(char* seq_a, char* seq_b, int size){
    int i;
    for (i = 0; i < size; i++){
        srand(time(NULL));
        int random = rand() % 5;
        switch(random) {
            case 0:
                //G
                seq_a[i] = 'G';
                break;
            case 1:
                //A
                seq_a[i] = 'A';
                break;
            case 2:
                //T
                seq_a[i] = 'T';
                break;
            case 3:
                //C
                seq_a[i] = 'C';
                break;
            case 4:
                //U
                seq_a[i] = 'U';
                break;
            default:
                printf("ERROR\n");
        }

        random = rand() % 5;
        switch(random) {
            case 0:
                //G
                seq_b[i] = 'G';
                break;
            case 1:
                //A
                seq_b[i] = 'A';
                break;
            case 2:
                //T
                seq_b[i] = 'T';
                break;
            case 3:
                //C
                seq_b[i] = 'C';
                break;
            case 4:
                //U
                seq_b[i] = 'U';
                break;
            default:
                printf("ERROR\n");
        }
    }
}

double CPUtime(){
       struct timeval tp;
       gettimeofday (&tp, NULL);
       return ((double)tp.tv_sec + (double)tp.tv_usec * 1e-6);
}

int main(void){
    int size;
    for (size = MIN_SEQ_LEN; size <= MAX_SEQ_LEN; size *= 2){
        //generate sequences
        printf("size is %d\n", size);
        char* h_seq_a = (char*) malloc(sizeof(char) * size);
        char* h_seq_b = (char*) malloc(sizeof(char) * size);
        generateSequence(h_seq_a, h_seq_b, size);

        //time()
        double time1 = CPUtime();
        //runing alignment on CPU
        alighment_cpu(h_seq_a, h_seq_b);
        //time()
        double time2 = CPUtime();
        printf("CPU alightment took %f\n", time2 - time1);    



        //run alignment on GPU
        alighment_gpu(h_seq_a, h_seq_b, size);

        printf("***************************************\n");
    }
}

