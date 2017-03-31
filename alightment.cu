#define DEBUG 1

void alighment_cpu (char* seq_a, char* seq_b,
                    char* alignet_a, chat* alignment_b)
{

  size_t length_a = strlen(seq_a);
  size_t length_b = strlen(seq_b);
  
  // Calculate largest amount of mem needed
  unsigned int longest_alignment = (unsigned int)length_a +
                                   (unsigned int)length_b;
  
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
  printf("Malloc'd score matrices!\n");
  #endif

  // Fill in traceback matrix

  unsigned int i, j;
  
  // [0][0]
  match_score[0] = 0;
    int mismatch = -1;
    int match = 1;
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
        int current_score = seq_a[seq_i] == seq_b[seq_j] ? 1 : -1;
                                     
      match_score = max(match_score[old_index1]
                        , match_score[old_index2]
                        , match_score[old_index3])
                    + current_score;

    }
  }
}

__device__
void alighment_gpu(char* seq_a, char* seq_b,
                    char* alignet_a, chat* alignment_b)
{


}
int main(void){
    //generate sequences


    //time()
    //runing alignment on CPU
    //time()    


    //run alignment on GPU




}
