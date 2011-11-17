/*
 needleman_wunsch.c
 project: NeedlemanWunsch
 author: Isaac Turner <isaac.turner@dtc.ox.ac.uk>
 Copyright (C) 25-May-2011
 
 To test/compile, see nw_test.c
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 ToDo:
 - Testing...
 - inversion alignment?
 - Gap end penalty??
 */

// Turn on debugging output by defining DEBUG
//#define DEBUG

#define arr_lookup(arr,width,i,j) arr[((j)*(width)) + (i)]
#define ARR_INDEX(width,i,j) ((j)*(width)) + (i)
#define MATCH_PENALTY(sim_mat,seq1,seq2,i,j) sim_mat[_get_base_index(seq1[i])][_get_base_index(seq2[j])]
#define MAX_3(x,y,z) ((x) >= (y) && (x) >= (z) ? (x) : ((y) >= (z) ? (y) : (z)))

#define MATCH 0
#define GAP_A 1
#define GAP_B 2

#include <stdlib.h>
#include <stdio.h>
#include <strings.h>

// DEBUG
#ifdef DEBUG
#define INT_MIN -99
#else
#include <limits.h>
#endif

#include "needleman_wunsch.h"

// Default penalty scores
const int nw_match_penalty = 1;
const int nw_mismatch_penalty = -1;

const int nw_gap_open_penalty = -10;
const int nw_gap_extend_penalty = -2;

// This is consistent with nw_match_penalty and nw_mismatch_penalty
// should be equivalent to the result from:
// nw_create_sub_matrix(nw_match_penalty, nw_mismatch_penalty, &matrix)
int nw_simple_sub_matrix[4][4] = {{ 1,-1,-1,-1},
                                  {-1, 1,-1,-1},
                                  {-1,-1, 1,-1},
                                  {-1,-1,-1, 1}};

int _get_base_index(const char b)
{
  switch(b)
  {
    case 'A':
    case 'a':
      return 0;
    case 'C':
    case 'c':
      return 1;
    case 'G':
    case 'g':
      return 2;
    case 'T':
    case 't':
      return 3;
    default:
      fprintf(stderr, "Error: Invalid base in _get_base_index (%i) '%c'\n",
              (int)b, b);
      exit(EXIT_FAILURE);
  }
}

// construct substitution matrix with match/mismatch scores
void nw_create_sub_matrix(int match, int mismatch, int (*matrix)[4][4])
//int[4][4] nw_create_sub_matrix(int match, int mismatch)
{
  int i,j;
  
  for(j = 0; j < 4; j++)
  {
    for(i = 0; i < 4; i++)
    {
      (*matrix)[i][j] = (i == j ? match : mismatch);
    }
  }
}

/*
 * - both strings should be of the same length and end with \0 char
 * - characters should be one of: aAcCgGtT
 * - gaps should be '-'
 */
int score_alignment(const char* alignment_a, const char* alignment_b,
                    int match_penalties[4][4], const int gap_penalty)
{
  int i;
  int score = 0;

  for(i = 0; alignment_a[i] != '\0'; i++)
  {
    if(alignment_a[i] == '-' || alignment_b[i] == '-')
    {
      // gap
      score += gap_penalty;
    }
    else
    {
      // both bases one of: ACGT
      score += MATCH_PENALTY(match_penalties, alignment_a, alignment_b, i, i);
      
#ifdef DEBUG
      int a = _get_base_index(alignment_a[i]);
      int b = _get_base_index(alignment_b[i]);
      printf("%c %i, %c %i = %i\n", alignment_a[i], a, alignment_b[i], b,
             match_penalties[a][b]);
#endif DEBUG
    }
  }

  return score;
}

/*
 * - both strings should be of the same length and end with \0 char
 * - characters should be one of: aAcCgGtT
 * - gaps should be '-'
 */
int score_alignment_affine(const char* alignment_a, const char* alignment_b,
                           int match_penalties[4][4],
                           const int gap_penalty_start,
                           const int gap_penalty_cont)
{
  int score = 0;
  char in_gap_a = 0;
  char in_gap_b = 0;
  
  int i;
  for(i = 0; alignment_a[i] != '\0'; i++)
  {
    if(alignment_a[i] == '-')
    {
      if(in_gap_a) {
        score += gap_penalty_cont;
      }
      else {
        in_gap_a = 1;
        score += gap_penalty_start;
      }
      in_gap_b = 0;
    }
    else if(alignment_b[i] == '-')
    {
      if(in_gap_b) {
        score += gap_penalty_cont;
      }
      else {
        in_gap_b = 1;
        score += gap_penalty_start;
      }
      in_gap_a = 0;
    }
    else
    {
      in_gap_a = 0;
      in_gap_b = 0;
      
      score += MATCH_PENALTY(match_penalties, alignment_a, alignment_b, i, i);
    }
  }
  
  return score;
}


/* Allocate memory for alignment results */

int nw_alloc_mem(const char* seq_a, const char* seq_b,
                 char** alignment_a, char** alignment_b)
{
  int length_a = strlen(seq_a);
  int length_b = strlen(seq_b);
  
  // Calculate largest amount of mem needed
  int longest_alignment = length_a + length_b;
  
  // longest_alignment + 1 to allow for \0
  *alignment_a = (char*) malloc((longest_alignment+1) * sizeof(char));
  *alignment_b = (char*) malloc((longest_alignment+1) * sizeof(char));

  return longest_alignment;
}

// length is = length_a + length_b
int nw_realloc_mem(unsigned int length, char** alignment_a, char** alignment_b)
{
  // longest_alignment + 1 to allow for \0
  *alignment_a = realloc(*alignment_a, (length+1) * sizeof(char));
  *alignment_b = realloc(*alignment_b, (length+1) * sizeof(char));

  return length;
}


/* Needleman-Wunsch alignment WITHOUT affine gap penalty scores */

int needleman_wunsch(const char* seq_a, const char* seq_b,
                     char* alignment_a, char* alignment_b,
                     int match_penalties[4][4], const int gap_penalty)
{
  int length_a = strlen(seq_a);
  int length_b = strlen(seq_b);

  // Calculate largest amount of mem needed
  int longest_alignment = length_a + length_b;
  
  // Store result in positions passed
  //char *alignment_a = *result_a;
  //char *alignment_b = *result_b;

  int score_width = length_a+1;
  int score_height = length_b+1;
  
  int arr_size = score_width * score_height;
  
  // 2d array (length_a x length_b)
  int* max_score = (int*) malloc(arr_size * sizeof(int));
  
  int i, j;

  for(i = 0; i < score_width; i++)
  {
    // [i][0]
    max_score[i] = gap_penalty*i;
  }
  
  for(j = 1; j < score_height; j++)
  {
    // [0][j]
    max_score[j*score_width] = gap_penalty*j;
  }
  
  for(i = 0; i < length_a; i++)
  {
    for(j = 0; j < length_b; j++)
    {
      int match = arr_lookup(max_score, score_width, i, j) +
                  MATCH_PENALTY(match_penalties, seq_a, seq_b, i, j);
    
      // Add a base from the first sequence
      int delete = arr_lookup(max_score, score_width, i, j+1) + gap_penalty;
      // Add a base from the second sequence
      int insert = arr_lookup(max_score, score_width, i+1, j) + gap_penalty;
      
      arr_lookup(max_score, score_width, i+1, j+1) = MAX_3(match, delete, insert);
    }
  }
  
#ifdef DEBUG
  // Debug: print max_score matrix
  for(j = 0; j <= length_b; j++)
  {
    printf("%i:", j);
    for(i = 0; i <= length_a; i++)
    {
      printf(" %3i", arr_lookup(max_score, length_a+1, i, j));
    }
    printf("\n");
  }
#endif
  
  // (backwards, then shift into place)
  int next_char;
  
  i = length_a - 1;
  j = length_b - 1;
  
  for(next_char = longest_alignment-1; i >= 0 && j >= 0; next_char--)
  {
    int curr_score = arr_lookup(max_score, score_width, i+1, j+1);

    int score_up = arr_lookup(max_score, score_width, i+1, j);
    int score_diag = arr_lookup(max_score, score_width, i, j);
    int score_left = arr_lookup(max_score, score_width, i, j+1);
  
    if(curr_score == score_diag + MATCH_PENALTY(match_penalties, seq_a, seq_b, i, j))
    {
      // Both sequences have a base
      alignment_a[next_char] = seq_a[i];
      alignment_b[next_char] = seq_b[j];
      i--;
      j--;
    }
    else if(curr_score == score_left + gap_penalty)
    {
      // Gap in seq_b
      alignment_a[next_char] = seq_a[i];
      alignment_b[next_char] = '-';
      i--;
    }
    else if(curr_score == score_up + gap_penalty)
    {
      // Gap in seq_a
      alignment_a[next_char] = '-';
      alignment_b[next_char] = seq_b[j];
      j--;
    }
    else
    {
      fprintf(stderr, "Something went wrong\n");
      exit(EXIT_FAILURE);
    }
  }
  
  while(i >= 0)
  {
    alignment_a[next_char] = seq_a[i];
    alignment_b[next_char] = '-';
    next_char--;
    i--;
  }
  
  while(j >= 0)
  {
    alignment_a[next_char] = '-';
    alignment_b[next_char] = seq_b[j];
    next_char--;
    j--;
  }
  
  // length of alignment
  int first_char = next_char+1;
  int alignment_len = longest_alignment - first_char;
  
#ifdef DEBUG
  printf("first_char %i; longest_alignment %i; length %i;\n",
         first_char, longest_alignment, alignment_len);
#endif

  // shift back into 0th position in char arrays
  
  int pos;
  for(pos = 0; pos < alignment_len; pos++)
  {
    alignment_a[pos] = alignment_a[pos+first_char];
    alignment_b[pos] = alignment_b[pos+first_char];
  }
  
  alignment_a[pos] = '\0';
  alignment_b[pos] = '\0';
  
  int max_alignment_score = max_score[arr_size-1];

  // free memory
  free(max_score);
  
  // Highest score (botttom right)
  return max_alignment_score;
}


/**
 * Align with gap start and continue penalties
 */
int needleman_wunsch_affine(const char* seq_a, const char* seq_b,
                            char* alignment_a, char* alignment_b,
                            int match_penalties[4][4], // subsitution penalty matrix
                            const int gap_penalty_start,
                            const int gap_penalty_cont)
{
  if(gap_penalty_start == gap_penalty_cont)
  {
    // Speed up - if gap start and gap extend are the same,
    // no need to call _affine
    return needleman_wunsch(seq_a, seq_b,
                            alignment_a, alignment_b,
                            match_penalties,
                            gap_penalty_cont);
  }

#ifdef DEBUG
   // DEBUG
  printf("sim matrix:\n");
  printf("     A    C    G    T\n");
  int x,y;
  char* str = "ACGT";
  for(y = 0; y<4; y++) {
    printf("%c", str[y]);
    for(x = 0; x<4; x++) {
      printf(" %4i", match_penalties[x][y]);
    }
    printf("\n");
  }

  printf("gap_penalty_start: %i\n", gap_penalty_start);
  printf("gap_penalty_cont: %i\n", gap_penalty_cont);
  printf("\n");
#endif
  
  int length_a = strlen(seq_a);
  int length_b = strlen(seq_b);
  
  // Calculate largest amount of mem needed
  int longest_alignment = length_a + length_b;
  
  // Store result in positions passed
  //char *alignment_a = *result_a;
  //char *alignment_b = *result_b;
  
  int score_width = length_a+1;
  int score_height = length_b+1;
  
  int arr_size = score_width * score_height;
  
  // 2d array (length_a x length_b)
  int* match_score = (int*) malloc(arr_size * sizeof(int));
  int* gap_a_score = (int*) malloc(arr_size * sizeof(int)); // aka delete wrt seq_a
  int* gap_b_score = (int*) malloc(arr_size * sizeof(int)); // aka insert wrt seq_a
  
  int i, j, index;
  
  // [0][0]
  match_score[0] = 0;
  gap_a_score[0] = 0;
  gap_b_score[0] = 0;
  
  // work along first row -> [i][0]
  for(i = 1; i < score_width; i++)
  {
    match_score[i] = INT_MIN;
    
    // Think carefully about which way round these are
    // -> can't have a gap in B at first row
    // x(i - 1) because first gap has gap_penalty_start
    gap_a_score[i] = INT_MIN;
    gap_b_score[i] = gap_penalty_start + (i - 1) * gap_penalty_cont;
  }
  
  // work down first column -> [0][j]
  for(j = 1; j < score_height; j++)
  {
    index = j*score_width;
    match_score[index] = INT_MIN;
    
    // Think carefully about which way round these are
    // -> can't have a gap in A at first column
    // x(i - 1) because first gap has gap_penalty_start
    gap_a_score[index] = gap_penalty_start + (j - 1) * gap_penalty_cont;
    gap_b_score[index] = INT_MIN;
  }
  
  // Update Dynamic Programming arrays
  int seq_i, seq_j, sub_penalty;
  int old_index, new_index;

  for (i = 1; i < score_width; i++)
  {
    for (j = 1; j < score_height; j++)
    {
      // It's an indexing thing...
      seq_i = i - 1;
      seq_j = j - 1;
      
      sub_penalty = MATCH_PENALTY(match_penalties, seq_a, seq_b, seq_i, seq_j);
      
      // Update match_score[i][j] with position [i-1][j-1]
      new_index = j*score_width + i;
      old_index = (j-1)*score_width + (i-1);
      
      // substitution
      match_score[new_index] = MAX_3(match_score[old_index], // continue alignment
                                     gap_a_score[old_index], // close gap in seq_a
                                     gap_b_score[old_index]) + sub_penalty; // close gap in seq_b
      
      // Update gap_a_score[i][j] with position [i][j-1]
      old_index = (j-1)*score_width + i;
      
      // Long arithmetic since some INTs are set to INT_MIN and penalty is -ve
      // (adding as ints would cause an integer overflow)
      gap_a_score[new_index] = MAX_3((long)match_score[old_index] + gap_penalty_start, // opening gap
                                     (long)gap_a_score[old_index] + gap_penalty_cont, // continuing gap
                                     (long)gap_b_score[old_index] + gap_penalty_start); // opening gap
    
      // Update gap_b_score[i][j] with position [i-1][j]
      old_index = j*score_width + (i-1);
      
      gap_b_score[new_index] = MAX_3((long)match_score[old_index] + gap_penalty_start, // opening gap
                                     (long)gap_a_score[old_index] + gap_penalty_start, // opening gap
                                     (long)gap_b_score[old_index] + gap_penalty_cont); // continuing gap
    }
  }
  
  
#ifdef DEBUG
  // Debug: print score matrices
  printf("match_score:\n");
  for(j = 0; j <= length_b; j++)
  {
    printf("%i:", j);
    for(i = 0; i <= length_a; i++)
    {
      printf(" %3i", arr_lookup(match_score, score_width, i, j));
    }
    printf("\n");
  }
  printf("gap_a_score:\n");
  for(j = 0; j <= length_b; j++)
  {
    printf("%i:", j);
    for(i = 0; i <= length_a; i++)
    {
      printf(" %3i", arr_lookup(gap_a_score, score_width, i, j));
    }
    printf("\n");
  }
  printf("gap_b_score:\n");
  for(j = 0; j <= length_b; j++)
  {
    printf("%i:", j);
    for(i = 0; i <= length_a; i++)
    {
      printf(" %3i", arr_lookup(gap_b_score, score_width, i, j));
    }
    printf("\n");
  }
#endif
   
  //
  // Trace back now (score matrices all calculated)
  //
  
  // work backwards re-tracing optimal alignment, then shift sequences into place
  
  // Get max score (and therefore current matrix)
  char curr_matrix = MATCH;
  int curr_score = match_score[arr_size-1];
  
  if(gap_a_score[arr_size-1] > curr_score)
  {
    curr_matrix = GAP_A;
    curr_score = gap_a_score[arr_size-1];
  }
  
  if(gap_b_score[arr_size-1] > curr_score)
  {
    curr_matrix = GAP_B;
    curr_score = gap_b_score[arr_size-1];
  }

  long max_alignment_score = curr_score;
  
  int score_i, score_j;
  
  seq_i = length_a - 1;
  seq_j = length_b - 1;
  
  int next_char;
  
  // Previous scores on each matrix
  int prev_match_score, prev_gap_a_score, prev_gap_b_score;
  
  // penalties if coming from each of the prev matrices
  int prev_match_penalty, prev_gap_a_penalty, prev_gap_b_penalty;
  
  // longest_alignment = strlen(seq_a) + strlen(seq_b)
  for(next_char = longest_alignment-1; seq_i >= 0 && seq_j >= 0; next_char--)
  {
    switch (curr_matrix)
    {
      case MATCH:
        #ifdef DEBUG
          printf("MATCH\n");
        #endif

        alignment_a[next_char] = seq_a[seq_i];
        alignment_b[next_char] = seq_b[seq_j];
        
        // Match
        prev_match_penalty = MATCH_PENALTY(match_penalties,
                                           seq_a, seq_b,
                                           seq_i, seq_j);
        
        prev_gap_a_penalty = prev_match_penalty; // match
        prev_gap_b_penalty = prev_match_penalty; // match
        
        // Moving back on i and j
        seq_i--;
        seq_j--;
        break;
      case GAP_A:
        #ifdef DEBUG
          printf("GAP_A\n");
        #endif

        alignment_a[next_char] = '-';
        alignment_b[next_char] = seq_b[seq_j];
        
        prev_match_penalty = gap_penalty_start; // match
        prev_gap_a_penalty = gap_penalty_cont; // starting gap on a
        prev_gap_b_penalty = gap_penalty_start; // continuing gap on b
        
        // Moving back on j
        seq_j--;
        break;
      case GAP_B:
        #ifdef DEBUG
          printf("GAP_B\n");
        #endif
        
        alignment_a[next_char] = seq_a[seq_i];
        alignment_b[next_char] = '-';
        
        prev_match_penalty = gap_penalty_start;
        prev_gap_a_penalty = gap_penalty_start;
        prev_gap_b_penalty = gap_penalty_cont;
        
        // Moving back on i
        seq_i--;
        break;
      default:
        fprintf(stderr, "Err: invalid matrix number\n");
        exit(EXIT_FAILURE);
    }
    
    // Current score matrix position is [seq_i+1][seq_j+1]
    score_i = seq_i + 1;
    score_j = seq_j + 1;
    
    // [score_i][score_j] is the next position in the score matrices
    
    prev_match_score = arr_lookup(match_score, score_width, score_i, score_j);
    prev_gap_a_score = arr_lookup(gap_a_score, score_width, score_i, score_j);
    prev_gap_b_score = arr_lookup(gap_b_score, score_width, score_i, score_j);
    
#ifdef DEBUG
     // DEBUG
    printf("(%i, %i) curr_score: %i\n", seq_i, seq_j, curr_score);
    printf("prev_match_score: %i\n", prev_match_score);
    printf("prev_gap_a_score: %i\n", prev_gap_a_score);
    printf("prev_gap_b_score: %i\n", prev_gap_b_score);
#endif
    
    // Now figure out which matrix we came from
    if(prev_match_score + prev_match_penalty == curr_score)
    {
      // Both sequences have a base
      curr_matrix = MATCH;
      curr_score = prev_match_score;
    }
    else if(prev_gap_a_score + prev_gap_a_penalty == curr_score)
    {
      // Gap in seq_a
      curr_matrix = GAP_A;
      curr_score = prev_gap_a_score;
    }
    else if(prev_gap_b_score + prev_gap_b_penalty == curr_score)
    {
      // Gap in seq_b
      curr_matrix = GAP_B;
      curr_score = prev_gap_b_score;
    }
    else {
      fprintf(stderr, "Fail\n");
      exit(EXIT_FAILURE);
    }
  }
  
  // Free memory
  free(match_score);
  free(gap_a_score);
  free(gap_b_score);
  
  /*
  if(curr_matrix == MATCH && seq_i >= 0 && seq_j >= 0)
  {
    // Needed to check seq_i and seq_j incase length was < 1
    // starts with a match
    alignment_a[next_char] = seq_a[seq_i];
    alignment_b[next_char] = seq_b[seq_j];
    next_char--;
    seq_i--;
    seq_j--;
  }
  else */
  
  if(curr_matrix == GAP_A)
  {
    // Gap in A
    while(seq_j >= 0)
    {
      alignment_a[next_char] = '-';
      alignment_b[next_char] = seq_b[seq_j];
      next_char--;
      seq_j--;
    }
  }
  else if(curr_matrix == GAP_B)
  {
    // Gap in B
    while(seq_i >= 0)
    {
      alignment_a[next_char] = seq_a[seq_i];
      alignment_b[next_char] = '-';
      next_char--;
      seq_i--;
    }
  }
  
  // Shift alignment strings back into 0th position in char arrays
  int first_char = next_char+1;
  int alignment_len = longest_alignment - first_char;

#ifdef DEBUG
  printf("first_char %i; longest_alignment %i; length %i;\n",
         first_char, longest_alignment, alignment_len);
#endif DEBUG

  int pos;
  for(pos = 0; pos < alignment_len; pos++)
  {
    alignment_a[pos] = alignment_a[pos+first_char];
    alignment_b[pos] = alignment_b[pos+first_char];
  }
  
  alignment_a[pos] = '\0';
  alignment_b[pos] = '\0';

  return max_alignment_score;
}
