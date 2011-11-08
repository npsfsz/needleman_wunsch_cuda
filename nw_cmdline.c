/*
 nw_cmdline.c
 project: NeedlemanWunsch
 author: Isaac Turner <isaac.turner@dtc.ox.ac.uk>
 Copyright (C) 6-Nov-2011
 
 == Build
 cp ../utility_lib/utility_lib.* .
 cp ../string_buffer/string_buffer.* .
 gcc -o NeedlemanWunsch -Wall -lz nw_cmdline.c needleman_wunsch.c utility_lib.c string_buffer.c

 == Development
 - add non-affine gap alignment

 == License
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <zlib.h>

#include "utility_lib.h"
#include "string_buffer.h"

#include "needleman_wunsch.h"

char* cmd;

void print_usage(char* err_msg)
{
  if(err_msg != NULL)
  {
    fprintf(stderr, "Error: %s\n", err_msg);
  }

  fprintf(stderr, "usage: %s [OPTIONS] [seq1 seq2]\n", cmd);
  fprintf(stderr, "  Needleman-Wunsch optimal global alignment (maximises score).  \n");
  fprintf(stderr, "  Takes a pair of sequences on the command line, reads from a\n");
  fprintf(stderr, "  file and from sequence piped in.  Can read gzip files and FASTA.  \n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  OPTIONS:\n");
  fprintf(stderr, "    --match <score>      default: %i\n", nw_match_penalty);
  fprintf(stderr, "    --mismatch <score>   default: %i\n", nw_mismatch_penalty);
  fprintf(stderr, "    --gapopen <score>    default: %i\n", nw_gap_open_penalty);
  fprintf(stderr, "    --gapextend <score>  default: %i\n", nw_gap_extend_penalty);
  fprintf(stderr, "    --file <file>\n");
  fprintf(stderr, "    --printscores\n");
  fprintf(stderr, "    --stdin              Read from STDIN or piped input\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  To do alignment without affine gap, set `gapopen' to match `gapextend'\n");
  fprintf(stderr, "                       <turner.isaac@gmail.com>\n");

  exit(EXIT_FAILURE);
}

// Reads the sequence after the FASTA header line
// Returns: 1 if read successful, 0 if not
char read_fasta_sequence(STRING_BUFFER* sbuf, gzFile* file)
{
  // By this point have already read FASTA header line ('>..')
  // Just need to read sequence

  int line_first_char = sbuf->len;
  int read_len = string_buff_readline(sbuf, file);

  if(read_len == -1) {
    return 0;
  }

  while(read_len > 0 && sbuf->buff[line_first_char] != '>')
  {
    string_buff_chomp(sbuf);
    line_first_char = sbuf->len;
    read_len = string_buff_readline(sbuf, file);
  }

  if(sbuf->buff[line_first_char] == '>')
  {
    // Remove '>...' line from buffer
    string_buff_shrink(sbuf, line_first_char);
  }
  
  return 1;
}

void align_from_file(gzFile* file,
                     int sub_matrix[4][4],
                     const int gap_open, const int gap_extend,
                     const char print_score,
                     char **alignment_a, char **alignment_b, int *max_alignment)
{
  STRING_BUFFER* line1 = string_buff_init(200);
  STRING_BUFFER* line2 = string_buff_init(200);
  
  // Read first line
  char is_fasta = 0;
  int first_char = gzgetc(file);
  
  if(first_char == -1) {
    return;
  }
  else if(first_char == '>')
  {
    // Reading FASTA
    is_fasta = 1;
    string_buff_readline(line1, file);
  }
  else {
    is_fasta = 0;
    // Put char back
    gzungetc(first_char, file);
  }

  char reading = 1;

  while(reading)
  {
    if(is_fasta)
    {
      string_buff_reset(line1);
      string_buff_reset(line2);

      reading = read_fasta_sequence(line1, file);
      
      if(!reading)
      {
        break;
      }
      
      reading = read_fasta_sequence(line2, file);

      if(!reading)
      {
        print_usage("Odd number of FASTA sequences - I read in pairs!");
      }
    }
    else
    {
      int bases_read = string_buff_reset_readline(line1, file);

      if(bases_read == -1)
      {
        reading = 0;
        break;
      }
      
      bases_read = string_buff_reset_readline(line2, file);
      
      if(bases_read == -1)
      {
        print_usage("Odd number of sequences - I read in pairs!");
        reading = 0;
        break;
      }
      
      string_buff_chomp(line1);
      string_buff_chomp(line2);
    }

    // Check memory
    int new_max_alignment = line1->len + line2->len;

    if(new_max_alignment > *max_alignment)
    {
      // Expand memory used for storing result
      *max_alignment = new_max_alignment;
      nw_realloc_mem(new_max_alignment, alignment_a, alignment_b);
    }

    // Align
    int score = needleman_wunsch_affine(line1->buff, line2->buff,
                                        alignment_a, alignment_b,
                                        sub_matrix,
                                        gap_open, gap_extend);

    if(print_score)
    {
      printf("%s\n%s\nscore: %i\n", *alignment_a, *alignment_b, score);
    }
    else
    {
      printf("%s\n%s\n", *alignment_a, *alignment_b);
    }
  }
}

void unknown_option(char *option)
{
  char *format = "Unkown option '%s'";
  int length = strlen(format) + strlen(option) + 1;
  char *err_msg = (char*) malloc(length*sizeof(char));
  sprintf(err_msg, format, option);
  print_usage(err_msg);
}

int main(int argc, char* argv[])
{
  if(argc == 3 && argv[1][0] != '-' && argv[2][0] != '-')
  {
    // default needleman wunsch
    // > nw ACGT AAGT
    char *alignment_a, *alignment_b;

    nw_alloc_mem(argv[1], argv[2], &alignment_a, &alignment_b);

    needleman_wunsch_affine(argv[1], argv[2],
                            &alignment_a, &alignment_b,
                            nw_simple_sub_matrix,
                            nw_gap_open_penalty, nw_gap_extend_penalty);

    printf("%s\n%s\n", alignment_a, alignment_b);
    return 0;
  }

  cmd = argv[0];
  
  char *seq1 = NULL, *seq2 = NULL;
  
  // Set penalty defaults
  int match = nw_match_penalty;
  int mismatch = nw_mismatch_penalty;
  int gapopen = nw_gap_open_penalty;
  int gapextend = nw_gap_extend_penalty;

  gzFile* file = NULL;

  char print_score = 0;
  char read_stdin = 0;

  int argi;
  for(argi = 1; argi < argc; argi++)
  {
    if(argv[argi][0] == '-')
    {
      // strcasecmp does case insensitive comparison
      if(strcasecmp(argv[argi], "--printscores") == 0)
      {
        print_score = 1;
      }
      else if(strcasecmp(argv[argi], "--stdin") == 0)
      {
        read_stdin = 1;
      }
      else if(argi == argc-1)
      {
        // All the remaining options take an extra argument
        unknown_option(argv[argi]);
      }
      else if(strcasecmp(argv[argi], "--match") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &match)) {
          print_usage("Invalid match -- must be int");
        }
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--mismatch") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &mismatch)) {
          print_usage("Invalid mismatch score -- must be int");
        }
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--gapopen") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &gapopen)) {
          print_usage("Invalid gap open score -- must be int");
        }
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--gapextend") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &gapextend)) {
          print_usage("Invalid gap extend score -- must be int");
        }
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--file") == 0)
      {
        file = gzopen(argv[argi+1], "r");
        argi++; // took an argument
      }
      else
      {
        // Error - unknown option
        unknown_option(argv[argi]);
      }
    }
    else {
      break;
    }
  }
  
  // Check for extra unused arguments
  // and set seq1 and seq2 if they have been passed
  if(argi < argc)
  {
    if(argc - argi != 2)
    {
      print_usage("Unknown options");
    }
    else {
      seq1 = argv[argi];
      seq2 = argv[argi+1];
    }
  }
  
  // adjust nw_simple_sub_matrix
  int substitution_matrix[4][4];
  nw_create_sub_matrix(match, mismatch, &substitution_matrix);

  //int substitution_matrix[4][4] = nw_create_sub_matrix(match, mismatch);

  char *alignment_a = NULL, *alignment_b = NULL;
  int alignment_max_length;
  int alignment_score;

  if(seq1 != NULL)
  {
    // Align seq1 and seq2
    alignment_max_length = nw_alloc_mem(seq1, seq2, &alignment_a, &alignment_b);
    
    alignment_score = needleman_wunsch_affine(seq1, seq2,
                                              &alignment_a, &alignment_b,
                                              substitution_matrix,
                                              gapopen, gapextend);
    
    printf("%s\n%s\n", alignment_a, alignment_b);
    
    if(print_score) {
      printf("%i\n", alignment_score);
    }
  }
  else
  {
    alignment_max_length = 1000; // equivalent to two strings of 500bp
    alignment_a = (char*) malloc((alignment_max_length+1) * sizeof(char));
    alignment_b = (char*) malloc((alignment_max_length+1) * sizeof(char));
  }


  if(file != NULL)
  {
    // read file
    align_from_file(file, substitution_matrix, gapopen, gapextend, print_score,
                    &alignment_a, &alignment_b, &alignment_max_length);
    
    gzclose(file);
  }

  //char piped_data = stdin_is_ready();
  char piped_data = 0;

  if(read_stdin || piped_data)
  {
    // Read from STDIN
    file = gzdopen(fileno(stdin), "r");
    
    // read file
    align_from_file(file, substitution_matrix, gapopen, gapextend, print_score,
                    &alignment_a, &alignment_b, &alignment_max_length);
    
    gzclose(file);
  }
  else if(argc == 1)
  {
    print_usage(NULL);
  }

  // Free memory
  free(alignment_a);
  free(alignment_b);

  return 0;
}
