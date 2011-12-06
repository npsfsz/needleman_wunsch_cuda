/*
 nw_cmdline.c
 project: NeedlemanWunsch
 author: Isaac Turner <turner.isaac@gmail.com>
 url: http://sourceforge.net/projects/needlemanwunsch
 Copyright (C) 06-Dec-2011
 
 see: README

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
#include <ctype.h> // tolower isspace

#include <zlib.h>

#include "utility_lib.h"
#include "string_buffer.h"

#include "nw_load_scores.h"
#include "needleman_wunsch.h"

// Printing
char* mismatch_start_colour = "\033[92m"; // Mismatch (GREEN)
char* indel_start_colour = "\033[91m"; // Insertion / deletion (RED)
char* colour_stop = "\033[0m";

// Defaults
int nw_match_default = 1;
int nw_mismatch_default = -2;
int nw_gap_open_default = -4;
int nw_gap_extend_default = -1;

// For this run
char* cmd;
char print_colour = 0, print_pretty = 0, print_scores = 0,
     print_fasta = 0, print_zam = 0;

NW_SCORING* scoring;

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
  fprintf(stderr, "    --file <file>        Sequence file reading with gzip support\n");
  fprintf(stderr, "    --stdin              Read from STDIN (same as '--file -')\n");
  fprintf(stderr, "    --case_sensitive     Case sensitive character comparison\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "    --scoring <PAM30|PAM70|BLOSUM80|BLOSUM62>\n");
  fprintf(stderr, "    --substitution_matrix <file>  see details for formatting\n");
  fprintf(stderr, "    --substitution_pairs <file>   see details for formatting\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "    --match <score>      default: %i\n", nw_match_default);
  fprintf(stderr, "    --mismatch <score>   default: %i\n", nw_mismatch_default);
  fprintf(stderr, "    --gapopen <score>    default: %i\n", nw_gap_open_default);
  fprintf(stderr, "    --gapextend <score>  default: %i\n", nw_gap_extend_default);
  fprintf(stderr, "\n");
  fprintf(stderr, "    --freestartgap       No penalty for gap at start of alignment\n");
  fprintf(stderr, "    --freeendgap         No penalty for gap at end of alignment\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "    --printscores        Print optimal alignment scores\n");
  fprintf(stderr, "    --printfasta         Print fasta header lines\n");
  fprintf(stderr, "    --pretty             Print with a descriptor line\n");
  fprintf(stderr, "    --colour             Print with in colour\n");
  fprintf(stderr, "    --zam                A funky type of output\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " DETAILS:\n");
  fprintf(stderr, "  * For help choosing scoring, see the README file. \n");
  fprintf(stderr, "  * Gap (of length N) penalty is: (open+N*extend)\n");
  fprintf(stderr, "  * To do alignment without affine gap, set '--gapopen 0'.\n");
  fprintf(stderr, "  * Scoring files should be matrices, with entries separated\n"
                  "    by a single character or whitespace.  See files in the\n"
                  "    'scores' directory for examples.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  turner.isaac@gmail.com  06 Dec 2011\n");

  exit(EXIT_FAILURE);
}

char read_fasta_entry(STRING_BUFFER* header, STRING_BUFFER* sequence,
                      gzFile* file)
{
  string_buff_reset(header);
  string_buff_reset(sequence);

  // Read until we have a header line
  int read_len;

  do
  {
    read_len = string_buff_readline(header, file);
    string_buff_chomp(header);
  }
  while (read_len > 0 && header->len == 0);

  if(read_len <= 0)
  {
    return 0;
  }
  else if(header->buff[0] != '>')
  {
    fprintf(stderr, "Warning: FASTA header does not begin with '>' ('%s')\n",
            header->buff);
  }

  string_buff_chomp(header);

  char success = 0;

  while(1)
  {
    int first_char = gzgetc(file);
  
    if(first_char == -1)
    {
      break;
    }
    else if(first_char == '>')
    {
      // Push char back onto buffer
      gzungetc(first_char, file);
      break;
    }

    // Push char onto string
    string_buff_add_char(sequence, first_char);
    
    string_buff_chomp(sequence);
    success = 1;

    if(first_char != '\n' && first_char != '\r')
    {
      // Read the rest of the line
      if(string_buff_readline(sequence, file) < 0) {
        break;
      }

      string_buff_chomp(sequence);
    }
  }

  return success;
}

void colour_print_alignment_against(char *alignment_a, char *alignment_b)
{
  int i;
  char red = 0, green = 0;

  for(i = 0; alignment_a[i] != '\0'; i++)
  {
    if(alignment_b[i] == '-')
    {
      if(!red)
      {
        printf("%s", indel_start_colour);
        red = 1;
      }
    }
    else if(red)
    {
      red = 0;
      printf("%s", colour_stop);
    }
    
    
    if(((scoring->case_sensitive && alignment_a[i] != alignment_b[i]) ||
        tolower(alignment_a[i]) != tolower(alignment_b[i])) &&
       alignment_a[i] != '-' && alignment_b[i] != '-')
    {
      if(!green)
      {
        printf("%s", mismatch_start_colour);
        green = 1;
      }
    }
    else if(green)
    {
      green = 0;
      printf("%s", colour_stop);
    }

    printf("%c", alignment_a[i]);
  }

  if(green || red)
  {
    // Stop all colours
    printf("%s", colour_stop);
  }

  printf("\n");
}

void align_zam(char *seq_a, char *seq_b, char *alignment_a, char *alignment_b)
{
  needleman_wunsch(seq_a, seq_b, alignment_a, alignment_b, scoring);

  int i;
  for(i = 0; alignment_a[i] != '\0'; i++)
  {
    if(alignment_a[i] == '-')
    {
      alignment_a[i] = '_';
    }

    if(alignment_b[i] == '-')
    {
      alignment_b[i] = '_';
    }
  }

  int num_of_mismatches = 0;
  int num_of_indels = 0;

  // Print branch 1
  printf("Br1:%s\n", alignment_a);
  
  // Print spacer
  printf("    ");

  for(i = 0; alignment_a[i] != '\0'; i++)
  {
    if(alignment_a[i] == '_' || alignment_b[i] == '_')
    {
      printf(" ");
      num_of_indels++;
    }
    else if((scoring->case_sensitive && alignment_a[i] != alignment_b[i]) ||
            tolower(alignment_a[i]) != tolower(alignment_b[i]))
    {
      printf("*");
      num_of_mismatches++;
    }
    else
    {
      printf("|");
    }
  }

  printf("\n");

  // Print branch 2
  printf("Br2:%s\n", alignment_b);

  // print mismatch indel numbers
  printf("%i %i\n\n", num_of_mismatches, num_of_indels);
}

void align(char *seq_a, char *seq_b, char *alignment_a, char *alignment_b,
           char *seq_a_name, char *seq_b_name)
{
  if(print_zam)
  {
    return align_zam(seq_a, seq_b, alignment_a, alignment_b);
  }

  int score = needleman_wunsch(seq_a, seq_b, alignment_a, alignment_b, scoring);

  if(print_fasta && seq_a_name != NULL)
  {
    printf("%s\n", seq_a_name);
  }

  if(print_fasta && print_pretty && seq_b_name != NULL)
  {
    printf("%s\n", seq_b_name);
  }

  if(print_colour)
  {
    // Print alignment line 1
    colour_print_alignment_against(alignment_a, alignment_b);
  }
  else
  {
    printf("%s\n", alignment_a);
  }
  
  if(print_pretty)
  {
    // Print spacer
    int i;

    for(i = 0; alignment_a[i] != '\0'; i++)
    {
      if(alignment_a[i] == '-' || alignment_b[i] == '-')
      {
        printf(" ");
      }
      else if((scoring->case_sensitive && alignment_a[i] == alignment_b[i]) ||
              tolower(alignment_a[i]) == tolower(alignment_b[i]))
      {
        printf("|");
      }
      else
      {
        printf("*");
      }
    }
    
    printf("\n");
  }
  else if(print_fasta && seq_b_name != NULL)
  {
    printf("%s\n", seq_b_name);
  }

  if(print_colour)
  {
    // Print alignment line 2
    colour_print_alignment_against(alignment_b, alignment_a);
  }
  else
  {
    printf("%s\n", alignment_b);
  }

  if(print_scores)
  {
    printf("score: %i\n", score);
  }
  
  printf("\n");
}

void align_from_file(gzFile* file, char **alignment_a, char **alignment_b,
                     int *max_alignment)
{

  STRING_BUFFER *entry1_header = NULL, *entry2_header = NULL;
  STRING_BUFFER* entry1_seq = string_buff_init(200);
  STRING_BUFFER* entry2_seq = string_buff_init(200);
  
  // Read first line
  char is_fasta = 0;
  int first_char = gzgetc(file);
  
  if(first_char == -1)
  {
    fprintf(stderr, "Warning: empty sequence file\n");
    return;
  }
  else if(first_char == '>')
  {
    // Reading FASTA
    is_fasta = 1;
    
    entry1_header = string_buff_init(200);
    entry2_header = string_buff_init(200);
  }
  else {
    is_fasta = 0;
  }

  // Put char back
  gzungetc(first_char, file);

  char reading = 1;

  while(reading)
  {
    if(is_fasta)
    {
      reading = read_fasta_entry(entry1_header, entry1_seq, file);
      
      if(!reading)
      {
        break;
      }
      
      reading = read_fasta_entry(entry2_header, entry2_seq, file);

      if(!reading)
      {
        print_usage("Odd number of FASTA sequences - I read in pairs!");
      }
    }
    else
    {
      int bases_read = string_buff_reset_readline(entry1_seq, file);

      if(bases_read == -1)
      {
        reading = 0;
        break;
      }
      
      bases_read = string_buff_reset_readline(entry2_seq, file);
      
      if(bases_read == -1)
      {
        print_usage("Odd number of sequences - I read in pairs!");
        reading = 0;
        break;
      }
      
      string_buff_chomp(entry1_seq);
      string_buff_chomp(entry2_seq);
    }

    // Check memory
    int new_max_alignment = entry1_seq->len + entry2_seq->len;

    if(new_max_alignment > *max_alignment)
    {
      // Expand memory used for storing result
      *max_alignment = new_max_alignment;
      nw_realloc_mem(new_max_alignment, alignment_a, alignment_b);
    }

    // Align
    if(is_fasta)
    {
      align(entry1_seq->buff, entry2_seq->buff, *alignment_a, *alignment_b,
            entry1_header->buff, entry2_header->buff);
    }
    else
    {
      align(entry1_seq->buff, entry2_seq->buff, *alignment_a, *alignment_b,
            NULL, NULL);
    }
  }

  // Free memory
  string_buff_free(entry1_seq);
  string_buff_free(entry2_seq);

  if(is_fasta)
  {
    string_buff_free(entry1_header);
    string_buff_free(entry2_header);
  }
}

void unknown_option(char *option)
{
  char *format = "Unkown option or argument required '%s'";
  int length = strlen(format) + strlen(option) + 1;
  char *err_msg = (char*) malloc(length*sizeof(char));
  sprintf(err_msg, format, option);
  print_usage(err_msg);
}

int main(int argc, char* argv[])
{
  cmd = argv[0];

  #ifdef DEBUG
  printf("DEBUG: on\n");
  #endif

  if(argc == 1)
  {
    print_usage(NULL);
  }
  else if(argc == 3 && argv[1][0] != '-' && argv[2][0] != '-')
  {
    // default needleman wunsch
    // > nw ACGT AAGT
    char *alignment_a, *alignment_b;

    nw_alloc_mem(argv[1], argv[2], &alignment_a, &alignment_b);

    // set up scoring (zeros are to penalise gaps everwhere)
    scoring = simple_scoring(nw_match_default, nw_mismatch_default,
                             nw_gap_open_default, nw_gap_extend_default,
                             0, 0, 0);

    align(argv[1], argv[2], alignment_a, alignment_b, NULL, NULL);

    return 0;
  }

  char *seq1 = NULL, *seq2 = NULL;
  
  scoring = NULL;
  
  // Set penalty defaults
  int match = nw_match_default;
  int mismatch = nw_mismatch_default;
  int gap_open = nw_gap_open_default;
  int gap_extend = nw_gap_extend_default;

  char no_start_gap_penalty = 0, no_end_gap_penalty = 0;

  // Indicates which of the above parameters have been set on the command line
  // This is required because --scoring may be specified after and
  // you don't want to lose the specified values
  char use_match = 0, use_mismatch = 0, use_gap_open = 0, use_gap_extend = 0;
  char use_no_start_gap_penalty = 0, use_no_end_gap_penalty = 0;

  char* file_path = NULL;
  char read_stdin = 0;

  char case_sensitive = 0;

  // case sensitive needs to be dealt with first
  // (is is used to construct hash table for swap_table)
  int argi;
  for(argi = 1; argi < argc; argi++)
  {
    if(strcasecmp(argv[argi], "--case_sensitive") == 0)
    {
      case_sensitive = 1;
      break;
    }
  }

  for(argi = 1; argi < argc; argi++)
  {
    if(argv[argi][0] == '-')
    {
      // strcasecmp does case insensitive comparison
      if(strcasecmp(argv[argi], "--freestartgap") == 0)
      {
        no_start_gap_penalty = 1;
        use_no_start_gap_penalty = 1;
      }
      else if(strcasecmp(argv[argi], "--freeendgap") == 0)
      {
        no_end_gap_penalty = 1;
        use_no_end_gap_penalty = 1;
      }
      else if(strcasecmp(argv[argi], "--case_sensitive") == 0)
      {
        // Already dealt with
        //case_sensitive = 1;
      }
      else if(strcasecmp(argv[argi], "--printscores") == 0)
      {
        print_scores = 1;
      }
      else if(strcasecmp(argv[argi], "--printfasta") == 0)
      {
        print_fasta = 1;
      }
      else if(strcasecmp(argv[argi], "--pretty") == 0)
      {
        print_pretty = 1;
      }
      else if(strcasecmp(argv[argi], "--colour") == 0)
      {
        print_colour = 1;
      }
      else if(strcasecmp(argv[argi], "--zam") == 0)
      {
        print_zam = 1;
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
      else if(strcasecmp(argv[argi], "--scoring") == 0)
      {
        if(strcasecmp(argv[argi+1], "PAM30") == 0)
        {
          scoring = scoring_system_PAM30();
        }
        else if(strcasecmp(argv[argi+1], "PAM70") == 0)
        {
          scoring = scoring_system_PAM70();
        }
        else if(strcasecmp(argv[argi+1], "BLOSUM80") == 0)
        {
          scoring = scoring_system_BLOSUM80();
        }
        else if(strcasecmp(argv[argi+1], "BLOSUM62") == 0)
        {
          scoring = scoring_system_BLOSUM62();
        }
        else if(strcasecmp(argv[argi+1], "DNA_HYBRIDIZATION") == 0)
        {
          scoring = scoring_system_DNA_hybridization();
        }
        else {
          print_usage("Unknown --scoring choice, not one of PAM30|PAM70|BLOSUM80|BLOSUM62");
        }

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--substitution_matrix") == 0)
      {
        if(scoring == NULL)
        {
          scoring = simple_scoring(match, mismatch, gap_open, gap_extend,
                                   no_start_gap_penalty, no_end_gap_penalty,
                                   case_sensitive);
        
          scoring->use_match_mismatch = 0;
        }
      
        gzFile* sub_matrix_file = gzopen(argv[argi+1], "r");
        load_matrix_scores(sub_matrix_file, scoring, case_sensitive, argv[argi+1]);
        gzclose(sub_matrix_file);

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--substitution_pairs") == 0)
      {
        if(scoring == NULL)
        {
          scoring = simple_scoring(match, mismatch, gap_open, gap_extend,
                                   no_start_gap_penalty, no_end_gap_penalty,
                                   case_sensitive);
        
          scoring->use_match_mismatch = 0;
        }
      
        gzFile* sub_pairs_file = gzopen(argv[argi+1], "r");
        load_pairwise_scores(sub_pairs_file, scoring, case_sensitive, argv[argi+1]);
        gzclose(sub_pairs_file);

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--match") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &match)) {
          print_usage("Invalid match -- must be int");
        }
        use_match = 1;
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--mismatch") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &mismatch)) {
          print_usage("Invalid mismatch score -- must be int");
        }
        use_mismatch = 1;
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--gapopen") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &gap_open)) {
          print_usage("Invalid gap open score -- must be int");
        }
        use_gap_open = 1;
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--gapextend") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &gap_extend)) {
          print_usage("Invalid gap extend score -- must be int");
        }
        use_gap_extend = 1;
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--file") == 0)
      {
        //file = gzopen(argv[argi+1], "r");
        file_path = argv[argi+1];
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
    else
    {
      seq1 = argv[argi];
      seq2 = argv[argi+1];
    }
  }

  if(seq1 == NULL && file_path == NULL && !read_stdin)
  {
    print_usage("No input specified");
  }

  if(print_zam && (print_pretty || print_scores || print_colour || print_fasta))
  {
    print_usage("Cannot use --printscore, --printfasta, --pretty or --colour with --zam");
  }

  // Set up scoring now
  if(scoring == NULL)
  {
    scoring = simple_scoring(match, mismatch, gap_open, gap_extend,
                             no_start_gap_penalty, no_end_gap_penalty,
                             case_sensitive);
  }
  else
  {
    // Update existing chosen scoring scheme
    if(use_match)
    {
      scoring->match = match;
    }
    
    if(use_mismatch)
    {
      scoring->mismatch = mismatch;
    }
    
    if(use_match != use_mismatch)
    {
      print_usage("--match --mismatch must both be set or neither set");
    }
    else if(use_match && use_mismatch)
    {
      scoring->use_match_mismatch = 1;
    }
    else {
      scoring->use_match_mismatch = 0;
    }

    
    if(use_gap_open)
    {
      scoring->gap_open = gap_open;
    }
    
    if(use_gap_extend)
    {
      scoring->gap_extend = gap_extend;
    }
    
    if(use_no_start_gap_penalty)
    {
      scoring->no_start_gap_penalty = no_start_gap_penalty;
    }
    
    if(use_no_end_gap_penalty)
    {
      scoring->no_end_gap_penalty = no_end_gap_penalty;
    }
    
    scoring->case_sensitive = case_sensitive;
  }

  // End of set up

  // Align!

  char *alignment_a = NULL, *alignment_b = NULL;
  int alignment_max_length;

  if(seq1 != NULL)
  {
    // Align seq1 and seq2
    alignment_max_length = nw_alloc_mem(seq1, seq2, &alignment_a, &alignment_b);
    
    align(seq1, seq2, alignment_a, alignment_b, NULL, NULL);
  }
  else
  {
    alignment_max_length = 1000; // equivalent to two strings of 500bp
    alignment_a = (char*) malloc((alignment_max_length+1) * sizeof(char));
    alignment_b = (char*) malloc((alignment_max_length+1) * sizeof(char));
  }

  if(file_path != NULL && strcmp(file_path, "-") != 0)
  {
    // read file
    gzFile* file = gzopen(file_path, "r");

    align_from_file(file, &alignment_a, &alignment_b, &alignment_max_length);
    
    gzclose(file);
  }

  if((file_path != NULL && strcmp(file_path, "-") == 0) || read_stdin)
  {
    // Read from STDIN
    gzFile* file = gzdopen(fileno(stdin), "r");
    
    // read file
    align_from_file(file, &alignment_a, &alignment_b, &alignment_max_length);
    
    gzclose(file);
  }
  else if(argc == 1)
  {
    print_usage(NULL);
  }

  // Free memory
  free(alignment_a);
  free(alignment_b);

  // Scoring also needs to be freed
  free_nw_scoring(scoring);

  return 0;
}
