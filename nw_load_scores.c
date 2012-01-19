/*
 nw_load_scores.c
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

void loading_error(char* err_msg, char* file_path, int line_num, char is_matrix)
{
  if(is_matrix)
  {
    fprintf(stderr, "--substitution_matrix Error: %s\n", err_msg);
  }
  else
  {
    fprintf(stderr, "--substitution_pairs Error: %s\n", err_msg);
  }
  
  if(file_path != NULL)
  {
    fprintf(stderr, "File: %s\n", file_path);
  }
  
  if(line_num != -1)
  {
    fprintf(stderr, "Line: %s\n", file_path);
  }

  exit(EXIT_FAILURE);
}

void load_matrix_scores(gzFile* file, NW_SCORING* scoring, char case_sensitive,
                        char* file_path)
{
  STRING_BUFFER* sbuf = string_buff_init(500);
  int read_length;
  int line_num = 0;

  // Read first line (column headings)
  while((read_length = string_buff_reset_zreadline(sbuf, file)) > 0)
  {
    string_buff_chomp(sbuf);

    if(sbuf->len > 0 && sbuf->buff[0] != '#' && // line is not empty, not comment
       !string_is_all_whitespace(sbuf->buff)) // and not whitespace
    {
      // Read first line

      if(sbuf->len < 2)
      {
        loading_error("Too few column headings", file_path, line_num, 0);
      }

      break;
    }

    line_num++;
  }

  if(line_num == 0 && sbuf->len <= 0)
  {
    loading_error("Empty file", file_path, -1, 0);
  }

  // If the separator character is whitespace,
  // the set of whitespace characters is used
  char sep = sbuf->buff[0];

  if((sep >= (int)'0' && sep <= (int)'9') || sep == '-')
  {
    loading_error("Numbers (0-9) and dashes (-) do not make good separators",
                  file_path, line_num, 0);
  }

  char* characters = (char*)malloc(sbuf->len);
  int num_of_chars = 0;
  
  if(isspace(sep))
  {
    char* next = sbuf->buff;

    while((next = next_nonwhitespace(next+1)) != NULL)
    {
      characters[num_of_chars++] = case_sensitive ? *next : tolower(*next);
    }

    // Now read lines below
    while((read_length = string_buff_reset_zreadline(sbuf, file)) > 0)
    {
      string_buff_chomp(sbuf);

      char* from_char_pos = next_nonwhitespace(sbuf->buff);

      if(from_char_pos == NULL || sbuf->buff[0] == '#')
      {
        // skip this line
        continue;
      }

      char from_char = case_sensitive ? *from_char_pos : tolower(*from_char_pos);
      char to_char;

      char* score_txt = sbuf->buff+1;
      int score;

      int i;
      for(i = 0; i < num_of_chars; i++)
      {
        to_char = characters[i];

        if(!isspace(*score_txt))
        {
          loading_error("Expected whitespace between elements - found character",
                        file_path, line_num, 0);
        }
      
        score_txt = next_nonwhitespace(score_txt+1);
        
        char* strtol_last_char_ptr = score_txt;
        score = (int)strtol(strtol_last_char_ptr, &strtol_last_char_ptr, 10);
  
        // If pointer to end of number string hasn't moved -> error
        if(strtol_last_char_ptr == score_txt)
        {
          loading_error("Missing number value on line", file_path, line_num, 0);
        }
        
        add_pair_to_scoring(scoring, from_char, to_char, score);
        
        score_txt = strtol_last_char_ptr;
      }
      
      if(*score_txt != '\0' && !string_is_all_whitespace(score_txt))
      {
        loading_error("Too many columns on row", file_path, line_num, 0);
      }

      line_num++;
    }
  }
  else
  {
    int i;
    for(i = 0; i < sbuf->len; i += 2)
    {
      if(sbuf->buff[i] != sep)
      {
        loading_error("Separator missing from line", file_path, line_num, 0);
      }

      char c = case_sensitive ? sbuf->buff[i+1] : tolower(sbuf->buff[i+1]);
      characters[num_of_chars++] = c;
    }
    
    int score;
    
    // Read rows
    while((read_length = string_buff_reset_zreadline(sbuf, file)) > 0)
    {
      string_buff_chomp(sbuf);

      char from_char = case_sensitive ? sbuf->buff[0] : tolower(sbuf->buff[0]);

      if(from_char == '#' || string_is_all_whitespace(sbuf->buff))
      {
        // skip this line
        continue;
      }
      
      char* str_pos = sbuf->buff;

      int to_char_index = 0;
      char to_char;

      while(*str_pos != '\0')
      {
        to_char = characters[to_char_index++];

        if(*str_pos != sep)
        {
          loading_error("Separator missing from line", file_path, line_num, 0);
        }
        
        // Move past separator
        str_pos++;

        char* after_num_str = str_pos;
        score = (int)strtol(str_pos, &after_num_str, 10);
  
        // If pointer to end of number string hasn't moved -> error
        if(str_pos == after_num_str)
        {
          loading_error("Missing number value on line", file_path, line_num, 0);
        }

        if(to_char_index >= num_of_chars)
        {
          loading_error("Too many columns on row", file_path, line_num, 0);
        }

        add_pair_to_scoring(scoring, from_char, to_char, score);

        str_pos = after_num_str;
      }
      
      line_num++;
    }
  }
}


void load_pairwise_scores(gzFile* file, NW_SCORING* scoring, char case_sensitive,
                          char* file_path)
{
  // Adds to hash table in scoring->swap_table (it needs to be already malloc'ed)

  STRING_BUFFER* sbuf = string_buff_init(200);
  int read_length;
  int line_num = 0;

  char a, b;
  int score;
  
  int num_pairs_added = 0;

  while((read_length = string_buff_reset_zreadline(sbuf, file)) > 0)
  {
    string_buff_chomp(sbuf);

    if(sbuf->len > 0 && sbuf->buff[0] != '#' && // line is not empty, not comment
       !string_is_all_whitespace(sbuf->buff)) // and not whitespace
    {
      if(read_length < 5)
      {
        loading_error("Too few column headings", file_path, line_num, 1);
      }
      
      if(isspace(sbuf->buff[1]))
      {
        // split by whitespace
        a = sbuf->buff[0];

        int char2_pos;

        for(char2_pos = 1;
            sbuf->buff[char2_pos] != '\0' && isspace(sbuf->buff[char2_pos]);
            char2_pos++);

        if(char2_pos+2 >= sbuf->len || !isspace(sbuf->buff[char2_pos+1]))
        {
          loading_error("Line too short", file_path, line_num, 0);
        }
        
        b = sbuf->buff[char2_pos];

        parse_entire_int(sbuf->buff+char2_pos+2, &score);
      }
      else
      {
        if(sbuf->buff[1] != sbuf->buff[3])
        {
          loading_error("Inconsistent separators used", file_path, line_num, 0);
        }
        
        a = sbuf->buff[0];
        b = sbuf->buff[2];
        parse_entire_int(sbuf->buff + 4, &score);
      }

      if(!case_sensitive)
      {
        a = tolower(a);
        b = tolower(b);
      }

      add_pair_to_scoring(scoring, a, b, score);
      num_pairs_added++;
    }
    
    line_num++;
  }

  if(num_pairs_added == 0)
  {
    loading_error("No pairs added from file (file empty?)",
                  file_path, line_num, 0);
  }
}
