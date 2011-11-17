/*
 needleman_wunsch.h
 project: NeedlemanWunsch
 author: Isaac Turner <isaac.turner@dtc.ox.ac.uk>
 Copyright (C) 25-May-2011
 
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


#ifndef NEEDLEMAN_WUNSCH_HEADER_SEEN
#define NEEDLEMAN_WUNSCH_HEADER_SEEN

const int nw_match_penalty, nw_mismatch_penalty,
          nw_gap_open_penalty, nw_gap_extend_penalty;

int nw_simple_sub_matrix[4][4];

// construct substitution matrix with match/mismatch scores
void nw_create_sub_matrix(int match, int mismatch, int (*matrix)[4][4]);
//int[4][4] nw_create_sub_matrix(int match, int mismatch);

/* Allocate memory for result */

// alloc memory for result (returns length of seq_a + seq_b)
int nw_alloc_mem(const char* seq_a, const char* seq_b,
                 char** alignment_a, char** alignment_b);

// length is = length_a + length_b
int nw_realloc_mem(unsigned int length, char** alignment_a, char** alignment_b);

/* Vanilla Needleman Wunsch */

int needleman_wunsch(const char* seq_a, const char* seq_b,
                     char* result_a, char* result_b,
                     int similarity_matrix[4][4],
                     const int gap_penalty);

int score_alignment(const char* alignment_a, const char* alignment_b,
                    int similarity_matrix[4][4],
                    const int gap_penalty);

/* Affine Gap Penalty */

int needleman_wunsch_affine(const char* seq_a, const char* seq_b,
                            char* result_a, char* result_b,
                            int similarity_matrix[4][4],
                            const int gap_penalty_start,
                            const int gap_penalty_cont);

int score_alignment_affine(const char* alignment_a, const char* alignment_b,
                           int similarity_matrix[4][4],
                           const int gap_penalty_start,
                           const int gap_penalty_cont);

// Private

int _get_base_index(const char b);

#endif
