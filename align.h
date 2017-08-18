//
// Created by john on 7/29/2017.
//

#ifndef QUICKTRIM_ALIGN_H
#define QUICKTRIM_ALIGN_H

#include "mem.h"
#include <string.h>

#define N 2 //controls size of matrix, N * adapterSize - adapterLength + 1 in one dimension
#define GAP_PENALTY 3
#define MISMATCH_PENALTY 2
#define MATCH_BONUS 1

struct AlignmentMatrix
{
    char* seq;
    unsigned int seqLength;
    char* adapter;
    unsigned int adapterLength;
    unsigned int minOverlap;
    int* matrix;
    unsigned int start;
    unsigned int end;
    unsigned int score;
};

unsigned int get_3_adapter_start_position(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap, unsigned int maxErrors);

unsigned int build_end_alignment_matrix(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap, struct AlignmentMatrix* mat);

unsigned int build_start_alignment_matrix(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap, struct AlignmentMatrix* mat);

unsigned int get_5_adapter_end_position(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap, unsigned int minScore);

/*
inline unsigned int check_if_adapter_at_end(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap);

inline unsigned int find_matrix_start(unsigned int i, struct AlignmentMatrix* mat);

inline int max4(int a, int b, int c, int d);
 */

void free_matrix(struct AlignmentMatrix* mat);

#endif //QUICKTRIM_ALIGN_H
