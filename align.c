//
// Created by john on 7/29/2017.
//
#include <stdio.h>
#include "align.h"

unsigned int check_if_adapter_at_end(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap);

unsigned int check_if_adapter_at_start(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap);

unsigned int find_matrix_start(unsigned int i, struct AlignmentMatrix* mat);

int max4(int a, int b, int c, int d);

unsigned int get_3_adapter_start_position(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap, unsigned int minScore)
{
    //work out return types here... maybe take a record and set seqLength.
    struct AlignmentMatrix mat;

    //check if seq present at end
    unsigned int res = check_if_adapter_at_end(seq, seqLength, adapter, adapterLength, minOverlap);
    if(res)
    {
        return res;
    }

    res = build_end_alignment_matrix(seq, seqLength, adapter, adapterLength, minOverlap, &mat);
    if(res || mat.score < minScore)
    {
        return 0;
    }

    unsigned int start_pos = mat.start;
    free_matrix(&mat);
    return start_pos;

}

unsigned int get_5_adapter_end_position(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap, unsigned int minScore)
{
    //work out return types here... maybe take a record and set seqLength.
    struct AlignmentMatrix mat;

    //check if seq present at end
    unsigned int res = check_if_adapter_at_start(seq, seqLength, adapter, adapterLength, minOverlap);
    if(res)
    {
        return res;
    }

    res = build_start_alignment_matrix(seq, seqLength, adapter, adapterLength, minOverlap, &mat);
    if(res || mat.score < minScore)
    {
        return 0;
    }

    unsigned int end_pos = mat.end;
    free_matrix(&mat);
    return end_pos;

}

unsigned int build_end_alignment_matrix(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap, struct AlignmentMatrix* mat)
{
    //Modification of the Smith-Waterman Algorithm to match the adapter to the end of the read

    //build the matrix such that the adapter sequence can overhang (but also such that minOverlap bases must align)
    int* matrix = aligned_malloc(sizeof(int) * (N*adapterLength - minOverlap + 1) * sizeof(int) * (adapterLength + 1), sizeof(int));
    //matrix = memset(matrix, 0, sizeof(int) * (N*adapterLength - minOverlap + 1) * sizeof(int) * (adapterLength + 1));
    if(matrix == NULL)
    {
        return 1;
    }

    //set the first row and column to zero. In the case of full alignment this will be our stopping point
    for(unsigned int i = 0; i < (N*adapterLength - minOverlap + 1); i++)
    {
        matrix[i] = 0;
    }
    for(unsigned int i = 0; i < (N*adapterLength - minOverlap + 1) * (adapterLength + 1); i+=(N*adapterLength - minOverlap + 1))
    {
        matrix[i] = 0;
    }

    //track the position of the max score so we can get it later
    int max = 0;
    unsigned int max_position = 0;

    //fill in the rest of the matrix using linear gap penalties and a mismatch penalty of one (for speed --
    //  getting the best alignment isn't as important here)
    for(unsigned int i = 1; i < adapterLength + 1; i++) {
        for (unsigned int j = 1; j < (N*adapterLength - minOverlap + 1); j++) {
            int score = matrix[((i - 1) * (N*adapterLength - minOverlap + 1)) + (j - 1)];
            if(j + (seqLength - (N*adapterLength - minOverlap + 1)) < seqLength)
            {
                score += seq[j + (seqLength - (N*adapterLength - minOverlap + 1))] == adapter[i - 1] ? MATCH_BONUS : -1*MISMATCH_PENALTY;
                //printf("%c", seq[j + (seqLength - (2 * adapterLength - minOverlap + 1)) - 1]);
            }
            int gap1 = matrix[((i) * (N*adapterLength - minOverlap + 1)) + (j - 1)] - GAP_PENALTY;
            int gap2 = matrix[((i - 1) * (N*adapterLength - minOverlap + 1)) + (j)] - GAP_PENALTY;

            matrix[((i) * (N*adapterLength - minOverlap + 1)) + (j)] = max4(score, gap1, gap2, 0);

            if(matrix[((i) * (N*adapterLength - minOverlap + 1)) + (j)] > max)
            {
                max = matrix[((i) * (N*adapterLength - minOverlap + 1)) + (j)];
                max_position = ((i) * (N*adapterLength - minOverlap + 1)) + (j);
            }
        }
        //printf("\n");
    }

    /*
      //print matrix for debugging
    printf("\n");
    for(unsigned int i = 0; i < adapterLength + 1; i++) {
        for (unsigned int j = 0; j < (N*adapterLength - minOverlap + 1); j++) {
            printf("%d|", matrix[((i) * (N*adapterLength - minOverlap + 1)) + (j)]);
        }
        printf("\n");
    }
    printf("\n");
    */

    mat->seqLength = seqLength;
    mat->seq = seq;
    mat->adapterLength = adapterLength;
    mat->adapter = adapter;
    mat->matrix = matrix;
    mat->minOverlap = minOverlap;
    mat->end = max_position % (N*adapterLength - minOverlap + 1) + (seqLength - (N*adapterLength - minOverlap));
    mat->start = find_matrix_start(max_position, mat) % (N*adapterLength - minOverlap + 1) + (seqLength - (N*adapterLength - minOverlap)); //start position in read
    mat->score = (unsigned int)max;

    //printf("%d\t%d\t%d\n", mat->start, mat->end, mat->score);

    return 0;
}

unsigned int build_start_alignment_matrix(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap, struct AlignmentMatrix* mat)
{
    //Modification of the Smith-Waterman Algorithm to match the adapter to the start of the read

    //build the matrix such that the adapter sequence can overhang (but also such that minOverlap bases must align)
    int* matrix = aligned_malloc(sizeof(int) * (N*adapterLength - minOverlap + 1) * sizeof(int) * (adapterLength + 1), sizeof(int));
    //matrix = memset(matrix, 0, sizeof(int) * (N*adapterLength - minOverlap + 1) * sizeof(int) * (adapterLength + 1));
    if(matrix == NULL)
    {
        return 1;
    }

    //set the first row and column to zero. In the case of full alignment this will be our stopping point
    for(unsigned int i = 0; i < (N*adapterLength - minOverlap + 1); i++)
    {
        matrix[i] = 0;
    }
    for(unsigned int i = 0; i < (N*adapterLength - minOverlap + 1) * (adapterLength + 1); i+=(N*adapterLength - minOverlap + 1))
    {
        matrix[i] = 0;
    }

    //track the position of the max score so we can get it later
    int max = 0;
    unsigned int max_position = 0;

    //fill in the rest of the matrix using linear gap penalties and a mismatch penalty of one (for speed --
    //  getting the best alignment isn't as important here)
    for(unsigned int i = 1; i < adapterLength + 1; i++) {
        for (unsigned int j = 1; j < (N*adapterLength - minOverlap + 1); j++) {
            int score = matrix[((i - 1) * (N*adapterLength - minOverlap + 1)) + (j - 1)];
            if((adapterLength - minOverlap + 1) - j > 0)
            {
                score += seq[j - (adapterLength - minOverlap + 1) - 1] == adapter[i - 1] ? MATCH_BONUS : -1*MISMATCH_PENALTY;
            }
            int gap1 = matrix[((i) * (N*adapterLength - minOverlap + 1)) + (j - 1)] - GAP_PENALTY;
            int gap2 = matrix[((i - 1) * (N*adapterLength - minOverlap + 1)) + (j)] - GAP_PENALTY;

            matrix[((i) * (N*adapterLength - minOverlap + 1)) + (j)] = max4(score, gap1, gap2, 0);

            if(matrix[((i) * (N*adapterLength - minOverlap + 1)) + (j)] > max)
            {
                max = matrix[((i) * (N*adapterLength - minOverlap + 1)) + (j)];
                max_position = ((i) * (N*adapterLength - minOverlap + 1)) + (j);
            }
        }
        //printf("\n");
    }

    /*
      //print matrix for debugging
    printf("\n");
    for(unsigned int i = 0; i < adapterLength + 1; i++) {
        for (unsigned int j = 0; j < (N*adapterLength - minOverlap + 1); j++) {
            printf("%d|", matrix[((i) * (N*adapterLength - minOverlap + 1)) + (j)]);
        }
        printf("\n");
    }
    printf("\n");
     */


    mat->seqLength = seqLength;
    mat->seq = seq;
    mat->adapterLength = adapterLength;
    mat->adapter = adapter;
    mat->matrix = matrix;
    mat->minOverlap = minOverlap;
    mat->end = (max_position % (N*adapterLength - minOverlap + 1)) - (adapterLength - minOverlap); //end noninclusive; this is start of the rest of the sequence
    //start not needed
    //mat->start = find_matrix_start(max_position, mat) % (N*adapterLength - minOverlap + 1) + (seqLength - (N*adapterLength - minOverlap + 1)); //start position in read
    mat->score = (unsigned int)max;

    //printf("%d\t%d\t%d\n", 0, mat->end, mat->score);

    return 0;
}

inline unsigned int check_if_adapter_at_end(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap)
{
    int found;
    for(unsigned int i = seqLength - N*adapterLength; i <= seqLength - minOverlap; i++)
    {
        found = 1;
        for(unsigned int j = 0; j < adapterLength && i + j < seqLength; j++)
        {
            if(adapter[j] != seq[i + j])
            {
                found = 0;
                break;
            }
        }
        if(found)
        {
            return i;
        }
    }
    return 0;
}

inline unsigned int check_if_adapter_at_start(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap)
{
    int found;
    for(unsigned int i = adapterLength - minOverlap; i > 0; i--)
    {
        found = 1;
        for(unsigned int j = 0; (i-1)+j < adapterLength && j < seqLength; j++)
        {
            if(adapter[(i-1)+j] != seq[j])
            {
                found = 0;
                break;
            }
        }
        if(found)
        {
            return adapterLength - (i-1);
        }
    }
    for(unsigned int i = 0; i < (N-1)*adapterLength; i++)
    {
        found = 1;
        for(unsigned int j = 0; j < adapterLength && i+j < seqLength; j++)
        {
            if(adapter[j] != seq[i+j])
            {
                found = 0;
                break;
            }
        }
        if(found)
        {
            return i+adapterLength; //next nt is start
        }
    }

    return 0;
}

inline unsigned int find_matrix_start(unsigned int i, struct AlignmentMatrix* mat)
{
    if(mat->matrix[i] == 0)
    {
        return i;
    }
    //check the cells to the left, upper left, and top of the current cell in the Smith-Waterman matrix
    int max = mat->matrix[i-1];
    unsigned int max_position = i-1;
    if(mat->matrix[i-(N*mat->adapterLength - mat->minOverlap + 1)] > max)
    {
        max = mat->matrix[i-(N*mat->adapterLength - mat->minOverlap + 1)];
        max_position = i-(N*mat->adapterLength - mat->minOverlap + 1);
    }
    if(mat->matrix[i-(N*mat->adapterLength - mat->minOverlap + 1)-1] > max)
    {
        //max = mat->matrix[i-(2*mat->adapterLength - mat->minOverlap + 1)-1]; //unneeded
        max_position = i-(N*mat->adapterLength - mat->minOverlap + 1)-1;
    }

    //recurse to find the start point
    return find_matrix_start(max_position, mat);
}

inline int max4(int a, int b, int c, int d)
{
    int max = a;
    max = b > max ? b : max;
    max = c > max ? c : max;
    max = d > max ? d : max;
    return max;
}

void free_matrix(struct AlignmentMatrix* mat)
{
    aligned_free(mat->matrix);
    free(mat);
}

/*
 * Prefix Array/tree using Kmers to infer best position, then checking any mismatches?
 */

