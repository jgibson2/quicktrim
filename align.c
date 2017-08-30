//
// Created by john on 7/29/2017.
//
#include <stdio.h>
#include "align.h"

inline unsigned int check_if_adapter_at_end(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap, int* delta1, int* delta2);

inline unsigned int check_if_adapter_at_start(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap, int* delta1, int* delta2);

unsigned int find_matrix_start(unsigned int i, struct AlignmentMatrix* mat);

int max4(int a, int b, int c, int d);

unsigned int get_3_adapter_start_position(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap, unsigned int minScore, struct deltas dlt)
/*
 * Gets the 3' adapter start position for trimming
 *  seq: sequence
 *  seqLength: sequence length
 *  adapter: 3' adapter
 *  adapterLength: adapter length
 *  minOverlap: minimum overlap
 *  minScore: minimum alignment score
 */
{
    //work out return types here... maybe take a record and set seqLength.
    struct AlignmentMatrix mat;

    //check if seq present at end
    unsigned int res = check_if_adapter_at_end(seq, seqLength, adapter, adapterLength, minOverlap, dlt.delta1, dlt.delta2);
    if(res)
    {
        //printf("Here!\n");
        return res;
    }

    res = build_end_alignment_matrix(seq, seqLength, adapter, adapterLength, minOverlap, &mat);
    if(res || mat.score < minScore)
    {
        return 0;
    }

    unsigned int start_pos = mat.start;
    //free_matrix(&mat); mat is on stack, so it does not need to be freed -- necessarily -- though the matrix does
    aligned_free(mat.matrix);
    return start_pos;

}

unsigned int get_5_adapter_end_position(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap, unsigned int minScore, struct deltas dlt)
/*
 * Gets the 5' adapter end position for trimming
 *  seq: sequence
 *  seqLength: sequence length
 *  adapter: 5' adapter
 *  adapterLength: adapter length
 *  minOverlap: minimum overlap
 *  minScore: minimum alignment score
 */
{
    //work out return types here... maybe take a record and set seqLength.
    struct AlignmentMatrix mat;

    //check if seq present at end
    unsigned int res = check_if_adapter_at_start(seq, seqLength, adapter, adapterLength, minOverlap, dlt.delta1, dlt.delta2);
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
    //free_matrix(&mat); mat is on stack, so it does not need to be freed -- necessarily -- though the matrix does
    aligned_free(mat.matrix);
    return end_pos;

}

unsigned int build_end_alignment_matrix(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap, struct AlignmentMatrix* mat)
/*
 * Builds end alignment matrix
 *  seq: sequence
 *  seqLength: sequence length
 *  adapter: 5' adapter
 *  adapterLength: adapter length
 *  minOverlap: minimum overlap
 *  mat: AlignmentMatrix struct
 */
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
    mat->start = find_matrix_start(max_position, mat) % (N*adapterLength - minOverlap + 1) + (seqLength - (N*adapterLength - minOverlap)) + 1; //start position in read
    mat->score = (unsigned int)max;

    //printf("%d\t%d\t%d\n", mat->start, mat->end, mat->score);

    return 0;
}

unsigned int build_start_alignment_matrix(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap, struct AlignmentMatrix* mat)
/*
 * Builds start alignment matrix
 *  seq: sequence
 *  seqLength: sequence length
 *  adapter: 5' adapter
 *  adapterLength: adapter length
 *  minOverlap: minimum overlap
 *  mat: AlignmentMatrix struct
 */
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

inline unsigned int find_matrix_start(unsigned int i, struct AlignmentMatrix* mat)
/*
 * Finds the matrix start
 *  i: position
 *  mat: struct AlignmentMatrix
 */
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
 * Boyer-Moore inspired algorithms to check for adapter at start or end
 * https://en.wikipedia.org/wiki/Boyer%E2%80%93Moore_string_search_algorithm
 */
void make_delta1(int *delta1, char *pat, int patlen);
void make_delta2(int *delta2, char *pat, int patlen);

void make_deltas(struct deltas* dlt, char* pat, int patlen)
{
    make_delta1(dlt->delta1, pat, patlen);
    make_delta2(dlt->delta2, pat, patlen);
}

// delta1 table: delta1[c] contains the distance between the last
// character of pat and the rightmost occurrence of c in pat.
// If c does not occur in pat, then delta1[c] = patlen.
// If c is at string[i] and c != pat[patlen-1], we can
// safely shift i over by delta1[c], which is the minimum distance
// needed to shift pat forward to get string[i] lined up
// with some character in pat.
// this algorithm runs in alphabet_len+patlen time.
void make_delta1(int *delta1, char *pat, int patlen) {
    int i;
    for (i=0; i < ALPHABET_LEN; i++) {
        delta1[i] = NOT_FOUND;
    }
    for (i=0; i < patlen-1; i++) {
        delta1[(unsigned int)pat[i]] = patlen-1 - i;
    }
}

// true if the suffix of word starting from word[pos] is a prefix
// of word
int is_prefix(char *word, int wordlen, int pos) {
    int i;
    int suffixlen = wordlen - pos;
    // could also use the strncmp() library function here
    for (i = 0; i < suffixlen; i++) {
        if (word[i] != word[pos+i]) {
            return 0;
        }
    }
    return 1;
}

// length of the longest suffix of word ending on word[pos].
// suffix_length("dddbcabc", 8, 4) = 2
int suffix_length(char *word, int wordlen, int pos) {
    int i;
    // increment suffix length i to the first mismatch or beginning
    // of the word
    for (i = 0; (word[pos-i] == word[wordlen-1-i]) && (i < pos); i++);
    return i;
}

// delta2 table: given a mismatch at pat[pos], we want to align
// with the next possible full match could be based on what we
// know about pat[pos+1] to pat[patlen-1].
//
// In case 1:
// pat[pos+1] to pat[patlen-1] does not occur elsewhere in pat,
// the next plausible match starts at or after the mismatch.
// If, within the substring pat[pos+1 .. patlen-1], lies a prefix
// of pat, the next plausible match is here (if there are multiple
// prefixes in the substring, pick the longest). Otherwise, the
// next plausible match starts past the character aligned with
// pat[patlen-1].
//
// In case 2:
// pat[pos+1] to pat[patlen-1] does occur elsewhere in pat. The
// mismatch tells us that we are not looking at the end of a match.
// We may, however, be looking at the middle of a match.
//
// The first loop, which takes care of case 1, is analogous to
// the KMP table, adapted for a 'backwards' scan order with the
// additional restriction that the substrings it considers as
// potential prefixes are all suffixes. In the worst case scenario
// pat consists of the same letter repeated, so every suffix is
// a prefix. This loop alone is not sufficient, however:
// Suppose that pat is "ABYXCDBYX", and text is ".....ABYXCDEYX".
// We will match X, Y, and find B != E. There is no prefix of pat
// in the suffix "YX", so the first loop tells us to skip forward
// by 9 characters.
// Although superficially similar to the KMP table, the KMP table
// relies on information about the beginning of the partial match
// that the BM algorithm does not have.
//
// The second loop addresses case 2. Since suffix_length may not be
// unique, we want to take the minimum value, which will tell us
// how far away the closest potential match is.
void make_delta2(int *delta2, char *pat, int patlen) {
    int p;
    int last_prefix_index = patlen-1;

    // first loop
    for (p=patlen-1; p>=0; p--) {
        if (is_prefix(pat, patlen, p+1)) {
            last_prefix_index = p+1;
        }
        delta2[p] = last_prefix_index + (patlen-1 - p);
    }

    // second loop
    for (p=0; p < patlen-1; p++) {
        int slen = suffix_length(pat, patlen, p);
        if (pat[p - slen] != pat[patlen-1 - slen]) {
            delta2[patlen-1 - slen] = patlen-1 - p + slen;
        }
    }
}

inline unsigned int check_if_adapter_at_end(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap, int* delta1, int* delta2)
/*
 * Checks for adapter at the end of the read
 *  seq: sequence
 *  seqLength: sequence length
 *  adapter: 5' adapter
 *  adapterLength: adapter length
 *  minOverlap: minimum overlap
 */
{
    //TODO: not working (da fuq?)
    int found;
    for(unsigned int i = seqLength - N*adapterLength; i <= seqLength - minOverlap;)
    {
        unsigned int j;
        found = 1;
        for(j = 0; j < adapterLength && i + j < seqLength; j++)
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

        //wtf? BM not working...
        /*
        printf("%d\n", i);
        printf("%d\t%d\n",delta1[(unsigned int)seq[i+j]], delta2[j]);
        i += max(delta1[(unsigned int)seq[i+j]], delta2[adapterLength-j-1]) - j;
        printf("New: %d\n", i);
        */

        i += 1;
    }
    return 0;
}

inline unsigned int check_if_adapter_at_start(char* seq, unsigned int seqLength, char* adapter, unsigned int adapterLength, unsigned int minOverlap, int* delta1, int* delta2)
/*
 * Checks for adapter at the start of the read
 *  seq: sequence
 *  seqLength: sequence length
 *  adapter: 5' adapter
 *  adapterLength: adapter length
 *  minOverlap: minimum overlap
 */
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
    for(unsigned int i = 0; i < (N-1)*adapterLength;)
    {
        found = 1;
        unsigned int j;
        for(j = 0; j < adapterLength && i+j < seqLength; j++)
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
        //i += max(delta1[(unsigned int)seq[i+j]], delta2[j]);
        i += 1;
    }

    return 0;
}