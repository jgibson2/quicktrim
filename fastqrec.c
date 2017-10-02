//
// Created by john on 6/29/2017.
//
#include <stdio.h>
#include "fastqrec.h"

//TODO: do NOT copy memory into other buffer!

int allocatefqrec(char* bytes, unsigned long name_start_offset, unsigned long name_end_offset,
    unsigned long seq_start_offset, unsigned long seq_end_offset, unsigned long qual_start_offset,
    unsigned long qual_end_offset, struct fqrec* rec)
    /*
     * offsets should be in terms of the bytes array and in [start, end) format
     * Returns:
     *      0 if success
     *      1 if qual and seq are not the same length
     *      2 if memory allocation fails
     */
{

    if(qual_end_offset - qual_start_offset != seq_end_offset - seq_start_offset)
    {
        return 1;
    }

    //allocate space for name and copy it into the name buffer
    rec -> nameLength = name_end_offset - name_start_offset + 1;
    rec->name = bytes + name_start_offset;

    //calculate offset (value to add in order to get to corresponding qual score from base
    rec->seqLength = seq_end_offset-seq_start_offset;

    //use aligned_alloc to allocate memory for sequence and qual in the same buffer
    rec->seq = bytes + seq_start_offset;
    rec->originalSeq = rec->seq;

    rec->qual = bytes + qual_start_offset;
    rec->originalQual = rec->qual;

    return 0;
}

int reallocatefqrec(char* bytes, unsigned long name_start_offset, unsigned long name_end_offset,
    unsigned long seq_start_offset, unsigned long seq_end_offset, unsigned long qual_start_offset,
    unsigned long qual_end_offset, struct fqrec* rec)
/*
 * offsets should be in terms of the bytes array and in [start, end) format
 * Returns:
 *      0 if success
 *      1 if qual and seq are not the same length
 *      2 if memory allocation fails
 */
{
    if(qual_end_offset - qual_start_offset != seq_end_offset - seq_start_offset)
    {
        return 1;
    }

    //allocate space for name and copy it into the name buffer
    rec->nameLength = name_end_offset - name_start_offset + 1;
    rec->name = bytes + name_start_offset;

    //calculate offset (value to add in order to get to corresponding qual score from base
    rec->seqLength = seq_end_offset-seq_start_offset;

    //use aligned_alloc to allocate memory for sequence and qual in the same buffer
    rec->seq = bytes + seq_start_offset;
    rec->originalSeq = rec->seq;

    rec->qual = bytes + qual_start_offset;
    rec->originalQual = rec->qual;

    return 0;
}

void freefqrec(struct fqrec* rec)
{
    aligned_free(rec->originalSeq);
    aligned_free(rec->originalQual);
}