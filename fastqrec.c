//
// Created by john on 6/29/2017.
//
#include <stdio.h>
#include "fastqrec.h"

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
    rec->name = (char*)calloc((rec->nameLength), sizeof(char));
    if(rec->name == 0)
    {
        return 2;
    }
    strncpy(rec->name, bytes + name_start_offset, name_end_offset - name_start_offset);
    rec->name[name_end_offset - name_start_offset] = '\0';

    //calculate offset (value to add in order to get to corresponding qual score from base
    rec->offset = seq_end_offset-seq_start_offset;

    //use aligned_alloc to allocate memory for sequence and qual in the same buffer
    rec->seqAndQualBuf = (char*)aligned_malloc(sizeof(char) * (2*(rec->offset) + 1), sizeof(char));
    if(rec->seqAndQualBuf == 0)
    {
        return 2;
    }
    strncpy(rec->seqAndQualBuf, bytes + seq_start_offset, seq_end_offset - seq_start_offset);
    strncpy(rec->seqAndQualBuf + rec->offset, bytes + qual_start_offset, qual_end_offset - qual_start_offset);
    rec->seqAndQualBuf[2 * (rec->offset)] = '\0';
    rec->originalBuf = rec->seqAndQualBuf;
    rec->seqLength = rec->offset;

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
    if(name_end_offset - name_start_offset + 1 != rec->nameLength) {
        rec->nameLength = name_end_offset - name_start_offset + 1;
        free(rec->name);
        rec->name = (char *) calloc((rec->nameLength), sizeof(char));
        if(rec->name == 0)
        {
            return 2;
        }
    }
    strncpy(rec->name, bytes + name_start_offset, name_end_offset - name_start_offset);
    rec->name[name_end_offset - name_start_offset] = '\0';

    //calculate offset (value to add in order to get to corresponding qual score from base
    if(seq_end_offset-seq_start_offset != rec->offset || rec->originalBuf != rec->seqAndQualBuf) {
        rec->offset = seq_end_offset - seq_start_offset;

        //use aligned_alloc to allocate memory for sequence and qual in the same buffer
        aligned_free(rec->originalBuf);
        rec->seqAndQualBuf = (char*)aligned_malloc(sizeof(char) * (2*(rec->offset) + 1), sizeof(char));
        if (rec->seqAndQualBuf == 0) {
            return 2;
        }
    }
    strncpy(rec->seqAndQualBuf, bytes + seq_start_offset, seq_end_offset - seq_start_offset);
    strncpy(rec->seqAndQualBuf + rec->offset, bytes + qual_start_offset, qual_end_offset - qual_start_offset);
    rec->seqAndQualBuf[2 * (rec->offset)] = '\0';
    rec->originalBuf = rec->seqAndQualBuf;
    rec->seqLength = rec->offset;

    return 0;
}


int freefqrec(struct fqrec* rec)
    /*
     * Frees memory associated with rec
     */
{

    free(rec->name);
    //free(&(rec->offset));
    aligned_free(rec->originalBuf);
    //free(&(rec->seqLength));
    //free(&(rec->nameLength));
    free(rec);
    return 0;
}

/*
 * https://stackoverflow.com/questions/16376942/best-cross-platform-method-to-get-aligned-memory
 */

