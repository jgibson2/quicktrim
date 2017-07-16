//
// Created by john on 6/29/2017.
//

#ifndef QUICKTRIM_FASTQREC_H
#define QUICKTRIM_FASTQREC_H

#include <stdlib.h>
#include <string.h>
#include "mem.h"


//TODO: offsets based on reference to char buf instead of memory? But then it gets rid of aligned memory... hmm...
struct fqrec
{
    char* name;
    char* seqAndQualBuf;
    unsigned int offset;
    unsigned int nameLength;
    unsigned int seqLength;
};

int allocatefqrec(char* bytes, unsigned long name_start_offset, unsigned long name_end_offset,
    unsigned long seq_start_offset, unsigned long seq_end_offset, unsigned long qual_start_offset,
    unsigned long qual_end_offset, struct fqrec* rec);

int reallocatefqrec(char* bytes, unsigned long name_start_offset, unsigned long name_end_offset,
    unsigned long seq_start_offset, unsigned long seq_end_offset, unsigned long qual_start_offset,
    unsigned long qual_end_offset, struct fqrec* rec);

int freefqrec(struct fqrec* rec);

#endif //QUICKTRIM_FASTQREC_H
