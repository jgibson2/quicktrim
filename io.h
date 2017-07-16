//
// Created by john on 6/30/2017.
//

#ifndef QUICKTRIM_IO_H
#define QUICKTRIM_IO_H

#include <stdio.h>
#include "fastqrec.h"

#define RECORD_LENGTH_DIVIDER 20

/*
struct fqfiledata {
    FILE* file;
    char* buf;
    unsigned long old_start;
    unsigned long old_end;
    unsigned long new_start;
    unsigned long new_end;
    unsigned long current_offset;
    unsigned long buf_size;
    int eof;
    char next();
};
*/

struct fqfiledata {
    FILE* file;
    char* buf;
    unsigned long current_offset;
    unsigned long buf_size;
    unsigned long current_size;
    int eof;

};
int init_file_data(unsigned long size, const char* filename, struct fqfiledata* data);
int writeRecord(struct fqrec* rec, FILE* file);
int writePairedRecords(struct fqrec* rec1, FILE* file1, struct fqrec* rec2, FILE* file2);

int getNextRecord(struct fqfiledata* data, struct fqrec* rec);

#endif //QUICKTRIM_IO_H
