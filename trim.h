//
// Created by john on 6/29/2017.
//

#ifndef QUICKTRIM_TRIM_H
#define QUICKTRIM_TRIM_H

#include "fastqrec.h"

int trim_se(struct fqrec* rec, unsigned int qual_cutoff, unsigned int length_cutoff, unsigned int in_a_row, int phred);

int trim_pe(struct fqrec* rec1, struct fqrec* rec2, unsigned int qual_cutoff, unsigned int length_cutoff,
            unsigned int in_a_row, int phred);

#endif //QUICKTRIM_TRIM_H
