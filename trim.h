//
// Created by john on 6/29/2017.
//

#ifndef QUICKTRIM_TRIM_H
#define QUICKTRIM_TRIM_H

#include "fastqrec.h"
#include "align.h"

int trim_se(struct fqrec* rec, unsigned int qual_cutoff, unsigned int length_cutoff, unsigned int in_a_row, unsigned char phred, unsigned char method);

int trim_pe(struct fqrec* rec1, struct fqrec* rec2, unsigned int qual_cutoff, unsigned int length_cutoff, unsigned int in_a_row, unsigned char phred, unsigned char method);

int trim_3_adapter_se(struct fqrec* rec, char* adapter, unsigned int adapterLength, unsigned int minOverlap, unsigned int minScore, struct deltas dlt);

int trim_3_adapter_pe(struct fqrec* rec1, struct fqrec* rec2, char* adapter, unsigned int adapterLength, unsigned int minOverlap, unsigned int minScore, struct deltas dlt);

int trim_5_adapter_se(struct fqrec* rec, char* adapter, unsigned int adapterLength, unsigned int minOverlap, unsigned int minScore, struct deltas dlt);

int trim_5_adapter_pe(struct fqrec* rec1, struct fqrec* rec2, char* adapter, unsigned int adapterLength, unsigned int minOverlap, unsigned int minScore, struct deltas dlt);

int trim_rev_5_adapter_pe(struct fqrec* rec1, struct fqrec* rec2, char* adapter, unsigned int adapterLength, char* adapter_rev, unsigned int adapter_revLength, unsigned int minOverlap, unsigned int minScore, struct deltas dlt, struct deltas dlt_rev);

int trim_rev_3_adapter_pe(struct fqrec* rec1, struct fqrec* rec2, char* adapter, unsigned int adapterLength, char* adapter_rev, unsigned int adapter_revLength, unsigned int minOverlap, unsigned int minScore, struct deltas dlt, struct deltas dlt_rev);

#endif //QUICKTRIM_TRIM_H
