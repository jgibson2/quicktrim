//
// Created by john on 6/29/2017.
//

#include "trim.h"


int trim_se(struct fqrec *rec, unsigned int qual_cutoff, unsigned int length_cutoff, unsigned int in_a_row,
            unsigned char phred, unsigned char method)
/*
 * Single end trimming
 *  Parameters:
 *      rec: record to trim
 *      qual_cutoff: quality cutoff for bases
 *      length_cutoff: length cutoff for the read. If the length goes below this number, seqLength drops to 0
 *      in_a_row: number of high-quality bases in a row to stop filtering
 *      phred : phred version
 */
{
    if (method == 0) {
        //algorithm based on minimum cumulative quality score
        int sum = 0;
        int min = 0;
        unsigned int min_position = rec->seqLength - 1;
        for (unsigned int i = rec->seqLength - 1; i > 0; i--) {
            if (rec->seqLength < length_cutoff) {
                rec->seqLength = 0;
                return 1;
            }
            sum += (rec->seqAndQualBuf[i + rec->offset] - phred) - qual_cutoff;
            if (sum < min) {
                min_position = i;
                min = sum;
            } else if (sum > 0) {
                break;
            }
        }
        rec->seqLength = min_position;
        return 0;
    } else
        //algorithm based on # of HQ bases in a row
    {
        unsigned int hq_bases = 0;
        for (unsigned int i = rec->seqLength - 1; i > 0; i--) {
            if (rec->seqLength < length_cutoff) {
                rec->seqLength = 0;
                return 1;
            }
            if ((unsigned char)rec->seqAndQualBuf[i + rec->offset] > qual_cutoff + phred) { //equivalent form to subtracting phred from left-hand side but protects signed-ness
                if (hq_bases == in_a_row) {
                    rec->seqLength += in_a_row;
                    return 0;
                }
                hq_bases++;
            } else if (hq_bases > 0) {
                hq_bases = 0;
            }
            rec->seqLength--;
        }
        return 2;
    }
}

int trim_pe(struct fqrec* rec1, struct fqrec* rec2, unsigned int qual_cutoff, unsigned int length_cutoff, unsigned int in_a_row, unsigned char phred, unsigned char method)
/*
 *  Paired-end trimming
 *  Parameters:
 *      rec1, rec2: records to trim
 *      qual_cutoff: quality cutoff for bases
 *      length_cutoff: length cutoff for the read. If the length goes below this number, seqLength drops to 0
 *      in_a_row: number of high-quality bases in a row to stop filtering
 *      phred : phred version
 */
{
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            trim_se(rec1, qual_cutoff, length_cutoff, in_a_row, phred, method);
        }

        #pragma omp section
        {
            trim_se(rec2, qual_cutoff, length_cutoff, in_a_row, phred, method);
        }
    }
    if (rec1->seqLength == 0 || rec2->seqLength == 0) {
        return 1;
    }
    return 0;
}

int trim_3_adapter_se(struct fqrec* rec, char* adapter, unsigned int adapterLength, unsigned int minOverlap, unsigned int minScore, struct deltas dlt)
{
    unsigned int start = get_3_adapter_start_position(rec->seqAndQualBuf, rec->seqLength, adapter, adapterLength, minOverlap, minScore, dlt);
    if(start == 0)
    {
        return 1;
    }
    rec->seqLength = start; //modify length to not include adapter
    return 0;
}

int trim_3_adapter_pe(struct fqrec* rec1, struct fqrec* rec2, char* adapter, unsigned int adapterLength, unsigned int minOverlap, unsigned int minScore, struct deltas dlt)
/*
 *  Paired-end trimming
 *  Parameters:
 *      most self-explanatory
 *      minScore: minimum alignment score in Smith-Waterman matrix for adapter to be considered found.
 */
{
    int res1, res2;
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            res1 = trim_3_adapter_se(rec1, adapter, adapterLength, minOverlap, minScore, dlt);
        }

        #pragma omp section
        {
            res2 = trim_3_adapter_se(rec2, adapter, adapterLength, minOverlap, minScore, dlt);
        }
    }
    return res1 + res2;
}

int trim_5_adapter_se(struct fqrec* rec, char* adapter, unsigned int adapterLength, unsigned int minOverlap, unsigned int minScore, struct deltas dlt)
{
    unsigned int end = get_5_adapter_end_position(rec->seqAndQualBuf, rec->seqLength, adapter, adapterLength, minOverlap, minScore, dlt);
    if(end == 0)
    {
        return 1;
    }
    rec->seqAndQualBuf = rec->seqAndQualBuf + end; //modify seq buf to not include adapter
    rec->seqLength = rec->seqLength - end; //modify seqLength to not include buffer

    //printf("%s\t%d\n", rec->name, end);

    return 0;
}

int trim_5_adapter_pe(struct fqrec* rec1, struct fqrec* rec2, char* adapter, unsigned int adapterLength, unsigned int minOverlap, unsigned int minScore, struct deltas dlt)
/*
 *  Paired-end trimming
 *  Parameters:
 *
 */
{
    int res1, res2;
#pragma omp parallel sections
    {
#pragma omp section
        {
            res1 = trim_5_adapter_se(rec1, adapter, adapterLength, minOverlap, minScore, dlt);
        }

#pragma omp section
        {
            res2 = trim_5_adapter_se(rec2, adapter, adapterLength, minOverlap, minScore, dlt);
        }
    }
    return res1 + res2;
}


int trim_rev_5_adapter_pe(struct fqrec* rec1, struct fqrec* rec2, char* adapter, unsigned int adapterLength, char* adapter_rev, unsigned int adapter_revLength, unsigned int minOverlap, unsigned int minScore, struct deltas dlt, struct deltas dlt_rev)
/*
 *  Paired-end trimming
 *  Parameters:
 *
 */
{
    int res1, res2;
#pragma omp parallel sections
    {
#pragma omp section
        {
            res1 = trim_5_adapter_se(rec1, adapter, adapterLength, minOverlap, minScore, dlt);
        }

#pragma omp section
        {
            res2 = trim_5_adapter_se(rec2, adapter_rev, adapter_revLength, minOverlap, minScore, dlt_rev);
        }
    }
    return res1 + res2;
}

int trim_rev_3_adapter_pe(struct fqrec* rec1, struct fqrec* rec2, char* adapter, unsigned int adapterLength, char* adapter_rev, unsigned int adapter_revLength, unsigned int minOverlap, unsigned int minScore, struct deltas dlt, struct deltas dlt_rev)
/*
 *  Paired-end trimming
 *  Parameters:
 *      most self-explanatory
 *      minScore: minimum alignment score in Smith-Waterman matrix for adapter to be considered found.
 */
{
    int res1, res2;
#pragma omp parallel sections
    {
#pragma omp section
        {
            res1 = trim_3_adapter_se(rec1, adapter, adapterLength, minOverlap, minScore, dlt);
        }

#pragma omp section
        {
            res2 = trim_3_adapter_se(rec2, adapter_rev, adapter_revLength, minOverlap, minScore, dlt_rev);
        }
    }
    return res1 + res2;
}