//
// Created by john on 6/29/2017.
//

//TODO: Possibly implement average quality over window?

#include "fastqrec.h"

int trim_se(struct fqrec* rec, unsigned int qual_cutoff, unsigned int length_cutoff, unsigned int in_a_row, unsigned int phred)
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
    unsigned int hq_bases = 0;
    for(unsigned int i = rec->seqLength - 1; i > 0; i--)
    {
        if(rec->seqLength < length_cutoff)
        {
            rec->seqLength = 0;
            return 1;
        }
        if(rec->seqAndQualBuf[i + rec->offset] - phred > qual_cutoff)
        {
            if(hq_bases >= in_a_row) {
                rec->seqLength += in_a_row;
                return 0;
            }
            hq_bases++;
        }
        else if(hq_bases > 0)
        {
            hq_bases = 0;
        }
        rec->seqLength--;
    }
    return 2;
}

int trim_pe(struct fqrec* rec1, struct fqrec* rec2, unsigned int qual_cutoff, unsigned int length_cutoff, unsigned int in_a_row, unsigned int phred)
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
            unsigned int hq_bases_rec1 = 0;
            for (unsigned int i = rec1->seqLength - 1; i > 0; i--) {
                if (rec1->seqLength < length_cutoff) {
                    rec1->seqLength = 0;
                    break;
                }
                if (rec1->seqAndQualBuf[i + rec1->offset] - phred > qual_cutoff) {
                    if (hq_bases_rec1 >= in_a_row) {
                        rec1->seqLength += in_a_row;
                        break;
                    }
                    hq_bases_rec1++;
                } else if (hq_bases_rec1 > 0) {
                    hq_bases_rec1 = 0;
                }
                rec1->seqLength--;
            }
        }

        #pragma omp section
        {
            unsigned int hq_bases_rec2 = 0;
            for (unsigned int i = rec2->seqLength - 1; i > 0; i--) {
                if (rec2->seqLength < length_cutoff) {
                    rec2->seqLength = 0;
                    break;
                }
                if (rec2->seqAndQualBuf[i + rec2->offset] - phred > qual_cutoff) {
                    if (hq_bases_rec2 >= in_a_row) {
                        rec2->seqLength += in_a_row;
                        break;
                    }
                    hq_bases_rec2++;
                } else if (hq_bases_rec2 > 0) {
                    hq_bases_rec2 = 0;
                }
                rec2->seqLength--;
            }
        }
    }
    if (rec1->seqLength == 0 || rec2->seqLength == 0) {
        return 1;
    }
    return 0;
}