//
// Created by john on 8/2/2017.
//

#ifndef QUICKTRIM_PIPELINES_H
#define QUICKTRIM_PIPELINES_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "trim.h"
#include "io.h"


void single_end_pipeline(unsigned long buf_size, unsigned char qual_cutoff, unsigned int length_cutoff,
                         unsigned int in_a_row, unsigned char phred, unsigned char method, char* input_file, char* output_base_name,
                         char* adapter_3, unsigned int adapter_3_length, unsigned int min_3_overlap, unsigned int min_3_score,
                         int trim_3_adapters, char* adapter_5, unsigned int adapter_5_length, unsigned int min_5_overlap, unsigned int min_5_score,
                         int trim_5_adapters);

void paired_end_pipeline(unsigned long buf_size, unsigned char qual_cutoff, unsigned int length_cutoff,
                         unsigned int in_a_row, unsigned char phred, unsigned char method, char* forward_file, char* reverse_file,
                         char* output_base_name, char* adapter_3, unsigned int adapter_3_length, unsigned int min_3_overlap, unsigned int min_3_score,
                         int trim_3_adapters, char* adapter_5, unsigned int adapter_5_length, unsigned int min_5_overlap, unsigned int min_5_score,
                         int trim_5_adapters, char* adapter_rev_5, unsigned int adapter_rev_5_length, int use_5_rev_adapters, char* adapter_rev_3,
                         unsigned int adapter_rev_3_length, int use_3_rev_adapters);

#endif //QUICKTRIM_PIPELINES_H
