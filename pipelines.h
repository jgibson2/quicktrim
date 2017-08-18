//
// Created by john on 8/2/2017.
//

#ifndef QUICKTRIM_PIPELINES_H
#define QUICKTRIM_PIPELINES_H

#include <stdio.h>
#include <stdlib.h>
#include "trim.h"
#include "io.h"


void single_end_pipeline(unsigned long buf_size, unsigned int qual_cutoff, unsigned int length_cutoff,
                         unsigned int in_a_row, unsigned int phred, unsigned int method, char* input_file, char* output_base_name,
                         char* adapter_3, unsigned int adapter_3_length, unsigned int min_3_overlap, unsigned int min_3_score,
                         int trim_adapters_3, char* adapter_5, unsigned int adapter_5_length, unsigned int min_5_overlap, unsigned int min_5_score,
                         int trim_adapters_5);

void paired_end_pipeline(unsigned long buf_size, unsigned int qual_cutoff, unsigned int length_cutoff,
                         unsigned int in_a_row, unsigned int phred, unsigned int method, char* forward_file, char* reverse_file,
                         char* output_base_name, char* adapter_3, unsigned int adapter_3_length, unsigned int min_3_overlap, unsigned int min_3_score,
                         int trim_adapters_3, char* adapter_5, unsigned int adapter_5_length, unsigned int min_5_overlap, unsigned int min_5_score,
                         int trim_adapters_5);

#endif //QUICKTRIM_PIPELINES_H
