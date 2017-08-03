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
                         unsigned int in_a_row, unsigned int phred, char* input_file, char* output_base_name,
                         char* adapter, unsigned int adapter_length, unsigned int min_overlap, unsigned int min_score,
                         int trim_adapters);

void paired_end_pipeline(unsigned long buf_size, unsigned int qual_cutoff, unsigned int length_cutoff,
                         unsigned int in_a_row, unsigned int phred, char* forward_file, char* reverse_file,
                         char* output_base_name, char* adapter, unsigned int adapter_length, unsigned int min_overlap,
                         unsigned int min_score, int trim_adapters);

#endif //QUICKTRIM_PIPELINES_H
