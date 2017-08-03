//
// Created by john on 8/2/2017.
//

#include "pipelines.h"

void single_end_pipeline(unsigned long buf_size, unsigned int qual_cutoff, unsigned int length_cutoff,
                         unsigned int in_a_row, unsigned int phred, char* input_file, char* output_base_name,
                         char* adapter, unsigned int adapter_length, unsigned int min_overlap, unsigned int min_score,
                         int trim_adapters)
{
    struct fqfiledata data;
    int res = init_file_data(buf_size, input_file, &data);
    if(res == 1) {
        printf("Error allocating memory.\n");
        return;
    } else if (res == 2) {
        printf("Error opening file.\n");
    }
    FILE* out = fopen(strcat(output_base_name, ".fastq"), "wb");

    struct fqrec rec;
    rec.nameLength = 0;
    rec.seqLength = 0;
    rec.offset = 0;

    int err = 0;

    while(1) {
        err = getNextRecord(&data, &rec);
        if(err) {
            printf("Error getting record: %d.\n", err);
            break;
        }
        if(trim_adapters)
        {
            trim_adapter_se(&rec, adapter, adapter_length, min_overlap, min_score);
        }
        err = trim_se(&rec, qual_cutoff, length_cutoff, in_a_row, phred);
        if(err) {
            //printf("Error trimming: %d.\n", err);
            continue;
        }
        err = writeRecord(&rec, out);
        if(err) {
            printf("Error writing: %d.\n", err);
            break;
        }
    }
    printf("Done reading file or error reading.\n");
    fclose(out);
    fclose(data.file);
    freefqrec(&rec);
    printf("Done.");
}

void paired_end_pipeline(unsigned long buf_size, unsigned int qual_cutoff, unsigned int length_cutoff,
                         unsigned int in_a_row, unsigned int phred, char* forward_file, char* reverse_file,
                         char* output_base_name, char* adapter, unsigned int adapter_length, unsigned int min_overlap,
                         unsigned int min_score, int trim_adapters)
{
    char* output_base_name_copy = (char*)calloc(strlen(output_base_name) + 50, sizeof(char));
    strncpy(output_base_name_copy, output_base_name, strlen(output_base_name));
    output_base_name[strlen(output_base_name)] = '\0';
    struct fqfiledata forward_data;
    int res = init_file_data(buf_size, forward_file, &forward_data);
    if(res == 1) {
        printf("Error allocating memory.\n");
        return;
    } else if (res == 2) {
        printf("Error opening file.\n");
    }
    FILE* forward_out = fopen(strcat(output_base_name, ".forward.fastq"), "wb");

    struct fqfiledata reverse_data;
    res = init_file_data(buf_size, reverse_file, &reverse_data);
    if(res == 1) {
        printf("Error allocating memory.\n");
        return;
    } else if (res == 2) {
        printf("Error opening file.\n");
    }
    FILE* reverse_out = fopen(strcat(output_base_name_copy, ".reverse.fastq"), "wb");
    struct fqrec rec1;
    rec1.nameLength = 0;
    rec1.seqLength = 0;
    rec1.offset = 0;

    struct fqrec rec2;
    rec2.nameLength = 0;
    rec2.seqLength = 0;
    rec2.offset = 0;

    int err = 0;
    int readerr = 0;

    allocatefqrec(forward_data.buf, 0,0,0,0,0,0, &rec1);
    allocatefqrec(reverse_data.buf, 0,0,0,0,0,0, &rec2);


    while(1) {
        err = getNextRecord(&forward_data, &rec1);
        if(err) {
            printf("Error getting fwd record: %d.\n", err);
            readerr = 1;
        }
        err = getNextRecord(&reverse_data, &rec2);
        if(err) {
            printf("Error getting rev record: %d.\n", err);
            readerr = 1;
        }
        if(readerr)
        {
            break;
        }
        if(trim_adapters)
        {
            trim_adapter_pe(&rec1, &rec2, adapter, adapter_length, min_overlap, min_score);
        }
        err = trim_pe(&rec1, &rec2, qual_cutoff, length_cutoff, in_a_row, phred);
        if(err) {
            //printf("Error trimming: %d.\n", err);
            continue;
        }
        err = writePairedRecords(&rec1, forward_out, &rec2, reverse_out);
        if(err) {
            printf("Error writing: %d.\n", err);
            break;
        }
    }
    printf("Done reading files or error reading.\n");
    fclose(forward_out);
    fclose(reverse_out);
    fclose(forward_data.file);
    fclose(reverse_data.file);
    freefqrec(&rec1);
    freefqrec(&rec2);
    printf("Done.");
}