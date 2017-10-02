//
// Created by john on 8/2/2017.
//

#include "pipelines.h"

void single_end_pipeline(unsigned long buf_size, unsigned char qual_cutoff, unsigned int length_cutoff,
                         unsigned int in_a_row, unsigned char phred, unsigned char method, char* input_file, char* output_base_name,
                         char* adapter_3, unsigned int adapter_3_length, unsigned int min_3_overlap, unsigned int min_3_score,
                         int trim_3_adapters, char* adapter_5, unsigned int adapter_5_length, unsigned int min_5_overlap, unsigned int min_5_score,
                         int trim_5_adapters)
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
    allocatefqrec(data.buf, 0,0,0,0,0,0, &rec);


    int err = 0;

    struct deltas dlt_3; //from align.h
    dlt_3.delta2 = aligned_malloc(adapter_3_length * sizeof(int), sizeof(int));
    struct deltas dlt_5;
    dlt_5.delta2 = aligned_malloc(adapter_5_length * sizeof(int), sizeof(int));

    if(trim_3_adapters) {make_deltas(&dlt_3, adapter_3, adapter_3_length);}
    if(trim_5_adapters) {make_deltas(&dlt_5, adapter_5, adapter_5_length);}

    printf("Starting run at ");
    time_t now = time(NULL);
    printf(ctime(&now));
    printf("\n");

    unsigned long int num_records = 0;

    while(1) {
        err = getNextRecord(&data, &rec);
        if(err) {
            printf("Error getting record: %d.\n", err);
            break;
        }
        if(trim_3_adapters)
        {
            trim_3_adapter_se(&rec, adapter_3, adapter_3_length, min_3_overlap, min_3_score, dlt_3);
        }
        if(trim_5_adapters)
        {
            trim_5_adapter_se(&rec, adapter_5, adapter_5_length, min_5_overlap, min_5_score, dlt_5);
        }
        err = trim_se(&rec, qual_cutoff, length_cutoff, in_a_row, phred, method);
        if(err) {
            //printf("Error trimming: %d.\n", err);
            continue;
        }
        err = writeRecord(&rec, out);
        if(err) {
            printf("Error writing: %d.\n", err);
            break;
        }
        num_records += 1;
    }
    printf("Done reading file or error reading.\n");
    fclose(out);
    fclose(data.file);
    printf("Done.");

    printf("Run ended at ");
    now = time(NULL);
    printf(ctime(&now));
    printf("\n");

    printf("Read %lu pairs of reads.", num_records);
}

void paired_end_pipeline(unsigned long buf_size, unsigned char qual_cutoff, unsigned int length_cutoff,
                         unsigned int in_a_row, unsigned char phred, unsigned char method, char* forward_file, char* reverse_file,
                         char* output_base_name, char* adapter_3, unsigned int adapter_3_length, unsigned int min_3_overlap, unsigned int min_3_score,
                         int trim_3_adapters, char* adapter_5, unsigned int adapter_5_length, unsigned int min_5_overlap, unsigned int min_5_score,
                         int trim_5_adapters, char* adapter_rev_5, unsigned int adapter_rev_5_length, int use_5_rev_adapters, char* adapter_rev_3,
                         unsigned int adapter_rev_3_length, int use_3_rev_adapters)
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

    struct fqrec rec2;
    rec2.nameLength = 0;
    rec2.seqLength = 0;

    int err = 0;
    int readerr = 0;

    struct deltas dlt_3; //from align.h
    dlt_3.delta2 = aligned_malloc(adapter_3_length * sizeof(int), sizeof(int));
    struct deltas dlt_5;
    dlt_5.delta2 = aligned_malloc(adapter_5_length * sizeof(int), sizeof(int));

    if(trim_3_adapters) {make_deltas(&dlt_3, adapter_3, adapter_3_length);}
    if(trim_5_adapters) {make_deltas(&dlt_5, adapter_5, adapter_5_length);}

    struct deltas dlt_rev_3; //from align.h
    dlt_rev_3.delta2 = aligned_malloc(adapter_rev_3_length * sizeof(int), sizeof(int));
    struct deltas dlt_rev_5;
    dlt_rev_5.delta2 = aligned_malloc(adapter_rev_5_length * sizeof(int), sizeof(int));

    if(trim_3_adapters && use_3_rev_adapters) {make_deltas(&dlt_rev_3, adapter_rev_3, adapter_rev_3_length);}
    if(trim_5_adapters && use_5_rev_adapters) {make_deltas(&dlt_rev_5, adapter_rev_5, adapter_rev_5_length);}


    allocatefqrec(forward_data.buf, 0,0,0,0,0,0, &rec1);
    allocatefqrec(reverse_data.buf, 0,0,0,0,0,0, &rec2);

    printf("Starting run at ");
    time_t now = time(NULL);
    printf(ctime(&now));
    printf("\n");

    unsigned long int num_records = 0;

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
        if(trim_3_adapters)
        {
            if(use_3_rev_adapters) {
                trim_rev_3_adapter_pe(&rec1, &rec2, adapter_3, adapter_3_length, adapter_rev_3, adapter_rev_3_length, min_3_overlap, min_3_score, dlt_3, dlt_rev_3);
            } else {
                trim_3_adapter_pe(&rec1, &rec2, adapter_3, adapter_3_length, min_3_overlap, min_3_score, dlt_3);
            }
        }
        if(trim_5_adapters)
        {
            if(use_5_rev_adapters) {
                trim_rev_5_adapter_pe(&rec1, &rec2, adapter_5, adapter_5_length, adapter_rev_5, adapter_rev_5_length, min_5_overlap, min_5_score, dlt_5, dlt_rev_5);
            } else {
                trim_5_adapter_pe(&rec1, &rec2, adapter_5, adapter_5_length, min_5_overlap, min_5_score, dlt_5);
            }
        }
        err = trim_pe(&rec1, &rec2, qual_cutoff, length_cutoff, in_a_row, phred, method);
        if(err) {
            //printf("Error trimming: %d.\n", err);
            continue;
        }
        err = writePairedRecords(&rec1, forward_out, &rec2, reverse_out);
        if(err) {
            printf("Error writing: %d.\n", err);
            break;
        }
        num_records += 1;
    }
    printf("Done reading files or error reading.\n");
    fclose(forward_out);
    fclose(reverse_out);
    fclose(forward_data.file);
    aligned_free(forward_data.buf);
    fclose(reverse_data.file);
    aligned_free(reverse_data.buf);
    printf("Done.\n");

    printf("Run ended at ");
    now = time(NULL);
    printf(ctime(&now));
    printf("\n");

    printf("Read %lu pairs of reads.", num_records);
}