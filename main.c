#include <stdio.h>
#include "io.h"
#include "trim.h"
#include <unistd.h>
#include <ctype.h>

void single_end_pipeline(unsigned long buf_size, unsigned int qual_cutoff, unsigned int length_cutoff,
                         unsigned int in_a_row, unsigned int phred, char* input_file, char* output_base_name);
void paired_end_pipeline(unsigned long buf_size, unsigned int qual_cutoff, unsigned int length_cutoff,
                         unsigned int in_a_row, unsigned int phred, char* forward_file, char* reverse_file,
                         char* output_base_name);

int main(int argc, char** argv) {
    unsigned long buf_size = 100000;
    unsigned int qual_cutoff = 30, length_cutoff = 20, in_a_row = 5, phred = 33;
    char* input_file;
    char* forward_file;
    char* reverse_file;
    char* output_base_name;
    int paired_end = 0;

    char* optstring = "q:l:r:p:i:o:1:2:";
    int c;
    opterr = 0;

    while ((c = getopt (argc, argv, optstring)) != -1) {
        char *end;
        switch (c) {
            case 'b':
                buf_size = strtoul(optarg, &end, 10);
                break;
            case 'q':
                qual_cutoff = strtoul(optarg, &end, 10);
                break;
            case 'l':
                length_cutoff = strtoul(optarg, &end, 10);
                break;
            case 'r':
                in_a_row = strtoul(optarg, &end, 10);
                break;
            case 'p':
                phred = strtoul(optarg, &end, 10);
                break;
            case 'i':
                input_file = (char*)calloc(strlen(optarg) + 10, sizeof(char));
                strncpy(input_file, optarg, strlen(optarg));
                input_file[strlen(optarg)] = '\0';
                break;
            case 'o':
                output_base_name = (char*)calloc(strlen(optarg) + 50, sizeof(char));
                strncpy(output_base_name, optarg, strlen(optarg));
                output_base_name[strlen(optarg)] = '\0';
                break;
            case '1':
                forward_file = (char*)calloc(strlen(optarg) + 10, sizeof(char));
                strncpy(forward_file, optarg, strlen(optarg));
                forward_file[strlen(optarg)] = '\0';
                paired_end = 1;
                break;
            case '2':
                reverse_file = (char*)calloc(strlen(optarg) + 10, sizeof(char));
                strncpy(reverse_file, optarg, strlen(optarg));
                reverse_file[strlen(optarg)] = '\0';
                paired_end = 1;
                break;
            case '?':
                if (strchr(optstring, optopt) != NULL)
                    fprintf(stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint(optopt))
                    fprintf(stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf(stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                return 1;
            default:
                abort();
        }
    }

    if(paired_end) {
        paired_end_pipeline(buf_size,  qual_cutoff,  length_cutoff, in_a_row,  phred,  forward_file,  reverse_file,
                                  output_base_name);
    }
    else
    {
        single_end_pipeline(buf_size,  qual_cutoff,  length_cutoff, in_a_row,  phred,  input_file,
                            output_base_name);
    }


    return 0;

}

void single_end_pipeline(unsigned long buf_size, unsigned int qual_cutoff, unsigned int length_cutoff,
                        unsigned int in_a_row, unsigned int phred, char* input_file, char* output_base_name)
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
        err = trim_se(&rec, qual_cutoff, length_cutoff, in_a_row, phred);
        if(err) {
            printf("Error trimming: %d.\n", err);
            break;
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
                         char* output_base_name)
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
        err = trim_pe(&rec1, &rec2, qual_cutoff, length_cutoff, in_a_row, phred);
        if(err) {
            printf("Error trimming: %d.\n", err);
            break;
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