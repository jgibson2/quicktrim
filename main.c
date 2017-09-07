#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include "pipelines.h"

int main(int argc, char** argv) {
    unsigned long buf_size = 1000000;
    unsigned char qual_cutoff = 30, phred = 33;
    unsigned int length_cutoff = 20, in_a_row = 5;
    char* input_file = "";
    char* forward_file = "";
    char* reverse_file = "";
    char* output_base_name = "out";
    int paired_end = 0;
    char* adapter_3 = "";
    char* adapter_5 = "";
    char* adapter_rev_3 = "";
    char* adapter_rev_5 = "";
    unsigned int adapter_3_length = 0, adapter_5_length = 0, adapter_rev_3_length = 0, adapter_rev_5_length = 0;
    unsigned int min_3_overlap = 10, min_5_overlap = 10;
    unsigned int min_3_score = 8, min_5_score = 8;
    int trim_3_adapters = 0, trim_5_adapters = 0;
    int use_3_rev_adapters = 0, use_5_rev_adapters = 0;
    unsigned int method = 0;

    char* optstring = "q:l:r:p:z:i:o:1:2:a:b:v:s:A:B:V:S:mh";
    int c;
    opterr = 0;

    while ((c = getopt (argc, argv, optstring)) != -1) {
        char *end;
        switch (c) {
            case 'h':
                printf("quicktrim: a short-read trimming application written in C.\n");
                printf("Options:\n");
                printf("\to: Output base filename (default out)\n");
                printf("\ti: Single-ended input file\n");
                printf("\t1: Paired-end input file 1\n");
                printf("\t2: Paired-end input file 2\n");
                printf("\ta: 3' fwd adapter\n");
                printf("\tA: 5' fwd adapter\n");
                printf("\ta: 3' rev adapter\n");
                printf("\tA: 5' rev adapter\n");
                printf("\tq: Set quality cutoff (default 30)\n");
                printf("\tl: Set length cutoff (default 20)\n");
                printf("\tr: Set number of high-quality bases in a row for use with fast trimming (default 5)\n");
                printf("\tp: Phred score base (default 33)\n");
                printf("\tb: Change buffer size in bytes (default 1000000) {NOTE: must be larger than record length / 10}\n");
                printf("\tv: Minimum 3' adapter overlap (default 10)\n");
                printf("\ts: Minimum 3' alignment score (default 8)\n");
                printf("\tV: Minimum 5' adapter overlap (default 10)\n");
                printf("\tS: Minimum 5' alignment score (default 8)\n");
                printf("\tm: Switch to faster, but less accurate mode of triming (HQ bases in a row)\n");
                printf("\th: Print this menu");
                return 1;
            case 'm':
                method = 1;
            case 'z':
                buf_size = strtoul(optarg, &end, 10);
                break;
            case 'q':
                qual_cutoff = (char)strtoul(optarg, &end, 10);
                break;
            case 'l':
                length_cutoff = strtoul(optarg, &end, 10);
                break;
            case 'r':
                in_a_row = strtoul(optarg, &end, 10);
                break;
            case 'p':
                phred = (char)strtoul(optarg, &end, 10);
                break;
            case 'v':
                min_3_overlap = strtoul(optarg, &end, 10);
                break;
            case 's':
                min_3_score = strtoul(optarg, &end, 10);
                break;
            case 'V':
                min_5_overlap = strtoul(optarg, &end, 10);
                break;
            case 'S':
                min_5_score = strtoul(optarg, &end, 10);
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
            case 'a':
                adapter_3 = (char*)calloc(strlen(optarg) + 10, sizeof(char));
                strncpy(adapter_3, optarg, strlen(optarg));
                adapter_3[strlen(optarg)] = '\0';
                trim_3_adapters = 1;
                adapter_3_length = strlen(optarg);
                break;
            case 'A':
                adapter_5 = (char*)calloc(strlen(optarg) + 10, sizeof(char));
                strncpy(adapter_5, optarg, strlen(optarg));
                adapter_5[strlen(optarg)] = '\0';
                trim_5_adapters = 1;
                adapter_5_length = strlen(optarg);
                break;
            case 'b':
                adapter_rev_3 = (char*)calloc(strlen(optarg) + 10, sizeof(char));
                strncpy(adapter_rev_3, optarg, strlen(optarg));
                adapter_rev_3[strlen(optarg)] = '\0';
                trim_3_adapters = 1;
                use_3_rev_adapters = 1;
                adapter_rev_3_length = strlen(optarg);
                break;
            case 'B':
                adapter_rev_5 = (char*)calloc(strlen(optarg) + 10, sizeof(char));
                strncpy(adapter_rev_5, optarg, strlen(optarg));
                adapter_rev_5[strlen(optarg)] = '\0';
                trim_5_adapters = 1;
                use_5_rev_adapters = 1;
                adapter_rev_5_length = strlen(optarg);
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
        paired_end_pipeline(buf_size,  qual_cutoff,  length_cutoff, in_a_row, phred, method, forward_file,  reverse_file,
                                  output_base_name, adapter_3, adapter_3_length, min_3_overlap, min_3_score, 
                            trim_3_adapters, adapter_5, adapter_5_length, min_5_overlap, min_5_score, trim_5_adapters, adapter_rev_5,
                            adapter_rev_5_length, use_5_rev_adapters, adapter_rev_3, adapter_rev_3_length, use_3_rev_adapters);
    }
    else
    {
        single_end_pipeline(buf_size,  qual_cutoff,  length_cutoff, in_a_row,  phred, method, input_file,
                            output_base_name, adapter_3, adapter_3_length, min_3_overlap, min_3_score,
                            trim_3_adapters, adapter_5, adapter_5_length, min_5_overlap, min_5_score, trim_5_adapters);
    }

    return 0;

}

