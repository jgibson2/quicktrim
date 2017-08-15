#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include "pipelines.h"

int main(int argc, char** argv) {
    unsigned long buf_size = 100000;
    unsigned int qual_cutoff = 30, length_cutoff = 20, in_a_row = 5, phred = 33;
    char* input_file;
    char* forward_file;
    char* reverse_file;
    char* output_base_name;
    int paired_end = 0;
    char* adapter_3;
    char* adapter_5;
    unsigned int adapter_3_length = 0, adapter_5_length = 0;
    unsigned int min_3_overlap = 10, min_5_overlap = 10;
    unsigned int min_3_score = 8, min_5_score = 8;
    int trim_3_adapters = 0, trim_5_adapters = 0;

    char* optstring = "q:l:r:p:i:o:1:2:a:v:s:A:V:S:";
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
                                  output_base_name, adapter_3, adapter_3_length, min_3_overlap, min_3_score, 
                            trim_3_adapters, adapter_5, adapter_5_length, min_5_overlap, min_5_score, trim_5_adapters);
    }
    else
    {
        single_end_pipeline(buf_size,  qual_cutoff,  length_cutoff, in_a_row,  phred,  input_file,
                            output_base_name, adapter_3, adapter_3_length, min_3_overlap, min_3_score,
                            trim_3_adapters, adapter_5, adapter_5_length, min_5_overlap, min_5_score, trim_5_adapters);
    }

    return 0;

}

