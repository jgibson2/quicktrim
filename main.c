#include <stdio.h>
#include "io.h"
#include "trim.h"
int main() {
    struct fqfiledata data;
    int res = init_file_data(100000, "C:\\Users\\john\\CLionProjects\\quicktrim\\test.fq", &data);
    if(res == 1) {
        printf("Error allocating memory.\n");
        return 1;
    } else if (res == 2) {
        printf("Error opening file.\n");
    }
    FILE* out = fopen("C:\\Users\\john\\CLionProjects\\quicktrim\\test.trimmed.fq", "wb");

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
        err = trim_se(&rec, 30, 50, 5, 33);
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
    return 0;
}