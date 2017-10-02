//
// Created by john on 6/30/2017.
//

#include "io.h"

int writeRecord(struct fqrec* rec, FILE* file)
    /*
     * writes rec to file
     */
{
    if(rec->seqLength == 0)
    {
        return 0;
    }

    //write name
    size_t written = fwrite(rec->name, sizeof(char), rec->nameLength - 1, file);
    if(written != rec->nameLength - 1)
    {
        return 1;
    }

    if(fputc('\n', file) != '\n')
    {
        return 1;
    }

    //write seq
    written = fwrite(rec->seq, sizeof(char), rec->seqLength, file);
    if(written != rec->seqLength * sizeof(char))
    {
        return 1;
    }

    //write separator
    if(fputc('\n', file) != '\n')
    {
        return 1;
    }

    if(fputc('+', file) != '+')
    {
        return 1;
    }

    if(fputc('\n', file) != '\n')
    {
        return 1;
    }

    //write qual
    written = fwrite(rec->qual, sizeof(char), rec->seqLength, file);
    if(written != rec->seqLength * sizeof(char))
    {
        return 1;
    }

    if(fputc('\n', file) != '\n')
    {
        return 1;
    }

    return 0;
}


int writePairedRecords(struct fqrec* rec1, FILE* file1, struct fqrec* rec2, FILE* file2)
{
    /*
     * writes rec1 and rec2 to file1 and file2, respectively
     */
    int m_res1 = 0, m_res2 = 0;
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            int res1 = writeRecord(rec1, file1);
            if(res1 > 0)
            {
                m_res1 = res1;
            }
        }

        #pragma omp section
        {
            int res2 = writeRecord(rec2, file2);
            if(res2 > 0)
            {
                m_res2 = res2;
            }
        }
    }
    return m_res1 + m_res2;
}

int init_file_data(unsigned long size, const char* filename, struct fqfiledata* data)
{
    data->buf_size = size;
    data->current_offset = 0;
    data->buf = aligned_malloc(data->buf_size * sizeof(char), sizeof(char));
    data->file = fopen(filename, "rb");
    data->current_size = 0;
    data->eof = 0;
    if(!data->buf)
    {
        //could not allocate memory
        return 1;
    }
    if(!data->file)
    {
        //could not open file
        return 2;
    }
    return 0;
}

int getNextRecord(struct fqfiledata* data, struct fqrec* rec)
    /*
     * gets the next record from the file and assigns it to rec
     */
{
    //check if we have enough data to read a record. If not, copy the remaining data to the start of buf and get more.

    if(!data->eof && data->current_size - data->current_offset < data->buf_size / RECORD_LENGTH_DIVIDER)
    {
        //move the data to the start of the array
        memcpy(data->buf, data->buf + data->current_offset, data->current_size - data->current_offset);
        size_t read_bytes = fread(data->buf + (data->current_size - data->current_offset), sizeof(char),
                                  data->buf_size - (data->current_size - data->current_offset), data->file);
        if(read_bytes == 0)
        {
            data->eof = 1;
        }
        data->current_size = (data->current_size - data->current_offset) + read_bytes;
        data->current_offset = 0;
    }


    //go through and mark offsets
    char line = 1; char brk = 0; char c = 0;
    unsigned long name_start_offset = data->current_offset; unsigned long name_end_offset = 0;
    unsigned long seq_start_offset = 0; unsigned long seq_end_offset = 0;
    unsigned long qual_start_offset = 0; unsigned long qual_end_offset = 0;
    for(unsigned long i = data->current_offset; i < data->current_size; i++)
    {
        c = data->buf[i];
        if(c == '\n')
        {
            switch(line)
            {
                case 1:
                    name_end_offset = i;
                    seq_start_offset = i+1;
                    line++;
                    break;
                case 2:
                    seq_end_offset = i;
                    line++;
                    break;
                case 3:
                    qual_start_offset = i+1;
                    line++;
                    break;
                case 4:
                    qual_end_offset = i;
                    brk = 1;
                    break;
                default:
                    break;
            }
        }
        else if(data->eof && i == data->current_size - 1)
        {
            qual_end_offset = i+1;
        }

        if(brk)
        {
            break;
        }

    }
    if(line < 4)
    {
        allocatefqrec(data->buf, 0,0,0,0,0,0, rec);
        return 2;
    }
    data->current_offset = qual_end_offset + 1;

    //printf("%lu %lu %lu %lu %lu %lu\n", name_start_offset, name_end_offset, seq_start_offset, seq_end_offset, qual_start_offset,
    //     qual_end_offset);

    return reallocatefqrec(data->buf, name_start_offset, name_end_offset, seq_start_offset, seq_end_offset, qual_start_offset, qual_end_offset, rec);
}