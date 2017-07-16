//
// Created by john on 7/8/2017.
//

#ifndef QUICKTRIM_MEM_H
#define QUICKTRIM_MEM_H

#include <stdlib.h>

void* aligned_malloc(size_t size, size_t align);
void aligned_free(void *ptr);

#endif //QUICKTRIM_MEM_H
