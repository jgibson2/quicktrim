//
// Created by john on 7/9/2017.
//

#include "mem.h"

void* aligned_malloc(size_t size, size_t align) {
    void *result;
#ifdef __INTEL_COMPILER
    result = _mm_alloc(size, align);
#elif _MSC_VER
    result = _aligned_malloc(size, align);
#elif __STDC_VERSION__ >= 201112L
    result = aligned_alloc(align, size);
#elif _POSIX_VERSION >= 200112L
    if(posix_memalign(&result, align, size))
    {
        result = 0;
    }
#else
    result = malloc(size);
#endif
    return result;
}

void aligned_free(void *ptr) {
#ifdef __INTEL_COMPILER
    _mm_free(ptr);
#elif _MSC_VER
    _aligned_free(ptr);
#else
    free(ptr);
#endif

}
