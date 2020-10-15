#ifndef _MEMORY_H
#define _MEMORY_H

#include <stdlib.h>

int posix_memalign_wrapper(void **memptr,size_t alignment,size_t size);
void free_wrapper(void *memptr);

extern size_t memory_used;
extern size_t max_memory_used;

#endif
