#ifndef _MEMORY_H
#define _MEMORY_H

#include <stdlib.h>

int posix_memalign_wrapper(void **memptr,size_t alignment,size_t size,const char* varname);
void free_wrapper(void *memptr);

#define POSIX_MEMALIGN_WRAPPER(MEMPTR,ALIGNMENT,SIZE)		\
  posix_memalign_wrapper(MEMPTR,ALIGNMENT,SIZE,#MEMPTR)

extern size_t memory_used;
extern size_t max_memory_used;

struct memory_allocated_t
{
  void* ptr;
  const char* varname;
  size_t size;
  struct memory_allocated_t *next;
  
};

extern struct memory_allocated_t *memory_allocated_base;

#endif
