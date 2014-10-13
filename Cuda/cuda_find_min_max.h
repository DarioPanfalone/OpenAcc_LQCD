#ifndef CUDA_FIND_MIN_MAX_H
#define CUDA_FIND_MIN_MAX_H

extern "C" void cuda_find_max(REAL *max);
extern "C" void cuda_find_min(REAL *min, const REAL max);

#endif
