#ifndef GLOBAL_MACRO_CC_
#define GLOBAL_MACRO_CC_



#define NITER 100
#define ALIGN 128

#define REAL double

// number of threads if GPU is used
#define NUM_THREADS 128

#define CONF_FILE config
#define ERROR_FILE data.err

// if #define val1 val2 then QUOTEME(val1) gives the string "val2"
#define _QUOTEME(x) #x
#define QUOTEME(x) _QUOTEME(x)

// to be used in cuda_err.cu
#define AT __FILE__ " : " QUOTEME(__LINE__)

#endif


