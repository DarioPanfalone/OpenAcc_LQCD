#ifndef RANDOM_H_
#define RANDOM_H_

extern int rng_fake_gl_index; 
double casuale(void);
void initrand(unsigned long s);
void initrand_fromfile(const char * filename, unsigned long seed_default);
void saverand_tofile(const char * filename);
void su2_rand(double *pp);

#endif
