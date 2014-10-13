
#ifndef RANDOM_CC_
#define RANDOM_CC_

#include <cmath>
#include <cstdlib>


#include"./RANDOM/dSFMT.c"

dsfmt_t dsfmt;

// random number generator in (0,1)
REAL casuale(void)
   {   
    return dsfmt_genrand_open_open(&dsfmt);
   }


// random number initialization
void initrand(unsigned long s)
  {
   if(s==0)
    {
    dsfmt_init_gen_rand(&dsfmt, time(NULL));
    }
  else
    {
    dsfmt_init_gen_rand(&dsfmt, s);
    }
  }


// 4 parameters for random SU(2) matrix
void su2_rand(REAL &p0, REAL &p1, REAL &p2, REAL &p3)
  { 
  REAL p=2.0;
  while(p>1.0)
       {
       p0=1.0-2.0*casuale();
       p1=1.0-2.0*casuale();
       p2=1.0-2.0*casuale();
       p3=1.0-2.0*casuale();
       p=sqrt(p0*p0+p1*p1+p2*p2+p3*p3);
       }

  p0/=p;
  p1/=p;
  p2/=p;
  p3/=p;
  }


#endif
