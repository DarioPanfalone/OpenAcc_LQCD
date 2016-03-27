#ifndef RANDOM_C_
#define RANDOM_C_

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include"./RANDOM/dSFMT.c"
#include"./random.h"

dsfmt_t dsfmt;

// random number generator in (0,1)
double casuale(void)
   {   
     return dsfmt_genrand_open_open(&dsfmt);
     //return (double)rand()/((double)(RAND_MAX));
   }


// random number initialization
void initrand(unsigned long s)
  {
   if(s==0)
    {
      time_t t;
      //      srand((unsigned) time(&t));
    
      dsfmt_init_gen_rand(&dsfmt, time(NULL));
    }
  else
    {
      //      srand(s);
      dsfmt_init_gen_rand(&dsfmt, s);
    }
  }



void initrand_fromfile(const char * filename, unsigned long seed_default){

#ifndef NORANDOM
    FILE * seedfile = fopen(filename,"r");
    int reads;
    int problems_in_read = 0;
    if(!seedfile){
        problems_in_read = 1;
    }else{
        reads = fread(&dsfmt,1,sizeof(dsfmt_t), seedfile);
        if(reads!=sizeof(dsfmt_t))
            problems_in_read = 1;
        fclose(seedfile);
    }
    if(problems_in_read)
        initrand(seed_default);
    if(verbosity_lv > 0){
        if(problems_in_read){
            printf("Problems in reading file %s, ",filename );
            printf("Initialising RNG with seed %lu\n",seed_default);
        }
        else{
            printf("RNG initialised from file %s.\n",filename );

        }
    }
#else
    initrand(seed_default);
#endif


}

void saverand_tofile(const char * filename){

    FILE * seedfile = fopen(filename,"w");

    fwrite(&dsfmt, 1,sizeof(dsfmt_t), seedfile);
    fclose(seedfile);
    if(verbosity_lv > 0)
            printf("RNG status saved in file %s.\n",filename );

}






// 4 parameters for random SU(2) matrix
//void su2_rand(double &p0, double &p1, double &p2, double &p3)
void su2_rand(double *pp)
  { 
  double p=2.0;
  while(p>1.0)
       {
       pp[0]=1.0-2.0*casuale();
       pp[1]=1.0-2.0*casuale();
       pp[2]=1.0-2.0*casuale();
       pp[3]=1.0-2.0*casuale();
       p=sqrt(pp[0]*pp[0]+pp[1]*pp[1]+pp[2]*pp[2]+pp[3]*pp[3]);
       }

  pp[0]/=p;
  pp[1]/=p;
  pp[2]/=p;
  pp[3]/=p;
  }





#endif
