
#ifndef RANDOM_CC_
#define RANDOM_CC_

// random number generator in (0,1)
extern "C" {
  double casuale(void);
}

// random number initialization
extern "C" {
  void initrand(unsigned long s);
}

// 4 parameters for random SU(2) matrix
extern "C" {
  //  void su2_rand(double &p0, double &p1, double &p2, double &p3);
  void su2_rand(double *pp);
}

#endif
