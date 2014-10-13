#ifndef GAUSS_CC_
#define GAUSS_CC_

#include <cmath>
#include <complex>
#include "../Rand/random.cc"

using namespace std;


// gaussian random number generator
REAL real_gauss(void)
   {
   REAL phi, temp, radius, ris;
   
   phi=pi2*casuale();
   temp=-1.0*log(casuale());
   radius=sqrt(temp);

   ris=radius*cos(phi);

   return ris;
   }


// gaussian random number generator
void real_gauss(REAL &r1, REAL &r2)
   {
   REAL phi, temp, radius;
   
   phi=pi2*casuale();
   temp=-1.0*log(casuale());
   radius=sqrt(temp);

   r1=radius*cos(phi);
   r2=radius*sin(phi);
   }


// complex random gaussian number
complex<REAL> complex_gauss(void)
   {
   REAL re, im;

   real_gauss(re, im);

   return complex<REAL>(re, im);
   }


// gaussian disrtributed matrix in Su3 algebra:  a_{i-1}*lambda_i , a_{i-1} with normal distribution
//
//            | 0  1  0 |            | 0 -i  0 |            | 1  0  0 |
//  lambda_1= | 1  0  0 |  lambda_2= | i  0  0 |  lambda_3= | 0 -1  0 |
//            | 0  0  0 |            | 0  0  0 |            | 0  0  0 |
//
//            | 0  0  1 |            | 0  0 -i |            | 0  0  0 |
//  lambda_4= | 0  0  0 |  lambda_5= | 0  0  0 |  lambda_6= | 0  0  1 |
//            | 1  0  0 |            | i  0  0 |            | 0  1  0 |
//
//            | 0  0  0 |            | 1  0  0 |
//  lambda_7= | 0  0 -i |  lambda_8= | 0  1  0 |*(1/sqrt(3))
//            | 0  i  0 |            | 0  0 -2 |
//




#endif
