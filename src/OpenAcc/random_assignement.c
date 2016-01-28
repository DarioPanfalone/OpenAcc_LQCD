#ifndef RANDOM_ASSIGNEMENT_C_
#define RANDOM_ASSIGNEMENT_C_
// random number generator in (0,1)
double casuale(void);

/*
// random number initialization
extern "C" {
  void initrand(unsigned long s);
}

// 4 parameters for random SU(2) matrix
extern "C" {
  void su2_rand(double *pp);
}
*/

#include "./geometry.h"
#include "./su3_utilities.h"
#include "./alloc_vars.h"
#include "./random_assignement.h"
#include "../DbgTools/debug_macros_glvarcheck.h"
#include "./single_types.h"

#define acc_twopi 2*3.14159265358979323846

#ifdef __GNUC__
 #include "math.h"
 #ifndef M_PI
  #define M_PI 3.14159265358979323846
 #endif
#endif





// gaussian random number generator
double double_gauss(void)
{
  double phi, temp, radius, ris;

  phi=acc_twopi*casuale();
  temp=-1.0*log(casuale());
  radius=sqrt(temp);

  ris=radius*cos(phi);

  return ris;
}


// gaussian random number generator
void two_double_gauss(double *r)
{
  double phi, temp, radius;

  phi=acc_twopi*casuale();
  temp=-1.0*log(casuale());
  radius=sqrt(temp);

  r[0]=radius*cos(phi);
  r[1]=radius*sin(phi);
}


// complex random gaussian number
d_complex d_complex_gauss(void)
{
  double num[2];

  two_double_gauss(num);

  return num[0] + I*num[1];
}

void generate_vec3_soa_gauss(__restrict vec3_soa * const vect){
  SETINUSE(vect);
  int t;
  for(t=0; t<sizeh; t++) {
    vect->c0[t]=d_complex_gauss();
    vect->c1[t]=d_complex_gauss();
    vect->c2[t]=d_complex_gauss();
  }
}


inline double norm3(d_complex a,d_complex b, d_complex c){

    double norm ;
    norm  = creal(a)*creal(a) +  cimag(a)*cimag(a);
    norm += creal(b)*creal(b) +  cimag(b)*cimag(b);
    norm += creal(c)*creal(c) +  cimag(c)*cimag(c);
     return norm;

}

void generate_random_su3(single_su3 * m, double f){

//////////GENERATE RANDOM MATRIX ///////////////////////////
// only first two rows
   m->comp[0][0] =1+(1.0-2.0*casuale()+I*(1.0-2.0*casuale()))*f;
   m->comp[0][1] =  (1.0-2.0*casuale()+I*(1.0-2.0*casuale()))*f;
   m->comp[0][2] =  (1.0-2.0*casuale()+I*(1.0-2.0*casuale()))*f;
   m->comp[1][0] =  (1.0-2.0*casuale()+I*(1.0-2.0*casuale()))*f;
   m->comp[1][1] =1+(1.0-2.0*casuale()+I*(1.0-2.0*casuale()))*f;
   m->comp[1][2] =  (1.0-2.0*casuale()+I*(1.0-2.0*casuale()))*f;

////////// Unitarize it ///////////////////////////////////////////
   double norm = norm3(m->comp[0][0],m->comp[0][1],m->comp[0][2]);
   norm=1.0/sqrt(norm);
    
   m->comp[0][0] *= norm;
   m->comp[0][1] *= norm;
   m->comp[0][2] *= norm;
    
   d_complex prod  = conj(m->comp[0][0]) * m->comp[1][0] + conj(m->comp[0][1]) * m->comp[1][1] + conj(m->comp[0][2]) * m->comp[1][2];
    
    m->comp[1][0] -= prod * m->comp[0][0];
    m->comp[1][1] -= prod * m->comp[0][1];
    m->comp[1][2] -= prod * m->comp[0][2];
    
    norm = norm3(m->comp[1][0],m->comp[1][1],m->comp[1][2]);
    norm=1.0/sqrt(norm);
    


    m->comp[1][0] *= norm;
    m->comp[1][1] *= norm;
    m->comp[1][2] *= norm;
    

    rebuild3row(m);

}

// iterations of random assignements over all the 8 components
void generate_Momenta_gauss(__restrict thmat_soa * const mom8){
    SETINUSE(mom);
    int mu;
    for(mu=0; mu<8; mu++){
        thmat_soa * mom = &mom8[mu];

        int i,t;
        double casuali[8], aux[2];
        double uno_su_radice_di_tre = 1.0/sqrt(3.0);
        for(t=0; t<sizeh; t++) {

            for(i=0; i<4; i++)
            {
                two_double_gauss(aux);
                casuali[2*i]   = aux[0];
                casuali[2*i+1] = aux[1];
            }
            mom->rc00[t] =  casuali[2] + casuali[7] * uno_su_radice_di_tre;
            mom->rc11[t] = -casuali[2] + casuali[7] * uno_su_radice_di_tre;
            mom->c01[t]  =  casuali[0] - casuali[1] * I;
            mom->c02[t]  =  casuali[3] - casuali[4] * I;
            mom->c12[t]  =  casuali[5] - casuali[6] * I;
        }  // t  
    }
}

// iterations of random assignements over all the 8 components
void generate_Conf_cold(__restrict su3_soa * const conf,double factor){
  SETINUSE(conf);
  for(int mu=0; mu<8; mu++)
      for(int idx=0;idx<sizeh;idx++){
         single_su3 aux;
         generate_random_su3(&aux,factor);
         single_gl3_into_su3_soa(&conf[mu],idx,&aux);
  }
}


void generate_vec3_soa_z2noise(__restrict vec3_soa * const vect){
  SETINUSE(vect);
  int t;
  double p;
  for(t=0; t<sizeh; t++) {
    p = casuale();
    if(p<0.5){
      vect->c0[t]=  1.0 + 0.0*I;
    }else{
      vect->c0[t]= -1.0 + 0.0*I;
    }
    p = casuale();
    if(p<0.5){
      vect->c1[t]=  1.0 + 0.0*I;
    }else{
      vect->c1[t]= -1.0 + 0.0*I;
    }
    p = casuale();
    if(p<0.5){
      vect->c2[t]=  1.0 + 0.0*I;
    }else{
      vect->c2[t]= -1.0 + 0.0*I;
    }
  }
}



#endif


