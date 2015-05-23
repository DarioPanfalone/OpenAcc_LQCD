#ifndef RANDOM_ASSIGNEMENT_C_
#define RANDOM_ASSIGNEMENT_C_
/*
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
  void su2_rand(double *pp);
}
*/


// gaussian random number generator
double double_gauss(void)
{
  double phi, temp, radius, ris;

  phi=acc_pi2*casuale();
  temp=-1.0*log(casuale());
  radius=sqrt(temp);

  ris=radius*cos(phi);

  return ris;
}


// gaussian random number generator
void two_double_gauss(double *r)
{
  double phi, temp, radius;

  phi=acc_pi2*casuale();
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
  int t;
  for(t=0; t<sizeh; t++) {
    vect->c0[t]=d_complex_gauss();
    vect->c1[t]=d_complex_gauss();
    vect->c2[t]=d_complex_gauss();
  }
}

void generate_MultiFermion_gauss( __restrict ACC_MultiFermion * const vect){
  int t;
  int ips;
  
  for(ips=0;ips<no_ps;ips++){
    for(t=0; t<sizeh; t++) {
      vect->multi[ips].c0[t] = d_complex_gauss();
      vect->multi[ips].c1[t] = d_complex_gauss();
      vect->multi[ips].c2[t] = d_complex_gauss();
    }
  }
}

void generate_Momenta_gauss_comp(__restrict thmat_soa * const mom){
  int i,t;
  double casuali[8], aux[2];
  double uno_su_radice_di_tre = 1.0/sqrt(3.0);
  for(i=0; i<4; i++)
    {
      two_double_gauss(aux);
      casuali[2*i]   = aux[0];
      casuali[2*i+1] = aux[1];
    }

  for(t=0; t<sizeh; t++) {
  mom->rc00[t] =  casuali[2] + casuali[7] * uno_su_radice_di_tre;
  mom->rc11[t] = -casuali[2] + casuali[7] * uno_su_radice_di_tre;
  mom->c01[t]  =  casuali[0] - casuali[1] * I;
  mom->c02[t]  =  casuali[3] - casuali[4] * I;
  mom->c12[t]  =  casuali[5] - casuali[6] * I;
  }  // t  

}

void generate_Momenta_gauss(__restrict thmat_soa * const mom){
  int mu;
  for(mu=0; mu<8; mu++){
    generate_Momenta_gauss_comp(&mom[mu]);
  }
}



void generate_vec3_soa_z2noise(__restrict vec3_soa * const vect){
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


