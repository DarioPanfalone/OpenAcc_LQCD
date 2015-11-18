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

#include "../OpenAcc/su3_utilities.c"

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

void generate_Momenta_gauss_comp(__restrict thmat_soa * const mom){
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

void generate_Conf_cold_comp(__restrict su3_soa * const conf,double factor){
  d_complex mat00,mat01,mat02 ;
  d_complex mat10,mat11,mat12 ;
  d_complex mat20,mat21,mat22 ;
  double norm;
  d_complex prod;
  int t;

  ///////////////////////////////////////////////////////////
  // LOGICAL STEPS FOR THE GENERATION OF THE COLD CONF  /////
  ///////////////////////////////////////////////////////////
  // A) extract randomly the first row and normalize it
  // B) extract randomly the second row, compute the scalar prod with the first
  // C) subtract from the second row the projection of the first on the second
  // D) normalize the second row
  // E) compute the third row as the vector product of the first two
  // F) add the identity to the resulting matrix (times a small number) and re-unitarize
  // G) multiply for the staggered phase
  /////////////////////////////////////////////////////////////
  for(t=0; t<sizeh; t++) {
    ///////////GENERATE RANDOM MATRIX //////////////////////////////////////////////////////
    mat00 = 1.0-2.0*casuale()+I*(1.0-2.0*casuale());
    mat01 = 1.0-2.0*casuale()+I*(1.0-2.0*casuale());
    mat02 = 1.0-2.0*casuale()+I*(1.0-2.0*casuale());
    mat10 = 1.0-2.0*casuale()+I*(1.0-2.0*casuale());
    mat11 = 1.0-2.0*casuale()+I*(1.0-2.0*casuale());
    mat12 = 1.0-2.0*casuale()+I*(1.0-2.0*casuale());
    
    ////////// Unitarize it ////////////////////////////////////////////////////////////////
    norm  = creal(mat00)*creal(mat00) +  cimag(mat00)*cimag(mat00);
    norm += creal(mat01)*creal(mat01) +  cimag(mat01)*cimag(mat01);
    norm += creal(mat02)*creal(mat02) +  cimag(mat02)*cimag(mat02);
    norm=1.0/sqrt(norm);
    
    mat00 *= norm;
    mat01 *= norm;
    mat02 *= norm;
    
    prod  = conj(mat00) * mat10 + conj(mat01) * mat11 + conj(mat02) * mat12;
    
    mat10 -= prod * mat00;
    mat11 -= prod * mat01;
    mat12 -= prod * mat02;
    
    norm  = creal(mat10)*creal(mat10) +  cimag(mat10)*cimag(mat10);
    norm += creal(mat11)*creal(mat11) +  cimag(mat11)*cimag(mat11);
    norm += creal(mat12)*creal(mat12) +  cimag(mat12)*cimag(mat12);
    norm=1.0/sqrt(norm);
    
    mat10 *= norm;
    mat11 *= norm;
    mat22 *= norm;
    
    mat20 = conj(mat01*mat12-mat02*mat11);
    mat21 = conj(mat02*mat10-mat00*mat12);
    mat22 = conj(mat00*mat11-mat01*mat10);
    
    ////////multiply the matrix times a factor and add the identity///////////////////////
    mat00 = mat00 * factor + 1.0;
    mat01 = mat01 * factor ;
    mat02 = mat02 * factor ;
    mat10 = mat10 * factor ;
    mat11 = mat11 * factor + 1.0;
    mat12 = mat12 * factor ;
    mat20 = mat20 * factor ;
    mat21 = mat21 * factor ;
    mat22 = mat22 * factor + 1.0;
    
    /////////// re-unitarize the matrix //////////////////////////////////////////////////
    norm  = creal(mat00)*creal(mat00) +  cimag(mat00)*cimag(mat00);
    norm += creal(mat01)*creal(mat01) +  cimag(mat01)*cimag(mat01);
    norm += creal(mat02)*creal(mat02) +  cimag(mat02)*cimag(mat02);
    norm=1.0/sqrt(norm);
    
    mat00 *= norm;
    mat01 *= norm;
    mat02 *= norm;
    
    prod  = conj(mat00) * mat10 + conj(mat01) * mat11 + conj(mat02) * mat12;
    
    mat10 -= prod * mat00;
    mat11 -= prod * mat01;
    mat12 -= prod * mat02;
    
    norm  = creal(mat10)*creal(mat10) +  cimag(mat10)*cimag(mat10);
    norm += creal(mat11)*creal(mat11) +  cimag(mat11)*cimag(mat11);
    norm += creal(mat12)*creal(mat12) +  cimag(mat12)*cimag(mat12);
    norm=1.0/sqrt(norm);
    
    mat10 *= norm;
    mat11 *= norm;
    mat22 *= norm;
    
    mat20 = conj(mat01*mat12-mat02*mat11);
    mat21 = conj(mat02*mat10-mat00*mat12);
    mat22 = conj(mat00*mat11-mat01*mat10);
    
    /////assign the result to the array  /////////////////////
    conf->r0.c0[t] = mat00;
    conf->r0.c1[t] = mat01;
    conf->r0.c2[t] = mat02;
    conf->r1.c0[t] = mat10;
    conf->r1.c1[t] = mat11;
    conf->r1.c2[t] = mat12;
    conf->r2.c0[t] = mat20;
    conf->r2.c1[t] = mat21;
    conf->r2.c2[t] = mat22;
  }
}

// iterations of random assignements over all the 8 components
void generate_Momenta_gauss(__restrict thmat_soa * const mom){
  SETINUSE(mom);
  int mu;
  for(mu=0; mu<8; mu++){
    generate_Momenta_gauss_comp(&mom[mu]);
  }
}

// iterations of random assignements over all the 8 components
void generate_Conf_cold(__restrict su3_soa * const conf,double factor){
  SETINUSE(conf);
  int mu;
  for(mu=0; mu<8; mu++){
    generate_Conf_cold_comp(&conf[mu],factor);
  }
  mult_conf_times_stag_phases_nodev(conf);
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

/*
void init_backfield(){
  int X,T,Z,Y;
  double  arg;

  int hx, y, z, t, parity,evod;
    for(t=0; t<nt; t++) {
      for(z=0; z<nz; z++) {
	for(y=0; y<ny; y++) {
	  for(hx=0; hx < nxh; hx++) {
	    for(evod=0; evod < 2; evod++) {
	      int x, idxh, matdir;
	      if(evod==0)	    x = 2*hx + ((y+z+t) & 0x1);
	      if(evod==1)	    x = 2*hx + ((y+z+t+1) & 0x1);
	      parity = 1; //dispari
	      if((x+y+z+t)%2==0) parity = 0; //pari
	      idxh = snum_acc(x,y,z,t);
	      
	      X = x + 1;
	      Y = y + 1;
	      Z = z + 1;
	      
	      ////////X-oriented////////
	      // X even --> dir=0
	      // X odd  --> dir=1
	      if(X == nx){             // x-oriented links on the boundary
		arg = ((double)(y+1))*((double)nx)*((double)bz_quantum)/(((double)nx)*((double)ny));
		arg += ((double)(t+1))*((double)nx)*((double)ex_quantum)/(((double)nx)*((double)nt));
		arg -= ((double)(z+1))*((double)by_quantum)/(((double)nz)*((double)nx));
		u1_back_field_phases[parity].d[idxh]= arg; 
	      }
	      else {
		arg = -((double)(z+1))*((double)by_quantum)/(((double)nz)*((double)nx));
		u1_back_field_phases[parity].d[idxh]= arg; 
	      }
	      
	      ////////Y-oriented/////////
	      // Y even --> dir=2
	      // Y odd  --> dir=3
	      if(Y == ny){             // y-oriented links on the boundary
		arg = ((double)(z+1))*((double)ny)*((double)bx_quantum)/(((double)ny)*((double)nz));
		arg += ((double)(t+1))*((double)ny)*((double)ey_quantum)/(((double)ny)*((double)nt));
		arg -= ((double)(x+1))*((double)bz_quantum)/(((double)nx)*((double)ny));
		u1_back_field_phases[2+parity].d[idxh]= arg;
	      }
	      else {
		arg = -((double)(x+1))*((double)bz_quantum)/(((double)nx)*((double)ny));
		u1_back_field_phases[2+parity].d[idxh]= arg;
	      }

	      ////////Z-oriented////////
	      // Z even --> dir=4
	      // Z odd  --> dir=5
	      if(Z == nz){          // z-oriented links on the boundary
		arg = ((double)(t+1))*((double)nz)*((double)ez_quantum)/(((double)nz)*((double)nt));
		arg += ((double)(x+1))*((double)nz)*((double)by_quantum)/(((double)nz)*((double)nx));
		arg -= ((double)(y+1))*((double)bx_quantum)/(((double)ny)*((double)nz));
		u1_back_field_phases[4+parity].d[idxh]= arg;
	      }
	      else{
		arg = -((double)(y+1))*((double)bx_quantum)/(((double)ny)*((double)nz));
		u1_back_field_phases[4+parity].d[idxh]= arg;
	      }

	      ///////T-oriented////////
	      // T even --> dir=6
	      // T odd  --> dir=7
	      arg = -((double)(z+1))*((double)ez_quantum)/(((double)nz)*((double)nt));
	      arg -= ((double)(y+1))*((double)ey_quantum)/(((double)ny)*((double)nt));
	      arg -= ((double)(x+1))*((double)ex_quantum)/(((double)nx)*((double)nt));
	      u1_back_field_phases[6+parity].d[idxh]= arg;
	    } // evod loop
	  } // hx loop
	} // y loop
      } // z loop
    } // t loop
}
*/



void init_backfield(){
  int X,T,Z,Y;
  double  arg;

  int x, y, z, t, parity,idxh;
    for(t=0; t < nt; t++){
      for(z=0; z < nz; z++){
	for(y=0; y < ny; y++){
	  for(x=0; x < nx; x++){
	    
	      if((x+y+z+t)%2==0){
		parity = 0; //pari
	      }else{
		parity = 1; //dispari
	      }
	      idxh = snum_acc(x,y,z,t);
	      
	      X = x + 1;
	      Y = y + 1;
	      Z = z + 1;
	      
	      ////////X-oriented////////
	      // X even --> dir=0
	      // X odd  --> dir=1
	      if(X == nx){             // x-oriented links on the boundary
		arg = ((double)(y+1))*((double)nx)*((double)bz_quantum)/(((double)nx)*((double)ny));
		arg += ((double)(t+1))*((double)nx)*((double)ex_quantum)/(((double)nx)*((double)nt));
		arg -= ((double)(z+1))*((double)by_quantum)/(((double)nz)*((double)nx));
		u1_back_field_phases[parity].d[idxh]= acc_twopi*arg; 
	      }
	      else {
		arg = -((double)(z+1))*((double)by_quantum)/(((double)nz)*((double)nx));
		u1_back_field_phases[parity].d[idxh]= acc_twopi*arg; 
	      }
	      
	      ////////Y-oriented/////////
	      // Y even --> dir=2
	      // Y odd  --> dir=3
	      if(Y == ny){             // y-oriented links on the boundary
		arg = ((double)(z+1))*((double)ny)*((double)bx_quantum)/(((double)ny)*((double)nz));
		arg += ((double)(t+1))*((double)ny)*((double)ey_quantum)/(((double)ny)*((double)nt));
		arg -= ((double)(x+1))*((double)bz_quantum)/(((double)nx)*((double)ny));
		u1_back_field_phases[2+parity].d[idxh]= acc_twopi*arg;
	      }
	      else {
		arg = -((double)(x+1))*((double)bz_quantum)/(((double)nx)*((double)ny));
		u1_back_field_phases[2+parity].d[idxh]= acc_twopi*arg;
	      }

	      ////////Z-oriented////////
	      // Z even --> dir=4
	      // Z odd  --> dir=5
	      if(Z == nz){          // z-oriented links on the boundary
		arg = ((double)(t+1))*((double)nz)*((double)ez_quantum)/(((double)nz)*((double)nt));
		arg += ((double)(x+1))*((double)nz)*((double)by_quantum)/(((double)nz)*((double)nx));
		arg -= ((double)(y+1))*((double)bx_quantum)/(((double)ny)*((double)nz));
		u1_back_field_phases[4+parity].d[idxh]=acc_twopi* arg;
	      }
	      else{
		arg = -((double)(y+1))*((double)bx_quantum)/(((double)ny)*((double)nz));
		u1_back_field_phases[4+parity].d[idxh]= acc_twopi*arg;
	      }

	      ///////T-oriented////////
	      // T even --> dir=6
	      // T odd  --> dir=7
	      arg = -((double)(z+1))*((double)ez_quantum)/(((double)nz)*((double)nt));
	      arg -= ((double)(y+1))*((double)ey_quantum)/(((double)ny)*((double)nt));
	      arg -= ((double)(x+1))*((double)ex_quantum)/(((double)nx)*((double)nt));
	      u1_back_field_phases[6+parity].d[idxh]= acc_twopi*arg;
	  } // x loop
	} // y loop
      } // z loop
    } // t loop
}


#endif


