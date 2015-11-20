#ifndef BACKFIELD_C_
#define BACKFIELD_C_

#include "./geometry.h"
#include "./backfield.h"
#ifdef __GNUC__
 #define acc_twopi 2*3.14159265358979323846
#endif

void init_backfield(double_soa * tu1_back_field_phases){
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
		tu1_back_field_phases[parity].d[idxh]= acc_twopi*arg; 
	      }
	      else {
		arg = -((double)(z+1))*((double)by_quantum)/(((double)nz)*((double)nx));
		tu1_back_field_phases[parity].d[idxh]= acc_twopi*arg; 
	      }
	      
	      ////////Y-oriented/////////
	      // Y even --> dir=2
	      // Y odd  --> dir=3
	      if(Y == ny){             // y-oriented links on the boundary
		arg = ((double)(z+1))*((double)ny)*((double)bx_quantum)/(((double)ny)*((double)nz));
		arg += ((double)(t+1))*((double)ny)*((double)ey_quantum)/(((double)ny)*((double)nt));
		arg -= ((double)(x+1))*((double)bz_quantum)/(((double)nx)*((double)ny));
		tu1_back_field_phases[2+parity].d[idxh]= acc_twopi*arg;
	      }
	      else {
		arg = -((double)(x+1))*((double)bz_quantum)/(((double)nx)*((double)ny));
		tu1_back_field_phases[2+parity].d[idxh]= acc_twopi*arg;
	      }

	      ////////Z-oriented////////
	      // Z even --> dir=4
	      // Z odd  --> dir=5
	      if(Z == nz){          // z-oriented links on the boundary
		arg = ((double)(t+1))*((double)nz)*((double)ez_quantum)/(((double)nz)*((double)nt));
		arg += ((double)(x+1))*((double)nz)*((double)by_quantum)/(((double)nz)*((double)nx));
		arg -= ((double)(y+1))*((double)bx_quantum)/(((double)ny)*((double)nz));
		tu1_back_field_phases[4+parity].d[idxh]=acc_twopi* arg;
	      }
	      else{
		arg = -((double)(y+1))*((double)bx_quantum)/(((double)ny)*((double)nz));
		tu1_back_field_phases[4+parity].d[idxh]= acc_twopi*arg;
	      }

	      ///////T-oriented////////
	      // T even --> dir=6
	      // T odd  --> dir=7
	      arg = -((double)(z+1))*((double)ez_quantum)/(((double)nz)*((double)nt));
	      arg -= ((double)(y+1))*((double)ey_quantum)/(((double)ny)*((double)nt));
	      arg -= ((double)(x+1))*((double)ex_quantum)/(((double)nx)*((double)nt));
	      tu1_back_field_phases[6+parity].d[idxh]= acc_twopi*arg;
	  } // x loop
	} // y loop
      } // z loop
    } // t loop
}


#endif
