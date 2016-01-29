#ifndef SU3_UTILITIES_C_
#define SU3_UTILITIES_C_

#include "../Include/common_defines.h"
#include "./geometry.h"
#ifdef __GNUC__
 #include <math.h>
#else // assuming PGI is used to compile on accelerators
 #include <accelmath.h>
#endif
#include "./su3_utilities.h"
#include "./single_types.h"
#include "../DbgTools/debug_macros_glvarcheck.h"

/*
// multiply the whole configuration for the staggered phases field
// (only the first two lines)
void mult_conf_times_stag_phases( __restrict su3_soa * const u){
  int hx,y,z,t,idxh;
#pragma acc kernels present(u)
#pragma acc loop independent gang // gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector // gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector // //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector // vector(DIM_BLOCK_X)
        for(hx=0; hx < nxh; hx++) {
          int x,eta;

	  //even sites
	  x = 2*hx + ((y+z+t) & 0x1);
	  idxh = snum_acc(x,y,z,t);
	  // dir  0  =  x even   --> eta = 1 , no multiplication needed
	  // dir  2  =  y even
	  eta = 1 - ( 2*(x & 0x1) );
	  mat1_times_int_factor(&u[2], idxh, eta);
	  // dir  4  =  z even
	  eta = 1 - ( 2*((x+y) & 0x1) );
	  mat1_times_int_factor(&u[4], idxh, eta);
	  // dir  6  =  t even
	  eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
	  eta *= (1- 2*(int)(t/(nt-1)));
#endif
	  mat1_times_int_factor(&u[6], idxh, eta);

	  //odd sites
	  x = 2*hx + ((y+z+t+1) & 0x1);
	  idxh = snum_acc(x,y,z,t);	  
	  // dir  1  =  x odd    --> eta = 1 , no multiplication needed
	  // dir  3  =  y odd
	  eta = 1 - ( 2*(x & 0x1) );
	  mat1_times_int_factor(&u[3], idxh, eta);
	  // dir  5  =  z odd
	  eta = 1 - ( 2*((x+y) & 0x1) );
	  mat1_times_int_factor(&u[5], idxh, eta);
	  // dir  7  =  t odd
	  eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
	  eta *= (1- 2*(int)(t/(nt-1)));
#endif
	  mat1_times_int_factor(&u[7], idxh, eta);	  

	  
	}
      }
    }
  }

}
// multiply the whole configuration for the staggered phases field
// (all three lines)
void mult_gl3_soa_times_stag_phases( __restrict su3_soa * const u){
  int hx,y,z,t,idxh;
#pragma acc kernels present(u)
#pragma acc loop independent gang // gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector // gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector // //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector // vector(DIM_BLOCK_X)
        for(hx=0; hx < nxh; hx++) {
          int x,eta;

	  //even sites
	  x = 2*hx + ((y+z+t) & 0x1);
	  idxh = snum_acc(x,y,z,t);
	  // dir  0  =  x even   --> eta = 1 , no multiplication needed
	  // dir  2  =  y even
	  eta = 1 - ( 2*(x & 0x1) );
	  gl3_times_int_factor(&u[2], idxh, eta);
	  // dir  4  =  z even
	  eta = 1 - ( 2*((x+y) & 0x1) );
	  gl3_times_int_factor(&u[4], idxh, eta);
	  // dir  6  =  t even
	  eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
	  eta *= (1- 2*(int)(t/(nt-1)));
#endif
	  gl3_times_int_factor(&u[6], idxh, eta);

	  //odd sites
	  x = 2*hx + ((y+z+t+1) & 0x1);
	  idxh = snum_acc(x,y,z,t);
	  // dir  1  =  x odd    --> eta = 1 , no multiplication needed
	  // dir  3  =  y odd
	  eta = 1 - ( 2*(x & 0x1) );
	  gl3_times_int_factor(&u[3], idxh, eta);
	  // dir  5  =  z odd
	  eta = 1 - ( 2*((x+y) & 0x1) );
	  gl3_times_int_factor(&u[5], idxh, eta);
	  // dir  7  =  t odd
	  eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
	  eta *= (1- 2*(int)(t/(nt-1)));
#endif
	  gl3_times_int_factor(&u[7], idxh, eta);	  

	  
	}
      }
    }
  }

}
// multiply the whole configuration for the staggered phases field
void mult_conf_times_stag_phases_nodev( __restrict su3_soa * const u){
  int hx,y,z,t,idxh;
  for(t=0; t<nt; t++) {
    for(z=0; z<nz; z++) {
      for(y=0; y<ny; y++) {
        for(hx=0; hx < nxh; hx++) {
          int x,eta;

	  //even sites
	  x = 2*hx + ((y+z+t) & 0x1);
	  idxh = snum_acc(x,y,z,t);
	  // dir  0  =  x even   --> eta = 1 , no multiplication needed
	  // dir  2  =  y even
	  eta = 1 - ( 2*(x & 0x1) );
	  mat1_times_int_factor(&u[2], idxh, eta);
	  // dir  4  =  z even
	  eta = 1 - ( 2*((x+y) & 0x1) );
	  mat1_times_int_factor(&u[4], idxh, eta);
	  // dir  6  =  t even
	  eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
	  eta *= (1- 2*(int)(t/(nt-1)));
#endif
	  mat1_times_int_factor(&u[6], idxh, eta);

	  //odd sites
	  x = 2*hx + ((y+z+t+1) & 0x1);
	  idxh = snum_acc(x,y,z,t);
	  // dir  1  =  x odd    --> eta = 1 , no multiplication needed
	  // dir  3  =  y odd
	  eta = 1 - ( 2*(x & 0x1) );
	  mat1_times_int_factor(&u[3], idxh, eta);
	  // dir  5  =  z odd
	  eta = 1 - ( 2*((x+y) & 0x1) );
	  mat1_times_int_factor(&u[5], idxh, eta);
	  // dir  7  =  t odd
	  eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
	  eta *= (1- 2*(int)(t/(nt-1)));
#endif
	  mat1_times_int_factor(&u[7], idxh, eta);
	  
	}
      }
    }
  }

}
*/


// reunitarize the conf by brute force
void unitarize_conf( __restrict su3_soa * const u){
  int idxh, dirindex;
#pragma acc kernels present(u)
#pragma acc loop independent
  for(dirindex = 0 ; dirindex < 8 ; dirindex++){
#pragma acc loop independent
    for( idxh = 0 ; idxh < sizeh; idxh++){
      loc_unitarize_conf(&u[dirindex],idxh);
    }
  }
}
void set_su3_soa_to_zero( __restrict su3_soa * const matrix){
  SETINUSE(matrix);
  int hx, y, z, t;
#pragma acc kernels present(matrix)
#pragma acc loop independent gang //gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector //vector(DIM_BLOCK_X)
	for(hx=0; hx < nxh; hx++) {
	  int x,idxh;
          x = 2*hx + ((y+z+t) & 0x1);
          idxh = snum_acc(x,y,z,t);
	  assign_zero_to_su3_soa_component(&matrix[0],idxh);
	  assign_zero_to_su3_soa_component(&matrix[1],idxh);
	  assign_zero_to_su3_soa_component(&matrix[2],idxh);
	  assign_zero_to_su3_soa_component(&matrix[3],idxh);
	  assign_zero_to_su3_soa_component(&matrix[4],idxh);
	  assign_zero_to_su3_soa_component(&matrix[5],idxh);
	  assign_zero_to_su3_soa_component(&matrix[6],idxh);
	  assign_zero_to_su3_soa_component(&matrix[7],idxh);
	}
      }
    }
  }
}
//copy matrix in into matrix out, this has to happen on the host
void set_su3_soa_to_su3_soa( __restrict su3_soa * const matrix_in,
			     __restrict su3_soa * const matrix_out){
  int hx, y, z, t;
  for(t=0; t<nt; t++) {
    for(z=0; z<nz; z++) {
      for(y=0; y<ny; y++) {
	for(hx=0; hx < nxh; hx++) {
	  int x,idxh;
          x = 2*hx + ((y+z+t) & 0x1);
          idxh = snum_acc(x,y,z,t);
	  assign_su3_soa_to_su3_soa_component(&matrix_in[0],&matrix_out[0],idxh);
	  assign_su3_soa_to_su3_soa_component(&matrix_in[1],&matrix_out[1],idxh);
	  assign_su3_soa_to_su3_soa_component(&matrix_in[2],&matrix_out[2],idxh);
	  assign_su3_soa_to_su3_soa_component(&matrix_in[3],&matrix_out[3],idxh);
	  assign_su3_soa_to_su3_soa_component(&matrix_in[4],&matrix_out[4],idxh);
	  assign_su3_soa_to_su3_soa_component(&matrix_in[5],&matrix_out[5],idxh);
	  assign_su3_soa_to_su3_soa_component(&matrix_in[6],&matrix_out[6],idxh);
	  assign_su3_soa_to_su3_soa_component(&matrix_in[7],&matrix_out[7],idxh);
	}
      }
    }
  }
}
// tamattamat
void conf_times_staples_ta_part(__restrict su3_soa * const u,        // constant --> is not updated
			        __restrict su3_soa * const loc_stap, // constant --> is not updated
				__restrict tamat_soa * const tipdot)
{

  int x, y, z, t;
#pragma acc kernels present(u) present(loc_stap) present(tipdot)
#pragma acc loop independent gang //gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector //vector(DIM_BLOCK_X)
	for(x=0; x < nx; x++) {
	  int idxh;
	  int parity;
	  int dir_link;
	  int mu;
	  idxh = snum_acc(x,y,z,t);  // r 
	  parity = (x+y+z+t) % 2;
	  for(mu=0;mu<4;mu++){ 
	    dir_link = 2*mu + parity;
	    mat1_times_mat2_into_tamat3(&u[dir_link],idxh,&loc_stap[dir_link],idxh,&tipdot[dir_link],idxh);

	  }

	}  // x
      }  // y
    }  // z
  }  // t

}// closes routine
// tamattamat
void RHO_times_conf_times_staples_ta_part(__restrict su3_soa * const u,        // constant --> is not updated
					  __restrict su3_soa * const loc_stap, // constant --> is not updated
					  __restrict tamat_soa * const tipdot)
{

    SETINUSE(tipdot);
  int x, y, z, t;
#pragma acc kernels present(u) present(loc_stap) present(tipdot)
#pragma acc loop independent gang //gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector //vector(DIM_BLOCK_X)
	for(x=0; x < nx; x++) {
	  int idxh;
	  int parity;
	  int dir_link;
	  int mu;
	  idxh = snum_acc(x,y,z,t);  // r 
	  parity = (x+y+z+t) % 2;
	  for(mu=0;mu<4;mu++){ 
	    dir_link = 2*mu + parity;
	    RHO_times_mat1_times_mat2_into_tamat3(&u[dir_link],idxh,&loc_stap[dir_link],idxh,&tipdot[dir_link],idxh);

	  }

	}  // x
      }  // y
    }  // z
  }  // t

}// closes routine
void mom_sum_mult( __restrict thmat_soa * const mom,
		   const __restrict tamat_soa * const ipdot,
		   double * factor,
		   int id_factor)
{
  // !!!!!!!!!!!!!!!  factor  is  equal to   -beta/3.0*timestep !!!!!!!!!!!!!!!!!!!!
  int x, y, z, t;
#pragma acc kernels present(mom) present(ipdot) present(factor)
#pragma acc loop independent //gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent //vector(DIM_BLOCK_X)
	for(x=0; x < nx; x++) {
	  int idxh;
	  int parity;
	  int dir_link;
	  int mu;
	  idxh = snum_acc(x,y,z,t);  // r 
	  parity = (x+y+z+t) % 2;
	  for(mu=0;mu<4;mu++){ 
	    dir_link = 2*mu + parity;
            thmat1_plus_tamat2_times_factor_into_thmat1(&mom[dir_link],&ipdot[dir_link],idxh,factor[id_factor]);
	  }
	}  // x
      }  // y
    }  // z
  }  // t
}// closes routine
void kernel_acc_mom_exp_times_conf( __restrict su3_soa * const conf,
				    thmat_soa * const mom, // e' costante e qui dentro non viene modificato
				    double * factor, // questo e' il vettore delta dove sono contenuti tutti i dt richiesti nell'omelyan
				    int id_factor)
{

  int x, y, z, t;
#pragma acc kernels present(mom) present(conf) present(factor)
#pragma acc loop independent gang //gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector //vector(DIM_BLOCK_X)
        for(x=0; x < nx; x++) {
          int idxh;
          int parity;
          int dir_link;
          int mu;
	  single_su3 mom_aux[1];
	  single_su3 expo[1];
	  single_su3 aux[1];
          idxh = snum_acc(x,y,z,t);  // r
          parity = (x+y+z+t) % 2;
	  for(mu=0;mu<4;mu++){
	    dir_link = 2*mu + parity;

	    extract_mom(&mom[dir_link],idxh,factor[id_factor],&mom_aux[0]);
	    matrix_exp_openacc(&mom_aux[0],&aux[0],&expo[0]);
	    conf_left_exp_multiply(&conf[dir_link],idxh,&expo[0],&aux[0],&mom_aux[0]);
	    project_on_su3(&conf[dir_link],idxh,&mom_aux[0]);
	  }
	  
        }  // x
      }  // y
    }  // z
  }  // t

}
void mom_exp_times_conf_soloopenacc(
        __restrict  su3_soa * const tconf_acc,
 thmat_soa * const tmomenta, // e' costante e qui dentro non viene modificata
 double * tdelta,  int id_delta)
{
  kernel_acc_mom_exp_times_conf(tconf_acc,tmomenta,tdelta, id_delta);
}
#endif





