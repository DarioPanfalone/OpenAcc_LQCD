#ifndef SU3_UTILITIES_C_
#define SU3_UTILITIES_C_

#include "../Include/common_defines.h"
#include "./geometry.h"
#ifdef __GNUC__
 #include <math.h>
#else // assuming PGI is used to compile on accelerators
 #include <openacc.h>
 #include <accelmath.h>
#endif
#include "./su3_utilities.h"
#include "./struct_c_def.h"
#include "./single_types.h"

extern int TOPO_GLOBAL_DONT_TOUCH;

// reunitarize the conf by brute force
void unitarize_conf( __restrict su3_soa * const u)
{
    int dirindex;
    int d0h, d1, d2, d3;
#pragma acc kernels present(u)
#pragma acc loop independent gang
    for(dirindex = 0 ; dirindex < 8 ; dirindex++){
#pragma acc loop independent gang(STAPGANG3)
        for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) { 
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
            for(d2=0; d2<nd2; d2++) {
                for(d1=0; d1<nd1; d1++) {
                    for(d0h=0; d0h < nd0h; d0h++) {
                        // I take the size to be even, but it's the same
                        int d0 = 2*d0h + ((d1+d2+d3) & 0x1);
                        int t = snum_acc(d0,d1,d2,d3);  
                        loc_unitarize_conf(&u[dirindex],t);

                    }
                }
            }
        }
    }
}


void set_su3_soa_to_zero( __restrict su3_soa * const matrix)
{
  int hd0, d1, d2, d3;
#pragma acc kernels present(matrix)
#pragma acc loop independent gang(STAPGANG3)
  for(d3=0; d3<nd3; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
	for(hd0=0; hd0 < nd0h; hd0++) {
	  int d0,idxh;
          d0 = 2*hd0 + ((d1+d2+d3) & 0x1);
          idxh = snum_acc(d0,d1,d2,d3);
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
void set_su3_soa_to_su3_soa( __restrict const su3_soa * const matrix_in,
			     __restrict su3_soa * const matrix_out)
{
  int hd0, d1, d2, d3;
  for(d3=0; d3<nd3; d3++) {
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
	for(hd0=0; hd0 < nd0h; hd0++) {
	  int d0,idxh;
          d0 = 2*hd0 + ((d1+d2+d3) & 0x1);
          idxh = snum_acc(d0,d1,d2,d3);
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

void set_su3_soa_to_su3_soa_sp( __restrict const su3_soa * const matrix_in,
                            __restrict su3_soa * const matrix_out)
{
    int hd0, d1, d2, d3;
    for(d3=0; d3<nd3; d3++) {
        for(d2=0; d2<nd2; d2++) {
            for(d1=0; d1<nd1; d1++) {
                for(hd0=0; hd0 < nd0h; hd0++) {
                    int d0,idxh;
                    d0 = 2*hd0 + ((d1+d2+d3) & 0x1);
                    idxh = snum_acc(d0,d1,d2,d3);
                    assign_su3_soa_to_su3_soa_component_sp(&matrix_in[0],&matrix_out[0],idxh);
                    assign_su3_soa_to_su3_soa_component_sp(&matrix_in[1],&matrix_out[1],idxh);
                    assign_su3_soa_to_su3_soa_component_sp(&matrix_in[2],&matrix_out[2],idxh);
                    assign_su3_soa_to_su3_soa_component_sp(&matrix_in[3],&matrix_out[3],idxh);
                    assign_su3_soa_to_su3_soa_component_sp(&matrix_in[4],&matrix_out[4],idxh);
                    assign_su3_soa_to_su3_soa_component_sp(&matrix_in[5],&matrix_out[5],idxh);
                    assign_su3_soa_to_su3_soa_component_sp(&matrix_in[6],&matrix_out[6],idxh);
                    assign_su3_soa_to_su3_soa_component_sp(&matrix_in[7],&matrix_out[7],idxh);
                }
            }
        }
    }
}




void set_su3_soa_to_su3_soa_trasl( __restrict const su3_soa * const matrix_in,
                            __restrict su3_soa * const matrix_out,int dir)
{
    int hd0, d1, d2, d3;
    for(d3=0; d3<nd3; d3++) {
        for(d2=0; d2<nd2; d2++) {
            for(d1=0; d1<nd1; d1++) {
                for(hd0=0; hd0 < nd0h; hd0++) {
                    int d0,idxh,parity,idxpmu;
                    d0 = 2*hd0 + ((d1+d2+d3) & 0x1);
                    idxh = snum_acc(d0,d1,d2,d3);
                     parity=(d0+d1+d2+d3)%2;
                    idxpmu=nnp_openacc[idxh][dir][parity];
                    assign_su3_soa_to_su3_soa_component_trasl(&matrix_in[0],&matrix_out[1],idxh,idxpmu);
                    assign_su3_soa_to_su3_soa_component_trasl(&matrix_in[1],&matrix_out[0],idxh,idxpmu);
                    assign_su3_soa_to_su3_soa_component_trasl(&matrix_in[2],&matrix_out[3],idxh,idxpmu);
                    assign_su3_soa_to_su3_soa_component_trasl(&matrix_in[3],&matrix_out[2],idxh,idxpmu);
                    assign_su3_soa_to_su3_soa_component_trasl(&matrix_in[4],&matrix_out[5],idxh,idxpmu);
                    assign_su3_soa_to_su3_soa_component_trasl(&matrix_in[5],&matrix_out[4],idxh,idxpmu);
                    assign_su3_soa_to_su3_soa_component_trasl(&matrix_in[6],&matrix_out[7],idxh,idxpmu);
                    assign_su3_soa_to_su3_soa_component_trasl(&matrix_in[7],&matrix_out[6],idxh,idxpmu);
                }
            }
        }
    }
}

void set_su3_soa_to_su3_soa_device(__restrict const su3_soa * const matrix_in,
			     	   __restrict su3_soa * const matrix_out)
{
  int hd0, d1, d2, d3;
#pragma acc kernels present(matrix_in) present(matrix_out)
#pragma acc loop independent gang(STAPGANG3)
  for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
	for(hd0=0; hd0 < nd0h; hd0++) {
	  int d0,idxh;
          d0 = 2*hd0 + ((d1+d2+d3) & 0x1);
          idxh = snum_acc(d0,d1,d2,d3);
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
void conf_times_staples_ta_part(
        __restrict const su3_soa * const u,
        __restrict const su3_soa * const loc_stap,
        __restrict tamat_soa * const tipdot)
{

  int d0, d1, d2, d3;
#pragma acc kernels present(u) present(loc_stap) present(tipdot)
#pragma acc loop independent gang(STAPGANG3)
  for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
	for(d0=0; d0 < nd0; d0++) {
	  int idxh;
	  int parity;
	  int dir_link;
	  int mu;
	  idxh = snum_acc(d0,d1,d2,d3);  // r 
	  parity = (d0+d1+d2+d3) % 2;
	  for(mu=0;mu<4;mu++){ 
	    dir_link = 2*mu + parity;
	    mat1_times_mat2_into_tamat3(&u[dir_link],idxh,&loc_stap[dir_link],idxh,&tipdot[dir_link],idxh);

	  }

	}  // d0
      }  // d1
    }  // d2
  }  // d3

}// closes routine

// tamattamat
void conf_times_staples_ta_part_addto_tamat(
        __restrict const su3_soa * const u,
        __restrict const su3_soa * const loc_stap,
        __restrict tamat_soa * const tipdot)
{

  int d0, d1, d2, d3;
#pragma acc kernels present(u) present(loc_stap) present(tipdot)
#pragma acc loop independent gang(STAPGANG3)
  for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
	for(d0=0; d0 < nd0; d0++) {
	  int idxh;
	  int parity;
	  int dir_link;
	  int mu;
	  idxh = snum_acc(d0,d1,d2,d3);  // r 
	  parity = (d0+d1+d2+d3) % 2;
	  for(mu=0;mu<4;mu++){ 
	    dir_link = 2*mu + parity;
	    mat1_times_mat2_addto_tamat3(&u[dir_link],idxh,&loc_stap[dir_link],idxh,&tipdot[dir_link],idxh);

	  }

	}  // d0
      }  // d1
    }  // d2
  }  // d3

}// closes routine
// tamattamat
void RHO_times_conf_times_staples_ta_part(
        __restrict const su3_soa * const u,
        __restrict const su3_soa * const loc_stap,
        __restrict tamat_soa * const tipdot)
{

  int d0, d1, d2, d3;
#pragma acc kernels present(u) present(loc_stap) present(tipdot)
#pragma acc loop independent gang(STAPGANG3)
  for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
	for(d0=0; d0 < nd0; d0++) {
	  int idxh;
	  int parity;
	  int dir_link;
	  int mu;
	  idxh = snum_acc(d0,d1,d2,d3);  // r 
	  parity = (d0+d1+d2+d3) % 2;
	  for(mu=0;mu<4;mu++){ 
	    dir_link = 2*mu + parity;
	    RHO_times_mat1_times_mat2_into_tamat3(&u[dir_link],idxh,&loc_stap[dir_link],idxh,&tipdot[dir_link],idxh);

	  }

	}  // d0
      }  // d1
    }  // d2
  }  // d3

}// closes routine

void mom_sum_mult( __restrict thmat_soa * const mom,
		   const __restrict tamat_soa * const ipdot,
		   const double * factor,
		   int id_factor)
{
  // !!!!!!!!!!!!!!!  factor  is  equal to   -beta/3.0*timestep !!!!!!!!!!!!!!!!!!!!
  int d0, d1, d2, d3;
#pragma acc kernels present(mom) present(ipdot) present(factor)
#pragma acc loop independent gang(STAPGANG3)
  for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
	for(d0=0; d0 < nd0; d0++) {
	  int idxh;
	  int parity;
	  int dir_link;
	  int mu;
	  idxh = snum_acc(d0,d1,d2,d3);  // r 
	  parity = (d0+d1+d2+d3) % 2;
	  for(mu=0;mu<4;mu++){ 
	    dir_link = 2*mu + parity;
            thmat1_plus_tamat2_times_factor_into_thmat1(&mom[dir_link],&ipdot[dir_link],idxh,factor[id_factor]);
	  }
	}  // d0
      }  // d1
    }  // d2
  }  // d3
}// closes routine

 
void mom_exp_times_conf_soloopenacc( 
        __restrict su3_soa * const conf,
        __restrict const thmat_soa * const mom,
        const double * factor, 
        // questo e' il vettore delta dove sono contenuti 
        // tutti i dt richiesti nell'omelyan
        int id_factor)
{

  int d0, d1, d2, d3;
#pragma acc kernels present(mom) present(conf) present(factor)
#pragma acc loop independent gang(STAPGANG3)
  for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
        for(d0=0; d0 < nd0; d0++) {
          int idxh;
          int parity;
          int dir_link;
          int mu;
          //single_su3 mom_aux, expo, aux;
          idxh = snum_acc(d0,d1,d2,d3);  // r
          parity = (d0+d1+d2+d3) % 2;
	  for(mu=0;mu<4;mu++){
	    dir_link = 2*mu + parity;

        mom_exp_times_conf_soloopenacc_loc(&mom[dir_link],&conf[dir_link],idxh,factor[id_factor]);
	    //extract_mom(&mom[dir_link],idxh,factor[id_factor],&mom_aux);
	    //matrix_exp_openacc(&mom_aux,&aux,&expo);
	    //conf_left_exp_multiply(&conf[dir_link],idxh,&expo,&aux,&mom_aux);
	    //project_on_su3(&conf[dir_link],idxh,&mom_aux);
	  }
	  
        }  // d0
      }  // d1
    }  // d2
  }  // d3

}

#ifdef MULTIDEVICE

void set_su3_soa_to_zero_bulk( __restrict su3_soa * const matrix)
{
  int hd0, d1, d2, d3;
#pragma acc kernels present(matrix)
#pragma acc loop independent gang(STAPGANG3)
  for(d3=D3_HALO+GAUGE_HALO; d3<nd3-D3_HALO-GAUGE_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
	for(hd0=0; hd0 < nd0h; hd0++) {
	  int d0,idxh;
          d0 = 2*hd0 + ((d1+d2+d3) & 0x1);
          idxh = snum_acc(d0,d1,d2,d3);
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

void conf_times_staples_ta_part_bulk(
        __restrict const su3_soa * const u,
        __restrict const su3_soa * const loc_stap,
        __restrict tamat_soa * const tipdot)
{

  int d0, d1, d2, d3;
#pragma acc kernels present(u) present(loc_stap) present(tipdot)
#pragma acc loop independent gang(STAPGANG3)
  for(d3=D3_HALO+GAUGE_HALO; d3<nd3-D3_HALO-GAUGE_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
	for(d0=0; d0 < nd0; d0++) {
	  int idxh;
	  int parity;
	  int dir_link;
	  int mu;
	  idxh = snum_acc(d0,d1,d2,d3);  // r 
	  parity = (d0+d1+d2+d3) % 2;
	  for(mu=0;mu<4;mu++){ 
	    dir_link = 2*mu + parity;
	    mat1_times_mat2_into_tamat3(&u[dir_link],idxh,&loc_stap[dir_link],idxh,&tipdot[dir_link],idxh);

	  }

	}  // d0
      }  // d1
    }  // d2
  }  // d3

}// closes routine

void mom_sum_mult_bulk( __restrict thmat_soa * const mom,
		   const __restrict tamat_soa * const ipdot,
		   const double * factor,
		   int id_factor)
{
  // !!!!!!!!!!!!!!!  factor  is  equal to   -beta/3.0*timestep !!!!!!!!!!!!!!!!!!!!
  int d0, d1, d2, d3;
#pragma acc kernels present(mom) present(ipdot) present(factor)
#pragma acc loop independent gang(STAPGANG3)
  for(d3=D3_HALO+GAUGE_HALO; d3<nd3-D3_HALO-GAUGE_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
	for(d0=0; d0 < nd0; d0++) {
	  int idxh;
	  int parity;
	  int dir_link;
	  int mu;
	  idxh = snum_acc(d0,d1,d2,d3);  // r 
	  parity = (d0+d1+d2+d3) % 2;
	  for(mu=0;mu<4;mu++){ 
	    dir_link = 2*mu + parity;
            thmat1_plus_tamat2_times_factor_into_thmat1(&mom[dir_link],&ipdot[dir_link],idxh,factor[id_factor]);
	  }
	}  // d0
      }  // d1
    }  // d2
  }  // d3
}// closes routine

 
void mom_exp_times_conf_soloopenacc_bulk( 
        __restrict const su3_soa * const conf_old,
        __restrict su3_soa * const conf_new,
        __restrict const thmat_soa * const mom,
        const double * factor, 
        // questo e' il vettore delta dove sono contenuti 
        // tutti i dt richiesti nell'omelyan
        int id_factor)
{

  int d0, d1, d2, d3;
#pragma acc kernels present(mom) present(conf_old) present(conf_new) present(factor)
#pragma acc loop independent gang(STAPGANG3)
  for(d3=D3_HALO+GAUGE_HALO; d3<nd3-D3_HALO-GAUGE_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
        for(d0=0; d0 < nd0; d0++) {
          int idxh;
          int parity;
          int dir_link;
          int mu;
	  //single_su3 mom_aux, expo, aux;
          idxh = snum_acc(d0,d1,d2,d3);  // r
          parity = (d0+d1+d2+d3) % 2;
	  for(mu=0;mu<4;mu++){
	    dir_link = 2*mu + parity;

        mom_exp_times_conf_soloopenacc_loc_split(&mom[dir_link],
                &conf_old[dir_link],
                &conf_new[dir_link],
                idxh,factor[id_factor]);
	    //extract_mom(&mom[dir_link],idxh,factor[id_factor],&mom_aux);
	    //matrix_exp_openacc(&mom_aux,&aux,&expo);
	    //conf_left_exp_multiply(&conf_old[dir_link],idxh,&expo,&aux,&mom_aux);
	    //project_on_su3(&conf_new[dir_link],idxh,&mom_aux);
	  }
	  
        }  // d0
      }  // d1
    }  // d2
  }  // d3

}

void set_su3_soa_to_zero_d3c( __restrict su3_soa * const matrix,
        int offset3, int thickness3)
{
  int hd0, d1, d2, d3;
#pragma acc kernels present(matrix)
#pragma acc loop independent gang
  for(d3=offset3; d3<offset3+thickness3; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
	for(hd0=0; hd0 < nd0h; hd0++) {
	  int d0,idxh;
          d0 = 2*hd0 + ((d1+d2+d3) & 0x1);
          idxh = snum_acc(d0,d1,d2,d3);
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

void conf_times_staples_ta_part_d3c(
        __restrict const su3_soa * const u,
        __restrict const su3_soa * const loc_stap,
        __restrict tamat_soa * const tipdot,
        int offset3, int thickness3)
{

  int d0, d1, d2, d3;
#pragma acc kernels present(u) present(loc_stap) present(tipdot)
#pragma acc loop independent gang
  for(d3=offset3; d3<offset3+thickness3; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
	for(d0=0; d0 < nd0; d0++) {
	  int idxh;
	  int parity;
	  int dir_link;
	  int mu;
	  idxh = snum_acc(d0,d1,d2,d3);  // r 
	  parity = (d0+d1+d2+d3) % 2;
	  for(mu=0;mu<4;mu++){ 
	    dir_link = 2*mu + parity;
	    mat1_times_mat2_into_tamat3(&u[dir_link],idxh,&loc_stap[dir_link],idxh,&tipdot[dir_link],idxh);

	  }

	}  // d0
      }  // d1
    }  // d2
  }  // d3

}// closes routine

void mom_sum_mult_d3c( __restrict thmat_soa * const mom,
		   const __restrict tamat_soa * const ipdot,
		   const double * factor,
		   int id_factor,
           int offset3, int thickness3)
{
  // !!!!!!!  factor  is  equal to   -beta/3.0*timestep !!!!!!!!!!!!!!!!
  int d0, d1, d2, d3;
#pragma acc kernels present(mom) present(ipdot) present(factor)
#pragma acc loop independent gang
  for(d3=offset3; d3<offset3+thickness3; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
	for(d0=0; d0 < nd0; d0++) {
	  int idxh;
	  int parity;
	  int dir_link;
	  int mu;
	  idxh = snum_acc(d0,d1,d2,d3);  // r 
	  parity = (d0+d1+d2+d3) % 2;
	  for(mu=0;mu<4;mu++){ 
	    dir_link = 2*mu + parity;
            thmat1_plus_tamat2_times_factor_into_thmat1(&mom[dir_link],&ipdot[dir_link],idxh,factor[id_factor]);
	  }
	}  // d0
      }  // d1
    }  // d2
  }  // d3
}// closes routine

void mom_exp_times_conf_soloopenacc_d3c( 
        __restrict const su3_soa * const conf_old,
        __restrict su3_soa * const conf_new,
        __restrict const thmat_soa * const mom,
        const double * factor, 
        // questo e' il vettore delta dove sono contenuti 
        // tutti i dt richiesti nell'omelyan
        int id_factor,
        int offset3, int thickness3)
{

  int d0, d1, d2, d3;
#pragma acc kernels present(mom) present(conf_old) present(conf_new) present(factor)
#pragma acc loop independent gang 
  for(d3=offset3; d3<offset3+thickness3; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
        for(d0=0; d0 < nd0; d0++) {
          int idxh;
          int parity;
          int dir_link;
          int mu;
	  //single_su3 mom_aux, expo, aux;
          idxh = snum_acc(d0,d1,d2,d3);  // r
          parity = (d0+d1+d2+d3) % 2;
	  for(mu=0;mu<4;mu++){
	    dir_link = 2*mu + parity;

        mom_exp_times_conf_soloopenacc_loc_split(&mom[dir_link],
                &conf_old[dir_link],
                &conf_new[dir_link],
                idxh,factor[id_factor]);
	    //extract_mom(&mom[dir_link],idxh,factor[id_factor],&mom_aux);
	    //matrix_exp_openacc(&mom_aux,&aux,&expo);
	    //conf_left_exp_multiply(&conf_old[dir_link],idxh,&expo,&aux,&mom_aux);
	    //project_on_su3(&conf_new[dir_link],idxh,&mom_aux);
	  }
	  
        }  // d0
      }  // d1
    }  // d2
  }  // d3

}

#endif

#endif





