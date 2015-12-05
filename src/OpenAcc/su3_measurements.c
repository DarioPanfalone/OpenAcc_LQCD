#ifndef SU3_MEASUREMENTS_C_
#define SU3_MEASUREMENTS_C_

#include "../Include/common_defines.h"
#include "../OpenAcc/struct_c_def.h"
#include "./su3_measurements.h"
#include "./su3_utilities.h"
#include "./single_types.h"
#ifndef __GNUC__
#include <openacc.h>
#endif

void check_unitarity_device( __restrict su3_soa * const u, double * max_unitarity_deviation, double *avg_unitarity_deviation){


  // removing stag phases
  mult_conf_times_stag_phases(u);


    printf("HEE\n");
    double r = 0;
    double rmax = 0;

#pragma acc kernels present(u)
#pragma acc loop reduction(+:r) 
  for(int dir = 0; dir < 8 ; dir++){
  for(int idx = 0; idx < sizeh ; idx++){
     single_su3 m;
     single_su3_from_su3_soa(&u[dir],idx,&m);
     rebuild3row(&m);

     d_complex err = 1 - detSu3(&m);
     r += creal(err * conj(err));
//     rmax = fmax(rmax,creal(err * conj(err)));

    }
  }
  //adding them again
  mult_conf_times_stag_phases(u);

  *avg_unitarity_deviation = r/(sizeh*8);
  *max_unitarity_deviation = rmax;


}
void check_unitarity_host( __restrict su3_soa * const u, double * max_unitarity_deviation, double *avg_unitarity_deviation){


  // removing stag phases
  mult_conf_times_stag_phases_nodev(u);


    double r = 0;
    double rmax = 0;

  for(int dir = 0; dir < 8 ; dir++){
  for(int idx = 0; idx < sizeh ; idx++){
     single_su3 m;
     single_su3_from_su3_soa(&u[dir],idx,&m);
     rebuild3row(&m);

     d_complex err = 1 - detSu3(&m);
     r += creal(err * conj(err));
     rmax = fmax(rmax,creal(err * conj(err)));

    }
  }
  //adding them again
  mult_conf_times_stag_phases_nodev(u);

  *avg_unitarity_deviation = r/(sizeh*8);
  *max_unitarity_deviation = rmax;


}


double calc_momenta_action( const __restrict thmat_soa * const mom,
			    double_soa * tr_local,
			    const int mu){
  int t;

#pragma acc kernels present(mom) present(tr_local)
#pragma acc loop independent //gang(nt)
  for(t=0; t<sizeh; t++) {
	  tr_local[0].d[t] = half_tr_thmat_squared(&mom[mu],t);
  }  // t


  double result=0.0;
#pragma acc kernels present(tr_local)
#pragma acc loop reduction(+:result)
  for(t=0; t<sizeh; t++) {
    result += tr_local[0].d[t];
  }
  
  return result;


}// closes routine


double  calc_plaquette_soloopenacc( __restrict  su3_soa * const tconf_acc, __restrict su3_soa * const local_plaqs, dcomplex_soa * const tr_local_plaqs){

  double tempo=0.0;
  // tolgo le fasi staggered
  mult_conf_times_stag_phases(tconf_acc);
  // calcolo il valore della plaquette sommata su tutti i siti a fissato piano mu-nu (6 possibili piani)
  for(int mu=0;mu<3;mu++){
    for(int nu=mu+1;nu<4;nu++){
      // sommo i 6 risultati in tempo
      tempo  += calc_loc_plaquettes_removing_stag_phases_nnptrick(tconf_acc,local_plaqs,tr_local_plaqs,mu,nu);
    }
  }
  // rimetto le fasi staggered
  mult_conf_times_stag_phases(tconf_acc);

  return tempo;

  }

#endif
