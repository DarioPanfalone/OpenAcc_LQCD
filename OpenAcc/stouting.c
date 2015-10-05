#ifndef STOUTING_C
#define STOUTING_C

#include "./cayley_hamilton.c"

#pragma acc routine seq
static inline void su3_soa_to_single_su3(__restrict su3_soa * const in,
					int idx,
					single_su3 * out){
  out->comp[0][0] = in->r0.c0[idx];
  out->comp[0][1] = in->r0.c1[idx];
  out->comp[0][2] = in->r0.c2[idx];
  out->comp[1][0] = in->r1.c0[idx];
  out->comp[1][1] = in->r1.c1[idx];
  out->comp[1][2] = in->r1.c2[idx];
}

#pragma acc routine seq
static inline void tamat_soa_to_single_tamat( __restrict tamat_soa * const in,
					      int idx,
					      single_tamat * out){
  out->c01  = in->c01[idx];
  out->c02  = in->c02[idx];
  out->c12  = in->c12[idx];
  out->rc00 = in->rc00[idx];
  out->rc11 = in->rc11[idx];

}


static inline void conf_left_exp_multiply_to_su3_soa(__restrict su3_soa * const cnf,
						     const int idx,
						     __restrict su3_soa * const  EXP,
						     __restrict su3_soa * const cnf_out){

  single_su3 AUX;
  //Multiply: U_new = EXP * U_old
  AUX.comp[0][0] = cnf->r0.c0[idx];
  AUX.comp[0][1] = cnf->r0.c1[idx];
  AUX.comp[0][2] = cnf->r0.c2[idx];
  AUX.comp[1][0] = cnf->r1.c0[idx];
  AUX.comp[1][1] = cnf->r1.c1[idx];
  AUX.comp[1][2] = cnf->r1.c2[idx];
  // ricostruisco la terza
  AUX.comp[2][0] = conj(AUX.comp[0][1] * AUX.comp[1][2] - AUX.comp[0][2] * AUX.comp[1][1]);
  AUX.comp[2][1] = conj(AUX.comp[0][2] * AUX.comp[1][0] - AUX.comp[0][0] * AUX.comp[1][2]);
  AUX.comp[2][2] = conj(AUX.comp[0][0] * AUX.comp[1][1] - AUX.comp[0][1] * AUX.comp[1][0]);

  // moltiplica
  cnf_out->r0.c0[idx] = EXP->r0.c0[idx] * AUX.comp[0][0] + EXP->r0.c1[idx] * AUX.comp[1][0] + EXP->r0.c2[idx] * AUX.comp[2][0];
  cnf_out->r0.c1[idx] = EXP->r0.c0[idx] * AUX.comp[0][1] + EXP->r0.c1[idx] * AUX.comp[1][1] + EXP->r0.c2[idx] * AUX.comp[2][1];
  cnf_out->r0.c2[idx] = EXP->r0.c0[idx] * AUX.comp[0][2] + EXP->r0.c1[idx] * AUX.comp[1][2] + EXP->r0.c2[idx] * AUX.comp[2][2];

  cnf_out->r1.c0[idx] = EXP->r1.c0[idx] * AUX.comp[0][0] + EXP->r1.c1[idx] * AUX.comp[1][0] + EXP->r1.c2[idx] * AUX.comp[2][0];
  cnf_out->r1.c1[idx] = EXP->r1.c0[idx] * AUX.comp[0][1] + EXP->r1.c1[idx] * AUX.comp[1][1] + EXP->r1.c2[idx] * AUX.comp[2][1];
  cnf_out->r1.c2[idx] = EXP->r1.c0[idx] * AUX.comp[0][2] + EXP->r1.c1[idx] * AUX.comp[1][2] + EXP->r1.c2[idx] * AUX.comp[2][2];

}


void exp_minus_QA_times_conf(__restrict su3_soa * const tu,
			     __restrict tamat_soa * const QA,
			     __restrict su3_soa * const tu_out,
			     __restrict su3_soa * const exp_aux){
  
  int x, y, z, t;
#pragma acc kernels present(tu) present(QA) present(tu_out) present(exp_aux)
#pragma acc loop independent gang
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector
        for(x=0; x < nx; x++) {
	  const  int idxh   = snum_acc(x,y,z,t);  // r
          const  int parity = (x+y+z+t) % 2;
          int dir_link;
          int mu;

          for(mu=0;mu<4;mu++){

            dir_link = 2*mu + parity;
	    //	    single_su3 tempU;
	    //	    single_su3 EXP;
	    //	    single_tamat tempQA;
	    //	    su3_soa_to_single_su3(&tu[dir_link],idxh,&tempU);
	    //	    tamat_soa_to_single_tamat(&QA[dir_link],idxh,&tempQA);

	    CH_exponential_antihermitian_soa(&exp_aux[dir_link],&QA[dir_link],idxh); // qui la variabile tu_out viene usata come output temporaneo
	    //	    CH_exponential_antihermitian(&EXP,&tempQA);
	    conf_left_exp_multiply_to_su3_soa(&tu[dir_link],idxh,&exp_aux[dir_link],&tu_out[dir_link]);
          }

        }  // x
      }  // y
    }  // z
  }  // t

}// closes routine



void stout_isotropic( __restrict su3_soa * const u,
		      __restrict su3_soa * const uprime,
		      __restrict su3_soa * const local_staples,
		      __restrict su3_soa * const auxiliary,
		      __restrict tamat_soa * const tipdot){

  set_su3_soa_to_zero(local_staples);
  mult_conf_times_stag_phases(u);


  calc_loc_staples_removing_stag_phases_nnptrick_all(u,local_staples);

  RHO_times_conf_times_staples_ta_part(u,local_staples,tipdot);

  exp_minus_QA_times_conf(u,tipdot,uprime,auxiliary);

  mult_conf_times_stag_phases(u);
  mult_conf_times_stag_phases(uprime);

}



#endif
