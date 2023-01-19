#ifndef FIELD_CORR_C_
#define FIELD_CORR_C_


#include "./field_corr.h"



void calc_field_corr(
										 __restrict su3_soa * const u,
										 __restrict su3_soa * const field_corr,
										 __restrict su3_soa * const field_corr_aux,
										 __restrict su3_soa * const loc_plaq,	
										 dcomplex_soa * const trace_local,
										 d_complex * const corr,
										 __restrict single_su3 * const closed_corr,
										 const int mu, const int nu, const int ro)
{
	
  //calcolo le placchette nel piano mu, nu in ogni sito
#pragma acc kernels present(u) present(field_corr) present(loc_plaq) 
#pragma acc loop independent gang(STAPGANG3)
	//d3=tempo
  for(int d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(int d2=0; d2<nd2; d2++) {
      for(int d1=0; d1<nd1; d1++) {
				for(int d0=0; d0 < nd0; d0++) {
	  int idxh,idxpmu,idxpnu;
	  int parity;
	  int dir_muA,dir_nuB;
	  int dir_muC,dir_nuD;
	  idxh = snum_acc(d0,d1,d2,d3);  // r  
	  parity = (d0+d1+d2+d3) % 2;  
	  
	  dir_muA = 2*mu +  parity; 
	  dir_muC = 2*mu + !parity;  
	  idxpmu = nnp_openacc[idxh][mu][parity];// r+mu
	    
	  dir_nuB = 2*nu + !parity;
	  dir_nuD = 2*nu +  parity;
	  idxpnu = nnp_openacc[idxh][nu][parity];// r+nu

		//      r+nu (C) r+mu+nu
	  //          +<---+
	  // nu       |    ^
	  // ^    (D) v    | (B)
	  // |        +--->+
	  // |       r (A) r+mu
	  // +---> mu

		mat1_times_mat2_into_mat3_absent_stag_phases_nc(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&field_corr[parity],idxh);   // field_corr = A * B
		mat1_times_conj_mat2_into_mat1_absent_stag_phases_nc(&field_corr[parity],idxh,&u[dir_muC],idxpnu);              // field_corr = field_corr * C 
		mat1_times_conj_mat2_into_mat1_absent_stag_phases_nc(&field_corr[parity],idxh,&u[dir_nuD],idxh);                // field_corr = field_corr * D
		
		//copia delle placchette che serve per poi chiudere i correlatori 
	  assign_su3_soa_to_su3_soa_component_nc(&field_corr[parity], &loc_plaq[parity], idxh); 
	 
				}  // d0
      }  // d1
    }  // d2
  }  // d3
  
  //calcolo dei correlatori al variare della lunghezza 

	for(int L=0; L<nd0/2; L++){

	//d3 tempo
#pragma acc kernels present(u) present(field_corr) present(loc_plaq) present(field_corr_aux) present(closed_corr) present(trace_local)
#pragma acc loop independent gang(STAPGANG3)
		for(int d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
	  	for(int d2=0 ; d2<nd2; d2++) {
				for(int d1=0; d1<nd1; d1++) {
						for(int d0=0; d0 < nd0; d0++) {
							
							int idxh = snum_acc(d0,d1,d2,d3);  // r  
							int parity = (d0+d1+d2+d3) % 2; 
	 
	  // moltiplico fieldcorr in ogni sito r per la linea di Schwinger lungo ro	
	  //
	  //                            r+nu       r+nu+mu
	  // 	nu	                          +<------+
	  //      ^                         |       ^
	  //      |                         |	      |  	
	  //      |                         v       |   	
	  //      +---->mu                r + ----->+ r+mu   
	  //     /                          |^
	  //    v             (E(^dagger))  v| (E)
	  //    ro	 	                	    +
	  //                     		       r+ro 
	     
							int dir_roE = 2*ro + parity; 
							int idxpro = nnp_openacc[idxh][ro][parity];  // r+ro
							//int idxmro = nnm_openacc[idxh][ro][parity];  // r-ro

// FxU (F=fieldcorr, U=link E)     
							mat1_times_mat2_into_mat1_absent_stag_phases_nc(&field_corr[parity], idxh, &u[dir_roE], idxh);
// U^(dagger)xFxU
							conj_mat1_times_mat2_into_mat2_absent_stag_phases_nc(&u[dir_roE],idxh,&field_corr[parity],idxh);

//Per fare la prova con la configuazione di identità: UxFxU^(dagger)xG^(dagger)
//						mat1_times_conj_mat2_into_single_mat3_absent_stag_phases_nc(&field_corr[parity], idxh, &loc_plaq[!parity], idxpro, closed_corr);

//copia field_corr in field_corr_aux 
							assign_su3_soa_to_su3_soa_diff_idx_component(&field_corr[parity], idxh, &field_corr_aux[parity], idxh);

//chiudo moltiplicando per il complesso coniugato della placchetta meno la placchetta nel sito r+ro: UxFxU^(dagger)x(G^(dagger)-G-1/3trace(G^(dagger)-G)) [G=placchetta]
				 			mat1_times_conj_mat2_minus_mat2_into_single_mat3_absent_stag_phases_nc(&field_corr[parity], idxh, &loc_plaq[!parity], idxpro, closed_corr);

						  d_complex tmp_aux =  single_matrix_trace_absent_stag_phase(closed_corr);
							trace_local[parity].c[idxh] = creal(tmp_aux) + cimag(tmp_aux)*I ;
 	
						}  // d0
				}  // d1
			}  // d2
		}  // d3

		double res_R_p = 0.0;
		double res_I_p = 0.0;
		double resR = 0.0;
		int t; 
		
#pragma acc kernels present(trace_local)
#pragma acc loop reduction(+:res_R_p) reduction(+:res_I_p)
		for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
			res_R_p += creal(trace_local[0].c[t]); // even sites correlator
			res_R_p += creal(trace_local[1].c[t]); // odd sites correlator
		}
		//	printf("\n%d\t%d\t%d\t%d\t%lf", mu, nu, ro, L, res_R_p);
		corr[L] = res_R_p;
				
		//traslazione 
//d3 tempo
#pragma acc kernels present(field_corr) present(field_corr_aux)
#pragma acc loop independent gang(STAPGANG3)
		for(int d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)	
			for(int d2=0; d2<nd2; d2++) {
				for(int d1=0; d1<nd1; d1++) {
					for(int d0=0; d0 < nd0; d0++) {
						
						int idxh = snum_acc(d0,d1,d2,d3);  // r  
						int parity = (d0+d1+d2+d3) % 2; 
						//int idxmro = nnm_openacc[idxh][ro][parity];// r-ro
						int idxpro = nnp_openacc[idxh][ro][parity];  // r+ro
						
						assign_su3_soa_to_su3_soa_diff_idx_component(&field_corr_aux[!parity], idxpro, &field_corr[parity], idxh);	
	
		    	}  // d0
        }  // d1
      }  // d2
    }  // d3
	  

	}//l

}// closes routine 


void random_gauge_transformation(
																 __restrict su3_soa * const u,
																 __restrict su3_soa * const u_new,
																 __restrict su3_soa * const m_soa)
{
	/*
//Genero per ogni sito del reticolo un elemento G di SU(3) random
#pragma acc kernels present(m_soa)
#pragma acc loop independent gang(STAPGANG3)
	for(int d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
		for(int d2=0; d2<nd2; d2++) {
			for(int d1=0; d1<nd1; d1++) {
				for(int d0=0; d0 < nd0; d0++) {
					
					int idxh = snum_acc(d0,d1,d2,d3);  // r
					int parity = (d0+d1+d2+d3) % 2;

					int factor=2.500000e-01;
					single_su3 aux;
					generate_random_su3(&aux, factor);
					single_gl3_into_su3_soa(&m_soa[parity], idxh, &aux);
					
				}  // d0
			}  // d1
		}  // d2
	}  // d3
*/
	//Faccio la trasformazione di gauge: G(n)xU(n)xG*(n+mu) dove U è il link n-->n+mu
#pragma acc kernels present(m_soa) present(u) present(nnp_openacc)	
#pragma acc loop independent gang(STAPGANG3)
		for(int d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
			for(int d2=0; d2<nd2; d2++) {
				for(int d1=0; d1<nd1; d1++) {
					for(int d0=0; d0 < nd0; d0++) {
					
						int idxh = snum_acc(d0,d1,d2,d3);  // r
						int parity = (d0+d1+d2+d3) % 2;
						int dir_link, idxpmu;
						
						for(int mu=0; mu<4; mu++){
					
						dir_link = 2*mu + parity;
						idxpmu = nnp_openacc[idxh][mu][parity]; // r+mu
						//U_(new) = GxU
						mat1_times_mat2_into_mat3_absent_stag_phases_nc(&m_soa[parity], idxh, &u[dir_link], idxh, &u_new[dir_link], idxh);
						//U_(new) = U_(new)xG*
						mat1_times_conj_mat2_into_mat1_absent_stag_phases_nc(&u_new[dir_link], idxh, &m_soa[!parity], idxpmu);

						}
					
					}  // d0
				}  // d1
			}  // d2
		}  // d3
} 





#endif
