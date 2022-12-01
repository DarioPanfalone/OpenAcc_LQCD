#ifndef FIELD_CORR_C_
#define FIELD_CORR_C_


#include "./field_corr.h"



void calc_field_corr(
    __restrict const su3_soa * const u,
		__restrict su3_soa * const field_corr,
		__restrict su3_soa * const field_corr_aux,
		__restrict su3_soa * const loc_plaq,	
		__restrict d_complex * const trace,
		__restrict single_su3 * const closed_corr,
        const int mu, const int nu, const int ro)
{

  //calcolo le placchette nel piano mu, nu in ogni sito
#pragma acc kernels present(u) present(field_corr) present(loc_plaq) present(nnp_openacc)
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
	  
	  dir_muA = 2*nu +  parity; 
	  dir_muC = 2*nu + !parity;  
	  idxpnu = nnp_openacc[idxh][mu][parity];// r+nu
	    
	  dir_nuB = 2*mu + !parity;
	  dir_nuD = 2*mu +  parity;
	  idxpmu = nnp_openacc[idxh][nu][parity];// r+mu
	  //      r+nu (B) r+mu+nu
	  //          +--->+
	  // nu       ^    |
	  // ^    (A) |    v (C)
	  // |        +<---+
	  // |       r (D) r+mu
	  // +---> mu

		mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpnu,&field_corr[parity],idxh);   // field_corr = A * B
		mat1_times_conj_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_muC],idxpmu);              // field_corr = field_corr * C 
		mat1_times_conj_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_nuD],idxh);                // field_corr = field_corr * D

		//copia delle placchette che serve per poi chiudere i correlatori 
	  assign_su3_soa_to_su3_soa_component_nc(&field_corr[parity], &loc_plaq[parity], idxh); 
	 
	    }  // d0
      }  // d1
    }  // d2
  }  // d3
  
  //calcolo dei correlatori al variare della lunghezza 
  int idxh, parity, idxpro,idxmro, dir_roE;    

 for(int L=1; L<=nd0/2; L++){
	//trace[L]=0;	
	//d3 tempo
#pragma acc kernels present(u) present(field_corr) present(loc_plaq) present(field_corr_aux) present(closed_corr) present(trace) present(nnp_openacc)
#pragma acc loop independent gang(STAPGANG3)
  	for(int d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
	  	for(int d2=0 ; d2<nd2; d2++) {
      		for(int d1=0; d1<nd1; d1++) {
				     for(int d0=0; d0 < nd0; d0++) {
	
	idxh = snum_acc(d0,d1,d2,d3);  // r  
	parity = (d0+d1+d2+d3) % 2; 
	 
	  // moltiplico fieldcorr in ogni sito r per la linea di Schwinger lungo ro	
	  //
	  //                            r+nu   (B)   r+nu+mu
	  // 	nu	                          +------>+
	  //      ^                         ^       |
	  //      |                     (A) |	      | (C) 	
	  //      |                         |  (D)  v   	
	  //      +---->mu                r + <-----+ r+mu   
	  //     /                          ^|
	  //    v                      (E*) |v (E)
	  //    ro	 	                	+
	  //                     		 r+ro 
	     
	dir_roE = 2*ro + !parity; 
	idxpro = nnp_openacc[idxh][ro][parity];  // r+ro
	//idxmro = nnm_openacc[idxh][ro][parity];  // r-ro                       	  

	// FxU (F=fieldcorr, U=link E)     
	mat1_times_mat2_into_mat3_absent_stag_phases_nc(&field_corr[parity], idxh, &u[dir_roE], idxh, &field_corr[parity], idxh);
	// U*FxU 
	conj_mat1_times_mat2_into_mat2_absent_stag_phases_nc(&u[dir_roE],idxh,&field_corr[parity],idxh);
  // chiudo moltiplicando per il complesso coniugato della placchetta meno la placchetta nel sito r+ro: U*FxUx(G*-G-1/3trace(G*-G))  [G=placchetta]
	mat1_times_conj_mat2_minus_mat2_into_single_mat3_absent_stag_phases_nc(&field_corr[parity], idxh, &loc_plaq[!parity], idxpro, closed_corr);
		 
	//assign_su3_soa_to_su3_soa_diff_idx_component(&field_corr[parity], idxh, &field_corr2[!parity], idxmro);	
	//copia fieldcorr nel sito r in fieldcorr2
	assign_su3_soa_to_su3_soa_component_nc(&field_corr[parity], &field_corr_aux[parity], idxh);

	
	trace[L] = trace[L] + single_matrix_trace_absent_stag_phase(closed_corr);
 	
 			 }  // d0
      }  // d1
    }  // d2
  }  // d3

 // scambio fieldcorr nel sito r con fieldcorr2 
 //d3 tempo
#pragma acc kernels present(field_corr) present(field_corr_aux)
#pragma acc loop independent gang(STAPGANG3)
 for(int d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)	
 	  for(int d2=0; d2<nd2; d2++) {
     		for(int d1=0; d1<nd1; d1++) {
						for(int d0=0; d0 < nd0; d0++) {
	
	idxh = snum_acc(d0,d1,d2,d3);  // r  
	parity = (d0+d1+d2+d3) % 2; 
	//idxmro = nnm_openacc[idxh][ro][parity];// r-ro 
	
	//assign_su3_soa_to_su3_soa_diff_idx_component(&field_corr2[!parity], idxmro, &field_corr[parity], idxh);	
	assign_su3_soa_to_su3_soa_component_nc(&field_corr_aux[parity], &field_corr[parity],idxh);
	
		    	}  // d0
        }  // d1
      }  // d2
    }  // d3


  
 } //l

}// closes routine 

/*





void calc_field_corr(
        __restrict const su3_soa * const u,
		__restrict su3_soa * const field_corr,
		__restrict d_complex * const traccia,
        const int mu, const int nu, const int ro)
{

  int d0, d1, d2, d3;
  su3_soa *loc_plaq;
  //calcolo le placchette nel piano mu, nu in ogni sito
#pragma acc kernels present(u) present(field_corr)
#pragma acc loop independent gang(STAPGANG3)
 //d3=tempo
  for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
          for(d0=0; d0 < nd0; d0++) {
	  int idxh,idxpmu,idxpnu;
	  int parity;
	  int dir_muA,dir_nuB;
	  int dir_muC,dir_nuD;
	  idxh = snum_acc(d0,d1,d2,d3);  // r  
	  parity = (d0+d1+d2+d3) % 2;  
	  
	  dir_muA = 2*nu +  parity; 
	  dir_muC = 2*nu + !parity;  
	  idxpnu = nnp_openacc[idxh][mu][parity];// r+nu
	    
	  dir_nuB = 2*mu + !parity;
	  dir_nuD = 2*mu +  parity;
	  idxpmu = nnp_openacc[idxh][nu][parity];// r+mu
	  //      r+nu (B) r+mu+nu
	  //          +--->+
	  // nu       ^    |
	  // ^    (A) |    v (C)
	  // |        +<---+
	  // |       r (D) r+mu
	  // +---> mu

	  mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpnu,&field_corr[parity],idxh);   // field_corr = A * B
	  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_muC],idxpmu);              // field_corr = field_corr * C 
	  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_nuD],idxh);                // field_corr = field_corr * D

      //copia delle placchette che serve per poi chiudere i correlatori 
	  assign_su3_soa_to_su3_soa_component(&field_corr[parity], &loc_plaq[parity]_, idx) 
	 
	    }  // d0
      }  // d1
    }  // d2
  }  // d3
  
  //calcolo dei correlatori al variare della lunghezza 
  int idxh, parity, idxpro,idxmro, dir_roE;
  su3_soa *loc_plaq;
  su3_soa *field_corr2;
  single_su3 *closed_corr;
  traccia[nd0/2];
  
  for(int L=1; L<=nd0/2; L++){
	traccia[L]=0;	
	//d3 tempo
	for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
		for(d2=0; d2<nd2; d2++) {
      		for(d1=0; d1<nd1; d1++) {
				for(d0=0; d0 < nd0; d0++) {
	
	idxh = snum_acc(d0,d1,d2,d3);  // r  
	parity = (d0+d1+d2+d3) % 2; 
	 
	  // moltiplico fieldcorr in ogni sito r per la linea di Schwinger lungo ro	
	  //
	  //                            r+nu   (B)   r+nu+mu
	  // 	nu	                        +------>+
	  //      ^                         ^       |
	  //      |                     (A) |	    | (C) 	
	  //      |                         |  (D)  v   	
	  //      +---->mu                r + <-----+ r+mu   
	  //     /                          ^|
	  //    v                      (E*) |v (E)
	  //    ro	 	                	+
	  //                     		 r+ro 
	     
	dir_roE = 2*ro + !parity; 
	idxpro = nnp_openacc[idxh][ro][parity];  // r+ro
	idxmro = nnm_openacc[idxh][ro][parity];  // r-ro                       	  

	// FxU (F=fieldcorr, U=link E)     
	mat1_times_mat2_into_mat3_absent_stag_phases(&field_corr[parity], idxh,&u[dir_roE], idxh, &field_corr[parity], idxh);         
	// U*FxU 
	conj_mat1_times_mat2_into_mat2_absent_stag_phases(&u[dir_roE],idxh,&field_corr[parity],idxh);	   	    
    // chiudo moltiplicando per il complesso coniugato della placchetta meno la placchetta nel sito r+ro: U*FxUx(G*-G)  [G=placchetta]
	mat1_times_conj_mat2_minus_mat2_into_single_mat3_absent_stag_phases(&field_corr[parity], idxh, &loc_plaq[!parity], idxpro, closed_corr); 
	
	// copia fieldcorr nel sito r in fieldcorr2 in r-ro
	assign_su3_soa_to_su3_soa_diff_idx_component(&field_corr[parity], idxh, &field_corr2[!parity], idxmro);	
	
	traccia[L] = traccia[L] + single_matrix_trace_absent_stag_phase(&closed_corr);
	
			}  // d0
      }  // d1
    }  // d2
  }  // d3

 // scambio fieldcorr nel sito r con fieldcorr2 
 //d3 tempo
	for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
		for(d2=0; d2<nd2; d2++) {
      		for(d1=0; d1<nd1; d1++) {
				for(d0=0; d0 < nd0; d0++) {
	
	idxh = snum_acc(d0,d1,d2,d3);  // r  
	parity = (d0+d1+d2+d3) % 2; 
	idxmro = nnm_openacc[idxh][ro][parity];  // r-ro 
	
	assign_su3_soa_to_su3_soa_diff_idx_component(&field_corr2[!parity], idxmro, &field_corr[parity], idxh);		
	
			}  // d0
      }  // d1
    }  // d2
  }  // d3


  
} //l

}// closes routine 



























void calc_field_corr(
        __restrict const su3_soa * const u,
		__restrict su3_soa * const field_corr,
		__restrict d_complex * const traccia,
        const int mu, const int nu, const int ro)
{

  int d0, d1, d2, d3;
  //calcolo le placchette nel piano mu, nu in ogni sito
#pragma acc kernels present(u) present(field_corr)
#pragma acc loop independent gang(STAPGANG3)
 //d3=tempo
  for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
          for(d0=0; d0 < nd0; d0++) {
	  int idxh,idxpmu,idxpnu;
	  int parity;
	  int dir_muA,dir_nuB;
	  int dir_muC,dir_nuD;
	  idxh = snum_acc(d0,d1,d2,d3);  // r  
	  parity = (d0+d1+d2+d3) % 2; //  
	  
	  dir_muA = 2*nu +  parity; 
	  dir_muC = 2*nu + !parity; // 
	  idxpnu = nnp_openacc[idxh][mu][parity];// r+nu
	    
	  dir_nuB = 2*mu + !parity;
	  dir_nuD = 2*mu +  parity;
	  idxpmu = nnp_openacc[idxh][nu][parity];// r+mu
	  //      r+nu (B) r+mu+nu
	  //          +--->+
	  // nu       ^    |
	  // ^    (A) |    v (C)
	  // |        +<---+
	  // |       r (D) r+mu
	  // +---> mu

	  mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpnu,&field_corr[parity],idxh);   // LOC_PLAQ = A * B
	  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_muC],idxpmu);              // LOC_PLAQ = LOC_PLAQ * C 
	  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_nuD],idxh);                // LOC_PLAQ = LOC_PLAQ * D
	
	  
	    }  // d0
      }  // d1
    }  // d2
  }  // d3
  
  //calcolo dei correlatori al variare della lunghezza 
  int idxh, parity, idxpro,idxmro, dir_roE;
  su3_soa *loc_plaq;
  set_su3_soa_to_su3_soa(&field_corr, &loc_plaq);
  su3_soa *field_corr2;
  su3_soa *field_aux;
  single_su3 *closed_corr;
  traccia[nd0/2];
  
  for(int L=1; L<=nd0/2; L++){
	traccia[L]=0;	
	//d3 tempo
	for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
		for(d2=0; d2<nd2; d2++) {
      		for(d1=0; d1<nd1; d1++) {
				for(d0=0; d0 < nd0; d0++) {
	
	idxh = snum_acc(d0,d1,d2,d3);  // r  
	parity = (d0+d1+d2+d3) % 2; 
	 
	  // moltiplico la placchetta in ogni sito r per la linea di Schwinger lungo ro	
	  //
	  //                            r+nu   (B)   r+nu+mu
	  // 	nu	                        +------>+
	  //      ^                         ^       |
	  //      |                     (A) |	    | (C) 	
	  //      |                         |  (D)  v   	
	  //      +---->mu                r + <-----+ r+mu   
	  //     /                          ^|
	  //    v                      (E*) |v (E)
	  //    ro	 	                	+
	  //                     		 r+ro 
	     
	dir_roE = 2*ro + !parity; 
	idxpro = nnp_openacc[idxh][ro][parity];  // r+ro
	idxmro = nnm_openacc[idxh][ro][parity];  // r-ro                       	  
	su3_soa * conf_to_use;
	su3_soa * conf_to_copy;
	su3_soa * conf_aux;
	conf_to_use = &field_corr[parity];
	conf_to_copy = &field_corr2[!parity];
	conf_aux = &field_aux; //?	

	//GxU (G=placchetta, U=link E)     
	mat1_times_mat2_into_mat3_absent_stag_phases(&field_corr[parity], idxh,&u[dir_roE], idxh, &field_corr[parity], idxh);         
	// U*GxU 
	conj_mat1_times_mat2_into_mat2_absent_stag_phases(&u[dir_roE],idxh,&field_corr[parity],idxh);
		   	    
    //chiudo moltiplicando per il complesso coniugato della placchetta meno la placchetta nel sito r+ro: U*GxUx(G*-G)
	mat1_times_conj_mat2_minus_mat2_into_single_mat3_absent_stag_phases(&field_corr[parity], idxh, &loc_plaq[!parity], idxpro, closed_corr); 
	
	//copia fieldcorr nel sito r in fieldcorr2 in r-ro
	assign_su3_soa_to_su3_soa_diff_idx_component(&field_corr[parity], idxh, &field_corr2[!parity], idxmro);
	//scambio fieldcorr2 e fieldcorr
	*conf_aux = *conf_to_use; 
	*conf_to_use = *conf_to_copy;
	*conf_to_copy = *conf_aux;
	
	traccia[L] = traccia[L] + single_matrix_trace_absent_stag_phase(&closed_corr);
	
			}  // d0
      }  // d1
    }  // d2
  }  // d3



} //l


}// closes routine corr















double calc_loc_plaquettes_nnptrick(
        __restrict const su3_soa * const u,
        __restrict su3_soa * const field_corr,
        dcomplex_soa * const tr_field_corr,
        const int mu, const int nu, const int ro)
{

  int d0, d1, d2, d3, l;
#pragma acc kernels present(u) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent gang(STAPGANG3)
 //d3=tempo
  for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
          for(d0=0; d0 < nd0; d0++) {
	  int i,idxh,idxpmu,idxpimu, idxpnu,idxpro,idxplro,idxplnu;
	  int idxhl,idxhi;
	  int parity,parityl,parityi;
	  int dir_nuA, dir_roB, dir_nuC, dir_roD, dir_muE, dir_muF, dir_roG, dir_nuH, dir_roI, dir_nuL;

	  idxh = snum_acc(d0,d1,d2,d3);  // r  
	  parity = (d0+d1+d2+d3) % 2; 
	
	  dir_nuA = 2*nu +  parity; 
	  dir_nuC = 2*nu + !parity;  
	  idxpnu = nnp_openacc[idxh][nu][parity]; // r+nu
	
	  dir_roB = 2*ro + !parity;
	  dir_roD = 2*ro +  parity;
	  idxpro = nnp_openacc[idxh][ro][parity]; // r+ro
                    
	  //                     
	  //				             r+nu  (B)   r+nu+ro
	  // 		                        +------->+
	  //                                ^        |
	  //                            (A) |		 | (C) 	
	  //                                |  (D)   v   	
	  //                              r + <------+ r+ro   
	  //                                ^|
	  //     nu                     (E) || (F)
	  //     ^		                    |v 	
	  //     |                r+l(mu)   + ------>+ r+l(mu)+ro
	  //     |                     	    ^  (G)   |
	  //     +----> mu              (L) |        | (H)
	  //    /                           |        v
 	  //   v 				            +<-------+
	  //    ro		          r+l(mu)+nu   (I)   r+l(mu)+nu+ro
	  //

	  mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_nuA],idxh,&u[dir_roB],idxpnu,&field_corr[parity],idxh);  // field_corr = A * B
	  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_nuC],idxpro);              // field_corr = field_corr * C 
	  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_roD],idxh);                // field_corr = field_corr * D
	 
	  for(l=d0+1 ;l<nd0 ;l++){
     		
			idxhl = snum_acc(l,d1,d2,d3); // r+l(mu)
			parityl = (l+d1+d2+d3) % 2;
			
			for(i=d0+1; i<l; i++){
					 	
						idxhi = snum_acc(i,d1,d2,d3);	//r+i
						parityi = (i+d1+d2+d3) % 2;

						dir_muF = 2*mu +  parityi;  
	  					dir_muE = 2*mu + !parityi; 
						idxpimu = nnp_openacc[idxhi][mu][parityi]; // r+i+mu

			mat1_times_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_muF],idxh);              // field_corr = field_corr * F
	  		mat1_times_conj_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_muE],idxpimu);      // field_corr = field_corr * E
			}
			
		
		dir_roG = 2*ro +  parityl; 
	  	dir_roI = 2*ro + !parityl;  
	  	idxplro = nnp_openacc[idxhl][ro][parityl]; // r+l(mu)+ro
	
	 	dir_nuH = 2*nu + !parityl;
	  	dir_nuL = 2*nu +  parityl;
	 	idxplnu = nnp_openacc[idxhl][nu][parityl]; // r+l(mu)+nu
	
	  	mat1_times_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_roG],idxhl);  				    // field_corr = field_corr * G
		mat1_times_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_nuH],idxplro);  			    // field_corr = field_corr * H
		mat1_times_conj_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_roI],idxplnu);            // field_corr = field_corr * I
		mat1_times_conj_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_nuC],idxhl);              // field_corr = field_corr * L

		d_complex traccia = matrix_trace_absent_stag_phase(&field_corr[parity],idxh);
		tr_field_corr[parity].c[idxh] = creal(traccia)+cimag(traccia)*I;
	  }
	   
	  }  // d0
      }  // d1
    }  // d2
  }  // d3


  }
// closes routine


                                               //////////////////////////////////////


double calc_loc_plaquettes_nnptrick(
        __restrict const su3_soa * const u,
        __restrict su3_soa * const field_corr,
        dcomplex_soa * const tr_field_corr,
        const int mu, const int nu, const int ro)
{

  int d0, d1, d2, d3, l;
#pragma acc kernels present(u) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent gang(STAPGANG3)
 //d3=tempo
  for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
          for(d0=0; d0 < nd0; d0++) {
	  int i,idxh,idxpmu,idxpimu, idxpnu,idxpro,idxplro,idxplnu;
	  int idxhl,idxhi;
	  int parity,parityl,parityi;
	  int dir_nuA, dir_roB, dir_nuC, dir_roD, dir_muE, dir_muF, dir_roG, dir_nuH, dir_roI, dir_nuL;

	  idxh = snum_acc(d0,d1,d2,d3);  // r  
	  parity = (d0+d1+d2+d3) % 2; 
	
	  dir_nuA = 2*nu +  parity; 
	  dir_nuC = 2*nu + !parity;  
	  idxpnu = nnp_openacc[idxh][nu][parity]; // r+nu
	
	  dir_roB = 2*ro + !parity;
	  dir_roD = 2*ro +  parity;
	  idxpro = nnp_openacc[idxh][ro][parity]; // r+ro
                    
	  //                     
	  //				             r+nu  (B)   r+nu+ro
	  // 		                        +------->+
	  //                                ^        |
	  //                            (A) |		 | (C) 	
	  //                                |  (D)   v   	
	  //                              r + <------+ r+ro   
	  //                                ^|
	  //     nu                     (E) || (F)
	  //     ^		                    |v 	
	  //     |                r+l(mu)   + ------>+ r+l(mu)+ro
	  //     |                     	    ^  (G)   |
	  //     +----> mu              (L) |        | (H)
	  //    /                           |        v
 	  //   v 				            +<-------+
	  //    ro		          r+l(mu)+nu   (I)   r+l(mu)+nu+ro
	  //

	  mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_nuA],idxh,&u[dir_roB],idxpnu,&field_corr[parity],idxh);  // field_corr = A * B
	  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_nuC],idxpro);              // field_corr = field_corr * C 
	  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_roD],idxh);                // field_corr = field_corr * D
	 idxhi=idxh;
	  for(l=d0+1 ;l<nd0 ;l++){
     		//calcolo la parità del sito r+l*mu
			if(l%2==1){				
				parityl=!parity;
			}
			 else {
				parityl=parity;
				}
			
			for(i=0; i<l; i++){
					if(i%2==1){
						idxhi = nnp_openacc[idxhi][mu][!parity]; //r+i*mu
						dir_muF = 2*mu +  parity;  
	  				  	dir_muE = 2*mu + !parity;
			 mat1_times_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_muF],idxh);              // field_corr = field_corr * F
	  		 mat1_times_conj_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_muE],idxhi);      // field_corr = field_corr * E
					}
				else{
					idxhi = nnp_openacc[idxhi][mu][parity]; //r+i*mu
					dir_muF = 2*mu + !parity;  
	  			  	dir_muE = 2*mu +  parity;
			 mat1_times_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_muF],idxh);              // field_corr = field_corr * F
	  		 mat1_times_conj_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_muE],idxhi);      // field_corr = field_corr * E
				}  
			} //alla fine sono nel sito idxhi=r+l*mu
						
		dir_roG = 2*ro +  parityl; 
	  	dir_roI = 2*ro + !parityl;  
	  	idxplro = nnp_openacc[idxhl][ro][parityl]; // r+l*mu+ro
	
	 	dir_nuH = 2*nu + !parityl;
	  	dir_nuL = 2*nu +  parityl;
	 	idxplnu = nnp_openacc[idxhl][nu][parityl]; // r+l*mu+nu
	
	  	mat1_times_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_roG],idxhl);  				    // field_corr = field_corr * G
		mat1_times_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_nuH],idxplro);  			    // field_corr = field_corr * H
		mat1_times_conj_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_roI],idxplnu);            // field_corr = field_corr * I
		mat1_times_conj_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_nuC],idxhl);              // field_corr = field_corr * L

		d_complex traccia = matrix_trace_absent_stag_phase(&field_corr[parity],idxh);
		tr_field_corr[parity].c[idxh] = creal(traccia)+cimag(traccia)*I;
	  }
	   
	  }  // d0
      }  // d1
    }  // d2
  }  // d3


  }
// closes routine



                                  //////////////////////////////

double calc_loc_plaquettes_nnptrick(
        __restrict const su3_soa * const u,
        __restrict su3_soa * const field_corr,
        dcomplex_soa * const tr_field_corr,
        const int mu, const int nu)
{

  int d0, d1, d2, d3;
#pragma acc kernels present(u) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent gang(STAPGANG3)
 //d3=tempo
  for(d3=D3_HALO; d3<nd3-D3_HALO; d3++) {
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
    for(d2=0; d2<nd2; d2++) {
      for(d1=0; d1<nd1; d1++) {
          for(d0=0; d0 < nd0; d0++) {
	  int idxh,idxpmu,idxpnu;
	  int parity;
	  int dir_muA,dir_nuB;
	  int dir_muC,dir_nuD;
	  int dir_muE,dir_muF;
	  idxh = snum_acc(d0,d1,d2,d3);  // r  
	  parity = (d0+d1+d2+d3) % 2; 
	
	  dir_muA = 2*nu +  parity; 
	  dir_muC = 2*nu + !parity; // 
	  idxpmu = nnp_openacc[idxh][nu][parity];// r+mu
	    
	  dir_nuB = 2*mu + !parity;
	  dir_nuD = 2*mu +  parity;
	  idxpnu = nnp_openacc[idxh][mu][parity];// r+nu
	  //                    
	  //                     
	  //				  r+nu  (B)  r+nu+mu
	  // 		             +------->+
	  //                     ^        |
	  //                 (A) |		  | (C) 	
	  //                     |  (D)   v   	
	  //                   r + <------+ r+mu 
	  //		        	    (E)     
	  //     nu              +------->+
	  //     ^		         +<-------+
	  //     |                  (F) 
	  //     |          	   
	  //     +---> mu
	    

	  mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&field_corr[parity],idxh);  // field_corr = A * B
	  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_muC],idxpmu);              // field_corr = field_corr * C 
	  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_nuD],idxh);                // field_corr = field_corr * D
	  
	  mat1_times_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_nuB],idxpmu);               	   // field_corr = field_corr * E
	  mat1_times_conj_mat2_into_mat1_absent_stag_phases(&field_corr[parity],idxh,&u[dir_nuD],idxh);                // field_corr = field_corr * F
	  
	
	  

	}  // d0
      }  // d1
    }  // d2
  }  // d3

  double res_R_p = 0.0;
  double res_I_p = 0.0;
  double resR = 0.0;
  int t;  // ONLY GOOD FOR 1D CUT
	#pragma acc kernels present(tr_local_plaqs)
	#pragma acc loop reduction(+:res_R_p) reduction(+:res_I_p)
		for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
			res_R_p += creal(tr_local_plaqs[0].c[t]);
			res_R_p += creal(tr_local_plaqs[1].c[t]);
  }

  //printf("res_R_p %e , mu %d  nu %d\n", res_R_p, mu ,nu);
  return res_R_p;
}// closes routine
*/
#endif
