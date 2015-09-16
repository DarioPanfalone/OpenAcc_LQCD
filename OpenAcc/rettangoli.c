#ifndef RETTANGOLI_C
#define RETTANGOLI_C

double calc_loc_rectangles_2x1_removing_stag_phases_nnptrick(   __restrict su3_soa * const u,
								__restrict su3_soa * const loc_plaq,
								dcomplex_soa * const tr_local_plaqs,
								const int mu,
								const int nu){
  
  int x, y, z, t;
#pragma acc kernels present(u) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent gang
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector
        for(x=0; x < nx; x++) {
          int idxh,idxpmu,idxpmupmu,idxpmupnu,idxpnu;
          int parity;
          int dir_muA,dir_muB,dir_nuC,dir_muD,dir_muE,dir_nuF;

          idxh = snum_acc(x,y,z,t);  // r
          parity = (x+y+z+t) % 2;

          dir_muA = 2*mu +  parity;
          dir_muB = 2*mu + !parity;
          dir_nuC = 2*nu +  parity;
          dir_muD = 2*mu +  parity;
          dir_muE = 2*mu + !parity;
          dir_nuF = 2*nu +  parity;
	  idxpmu = nnp_openacc[idxh][mu][parity];// r+mu 
          idxpmupmu = nnp_openacc[idxpmu][mu][!parity];// r+2mu
          idxpmupnu = nnp_openacc[idxpmu][nu][!parity];// r+mu+nu
          idxpnu = nnp_openacc[idxh][nu][parity];// r+nu

          //       r+nu r+mu+nu r+2mu+nu   
          //          +<---+<---+         
          // nu       | (E) (D) ^         
          // ^    (F) V (A) (B) |  (C)     
          // |        +--->+--->+         
          // |       r   r+mu r+2mu       
          // +---> mu    

          mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_muB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
          mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuC],idxpmupmu);                 // LOC_RECT = LOC_RECT * C
          mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muD],idxpmupnu);            // LOC_RECT = LOC_RECT * D 
          mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muE],idxpnu);               // LOC_RECT = LOC_RECT * E 
          mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F
          tr_local_plaqs[parity].c[idxh] = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);

        }  // x          
      }  // y            
    }  // z              
  }  // t                

  double res_R_p = 0.0;
  double res_I_p = 0.0;
  double resR = 0.0;
#pragma acc kernels present(tr_local_plaqs)
#pragma acc loop reduction(+:res_R_p) reduction(+:res_I_p)
  for(t=0; t<sizeh; t++) {
    res_R_p += creal(tr_local_plaqs[0].c[t]);
    res_R_p += creal(tr_local_plaqs[1].c[t]);
    res_I_p += cimag(tr_local_plaqs[0].c[t]);
    res_I_p += cimag(tr_local_plaqs[1].c[t]);
  }

  return res_R_p;
}// closes routine 



double calc_loc_rectangles_1x2_removing_stag_phases_nnptrick(   __restrict su3_soa * const u,
								__restrict su3_soa * const loc_plaq,
								dcomplex_soa * const tr_local_plaqs,
								const int mu,
								const int nu){
  
  int x, y, z, t;
#pragma acc kernels present(u) present(loc_plaq) present(tr_local_plaqs)
#pragma acc loop independent gang
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector
        for(x=0; x < nx; x++) {
          int idxh,idxpmu,idxpnupnu,idxpmupnu,idxpnu;
          int parity;
          int dir_muA,dir_nuB,dir_nuC,dir_muD,dir_nuE,dir_nuF;

          idxh = snum_acc(x,y,z,t);  // r
          parity = (x+y+z+t) % 2;

          dir_muA = 2*mu +  parity;
          dir_nuB = 2*nu + !parity;
          dir_nuC = 2*nu +  parity;
          dir_muD = 2*mu +  parity;
          dir_nuE = 2*nu + !parity;
          dir_nuF = 2*nu +  parity;

	  idxpmu = nnp_openacc[idxh][mu][parity];        // r+mu 
          idxpmupnu = nnp_openacc[idxpmu][nu][!parity];  // r+mu+nu
          idxpnu = nnp_openacc[idxh][nu][parity];        // r+nu
          idxpnupnu = nnp_openacc[idxpnu][nu][!parity];  // r+nu+nu
	  //            (D)
	  //    r+2nu +<---+ r+mu+2nu
	  //          |    ^
          //      (E) V    | (C)
          //     r+nu +    + r+mu+nu
          // nu       |    ^
          // ^    (F) V    | (B)
          // |        +--->+
          // |       r  (A)  r+mu
          // +---> mu    

          mat1_times_mat2_into_mat3_absent_stag_phases(&u[dir_muA],idxh,&u[dir_nuB],idxpmu,&loc_plaq[parity],idxh);   // LOC_RECT = A * B
          mat1_times_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuC],idxpmupnu);                 // LOC_RECT = LOC_RECT * C
          mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_muD],idxpnupnu);            // LOC_RECT = LOC_RECT * D 
          mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuE],idxpnu);               // LOC_RECT = LOC_RECT * E 
          mat1_times_conj_mat2_into_mat1_absent_stag_phases(&loc_plaq[parity],idxh,&u[dir_nuF],idxh);                 // LOC_RECT = LOC_RECT * F

          tr_local_plaqs[parity].c[idxh] = matrix_trace_absent_stag_phase(&loc_plaq[parity],idxh);

        }  // x          
      }  // y            
    }  // z              
  }  // t                

  double res_R_p = 0.0;
  double res_I_p = 0.0;
  double resR = 0.0;
#pragma acc kernels present(tr_local_plaqs)
#pragma acc loop reduction(+:res_R_p) reduction(+:res_I_p)
  for(t=0; t<sizeh; t++) {
    res_R_p += creal(tr_local_plaqs[0].c[t]);
    res_R_p += creal(tr_local_plaqs[1].c[t]);
    res_I_p += cimag(tr_local_plaqs[0].c[t]);
    res_I_p += cimag(tr_local_plaqs[1].c[t]);
  }

  return res_R_p;
}// closes routine 



// Routine for the computation of the 5 matrices which contributes to the staples A-Right and B-Right.
// mat6 += mat1 * mat2 * hermitian_conjugate(mat3) * hermitian_conjugate(mat4) * hermitian_conjugate(mat5)
//static inline 
#pragma acc routine seq
void    PPMMM_5mat_prod_addto_mat6_absent_stag_phases(  __restrict su3_soa * const mat1,
								      const int idx_mat1,
								      __restrict su3_soa * const mat2,
								      const int idx_mat2,
								      __restrict su3_soa * const mat3,
								      const int idx_mat3,
								      __restrict su3_soa * const mat4,
								      const int idx_mat4,
								      __restrict su3_soa * const mat5,
								      const int idx_mat5,
								      __restrict su3_soa * const mat6,
								      const int idx_mat6){

  // La logica della operazioni che seguono:
  // M1 = m1
  // M2 = m2
  //                 A  = M1 * M2
  // M1 = hc(m3)
  //                 M2 = A  * M1
  // M1 = hc(m4)
  //                 A  = M2 * M1
  // M1 = hc(m5)
  //                 M2 = A  * M1
  //                 m6+= C1 * M2
  // 
  // C1 is a coeff defined in struct_c_def.c

  single_su3 AUX,MAT1,MAT2;
  int i,j;
 
  // Load rows 1 and 2 of mat1
  MAT1.comp[0][0] = mat1->r0.c0[idx_mat1];
  MAT1.comp[0][1] = mat1->r0.c1[idx_mat1];
  MAT1.comp[0][2] = mat1->r0.c2[idx_mat1];

  MAT1.comp[1][0] = mat1->r1.c0[idx_mat1];
  MAT1.comp[1][1] = mat1->r1.c1[idx_mat1];
  MAT1.comp[1][2] = mat1->r1.c2[idx_mat1];

  // Load rows 1 and 2 of mat2
  MAT2.comp[0][0] = mat2->r0.c0[idx_mat2];
  MAT2.comp[0][1] = mat2->r0.c1[idx_mat2];
  MAT2.comp[0][2] = mat2->r0.c2[idx_mat2];

  MAT2.comp[1][0] = mat2->r1.c0[idx_mat2];
  MAT2.comp[1][1] = mat2->r1.c1[idx_mat2];
  MAT2.comp[1][2] = mat2->r1.c2[idx_mat2];

  //Compute 3rd mat2 column from the first two
  MAT2.comp[2][0] = conj(MAT2.comp[0][1] * MAT2.comp[1][2] - MAT2.comp[0][2] * MAT2.comp[1][1]);
  MAT2.comp[2][1] = conj(MAT2.comp[0][2] * MAT2.comp[1][0] - MAT2.comp[0][0] * MAT2.comp[1][2]);
  MAT2.comp[2][2] = conj(MAT2.comp[0][0] * MAT2.comp[1][1] - MAT2.comp[0][1] * MAT2.comp[1][0]);


// #pragma acc loop seq
  for(i=0;i<2;i++){
// #pragma acc loop seq
    for(j=0;j<3;j++){
      AUX.comp[i][j] = (MAT1.comp[i][0] * MAT2.comp[0][j] + MAT1.comp[i][1] * MAT2.comp[1][j] + MAT1.comp[i][2] * MAT2.comp[2][j]);
    }
  }
  
  // Load columns 1 and 2 of h.c.(mat3)
  MAT1.comp[0][0] = conj(mat3->r0.c0[idx_mat3]);
  MAT1.comp[1][0] = conj(mat3->r0.c1[idx_mat3]);
  MAT1.comp[2][0] = conj(mat3->r0.c2[idx_mat3]);

  MAT1.comp[0][1] = conj(mat3->r1.c0[idx_mat3]);
  MAT1.comp[1][1] = conj(mat3->r1.c1[idx_mat3]);
  MAT1.comp[2][1] = conj(mat3->r1.c2[idx_mat3]);

  //Compute 3rd mat3 column from the first two
  MAT1.comp[0][2] = conj(MAT1.comp[1][0] * MAT1.comp[2][1] - MAT1.comp[2][0] * MAT1.comp[1][1]);
  MAT1.comp[1][2] = conj(MAT1.comp[2][0] * MAT1.comp[0][1] - MAT1.comp[0][0] * MAT1.comp[2][1]);
  MAT1.comp[2][2] = conj(MAT1.comp[0][0] * MAT1.comp[1][1] - MAT1.comp[1][0] * MAT1.comp[0][1]);

// #pragma acc loop seq
  for(i=0;i<2;i++){
// #pragma acc loop seq
    for(j=0;j<3;j++){
      MAT2.comp[i][j] = (AUX.comp[i][0] * MAT1.comp[0][j] + AUX.comp[i][1] * MAT1.comp[1][j] + AUX.comp[i][2] * MAT1.comp[2][j]);
    }
  }

  // Load columns 1 and 2 of h.c.(mat4)
  MAT1.comp[0][0] = conj(mat4->r0.c0[idx_mat4]);
  MAT1.comp[1][0] = conj(mat4->r0.c1[idx_mat4]);
  MAT1.comp[2][0] = conj(mat4->r0.c2[idx_mat4]);

  MAT1.comp[0][1] = conj(mat4->r1.c0[idx_mat4]);
  MAT1.comp[1][1] = conj(mat4->r1.c1[idx_mat4]);
  MAT1.comp[2][1] = conj(mat4->r1.c2[idx_mat4]);

  //Compute 3rd mat4 column from the first two
  MAT1.comp[0][2] = conj(MAT1.comp[1][0] * MAT1.comp[2][1] - MAT1.comp[2][0] * MAT1.comp[1][1]);
  MAT1.comp[1][2] = conj(MAT1.comp[2][0] * MAT1.comp[0][1] - MAT1.comp[0][0] * MAT1.comp[2][1]);
  MAT1.comp[2][2] = conj(MAT1.comp[0][0] * MAT1.comp[1][1] - MAT1.comp[1][0] * MAT1.comp[0][1]);

// #pragma acc loop seq
  for(i=0;i<2;i++){
// #pragma acc loop seq
    for(j=0;j<3;j++){
      AUX.comp[i][j] = (MAT2.comp[i][0] * MAT1.comp[0][j] + MAT2.comp[i][1] * MAT1.comp[1][j] + MAT2.comp[i][2] * MAT1.comp[2][j]);
    }
  }

  // Load columns 1 and 2 of h.c.(mat5)
  MAT1.comp[0][0] = conj(mat5->r0.c0[idx_mat5]);
  MAT1.comp[1][0] = conj(mat5->r0.c1[idx_mat5]);
  MAT1.comp[2][0] = conj(mat5->r0.c2[idx_mat5]);

  MAT1.comp[0][1] = conj(mat5->r1.c0[idx_mat5]);
  MAT1.comp[1][1] = conj(mat5->r1.c1[idx_mat5]);
  MAT1.comp[2][1] = conj(mat5->r1.c2[idx_mat5]);

  //Compute 3rd mat5 column from the first two
  MAT1.comp[0][2] = conj(MAT1.comp[1][0] * MAT1.comp[2][1] - MAT1.comp[2][0] * MAT1.comp[1][1]);
  MAT1.comp[1][2] = conj(MAT1.comp[2][0] * MAT1.comp[0][1] - MAT1.comp[0][0] * MAT1.comp[2][1]);
  MAT1.comp[2][2] = conj(MAT1.comp[0][0] * MAT1.comp[1][1] - MAT1.comp[1][0] * MAT1.comp[0][1]);

// #pragma acc loop seq
  for(i=0;i<2;i++){
// #pragma acc loop seq
    for(j=0;j<3;j++){
      MAT2.comp[i][j] = (AUX.comp[i][0] * MAT1.comp[0][j] + AUX.comp[i][1] * MAT1.comp[1][j] + AUX.comp[i][2] * MAT1.comp[2][j]);
    }
  }

  mat6->r0.c0[idx_mat6] += C_ONE * MAT2.comp[0][0];
  mat6->r0.c1[idx_mat6] += C_ONE * MAT2.comp[0][1];
  mat6->r0.c2[idx_mat6] += C_ONE * MAT2.comp[0][2];

  mat6->r1.c0[idx_mat6] += C_ONE * MAT2.comp[1][0];
  mat6->r1.c1[idx_mat6] += C_ONE * MAT2.comp[1][1];
  mat6->r1.c2[idx_mat6] += C_ONE * MAT2.comp[1][2];

  //Compute 3rd result column from the first two
  mat6->r2.c0[idx_mat6] += C_ONE * conj(MAT2.comp[0][1] * MAT2.comp[1][2] - MAT2.comp[0][2] * MAT2.comp[1][1]);
  mat6->r2.c1[idx_mat6] += C_ONE * conj(MAT2.comp[0][2] * MAT2.comp[1][0] - MAT2.comp[0][0] * MAT2.comp[1][2]);
  mat6->r2.c2[idx_mat6] += C_ONE * conj(MAT2.comp[0][0] * MAT2.comp[1][1] - MAT2.comp[0][1] * MAT2.comp[1][0]);

}  // closes PPMMM




// Routine for the computation of the 5 matrices which contributes to the staples C-Right and B-Left.
// mat6 += mat1 * hermitian_conjugate(mat2) * hermitian_conjugate(mat3) * hermitian_conjugate(mat4) * mat5
//static inline 
#pragma acc routine seq
void    PMMMP_5mat_prod_addto_mat6_absent_stag_phases(  __restrict su3_soa * const mat1,
								      const int idx_mat1,
								      __restrict su3_soa * const mat2,
								      const int idx_mat2,
								      __restrict su3_soa * const mat3,
								      const int idx_mat3,
								      __restrict su3_soa * const mat4,
								      const int idx_mat4,
								      __restrict su3_soa * const mat5,
								      const int idx_mat5,
								      __restrict su3_soa * const mat6,
								      const int idx_mat6){

  // La logica della operazioni che seguono:
  // M1 = m1
  // M2 = hc(m2)
  //                 A  = M1 * M2
  // M1 = hc(m3)
  //                 M2 = A  * M1
  // M1 = hc(m4)
  //                 A  = M2 * M1
  // M1 = m5
  //                 M2 = A  * M1
  //                 m6+= C1 * M2
  // 
  // C1 is a coeff defined in struct_c_def.c

  single_su3 AUX,MAT1,MAT2;
  int i,j;
 
  // Load rows 1 and 2 of mat1
  MAT1.comp[0][0] = mat1->r0.c0[idx_mat1];
  MAT1.comp[0][1] = mat1->r0.c1[idx_mat1];
  MAT1.comp[0][2] = mat1->r0.c2[idx_mat1];

  MAT1.comp[1][0] = mat1->r1.c0[idx_mat1];
  MAT1.comp[1][1] = mat1->r1.c1[idx_mat1];
  MAT1.comp[1][2] = mat1->r1.c2[idx_mat1];

  // Load columns 1 and 2 of h.c.(mat2)
  MAT2.comp[0][0] = conj(mat2->r0.c0[idx_mat2]);
  MAT2.comp[1][0] = conj(mat2->r0.c1[idx_mat2]);
  MAT2.comp[2][0] = conj(mat2->r0.c2[idx_mat2]);

  MAT2.comp[0][1] = conj(mat2->r1.c0[idx_mat2]);
  MAT2.comp[1][1] = conj(mat2->r1.c1[idx_mat2]);
  MAT2.comp[2][1] = conj(mat2->r1.c2[idx_mat2]);

  //Compute 3rd mat2 column from the first two
  MAT2.comp[0][2] = conj(MAT2.comp[1][0] * MAT2.comp[2][1] - MAT2.comp[2][0] * MAT2.comp[1][1]);
  MAT2.comp[1][2] = conj(MAT2.comp[2][0] * MAT2.comp[0][1] - MAT2.comp[0][0] * MAT2.comp[2][1]);
  MAT2.comp[2][2] = conj(MAT2.comp[0][0] * MAT2.comp[1][1] - MAT2.comp[1][0] * MAT2.comp[0][1]);

// #pragma acc loop seq
  for(i=0;i<2;i++){
// #pragma acc loop seq
    for(j=0;j<3;j++){
      AUX.comp[i][j] = (MAT1.comp[i][0] * MAT2.comp[0][j] + MAT1.comp[i][1] * MAT2.comp[1][j] + MAT1.comp[i][2] * MAT2.comp[2][j]);
    }
  }
  
  // Load columns 1 and 2 of h.c.(mat3)
  MAT1.comp[0][0] = conj(mat3->r0.c0[idx_mat3]);
  MAT1.comp[1][0] = conj(mat3->r0.c1[idx_mat3]);
  MAT1.comp[2][0] = conj(mat3->r0.c2[idx_mat3]);

  MAT1.comp[0][1] = conj(mat3->r1.c0[idx_mat3]);
  MAT1.comp[1][1] = conj(mat3->r1.c1[idx_mat3]);
  MAT1.comp[2][1] = conj(mat3->r1.c2[idx_mat3]);

  //Compute 3rd mat3 column from the first two
  MAT1.comp[0][2] = conj(MAT1.comp[1][0] * MAT1.comp[2][1] - MAT1.comp[2][0] * MAT1.comp[1][1]);
  MAT1.comp[1][2] = conj(MAT1.comp[2][0] * MAT1.comp[0][1] - MAT1.comp[0][0] * MAT1.comp[2][1]);
  MAT1.comp[2][2] = conj(MAT1.comp[0][0] * MAT1.comp[1][1] - MAT1.comp[1][0] * MAT1.comp[0][1]);

// #pragma acc loop seq
  for(i=0;i<2;i++){
// #pragma acc loop seq
    for(j=0;j<3;j++){
      MAT2.comp[i][j] = (AUX.comp[i][0] * MAT1.comp[0][j] + AUX.comp[i][1] * MAT1.comp[1][j] + AUX.comp[i][2] * MAT1.comp[2][j]);
    }
  }

  // Load columns 1 and 2 of h.c.(mat4)
  MAT1.comp[0][0] = conj(mat4->r0.c0[idx_mat4]);
  MAT1.comp[1][0] = conj(mat4->r0.c1[idx_mat4]);
  MAT1.comp[2][0] = conj(mat4->r0.c2[idx_mat4]);

  MAT1.comp[0][1] = conj(mat4->r1.c0[idx_mat4]);
  MAT1.comp[1][1] = conj(mat4->r1.c1[idx_mat4]);
  MAT1.comp[2][1] = conj(mat4->r1.c2[idx_mat4]);

  //Compute 3rd mat3 column from the first two
  MAT1.comp[0][2] = conj(MAT1.comp[1][0] * MAT1.comp[2][1] - MAT1.comp[2][0] * MAT1.comp[1][1]);
  MAT1.comp[1][2] = conj(MAT1.comp[2][0] * MAT1.comp[0][1] - MAT1.comp[0][0] * MAT1.comp[2][1]);
  MAT1.comp[2][2] = conj(MAT1.comp[0][0] * MAT1.comp[1][1] - MAT1.comp[1][0] * MAT1.comp[0][1]);

// #pragma acc loop seq
  for(i=0;i<2;i++){
// #pragma acc loop seq
    for(j=0;j<3;j++){
      AUX.comp[i][j] = (MAT2.comp[i][0] * MAT1.comp[0][j] + MAT2.comp[i][1] * MAT1.comp[1][j] + MAT2.comp[i][2] * MAT1.comp[2][j]);
    }
  }

  // Load rows 1 and 2 of mat5
  MAT1.comp[0][0] = mat5->r0.c0[idx_mat5];
  MAT1.comp[0][1] = mat5->r0.c1[idx_mat5];
  MAT1.comp[0][2] = mat5->r0.c2[idx_mat5];

  MAT1.comp[1][0] = mat5->r1.c0[idx_mat5];
  MAT1.comp[1][1] = mat5->r1.c1[idx_mat5];
  MAT1.comp[1][2] = mat5->r1.c2[idx_mat5];

  //Compute 3rd mat2 column from the first two
  MAT1.comp[2][0] = conj(MAT1.comp[0][1] * MAT1.comp[1][2] - MAT1.comp[0][2] * MAT1.comp[1][1]);
  MAT1.comp[2][1] = conj(MAT1.comp[0][2] * MAT1.comp[1][0] - MAT1.comp[0][0] * MAT1.comp[1][2]);
  MAT1.comp[2][2] = conj(MAT1.comp[0][0] * MAT1.comp[1][1] - MAT1.comp[0][1] * MAT1.comp[1][0]);

// #pragma acc loop seq
  for(i=0;i<2;i++){
// #pragma acc loop seq
    for(j=0;j<3;j++){
      MAT2.comp[i][j] = (AUX.comp[i][0] * MAT1.comp[0][j] + AUX.comp[i][1] * MAT1.comp[1][j] + AUX.comp[i][2] * MAT1.comp[2][j]);
    }
  }

  mat6->r0.c0[idx_mat6] += C_ONE * MAT2.comp[0][0];
  mat6->r0.c1[idx_mat6] += C_ONE * MAT2.comp[0][1];
  mat6->r0.c2[idx_mat6] += C_ONE * MAT2.comp[0][2];

  mat6->r1.c0[idx_mat6] += C_ONE * MAT2.comp[1][0];
  mat6->r1.c1[idx_mat6] += C_ONE * MAT2.comp[1][1];
  mat6->r1.c2[idx_mat6] += C_ONE * MAT2.comp[1][2];

  //Compute 3rd result column from the first two
  mat6->r2.c0[idx_mat6] += C_ONE * conj(MAT2.comp[0][1] * MAT2.comp[1][2] - MAT2.comp[0][2] * MAT2.comp[1][1]);
  mat6->r2.c1[idx_mat6] += C_ONE * conj(MAT2.comp[0][2] * MAT2.comp[1][0] - MAT2.comp[0][0] * MAT2.comp[1][2]);
  mat6->r2.c2[idx_mat6] += C_ONE * conj(MAT2.comp[0][0] * MAT2.comp[1][1] - MAT2.comp[0][1] * MAT2.comp[1][0]);
}  // closes PMMMP


// Routine for the computation of the 5 matrices which contributes to the staples C-Right and B-Left.
// mat6 += hermitian_conjugate(mat1) * hermitian_conjugate(mat2) * hermitian_conjugate(mat3) * mat4 * mat5
//static inline 
#pragma acc routine seq
void    MMMPP_5mat_prod_addto_mat6_absent_stag_phases(  __restrict su3_soa * const mat1,
								      const int idx_mat1,
								      __restrict su3_soa * const mat2,
								      const int idx_mat2,
								      __restrict su3_soa * const mat3,
								      const int idx_mat3,
								      __restrict su3_soa * const mat4,
								      const int idx_mat4,
								      __restrict su3_soa * const mat5,
								      const int idx_mat5,
								      __restrict su3_soa * const mat6,
								      const int idx_mat6){

  // La logica della operazioni che seguono:
  // M1 = hc(m1)
  // M2 = hc(m2)
  //                 A  = M1 * M2
  // M1 = hc(m3)
  //                 M2 = A  * M1
  // M1 = m4
  //                 A  = M2 * M1
  // M1 = m5
  //                 M2 = A  * M1
  //                 m6+= C1 * M2
  // 
  // C1 is a coeff defined in struct_c_def.c

  single_su3 AUX,MAT1,MAT2;
  int i,j;
 
  // Load columns 1 and 2 of hc(mat1)
  MAT1.comp[0][0] = conj(mat1->r0.c0[idx_mat1]);
  MAT1.comp[1][0] = conj(mat1->r0.c1[idx_mat1]);
  MAT1.comp[2][0] = conj(mat1->r0.c2[idx_mat1]);

  MAT1.comp[0][1] = conj(mat1->r1.c0[idx_mat1]);
  MAT1.comp[1][1] = conj(mat1->r1.c1[idx_mat1]);
  MAT1.comp[2][1] = conj(mat1->r1.c2[idx_mat1]);

  //Compute 3rd hc(mat1) column from the first two
  MAT1.comp[0][2] = conj(MAT1.comp[1][0] * MAT1.comp[2][1] - MAT1.comp[2][0] * MAT1.comp[1][1]);
  MAT1.comp[1][2] = conj(MAT1.comp[2][0] * MAT1.comp[0][1] - MAT1.comp[0][0] * MAT1.comp[2][1]);
  MAT1.comp[2][2] = conj(MAT1.comp[0][0] * MAT1.comp[1][1] - MAT1.comp[1][0] * MAT1.comp[0][1]); // quest'ultimo forse non serve

  // Load columns 1 and 2 of h.c.(mat2)
  MAT2.comp[0][0] = conj(mat2->r0.c0[idx_mat2]);
  MAT2.comp[1][0] = conj(mat2->r0.c1[idx_mat2]);
  MAT2.comp[2][0] = conj(mat2->r0.c2[idx_mat2]);

  MAT2.comp[0][1] = conj(mat2->r1.c0[idx_mat2]);
  MAT2.comp[1][1] = conj(mat2->r1.c1[idx_mat2]);
  MAT2.comp[2][1] = conj(mat2->r1.c2[idx_mat2]);

  //Compute 3rd mat2 column from the first two
  MAT2.comp[0][2] = conj(MAT2.comp[1][0] * MAT2.comp[2][1] - MAT2.comp[2][0] * MAT2.comp[1][1]);
  MAT2.comp[1][2] = conj(MAT2.comp[2][0] * MAT2.comp[0][1] - MAT2.comp[0][0] * MAT2.comp[2][1]);
  MAT2.comp[2][2] = conj(MAT2.comp[0][0] * MAT2.comp[1][1] - MAT2.comp[1][0] * MAT2.comp[0][1]);

// #pragma acc loop seq
  for(i=0;i<2;i++){
// #pragma acc loop seq
    for(j=0;j<3;j++){
      AUX.comp[i][j] = (MAT1.comp[i][0] * MAT2.comp[0][j] + MAT1.comp[i][1] * MAT2.comp[1][j] + MAT1.comp[i][2] * MAT2.comp[2][j]);
    }
  }
  
  // Load columns 1 and 2 of h.c.(mat3)
  MAT1.comp[0][0] = conj(mat3->r0.c0[idx_mat3]);
  MAT1.comp[1][0] = conj(mat3->r0.c1[idx_mat3]);
  MAT1.comp[2][0] = conj(mat3->r0.c2[idx_mat3]);

  MAT1.comp[0][1] = conj(mat3->r1.c0[idx_mat3]);
  MAT1.comp[1][1] = conj(mat3->r1.c1[idx_mat3]);
  MAT1.comp[2][1] = conj(mat3->r1.c2[idx_mat3]);

  //Compute 3rd mat3 column from the first two
  MAT1.comp[0][2] = conj(MAT1.comp[1][0] * MAT1.comp[2][1] - MAT1.comp[2][0] * MAT1.comp[1][1]);
  MAT1.comp[1][2] = conj(MAT1.comp[2][0] * MAT1.comp[0][1] - MAT1.comp[0][0] * MAT1.comp[2][1]);
  MAT1.comp[2][2] = conj(MAT1.comp[0][0] * MAT1.comp[1][1] - MAT1.comp[1][0] * MAT1.comp[0][1]);

// #pragma acc loop seq
  for(i=0;i<2;i++){
// #pragma acc loop seq
    for(j=0;j<3;j++){
      MAT2.comp[i][j] = (AUX.comp[i][0] * MAT1.comp[0][j] + AUX.comp[i][1] * MAT1.comp[1][j] + AUX.comp[i][2] * MAT1.comp[2][j]);
    }
  }

  // Load rows 1 and 2 of mat4
  MAT1.comp[0][0] = mat4->r0.c0[idx_mat4];
  MAT1.comp[0][1] = mat4->r0.c1[idx_mat4];
  MAT1.comp[0][2] = mat4->r0.c2[idx_mat4];

  MAT1.comp[1][0] = mat4->r1.c0[idx_mat4];
  MAT1.comp[1][1] = mat4->r1.c1[idx_mat4];
  MAT1.comp[1][2] = mat4->r1.c2[idx_mat4];

  //Compute 3rd mat4 column from the first two
  MAT1.comp[2][0] = conj(MAT1.comp[0][1] * MAT1.comp[1][2] - MAT1.comp[0][2] * MAT1.comp[1][1]);
  MAT1.comp[2][1] = conj(MAT1.comp[0][2] * MAT1.comp[1][0] - MAT1.comp[0][0] * MAT1.comp[1][2]);
  MAT1.comp[2][2] = conj(MAT1.comp[0][0] * MAT1.comp[1][1] - MAT1.comp[0][1] * MAT1.comp[1][0]);

// #pragma acc loop seq
  for(i=0;i<2;i++){
// #pragma acc loop seq
    for(j=0;j<3;j++){
      AUX.comp[i][j] = (MAT2.comp[i][0] * MAT1.comp[0][j] + MAT2.comp[i][1] * MAT1.comp[1][j] + MAT2.comp[i][2] * MAT1.comp[2][j]);
    }
  }

  // Load rows 1 and 2 of mat5
  MAT1.comp[0][0] = mat5->r0.c0[idx_mat5];
  MAT1.comp[0][1] = mat5->r0.c1[idx_mat5];
  MAT1.comp[0][2] = mat5->r0.c2[idx_mat5];

  MAT1.comp[1][0] = mat5->r1.c0[idx_mat5];
  MAT1.comp[1][1] = mat5->r1.c1[idx_mat5];
  MAT1.comp[1][2] = mat5->r1.c2[idx_mat5];

  //Compute 3rd mat2 column from the first two
  MAT1.comp[2][0] = conj(MAT1.comp[0][1] * MAT1.comp[1][2] - MAT1.comp[0][2] * MAT1.comp[1][1]);
  MAT1.comp[2][1] = conj(MAT1.comp[0][2] * MAT1.comp[1][0] - MAT1.comp[0][0] * MAT1.comp[1][2]);
  MAT1.comp[2][2] = conj(MAT1.comp[0][0] * MAT1.comp[1][1] - MAT1.comp[0][1] * MAT1.comp[1][0]);

// #pragma acc loop seq
  for(i=0;i<2;i++){
// #pragma acc loop seq
    for(j=0;j<3;j++){
      MAT2.comp[i][j] = (AUX.comp[i][0] * MAT1.comp[0][j] + AUX.comp[i][1] * MAT1.comp[1][j] + AUX.comp[i][2] * MAT1.comp[2][j]);
    }
  }

  mat6->r0.c0[idx_mat6] += C_ONE * MAT2.comp[0][0];
  mat6->r0.c1[idx_mat6] += C_ONE * MAT2.comp[0][1];
  mat6->r0.c2[idx_mat6] += C_ONE * MAT2.comp[0][2];

  mat6->r1.c0[idx_mat6] += C_ONE * MAT2.comp[1][0];
  mat6->r1.c1[idx_mat6] += C_ONE * MAT2.comp[1][1];
  mat6->r1.c2[idx_mat6] += C_ONE * MAT2.comp[1][2];

  //Compute 3rd result column from the first two
  mat6->r2.c0[idx_mat6] += C_ONE * conj(MAT2.comp[0][1] * MAT2.comp[1][2] - MAT2.comp[0][2] * MAT2.comp[1][1]);
  mat6->r2.c1[idx_mat6] += C_ONE * conj(MAT2.comp[0][2] * MAT2.comp[1][0] - MAT2.comp[0][0] * MAT2.comp[1][2]);
  mat6->r2.c2[idx_mat6] += C_ONE * conj(MAT2.comp[0][0] * MAT2.comp[1][1] - MAT2.comp[0][1] * MAT2.comp[1][0]);
}  // closes MMMPP









void calc_loc_improved_staples_typeA_removing_stag_phases_nnptrick_all(  __restrict su3_soa * const u,
									 __restrict su3_soa * const loc_stap ){
  int x, y, z, t, mu, iter;
#pragma acc kernels present(u) present(loc_stap) present(nnp_openacc) present(nnm_openacc)
 #pragma acc loop independent gang
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector(4)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector(4)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(32)
        for(x=0; x < nx; x++) {
     #pragma acc loop seq
          for(mu=0; mu<4; mu++){
      #pragma acc loop seq
            for(iter=0; iter<3; iter++){

              int nu;
              if (mu==0) { nu = iter + 1; }
              else if (mu==1) { nu = iter + (iter & 1) + (iter >> 1); }
              else if (mu==2) { nu = iter + (iter >> 1); }
              else if (mu==3) { nu = iter; }
              else { //error
              }

              const int idxh = snum_acc(x,y,z,t);  // r
              const int parity = (x+y+z+t) % 2;

              const int dir_link = 2*mu + parity;
              const int dir_nu_1L = 2*nu +  parity;
              const int dir_nu_2L = 2*nu + !parity;
              const int dir_mu_3L = 2*mu +  parity;
              const int dir_nu_4L = 2*nu +  parity;
              const int dir_nu_5L = 2*nu + !parity;
              const int dir_nu_1R = 2*nu + !parity;
              const int dir_nu_2R = 2*nu +  parity;
              const int dir_mu_3R = 2*mu +  parity;
              const int dir_nu_4R = 2*nu + !parity;
              const int dir_nu_5R = 2*nu +  parity;

#pragma acc cache (nnp_openacc[idxh:8])
              const int idx_pmu = nnp_openacc[idxh][mu][parity];          // r+mu
#pragma acc cache (nnm_openacc[idx_pmu:8])
              const int idx_pnu = nnp_openacc[idxh][nu][parity];          // r+nu
              const int idx_mnu = nnm_openacc[idxh][nu][parity] ;         // r-nu
              const int idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu
              const int idx_pmu_pnu = nnp_openacc[idx_pmu][nu][!parity];  // r+mu+nu
              const int idx_2pnu = nnp_openacc[idx_pnu][nu][!parity];     // r+2nu
              const int idx_pmu_2mnu = nnm_openacc[idx_pmu_mnu][nu][parity];  // r+mu-2nu
              const int idx_2mnu = nnm_openacc[idx_mnu][nu][!parity] ;         // r-2nu

	      // IMPROVEMENT TYPE A
	      //   AR is of type PPMMM
	      //   AL is of type MMMPP
	      //    r+mu-2nu  r+mu-nu  r+mu   r+mu+nu  r+mu+2nu
	      //          +<-----+<-----+----->+----->+
	      //          |  2L     1L  ^  1R     2R  |
	      // mu    3L |             |             | 3R
	      // ^        V  4L     5L  |  5R     4R  V
	      // |        +----->+----->+<-----+<-----+ 
	      // |     r-2nu    r-nu    r     r+nu     r+2nu
	      // +---> nu
	      //            r is idxh in the following

	      //computation of the Right part of the staple
              PPMMM_5mat_prod_addto_mat6_absent_stag_phases(&u[dir_nu_1R],        idx_pmu,
							    &u[dir_nu_2R],        idx_pmu_pnu,
							    &u[dir_mu_3R],        idx_2pnu,
							    &u[dir_nu_4R],        idx_pnu,
							    &u[dir_nu_5R],        idxh,
							    &loc_stap[dir_link],  idxh);							    
	      
              //computation of the Left  part of the staple
              MMMPP_5mat_prod_addto_mat6_absent_stag_phases(&u[dir_nu_1L],        idx_pmu_mnu,
							    &u[dir_nu_2L],        idx_pmu_2mnu,
							    &u[dir_mu_3L],        idx_2mnu,
							    &u[dir_nu_4L],        idx_2mnu,
							    &u[dir_nu_5L],        idx_mnu,
							    &loc_stap[dir_link],  idxh);							    
	      
            }  // mu
          }  // iter
        }  // x
      }  // y
    }  // z
  }  // t

}// closes routine





void calc_loc_improved_staples_typeB_removing_stag_phases_nnptrick_all(  __restrict su3_soa * const u,
									 __restrict su3_soa * const loc_stap ){
  int x, y, z, t, mu, iter;
#pragma acc kernels present(u) present(loc_stap) present(nnp_openacc) present(nnm_openacc)
 #pragma acc loop independent gang
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector(4)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector(4)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(32)
        for(x=0; x < nx; x++) {
     #pragma acc loop seq
          for(mu=0; mu<4; mu++){
      #pragma acc loop seq
            for(iter=0; iter<3; iter++){

              int nu;
              if (mu==0) { nu = iter + 1; }
              else if (mu==1) { nu = iter + (iter & 1) + (iter >> 1); }
              else if (mu==2) { nu = iter + (iter >> 1); }
              else if (mu==3) { nu = iter; }
              else { //error
              }

              const int idxh = snum_acc(x,y,z,t);  // r
              const int parity = (x+y+z+t) % 2;

              const int dir_link = 2*mu + parity;
              const int dir_mu_1L = 2*mu + !parity;
              const int dir_nu_2L = 2*nu + !parity;
              const int dir_mu_3L = 2*mu +  parity;
              const int dir_mu_4L = 2*mu + !parity;
              const int dir_nu_5L = 2*nu + !parity;
              const int dir_mu_1R = 2*mu + !parity;
              const int dir_nu_2R = 2*nu +  parity;
              const int dir_mu_3R = 2*mu +  parity;
              const int dir_mu_4R = 2*mu + !parity;
              const int dir_nu_5R = 2*nu +  parity;

#pragma acc cache (nnp_openacc[idxh:8])
              const int idx_pmu = nnp_openacc[idxh][mu][parity];             // r+mu
#pragma acc cache (nnm_openacc[idx_pmu:8])
              const int idx_pnu = nnp_openacc[idxh][nu][parity];             // r+nu
              const int idx_mnu = nnm_openacc[idxh][nu][parity] ;            // r-nu
              const int idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];     // r+mu-nu
              const int idx_pmu_pnu = nnp_openacc[idx_pmu][nu][!parity];     // r+mu+nu
              const int idx_2pmu = nnp_openacc[idx_pmu][mu][!parity];        // r+2mu
              const int idx_2pmu_mnu = nnp_openacc[idx_pmu_mnu][mu][parity]; // r+2mu-nu


	      // IMPROVEMENT TYPE B      
	      //   BR is of type PPMMM
	      //   BL is of type PMMMP
	      //      r+2mu-nu r+2mu   r+2mu+nu
	      //          +<-----+----->+
	      //          |  2L  ^  2R  |
	      //       3L |    1L|1R    | 3R
	      //          V      |      V
	      //  r+mu-nu +     r+mu    + r+mu+nu
	      //          |      ^      |
	      // mu    4L |      |      | 4R
	      // ^        V  5L  |  5R  V
	      // |        +----->+<-----+
	      // |       r-nu    r     r+nu
	      // +---> nu
	      //            r is idxh in the following

	      //computation of the Right part of the staple
              PPMMM_5mat_prod_addto_mat6_absent_stag_phases(&u[dir_mu_1R],        idx_pmu,
							    &u[dir_nu_2R],        idx_2pmu,
							    &u[dir_mu_3R],        idx_pmu_pnu,
							    &u[dir_mu_4R],        idx_pnu,
							    &u[dir_nu_5R],        idxh,
							    &loc_stap[dir_link],  idxh);
	      
              //computation of the Left  part of the staple
              PMMMP_5mat_prod_addto_mat6_absent_stag_phases(&u[dir_mu_1L],        idx_pmu,
							    &u[dir_nu_2L],        idx_2pmu_mnu,
							    &u[dir_mu_3L],        idx_pmu_mnu,
							    &u[dir_mu_4L],        idx_mnu,
							    &u[dir_nu_5L],        idx_mnu,
							    &loc_stap[dir_link],  idxh);
	      
            }  // mu
          }  // iter
        }  // x
      }  // y
    }  // z
  }  // t

}// closes routine




void calc_loc_improved_staples_typeC_removing_stag_phases_nnptrick_all(  __restrict su3_soa * const u,
									 __restrict su3_soa * const loc_stap ){
  int x, y, z, t, mu, iter;
#pragma acc kernels present(u) present(loc_stap) present(nnp_openacc) present(nnm_openacc)
 #pragma acc loop independent gang
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector(4)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector(4)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(32)
        for(x=0; x < nx; x++) {
     #pragma acc loop seq
          for(mu=0; mu<4; mu++){
      #pragma acc loop seq
            for(iter=0; iter<3; iter++){

              int nu;
              if (mu==0) { nu = iter + 1; }
              else if (mu==1) { nu = iter + (iter & 1) + (iter >> 1); }
              else if (mu==2) { nu = iter + (iter >> 1); }
              else if (mu==3) { nu = iter; }
              else { //error
              }

              const int idxh = snum_acc(x,y,z,t);  // r
              const int parity = (x+y+z+t) % 2;

              const int dir_link = 2*mu + parity;
              const int dir_nu_1L = 2*nu +  parity;
              const int dir_mu_2L = 2*mu + !parity;
              const int dir_mu_3L = 2*mu +  parity;
              const int dir_nu_4L = 2*nu +  parity;
              const int dir_mu_5L = 2*mu + !parity;
              const int dir_nu_1R = 2*nu + !parity;
              const int dir_mu_2R = 2*mu + !parity;
              const int dir_mu_3R = 2*mu +  parity;
              const int dir_nu_4R = 2*nu + !parity;
              const int dir_mu_5R = 2*mu + !parity;

#pragma acc cache (nnp_openacc[idxh:8])
              const int idx_pmu = nnp_openacc[idxh][mu][parity];          // r+mu
#pragma acc cache (nnm_openacc[idx_pmu:8])
              const int idx_pnu = nnp_openacc[idxh][nu][parity];          // r+nu
              const int idx_mnu = nnm_openacc[idxh][nu][parity] ;         // r-nu
              const int idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu
              const int idx_mmu = nnm_openacc[idxh][mu][parity];          // r-mu
              const int idx_mmu_pnu = nnp_openacc[idx_mmu][nu][!parity];  // r-mu+nu
              const int idx_mmu_mnu = nnm_openacc[idx_mmu][nu][!parity];  // r-mu-nu


	      // IMPROVEMENT TYPE C
	      //   CR is of type PMMMP
	      //   CL is of type MMMPP
	      //       r+mu-nu  r+mu   r+mu+nu
	      //          +<-----+----->+
	      //          |  1L  ^  1R  |
	      //       2L |      |      | 2R
	      //          V      |      V
	      //     r-nu +      r      + r+nu
	      //          |      ^      |
	      // mu    3L |    5L|5R    | 3R
	      // ^        V  4L  |  4R  V
	      // |        +----->+<-----+
	      // |    r-mu-nu   r-mu   r-mu+nu
	      // +---> nu
	      //            r is idxh in the following
	      
	      //computation of the Right part of the staple
              PMMMP_5mat_prod_addto_mat6_absent_stag_phases(&u[dir_nu_1R],        idx_pmu,
							    &u[dir_mu_2R],        idx_pnu,
							    &u[dir_mu_3R],        idx_mmu_pnu,
							    &u[dir_nu_4R],        idx_mmu,
							    &u[dir_mu_5R],        idx_mmu,
							    &loc_stap[dir_link],  idxh);							    
	      
              //computation of the Left  part of the staple
              MMMPP_5mat_prod_addto_mat6_absent_stag_phases(&u[dir_nu_1L],        idx_pmu_mnu,
							    &u[dir_mu_2L],        idx_mnu,
							    &u[dir_mu_3L],        idx_mmu_mnu,
							    &u[dir_nu_4L],        idx_mmu_mnu,
							    &u[dir_mu_5L],        idx_mmu,
							    &loc_stap[dir_link],  idxh);							    
	      
            }  // mu
          }  // iter
        }  // x
      }  // y
    }  // z
  }  // t

}// closes routine



void calc_loc_improved_staples_typeABC_removing_stag_phases_nnptrick_all(  __restrict su3_soa * const u,
									   __restrict su3_soa * const loc_stap ){
  int x, y, z, t, mu, iter;
#pragma acc kernels present(u) present(loc_stap) present(nnp_openacc) present(nnm_openacc)
 #pragma acc loop independent gang
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector(4)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector(4)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(32)
        for(x=0; x < nx; x++) {
     #pragma acc loop seq
          for(mu=0; mu<4; mu++){
      #pragma acc loop seq
            for(iter=0; iter<3; iter++){

              int nu;
              if (mu==0) { nu = iter + 1; }
              else if (mu==1) { nu = iter + (iter & 1) + (iter >> 1); }
              else if (mu==2) { nu = iter + (iter >> 1); }
              else if (mu==3) { nu = iter; }
              else { //error
              }

              const int idxh = snum_acc(x,y,z,t);  // r
              const int parity = (x+y+z+t) % 2;

              const int dir_link = 2*mu + parity;

              const int Adir_nu_1L = 2*nu +  parity;
              const int Adir_nu_2L = 2*nu + !parity;
              const int Adir_mu_3L = 2*mu +  parity;
              const int Adir_nu_4L = 2*nu +  parity;
              const int Adir_nu_5L = 2*nu + !parity;
              const int Adir_nu_1R = 2*nu + !parity;
              const int Adir_nu_2R = 2*nu +  parity;
              const int Adir_mu_3R = 2*mu +  parity;
              const int Adir_nu_4R = 2*nu + !parity;
              const int Adir_nu_5R = 2*nu +  parity;

              const int Bdir_mu_1L = 2*mu + !parity;
              const int Bdir_nu_2L = 2*nu + !parity;
              const int Bdir_mu_3L = 2*mu +  parity;
              const int Bdir_mu_4L = 2*mu + !parity;
              const int Bdir_nu_5L = 2*nu + !parity;
              const int Bdir_mu_1R = 2*mu + !parity;
              const int Bdir_nu_2R = 2*nu +  parity;
              const int Bdir_mu_3R = 2*mu +  parity;
              const int Bdir_mu_4R = 2*mu + !parity;
              const int Bdir_nu_5R = 2*nu +  parity;


              const int Cdir_nu_1L = 2*nu +  parity;
              const int Cdir_mu_2L = 2*mu + !parity;
              const int Cdir_mu_3L = 2*mu +  parity;
              const int Cdir_nu_4L = 2*nu +  parity;
              const int Cdir_mu_5L = 2*mu + !parity;
              const int Cdir_nu_1R = 2*nu + !parity;
              const int Cdir_mu_2R = 2*mu + !parity;
              const int Cdir_mu_3R = 2*mu +  parity;
              const int Cdir_nu_4R = 2*nu + !parity;
              const int Cdir_mu_5R = 2*mu + !parity;

#pragma acc cache (nnp_openacc[idxh:8])
              const int idx_pmu = nnp_openacc[idxh][mu][parity];          // r+mu
#pragma acc cache (nnm_openacc[idx_pmu:8])
              const int idx_pnu = nnp_openacc[idxh][nu][parity];          // r+nu
              const int idx_mnu = nnm_openacc[idxh][nu][parity] ;         // r-nu
              const int idx_pmu_mnu = nnm_openacc[idx_pmu][nu][!parity];  // r+mu-nu
              const int idx_mmu = nnm_openacc[idxh][mu][parity];          // r-mu
              const int idx_mmu_pnu = nnp_openacc[idx_mmu][nu][!parity];  // r-mu+nu
              const int idx_mmu_mnu = nnm_openacc[idx_mmu][nu][!parity];  // r-mu-nu
              const int idx_pmu_pnu = nnp_openacc[idx_pmu][nu][!parity];     // r+mu+nu
              const int idx_2pmu = nnp_openacc[idx_pmu][mu][!parity];        // r+2mu
              const int idx_2pmu_mnu = nnp_openacc[idx_pmu_mnu][mu][parity]; // r+2mu-nu
              const int idx_2pnu = nnp_openacc[idx_pnu][nu][!parity];     // r+2nu
              const int idx_pmu_2mnu = nnm_openacc[idx_pmu_mnu][nu][parity];  // r+mu-2nu
              const int idx_2mnu = nnm_openacc[idx_mnu][nu][!parity] ;         // r-2nu              



	      // IMPROVEMENT TYPE ABC
              //computation of the Right part of the A staple
              PPMMM_5mat_prod_addto_mat6_absent_stag_phases(&u[Adir_nu_1R],        idx_pmu,
                                                            &u[Adir_nu_2R],        idx_pmu_pnu,
                                                            &u[Adir_mu_3R],        idx_2pnu,
                                                            &u[Adir_nu_4R],        idx_pnu,
                                                            &u[Adir_nu_5R],        idxh,
                                                            &loc_stap[dir_link],  idxh);

              //computation of the Left  part of the A staple
              MMMPP_5mat_prod_addto_mat6_absent_stag_phases(&u[Adir_nu_1L],        idx_pmu_mnu,
                                                            &u[Adir_nu_2L],        idx_pmu_2mnu,
                                                            &u[Adir_mu_3L],        idx_2mnu,
                                                            &u[Adir_nu_4L],        idx_2mnu,
                                                            &u[Adir_nu_5L],        idx_mnu,
                                                            &loc_stap[dir_link],  idxh);



              //computation of the Right part of the B staple
              PPMMM_5mat_prod_addto_mat6_absent_stag_phases(&u[Bdir_mu_1R],        idx_pmu,
                                                            &u[Bdir_nu_2R],        idx_2pmu,
                                                            &u[Bdir_mu_3R],        idx_pmu_pnu,
                                                            &u[Bdir_mu_4R],        idx_pnu,
                                                            &u[Bdir_nu_5R],        idxh,
                                                            &loc_stap[dir_link],  idxh);

              //computation of the Left  part of the B staple
              PMMMP_5mat_prod_addto_mat6_absent_stag_phases(&u[Bdir_mu_1L],        idx_pmu,
                                                            &u[Bdir_nu_2L],        idx_2pmu_mnu,
                                                            &u[Bdir_mu_3L],        idx_pmu_mnu,
                                                            &u[Bdir_mu_4L],        idx_mnu,
                                                            &u[Bdir_nu_5L],        idx_mnu,
                                                            &loc_stap[dir_link],  idxh);


	      
	      //computation of the Right part of the C staple
              PMMMP_5mat_prod_addto_mat6_absent_stag_phases(&u[Cdir_nu_1R],        idx_pmu,
							    &u[Cdir_mu_2R],        idx_pnu,
							    &u[Cdir_mu_3R],        idx_mmu_pnu,
							    &u[Cdir_nu_4R],        idx_mmu,
							    &u[Cdir_mu_5R],        idx_mmu,
							    &loc_stap[dir_link],  idxh);							    
	      
              //computation of the Left  part of the C staple
              MMMPP_5mat_prod_addto_mat6_absent_stag_phases(&u[Cdir_nu_1L],        idx_pmu_mnu,
							    &u[Cdir_mu_2L],        idx_mnu,
							    &u[Cdir_mu_3L],        idx_mmu_mnu,
							    &u[Cdir_nu_4L],        idx_mmu_mnu,
							    &u[Cdir_mu_5L],        idx_mmu,
							    &loc_stap[dir_link],  idxh);							    
	      
            }  // mu
          }  // iter
        }  // x
      }  // y
    }  // z
  }  // t

}// closes routine






double  calc_rettangolo_soloopenacc( __restrict  su3_soa * const tconf_acc, __restrict su3_soa * const local_plaqs, dcomplex_soa * const tr_local_plaqs){

  double tempo=0.0;
  // tolgo le fasi staggered
  mult_conf_times_stag_phases(tconf_acc);
  // calcolo il valore della plaquette sommata su tutti i siti a fissato piano mu-nu (6 possibili piani)
  for(int mu=0;mu<3;mu++){
    for(int nu=mu+1;nu<4;nu++){
      // sommo i 12 risultati in tempo
      tempo  += calc_loc_rectangles_2x1_removing_stag_phases_nnptrick(tconf_acc,local_plaqs,tr_local_plaqs,mu,nu);
      tempo  += calc_loc_rectangles_1x2_removing_stag_phases_nnptrick(tconf_acc,local_plaqs,tr_local_plaqs,mu,nu);
    }
  }
  // rimetto le fasi staggered
  mult_conf_times_stag_phases(tconf_acc);

  return tempo;

}



#endif
