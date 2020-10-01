#ifndef TOPOLOGICAL_FORCE_C
#define TOPOLOGICAL_FORCE_C
#include "../Include/common_defines.h"
#include "./topological_force.h"
#include <math.h>
#include "../Mpi/multidev.h"
#include "../Meas/gauge_meas.h"
#include "./fermion_force_utilities.h"
#include "./alloc_vars.h"
#include "./action.h"
#include "./topological_action.h"

#define ALLOLLOCHECK(control_int,var)  if(control_int != 0 ) \
    printf("MPI%02d: \tError in  allocation of %s . \n",devinfo.myrank, #var);\
    else if(verbosity_lv > 2) printf("MPI%02d: \tAllocation of %s : OK , %p\n",\
         devinfo.myrank, #var, var );\

#define FREELOLLOCHECK(var) if(verbosity_lv >2) \
    printf("\tFreed %s, %p ...", #var,var);\
    free(var); if(verbosity_lv > 2)  printf(" done.\n");\


double compute_topodynamical_potential_der(const double Q)
{
	double barr=act_params.barrier, width=act_params.width;
	int ngrid= (int) floor((2*barr+width/2)/width);
	int kgrid= (int) floor((Q+barr)/width);
	
	if(kgrid>=0 && kgrid<=ngrid)
		{
			double grid[ngrid];
			if(verbosity_lv>3)
				printf("\t\t\tMPI%02d - load_topo(path,barr,width,grid,ngrid)\n",devinfo.myrank);
			
			load_topo(act_params.topo_file_path,barr,width,grid,ngrid);
			return (grid[kgrid+1]-grid[kgrid])/width;
		}
	else
		return 0;
}

void four_leaves(su3_soa * const leaves, su3_soa * const u)
{
  int mu, nu;
    int d0, d1, d2, d3;
	
#pragma acc kernels present(u) present(leaves) present(nnp_openacc) present(nnm_openacc)
 for(mu=0;mu<3;mu++)
   for(nu=mu+1;nu<4;nu++){
#pragma acc loop independent gang
     for(d3 = D3_HALO; d3<nd3-D3_HALO; d3++){
#pragma acc loop independent gang vector
       for(d2=0; d2<nd2; d2++){
#pragma acc loop independent gang vector
	 for(d1=0; d1<nd1; d1++){
#pragma acc loop independent vector 
	   for(d0=0; d0 < nd0; d0++){
	     const int idxh = snum_acc(d0,d1,d2,d3);
	     const int parity = (d0+d1+d2+d3)%2;
	     const int idx_plane = 2*mu + nu - 1 - ((int)mu/2);//mu=0 nu=1 --> idx_plane=0;... mu=1 nu=2 --> idx_plane=3;... mu=2 nu=3 --> idx_plane=5.
	     //leaves points
	     int idx_p_mu = nnp_openacc[idxh][mu][parity];
	     int idx_p_nu = nnp_openacc[idxh][nu][parity];
	     int idx_m_mu = nnm_openacc[idxh][mu][parity];
	     int idx_m_nu = nnm_openacc[idxh][nu][parity];

	     int idx_p_nu_m_mu = nnm_openacc[idx_p_nu][mu][!parity];
	     int idx_m_mu_m_nu = nnm_openacc[idx_m_mu][nu][!parity];
	     int idx_m_nu_p_mu = nnp_openacc[idx_m_nu][mu][!parity];
							
	     //first leave
	     comp_U_U_Udag_Udag(&u[2*mu+parity] , idxh,
				&u[2*nu+!parity], idx_p_mu,
				&u[2*mu+!parity], idx_p_nu,
				&u[2*nu+parity] , idxh,
				&leaves[idx_plane+parity], idxh);

	     comp_and_add_U_Udag_Udag_U(&u[2*nu+parity] , idxh,
					&u[2*mu+parity] , idx_p_nu_m_mu,
					&u[2*nu+!parity], idx_m_mu,
					&u[2*mu+!parity], idx_m_mu,
					&leaves[idx_plane+parity], idxh);

	     comp_and_add_Udag_Udag_U_U(&u[2*mu+!parity], idx_m_mu,           
					&u[2*nu+parity] , idx_m_mu_m_nu,  
					&u[2*mu+parity] , idx_m_mu_m_nu,      
					&u[2*nu+!parity], idx_m_nu,      
					&leaves[idx_plane+parity], idxh);

	     comp_and_add_Udag_U_U_Udag(&u[2*nu+!parity], idx_m_nu,           
					&u[2*mu+!parity], idx_m_nu,  
					&u[2*nu+parity] , idx_m_nu_p_mu,      
					&u[2*mu+parity] , idxh,      
					&leaves[idx_plane+parity], idxh);
							
	   }//d0
	 }//d1
       }//d2
     }//d3
   }//nu
 
}

void antihermatize_unsafe(su3_soa * const leaves)
{
	int d0,d1,d2,d3,plane;
#pragma acc kernels present(leaves)
	for(plane=0;plane<6;plane++)
#pragma acc loop independent gang
		for(d3 = D3_HALO; d3<nd3-D3_HALO; d3++)
#pragma acc loop independent gang vector
			for(d2=0; d2<nd2; d2++)
#pragma acc loop independent gang vector
				for(d1=0; d1<nd1; d1++){
					const int idxh = snum_acc(d0,d1,d2,d3);
					const int parity = (d0+d1+d2+d3)%2;

					leaves[plane+parity].r0.c0[idxh] = cimag(leaves[plane+parity].r0.c0[idxh]) * _Complex_I;
					leaves[plane+parity].r1.c1[idxh] = cimag(leaves[plane+parity].r1.c1[idxh]) * _Complex_I;
					leaves[plane+parity].r2.c2[idxh] = cimag(leaves[plane+parity].r2.c2[idxh]) * _Complex_I;

					d_complex tmp01 = conj(leaves[plane+parity].r1.c0[idxh]);
					d_complex tmp02 = conj(leaves[plane+parity].r2.c0[idxh]);
					d_complex tmp12 = conj(leaves[plane+parity].r2.c1[idxh]);
					
					leaves[plane+parity].r1.c0[idxh] = 0.5 * (leaves[plane+parity].r1.c0[idxh] - conj(leaves[plane+parity].r0.c1[idxh]));
					leaves[plane+parity].r2.c0[idxh] = 0.5 * (leaves[plane+parity].r2.c0[idxh] - conj(leaves[plane+parity].r0.c2[idxh]));
					leaves[plane+parity].r2.c1[idxh] = 0.5 * (leaves[plane+parity].r2.c1[idxh] - conj(leaves[plane+parity].r1.c2[idxh]));

					leaves[plane+parity].r0.c1[idxh] = 0.5 * (leaves[plane+parity].r0.c1[idxh] - tmp01);
					leaves[plane+parity].r0.c2[idxh] = 0.5 * (leaves[plane+parity].r0.c2[idxh] - tmp02);
					leaves[plane+parity].r1.c2[idxh] = 0.5 * (leaves[plane+parity].r1.c2[idxh] - tmp12);
					
				}//closing all the loops at once

}

void topo_staples(__restrict su3_soa * const u,__restrict su3_soa * const staples, double norm)
{
  //compute leaves
  su3_soa * leaves;
  posix_memalign((void **)&leaves,128,12*sizeof(su3_soa));
#pragma acc enter data create(leaves[0:12]) copyin(norm[0:1])

  if(verbosity_lv>3)
    printf("\t\t\tMPI%02d - four_leaves(leaves,u)\n",devinfo.myrank);
  
  //compute leaves
  
  four_leaves(leaves, u);
  
  if(verbosity_lv>3)
    printf("\t\t\tMPI%02d - antihermatize_usafe(leaves)\n",devinfo.myrank);
  
  //antihermatize leaves
  
  antihermatize_unsafe(leaves);
  int d0, d1, d2, d3;
  if(verbosity_lv>3)
    {
      printf("MPI%d - computing staples\n",devinfo.myrank);
    }
  
  su3_soa *ABC, *BCF, *ABCF, *ADE, *DEF, *ADEF, *temp_r, *temp_l;
  int allocation_check;
  

  allocation_check = posix_memalign((void **)&ABC,128,8*sizeof(su3_soa));
  ALLOLLOCHECK(allocation_check,ABC);

#pragma acc enter data create(ABC[0:8])
  
  allocation_check = posix_memalign((void **)&BCF,128,8*sizeof(su3_soa));
  ALLOLLOCHECK(allocation_check,BCF);
  
#pragma acc enter data create(BCF[0:8])
  
  allocation_check = posix_memalign((void **)&ABCF,128,8*sizeof(su3_soa));
  ALLOLLOCHECK(allocation_check,ABCF);
  
#pragma acc enter data create(ABCF[0:8])
  
  allocation_check = posix_memalign((void **)&ADE,128,8*sizeof(su3_soa));
  ALLOLLOCHECK(allocation_check,ADE);
  
#pragma acc enter data create(ADE[0:8])
  
  allocation_check = posix_memalign((void **)&DEF,128,8*sizeof(su3_soa));
  ALLOLLOCHECK(allocation_check,DEF);
  
#pragma acc enter data create(DEF[0:8])
  
  allocation_check = posix_memalign((void **)&ADEF,128,8*sizeof(su3_soa));
  ALLOLLOCHECK(allocation_check,ADEF);
  
#pragma acc enter data create(ADEF[0:8])
  
  allocation_check = posix_memalign((void **)&temp_r,128,8*sizeof(su3_soa));
  ALLOLLOCHECK(allocation_check,temp_r);
  
#pragma acc enter data create(temp_r[0:8])
  
  allocation_check = posix_memalign((void **)&temp_l,128,8*sizeof(su3_soa));
  ALLOLLOCHECK(allocation_check,temp_l);
  
#pragma acc enter data create(temp_l[0:8])
  

  for(int mu=0; mu<4; mu++)
    for(int inu=0; inu<3; inu++){
#pragma acc data present(u[0:8]) present(nnp_openacc) present(nnm_openacc) present(staples) \
  present(ABC[0:8]) present(leaves[0:12]) present(BCF[0:8]) present(ABCF[0:8]) \
  present(ADE[0:8]) present(DEF[0:8]) present(ADEF[0:8])		\
  present(temp_l[0:8]) present(temp_r[0:8])
      
#pragma acc parallel loop independent gang
      for(d3 = D3_HALO; d3<nd3-D3_HALO; d3++){
#pragma acc loop independent gang vector
	for(d2=0; d2<nd2; d2++){
#pragma acc loop independent gang vector
	  for(d1=0; d1<nd1; d1++){
#pragma acc loop independent gang vector
	    for(d0=0; d0<nd0; d0++){


	      int perp_dir[4][3] = {{ 1, 2, 3}, { 0, 2, 3}, { 0, 1, 3}, { 0, 1, 2}};
	      int plan_perp[4][3]= {{ 5, 4, 3}, { 5, 2, 1}, { 4, 2, 0}, { 3, 1, 0}};
	      int plan_sign[4][3]= {{+1,-1,+1}, {-1,+1,-1}, {+1,-1,+1}, {-1,+1,-1}};

	      const int nu=perp_dir[mu][inu];                  //  E---F---C   
	      const int idxA = snum_acc(d0,d1,d2,d3);          //  |   |   | mu
	      const int parity = (d0+d1+d2+d3)%2;	       //  D---A---B
	      const int iplan=plan_perp[mu][inu];              //        nu
	      
	      
	      const int idxF = nnp_openacc[idxA][mu][parity];
	      const int idxB = nnp_openacc[idxA][nu][parity];
	      const int idxC = nnp_openacc[idxB][mu][!parity];
	      const int idxD = nnm_openacc[idxA][nu][parity];
	      const int idxE = nnp_openacc[idxD][mu][!parity];
						
	      //compute ABC BCF ABCF
	      mat1_times_mat2_into_mat3_absent_stag_phases(&u[nu+parity], idxA,
							   &u[mu+!parity],idxB,
							   &ABC[mu+parity],idxA);

	      
		mat1_times_conj_mat2_into_mat3_absent_stag_phases(&u[mu+!parity], idxB,
								  &u[nu+!parity],idxF,
								  &BCF[mu+parity],idxA);
		mat1_times_conj_mat2_into_mat3_absent_stag_phases(&ABC[mu+parity], idxA,
								  &u[nu+!parity],idxF,
								  &ABCF[mu+parity],idxA);

		//compute ADE DEF ADEF
		conj_mat1_times_mat2_into_mat3_absent_stag_phases(&u[nu+!parity], idxD,
								  &u[mu+!parity],idxD,
								  &ADE[mu+parity],idxA);
		mat1_times_mat2_into_mat3_absent_stag_phases(&u[nu+!parity], idxD,
							     &u[mu+parity],idxE,
							     &DEF[mu+parity],idxA);
		conj_mat1_times_mat2_into_mat3_absent_stag_phases(&u[nu+!parity], idxD,
								  &DEF[mu+parity],idxA,
								  &ADEF[mu+parity],idxA);

		//compute local staples. Beware, this step should differ from nissa's
		mat1_times_mat2_times_fact_addto_mat3_absent_stag_phases(&leaves[iplan+parity], idxA, //--
									 &ABCF[iplan+parity], idxA,   // |
									 &staples[mu+parity],idxA,    // |
									 norm*plan_sign[mu][inu]);    // |
		mat1_times_mat2_into_mat3_absent_stag_phases(&u[nu+parity], idxA,		      // |
							     &leaves[iplan+!parity], idxB,	      // |
							     &temp_r[nu+parity], idxA);               // |
		mat1_times_mat2_times_fact_addto_mat3_absent_stag_phases(&temp_r[nu+parity], idxA,    // |
									 &BCF[iplan+parity], idxA,    // |
									 &staples[mu+parity],idxA,    // |
									 norm*plan_sign[mu][inu]);    // |
		mat1_times_conj_mat2_into_mat3_absent_stag_phases(&leaves[iplan+parity],idxC,         //   RIGHT STAPLE    
								  &u[nu+!parity], idxF,               // |
								  &temp_r[nu+!parity], idxF);	      // |
		mat1_times_mat2_times_fact_addto_mat3_absent_stag_phases(&ABC[iplan+parity], idxA,    // |
									 &temp_r[nu+!parity], idxF,   // |
									 &staples[mu+parity],idxA,    // |
									 norm*plan_sign[mu][inu]);    // |
		mat1_times_mat2_times_fact_addto_mat3_absent_stag_phases(&ABCF[iplan+parity], idxA,   // |
									 &leaves[iplan+parity], idxF, // |
									 &staples[mu+parity],idxA,    // |
									 norm*plan_sign[mu][inu]);    //--

		conj_mat1_times_mat2_times_fact_addto_mat3_absent_stag_phases(&leaves[iplan+parity], idxA, //--
									      &ADEF[iplan+parity], idxA,   // |
									      &staples[mu+parity],idxA,    // |
									      norm*plan_sign[mu][inu]);    // |
		conj_mat1_times_conj_mat2_into_mat3_absent_stag_phases(&u[nu+parity], idxD,		   // |
								       &leaves[iplan+!parity], idxD,       // |
								       &temp_l[nu+parity], idxD);          // |
		mat1_times_mat2_times_fact_addto_mat3_absent_stag_phases(&temp_l[nu+parity], idxD,	   // |
									 &DEF[iplan+parity], idxA,	   // |
									 &staples[mu+parity],idxA,	   // |
									 norm*plan_sign[mu][inu]);	   // |
		conj_mat1_times_mat2_into_mat3_absent_stag_phases(&leaves[iplan+parity],idxE,              //   LEFT STAPLE    
								  &u[nu+!parity], idxF,                    // |
								  &temp_l[nu+!parity], idxF);	           // |
		mat1_times_mat2_times_fact_addto_mat3_absent_stag_phases(&ADE[iplan+parity], idxA,	   // |
									 &temp_l[nu+!parity], idxF,	   // |
									 &staples[mu+parity],idxA,         // |
									 norm*plan_sign[mu][inu]);	   // |
		mat1_times_conj_mat2_times_fact_addto_mat3_absent_stag_phases(&ADEF[iplan+parity], idxA,   // |
									      &leaves[iplan+parity], idxF, // |
									      &staples[mu+parity],idxA,    // |
									      norm*plan_sign[mu][inu]);    //--

	    }//d0
	  }//d1
	}//d2
      }//d3
      printf("this should appear 12 times\n");
    }//inu

  FREELOLLOCHECK(ABC);
#pragma acc exit data delete(ABC[:8])
  FREELOLLOCHECK(BCF);
#pragma acc exit data delete(BCF[:8])
  FREELOLLOCHECK(ABCF);
#pragma acc exit data delete(ABCF[:8])
  FREELOLLOCHECK(ADE);
#pragma acc exit data delete(ADE[:8])
  FREELOLLOCHECK(DEF);
#pragma acc exit data delete(DEF[:8])
  FREELOLLOCHECK(ADEF);
#pragma acc exit data delete(ADEF[:8])
  FREELOLLOCHECK(temp_r);
#pragma acc exit data delete(temp_r[:8])
  FREELOLLOCHECK(temp_l);
#pragma acc exit data delete(temp_l[:8])
  FREELOLLOCHECK(leaves);
#pragma acc exit data delete(leaves[:12]) 

#pragma acc exit data delete(norm)

	
}



void calc_loc_topo_staples(__restrict su3_soa * const u, __restrict su3_soa * const staples)
{
	su3_soa tstout_conf_acc_arr, *quadri;
	double_soa * loc_q;

	if(verbosity_lv>3)
	  printf("\t\tMPI%02d - compute_topological_charge(u,quadri,loc_q)\n",devinfo.myrank);
	
	posix_memalign((void **)&quadri,128,8*sizeof(su3_soa));
#pragma acc enter data create(quadri[0:8])
	posix_memalign((void **)&loc_q,128,2*sizeof(double_soa));
#pragma acc enter data create(loc_q[0:2])

	double Q = compute_topological_charge(u, quadri, loc_q);

	free(quadri);
#pragma acc exit data delete(quadri)
	free(loc_q);
#pragma acc exit data delete(loc_q)

	if(verbosity_lv>3)
	  printf("Topological Charge: %lf\n",Q);


	if(verbosity_lv>3)
	  printf("\t\tMPI%02d - compute_topodynamical_potential_der(Q)\n",devinfo.myrank);
	
	double pot_der=compute_topodynamical_potential_der(Q);
	double norm=pot_der/(M_PI*M_PI*64);
	
	if(verbosity_lv>3)
	  printf("\t\tMPI%02d - topo_staples(u, staples, norm)\n",devinfo.myrank);
	
	topo_staples(u, staples, norm);
	
#ifdef STOUT_TOPO
	stout_wrapper(u,tstout_conf_acc_arr); //INSERIRE PARAMETRI STOUTING TOPOLOGICO
#else
	tstout_conf_acc_arr=*u;
#endif
	int d0, d1, d2, d3;
	for(int mu=0; mu<4; mu++)
#pragma acc loop independent gang
		for(d3 = D3_HALO; d3<nd3-D3_HALO; d3++)
#pragma acc loop independent gang vector
			for(d2=0; d2<nd2; d2++)
#pragma acc loop independent gang vector
				for(d2=0; d2<nd2; d2++){
					const int idxm = snum_acc(d0,d1,d2,d3);
					const int parity = (d0+d1+d2+d3)%2;
					mat1_times_double_factor(&staples[mu+parity], idxm, norm);
		}

}


#endif
