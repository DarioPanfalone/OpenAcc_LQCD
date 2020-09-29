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

static inline void where_error(int a){printf("Where is the error? %d\n",a);}

void four_leaves(su3_soa * const leaves, su3_soa * const u)
{
	int mu, nu;
	int d0, d1, d2, d3;
	
where_error(0);	
#pragma acc kernels present(u) present(leaves) present(nnp_openacc) present(nnm_openacc)
	for(mu=0;mu<3;mu++)
		for(nu=mu+1;nu<4;nu++)
			{
#pragma acc loop independent gang
				for(d3 = D3_HALO; d3<nd3-D3_HALO; d3++);
#pragma acc loop independent gang vector
					for(d2=0; d2<nd2; d2++)
#pragma acc loop independent gang vector
					for(d1=0; d1<nd1; d1++)
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
							where_error(1);
							
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
							
						}//closing all loops at once
					
					
			}
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

#pragma acc enter data create(leaves[0:12]) create(norm[0:1])

#pragma kernels present(u[:8]) present(nnp_openacc) present(nnm_openacc) present(staples)
	

	if(verbosity_lv>3)
		printf("\t\t\tMPI%02d - four_leaves(leaves,u)\n",devinfo.myrank);
		
	//compute leaves

	four_leaves(leaves, u);

	if(verbosity_lv>3)
		printf("\t\t\tMPI%02d - antihermatize_usafe(leaves)\n",devinfo.myrank);

	//antihermatize leaves

	antihermatize_unsafe(leaves);
	int d0, d1, d2, d3;
	int perp_dir[4][3] = {{ 1, 2, 3}, { 0, 2, 3}, { 0, 1, 3}, { 0, 1, 2}};
	int plan_perp[4][3]= {{ 5, 4, 3}, { 5, 2, 1}, { 4, 2, 0}, { 3, 1, 0}};
	int plan_sign[4][3]= {{+1,-1,+1}, {-1,+1,-1}, {+1,-1,+1}, {-1,+1,-1}};
#pragma acc enter data create(perp_dir[:][:]) create(plan_perp[:][:]) create(plan_sign[:][:])
	if(verbosity_lv>3)
		{
			printf("MPI%d - computing staples\n",devinfo.myrank);
		}
	
	su3_soa * ABC, * BCF, * ABCF, * ADE, * DEF, * ADEF, * temp_r, * temp_l;


	posix_memalign((void **)&BCF,128,8*sizeof(su3_soa));
#pragma acc enter data create(BCF[:8])

	posix_memalign((void **)&ABCF,128,8*sizeof(su3_soa));
#pragma acc enter data create(ABCF[:8])

	posix_memalign((void **)&ADE,128,8*sizeof(su3_soa));
#pragma acc enter data create(ADE[:8])

	posix_memalign((void **)&DEF,128,8*sizeof(su3_soa));
#pragma acc enter data create(DEF[:8])

	posix_memalign((void **)&ADEF,128,8*sizeof(su3_soa));
#pragma acc enter data create(ADEF[:8])

	posix_memalign((void **)&temp_r,128,8*sizeof(su3_soa));
#pragma acc enter data create(temp_r[:8])

	posix_memalign((void **)&temp_l,128,8*sizeof(su3_soa));
#pragma acc enter data create(temp_l[:8])
	
	for(int mu=0; mu<4; mu++)
		for(int inu=0; inu<3; inu++)
#pragma acc loop independent gang
			for(d3 = D3_HALO; d3<nd3-D3_HALO; d3++)
#pragma acc loop independent gang vector
				for(d2=0; d2<nd2; d2++)
#pragma acc loop independent gang vector
					for(d2=0; d2<nd2; d2++){
						const int idxA = snum_acc(d0,d1,d2,d3);          //  E---F---C   
						const int parity = (d0+d1+d2+d3)%2;				 //  |   |   | mu
																		 //  D---A---B   
																		 //        nu    


						int nu=perp_dir[mu][inu];
						int iplan=plan_perp[mu][inu];
						int idxF = nnp_openacc[idxA][mu][parity];
						int idxB = nnp_openacc[idxA][nu][parity];
						int idxC = nnp_openacc[idxB][mu][!parity];
						int idxD = nnm_openacc[idxA][nu][parity];
						int idxE = nnp_openacc[idxD][mu][!parity];
						
						printf("mu = %d, nu = %d, idxA = %d, parity = %d, idxB = %d, idxC = %d, idxD = %d, idxE = %d, idxF = %d\n",mu,nu,idxA,parity,idxB,idxC,idxD,idxE,idxF);
						printf("mu+parity = %d, nu+parity =%d, mu+!parity = %d, nu+!parity = %d\n",mu+parity,nu+parity,mu+!parity,nu+!parity);
						posix_memalign((void **)&ABC,128,8*sizeof(su3_soa));
#pragma acc enter data create(ABC[:8])

						printf("I can access these elements: ABC[%d],(idx: %d):\n%lf, %lf, %lf\n%lf, %lf, %lf\n%lf, %lf, %lf\n",mu+parity,idxA,creal(ABC[mu+parity].r0.c0[idxA]),creal(ABC[mu+parity].r0.c1[idxA]),creal(ABC[mu+parity].r0.c2[idxA]),creal(ABC[mu+parity].r1.c0[idxA]),creal(ABC[mu+parity].r1.c1[idxA]),creal(ABC[mu+parity].r1.c2[idxA]),creal(ABC[mu+parity].r2.c0[idxA]),creal(ABC[mu+parity].r2.c1[idxA]),creal(ABC[mu+parity].r2.c2[idxA]));

						printf("I'm computing right part\n");
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
						printf("I'm computing left part\n");
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

						printf("I'm computing right staple\n");
						//compute local staples. Beware, this step should differ from nissa's
						mat1_times_mat2_times_fact_addto_mat3_absent_stag_phases(&leaves[iplan+parity], idxA, //--
																				 &ABCF[iplan+parity], idxA,   // |
																				 &staples[mu+parity],idxA,    // |
																				 norm*plan_sign[mu][nu]);    // |
						mat1_times_mat2_into_mat3_absent_stag_phases(&u[nu+parity], idxA,					 // |
																	 &leaves[iplan+!parity], idxB,			 // |
																	 &temp_r[nu+parity], idxA);               // |
						mat1_times_mat2_times_fact_addto_mat3_absent_stag_phases(&temp_r[nu+parity], idxA,	 // |
																				 &BCF[iplan+parity], idxA,	 // |
																				 &staples[mu+parity],idxA,	 // |
																				 norm*plan_sign[mu][nu]);    // |
						mat1_times_conj_mat2_into_mat3_absent_stag_phases(&leaves[iplan+parity],idxC,         //   RIGHT STAPLE    
																		  &u[nu+!parity], idxF,               // |
																		  &temp_r[nu+!parity], idxF);		 //	|
						mat1_times_mat2_times_fact_addto_mat3_absent_stag_phases(&ABC[iplan+parity], idxA,	 //	|
																				 &temp_r[nu+!parity], idxF,	 //	|
																				 &staples[mu+parity],idxA,	 //	|
																				 norm*plan_sign[mu][nu]);	 //	|
						mat1_times_mat2_times_fact_addto_mat3_absent_stag_phases(&ABCF[iplan+parity], idxA,   //	|
																				 &leaves[iplan+parity], idxF, //	|
																				 &staples[mu+parity],idxA,	 //	|
																				 norm*plan_sign[mu][nu]);    //--

						printf("I'm computing left staple\n");						
						conj_mat1_times_mat2_times_fact_addto_mat3_absent_stag_phases(&leaves[iplan+parity], idxA, //--
																					  &ADEF[iplan+parity], idxA,   // |
																					  &staples[mu+parity],idxA,    // |
																					  norm*plan_sign[mu][nu]);    // |
						conj_mat1_times_conj_mat2_into_mat3_absent_stag_phases(&u[nu+parity], idxD,			      // |
																			   &leaves[iplan+!parity], idxD,		  // |
																			   &temp_l[nu+parity], idxD);          // |
						mat1_times_mat2_times_fact_addto_mat3_absent_stag_phases(&temp_l[nu+parity], idxD,	      // |
																				 &DEF[iplan+parity], idxA,	      // |
																				 &staples[mu+parity],idxA,	      // |
																				 norm*plan_sign[mu][nu]);	      // |
						conj_mat1_times_mat2_into_mat3_absent_stag_phases(&leaves[iplan+parity],idxE,              //   LEFT STAPLE    
																		  &u[nu+!parity], idxF,                    // |
																		  &temp_l[nu+!parity], idxF);		      // |
						mat1_times_mat2_times_fact_addto_mat3_absent_stag_phases(&ADE[iplan+parity], idxA,	      // |
																				 &temp_l[nu+!parity], idxF,	      // |
																				 &staples[mu+parity],idxA,         // |
																				 norm*plan_sign[mu][nu]);	      // |
						mat1_times_conj_mat2_times_fact_addto_mat3_absent_stag_phases(&ADEF[iplan+parity], idxA,   // |
																					  &leaves[iplan+parity], idxF, // |
																					  &staples[mu+parity],idxA,	  // |
																					  norm*plan_sign[mu][nu]);    //--
												
					}
	
}



void calc_loc_topo_staples(__restrict su3_soa * const u, __restrict su3_soa * const staples)
{
	su3_soa tstout_conf_acc_arr, *quadri;
	double_soa * loc_q;
	posix_memalign((void **)&quadri,128,8*sizeof(su3_soa));
	posix_memalign((void **)&loc_q,128,2*sizeof(double_soa));
#pragma acc enter data create(quadri[0:8]) create(loc_q[0:2])
	
	if(verbosity_lv>3)
			printf("\t\tMPI%02d - compute_topological_charge(u,quadri,loc_q)\n",devinfo.myrank);
	

	double Q = compute_topological_charge(u, quadri, loc_q);

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
	free(quadri);
	free(loc_q);
}


#endif
