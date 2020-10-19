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
#include "./single_types.h"

//#define DEBUG_LOLLO

#ifdef DEBUG_LOLLO
 #include "./cayley_hamilton.h"
 #include "./alloc_vars.h"
 #include "../DbgTools/dbgtools.h"


#define ALLOLLOCHECK(control_int,var)  if(control_int != 0 ) \
    printf("MPI%02d: \tError in  allocation of %s . \n",devinfo.myrank, #var);\
    else if(verbosity_lv > 2) printf("MPI%02d: \tAllocation of %s : OK , %p\n",\
         devinfo.myrank, #var, var );

#define FREELOLLOCHECK(var) if(verbosity_lv >2) \
    printf("\tFreed %s, %p ...", #var,var);\
    free(var); if(verbosity_lv > 2)  printf(" done.\n");

#endif

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

void four_leaves(__restrict su3_soa * const leaves,__restrict su3_soa * const u)
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
				 &leaves[2*idx_plane+parity], idxh);
	      comp_and_add_U_Udag_Udag_U(&u[2*nu+parity] , idxh,
					 &u[2*mu+parity] , idx_p_nu_m_mu,
					 &u[2*nu+!parity], idx_m_mu,
					 &u[2*mu+!parity], idx_m_mu,
					 &leaves[2*idx_plane+parity], idxh);
	      comp_and_add_Udag_Udag_U_U(&u[2*mu+!parity], idx_m_mu,           
					 &u[2*nu+parity] , idx_m_mu_m_nu,  
					 &u[2*mu+parity] , idx_m_mu_m_nu,      
					 &u[2*nu+!parity], idx_m_nu,      
					 &leaves[2*idx_plane+parity], idxh);
	      comp_and_add_Udag_U_U_Udag(&u[2*nu+!parity], idx_m_nu,           
					 &u[2*mu+!parity], idx_m_nu,  
					 &u[2*nu+parity] , idx_m_nu_p_mu,      
					 &u[2*mu+parity] , idxh,      
					 &leaves[2*idx_plane+parity], idxh);
	    }//d0
	  }//d1
	}//d2
      }//d3
    }//nu
}

void antihermatize_unsafe(__restrict su3_soa * const leaves)
{
  int d0,d1,d2,d3,plane;

#pragma acc kernels present(leaves)
#pragma acc loop independent vector
  for(plane=0; plane<12; plane++)
#pragma acc loop independent gang(STAPGANG3)
    for(d3 = D3_HALO; d3<nd3-D3_HALO; d3++)
#pragma acc loop independent tile(STAPTILE0,STAPTILE1,STAPTILE2)
      for(d2=0; d2<nd2; d2++)
	for(d1=0; d1<nd1; d1++)
	  for(d0=0; d0<nd0; d0++)
	    {
	      const int idx = snum_acc(d0,d1,d2,d3);
	      
	      leaves[plane].r0.c0[idx] = 2* cimag(leaves[plane].r0.c0[idx]) * _Complex_I;
	      leaves[plane].r1.c1[idx] = 2* cimag(leaves[plane].r1.c1[idx]) * _Complex_I;
	      leaves[plane].r2.c2[idx] = 2* cimag(leaves[plane].r2.c2[idx]) * _Complex_I;
	      
	      d_complex tmp01 = conj(leaves[plane].r1.c0[idx]);
	      d_complex tmp02 = conj(leaves[plane].r2.c0[idx]);
	      d_complex tmp12 = conj(leaves[plane].r2.c1[idx]);
	      
	      leaves[plane].r1.c0[idx] = (leaves[plane].r1.c0[idx] - conj(leaves[plane].r0.c1[idx]));
	      leaves[plane].r2.c0[idx] = (leaves[plane].r2.c0[idx] - conj(leaves[plane].r0.c2[idx]));
	      leaves[plane].r2.c1[idx] = (leaves[plane].r2.c1[idx] - conj(leaves[plane].r1.c2[idx]));
	      
	      leaves[plane].r0.c1[idx] = (leaves[plane].r0.c1[idx] - tmp01);
	      leaves[plane].r0.c2[idx] = (leaves[plane].r0.c2[idx] - tmp02);
	      leaves[plane].r1.c2[idx] = (leaves[plane].r1.c2[idx] - tmp12);
	      
	    }//closing all the loops at once
}

void topo_staples(__restrict su3_soa * const u,__restrict su3_soa * const staples, double norm)
{
  //compute leaves
  su3_soa * leaves;
  posix_memalign((void **)&leaves,128,12*sizeof(su3_soa));
#pragma acc enter data create(leaves[0:12]) copyin(norm)

  if(verbosity_lv>3)
    printf("\t\t\tMPI%02d - four_leaves(leaves,u)\n",devinfo.myrank);
  
  //compute leaves
  
  four_leaves(leaves,u);

  
  if(verbosity_lv>3)
    printf("\t\t\tMPI%02d - antihermatize_usafe(leaves)\n",devinfo.myrank);
  
  //antihermatize leaves
  
  antihermatize_unsafe(leaves);
  
  int d0, d1, d2, d3;
  if(verbosity_lv>3)
    {
      printf("MPI%d - computing staples\n",devinfo.myrank);
    }  
  
  for(int mu=0; mu<4; mu++)
    for(int inu=0; inu<3; inu++){
#pragma acc parallel loop independent gang present(leaves) present(u) present(nnp_openacc) present(nnm_openacc) present(staples)
      for(d3 = D3_HALO; d3<nd3-D3_HALO; d3++){
#pragma acc loop independent gang vector
	for(d2=0; d2<nd2; d2++){
#pragma acc loop independent gang vector
	  for(d1=0; d1<nd1; d1++){
#pragma acc loop independent gang vector
	    for(d0=0; d0<nd0; d0++){
	      single_su3 ABC, BCF, ABCF;
	      single_su3 ADE, DEF, ADEF;
	      single_su3 temp;
	      single_su3 loc_stap;
	      set_to_zero_single_su3(&loc_stap);
	      
	      int perp_dir[4][3] = {{ 1, 2, 3}, { 0, 2, 3}, { 0, 1, 3}, { 0, 1, 2}};
	      int plan_perp[4][3]= {{ 5, 4, 3}, { 5, 2, 1}, { 4, 2, 0}, { 3, 1, 0}};
	      int plan_sign[4][3]= {{+1,-1,+1}, {-1,+1,-1}, {+1,-1,+1}, {-1,+1,-1}};
	      
	      const int nu = perp_dir[mu][inu];                  //  E---F---C   
	      const int idxA = snum_acc(d0,d1,d2,d3);            //  |   |   | mu
	      const int parity = (d0+d1+d2+d3)%2;                //  D---A---B
	      const int iplan = plan_perp[mu][inu];              //        nu
	      
	      const int idxF = nnp_openacc[idxA][mu][parity];
	      const int idxB = nnp_openacc[idxA][nu][parity];
	      const int idxC = nnp_openacc[idxB][mu][!parity];
	      const int idxD = nnm_openacc[idxA][nu][parity];
	      const int idxE = nnp_openacc[idxD][mu][!parity];

	      //compute ABC BCF ABCF
	      su3_soa_times_su3_soa_into_single_su3(&u[2*nu+parity],  idxA,
	      					    &u[2*mu+!parity], idxB,
	      					    &ABC);
	      su3_soa_times_su3_soa_dag_into_single_su3(&u[2*mu+!parity], idxB,
	      						&u[2*nu+!parity], idxF,
							&BCF);
	      su3_soa_times_single_su3_into_single_su3(&u[2*nu+parity],idxA,
	      					       &BCF, &ABCF);
	      //compute ADE DEF ADEF
	      su3_soa_dag_times_su3_soa_into_single_su3(&u[2*nu+!parity], idxD,
	      						&u[2*mu+!parity], idxD,
	      						&ADE);
	      su3_soa_times_su3_soa_into_single_su3(&u[2*mu+!parity], idxD,
	      					    &u[2*nu+parity],  idxE,
	      					    &DEF);
	      single_su3_times_su3_soa_into_single_su3(&ADE,
	      					       &u[2*nu+parity], idxE,
	      					       &ADEF);
	      //compute local staples.
	      gl3_soa_times_single_su3_addto_gl3(&leaves[2*iplan+parity],idxA,
	      					 &ABCF,&loc_stap);
	      su3_soa_times_gl3_soa_into_gl3(&u[2*nu+parity],idxA,
	      				     &leaves[2*iplan+!parity],idxB,
	      				     &temp);
	      single_gl3xsu3_add_to_out(&loc_stap,&temp,&BCF);
	      gl3_soa_times_su3_soa_dag_into_gl3(&leaves[2*iplan+parity],idxC,
	      					 &u[2*nu+!parity],       idxF,
	      					 &temp);
              single_su3xgl3_add_to_out(&loc_stap,&ABC,&temp);
	      single_su3_times_gl3_soa_addto_gl3(&ABCF,
	      					 &leaves[2*iplan+!parity],idxF,
	      					 &loc_stap);
	      //^--Right staple v--Left staple
	      gl3_soa_dag_times_single_su3_addto_gl3(&leaves[2*iplan+parity],idxA,
	      					     &ADEF,&loc_stap);
	      su3_soa_dag_times_gl3_soa_dag_into_gl3(&u[2*nu+!parity],        idxD,
	      					     &leaves[2*iplan+!parity],idxD,
	      					     &temp);
	      single_gl3xsu3_add_to_out(&loc_stap,&temp,&DEF);
	      gl3_soa_dag_times_su3_soa_into_gl3(&leaves[2*iplan+parity],idxE,
	      					 &u[2*nu+parity],        idxE,
	      					 &temp);
	      single_su3xgl3_add_to_out(&loc_stap,&ADE,&temp);
	      single_su3_times_gl3_soa_dag_addto_gl3(&ADEF,
	      					 &leaves[2*iplan+!parity],idxF,
	      					 &loc_stap);
	      //finalizing
	      if(plan_sign[mu][inu] >= 0)
		single_gl3_addinto_su3_soa(&staples[2*mu+parity], idxA,&loc_stap);
	      else
		single_gl3_subtinto_su3_soa(&staples[2*mu+parity], idxA,&loc_stap);		
	    }//d0
	  }//d1
	}//d2
      }//d3
    }//inu

#pragma acc exit data delete(leaves)
  free(leaves);

#pragma acc parallel loop independent gang present(staples) present(norm)
  for(d3 = D3_HALO; d3<nd3-D3_HALO; d3++){
#pragma acc loop independent gang vector
    for(d2=0; d2<nd2; d2++){
#pragma acc loop independent gang vector
      for(d1=0; d1<nd1; d1++){
#pragma acc loop independent gang vector  
	for(d0=0; d0<nd0; d0++){
	
	const int idxA = snum_acc(d0,d1,d2,d3);
	const int parity = (d0+d1+d2+d3)%2;
	
	for(int mu = 0; mu<4; mu++)
	  gl3_dag_times_double_factor(&staples[2*mu+parity],idxA,norm);
	
	}
      }
    }	    
  }
}



void calc_loc_topo_staples(__restrict const su3_soa * const u, __restrict su3_soa * const staples)
{    
#ifdef DEBUG_LOLLO
 set_su3_soa_to_zero(staples);
#endif


  su3_soa tstout_conf_acc_arr, *quadri;
  double_soa * loc_q;
  
  if(verbosity_lv>3)
    printf("\t\tMPI%02d - compute_topological_charge(u,quadri,loc_q)\n",devinfo.myrank);
  
  posix_memalign((void **)&quadri,128,8*sizeof(su3_soa));
#pragma acc enter data create(quadri[0:8])
  posix_memalign((void **)&loc_q,128,2*sizeof(double_soa));
#pragma acc enter data create(loc_q[0:2])
  
  double Q = compute_topological_charge(u, quadri, loc_q);
  
#pragma acc exit data delete(quadri)  
  free(quadri);
#pragma acc exit data delete(loc_q)
  free(loc_q);
  
  if(verbosity_lv>4)
    printf("Topological Charge: %lf\n",Q);
  
  
  if(verbosity_lv>3)
    printf("\t\tMPI%02d - compute_topodynamical_potential_der(Q)\n",devinfo.myrank);
	
  double pot_der=compute_topodynamical_potential_der(Q);
  double norm=pot_der/(M_PI*M_PI*128);

  if(verbosity_lv>3)
    printf("\t\tMPI%02d - topo_staples(u, staples, norm)\n",devinfo.myrank);
  
  topo_staples(u,staples,norm);
  
#ifdef STOUT_TOPO
  stout_wrapper(u,tstout_conf_acc_arr); //INSERIRE PARAMETRI STOUTING TOPOLOGICO
#else
  tstout_conf_acc_arr=*u;
#endif
  
  
#ifdef DEBUG_LOLLO
  tamat_soa *tipdot;
  int allocation_check =  posix_memalign((void **)&tipdot, 128, 8*sizeof(tamat_soa)); 
  ALLOLLOCHECK(allocation_check, tipdot) ;
#pragma acc enter data create(tipdot[0:8])
  set_tamat_soa_to_zero(tipdot);
#pragma acc update self(tipdot[0:8]) self(u[0:8]) self(staples[0:8])
  printf("check multiplication into tamat:\n");
  STAMPA_DEBUG_SU3_SOA(u,0,0);
  STAMPA_DEBUG_SU3_SOA(staples,0,0);
  
  conf_times_staples_ta_part(u,staples,tipdot);
  
#pragma acc update self(tipdot[0:8])
  STAMPA_DEBUG_TAMAT_SOA(tipdot,0,0);
  
#define SQRT_3 1.732050807568877
  //i*Gell-mann matrices as from eq.A.10 of Gattringer - note that T=lambda/2
  single_tamat i_gell_mann_matr[8]={ 
    { 0+1*I, 0+0*I, 0+0*I, 0, 0 },
    { 1+0*I, 0+0*I, 0+0*I, 0, 0 },
    { 0+0*I, 0+0*I, 0+0*I, 1, -1 },
    { 0+0*I, 0+1*I, 0+0*I, 0, 0 },
    { 0+0*I, 1+0*I, 0+0*I, 0, 0 },
    { 0+0*I, 0+0*I, 0+1*I, 0, 0 },
    { 0+0*I, 0+0*I, 1+0*I, 0, 0 },
    { 0+0*I, 0+0*I, 0+0*I, 1/SQRT_3, 1/SQRT_3 }
  };
  double eps=1e-4;
  
  //store initial link and comp action
  single_su3 sto;
#pragma acc update self(u[0:8])
  single_su3_from_su3_soa(&u[0],0,&sto);
  rebuild3row(&sto);
  
  double ori_act = compute_topo_action(u);
  //store derivative
  single_tamat posi={ 0+0*I, 0+0*I, 0+0*I, 0, 0 };
  single_tamat nega={ 0+0*I, 0+0*I, 0+0*I, 0, 0 };
  
  for(int igen=0;igen<8;igen++)
    {	    
      //prepare increment and change
      single_tamat ba;
      single_tamat_times_scalar_into_tamat(&ba,&i_gell_mann_matr[igen],eps/2+0*I);
      
      single_su3 exp_mod;
      CH_exponential_antihermitian_nissalike(&exp_mod,&ba);
      rebuild3row(&exp_mod);
      
      //change +, compute action
      single_su3 shilink;
      single_su3xsu3(&shilink, &exp_mod, &sto);
      single_su3_into_su3_soa(&u[0], 0, &shilink);
      
#pragma acc update device(u[0:8])
      
      double act_minus = compute_topo_action(u);
      
      //change -, compute action
      gl3_dagger(&exp_mod);
      single_su3xsu3(&shilink, &exp_mod, &sto);
      single_su3_into_su3_soa(&u[0], 0, &shilink);
      
#pragma acc update device(u[0:8])
      
      double act_plus = compute_topo_action(u);

      //set back everything
      single_su3_into_su3_soa(&u[0], 0, &sto);
	    
#pragma acc update device(u[0:8])

      double check_ori_act = compute_topo_action(u);
      /* printf("plus = %+016.016le ",act_plus); */
      /* printf("ori = %+016.016le ",check_ori_act); */
      /* printf("minus = %+016.016le\n",act_minus);	     */
	    
      double gr_plus = -(act_plus-ori_act)/eps;
      double gr_minus = -(ori_act-act_minus)/eps;
      single_tamat_times_scalar_add_to_tamat(&posi,&i_gell_mann_matr[igen],gr_plus+0*I);
      single_tamat_times_scalar_add_to_tamat(&nega,&i_gell_mann_matr[igen],gr_minus+0*I);

      /* printf("gr_plus = %+016.016le\ngr_minus = %+016.016le\n", gr_plus, gr_minus); */
    }
  //take the average
  single_tamat Numerical_derivative;
  summ_single_tamats_times_scalar(&Numerical_derivative, &posi, &nega, 0.5+0*I);
  STAMPA_DEBUG_SINGLE_TAMAT(Numerical_derivative);
	
  printf("Ringrazia Iddio se si assomigliano almeno un po'\n");
  FREELOLLOCHECK(tipdot);
#pragma acc exit data delete(tipdot)
  mem_free_core();
  mem_free_extended();
  exit(0);
#endif


}

#endif
