double casuale(void);
void initrand(unsigned long s);
void su2_rand(double *pp);
#define PRINT_DETAILS_INSIDE_UPDATE


#include "openacc.h"
#include "../RationalApprox/rationalapprox.c"
#include "./struct_c_def.c"
#include "./alloc_vars.c"
#include "./dbgtools.c"
#include "./fermionic_utilities.c"
#include "./single_types.c"
#include "./su3_utilities.c"
#include "./random_assignement.c"
#include "./fermion_matrix.c"
#include "./inverter_full.c"
#include "./find_min_max.c"
#include "./inverter_multishift_full.c"
#include "./rettangoli.c"
#include "./ipdot_gauge.c"
#include "../Meas/ferm_meas.c"
#include "./homebrew_acos.c"
#include "./stouting.c"
#include "./stouting_deottimizzato.c"
#include "./fermion_force.c"
#include "./md_integrator.c"
#include "./update_versatile.c"
int main(){

  initrand(111);
  fflush(stdout);
  printf("INIZIO DEL PROGRAMMA \n");
  su3_soa  * conf_acc;
  int  allocation_check =  posix_memalign((void **)&conf_acc, ALIGN, 8*sizeof(su3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di conf_acc \n");
  SETREQUESTED(conf_acc);
  printf("Allocazione della configurazione : OK \n");

  // INIT FERM PARAMS AND READ RATIONAL APPROX COEFFS
  init_ferm_params();

  mem_alloc();
  printf("Allocazione della memoria : OK \n");
  initialize_global_variables();
  printf("init vars : OK \n");
  compute_nnp_and_nnm_openacc();
  printf("nn computation : OK \n");
#ifdef BACKFIELD
  init_backfield();
  print_double_soa(u1_back_field_phases,"backfield");
  printf("u1_backfield initialization : OK \n");
#endif

  generate_Conf_cold(conf_acc,0.0);
  //  generate_Conf_cold(conf_acc,0.1);
  printf("Cold Gauge Conf Generated : OK \n");
  conf_id_iter=0;


#ifdef STOUT_FERMIONS        
    su3_soa *tstout_conf_acc_arr = gstout_conf_acc_arr;
#endif

    double minmaxeig[2]; 
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
	int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;
	generate_vec3_soa_gauss(&ferm_phi_acc[ps_index]);
      }
    }// end for iflav


#ifdef STOUT_FERMIONS 
    // USO DELLA VERSIONE STOUTATA GIA' PER LO STIRACCHIAMENTO
    // STOUTING...(ALREADY ON DEVICE)
    stout_wrapper(conf_acc,tstout_conf_acc_arr,0);
    gconf_as_fermionmatrix = &tstout_conf_acc_arr[8*(STOUT_STEPS-1)];
#else
    gconf_as_fermionmatrix = conf_acc;
#endif
    
    
    // STIRACCHIAMENTO DELL'APPROX RAZIONALE FIRST_INV
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      generate_vec3_soa_gauss(kloc_p);
      generate_vec3_soa_gauss(kloc_s);
      find_min_max_eigenvalue_soloopenacc(gconf_as_fermionmatrix,u1_back_field_phases,&(fermions_parameters[iflav]),kloc_r,kloc_h,kloc_p,kloc_s,minmaxeig);
      RationalApprox *approx_fi = &(fermions_parameters[iflav].approx_fi);
      RationalApprox *approx_fi_mother = &(fermions_parameters[iflav].approx_fi_mother);
      rescale_rational_approximation(approx_fi_mother,approx_fi,minmaxeig);
    }//end for iflav


    /////////////// INITIAL ACTION COMPUTATION ////////////////////////////////////////////
    double    action_ferm_in=0;
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
	int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;
	action_ferm_in += real_scal_prod_global(&ferm_phi_acc[ps_index],&ferm_phi_acc[ps_index]);
      }
    }// end for iflav
    ///////////////////////////////////////////////////////////////////////////////////////

    
    // FIRST INV APPROX CALC --> calcolo del fermione CHI
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
	int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;
	multishift_invert(gconf_as_fermionmatrix, &fermions_parameters[iflav], &(fermions_parameters[iflav].approx_fi), u1_back_field_phases, ferm_shiftmulti_acc, &(ferm_phi_acc[ps_index]), residue_metro, kloc_r, kloc_h, kloc_s, kloc_p, k_p_shiftferm);
        recombine_shifted_vec3_to_vec3(ferm_shiftmulti_acc, &(ferm_phi_acc[ps_index]), &(ferm_chi_acc[ps_index]),&(fermions_parameters[iflav].approx_fi));
      }
    }// end for iflav


    // STIRACCHIAMENTO DELL'APPROX RAZIONALE MD
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      minmaxeig[0] = (fermions_parameters[iflav].approx_fi.lambda_min);
      minmaxeig[1] = (fermions_parameters[iflav].approx_fi.lambda_max);
      minmaxeig[0] = minmaxeig[0] / minmaxeig[1];
      minmaxeig[0] = minmaxeig[0] / 0.95;
      minmaxeig[1] = minmaxeig[1] / 1.05;
      RationalApprox *approx_md = &(fermions_parameters[iflav].approx_md);
      RationalApprox *approx_md_mother = &(fermions_parameters[iflav].approx_md_mother);
      rescale_rational_approximation(approx_md_mother,approx_md,minmaxeig);
    }

    DEOTT_fermion_force_soloopenacc( conf_acc,
#ifdef STOUT_FERMIONS
			      gstout_conf_acc_arr, auxbis_conf_acc, // parkeggio
#endif
			      u1_back_field_phases, ipdot_acc, fermions_parameters, NDiffFlavs,
			      ferm_chi_acc, residue_metro, aux_conf_acc, ferm_shiftmulti_acc, kloc_r,
			      kloc_h, kloc_s, kloc_p, k_p_shiftferm);


    // CAMBIARE UN LINK!!!!!!!!!!
    //    conf_acc[0].r0.c0[0] += 0.001 + 0.0*I;
    double eps=0.00001;
    set_tamat_soa_to_zero(ipdot_acc);
    int i=0;
    int s=0;
    //    ipdot_acc[i].c01[s]= 0.0 + eps*I;
    ipdot_acc[i].c01[s]= 0.0 + eps*I;
    //    ipdot_acc[i].c02[s]= eps + 0.0*I;
    //    ipdot_acc[i].rc00[s]= eps/sqrt(3.0) ;
    //    ipdot_acc[i].rc11[s]= eps/sqrt(3.0) ;
    CH_exponential_antihermitian_soa_nissalike(&(aux_conf_acc[i]),&(ipdot_acc[i]),s);
    single_su3 M;
    M.comp[0][0] =aux_conf_acc[i].r0.c0[s]*conf_acc[i].r0.c0[s]+aux_conf_acc[i].r0.c1[s]*conf_acc[i].r1.c0[s]+aux_conf_acc[i].r0.c2[s]*conf_acc[i].r2.c0[s];
    M.comp[0][1] =aux_conf_acc[i].r0.c0[s]*conf_acc[i].r0.c1[s]+aux_conf_acc[i].r0.c1[s]*conf_acc[i].r1.c1[s]+aux_conf_acc[i].r0.c2[s]*conf_acc[i].r2.c1[s];
    M.comp[0][2] =aux_conf_acc[i].r0.c0[s]*conf_acc[i].r0.c2[s]+aux_conf_acc[i].r0.c1[s]*conf_acc[i].r1.c2[s]+aux_conf_acc[i].r0.c2[s]*conf_acc[i].r2.c2[s];
    M.comp[1][0] =aux_conf_acc[i].r1.c0[s]*conf_acc[i].r0.c0[s]+aux_conf_acc[i].r1.c1[s]*conf_acc[i].r1.c0[s]+aux_conf_acc[i].r1.c2[s]*conf_acc[i].r2.c0[s];
    M.comp[1][1] =aux_conf_acc[i].r1.c0[s]*conf_acc[i].r0.c1[s]+aux_conf_acc[i].r1.c1[s]*conf_acc[i].r1.c1[s]+aux_conf_acc[i].r1.c2[s]*conf_acc[i].r2.c1[s];
    M.comp[1][2] =aux_conf_acc[i].r1.c0[s]*conf_acc[i].r0.c2[s]+aux_conf_acc[i].r1.c1[s]*conf_acc[i].r1.c2[s]+aux_conf_acc[i].r1.c2[s]*conf_acc[i].r2.c2[s];
    M.comp[2][0] =aux_conf_acc[i].r2.c0[s]*conf_acc[i].r0.c0[s]+aux_conf_acc[i].r2.c1[s]*conf_acc[i].r1.c0[s]+aux_conf_acc[i].r2.c2[s]*conf_acc[i].r2.c0[s];
    M.comp[2][1] =aux_conf_acc[i].r2.c0[s]*conf_acc[i].r0.c1[s]+aux_conf_acc[i].r2.c1[s]*conf_acc[i].r1.c1[s]+aux_conf_acc[i].r2.c2[s]*conf_acc[i].r2.c1[s];
    M.comp[2][2] =aux_conf_acc[i].r2.c0[s]*conf_acc[i].r0.c2[s]+aux_conf_acc[i].r2.c1[s]*conf_acc[i].r1.c2[s]+aux_conf_acc[i].r2.c2[s]*conf_acc[i].r2.c2[s];

    conf_acc[i].r0.c0[s] = M.comp[0][0];
    conf_acc[i].r0.c1[s] = M.comp[0][1];
    conf_acc[i].r0.c2[s] = M.comp[0][2];
    conf_acc[i].r1.c0[s] = M.comp[1][0];
    conf_acc[i].r1.c1[s] = M.comp[1][1];
    conf_acc[i].r1.c2[s] = M.comp[1][2];
    conf_acc[i].r2.c0[s] = M.comp[2][0];
    conf_acc[i].r2.c1[s] = M.comp[2][1];
    conf_acc[i].r2.c2[s] = M.comp[2][2];

    print_su3_soa(conf_acc,"conf_mod");




#ifdef STOUT_FERMIONS
    // STOUTING...(ALREADY ON DEVICE)
    stout_wrapper(conf_acc,tstout_conf_acc_arr,0);
    gconf_as_fermionmatrix = &(tstout_conf_acc_arr[8*(STOUT_STEPS-1)]);
#else
    gconf_as_fermionmatrix = conf_acc;
#endif



    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      generate_vec3_soa_gauss(kloc_p);
      generate_vec3_soa_gauss(kloc_s);
      find_min_max_eigenvalue_soloopenacc(gconf_as_fermionmatrix,u1_back_field_phases,&(fermions_parameters[iflav]),kloc_r,kloc_h,kloc_p,kloc_s,minmaxeig);
      RationalApprox *approx_li = &(fermions_parameters[iflav].approx_li);
      RationalApprox *approx_li_mother = &(fermions_parameters[iflav].approx_li_mother);
      rescale_rational_approximation(approx_li_mother,approx_li,minmaxeig);
    }


    // LAST INV APPROX CALC 
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
        int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;
        // USING STOUTED CONF
        multishift_invert(gconf_as_fermionmatrix, &fermions_parameters[iflav], 
			  &(fermions_parameters[iflav].approx_li), u1_back_field_phases,
			  ferm_shiftmulti_acc, &(ferm_chi_acc[ps_index]), residue_metro, 
			  kloc_r, kloc_h, kloc_s, kloc_p, k_p_shiftferm);
        recombine_shifted_vec3_to_vec3(ferm_shiftmulti_acc, &(ferm_chi_acc[ps_index]), &(ferm_phi_acc[ps_index]),&(fermions_parameters[iflav].approx_li));
      }
    }

    ///////////////   FINAL ACTION COMPUTATION  ////////////////////////////////////////////
    double  action_ferm_fin=0;
    for(int iflav = 0 ; iflav < NDiffFlavs ; iflav++){
      for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
	int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;
	action_ferm_fin += real_scal_prod_global(&ferm_chi_acc[ps_index],&ferm_phi_acc[ps_index]);
      }
    } // end for iflav


    printf("ACTION_IN= %.18lf\n",action_ferm_in);    
    printf("ACTION_AF= %.18lf\n",action_ferm_fin);    
    printf("eps= %.18lf  DELTA_ACTION= %.18lf\n",eps,action_ferm_fin-action_ferm_in);    


  free(conf_acc);
  mem_free();

  return 0;
}


