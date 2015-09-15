double casuale(void);
void initrand(unsigned long s);
void su2_rand(double *pp);


#include "openacc.h"
#include "../RationalApprox/rationalapprox.c"
#include "./struct_c_def.c"
#include "./alloc_vars.c"
#include "./dbgtools.c"
#include "./fermionic_utilities.c"
#include "./su3_utilities.c"
#include "./random_assignement.c"
#include "./fermion_matrix.c"
#include "./inverter_full.c"
#include "./find_min_max.c"
#include "./inverter_multishift_full.c"
#include "./md_integrator.c"
#include "./update_standard_action.c"
#include "./rettangoli.c"
#include "../Meas/ferm_meas.c"


int main(){

  initrand(0);
  fflush(stdout);
  printf("INIZIO DEL PROGRAMMA \n");
  su3_soa  * conf_acc;
  int  allocation_check =  posix_memalign((void **)&conf_acc, ALIGN, 8*sizeof(su3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di conf_acc \n");
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

  //////  OPENACC CONTEXT INITIALIZATION    //////////////////////////////////////////////////////
  // NVIDIA GPUs
  acc_device_t my_device_type = acc_device_nvidia;
  // AMD GPUs
  // acc_device_t my_device_type = acc_device_radeon;
  // Intel XeonPhi
  //acc_device_t my_device_type = acc_device_xeonphi;
  // Select device ID
  int dev_index = 0;
  SELECT_INIT_ACC_DEVICE(my_device_type, dev_index);
  printf("Device Selected : OK \n");


  //###################### INIZIALIZZAZIONE DELLA CONFIGURAZIONE #################################
  // cold start
  if(start_opt==0){ 
    generate_Conf_cold(conf_acc,0.0);
    printf("Cold Gauge Conf Generated : OK \n");
    conf_id_iter=0;
  }
 // start from saved conf
  if(start_opt==1){
    read_su3_soa(conf_acc,"stored_config"); // READS ALSO THE conf_id_iter
    printf("Stored Gauge Conf Read : OK \n");
  }
  //###############################################################################################  


#pragma acc data   copy(conf_acc[0:8]) copyin(u1_back_field_phases[0:8]) create(ipdot_acc[0:8]) create(aux_conf_acc[0:8]) create(ferm_chi_acc[0:NPS_tot]) create(ferm_phi_acc[0:NPS_tot])  create(ferm_out_acc[0:NPS_tot]) create(ferm_shiftmulti_acc[0:max_ps*max_approx_order]) create(kloc_r[0:1])  create(kloc_h[0:1])  create(kloc_s[0:1])  create(kloc_p[0:1])  create(k_p_shiftferm[0:max_approx_order]) create(momenta[0:8]) copyin(nnp_openacc) copyin(nnm_openacc) create(local_sums[0:2]) create(d_local_sums[0:1])  copyin(fermions_parameters[0:NDiffFlavs])
    {

      double plq,rect;
    /*
    plq = calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
    printf("ALL'INIZIO   Placchetta=%.18lf \n",plq/size/6.0/3.0);
    */

    int accettate_therm=0;
    int accettate_metro=0;
    int id_iter_offset=conf_id_iter;
    //################### THERMALIZATION & METRO    ----   UPDATES ####################//
    for(int id_iter=id_iter_offset;id_iter<(ITERATIONS+id_iter_offset);id_iter++){
               conf_id_iter++;
	       printf("\n#################################################\n");
               printf(  "            GENERATING CONF %d of %d\n",conf_id_iter,ITERATIONS+id_iter_offset);
	       printf(  "#################################################\n\n");
	       //--------- CONF UPDATE ----------------//
	       if(id_iter<therm_ITERATIONS){
		 accettate_therm = UPDATE_SOLOACC_UNOSTEP_VERSATILE(conf_acc,residue_metro,residue_md,id_iter-id_iter_offset,accettate_therm,0);
	       }else{
		 accettate_metro = UPDATE_SOLOACC_UNOSTEP_VERSATILE(conf_acc,residue_metro,residue_md,id_iter-id_iter_offset-accettate_therm,accettate_metro,1);
	       }
#pragma acc update host(conf_acc[0:8])
	       //---------------------------------------//

	       //--------- MISURA ROBA FERMIONICA ----------------//
	       FILE *outfile = fopen(nome_file_ferm_output,"at");
	       if(!outfile) outfile = fopen(nome_file_ferm_output,"wt");
	       if(outfile){
		 fprintf(outfile,"%d\t",conf_id_iter);
		 for(int iflv=0;iflv<NDiffFlavs;iflv++) perform_chiral_measures(conf_acc,u1_back_field_phases,&(fermions_parameters[iflv]),residue_metro,outfile);
		 fprintf(outfile,"\n");
	       }
	       fclose(outfile);
	       //-------------------------------------------------// 
	       //--------- MISURA ROBA DI GAUGE ------------------//
	       plq = calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
	       rect = calc_rettangolo_soloopenacc(conf_acc,aux_conf_acc,local_sums);
	       if(id_iter<therm_ITERATIONS){
		 printf("Therm_iter %d   Placchetta=%.18lf \n",conf_id_iter,plq/size/6.0/3.0);
		 printf("Therm_iter %d   Rettangolo=%.18lf \n",conf_id_iter,rect/size/6.0/3.0/2.0);
	       }else{
		 printf("Metro_iter %d   Placchetta=%.18lf \n",conf_id_iter,plq/size/6.0/3.0);
		 printf("Metro_iter %d   Rettangolo=%.18lf \n",conf_id_iter,rect/size/6.0/3.0/2.0);
	       }
	       //-------------------------------------------------//
	       
	       //--------- SALVA LA CONF SU FILE ------------------//
	       if(conf_id_iter%save_conf_every==0)	print_su3_soa(conf_acc,"stored_config");
	       //-------------------------------------------------//
	       
    }// id_iter loop ends here
    

    //--------- SALVA LA CONF SU FILE ------------------//
    print_su3_soa(conf_acc,"stored_config");
    //-------------------------------------------------//


    }// end pragma acc data


  //////  OPENACC CONTEXT CLOSING    //////////////////////////////////////////////////////////////
  SHUTDOWN_ACC_DEVICE(my_device_type);
  /////////////////////////////////////////////////////////////////////////////////////////////////


  free(conf_acc);
  mem_free();


  return 0;
}

