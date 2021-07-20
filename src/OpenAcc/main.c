//here macros are defined
#define PRINT_DETAILS_INSIDE_UPDATE
#define ALIGN 128



// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
#define __restrict
#endif

#ifndef __GNUC__
#include "openacc.h"
#endif

#ifdef ONE_FILE_COMPILATION
#include "../Include/all_include.h"
#endif

//######################################################################################################################################//
//############### INCLUTIONS #########################################################################################################//
//######################################################################################################################################//


#include "../DbgTools/dbgtools.h" // DEBUG
#include "../Include/debug.h"
#include "../Include/fermion_parameters.h"
#include "../Include/montecarlo_parameters.h"
#include "../Include/inverter_tricks.h"
#include "../Include/memory_wrapper.h"
#include "../Include/setting_file_parser.h"
#include "../Include/tell_geom_defines.h"
#include "../Include/rep_info.h"
#include "../Meas/ferm_meas.h"
#include "../Meas/gauge_meas.h"
#include "../Meas/polyakov.h"
#include "../Meas/measure_topo.h"
#include "../Mpi/communications.h"
#include "../Mpi/multidev.h"
#include "../Rand/random.h"
#include "../RationalApprox/rationalapprox.h"
#include "./action.h"
#include "./alloc_vars.h"
#include "./backfield_parameters.h"
#include "./deviceinit.h"
#include "./fermion_matrix.h"
#include "./fermionic_utilities.h"
#include "./find_min_max.h"
#include "./inverter_full.h"
#include "./inverter_multishift_full.h"
#include "./io.h"
#include "./ipdot_gauge.h"
#include "./md_integrator.h"
#include "./md_parameters.h"
#include "./random_assignement.h"
#include "./rettangoli.h"
#include "./sp_alloc_vars.h"
#include "./stouting.h"
#include "./struct_c_def.h"
#include "./su3_measurements.h"
#include "./su3_utilities.h"
#include "./update_versatile.h"
#include "./alloc_settings.h"
#ifdef __GNUC__
#include "sys/time.h"
#endif
#include "./cooling.h"
#include <unistd.h>
#include <mpi.h>
#include "../Include/stringify.h"

#include <errno.h>
#include "./Mod_functions.h" //funzioni moddate.

//######################################################################################################################################//
//######################################################################################################################################//


// double level macro, necessary to stringify
// https://gcc.gnu.org/onlinedocs/cpp/Stringification.html
#define xstr(s) str(s) 
#define str(s) #s


//definitions outside the main.

int conf_id_iter;
int verbosity_lv;

extern int TOPO_GLOBAL_DONT_TOUCH; 





int main(int argc, char* argv[]){

  
    
    gettimeofday ( &(mc_params.start_time), NULL );
    FILE *file_label;
     file_label=fopen("./file_label.txt","w");
//######################################################################################################################################//
//############### FILE READING #########################################################################################################//
//######################################################################################################################################//

    // READ input file.
#ifdef MULTIDEVICE
    pre_init_multidev1D(&devinfo);
    {
        volatile int flag = 0;
        if(0==devinfo.myrank && getenv("GDBHOOK")){
            printf("MPI%02d: Thank you for using the GDB HOOK!\n",devinfo.myrank);
            printf("         Waiting for user intervention.Please attach to process %d,\n",
                    getpid());
            printf("         [hint: gdb -p %d ] \n",  getpid());
            printf("         and set the value of 'flag' to 1.\n");
            printf("         You may have to use 'finish' a couple of times to go down the call stack.\n");
            printf("         HAPPY HUNTING AND GOOD LUCK!!!!\n");
            while(1 != flag)
                sleep(1);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    if(0==devinfo.myrank){
        printf("****************************************************\n");
        printf("          PRE INIT - READING SETTING  FILE          \n");
        printf("     check which parameter corresponds to what! \n");
        printf("commit: %s\n", xstr(COMMIT_HASH) );
        printf("****************************************************\n");
    }


    int input_file_read_check = set_global_vars_and_fermions_from_input_file(argv[1]); //here the reading happens.
    


   
    printf("ECCO IL NUMERO DI REPLICHE %d  %d!!!\n",rep->replicas_total_number,alloc_info.num_replicas); // further controls.
    int i3=0;
    
    printf("il defect è sulla linea %d\n",rep->defect_boundary);
    
    printf("Ecco i valori del defect:\n"); //defect values.
    for(i3=0;i3<3;i3++){
        
        printf("def%d: %d\n",i3,rep->defect_coordinates[i3]);
    }
    
    
    printf("Ecco i valori dei cr:\n");
    for(i3=0;i3<alloc_info.num_replicas;i3++){
        
        printf("c%d: %f\n",i3,rep->cr_vet[i3]);
    }
    
    
    /*counter_size_function(nd0,nd1,nd2,nd3);*/
    

    
    
    
#ifdef MULTIDEVICE
    if(input_file_read_check){
        printf("MPI%02d: input file reading failed, Aborting...\n",devinfo.myrank);
        MPI_Abort(MPI_COMM_WORLD,1);
    }else init_multidev1D(&devinfo);
#else
    devinfo.myrank = 0;
    devinfo.nranks = 1;
#endif

    if(input_file_read_check){
        printf("MPI%02d: input file reading failed, aborting...\n",devinfo.myrank);
        exit(1);
    }
    if(0==devinfo.myrank) print_geom_defines();
    verbosity_lv = debug_settings.input_vbl;

    //######################################################################################################################################//
    //######################################################################################################################################//
    
    
    //######################################################################################################################################//
    //############### INITIALIZATION #########################################################################################################//
    //######################################################################################################################################//
  
    
    //just printing headtitles
    if(0==devinfo.myrank){
        if(0 != mc_params.JarzynskiMode){
            printf("****************************************************\n");
            printf("                   JARZYNSKI MODE              \n");
            printf("     check which parameter corresponds to what! \n");
            printf("****************************************************\n");


        }else {
            printf("****************************************************\n");
            printf("                    NORMAL MODE                \n");
            printf("****************************************************\n");
        }
        if(debug_settings.do_norandom_test){
            printf("****************************************************\n");
            printf("         WELCOME. This is a NORANDOM test.    \n");
            printf("     MOST things will not be random generated,\n");
            printf("            but read from memory instead.     \n");
            printf("                  CHECK THE CODE!!            \n");
            printf("   ALSO: setting the number of trajectories to 1.\n");
            printf("****************************************************\n");
            mc_params.ntraj = 1;

        }
    }




    if(verbosity_lv > 2) 
        printf("MPI%02d, Input file read and initialized multidev1D...\n",
                devinfo.myrank);



#ifndef __GNUC__
    //////  OPENACC CONTEXT INITIALIZATION    //////////////////////////////////////////////////////
    // NVIDIA GPUs
    acc_device_t my_device_type = acc_device_nvidia;
    // AMD GPUs
    // acc_device_t my_device_type = acc_device_radeon;
    // Intel XeonPhi
    //acc_device_t my_device_type = acc_device_xeonphi;
    // Select device ID
    printf("MPI%02d: Selecting device.\n", devinfo.myrank);
#ifdef MULTIDEVICE
    select_init_acc_device(my_device_type, (devinfo.single_dev_choice + devinfo.myrank)%devinfo.proc_per_node);
#else
    select_init_acc_device(my_device_type, devinfo.single_dev_choice);
#endif
    printf("Device Selected : OK \n"); //checking printing.
#endif

    unsigned int myseed_default =  (unsigned int) mc_params.seed; 

#ifdef MULTIDEVICE
    myseed_default =  (unsigned int) (myseed_default + devinfo.myrank) ;
    char myrank_string[6];
    sprintf(myrank_string,".R%d",devinfo.myrank);
    strcat(mc_params.RandGenStatusFilename,myrank_string);
#endif




    initrand_fromfile(mc_params.RandGenStatusFilename,myseed_default);


    // INIT FERM PARAMS AND READ RATIONAL APPROX COEFFS
    if(init_ferm_params(fermions_parameters)){
        printf("MPI%02d - Finalizing...\n",devinfo.myrank); //cp
#ifdef MULTIDEVICE
        MPI_Finalize();
#endif
        exit(1);
    }
#pragma acc enter data copyin(fermions_parameters[0:alloc_info.NDiffFlavs])
    
//********************* ALLOCATION OF CONF *******************************************//
    
    mem_alloc_core(); // Allocation has been done here.
    mem_alloc_extended(); //extend alloc.
 
     //************NOW HASENBUSCH CONF EXISTS!!*****************//
    
    //////direi che il for per la conf va messo qui./////
    
    //eventual single precision allocation
    printf("\n   MPI%02d - Allocazione della memoria (double) : OK \n\n\n",devinfo.myrank);
    if(inverter_tricks.useMixedPrecision || md_parameters.singlePrecMD){
        mem_alloc_core_f();
        printf("\n  MPI%02d - Allocazione della memoria (float) [CORE]: OK \n\n\n",devinfo.myrank);
    }

    if( md_parameters.singlePrecMD){
        mem_alloc_extended_f();
        printf("\n  MPI%02d - Allocazione della memoria (float) [EXTENDED]: OK \n\n\n",devinfo.myrank);
    }
   
    
    //possible conf for cycle position //
    
     //volendo anche dopo la definzione della tavola dei nn, il ciclo, non serve riderfiniral ogni volta.
       //direi che il for per la conf va messo qui.///// da qui in avanti -->conf_hasenbusch[i].
    
    printf("\n  MPI%02d - Allocazione della memoria totale: %zu \n\n\n",devinfo.myrank,max_memory_used);
 
    
   //*****************NEAREST NEIGHBOURS DEFINITION!!******************************//
    
    compute_nnp_and_nnm_openacc(); //nn def. Here 2 matrices nn(p/m)_openacc[sizeh][4][2] are used to contain nn information. The first index refers to the half-lattice site, the second to the direction, the last to the parity(even or odd half-lattice.).
#pragma acc enter data copyin(nnp_openacc) //copied on GPU/CPU.
#pragma acc enter data copyin(nnm_openacc)
    printf("MPI%02d - nn computation : OK \n",devinfo.myrank);
    init_all_u1_phases(backfield_parameters,fermions_parameters);
#pragma acc update device(u1_back_phases[0:8*alloc_info.NDiffFlavs])
#pragma acc update device(mag_obs_re[0:8*alloc_info.NDiffFlavs])
#pragma acc update device(mag_obs_im[0:8*alloc_info.NDiffFlavs])

    if(inverter_tricks.useMixedPrecision || md_parameters.singlePrecMD){
#pragma acc update device(u1_back_phases_f[0:8*alloc_info.NDiffFlavs])
    }

    printf("MPI%02d - u1_backfield initialization (float & double): OK \n",devinfo.myrank);

    initialize_md_global_variables(md_parameters);//md_paramaters init.
    printf("MPI%02d - init md vars : OK \n",devinfo.myrank);

//anche da qui il ciclo.

    //***************** INIZIALIZZAZIONE DELLA CONFIGURAZIONE *************//
    
    
    
    int replicas_counter;
    char rep_str [10];
    char aux_name_file[200];
    strcpy(aux_name_file,mc_params.save_conf_name); //salva il nome originale della conf in aux_name_file
    //init for starting
    
    
    //for per aprire o creare i file delle configurazioni.
    for(replicas_counter=0;replicas_counter<rep->replicas_total_number;replicas_counter++){
        printf("replicas counter %d\n",replicas_counter);
        snprintf(rep_str,10,"replica_%d",replicas_counter);//inizializza rep_str
        strcat(mc_params.save_conf_name,rep_str); //appiccica rep_str in fondo.
        
        
     
    if(debug_settings.do_norandom_test){
        if(!read_conf_wrapper(conf_hasenbusch[replicas_counter],"conf_norndtest",&conf_id_iter,debug_settings.use_ildg)){
            // READS ALSO THE conf_id_iter
            printf("MPI%02d - Stored Gauge Conf conf_norndtest Read : OK\n",devinfo.myrank);
        }
        else{
            // cold start
            printf("MPI%02d - COMPILED IN NORANDOM MODE. A CONFIGURATION FILE NAMED\
                    \"conf_norndtest\" MUST BE PRESENT\n",devinfo.myrank);
            exit(1);
        }
    }
    else{
        if(!read_conf_wrapper(conf_hasenbusch[replicas_counter],mc_params.save_conf_name,
                    &conf_id_iter,debug_settings.use_ildg)){
            // READS ALSO THE conf_id_iter
            printf("MPI%02d - Stored Gauge Conf \"%s\" Read : OK \n",
                    devinfo.myrank, mc_params.save_conf_name);

        }
        else{
            
            
            
            generate_Conf_cold(conf_hasenbusch[replicas_counter],mc_params.eps_gen);
            printf("MPI%02d - Cold Gauge Conf Generated : OK \n",
                    devinfo.myrank);
            conf_id_iter=0;
        }
    }
        
    
        strcpy(mc_params.save_conf_name,aux_name_file);
        
        
    }//for  replicas closing
    
       //*****************inizializzazione Ku della conf****************//
  
    //Here. I include the K_mu inizialization for each conference.
    

          printf("K_mu Initialization starts\n");
    
    for(replicas_counter=0;replicas_counter<rep->replicas_total_number;replicas_counter++){
        
     printf("Initializing K_mu & label Replica %d\n",replicas_counter);
        
        int mu1;
        for (mu1=0;mu1<8;mu1++){
            conf_hasenbusch[replicas_counter][mu1].label=replicas_counter;
        }
    int init_result=init_k(conf_hasenbusch[replicas_counter],rep->cr_vet[replicas_counter],rep->defect_boundary,rep->defect_coordinates);
   #pragma acc update device(conf_hasenbusch[0:rep->replicas_total_number][0:8])
        
    if(init_result!=0){printf("Error in Initialization Replica %d!\n",replicas_counter); return 1;};
    
    }
    
        printf(" Initialization K_mu success\n");
   
    


    
    
    
    //#################################################################################  


//Here again, the job done for 1 replica simulations it has been turned into a multi conf process.
    double max_unitarity_deviation,avg_unitarity_deviation;
    
     for(replicas_counter=0;replicas_counter<rep->replicas_total_number;replicas_counter++){ //for start
         
    check_unitarity_host(conf_hasenbusch[replicas_counter],&max_unitarity_deviation,&avg_unitarity_deviation);
    printf("\tMPI%02d: Avg_unitarity_deviation on host: %e\n", devinfo.myrank, 
            avg_unitarity_deviation);
    printf("\tMPI%02d: Max_unitarity_deviation on host: %e\n", devinfo.myrank,
            max_unitarity_deviation);

     }// for end

 
    //######################################################################################################################################//
    //######################################################################################################################################//

    
    //TEST SWAP!!
    double Delta_S_SWAP_0=0.0;
     double Delta_S_SWAP_1=0.0;
      double Delta_S_SWAP_2=0.0;
    double S_0_0, S_0_2,S_2_2,S_2_0;
    double S_0_0_RET, S_0_2_RET,S_2_2_RET,S_2_0_RET;
    double Delta_S_SWAPS=0.0;
    printf("ECCO IL TEST SWAP!!!\n");
    
    #pragma acc update device(conf_hasenbusch[0:rep->replicas_total_number][0:8])
   // printing_k_mu(conf_hasenbusch[2]);
    S_0_0=BETA_BY_THREE*calc_plaquette_soloopenacc(conf_hasenbusch[0],aux_conf_acc,local_sums);
    
    S_0_0_RET=C_ZERO*S_0_0+C_ONE*BETA_BY_THREE*calc_rettangolo_soloopenacc(conf_hasenbusch[0],aux_conf_acc,local_sums);
    
    
    
    
    S_2_2=BETA_BY_THREE*calc_plaquette_soloopenacc(conf_hasenbusch[2],aux_conf_acc,local_sums);
    
    S_2_2_RET=C_ZERO*S_2_2+C_ONE*BETA_BY_THREE*calc_rettangolo_soloopenacc(conf_hasenbusch[2],aux_conf_acc,local_sums);
    
    
    
    
   
    
    
replicas_swap(conf_hasenbusch[0],conf_hasenbusch[2],rep->defect_boundary,rep->defect_coordinates);
#pragma acc update device(conf_hasenbusch[0:rep->replicas_total_number][0:8])
    
    S_0_2=BETA_BY_THREE*calc_plaquette_soloopenacc(conf_hasenbusch[0],aux_conf_acc,local_sums);
    S_0_2_RET=C_ZERO*S_0_2+C_ONE*BETA_BY_THREE*calc_rettangolo_soloopenacc(conf_hasenbusch[0],aux_conf_acc,local_sums);
    
    
    S_2_0=BETA_BY_THREE*calc_plaquette_soloopenacc(conf_hasenbusch[2],aux_conf_acc,local_sums);
    S_2_0_RET=C_ZERO*S_2_0+C_ONE*BETA_BY_THREE*calc_rettangolo_soloopenacc(conf_hasenbusch[2],aux_conf_acc,local_sums);
    
    Delta_S_SWAPS=(S_2_0+S_0_2)-(S_0_0+S_2_2);
    
    Delta_S_SWAPS=C_ZERO*Delta_S_SWAPS;
 
    printf("S_2_0 %f  S_0_2 %f  S_0_0 %f  S_2_2 %f\n",S_2_0,S_0_2,S_0_0,S_2_2);
   printf("RET S_2_0 %f  S_0_2 %f  S_0_0 %f  S_2_2 %f\n",S_2_0_RET,S_0_2_RET,S_0_0_RET,S_2_2_RET);

   
    Delta_S_SWAP_0=calc_plaquette_soloopenacc_SWAP(conf_hasenbusch[0],conf_hasenbusch[2],aux_conf_acc,local_sums,rep->defect_boundary,rep->defect_coordinates,0);
    int vet_prova[3]={nd1,nd2,nd3};
    

    int vet_prova_sub[3]={nd1/2,nd2/2,nd3/2};
    
  //  Delta_S_SWAP_1=calc_plaquette_soloopenacc_SWAP(conf_hasenbusch[0],conf_hasenbusch[2],aux_conf_acc,local_sums,rep->defect_boundary,vet_prova,0);

    
    Delta_S_SWAP_2=calc_plaquette_soloopenacc_SWAP(conf_hasenbusch[0],conf_hasenbusch[2],aux_conf_acc,local_sums,rep->defect_boundary,vet_prova_sub,0);

    printf("CONFRONTO DELTA_SWAP\n");
    printf("%lf (S_2_0+S_0_2)-(S_0_0+S_2_2)  ||    %lf (DELTA_S)  \n",Delta_S_SWAPS,Delta_S_SWAP_0);
    
    printf("DELTA_S_SWAP: %f(optimus)|| %f(casual 8th) || %f(full)\n",Delta_S_SWAP_0,Delta_S_SWAP_2,Delta_S_SWAP_1);

    Delta_S_SWAPS=(S_2_0_RET+S_0_2_RET)-(S_0_0_RET+S_2_2_RET);
    Delta_S_SWAP_0=calc_plaquette_soloopenacc_SWAP(conf_hasenbusch[0],conf_hasenbusch[2],aux_conf_acc,local_sums,rep->defect_boundary,rep->defect_coordinates,1);
    
    printf("CONFRONTO DELTA_SWAP 2\n");
    printf("%lf RET (S_2_0+S_0_2)-(S_0_0+S_2_2)  ||    %lf (DELTA_S)  \n",Delta_S_SWAPS,Delta_S_SWAP_0);
    
    

    printf("DELTA_S_SWAP: %f(optimus)|| %f(casual 8th) || %f(full)\n",Delta_S_SWAP_0,Delta_S_SWAP_2,Delta_S_SWAP_1);
    
    //NEW TEST SWAP:
    
    int  number_accept; //acceptance number
    int swap_number=0;
   
    label_print(conf_hasenbusch, rep->replicas_total_number,file_label,swap_number);

    //acceptance vector
    int all_swap_vector[rep->replicas_total_number];
    int acceptance_vector[rep->replicas_total_number];
    int mu1,mu2;
    //intialization
    for(mu1=0;mu1<rep->replicas_total_number;mu1++){
        all_swap_vector[mu1]=0;
        acceptance_vector[mu1]=0;
    }
    
    
    for(mu2=0;mu2<5;mu2++){
    
    for(mu1=0;mu1<8;mu1++){
     printf("beforrre (%d) %.18lf %.18lf\n",mu1,creal(conf_hasenbusch[0][mu1].r1.c1[snum_acc(31,6,6,6)]),creal(conf_hasenbusch[2][mu1].r1.c1[snum_acc(31,6,6,6)]));
    }
    
    
    int accettata=0;
        
        accettata=metro_SWAP( conf_hasenbusch, 0, 2,rep->defect_boundary,rep->defect_coordinates);
        all_swap_vector[0]++;
        all_swap_vector[2]++;
        
        if(accettata==1){
            acceptance_vector[0]++;
            acceptance_vector[2]++;
        }
    #pragma acc update device(conf_hasenbusch[0:rep->replicas_total_number][0:8])
        swap_number++;
        label_print(conf_hasenbusch, rep->replicas_total_number,file_label,swap_number);


    printf("acpt :%d\n",accettata);
        double mean_accept;
        mean_accept=acceptance_vector[0]/all_swap_vector[0];
        printf("mean accept:%f\n",mean_accept);
    
     for(mu1=0;mu1<8;mu1++){
     printf("aftermath (%d) %.18lf %.18lf\n",mu1,creal(conf_hasenbusch[0][mu1].r1.c1[snum_acc(31,6,6,6)]),creal(conf_hasenbusch[2][mu1].r1.c1[snum_acc(31,6,6,6)]));
     }
    
    
    }

    
    
    //######################################################################################################################################//
    //############### MEASURES #########################################################################################################//
    //######################################################################################################################################//
    
    
    
    double plq,rect;
    double cool_topo_ch[meastopo_params.coolmeasstep/meastopo_params.cool_measinterval+1];
    double stout_topo_ch[meastopo_params.stoutmeasstep/meastopo_params.stout_measinterval+1];
    d_complex poly;

    int accettate_therm=0;
    int accettate_metro=0;
    int accettate_therm_old=0;
    int accettate_metro_old=0;
    int id_iter_offset=conf_id_iter;
    
   
    //Plaquette measures and polyakov loop measures.
    printf("PLAQUETTE START\n");
    //TEST MOD
    printf("Plaquette  MOD   : %.18lf\n" , S_0_2/GL_SIZE/BETA_BY_THREE/6.0/3.0);
    

    
    
    plq = calc_plaquette_soloopenacc(conf_hasenbusch[0],aux_conf_acc,local_sums);
    printf("\tMPI%02d: Therm_iter %d Placchetta    = %.18lf \n",
            devinfo.myrank, conf_id_iter,plq/GL_SIZE/6.0/3.0);
    
    printf("PLAQUETTE END\n");
    
    rect = calc_rettangolo_soloopenacc(conf_hasenbusch[0],aux_conf_acc,local_sums);
    printf("\tMPI%02d: Therm_iter %d Rettangolo    = %.18lf \n",
            devinfo.myrank, conf_id_iter,rect/GL_SIZE/6.0/3.0/2.0);

    poly =  (*polyakov_loop[geom_par.tmap])(conf_hasenbusch[0]);//misura polyakov loop
    printf("\tMPI%02d: Therm_iter %d Polyakov Loop = (%.18lf, %.18lf)  \n",
            devinfo.myrank, conf_id_iter,creal(poly),cimag(poly));

    
    //char confile_dbg[50];
    //sprintf(confile_dbg,"conf_ascii_test_%s", devinfo.myrankstr);
    //dbg_print_su3_soa(conf_hasenbusch[0],confile_dbg,0);


    //            MPI_Finalize(); // DEBUG
    //            return 0 ;      // DEBUG
    
    

//Here we are in Jarzynski mode.
    
    
    if(0 == mc_params.ntraj && 0 == mc_params.JarzynskiMode ){ // MEASURES ONLY
      
        printf("\n#################################################\n");
        printf("\tMEASUREMENTS ONLY ON FILE %s\n", mc_params.save_conf_name);
        printf("\n#################################################\n");

        //-------------------------------------------------// 
        //--------- MISURA ROBA DI GAUGE ------------------//
        if(0 == devinfo.myrank ) printf("Misure di Gauge:\n");
        plq = calc_plaquette_soloopenacc(conf_hasenbusch[0],aux_conf_acc,local_sums);
        rect = calc_rettangolo_soloopenacc(conf_hasenbusch[0],aux_conf_acc,local_sums);
        poly =  (*polyakov_loop[geom_par.tmap])(conf_hasenbusch[0]);//misura polyakov loop
    

        
        
	    printf("Plaquette     : %.18lf\n" ,plq/GL_SIZE/3.0/6.0);
        printf("Rectangle     : %.18lf\n" ,rect/GL_SIZE/3.0/6.0/2.0);
        printf("Polyakov Loop : (%.18lf,%.18lf) \n",creal(poly),cimag(poly));

        //--------- MISURA ROBA FERMIONICA ----------------//
        //
        if(0 == devinfo.myrank)  printf("Fermion Measurements: see file %s\n",
                fm_par.fermionic_outfilename);
        fermion_measures(conf_hasenbusch[0],fermions_parameters,
                &fm_par, md_parameters.residue_metro, 
                md_parameters.max_cg_iterations, id_iter_offset,
                plq/GL_SIZE/3.0/6.0,
                rect/GL_SIZE/3.0/6.0/2.0);   


      


     }else printf("MPI%02d: Starting generation of Configurations.\n",
            devinfo.myrank);
    

    
    
    
    
    
    
     
    // THERMALIZATION & METRO    ----   UPDATES //

    int id_iter=id_iter_offset;
    
    init_global_program_status(); 

    printf("run_condition: %d\n",mc_params.run_condition) ;
  
    while ( RUN_CONDITION_TERMINATE != mc_params.run_condition)
    {

        if(GPSTATUS_UPDATE == mc_params.next_gps){
            struct timeval tstart_cycle,tend_cycle;
            gettimeofday(&tstart_cycle, NULL);

            if(0 != mc_params.JarzynskiMode ){

                bf_param new_backfield_parameters = backfield_parameters;

                // DIRECT MODE 
                if(1 == mc_params.JarzynskiMode)
                    new_backfield_parameters.bz = backfield_parameters.bz + 
                        (double) id_iter/mc_params.MaxConfIdIter;
                // REVERSE MODE
                if(-1 == mc_params.JarzynskiMode)
                    new_backfield_parameters.bz = backfield_parameters.bz -
                        (double) id_iter/mc_params.MaxConfIdIter;




                if(0==devinfo.myrank){

                    if(1 == mc_params.JarzynskiMode)
                        printf("\n\nJarzynskiMode - DIRECT - From bz=%f to bz=%f+1 in %d steps.\n",
                                backfield_parameters.bz , backfield_parameters.bz, 
                                mc_params.MaxConfIdIter);
                    if(-1 == mc_params.JarzynskiMode)
                        printf("\n\nJarzynskiMode - REVERSE - From bz=%f to bz=%f-1 in %d steps.\n",
                                backfield_parameters.bz , backfield_parameters.bz, 
                                mc_params.MaxConfIdIter);


                    printf("JarzynskiMode, iteration %d/%d (%d max for this run)\n",
                            id_iter,mc_params.MaxConfIdIter,mc_params.ntraj);
                    printf("JarzynskiMode - current bz value : %f\n", new_backfield_parameters.bz);
                }

                init_all_u1_phases(new_backfield_parameters,fermions_parameters);
#pragma acc update device(u1_back_phases[0:8*alloc_info.NDiffFlavs])
#pragma acc update device(u1_back_phases_f[0:8*alloc_info.NDiffFlavs])
                printf("Jarzynski mode's end\n");
            }


            check_unitarity_device(conf_hasenbusch[0],&max_unitarity_deviation,
                    &avg_unitarity_deviation);
            printf("\tMPI%02d: Avg/Max unitarity deviation on device: %e / %e\n", 
                    devinfo.myrank,avg_unitarity_deviation,max_unitarity_deviation);
            accettate_therm_old = accettate_therm;
            accettate_metro_old = accettate_metro;


            if(devinfo.myrank ==0 ){
                printf("\n#################################################\n");
                printf(  "   GENERATING CONF %d of %d, %dx%dx%dx%d,%1.3f \n",
                        conf_id_iter,mc_params.ntraj+id_iter_offset,
                        geom_par.gnx,geom_par.gny,
                        geom_par.gnz,geom_par.gnt,
                        act_params.beta);
                printf(  "#################################################\n\n");
            }
            
            
            

            //--------- CONF UPDATE ----------------//
        
            //If metro, the second-last parameter, is 0 the update function thermalizes the conf, otherwise( metro=1) it does metropolis.
            
            for(replicas_counter=0;replicas_counter<rep->replicas_total_number;replicas_counter++){
                 printf("\n#################################################\n");
                printf("REPLICA %d:\n",replicas_counter);
                printf(  "#################################################\n\n");
                
            if(id_iter<mc_params.therm_ntraj){
                accettate_therm = UPDATE_SOLOACC_UNOSTEP_VERSATILE(conf_hasenbusch[replicas_counter],
                        md_parameters.residue_metro,md_parameters.residue_md,
                        id_iter-id_iter_offset,
                        accettate_therm,0,md_parameters.max_cg_iterations);
            }else{
                accettate_metro = UPDATE_SOLOACC_UNOSTEP_VERSATILE(conf_hasenbusch[replicas_counter],
                        md_parameters.residue_metro,md_parameters.residue_md,
                        id_iter-id_iter_offset-accettate_therm,accettate_metro,1,
                        md_parameters.max_cg_iterations);
                if(0==devinfo.myrank){ //multidevice control
                    int iterations = id_iter-id_iter_offset-accettate_therm +1;
                    double acceptance = (double) accettate_metro / iterations;
                    double acc_err = 
                        sqrt((double)accettate_metro*(iterations-accettate_metro)/iterations)
                        /iterations;
                    printf("Estimated acceptance for this run: %f +- %f\n",acceptance,
                            acc_err);
                    
                }
                //QUI VA IL CONF SWAP
                printf("CONF SWAP HERE!\n ");
                All_Conf_SWAP(conf_hasenbusch, rep->replicas_total_number,rep->defect_boundary, rep->defect_coordinates,file_label, swap_number);
                printf("swap_number %d", swap_number);
            }
            }
#pragma acc update self(conf_hasenbusch[0:rep->replicas_total_number][0:8]) //updating conf sul device
             //-----------------------------------------------//
             //---------------CONF SWAP---------------------------//
            
            //FILE *file_label;
           // file_label=fopen("./file_label.txt","w");
            //int swap_number=0;
            double Delta_S_SWAP=0.0;
            //label_print(conf_hasenbusch, rep->replicas_total_number,file_label,swap_number);
            printf("ECCO LO SWAP\n");
            
            
            
            Delta_S_SWAP=calc_plaquette_soloopenacc_SWAP(conf_hasenbusch[0],conf_hasenbusch[2],aux_conf_acc,local_sums,rep->defect_boundary,rep->defect_coordinates,0);
            printf("DELTA_S_SWAP: %f\n",Delta_S_SWAP);
            
            
            replicas_swap(conf_hasenbusch[0],conf_hasenbusch[2],rep->defect_boundary,rep->defect_coordinates);
            //swap_number++;
            
           // label_print(conf_hasenbusch, rep->replicas_total_number,file_label,swap_number);
            
        //-----------------------------------------------//
            
            id_iter++;
            conf_id_iter++;
            //-------------------------------------------------// 

            //--------- MISURA ROBA DI GAUGE ------------------//
            plq  = calc_plaquette_soloopenacc(conf_hasenbusch[0],aux_conf_acc,local_sums);
            rect = calc_rettangolo_soloopenacc(conf_hasenbusch[0],aux_conf_acc,local_sums);
            poly =  (*polyakov_loop[geom_par.tmap])(conf_hasenbusch[0]);
	    
	    if(meastopo_params.meascool && conf_id_iter%meastopo_params.cooleach==0){
	      su3_soa *conf_to_use;
	      cool_topo_ch[0]=compute_topological_charge(conf_hasenbusch[0],auxbis_conf_acc,topo_loc);
		    for(int cs = 1; cs <= meastopo_params.coolmeasstep; cs++){
		      if(cs==1)
			conf_to_use=(su3_soa*)conf_hasenbusch[0];
		      else
			conf_to_use=(su3_soa*)aux_conf_acc;
		      cool_conf(conf_to_use,aux_conf_acc,auxbis_conf_acc);
		      if(cs%meastopo_params.cool_measinterval==0)
			cool_topo_ch[cs/meastopo_params.cool_measinterval]=compute_topological_charge(aux_conf_acc,auxbis_conf_acc,topo_loc);
		    }
	            printf("MPI%02d - Printing cooled charge - only by master rank...\n",
                    devinfo.myrank);
        	    if(devinfo.myrank ==0){
                	FILE *cooloutfile = fopen(meastopo_params.pathcool,"at");
	                if(!cooloutfile){
        	            cooloutfile = fopen(meastopo_params.pathcool,"wt");
	                    char coolheader[35];
			    strcpy(coolheader,"#conf_id\tCoolStp\tTopoChCool\n");
                	    fprintf(cooloutfile,"%s",coolheader);
	                }
        	        if(cooloutfile){
				for(int i = 0; i <= meastopo_params.coolmeasstep/meastopo_params.cool_measinterval;i++)
					fprintf(cooloutfile,"%d\t%d\t%18.18lf\n",conf_id_iter,
						i*meastopo_params.cool_measinterval,
						cool_topo_ch[i]);
			}
        	        fclose(cooloutfile);
            	    }
	    }
	    if(meastopo_params.measstout && conf_id_iter%meastopo_params.stouteach==0){
		    TOPO_GLOBAL_DONT_TOUCH=1;
		    stout_wrapper(conf_hasenbusch[0],gstout_conf_acc_arr);
		    TOPO_GLOBAL_DONT_TOUCH=0;
		    stout_topo_ch[0]=compute_topological_charge(conf_hasenbusch[0],auxbis_conf_acc,topo_loc);
		    for(int ss = 0; ss < meastopo_params.stoutmeasstep; ss+=meastopo_params.stout_measinterval){
	    	    	/* aux_conf_acc = &(gstout_conf_acc_arr[8*ss]); */
			int topoindx =1+ss/meastopo_params.stout_measinterval; 
		    	stout_topo_ch[topoindx]=compute_topological_charge(&gstout_conf_acc_arr[8*ss],auxbis_conf_acc,topo_loc);
		    }
	    
	            printf("MPI%02d - Printing stouted charge - only by master rank...\n",
                    devinfo.myrank);
        	    if(devinfo.myrank ==0){
                	FILE *stoutoutfile = fopen(meastopo_params.pathstout,"at");
	                if(!stoutoutfile){
        	            stoutoutfile = fopen(meastopo_params.pathstout,"wt");
	                    char stoutheader[35];
			    strcpy(stoutheader,"#conf_id\tStoutStp\tTopoChStout\n");
                	    fprintf(stoutoutfile,"%s",stoutheader);
	                }
        	        if(stoutoutfile){
				for(int i = 0; i <= meastopo_params.stoutmeasstep/meastopo_params.stout_measinterval;i++)
					fprintf(stoutoutfile,"%d\t%d\t%18.18lf\n",conf_id_iter,
						i*meastopo_params.stout_measinterval,
						stout_topo_ch[i]);
			}
        	        fclose(stoutoutfile);
            	    }
	    }

            printf("MPI%02d - Printing gauge obs - only by master rank...\n",
                    devinfo.myrank);
            if(devinfo.myrank ==0){
                FILE *goutfile = fopen(gauge_outfilename,"at");
                if(!goutfile){
                    goutfile = fopen(gauge_outfilename,"wt");
                    strcpy(gauge_outfile_header,"#conf_id\tacc\tplq\trect\tReP\tImP\n");
                    fprintf(goutfile,"%s",gauge_outfile_header);
                }
                if(goutfile){
                    if(id_iter<mc_params.therm_ntraj){
                        printf("Therm_iter %d",conf_id_iter );
                        printf("Placchetta= %.18lf    ", plq/GL_SIZE/6.0/3.0);
                        printf("Rettangolo= %.18lf\n",rect/GL_SIZE/6.0/3.0/2.0);
                    }else printf("Metro_iter %d   Placchetta= %.18lf    Rettangolo= %.18lf\n",conf_id_iter,plq/GL_SIZE/6.0/3.0,rect/GL_SIZE/6.0/3.0/2.0);


                    fprintf(goutfile,"%d\t%d\t",conf_id_iter,
                            accettate_therm+accettate_metro
                            -accettate_therm_old-accettate_metro_old);
                    fprintf(goutfile,"%.18lf\t%.18lf\t%.18lf\t%.18lf\n",
                            plq/GL_SIZE/6.0/3.0,
                            rect/GL_SIZE/6.0/3.0/2.0, 
                            creal(poly), cimag(poly));
                }
                fclose(goutfile);
            }
            //-------------------------------------------------// 
            fclose(file_label);
       
            //######################################################################################################################################//
            //############### SAVING #########################################################################################################//
            //######################################################################################################################################//
            
            
             //*****************SAVING CONF_STORE****************//
            
        //Here Store_conf is saved.
         //conf[0] is the store_conf.
            
            //---- SAVES GAUGE CONF AND RNG STATUS TO FILE ----//
            if(conf_id_iter%mc_params.storeconfinterval==0){ 
                char tempname[50];
                char serial[10];
                strcpy(tempname,mc_params.store_conf_name);
                sprintf(serial,".%05d",conf_id_iter);
                strcat(tempname,serial);
                printf("MPI%02d - Storing conf %s.\n",
                        devinfo.myrank, tempname);
                save_conf_wrapper(conf_hasenbusch[0],tempname,conf_id_iter,
                        debug_settings.use_ildg);
                strcpy(tempname,mc_params.RandGenStatusFilename);
                sprintf(serial,".%05d",conf_id_iter);
                strcat(tempname,serial);
                printf("MPI%02d - Storing rng status in %s.\n", 
                        devinfo.myrank , tempname);
                saverand_tofile(tempname);
            }

            if(conf_id_iter%mc_params.saveconfinterval==0){
                if (debug_settings.SaveAllAtEnd){
                    printf("MPI%02d - Saving conf %s.\n", devinfo.myrank,
                            mc_params.save_conf_name);
                    save_conf_wrapper(conf_hasenbusch[0],mc_params.save_conf_name, conf_id_iter,
                            debug_settings.use_ildg);
                    printf("MPI%02d - Saving rng status in %s.\n", devinfo.myrank, 
                            mc_params.RandGenStatusFilename);
                    saverand_tofile(mc_params.RandGenStatusFilename);
                }else printf(
                        "\n\nMPI%02d: WARNING, \'SaveAllAtEnd\'=0,NOT SAVING/OVERWRITING CONF AND RNG STATUS.\n\n\n", devinfo.myrank);
            }
            //-------------------------------------------------//

            gettimeofday(&tend_cycle, NULL);

            double update_time = (double) 
                (tend_cycle.tv_sec - tstart_cycle.tv_sec)+
                (double)(tend_cycle.tv_usec - tstart_cycle.tv_usec)/1.0e6;
            
            mc_params.max_update_time = (update_time > mc_params.max_update_time)?
                update_time :mc_params.max_update_time;


            if(0==devinfo.myrank){
                printf("Tot time : %f sec (with measurements)\n", update_time);
                if(debug_settings.save_diagnostics == 1){
                    FILE *foutfile = fopen(debug_settings.diagnostics_filename,"at");
                    fprintf(foutfile,"TOTTIME  %f \n",update_time);
                    fclose(foutfile);
                }

            }
  
        }


        if (GPSTATUS_FERMION_MEASURES == mc_params.next_gps){
    
            
            //--------- MISURA ROBA FERMIONICA ----------------//
            
            if(0 != mc_params.JarzynskiMode ){ // HALFWAY MEASUREMENTS FOR JARZYNSKI

                bf_param new_backfield_parameters = backfield_parameters;

                // DIRECT MODE 
                if(1 == mc_params.JarzynskiMode)
                    new_backfield_parameters.bz = backfield_parameters.bz + 
                        (double) (id_iter+0.5)/mc_params.MaxConfIdIter;
                // REVERSE MODE
                if(-1 == mc_params.JarzynskiMode)
                    new_backfield_parameters.bz = backfield_parameters.bz -
                        (double) (id_iter+0.5)/mc_params.MaxConfIdIter;


                if(0==devinfo.myrank){

                    printf("JarzynskiMode, iteration %d/%d (%d max for this run) - MEASUREMENTS AT HALFWAY \n",
                            id_iter,mc_params.MaxConfIdIter,mc_params.ntraj);
                    printf("JarzynskiMode - current bz value : %f (HALFWAY)\n", new_backfield_parameters.bz);
                }

                init_all_u1_phases(new_backfield_parameters,fermions_parameters);
#pragma acc update device(u1_back_phases[0:8*alloc_info.NDiffFlavs])
#pragma acc update device(u1_back_phases_f[0:8*alloc_info.NDiffFlavs])

            }

            check_unitarity_device(conf_hasenbusch[0],&max_unitarity_deviation,
                    &avg_unitarity_deviation);
            printf("\tMPI%02d: Avg/Max unitarity deviation on device: %e / %e\n", 
                    devinfo.myrank,avg_unitarity_deviation,max_unitarity_deviation);

            if(conf_id_iter % fm_par.measEvery == 0 )
                mc_params.next_gps = GPSTATUS_FERMION_MEASURES;

            struct timeval tf0, tf1;
            gettimeofday(&tf0, NULL);
            fermion_measures(conf_hasenbusch[0],fermions_parameters,
                    &fm_par, md_parameters.residue_metro,
                    md_parameters.max_cg_iterations,conf_id_iter,
                    plq/GL_SIZE/3.0/6.0,
                    rect/GL_SIZE/3.0/6.0/2.0);   

            gettimeofday(&tf1, NULL);

            double fermionMeasureTiming =
                (double) (tf1.tv_sec - tf0.tv_sec)+
                (double)(tf1.tv_usec - tf0.tv_usec)/1.0e6;

            if(debug_settings.save_diagnostics == 1){
                FILE *foutfile = 
                    fopen(debug_settings.diagnostics_filename,"at");

                if(conf_id_iter % fm_par.measEvery == 0 )
                    fprintf(foutfile,"FERMMEASTIME  %f \n",fermionMeasureTiming);
                fclose(foutfile);
            }
      
            
            
            //---- SAVES RNG STATUS TO FILE ----//
            
            //RNG STATUS
            if(conf_id_iter%mc_params.storeconfinterval==0){
                char tempname[50];
                char serial[10];
                strcpy(tempname,mc_params.RandGenStatusFilename);
                sprintf(serial,".%05d",conf_id_iter);
                strcat(tempname,serial);
                printf("MPI%02d - Storing rng status in %s.\n", 
                        devinfo.myrank , tempname);
                saverand_tofile(tempname);
            } 

            if(conf_id_iter%mc_params.saveconfinterval==0){
               if( debug_settings.SaveAllAtEnd){
                printf("MPI%02d - Saving rng status in %s.\n", devinfo.myrank, 
                        mc_params.RandGenStatusFilename);
                saverand_tofile(mc_params.RandGenStatusFilename);
               }
               else printf(
                    "\n\nMPI%02d: WARNING, \'SaveAllAtEnd\'=0,NOT SAVING/OVERWRITING RNG STATUS.\n\n\n", devinfo.myrank);
            }
            //-------------------------------------------------//

        }

        // determining next thing to do
        if(0 == conf_id_iter % fm_par.measEvery)
            mc_params.next_gps = GPSTATUS_FERMION_MEASURES;
        if(mc_params.measures_done == fm_par.SingleInvNVectors){
            mc_params.next_gps = GPSTATUS_UPDATE;
            mc_params.measures_done = 0;
        }

        // determining run condition
        if(0 == devinfo.myrank && RUN_CONDITION_TERMINATE != mc_params.run_condition){
        
            // program exits if it finds a file called "stop"

            FILE * test_stop = fopen("stop","r");
            if(test_stop){
                fclose(test_stop);
                printf("File  \'stop\' found, stopping cycle now.\n");
                mc_params.run_condition = RUN_CONDITION_TERMINATE;
            }

            // program exits if it time is running out
            
            struct timeval now;
            gettimeofday(&now,NULL);
            double total_duration = (double) 
                (now.tv_sec - mc_params.start_time.tv_sec)+
                (double)(now.tv_usec - mc_params.start_time.tv_usec)/1.0e6;

            double max_expected_duration_with_another_cycle;
            if(GPSTATUS_UPDATE == mc_params.next_gps){
               max_expected_duration_with_another_cycle = 
                    total_duration + 1.3*mc_params.max_update_time;
               printf("Next step, update : %ds\n",(int) mc_params.max_update_time);
            }
            if(GPSTATUS_FERMION_MEASURES == mc_params.next_gps){
               max_expected_duration_with_another_cycle = 
                    total_duration + 2*mc_params.max_flavour_cycle_time;
               printf("Next step, flavour measure cycle : %ds\n",
                       (int) mc_params.max_flavour_cycle_time);
            }

            if(max_expected_duration_with_another_cycle > mc_params.MaxRunTimeS){
                printf("Time is running out (%d of %d seconds elapsed),",
                        (int) total_duration, (int) mc_params.MaxRunTimeS);
                printf(" shutting down now.\n");
                printf("Total max expected duration: %d seconds",
                        (int) max_expected_duration_with_another_cycle);
                printf("(%d elapsed now)\n",(int) total_duration);
                //https://www.youtube.com/watch?v=MfGhlVcrc8U
                // but without that much pathos
                mc_params.run_condition = RUN_CONDITION_TERMINATE;
            }

            // program exits if MaxConfIdIter is reached
            if(conf_id_iter > mc_params.MaxConfIdIter ){

                printf("%s - MaxConfIdIter=%d reached, job done!",
                        devinfo.myrankstr, mc_params.MaxConfIdIter);
                printf("%s - shutting down now.\n", devinfo.myrankstr);
                mc_params.run_condition = RUN_CONDITION_TERMINATE;
            }
            // program exits if MTraj is reached
            if( id_iter > (mc_params.ntraj+id_iter_offset)){
                printf("%s - NTraj=%d reached, job done!",
                        devinfo.myrankstr, mc_params.ntraj);
                printf("%s - shutting down now.\n", devinfo.myrankstr);
                mc_params.run_condition = RUN_CONDITION_TERMINATE;
            }

        }
#ifdef MULTIDEVICE

        MPI_Bcast((void*)&(mc_params.run_condition),1,MPI_INT,0,MPI_COMM_WORLD);
        printf("MPI%02d - Broadcast of run condition %d from master...\n",
                devinfo.myrank, mc_params.run_condition);
        MPI_Bcast((void*)&(mc_params.next_gps),1,MPI_INT,0,MPI_COMM_WORLD);
        printf("MPI%02d - Broadcast of next global program status %d from master...\n",
                devinfo.myrank, mc_params.next_gps);

#endif

    }// while id_iter loop ends here             

    //QUI SALVA Le CONF_REPLICAS!!
    
     //*****************SAVING CONFS****************//
    
    //---- SAVES GAUGE CONF AND RNG STATUS TO FILE ----//
    for(replicas_counter=0;replicas_counter<rep->replicas_total_number;replicas_counter++){
      
        snprintf(rep_str,10,"replica_%d",replicas_counter);//inizializza rep_str
        strcat(mc_params.save_conf_name,rep_str); //appiccica rep_str in fondo.
        
        
        
    if (debug_settings.SaveAllAtEnd){
        printf("MPI%02d - Saving conf %s.\n", devinfo.myrank,
                mc_params.save_conf_name);
        save_conf_wrapper(conf_hasenbusch[replicas_counter],mc_params.save_conf_name, conf_id_iter,
                debug_settings.use_ildg);
        printf("MPI%02d - Saving rng status in %s.\n", devinfo.myrank, 
                mc_params.RandGenStatusFilename);
        saverand_tofile(mc_params.RandGenStatusFilename);
    }else printf(
            "\n\nMPI%02d: WARNING, \'SaveAllAtEnd\'=0,NOT SAVING/OVERWRITING CONF AND RNG STATUS.\n\n\n", devinfo.myrank);
    //-------------------------------------------------//
        
        strcpy(mc_params.save_conf_name,aux_name_file);
        
    }//end for replicas.


    if(0 == devinfo.myrank && debug_settings.SaveAllAtEnd){
/*
	    printf("Saving global program status...\n");
	    printf("%d %f %f %d\n",
			    mc_params.next_gps,
			    mc_params.max_flavour_cycle_time,
			    mc_params.max_update_time,
			    mc_params.measures_done);

	    printf("#mc_params.next_gps,mc_params.max_flavour_cycle_time,\n#mc_params.max_update_time,mc_params.measures_done\n");

	    FILE * gps_file = fopen(mc_params.statusFileName, "w");  
	    fprintf(gps_file,"%d %f %f %d\n",
			    mc_params.next_gps,
			    mc_params.max_flavour_cycle_time,
			    mc_params.max_update_time,
			    mc_params.measures_done);
	    fprintf(gps_file,"#mc_params.next_gps,mc_params.max_flavour_cycle_time,\n#mc_params.max_update_time,mc_params.measures_done\n");
	    fclose(gps_file);
*/
	    save_global_program_status(mc_params); // THIS FUNCTION IN SOME CASES DOES NOT WORK
     }




    printf("MPI%02d: Double precision free [CORE]\n", devinfo.myrank);
    mem_free_core();
    
 
    
    printf("MPI%02d: Double precision free [EXTENDED]\n", devinfo.myrank);
    mem_free_extended();


    

    if(inverter_tricks.useMixedPrecision || md_parameters.singlePrecMD){
        printf("MPI%02d: Single precision free [CORE]\n", devinfo.myrank);
        mem_free_core_f();
    }
    if( md_parameters.singlePrecMD){
        printf("MPI%02d: Signle precision free [EXTENDED]\n", devinfo.myrank);
        mem_free_extended_f();
    }

    printf("MPI%02d: freeing device nnp and nnm\n", devinfo.myrank);
#pragma acc exit data delete(nnp_openacc)
#pragma acc exit data delete(nnm_openacc)

    printf("\n  MPI%02d - Prima dello shutdown, memoria allocata: %zu \n\n\n",devinfo.myrank,memory_used);
    
    struct memory_allocated_t *all=memory_allocated_base;
    
    while(all!=NULL)
      {
	printf("\n  MPI%02d - Va disallocato: %s di taglia %zu (o va conteggiato)\n\n\n",devinfo.myrank,all->varname,all->size);
	//free_wrapper(all->ptr);
	all=all->next;
    };

#ifndef __GNUC__
    //////  OPENACC CONTEXT CLOSING    //////////////////////////////////////////////////////////////
    shutdown_acc_device(my_device_type);
    /////////////////////////////////////////////////////////////////////////////////////////////////
#endif


#ifdef MULTIDEVICE
    shutdown_multidev();
#endif
    printf("The End\n");
    return 0;
}
