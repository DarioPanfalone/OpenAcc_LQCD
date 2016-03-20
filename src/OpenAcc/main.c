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

#include "../Include/fermion_parameters.h"
#include "../Include/setting_file_parser.h"
#include "../Include/markowchain.h"
#include "../DbgTools/debug_macros_glvarcheck.h"
#include "../RationalApprox/rationalapprox.h"
#include "./struct_c_def.h"
#include "./alloc_vars.h"
#include "./io.h"
#include "./fermionic_utilities.h"
#include "./su3_utilities.h"
#include "./su3_measurements.h"
#include "./random_assignement.h"
#include "./fermion_matrix.h"
#include "./inverter_full.h"
#include "./find_min_max.h"
#include "./inverter_multishift_full.h"
#include "./rettangoli.h"
#include "./ipdot_gauge.h"
#include "../Meas/gauge_meas.h"
#include "../Meas/polyakov.h"
#include "../Meas/ferm_meas.h"
#include "./stouting.h"
#include "./fermion_force.h"
#include "./md_integrator.h"
#include "./update_versatile.h"
#include "./cooling.h"
#include "./backfield.h"
#include "./action.h"
#include "../Rand/random.h"

#include "../Mpi/multidev.h"
#include "../Mpi/communications.h"
#include "./deviceinit.h"

#ifdef __GNUC__
#include "sys/time.h"
#endif

#include "../DbgTools/dbgtools.h" // DEBUG




//#define NORANDOM  // FOR debug, check also update_versatile.c 

int conf_id_iter;
int verbosity_lv = 5;// 5 should print everything.

int main(int argc, char* argv[]){

#define  start_opt 0 // 0 --> COLD START; 1 --> START FROM SAVED CONF
    struct timeval tinit;
    gettimeofday ( &tinit, NULL );

#ifdef NORANDOM
    printf("WELCOME! NORANDOM MODE. (main()) \n" );
#endif

    printf("WELCOME! \n");
    // READ input file.
#ifdef MULTIDEVICE
    pre_init_multidev1D(&devinfo);
#endif

    set_global_vars_and_fermions_from_input_file(argv[1]);
    verbosity_lv = mkwch_pars.input_vbl;

#ifdef MULTIDEVICE
    init_multidev1D(&devinfo);
#else
    devinfo.myrank = 0;
    devinfo.nranks = 1;
#endif


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
    select_init_acc_device(my_device_type, devinfo.single_dev_choice);
    printf("Device Selected : OK \n");
#endif

    unsigned int myseed_default =  (unsigned int) mkwch_pars.seed; 

#ifdef MULTIDEVICE
    myseed_default =  (unsigned int) (myseed_default + devinfo.myrank) ;
    char myrank_string[6];
    sprintf(myrank_string,".R%d",devinfo.myrank);
    strcat(mkwch_pars.RandGenStatusFilename,myrank_string);
#endif




    initrand_fromfile(mkwch_pars.RandGenStatusFilename,myseed_default);


    // INIT FERM PARAMS AND READ RATIONAL APPROX COEFFS
    if(init_ferm_params(fermions_parameters)){
        printf("MPI%02d - Finalizing...\n",devinfo.myrank);
#ifdef MULTIDEVICE
    MPI_Finalize();
#endif
        exit(1);
    }



    mem_alloc();
    printf("MPI%02d - Allocazione della memoria : OK \n",devinfo.myrank);
    compute_nnp_and_nnm_openacc();
    printf("MPI%02d - nn computation : OK \n",devinfo.myrank);
    init_all_u1_phases(backfield_parameters,fermions_parameters);

    printf("MPI%02d - u1_backfield initialization : OK \n",devinfo.myrank);

    initialize_md_global_variables(md_parameters);
    printf("MPI%02d - init md vars : OK \n",devinfo.myrank);



    //################## INIZIALIZZAZIONE DELLA CONFIGURAZIONE #######################
    // start from saved conf

#ifdef NORANDOM
    if(!read_conf_wrapper(conf_acc,"conf_norndtest",&conf_id_iter,mkwch_pars.use_ildg)){
        // READS ALSO THE conf_id_iter
        printf("MPI%02d - Stored Gauge Conf conf_norndtest Read : OK\n",devinfo.myrank);
    }
    else{
        // cold start
        printf("MPI%02d - COMPILED IN NORANDOM MODE. A CONFIGURATION FILE NAMED\
                \"conf_norndtest\" MUST BE PRESENT\n",devinfo.myrank);
        exit(1);
    }

#else
    if(!read_conf_wrapper(conf_acc,mkwch_pars.save_conf_name,
                &conf_id_iter,mkwch_pars.use_ildg)){
        // READS ALSO THE conf_id_iter
        printf("MPI%02d - Stored Gauge Conf \"%s\" Read : OK \n",
                devinfo.myrank, mkwch_pars.save_conf_name);

    }
    else{
        generate_Conf_cold(conf_acc,mkwch_pars.eps_gen);
        printf("MPI%02d - Cold Gauge Conf Generated : OK \n",
                devinfo.myrank);
        conf_id_iter=0;
    }
#endif
    //#################################################################################  



    double max_unitarity_deviation,avg_unitarity_deviation;
    check_unitarity_host(conf_acc,&max_unitarity_deviation,&avg_unitarity_deviation);
    printf("\tMPI%02d: Avg_unitarity_deviation on host: %e\n", devinfo.myrank, 
            avg_unitarity_deviation);
    printf("\tMPI%02d: Max_unitarity_deviation on host: %e\n", devinfo.myrank,
            max_unitarity_deviation);



#pragma acc data   copy(conf_acc[0:8]) \
    create(ipdot_acc[0:8]) create(aux_conf_acc[0:8])\
    create(auxbis_conf_acc[0:8]) create(ferm_chi_acc[0:NPS_tot])\
    create(ferm_phi_acc[0:NPS_tot])  create(ferm_out_acc[0:NPS_tot])\
    create(ferm_shiftmulti_acc[0:max_ps*MAX_APPROX_ORDER])\
    create(kloc_r[0:1])  create(kloc_h[0:1])  create(kloc_s[0:1])\
    create(kloc_p[0:1])  create(k_p_shiftferm[0:MAX_APPROX_ORDER])\
    create(momenta[0:8]) copyin(nnp_openacc) copyin(nnm_openacc)\
    create(local_sums[0:2]) create(d_local_sums[0:2])\
    copyin(fermions_parameters[0:NDiffFlavs])\
    copyin(deltas_Omelyan[0:7]) \
    copyin(u1_back_phases[0:8*NDiffFlavs])\
    create(ipdot_g_old[0:8]) create(ipdot_f_old[0:8])
    {
#ifdef STOUT_FERMIONS
#pragma acc data create(aux_th[0:8]) create(aux_ta[0:8])\
        create(gstout_conf_acc_arr[0:(8*act_params.stout_steps)])\
        create(glocal_staples[0:8]) create(gipdot[0:8]) 
        {
#endif

            double plq,rect,topoch;
            d_complex poly;

            int accettate_therm=0;
            int accettate_metro=0;
            int accettate_therm_old=0;
            int accettate_metro_old=0;
            int id_iter_offset=conf_id_iter;
            plq = calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
            rect = calc_rettangolo_soloopenacc(conf_acc,aux_conf_acc,local_sums);
            poly =  (*polyakov_loop[geom_par.tmap])(conf_acc);//misura polyakov loop
            printf("\tMPI%02d: Therm_iter %d Placchetta    = %.18lf \n",
                    devinfo.myrank, conf_id_iter,plq/NSITES/6.0/3.0);
            printf("\tMPI%02d: Therm_iter %d Rettangolo    = %.18lf \n",
                    devinfo.myrank, conf_id_iter,rect/NSITES/6.0/3.0/2.0);
            printf("\tMPI%02d: Therm_iter %d Polyakov Loop = (%.18lf, %.18lf)  \n",
                    devinfo.myrank, conf_id_iter,creal(poly),cimag(poly));

            char confile_dbg[50];
            sprintf(confile_dbg,"conf_ascii_test_%s", devinfo.myrankstr);
            dbg_print_su3_soa(conf_acc,confile_dbg,0);

//            MPI_Finalize();
//            return 0;




            if(mkwch_pars.ntraj==0){ // MEASURES ONLY

                printf("\n#################################################\n");
                printf("\tMEASUREMENTS ONLY ON FILE %s\n", mkwch_pars.save_conf_name);
                printf("\n#################################################\n");

                //--------- MISURA ROBA FERMIONICA ----------------//
                //
                if(devinfo.myrank == 0)  printf("Fermion Measurements: see file %s\n",
                        fm_par.fermionic_outfilename);
                fermion_measures(conf_acc,fermions_parameters,
                        &fm_par, mkwch_pars.residue_metro, id_iter_offset) ;


                //-------------------------------------------------// 
                //--------- MISURA ROBA DI GAUGE ------------------//
                if(devinfo.myrank == 0 ) printf("Misure di Gauge:\n");
                plq = calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
                rect = calc_rettangolo_soloopenacc(conf_acc,aux_conf_acc,local_sums);
                poly =  (*polyakov_loop[geom_par.tmap])(conf_acc);//misura polyakov loop

                printf("Plaquette     : %.18lf\n" ,plq/NSITES/3.0/6.0);
                printf("Rectangle     : %.18lf\n" ,rect/NSITES/3.0/6.0/2.0);
                printf("Polyakov Loop : (%.18lf,%.18lf) \n",creal(poly),cimag(poly));


            }else printf("MPI%02d: Starting generation of Configurations.\n",
                    devinfo.myrank);

            // THERMALIZATION & METRO    ----   UPDATES //

            for(int id_iter=id_iter_offset;id_iter<(mkwch_pars.ntraj+id_iter_offset);
                    id_iter++){
                struct timeval tstart_cycle;
                gettimeofday(&tstart_cycle, NULL);

                check_unitarity_device(conf_acc,&max_unitarity_deviation,
                        &avg_unitarity_deviation);
                printf("\tMPI%02d: Avg/Max unitarity deviation on device: %e / %e\n", 
                        devinfo.myrank,avg_unitarity_deviation,max_unitarity_deviation);
                accettate_therm_old = accettate_therm;
                accettate_metro_old = accettate_metro;
                conf_id_iter++;


                if(devinfo.myrank ==0 ){
                    printf("\n#################################################\n");
                    printf(  "   GENERATING CONF %d of %d, %dx%dx%dx%d,%1.3f \n",
                            conf_id_iter,mkwch_pars.ntraj+id_iter_offset,
                            geom_par.gnx,geom_par.gny,
                            geom_par.gnz,geom_par.gnt,
                            act_params.beta);
                    printf(  "#################################################\n\n");
                }



                //--------- CONF UPDATE ----------------//
                if(id_iter<mkwch_pars.therm_ntraj){
                    accettate_therm = UPDATE_SOLOACC_UNOSTEP_VERSATILE(conf_acc,
                            mkwch_pars.residue_metro,md_parameters.residue_md,
                            id_iter-id_iter_offset,
                            accettate_therm,0);
                }else{
                    accettate_metro = UPDATE_SOLOACC_UNOSTEP_VERSATILE(conf_acc,
                            mkwch_pars.residue_metro,md_parameters.residue_md,
                            id_iter-id_iter_offset-accettate_therm,accettate_metro,1);
                }



#pragma acc update host(conf_acc[0:8])
                //---------------------------------------//

                //--------- MISURA ROBA FERMIONICA ----------------//
                //
                fermion_measures(conf_acc,fermions_parameters,
                        &fm_par, mkwch_pars.residue_metro,id_iter) ;


                //-------------------------------------------------// 
                //--------- MISURA ROBA DI GAUGE ------------------//
                plq  = calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
                rect = calc_rettangolo_soloopenacc(conf_acc,aux_conf_acc,local_sums);
                poly =  (*polyakov_loop[geom_par.tmap])(conf_acc);


                printf("MPI%02d - Printing gauge obs - only by master rank...\n",
                        devinfo.myrank);
                if(devinfo.myrank ==0 ){

                    FILE *goutfile = fopen(gauge_outfilename,"at");
                    if(!goutfile){
                        goutfile = fopen(gauge_outfilename,"wt");
                        strcpy(gauge_outfile_header,"#conf_id\tacc\tplq\trect\tReP\tImP\n");
                        fprintf(goutfile,"%s",gauge_outfile_header);
                    }
                    if(goutfile){
                        if(id_iter<mkwch_pars.therm_ntraj){
                            printf("Therm_iter %d",conf_id_iter );
                            printf("Placchetta= %.18lf    ", plq/NSITES/6.0/3.0);
                            printf("Rettangolo= %.18lf\n",rect/NSITES/6.0/3.0/2.0);


                        }else printf("Metro_iter %d   Placchetta= %.18lf    Rettangolo= %.18lf\n",conf_id_iter,plq/NSITES/6.0/3.0,rect/NSITES/6.0/3.0/2.0);


                        fprintf(goutfile,"%d\t%d\t",conf_id_iter,
                                accettate_therm+accettate_metro
                                -accettate_therm_old-accettate_metro_old);
                        fprintf(goutfile,"%.18lf\t%.18lf\t%.18lf\t%.18lf\n",
                                plq/NSITES/6.0/3.0,
                                rect/NSITES/6.0/3.0/2.0, 
                                creal(poly), cimag(poly));

                    }
                    fclose(goutfile);
                }
                //-------------------------------------------------//

                //---- SAVES GAUGE CONF AND RNG STATUS TO FILE ----//
                if(conf_id_iter%mkwch_pars.storeconfinterval==0){
                    char tempname[50];
                    char serial[10];
                    strcpy(tempname,mkwch_pars.store_conf_name);
                    sprintf(serial,".%05d",conf_id_iter);
                    strcat(tempname,serial);
                    printf("MPI%02d - Storing conf %s.\n",
                            devinfo.myrank, tempname);
                    save_conf_wrapper(conf_acc,tempname,conf_id_iter,
                            mkwch_pars.use_ildg);
                    strcpy(tempname,mkwch_pars.RandGenStatusFilename);
                    sprintf(serial,".%05d",conf_id_iter);
                    strcat(tempname,serial);
                    printf("MPI%02d - Storing rng status in %s.\n", 
                            devinfo.myrank , tempname);
                    saverand_tofile(tempname);
                }
                if(conf_id_iter%mkwch_pars.saveconfinterval==0){
                    printf("MPI%02d - Saving conf %s.\n", devinfo.myrank,
                            mkwch_pars.save_conf_name);
                    save_conf_wrapper(conf_acc,mkwch_pars.save_conf_name, conf_id_iter,
                            mkwch_pars.use_ildg);
                    printf("MPI%02d - Saving rng status in %s.\n", devinfo.myrank, 
                            mkwch_pars.RandGenStatusFilename);
                    saverand_tofile(mkwch_pars.RandGenStatusFilename);
                }

                //-------------------------------------------------//
                // program exits if it finds a file called "stop"

                int run_condition = 1;

                
                if(devinfo.myrank ==0 ){
                    FILE * test_stop = fopen("stop","r");
                    if(test_stop){
                        fclose(test_stop);
                        printf("File  \'stop\' found, stopping cycle now.\n");
                        run_condition = 0;
                    }

                    // program exits if it time is running out
                    struct timeval tend_cycle;
                    gettimeofday(&tend_cycle, NULL);

                    double cycle_duration = (double) 
                        (tend_cycle.tv_sec - tstart_cycle.tv_sec)+
                        (double)(tend_cycle.tv_usec - tstart_cycle.tv_usec)/1.0e6;
                    double total_duration = (double) 
                        (tend_cycle.tv_sec - tinit.tv_sec)+
                        (double)(tend_cycle.tv_usec - tinit.tv_usec)/1.0e6;
                    double max_expected_duration_with_another_cycle = 
                        total_duration + 2*cycle_duration ; 

                    if(max_expected_duration_with_another_cycle > mkwch_pars.MaxRunTimeS){
                        printf("Time is running out (%d of %d seconds elapsed),",
                                (int) total_duration, (int) mkwch_pars.MaxRunTimeS);
                        printf(" shutting down now.\n");
                        //https://www.youtube.com/watch?v=MfGhlVcrc8U
                        // but without that much pathos
                        run_condition = 0;
                    }

                    // program exits if MaxConfIdIter is reached
                    if(conf_id_iter >= mkwch_pars.MaxConfIdIter ){

                        printf("%s - MaxConfIdIter=%d reached, job done!",
                                devinfo.myrankstr, mkwch_pars.MaxConfIdIter);
                        printf("%s - shutting down now.\n", devinfo.myrankstr);
                        run_condition = 0;
                    }
                }
#ifdef MULTIDEVICE

                MPI_Bcast((void*)&run_condition,1,MPI_INT,0,MPI_COMM_WORLD);
                printf("MPI%02d - Broadcast of run %d condition from master...\n",
                        devinfo.myrank, run_condition);
#endif
                if(run_condition == 0) break;

            }// id_iter loop ends here             

            //---- SAVES GAUGE CONF AND RNG STATUS TO FILE ----//

            if(mkwch_pars.ntraj > 0){
                save_conf_wrapper(conf_acc,mkwch_pars.save_conf_name, conf_id_iter,
                        mkwch_pars.use_ildg );
                saverand_tofile(mkwch_pars.RandGenStatusFilename);
            }
            //-------------------------------------------------//


            plq = calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
            topoch = compute_topological_charge(conf_acc,aux_conf_acc,d_local_sums);
            printf("COOL 0  Placchetta= %.18lf  TopCh= %.18lf \n",plq/NSITES/6.0/3.0,topoch);

            //               // You might want to put this inside the loop
            //               for(int icool=0;icool<5000;icool++){
            //               cool_conf(conf_acc,aux_conf_acc);
            //               plq = calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
            //               topoch = compute_topological_charge(conf_acc,aux_conf_acc,d_local_sums);
            //               printf("COOL %d  Placchetta= %.18lf  TopCh= %.18lf \n",icool+1,plq/NSITES/6.0/3.0,topoch);
            //               }
            //               
            //

#ifdef STOUT_FERMIONS
        } // end pragma acc data (le cose del caso stout)
#endif

    }// end pragma acc data

    CHECKSTATUS(conf_acc);


#ifndef __GNUC__
    //////  OPENACC CONTEXT CLOSING    //////////////////////////////////////////////////////////////
    shutdown_acc_device(my_device_type);
    /////////////////////////////////////////////////////////////////////////////////////////////////
#endif

#ifdef MULTIDEVICE
    shutdown_multidev();
#endif

    mem_free();

    return 0;
}

