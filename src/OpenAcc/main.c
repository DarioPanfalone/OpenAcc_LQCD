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

#include "../DbgTools/dbgtools.h" // DEBUG
#include "../Include/debug.h"
#include "../Include/fermion_parameters.h"
#include "../Include/montecarlo_parameters.h"
#include "../Include/inverter_tricks.h"
#include "../Include/setting_file_parser.h"
#include "../Meas/ferm_meas.h"
#include "../Meas/gauge_meas.h"
#include "../Meas/polyakov.h"
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
#ifdef __GNUC__
#include "sys/time.h"
#endif






int conf_id_iter;
int verbosity_lv;

int main(int argc, char* argv[]){

#define  start_opt 0 // 0 --> COLD START; 1 --> START FROM SAVED CONF
    struct timeval tinit;
    gettimeofday ( &tinit, NULL );



    // READ input file.
#ifdef MULTIDEVICE
    pre_init_multidev1D(&devinfo);
#endif

    set_global_vars_and_fermions_from_input_file(argv[1]);
    verbosity_lv = debug_settings.input_vbl;

#ifdef MULTIDEVICE
    init_multidev1D(&devinfo);
#else
    devinfo.myrank = 0;
    devinfo.nranks = 1;
#endif


    if(0==devinfo.myrank){
        if(0 != mc_params.JarzynskiMode){
            printf("********************************************\n");
            printf("                JARZYNSKI MODE              \n");
            printf(" check which parameter corresponds to what! \n");
            printf("********************************************\n");


        }
        if(debug_settings.do_norandom_test){
            printf("*******************************************\n");
            printf("      WELCOME. This is a NORANDOM test.    \n");
            printf("  MOST things will not be random generated,\n");
            printf("         but read from memory instead.     \n");
            printf("               CHECK THE CODE!!            \n");
            printf("ALSO: setting the number of trajectories to 1.\n");
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
    select_init_acc_device(my_device_type, devinfo.myrank%devinfo.proc_per_node);
#else
    select_init_acc_device(my_device_type, devinfo.single_dev_choice);
#endif
    printf("Device Selected : OK \n");
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
        printf("MPI%02d - Finalizing...\n",devinfo.myrank);
#ifdef MULTIDEVICE
        MPI_Finalize();
#endif
        exit(1);
    }



    mem_alloc();
    printf("MPI%02d - Allocazione della memoria (double) : OK \n",devinfo.myrank);
    mem_alloc_f();
    printf("MPI%02d - Allocazione della memoria (float) : OK \n",devinfo.myrank);
    compute_nnp_and_nnm_openacc();
    printf("MPI%02d - nn computation : OK \n",devinfo.myrank);
    init_all_u1_phases(backfield_parameters,fermions_parameters);

    printf("MPI%02d - u1_backfield initialization (float & double): OK \n",devinfo.myrank);

    initialize_md_global_variables(md_parameters);
    printf("MPI%02d - init md vars : OK \n",devinfo.myrank);



    //################## INIZIALIZZAZIONE DELLA CONFIGURAZIONE #######################
    // start from saved conf

    if(debug_settings.do_norandom_test){
        if(!read_conf_wrapper(conf_acc,"conf_norndtest",&conf_id_iter,debug_settings.use_ildg)){
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
        if(!read_conf_wrapper(conf_acc,mc_params.save_conf_name,
                    &conf_id_iter,debug_settings.use_ildg)){
            // READS ALSO THE conf_id_iter
            printf("MPI%02d - Stored Gauge Conf \"%s\" Read : OK \n",
                    devinfo.myrank, mc_params.save_conf_name);

        }
        else{
            generate_Conf_cold(conf_acc,mc_params.eps_gen);
            printf("MPI%02d - Cold Gauge Conf Generated : OK \n",
                    devinfo.myrank);
            conf_id_iter=0;
        }
    }

    //#################################################################################  



    double max_unitarity_deviation,avg_unitarity_deviation;
    check_unitarity_host(conf_acc,&max_unitarity_deviation,&avg_unitarity_deviation);
    printf("\tMPI%02d: Avg_unitarity_deviation on host: %e\n", devinfo.myrank, 
            avg_unitarity_deviation);
    printf("\tMPI%02d: Max_unitarity_deviation on host: %e\n", devinfo.myrank,
            max_unitarity_deviation);



#pragma acc data  copy(conf_acc[0:conf_acc_size]) \
    create(ipdot_acc[0:8]) create(aux_conf_acc[0:8])\
    create(auxbis_conf_acc[0:8]) create(ferm_chi_acc[0:NPS_tot])\
    create(ferm_phi_acc[0:NPS_tot])  create(ferm_out_acc[0:NPS_tot])\
    create(ferm_shiftmulti_acc[0:maxNeededShifts])\
    create(kloc_r[0:1])  create(kloc_h[0:1])  create(kloc_s[0:1])\
    create(kloc_p[0:1])  create(k_p_shiftferm[0:maxApproxOrder])\
    create(aux1[0:1])\
    create(momenta[0:8])\
    create(momenta_backup[0:8*debug_settings.do_reversibility_test])\
    create(  conf_acc_bkp[0:8*debug_settings.do_reversibility_test])\
    copyin(nnp_openacc) copyin(nnm_openacc)\
    create(local_sums[0:2]) create(d_local_sums[0:2])\
    copyin(fermions_parameters[0:NDiffFlavs])\
    copyin(deltas_Omelyan[0:7]) \
    copyin(u1_back_phases[0:8*NDiffFlavs])\
    create(conf_acc_f[0:8])\
    create(ipdot_acc_f[0:8]) create(aux_conf_acc_f[0:8])\
    create(auxbis_conf_acc_f[0:8]) create(ferm_chi_acc_f[0:NPS_tot])\
    create(ferm_phi_acc_f[0:NPS_tot])  create(ferm_out_acc_f[0:NPS_tot])\
    create(ferm_shiftmulti_acc_f[0:maxNeededShifts])\
    create(kloc_r_f[0:1])  create(kloc_h_f[0:1])  create(kloc_s_f[0:1])\
    create(kloc_p_f[0:1])  create(k_p_shiftferm_f[0:maxApproxOrder])\
    create(aux1_f[0:1])\
    create(momenta_f[0:8]) \
    create(local_sums_f[0:2]) create(d_local_sums_f[0:2])\
    copyin(deltas_Omelyan_f[0:7]) \
    copyin(u1_back_phases_f[0:8*NDiffFlavs])\
    create(ipdot_g_old_f[0:8]) create(ipdot_f_old_f[0:8])\
    copyin(mag_obs_re[0:8*NDiffFlavs])\
    copyin(mag_obs_im[0:8*NDiffFlavs])\
    create(ipdot_g_old[0:8*debug_settings.save_diagnostics])\
    create(ipdot_f_old[0:8*debug_settings.save_diagnostics])
    {
#ifdef STOUT_FERMIONS
#pragma acc data create(aux_th[0:8]) create(aux_ta[0:8])\
        create(gstout_conf_acc_arr[0:(8*act_params.stout_steps)])\
        create(glocal_staples[0:8]) create(gipdot[0:8])\
        create(aux_th_f[0:8]) create(aux_ta_f[0:8])\
        create(gstout_conf_acc_arr_f[0:(8*act_params.stout_steps)])\
        create(glocal_staples_f[0:8]) create(gipdot_f[0:8]) 
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
            printf("\tMPI%02d: Therm_iter %d Placchetta    = %.18lf \n",
                    devinfo.myrank, conf_id_iter,plq/GL_SIZE/6.0/3.0);
            rect = calc_rettangolo_soloopenacc(conf_acc,aux_conf_acc,local_sums);


            printf("\tMPI%02d: Therm_iter %d Rettangolo    = %.18lf \n",
                    devinfo.myrank, conf_id_iter,rect/GL_SIZE/6.0/3.0/2.0);

            poly =  (*polyakov_loop[geom_par.tmap])(conf_acc);//misura polyakov loop
            printf("\tMPI%02d: Therm_iter %d Polyakov Loop = (%.18lf, %.18lf)  \n",
                    devinfo.myrank, conf_id_iter,creal(poly),cimag(poly));


            //char confile_dbg[50];
            //sprintf(confile_dbg,"conf_ascii_test_%s", devinfo.myrankstr);
            //dbg_print_su3_soa(conf_acc,confile_dbg,0);


            //            MPI_Finalize(); // DEBUG
            //            return 0 ;      // DEBUG


            if(0 == mc_params.ntraj && 0 == mc_params.JarzynskiMode ){ // MEASURES ONLY

                printf("\n#################################################\n");
                printf("\tMEASUREMENTS ONLY ON FILE %s\n", mc_params.save_conf_name);
                printf("\n#################################################\n");

                //--------- MISURA ROBA FERMIONICA ----------------//
                //
                if(0 == devinfo.myrank)  printf("Fermion Measurements: see file %s\n",
                        fm_par.fermionic_outfilename);
                fermion_measures(conf_acc,fermions_parameters,
                        &fm_par, md_parameters.residue_metro, 
                        md_parameters.max_cg_iterations, id_iter_offset);


                //-------------------------------------------------// 
                //--------- MISURA ROBA DI GAUGE ------------------//
                if(0 == devinfo.myrank ) printf("Misure di Gauge:\n");
                plq = calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
                rect = calc_rettangolo_soloopenacc(conf_acc,aux_conf_acc,local_sums);
                poly =  (*polyakov_loop[geom_par.tmap])(conf_acc);//misura polyakov loop

                printf("Plaquette     : %.18lf\n" ,plq/GL_SIZE/3.0/6.0);
                printf("Rectangle     : %.18lf\n" ,rect/GL_SIZE/3.0/6.0/2.0);
                printf("Polyakov Loop : (%.18lf,%.18lf) \n",creal(poly),cimag(poly));


            }else printf("MPI%02d: Starting generation of Configurations.\n",
                    devinfo.myrank);

            // THERMALIZATION & METRO    ----   UPDATES //

            for(int id_iter=id_iter_offset;id_iter<(mc_params.ntraj+id_iter_offset);
                    id_iter++){

                struct timeval tstart_cycle;
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
#pragma acc update device(u1_back_phases[0:8*NDiffFlavs])
#pragma acc update device(u1_back_phases_f[0:8*NDiffFlavs])

                }


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
                            conf_id_iter,mc_params.ntraj+id_iter_offset,
                            geom_par.gnx,geom_par.gny,
                            geom_par.gnz,geom_par.gnt,
                            act_params.beta);
                    printf(  "#################################################\n\n");
                }



                //--------- CONF UPDATE ----------------//
                if(id_iter<mc_params.therm_ntraj){
                    accettate_therm = UPDATE_SOLOACC_UNOSTEP_VERSATILE(conf_acc,
                            md_parameters.residue_metro,md_parameters.residue_md,
                            id_iter-id_iter_offset,
                            accettate_therm,0,md_parameters.max_cg_iterations);
                }else{
                    accettate_metro = UPDATE_SOLOACC_UNOSTEP_VERSATILE(conf_acc,
                            md_parameters.residue_metro,md_parameters.residue_md,
                            id_iter-id_iter_offset-accettate_therm,accettate_metro,1,
                            md_parameters.max_cg_iterations);
                    if(0==devinfo.myrank){
                        int iterations = id_iter-id_iter_offset-accettate_therm +1;
                        double acceptance = (double) accettate_metro / iterations;
                        double acc_err = 
                            sqrt(accettate_metro*(iterations-accettate_metro)/iterations)
                            /iterations;
                        printf("Estimated acceptance for this run: %f +- %f\n",acceptance,
                                acc_err);
                    }
                }



#pragma acc update host(conf_acc[0:8])
                //---------------------------------------//

                //--------- MISURA ROBA FERMIONICA ----------------//
                //
                check_unitarity_device(conf_acc,&max_unitarity_deviation,
                        &avg_unitarity_deviation);
                printf("\tMPI%02d: Avg/Max unitarity deviation on device: %e / %e\n", 
                        devinfo.myrank,avg_unitarity_deviation,max_unitarity_deviation);

                fermion_measures(conf_acc,fermions_parameters,
                        &fm_par, md_parameters.residue_metro,
                        md_parameters.max_cg_iterations,
                        id_iter) ;

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

                //---- SAVES GAUGE CONF AND RNG STATUS TO FILE ----//
                if(conf_id_iter%mc_params.storeconfinterval==0){
                    char tempname[50];
                    char serial[10];
                    strcpy(tempname,mc_params.store_conf_name);
                    sprintf(serial,".%05d",conf_id_iter);
                    strcat(tempname,serial);
                    printf("MPI%02d - Storing conf %s.\n",
                            devinfo.myrank, tempname);
                    save_conf_wrapper(conf_acc,tempname,conf_id_iter,
                            debug_settings.use_ildg);
                    strcpy(tempname,mc_params.RandGenStatusFilename);
                    sprintf(serial,".%05d",conf_id_iter);
                    strcat(tempname,serial);
                    printf("MPI%02d - Storing rng status in %s.\n", 
                            devinfo.myrank , tempname);
                    saverand_tofile(tempname);
                }
                if(conf_id_iter%mc_params.saveconfinterval==0){
                    printf("MPI%02d - Saving conf %s.\n", devinfo.myrank,
                            mc_params.save_conf_name);
                    save_conf_wrapper(conf_acc,mc_params.save_conf_name, conf_id_iter,
                            debug_settings.use_ildg);
                    printf("MPI%02d - Saving rng status in %s.\n", devinfo.myrank, 
                            mc_params.RandGenStatusFilename);
                    saverand_tofile(mc_params.RandGenStatusFilename);
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

                    if(0==devinfo.myrank) printf("Tot time : %f sec (with measurements)\n", cycle_duration);


                    double total_duration = (double) 
                        (tend_cycle.tv_sec - tinit.tv_sec)+
                        (double)(tend_cycle.tv_usec - tinit.tv_usec)/1.0e6;
                    double max_expected_duration_with_another_cycle = 
                        total_duration + 2*cycle_duration ; 

                    if(max_expected_duration_with_another_cycle > mc_params.MaxRunTimeS){
                        printf("Time is running out (%d of %d seconds elapsed),",
                                (int) total_duration, (int) mc_params.MaxRunTimeS);
                        printf(" shutting down now.\n");
                        //https://www.youtube.com/watch?v=MfGhlVcrc8U
                        // but without that much pathos
                        run_condition = 0;
                    }

                    // program exits if MaxConfIdIter is reached
                    if(conf_id_iter > mc_params.MaxConfIdIter ){

                        printf("%s - MaxConfIdIter=%d reached, job done!",
                                devinfo.myrankstr, mc_params.MaxConfIdIter);
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

            if(debug_settings.SaveAllAtEnd){
                if(mc_params.ntraj > 0 ){
                    save_conf_wrapper(conf_acc,mc_params.save_conf_name, conf_id_iter,
                            debug_settings.use_ildg );
                }
                saverand_tofile(mc_params.RandGenStatusFilename);
            }
            else 
                printf(
                        "\n\nMPI%02d: WARNING, \'SaveAllAtEnd\'=0,NOT SAVING/OVERWRITING CONF AND RNG STATUS.\n\n\n", devinfo.myrank);
            //-------------------------------------------------//


            plq = calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
            topoch = compute_topological_charge(conf_acc,aux_conf_acc,d_local_sums);
            printf("COOL 0  Placchetta= %.18lf  TopCh= %.18lf \n",plq/GL_SIZE/6.0/3.0,topoch);

            //               // You might want to put this inside the loop
            //               for(int icool=0;icool<5000;icool++){
            //               cool_conf(conf_acc,aux_conf_acc);
            //               plq = calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
            //               topoch = compute_topological_charge(conf_acc,aux_conf_acc,d_local_sums);
            //               printf("COOL %d  Placchetta= %.18lf  TopCh= %.18lf \n",icool+1,plq/GL_SIZE/6.0/3.0,topoch);
            //               }
            //               
            //

#ifdef STOUT_FERMIONS
        } // end pragma acc data (le cose del caso stout)
#endif

    }// end pragma acc data



#ifndef __GNUC__
    //////  OPENACC CONTEXT CLOSING    //////////////////////////////////////////////////////////////
    shutdown_acc_device(my_device_type);
    /////////////////////////////////////////////////////////////////////////////////////////////////
#endif


#ifdef MULTIDEVICE
    shutdown_multidev();
#endif

    mem_free();
    mem_free_f();

    return 0;
}

