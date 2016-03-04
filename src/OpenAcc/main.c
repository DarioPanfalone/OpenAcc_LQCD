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
#include "../Include/init.h"
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
#include "../Mpi/communications.h"


//#define NORANDOM  // FOR debug, check also update_versatile.c 

int conf_id_iter;
int verbosity_lv = 5;// 5 should print everything.

int main(int argc, char* argv[]){

#define  start_opt 0 // 0 --> COLD START; 1 --> START FROM SAVED CONF


#ifdef NORANDOM
    printf("WELCOME! NORANDOM MODE. (main()) \n" );
#endif

    printf("WELCOME! \n");
    // READ input file.
    set_global_vars_and_fermions_from_input_file(argv[1]);
    //


#ifndef __GNUC__
    //////  OPENACC CONTEXT INITIALIZATION    //////////////////////////////////////////////////////
    // NVIDIA GPUs
    acc_device_t my_device_type = acc_device_nvidia;
    // AMD GPUs
    // acc_device_t my_device_type = acc_device_radeon;
    // Intel XeonPhi
    //acc_device_t my_device_type = acc_device_xeonphi;
    // Select device ID
    SELECT_INIT_ACC_DEVICE(my_device_type, dev_settings.device_choice);
    printf("Device Selected : OK \n");
#endif

    initrand((unsigned int) mkwch_pars.seed);
    verbosity_lv = mkwch_pars.input_vbl;
    // INIT FERM PARAMS AND READ RATIONAL APPROX COEFFS
    if(init_ferm_params(fermions_parameters)) exit(1);



    mem_alloc();
    printf("Allocazione della memoria : OK \n");
    compute_nnp_and_nnm_openacc();
    printf("nn computation : OK \n");
    init_all_u1_phases(backfield_parameters,fermions_parameters);



    printf("u1_backfield initialization : OK \n");
    
    initialize_md_global_variables(md_parameters);
    printf("init md vars : OK \n");

    //################## INIZIALIZZAZIONE DELLA CONFIGURAZIONE #######################
    // start from saved conf
    

#ifdef NORANDOM
    if(!read_conf(conf_rw,"conf_norndtest",&conf_id_iter),mkwch_pars.use_ildg){
        // READS ALSO THE conf_id_iter
        printf("Stored Gauge Conf conf_norndtest Read : OK\n",mkwch_pars.save_conf_name);
        send_lnh_subconf_to_buffer(conf_rw,conf_acc,0);
    }
    else{
        // cold start
        printf("COMPILED IN NORANDOM MODE. A CONFIGURATION FILE NAMED\
\"conf_norndtest\" MUST BE PRESENT\n");
        exit(1);
    }
    
#else
    if(!read_conf(conf_rw,mkwch_pars.save_conf_name,
                &conf_id_iter,mkwch_pars.use_ildg)){
       // READS ALSO THE conf_id_iter
       printf("Stored Gauge Conf \"%s\" Read : OK \n", mkwch_pars.save_conf_name);
       send_lnh_subconf_to_buffer(conf_rw,conf_acc,0);

    }
    else{
        generate_Conf_cold(conf_acc,mkwch_pars.eps_gen);
        printf("Cold Gauge Conf Generated : OK \n");
        conf_id_iter=0;
    }
#endif
    //#################################################################################  



    double max_unitarity_deviation,avg_unitarity_deviation;
    check_unitarity_host(conf_acc,&max_unitarity_deviation,&avg_unitarity_deviation);
    printf("\tAvg_unitarity_deviation on host: %e\n", avg_unitarity_deviation);
    printf("\tMax_unitarity_deviation on host: %e\n", max_unitarity_deviation);



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
            printf("Therm_iter %d Placchetta    = %.18lf \n",
                    conf_id_iter,plq/size/6.0/3.0);
            printf("Therm_iter %d Rettangolo    = %.18lf \n",
                    conf_id_iter,rect/size/6.0/3.0/2.0);
            printf("Therm_iter %d Polyakov Loop = (%.18lf, %.18lf)  \n",conf_id_iter,
                    creal(poly),cimag(poly));

            if(mkwch_pars.ntraj==0){ // MEASURES ONLY

                printf("\n#################################################\n");
                printf("\tMEASUREMENTS ONLY ON FILE %s\n", mkwch_pars.save_conf_name);
                printf("\n#################################################\n");

                //--------- MISURA ROBA FERMIONICA ----------------//
                //
             printf("Fermion Measurements: see file %s\n",fm_par.fermionic_outfilename);
             fermion_measures(conf_acc,fermions_parameters,
                        &fm_par, mkwch_pars.residue_metro, id_iter_offset) ;


                //-------------------------------------------------// 
                //--------- MISURA ROBA DI GAUGE ------------------//
                printf("Misure di Gauge:\n");
                plq = calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
                rect = calc_rettangolo_soloopenacc(conf_acc,aux_conf_acc,local_sums);
                poly =  (*polyakov_loop[geom_par.tmap])(conf_acc);//misura polyakov loop

                printf("Plaquette     : %.18lf\n" ,plq/size/3.0/6.0);
                printf("Rectangle     : %.18lf\n" ,rect/size/3.0/6.0/2.0);
                printf("Polyakov Loop : (%.18lf,%.18lf) \n",creal(poly),cimag(poly));


            }else printf("Starting generation of Configurations.\n");

            // THERMALIZATION & METRO    ----   UPDATES //
    
            for(int id_iter=id_iter_offset;id_iter<(mkwch_pars.ntraj+id_iter_offset);id_iter++){

                check_unitarity_device(conf_acc,&max_unitarity_deviation,&avg_unitarity_deviation);
                printf("\tAvg/Max unitarity deviation on device: %e / %e\n", avg_unitarity_deviation, max_unitarity_deviation);
                accettate_therm_old = accettate_therm;
                accettate_metro_old = accettate_metro;
                conf_id_iter++;
                printf("\n#################################################\n");
                printf(  "   GENERATING CONF %d of %d, %dx%dx%dx%d,%1.3f \n",
                        conf_id_iter,mkwch_pars.ntraj+id_iter_offset,
                        geom_par.gnx,geom_par.gny,
                        geom_par.gnz,geom_par.gnt,
                        act_params.beta);
                printf(  "#################################################\n\n");
                //--------- CONF UPDATE ----------------//
                if(id_iter<mkwch_pars.therm_ntraj){
                    accettate_therm = UPDATE_SOLOACC_UNOSTEP_VERSATILE(conf_acc,
                          mkwch_pars.residue_metro,md_parameters.residue_md,id_iter-id_iter_offset,
                            accettate_therm,0);
                }else{
                   accettate_metro = UPDATE_SOLOACC_UNOSTEP_VERSATILE(conf_acc,mkwch_pars.residue_metro,md_parameters.residue_md,id_iter-id_iter_offset-accettate_therm,accettate_metro,1);
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

                FILE *goutfile = fopen(gauge_outfilename,"at");
                if(!goutfile){
                    goutfile = fopen(gauge_outfilename,"wt");
                    strcpy(gauge_outfile_header,"#conf_id\tacc\tplq\trect\tReP\tImP\n");
                    fprintf(goutfile,"%s",gauge_outfile_header);
                }
                if(goutfile){
                    if(id_iter<mkwch_pars.therm_ntraj){
                        printf("Therm_iter %d",conf_id_iter );
                        printf("Placchetta= %.18lf    ", plq/size/6.0/3.0);
                        printf("Rettangolo= %.18lf\n",rect/size/6.0/3.0/2.0);
                    

                    }else printf("Metro_iter %d   Placchetta= %.18lf    Rettangolo= %.18lf\n",conf_id_iter,plq/size/6.0/3.0,rect/size/6.0/3.0/2.0);


                    fprintf(goutfile,"%d\t%d\t",conf_id_iter,
                            accettate_therm+accettate_metro
                            -accettate_therm_old-accettate_metro_old);
                    fprintf(goutfile,"%.18lf\t%.18lf\t%.18lf\t%.18lf\n",
                            plq/size/6.0/3.0,
                            rect/size/6.0/3.0/2.0, 
                            creal(poly), cimag(poly));
                    
                }
                fclose(goutfile);
                //-------------------------------------------------//

                //--------- SALVA LA CONF SU FILE ------------------//
                if(conf_id_iter%mkwch_pars.storeconfinterval==0){
                    char tempname[50];char serial[10];
                    strcpy(tempname,mkwch_pars.store_conf_name);
                    sprintf(serial,".%05d",conf_id_iter);
                    strcat(tempname,serial);
                    printf("Storing conf %s.\n", tempname);
                    recv_loc_subconf_from_buffer(conf_rw,conf_acc,0);
                    save_conf(conf_rw,tempname,conf_id_iter,mkwch_pars.use_ildg);
                }
                if(conf_id_iter%mkwch_pars.saveconfinterval==0){
                    printf("Saving conf %s.\n", mkwch_pars.save_conf_name);
                    recv_loc_subconf_from_buffer(conf_rw,conf_acc,0);
                    save_conf(conf_rw,mkwch_pars.save_conf_name, conf_id_iter,
                            mkwch_pars.use_ildg);
                }

                //-------------------------------------------------//

            }// id_iter loop ends here


            //--------- SALVA LA CONF SU FILE ------------------//

            if(mkwch_pars.ntraj > 0){
                    recv_loc_subconf_from_buffer(conf_rw,conf_acc,0);
                    save_conf(conf_rw,mkwch_pars.save_conf_name, conf_id_iter,
                    mkwch_pars.use_ildg );
            }
            //-------------------------------------------------//


            plq = calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
            topoch = compute_topological_charge(conf_acc,aux_conf_acc,d_local_sums);
            printf("COOL 0  Placchetta= %.18lf  TopCh= %.18lf \n",plq/size/6.0/3.0,topoch);

//
//               for(int icool=0;icool<5000;icool++){
//               cool_conf(conf_acc,aux_conf_acc);
//               plq = calc_plaquette_soloopenacc(conf_acc,aux_conf_acc,local_sums);
//               topoch = compute_topological_charge(conf_acc,aux_conf_acc,d_local_sums);
//               printf("COOL %d  Placchetta= %.18lf  TopCh= %.18lf \n",icool+1,plq/size/6.0/3.0,topoch);
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
    SHUTDOWN_ACC_DEVICE(my_device_type);
    /////////////////////////////////////////////////////////////////////////////////////////////////
#endif

    mem_free();

    return 0;
}

