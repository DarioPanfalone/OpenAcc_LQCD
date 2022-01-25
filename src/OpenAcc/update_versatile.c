// se metro==1 allora fa il test di metropolis
// se metro==0 allora non fa il test di metropolis --> termalizzazione
#ifndef UPDATE_VERSATILE_C_
#define UPDATE_VERSATILE_C_


#include "../DbgTools/dbgtools.h" // useful only for debug
#include "../Include/common_defines.h"
#include "../Include/debug.h"
#include "../Mpi/multidev.h"
#include "../Rand/random.h"
#include "./action.h"
#include "./inverter_package.h"
#include "./alloc_vars.h"
#include "./sp_alloc_vars.h"
#include "./fermionic_utilities.h"
#include "./sp_fermionic_utilities.h"
#include "./find_min_max.h"
#include "./float_double_conv.h"
#include "./inverter_wrappers.h"
#include "./inverter_multishift_full.h"
#include "./io.h"
#include "./md_integrator.h"
#include "./sp_md_integrator.h"
#include "./random_assignement.h"
#include "./rettangoli.h"
#include "./stouting.h"
#include "./su3_measurements.h"
#include "./sp_su3_measurements.h"
#include "./su3_utilities.h"
#include "../tests_and_benchmarks/test_and_benchmarks.h"
#include "./alloc_settings.h"
#include "./update_versatile.h"
#include "./topological_action.h"

#ifdef __GNUC__
#include "sys/time.h"
#endif

#ifdef MULTIDEVICE
#include <mpi.h>
#endif 

action_param act_params; //parametri dell'azione action.h



int UPDATE_SOLOACC_UNOSTEP_VERSATILE(su3_soa *tconf_acc,
        double res_metro, double res_md, int id_iter,int acc,int metro, int max_cg){

    
    
#ifdef STOUT_FERMIONS //Per adesso va bene così, ma dopo aggiungerà la topossibilità       
    su3_soa *tstout_conf_acc_arr = gstout_conf_acc_arr;
    su3_soa_f *tstout_conf_acc_arr_f = gstout_conf_acc_arr_f;
#endif
#ifdef NORANDOM
    printf("MIP%02d: WELCOME! NORANDOM MODE. (UPDATE_SOLOACC_UNOSTEP_VERSATILE())\n",
            devinfo.myrank);
#endif



    printf("MPI%02d: UPDATE_SOLOACC_UNOSTEP_VERSATILE_TLSM_STDFERM: starting... \n",
            devinfo.myrank);
    // DEFINIZIONE DI TUTTI I dt NECESSARI PER L'INTEGRATORE OMELYAN
    
    int iterazioni = id_iter+1;
    double dt_tot;
    double dt_saveoldconf;
    double dt_mdsetup;
    double dt_metropolis;
    gauge_mdtimes = gauge_mdtimes0;
    
//Just Debugging diagnostics
    if(debug_settings.save_diagnostics == 1){
        FILE *foutfile = 
            fopen(debug_settings.diagnostics_filename,"at");
        fprintf(foutfile,"\nIteration %d \n",id_iter);
        fclose(foutfile);
    }


    int mu;
    double **minmaxeig;  //matrix of doubles(actually addresses ) Of Course it has to be allocated.
    minmaxeig = (double**)malloc(alloc_info.NDiffFlavs*sizeof(double*));
    for(int iflav = 0 ; iflav < alloc_info.NDiffFlavs ; iflav++)
        minmaxeig[iflav] = (double*) malloc(2*sizeof(double));
    

    double p1; //definition of p1, p2. Acceptance rates
    double p2;
    int accettata;
    double delta_S; //delta_S
    double action_in,action_fin,action_mom_in,action_mom_fin,action_ferm_in,action_ferm_fin;
	double action_topo_in,action_topo_fin;

    struct timeval t_start, t_saved,t_end;
    gettimeofday ( &t_start, NULL ); 
    if(metro==1){
        // store old conf   set_su3_soa_to_su3_soa(arg1,arg2) ===>   arg2=arg1;
        set_su3_soa_to_su3_soa(tconf_acc,conf_acc_bkp);
        if(verbosity_lv > 1) printf("Backup copy of the initial gauge conf : OK \n");

    }



    gettimeofday ( &t_saved, NULL );

    // ESTRAZIONI RANDOM //momenta are randomly generated.
    if(debug_settings.do_norandom_test){ // NORANDOM
        printf("NORANDOM mode, loading momenta from memory.\n");
        if(read_thmat_soa(momenta,"momenta_norndtest")){
            printf("GENERATING MOMENTA FILE FOR YOUR CONVENIENCE, RE-RUN THIS TEST\n");
            generate_Momenta_gauss(momenta);
            print_thmat_soa(momenta,"momenta_norndtest");
        }
    }
    else generate_Momenta_gauss(momenta); // NORMAL, RANDOM

    printf("MPI%02d - Momenta generated/read : OK \n", devinfo.myrank);

    //    read_thmat_soa(momenta,"momenta");
#pragma acc update device(momenta[0:8])
    if(debug_settings.do_reversibility_test)
        copy_momenta_into_old(momenta,momenta_backup);


    //HERE A FERMIONIC CYCLE:
    for(int iflav = 0 ; iflav < alloc_info.NDiffFlavs ; iflav++){//iflav start
        for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
            int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;

            if(debug_settings.do_norandom_test){ // NORANDOM
                char psferm_filename[20];

                char ps_index_str[5];
                sprintf(ps_index_str,"%d",ps_index);
                strcpy(psferm_filename,"fermion_norndtest");
                strcat(psferm_filename,ps_index_str);
                if(read_vec3_soa(&ferm_phi_acc[ps_index],psferm_filename)){
                    printf("GENERATING FERMION FILE FOR YOUR CONVENIENCE, RE-RUN THIS TEST\n");
                    generate_vec3_soa_gauss(&ferm_phi_acc[ps_index]);
                    print_vec3_soa(&ferm_phi_acc[ps_index],psferm_filename);

                }
            }
            else{    // NORMAL, RANDOM
                if(verbosity_lv > 3 )printf("Ferm generation (flav=%d,ps=%d,ps_index=%d) : OK \n",iflav,ips,ps_index);
                generate_vec3_soa_gauss(&ferm_phi_acc[ps_index]);
            }
        }
    }// end for iflav
#pragma acc update device(ferm_phi_acc[0:alloc_info.NPS_tot])

    gconf_as_fermionmatrix_f = conf_acc_f;
#ifdef STOUT_FERMIONS  //STOUT FERM DIRECTIVE.
    // DILATION USING STOUTED DIRAC OPERATOR
    // STOUTING...(ALREADY ON DEVICE)
    if ( (act_params.stout_steps > 0 ) && (alloc_info.NDiffFlavs > 0) ){
        stout_wrapper(tconf_acc,tstout_conf_acc_arr);
        gconf_as_fermionmatrix = 
            &(tstout_conf_acc_arr[8*(act_params.stout_steps-1)]);

        if(md_parameters.singlePrecMD)
            convert_double_to_float_su3_soa(gconf_as_fermionmatrix,gconf_as_fermionmatrix_f);
    }
    else gconf_as_fermionmatrix = tconf_acc;
#else
    gconf_as_fermionmatrix = tconf_acc;
    if(md_parameters.singlePrecMD)
        convert_double_to_float_su3_soa(gconf_as_fermionmatrix,gconf_as_fermionmatrix_f);
#endif

    // DILATION OF FIRST_INV RATIONAL APPROXIMATION
    for(int iflav = 0 ; iflav < alloc_info.NDiffFlavs ; iflav++){ //start iflav for
        if(verbosity_lv > 2 ) printf("Rat approx rescale (flav=%d)\n",iflav);

        if(debug_settings.do_norandom_test){ // NORANDOM
            if(read_vec3_soa(kloc_p,"kloc_p_norndtest")){
                generate_vec3_soa_gauss(kloc_p);
                print_vec3_soa(kloc_p,"kloc_p_norndtest");
                printf("GENERATED kloc_p_norndtest FOR NORANDOM TEST, RE-RUN THIS TEST\n");
            }
            if(read_vec3_soa(kloc_s,"kloc_s_norndtest")){
                generate_vec3_soa_gauss(kloc_s);
                print_vec3_soa(kloc_s,"kloc_s_norndtest");
                printf("GENERATED kloc_s_norndtest FOR NORANDOM TEST, RE-RUN THIS TEST\n");
            }
        }
        else{   // NORMAL, RANDOM
            // generate gauss-randomly the fermion kloc_p that will be used in the computation of the max eigenvalue
            generate_vec3_soa_gauss(kloc_p);
            generate_vec3_soa_gauss(kloc_s);
        }
        // update the fermion kloc_p copying it from the host to the device
#pragma acc update device(kloc_p[0:1])
#pragma acc update device(kloc_s[0:1])
        // USING STOUTED GAUGE MATRIX
        //printf("    before min and max eig comp : OK \n");
        find_min_max_eigenvalue_soloopenacc(gconf_as_fermionmatrix,&(fermions_parameters[iflav]),kloc_r,kloc_h,kloc_p,kloc_s,minmaxeig[iflav]);
        if(verbosity_lv > 3 ) printf("    find min and max eig : OK \n");
        RationalApprox *approx_fi = &(fermions_parameters[iflav].approx_fi);
        RationalApprox *approx_fi_mother = &(fermions_parameters[iflav].approx_fi_mother);
        rescale_rational_approximation(approx_fi_mother,approx_fi,minmaxeig[iflav]);
        if(verbosity_lv > 4 ) printf("    rat approx rescaled : OK \n\n");

#pragma acc update device(approx_fi[0:1])
    }//end for iflav

//it goes into this for if metro==1, that means that metropolis will be executed.
    if(metro==1){
        /////////////// INITIAL ACTION COMPUTATION ////////////////////////////////////////////
        
        if(GAUGE_ACTION == 0)// Standard gauge action 
            action_in = BETA_BY_THREE*calc_plaquette_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
        if(GAUGE_ACTION == 1){ //Tlsym gauge action
            action_in = C_ZERO * BETA_BY_THREE * calc_plaquette_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
            action_in += C_ONE * BETA_BY_THREE * calc_rettangolo_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
        }
		
	if(verbosity_lv>3 && act_params.topo_action==1)printf("\tComputing topological action\n");

#ifdef STOUT_TOPO
	action_topo_in = (act_params.topo_action==0)? 0.0 : compute_topo_action(tconf_acc,tstout_conf_acc_arr);

#ifdef STOUT_FERMIONS 
    // DILATION USING STOUTED DIRAC OPERATOR
    // STOUTING...(ALREADY ON DEVICE)
    if ( (act_params.stout_steps > 0) && (alloc_info.NDiffFlavs > 0) ){
        stout_wrapper(tconf_acc,tstout_conf_acc_arr);
        gconf_as_fermionmatrix = 
            &(tstout_conf_acc_arr[8*(act_params.stout_steps-1)]);

        if(md_parameters.singlePrecMD)
            convert_double_to_float_su3_soa(gconf_as_fermionmatrix,gconf_as_fermionmatrix_f);
    }
    else gconf_as_fermionmatrix = tconf_acc;
#endif
#else
	action_topo_in = (act_params.topo_action==0)? 0.0 : compute_topo_action(tconf_acc);
#endif

        action_mom_in = 0.0;
        for(mu =0;mu<8;mu++)  action_mom_in += calc_momenta_action(momenta,d_local_sums,mu);//calculating momenta action.
        
        
        action_ferm_in=0;
        for(int iflav = 0 ; iflav < alloc_info.NDiffFlavs ; iflav++){ //iflav start
	  for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
	    
	    int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;
	    action_ferm_in += real_scal_prod_global(&ferm_phi_acc[ps_index],&ferm_phi_acc[ps_index]);
	  }
        }// end for iflav
        ///////////////////////////////////////////////////////////////////////////////////////
        printf("MPI%02d - Initial Action Computed : OK \n", devinfo.myrank);
    }
    
    
    multishift_invert_iterations = 0 ; 
    // FIRST INV APPROX CALC --> calculation of CHI fermion
    
    inverter_package ip;
    setup_inverter_package_sp(&ip,gconf_as_fermionmatrix_f,k_p_shiftferm_f,alloc_info.maxApproxOrder,
			      kloc_r_f,kloc_h_f,kloc_s_f,kloc_p_f,aux1_f);  
    setup_inverter_package_dp(&ip,gconf_as_fermionmatrix,  k_p_shiftferm,  alloc_info.maxApproxOrder,
			      kloc_r,  kloc_h,  kloc_s,  kloc_p);  



    for(int iflav = 0 ; iflav < alloc_info.NDiffFlavs ; iflav++){ //for iflav.
        for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
            if(0==devinfo.myrank)
                printf("Calculation of chi for fermion %d, copy %d\n", iflav,ips);

            int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;
            // USING STOUTED GAUGE MATRIX

            inverter_multishift_wrapper(ip, &fermions_parameters[iflav], 
                    &(fermions_parameters[iflav].approx_fi), ferm_shiftmulti_acc,
                    &(ferm_phi_acc[ps_index]), res_metro,max_cg, CONVERGENCE_CRITICAL);
            if(0==devinfo.myrank) printf("Inversion performed for fermion %d, copy %d\n", iflav,ips);
            recombine_shifted_vec3_to_vec3(ferm_shiftmulti_acc, &(ferm_phi_acc[ps_index]), &(ferm_chi_acc[ps_index]),&(fermions_parameters[iflav].approx_fi));
            if(0==devinfo.myrank) printf("Calculated chi for fermion %d, copy %d\n", iflav,ips);


        }
    }// end for iflav
    if(verbosity_lv > 3) printf("MPI%02d: Computed the fermion CHI : OK \n", devinfo.myrank);

    // DILATION OF MOLECULAR DYNAMICS RATIONAL APPROXIMATION
    for(int iflav = 0 ; iflav < alloc_info.NDiffFlavs ; iflav++){
        // recovering eigenvalues from approx_fi
        RationalApprox *approx_md = &(fermions_parameters[iflav].approx_md);
        RationalApprox *approx_md_mother = &(fermions_parameters[iflav].approx_md_mother);
        rescale_rational_approximation(approx_md_mother,approx_md,minmaxeig[iflav]);
        printf("MPI%02d: Rescaled Rational Approximation for flavour %d\n", devinfo.myrank, iflav);
#pragma acc update device(approx_md[0:1])

    }//end for iflav

    // DINAMICA MOLECOLARE (stouting implicitamente usato in calcolo forza fermionica)
    struct timeval md_start, md_end;
    int multishift_iterations_before_md = multishift_invert_iterations ;  
    gettimeofday(&md_start,NULL);
    if(1 == md_parameters.singlePrecMD){
        if(verbosity_lv > 1) 
            printf("MPI%02d: SINGLE PRECISION MOLECULAR DYNAMICS...\n", devinfo.myrank);

        // conversion double to float

        su3_soa_f * tconf_acc_f = conf_acc_f;

        printf("Converting mommenta...\n");
        convert_double_to_float_thmat_soa(momenta,momenta_f);
        double act_mom_check_f = 0;
        for(mu =0;mu<8;mu++)  act_mom_check_f += 
            calc_momenta_action_f(momenta_f,d_local_sums_f,mu);
        double act_mom_check = 0; // should not be necessary
        for(mu =0;mu<8;mu++)  act_mom_check += 
            calc_momenta_action(momenta,d_local_sums,mu);

        printf("Converting conf...\n");
        convert_double_to_float_su3_soa(tconf_acc,tconf_acc_f);
        double act_links_check_f = BETA_BY_THREE* calc_plaquette_soloopenacc_f(
                tconf_acc_f, aux_conf_acc_f, local_sums_f);
        double act_links_check = BETA_BY_THREE*calc_plaquette_soloopenacc(tconf_acc,
                aux_conf_acc,local_sums); // should not be necessary

        double act_ferm_check_f = 0;
        double act_ferm_check = 0; // should not be necessary
        int ips;
        for(ips = 0; ips < alloc_info.NPS_tot;ips++){
            printf("Converting ferm_chi_acc[%d]...\n",ips);
            convert_double_to_float_vec3_soa(&ferm_chi_acc[ips],&ferm_chi_acc_f[ips]);
            act_ferm_check += real_scal_prod_global(&ferm_phi_acc[ips],
                    &ferm_phi_acc[ips]);
            act_ferm_check_f += real_scal_prod_global_f(&ferm_phi_acc_f[ips],
                    &ferm_phi_acc_f[ips]);
        }

        if(verbosity_lv>1){
            printf("MPI%02d: DOUBLE->SINGLE PRECISION conversion done.\n", 
                    devinfo.myrank);
            if(verbosity_lv>3){
                printf("MPI%02d: Mom action  (single/double precision): %lf / %lf \n",
                        devinfo.myrank,act_mom_check_f,act_mom_check );
                printf("MPI%02d: Link Action (single/double precision): %lf / %lf \n",
                        devinfo.myrank,act_links_check_f,act_links_check );
                printf("MPI%02d: Ferm Action (single/double precision): %lf / %lf \n",
                        devinfo.myrank,act_ferm_check_f,act_ferm_check );
            }
        }

        multistep_2MN_SOLOOPENACC_f(ipdot_acc_f,tconf_acc_f,
#if (defined STOUT_FERMIONS) || (defined STOUT_TOPO)
                tstout_conf_acc_arr_f,
#endif
                auxbis_conf_acc_f, // globale
                aux_conf_acc_f,fermions_parameters,alloc_info.NDiffFlavs,
                ferm_chi_acc_f,ferm_shiftmulti_acc_f,
                ip,
                momenta_f,local_sums_f,res_md,max_cg);

        if(verbosity_lv > 1) printf("MPI%02d: Single Precision Molecular Dynamics Completed \n",devinfo.myrank );

        convert_float_to_double_thmat_soa(momenta_f,momenta);
        convert_float_to_double_su3_soa(tconf_acc_f,tconf_acc);

        if(verbosity_lv > 1) printf("MPI%02d: SINGLE -> DOUBLE PRECISION conversion done.\n",devinfo.myrank );

    } 
    else{

        printf("DOUBLE PRECISION MOLECULAR DYNAMICS...\n");

        multistep_2MN_SOLOOPENACC(ipdot_acc,tconf_acc,
#if (defined STOUT_FERMIONS) || (defined STOUT_TOPO)
                tstout_conf_acc_arr,
#endif
                auxbis_conf_acc, // globale
                aux_conf_acc,fermions_parameters,alloc_info.NDiffFlavs,
                ferm_chi_acc,ferm_shiftmulti_acc,
                ip,
                momenta,local_sums,res_md,max_cg);


    }
    gettimeofday(&md_end,NULL);


    int multishift_iterations_after_md = multishift_invert_iterations ; 
    int multishift_iterations_md = multishift_iterations_after_md - 
        multishift_iterations_before_md;

    double mdtime = (double)(md_end.tv_sec - md_start.tv_sec) + 
        ((double)(md_end.tv_usec - md_start.tv_usec)/1.0e6);



    if(devinfo.myrank == 0){

        printf("MPI%02d:Time needed for moledular dynamics: %f s, ",devinfo.myrank,mdtime);
        printf("total CG-M iterations during MD: %d\n",multishift_iterations_md );
    }

    if(debug_settings.do_reversibility_test){

        printf("MPI%02d: PERFORMING REVERSIBILITY TEST, DOUBLE PRECISION.\n", devinfo.myrank);
        printf("MPI%02d: Inverting momenta.\n", devinfo.myrank);

        invert_momenta(momenta);

        multistep_2MN_SOLOOPENACC(ipdot_acc,tconf_acc,
#if (defined STOUT_FERMIONS) || (defined STOUT_TOPO)
                tstout_conf_acc_arr,
#endif
                auxbis_conf_acc, // globale
                aux_conf_acc,fermions_parameters,alloc_info.NDiffFlavs,
                ferm_chi_acc,ferm_shiftmulti_acc,
                ip,
                momenta,local_sums,res_md, max_cg);

#pragma acc update device(conf_acc_bkp[0:8])
        double conf_error =  calc_diff_su3_soa_norm(tconf_acc,conf_acc_bkp);
        // Note: sum instead of difference must be calculated because momenta 
        // were inverted.
        double momenta_error = calc_sum_momenta_norm(momenta,momenta_backup);

        printf("MPI%02d: Reversibility test: Conf_error: %e , momenta_error %e \n",
                devinfo.myrank, conf_error, momenta_error);


        printf("MPI%02d: Since a reversibility test is being performed, no need to  go on. Exiting now.\n", 
                devinfo.myrank);

        save_conf_wrapper(tconf_acc,"conf_REVTEST",0,0); // conf_id_iter =0, not using ILDG
        save_conf_wrapper(conf_acc_bkp,"conf_REVTEST_bkp",0,0);


#ifdef MULTIDEVICE
        MPI_Finalize();
#endif
        mem_free_core();
        mem_free_extended();
        exit(0);

    }



    if(verbosity_lv > 1) printf("MPI%02d - MOLECULAR DYNAMICS COMPLETED \n", 
            devinfo.myrank);


#ifdef STOUT_FERMIONS
    // STOUTING...(ALREADY ON DEVICE)
    if ( (act_params.stout_steps > 0) && (alloc_info.NDiffFlavs > 0) ){
        stout_wrapper(tconf_acc,tstout_conf_acc_arr);
        gconf_as_fermionmatrix = 
            &(tstout_conf_acc_arr[8*(act_params.stout_steps-1)]);
    }
    else gconf_as_fermionmatrix = tconf_acc;
#else
    gconf_as_fermionmatrix = tconf_acc;
#endif

    if(metro==1){
        // DILATION OF LAST_INV RATIONAL APPROX
        for(int iflav = 0 ; iflav < alloc_info.NDiffFlavs ; iflav++){
            // generate gauss-randomly the fermion kloc_p that will be used in the computation of the max eigenvalue
            if(debug_settings.do_norandom_test){ // NORANDOM
                if(read_vec3_soa(kloc_p,"kloc_p_norndtest")){
                    generate_vec3_soa_gauss(kloc_p);
                    print_vec3_soa(kloc_p,"kloc_p_norndtest");
                    printf("GENERATED kloc_p_norndtest FOR NORANDOM TEST, RE-RUN THIS TEST\n");
                }
                if(read_vec3_soa(kloc_s,"kloc_s_norndtest")){
                    generate_vec3_soa_gauss(kloc_s);
                    print_vec3_soa(kloc_s,"kloc_s_norndtest");
                    printf("GENERATED kloc_s_norndtest FOR NORANDOM TEST, RE-RUN THIS TEST\n");
                }
            }
            else{    // NORMAL,/ RANDOM
                generate_vec3_soa_gauss(kloc_p);
                generate_vec3_soa_gauss(kloc_s);
            }
            // update the fermion kloc_p copying it from the host to the device
#pragma acc update device(kloc_p[0:1])
#pragma acc update device(kloc_s[0:1])
            // USING STOUTED CONF
            find_min_max_eigenvalue_soloopenacc(gconf_as_fermionmatrix,&(fermions_parameters[iflav]),kloc_r,kloc_h,kloc_p,kloc_s,minmaxeig[iflav]);
            //#pragma acc update device(minmaxeig[0:2])
            RationalApprox *approx_li = &(fermions_parameters[iflav].approx_li);
            RationalApprox *approx_li_mother = &(fermions_parameters[iflav].approx_li_mother);
            rescale_rational_approximation(approx_li_mother,approx_li,minmaxeig[iflav]);
#pragma acc update device(approx_li[0:1])
        }//end iflav.

       // LAST INV APPROX CALC 
       
        setup_inverter_package_dp(&ip,gconf_as_fermionmatrix,  k_p_shiftferm,  alloc_info.maxApproxOrder,
                kloc_r,  kloc_h,  kloc_s,  kloc_p);  
        if(md_parameters.singlePrecMD){
            convert_double_to_float_su3_soa(gconf_as_fermionmatrix,
                    gconf_as_fermionmatrix_f);
            setup_inverter_package_sp(&ip,gconf_as_fermionmatrix_f,k_p_shiftferm_f,
                    alloc_info.maxApproxOrder,kloc_r_f,kloc_h_f,kloc_s_f,kloc_p_f,aux1_f);        }
        for(int iflav = 0 ; iflav < alloc_info.NDiffFlavs ; iflav++){ //iflav
            for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
                int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;
                // USING STOUTED CONF
                inverter_multishift_wrapper(ip, &fermions_parameters[iflav], 
                        &(fermions_parameters[iflav].approx_li),
                        ferm_shiftmulti_acc, &(ferm_chi_acc[ps_index]), res_metro,max_cg,
                        CONVERGENCE_CRITICAL);
                recombine_shifted_vec3_to_vec3(ferm_shiftmulti_acc, &(ferm_chi_acc[ps_index]), &(ferm_phi_acc[ps_index]),&(fermions_parameters[iflav].approx_li));
            }
        }//end iflav
        printf(" MPI%02d - Final Action Computed : OK \n", devinfo.myrank);

        //I should not modify the action computation...well actually yes. I Should modify the way action is computed.
        //I have only to change the gauge action
        
        ///////////////   FINAL ACTION COMPUTATION  ////////////////////////////////////////////
        if(GAUGE_ACTION == 0) // Standard gauge action
            action_fin = BETA_BY_THREE * calc_plaquette_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
        if(GAUGE_ACTION == 1){ // Tlsym gauge action
            action_fin = C_ZERO * BETA_BY_THREE * calc_plaquette_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
            action_fin += C_ONE * BETA_BY_THREE * calc_rettangolo_soloopenacc(tconf_acc,aux_conf_acc,local_sums);
        }
	
#ifdef STOUT_TOPO
	action_topo_fin = (act_params.topo_action==0)? 0.0 : compute_topo_action(tconf_acc,tstout_conf_acc_arr);
	
#ifdef STOUT_FERMIONS
    // STOUTING...(ALREADY ON DEVICE)
    if ( (act_params.stout_steps > 0) && (alloc_info.NDiffFlavs > 0) ){
        stout_wrapper(tconf_acc,tstout_conf_acc_arr);
        gconf_as_fermionmatrix = 
            &(tstout_conf_acc_arr[8*(act_params.stout_steps-1)]);
    }
    else gconf_as_fermionmatrix = tconf_acc;
#else
    gconf_as_fermionmatrix = tconf_acc;
#endif
#else
	action_topo_fin = (act_params.topo_action==0)? 0.0 : compute_topo_action(tconf_acc);
#endif
        action_mom_fin = 0.0;
        //action_mom doesn't change by adding the defect.
        for(mu =0;mu<8;mu++)    action_mom_fin += calc_momenta_action(momenta,d_local_sums,mu);

        action_ferm_fin=0;
        for(int iflav = 0 ; iflav < alloc_info.NDiffFlavs ; iflav++){
            for(int ips = 0 ; ips < fermions_parameters[iflav].number_of_ps ; ips++){
                int ps_index = fermions_parameters[iflav].index_of_the_first_ps + ips;
                action_ferm_fin += real_scal_prod_global(&ferm_chi_acc[ps_index],&ferm_phi_acc[ps_index]);
            }
        } // end for iflav


        ////////////////////////////////////////////////////////////////////////////////////////
        
        
        //******************************************************************************************//
        //THE METROPOLIS STEP:
        //******************************************************************************************//
        
        
        // delta_S = action_new - action_old
        delta_S  = - (-action_in+action_mom_in+action_ferm_in+action_topo_in) + (-action_fin+action_mom_fin+action_ferm_fin+action_topo_fin);
        if(verbosity_lv > 2 && 0 == devinfo.myrank ){
            printf("MPI%02d-iterazione %i:  Gauge_ACTION  (in and out) = %.18lf , %.18lf\n",
                    devinfo.myrank,iterazioni,-action_in,-action_fin);
			printf("MPI%02d-iterazione %i:  Topol_ACTION  (in and out) = %.18lf , %.18lf\n",
                    devinfo.myrank,iterazioni,+action_topo_in,+action_topo_fin);
            printf("MPI%02d-iterazione %i:  Momen_ACTION  (in and out) = %.18lf , %.18lf\n",
                    devinfo.myrank,iterazioni,action_mom_in,action_mom_fin);
            printf("MPI%02d-iterazione %i:  Fermi_ACTION  (in and out) = %.18lf , %.18lf\n",
                    devinfo.myrank,iterazioni,action_ferm_in,action_ferm_fin);
        }       
        printf("MPI%02d-iterazione %i:  DELTA_ACTION = %.18lf. ",
                devinfo.myrank,iterazioni,delta_S);


        if(debug_settings.do_norandom_test)  // NORANDOM
            printf("Always accept in NORANDOM MODE!!!\n");

        if(delta_S<0){
            accettata=1;
        }
        else
        {
            p1=exp(-delta_S);
            if(debug_settings.do_norandom_test) p2=0; // NORANDOM
            else{   // NORMAL, RANDOM
                if(0==devinfo.myrank)p2=casuale();
#ifdef MULTIDEVICE
                MPI_Bcast((void*) &p2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
                printf("MPI%02d p2 : %f, p1 %f \n",devinfo.myrank, p2,p1);
#endif
            }
            if(p2<p1)
            {
                accettata=1;
            }
            else
            {
                accettata=0;
                // configuration reject
            }
        }
    }


    if(metro==1){
        if(accettata==1){
            acc++;
            printf("MPI%02d,ACCEPTED   ---> [acc/iter] = [%i/%i] \n",
                    devinfo.myrank,acc,iterazioni);
            // configuration accepted   set_su3_soa_to_su3_soa(arg1,arg2) ===>   arg2=arg1;
            set_su3_soa_to_su3_soa(tconf_acc,conf_acc_bkp);
        }else{
            printf("MPI%02d,REJECTED   ---> [acc/iter] = [%i/%i] \n",
                    devinfo.myrank,acc,iterazioni);
            // configuration rejected   set_su3_soa_to_su3_soa(arg1,arg2) ===>   arg2=arg1;
            set_su3_soa_to_su3_soa(conf_acc_bkp,tconf_acc);
#pragma acc update device(tconf_acc[0:8])
            // sul device aggiorniamo la conf rimettendo quella del passo precedente
        }
    }
    gettimeofday ( &t_end, NULL );


    if(metro==0){// accetta sempre in fase di termalizzazione
        acc++;
    }


    dt_tot = (double)(t_end.tv_sec - t_start.tv_sec) + ((double)(t_end.tv_usec - t_start.tv_usec)/1.0e6);
    dt_saveoldconf = (double)(t_saved.tv_sec - t_start.tv_sec) + ((double)(t_saved.tv_usec - t_start.tv_usec)/1.0e6);
    dt_mdsetup = (double)(md_start.tv_sec - t_saved.tv_sec) + ((double)(md_start.tv_usec - t_saved.tv_usec)/1.0e6);
    dt_metropolis = (double)(t_end.tv_sec - md_end.tv_sec) + ((double)(t_end.tv_usec - md_end.tv_usec)/1.0e6);

    if(0==devinfo.myrank){

        printf("   FULL UPDATE COMPUTATION TIME ");
        if(metro==0)printf("NO");
        else printf("SI");

        printf("METRO - Tot Update time : %f sec \n",dt_tot);
        printf("\t\tTime to save old conf    : %f sec\n",dt_saveoldconf);
        printf("\t\tTime for MD setup        : %f sec\n",dt_mdsetup);
        printf("\t\tTime for MD              : %f sec\n",mdtime);
        printf("\t\tTime for Metropolis test : %f sec\n",dt_metropolis);
        printf("\t\tTotal CG-M iterations    : %d\n",multishift_invert_iterations);
        printf("\t\t      CG-M iterations[MD]: %d\n",multishift_iterations_md);
        printf("\t\t      CG-M iterations[FI]: %d\n",multishift_iterations_before_md);
        printf("\t\t      CG-M iterations[LI]: %d\n",multishift_invert_iterations-
                multishift_iterations_after_md);

    }
    if(1 == debug_settings.save_diagnostics && 0 == devinfo.myrank){
        FILE *foutfile = 
            fopen(debug_settings.diagnostics_filename,"at");

        if(metro==1){
            fprintf(foutfile,"GAS %.18lf GAF %.18lf \t",-action_in,-action_fin);
            fprintf(foutfile,"MAS %.18lf MAF %.18lf \t",action_mom_in,action_mom_fin);
            fprintf(foutfile,"FAS %.18lf FAF %.18lf \t",action_ferm_in,action_ferm_fin);
            fprintf(foutfile," D %.18lf\n",delta_S);
        }else{

            fprintf(foutfile," THERM_ITERATION ");

        }
        fprintf(foutfile,"TOTUPDATETIME  %f \n",dt_tot);
        fprintf(foutfile,"OLDCONFSAVETIME  %f \n",dt_saveoldconf);
        fprintf(foutfile,"MDSETUPTIME    %f \n",dt_mdsetup);
        fprintf(foutfile,"MDTIME %f (STEPS: %d , GAUGESTEPS: %d )\n",mdtime,md_parameters.no_md,md_parameters.gauge_scale);
        fprintf(foutfile,"METROTIME %f \n",dt_metropolis);
        fprintf(foutfile,"CG-M-iters[TOT] %d \n",multishift_invert_iterations);
        fprintf(foutfile,"CG-M-iters[MD]  %d \n",multishift_iterations_md);
        fprintf(foutfile,"CG-M-iters[FI]  %d \n",multishift_iterations_before_md);
        fprintf(foutfile,"CG-M-iters[LI]  %d \n",multishift_invert_iterations-
                multishift_iterations_after_md);
        if(metro==1) fprintf(foutfile,"ACCEPTED  %d\n",accettata); 

        fprintf(foutfile,"GIPDSTIME %e\n",gauge_mdtimes.calcIpdotTimeBorder      / gauge_mdtimes.count );
        fprintf(foutfile,"GIPDBTIME %e\n",gauge_mdtimes.calcIpdotTimeBulk        / gauge_mdtimes.count );
        fprintf(foutfile,"GMSMSTIME %e\n",gauge_mdtimes.momSumMultTimeBorder     / gauge_mdtimes.count );
        fprintf(foutfile,"GMSMBTIME %e\n",gauge_mdtimes.momSumMultTimeBulk       / gauge_mdtimes.count );
        fprintf(foutfile,"GMETSTIME %e\n",gauge_mdtimes.momExpTimesConfTimeBorder/ gauge_mdtimes.count );
        fprintf(foutfile,"GMETBTIME %e\n",gauge_mdtimes.momExpTimesConfTimeBulk  / gauge_mdtimes.count );
        fprintf(foutfile,"GCOSTTIME %e\n",gauge_mdtimes.communicationsStartTime  / gauge_mdtimes.count );
        fprintf(foutfile,"GCOMMTIME %e\n",gauge_mdtimes.communicationsTime       / gauge_mdtimes.count );
        gaugeMdCountersReset(&gauge_mdtimes);
#ifdef MULTIDEVICE
        if( ! devinfo.async_comm_fermion){
            fprintf(foutfile,"FCOMMTIME %e\n", dirac_times.totTransferTime       / dirac_times.count );

        diracCountersReset(&dirac_times);
        }
#endif



        fclose(foutfile);
    }
    for(int iflav = 0 ; iflav < alloc_info.NDiffFlavs ; iflav++)
        free(minmaxeig[iflav]);
    free(minmaxeig);

    return acc;

}

#endif
