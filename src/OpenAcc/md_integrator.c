// MULTISTEP VERSION of minimum_norm2_B

// 2nd order Minimum Norm integrator (2MN)
// See reference hep-lat/0505020 Takaishi, De Forcrand
//
// Scheme (needs three force calculations per step):
// 2MN(dt) = exp (-l * dt * dS/dq) exp (dt/2 * p) * ...
//           exp (-(1- 2l) * dt * dS/dq) exp (dt/2 * p)  exp (-l * dt * dS/dq)
// p       : momenta
// - dS/sq : force contibution
// dt      : integration step
// l       : lambda parameter (1/6 Sexton Weingarten, 0.193... Omelyan et al)
//    
// Total scheme:
// [2MN]^N
// N       : number of Molecular Dynamics steps
//
// See reference hep-lat/0506011 Urbach et al. for the multiple time scale 
//
// p->p+a dS/dq=p+ia(-i dS/dq)=p+ia*ipdot

// the various dt involved are stored into the delta array which is allocated into UPDATE_ACC(args)

#ifndef MD_INTEGRATOR_C
#define MD_INTEGRATOR_C

#include "../Include/common_defines.h"
#include "./struct_c_def.h"
#include "./fermion_force.h"
#include "./md_integrator.h"
#include "./alloc_vars.h"
#include "./ipdot_gauge.h"
#include "./su3_utilities.h"
#include "../Include/common_defines.h"
#include "../Include/fermion_parameters.h"
#include "./action.h"
#include "../Mpi/multidev.h"

#ifdef MULTIDEVICE
#include "../Mpi/communications.h"
#endif

md_param md_parameters;

double deltas_Omelyan[7];


#ifdef DEBUG_MD
#include "../DbgTools/dbgtools.h"
int already_printed_debug = 0;
#endif

void initialize_md_global_variables(md_param md_params )
{
    int no_md = md_params.no_md;// number of MD steps
    int gauge_scale = md_params.gauge_scale;   // Update fermions every gauge_scale gauge updates
    double t = md_params.t ;

    double epsilon_acc;
    d_complex ieps_acc,iepsh_acc;


    epsilon_acc = t/((double)(no_md));
    ieps_acc  = 0.0 + (epsilon_acc) * 1.0I;
    iepsh_acc = 0.0 + (epsilon_acc) * 0.5 * 1.0I;

    const double lambda=0.1931833275037836; // Omelyan Et Al.
    //const double lambda=0.1931833; // Omelyan Et Al.
    const double gs=t*0.5/(double) gauge_scale;

    deltas_Omelyan[0]= -cimag(ieps_acc) * lambda;
    deltas_Omelyan[1]= -cimag(ieps_acc) * (1.0-2.0*lambda);
    deltas_Omelyan[2]= -cimag(ieps_acc) * 2.0*lambda;
    deltas_Omelyan[3]= -cimag(ieps_acc) * gs*lambda * BETA_BY_THREE;
    deltas_Omelyan[4]=  cimag(iepsh_acc)* gs;
    deltas_Omelyan[5]= -cimag(ieps_acc) * gs*(1.0-2.0*lambda)*BETA_BY_THREE;
    deltas_Omelyan[6]= -cimag(ieps_acc) * gs*2.0*lambda*BETA_BY_THREE;

}



#ifdef MULTIDEVICE
void multistep_2MN_gauge_async_bloc(su3_soa *tconf_acc,su3_soa *local_staples,
        tamat_soa *tipdot,thmat_soa *tmomenta, int omelyan_index){
    MPI_Request send_border_requests[96]; 
    MPI_Request recv_border_requests[96];
    if(verbosity_lv > 3) printf("MPI%02d - In async bloc - Index %d\n",
            devinfo.myrank, omelyan_index);

/* // At present, useless!
    calc_ipdot_gauge_soloopenacc_d3c(tconf_acc,local_staples,tipdot,
            HALO_WIDTH,GAUGE_HALO); 
    calc_ipdot_gauge_soloopenacc_d3c(tconf_acc,local_staples,tipdot,
            nd3-HALO_WIDTH-GAUGE_HALO,GAUGE_HALO); 
    calc_ipdot_gauge_soloopenacc_bulk(tconf_acc,local_staples,tipdot);
*/
    // ISSUE : calc_ipdot_gauge_bulk depends on tconf_acc in the surface.
    // mom_exp_times_conf_soloopenacc_d3c() is goig to modify it. 
    // So, calc_ipdot_gauge_soloopenacc_bulk must be done before 
    // mom_exp_times_conf_soloopenacc_d3c(), unless we find another 
    // clever way to fix it, like using another conf in 
    // mom_exp_times_conf_soloopenacc().

    calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);

    mom_sum_mult_d3c(tmomenta,tipdot,deltas_Omelyan,omelyan_index,
            HALO_WIDTH,GAUGE_HALO);
    mom_sum_mult_d3c(tmomenta,tipdot,deltas_Omelyan,omelyan_index,
            nd3-HALO_WIDTH-GAUGE_HALO,GAUGE_HALO); 


    mom_exp_times_conf_soloopenacc_d3c(tconf_acc,tmomenta,
            deltas_Omelyan,4,
            HALO_WIDTH,GAUGE_HALO);
    mom_exp_times_conf_soloopenacc_d3c(tconf_acc,tmomenta,
            deltas_Omelyan,4,
            nd3-HALO_WIDTH-GAUGE_HALO,GAUGE_HALO); 

    communicate_su3_borders_async(tconf_acc,GAUGE_HALO,
            send_border_requests,recv_border_requests);

    mom_sum_mult_bulk(tmomenta,tipdot,deltas_Omelyan,omelyan_index);

    mom_exp_times_conf_soloopenacc_bulk(tconf_acc,tmomenta,
            deltas_Omelyan,4);


    MPI_Waitall(96,send_border_requests,MPI_STATUSES_IGNORE);
    MPI_Waitall(96,recv_border_requests,MPI_STATUSES_IGNORE);

}

void multistep_2MN_gauge_async(su3_soa *tconf_acc,su3_soa *local_staples,tamat_soa *tipdot,thmat_soa *tmomenta)
{
    int md;
    if(verbosity_lv>1) 
        printf("MPI%02d - Performing Async Gauge substeps\n",
                devinfo.myrank);

    multistep_2MN_gauge_async_bloc(tconf_acc,local_staples,
            tipdot,tmomenta,3);

    for(md=1; md<md_parameters.gauge_scale; md++){
        if(verbosity_lv > 2) printf("MPI%02d - Gauge step %d of %d...\n",
                devinfo.myrank,md,md_parameters.gauge_scale);

        multistep_2MN_gauge_async_bloc(tconf_acc,local_staples,
                tipdot,tmomenta,5);

        multistep_2MN_gauge_async_bloc(tconf_acc,local_staples,
                tipdot,tmomenta,6);

    }
    if(verbosity_lv > 2) printf("MPI%02d - Last Gauge step of %d...\n",
            devinfo.myrank,md_parameters.gauge_scale);

    multistep_2MN_gauge_async_bloc(tconf_acc,local_staples,
            tipdot,tmomenta,5);

    // LAST STEP IS NOT ASYNCED
    // Step for the P
    // P' = P - l*dt*dS/dq
    // deltas_Omelyan[3]=-cimag(ieps_acc)*lambda*scale;
    calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);
    mom_sum_mult(tmomenta,tipdot,deltas_Omelyan,3);

}
#endif


void multistep_2MN_gauge(su3_soa *tconf_acc,su3_soa *local_staples,tamat_soa *tipdot,thmat_soa *tmomenta)
{
    if(verbosity_lv>1) 
        printf("MPI%02d - Performing Gauge substeps\n",
                devinfo.myrank);
    int md;
    // Step for the P
    // P' = P - l*dt*dS/dq
    // deltas_Omelyan[3]=-cimag(ieps_acc)*scale*lambda;


    calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot); 
#ifdef DEBUG_MD
    char conffilename[50];
    char momfilename[50];
    char ipdotfilename[50];
    if(!already_printed_debug){
        sprintf(conffilename,"conf_md_0_%s",devinfo.myrankstr);
        sprintf(momfilename,"tmomenta_0_%s",devinfo.myrankstr);
        sprintf(ipdotfilename,"tipdot_0_%s",devinfo.myrankstr);
        dbg_print_su3_soa(tconf_acc,conffilename, 1);
        print_thmat_soa(tmomenta,momfilename);
        print_tamat_soa(tipdot,ipdotfilename);
    }
#endif

    mom_sum_mult(tmomenta,tipdot,deltas_Omelyan,3);
#ifdef DEBUG_MD
    if(!already_printed_debug){
        sprintf(momfilename,"tmomenta_1_%s",devinfo.myrankstr);
        print_thmat_soa(tmomenta,momfilename);
    }
#endif
    //---------------------------
    // Step for the Q
    // Q' = exp[dt/2 *i P] Q
    // deltas_Omelyan[4]=cimag(iepsh_acc)*scale;
    mom_exp_times_conf_soloopenacc(tconf_acc,tmomenta,deltas_Omelyan,4);
#ifdef MULTIDEVICE
    communicate_su3_borders(tconf_acc,GAUGE_HALO);  
#endif
    //-------------------------




    for(md=1; md<md_parameters.gauge_scale; md++){
        if(verbosity_lv > 2) printf("MPI%02d - Gauge step %d of %d...\n",
                devinfo.myrank,md,md_parameters.gauge_scale);
        // Step for the P
        // P' = P - (1-2l)*dt*dS/dq
        // deltas_Omelyan[5]=-cimag(ieps_acc)*(1.0-2.0*lambda)*scale;
        calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);

        mom_sum_mult(tmomenta,tipdot,deltas_Omelyan,5);
        // Step for the Q
        // Q' = exp[dt/2 *i P] Q
        // deltas_Omelyan[4]=cimag(iepsh_acc)*scale;
        mom_exp_times_conf_soloopenacc(tconf_acc,tmomenta,deltas_Omelyan,4);
#ifdef MULTIDEVICE
        communicate_su3_borders(tconf_acc, GAUGE_HALO);  
#endif
        // Step for the P
        // P' = P - 2l*dt*dS/dq
        // deltas_Omelyan[6]=-cimag(ieps_acc)*2.0*lambda*scale;
        calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);

        mom_sum_mult(tmomenta,tipdot,deltas_Omelyan,6);
        //--------------
        // Step for the Q
        // Q' = exp[dt/2 *i P] Q
        // deltas_Omelyan[4]=cimag(iepsh_acc)*scale;
        mom_exp_times_conf_soloopenacc(tconf_acc,tmomenta,deltas_Omelyan,4);
#ifdef MULTIDEVICE
        communicate_su3_borders(tconf_acc, GAUGE_HALO);  
#endif
        //---------------------



    }
    if(verbosity_lv > 2) printf("MPI%02d - Last Gauge step of %d...\n",
            devinfo.myrank,md_parameters.gauge_scale);


    // Step for the P
    // P' = P - (1-2l)*dt*dS/dq
    calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);

    // calc_ipdot_gauge();
    // deltas_Omelyan[5]=-cimag(ieps_acc)*(1.0-2.0*lambda)*scale;
    mom_sum_mult(tmomenta,tipdot,deltas_Omelyan,5);

#ifdef DEBUG_MD
    if(!already_printed_debug){
        sprintf(conffilename,"conf_md_1_%s",devinfo.myrankstr);
        sprintf(momfilename,"tmomenta_3_%s",devinfo.myrankstr);
        sprintf(ipdotfilename,"tipdot_2_%s",devinfo.myrankstr);
        dbg_print_su3_soa(tconf_acc,conffilename, 2);
        print_thmat_soa(tmomenta,momfilename);
        print_tamat_soa(tipdot,ipdotfilename);
    }
    already_printed_debug = 1;
#endif

    // Step for the Q
    // Q' = exp[dt/2 *i P] Q
    // deltas_Omelyan[4]=cimag(iepsh_acc)*scale;
    mom_exp_times_conf_soloopenacc(tconf_acc,tmomenta,deltas_Omelyan,4);
#ifdef MULTIDEVICE
    communicate_su3_borders(tconf_acc, GAUGE_HALO);  
#endif
    // Step for the P
    // P' = P - l*dt*dS/dq
    // deltas_Omelyan[3]=-cimag(ieps_acc)*lambda*scale;
    calc_ipdot_gauge_soloopenacc(tconf_acc,local_staples,tipdot);
    mom_sum_mult(tmomenta,tipdot,deltas_Omelyan,3);


}
void multistep_2MN_SOLOOPENACC( tamat_soa * tipdot_acc,
        su3_soa  * tconf_acc,
#ifdef STOUT_FERMIONS
        su3_soa  * tstout_conf_acc_arr, // huge parking for stouting
#endif
        su3_soa  * tauxbis_conf_acc, 
        su3_soa  * taux_conf_acc,
        ferm_param * tfermions_parameters,// [nflavs]
        int tNDiffFlavs,
        vec3_soa * ferm_in_acc, //[NPS_tot], will be ferm_chi_acc
        vec3_soa * tferm_shiftmulti_acc,// parking variable [max_ps*max_approx_order]
        vec3_soa * tkloc_r, // parking
        vec3_soa * tkloc_h, // parking
        vec3_soa * tkloc_s, // parking
        vec3_soa * tkloc_p, // parking
        vec3_soa * tk_p_shiftferm, // parking, [max_nshift]
        thmat_soa * tmomenta,
        dcomplex_soa * local_sums,
        double res)
{


    int md;

    // Step for the P
    // P' = P - l*dt*dS/dq
    //    deltas_Omelyan[0]=-cimag(ieps_acc)*lambda;
    //  DEOTT_fermion_force_soloopenacc(tconf_acc, 
    fermion_force_soloopenacc(tconf_acc, 
#ifdef STOUT_FERMIONS
            tstout_conf_acc_arr,
#endif
            tauxbis_conf_acc, // parkeggio
            tipdot_acc, tfermions_parameters, tNDiffFlavs, 
            ferm_in_acc, res, taux_conf_acc, tferm_shiftmulti_acc, tkloc_r,
            tkloc_h, tkloc_s, tkloc_p, tk_p_shiftferm);

    if(verbosity_lv > 4) printf("MPI%02d - Calculated first fermion force/n", 
            devinfo.myrank);


    mom_sum_mult(tmomenta,tipdot_acc,deltas_Omelyan,0);

    for(md=1; md<md_parameters.no_md; md++){

        printf("\n\nMPI%02d\t\tRUNNING MD STEP %d OF %d...\n",
                devinfo.myrank, md, md_parameters.no_md);
        // Step for the Q
        // Q' = exp[dt/2 *i P] Q
#ifdef MULTIDEVICE
        if(devinfo.async_comm_gauge)
            multistep_2MN_gauge_async(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta);
        else 
#endif
            multistep_2MN_gauge(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta);
        // Step for the P
        // P' = P - (1-2l)*dt*dS/dq
        // deltas_Omelyan[1]=-cimag(ieps_acc)*(1.0-2.0*lambda);
        //    DEOTT_fermion_force_soloopenacc(tconf_acc, 
        fermion_force_soloopenacc(tconf_acc, 
#ifdef STOUT_FERMIONS
                tstout_conf_acc_arr,
#endif
                tauxbis_conf_acc, // parkeggio
                tipdot_acc, tfermions_parameters, tNDiffFlavs,
                ferm_in_acc, res, taux_conf_acc, tferm_shiftmulti_acc,
                tkloc_r, tkloc_h, tkloc_s, tkloc_p, tk_p_shiftferm);



        mom_sum_mult(tmomenta,tipdot_acc,deltas_Omelyan,1);
        // Step for the Q
        // Q' = exp[dt/2 *i P] Q
#ifdef MULTIDEVICE
        if(devinfo.async_comm_gauge)
            multistep_2MN_gauge_async(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta);
        else
#endif           
            multistep_2MN_gauge(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta);



        // Step for the P
        // P' = P - 2l*dt*dS/dq
        // deltas_Omelyan[2]=-cimag(ieps_acc)*(2.0*lambda);
        //    DEOTT_fermion_force_soloopenacc(tconf_acc,
        fermion_force_soloopenacc(tconf_acc,
#ifdef STOUT_FERMIONS
                tstout_conf_acc_arr, 
#endif
                tauxbis_conf_acc, // parkeggio
                tipdot_acc, tfermions_parameters, tNDiffFlavs, ferm_in_acc, res, taux_conf_acc, tferm_shiftmulti_acc, tkloc_r, tkloc_h, tkloc_s, tkloc_p, tk_p_shiftferm);

        mom_sum_mult(tmomenta,tipdot_acc,deltas_Omelyan,2);
    }  

    printf("\n\nMPI%02d\t\tRUNNING LAST MD STEP OF %d...\n",
            devinfo.myrank,md_parameters.no_md);
    // Step for the Q
    // Q' = exp[dt/2 *i P] Q
#ifdef MULTIDEVICE
    if(devinfo.async_comm_gauge)
        multistep_2MN_gauge_async(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta);
    else
#endif
        multistep_2MN_gauge(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta);


    // Step for the P
    // P' = P - (1-2l)*dt*dS/dq
    // deltas_Omelyan[1]=-cimag(ieps_acc)*(1.0-2.0*lambda);
    //  DEOTT_fermion_force_soloopenacc(tconf_acc,
    fermion_force_soloopenacc(tconf_acc,
#ifdef STOUT_FERMIONS
            tstout_conf_acc_arr,
#endif
            tauxbis_conf_acc, // parkeggio
            tipdot_acc, tfermions_parameters, tNDiffFlavs, ferm_in_acc, res, taux_conf_acc, tferm_shiftmulti_acc, tkloc_r, tkloc_h, tkloc_s, tkloc_p, tk_p_shiftferm);



    mom_sum_mult(tmomenta,ipdot_acc,deltas_Omelyan,1);
    // Step for the Q
    // Q' = exp[dt/2 *i P] Q
#ifdef MULTIDEVICE
    if(devinfo.async_comm_gauge)
        multistep_2MN_gauge_async(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta);
    else 
#endif
        multistep_2MN_gauge(tconf_acc,taux_conf_acc,tipdot_acc,tmomenta);

    // Step for the P
    // P' = P - l*dt*dS/dq
    // deltas_Omelyan[0]=-cimag(ieps_acc)*lambda;
    //DEOTT_fermion_force_soloopenacc(tconf_acc,
    fermion_force_soloopenacc(tconf_acc,
#ifdef STOUT_FERMIONS
            tstout_conf_acc_arr, 
#endif
            tauxbis_conf_acc, // parkeggio
            tipdot_acc, tfermions_parameters, tNDiffFlavs, ferm_in_acc, res, taux_conf_acc, tferm_shiftmulti_acc, tkloc_r, tkloc_h, tkloc_s, tkloc_p, tk_p_shiftferm);
    mom_sum_mult(tmomenta,tipdot_acc,deltas_Omelyan,0);



}// end multistep_2MN_SOLOOPENACC()

#endif

