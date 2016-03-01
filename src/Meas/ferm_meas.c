// use z2 noise instead of gaussian noise (see hep-lat/9308015)
// use the global defined fermions loc_chi, loc_phi, rnd_o, rnd_e, chi_o and loc_h

#ifndef FERM_MEAS_C 
#define FERM_MEAS_C

#define ALIGN 128

#ifdef __GNUC__
#define _POSIX_C_SOURCE 200809L   // not to have warning on posix memalign
#endif




#include <stdio.h>
#include <stdlib.h>
#include "../OpenAcc/struct_c_def.h"
#include "../OpenAcc/inverter_full.h"
#include "../OpenAcc/alloc_vars.h"
#include "../OpenAcc/random_assignement.h"
#include "./ferm_meas.h"
#include "./baryon_number_utilities.h"
#include "../Include/fermion_parameters.h"
#include "../OpenAcc/fermion_matrix.h"
#include "../OpenAcc/fermionic_utilities.h"
#include "../DbgTools/debug_macros_glvarcheck.h"
#include "../Include/common_defines.h"
#ifdef STOUT_FERMIONS
#include "../OpenAcc/stouting.h"
#endif
#include "../OpenAcc/action.h"

ferm_meas_params  fm_par;

// vedi tesi LS F.Negro per ragguagli (Appendici)
void eo_inversion(su3_soa *tconf_acc,
        ferm_param * tfermions_parameters,
        double res,
        vec3_soa * in_e,     // z2 noise
        vec3_soa * in_o,     // z2 noise
        vec3_soa * out_e,
        vec3_soa * out_o,
        vec3_soa * phi_e,    // parking variable
        vec3_soa * phi_o,    // parking variable
        vec3_soa * trialSolution,       // initial vector for the inversion 
        vec3_soa * tloc_r,    // parking variable for the inverter      
        vec3_soa * tloc_h,    // parking variable for the inverter
        vec3_soa * tloc_s,    // parking variable for the inverter
        vec3_soa * tloc_p){   // parking variable for the inverter




    acc_Deo(tconf_acc, phi_e, in_o,tfermions_parameters->phases);
    combine_in1_x_fact1_minus_in2_back_into_in2(in_e, tfermions_parameters->ferm_mass , phi_e);
    ker_invert_openacc(tconf_acc,tfermions_parameters,
            out_e,phi_e,res,trialSolution,
            tloc_r,tloc_h,tloc_s,tloc_p);
    acc_Doe(tconf_acc, phi_o, out_e,tfermions_parameters->phases);
    combine_in1_minus_in2_allxfact(in_o,phi_o,(double)1/tfermions_parameters->ferm_mass,out_o);


}// end eo_inversion


d_complex eo_scal_prod_global(vec3_soa * rnd_e,
        vec3_soa * rnd_o,
        vec3_soa * chi_e,
        vec3_soa * chi_o){

    return scal_prod_global(rnd_o,chi_o) + scal_prod_global(rnd_e,chi_e);

}




void set_fermion_file_header(ferm_meas_params * fmpar, ferm_param * tferm_par){

    strcpy(fmpar->fermionic_outfile_header,"#conf\ticopy\t");
    for(int iflv=0;iflv<NDiffFlavs;iflv++){
        char strtocat[120];
        sprintf(strtocat, "Reff_%-19sImff_%-19s",tferm_par[iflv].name,tferm_par[iflv].name);
        strcat(fmpar->fermionic_outfile_header,strtocat);
        sprintf(strtocat, "ReN_%-20sImN_%-20s",tferm_par[iflv].name,tferm_par[iflv].name);
        strcat(fmpar->fermionic_outfile_header,strtocat);
         if (fmpar->DoubleInvNVectors>0){

        sprintf(strtocat, "ReTrM^-2_%-14sImTrM^-2_%-14s",tferm_par[iflv].name,tferm_par[iflv].name);
        strcat(fmpar->fermionic_outfile_header,strtocat);
         }

    }
    strcat(fmpar->fermionic_outfile_header,"\n");

}


void fermion_measures( su3_soa * tconf_acc,
        ferm_param * tfermions_parameters,
        ferm_meas_params * tfm_par,
        double res,
        int conf_id_iter  ){
    vec3_soa * rnd_e,* rnd_o;
    vec3_soa * chi_e,* chi_o; //results of eo_inversion
    vec3_soa *chi2_e,*chi2_o; //results of second eo_inversion, 
    vec3_soa * phi_e,* phi_o; // parking variables for eo_inverter
    // also results of eo inversion multiplied by dM/dmu
    vec3_soa * trial_sol;
    su3_soa * conf_to_use;

#ifdef STOUT_FERMIONS
    SETREQUESTED(gstout_conf_acc_arr);
    stout_wrapper(tconf_acc ,gstout_conf_acc_arr);
    conf_to_use = &gstout_conf_acc_arr[8*(act_params.stout_steps-1)];
#else
    conf_to_use = tconf_acc;
#endif

    int allocation_check;
#define ALLOCCHECK(control_int,var)  if(control_int != 0 ) \
    printf("\tError in  allocation of %s . \n", #var);\
    else if(verbosity_lv > 2) printf("\tAllocation of %s : OK , %p\n", #var, var );\

    allocation_check =  posix_memalign((void **)&rnd_e, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,rnd_e);     
    allocation_check =  posix_memalign((void **)&rnd_o, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,rnd_o);     
    allocation_check =  posix_memalign((void **)&phi_e, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,phi_e);     
    allocation_check =  posix_memalign((void **)&phi_o, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,phi_o);     
    allocation_check =  posix_memalign((void **)&chi_e, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,chi_e);     
    allocation_check =  posix_memalign((void **)&chi_o, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,chi_o);     
    allocation_check =  posix_memalign((void **)&chi2_e, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,chi2_e);     
    allocation_check =  posix_memalign((void **)&chi2_o, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,chi2_o);     

    allocation_check =  posix_memalign((void **)&trial_sol, ALIGN, sizeof(vec3_soa));
    ALLOCCHECK(allocation_check,trial_sol);

#undef ALLOCCHECK

    // FILE checks
    // find the actual file size
    FILE *foutfile = fopen(tfm_par->fermionic_outfilename,"ab");
    int fsize;
    if(foutfile){
        fseek(foutfile, 0L, SEEK_END);
        fsize = ftell(foutfile);
        fseek(foutfile, 0L, SEEK_SET);
        fsize -= ftell(foutfile);

    }else {
        printf("File %s can't be opened for writing. Exiting.\n", fm_par.fermionic_outfilename);
        exit(1);
    }
    fclose(foutfile);// found file size

    // cycle on copies
    for(int icopy = 0; icopy < tfm_par->SingleInvNVectors ; icopy++) {
        FILE *foutfile = fopen(tfm_par->fermionic_outfilename,"at");
        if(fsize == 0){
            set_fermion_file_header(tfm_par, tfermions_parameters);
            fprintf(foutfile,"%s",fm_par.fermionic_outfile_header);
            fsize++;
        }

        if(!foutfile) {
            printf("File %s can't be opened for writing. Exiting.\n", fm_par.fermionic_outfilename);
            exit(1);
        }


        fprintf(foutfile,"%d\t%d\t",conf_id_iter,icopy);

        // cycle on flavours
        for(int iflv=0; iflv < NDiffFlavs ; iflv++){

            if(verbosity_lv > 1){
                printf("Performing %d of %d chiral measures for quark %s",
                    icopy+1,tfm_par->SingleInvNVectors, tfermions_parameters[iflv].name);
    
                if(icopy < tfm_par->DoubleInvNVectors )
                printf(" ( %d of %d double inversions)",
                    icopy+1,tfm_par->DoubleInvNVectors);
            
                printf(".\n");

            
            }
           

            generate_vec3_soa_z2noise(rnd_e);
            generate_vec3_soa_z2noise(rnd_o);
            generate_vec3_soa_gauss(trial_sol);
            d_complex chircond_size = 0.0 + 0.0*I;
            d_complex barnum_size = 0.0 + 0.0*I; // https://en.wikipedia.org/wiki/P._T._Barnum
            d_complex trMinvSq_size = 0.0 + 0.0*I;
            double factor = tfermions_parameters[iflv].degeneracy*0.25/size;
            double factor2 = factor*factor*size;

#pragma acc data create(phi_e[0:1]) create(phi_o[0:1])\
            create(chi_e[0:1]) create(chi_o[0:1])   \
            create(chi2_e[0:1]) create(chi2_o[0:1])   \
            copyin(rnd_e[0:1]) copyin(rnd_o[0:1]) copyin(trial_sol[0:1])
            {
                // i fermioni ausiliari kloc_* sono quelli GLOBALI !!!
                eo_inversion(conf_to_use,&tfermions_parameters[iflv],res,rnd_e,rnd_o,chi_e,chi_o,phi_e,phi_o,trial_sol,kloc_r,kloc_h,kloc_s,kloc_p);
                chircond_size = scal_prod_global(rnd_e,chi_e)+
                    scal_prod_global(rnd_o,chi_o);

                dM_dmu_eo[geom_par.tmap](conf_to_use,phi_e,chi_o,tfermions_parameters[iflv].phases);
                dM_dmu_oe[geom_par.tmap](conf_to_use,phi_o,chi_e,tfermions_parameters[iflv].phases);
                barnum_size = scal_prod_global(rnd_e,phi_e)+
                    scal_prod_global(rnd_o,phi_o);

                fprintf(foutfile,"%.16lf\t%.16lf\t",
                        creal(chircond_size*factor),cimag(chircond_size*factor));
                fprintf(foutfile,"%.16lf\t%.16lf\t",
                        creal(barnum_size*factor),cimag(barnum_size*factor));

                if(icopy < tfm_par->DoubleInvNVectors ){
                eo_inversion(conf_to_use,&tfermions_parameters[iflv],res,chi_e,chi_o,chi2_e,chi2_o,phi_e,phi_o,trial_sol,kloc_r,kloc_h,kloc_s,kloc_p);

                trMinvSq_size = scal_prod_global(chi2_e,rnd_e)+
                  scal_prod_global(chi2_o,rnd_o); 
                fprintf(foutfile,"%.16lf\t%.16lf\t",
                        creal(trMinvSq_size*factor2),cimag(trMinvSq_size*factor2));

                } else if (tfm_par->DoubleInvNVectors>0)                  
                    fprintf(foutfile,"%-24s%-24s","none","none" );




            }
        }// end of cycle over flavours

        fprintf(foutfile,"\n");
        fclose(foutfile);
    }// end of cycle over copies

    free(rnd_e);
    free(rnd_o);
    free(phi_e);
    free(phi_o);
    free(chi_e);
    free(chi_o);
    free(chi2_e);
    free(chi2_o);
    free(trial_sol);


}


#endif
