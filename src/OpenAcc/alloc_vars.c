#ifndef ALLOC_VARS_C_
#define ALLOC_VARS_C_

#ifdef __GNUC__
#define _POSIX_C_SOURCE 200809L // not to have warning on posix memalign
#endif

#include <stdio.h>
#include <stdlib.h>

#include "../DbgTools/debug_macros_glvarcheck.h"
#include "../Include/debug.h"
#include "../Include/fermion_parameters.h"
#include "../Mpi/multidev.h"
#include "./action.h"
#include "./alloc_vars.h"
#include "./struct_c_def.h"


#define ALIGN 128
global_su3_soa  * conf_rw; // the gauge configuration, only for read-write
global_vec3_soa  * ferm_rw; // a global fermion, only for read-write
int conf_acc_size;
su3_soa  * conf_acc; // the gauge configuration.
su3_soa  * conf_acc_bkp; // the old stored conf that will be recovered 
// if the metro test fails.
su3_soa  * aux_conf_acc; // auxiliary 
su3_soa  * auxbis_conf_acc; // auxiliary 
double_soa * u1_back_phases; //Background,staggered,chempot phases
// 8 for each flavour

thmat_soa * momenta;// GAUGE FIELD EVOLUTION
int momenta_backupped; 
thmat_soa * momenta_backup;// GAUGE FIELD EVOLUTION - REVERSIBILITY TEST
tamat_soa * ipdot_acc;// GAUGE FIELD EVOLUTION
tamat_soa * ipdot_g_old;// for HMC diagnostics
tamat_soa * ipdot_f_old;// for HMC diagnostics

su3_soa * gconf_as_fermionmatrix; //(only a pointer) conf to use in either cases 
// in fermion related computation (with or without stouting)


// STOUTING 
#ifdef STOUT_FERMIONS
su3_soa * gstout_conf_acc; // max stouted conf, just pointer
su3_soa * gstout_conf_acc_arr; // all stouting steps except the zeroth
su3_soa * glocal_staples;
tamat_soa * gipdot;
tamat_soa * aux_ta; // aggiunta per il calcolo della forza stoutata
thmat_soa * aux_th; // aggiunta per il calcolo della forza stoutata
#endif

// FERMIONS

vec3_soa * ferm_chi_acc; // questo e' il chi [NPS_tot]
vec3_soa * ferm_phi_acc; // questo e' il phi [NPS_tot]
vec3_soa * ferm_out_acc; // questo e' uno ausiliario [NPS_tot]
vec3_soa * ferm_shiftmulti_acc; // ausiliario per l'invertitore multishift [max_ps*MAX_APPROX_ORDER]
vec3_soa * kloc_r;  // vettore ausiliario
vec3_soa * kloc_h;  // vettore ausiliario
vec3_soa * kloc_s;  // vettore ausiliario
vec3_soa * kloc_p;  // vettore ausiliario
vec3_soa * k_p_shiftferm; // ausiliario [max_nshift=MAX_APPROX_ORDER]


// LOCAL SUMS
dcomplex_soa * local_sums;
double_soa * d_local_sums;







void mem_alloc()
{
    printf("\n\nAllocating resources for NPS_tot=%d pseudofermions in total, with MAX_APPROX_ORDER=%d\n", NPS_tot, MAX_APPROX_ORDER);
    int allocation_check;  


#define ALLOCCHECK(control_int,var)  if(control_int != 0 ) \
    printf("MPI%02d: \tError in  allocation of %s . \n",devinfo.myrank, #var);\
    else if(verbosity_lv > 2) printf("MPI%02d: \tAllocation of %s : OK , %p\n",\
         devinfo.myrank, #var, var );\

    allocation_check =  posix_memalign((void **)&u1_back_phases, ALIGN,
            NDiffFlavs*8*sizeof(double_soa));   
    //  --> NDiffFlavs*4*NSITES phases (as many as links)
    SETFREE(u1_back_phases);
    ALLOCCHECK(allocation_check, u1_back_phases);



#ifdef MULTIDEVICE 
    if(devinfo.myrank == 0){
#endif
        allocation_check =  posix_memalign((void **)&conf_rw, ALIGN,8*sizeof(global_su3_soa));
        ALLOCCHECK(allocation_check, conf_rw);
        allocation_check =  posix_memalign((void **)&ferm_rw, ALIGN,sizeof(global_vec3_soa));
        ALLOCCHECK(allocation_check, ferm_rw);
#ifdef MULTIDEVICE
    }
#endif

    conf_acc_size = 8;
#ifdef MULTIDEVICE
    if(devinfo.async_comm_gauge) conf_acc_size *=2 ; 
#endif
    allocation_check =  posix_memalign((void **)&conf_acc, ALIGN, 
            conf_acc_size*sizeof(su3_soa));
    ALLOCCHECK(allocation_check, conf_acc);





    //the double bracket in the setfree macro MUST be there(because of operators precedence)
    allocation_check =  posix_memalign((void **)&aux_conf_acc, ALIGN, 8*sizeof(su3_soa)); for(int mu=0;mu<8;mu++) SETFREE((&aux_conf_acc[mu])); 
    allocation_check =  posix_memalign((void **)&auxbis_conf_acc, ALIGN, 8*sizeof(su3_soa)); for(int mu=0;mu<8;mu++) SETFREE((&auxbis_conf_acc[mu]));
    ALLOCCHECK(allocation_check, auxbis_conf_acc ) ;
    allocation_check =  posix_memalign((void **)&conf_acc_bkp, ALIGN, 8*sizeof(su3_soa)); for(int mu=0;mu<8;mu++) SETFREE((&conf_acc_bkp[mu])); 
    ALLOCCHECK(allocation_check, conf_acc_bkp) ;


    // GAUGE EVOLUTION
    allocation_check =  posix_memalign((void **)&momenta, ALIGN, 8*sizeof(thmat_soa));  for(int mu=0;mu<8;mu++) SETFREE((&momenta[mu]));  //  -->  4*NSITES
    ALLOCCHECK(allocation_check, momenta ) ;
    if(debug_settings.do_reversibility_test){
        momenta_backupped = 1;

        allocation_check =  posix_memalign((void **)&momenta_backup, ALIGN, 8*sizeof(thmat_soa));  for(int mu=0;mu<8;mu++) SETFREE((&momenta_backup[mu]));  //  -->  4*NSITES
        ALLOCCHECK(allocation_check, momenta_backup ) ;
    }
    else momenta_backupped = 0;


    allocation_check =  posix_memalign((void **)&ipdot_acc, ALIGN, 8*sizeof(tamat_soa)); for(int mu=0;mu<8;mu++) SETFREE((&ipdot_acc[mu])); 
    ALLOCCHECK(allocation_check, ipdot_acc) ;
    allocation_check =  posix_memalign((void **)&ipdot_g_old, ALIGN, 8*sizeof(tamat_soa)); for(int mu=0;mu<8;mu++) SETFREE((&ipdot_g_old[mu])); 
    ALLOCCHECK(allocation_check, ipdot_g_old) ;
    allocation_check =  posix_memalign((void **)&ipdot_f_old, ALIGN, 8*sizeof(tamat_soa)); for(int mu=0;mu<8;mu++) SETFREE((&ipdot_f_old[mu])); 
    ALLOCCHECK(allocation_check, ipdot_f_old) ;



#ifdef STOUT_FERMIONS
    // STOUTING
    allocation_check =  posix_memalign((void **)&gstout_conf_acc_arr, ALIGN, act_params.stout_steps*8*sizeof(su3_soa)); for(int mu=0;mu<8*act_params.stout_steps;mu++) SETFREE((&gstout_conf_acc_arr[mu]));
    gstout_conf_acc = &gstout_conf_acc_arr[8*(act_params.stout_steps-1)];
    ALLOCCHECK(allocation_check,gstout_conf_acc_arr ) ;
    allocation_check =  posix_memalign((void **)&glocal_staples, ALIGN, 8*sizeof(su3_soa)); for(int mu=0;mu<8;mu++) SETFREE((&glocal_staples[mu]));
    ALLOCCHECK(allocation_check, glocal_staples) ;
    allocation_check =  posix_memalign((void **)&gipdot, ALIGN, 8*sizeof(tamat_soa)); for(int mu=0;mu<8;mu++) SETFREE((&gipdot[mu])); 
    ALLOCCHECK(allocation_check, gipdot) ;

    allocation_check =  posix_memalign((void **)&aux_th, ALIGN, 8*sizeof(thmat_soa)); for(int mu=0;mu<8;mu++) SETFREE((&aux_th[mu]));   //  -->  4*NSITES
    ALLOCCHECK(allocation_check, aux_th ) ;
    allocation_check =  posix_memalign((void **)&aux_ta, ALIGN, 8*sizeof(tamat_soa)); for(int mu=0;mu<8;mu++)  SETFREE((&aux_ta[mu])); //  -->  4*NSITES
    ALLOCCHECK(allocation_check, aux_ta ) ;



#endif



    // FERMION ALLOCATIONS
    allocation_check =  posix_memalign((void **)&kloc_r, ALIGN, sizeof(vec3_soa)); SETFREE(kloc_r);
    ALLOCCHECK(allocation_check, kloc_r) ;
    allocation_check =  posix_memalign((void **)&kloc_h, ALIGN, sizeof(vec3_soa)); SETFREE(kloc_h);
    ALLOCCHECK(allocation_check, kloc_h ) ;
    allocation_check =  posix_memalign((void **)&kloc_s, ALIGN, sizeof(vec3_soa)); SETFREE(kloc_s);
    ALLOCCHECK(allocation_check, kloc_s ) ;
    allocation_check =  posix_memalign((void **)&kloc_p, ALIGN, sizeof(vec3_soa)); SETFREE(kloc_p);
    ALLOCCHECK(allocation_check, kloc_p) ;
    allocation_check =  posix_memalign((void **)&k_p_shiftferm, ALIGN, MAX_APPROX_ORDER* sizeof(vec3_soa)); for(int mu=0;mu<MAX_APPROX_ORDER;mu++) SETFREE((&k_p_shiftferm[mu]));
    ALLOCCHECK(allocation_check, k_p_shiftferm) ;


    allocation_check =  posix_memalign((void **)&ferm_chi_acc  , ALIGN, NPS_tot * sizeof(vec3_soa)); for(int mu=0;mu<NPS_tot;mu++) SETFREE((&ferm_chi_acc[mu]));
    ALLOCCHECK(allocation_check, ferm_chi_acc) ;
    allocation_check =  posix_memalign((void **)&ferm_phi_acc  , ALIGN, NPS_tot * sizeof(vec3_soa)); for(int mu=0;mu<NPS_tot;mu++) SETFREE((&ferm_phi_acc[mu]));
    ALLOCCHECK(allocation_check, ferm_phi_acc) ;
    allocation_check =  posix_memalign((void **)&ferm_out_acc  , ALIGN, NPS_tot * sizeof(vec3_soa)); for(int mu=0;mu<NPS_tot;mu++) SETFREE((&ferm_out_acc[mu]));
    ALLOCCHECK(allocation_check, ferm_out_acc) ;
    allocation_check =  posix_memalign((void **)&ferm_shiftmulti_acc, ALIGN, max_ps*MAX_APPROX_ORDER*sizeof(vec3_soa)); for(int mu=0;mu<max_ps*MAX_APPROX_ORDER;mu++) SETFREE((&ferm_shiftmulti_acc[mu]));
    ALLOCCHECK(allocation_check, ferm_shiftmulti_acc ) ;



    allocation_check =  posix_memalign((void **)&d_local_sums, ALIGN, 2*sizeof(double_soa)); SETFREE(d_local_sums);
    ALLOCCHECK(allocation_check, d_local_sums) ;
    allocation_check =  posix_memalign((void **)&local_sums, ALIGN, 2*sizeof(dcomplex_soa)); for(int mu=0;mu<2;mu++) SETFREE((&local_sums[mu]));
    ALLOCCHECK(allocation_check, local_sums) ;



#undef ALLOCCHECK

}

inline void mem_free()
{


#define FREECHECK(var) if(verbosity_lv >2) \
    printf("\tFreed %s, %p ...", #var,var);\
    free(var); if(verbosity_lv > 2)  printf(" done.\n");

    printf("Deallocation.\n");

#ifdef MULTIDEVICE 
    if(devinfo.myrank == 0){
#endif
        FREECHECK(conf_rw);
        //  FREECHECK(ferm_rw);
#ifdef MULTIDEVICE
    }
#endif

    //  FREECHECK(conf_acc);
    //  FREECHECK(u1_back_phases);        
    FREECHECK(momenta);               
    if(debug_settings.do_reversibility_test){
        FREECHECK(momenta_backup);               
    }
    FREECHECK(aux_conf_acc);          
    FREECHECK(auxbis_conf_acc);       

#ifdef STOUT_FERMIONS               
    FREECHECK(gstout_conf_acc_arr);   
    FREECHECK(glocal_staples);        
    FREECHECK(gipdot);              
    FREECHECK(aux_ta);                
    FREECHECK(aux_th);                
#endif                              


    FREECHECK(conf_acc_bkp);          
    FREECHECK(ipdot_acc);           
    FREECHECK(ipdot_g_old);           
    FREECHECK(ipdot_f_old);           

    //  FREECHECK(ferm_chi_acc);          
    //  FREECHECK(ferm_phi_acc);          
    //  FREECHECK(ferm_out_acc);          
    //    
    FREECHECK(ferm_shiftmulti_acc);   
    //                                
    FREECHECK(kloc_r);                
    FREECHECK(kloc_s);                
    FREECHECK(kloc_h);                
    FREECHECK(kloc_p);                
    FREECHECK(k_p_shiftferm);         
    //                                
    FREECHECK(local_sums);            
    FREECHECK(d_local_sums);          
}

#endif
