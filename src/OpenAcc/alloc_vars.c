#ifndef ALLOC_VARS_C_
#define ALLOC_VARS_C_

#ifdef __GNUC__
#define _POSIX_C_SOURCE 200809L // not to have warning on posix memalign
#endif

#include <stdio.h>
#include <stdlib.h>

#include "../Include/debug.h"
#include "../Include/fermion_parameters.h"
#include "../Mpi/multidev.h"
#include "./action.h"
#include "./alloc_vars.h"
#include "./struct_c_def.h"
#include "./alloc_settings.h"


#define ALIGN 128
global_su3_soa  * conf_rw; // the gauge configuration, only for read-write

//used in debugging/testing
global_vec3_soa  * ferm_rw; // a global fermion, only for read-write
global_tamat_soa  * tamat_rw; // a global tamat, only for read-write
global_thmat_soa  * thmat_rw; // a global thmat, only for read-write
// a global dcomplex_soa, only for read-write
global_dcomplex_soa *dcomplex_rw;
// a global double_soa, only for read-write
global_double_soa *double_rw; 




su3_soa  * conf_acc; // the gauge configuration.
su3_soa  * conf_acc_bkp; // the old stored conf that will be recovered 
// if the metro test fails.
su3_soa  * aux_conf_acc; // auxiliary 
su3_soa  * auxbis_conf_acc; // auxiliary 
double_soa * u1_back_phases; //Background,staggered,chempot phases
                             // 8 for each flavour

double_soa * mag_obs_re;     // Real part of the 'algebra-prefix'
                             // of magnetization observable 
                             // 8 for each flavour
                                                                       
double_soa * mag_obs_im;     // Imaginary part of the 'algebra-prefix'
                             // of magnetization observable 
                             // 8 for each flavour


thmat_soa * momenta;// GAUGE FIELD EVOLUTION
thmat_soa * momenta_backup;// GAUGE FIELD EVOLUTION - REVERSIBILITY TEST
tamat_soa * ipdot_acc;// GAUGE FIELD EVOLUTION
tamat_soa * ipdot_g_old;// for HMC diagnostics
tamat_soa * ipdot_f_old;// for HMC diagnostics


su3_soa * gconf_as_fermionmatrix; //(only a pointer) conf to use in either cases 
// in fermion related computation (with or without stouting)


// STOUTING 
su3_soa * gstout_conf_acc_arr; // all stouting steps except the zeroth
su3_soa * glocal_staples;
tamat_soa * gipdot;
tamat_soa * aux_ta; // aggiunta per il calcolo della forza stoutata
thmat_soa * aux_th; // aggiunta per il calcolo della forza stoutata

// FERMIONS

vec3_soa * ferm_chi_acc; // questo e' il chi [alloc_info.NPS_tot]
vec3_soa * ferm_phi_acc; // questo e' il phi [alloc_info.NPS_tot]
vec3_soa * ferm_out_acc; // questo e' uno ausiliario [alloc_info.NPS_tot]
vec3_soa * ferm_shiftmulti_acc; // ausiliario per l'invertitore multishift [alloc_info.maxNeededShifts]
vec3_soa * kloc_r;  // vettore ausiliario
vec3_soa * kloc_h;  // vettore ausiliario
vec3_soa * kloc_s;  // vettore ausiliario
vec3_soa * kloc_p;  // vettore ausiliario
vec3_soa * k_p_shiftferm; // ausiliario [alloc_info.maxApproxOrder]

vec3_soa * aux1; // used in fermion force calculation, 
                 // for single precision acceleration


// LOCAL SUMS
dcomplex_soa * local_sums;
double_soa * d_local_sums;


#define ALLOCCHECK(control_int,var)  if(control_int != 0 ) \
    printf("MPI%02d: \tError in  allocation of %s . \n",devinfo.myrank, #var);\
    else if(verbosity_lv > 2) printf("MPI%02d: \tAllocation of %s : OK , %p\n",\
         devinfo.myrank, #var, var );\


void mem_alloc_core(){

    int allocation_check;  

    printf("\n\n[CORE] Allocations..\n");
    
    allocation_check =  posix_memalign((void **)&kloc_r, ALIGN, sizeof(vec3_soa)); 
    ALLOCCHECK(allocation_check, kloc_r) ;
#pragma acc enter data create(kloc_r[0:1])
    allocation_check =  posix_memalign((void **)&kloc_h, ALIGN, sizeof(vec3_soa)); 
    ALLOCCHECK(allocation_check, kloc_h ) ;
#pragma acc enter data create(kloc_h[0:1])
    allocation_check =  posix_memalign((void **)&kloc_s, ALIGN, sizeof(vec3_soa)); 
    ALLOCCHECK(allocation_check, kloc_s ) ;
#pragma acc enter data create(kloc_s[0:1])
    allocation_check =  posix_memalign((void **)&kloc_p, ALIGN, sizeof(vec3_soa)); 
    ALLOCCHECK(allocation_check, kloc_p) ;
#pragma acc enter data create(kloc_p[0:1])

    allocation_check =  posix_memalign((void **)&aux1, ALIGN, sizeof(vec3_soa)); 
    ALLOCCHECK(allocation_check, aux1) ; // used in fermion force calculation, 
                                         // for single precision acceleration
#pragma acc enter data create(aux1[0:1])




    allocation_check =  posix_memalign((void **)&u1_back_phases, ALIGN,
            alloc_info.NDiffFlavs*8*sizeof(double_soa));   
    //  --> alloc_info.NDiffFlavs*4*NSITES phases (as many as links)
    ALLOCCHECK(allocation_check, u1_back_phases);
#pragma acc enter data create(u1_back_phases[0:alloc_info.NDiffFlavs*8])

    alloc_info.conf_acc_size = 8;
#ifdef MULTIDEVICE
    if(devinfo.async_comm_gauge) alloc_info.conf_acc_size *=2 ; 
#endif
    allocation_check =  posix_memalign((void **)&conf_acc, ALIGN, 
            alloc_info.conf_acc_size*sizeof(su3_soa));
    ALLOCCHECK(allocation_check, conf_acc);
#pragma acc enter data create(conf_acc[0:alloc_info.conf_acc_size])


}


void mem_alloc_extended()
{
    printf("\n\n[EXTENDED] Allocating resources for alloc_info.NPS_tot=%d pseudofermions in total, with alloc_info.maxApproxOrder=%d, alloc_info.maxNeededShifts=%d\n", 
            alloc_info.NPS_tot, alloc_info.maxApproxOrder,alloc_info.maxNeededShifts);
    int allocation_check;  

#ifdef MULTIDEVICE 
    if(devinfo.myrank == 0){
#endif 
        // These containers shall not be allocated on the device
        allocation_check =  posix_memalign((void **)&conf_rw, ALIGN,
                8*sizeof(global_su3_soa));
        ALLOCCHECK(allocation_check, conf_rw);

        // used in debugging/testing
        allocation_check =  posix_memalign((void **)&ferm_rw, ALIGN,
                sizeof(global_vec3_soa));
        ALLOCCHECK(allocation_check, ferm_rw);
        allocation_check =  posix_memalign((void **)&tamat_rw, ALIGN,
                8*sizeof(global_tamat_soa));
        ALLOCCHECK(allocation_check, tamat_rw);
        allocation_check =  posix_memalign((void **)&thmat_rw, ALIGN,
                8*sizeof(global_thmat_soa));
        ALLOCCHECK(allocation_check, thmat_rw);
        allocation_check =  posix_memalign((void **)&dcomplex_rw, ALIGN,
                8*sizeof(global_dcomplex_soa));
        ALLOCCHECK(allocation_check, dcomplex_rw);
        allocation_check =  posix_memalign((void **)&double_rw, ALIGN,
                8*sizeof(global_double_soa));
        ALLOCCHECK(allocation_check, double_rw);

#ifdef MULTIDEVICE
    }
#endif



    allocation_check =  posix_memalign((void **)&mag_obs_re, ALIGN,
            alloc_info.NDiffFlavs*8*sizeof(double_soa));   
    //  --> alloc_info.NDiffFlavs*4*NSITES phases (as many as links)
    ALLOCCHECK(allocation_check, mag_obs_re);
#pragma acc enter data create(mag_obs_re[0:alloc_info.NDiffFlavs*8])

    allocation_check =  posix_memalign((void **)&mag_obs_im, ALIGN,
            alloc_info.NDiffFlavs*8*sizeof(double_soa));   
    //  --> alloc_info.NDiffFlavs*4*NSITES phases (as many as links)
    ALLOCCHECK(allocation_check, mag_obs_im);
#pragma acc enter data create(mag_obs_im[0:alloc_info.NDiffFlavs*8])

    allocation_check =  posix_memalign((void **)&aux_conf_acc, ALIGN, 8*sizeof(su3_soa)); 
    ALLOCCHECK(allocation_check, aux_conf_acc );
#pragma acc enter data create(aux_conf_acc[0:8])

    allocation_check =  posix_memalign((void **)&auxbis_conf_acc, ALIGN, 8*sizeof(su3_soa));
    ALLOCCHECK(allocation_check, auxbis_conf_acc ) ;
#pragma acc enter data create(auxbis_conf_acc[0:8])

    // GAUGE EVOLUTION
    allocation_check =  posix_memalign((void **)&momenta, ALIGN, 8*sizeof(thmat_soa));  
    ALLOCCHECK(allocation_check, momenta ) ;
#pragma acc enter data create(momenta[0:8])
    allocation_check =  posix_memalign((void **)&conf_acc_bkp, ALIGN, 8*sizeof(su3_soa));
        ALLOCCHECK(allocation_check, conf_acc_bkp) ;
    alloc_info.revTestAllocations = debug_settings.do_reversibility_test;
    if(alloc_info.revTestAllocations){
        // we allocate it also on the device only if we have to compute
        // differences with the reverse-evolved version
#pragma acc enter data create(conf_acc_bkp[0:8])
        allocation_check =  posix_memalign((void **)&momenta_backup, ALIGN, 8*sizeof(thmat_soa));
        ALLOCCHECK(allocation_check, momenta_backup ) ;
#pragma acc enter data create(momenta_backup[0:8])
    }

    allocation_check =  posix_memalign((void **)&ipdot_acc, ALIGN, 8*sizeof(tamat_soa)); 
    ALLOCCHECK(allocation_check, ipdot_acc) ;
#pragma acc enter data create(ipdot_acc[0:8])

    alloc_info.diagnosticsAllocations = debug_settings.save_diagnostics;

    if(alloc_info.diagnosticsAllocations){
        allocation_check =  posix_memalign((void **)&ipdot_g_old, ALIGN, 8*sizeof(tamat_soa)); 
        ALLOCCHECK(allocation_check, ipdot_g_old) ;
#pragma acc enter data create(ipdot_g_old[0:8])
        allocation_check =  posix_memalign((void **)&ipdot_f_old, ALIGN, 8*sizeof(tamat_soa)); 
        ALLOCCHECK(allocation_check, ipdot_f_old) ;
#pragma acc enter data create(ipdot_f_old[0:8])
    }


    // STOUTING
    if(alloc_info.stoutAllocations){ // not always the same as act_params.stout_steps,
                                     // e.g., in benchmarks
    allocation_check =  posix_memalign((void **)&gstout_conf_acc_arr, ALIGN, act_params.stout_steps*8*sizeof(su3_soa)); 
    ALLOCCHECK(allocation_check,gstout_conf_acc_arr ) ;
#pragma acc enter data create(gstout_conf_acc_arr[0:8*act_params.stout_steps])

    allocation_check =  posix_memalign((void **)&glocal_staples, ALIGN, 8*sizeof(su3_soa)); 
    ALLOCCHECK(allocation_check, glocal_staples) ;
#pragma acc enter data create(glocal_staples[0:8])

    allocation_check =  posix_memalign((void **)&gipdot, ALIGN, 8*sizeof(tamat_soa)); 
    ALLOCCHECK(allocation_check, gipdot) ;
#pragma acc enter data create(gipdot[0:8])

    allocation_check =  posix_memalign((void **)&aux_th, ALIGN, 8*sizeof(thmat_soa)); 
    ALLOCCHECK(allocation_check, aux_th ) ;
#pragma acc enter data create(aux_th[0:8])

    allocation_check =  posix_memalign((void **)&aux_ta, ALIGN, 8*sizeof(tamat_soa)); 
    ALLOCCHECK(allocation_check, aux_ta ) ;
#pragma acc enter data create(aux_ta[0:8])
    }

    // FERMION ALLOCATIONS

    allocation_check =  posix_memalign((void **)&k_p_shiftferm, ALIGN, alloc_info.maxApproxOrder* sizeof(vec3_soa)); 
    ALLOCCHECK(allocation_check, k_p_shiftferm) ;
#pragma acc enter data create(k_p_shiftferm[0:alloc_info.maxApproxOrder])


    allocation_check =  posix_memalign((void **)&ferm_chi_acc  , ALIGN, alloc_info.NPS_tot * sizeof(vec3_soa)); 
    ALLOCCHECK(allocation_check, ferm_chi_acc) ;
#pragma acc enter data create(ferm_chi_acc[0:alloc_info.NPS_tot])

    allocation_check =  posix_memalign((void **)&ferm_phi_acc  , ALIGN, alloc_info.NPS_tot * sizeof(vec3_soa));
    ALLOCCHECK(allocation_check, ferm_phi_acc) ;
#pragma acc enter data create(ferm_phi_acc[0:alloc_info.NPS_tot])

    allocation_check =  posix_memalign((void **)&ferm_out_acc  , ALIGN, alloc_info.NPS_tot * sizeof(vec3_soa));
    ALLOCCHECK(allocation_check, ferm_out_acc) ;
#pragma acc enter data create(ferm_out_acc[0:alloc_info.NPS_tot])

    allocation_check =  posix_memalign((void **)&ferm_shiftmulti_acc, ALIGN, alloc_info.maxNeededShifts*sizeof(vec3_soa)); 
    ALLOCCHECK(allocation_check, ferm_shiftmulti_acc ) ;
#pragma acc enter data create(ferm_shiftmulti_acc[0:alloc_info.maxNeededShifts])


    // REDUCTIONS
    allocation_check =  posix_memalign((void **)&d_local_sums, ALIGN, 2*sizeof(double_soa)); 
    ALLOCCHECK(allocation_check, d_local_sums) ;
#pragma acc enter data create(d_local_sums[0:2])
    
    allocation_check =  posix_memalign((void **)&local_sums, ALIGN, 2*sizeof(dcomplex_soa));
    ALLOCCHECK(allocation_check, local_sums) ;
#pragma acc enter data create(local_sums[0:2])

}

#undef ALLOCCHECK

#define FREECHECK(var) if(verbosity_lv >2) \
    printf("\tFreed %s, %p ...", #var,var);\
    free(var); if(verbosity_lv > 2)  printf(" done.\n");



inline void mem_free_core()
{

    printf("[CORE] Deallocation.\n");
    FREECHECK(kloc_r);                
#pragma acc exit data delete(kloc_r)
    FREECHECK(kloc_h);                
#pragma acc exit data delete(kloc_h)
    FREECHECK(kloc_s);              
#pragma acc exit data delete(kloc_s)
    FREECHECK(kloc_p);                
#pragma acc exit data delete(kloc_p)
    FREECHECK(aux1);                
#pragma acc exit data delete(aux1)                




    FREECHECK(u1_back_phases);        
#pragma acc exit data delete(u1_back_phases)        

    FREECHECK(conf_acc);
#pragma acc exit data delete(conf_acc)


}



inline void mem_free_extended()
{


    printf("[EXTENDED] Deallocation.\n");

#ifdef MULTIDEVICE 
    if(devinfo.myrank == 0){
#endif
        // NOT ON DEVICE
        FREECHECK(conf_rw);
        FREECHECK(ferm_rw);
        FREECHECK(tamat_rw);
        FREECHECK(thmat_rw);
        FREECHECK(dcomplex_rw);
        FREECHECK(double_rw);
#ifdef MULTIDEVICE
    }
#endif

    FREECHECK(mag_obs_re);
#pragma acc exit data delete(mag_obs_re)
    FREECHECK(mag_obs_im);
#pragma acc exit data delete(mag_obs_im)
    FREECHECK(aux_conf_acc);          
#pragma acc exit data delete(aux_conf_acc)          
    FREECHECK(auxbis_conf_acc);       
#pragma acc exit data delete(auxbis_conf_acc)       
    // GAUGE EVOLUTION
    FREECHECK(momenta);               
#pragma acc exit data delete(momenta)               
    FREECHECK(conf_acc_bkp);          
    //alloc_info.revTestAllocations = debug_settings.do_reversibility_test;
    if(alloc_info.revTestAllocations){
        // we allocate it also on the device only if we have to compute
        // differences with the reverse-evolved version
#pragma acc exit data delete(conf_acc_bkp)  
        FREECHECK(momenta_backup);               
#pragma acc exit data delete(momenta_backup)               
    }
    FREECHECK(ipdot_acc);  
#pragma acc exit data delete(ipdot_acc)  

    //alloc_info.diagnosticsAllocations = debug_settings.save_diagnostics;
    if(alloc_info.diagnosticsAllocations){
        FREECHECK(ipdot_g_old);           
#pragma acc exit data delete(ipdot_g_old)           
        FREECHECK(ipdot_f_old);           
#pragma acc exit data delete(ipdot_f_old)           
    }

    // STOUTING
    if(alloc_info.stoutAllocations){
        FREECHECK(gstout_conf_acc_arr);   
#pragma acc exit data delete(gstout_conf_acc_arr)   
        FREECHECK(glocal_staples);        
#pragma acc exit data delete(glocal_staples)        
        FREECHECK(gipdot);              
#pragma acc exit data delete(gipdot)              
        FREECHECK(aux_th);                
#pragma acc exit data delete(aux_th)                
        FREECHECK(aux_ta);                
#pragma acc exit data delete(aux_ta)                
    }
    // FERMION ALLOCATIONS
    FREECHECK(k_p_shiftferm);         
#pragma acc exit data delete(k_p_shiftferm)         
    FREECHECK(ferm_chi_acc);          
#pragma acc exit data delete(ferm_chi_acc)          
    FREECHECK(ferm_phi_acc);          
#pragma acc exit data delete(ferm_phi_acc)          
    FREECHECK(ferm_out_acc);          
#pragma acc exit data delete(ferm_out_acc)          
    FREECHECK(ferm_shiftmulti_acc);   
#pragma acc exit data delete(ferm_shiftmulti_acc)   
    
    // REDUCTIONS
    FREECHECK(local_sums);            
#pragma acc exit data delete(local_sums)            
    FREECHECK(d_local_sums);          
#pragma acc exit data delete(d_local_sums)          

}

#undef FREECHECK



#endif
