
// here macros are defined
#define PRINT_DETAILS_INSIDE_UPDATE
#define ALIGN 128

// if using GCC, there are some problems with __restrict.
#ifdef _GNUC_
#define __restrict
#endif

#ifndef _GNUC_
#include "openacc.h"
#endif

#ifdef ONE_FILE_COMPILATION
#include "../Include/all_include.h"
#endif

#include "./field_corr.h"
#include "../Include/memory_wrapper.h"
#include "../Include/debug.h"
#include "../Include/memory_wrapper.h"
#include "../Mpi/multidev.h"
#include "../Include/setting_file_parser.h"
#include "../Include/tell_geom_defines.h"
#include "./alloc_vars.h"
#include "./struct_c_def.h"
#include "./alloc_settings.h"
#include "../Meas/measure_topo.h"
// definitions outside the main.
int conf_id_iter;
int verbosity_lv;


#define ALLOCCHECK(control_int,var)  if(control_int != 0 ) \
    printf("MPI%02d: \tError in  allocation of %s . \n",devinfo.myrank, #var);\
    else if(verbosity_lv > 2) printf("MPI%02d: \tAllocation of %s : OK , %p\n",\
         devinfo.myrank, #var, var );\


int main(int argc, char* argv[]){
  
    su3_soa *  u;
    su3_soa *  field_corr;
    d_complex * trace ;
    int mu, nu, ro, L;
    int allocation_check;
    
    allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&u, ALIGN, 8*sizeof(su3_soa)); 
    ALLOCCHECK(allocation_check, u);
#pragma acc enter data create(u[0:8])
    allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&field_corr, ALIGN, 8*sizeof(su3_soa)); 
    ALLOCCHECK(allocation_check, field_corr );
#pragma acc enter data create(field_corr[0:8])

#ifdef MULTIDEVICE
    pre_init_multidev1D(&devinfo);
    gdbhook();
#endif



    if(0==devinfo.myrank){
        printf("****************************************************\n");
        printf("          PRE INIT - READING SETTING  FILE          \n");
        printf("     check which parameter corresponds to what! \n");
        printf("commit: %s\n", xstr(COMMIT_HASH) );
        printf("****************************************************\n");
    }


    int input_file_read_check = set_global_vars_and_fermions_from_input_file(argv[1]);
		
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
compute_nnp_and_nnm_openacc();
#pragma acc enter data copyin(nnp_openacc)
#pragma acc enter data copyin(nnm_openacc)

for(ro=0; ro<4; ro++){
   for(mu=0 ; mu<3; mu++){
        for(nu=mu+1; nu<4; nu++){          

         calc_field_corr(u, field_corr, traccia, mu, nu, ro);

            for(L=1, L<=nd0/2; L++) {
                    fprintf(fp,"%d;%d;%d;%d;%lf\n", mu, nu, ro, L, traccia[L]);
                                        }
                            }
                    }   
        }

    return 0;
}
