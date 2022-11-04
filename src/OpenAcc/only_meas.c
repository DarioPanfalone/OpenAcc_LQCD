// here macros are defined
#define PRINT_DETAILS_INSIDE_UPDATE
#define ALIGN 128


#include <mpi.h>
#include "../Include/stringify.h"

#include <errno.h>
#include "./HPT_utilities.h"
#include <time.h>

#ifdef ONE_FILE_COMPILATION
#include "../Include/all_include.h"
#endif
#include "../DbgTools/dbgtools.h" // DEBUG
#include "../DbgTools/debugger_hook.h"
#include "../Include/debug.h"
#include "../Include/fermion_parameters.h"
#include "../Include/montecarlo_parameters.h"
#include "../Include/inverter_tricks.h"
#include "../Include/memory_wrapper.h"
#include "../Include/setting_file_parser.h"
#include "../Include/tell_geom_defines.h"
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
#include "./rectangles.h"
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

// definitions outside the main.
int conf_id_iter;
int verbosity_lv;
//
#define STAMPA_DEBUG_SU3_SOA(var,dir,idx)\
	printf("%s[%d], idx %d:\n(%le,%le)    (%le,%le)    (%le,%le)\n(%le,%le)    (%le,%le)    (%le,%le)\n(%le,%le)    (%le,%le)    (%le,%le)\n\n", #var, dir, idx, \
		creal(var[dir].r0.c0[idx]),cimag(var[dir].r0.c0[idx]),creal(var[dir].r0.c1[idx]),cimag(var[dir].r0.c1[idx]),creal(var[dir].r0.c2[idx]),cimag(var[dir].r0.c2[idx]), \
		creal(var[dir].r1.c0[idx]),cimag(var[dir].r1.c0[idx]),creal(var[dir].r1.c1[idx]),cimag(var[dir].r1.c1[idx]),creal(var[dir].r1.c2[idx]),cimag(var[dir].r1.c2[idx]), \
		creal(var[dir].r2.c0[idx]),cimag(var[dir].r2.c0[idx]),creal(var[dir].r2.c1[idx]),cimag(var[dir].r2.c1[idx]),creal(var[dir].r2.c2[idx]),cimag(var[dir].r2.c2[idx]));
//
int main(int argc, char **argv){
  
    su3_soa *  u;
    su3_soa *  field_corr;
		su3_soa * conf_au;
		double plq, rect;
		global_su3_soa * conf;
    dcomplex_soa * trace ;
		dcomplex_soa * local_sum;
		//    int mu, nu, ro, L;
    int allocation_check;
		//		verbosity_lv=6;		
#define ALLOCCHECK(control_int,var)  if(control_int != 0 )							\
    printf("MPI%02d: \tError in  allocation of %s . \n",devinfo.myrank, #var);\
    else if(verbosity_lv > 2) printf("MPI%02d: \tAllocation of %s : OK , %p\n",\
         devinfo.myrank, #var, var );\
	
	 	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&conf, ALIGN,8*sizeof(global_su3_soa));
		ALLOCCHECK(allocation_check, conf);
#pragma acc enter data create(conf[0:8])
    allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&u, ALIGN, 8*sizeof(su3_soa)); 
		ALLOCCHECK(allocation_check, u);
#pragma acc enter data create(u[0:8])
		allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&conf_au, ALIGN, 8*sizeof(su3_soa));
		ALLOCCHECK(allocation_check, conf_au);
#pragma acc enter data create(conf_au[0:8])
		allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&local_sum, ALIGN, 2*sizeof(dcomplex_soa));
		ALLOCCHECK(allocation_check, local_sum) ;
#pragma acc enter data create(local_sum[0:2])
		allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&field_corr, ALIGN, 8*sizeof(su3_soa)); 
    ALLOCCHECK(allocation_check, field_corr );
#pragma acc enter data create(field_corr[0:8])
		
		//		allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&trace, ALIGN, nd0*sizeof(d_complex));
		//		ALLOCCHECK(allocation_check, trace );
		//pragma acc enter data create(trace[0:nd0])
		
#ifdef MULTIDEVICE
    pre_init_multidev1D(&devinfo);
    gdbhook();
#endif
		
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
		// select_init_acc_device(my_device_type, devinfo.myrank%devinfo.proc_per_node);
		#else
		select_init_acc_device(my_device_type, devinfo.single_dev_choice);
		#endif
		printf("Device Selected : OK \n");
    #endif


		if(0==devinfo.myrank) print_geom_defines(); 
		compute_nnp_and_nnm_openacc();
#pragma acc enter data copyin(nnp_openacc) 
#pragma acc enter data copyin(nnm_openacc) 

		printf("HARDCODED LATTICE DIMENSIONS:\n");
		printf("GL_N0: %d\n", GL_N0) ;
		printf("GL_N1: %d\n", GL_N1) ;
		printf("GL_N2: %d\n", GL_N2) ;
		printf("GL_N3: %d\n", GL_N3) ;
		
		geom_par.gnx = GL_N0 ;
		geom_par.gny = GL_N1 ;
		geom_par.gnz = GL_N2 ;
		geom_par.gnt = GL_N3 ;

	 	geom_par.xmap = 3 ;
	 	geom_par.ymap = 2 ;
		geom_par.zmap = 1 ;
		geom_par.tmap = 0 ;

		set_geom_glv(&geom_par);
		
 	//	global_su3_soa * conf_rw = (global_su3_soa * ) malloc(8*sizeof(global_su3_soa));	
 		int conf_id;
		//argv[0]="save_conf";
		printf("Reading conf..\n");
		//int r=read_su3_soa_ASCII(&conf_rw,"save_conf",&conf_id);
 		int  r = read_su3_soa_ildg_binary(conf,"save_conf",&conf_id); 
		
	 	if(0==r){
		send_lnh_subconf_to_buffer(conf, u, 0);
#pragma acc update device(u[0:8])
		    } else {
			printf("Some error in reading occured\n");
		}
		//printf("fatto\n");
		//prova stampa link
		int d0=1, d1=11, d2=1, d3=D3_HALO;
		int mu=1;
		int idxh = snum_acc(d0, d1, d2, d3);
		int parity = (d0 + d1 + d2 + d3) % 2;
		int dir_link = 2*mu + parity;
		
	  
#pragma acc update self(u[0:8])
		STAMPA_DEBUG_SU3_SOA(u,dir_link,idxh);
		printf("%d\n",nnp_openacc[idxh][mu][parity]);

		plq = calc_plaquette_soloopenacc(u, conf_au, local_sum);
		//		rect = calc_rettangolo_soloopenacc(u, conf_au, local_sum);
		//	printf("rettangolo     : %.18lf\n" ,rect);
		printf("Plaquette     : %.18lf\n" ,plq/GL_SIZE/3.0/6.0);

	#pragma acc exit data delete(nnp_openacc)
  #pragma acc exit data delete(nnm_openacc)
		
/* for(ro=0; ro<4; ro++){ */
/*    for(mu=0 ; mu<3; mu++){ */
/*         for(nu=mu+1; nu<4; nu++){           */

/*          calc_field_corr(u, field_corr, traccia, mu, nu, ro); */

/*             for(L=1, L<=nd0/2; L++) { */
/*                     fprintf(fp,"%d;%d;%d;%d;%lf\n", mu, nu, ro, L, traccia[L]); */
/*                                         } */
/*                             } */
/*                     }    */
/*         } */

    return 0;
}
