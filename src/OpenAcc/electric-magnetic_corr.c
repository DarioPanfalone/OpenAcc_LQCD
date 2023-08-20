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
#define ALLOCCHECK(control_int,var)  if(control_int != 0 )              \
		printf("MPI%02d: \tError in  allocation of %s . \n",devinfo.myrank, #var);\
	else if(verbosity_lv > 2) printf("MPI%02d: \tAllocation of %s : OK , %p\n",\
																	 devinfo.myrank, #var, var );\
	
#define STAMPA_DEBUG_SU3_SOA(var,dir,idx)																\
	printf("%s[%d], idx %d:\n(%le,%le)    (%le,%le)    (%le,%le)\n(%le,%le)    (%le,%le)    (%le,%le)\n(%le,%le)    (%le,%le)    (%le,%le)\n\n", #var, dir, idx, \
		creal(var[dir].r0.c0[idx]),cimag(var[dir].r0.c0[idx]),creal(var[dir].r0.c1[idx]),cimag(var[dir].r0.c1[idx]),creal(var[dir].r0.c2[idx]),cimag(var[dir].r0.c2[idx]), \
		creal(var[dir].r1.c0[idx]),cimag(var[dir].r1.c0[idx]),creal(var[dir].r1.c1[idx]),cimag(var[dir].r1.c1[idx]),creal(var[dir].r1.c2[idx]),cimag(var[dir].r1.c2[idx]), \
		creal(var[dir].r2.c0[idx]),cimag(var[dir].r2.c0[idx]),creal(var[dir].r2.c1[idx]),cimag(var[dir].r2.c1[idx]),creal(var[dir].r2.c2[idx]),cimag(var[dir].r2.c2[idx]));
//
#define MAX 140
int lettura_parole(char file[], char matrice[][MAX]);

int main(int argc, char **argv){
	
	su3_soa *  conf_acc;
	su3_soa *  conf_random;
	su3_soa *  conf_to_use;
	su3_soa *  aux_conf_acc;
	su3_soa *  field_corr;
	su3_soa *  field_corr_aux;
	su3_soa * conf_au;
	su3_soa * aux_staple;
	su3_soa * m_soa;
	single_su3 * closed_corr;
	single_su3 * loc_plaq_aux;
	//single_su3 * m; //elemento random di SU(3)
	double plq, rect;
	global_su3_soa * conf;
	d_complex * corr ;
	dcomplex_soa * local_sum;
		//    int mu, nu, ro, L;
	int allocation_check;
		//		verbosity_lv=6;
//Variabili per lettura file conf
	char confs[MAX][MAX], file[MAX],  nome[MAX];
	FILE *fp;
	FILE *fd;
		
	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&conf, ALIGN,8*sizeof(global_su3_soa));
	ALLOCCHECK(allocation_check, conf);
#pragma acc enter data create(conf[0:8])
	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&conf_acc, ALIGN, 8*sizeof(su3_soa)); 
	ALLOCCHECK(allocation_check, conf_acc);
#pragma acc enter data create(conf_acc[0:8])
	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&conf_random, ALIGN, 8*sizeof(su3_soa));
  ALLOCCHECK(allocation_check, conf_random);
#pragma acc enter data create(conf_random[0:8])
	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&aux_conf_acc, ALIGN, 8*sizeof(su3_soa));
	ALLOCCHECK(allocation_check, aux_conf_acc);
#pragma acc enter data create(aux_conf_acc[0:8])
	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&conf_to_use, ALIGN, 8*sizeof(su3_soa));
	ALLOCCHECK(allocation_check, conf_to_use);
#pragma acc enter data create(conf_to_use[0:8])
	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&conf_au, ALIGN, 8*sizeof(su3_soa));
	ALLOCCHECK(allocation_check, conf_au);
#pragma acc enter data create(conf_au[0:8])
	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&aux_staple, ALIGN, 8*sizeof(su3_soa));
	ALLOCCHECK(allocation_check, aux_staple);
#pragma acc enter data create(aux_staple[0:8])
	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&local_sum, ALIGN, 2*sizeof(dcomplex_soa));
	ALLOCCHECK(allocation_check, local_sum) ;
#pragma acc enter data create(local_sum[0:2])
	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&field_corr, ALIGN, 8*sizeof(su3_soa)); 
	ALLOCCHECK(allocation_check, field_corr );
#pragma acc enter data create(field_corr[0:8])
	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&field_corr_aux, ALIGN, 8*sizeof(su3_soa));
	ALLOCCHECK(allocation_check, field_corr_aux );
#pragma acc enter data create(field_corr_aux[0:8])
	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&m_soa, ALIGN, 8*sizeof(su3_soa));
	ALLOCCHECK(allocation_check, m_soa );
#pragma acc enter data create(m_soa[0:8])
	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&loc_plaq_aux, ALIGN, 8*sizeof(su3_soa));
	ALLOCCHECK(allocation_check, loc_plaq_aux);
#pragma acc enter data create(loc_plaq_aux[0:8])
	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&corr, ALIGN, nd0*sizeof(d_complex));
	ALLOCCHECK(allocation_check, corr);
#pragma acc enter data create(corr[0:nd0])
	allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&closed_corr, ALIGN, sizeof(single_su3));
	ALLOCCHECK(allocation_check, closed_corr);
#pragma acc enter data create(closed_corr[0:1])
	/*  allocation_check =  POSIX_MEMALIGN_WRAPPER((void **)&m, ALIGN, sizeof(single_su3));
	ALLOCCHECK(allocation_check, m);
	#pragma acc enter data create(m[0:1])*/
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
		
	geom_par.xmap = 0 ;
	geom_par.ymap = 1 ;
	geom_par.zmap = 2 ;
	geom_par.tmap = 3 ;

	set_geom_glv(&geom_par);
		
 	//	global_su3_soa * conf_rw = (global_su3_soa * ) malloc(8*sizeof(global_su3_soa));	
	int conf_id;
		//argv[0]="save_conf";
		//printf("Reading conf..\n");
		//int r=read_su3_soa_ASCII(&conf_rw,"save_conf",&conf_id);
	int dim=lettura_parole(argv[1], confs); //leggo i nomi dei file delle conf da input
		//strcpy(nome, confs[0]);
	int maxstep=15, L, confmax=dim;
	fp=fopen(argv[2], "w"); //file dove scrivere le misure
	
	for(int conf_num=0; conf_num<confmax; conf_num++){
				
		int  r = read_su3_soa_ildg_binary(conf,confs[conf_num],&conf_id); 
		
		if(0==r){
			send_lnh_subconf_to_buffer(conf, conf_acc, 0);
#pragma acc update device(conf_acc[0:8])
		} else {
			printf("Some error in reading occured\n");
		}

		for(int coolstep=0; coolstep<maxstep; coolstep++){
			  	
			if(coolstep==0){
				aux_conf_acc=(su3_soa*)conf_acc;
			}
			else if(coolstep==1){
				conf_to_use=(su3_soa*)conf_acc;
				cool_conf(conf_to_use, aux_conf_acc, aux_staple);
			}else{
				conf_to_use=(su3_soa*)aux_conf_acc;
				cool_conf(conf_to_use, aux_conf_acc, aux_staple);		
			}
		  	
			//double  D_paral= 0.0, D_perp = 0.0;
			double D_paral[GL_N0/2], D_perp[GL_N0/2];
			for(int ro=0; ro<4; ro++){ 
				for(int mu=0 ; mu<3; mu++){ 
					for(int nu=mu+1; nu<4; nu++){           
						//								int mu=0, nu=1, ro=2;
						calc_field_corr(aux_conf_acc, field_corr, field_corr_aux, conf_au, local_sum, corr, closed_corr, mu, nu, ro); 
							
						//for(int L=1; L<=nd0/2; L++) { 
						//fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%.18lf\n", coolstep, L, ro, mu, nu, corr[L]/GL_SIZE);
						//ELECTRIC-ELECTRIC
						if(nu==3){
							if(ro==mu || ro==nu){
								for(L=0; L<GL_N0/2; L++){
									
									D_E_paral[L] = D_E_paral[L] + corr[L];
									//fprintf(fp,"%d;%d;%d;%d;%.18lf\n", coolstep, ro, mu, nu, 0.5*creal(trace[L])/GL_SIZE);
								}
							} else {
								for(L=0; L<GL_N0/2; L++){
										
									D_E_perp[L] = D_E_perp[L] + corr[L];
									//fprintf(fd,"%d;%d;%d;%d;%.18lf\n", coolstep, ro, mu, nu, 0.5*creal(trace[L])/GL_SIZE);
								}
							}
							//MAGNETIC-MAGNETIC
							else{
								if(ro==mu || ro==nu){
									for(L=0; L<GL_N0/2; L++){
												
										D_B_paral[L] = D_B_paral[L] + corr[L];
										//fprintf(fp,"%d;%d;%d;%d;%.18lf\n", coolstep, ro, mu, nu, 0.5*creal(trace[L])/GL_SIZE);
									}
								}
							} else {
								for(L=0; L<GL_N0/2; L++){
										
									D_B_perp[L] = D_B_perp[L] + corr[L];
									//fprintf(fd,"%d;%d;%d;%d;%.18lf\n", coolstep, ro, mu, nu, 0.5*creal(trace[L])/GL_SIZE);
								}
							}
						}

					}
				} 
			} 
		}
		//printf("%d\t%.18lf\t%.18lf\n", coolstep, D_paral/((double)4),  D_perp/((double)4));
		//for(L=0; L<GL_N0/2; L++){
		//double D_par, D_per;
		//D_par = D_paral[5], D_per = D_perp[11];
		fprintf(fp,"%d\t%.18lf\t%.18lf\n", coolstep, D_paral[5]/((double)24*GL_SIZE), D_perp[11]/((double)24*GL_SIZE));			
		//}
		// printf("\nfine step %d\n", coolstep);
	}
	printf("\nfine conf %d\n", conf_num);
	}
	fclose(fp);
	
#pragma acc exit data delete(nnp_openacc)
#pragma acc exit data delete(nnm_openacc)					 
		
	return 0;
}

int lettura_parole(char file[], char matrice[][MAX]){
	
	FILE *fp;
	int i=0;
	fp=fopen(file,"r");
	while(!feof(fp)){
		fscanf(fp, "%s", matrice[i]);
		i++;
	}
	fclose(fp);
	return i;

}
