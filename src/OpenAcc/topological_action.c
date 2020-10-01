#ifndef TOPOLOGICAL_ACTION_C
#define TOPOLOGICAL_ACTION_C

#include "./topological_action.h"
#include "../Mpi/multidev.h"
#include "../Meas/gauge_meas.h"
#include "./alloc_vars.h"
#include "./action.h"
#ifndef STOUTING_H
 #include "./stouting.h"
#endif

#ifdef MULTIDEVICE
#include <mpi.h>
#endif


void load_topo(const char *path,const double barr,const double width, double *grid,const int ngrid)
{
	if(verbosity_lv>3)
		printf("Searching file %s\n",path);
	double widthh = width/2;
	FILE *fin=fopen(path,"r");
	if(fin==NULL)
		{
			printf("ERROR: topological file not found");
			exit(1);
		}
	for(int igrid=0; igrid<=ngrid; igrid++)
		{
			double xread;
			int rc=fscanf(fin,"%lf %lf\n",&xread,&grid[igrid]);
			if(rc!=2)
				{
					printf("ERROR: topological file ill-formatted\n");
					exit(1);
				}
			
			int jgrid=(int)floor((xread+barr+widthh)/width);
			if(igrid!=jgrid)
				{
					printf("ERROR: found %d (%lf) when expecting %d - topological file ill-formatted",jgrid,xread,igrid);
					exit(1);
				}
		}
	int check_closure=fclose(fin);
	if(check_closure!=0)
		{
			printf("ERROR: file %s can't be closed",path);
			exit(1);
		}
}



double topodynamical_pot(double Q)
{
	double barr=act_params.barrier, width=act_params.width;
	int ngrid= (int) floor((2*barr+width/2)/width);
	int igrid= (int) floor((Q+barr)/width);

	if(igrid>=0 && igrid<=ngrid)
		{
			double grid[ngrid];
			if(verbosity_lv>3)
				printf("\t\t\tMPI%02d - load_topo(path,barr,width,grid,ngrid)\n",devinfo.myrank);
			
			load_topo(act_params.topo_file_path,barr,width,grid,ngrid);			
			//interpolate
			double x0=igrid*width-barr;
			double m=(grid[igrid+1]-grid[igrid])/width;
			double q=grid[igrid]-m*x0;
		    
			return q+m*Q;
		}
	else{
		printf("\tCarica topologica oltre la barriera\n");
		return 0;
	}
}




double compute_topo_action(su3_soa * const u)
{
	su3_soa *quadri, tstout_conf_acc_arr;
	double_soa *loc_q;
	
	posix_memalign((void **)&quadri,128,8*sizeof(su3_soa));

	posix_memalign((void **)&loc_q,128,2*sizeof(double_soa));
#pragma acc enter data create(quadri[0:8]) create(loc_q[0:2])

#ifdef STOUT_TOPO
	stout_wrapper(u,tstout_conf_acc_arr); //INSERIRE PARAMETRI STOUTING TOPOLOGICO
#else
	tstout_conf_acc_arr=*u;
#endif
	
	if(verbosity_lv>3)
		printf("\t\t\tMPI%02d - compute_topological_charge(u,quadri,loc_q)\n",devinfo.myrank);
		
	double Q = compute_topological_charge(u, quadri, loc_q);
	
	if(verbosity_lv>3)
		printf("\t\t\tMPI%02d - topodynamical_pot(Q)\n",devinfo.myrank);

	double topo_action = topodynamical_pot(Q);
	
	free(quadri);
#pragma acc exit data delete(quadri)

	free(loc_q);
#pragma acc exit data delete(loc_q)

	return topo_action;
}

#endif
