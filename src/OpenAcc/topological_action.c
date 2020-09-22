#ifndef TOPOLOGICAL_ACTION_C
#define TOPOLOGICAL_ACTION_C
#include "./topological_action.h"
#include "./action.h"


double topodynamical_pot(double Q)
{
	double barr=act_params.barrier, width=act_params.width;
	int ngrid=(2*barr+width/2)/width;
	int igrid=floor((Q+barr)/width);
  
	if(igrid>=0 && igrid<ngrid)
		{
			double *grid = (double *)malloc(ngrid*sizeof(double));
			int check_load=load_topo(act_params.topo_file_path,barr,width,grid,ngrid);

			if(check_load!=0)
				{
					printf("ERROR: topological file not found");
					exit(1);
				}
			
			//interpolate
			double x0=igrid*width-barr;
			double m=(grid[igrid+1]-grid[igrid])/width;
			double q=grid[igrid]-m*x0;
			free(grid);
			return q+m*Q;
		}
	else
		return 0;
}




double compute_topo_action(su3_soa * const u)
{
	su3_soa * const quadri, tstout_conf_acc_arr;
	double_soa * const loc_q;
	

#ifdef STOUT_TOPO
	stout_wrapper(u,tstout_conf_acc_arr); //INSERIRE PARAMETRI STOUTING TOPOLOGICO
#else
	tstout_conf_acc_arr=*u;
#endif
	
	double Q = compute_topological_charge(u,quadri,loc_q);

	double topo_action = topodynamical_pot(Q);

	//PROBABILI PROBLEMI DI ALLOCAZIONE GENERATI
  
	return topo_action;
}

#endif
