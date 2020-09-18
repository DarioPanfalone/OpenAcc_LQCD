#ifndef TOPO_ACTION_C
#define TOPO_ACTION_C
#include "./topological_action.h"
#include "./action.h"


int stout_passi; //MOMENTANEO

//read COPIATO DA NISSA \
SI DEVE CREARE UN FILE CHE CONTENGA UN INSIEME DI VALORI CHE DEFINISCONO IL POTENZIALE

int load(const char *path, double barr, double width, double* grid, int ngrid)
{//DOMANDA: COSÌ CI SCRIVE SU GRID O NO?
	//QUI C'ERA UN GET_THREAD_ID()

	FILE *fin=fopen(path,"r");
  
	if(!fin){
		printf("%s file not found\n",path);
		return 1; //QUESTO DOVREBBE INDICARE UN ERRORE MA NON MI È BEN CHIARO COME FUNZIONI LA CODIFICA
	}
  
	//if(IS_MASTER_THREAD && rank==0) NON AVENDO GET_THREAD_ID, NON HA SENSO STO IF
	for(int igrid=0;igrid<=ngrid;igrid++)
		{
			double xread;
			int rc=fscanf(fin,"%lg %lg",&xread,&grid[igrid]);
			if(rc!=2) {printf("reading line %d of \"%s\"",igrid,path); return 2;}
			int jgrid=floor((xread+barr+width/2)/width);
			if(igrid!=jgrid) {printf("found %d (%lg) when expecting %d",jgrid,xread,igrid);return 2;}
		}
	fclose(fin);
  
	//broadcast DEVO ANCORA CAPIRE COME FUNZIONA QUESTA COSA IN OPENACC.
	//for(int igrid=0;igrid<=ngrid;igrid++) MPI_Bcast(&grid[igrid],1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	return 0;
}

double topodynamical_pot(double Q)
{
	double barr=act_params.barrier, width=act_params.width;
	int ngrid=(2*barr+width/2)/width;
	int igrid=floor((Q+barr)/width);
  
	if(igrid>=0 && igrid<ngrid)
		{
			double *grid = (double *)malloc(ngrid*sizeof(double));
			load(act_params.topo_file_path,barr,width,grid,ngrid);
      
			//interpolate
			double x0=igrid*width-barr;
			double m=(grid[igrid+1]-grid[igrid])/width;
			double q=grid[igrid]-m*x0;
			free(grid);
			return q+m*Q;
		}
	else
		return 0; //IN NISSA force_out*(roba) (MA NEL MUTLICANONICO È SOLITAMENTE 0) 
}




double compute_topo_action(su3_soa * const u)
{
	su3_soa * const quadri, tstout_conf_acc_arr;
	double_soa * const loc_q;
  
	//metti if viene richiesto lo stouting:
#ifdef STOUT_TOPO
	stout_wrapper(u,tstout_conf_acc_arr); //INSERIRE PARAMETRI STOUTING TOPOLOGICO
#else
	tstout_conf_acc_arr=*u;
#endif
  
	double Q = compute_topological_charge(u,quadri,loc_q);

	double topo_action = topodynamical_pot(act_params.topo_file_path, Q, act_params.barrier, act_params.width);

	//PROBABILI PROBLEMI DI ALLOCAZIONE GENERATI
  
	return topo_action;
}



double compute_topo_force(su3_soa * const u)
{
	su3_soa * const quadri, tstout_conf_acc_arr;
	double_soa * const loc_q;
  
#ifdef STOUT_TOPO
	stout_wrapper(u,tstout_conf_acc_arr); //INSERIRE PARAMETRI STOUTING TOPOLOGICO
#else
	tstout_conf_acc_arr=*u;
#endif

	double Q = compute_topological_charge(u, quadri, loc_q);
	double pot_der=topodynamical_pot_der(act_params.topo_file_path, Q);
  
  
}
#endif
