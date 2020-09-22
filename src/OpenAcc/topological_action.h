#ifndef TOPOLOGICAL_ACTION_H
#define TOPOLOGICAL_ACTION_H

#include <stdlib.h>
#include "../Meas/gauge_meas.h"
#include "./alloc_vars.h"
#include "./action.h"


#ifndef STOUTING_H
 #include "./stouting.h"
#endif

//read COPIATO DA NISSA													\
SI DEVE CREARE UN FILE CHE CONTENGA UN INSIEME DI VALORI CHE DEFINISCONO IL POTENZIALE

static int load_topo(const char *path, double barr, double width, double* grid, int ngrid)
{
	FILE *fin=fopen(path,"r");
	
	if(!fin)
		{
			printf("%s file not found\n",path);
			return 1;
		}
	
	for(int igrid=0;igrid<=ngrid;igrid++)
		{
			double xread;
			int rc=fscanf(fin,"%lg %lg",&xread,&grid[igrid]);
			if(rc!=2) {printf("reading line %d of \"%s\"",igrid,path); return 2;}
			int jgrid=floor((xread+barr+width/2)/width);
			if(igrid!=jgrid) {printf("found %d (%lg) when expecting %d",jgrid,xread,igrid);return 2;}
		}
	fclose(fin);
  
	return 0;
}


static double topodynamical_pot(double Q);


static double compute_topo_action(su3_soa * const u);




#endif
