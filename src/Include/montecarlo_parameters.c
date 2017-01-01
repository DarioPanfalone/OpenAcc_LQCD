#ifndef _MONTECARLO_PARAMS_C_
#define _MONTECARLO_PARAMS_C_

#include <stdio.h>

#include "./montecarlo_parameters.h"
#include "./common_defines.h"
#include "../Mpi/multidev.h"

mc_params_t mc_params;

void init_global_program_status(){
    
    printf("MPI%02d: reading global program status  from file %s\n",
            devinfo.myrank, mc_params.statusFileName);

    FILE * gps_file = fopen(mc_params.statusFileName, "r");  
    if(gps_file){
        fscanf(gps_file,"%d %lf %lf %d",
                &mc_params.next_gps,
                &mc_params.max_flavour_cycle_time,
                &mc_params.max_update_time,
                &mc_params.measures_done);

        fclose(gps_file);
    }
    else 
    {
        printf("MPI%02d: file %s not readable\n",
            devinfo.myrank, mc_params.statusFileName);

        mc_params.next_gps  = GPSTATUS_UPDATE;
        mc_params.max_flavour_cycle_time = 1.0;
        mc_params.max_update_time = 1.0;
        mc_params.measures_done = 0 ;
    }

    mc_params.run_condition = RUN_CONDITION_GO; 

    printf("%d %lf %lf %d\n",
            mc_params.next_gps,
            mc_params.max_flavour_cycle_time,
            mc_params.max_update_time,
            mc_params.measures_done);

    printf("#mc_params.next_gps,mc_params.max_flavour_cycle_time,\n\
#mc_params.max_update_time,mc_params.measures_done\n");


}

void save_global_program_status(){

    printf("Saving global program status...\n");
    printf("%d %lf %lf %d\n",
            mc_params.next_gps,
            mc_params.max_flavour_cycle_time,
            mc_params.max_update_time,
            mc_params.measures_done);

    printf("#mc_params.next_gps,mc_params.max_flavour_cycle_time,\n\
#mc_params.max_update_time,mc_params.measures_done\n");




    FILE * gps_file = fopen(mc_params.statusFileName, "w");  
    fprintf(gps_file,"%d %lf %lf %d\n",
            mc_params.next_gps,
            mc_params.max_flavour_cycle_time,
            mc_params.max_update_time,
            mc_params.measures_done);

    fprintf(gps_file,"#mc_params.next_gps,mc_params.max_flavour_cycle_time,\n\
#mc_params.max_update_time,mc_params.measures_done\n");




    fclose(gps_file);


}




#endif
