#ifndef BACKFIELD_C_
#define BACKFIELD_C_

#include "./geometry.h"
#include "./backfield.h"
#define acc_twopi 2*3.14159265358979323846
#ifndef __GNUC__
#include <accelmath.h>
#endif
#include "./struct_c_def.h"
#include "../Include/fermion_parameters.h"
#include "./geometry.h"

bf_param backfield_parameters;



void init_backfield(double_soa * tu1_back_field_phases,bf_param bf_pars ){

    double  ex_quantum = bf_pars.ex;
    double  ey_quantum = bf_pars.ey;
    double  ez_quantum = bf_pars.ez;
    double  bx_quantum = bf_pars.bx;
    double  by_quantum = bf_pars.by;
    double  bz_quantum = bf_pars.bz;

    int X,Z,Y;
    double  arg;

    int x, y, z, t, parity,idxh;
    for(t=0; t < nt; t++){
        for(z=0; z < nz; z++){
            for(y=0; y < ny; y++){
                for(x=0; x < nx; x++){

                    if((x+y+z+t)%2==0){
                        parity = 0; //pari
                    }else{
                        parity = 1; //dispari
                    }
                    idxh = snum_acc(x,y,z,t);

                    X = x + 1;
                    Y = y + 1;
                    Z = z + 1;

                    ////////X-oriented////////
                    // X even --> dir=0
                    // X odd  --> dir=1
                    if(X == nx){             // x-oriented links on the boundary
                        arg = ((double)(y+1))*((double)nx)*((double)bz_quantum)/(((double)nx)*((double)ny));
                        arg += ((double)(t+1))*((double)nx)*((double)ex_quantum)/(((double)nx)*((double)nt));
                        arg -= ((double)(z+1))*((double)by_quantum)/(((double)nz)*((double)nx));
                        tu1_back_field_phases[parity].d[idxh]= acc_twopi*arg; 
                    }
                    else {
                        arg = -((double)(z+1))*((double)by_quantum)/(((double)nz)*((double)nx));
                        tu1_back_field_phases[parity].d[idxh]= acc_twopi*arg; 
                    }

                    ////////Y-oriented/////////
                    // Y even --> dir=2
                    // Y odd  --> dir=3
                    if(Y == ny){             // y-oriented links on the boundary
                        arg = ((double)(z+1))*((double)ny)*((double)bx_quantum)/(((double)ny)*((double)nz));
                        arg += ((double)(t+1))*((double)ny)*((double)ey_quantum)/(((double)ny)*((double)nt));
                        arg -= ((double)(x+1))*((double)bz_quantum)/(((double)nx)*((double)ny));
                        tu1_back_field_phases[2+parity].d[idxh]= acc_twopi*arg;
                    }
                    else {
                        arg = -((double)(x+1))*((double)bz_quantum)/(((double)nx)*((double)ny));
                        tu1_back_field_phases[2+parity].d[idxh]= acc_twopi*arg;
                    }

                    ////////Z-oriented////////
                    // Z even --> dir=4
                    // Z odd  --> dir=5
                    if(Z == nz){          // z-oriented links on the boundary
                        arg = ((double)(t+1))*((double)nz)*((double)ez_quantum)/(((double)nz)*((double)nt));
                        arg += ((double)(x+1))*((double)nz)*((double)by_quantum)/(((double)nz)*((double)nx));
                        arg -= ((double)(y+1))*((double)bx_quantum)/(((double)ny)*((double)nz));
                        tu1_back_field_phases[4+parity].d[idxh]=acc_twopi* arg;
                    }
                    else{
                        arg = -((double)(y+1))*((double)bx_quantum)/(((double)ny)*((double)nz));
                        tu1_back_field_phases[4+parity].d[idxh]= acc_twopi*arg;
                    }

                    ///////T-oriented////////
                    // T even --> dir=6
                    // T odd  --> dir=7
                    arg = -((double)(z+1))*((double)ez_quantum)/(((double)nz)*((double)nt));
                    arg -= ((double)(y+1))*((double)ey_quantum)/(((double)ny)*((double)nt));
                    arg -= ((double)(x+1))*((double)ex_quantum)/(((double)nx)*((double)nt));
                    tu1_back_field_phases[6+parity].d[idxh]= acc_twopi*arg;
                } // x loop
            } // y loop
        } // z loop
    } // t loop
}



void mult_backfield_by_stag_phases(double_soa * phases){

    int X,Y,Z,T;
    double  arg2pi;
    const double pi = 0.5*acc_twopi;

    int x, y, z, t, parity;
    int d[4], nd[4], idxh;
    nd[0] = nd0;
    nd[1] = nd1;
    nd[2] = nd2;
    nd[3] = nd3;

    for(d[3]=0; d[3] < nd3; d[3]++){
        for(d[2]=0; d[2] < nd2; d[2]++){
            for(d[1]=0; d[1] < nd1; d[1]++){
                for(d[0]=0; d[0] < nd0; d[0]++){

                    idxh = snum_acc(d[0],d[1],d[2],d[3]);
                   
                    x = d[geom_par.xmap];
                    y = d[geom_par.ymap];
                    z = d[geom_par.zmap];
                    t = d[geom_par.tmap];

                    if((x+y+z+t)%2==0) parity = 0; //pari
                    else parity = 1; //dispari
                

                    X = x + 1;
                    Y = y + 1;
                    Z = z + 1;
                    T = t + 1;

                    ////////X-oriented////////
                    // X even --> dir=0
                    // X odd  --> dir=1
                    arg2pi = phases[parity].d[idxh]; 
                    if(KSphaseX(x,y,z,t) == -1) arg2pi += pi;
                    // ^this should be false, UNLESS one chooses
                    // a different setup for the staggered phases....
                    while(arg2pi >  pi) arg2pi -= 2*pi;
                    while(arg2pi < -pi) arg2pi += 2*pi;

                    phases[parity].d[idxh]= arg2pi; 

                    ////////Y-oriented/////////
                    // Y even --> dir=2
                    // Y odd  --> dir=3

                    arg2pi = phases[2+parity].d[idxh]; 
                    if(KSphaseY(x,y,z,t) == -1) arg2pi += pi;
                    while(arg2pi >  pi) arg2pi -= 2*pi;
                    while(arg2pi < -pi) arg2pi += 2*pi;

                    phases[2+parity].d[idxh]= arg2pi; 
                    
                    ////////Z-oriented////////
                    // Z even --> dir=4
                    // Z odd  --> dir=5

                    arg2pi = phases[4+parity].d[idxh]; 
                    if(KSphaseZ(x,y,z,t) == -1) arg2pi += pi;
                    while(arg2pi >  pi) arg2pi -= 2*pi;
                    while(arg2pi < -pi) arg2pi += 2*pi;

                    phases[4+parity].d[idxh]= arg2pi; 


                    ///////T-oriented////////
                    // T even --> dir=6
                    // T odd  --> dir=7
                   
                    if(KSphaseZ(x,y,z,t) == -1) arg2pi += pi;
                    if(T == nt) arg2pi += pi; // antiperiodic boundary conds
                    while(arg2pi >  pi) arg2pi -= 2*pi;
                    while(arg2pi < -pi) arg2pi += 2*pi;

                    phases[6+parity].d[idxh]=arg2pi;
                } // x loop
            } // y loop
        } // z loop
    } // t loop
}


// DEBUG
#include <stdio.h>
void print_double_soa(double_soa * const backfield, const char* nomefile){

    printf("Printing backfield on file %s", nomefile);

    FILE *fp;
    fp = fopen(nomefile,"w");
    for(int q = 0 ; q < 8 ; q++){
        for(int i = 0 ; i < sizeh ; i++){
            fprintf(fp, "%.18lf\n",backfield[q].d[i]);
        }
    }
    fclose(fp);
}


#endif
