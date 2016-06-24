#ifndef BACKFIELD_C_
#define BACKFIELD_C_

#include "./geometry.h"
#include "./backfield.h"
#include "./struct_c_def.h"
#include "../Mpi/multidev.h"
#include <stdio.h>


#define acc_twopi 2*3.14159265358979323846
#ifndef __GNUC__
#include <accelmath.h>
#endif


bf_param backfield_parameters;

void calc_u1_phases(double_soa * phases,bf_param bf_pars,
        double im_chem_pot, double ferm_charge)
{
    double  ex_quantum = bf_pars.ex;
    double  ey_quantum = bf_pars.ey;
    double  ez_quantum = bf_pars.ez;
    double  bx_quantum = bf_pars.bx;
    double  by_quantum = bf_pars.by;
    double  bz_quantum = bf_pars.bz;


    int X,Y,Z,T;
    double  arg;
    double chpotphase = im_chem_pot/geom_par.gnt; 

    int x, y, z, t, parity;
    int d[4], idxh;

    if(verbosity_lv > 2 && 0 == devinfo.myrank ) { 
        printf("Direction mapping  x y z t: %d %d %d %d\n",
                geom_par.xmap,geom_par.ymap,geom_par.zmap,geom_par.tmap);
        printf("EM field quanta: \n\tB: %f %f %f \n\tE: %f %f %f\n",
                bx_quantum, by_quantum,bz_quantum,ex_quantum,ey_quantum,ez_quantum);

    }

#ifdef MULTIDEVICE

    printf("MPI%02d: Origin coordinates: (x=%02d y=%02d z=%02d t=%02d)\n",
            devinfo.myrank ,  devinfo.origin_0123[geom_par.xmap],
            devinfo.origin_0123[geom_par.ymap], devinfo.origin_0123[geom_par.zmap],
            devinfo.origin_0123[geom_par.tmap]);


    if(0==devinfo.myrank) printf("Halo Widths: x: %d y: %d z: %d t :%d\n",
            devinfo.halo_widths0123[geom_par.xmap],
            devinfo.halo_widths0123[geom_par.ymap],
            devinfo.halo_widths0123[geom_par.zmap],
            devinfo.halo_widths0123[geom_par.tmap]);

#endif                  


    for(d[3]=0; d[3] < nd3; d[3]++) for(d[2]=0; d[2] < nd2; d[2]++)
        for(d[1]=0; d[1] < nd1; d[1]++) for(d[0]=0; d[0] < nd0; d[0]++){

            idxh = snum_acc(d[0],d[1],d[2],d[3]);

            x = d[geom_par.xmap];int tnx = geom_par.gnx;
            y = d[geom_par.ymap];int tny = geom_par.gny;
            z = d[geom_par.zmap];int tnz = geom_par.gnz;
            t = d[geom_par.tmap];int tnt = geom_par.gnt;

#ifdef MULTIDEVICE
            x+= devinfo.origin_0123[geom_par.xmap]
                - devinfo.halo_widths0123[geom_par.xmap];  
            y+= devinfo.origin_0123[geom_par.ymap]
                - devinfo.halo_widths0123[geom_par.ymap];
            z+= devinfo.origin_0123[geom_par.zmap]
                - devinfo.halo_widths0123[geom_par.zmap];
            t+= devinfo.origin_0123[geom_par.tmap]
                - devinfo.halo_widths0123[geom_par.tmap];

            if(x>tnx-1)  x-= tnx ; if(x<0)  x+= tnx ; 
            if(y>tny-1)  y-= tny ; if(y<0)  y+= tny ; 
            if(z>tnz-1)  z-= tnz ; if(z<0)  z+= tnz ; 
            if(t>tnt-1)  t-= tnt ; if(t<0)  t+= tnt ; 

#endif                  

            parity = (x+y+z+t)%2; 
            // NOTICE: (x+y+z+t)%2 =/= (nd0+nd1+nd2+nd3)%2


            X = x + 1;
            Y = y + 1;
            Z = z + 1;
            T = t + 1;

            ////////X-oriented////////
            if(X == tnx){
                // x-oriented links on the boundary
                arg = (y+1)*tnx*bz_quantum/(tnx*tny);
                arg+= (t+1)*tnx*ex_quantum/(tnx*tnt);
                arg-= (z+1)*by_quantum/(tnz*tnx);
            }
            else arg = -(z+1)*by_quantum/(tnz*tnx);

            arg *= ferm_charge;// only em phase so far
            if(KSphaseX(x,y,z,t) == -1) arg += 0.5;
            // ^this should be false, UNLESS one chooses
            // a different setup for the staggered phases....
            while(arg > 0.5) arg -= 1.0;
            while(arg < -0.5) arg += 1.0;

            phases[geom_par.xmap*2+parity].d[idxh]= acc_twopi*arg; 


            ////////Y-oriented/////////
            if(Y == tny){
                // y-oriented links on the boundary
                arg = (z+1)*tny*bx_quantum/(tny*tnz);
                arg+= (t+1)*tny*ey_quantum/(tny*tnt);
                arg-= (x+1)*bz_quantum/(tnx*tny);
            }
            else arg = -(x+1)*bz_quantum/(tnx*tny);

            arg *= ferm_charge;// only am phase so far
            if(KSphaseY(x,y,z,t) == -1) arg += 0.5;
            while(arg > 0.5) arg -= 1.0;
            while(arg < -0.5) arg += 1.0;

            phases[geom_par.ymap*2+parity].d[idxh]=acc_twopi*arg;


            ////////Z-oriented////////
            if(Z == tnz){
                // z-oriented links on the boundary
                arg = (t+1)*tnz*ez_quantum/(tnz*tnt);
                arg += (x+1)*tnz*by_quantum/(tnz*tnx);
                arg -= (y+1)*bx_quantum/(tny*tnz);
            }
            else arg = -(y+1)*bx_quantum/(tny*tnz);

            arg *= ferm_charge;// only am phase so far
            if(KSphaseZ(x,y,z,t) == -1) arg += 0.5;
            while(arg > 0.5) arg -= 1.0;
            while(arg < -0.5) arg += 1.0;

            phases[geom_par.zmap*2+parity].d[idxh]=acc_twopi*arg;


            ///////T-oriented////////
            arg = -(z+1)*ez_quantum/(tnz*tnt);
            arg -= (y+1)*ey_quantum/(tny*tnt);
            arg -= (x+1)*ex_quantum/(tnx*tnt);

            arg *= ferm_charge;// only am phase so far
            if(KSphaseT(x,y,z,t) == -1) arg += 0.5;
            arg += chpotphase*0.5; // it must be multiplied by pi, not 2pi
            if(T == tnt) arg += 0.5;
            // antiperiodic boundary conds
            // notice T = t+1 !!

            while(arg > 0.5) arg -= 1.0;
            while(arg < -0.5) arg += 1.0;

            phases[geom_par.tmap*2+parity].d[idxh]=acc_twopi*arg;
        } // d3,d2,d1,d0 loops







}


void phase_diff_in_place(double_soa * inout, double_soa * sottraendo){
    int dir,idxh;

    for(dir=0;dir<8;dir++)
        for(idxh=0;idxh<sizeh;idxh++)
            inout[dir].d[idxh] -= sottraendo[dir].d[idxh];

}

void u1_diff(double_soa * out_re, double_soa * out_im,
        double_soa * phase_p,double_soa * phase_m ){

     int dir,idxh;

    for(dir=0;dir<8;dir++)
        for(idxh=0;idxh<sizeh;idxh++){
            out_re[dir].d[idxh] = cos(phase_p[dir].d[idxh]) - cos(phase_m[dir].d[idxh]);
            out_im[dir].d[idxh] = sin(phase_p[dir].d[idxh]) - sin(phase_m[dir].d[idxh]);
        }
    
}


void set_double_soa_to_zero(double_soa * p){
    int dir,idxh;

    for(dir=0;dir<8;dir++)
        for(idxh=0;idxh<sizeh;idxh++)
            p[dir].d[idxh]=0;
    
}



/*
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
}*/


#endif

