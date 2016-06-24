#ifndef FERMION_PARAMETERS_C_
#define FERMION_PARAMETERS_C_


#include "../DbgTools/dbgtools.h"
#include "../Mpi/multidev.h"
#include "../OpenAcc/alloc_vars.h"
#include "../OpenAcc/backfield.h"
#include "../OpenAcc/geometry.h"
#include "../OpenAcc/md_integrator.h"
#include "../RationalApprox/rationalapprox.h"
#include "./fermion_parameters.h"
#include "./montecarlo_parameters.h"
#include "../Include/debug.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

#define ALIGN 128
#define acc_twopi 2*3.14159265358979323846

int NDiffFlavs;// set in init.c, from input file
int NPS_tot;
int max_ps;
ferm_param *fermions_parameters;// set in init.c, from input file

int init_ferm_params(ferm_param *fermion_settings){

    int errorstatus = 0;


    printf("Initializing fermions...\n");

    NPS_tot = 0;
    max_ps = 0;

    // calculation of NPS_tot, max_ps,index_of_the_first_ps; 
    for(int i=0;i<NDiffFlavs;i++){
        // compute the total number of ps
        NPS_tot += fermion_settings[i].number_of_ps;
        // compute the max number of ps among the various flavs
        if(fermion_settings[i].number_of_ps>=max_ps) max_ps = fermion_settings[i].number_of_ps;
        // determine the offset (where does the ps of the flavour i starts?)
        if(i==0){
            fermion_settings[i].index_of_the_first_ps=0;
        }else{
            fermion_settings[i].index_of_the_first_ps = fermion_settings[i-1].index_of_the_first_ps + fermion_settings[i-1].number_of_ps;
        }
    }

    printf("NPS_tot = %d \n",NPS_tot);
    printf("max_ps = %d \n",max_ps);

    // Rational Approximation related stuff
    for(int i=0;i<NDiffFlavs;i++){
        ferm_param *quark = &fermion_settings[i];
        quark->approx_fi_mother.exponent_num =  +quark->degeneracy;
        quark->approx_md_mother.exponent_num =  -quark->degeneracy;
        quark->approx_li_mother.exponent_num =  -quark->degeneracy;

        quark->approx_fi_mother.exponent_den =   quark->number_of_ps*8;
        quark->approx_md_mother.exponent_den =   quark->number_of_ps*4;
        quark->approx_li_mother.exponent_den =   quark->number_of_ps*4;

        quark->approx_fi_mother.lambda_min = quark->ferm_mass*quark->ferm_mass/md_parameters.expected_max_eigenvalue;
        quark->approx_md_mother.lambda_min = quark->ferm_mass*quark->ferm_mass/md_parameters.expected_max_eigenvalue;
        quark->approx_li_mother.lambda_min = quark->ferm_mass*quark->ferm_mass/md_parameters.expected_max_eigenvalue;

        quark->approx_fi_mother.lambda_max =  1.0;
        quark->approx_md_mother.lambda_max =  1.0;
        quark->approx_li_mother.lambda_max =  1.0;

        quark->approx_fi_mother.error =  md_parameters.residue_metro/
            pow(md_parameters.expected_max_eigenvalue,
                    (double) quark->approx_fi_mother.exponent_num/
                    quark->approx_fi_mother.exponent_den );
        quark->approx_md_mother.error =  md_parameters.residue_md/
            pow(md_parameters.expected_max_eigenvalue,
                    (double) quark->approx_md_mother.exponent_num/
                    quark->approx_md_mother.exponent_den );

        quark->approx_li_mother.error =  md_parameters.residue_metro/
            pow(md_parameters.expected_max_eigenvalue,
                    (double) quark->approx_li_mother.exponent_num/
                    quark->approx_li_mother.exponent_den );

        quark->approx_fi_mother.gmp_remez_precision = 100;
        quark->approx_md_mother.gmp_remez_precision = 100;
        quark->approx_li_mother.gmp_remez_precision = 100;

        // copy everything also in the daughter approxs
        quark->approx_fi.exponent_num =   quark->approx_fi_mother.exponent_num;
        quark->approx_md.exponent_num =   quark->approx_md_mother.exponent_num;
        quark->approx_li.exponent_num =   quark->approx_li_mother.exponent_num;
        quark->approx_fi.exponent_den =   quark->approx_fi_mother.exponent_den;
        quark->approx_md.exponent_den =   quark->approx_md_mother.exponent_den;
        quark->approx_li.exponent_den =   quark->approx_li_mother.exponent_den;
        quark->approx_fi.approx_order =   quark->approx_fi_mother.approx_order;
        quark->approx_md.approx_order =   quark->approx_md_mother.approx_order;
        quark->approx_li.approx_order =   quark->approx_li_mother.approx_order;
        quark->approx_fi.gmp_remez_precision =
            quark->approx_fi_mother.gmp_remez_precision;
        quark->approx_md.gmp_remez_precision =   
            quark->approx_md_mother.gmp_remez_precision;
        quark->approx_li.gmp_remez_precision =   
            quark->approx_li_mother.gmp_remez_precision;

        // READ THE RAT APPROXS FROM THE FILES
        errorstatus += rationalapprox_read(&(quark->approx_fi_mother));
        errorstatus += rationalapprox_read(&(quark->approx_md_mother));
        errorstatus += rationalapprox_read(&(quark->approx_li_mother));

    }

    fermion_settings->printed_bf_dbg_info = 0;
    return errorstatus;

}


void init_all_u1_phases(bf_param bfpars, ferm_param *fpar  )
{

    for(int i=0;i<NDiffFlavs;i++){
        fpar[i].phases = &u1_back_phases[i*8];
        init_fermion_backfield(bfpars,&(fpar[i]));

        // PRINTING DEBUG INFO
        if(debug_settings.print_bfield_dbginfo){
            char tempname[50];                           
            // phases
            sprintf(tempname,"backfield_%s_c%d",fpar[i].name,fpar[i].printed_bf_dbg_info);
#ifdef MULTIDEVICE      
            strcat(tempname,devinfo.myrankstr);          
#endif
            print_double_soa(fpar[i].phases,tempname);   

            // plaquettes
            sprintf(tempname,"abelian_plq_%s_c%d",fpar[i].name,
                    fpar[i].printed_bf_dbg_info);
#ifdef MULTIDEVICE      
            strcat(tempname,devinfo.myrankstr);          
#endif

            print_all_abelian_plaquettes(fpar[i].phases,tempname);
            fpar->printed_bf_dbg_info++; 
        }
    }


}

void init_fermion_backfield(bf_param bf_pars, ferm_param *fermion_parameters)
{

    if(verbosity_lv > 2 && 0 == devinfo.myrank ) { 
        printf("Generating external field (containing staggered phases) ");
        printf("for flavour %s\n",fermion_parameters->name);
    }


    calc_u1_phases(fermion_parameters->phases, bf_pars, 
            fermion_parameters->ferm_im_chem_pot, fermion_parameters->ferm_charge);



}


/*

   
   void init_fermion_backfield(bf_param bf_pars, ferm_param *fermion_parameters){

    double  ex_quantum = bf_pars.ex;
    double  ey_quantum = bf_pars.ey;
    double  ez_quantum = bf_pars.ez;
    double  bx_quantum = bf_pars.bx;
    double  by_quantum = bf_pars.by;
    double  bz_quantum = bf_pars.bz;

    int X,Y,Z,T;
    double  arg;
    double ferm_charge = fermion_parameters->ferm_charge;
    double chpotphase = fermion_parameters->ferm_im_chem_pot/geom_par.gnt; 

    double_soa * phases = fermion_parameters->phases;

    int x, y, z, t, parity;
    int d[4], idxh;

    if(verbosity_lv > 2 && 0 == devinfo.myrank ) { 
        printf("Generating external field (containing staggered phases) ");
        printf("for flavour %s\n",fermion_parameters->name);
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


*/


#endif
