#ifndef FERMION_PARAMETERS_C_
#define FERMION_PARAMETERS_C_

#include "./fermion_parameters.h"
#include "../OpenAcc/backfield.h"
#include "../OpenAcc/alloc_vars.h"
#include "./markowchain.h"
#include <string.h>
#include <math.h>


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

    quark->approx_fi_mother.lambda_min = quark->ferm_mass*quark->ferm_mass/mkwch_pars.expected_max_eigenvalue;
    quark->approx_md_mother.lambda_min = quark->ferm_mass*quark->ferm_mass/mkwch_pars.expected_max_eigenvalue;
    quark->approx_li_mother.lambda_min = quark->ferm_mass*quark->ferm_mass/mkwch_pars.expected_max_eigenvalue;

    quark->approx_fi_mother.lambda_max =  1.0;
    quark->approx_md_mother.lambda_max =  1.0;
    quark->approx_li_mother.lambda_max =  1.0;

    quark->approx_fi_mother.error =  mkwch_pars.residue_metro/
        pow(mkwch_pars.expected_max_eigenvalue,
                (double) quark->approx_fi_mother.exponent_num/
                quark->approx_fi_mother.exponent_den );
    quark->approx_md_mother.error =  mkwch_pars.residue_md/
        pow(mkwch_pars.expected_max_eigenvalue,
                (double) quark->approx_md_mother.exponent_num/
                quark->approx_md_mother.exponent_den );
        
    quark->approx_li_mother.error =  mkwch_pars.residue_metro/
        pow(mkwch_pars.expected_max_eigenvalue,
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

  return errorstatus;

}


void init_all_u1_phases(bf_param bfpars, ferm_param *fpar  )
{
  for(int i=0;i<NDiffFlavs;i++){
      fpar[i].phases = &u1_back_phases[i*8];
      init_fermion_backfield(bfpars,&(fpar[i]));
  }
}


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
    double chpotphase = fermion_parameters->ferm_im_chem_pot/nt; 

    double_soa * phases = fermion_parameters->phases;

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
                    if(X == nx){// x-oriented links on the boundary
                        arg = (y+1)*nx*bz_quantum/(nx*ny);
                        arg+= (t+1)*nx*ex_quantum/(nx*nt);
                        arg-= (z+1)*by_quantum/(nz*nx);
                    }
                    else arg = -(z+1)*by_quantum/(nz*nx);

                    arg *= ferm_charge;// only em phase so far
                    if(KSphaseX(x,y,z,t) == -1) arg += 0.5;
                    // ^this should be false, UNLESS one chooses
                    // a different setup for the staggered phases....
                    while(arg > 0.5) arg -= 1.0;
                    while(arg < -0.5) arg += 1.0;

                    phases[parity].d[idxh]= acc_twopi*arg; 

                    ////////Y-oriented/////////
                    // Y even --> dir=2
                    // Y odd  --> dir=3
                    if(Y == ny){// y-oriented links on the boundary
                        arg = (z+1)*ny*bx_quantum/(ny*nz);
                        arg+= (t+1)*ny*ey_quantum/(ny*nt);
                        arg-= (x+1)*bz_quantum/(nx*ny);
                    }
                    else arg = -(x+1)*bz_quantum/(nx*ny);

                    arg *= ferm_charge;// only am phase so far
                    if(KSphaseY(x,y,z,t) == -1) arg += 0.5;
                    while(arg > 0.5) arg -= 1.0;
                    while(arg < -0.5) arg += 1.0;

                    phases[2+parity].d[idxh]=acc_twopi*arg;
                    
                    ////////Z-oriented////////
                    // Z even --> dir=4
                    // Z odd  --> dir=5
                    if(Z == nz){// z-oriented links on the boundary
                        arg = (t+1)*nz*ez_quantum/(nz*nt);
                        arg += (x+1)*nz*by_quantum/(nz*nx);
                        arg -= (y+1)*bx_quantum/(ny*nz);
                    }
                    else arg = -(y+1)*bx_quantum/(ny*nz);

                    arg *= ferm_charge;// only am phase so far
                    if(KSphaseZ(x,y,z,t) == -1) arg += 0.5;
                    while(arg > 0.5) arg -= 1.0;
                    while(arg < -0.5) arg += 1.0;

                    phases[4+parity].d[idxh]=acc_twopi*arg;

                    ///////T-oriented////////
                    // T even --> dir=6
                    // T odd  --> dir=7
                    arg = -(z+1)*ez_quantum/(nz*nt);
                    arg -= (y+1)*ey_quantum/(ny*nt);
                    arg -= (x+1)*ex_quantum/(nx*nt);
                   
                    arg *= ferm_charge;// only am phase so far
                    if(KSphaseZ(x,y,z,t) == -1) arg += 0.5;
                    arg += chpotphase; 
                    if(T == nt) arg += 0.5;// antiperiodic boundary conds

                    while(arg > 0.5) arg -= 1.0;
                    while(arg < -0.5) arg += 1.0;

                    phases[6+parity].d[idxh]=acc_twopi*arg;
                } // x loop
            } // y loop
        } // z loop
    } // t loop

    int dir;
    for(dir = 0; dir < 8; dir++)
        phases[dir].staggered_status = 1;

    
}
#endif
