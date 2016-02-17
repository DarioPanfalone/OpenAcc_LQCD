#ifndef POLYAKOV_C_
#define POLYAKOV_C_

#include "../OpenAcc/su3_utilities.h"
#include "./polyakov.h"
#include "../OpenAcc/su3_measurements.h"

#define dir_tempo 3
#define dir_tempo_link 6

double polyakov_loop(	__restrict su3_soa * const u,
        __restrict su3_soa * const loopplk			)
{	
    // removing stag phases
    mult_conf_times_stag_phases(u);

    //set to identity
#pragma acc kernels present(loopplk)
#pragma acc loop independent gang  
    for(int j=0;j<2;j++){
#pragma acc loop independent vector
        for(int i=0;i<sizeh;i++){
            loopplk[j].r0.c0[i]=1;
            loopplk[j].r0.c1[i]=0;
            loopplk[j].r0.c2[i]=0;
            loopplk[j].r1.c0[i]=0;
            loopplk[j].r1.c1[i]=1;
            loopplk[j].r1.c2[i]=0;
            loopplk[j].r2.c0[i]=0;
            loopplk[j].r2.c1[i]=0;
            loopplk[j].r2.c2[i]=1;
        }}

    d_complex plkv[vol3];
    int x, y, z, t,h,idxh,id_zero,parity,parity_zero;
    double r = 0;
#pragma acc kernels present(u) present(loopplk)
#pragma acc loop independent gang //gang(nt)
    for(t=0; t<nt; t++) {
#pragma acc loop independent gang vector //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
        for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
            for(y=0; y<ny; y++) {
#pragma acc loop independent vector //vector(DIM_BLOCK_X)
                for(x=0; x < nx; x++) {	     
                    parity_zero = (x+y+z+0) % 2; 
                    id_zero = snum_acc(x,y,z,0);  // first 			        
                    parity = (x+y+z+t) % 2;
                    idxh = snum_acc(x,y,z,t);  	
                    mat1_times_mat2_into_mat1_absent_stag_phases(&loopplk[parity_zero],id_zero,&u[dir_tempo_link+parity],idxh);              // LOOP = LOOP * C
                }}}}


#pragma acc kernels present(u) present(loopplk) copyout(plkv)
#pragma acc loop independent gang //reduction(+:r)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang vector 
        for(y=0; y<ny; y++) {
#pragma acc loop independent vector 
            for(x=0; x < nx; x++) {

                h=x+(y*vol1)+(z*vol2); 			
                parity_zero = (x+y+z+0) % 2; 
                id_zero = snum_acc(x,y,z,0);
                // calcolo traccia
                plkv[h] = matrix_trace_absent_stag_phase(&loopplk[parity_zero],id_zero);    
                //controllo unitarietà
                /*   single_su3 m;
                     single_su3_from_su3_soa(&loopplk[parity_zero],id_zero,&m);
                     rebuild3row(&m);
                     id_zero = snum_acc(x,y,z,0);                         
                     d_complex err = 1 - detSu3(&m);
                     r += creal(err * conj(err));*/	
            }}}								  


    double res_L = 0.0;
    for (t=0; t<vol3; t++){
        res_L += creal(plkv[t]); //sum all polyakov loops			                  
    }
    res_L=res_L/vol3;
    res_L=res_L/3; //divide by 3 color

    //  printf("\t Errore su unitarieta polyakov: %f\n",r);


    // add stag phases
    mult_conf_times_stag_phases(u);

    return res_L;
}
#endif
