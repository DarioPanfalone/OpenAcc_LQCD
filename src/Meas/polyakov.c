#ifndef POLYAKOV_C_
#define POLYAKOV_C_

#include "../OpenAcc/su3_utilities.h"
#include "./polyakov.h"
#include "../OpenAcc/su3_measurements.h"


double polyakov_loop(	__restrict su3_soa * const u,
        __restrict su3_soa * const loopplk			)
{	

    if(geom_par.tmap != 3)
    printf("WARNING: tmap != 3 , this function won't work!\n");
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
    int size_plkv = vol3s[geom_par.tmap];
    d_complex* plkv = (d_complex*) 
        malloc(size_plkv*2*sizeof(double));
    int d0, d1, d2, d3,h,idxh,id_zero,parity,parity_zero;
    double r = 0;
#pragma acc kernels present(u) present(loopplk)
#pragma acc loop independent gang //gang(nd3)
    for(d3=0; d3<nd3; d3++) {
#pragma acc loop independent gang vector //gang(nd2/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
        for(d2=0; d2<nd2; d2++) {
#pragma acc loop independent gang vector //gang(nd1/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
            for(d1=0; d1<nd1; d1++) {
#pragma acc loop independent vector //vector(DIM_BLOCK_X)
                for(d0=0; d0 < nd0; d0++) {	     
                    parity_zero = (d0+d1+d2+d3 - nd[geom_par.tmap]) % 2; 
                    id_zero = snum_acc(d0,d1,d2,0);  // first //PROBLEM 
                    parity = (d0+d1+d2+d3) % 2;
                    idxh = snum_acc(d0,d1,d2,d3);  	
                    mat1_times_mat2_into_mat1_absent_stag_phases(&loopplk[parity_zero],id_zero,&u[geom_par.tmap*2+parity],idxh);              // LOOP = LOOP * C
                }}}}


    int x,y,z;
#pragma acc kernels present(u) present(loopplk) copyout(plkv[0:size_plkv])
#pragma acc loop independent gang //reduction(+:r)
    for(z=0; z<nd[geom_par.zmap]; z++) {
#pragma acc loop independent gang vector 
        for(y=0; y<nd[geom_par.ymap]; y++) {
#pragma acc loop independent vector 
            for(x=0; x < nd[geom_par.xmap]; x++) {

                h=x+(y*vol1)+(z*vol2); 			// PROBLEM
                parity_zero = (d0+d1+d2+0) % 2; // PROBLEM
                id_zero = snum_acc(d0,d1,d2,0); // PROBLEM
                // calcolo traccia
                plkv[h] = matrix_trace_absent_stag_phase(&loopplk[parity_zero],id_zero);    
                //controllo unitarietà
                /*   single_su3 m;
                     single_su3_from_su3_soa(&loopplk[parity_zero],id_zero,&m);
                     rebuild3row(&m);
                     id_zero = snum_acc(d0,d1,d2,0);                         
                     d_complex err = 1 - detSu3(&m);
                     r += creal(err * conj(err));*/	
            }
        }
    }								  


    double res_L = 0.0;
    for (d3=0; d3<vol3; d3++){
        res_L += creal(plkv[d3]); //sum all polyakov loops			                  
    }
    res_L=res_L/vol3s[geom_par.tmap];
    res_L=res_L/3; //divide by 3 color

    //  printf("\d3 Errore su unitarieta polyakov: %f\n",r);


    return res_L;
}
#endif
