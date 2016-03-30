#ifndef SU3_MEASUREMENTS_C_
#define SU3_MEASUREMENTS_C_

#include "../Include/common_defines.h"
#include "../OpenAcc/struct_c_def.h"
#include "./su3_measurements.h"
#include "./su3_utilities.h"
#include "./plaquettes.h"
#include "./single_types.h"
#include "math.h"
#ifndef __GNUC__
#include <openacc.h>
#endif

#ifdef MULTIDEVICE
#include <mpi.h>
#endif



void check_unitarity_device( __restrict su3_soa * const u, double * max_unitarity_deviation, double *avg_unitarity_deviation){


    // removing stag phases


    double r = 0;
    double rmax = 0;

#pragma acc kernels present(u)
#pragma acc loop reduction(+:r) 
    for(int idx = 0; idx < sizeh ; idx++){
        for(int dir = 0; dir < 8 ; dir++){
            single_su3 m;
            single_su3_from_su3_soa(&u[dir],idx,&m);
            rebuild3row(&m);

            d_complex err = 1 - detSu3(&m);
            r += creal(err * conj(err));
            //     rmax = fmax(rmax,creal(err * conj(err)));

        }
    }
    printf("Warning: max unitarity deviation NOT calculated (%s:%d)\n",
            __FILE__,__LINE__);
    //adding them again
    *avg_unitarity_deviation = r/(sizeh*8);
    *max_unitarity_deviation = rmax;


}
void check_unitarity_host( __restrict su3_soa * const u, double * max_unitarity_deviation, double *avg_unitarity_deviation){


    // removing stag phases


    double r = 0;
    double rmax = 0;

    for(int dir = 0; dir < 8 ; dir++){
        for(int idx = 0; idx < sizeh ; idx++){
            single_su3 m;
            single_su3_from_su3_soa(&u[dir],idx,&m);
            rebuild3row(&m);

            d_complex err = 1 - detSu3(&m);
            r += creal(err * conj(err));
            rmax = fmax(rmax,creal(err * conj(err)));

        }
    }
    //adding them again

    *avg_unitarity_deviation = r/(sizeh*8);
    *max_unitarity_deviation = rmax;


}


double calc_momenta_action( const __restrict thmat_soa * const mom,
        double_soa * tr_local,
        const int mu){
    int t;

#pragma acc kernels present(mom) present(tr_local)
#pragma acc loop independent //gang(nt)
    for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
        tr_local[0].d[t] = half_tr_thmat_squared(&mom[mu],t);
    }  // t


    double result=0.0;
    double total_result=0.0;
#pragma acc kernels present(tr_local)
#pragma acc loop reduction(+:result)
    for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
        result += tr_local[0].d[t];
    }


#ifdef MULTIDEVICE
     MPI_Allreduce((void*)&result,(void*)&total_result,
             1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
     total_result = result;
#endif
    return total_result;

}// closes routine


double  calc_plaquette_soloopenacc( __restrict  su3_soa * const tconf_acc, __restrict su3_soa * const local_plaqs, dcomplex_soa * const tr_local_plaqs){

    double result=0.0;
    double total_result=0.0;

    // calcolo il valore della plaquette sommata su tutti i siti a fissato piano mu-nu (6 possibili piani)
    for(int mu=0;mu<3;mu++){
        for(int nu=mu+1;nu<4;nu++){
            // sommo i 6 risultati in tempo
            result  += calc_loc_plaquettes_nnptrick(tconf_acc,local_plaqs,tr_local_plaqs,mu,nu);
        }
    }


#ifdef MULTIDEVICE
     MPI_Allreduce((void*)&result,(void*)&total_result,
             1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
     total_result = result;
#endif
    return total_result;

}



double calc_force_norm(const __restrict tamat_soa * tipdot){

    int t,mu;
    double result=0.0;
    double total_result=0.0;



#pragma acc kernels present(tipdot)
#pragma acc loop reduction(+:result)
    for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
#pragma acc loop reduction(+:result)
        for(mu=0; mu < 8; mu++) // CHECK IF WORKS
        {
            result += half_tr_tamat_squared(&tipdot[mu],t);
        }
    }  
#ifdef MULTIDEVICE
     MPI_Allreduce((void*)&result,(void*)&total_result,
             1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
     total_result = result;
#endif
    return sqrt(total_result);

};

double calc_diff_force_norm(const __restrict tamat_soa * tipdot,  
        const __restrict tamat_soa * tipdot_old){

    int t,mu;
    double result=0.0;
    double total_result=0.0;


#pragma acc kernels present(tipdot) present(tipdot_old)
#pragma acc loop reduction(+:result)
    for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
#pragma acc loop reduction(+:result)
        for(mu=0; mu < 8; mu++) // CHECK IF WORKS
        {
            double A = tipdot[mu].ic00[t]-tipdot_old[mu].ic00[t];
            double B = tipdot[mu].ic11[t]-tipdot_old[mu].ic11[t];
            d_complex C = tipdot[mu].c01[t]-tipdot_old[mu].c01[t]; 
            d_complex D = tipdot[mu].c02[t]-tipdot_old[mu].c02[t];
            d_complex E = tipdot[mu].c12[t]-tipdot_old[mu].c12[t];
            result +=  A*A + B*B + A*B 
                + creal(C)*creal(C) + cimag(C)*cimag(C) 
                + creal(D)*creal(D) + cimag(D)*cimag(D) 
                + creal(E)*creal(E) + cimag(E)*cimag(E);

        }
    }  

#ifdef MULTIDEVICE
     MPI_Allreduce((void*)&result,(void*)&total_result,
             1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
     total_result = result;
#endif
    return sqrt(total_result);

};

void copy_ipdot_into_old(
        const __restrict tamat_soa * tipdot,  
        __restrict tamat_soa * tipdot_old){

    // POSSIBLY SOME PROBLEMS HERE
    int t,mu;

#pragma acc kernels present(tipdot) present(tipdot_old)
#pragma acc loop
    for(t=(LNH_SIZEH-LOC_SIZEH)/2; t  < (LNH_SIZEH+LOC_SIZEH)/2; t++) {
#pragma acc loop
        for(mu=0; mu < 8; mu++) // CHECK IF WORKS
        {
            tipdot_old[mu].ic00[t] = tipdot[mu].ic00[t];
            tipdot_old[mu].ic11[t] = tipdot[mu].ic11[t];
            tipdot_old[mu].c01[t]  = tipdot[mu].c01[t] ; 
            tipdot_old[mu].c02[t]  = tipdot[mu].c02[t] ;
            tipdot_old[mu].c12[t]  = tipdot[mu].c12[t] ;

        }
    }  

}


#endif
