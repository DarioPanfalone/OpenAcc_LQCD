#ifndef FERMION_FORCE_UTILITIES_C
#define FERMION_FORCE_UTILITIES_C

#include "./geometry.h"
#include "./fermion_force_utilities.h"
#include "./struct_c_def.h"
#include "../Include/fermion_parameters.h"
#include "./fermionic_utilities.h"
#include "./fermion_matrix.h"

void set_tamat_soa_to_zero( __restrict tamat_soa * const matrix)
{
  int hx, y, z, t;
  int mu;
  SETINUSE(matrix);
#pragma acc kernels present(matrix)
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
	for(hx=0; hx < nxh; hx++) {
	  int x,idxh;
	  x = 2*hx + ((y+z+t) & 0x1);
	  idxh = snum_acc(x,y,z,t);
	  for(mu=0; mu<8; mu++) {
	    assign_zero_to_tamat_soa_component(&matrix[mu],idxh);
	  }
	}  // x
      }  // y
    }  // z
  }  // t
}






void direct_product_of_fermions_into_auxmat(__restrict vec3_soa  * const loc_s, // questo fermione e' costante e non viene modificato qui dentro
					    __restrict vec3_soa  * const loc_h, // questo fermione e' costante e non viene modificato qui dentro
					    __restrict su3_soa * const aux_u,
					    const RationalApprox * const approx,
					    int iter){
    SETINUSE(aux_u);
  
  //   ////////////////////////////////////////////////////////////////////////   //
  //    Riflettere se conviene tenere i loop sui siti pari e su quelli dispari    //
  //    separati come sono adesso o se invece conviene fare un unico kernel!!!    //
  //   ////////////////////////////////////////////////////////////////////////   //

  //LOOP SUI SITI PARI
  int xh, y, z, t;
#pragma acc kernels present(loc_s) present(loc_h)  present(approx) present(aux_u)
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
        for(xh=0; xh < nxh; xh++) {
          int idxh,idxpmu,x;
          int parity;
	  int dir_mu;
	  int mu;
	  x = 2*xh + ((y+z+t) & 0x1);
          idxh = snum_acc(x,y,z,t);  // r
	  //  parity = (x+y+z+t) % 2;
	  parity = 0; // la fisso cosi' perche' sto prendendo il sito pari

	  for(mu=0;mu<4;mu++){
	    idxpmu = nnp_openacc[idxh][mu][parity];// r+mu        
	    dir_mu = 2*mu +  parity;
	    vec1_directprod_conj_vec2_into_mat1(&aux_u[dir_mu],idxh,loc_h,idxpmu,loc_s,idxh,approx->RA_a[iter]);
	  }//mu

        }  // x     
      }  // y       
    }  // z         
  }  // t

  //LOOP SUI SITI DISPARI
#pragma acc kernels present(loc_s) present(loc_h)  present(approx) present(aux_u)
#pragma acc loop independent gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent vector(DIM_BLOCK_X)
        for(xh=0; xh < nxh; xh++) {
          int idxh,idxpmu,x;
          int parity;
	  int dir_mu;
	  int mu;
	  x = 2*xh + ((y+z+t+1) & 0x1);
          idxh = snum_acc(x,y,z,t);  // r
	  //  parity = (x+y+z+t) % 2;
	  parity = 1; // la fisso cosi' perche' sto prendendo il sito dispari
#pragma acc loop independent
	  for(mu=0;mu<4;mu++){
	    idxpmu = nnp_openacc[idxh][mu][parity];// r+mu        
	    dir_mu = 2*mu +  parity;
	    vec1_directprod_conj_vec2_into_mat1(&aux_u[dir_mu],idxh,loc_s,idxpmu,loc_h,idxh,-approx->RA_a[iter]);
	  }//mu
        }  // x     
      }  // y       
    }  // z         
  }  // t
}// closes routine     



void multiply_conf_times_force_and_take_ta_even(__restrict su3_soa * const u, // la conf e' costante e non viene modificata
						__restrict ferm_param * const tpars,
						__restrict double_soa * const backfield,
						__restrict su3_soa * const auxmat, // anche questa conf ausiliaria e' costante e non viene modificata
						__restrict tamat_soa * const ipdot){

  SETINUSE(ipdot);
  int hx,y,z,t,idxh;
#ifdef BACKFIELD
#pragma acc kernels present(u) present(auxmat) present(ipdot) present(tpars) present(backfield)
#else
#pragma acc kernels present(u) present(auxmat) present(ipdot) present(tpars)
#endif
#pragma acc loop independent //gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent //vector(DIM_BLOCK_X)
        for(hx=0; hx < nxh; hx++) {
          int x,eta;

          double arg;
          d_complex phase;
#ifdef IMCHEMPOT
          double imchempot = tpars->ferm_im_chem_pot/((double)(nt));
#endif

#ifdef BACKFIELD
          double charge = (double)(tpars->ferm_charge);
#endif

          //even sites
          x = 2*hx + ((y+z+t) & 0x1);
          idxh = snum_acc(x,y,z,t);

          // dir  0  =  x even   --> eta = 1 , no multiplication needed
	  eta = 1;
#ifdef BACKFIELD
	  arg = backfield[0].d[idxh] * charge;
          phase = cos(arg) + I * sin(arg);
	  mat1_times_auxmat_into_tamat(&u[0],idxh,eta,&auxmat[0],idxh,&ipdot[0],idxh,phase);
#else
	  mat1_times_auxmat_into_tamat(&u[0],idxh,eta,&auxmat[0],idxh,&ipdot[0],idxh);
#endif


          // dir  2  =  y even
          eta = 1 - ( 2*(x & 0x1) );
#ifdef BACKFIELD
	  arg = backfield[2].d[idxh] * charge;
          phase = cos(arg) + I * sin(arg);
	  mat1_times_auxmat_into_tamat(&u[2],idxh,eta,&auxmat[2],idxh,&ipdot[2],idxh,phase);
#else
	  mat1_times_auxmat_into_tamat(&u[2],idxh,eta,&auxmat[2],idxh,&ipdot[2],idxh);
#endif


          // dir  4  =  z even
          eta = 1 - ( 2*((x+y) & 0x1) );
#ifdef BACKFIELD
	  arg = backfield[4].d[idxh] * charge;
          phase = cos(arg) + I * sin(arg);
	  mat1_times_auxmat_into_tamat(&u[4],idxh,eta,&auxmat[4],idxh,&ipdot[4],idxh,phase);
#else
	  mat1_times_auxmat_into_tamat(&u[4],idxh,eta,&auxmat[4],idxh,&ipdot[4],idxh);
#endif

          // dir  6  =  t even
          eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
          eta *= (1- 2*(int)(t/(nt-1)));
#endif
	  arg = 0;
#ifdef BACKFIELD
	  arg += backfield[6].d[idxh] * charge;
#endif
#ifdef IMCHEMPOT
          arg += imchempot;
#endif
#ifdef PHASE_MAT_VEC_MULT
          phase = cos(arg) + I * sin(arg);
	  mat1_times_auxmat_into_tamat(&u[6],idxh,eta,&auxmat[6],idxh,&ipdot[6],idxh,phase);
#else
	  mat1_times_auxmat_into_tamat(&u[6],idxh,eta,&auxmat[6],idxh,&ipdot[6],idxh);
#endif

        }
      }
    }
  }

}



void multiply_conf_times_force_and_take_ta_odd(  __restrict su3_soa * const u, // e' costante e non viene modificata
						 __restrict ferm_param * const tpars,
						 __restrict double_soa * const backfield,
					         __restrict su3_soa * const auxmat, // e' costante e non viene modificata
					         __restrict tamat_soa * const ipdot){ 
    SETINUSE(ipdot);
  int hx,y,z,t,idxh;
#ifdef BACKFIELD
#pragma acc kernels present(u) present(auxmat) present(ipdot) present(tpars) present(backfield)
#else
#pragma acc kernels present(u) present(auxmat) present(ipdot) present(tpars)
#endif
#pragma acc loop independent //gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent //vector(DIM_BLOCK_X)
        for(hx=0; hx < nxh; hx++) {
          int x,eta;
	  double arg;
	  d_complex phase;
#ifdef IMCHEMPOT
          double imchempot = tpars->ferm_im_chem_pot/((double)(nt));
#endif
#ifdef BACKFIELD
          double charge = (double)(tpars->ferm_charge);
#endif

          //odd sites
          x = 2*hx + ((y+z+t+1) & 0x1);
          idxh = snum_acc(x,y,z,t);
          // dir  1  =  x odd    --> eta = 1 , no multiplication needed
	  eta = 1;
#ifdef BACKFIELD
          arg = backfield[1].d[idxh] * charge;
          phase = cos(arg) + I * sin(arg);
          mat1_times_auxmat_into_tamat(&u[1],idxh,eta,&auxmat[1],idxh,&ipdot[1],idxh,phase);
#else
          mat1_times_auxmat_into_tamat(&u[1],idxh,eta,&auxmat[1],idxh,&ipdot[1],idxh);
#endif

          // dir  3  =  y odd
          eta = 1 - ( 2*(x & 0x1) );
#ifdef BACKFIELD
          arg = backfield[3].d[idxh] * charge;
          phase = cos(arg) + I * sin(arg);
          mat1_times_auxmat_into_tamat(&u[3],idxh,eta,&auxmat[3],idxh,&ipdot[3],idxh,phase);
#else
          mat1_times_auxmat_into_tamat(&u[3],idxh,eta,&auxmat[3],idxh,&ipdot[3],idxh);
#endif

          // dir  5  =  z odd
	  eta = 1 - ( 2*((x+y) & 0x1) );
#ifdef BACKFIELD
          arg = backfield[5].d[idxh] * charge;
          phase = cos(arg) + I * sin(arg);
          mat1_times_auxmat_into_tamat(&u[5],idxh,eta,&auxmat[5],idxh,&ipdot[5],idxh,phase);
#else
	  mat1_times_auxmat_into_tamat(&u[5],idxh,eta,&auxmat[5],idxh,&ipdot[5],idxh);
#endif

          // dir  7  =  t odd
          eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
          eta *= (1- 2*(int)(t/(nt-1)));
#endif
          arg = 0;
#ifdef BACKFIELD
          arg += backfield[7].d[idxh] * charge;
#endif
#ifdef IMCHEMPOT
          arg += imchempot;
#endif
#ifdef PHASE_MAT_VEC_MULT
          phase = cos(arg) + I * sin(arg);
          mat1_times_auxmat_into_tamat(&u[7],idxh,eta,&auxmat[7],idxh,&ipdot[7],idxh,phase);
#else
          mat1_times_auxmat_into_tamat(&u[7],idxh,eta,&auxmat[7],idxh,&ipdot[7],idxh);
#endif


        } //hx
      } //y
    } // z
  } // t
} // end  multiply_conf_times_force_and_take_ta_odd()


#ifdef STOUT_FERMIONS
void multiply_conf_times_force_and_take_ta_even_nophase(__restrict su3_soa * const u, // la conf e' costante e non viene modificata
							__restrict su3_soa * const auxmat, // anche questa conf ausiliaria e' costante; non viene modificata
							__restrict tamat_soa * const ipdot){
SETINUSE(ipdot);
  int hx,y,z,t,idxh;
#pragma acc kernels present(u) present(auxmat) present(ipdot) 
#pragma acc loop independent //gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent //vector(DIM_BLOCK_X)
          for(hx=0; hx < nxh; hx++) {
              int x,eta;

              //even sites
              x = 2*hx + ((y+z+t) & 0x1);
              idxh = snum_acc(x,y,z,t);

              // dir  0  =  x even   --> eta = 1 , no multiplication needed
              eta = 1;
              mat1_times_auxmat_into_tamat_nophase(&u[0],idxh,eta,&auxmat[0],idxh,&ipdot[0],idxh); 
              // dir  2  =  y even
              eta = 1 - ( 2*(x & 0x1) );
              mat1_times_auxmat_into_tamat_nophase(&u[2],idxh,eta,&auxmat[2],idxh,&ipdot[2],idxh);
              // dir  4  =  z even
              eta = 1 - ( 2*((x+y) & 0x1) );
              mat1_times_auxmat_into_tamat_nophase(&u[4],idxh,eta,&auxmat[4],idxh,&ipdot[4],idxh);

              // dir  6  =  t even
              eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
              eta *= (1- 2*(int)(t/(nt-1)));
#endif
              mat1_times_auxmat_into_tamat_nophase(&u[6],idxh,eta,&auxmat[6],idxh,&ipdot[6],idxh);

          } // hx
      } // y
    } // z
  } // t
}


void multiply_conf_times_force_and_take_ta_odd_nophase(  __restrict su3_soa * const u, // e' costante e non viene modificata
							 __restrict su3_soa * const auxmat, // e' costante e non viene modificata
							 __restrict tamat_soa * const ipdot){
    SETINUSE(ipdot);
  int hx,y,z,t,idxh;
#pragma acc kernels present(u) present(auxmat) present(ipdot)
#pragma acc loop independent //gang(nt)
  for(t=0; t<nt; t++) {
#pragma acc loop independent //gang(nz/DIM_BLOCK_Z) vector(DIM_BLOCK_Z)
    for(z=0; z<nz; z++) {
#pragma acc loop independent //gang(ny/DIM_BLOCK_Y) vector(DIM_BLOCK_Y)
      for(y=0; y<ny; y++) {
#pragma acc loop independent //vector(DIM_BLOCK_X)
        for(hx=0; hx < nxh; hx++) {
          int x,eta;
          //odd sites
          x = 2*hx + ((y+z+t+1) & 0x1);
          idxh = snum_acc(x,y,z,t);
          // dir  1  =  x odd    --> eta = 1 , no multiplication needed
          eta = 1;
          mat1_times_auxmat_into_tamat_nophase(&u[1],idxh,eta,&auxmat[1],idxh,&ipdot[1],idxh);

          // dir  3  =  y odd
          eta = 1 - ( 2*(x & 0x1) );
          mat1_times_auxmat_into_tamat_nophase(&u[3],idxh,eta,&auxmat[3],idxh,&ipdot[3],idxh);

          // dir  5  =  z odd
	  eta = 1 - ( 2*((x+y) & 0x1) );
	  mat1_times_auxmat_into_tamat_nophase(&u[5],idxh,eta,&auxmat[5],idxh,&ipdot[5],idxh);

          // dir  7  =  t odd
          eta = 1 - ( 2*((x+y+z) & 0x1) );
#ifdef ANTIPERIODIC_T_BC
          eta *= (1- 2*(int)(t/(nt-1)));
#endif
          mat1_times_auxmat_into_tamat_nophase(&u[7],idxh,eta,&auxmat[7],idxh,&ipdot[7],idxh);

        } //hx
      } //y
    } // z
  } // t
} // end  multiply_conf_times_force_and_take_ta_odd()
#endif //ifdef STOUT_FERMIONS


#if defined(IMCHEMPOT) || defined(BACKFIELD)
void multiply_backfield_times_force(__restrict ferm_param * const tpars,
        __restrict double_soa * const backfield,
        __restrict su3_soa * const auxmat, // anche questa conf ausiliaria e' costante e non viene modificata
        __restrict su3_soa * const pseudo_ipdot){
    SETINUSE(pseudo_ipdot);

  double arg;
  d_complex phase;
#ifdef IMCHEMPOT
  double imchempot = tpars->ferm_im_chem_pot/((double)(nt));
#endif
  double charge = (double)(tpars->ferm_charge);
  int idxh;
#pragma acc data present(backfield) present(auxmat) present(pseudo_ipdot)
#pragma acc loop independent
  for(int dirindex = 0 ; dirindex < 8 ; dirindex++){
#pragma acc loop independent
    for( idxh = 0 ; idxh < sizeh; idxh++){
      
#ifdef BACKFIELD
      arg = backfield[dirindex].d[idxh] * charge;
#else
      arg = 0.0;
#endif

#ifdef IMCHEMPOT
      if(dirindex>6) arg += imchempot;
#endif
      phase = cos(arg) + I * sin(arg);
      phase_times_auxmat_into_auxmat(&auxmat[dirindex],&pseudo_ipdot[dirindex],idxh,phase);
    }
  }
} // end multiply_backfield_times_force()

#else 
void accumulate_gl3soa_into_gl3soa(
        __restrict su3_soa * const auxmat, // anche questa conf ausiliaria e' costante e non viene modificata
        __restrict su3_soa * const pseudo_ipdot){

    SETINUSE(pseudo_ipdot);

    int idxh, dirindex;
//
#pragma acc kernels present(auxmat) present(pseudo_ipdot)
#pragma acc loop independent
    for(dirindex = 0 ; dirindex < 8 ; dirindex++){
#pragma acc loop independent
        for( idxh = 0 ; idxh < sizeh; idxh++){
            accumulate_auxmat1_into_auxmat2(&auxmat[dirindex],&pseudo_ipdot[dirindex],idxh);
        }
    }
}

#endif

void ker_openacc_compute_fermion_force( __restrict su3_soa * const u, // e' costante e non viene mai modificato qui dentro
					double_soa * backfield,
					__restrict su3_soa * const aux_u,
					__restrict vec3_soa * const in_shiftmulti,  // e' costante e non viene mai modificato qui dentro
					__restrict vec3_soa  * const loc_s,
					__restrict vec3_soa  * const loc_h,
					ferm_param  *  tpars
					){
  int ih;
  int iter=0;

  for(iter=0; iter<tpars->approx_md.approx_order; iter++){
    assign_in_to_out(&in_shiftmulti[iter],loc_s);
    acc_Doe(u,loc_h,loc_s,tpars,backfield);
    direct_product_of_fermions_into_auxmat(loc_s,loc_h,aux_u,&(tpars->approx_md),iter);
  }
}



#endif
