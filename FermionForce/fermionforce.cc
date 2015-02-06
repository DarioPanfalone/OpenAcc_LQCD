#ifdef USE_GPU
#include"../Cuda/cuda_fermion_force.h"
#include"../Inverter/cu_inverter.h"
#include"../Packer/packer.h"
#endif

void fermionforce(int moltiplico)
  {
  #ifdef DEBUG_MODE
  cout << "DEBUG: inside fermionforce ..."<<endl;
  #endif

 #ifdef IM_CHEM_POT
  const complex<REAL> eim=complex<REAL>(eim_cos, eim_sin);
  const complex<REAL> emim=complex<REAL>(eim_cos, -eim_sin);
 #endif

  int pseudofermion, iter, mu;
  long int even, odd, x, y;
  RationalApprox approx;
  Vec3 vr_1, vr_2;
  Su3 aux, aux1;

  // approximation of x^{-beta}    beta=no_flavours/4/no_ps
  approx.md_inv_approx_coeff();

  // shift inverter
  #ifdef USE_GPU
    int order;
    get_order(order, approx);

    cu_multips_shifted_invert(residue_md, approx);
  #else
    multips_shifted_invert(fermion_shiftmulti, fermion_chi, residue_md, approx);   //////////////////////////////     che succede se B =/= 0 ??
  #endif

  // force reconstruction
  #ifdef USE_GPU 
    float num[max_approx_order];
    get_numerators(num, approx);
    cuda_fermion_force(order, num);  // here fermion force is also multiplied by u_work and taken TA

  #else
    for(pseudofermion=0; pseudofermion<no_ps; pseudofermion++)
      {
	for(iter=0; iter<(approx.approx_order); iter++)
          {
	    for(even=0; even<sizeh; even++)
	      {
		(loc_s->fermion[even])=(fermion_shiftmulti->fermion[pseudofermion][iter][even]);
	      }

#ifdef MAGN	    
	    Doe(loc_h, loc_s, pseudofermion);
#else
	    Doe(loc_h, loc_s);
#endif

	    for(even=0; even<sizeh; even++)
	      {
		odd=even+sizeh;
		for(mu=0; mu<3; mu++)
		  {
		    x=even;
		    y=nnp[even][mu]-sizeh;     // sizeh<=nnp[even][mu]<size
		    vr_1=(loc_h->fermion[y]);
		    vr_2=~(loc_s->fermion[x]);
		    vr_1*=(approx.RA_a[iter]);
		    aux=(vr_1^vr_2);
		    if(moltiplico==0){
		      (gauge_ipdot->ipdot[mu*size+even])+=aux;
		    }
		    if(moltiplico==1){
		      aux1=(gauge_conf->u_work[mu*size+even])*aux;
#ifdef MAGN
		      aux1*=b_field[pseudofermion][mu*size+even];
#endif 
		      (gauge_ipdot->ipdot[mu*size+even])+=aux1;
		    }
		    
		    y=nnp[odd][mu];
		    vr_1=(loc_s->fermion[y]);
		    vr_2=~(loc_h->fermion[x]);
		    vr_1*=(-approx.RA_a[iter]);
		    aux=(vr_1^vr_2);
		    if(moltiplico==0){
		      (gauge_ipdot->ipdot[mu*size+odd])+=aux;
		    }
		    if(moltiplico==1){
		      aux1=(gauge_conf->u_work[mu*size+odd])*aux;
#ifdef MAGN
		      aux1*=b_field[pseudofermion][mu*size+odd];
#endif
		      (gauge_ipdot->ipdot[mu*size+odd])+=aux1;
		    }



		  }
		for(mu=3; mu<4; mu++)
		  {
		    x=even;
		    y=nnp[even][mu]-sizeh;     // sizeh<=nnp[even][mu]<size
		    vr_1=(loc_h->fermion[y]);
		    vr_2=~(loc_s->fermion[x]);
#ifdef IM_CHEM_POT
		    vr_1*=(approx.RA_a[iter]*eim);
#else
		    vr_1*=(approx.RA_a[iter]);
#endif
		    aux=(vr_1^vr_2);
		    if(moltiplico==0){
		      (gauge_ipdot->ipdot[mu*size+even])+=aux;
		    }
		    if(moltiplico==1){
		      aux1=(gauge_conf->u_work[mu*size+even])*aux;
#ifdef MAGN
		      aux1*=b_field[pseudofermion][mu*size+even];
#endif
		      (gauge_ipdot->ipdot[mu*size+even])+=aux1;
		    }

		    y=nnp[odd][mu];
		    vr_1=(loc_s->fermion[y]);
		    vr_2=~(loc_h->fermion[x]);
#ifdef IM_CHEM_POT
		    vr_1*=(-approx.RA_a[iter]*eim);
#else
		    vr_1*=(-approx.RA_a[iter]);
#endif
		    aux=(vr_1^vr_2);
		    if(moltiplico==0){
		      (gauge_ipdot->ipdot[mu*size+odd])+=aux;
		    }
		    if(moltiplico==1){
		      aux1=(gauge_conf->u_work[mu*size+odd])*aux;
#ifdef MAGN
		      aux1*=b_field[pseudofermion][mu*size+odd];
#endif
		      (gauge_ipdot->ipdot[mu*size+odd])+=aux1;
		    }
		  }
	      } 
          }
      }
#endif

  #ifdef DEBUG_MODE
  cout << "\tterminated fermionforce"<<endl;
  #endif
 }
