
#ifndef FERMIONMATRIX_CC_
#define FERMIONMATRIX_CC_


#include "../Fermions/fermions.cc"
#include "../Include/global_var.cc"
#include "../Include/parameters.cc"

// Hopping term for links starting at EVEN sites and 
// ending at ODD sites
//
// (vector_even)_out=Deo*(vector_odd)_in
void Deo(Fermion *out, const Fermion *in)
 {
 long int even;
 Vec3 aux;
 Vec3 aux1;
 for(even=0; even<sizeh; even++)
    {
    // nnp[even,i]=odd  ---> sizeh<=nnp[even,i]<size   but fermion[j] 
    // defined for 0<=j<sizeh 
 
    aux =(gauge_conf->u_work[even      ])*(in->fermion[nnp[even][0]-sizeh]);
    aux+=(gauge_conf->u_work[even+size ])*(in->fermion[nnp[even][1]-sizeh]);
    aux+=(gauge_conf->u_work[even+size2])*(in->fermion[nnp[even][2]-sizeh]);
    
    aux+=(gauge_conf->u_work[even+size3])*(in->fermion[nnp[even][3]-sizeh]);
 

    aux-=(~(gauge_conf->u_work[nnm[even][0]      ]))*(in->fermion[nnm[even][0]-sizeh]);
    aux-=(~(gauge_conf->u_work[nnm[even][1]+size ]))*(in->fermion[nnm[even][1]-sizeh]);
    aux-=(~(gauge_conf->u_work[nnm[even][2]+size2]))*(in->fermion[nnm[even][2]-sizeh]);
    aux-=(~(gauge_conf->u_work[nnm[even][3]+size3]))*(in->fermion[nnm[even][3]-sizeh]);

    aux*=0.5;

    (out->fermion[even])=aux;
    } 
 }



// Hopping term for links starting at ODD sites and 
// ending at EVEN sites
//
// (vector_odd)_out=Deo*(vector_even)_in
void Doe(Fermion *out, const Fermion *in)
 {
 long int odd;
 Vec3 aux;
 Vec3 aux1;

 for(odd=sizeh; odd<size; odd++)
    {
    // nnp[odd,i]=even  ---> 0<=nnp[odd,i]<sizeh 
    aux =(gauge_conf->u_work[odd]      )*(in->fermion[nnp[odd][0]]);
    aux+=(gauge_conf->u_work[odd+size] )*(in->fermion[nnp[odd][1]]);
    aux+=(gauge_conf->u_work[odd+size2])*(in->fermion[nnp[odd][2]]);
    aux+=(gauge_conf->u_work[odd+size3])*(in->fermion[nnp[odd][3]]);

    aux-=(~(gauge_conf->u_work[nnm[odd][0]      ]))*(in->fermion[nnm[odd][0]]);
    aux-=(~(gauge_conf->u_work[nnm[odd][1]+size ]))*(in->fermion[nnm[odd][1]]);
    aux-=(~(gauge_conf->u_work[nnm[odd][2]+size2]))*(in->fermion[nnm[odd][2]]);
    
    aux-=(~(gauge_conf->u_work[nnm[odd][3]+size3]))*(in->fermion[nnm[odd][3]]);
    
    aux*=0.5;

    (out->fermion[odd-sizeh])=aux;
    } 
 }


#endif
