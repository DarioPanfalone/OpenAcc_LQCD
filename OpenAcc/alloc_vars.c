#ifndef ALLOC_DEF_
#define ALLOC_DEF_

#include "struct_c_def.c"

double_soa * u1_back_field_phases;
tamat_soa * ipdot_acc;

tamat_soa * aux_ta; // aggiunta per il calcolo della forza stoutata
thmat_soa * aux_th; // aggiunta per il calcolo della forza stoutata

su3_soa  * conf_acc_bkp; // the old stored conf that will be recovered if the metro test fails.
su3_soa  * aux_conf_acc; // auxiliary 
su3_soa  * auxbis_conf_acc; // auxiliary 


// BACKRGROUND EM FIELD
double_soa * u1_back_field_phases; 

// GAUGE FIELD EVOLUTION
thmat_soa * momenta;
tamat_soa * ipdot_acc;



su3_soa * gconf_as_fermionmatrix; // conf to use in either cases 
                                  //in fermion related
                                  // computation (with or without stouting)
// STOUTING 
#ifdef STOUT_FERMIONS
su3_soa  * gstout_conf_acc; // max stouted conf,
                           // just pointer
su3_soa  * gstout_conf_acc_arr; // all stouting steps 
                               // except the zeroth
su3_soa * glocal_staples;
tamat_soa * gtipdot;
#endif

// FERMIONS

vec3_soa * ferm_chi_acc; // questo e' il chi [NPS_tot]
vec3_soa * ferm_phi_acc; // questo e' il phi [NPS_tot]
vec3_soa * ferm_out_acc; // questo e' uno ausiliario [NPS_tot]
vec3_soa * ferm_shiftmulti_acc; // ausiliario per l'invertitore multishift [max_ps*max_approx_order]
vec3_soa * kloc_r;  // vettore ausiliario
vec3_soa * kloc_h;  // vettore ausiliario
vec3_soa * kloc_s;  // vettore ausiliario
vec3_soa * kloc_p;  // vettore ausiliario
vec3_soa * k_p_shiftferm; // ausiliario [max_nshift=max_approx_order]


// LOCAL SUMS
dcomplex_soa * local_sums;
double_soa * d_local_sums;

void mem_alloc(){
  printf("Allocating resources for NPS_tot=%d pseudofermions in total, with max_approx_order=%d\n", NPS_tot, max_approx_order);
  int allocation_check;  
#ifdef BACKFIELD
  allocation_check =  posix_memalign((void **)&u1_back_field_phases, ALIGN, 8*sizeof(double_soa));   //  -->  4*size phases (as many as links)
  if(allocation_check != 0)  printf("Errore nella allocazione di u1_back_field_phases \n");
#else
  u1_back_field_phases=NULL;
#endif

  allocation_check =  posix_memalign((void **)&aux_conf_acc, ALIGN, 8*sizeof(su3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di aux_conf_acc \n");
  allocation_check =  posix_memalign((void **)&auxbis_conf_acc, ALIGN, 8*sizeof(su3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di auxbis_conf_acc \n");
  allocation_check =  posix_memalign((void **)&conf_acc_bkp, ALIGN, 8*sizeof(su3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di aux_conf_bkp \n");


  // GAUGE EVOLUTION
  allocation_check =  posix_memalign((void **)&momenta, ALIGN, 8*sizeof(thmat_soa));   //  -->  4*size
  if(allocation_check != 0)  printf("Errore nella allocazione di momenta \n");
  allocation_check =  posix_memalign((void **)&ipdot_acc, ALIGN, 8*sizeof(tamat_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di ipdot_acc \n");


#ifdef STOUT_FERMIONS
  // STOUTING
  allocation_check =  posix_memalign((void **)&gstout_conf_acc_arr, ALIGN, STOUT_STEPS*8*sizeof(su3_soa));
  gstout_conf_acc = &gstout_conf_acc_arr[8*(STOUT_STEPS-1)];
  if(allocation_check != 0)  printf("Errore nella allocazione di stout_conf_acc_arr \n");
  allocation_check =  posix_memalign((void **)&glocal_staples, ALIGN, 8*sizeof(su3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di glocal_staples \n");
  allocation_check =  posix_memalign((void **)&gtipdot, ALIGN, 8*sizeof(tamat_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di gtipdot \n");

  allocation_check =  posix_memalign((void **)&aux_th, ALIGN, 8*sizeof(thmat_soa));   //  -->  4*size
  if(allocation_check != 0)  printf("Errore nella allocazione di aux_th \n");
  allocation_check =  posix_memalign((void **)&aux_ta, ALIGN, 8*sizeof(tamat_soa));   //  -->  4*size
  if(allocation_check != 0)  printf("Errore nella allocazione di aux_ta \n");



#endif


  // FERMION ALLOCATIONS
  allocation_check =  posix_memalign((void **)&kloc_r, ALIGN, sizeof(vec3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di kloc_r \n");
  allocation_check =  posix_memalign((void **)&kloc_h, ALIGN, sizeof(vec3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di kloc_h \n");
  allocation_check =  posix_memalign((void **)&kloc_s, ALIGN, sizeof(vec3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di kloc_s \n");
  allocation_check =  posix_memalign((void **)&kloc_p, ALIGN, sizeof(vec3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di kloc_p \n");
  allocation_check =  posix_memalign((void **)&k_p_shiftferm, ALIGN, max_approx_order* sizeof(vec3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di k_p_shiftferm \n");
  allocation_check =  posix_memalign((void **)&ferm_chi_acc  , ALIGN, NPS_tot * sizeof(vec3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di ferm_chi_acc \n");
  allocation_check =  posix_memalign((void **)&ferm_phi_acc  , ALIGN, NPS_tot * sizeof(vec3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di ferm_phi_acc \n");
   allocation_check =  posix_memalign((void **)&ferm_out_acc  , ALIGN, NPS_tot * sizeof(vec3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di ferm_out_acc \n");
  allocation_check =  posix_memalign((void **)&ferm_shiftmulti_acc, ALIGN, max_ps*max_approx_order*sizeof(vec3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di ferm_shiftmulti_acc \n");


  allocation_check =  posix_memalign((void **)&d_local_sums, ALIGN, sizeof(double_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di d_local_sums \n");
  allocation_check =  posix_memalign((void **)&local_sums, ALIGN, 2*sizeof(dcomplex_soa));  // --> size complessi --> vettore per sommare cose locali
  if(allocation_check != 0)  printf("Errore nella allocazione di local_sums \n");
}

void mem_free(){
#ifdef BACKFIELD
  free(u1_back_field_phases);
#endif
  free(momenta);
  free(aux_conf_acc);
  free(aux_ta);
  free(aux_th);
  free(stout_conf_acc);
  free(conf_acc_bkp);
  free(ipdot_acc);

  free(ferm_chi_acc);
  free(ferm_phi_acc);
  free(ferm_out_acc);

  free(ferm_shiftmulti_acc);

  free(kloc_r);
  free(kloc_s);
  free(kloc_h);
  free(kloc_p);
  free(k_p_shiftferm);

  free(local_sums);
  free(d_local_sums);
}

#endif
