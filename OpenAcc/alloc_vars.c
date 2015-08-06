#ifndef ALLOC_DEF_
#define ALLOC_DEF_

double_soa * u1_back_field_phases;
tamat_soa * ipdot_acc;
su3_soa  * conf_acc_bkp; // the old stored conf that will be recovered if the metro test fails.
su3_soa  * aux_conf_acc; // auxiliary 
ACC_MultiFermion * ferm_chi_acc; // questo e' il chi
ACC_MultiFermion * ferm_phi_acc; // questo e' il phi
ACC_MultiFermion * ferm_out_acc; // questo e' uno ausiliario
ACC_ShiftMultiFermion * ferm_shiftmulti_acc; // ausiliario per l'invertitore multishift
vec3_soa * kloc_r;  // vettore ausiliario
vec3_soa * kloc_h;  // vettore ausiliario
vec3_soa * kloc_s;  // vettore ausiliario
vec3_soa * kloc_p;  // vettore ausiliario
ACC_ShiftFermion *k_p_shiftferm;
thmat_soa * momenta;
dcomplex_soa * local_sums;
double_soa * d_local_sums;
// RATIONAL APPROX INITIALIZATION ///////////////////////////////////////////////////////////////                                                          
// la 1 e' la "first_inv_approx_norm_coeff"                                                                                                                
// la 2 e' la "md_inv_approx_norm_coeff"                                                                                                                   
// la 3 e' la "last_inv_approx_norm_coeff"                                                                                                                 
COM_RationalApprox *approx1;
COM_RationalApprox *approx2;
COM_RationalApprox *approx3;
COM_RationalApprox *approx_mother1;
COM_RationalApprox *approx_mother2;
COM_RationalApprox *approx_mother3;
//////////////////////////////////////////////////////////////////////////////////////////////////
double minmaxeig[2];




void mem_alloc(){
  int allocation_check;  
#ifdef BACKFIELD
  allocation_check =  posix_memalign((void **)&u1_back_field_phases, ALIGN, 8*sizeof(double_soa));   //  -->  4*size phases (as many as links)
  if(allocation_check != 0)  printf("Errore nella allocazione di u1_back_field_phases \n");
#else
  u1_back_field_phases=NULL;
#endif
  allocation_check =  posix_memalign((void **)&momenta, ALIGN, 8*sizeof(thmat_soa));   //  -->  4*size
  if(allocation_check != 0)  printf("Errore nella allocazione di momenta \n");
  allocation_check =  posix_memalign((void **)&kloc_r, ALIGN, sizeof(vec3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di kloc_r \n");
  allocation_check =  posix_memalign((void **)&kloc_h, ALIGN, sizeof(vec3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di kloc_h \n");
  allocation_check =  posix_memalign((void **)&kloc_s, ALIGN, sizeof(vec3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di kloc_s \n");
  allocation_check =  posix_memalign((void **)&kloc_p, ALIGN, sizeof(vec3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di kloc_p \n");
  allocation_check =  posix_memalign((void **)&k_p_shiftferm, ALIGN, sizeof(ACC_ShiftFermion));
  if(allocation_check != 0)  printf("Errore nella allocazione di k_p_shiftferm \n");
  allocation_check =  posix_memalign((void **)&aux_conf_acc, ALIGN, 8*sizeof(su3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di aux_conf_acc \n");
  allocation_check =  posix_memalign((void **)&conf_acc_bkp, ALIGN, 8*sizeof(su3_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di aux_conf_bkp \n");
  allocation_check =  posix_memalign((void **)&ipdot_acc, ALIGN, 8*sizeof(tamat_soa));
  if(allocation_check != 0)  printf("Errore nella allocazione di ipdot_acc \n");
  allocation_check =  posix_memalign((void **)&ferm_chi_acc  , ALIGN, sizeof(ACC_MultiFermion));
  if(allocation_check != 0)  printf("Errore nella allocazione di ferm_chi_acc \n");
  allocation_check =  posix_memalign((void **)&ferm_phi_acc  , ALIGN, sizeof(ACC_MultiFermion));
  if(allocation_check != 0)  printf("Errore nella allocazione di ferm_phi_acc \n");
  allocation_check =  posix_memalign((void **)&ferm_out_acc , ALIGN, sizeof(ACC_MultiFermion));
  if(allocation_check != 0)  printf("Errore nella allocazione di ferm_out_acc \n");
  allocation_check =  posix_memalign((void **)&ferm_shiftmulti_acc, ALIGN, sizeof(ACC_ShiftMultiFermion));
  if(allocation_check != 0)  printf("Errore nella allocazione di ferm_shiftmulti_acc \n");
  allocation_check =  posix_memalign((void **)&approx1, ALIGN, sizeof(COM_RationalApprox));
  if(allocation_check != 0)  printf("Errore nella allocazione di approx1 \n");
  allocation_check =  posix_memalign((void **)&approx2, ALIGN, sizeof(COM_RationalApprox));
  if(allocation_check != 0)  printf("Errore nella allocazione di approx2 \n");
  allocation_check =  posix_memalign((void **)&approx3, ALIGN, sizeof(COM_RationalApprox));
  if(allocation_check != 0)  printf("Errore nella allocazione di approx3 \n");
  allocation_check =  posix_memalign((void **)&approx_mother1, ALIGN, sizeof(COM_RationalApprox));
  if(allocation_check != 0)  printf("Errore nella allocazione di approx_mother1 \n");
  allocation_check =  posix_memalign((void **)&approx_mother2, ALIGN, sizeof(COM_RationalApprox));
  if(allocation_check != 0)  printf("Errore nella allocazione di approx_mother2 \n");
  allocation_check =  posix_memalign((void **)&approx_mother3, ALIGN, sizeof(COM_RationalApprox));
  if(allocation_check != 0)  printf("Errore nella allocazione di approx_mother3 \n");
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
  free(approx1);
  free(approx2);
  free(approx3);
  free(approx_mother1);
  free(approx_mother2);
  free(approx_mother3);
}

#endif
