#ifdef USE_GPU
#include"../Cuda/cuda_exponentiate.h"
#include"../Cuda/cuda_mom_sum_mul.h"
#include"../Cuda/cuda_init.h"
#include"../Packer/packer.h"
#endif

class Momenta {
  //private:
public:
 Su3 momenta[no_links];
  
  Momenta(void);

  friend void create_momenta(void);
  friend void conf_left_exp_multiply(const complex<REAL> p);
  friend void momenta_sum_multiply(const complex<REAL> p);
  
  void mom_aos_to_soaCOM(thmatCOM_soa *out,int matdir) const;
  void mom_soaCOM_to_aos(thmatCOM_soa const * const in,int matdir);

  // defined in Packer/packer.cc
  friend void smartpack_thmatrix(float out[8*no_links], Momenta *in);
  friend void smartunpack_thmatrix(Momenta *out, float in[8*no_links]);

  // defined in Update/action.cc
  friend void calc_action(double *value, int init);
  
#ifdef REVERSIBILITY_TEST
  friend void reverse_momenta(void);
  friend void save_momenta(void);
  friend void rev_update(void);
#endif
};

void Momenta::mom_aos_to_soaCOM( thmatCOM_soa *out,int matdir) const{
  int indice;
  for( int i =0 ; i < sizeh ; i++){
    indice = i + sizeh * matdir;
    out->c01[i].Re = (momenta[indice].comp[0][1]).real();
    out->c01[i].Im = (momenta[indice].comp[0][1]).imag();
    out->c02[i].Re = (momenta[indice].comp[0][2]).real();
    out->c02[i].Im = (momenta[indice].comp[0][2]).imag();
    out->c12[i].Re = (momenta[indice].comp[1][2]).real();
    out->c12[i].Im = (momenta[indice].comp[1][2]).imag();
    out->rc00[i] = (momenta[indice].comp[0][0]).real();
    out->rc11[i] = (momenta[indice].comp[1][1]).real();
  }
}


void Momenta::mom_soaCOM_to_aos( thmatCOM_soa const * const in,int matdir){
  int indice;
  for( int i =0 ; i < sizeh ; i++){
    indice = i + sizeh * matdir;
    // assign
    momenta[indice].comp[0][1] = complex<double>(in->c01[i].Re,in->c01[i].Im);
    momenta[indice].comp[0][2] = complex<double>(in->c02[i].Re,in->c02[i].Im);
    momenta[indice].comp[1][2] = complex<double>(in->c12[i].Re,in->c12[i].Im);
    momenta[indice].comp[0][0] = complex<double>(in->rc00[i],0.0);
    momenta[indice].comp[1][1] = complex<double>(in->rc11[i],0.0);
    //reconstruct the other elements
    momenta[indice].comp[1][0] = conj(momenta[indice].comp[0][1]);
    momenta[indice].comp[2][0] = conj(momenta[indice].comp[0][2]);
    momenta[indice].comp[2][1] = conj(momenta[indice].comp[1][2]);
    momenta[indice].comp[2][2] = -momenta[indice].comp[0][0] -momenta[indice].comp[1][1];
  }
}



// base constructor
Momenta::Momenta(void)
 {
 for(long int i=0; i<no_links; i++)
    {
    momenta[i].zero();
    }
 }


// create gaussian distributed momenta
void create_momenta(void)
  {
  #ifdef DEBUG_MODE
  cout << "DEBUG: inside create_momenta ..."<<endl;
  #endif

  for(long int i=0; i<no_links; i++)
     {
     (gauge_momenta->momenta[i]).gauss();
     }
  #ifdef DEBUG_MODE
  cout << "\tterminated create_momenta"<<endl;
  #endif
  } 


// u_work[i] -> exp(p*momenta[i])*u_work[i]
void conf_left_exp_multiply(const complex<REAL> p)
  {
  #ifdef DEBUG_MODE
  cout << "DEBUG: inside conf_left_exp_multiply ..."<<endl;
  #endif

  #ifdef TIMING_CUDA_CPP
  clock_t time_start, time_finish;
  time_start=clock();
  #endif

  #ifdef USE_GPU
    cuda_exp_su3(imag(p));
  #else
    Su3 aux;
    long int i;

    for(i=0; i<no_links; i++)
       {
       aux=(gauge_momenta->momenta[i]);
       aux*=p;
       aux.exp();
       aux*=(gauge_conf->u_work[i]);
       (gauge_conf->u_work[i])=aux;
       }
       gauge_conf->unitarize_with_eta(); 
  #endif

  #ifdef TIMING_CUDA_CPP
  time_finish=clock();
  cout << "time for conf_left_exp_multiply = " << ((REAL)(time_finish)-(REAL)(time_start))/CLOCKS_PER_SEC << " sec.\n";
  #endif

  #ifdef DEBUG_MODE
  cout << "\tterminated conf_left_exp_multiply"<<endl;
  #endif
  }


// momenta[i]+=p*ipdot_loc[i]
void momenta_sum_multiply(const complex<REAL> p)
  {
  #ifdef DEBUG_MODE
  cout << "DEBUG: inside momenta_sum_multiply ..."<<endl;
  #endif

  #ifdef TIMING_CUDA_CPP
  clock_t time_start, time_finish;
  time_start=clock();
  #endif

  #ifdef USE_GPU
    cuda_momenta_sum_multiply(imag(p));
  #else
    long int i;
    Su3 aux;

    for(i=0; i<no_links; i++)
       {
       aux=(gauge_ipdot->ipdot[i]);
       aux*=p;
       (gauge_momenta->momenta[i])+=aux;
       }
  #endif

  #ifdef TIMING_CUDA_CPP
  time_finish=clock();
  cout << "time for momenta_sum_multiply = " << ((REAL)(time_finish)-(REAL)(time_start))/CLOCKS_PER_SEC << " sec.\n";
  #endif

  #ifdef DEBUG_MODE
  cout << "\tterminated momenta_sum_multiply"<<endl;
  #endif
  }


#ifdef REVERSIBILITY_TEST
// reverse momenta
// to be used in reversibility tests
void reverse_momenta(void)
  {
  #ifndef USE_GPU
  for(long int i=0; i<no_links; i++)
     {
     (gauge_momenta->momenta[i])*=(-1.0);
     }
  #else
  cuda_get_momenta();
  smartunpack_thmatrix(gauge_momenta, momenta_packed);
  for(long int i=0; i<no_links; i++)
     {
     (gauge_momenta->momenta[i])*=(-1.0);
     }
  smartpack_thmatrix(momenta_packed, gauge_momenta);
  cuda_put_momenta();
  #endif
  } 


// to be used in reversibility tests
void save_momenta(void)
  {
  for(long int i=0; i<no_links; i++)
     {
     (gauge_momenta_save->momenta[i])=(gauge_momenta->momenta[i]);
     }
  } 
#endif

