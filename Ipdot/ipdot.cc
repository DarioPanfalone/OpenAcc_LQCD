#ifdef USE_GPU
#include"../Cuda/cuda_init.h"
#include"../Cuda/cuda_deriv.h"
#include"../Packer/packer.h"
#endif

// defined in FermionForce/fermionforce.cc
void fermionforce(int moltiplico); // 0 --> la vecchia forza fermionica
                                   // 1 --> restituisce la forza già moltiplicata per gli u_work (serve quando b_quantum =/= 0 perchè gli u_work sono diversi per i 2 pseudofermioni)
class Ipdot {
public:
  Su3 ipdot[no_links];

  Ipdot(void);

  friend void calc_ipdot_gauge(void);
  friend void calc_ipdot_fermion(void);
  friend void calc_ipdot(void);

  // defined in FermionForce/fermionforce.cc
  friend void fermionforce(int moltiplico);

  // defined in Momenta/momenta.cc
  friend void momenta_sum_multiply(const complex<REAL> p);

  // defined in Packer/packer.cc
  friend void smartpack_tamatrix(float out[8*no_links], Ipdot *in);
  friend void smartunpack_tamatrix(Ipdot *out, float in[8*no_links]);

  void ipdot_aos_to_soaCOM(tamatCOM_soa *out,int matdir) const;
  void ipdot_soaCOM_to_aos(tamatCOM_soa const * const in,int matdir);

};


void Ipdot::ipdot_aos_to_soaCOM(tamatCOM_soa *out,int matdir) const{
  int indice;
  for( int i =0 ; i < sizeh ; i++){
    indice = i + sizeh * matdir;
    out->c01[i].Re = (ipdot[indice].comp[0][1]).real();
    out->c01[i].Im = (ipdot[indice].comp[0][1]).imag();
    out->c02[i].Re = (ipdot[indice].comp[0][2]).real();
    out->c02[i].Im = (ipdot[indice].comp[0][2]).imag();
    out->c12[i].Re = (ipdot[indice].comp[1][2]).real();
    out->c12[i].Im = (ipdot[indice].comp[1][2]).imag();
    out->rc00[i] = (ipdot[indice].comp[0][0]).imag();
    out->rc11[i] = (ipdot[indice].comp[1][1]).imag();
  }
}


void Ipdot::ipdot_soaCOM_to_aos(tamatCOM_soa const * const in,int matdir){
  int indice;
  for( int i =0 ; i < sizeh ; i++){
    indice = i + sizeh * matdir;
    // assign
    ipdot[indice].comp[0][1] = complex<double>(in->c01[i].Re,in->c01[i].Im);
    ipdot[indice].comp[0][2] = complex<double>(in->c02[i].Re,in->c02[i].Im);
    ipdot[indice].comp[1][2] = complex<double>(in->c12[i].Re,in->c12[i].Im);
    ipdot[indice].comp[0][0] = complex<double>(0.0,in->rc00[i]);
    ipdot[indice].comp[1][1] = complex<double>(0.0,in->rc11[i]);
    //reconstruct the other elements
    ipdot[indice].comp[1][0] = -conj(ipdot[indice].comp[0][1]);
    ipdot[indice].comp[2][0] = -conj(ipdot[indice].comp[0][2]);
    ipdot[indice].comp[2][1] = -conj(ipdot[indice].comp[1][2]);
    ipdot[indice].comp[2][2] = -ipdot[indice].comp[0][0] -ipdot[indice].comp[1][1];
  }
}


// base constructor
Ipdot::Ipdot(void)
 {
 for(long int i=0; i<no_links; i++)
    {
    ipdot[i].zero();
    }
 }


// calculate the gauge part of ipdot
void calc_ipdot_gauge(void)
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside calc_ipdot_gauge ..."<<endl;
 #endif

 #ifdef TIMING_CUDA_CPP
 clock_t time_start, time_finish;
 time_start=clock();
 #endif

 #ifdef USE_GPU 
   cuda_gauge_deriv(0);
 #else
   long int r;
   Su3 aux;

   calc_staples(0); // 0 ---> non si aggiunge il contributo topologico alle staples ; !0 ---> ci mette anche il pezzo topo dei quadrifogli
   for(r=0; r<no_links; r++)
      {
      // ipdot = traceless anti-hermitian part of beta_by_three*(u_work*staple) 
      aux=(gauge_conf->u_work[r]);
      aux*=(beta_by_three);         // PLUS SIGN DUE TO STAGGERED PHASES
      aux*=(gauge_staples->staples[r]);
      aux.ta();
      (gauge_ipdot->ipdot[r])=aux;
      }
 #endif

 #ifdef TIMING_CUDA_CPP
 time_finish=clock();
 cout << "time for calc_ipdot_gauge = " << ((REAL)(time_finish)-(REAL)(time_start))/CLOCKS_PER_SEC << " sec.\n";
 #endif

 #ifndef USE_GPU
   #ifdef PARAMETER_TEST
   int mu;
   for(r=0; r<size; r++)
      {
      d_vector1[r]=0.0;
      for(mu=0; mu<4; mu++)
         {
         d_vector1[r]+=(gauge_ipdot->ipdot[r+mu*size]).l2norm2();
         }
      }
   global_sum(d_vector1, size);

   cout << "L2 norm of GAUGE force = " << sqrt(d_vector1[0])<<endl;
   #endif
 #endif

 #ifdef DEBUG_MODE
 cout << "\tterminated calc_ipdot_gauge"<<endl;
 #endif
 }


// calculate the fermionic part of ipdot
void calc_ipdot_fermion(void)
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside calc_ipdot_fermion ..."<<endl;
 #endif

 #ifdef TIMING_CUDA_CPP
 clock_t time_start, time_finish;
 time_start=clock();
 #endif

 #ifdef USE_GPU 
   fermionforce(0);
 #else
   long int r;
   Su3 aux;

   // initialize to zero
   for(r=0; r<no_links; r++)
      {
      (gauge_ipdot->ipdot[r]).zero();
      }
cout << "Ipdot messo a zero " << endl;

   // add fermionic term
   // Se si fa campo magnetico ci sono 2 contributi a ipdot
   // che vengono dai due pseudofermioni.
   // Gli u_work che sono coinvolti per i 2 pseudofermioni sono diversi
   // a causa delle fasi magnetiche
   // Quindi se faccio MAGN mi faccio ritornare da fermionforce il valore
   // dei due contributi di ipdot opportunamente moltiplicati per il giusto link e poi sommati

  #ifdef MAGN            
     fermionforce(1);
  #endif

  #ifndef MAGN
     fermionforce(0);
  #endif
     int index=19;
cout << "Forza fermionica calcolata " << endl;
 cout << "Link" << endl;
 cout << (gauge_conf->u_work[19]) << endl;
 cout << "Force" << endl;
 cout << (gauge_ipdot->ipdot[19]) << endl;

   // multiply by u_work and take Traceless-Antihermitean part
   for(r=0; r<no_links; r++)
     {
       #ifndef MAGN
         aux=(gauge_conf->u_work[r]);
         aux*=(gauge_ipdot->ipdot[r]);
         aux.ta();
         (gauge_ipdot->ipdot[r])=aux;
       #endif

       #ifdef MAGN
         aux=(gauge_ipdot->ipdot[r]);
         aux.ta();
         (gauge_ipdot->ipdot[r])=aux;
       #endif

     }
cout << "Moltiplicato per il link e presa la parte TA" << endl;

 #endif

 #ifdef TIMING_CUDA_CPP
 time_finish=clock();
 cout << "time for calc_ipdot_fermion = " << ((REAL)(time_finish)-(REAL)(time_start))/CLOCKS_PER_SEC << " sec.\n";
 #endif

 #ifdef PARAMETER_TEST
  #ifndef USE_GPU
    int mu;
    for(r=0; r<size; r++)
       {
       d_vector1[r]=0.0;
       for(mu=0; mu<4; mu++)
          {
          d_vector1[r]+=(gauge_ipdot->ipdot[r+mu*size]).l2norm2();
          }
       }
    global_sum(d_vector1, size);

    cout << "L2 norm of FERMION force = " << sqrt(d_vector1[0])<<endl;
    #endif
 #endif

 #ifdef DEBUG_MODE
 cout << "\tterminated calc_ipdot_fermion"<<endl;
 #endif
 }


// questa funzione non viene mai chiamata se si fa il metodo che separa l'evoluzione per la forza fermionica da quella di gauge
// quindi per ora non la modifico per quanto riguarda le parti di campo magnetico.
// calculate the complete ipdot
void calc_ipdot(void)
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside calc_ipdot ..."<<endl;
 #endif

 #ifdef TIMING_CUDA_CPP
 clock_t time_start, time_finish;
 time_start=clock();
 #endif

 #ifdef USE_GPU 
   fermionforce(0);
   cuda_gauge_deriv(1);  // add gauge part
 #else
   long int r;
   Su3 aux;

   calc_staples(0);  // 0 ---> non si aggiunge il contributo topologico alle staples ; !0 ---> ci mette anche il pezzo topo dei quadrifogli

   // add gauge_part
   for(r=0; r<no_links; r++)
      {
      (gauge_ipdot->ipdot[r])=beta_by_three*(gauge_staples->staples[r]);
      }

   // add fermionic term   
   fermionforce(0);

   // multiply by u_work and take TA
   for(r=0; r<no_links; r++)
      {
      aux=(gauge_conf->u_work[r]);
      aux*=(gauge_ipdot->ipdot[r]);
      aux.ta();
      (gauge_ipdot->ipdot[r])=aux;
      }
 #endif

 #ifdef TIMING_CUDA_CPP
 time_finish=clock();
 cout << "time for calc_ipdot = " << ((REAL)(time_finish)-(REAL)(time_start))/CLOCKS_PER_SEC << " sec.\n";
 #endif

 #ifdef DEBUG_MODE
 cout << "\tterminated calc_ipdot"<<endl;
 #endif
 }
