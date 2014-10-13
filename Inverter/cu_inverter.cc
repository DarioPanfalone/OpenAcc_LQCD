// shifted inverter for multiple pseudofermions: for each pseudofermion component of "in" 
// solve the equation {(M^dag M) + RA_b[i]}out=in   
// Jegerlehner hep-lat/9612014 "Krylov space solvers for shifted linear systems

#ifdef USE_GPU
#include"../Packer/packer.h"
#include"../Cuda/cuda_inverter.h"
#endif

void cu_invert(REAL res)
  {
  #if ((defined DEBUG_MODE) || (defined DEBUG_INVERTER))
  cout <<"DEBUG: inside cu_invert ..."<<endl;

  #endif

  #ifdef TIMING_CUDA_CPP
  clock_t time_start, time_finish;
  time_start=clock();
  #endif

  long int i;
  int iter, cg;

  // start loop on pseudofermions

     if(res<inv_single_double_prec)cuda_inverter_d(res, &cg);
     else  cuda_inverter(res, &cg);

     if(cg==max_cg)
       {
       ofstream err_file;
       err_file.open(QUOTEME(ERROR_FILE), ios::app);   
       err_file  << "WARNING: maximum number of iterations reached in cuda_multips_shift_invert\n";
       err_file.close();
       }

  #ifdef TIMING_CUDA_CPP
  time_finish=clock();
  cout << "time for cu_invert = " << ((REAL)(time_finish)-(REAL)(time_start))/CLOCKS_PER_SEC << " sec.\n";
  #endif


  }

