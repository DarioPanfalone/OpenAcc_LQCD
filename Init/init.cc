#ifndef INIT_CC_
#define INIT_CC_

#include <iostream>

#include "../Include/global_var.cc"
#include "../Include/parameters.cc"
#include "../Geometry/geometry.cc"
#include "../Rand/random.cc"
#include "../Su3/su3.cc"
#include "../Fermions/fermions.cc"
#include "../Conf/conf.cc"
#include "../FermionMatrix/fermionmatrix.cc"  
#include "../Extern/extern_c_func.cc" 
#include "../Findminmax/findminmax.cc"
#include "../SimplifiedRationalApprox/rationalapprox.cc"      
#include "../SimplifiedRationalApprox/rationalapprox_calc.cc" 
#include "../Staples/staples.cc"
#include "../Rand/random.cc"
#include "../Geometry/geometry.cc"
#include "../Include/global_var.cc"
#include "../Exception/exception.cc"
#include "../Ipdot/ipdot.cc"
#include "../FermionForce/fermionforce.cc"
#include "../Momenta/momenta.cc"
#include "../MD_integrators/leapfrog.cc"
#include "../MD_integrators/minimum_norm2.cc"
#include "../MD_integrators/multistep_2MN.cc"
#include "../MD_integrators/multistep_4MN.cc"
#include "../Action/action.cc"
#include "../Update/update.cc"

#include "../Inverter/inverter.cc"    
#include "../OpenAcc/inverter_simple.cc"  

using namespace std;

int init(int startMode = 1)
  {
  #ifdef DEBUG_MODE
  cout << "DEBUG: inside init ..."<<endl;
  #endif
  
  try{
     // initialize random number generator
     initrand(rand_seed);

     // initialize geometry
     // defined in Geometry/geometry.cc
     init_geo();


     cerr << " Random number inside init and before gauge conf generation :    " << casuale() << endl;
     // allocate gauge configuration
     gauge_conf=new Conf(startMode);
     cerr << " Random number inside init and after gauge conf generation :    " << casuale() << endl;

     // allocate staples
     gauge_staples=new Staples();
     // allocate ipdot
     gauge_ipdot=new Ipdot();
     // allocate momenta
     gauge_momenta=new Momenta();

     // allocate auxiliary global fermion 
     loc_r=new Fermion;
     loc_h=new Fermion;
     loc_s=new Fermion;
     loc_p=new Fermion;
     fermion_phi=new MultiFermion;
     fermion_chi=new MultiFermion;
     p_shiftferm=new ShiftFermion;
     fermion_shiftmulti=new ShiftMultiFermion;

     // auxiliary vectors for global sums
     d_vector1=new double[size];
     d_vector2=new double[size];
     plaq_test=new double[size];

#ifdef USEGPU
     gauge_field_packed=new float[12*no_links*2];            // 2 since 1double~2float
     shift_table = new int[8*size];//the new nnp and nnm, but to copy on the device
#endif
     // test if the simulation details are applicable
     test_param(); 

     }
  catch (exception& e)
     {
     cout << e.what() << endl;
     #ifdef DEBUG_MODE
     cout << "\tterminated init with exceptions"<< endl;
     #endif
     return 1;
     }

  #ifdef DEBUG_MODE
  cout << "\tterminated init without exceptions"<< endl;
  #endif
  return 0;
  }



void end(void)
  {
  #ifdef DEBUG_MODE
  cout << "DEBUG: inside end ..."<< endl;
  #endif


  // clear auxiliary vectors
  delete [] d_vector1;
  delete [] d_vector2;
  delete [] plaq_test;

  // clear auxiliary global fermions
  delete loc_r;
  delete loc_h;
  delete loc_s;
  delete loc_p;

  // clear fermions
  delete fermion_phi;
  delete fermion_chi;

  // clear shifted fermions
  delete p_shiftferm;
  delete fermion_shiftmulti;

  // clear ipdot 
  delete gauge_ipdot;

  // clear momenta 
  delete gauge_momenta;


   // clear gauge configuration
  delete gauge_conf;
  // clear staples
  delete gauge_staples;

  // clear geometry variables
  end_geo();

  #ifdef DEBUG_MODE
  cout << "\tterminated end"<< endl;
  #endif
  }

#endif


