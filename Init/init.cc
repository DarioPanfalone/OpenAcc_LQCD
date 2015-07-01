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


     // allocate gauge configuration
     gauge_conf=new Conf(startMode);

     // auxiliary vectors for global sums
     d_vector1=new double[size];
     d_vector2=new double[size];
     plaq_test=new double[size];

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

  // clear gauge configuration
  delete gauge_conf;

  // clear geometry variables
  end_geo();

  #ifdef DEBUG_MODE
  cout << "\tterminated end"<< endl;
  #endif
  }

#endif


