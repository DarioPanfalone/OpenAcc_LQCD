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
#include "../Rand/random.cc"
#include "../Geometry/geometry.cc"
#include "../Include/global_var.cc"
#include "../Exception/exception.cc"

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

     // allocate auxiliary global fermion 
     loc_r=new Fermion;
     loc_h=new Fermion;
     loc_s=new Fermion;
     loc_p=new Fermion;

     // auxiliary vectors for global sums
     d_vector1=new double[size];

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


   // clear gauge configuration
  delete gauge_conf;

  // clear geometry variables
  end_geo();

  #ifdef DEBUG_MODE
  cout << "\tterminated end"<< endl;
  #endif
  }

#endif


