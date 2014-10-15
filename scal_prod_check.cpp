#include <iostream>
#include <ctime>
//constants, global variables, geometric parameters et cetera 
#include "./Include/global_const.cc"
#include "./Include/global_macro.cc"
#include "./Include/global_var.cc"
#include "./Include/struct_def.cc"

//#define DEBUG_MODE  //useful in some cases
#include "./Init/init.cc"
#include "./Fermions/fermions.cc"
#include "./FermionMatrix/fermionmatrix.cc"
#include "./Extern/extern_c_func.cc"

using namespace std;

int main(){


   init(1);
   cout << "Initialized Random Gauge Matrix. (but it is useless here)\n ";
   Fermion* tempFermion1 = new Fermion();//initialized to 0
   Fermion* tempFermion2 = new Fermion();//initialized to 0

   tempFermion1->gauss();
   tempFermion2->gauss();

   cout << "Initialized Random Fermion Vectors.\n";

   vec3COM_soa soa1COM;
   vec3COM_soa soa2COM;

   tempFermion1->ferm_aos_to_soaCOM(&soa1COM);
   tempFermion2->ferm_aos_to_soaCOM(&soa2COM);


   // COMPUTATION OF THE SCALAR PRODUCT IN THE CPU C++ CODE
   complex<double> *cc_vector1;
   Vec3 vt1,vt2;
   cc_vector1=new complex<double>[sizeh];
   for(int ish=0;ish<sizeh;ish++){
     vt1=(tempFermion1->fermion[ish]);
     vt2=(tempFermion2->fermion[ish]);
     cc_vector1[ish]=c_scalprod(vt1,vt2);
   }
   global_sum(cc_vector1,sizeh);
   cout << "SCAL PROD CPU      = "  <<  cc_vector1[0]  << "\n";


   // COMPUTATION OF THE SCALAR PRODUCT IN THE OPENACC C CODE
   double p_re,p_im;
   scal_prod_openacc(&soa1COM,&soa2COM,&p_re,&p_im);
   cout << "SCAL PROD OPENACC  = ("  <<  p_re << "," << p_im << ")"  << "\n";





   delete tempFermion1;
   delete tempFermion2; 
   end();//deallocates all

   return 0;

}

