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

#include "./Inverter/inverter.cc"
#include "./OpenAcc/inverter_simple.cc"

using namespace std;

#define INV_CPU_ONLY
//#define INV_SIMPLE_OPENACC
//#define INV_FULL_OPENACC


int main(){


//All SU(3) links set to identity
   init(0);
   cout << "Initialized Random Gauge Matrix.\n";

   Fermion* tempFermion1 = new Fermion();//initialized to 0
   Fermion* tempFermion2 = new Fermion();//initialized to 0
   Fermion* tempFermion3 = new Fermion();//initialized to 0

   tempFermion1->gauss();
   tempFermion3->gauss();
   cout << "Initialized Random Fermion Vectors.\n";

   // CASO CPU ONLY
#ifdef INV_CPU_ONLY
   cout << "CPU ONLY INVERSION" << endl;
    invert(tempFermion2,tempFermion1,inv_single_double_prec,tempFermion3);
#endif

   //CASO SIMPLE OPENACC
#ifdef INV_SIMPLE_OPENACC
   cout << "SIMPLE OPENACC INVERSION" << endl;
   invert_openacc(tempFermion2,tempFermion1,inv_single_double_prec,tempFermion3);
#endif

   //CASO FULL OPENACC
#ifdef INV_FULL_OPENACC
   cout << "FULL OPENACC INVERSION" << endl;
   vec3COM_soa soa1COM;
   vec3COM_soa soa2COM;
   vec3COM_soa soa3COM;
   su3COM_soa conf_soaCOM[8];
   tempFermion1->ferm_aos_to_soaCOM(&soa1COM);
   tempFermion3->ferm_aos_to_soaCOM(&soa3COM);
   for(int index=0;index<8;index++)   gauge_conf->conf_aos_to_soaCOM(&conf_soaCOM[index],index);
   cout << inv_single_double_prec << endl;
   invert_openacc_full(conf_soaCOM,&soa2COM,&soa1COM,inv_single_double_prec,&soa3COM);
   tempFermion2->ferm_soaCOM_to_aos(&soa2COM);
#endif

   tempFermion2->saveToFile("invertedCPU.fer");
 

   /*
   vec3COM_soa soa1COM;
   vec3COM_soa soa2COM;
   su3COM_soa conf_soaCOM[8];

   tempFermion1->ferm_aos_to_soaCOM(&soa1COM);
   for(int index=0;index<8;index++)   gauge_conf->conf_aos_to_soaCOM(&conf_soaCOM[index],index);

   //   apply_Doe_openacc(conf_soaCOM,&soa1COM,&soa2COM);
   //   apply_Deo_openacc(conf_soaCOM,&soa1COM,&soa2COM);

   tempFermion2->ferm_soaCOM_to_aos(&soa2COM);
   for(int index=0;index<8;index++)   gauge_conf->conf_soaCOM_to_aos(&conf_soaCOM[index],index);

   //   tempFermion2->saveToFile("pippo2.fer");
   //   gauge_conf->saveToFile("TestConf_biconv.cnf");



   cout << "Saving fermion in StartFermion.fer\n";
   tempFermion1->saveToFile("StartFermion.fer");
 
   //   Fermion* tempFermion2 = new Fermion();

   cout << "Multiplying by Doe.\n";
   Doe(tempFermion2,tempFermion1);
   //   cout << "Multiplying by Deo.\n";
   //   Deo(tempFermion1,tempFermion2);
 
   cout << "Saving fermion in EndFermionCPU.fer\n";
   tempFermion2->saveToFile("EndFermionCPU.fer");
   */

   delete tempFermion1;
   delete tempFermion2; 
   delete tempFermion3; 
   end();//deallocates all

   return 0;

}

