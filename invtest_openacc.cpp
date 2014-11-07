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



int main(){


//All SU(3) links set to identity
   init(0);
   //   gauge_conf->saveToFile("TestConf_32.cnf");
   cout << "Initialized Random Gauge Matrix.\n";

   //   Fermion() initializes to 0
 
   Fermion* tempFermion1 = new Fermion(); // Fermion to invert
   Fermion* tempFermion2_cpu = new Fermion(); // Result fermion for cpu
   Fermion* tempFermion2_openacc_simple = new Fermion();// result fermions for openacc_simple
   Fermion* tempFermion2_openacc_full = new Fermion();// result fermion for openacc_full
   Fermion* tempFermion3 = new Fermion();// Trial solution

   //   tempFermion1->gauss();
   //   tempFermion3->gauss();
   tempFermion1->z2noise();
   tempFermion3->z2noise();
   cout << "Initialized Random Fermion Vectors.\n";
   //   tempFermion1->saveToFile("fer_32.fer");

   // CASO CPU ONLY
   cout << "CPU ONLY INVERSION" << endl;
   invert(tempFermion2_cpu,tempFermion1,inv_single_double_prec,tempFermion3);
   
   //CASO SIMPLE OPENACC
   cout << "SIMPLE OPENACC INVERSION" << endl;
   invert_openacc(tempFermion2_openacc_simple,tempFermion1,inv_single_double_prec,tempFermion3);
   
   //CASO FULL OPENACC
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
   tempFermion2_openacc_full->ferm_soaCOM_to_aos(&soa2COM);

   double diff_CPU_SIMPLE = difference(tempFermion2_cpu,tempFermion2_openacc_simple);
   double diff_CPU_FULL   = difference(tempFermion2_cpu,tempFermion2_openacc_full);

   cout << "Differenze: \n  |CPU - SIMPLE| = " << diff_CPU_SIMPLE << "\n  |CPU - FULL| = " << diff_CPU_FULL << endl;

   tempFermion2_cpu->saveToFile("invertedCPU.fer");
   tempFermion2_openacc_simple->saveToFile("invertedOPENACC_SIMPLE.fer");
   tempFermion2_openacc_full->saveToFile("invertedOPENACC_FULL.fer");

 

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
   delete tempFermion2_cpu; 
   delete tempFermion2_openacc_simple; 
   delete tempFermion2_openacc_full; 
   delete tempFermion3; 
   end();//deallocates all

   return 0;

}

