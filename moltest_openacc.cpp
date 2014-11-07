#include <iostream>
#include <ctime>
//constants, global variables, geometric parameters et cetera 
#include "./Include/global_const.cc"
#include "./Include/global_macro.cc"
#include "./Include/global_var.cc"

//#define DEBUG_MODE  //useful in some cases
#include "./Init/init.cc"
#include "./Fermions/fermions.cc"
#include "./FermionMatrix/fermionmatrix.cc"
#include "./Extern/extern_c_func.cc"

using namespace std;

int main(){


//All SU(3) links set to identity
   init(1);
   cout << "Initialized Random Gauge Matrix.\n";
//If needed the gauge field can be loaded from file
//using init(2), in that case a 
//file named 'config' will be looked for.
//The gauge field can also be saved :
   gauge_conf->saveToFile("TestConf.cnf");
   Fermion* tempFermion1 = new Fermion();//initialized to 0
   Fermion* tempFermion2 = new Fermion();//initialized to 0

   tempFermion1->gauss();
   tempFermion1->saveToFile("pippo.fer");
   cout << "Initialized Random Fermion Vector.\n";

   vec3COM_soa soa1COM;
   vec3COM_soa soa2COM;
   su3COM_soa conf_soaCOM[8];

   tempFermion1->ferm_aos_to_soaCOM(&soa1COM);
   for(int index=0;index<8;index++)   gauge_conf->conf_aos_to_soaCOM(&conf_soaCOM[index],index);

   apply_Doe_openacc(conf_soaCOM,&soa1COM,&soa2COM);
   //   apply_Deo_openacc(conf_soaCOM,&soa1COM,&soa2COM);

   tempFermion2->ferm_soaCOM_to_aos(&soa2COM);
   for(int index=0;index<8;index++)   gauge_conf->conf_soaCOM_to_aos(&conf_soaCOM[index],index);

   tempFermion2->saveToFile("pippo2.fer");
   gauge_conf->saveToFile("TestConf_biconv.cnf");



   cout << "Saving fermion in StartFermion.fer\n";
   tempFermion1->saveToFile("StartFermion.fer");
 
   //   Fermion* tempFermion2 = new Fermion();

   cout << "Multiplying by Doe.\n";
   Doe(tempFermion2,tempFermion1);
   //   cout << "Multiplying by Deo.\n";
   //   Deo(tempFermion1,tempFermion2);
 
   cout << "Saving fermion in EndFermionCPU.fer\n";
   tempFermion2->saveToFile("EndFermionCPU.fer");

   delete tempFermion1;
   delete tempFermion2; 
   end();//deallocates all

   return 0;

}

