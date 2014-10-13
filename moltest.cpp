#include <iostream>

//constants, global variables, geometric parameters et cetera 
#include "./Include/global_const.cc"
#include "./Include/global_macro.cc"
#include "./Include/global_var.cc"

//#define DEBUG_MODE  //useful in some cases

#include "./Init/init.cc"
#include "./Fermions/fermions.cc"
#include "./FermionMatrix/fermionmatrix.cc"

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

   tempFermion1->gauss();
   tempFermion1->saveToFile("pippo.fer");
   cout << "Initialized Random Fermion Vector.\n";

   Vec3_soa soa1;
   tempFermion1->aos_to_soa(&soa1);
   tempFermion1->soa_to_aos(&soa1);

   tempFermion1->saveToFile("pippo2.fer");

   cout << "Saving fermion in StartFermion.fer\n";
   tempFermion1->saveToFile("StartFermion.fer");
 
   Fermion* tempFermion2 = new Fermion();

   cout << "Multiplying by Doe.\n";
   Doe(tempFermion2,tempFermion1);
   cout << "Multiplying by Deo.\n";
   Deo(tempFermion1,tempFermion2);
 
   cout << "Saving fermion in EndFermionCPU.fer\n";
   tempFermion1->saveToFile("EndFermionCPU.fer");

   delete tempFermion1;
   delete tempFermion2; 
   end();//deallocates all

   return 0;

}

