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

#include "./Findminmax/findminmax.cc"
#include "./SimplifiedRationalApprox/rationalapprox.cc"
#include "./SimplifiedRationalApprox/rationalapprox_calc.cc"

#include "./Inverter/inverter.cc"
#include "./OpenAcc/inverter_simple.cc"

#include "./Meas/gaugemeas.cc"

using namespace std;



int main(){

  cout.precision(18);
  cerr.precision(18);
  rationalapprox_calc();

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

   tempFermion1->gauss();
   tempFermion3->gauss();
      //tempFermion1->z2noise();
      //tempFermion3->z2noise();
   cout << "Initialized Random Fermion Vectors.\n";
   //   tempFermion1->saveToFile("fer_32.fer");




   //////////////////////////////////////////////////////////////////////////////////////////////////
   /////////////////////// CONFRONTO INVERTITORE SINGOLO     ////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////////////////////////
   // CASO CPU ONLY
   //      cout << "CPU ONLY INVERSION" << endl;
   //      invert(tempFermion2_cpu,tempFermion1,inv_single_double_prec,tempFermion3);
   
   //CASO SIMPLE OPENACC
   //   cout << "SIMPLE OPENACC INVERSION" << endl;
   //   invert_openacc(tempFermion2_openacc_simple,tempFermion1,inv_single_double_prec,tempFermion3);
   
   //CASO FULL OPENACC
   /*
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
   */




   /*
   //////////////////////////////////////////////////////////////////////////////////////////////////
   /////////////////////// CONFRONTO FUNZIONI FIND MIN E FIND MAX  //////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////////////////////////
   double min, max, epsilon;
   double *minmax;
   minmax = new REAL [2];
   //TEST CALCOLO AUTOVALORE MASSIMO ED AUTOVALORE MINIMO
   //   findminmax(min, max);                                                                                                    
   cout << "CPU:      min = " << min << "    max = " << max << endl;

   findminmax_con_openacc(minmax);
   min = minmax[0];
   max = minmax[1];
   cout << "OPENACC:  min = " << min << "    max = " << max << endl;
   */




   //////////////////////////////////////////////////////////////////////////////////////////////////
   /////////////////////// CONFRONTO INVERTITORE MULTISHIFT    //////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////////////////////////
   create_phi();
   cout << "Created phi\n";

   su3COM_soa conf_soaCOM[8];
   for(int index=0;index<8;index++)   gauge_conf->conf_aos_to_soaCOM(&conf_soaCOM[index],index);

   calc_plaquette_openacc(conf_soaCOM);
   double ps,pt;
   gauge_conf->calc_plaq(ps,pt);
   cout << "PlaquetteCPU  " << (ps+pt)*0.5 << endl;

   vec3COM_soa soa1COM;
   vec3COM_soa soa2COM;
   tempFermion1->ferm_aos_to_soaCOM(&soa1COM);

   RationalApprox approx;
   COM_RationalApprox *COM_approx;
   COM_approx = new COM_RationalApprox[1];

   //   approx.md_inv_approx_coeff();
   approx.first_inv_approx_coeff();
   cerr << approx << endl;

   convert_RationalApprox_to_COM_RationalApprox(&COM_approx[0],approx);

   int ips=(int)(2.0*casuale());
   int ish=(int)((sizeh+1)*casuale());
   cout << "ECCOLIIIII   ----------->  " << ips << "    " << ish << endl;
   cout << COM_approx[0].COM_approx_order << endl ;

   cout << "Rat approx done\n";

   COM_ShiftFermion    COMMON_shift;
   COM_MultiFermion    COMMON_multi;
   COM_MultiFermion    COMMON_multi_out;
   COM_ShiftMultiFermion   COMMON_shiftmulti;


   MultiFermion *fermion_app;
   fermion_app = new MultiFermion;




   cout << "Chosen directions\n";

   for(int icomp=0;icomp<3;icomp++)
     cout << (fermion_phi->fermion[ips][ish].comp[icomp]).real() << "   " << (fermion_phi->fermion[ips][ish].comp[icomp]).imag() << endl;


   fermion_shiftmulti->ferm_ShiftMulti_to_ShiftMultiCOM(&COMMON_shiftmulti);
   fermion_phi->ferm_Multi_to_MultiCOM(&COMMON_multi);

   first_inv_approx_calc_openacc(conf_soaCOM,&COMMON_multi_out,&COMMON_multi,residue_metro,COM_approx);

   fermion_shiftmulti->ferm_ShiftMultiCOM_to_ShiftMulti(&COMMON_shiftmulti);
   fermion_app->ferm_MultiCOM_to_Multi(&COMMON_multi_out);
   for(int icomp=0;icomp<3;icomp++)
     cout << (fermion_app->fermion[ips][ish].comp[icomp]).real() << "   " << (fermion_app->fermion[ips][ish].comp[icomp]).imag() << endl;

   cout << "Conversions done!" << endl;



   //   multips_shifted_invert(fermion_shiftmulti, fermion_phi,residue_metro,approx);
   int iter, pseudofermion,i;
   Vec3 vr_1;

   for(pseudofermion=0; pseudofermion<no_ps; pseudofermion++)
     {
       for(i=0; i<sizeh; i++)
	 {
	   vr_1=(approx.RA_a0)*(fermion_phi->fermion[pseudofermion][i]);
	   for(iter=0; iter<(approx.approx_order); iter++)
	     {
	       vr_1+=(approx.RA_a[iter])*(fermion_shiftmulti->fermion[pseudofermion][iter][i]);
	     }
	   fermion_chi->fermion[pseudofermion][i]=vr_1;
	 }
     }
   cerr << "phi  " << fermion_phi->fermion[0][i] << endl ;
   cerr << "shi  " << fermion_shiftmulti->fermion[0][0][i] << endl ;

   cout << "--- INVERTITORE CPU --- --- --- RISULTATI DI MULTIPS" << endl << endl;

   /*
   ofstream oferm;
   string ofermname="inverted_openacc.fer";
   oferm.open(ofermname.c_str());
   oferm.precision(18);
   for(int ips=0;ips<2;ips++){
     for(int ish=0;ish<sizeh;ish++){
       oferm <<   (fermion_app->fermion[ips][ish].comp[0]).real() << "     " <<    (fermion_app->fermion[ips][ish].comp[0]).imag() << endl;
       oferm <<   (fermion_app->fermion[ips][ish].comp[1]).real() << "     " <<    (fermion_app->fermion[ips][ish].comp[1]).imag() << endl;
       oferm <<   (fermion_app->fermion[ips][ish].comp[2]).real() << "     " <<    (fermion_app->fermion[ips][ish].comp[2]).imag() << endl;
     }
   }
   oferm.close();
   */

   double diff_CPU_FULL   = difference_multi(fermion_chi,fermion_app);
   cout << "--- multinv differences: openacc - cpu      --------" << endl << endl;
   cout << diff_CPU_FULL  << endl;


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

