#ifndef CONF_CC_
#define CONF_CC_

#include <fstream>

#include "../Su3/su3.cc"
#include "../Include/global_var.cc"
#include "../Include/global_macro.cc"
#include "../Exception/exception.cc"

using namespace std;


class Conf {
  // base constructor private!
  Conf(void);
public:
  Conf(int initMode);

  void allIdentities();
  void allIdentitiesExcept(int ii,int sub);

  Su3 u_save[no_links];
  Su3 u_work[no_links];

  void unitarize_with_eta(void);
  void saveToFile(const char*);

  void conf_aos_to_soaCOM(su3COM_soa *out,int matdir) const;
  void conf_soaCOM_to_aos(su3COM_soa const* const in,int matdir);


};

void Conf::conf_aos_to_soaCOM(su3COM_soa *out,int matdir) const{
  int indice;
  for( int i =0 ; i < sizeh ; i++){
    indice = i + sizeh * matdir;
    out->r0.c0[i].Re = (u_work[indice].comp[0][0]).real();
    out->r1.c0[i].Re = (u_work[indice].comp[1][0]).real();
    out->r2.c0[i].Re = (u_work[indice].comp[2][0]).real();
    out->r0.c1[i].Re = (u_work[indice].comp[0][1]).real();
    out->r1.c1[i].Re = (u_work[indice].comp[1][1]).real();
    out->r2.c1[i].Re = (u_work[indice].comp[2][1]).real();
    out->r0.c2[i].Re = (u_work[indice].comp[0][2]).real();
    out->r1.c2[i].Re = (u_work[indice].comp[1][2]).real();
    out->r2.c2[i].Re = (u_work[indice].comp[2][2]).real();

    out->r0.c0[i].Im = (u_work[indice].comp[0][0]).imag();
    out->r1.c0[i].Im = (u_work[indice].comp[1][0]).imag();
    out->r2.c0[i].Im = (u_work[indice].comp[2][0]).imag();
    out->r0.c1[i].Im = (u_work[indice].comp[0][1]).imag();
    out->r1.c1[i].Im = (u_work[indice].comp[1][1]).imag();
    out->r2.c1[i].Im = (u_work[indice].comp[2][1]).imag();
    out->r0.c2[i].Im = (u_work[indice].comp[0][2]).imag();
    out->r1.c2[i].Im = (u_work[indice].comp[1][2]).imag();
    out->r2.c2[i].Im = (u_work[indice].comp[2][2]).imag();

  }
}

void Conf::conf_soaCOM_to_aos(su3COM_soa const* const in,int matdir){
  int indice;
  for( int i =0 ; i < sizeh ; i++){
    indice = i + sizeh * matdir;
    u_work[indice].comp[0][0] = complex<double>(in->r0.c0[i].Re,in->r0.c0[i].Im);
    u_work[indice].comp[1][0] = complex<double>(in->r1.c0[i].Re,in->r1.c0[i].Im);
    u_work[indice].comp[2][0] = complex<double>(in->r2.c0[i].Re,in->r2.c0[i].Im);
    u_work[indice].comp[0][1] = complex<double>(in->r0.c1[i].Re,in->r0.c1[i].Im);
    u_work[indice].comp[1][1] = complex<double>(in->r1.c1[i].Re,in->r1.c1[i].Im);
    u_work[indice].comp[2][1] = complex<double>(in->r2.c1[i].Re,in->r2.c1[i].Im);
    u_work[indice].comp[0][2] = complex<double>(in->r0.c2[i].Re,in->r0.c2[i].Im);
    u_work[indice].comp[1][2] = complex<double>(in->r1.c2[i].Re,in->r1.c2[i].Im);
    u_work[indice].comp[2][2] = complex<double>(in->r2.c2[i].Re,in->r2.c2[i].Im);
    u_save[indice] = u_work[indice];
  }
}


void Conf::allIdentities(){

   for(int i=0; i<no_links; i++){

       u_save[i].one();
       u_save[i]*=eta[i];
       u_work[i]=u_save[i];

   }

}


void Conf::allIdentitiesExcept(int ii, int sub){

   this->allIdentities();

   u_save[ii].isigmaZ(sub); 
   u_work[ii]=u_save[ii];

   cout <<   u_save[ii] ; 

}




// constructor & initialization
Conf::Conf(int initMode)
 {
  #ifdef DEBUG_MODE
  cout << "DEBUG: inside Conf::Conf ..."<<endl;
  #endif

 long int i;

 if(initMode==0) // start from ordered conf
   {
   #ifdef DEBUG_MODE
   cout << "\tstart from ordered conf ..."<<endl;
   #endif
   for(i=0; i<no_links; i++)
      {
      Su3 aux;               // aux is a small perturbation
      aux.rand_matrix();
      aux*=0.01;
      u_save[i].one();
      u_save[i]+=aux;
      u_save[i].sunitarize();
      u_save[i]*=eta[i];
      u_work[i]=u_save[i];
      }
   }

 if(initMode==1) // start from random conf
   {
   #ifdef DEBUG_MODE
   cout << "\tstart from random conf ..."<<endl;
   #endif
   for(i=0; i<no_links; i++)
      {
      u_save[i].rand_matrix();
      u_save[i]*=eta[i];
      u_work[i]=u_save[i];
      }
   update_iteration=0;
   }

 if(initMode==2) // start from saved conf.
   {
   #ifdef DEBUG_MODE
   cout << "\tstart from stored conf ..."<<endl;
   #endif
   int nx_l, ny_l, nz_l, nt_l;
   REAL beta_l, mass_l, no_flavours_l;
   long int r;
   ifstream file;

   file.open(QUOTEME(CONF_FILE), ios::in);

   file >> nx_l;
   file >> ny_l;
   file >> nz_l;
   file >> nt_l;
   file >> beta_l;
   file >> mass_l;
   file >> no_flavours_l;
   file >> update_iteration;

   if(nx!=nx_l || ny!=ny_l || nz!=nz_l || nt !=nt_l)
     {
     throw stored_conf_not_fit;
     }

   for(r=0; r<no_links; r++)
      {
      file >> u_save[r];

      u_work[r].sunitarize();
      u_work[r].row_multiply(2,eta[r]);

      u_work[r]=u_save[r];
      }
   file.close();
   }

 #ifdef DEBUG_MODE
 cout << "\tterminated Conf::Conf"<<endl;
 #endif
 }



// unitarize and restore staggered phases
void Conf::unitarize_with_eta(void)
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside Conf::unitarize_with_eta ..."<<endl;
 #endif

 long int r;

 for(r=0; r<no_links; r++)
    {
    u_work[r].sunitarize();
    u_work[r].row_multiply(2,eta[r]);
    }
 #ifdef DEBUG_MODE
 cout << "\tterminated Conf::unitarize_with_eta"<<endl;
 #endif
 }

void Conf::saveToFile(const char *filename){

 ofstream file;

 file.open(filename, ios::out);
 file.precision(16);

 file << nx << " " << ny << " " << nz << " " << nt << " " << beta << " " << mass << " " <<no_flavours << " " <<  update_iteration <<"\n";

 for(int r=0; r<no_links; r++)
    {
    file << u_save[r];
    }

 file.close();

}


#endif
