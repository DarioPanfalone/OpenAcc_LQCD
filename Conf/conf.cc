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

  Su3 u_save[no_links];
  Su3 u_work[no_links];

  Conf(int initMode);

  // defined in Action/action.cc
  void loc_action(double *value, int init) const;

  void  calc_plaq(REAL &pls, REAL &plt) const;
  void  calc_plaq_uwork(REAL &pls, REAL &plt) const;
  //  void  calc_all_plaq(REAL &pxyr,REAL &pxyi,REAL &pxzr,REAL &pxzi,REAL &pxtr,REAL &pxti,REAL &pyzr,REAL &pyzi,REAL &pytr,REAL &pyti,REAL &pztr,REAL &pzti);
  void  calc_poly(REAL &re, REAL &im) const;
  Su3 get_staple(int index, int mu);
  void  cooling();

  void allIdentities();
  void allIdentitiesExcept(int ii,int sub);


  void unitarize_with_eta(void);
  void saveToFile(const char*);

  void conf_aos_to_soaCOM(su3COM_soa *out,int matdir) const;
  void conf_soaCOM_to_aos(su3COM_soa const* const in,int matdir);

  void save(void);
  void copy_saved(void);
  void write(void);
  void write_last(void);
  void print(void);


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
      aux*=0.1;
      //aux*=0.0;
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

//Calcola la staple associata al link che parte dal sito index in direzione mu               
Su3 Conf::get_staple(int index, int mu)
{
  Su3 aux,staple;
  int nu, pos;
  long int index_mu, index_nu, helper;

  pos = index;
  index_mu=mu*size;
  aux.zero();
  staple.zero();
  for(nu=0; nu<4; nu++)
    {
      if(nu!=mu)
        {
          index_nu=nu*size;
          aux = (gauge_conf->u_work[index_nu + nnp[pos][mu]]);  // 1        \
          aux*=~(gauge_conf->u_work[index_mu + nnp[pos][nu]]);  // 2        \
          aux*=~(gauge_conf->u_work[index_nu + pos]);           // 3        \
          staple += aux;

          helper=nnm[pos][nu];
          aux =~(gauge_conf->u_work[index_nu + nnp[helper][mu]]); // 1      \
          aux*=~(gauge_conf->u_work[index_mu + helper]);          // 2      \
          aux*= (gauge_conf->u_work[index_nu + helper]);          // 3      \             
          staple += aux;
        }
    }
  return staple;
}


// copy u_work to u_save
void Conf::save(void)
{
 #ifdef DEBUG_MODE
  cout << "DEBUG: inside Conf::save ..."<<endl;
 #endif
  long int r;

  for(r=0; r<no_links; r++)
    {
      u_save[r]=u_work[r];
    }
 #ifdef DEBUG_MODE
  cout << "\tterminated Conf::save"<<endl;
 #endif
}

// copy u_save to u_work
void Conf::copy_saved(void)
{
 #ifdef DEBUG_MODE
  cout << "DEBUG: inside Conf::copy_saved ..."<<endl;
 #endif

  long int r;

  for(r=0; r<no_links; r++)
    {
      u_work[r]=u_save[r];
    }
 #ifdef DEBUG_MODE
  cout << "\tterminated Conf::copy_saved"<<endl;
 #endif
}

// write configuration to conf_file_0 or conf_file_1 (defined in Include/global_macro.cc) 
void Conf::write(void)
{
 #ifdef DEBUG_MODE
  cout << "DEBUG: inside Conf::write ..."<<endl;
 #endif

  static int conf_count=0;

  long int r;
  ofstream file;
  char conf_name[50], aux[10];

  strcpy(conf_name, QUOTEME(CONF_FILE));
  sprintf(aux, "_%d", conf_count);
  strcat(conf_name, aux);

  file.open(conf_name, ios::out);
  file.precision(16);

  file << nx << " " << ny << " " << nz << " " << nt << " " << beta << " " << mass << " " <<no_flavours << " " <<  update_iteration <<"\n";

  for(r=0; r<no_links; r++)
    {
      file << u_save[r];
    }

  file.close();

  conf_count=1-conf_count;

 #ifdef DEBUG_MODE
  cout << "\tterminated Conf::write"<<endl;
 #endif
}


// write configuration to conf_file (defined in Include/global_macro.cc) 
void Conf::write_last(void)
{
 #ifdef DEBUG_MODE
  cout << "DEBUG: inside Conf::write_last ..."<<endl;
 #endif

  long int r;
  ofstream file;

  file.open(QUOTEME(CONF_FILE), ios::out);
  file.precision(16);

  file << nx << " " << ny << " " << nz << " " << nt << " " << beta << " " << mass << " " <<no_flavours << " " <<  update_iteration <<"\n";

  for(r=0; r<no_links; r++)
    {
      file << u_save[r];
    }

  file.close();

 #ifdef DEBUG_MODE
  cout << "\tterminated Conf::write_last"<<endl;
 #endif
}

// print configuration 
void Conf::print(void)
{
 #ifdef DEBUG_MODE
  cout << "DEBUG: inside Conf::print ..."<<endl;
 #endif

  long int r;

  cout<<"U_SAVE\n";
  for(r=0; r<no_links; r++)
    {
      cout << u_save[r];
    }

  cout<<"U_WORK\n";
  for(r=0; r<no_links; r++)
    {
      cout << u_work[r];
    }

 #ifdef DEBUG_MODE
  cout << "\tterminated Conf::print"<<endl;
 #endif
}


#endif
