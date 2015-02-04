#ifndef FERMIONS_CC_
#define FERMIONS_CC_


#include "../Vec3/vec3.cc"
#include "../Global_sum/global_sum.cc"
#include "../Include/parameters.cc"
#include "../Include/global_macro.cc"
#include <fstream>
/*
typedef struct{

    complex <double> c0[sizeh];
    complex <double> c1[sizeh];
    complex <double> c2[sizeh];

} Vec3_soa ;
*/

							// SINGLE FERMION
class Fermion {
public:

  Vec3 fermion[sizeh];

  Fermion(void);
  Fermion(const Fermion&); //copy constructor
  
  void gauss(void);
  void gauss(int);
  void z2noise(void);
  void allOnes(void);
  void Growing(void);
  void allOnesI(void);
  void all123(void);
  double l2norm2(void);
  void allabc(REAL,REAL,REAL);
  void allOnesExcept(int);
  
  void ferm_aos_to_soaCOM(vec3COM_soa *out) const;
  void ferm_soaCOM_to_aos(vec3COM_soa const* const in);

  void saveToFile(const char*);
  void loadFromFile(const char*);

  double difference(Fermion a,Fermion b); 
 
};

// base constructor
Fermion::Fermion(void)
 {
 for(long int i=0; i<sizeh; i++)
    {
    fermion[i].zero();
    }
 }

//copy constructor
Fermion::Fermion(const Fermion& frm){

  for(int i = 0 ; i < sizeh ; i++)
     fermion[i] = Vec3(frm.fermion[i]);

}



// Initialize with random gaussian complex numbers
void Fermion::gauss(void)
 {
 for(long int i=0; i<sizeh; i++)
    {
    fermion[i].gauss();
    }
 }

void Fermion::gauss(int d)
 {
 for(long int i=0; i<sizeh; i++)
    {
    fermion[i].gauss(d);
    }
 }

void Fermion::Growing(void)
{
  for(long int i=0; i<sizeh; i++)
    {
      fermion[i].grow_vec3(i*6.0);
    }
}


// Initialize with z2 noise
void Fermion::z2noise(void)
 {
 for(long int i=0; i<sizeh; i++)
    {
    fermion[i].z2noise();
    }
 }

// Initialize with all ones
void Fermion::allOnes(void)
 {
 for(long int i=0; i<sizeh; i++)
    {
    fermion[i].ones();
    }
 }


void Fermion::allOnesI(void)
 {
 for(long int i=0; i<sizeh; i++)
    {
    fermion[i].onesI();
    }
 }



// Initialize with all {1,2,3}
void Fermion::all123(void)
 {
 for(long int i=0; i<sizeh; i++)
    {
    fermion[i].o123();
    }
 }

void Fermion::allabc(REAL a, REAL b, REAL c)
 {
  for(long int i=0; i<sizeh; i++)
    {
    fermion[i].oabc(a,b,c);
    }
}

void Fermion::allOnesExcept(int ii){

 for(long int i=0; i<sizeh; i++)
    {
    fermion[i].ones();
    }

 fermion[ii].comp[0] = 0;
 fermion[ii].comp[1] = 0;
 fermion[ii].comp[2] = 0;


}

// squared L2 norm of the fermion
double Fermion::l2norm2(void)
 {
 for(long int i=0; i<sizeh; i++)
    {
    d_vector1[i]=(fermion[i].l2norm2());
    }
 global_sum(d_vector1, sizeh);

 return d_vector1[0];
 }

void Fermion::saveToFile(const char *filename){

    ofstream ofile;
    ofile.open(filename);
    ofile.precision(16);
     
    for(int i = 0 ; i < sizeh ; i++)
         ofile << fermion[i];

    ofile.close();
}

void Fermion::loadFromFile(const char *filename){

    ifstream ifile;
    ifile.open(filename);
     
    for(int i = 0 ; i < sizeh ; i++)
         ifile >> fermion[i];

    ifile.close();
}


void Fermion::ferm_aos_to_soaCOM(vec3COM_soa *out) const{
    for( int i =0 ; i < sizeh ; i++){
      out->c0[i].Re = (fermion[i].comp[0]).real();
      out->c1[i].Re = (fermion[i].comp[1]).real();
      out->c2[i].Re = (fermion[i].comp[2]).real();
      out->c0[i].Im = (fermion[i].comp[0]).imag();
      out->c1[i].Im = (fermion[i].comp[1]).imag();
      out->c2[i].Im = (fermion[i].comp[2]).imag();
    }
}
 

void Fermion::ferm_soaCOM_to_aos(vec3COM_soa const* const in){
    for( int i =0 ; i < sizeh ; i++){
      fermion[i].comp[0] = complex<double>(in->c0[i].Re,in->c0[i].Im);
      fermion[i].comp[1] = complex<double>(in->c1[i].Re,in->c1[i].Im);
      fermion[i].comp[2] = complex<double>(in->c2[i].Re,in->c2[i].Im);
    }     
}


double difference(Fermion *a,Fermion *b){
  complex<double> cd0,cd1,cd2;
  double d0,d1,d2;
  for( int i =0 ; i < sizeh ; i++){
    cd0 = (a->fermion[i].comp[0])-(b->fermion[i].comp[0]);
    cd1 = (a->fermion[i].comp[1])-(b->fermion[i].comp[1]);
    cd2 = (a->fermion[i].comp[2])-(b->fermion[i].comp[2]);
    d0=cd0.real()*cd0.real()+cd0.imag()*cd0.imag();
    d1=cd1.real()*cd1.real()+cd1.imag()*cd1.imag();
    d2=cd2.real()*cd2.real()+cd2.imag()*cd2.imag();
    d_vector1[i]=d0+d1+d2;
  }
  global_sum(d_vector1, sizeh);
  return sqrt(d_vector1[0]);
}
//SHIFT FERMIONS        
class ShiftFermion {
public:

  Vec3 fermion[max_approx_order][sizeh];
  ShiftFermion(void);

  // defined below       
  friend void extract_fermion(ShiftFermion *out, ShiftMultiFermion *in, int i);

  // defined in Inverter/inverter.cc 
  friend void multips_shifted_invert (ShiftMultiFermion *chi, MultiFermion *phi, REAL res, RationalApprox approx);

  //defined here         
  void ferm_Shift_to_ShiftCOM(COM_ShiftFermion *out) const;
  void ferm_ShiftCOM_to_Shift(COM_ShiftFermion const* const in);

};


// base construction     
ShiftFermion::ShiftFermion(void)
{
  for(int i=0; i<max_approx_order; i++)
    {
      for(long int j=0; j<sizeh; j++)
        {
          fermion[i][j].zero();
        }
    }
}

void ShiftFermion::ferm_Shift_to_ShiftCOM(COM_ShiftFermion *out) const{
  for( int ia =0 ; ia < max_approx_order ; ia++){
    for( int i =0 ; i < sizeh ; i++){
      out->shift[ia].c0[i].Re = (fermion[ia][i].comp[0]).real();
      out->shift[ia].c1[i].Re = (fermion[ia][i].comp[1]).real();
      out->shift[ia].c2[i].Re = (fermion[ia][i].comp[2]).real();
      out->shift[ia].c0[i].Im = (fermion[ia][i].comp[0]).imag();
      out->shift[ia].c1[i].Im = (fermion[ia][i].comp[1]).imag();
      out->shift[ia].c2[i].Im = (fermion[ia][i].comp[2]).imag();
    }
  }
}

void ShiftFermion::ferm_ShiftCOM_to_Shift(COM_ShiftFermion const* const in){
  for( int ia =0 ; ia < max_approx_order ; ia++){
    for( int i =0 ; i < sizeh ; i++){
      fermion[ia][i].comp[0] = complex<double>(in->shift[ia].c0[i].Re,in->shift[ia].c0[i].Im);
      fermion[ia][i].comp[1] = complex<double>(in->shift[ia].c1[i].Re,in->shift[ia].c1[i].Im);
      fermion[ia][i].comp[2] = complex<double>(in->shift[ia].c2[i].Re,in->shift[ia].c2[i].Im);
    }
  }
}


// MULTIPLE FERMIONS      
class MultiFermion {
public:
  Vec3 fermion[no_ps][sizeh];

  MultiFermion(void);

  friend void create_phi(void);
  friend void extract_fermion(Fermion *out, MultiFermion *in, int i);

  // defined in Action/action.cc     
  friend void fermion_action(double *value, int init);

  // defined in FermionForce/fermionforce.cc     
  friend void fermionforce(void);

  // defined in Inverter/inverter.cc 
  friend void multips_shifted_invert (ShiftMultiFermion *chi, MultiFermion *phi, REAL res, RationalApprox approx);

  // defined in Inverter/cu_inverter.cc          
  friend void cu_multips_shifted_invert (REAL res, RationalApprox approx);

  // defined in RationalApprox/rationalapprox.cc 
  friend void first_inv_approx_calc(REAL res);
  friend void last_inv_approx_calc(REAL res);

  // defined in Packer/packer.cc     
  friend void smartpack_multifermion(float out[6*sizeh*no_ps*2] , const MultiFermion *in);
  friend void smartunpack_multifermion(MultiFermion *out, const float in[6*sizeh*no_ps*2]);

  //defined here         
  void ferm_Multi_to_MultiCOM(COM_MultiFermion *out) const;
  void ferm_MultiCOM_to_Multi(COM_MultiFermion const* const in);

  double difference_multi(MultiFermion *a,MultiFermion *b);

};


// base construction     
MultiFermion::MultiFermion(void)
{
  for(int i=0; i<no_ps; i++)
    {
      for(long int j=0; j<sizeh; j++)
        {
          fermion[i][j].zero();
        }
    }
}

double difference_multi(MultiFermion *a,MultiFermion *b){
  complex<double> cd0,cd1,cd2;
  double d0,d1,d2;
  for(long int j=0; j<sizeh; j++){
    d_vector1[j]=0.0;
    for(int i=0; i<no_ps; i++){
      
      cd0 = (a->fermion[i][j].comp[0])-(b->fermion[i][j].comp[0]);
      cd1 = (a->fermion[i][j].comp[1])-(b->fermion[i][j].comp[1]);
      cd2 = (a->fermion[i][j].comp[2])-(b->fermion[i][j].comp[2]);
      d0=cd0.real()*cd0.real()+cd0.imag()*cd0.imag();
      d1=cd1.real()*cd1.real()+cd1.imag()*cd1.imag();
      d2=cd2.real()*cd2.real()+cd2.imag()*cd2.imag();
      d_vector1[j]+=d0+d1+d2;
    }
  }
  global_sum(d_vector1, sizeh);
  return sqrt(d_vector1[0]);
}


void MultiFermion::ferm_Multi_to_MultiCOM(COM_MultiFermion *out) const{
  for( int ips =0 ; ips < no_ps ; ips++){
    for( int i =0 ; i < sizeh ; i++){
      out->multi[ips].c0[i].Re = (fermion[ips][i].comp[0]).real();
      out->multi[ips].c1[i].Re = (fermion[ips][i].comp[1]).real();
      out->multi[ips].c2[i].Re = (fermion[ips][i].comp[2]).real();
      out->multi[ips].c0[i].Im = (fermion[ips][i].comp[0]).imag();
      out->multi[ips].c1[i].Im = (fermion[ips][i].comp[1]).imag();
      out->multi[ips].c2[i].Im = (fermion[ips][i].comp[2]).imag();
    }
  }
}


void MultiFermion::ferm_MultiCOM_to_Multi(COM_MultiFermion const* const in){
  for( int ips =0 ; ips < no_ps ; ips++){
    for( int i =0 ; i < sizeh ; i++){
      fermion[ips][i].comp[0] = complex<double>(in->multi[ips].c0[i].Re,in->multi[ips].c0[i].Im);
      fermion[ips][i].comp[1] = complex<double>(in->multi[ips].c1[i].Re,in->multi[ips].c1[i].Im);
      fermion[ips][i].comp[2] = complex<double>(in->multi[ips].c2[i].Re,in->multi[ips].c2[i].Im);
    }
  }
}


// Initialize fermion_phi with gaussian random numbers       
void create_phi(void)
{
 #ifdef DEBUG_MODE
  cout << "DEBUG: inside create_phi ..."<<endl;
 #endif

  for(int i=0; i<no_ps; i++)
    {
      for(long int j=0; j<sizeh; j++)
        {
          (fermion_phi->fermion[i][j]).gauss();
        }
    }
 #ifdef DEBUG_MODE
  cout << "\tterminated create_phi"<<endl;
 #endif
}


// Extract the i-th fermion from the MultiFermion "in"       
void extract_fermion(Fermion *out, MultiFermion *in, int i)
{
  for(long int j=0; j<sizeh; j++)
    {
      out->fermion[j]=in->fermion[i][j];
    }
}



// SHIFT MULTIPLE FERMIONS 
class ShiftMultiFermion {
public:
  Vec3 fermion[no_ps][max_approx_order][sizeh];

  ShiftMultiFermion(void);

  friend void extract_fermion(ShiftFermion *out, ShiftMultiFermion *in, int i);

  // defined in FermionForce/fermionforce.cc     
  friend void fermionforce(int moltiplico);

  // defined in Inverter/inverter.cc 
  friend void multips_shifted_invert (ShiftMultiFermion *chi, MultiFermion *phi, REAL res, RationalApprox approx);

  // defined in Inverter/cu_inverter.cc          
  friend void cu_multips_shifted_invert (REAL res, RationalApprox approx);

  // defined in RationalApprox/rationalapprox.cc 
  friend void first_inv_approx_calc(REAL res);
  friend void last_inv_approx_calc(REAL res);

  // defined in Packer/packer.cc     
  friend void smartunpack_multishiftfermion(ShiftMultiFermion *out, const float in[6*sizeh*max_approx_order*no_ps*2],  int order);

  //defined here         
  void ferm_ShiftMulti_to_ShiftMultiCOM(COM_ShiftMultiFermion *out) const;
  void ferm_ShiftMultiCOM_to_ShiftMulti(COM_ShiftMultiFermion const* const in);
};

void ShiftMultiFermion::ferm_ShiftMulti_to_ShiftMultiCOM(COM_ShiftMultiFermion *out) const{
  for( int ia =0 ; ia < max_approx_order ; ia++){
    for( int ips =0 ; ips < no_ps ; ips++){
      for( int i =0 ; i < sizeh ; i++){
        out->shiftmulti[ia][ips].c0[i].Re = (fermion[ips][ia][i].comp[0]).real();
        out->shiftmulti[ia][ips].c1[i].Re = (fermion[ips][ia][i].comp[1]).real();
        out->shiftmulti[ia][ips].c2[i].Re = (fermion[ips][ia][i].comp[2]).real();
        out->shiftmulti[ia][ips].c0[i].Im = (fermion[ips][ia][i].comp[0]).imag();
        out->shiftmulti[ia][ips].c1[i].Im = (fermion[ips][ia][i].comp[1]).imag();
        out->shiftmulti[ia][ips].c2[i].Im = (fermion[ips][ia][i].comp[2]).imag();
      }
    }
  }
}

void ShiftMultiFermion::ferm_ShiftMultiCOM_to_ShiftMulti(COM_ShiftMultiFermion const* const in){
  for( int ia =0 ; ia < max_approx_order ; ia++){
    for( int ips =0 ; ips < no_ps ; ips++){
      for( int i =0 ; i < sizeh ; i++){
        fermion[ips][ia][i].comp[0] = complex<double>(in->shiftmulti[ia][ips].c0[i].Re,in->shiftmulti[ia][ips].c0[i].Im);
        fermion[ips][ia][i].comp[1] = complex<double>(in->shiftmulti[ia][ips].c1[i].Re,in->shiftmulti[ia][ips].c1[i].Im);
        fermion[ips][ia][i].comp[2] = complex<double>(in->shiftmulti[ia][ips].c2[i].Re,in->shiftmulti[ia][ips].c2[i].Im);
      }
    }
  }
}


// base construction     
ShiftMultiFermion::ShiftMultiFermion(void)
{
  for(int i=0; i<no_ps; i++)
    {
      for(int k=0; k<max_approx_order; k++)
        {
          for(long int j=0; j<sizeh; j++)
            {
              fermion[i][k][j].zero();
            }
        }
    }
}

// extract the i-th pseudfermion from the ShiftMultiFermion "in"
void extract_fermion(ShiftFermion *out, ShiftMultiFermion *in, int i)
{
  int j;
  long int l;

  for(j=0; j<max_approx_order; j++)
    {
      for(l=0; l>sizeh; l++)
        {
          out->fermion[j][l]=in->fermion[i][j][l];
        }
    }
}

 
#endif
