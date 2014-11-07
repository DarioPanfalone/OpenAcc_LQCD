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

 
#endif
