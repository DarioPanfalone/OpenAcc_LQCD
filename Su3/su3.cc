
#ifndef SU3_CC_
#define SU3_CC_

#include <cmath>
#include <complex>

#include "../Vec3/vec3.cc"
#include "../Gauss/gauss.cc"
#include "../Rand/random.cc"

using namespace std;


class Su3 {
public:
 complex<REAL> comp[3][3];
 Su3(void);
 Su3(complex<REAL> input[3][3]);

 Su3 operator~(void) const; 
 Su3& operator=(const Su3 &rhs);
 void operator+=(const Su3 &rhs);
 void operator-=(const Su3 &rhs);
 void operator*=(int rhs);
 void operator*=(float rhs);
 void operator*=(double rhs);
 void operator*=(complex<REAL> rhs);
 void operator*=(const Su3 &rhs);

 void one(void);
 void zero(void);
 void isigmaX(int sub);
 void isigmaY(int sub);
 void isigmaZ(int sub);
 void rand_matrix(void);

 REAL retr(void) const;
 REAL imtr(void) const;
 complex<REAL> tr(void) const;
 complex<REAL> det(void) const;
 void sunitarize(void);
 void sub_su2unitarize(int index,int sito);
  void sub_su2unitarize_persmearing(int index,int sito);

 REAL l2norm(void) const;
 REAL l2norm2(void) const;
 void ta(void);
 void exp(void);
 void row_multiply(int i, int p);
 void row_multiply(int i, float p);
 void row_multiply(int i, double p);
 void row_multiply(int i, complex<REAL> p);

 void gauss(void); //produces a random matrix

 friend ostream& operator<<(ostream &os, const Su3 &M);
 friend istream& operator>>(istream &is, Su3 &M);
 friend Su3 operator+(Su3 lhs, Su3 rhs);
 friend Su3 operator-(Su3 lhs, Su3 rhs);
 friend Su3 operator*(Su3 lhs, float p);
 friend Su3 operator*(float p, Su3 rhs);
 friend Su3 operator*(Su3 lhs, double p);
 friend Su3 operator*(double p, Su3 rhs);
 friend Su3 operator*(Su3 lhs, complex<REAL> p);
 friend Su3 operator*(complex<REAL> p, Su3 rhs);
 friend Su3 operator*(Su3 lhs, Su3 rhs);

 friend Vec3 operator*(const Su3 &mat, const Vec3 &vec);
 friend Su3 operator^(const Vec3 &lhs, const Vec3 &rhs);

};
 
                                // MEMBER FUNCTIONS
// base constructor
inline Su3::Su3()
  {
  comp[0][0]=complex<REAL>(0.0,0.0);
  comp[0][1]=complex<REAL>(0.0,0.0);
  comp[0][2]=complex<REAL>(0.0,0.0);

  comp[1][0]=complex<REAL>(0.0,0.0);
  comp[1][1]=complex<REAL>(0.0,0.0);
  comp[1][2]=complex<REAL>(0.0,0.0);

  comp[2][0]=complex<REAL>(0.0,0.0);
  comp[2][1]=complex<REAL>(0.0,0.0);
  comp[2][2]=complex<REAL>(0.0,0.0);
  }


// costructor from imput
inline Su3::Su3(complex<REAL> input[3][3])
  {
  comp[0][0]=input[0][0];
  comp[0][1]=input[0][1];
  comp[0][2]=input[0][2];

  comp[1][0]=input[1][0];
  comp[1][1]=input[1][1];
  comp[1][2]=input[1][2];

  comp[2][0]=input[2][0];
  comp[2][1]=input[2][1];
  comp[2][2]=input[2][2];
  }


// hermitian conjugate
inline Su3 Su3::operator~(void) const
  {
  register Su3 ris;

  ris.comp[0][0]=conj(comp[0][0]);
  ris.comp[0][1]=conj(comp[1][0]);
  ris.comp[0][2]=conj(comp[2][0]);

  ris.comp[1][0]=conj(comp[0][1]);
  ris.comp[1][1]=conj(comp[1][1]);
  ris.comp[1][2]=conj(comp[2][1]);

  ris.comp[2][0]=conj(comp[0][2]);
  ris.comp[2][1]=conj(comp[1][2]);
  ris.comp[2][2]=conj(comp[2][2]);

  return ris;
  }


// assignement operator
inline Su3& Su3::operator=(const Su3 &rhs)
  {
  comp[0][0]=rhs.comp[0][0];
  comp[0][1]=rhs.comp[0][1];
  comp[0][2]=rhs.comp[0][2];

  comp[1][0]=rhs.comp[1][0];
  comp[1][1]=rhs.comp[1][1];
  comp[1][2]=rhs.comp[1][2];

  comp[2][0]=rhs.comp[2][0];
  comp[2][1]=rhs.comp[2][1];
  comp[2][2]=rhs.comp[2][2];

  return *this;
  }


// operator += between matrices
inline void Su3::operator+=(const Su3 &rhs)
  {
  comp[0][0]+=rhs.comp[0][0];
  comp[0][1]+=rhs.comp[0][1];
  comp[0][2]+=rhs.comp[0][2];

  comp[1][0]+=rhs.comp[1][0];
  comp[1][1]+=rhs.comp[1][1];
  comp[1][2]+=rhs.comp[1][2];

  comp[2][0]+=rhs.comp[2][0];
  comp[2][1]+=rhs.comp[2][1];
  comp[2][2]+=rhs.comp[2][2];
  }


// operator -= between matrices
inline void Su3::operator-=(const Su3 &rhs)
  {
  comp[0][0]-=rhs.comp[0][0];
  comp[0][1]-=rhs.comp[0][1];
  comp[0][2]-=rhs.comp[0][2];

  comp[1][0]-=rhs.comp[1][0];
  comp[1][1]-=rhs.comp[1][1];
  comp[1][2]-=rhs.comp[1][2];

  comp[2][0]-=rhs.comp[2][0];
  comp[2][1]-=rhs.comp[2][1];
  comp[2][2]-=rhs.comp[2][2];
  }


// operator *= with a int
inline void Su3::operator*=(int rhs)
  {
  comp[0][0]*=rhs;
  comp[0][1]*=rhs;
  comp[0][2]*=rhs;

  comp[1][0]*=rhs;
  comp[1][1]*=rhs;
  comp[1][2]*=rhs;

  comp[2][0]*=rhs;
  comp[2][1]*=rhs;
  comp[2][2]*=rhs;
  }


// operator *= with a float
inline void Su3::operator*=(float rhs)
  {
  comp[0][0]*=rhs;
  comp[0][1]*=rhs;
  comp[0][2]*=rhs;

  comp[1][0]*=rhs;
  comp[1][1]*=rhs;
  comp[1][2]*=rhs;

  comp[2][0]*=rhs;
  comp[2][1]*=rhs;
  comp[2][2]*=rhs;
  }


// operator *= with a double
inline void Su3::operator*=(double rhs)
  {
  comp[0][0]*=rhs;
  comp[0][1]*=rhs;
  comp[0][2]*=rhs;

  comp[1][0]*=rhs;
  comp[1][1]*=rhs;
  comp[1][2]*=rhs;

  comp[2][0]*=rhs;
  comp[2][1]*=rhs;
  comp[2][2]*=rhs;
  }


// operator *= with a complex
inline void Su3::operator*=(complex<REAL> rhs)
  {
  comp[0][0]*=rhs;
  comp[0][1]*=rhs;
  comp[0][2]*=rhs;

  comp[1][0]*=rhs;
  comp[1][1]*=rhs;
  comp[1][2]*=rhs;

  comp[2][0]*=rhs;
  comp[2][1]*=rhs;
  comp[2][2]*=rhs;
  }


// operator *= with a matrix
inline void Su3::operator*=(const Su3 &rhs)
  {
  register Su3 aux;

  aux=*this;

  register complex<REAL>p0=rhs.comp[0][0];
  register complex<REAL>p1=rhs.comp[1][0];
  register complex<REAL>p2=rhs.comp[2][0];

  comp[0][0]=p0*aux.comp[0][0];
  comp[1][0]=p0*aux.comp[1][0];
  comp[2][0]=p0*aux.comp[2][0];

  comp[0][0]+=p1*aux.comp[0][1];
  comp[1][0]+=p1*aux.comp[1][1];
  comp[2][0]+=p1*aux.comp[2][1];

  comp[0][0]+=p2*aux.comp[0][2];
  comp[1][0]+=p2*aux.comp[1][2];
  comp[2][0]+=p2*aux.comp[2][2];

  p0=rhs.comp[0][1];
  p1=rhs.comp[1][1];
  p2=rhs.comp[2][1];

  comp[0][1]=p0*aux.comp[0][0];
  comp[1][1]=p0*aux.comp[1][0];
  comp[2][1]=p0*aux.comp[2][0];

  comp[0][1]+=p1*aux.comp[0][1];
  comp[1][1]+=p1*aux.comp[1][1];
  comp[2][1]+=p1*aux.comp[2][1];

  comp[0][1]+=p2*aux.comp[0][2];
  comp[1][1]+=p2*aux.comp[1][2];
  comp[2][1]+=p2*aux.comp[2][2];

  p0=rhs.comp[0][2];
  p1=rhs.comp[1][2];
  p2=rhs.comp[2][2];

  comp[0][2]=p0*aux.comp[0][0];
  comp[1][2]=p0*aux.comp[1][0];
  comp[2][2]=p0*aux.comp[2][0];

  comp[0][2]+=p1*aux.comp[0][1];
  comp[1][2]+=p1*aux.comp[1][1];
  comp[2][2]+=p1*aux.comp[2][1];

  comp[0][2]+=p2*aux.comp[0][2];
  comp[1][2]+=p2*aux.comp[1][2];
  comp[2][2]+=p2*aux.comp[2][2];
  }


// zero operator
inline void Su3::zero(void)
  {
  comp[0][0]=complex<REAL>(0.0,0.0);
  comp[0][1]=complex<REAL>(0.0,0.0);
  comp[0][2]=complex<REAL>(0.0,0.0);

  comp[1][0]=complex<REAL>(0.0,0.0);
  comp[1][1]=complex<REAL>(0.0,0.0);
  comp[1][2]=complex<REAL>(0.0,0.0);

  comp[2][0]=complex<REAL>(0.0,0.0);
  comp[2][1]=complex<REAL>(0.0,0.0);
  comp[2][2]=complex<REAL>(0.0,0.0);
  }


// identity operator
inline void Su3::one(void)
  {
  comp[0][0]=complex<REAL>(1.0,0.0);
  comp[0][1]=complex<REAL>(0.0,0.0);
  comp[0][2]=complex<REAL>(0.0,0.0);

  comp[1][0]=complex<REAL>(0.0,0.0);
  comp[1][1]=complex<REAL>(1.0,0.0);
  comp[1][2]=complex<REAL>(0.0,0.0);

  comp[2][0]=complex<REAL>(0.0,0.0);
  comp[2][1]=complex<REAL>(0.0,0.0);
  comp[2][2]=complex<REAL>(1.0,0.0);
  }


// random SU(3) matrix assignement
void Su3::rand_matrix(void)
  {
  int i, j;
  complex<REAL> prod;
  REAL norm;
  Su3 ris;
  double pre,pim;

  for(i=0; i<2; i++)
     {
     for(j=0; j<3; j++)
        {
	  pre = 1.0-2.0*casuale();
	  pim = 1.0-2.0*casuale();
	  ris.comp[i][j]=complex<REAL>(pre,pim);
        }
     }
 
  norm=0.0;
  for(i=0; i<3; i++)
     {
     norm+=real(ris.comp[0][i])*real(ris.comp[0][i])+imag(ris.comp[0][i])*imag(ris.comp[0][i]);
     }
  norm=1.0/sqrt(norm);
  for(i=0; i<3; i++)
     {
     ris.comp[0][i]*=norm;
     }
  
  prod=complex<REAL>(0.0,0.0);
  for(i=0; i<3; i++)
     {
     prod+=conj(ris.comp[0][i])*ris.comp[1][i];
     }
  for(i=0; i<3; i++)
     {
     ris.comp[1][i]-=prod*ris.comp[0][i];
     }
  norm=0.0;
  for(i=0; i<3; i++)
     {
     norm+=real(ris.comp[1][i])*real(ris.comp[1][i])+imag(ris.comp[1][i])*imag(ris.comp[1][i]);
     }
  norm=1.0/sqrt(norm);
  for(i=0; i<3; i++)
     {
     ris.comp[1][i]*=norm;
     }

  prod=ris.comp[0][1]*ris.comp[1][2]-ris.comp[0][2]*ris.comp[1][1];
  ris.comp[2][0]=conj(prod);
  prod=ris.comp[0][2]*ris.comp[1][0]-ris.comp[0][0]*ris.comp[1][2];
  ris.comp[2][1]=conj(prod);
  prod=ris.comp[0][0]*ris.comp[1][1]-ris.comp[0][1]*ris.comp[1][0];
  ris.comp[2][2]=conj(prod);

  *this=ris;
  }


// real part of trace
inline REAL Su3::retr(void) const 
  {
  complex<REAL> aux;

  aux=comp[0][0]+comp[1][1]+comp[2][2];

  return real(aux);
  }
  

// immaginary part of trace
inline REAL Su3::imtr(void) const
  {
  complex<REAL> aux;

  aux=comp[0][0]+comp[1][1]+comp[2][2];

  return imag(aux);
  }


// trace
inline complex<REAL> Su3::tr(void) const
  {
  complex<REAL> aux;

  aux=comp[0][0]+comp[1][1]+comp[2][2];

  return aux;
  }


// determinant
inline complex<REAL> Su3::det(void) const
  {
  complex<REAL> aux;

  aux=comp[0][0]*comp[1][1]*comp[2][2]+comp[1][0]*comp[2][1]*comp[0][2]+comp[2][0]*comp[0][1]*comp[1][2];
  aux-=comp[2][0]*comp[1][1]*comp[0][2]+comp[1][0]*comp[0][1]*comp[2][2]+comp[0][0]*comp[2][1]*comp[1][2];

  return aux;
  }


// project on SU(3)
void Su3::sunitarize(void)
  {
  int i;
  complex<REAL> prod;
  REAL norm;
 
  norm=0.0;
  for(i=0; i<3; i++)
     {
     norm+=real(comp[0][i])*real(comp[0][i])+imag(comp[0][i])*imag(comp[0][i]);
     }
  norm=1.0/sqrt(norm);
  for(i=0; i<3; i++)
     {
     comp[0][i]*=norm;
     }
  
  prod=complex<REAL>(0.0,0.0);
  for(i=0; i<3; i++)
     {
     prod+=conj(comp[0][i])*comp[1][i];
     }
  for(i=0; i<3; i++)
     {
     comp[1][i]-=prod*comp[0][i];
     }
  norm=0.0;
  for(i=0; i<3; i++)
     {
     norm+=real(comp[1][i])*real(comp[1][i])+imag(comp[1][i])*imag(comp[1][i]);
     }
  norm=1.0/sqrt(norm);
  for(i=0; i<3; i++)
     {
     comp[1][i]*=norm;
     }


  prod=comp[0][1]*comp[1][2]-comp[0][2]*comp[1][1];
  comp[2][0]=conj(prod);
  prod=comp[0][2]*comp[1][0]-comp[0][0]*comp[1][2];
  comp[2][1]=conj(prod);
  prod=comp[0][0]*comp[1][1]-comp[0][1]*comp[1][0];
  comp[2][2]=conj(prod);
  }


// L2 norm of the matrix
REAL Su3::l2norm(void) const 
  {
  REAL aux;

  aux=0.0;
  aux+=real(comp[0][0])*real(comp[0][0])+imag(comp[0][0])*imag(comp[0][0]);
  aux+=real(comp[0][1])*real(comp[0][1])+imag(comp[0][1])*imag(comp[0][1]);
  aux+=real(comp[0][2])*real(comp[0][2])+imag(comp[0][2])*imag(comp[0][2]);

  aux+=real(comp[1][0])*real(comp[1][0])+imag(comp[1][0])*imag(comp[1][0]);
  aux+=real(comp[1][1])*real(comp[1][1])+imag(comp[1][1])*imag(comp[1][1]);
  aux+=real(comp[1][2])*real(comp[1][2])+imag(comp[1][2])*imag(comp[1][2]);

  aux+=real(comp[2][0])*real(comp[2][0])+imag(comp[2][0])*imag(comp[2][0]);
  aux+=real(comp[2][1])*real(comp[2][1])+imag(comp[2][1])*imag(comp[2][1]);
  aux+=real(comp[2][2])*real(comp[2][2])+imag(comp[2][2])*imag(comp[2][2]);

  return sqrt(aux);
  }


// L2 norm squared of the matrix
REAL Su3::l2norm2(void) const 
  {
  REAL aux;

  aux=0.0;
  aux+=real(comp[0][0])*real(comp[0][0])+imag(comp[0][0])*imag(comp[0][0]);
  aux+=real(comp[0][1])*real(comp[0][1])+imag(comp[0][1])*imag(comp[0][1]);
  aux+=real(comp[0][2])*real(comp[0][2])+imag(comp[0][2])*imag(comp[0][2]);

  aux+=real(comp[1][0])*real(comp[1][0])+imag(comp[1][0])*imag(comp[1][0]);
  aux+=real(comp[1][1])*real(comp[1][1])+imag(comp[1][1])*imag(comp[1][1]);
  aux+=real(comp[1][2])*real(comp[1][2])+imag(comp[1][2])*imag(comp[1][2]);

  aux+=real(comp[2][0])*real(comp[2][0])+imag(comp[2][0])*imag(comp[2][0]);
  aux+=real(comp[2][1])*real(comp[2][1])+imag(comp[2][1])*imag(comp[2][1]);
  aux+=real(comp[2][2])*real(comp[2][2])+imag(comp[2][2])*imag(comp[2][2]);

  return aux;
  }



// traceless anti-hermitian part
void Su3::ta(void)
  {
  int i, j;
  complex<REAL> trace=complex<REAL>(0.0,0.0);
  Su3 aux;

  aux=*this;

  for(i=0; i<3; i++)
     {
     for(j=0; j<3; j++)
        {
        comp[i][j]=(aux.comp[i][j]-conj(aux.comp[j][i]))*complex<REAL>(0.5,0.0);
        }
     trace+=comp[i][i];
     }

  trace*=one_by_three;
  for(i=0; i<3; i++)
     {
     comp[i][i]-=trace;
     }
  }


// exponential
void Su3::exp(void)
  {
    Su3 aux, ris,uno;
    uno.one();
    // exp x = 1+x*(1+x/2*(1+x/3*(1+x/4*(1+x/5))))

    aux=*this;
    aux*=0.2;
    aux+=uno;
    aux*=*this;
    aux*=0.25;
    aux+=uno;
    aux*=*this;
    aux*=one_by_three;
    aux+=uno;
    aux*=*this;
    aux*=0.5;
    aux+=uno;
    aux*=*this;
    aux+=uno;
    *this=aux;
 

  /*
  ris.one();
  aux=*this;

  ris+=aux;  // ris=1+x
  
  aux*=*this;
  aux*=0.5;
  ris+=aux;  // 2 order

  aux*=*this;
  aux*=one_by_three;
  ris+=aux;  // 3 order

  aux*=*this;
  aux*=0.25;
  ris+=aux;  // 4 order

  aux*=*this;
  aux*=0.2;
  ris+=aux;  // 5 order
  
  aux*=*this;
  aux*=0.166666666666666666;
  ris+=aux;  // 6 order

  ris.sunitarize();
  *this=ris;
  */

  }


// multiply i-th row by p
inline void Su3::row_multiply(int i, int p)
 {
 comp[i][0]*=p;
 comp[i][1]*=p;
 comp[i][2]*=p;
 }


// multiply i-th row by p
inline void Su3::row_multiply(int i, float p)
 {
 comp[i][0]*=p;
 comp[i][1]*=p;
 comp[i][2]*=p;
 }


// multiply i-th row by p
inline void Su3::row_multiply(int i, double p)
 {
 comp[i][0]*=p;
 comp[i][1]*=p;
 comp[i][2]*=p;
 }


// multiply i-th row by p
inline void Su3::row_multiply(int i, complex<REAL> p)
 {
 comp[i][0]*=p;
 comp[i][1]*=p;
 comp[i][2]*=p;
 }


					// FRIEND FUNCTIONS

// print a matrix on ostream
ostream& operator<<(ostream &os, const Su3 &M)
  {
  os << M.comp[0][0] << " " << M.comp[0][1] << " " << M.comp[0][2] <<"\n";
  os << M.comp[1][0] << " " << M.comp[1][1] << " " << M.comp[1][2] <<"\n";
  os << M.comp[2][0] << " " << M.comp[2][1] << " " << M.comp[2][2] <<"\n\n";

  return os;
  }


// copy matrix from istream
istream& operator>>(istream &is, Su3 &M)
  {
  is >> M.comp[0][0] >> M.comp[0][1] >> M.comp[0][2];
  is >> M.comp[1][0] >> M.comp[1][1] >> M.comp[1][2];
  is >> M.comp[2][0] >> M.comp[2][1] >> M.comp[2][2];
  
  return is;
  }


// sum of matrices
inline Su3 operator+(Su3 lhs, Su3 rhs)
  {
  lhs+=rhs;
  return lhs;
  }


// difference of matrices
inline Su3 operator-(Su3 lhs, Su3 rhs)
  {
  lhs-=rhs;
  return lhs;
  }


// product with float
inline Su3 operator*(Su3 lhs, float p)
  {
  lhs*=p;
  return lhs;
  }


// product with float
inline Su3 operator*(float p, Su3 rhs)
  {
  rhs*=p;
  return rhs;
  }


// product with double
inline Su3 operator*(Su3 lhs, double p)
  {
  lhs*=p;
  return lhs;
  }


// product with double
inline Su3 operator*(double p, Su3 rhs)
  {
  rhs*=p;
  return rhs;
  }


// product with complex
inline Su3 operator*(Su3 lhs, complex<REAL> p)
  {
  lhs*=p;
  return lhs;
  }


// product with complex
inline Su3 operator*(complex<REAL> p, Su3 rhs)
  {
  rhs*=p;
  return rhs;
  }


// product between matrices
Su3 operator*(Su3 lhs, Su3 rhs)
  {
  lhs*=rhs;
  return lhs;
  }


// matrix * vector product
inline Vec3 operator*(const Su3 &mat, const Vec3 &vec)
  {
  register Vec3 ris;
  register complex<REAL> p0=vec.comp[0];
  register complex<REAL> p1=vec.comp[1];
  register complex<REAL> p2=vec.comp[2];
 
  ris.comp[0]=p0*mat.comp[0][0];
  ris.comp[1]=p0*mat.comp[1][0];
  ris.comp[2]=p0*mat.comp[2][0];

  ris.comp[0]+=p1*mat.comp[0][1];
  ris.comp[1]+=p1*mat.comp[1][1];
  ris.comp[2]+=p1*mat.comp[2][1];

  ris.comp[0]+=p2*mat.comp[0][2];
  ris.comp[1]+=p2*mat.comp[1][2];
  ris.comp[2]+=p2*mat.comp[2][2];
 
  return ris;
  }


// direct product of vectors
inline Su3 operator^(const Vec3 &lhs, const Vec3 &rhs)
  {
  Su3 ris;

  ris.comp[0][0]=lhs.comp[0]*rhs.comp[0];
  ris.comp[0][1]=lhs.comp[0]*rhs.comp[1];
  ris.comp[0][2]=lhs.comp[0]*rhs.comp[2];

  ris.comp[1][0]=lhs.comp[1]*rhs.comp[0];
  ris.comp[1][1]=lhs.comp[1]*rhs.comp[1];
  ris.comp[1][2]=lhs.comp[1]*rhs.comp[2];

  ris.comp[2][0]=lhs.comp[2]*rhs.comp[0];
  ris.comp[2][1]=lhs.comp[2]*rhs.comp[1];
  ris.comp[2][2]=lhs.comp[2]*rhs.comp[2];

  return ris;
  }

// proietta la matrice sul sottogruppo su2 identificato da index                                                                                                                                                                                                                         
void Su3::sub_su2unitarize(int index,int sito)
{
  Su3 aux=*this;
  REAL norma=0.0;
  // aux.row_multiply(2,eta[sito]); 
  int ind1,ind2;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      if(i==index) aux.comp[i][j] = 0.0;
      if(j==index){
	aux.comp[i][j] = 0.0;
	if(i==index) aux.comp[i][j] = 1.0;
      }
    }
  }

  ind1=1;
  ind2=2;
  if(index==1){
    ind1=0;
    ind2=2;
  }
  if(index==2){
    ind1=0;
    ind2=1;
  }
  REAL c0,c1,c2,c3;
  c0=(real(aux.comp[ind1][ind1]+aux.comp[ind2][ind2]));
  //  cout << "C0  = " << c0 << endl;
  c1=(imag(aux.comp[ind1][ind2]+aux.comp[ind2][ind1]));
  c2=(real(aux.comp[ind1][ind2]-aux.comp[ind2][ind1]));
  c3=(imag(aux.comp[ind1][ind1]-aux.comp[ind2][ind2]));
  norma=-sqrt(c0*c0+c1*c1+c2*c2+c3*c3);   // minus sign due to staggered phases
  //  cout << "norma  = " << norma << endl;

  c0/=norma;
  c1/=norma;
  c2/=norma;
  c3/=norma;

  aux.comp[ind1][ind1]=  complex<REAL>(c0,c3);
  aux.comp[ind1][ind2]=  complex<REAL>(c2,c1);
  aux.comp[ind2][ind1]= -conj(aux.comp[ind1][ind2]);
  aux.comp[ind2][ind2]=  conj(aux.comp[ind1][ind1]);
  //  aux.row_multiply(2,eta[sito]);
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      comp[i][j]=aux.comp[i][j];
    }
  }


}


// proietta la matrice sul sottogruppo su2 identificato da index 
void Su3::sub_su2unitarize_persmearing(int index,int sito)
{
  Su3 aux=*this;
  REAL norm = 0.0;
  REAL norma=0.0;
  int ind1,ind2;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      if(i==index) aux.comp[i][j] = 0.0;
      if(j==index){
        aux.comp[i][j] = 0.0;
        if(i==index) aux.comp[i][j] = 1.0;
      }
    }
  }
  ind1=1;
  ind2=2;
  if(index==1){
    ind1=0;
    ind2=2;
  }
  if(index==2){
    ind1=0;
    ind2=1;
  }
  REAL c0,c1,c2,c3;
  c0=(real(aux.comp[ind1][ind1]+aux.comp[ind2][ind2]));
  c1=(imag(aux.comp[ind1][ind2]+aux.comp[ind2][ind1]));
  c2=(real(aux.comp[ind1][ind2]-aux.comp[ind2][ind1]));
  c3=(imag(aux.comp[ind1][ind1]-aux.comp[ind2][ind2]));
  norma=sqrt(c0*c0+c1*c1+c2*c2+c3*c3);
  c0/=norma;
  c1/=norma;
  c2/=norma;
  c3/=norma;
  aux.comp[ind1][ind1]=  complex<REAL>(c0,c3);
  aux.comp[ind1][ind2]=  complex<REAL>(c2,c1);
  aux.comp[ind2][ind1]= -conj(aux.comp[ind1][ind2]);
  aux.comp[ind2][ind2]=  conj(aux.comp[ind1][ind1]);
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      comp[i][j]=aux.comp[i][j];
    }
  }
}



void Su3::gauss(void)
  {
  int i;
  REAL rand[8], aux1, aux2;

  for(i=0; i<4; i++)
     {
     real_gauss(aux1, aux2);
     rand[2*i]=aux1;
     rand[2*i+1]=aux2;
     }

  comp[0][0]=complex<REAL>(rand[2] + one_by_sqrt_three*rand[7], 0.0);
  comp[0][1]=complex<REAL>(rand[0], -rand[1]);
  comp[0][2]=complex<REAL>(rand[3], -rand[4]);

  comp[1][0]=complex<REAL>(rand[0], rand[1]);
  comp[1][1]=complex<REAL>(-rand[2] + one_by_sqrt_three*rand[7], 0.0);
  comp[1][2]=complex<REAL>(rand[5], -rand[6]);

  comp[2][0]=complex<REAL>(rand[3], rand[4]);
  comp[2][1]=complex<REAL>(rand[5], rand[6]);
  comp[2][2]=complex<REAL>( -two_by_sqrt_three * rand[7], 0.0);
  }



void Su3::isigmaX(int sub){

    this->zero();
    int i,j;
    switch(sub){
         case 0 : i = 1 ; j = 2 ; break;
         case 1 : i = 0 ; j = 2 ; break;
         case 2 : i = 0 ; j = 1 ; break;
    }
    comp[sub][sub] = 1 ;
    comp[i][j] = complex<REAL>(0,1.0);
    comp[j][i] = complex<REAL>(0,1.0);

}


void Su3::isigmaY(int sub){

    this->zero();
    int i,j;
    switch(sub){
         case 0 : i = 1 ; j = 2 ; break;
         case 1 : i = 0 ; j = 2 ; break;
         case 2 : i = 0 ; j = 1 ; break;
    }
    comp[sub][sub] = 1 ;
    comp[i][j] = complex<REAL>(-1.0,0);
    comp[j][i] = complex<REAL>(1.0,0);

}

void Su3::isigmaZ(int sub){

    this->zero();
    int i,j;
    switch(sub){
         case 0 : i = 1 ; j = 2 ; break;
         case 1 : i = 0 ; j = 2 ; break;
         case 2 : i = 0 ; j = 1 ; break;
    }
    comp[sub][sub] = 1 ;
    comp[i][i] = complex<REAL>(0,1.0);
    comp[j][j] = complex<REAL>(0,-1.0);

}















#endif
