
#ifndef VEC3_CC_
#define VEC3_CC_

#include <cmath>
#include <complex>

#include "../Gauss/gauss.cc"
#include "../Rand/random.cc"


using namespace std;

//#include "../Su3/su3.cc"
//#include "../Fermions/fermions.cc"

class Su3;
//class Fermion;
//class MultiFermion;
//class ShiftMultiFermion;


class Vec3 {
 

public:

 complex<REAL> comp[3];
 Vec3(void);
 Vec3(complex<REAL> input[3]);
 Vec3(const Vec3&);


 Vec3 operator~(void) const;
 Vec3& operator=(const Vec3 &rhs);
 void operator+=(const Vec3 &rhs);
 void operator-=(const Vec3 &rhs);
 void operator*=(float rhs);
 void operator*=(double rhs);
 void operator*=(complex<REAL> rhs);


 void zero(void);
 void unitarize(void);
 REAL l2norm(void);
 REAL l2norm2(void);

 void gauss(void);
  void gauss(int);
  void grow_vec3(double site);

 void z2noise(void);
 void ones(void);
 void onesI(void);
 void o123(void);
 void oabc(REAL,REAL,REAL);

 friend ostream& operator<<(ostream &os, const Vec3 &M);
 friend istream& operator>>(istream &is, Vec3 &M);
 friend Vec3 operator+(Vec3 lhs, Vec3 rhs);
 friend Vec3 operator-(Vec3 lhs, Vec3 rhs);
 friend Vec3 operator*(Vec3 lhs, float p);
 friend Vec3 operator*(float p, Vec3 rhs);
 friend Vec3 operator*(Vec3 lhs, double p);
 friend Vec3 operator*(double p, Vec3 rhs);
 friend Vec3 operator*(Vec3 lhs, complex<REAL> p);
 friend Vec3 operator*(complex<REAL> p, Vec3 rhs);
 friend complex<REAL> c_scalprod(Vec3 lhs, Vec3 rhs);
 friend REAL r_scalprod(Vec3 lhs, Vec3 rhs);

 // defined in Su3/su3.cc
 friend Vec3 operator*(const Su3 &mat, const Vec3 &vec);
 friend Su3 operator^(const Vec3 &lhs, const Vec3 &rhs);

 };


				// MEMBER FUNCTIONS
// base constructor
inline Vec3::Vec3(void)
  {
  comp[0]=complex<REAL>(0.0,0.0);
  comp[1]=complex<REAL>(0.0,0.0);
  comp[2]=complex<REAL>(0.0,0.0);
  }


// constructo from input
inline Vec3::Vec3(complex<REAL> input[3])
  {
  comp[0]=input[0];
  comp[1]=input[1];
  comp[2]=input[2];
  }
//copy constructor
inline Vec3::Vec3(const Vec3& input)
  {
  comp[0]=input.comp[0];
  comp[1]=input.comp[1];
  comp[2]=input.comp[2];
  }



// complex conjugate
inline Vec3 Vec3::operator~(void) const
  {
  register Vec3 aux;
  aux.comp[0]=conj(comp[0]); 
  aux.comp[1]=conj(comp[1]); 
  aux.comp[2]=conj(comp[2]); 

  return aux;
  }


// assignement operator
inline Vec3& Vec3::operator=(const Vec3 &rhs)
  {
  comp[0]=rhs.comp[0];
  comp[1]=rhs.comp[1];
  comp[2]=rhs.comp[2];

  return *this;
  }


// operator +=
inline void Vec3::operator+=(const Vec3 &rhs)
  {
  comp[0]+=rhs.comp[0];
  comp[1]+=rhs.comp[1];
  comp[2]+=rhs.comp[2];
  }


// operator -=
inline void Vec3::operator-=(const Vec3 &rhs)
  {
  comp[0]-=rhs.comp[0];
  comp[1]-=rhs.comp[1];
  comp[2]-=rhs.comp[2];
  }


// operator *= with float
inline void Vec3::operator*=(float p)
  {
  comp[0]*=p;
  comp[1]*=p;
  comp[2]*=p;
  }


// operator *= with double
inline void Vec3::operator*=(double p)
  {
  comp[0]*=p;
  comp[1]*=p;
  comp[2]*=p;
  }


// operator *= with complex
inline void Vec3::operator*=(complex<REAL> p)
  {
  comp[0]*=p;
  comp[1]*=p;
  comp[2]*=p;
  }


// zero vector
inline void Vec3::zero(void)
  {
  comp[0]=complex<REAL>(0.0,0.0);
  comp[1]=complex<REAL>(0.0,0.0);
  comp[2]=complex<REAL>(0.0,0.0);
  }


// unitarize
void Vec3::unitarize(void)
  {
  REAL norm=0.0;

  norm+=real(comp[0])*real(comp[0])+imag(comp[0])*imag(comp[0]);
  norm+=real(comp[1])*real(comp[1])+imag(comp[1])*imag(comp[1]);
  norm+=real(comp[2])*real(comp[2])+imag(comp[2])*imag(comp[2]);

  norm=1.0/sqrt(norm);
  comp[0]*=norm;
  comp[1]*=norm;
  comp[2]*=norm;
  }


// return L2 norm 
REAL Vec3::l2norm(void)
  {
  REAL norm=0.0;

  norm+=real(comp[0])*real(comp[0])+imag(comp[0])*imag(comp[0]);
  norm+=real(comp[1])*real(comp[1])+imag(comp[1])*imag(comp[1]);
  norm+=real(comp[2])*real(comp[2])+imag(comp[2])*imag(comp[2]);

  return sqrt(norm);
  }


// return L2 norm squared 
REAL Vec3::l2norm2(void)
  {
  REAL norm2=0.0;

  norm2+=real(comp[0])*real(comp[0])+imag(comp[0])*imag(comp[0]);
  norm2+=real(comp[1])*real(comp[1])+imag(comp[1])*imag(comp[1]);
  norm2+=real(comp[2])*real(comp[2])+imag(comp[2])*imag(comp[2]);

  return norm2;
  }


					// FRIEND FUNCTIONS
// print a vector on a stream
ostream& operator<<(ostream &os, const Vec3 &v)
  {
  os << v.comp[0] << " " << v.comp[1] << " " << v.comp[2] << "\n";
  return os;
  }


// copy a vector froma a stream
istream& operator>>(istream &is, Vec3 &v)
  {
  is >> v.comp[0] >> v.comp[1] >> v.comp[2];
  return is;
  }


// operator +
inline Vec3 operator+(Vec3 lhs, Vec3 rhs)
  {
  lhs+=rhs;
  return lhs;
  }


// operator -
inline Vec3 operator-(Vec3 lhs, Vec3 rhs)
  {
  lhs-=rhs;
  return lhs;
  }


// operator * with a float
inline Vec3 operator*(Vec3 lhs, float p)
  {
  lhs*=p;
  return lhs;
  }


// operator * with a float
inline Vec3 operator*(float p, Vec3 rhs)
  {
  rhs*=p;
  return rhs;
  }


// operator * with a double
inline Vec3 operator*(Vec3 lhs, double p)
  {
  lhs*=p;
  return lhs;
  }


// operator * with a double
inline Vec3 operator*(double p, Vec3 rhs)
  {
  rhs*=p;
  return rhs;
  }


// operator * with a complex
inline Vec3 operator*(Vec3 lhs, complex<REAL> p)
  {
  lhs*=p;
  return lhs;
  }


// operator * with a complex
inline Vec3 operator*(complex<REAL> p, Vec3 rhs)
  {
  rhs*=p;
  return rhs;
  }


// scalar product of two vectors
inline complex<REAL> c_scalprod(Vec3 lhs, Vec3 rhs)
  {
  complex<REAL> ris=complex<REAL>(0.0,0.0);

  ris+=conj(lhs.comp[0])*rhs.comp[0];
  ris+=conj(lhs.comp[1])*rhs.comp[1];
  ris+=conj(lhs.comp[2])*rhs.comp[2];
  
  return ris;
  }


// scalar product of two vectors, real part
inline REAL r_scalprod(Vec3 lhs, Vec3 rhs)
  {
  complex<REAL> ris=complex<REAL>(0.0,0.0);

  ris+=conj(lhs.comp[0])*rhs.comp[0];
  ris+=conj(lhs.comp[1])*rhs.comp[1];
  ris+=conj(lhs.comp[2])*rhs.comp[2];
  
  return real(ris);
  }

// Complex gaussian vector
void Vec3::gauss(void)
  {
  comp[0]=complex_gauss();
  comp[1]=complex_gauss();
  comp[2]=complex_gauss();
  }

// Complex gaussian vector
void Vec3::grow_vec3(double site)
  {
    int v2_sit = (int) site/6.0;
    int sit = (int) site;
    int sit1=sit+1;
    int sit2=sit+2;
    int sit3=sit+3;
    int sit4=sit+4;
    int sit5=sit+5;

    sit = sit%1000;
    sit1 = sit1%1000;
    sit2 = sit2%1000;
    sit3 = sit3%1000;
    sit4 = sit4%1000;
    sit5 = sit5%1000;

    comp[0]=complex<REAL>(sit,sit1);
    comp[1]=complex<REAL>(sit2,sit3);
    comp[2]=complex<REAL>(sit4,sit5);

    comp[0]=complex<REAL>(v2_sit,0.0);
    comp[1]=complex<REAL>(v2_sit,0.0);
    comp[2]=complex<REAL>(v2_sit,0.0);
  }

void Vec3::gauss(int d)
  {
  comp[0]=0;
  comp[1]=2;
  comp[2]=0;
  comp[d]=complex_gauss();
  
  }





// Z2 noise vector
void Vec3::z2noise(void)
  {
  double p;

  p=casuale();
  if(p<0.5) comp[0]=complex<REAL>(1.0,0.0);
  else  comp[0]=complex<REAL>(-1.0,0.0);

  p=casuale();
  if(p<0.5) comp[1]=complex<REAL>(1.0,0.0);
  else  comp[1]=complex<REAL>(-1.0,0.0);

  p=casuale();
  if(p<0.5) comp[2]=complex<REAL>(1.0,0.0);
  else  comp[2]=complex<REAL>(-1.0,0.0);
  }


void Vec3::ones(void)
  {
  comp[0]=1;
  comp[1]=1;
  comp[2]=1;
  }

void Vec3::onesI(void)
  {
  comp[0]=complex<REAL>(1,1);
  comp[1]=complex<REAL>(1,1);
  comp[2]=complex<REAL>(1,1);
  }


void Vec3::o123(void)
  {
  comp[0]=1;
  comp[1]=2;
  comp[2]=3;
  }

void Vec3::oabc(REAL a, REAL b, REAL c)
  {
  comp[0]=a;
  comp[1]=b;
  comp[2]=c;
}


#endif
