#ifdef USE_GPU
#include"../Packer/packer.h"
#include"../Cuda/cuda_sol_sum.h"
#include"../Inverter/cu_inverter.h"
#endif


//                                                                    __ i<approx_order       RA_a[i]
//  approximation valid in [min_epsilon, 1] of the form  f(x)=RA_a0 + \                  --------------------
//                                                                    /_ i=0                x + RA_b[i]
class RationalApprox {
  //private:
public:
  REAL min_epsilon;
  int approx_order;
  REAL RA_a0;
  REAL RA_a[max_approx_order];
  REAL RA_b[max_approx_order];

  RationalApprox();
  RationalApprox(const REAL eps, const int order, const REAL a0, const REAL *a, const REAL *b);
  RationalApprox& operator=(const RationalApprox &rhs); 
  
  friend void get_order(int &order, RationalApprox approx);
  friend void get_shifts(float shifts[max_approx_order], RationalApprox approx);
  friend void get_shifts(double shifts[max_approx_order], RationalApprox approx);
  friend void get_const(float &a, RationalApprox approx);
  friend void get_const(double &a, RationalApprox approx);
  friend void get_numerators(float numerators[max_approx_order], RationalApprox approx);
  friend void get_numerators(double numerators[max_approx_order], RationalApprox approx);
  friend std::ostream& operator<<(std::ostream &os, const RationalApprox &approx); 

  friend void convert_RationalApprox_to_COM_RationalApprox(COM_RationalApprox *COMapprox,RationalApprox CPPapprox);
  
  // rescaled coefficients
  void first_inv_approx_coeff(void);
  void md_inv_approx_coeff(void);
  void last_inv_approx_coeff(void);
  
  
  friend void first_inv_approx_calc(REAL res);
  friend void last_inv_approx_calc(REAL res);
  // defined in Inverter/inverter.cc
  friend void multips_shifted_invert(ShiftMultiFermion *chi, MultiFermion *phi, REAL res, RationalApprox approx);
  
  // defined in FermionForce/fermionforce.cc
  friend void fermionforce(int moltiplico);
};


void convert_RationalApprox_to_COM_RationalApprox(COM_RationalApprox *COMapprox,RationalApprox CPPapprox)
{
  COMapprox->COM_min_epsilon=CPPapprox.min_epsilon;
  COMapprox->COM_approx_order=CPPapprox.approx_order;
  COMapprox->COM_RA_a0=CPPapprox.RA_a0;
  for(int i=0;i<max_approx_order;i++){
    COMapprox->COM_RA_a[i]=CPPapprox.RA_a[i];
    COMapprox->COM_RA_b[i]=CPPapprox.RA_b[i];
  }
}


// default constructor
RationalApprox::RationalApprox(void)
 {
 min_epsilon=1.0;
 approx_order=max_approx_order;
 RA_a0=0.0;
 for(int i=0; i<approx_order; i++)
    {
    RA_a[i]=0.0;
    RA_b[i]=0.0;
    }
 }


// constuctor with initialization
RationalApprox::RationalApprox(const REAL eps, const int order, const REAL a0, const REAL *a, const REAL *b)
 {
 min_epsilon=eps;
 approx_order=order;
 RA_a0=a0;
 for(int i=0; i<order; i++)
    {
    RA_a[i]=a[i];
    RA_b[i]=b[i];
    }
 }


// asignement operator
RationalApprox& RationalApprox::operator=(const RationalApprox& rhs)
 {
 min_epsilon=rhs.min_epsilon;
 approx_order=rhs.approx_order;
 RA_a0=rhs.RA_a0;
 for(int i=0; i<approx_order; i++)
    {
    RA_a[i]=rhs.RA_a[i];
    RA_b[i]=rhs.RA_b[i];
    }
 
 return *this;
 }


void get_order(int &order, RationalApprox approx)
 {
 order=(approx.approx_order);
 }


void get_shifts(float shifts[max_approx_order], RationalApprox approx)
 {
 for(int i=0; i<max_approx_order; i++)
    {
    shifts[i]=(float) (approx.RA_b[i]); 
    }
 }


void get_const(float &a, RationalApprox approx)
 {
 a=(approx.RA_a0);
 }


void get_const(double &a, RationalApprox approx)
 {
 a=(approx.RA_a0);
 }


void get_shifts(double shifts[max_approx_order], RationalApprox approx)
 {
 for(int i=0; i<max_approx_order; i++)
    {
    shifts[i]=(double) (approx.RA_b[i]); 
    }
 }


void get_numerators(float numerators[max_approx_order], RationalApprox approx)
 {
 for(int i=0; i<max_approx_order; i++)
    {
    numerators[i]=(float) (approx.RA_a[i]); 
    }
 }


void get_numerators(double numerators[max_approx_order], RationalApprox approx)
 {
 for(int i=0; i<max_approx_order; i++)
    {
    numerators[i]=(double) (approx.RA_a[i]); 
    }
 }

std::ostream& operator<<(std::ostream &os, const RationalApprox &approx) 
   {
   os << "\t+++++++++++RationalApprox++++++++++++++++\n";
   os << "\t++ approx_order="<< approx.approx_order<<"\n";
   os << "\t++ min_epsilon="<< approx.min_epsilon<<"\n";
   os << "\t++ RA_a0="<< approx.RA_a0 <<"\n";
   for(int i=0; i<approx.approx_order; i++)
      {
      os <<"\t++ "<<i<<"\t"<<"RA_a["<<i<<"]="<< approx.RA_a[i] << "\t" <<"RA_b["<<i<<"]="<< approx.RA_b[i]<<"\n";
      }
   os << "\t+++++++++++++++++++++++++++++++++++++++++"<< endl;
   return os;
   }


                                             // RESCALED COEFFICIENTS

void RationalApprox::first_inv_approx_coeff(void)
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside RationalApprox::first_inv_approx_coef ..."<<endl;
 #endif
 int i;
 REAL min, max, epsilon;
 REAL *minmax;
 minmax = new REAL [2];

 cout.precision(18);

 // normalized coefficients
 *this=*first_inv_approx_norm_coeff;

 // findminmax(min, max);
 //  cout << "CPU:      min = " << min << "    max = " << max << endl;
 findminmax_con_openacc(minmax);
 min = minmax[0];
 max = minmax[1];

 cout << "OPENACC:  min = " << min << "    max = " << max << endl;

 min*=0.95;
 max*=1.05;
 cout << "OPENACC times security factors:  min = " << min << "    max = " << max << endl;

 epsilon=min/max;
 cout << "Epsilon = " << epsilon << endl ;

 if(epsilon<min_epsilon)
   {
   ofstream err_file;
   err_file.open(QUOTEME(ERROR_FILE), ios::app);   
   err_file << "WARNING: in first_inv_approx_coef epsilon="<<epsilon<<" < min_epsilon="<<min_epsilon;
   err_file << " at update_iteration="<<update_iteration <<endl;
   err_file.close();
   } 

  
 // rescale coeff.
 double exponent = ((double) no_flavours) /8./ ((double)no_ps);
 epsilon=pow(max, exponent);
 cout << "RA_a0  = " << RA_a0 << endl ;
 RA_a0*=epsilon;
 cout << "RA_a0  = " << RA_a0 << endl ;

 for(i=0; i<approx_order; i++)
    {
    RA_a[i]*=(max*epsilon);
    RA_b[i]*=max;
    }

 #ifdef DEBUG_MODE
 cout << "\tterminated RationalApprox::first_inv_approx_coeff, epsilon="<< min/max<< " min_epsilon="<<min_epsilon<<endl;
 #endif
 }


void RationalApprox::md_inv_approx_coeff(void)
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside RationalApprox::md_inv_approx_coef ..."<<endl;
 #endif
 int i;
 REAL min, max, epsilon;

 REAL *minmax;
 minmax = new REAL [2];
 // normalized coefficients
 *this=*md_inv_approx_norm_coeff;

 // findminmax(min, max);
 findminmax_con_openacc(minmax);
 min = minmax[0];
 max = minmax[1];


 min*=0.95;
 max*=1.05;

 epsilon=min/max;

 if(epsilon<min_epsilon)
   {
   ofstream err_file;
   err_file.open(QUOTEME(ERROR_FILE), ios::app);   
   err_file << "WARNING: in first_md_approx_coef epsilon="<<epsilon<<" < min_epsilon="<<min_epsilon;
   err_file << " at update_iteration="<<update_iteration <<endl;
   err_file.close();
   } 

  
 // rescale coeff.
 double exponent = -((double) no_flavours)/4./((double)no_ps);
 epsilon=pow(max, exponent);
 RA_a0*=epsilon;
 for(i=0; i<approx_order; i++)
    {
    RA_a[i]*=(max*epsilon);
    RA_b[i]*=max;
    }

 #ifdef DEBUG_MODE
 cout << "\tterminated RationalApprox::md_inv_approx_coeff, epsilon="<< min/max<< " min_epsilon="<<min_epsilon<<endl;
 #endif
 }


void RationalApprox::last_inv_approx_coeff(void)
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside RationalApprox::last_inv_approx_coef ..."<<endl;
 #endif
 int i;
 REAL min, max, epsilon;
 REAL *minmax;
 minmax = new REAL [2];

 // normalized coefficients
 *this=*last_inv_approx_norm_coeff;

 // findminmax(min, max);
 findminmax_con_openacc(minmax);
 min = minmax[0];
 max = minmax[1];

 min*=0.95;
 max*=1.05;

 epsilon=min/max;

 if(epsilon<min_epsilon)
   {
   ofstream err_file;
   err_file.open(QUOTEME(ERROR_FILE), ios::app);   
   err_file << "WARNING: in first_last_approx_coef epsilon="<<epsilon<<" < min_epsilon="<<min_epsilon;
   err_file << " at update_iteration="<<update_iteration <<endl;
   err_file.close();
   } 

  
 // rescale coeff.
 double exponent = -((double) no_flavours)/4./((double)no_ps);
 epsilon=pow(max, exponent);
 RA_a0*=epsilon;
 for(i=0; i<approx_order; i++)
    {
    RA_a[i]*=(max*epsilon);
    RA_b[i]*=max;
    }

 #ifdef DEBUG_MODE
 cout << "\tterminated RationalApprox::last_inv_approx_coeff, epsilon="<< min/max<< " min_epsilon="<<min_epsilon<<endl;
 #endif
 }


                                     // EXPLICIT CALCULATIONS ON FERMIONS

void first_inv_approx_calc(REAL res)
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside first_inv_approx_calc ..."<<endl;
 #endif

 RationalApprox approx;
 approx.first_inv_approx_coeff(); 

 #ifdef USE_GPU
 int order;
 double constant, numerators[max_approx_order];

 cu_multips_shifted_invert (res, approx);

 get_order(order, approx);
 get_const(constant, approx);
 get_numerators(numerators, approx);

 cuda_sol_sum(order, constant, numerators);

 smartunpack_multifermion(fermion_chi, chi_packed);
 #endif

 #ifndef USE_GPU
 int iter, pseudofermion;
 long int i;
 Vec3 vr_1;

 multips_shifted_invert (fermion_shiftmulti, fermion_phi, res, approx);

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
 #endif

 #ifdef DEBUG_MODE
 cout << "\tterminated first_inv_approx_calc"<<endl;
 #endif
 }

/*

void last_inv_approx_calc(REAL res)
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside  last_inv_approx_calc  .."<<endl;
 #endif

 RationalApprox approx;
 approx.last_inv_approx_coeff(); 

 #ifdef USE_GPU
 int order;
 double constant, numerators[max_approx_order];

 cu_multips_shifted_invert (res, approx);

 get_order(order, approx);
 get_const(constant, approx);
 get_numerators(numerators, approx);

 cuda_sol_sum(order, constant, numerators);

 smartunpack_multifermion(fermion_phi, chi_packed);
 #endif

 #ifndef USE_GPU
 int iter, pseudofermion;
 long int i;
 Vec3 vr_1;

 multips_shifted_invert (fermion_shiftmulti, fermion_chi, res, approx);

 for(pseudofermion=0; pseudofermion<no_ps; pseudofermion++)
    {
    for(i=0; i<sizeh; i++)
       {
       vr_1=(approx.RA_a0)*(fermion_chi->fermion[pseudofermion][i]);
       for(iter=0; iter<(approx.approx_order); iter++)
          {
          vr_1+=(approx.RA_a[iter])*(fermion_shiftmulti->fermion[pseudofermion][iter][i]);
          }
       fermion_phi->fermion[pseudofermion][i]=vr_1;
       }
    }
 #endif

 #ifdef DEBUG_MODE
 cout << "\tterminated root_1_12_calc"<<endl;
 #endif
 }

*/
