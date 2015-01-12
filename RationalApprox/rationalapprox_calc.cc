#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include"../Remez/bigfloat.h"
#include"../Remez/alg_remez.h"
#include"../Remez/alg_remez.cc"


//using namespace std;
char* concat(char *s1, char *s2)
{
  size_t len1 = strlen(s1);
  size_t len2 = strlen(s2);
  char *result = (char *) malloc(len1+len2+1);//+1 for the zero-terminator
  //in real code you would check for errors in malloc here
  memcpy(result, s1, len1);
  memcpy(result+len1, s2, len2+1);//+1 to copy the null-terminator
  return result;
}

const int gmp_remez_precision=100; // The precision that gmp uses

int no_flavours;
int no_ps;
//Rational approximations for Metropolis
int approx_metro;//=19;
const double lambda_min_metro=4.0e-7;  // rational approx valid on [lambda_min_metro, 1.0]
const double residue_metro=1.0e-8;    // stopping residual for CG

//Rational approximations for MD
int approx_md;//=9;
const double lambda_min_md=4.0e-7;  // rational approx valid on [lambda_min_metro, 1.0]
const double residue_md=1.0e-3;    // stopping residual for CG


int main(int argc, char **argv){
  if(argc!=5){
    printf("Use as arguments:    approx_order_x_metro     approx_order_x_md      no_flavours      no_ps \n");
    return 0;
  }
  char * primop  = "REMEZ_SAVE/approx_plus_";
  char * primom  = "REMEZ_SAVE/approx_minus_";
  char * secondo = "_over_";
  char * terzo   = "_order_";
  char * quarto  = ".REMEZ";

  approx_metro = atoi(argv[1]);
  approx_md    = atoi(argv[2]);
  no_flavours  = atoi(argv[3]);
  no_ps        = atoi(argv[4]);

  // CALCULATION OF COEFFICIENTS FOR FIRST_INV_APPROX_NORM_COEFF

  int n=approx_metro; // The degree of the numerator polynomial
  int d=approx_metro; // The degree of the denominator polynomial
  int y1=no_flavours; // The numerator of the exponent
  int z1=8*no_ps;     // The denominator of the exponent
  int precision=gmp_remez_precision; // The precision that gmp uses
  double lambda_low=lambda_min_metro;
  double  lambda_high=1.0;              // The bounds of the approximation

  // The error from the approximation (the relative error is minimised
  // - if another error minimisation is requried, then line 398 in
  // alg_remez.C is where to change it)
  double error;

  // The partial fraction expansion takes the form 
  // r(x) = norm + sum_{k=1}^{n} res[k] / (x + pole[k])
  double norm;
  double *res1 = new double[n];
  double *pole1 = new double[d];
 
  printf("\nApproximation to f(x) = (x)^(%d/%d) on [%e, %e]\n\n", y1, z1, lambda_low, lambda_high);
  fflush(stdout);
  // Instantiate the Remez class
  AlgRemez remez1(lambda_low,lambda_high,precision);
  // Generate the required approximation
  error = remez1.generateApprox(n,d,y1,z1);
  printf("approximation error = %e\n\n", error);

  // Find the partial fraction expansion of the approximation 
  // to the function x^{y/z} (this only works currently for 
  // the special case that n = d)
  remez1.getPFE(res1,pole1,&norm);

  char Cn1[3];
  char Cy1[3];
  char Cz1[3];
  sprintf(Cn1, "%d", n);  // order
  sprintf(Cy1, "%d", y1); // num
  sprintf(Cz1, "%d", z1); // den
  char * name1 = concat(primop,Cy1);
  name1 = concat(name1,secondo);
  name1 = concat(name1,Cz1);
  name1 = concat(name1,terzo);
  name1 = concat(name1,Cn1);
  name1 = concat(name1,quarto);

  FILE *output = fopen(name1, "w");
  fprintf(output, "\nApproximation to f(x) = (x)^(%d/%d)\n\n", y1, z1);
  fprintf(output, "RA_a0 = %18.16e\n", norm);
  for(int i = 0; i < n; i++) 
     {
     fprintf(output, "RA_a[%d] = %18.16e, RA_b[%d] = %18.16e\n", i, res1[i], i, pole1[i]);
     } 
  fclose(output);



  printf("RA_a0 = %18.16e\n", norm);
  for(int i = 0; i < n; i++) 
     {
     printf("RA_a[%d] = %18.16e, RA_b[%d] = %18.16e\n", i, res1[i], i, pole1[i]);
     } 

  delete [] res1;
  delete [] pole1;



  // CALCULATION OF COEFFICIENTS FOR MD_INV_APPROX_NORM_COEFF

  n=approx_md; // The degree of the numerator polynomial
  d=approx_md; // The degree of the denominator polynomial
  y1=no_flavours; // The numerator of the exponent
  z1=4*no_ps;     // The denominator of the exponent
  precision=gmp_remez_precision; // The precision that gmp uses
  lambda_low=lambda_min_md;      // The lower bounds of the approximation

  double *res2 = new double[n];
  double *pole2 = new double[d];

  printf("\nApproximation to f(x) = (x)^(-%d/%d) on [%e, %e]\n\n", y1, z1, lambda_low, lambda_high);
  fflush(stdout);
  // Instantiate the Remez class
  AlgRemez remez2(lambda_low,lambda_high,precision);

  // Generate the required approximation
  error = remez2.generateApprox(n,d,y1,z1);

  printf("approximation error = %e\n\n", error);

  // Find pfe of inverse function
  remez2.getIPFE(res2,pole2,&norm);

  char Cn2[3];
  char Cy2[3];
  char Cz2[3];
  sprintf(Cn2, "%d", n);  // order
  sprintf(Cy2, "%d", y1); // num
  sprintf(Cz2, "%d", z1); // den
  char * name2 = concat(primom,Cy2);
  name2 = concat(name2,secondo);
  name2 = concat(name2,Cz2);
  name2 = concat(name2,terzo);
  name2 = concat(name2,Cn2);
  name2 = concat(name2,quarto);

  FILE *output2 = fopen(name2, "w");
  fprintf(output2, "\nApproximation to f(x) = (x)^(%d/%d)\n\n", y1, z1);
  fprintf(output2, "RA_a0 = %18.16e\n", norm);
  for(int i = 0; i < n; i++) 
     {
     fprintf(output2, "RA_a[%d] = %18.16e, RA_b[%d] = %18.16e\n", i, res2[i], i, pole2[i]);
     } 
  fclose(output2);

  printf("RA_a0 = %18.16e\n", norm);
  for(int i = 0; i < n; i++) 
     {
     printf("RA_a[%d] = %18.16e, RA_b[%d] = %18.16e\n", i, res2[i], i, pole2[i]);
     } 


  delete [] res2;
  delete [] pole2;


  // CALCULATION OF COEFFICIENTS FOR LAST_INV_APPROX_NORM_COEFF

  n=approx_metro; // The degree of the numerator polynomial
  d=approx_metro; // The degree of the denominator polynomial
  y1=no_flavours; // The numerator of the exponent
  z1=4*no_ps;     // The denominator of the exponent
  precision=gmp_remez_precision; // The precision that gmp uses
  lambda_low=lambda_min_metro;      // The lower bounds of the approximation

  double *res3 = new double[n];
  double *pole3 = new double[d];

  printf("\nApproximation to f(x) = (x)^(-%d/%d) on [%e, %e]\n\n", y1, z1, lambda_low, lambda_high);
  fflush(stdout);

  // Instantiate the Remez class
  AlgRemez remez3(lambda_low,lambda_high,precision);

  // Generate the required approximation
  error = remez3.generateApprox(n,d,y1,z1);

  printf("approximation error = %e\n\n", error);

  // Find pfe of inverse function
  remez3.getIPFE(res3,pole3,&norm);

  char Cn3[3];
  char Cy3[3];
  char Cz3[3];
  sprintf(Cn3, "%d", n);  // order
  sprintf(Cy3, "%d", y1); // num
  sprintf(Cz3, "%d", z1); // den
  char * name3 = concat(primom,Cy3);
  name3 = concat(name3,secondo);
  name3 = concat(name3,Cz3);
  name3 = concat(name3,terzo);
  name3 = concat(name3,Cn3);
  name3 = concat(name3,quarto);

  FILE *output3 = fopen(name3, "w");

  fprintf(output3, "\nApproximation to f(x) = (x)^(%d/%d)\n\n", y1, z1);
  fprintf(output3, "RA_a0 = %18.16e\n", norm);
  for(int i = 0; i < n; i++) 
     {
     fprintf(output3, "RA_a[%d] = %18.16e, RA_b[%d] = %18.16e\n", i, res3[i], i, pole3[i]);
     } 
  fclose(output3);

  printf("RA_a0 = %18.16e\n", norm);
  for(int i = 0; i < n; i++) 
     {
     printf("RA_a[%d] = %18.16e, RA_b[%d] = %18.16e\n", i, res3[i], i, pole3[i]);
     } 

  delete [] res3;
  delete [] pole3;


}
