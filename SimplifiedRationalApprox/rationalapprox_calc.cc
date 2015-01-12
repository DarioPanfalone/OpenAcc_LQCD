#include <stdlib.h>
#include <stdio.h>
#include <math.h>


char* concat(char *s1, char *s2)
{
  size_t len1 = strlen(s1);
  size_t len2 = strlen(s2);
  char *result = (char *) malloc(len1+len2+1);//+1 for the zero-terminator
  memcpy(result, s1, len1);
  memcpy(result+len1, s2, len2+1);//+1 to copy the null-terminator
  return result;
}


void rationalapprox_calc(void)
  {
    char * primop  = "./REMEZ_SAVE/approx_plus_";
    char * primom  = "./REMEZ_SAVE/approx_minus_";
    char * secondo = "_over_";
    char * terzo   = "_order_";
    char * quarto  = ".REMEZ";

  // CALCULATION OF COEFFICIENTS FOR FIRST_INV_APPROX_NORM_COEFF
  int n=approx_metro; // The degree of the numerator polynomial
  int d=approx_metro; // The degree of the denominator polynomial
  int y1=no_flavours; // The numerator of the exponent
  int z1=8*no_ps;     // The denominator of the exponent

  int precision=gmp_remez_precision; // The precision that gmp uses
  double lambda_low=lambda_min_metro;
  double  lambda_high=1.0;              // The bounds of the approximation

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

  FILE *input1 = fopen(name1, "rt");
  printf("%s\n", name1);
  if (input1 == NULL) {
    printf("Could not open file %s \n", name1);
    exit(-1);
  }

  // The partial fraction expansion takes the form 
  // r(x) = norm + sum_{k=1}^{n} res[k] / (x + pole[k])
  double norm;
  double *res1 = new double[n];
  double *pole1 = new double[d];
  y1=0;
  z1=0;
  fscanf(input1,"\nApproximation to f(x) = (x)^(%i/%i)\n\n", &y1, &z1);
  fscanf(input1, "RA_a0 = %lf\n", &norm);
  for(int i = 0; i < n; i++)
    {
      fscanf(input1, "RA_a[%i] = %lf, RA_b[%i] = %lf\n", &i, &res1[i], &i, &pole1[i]);
    }
  fclose(input1);

  printf("RA_a0 = %18.16e\n", norm);
  for(int i = 0; i < n; i++) 
     {
     printf("RA_a[%d] = %18.16e, RA_b[%d] = %18.16e\n", i, res1[i], i, pole1[i]);
     } 

  first_inv_approx_norm_coeff = new RationalApprox(lambda_low, n, norm, res1, pole1);  

  delete res1;
  delete pole1;



  // CALCULATION OF COEFFICIENTS FOR MD_INV_APPROX_NORM_COEFF
  n=approx_md; // The degree of the numerator polynomial
  d=approx_md; // The degree of the denominator polynomial
  y1=no_flavours; // The numerator of the exponent
  z1=4*no_ps;     // The denominator of the exponent
  precision=gmp_remez_precision; // The precision that gmp uses
  lambda_low=lambda_min_md;      // The lower bounds of the approximation

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

  printf("\nApproximation to f(x) = (x)^(-%d/%d) on [%e, %e]\n\n", y1, z1, lambda_low, lambda_high);
  fflush(stdout);

  double *res2 = new double[n];
  double *pole2 = new double[d];
  FILE *input2 = fopen(name2, "rt");
  printf("%s\n", name2);
  if (input2 == NULL) {
    printf("Could not open file %s \n", name2);
    exit(-1);
  }
  fscanf(input2,"\nApproximation to f(x) = (x)^(%i/%i)\n\n", &y1, &z1);
  fscanf(input2, "RA_a0 = %lf\n", &norm);
  for(int i = 0; i < n; i++)
    {
      fscanf(input2, "RA_a[%i] = %lf, RA_b[%i] = %lf\n", &i, &res2[i], &i, &pole2[i]);
    }
  fclose(input2);

  printf("RA_a0 = %18.16e\n", norm);
  for(int i = 0; i < n; i++) 
     {
     printf("RA_a[%d] = %18.16e, RA_b[%d] = %18.16e\n", i, res2[i], i, pole2[i]);
     } 

  md_inv_approx_norm_coeff=new RationalApprox(lambda_low, n, norm, res2, pole2);  

  delete res2;
  delete pole2;
  

  // CALCULATION OF COEFFICIENTS FOR LAST_INV_APPROX_NORM_COEFF

  n=approx_metro; // The degree of the numerator polynomial
  d=approx_metro; // The degree of the denominator polynomial
  y1=no_flavours; // The numerator of the exponent
  z1=4*no_ps;     // The denominator of the exponent
  precision=gmp_remez_precision; // The precision that gmp uses
  lambda_low=lambda_min_metro;      // The lower bounds of the approximation


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


  double *res3 = new double[n];
  double *pole3 = new double[d];

  printf("\nApproximation to f(x) = (x)^(-%d/%d) on [%e, %e]\n\n", y1, z1, lambda_low, lambda_high);
  fflush(stdout);

  FILE *input3 = fopen(name3, "rt");
  printf("%s\n", name3);
  if (input3 == NULL) {
    printf("Could not open file %s \n", name3);
    exit(-1);
  }

  fscanf(input3,"\nApproximation to f(x) = (x)^(%i/%i)\n\n", &y1, &z1);
  fscanf(input3, "RA_a0 = %lf\n", &norm);
  for(int i = 0; i < n; i++)
    {
      fscanf(input3, "RA_a[%i] = %lf, RA_b[%i] = %lf\n", &i, &res3[i], &i, &pole3[i]);
    }
  fclose(input3);

  printf("RA_a0 = %18.16e\n", norm);
  for(int i = 0; i < n; i++) 
     {
     printf("RA_a[%d] = %18.16e, RA_b[%d] = %18.16e\n", i, res3[i], i, pole3[i]);
     } 

  last_inv_approx_norm_coeff= new RationalApprox(lambda_low, n, norm, res3, pole3);  



  delete res3;
  delete pole3;
  }
