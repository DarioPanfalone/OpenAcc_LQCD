/***************************************************
 *  * The format for the rational approximation files *
 *  (file starts)
 *  Approximation to f(x) = (x)^(%d/%i)
 *  Order: %i
 *  Lambda min: %lf
 *  Lambda max: %lf
 *  GMP Remez Precision: %i
 *  error: %lf
 *  RA_a0 = %lf
 *  RA_a[%i] = %lf, RA_b[%i] = %lf\n#(Order times!)   
 *  (file ends)
 *  **********************************************/


#ifndef RATIONAL_APPROX_C_
#define RATIONAL_APPROX_C_

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define MAX_APPROX_ORDER 18


typedef struct RationalApprox_t{
    int exponent_num;
    int exponent_den;
    int approx_order;
    double lambda_min;
    double lambda_max; // usually equal to one
    int gmp_remez_precision;
    double error;
    double RA_a0;
    double RA_a[MAX_APPROX_ORDER];
    double RA_b[MAX_APPROX_ORDER];
} RationalApprox;


char* rational_approx_filename(int approx_order, int exponent_num, int exponent_den){

  char Cn1[3];
  char Cy1[3];
  char Cz1[3];
  sprintf(Cn1, "%d", approx_order);  // order
  sprintf(Cy1, "%d", exponent_num); // num
  sprintf(Cz1, "%d", exponent_den); // den

  char * nomefile = (char*)malloc(50*sizeof(char));
  strcpy(nomefile,"approx_");
  strcat(nomefile,Cy1);
  strcat(nomefile,"_over_");
  strcat(nomefile,Cz1);
  strcat(nomefile,"_order_");
  strcat(nomefile,Cn1);
  strcat(nomefile,".REMEZ");

  return nomefile;

}


void rationalapprox_read(const char* nomefile, RationalApprox* rational_approx )
{

    // CALCULATION OF COEFFICIENTS FOR FIRST_INV_APPROX_NORM_COEFF

    FILE *input = fopen(nomefile, "rt");
    printf("%s\n", nomefile );
    if (input == NULL) {
        printf("Could not open file %s \n",nomefile );
        exit(-1);
    }

    // The partial fraction expansion takes the form 
    // r(x) = norm + sum_{k=1}^{n} res[k] / (x + pole[k])
    fscanf(input,"\nApproximation to f(x) = (x)^(%i/%i)\n", &(rational_approx->exponent_num), &(rational_approx->exponent_den));
    fscanf(input,"Order: %i\n", &(rational_approx->approx_order));
    fscanf(input,"Lambda Min: %lf\n", &(rational_approx->lambda_min));
    fscanf(input,"Lambda Max: %lf\n", &(rational_approx->lambda_max));
    fscanf(input,"GMP Remez Precision: %i\n", &(rational_approx->gmp_remez_precision));
    fscanf(input,"Error: %lf\n", &(rational_approx->error));
    fscanf(input, "RA_a0 = %lf\n",&(rational_approx->RA_a0));
    for(int i = 0; i < rational_approx->approx_order; i++)
    {
        fscanf(input, "RA_a[%i] = %lf, RA_b[%i] = %lf\n", &i, &(rational_approx->RA_a[i]), &i, &(rational_approx->RA_b[i]));
    }
    fclose(input);

    printf("RA_a0 = %18.16e\n", rational_approx->RA_a0);
    for(int i = 0; i < rational_approx->approx_order; i++) 
    {
        printf("RA_a[%d] = %18.16e, RA_b[%d] = %18.16e\n", i, rational_approx->RA_a[i], i, rational_approx->RA_b[i]);
    } 
}

void rationalapprox_save(const char* nomefile, RationalApprox* rational_approx){

    FILE *output = fopen(nomefile, "w");
    printf("%s\n", nomefile );
    if (output == NULL) {
        printf("Could not open file %s \n",nomefile );
        exit(-1);
    }

    // The partial fraction expansion takes the form 
    // r(x) = norm + sum_{k=1}^{n} res[k] / (x + pole[k])
    fprintf(output,"\nApproximation to f(x) = (x)^(%i/%i)\n", rational_approx->exponent_num, rational_approx->exponent_den);
    fprintf(output,"Order: %i\n", rational_approx->approx_order);
    fprintf(output,"Lambda Min: %18.16e\n", rational_approx->lambda_min);
    fprintf(output,"Lambda Max: %18.16e\n", rational_approx->lambda_max);
    fprintf(output,"GMP Remez Precision: %i\n", rational_approx->gmp_remez_precision);
    fprintf(output,"Error: %18.16e\n", rational_approx->error);
    fprintf(output, "RA_a0 = %18.16e\n",rational_approx->RA_a0);
    for(int i = 0; i < rational_approx->approx_order; i++)
    {
        fprintf(output,"RA_a[%d] = %18.16e, RA_b[%d] = %18.16e\n", i, rational_approx->RA_a[i], i, rational_approx->RA_b[i]);
    }
    fclose(output);

}


#endif


