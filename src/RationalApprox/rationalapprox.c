#ifndef RATIONAL_APPROX_C_
#define RATIONAL_APPROX_C_

#include "./rationalapprox.h"


char* rational_approx_filename(int approx_order, int exponent_num, int exponent_den, double lambda_min){

  char Cn1[3];
  char Cy1[3];
  char Cz1[3];
  char mloglmin[5];

  sprintf(Cn1, "%d", approx_order);  // order
  sprintf(Cy1, "%d", exponent_num); // num
  sprintf(Cz1, "%d", exponent_den); // den
  sprintf(mloglmin, "%1.1f", -log(lambda_min)/log(10.0)); // lambda min

  char * nomefile = (char*)malloc(50*sizeof(char));
  strcpy(nomefile,"approx_");
  strcat(nomefile,Cy1);
  strcat(nomefile,"_over_");
  strcat(nomefile,Cz1);
  strcat(nomefile,"_order_");
  strcat(nomefile,Cn1);
  strcat(nomefile,"_mloglm_");
  strcat(nomefile,mloglmin);
  strcat(nomefile,".REMEZ");

  return nomefile;

}


void rationalapprox_read(RationalApprox* rational_approx)
{

    // CALCULATION OF COEFFICIENTS FOR FIRST_INV_APPROX_NORM_COEFF

    char * nomefile = rational_approx_filename(rational_approx->approx_order,rational_approx->exponent_num,rational_approx->exponent_den,rational_approx->lambda_min);

    rationalapprox_read_custom_nomefile(rational_approx,nomefile);

}


void rationalapprox_read_custom_nomefile(RationalApprox* rational_approx, char* nomefile)
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
//      printf("RA_a[%d] = %18.16e, RA_b[%d] = %18.16e\n", i, rational_approx->RA_a[i], i, rational_approx->RA_b[i]);
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

void rescale_rational_approximation(RationalApprox *in, RationalApprox *out, double *minmax){

   double power = (double) in->exponent_num/in->exponent_den;

   

   out->exponent_num        = in->exponent_num       ;
   out->exponent_den        = in->exponent_den       ;
   out->approx_order        = in->approx_order       ;
   out->gmp_remez_precision = in->gmp_remez_precision;              
   out->error               = in->error              ;
   
   // HERE THE ASSUMPTION IS THAT in->lambda_max  = 1
   
   double min =  minmax[0];
   double max =  minmax[1];
   min*=0.95;
   max*=1.05;
   double epsilon=pow(max, power);  
   out->RA_a0               = in->RA_a0       *     epsilon ;
   for(int order = 0; order < in->approx_order; order ++){
   out->RA_a[order] = in->RA_a[order]*max * epsilon;
   out->RA_b[order] = in->RA_b[order]*max ;
   }

   out->lambda_min  = in->lambda_min * max ;
   out->lambda_max  = max ;
   //pray
   if(out->lambda_min > minmax[0]){
       printf("Warning: mother rational approx does not cover the range!\n");
       printf("out->lambda_min: %e , minmax[0]: %e\n", out->lambda_min, minmax[0] );
   }

}


void renormalize_rational_approximation(RationalApprox *in, RationalApprox *out){
   double power = (double) in->exponent_num/in->exponent_den;

   

   out->exponent_num        = in->exponent_num       ;
   out->exponent_den        = in->exponent_den       ;
   out->approx_order        = in->approx_order       ;
   out->gmp_remez_precision = in->gmp_remez_precision;              
   out->error               = in->error              ;
   
   // HERE THE ASSUMPTION IS THAT in->lambda_max  = 1
   
   double rescale_ratio =  1/in->lambda_max;
   double epsilon=pow(rescale_ratio, power);  
   out->RA_a0               = in->RA_a0       *     epsilon ;
   for(int order = 0; order < in->approx_order; order ++){
   out->RA_a[order] = in->RA_a[order]*rescale_ratio * epsilon;
   out->RA_b[order] = in->RA_b[order]*rescale_ratio ;
   }

   out->lambda_min  = in->lambda_min * rescale_ratio ;
   out->lambda_max  = 1 ;

}




#endif


