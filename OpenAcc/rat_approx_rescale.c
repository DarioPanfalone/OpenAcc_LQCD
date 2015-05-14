void rescale_rational_approximation(COM_RationalApprox * mother,COM_RationalApprox * daughter, double *minmax, double power)
{

  daughter[0].COM_approx_order = mother[0].COM_approx_order; 
  daughter[0].COM_min_epsilon = mother[0].COM_min_epsilon;
  daughter[0].COM_RA_a0 = mother[0].COM_RA_a0;
  for(int ii=0; ii<daughter[0].COM_approx_order; ii++)
    {
      daughter[0].COM_RA_a[ii] = mother[0].COM_RA_a[ii] ;
      daughter[0].COM_RA_b[ii] = mother[0].COM_RA_b[ii] ;
    }
  double min =  minmax[0];
  double max =  minmax[1];
  min*=0.95;
  max*=1.05;
  double epsilon=min/max;
  if(epsilon<daughter[0].COM_min_epsilon)
    {
      printf("WARNING: in first_inv_approx_coef epsilon= %lg  while min_epsilon= %lg \n", epsilon, daughter[0].COM_min_epsilon);
    }
  double exponent = power;
  epsilon=pow(max, exponent);  
  daughter[0].COM_RA_a0*=epsilon;
  for(int ii=0; ii<daughter[0].COM_approx_order; ii++)
    {
      daughter[0].COM_RA_a[ii] = daughter[0].COM_RA_a[ii] * max * epsilon ;
      daughter[0].COM_RA_b[ii] = daughter[0].COM_RA_b[ii] * max ;
    }

}
