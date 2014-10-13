
extern "C" void cuda_init0(void)
  {
  #ifdef DEBUG_MODE
  printf("DEBUG: inside cuda_init0 ...\n");
  #endif

  size_t gauge_field_size_f = 3*sizeof(float4)*no_links;    //first two lines only of each SU(3) matrix

  // allocate & initialize gauge
  // 2 since 1double~2float

  cudaSafe(AT,cudaMalloc((void**)&gauge_field_device, 2*gauge_field_size_f), "cudaMalloc");
  cudaSafe(AT,cudaMemcpy(gauge_field_device, gauge_field_packed, 2*gauge_field_size_f, cudaMemcpyHostToDevice), 
                 "cudaMemcpy");


  // allocate & initialize device_table
  cudaSafe(AT,cudaMalloc((void**)&device_table, sizeof(int)*size*8), "cudaMalloc");
  cudaSafe(AT,cudaMemcpy(device_table, shift_table, sizeof(int)*size*8, cudaMemcpyHostToDevice), "cudaMemcpy"); 

  // allocate & initialize device_phases
  cudaSafe(AT,cudaMalloc((void**)&device_phases, sizeof(int)*size*4), "cudaMalloc");
  cudaSafe(AT,cudaMemcpy(device_phases, eta, sizeof(int)*size*4, cudaMemcpyHostToDevice), "cudaMemcpy"); 


  // initialize constants
  float mass_l=(float) mass;
  cudaSafe(AT,cudaMemcpyToSymbol(mass_dev, &mass_l, sizeof(float), 0, cudaMemcpyHostToDevice), "cudaMemcpyToSymbol");
  cudaSafe(AT,cudaMemcpyToSymbol(mass_d_dev, &mass, sizeof(double), 0, cudaMemcpyHostToDevice), "cudaMemcpyToSymbol");
  int size_l=(int) size;
  cudaSafe(AT,cudaMemcpyToSymbol(size_dev, &size_l, sizeof(int), 0, cudaMemcpyHostToDevice), "cudaMemcpyToSymbol");
  size_l=(int) sizeh;
  cudaSafe(AT,cudaMemcpyToSymbol(size_dev_h, &size_l, sizeof(int), 0, cudaMemcpyHostToDevice), "cudaMemcpyToSymbol");


  #ifdef DEBUG_MODE
  printf("\tterminated cuda_init0\n");
  #endif
  }


extern "C" void cuda_init1(void)
  {
  #ifdef DEBUG_MODE
  printf("DEBUG: inside cuda_init1 ...\n");
  #endif

  size_t vector_size_f   = sizeof(float2)*3*size;           // 2(complex)*3(su3_vector)

  // allocate & initialize mf_device
  // again 2 since 1double~2float
  cudaSafe(AT,cudaMalloc((void**)&mf_device, 2*vector_size_f), "cudaMalloc"); 
  cudaSafe(AT,cudaMemset(mf_device, 0, 2*vector_size_f), "cudaMemset");  // initialize even and odd to 0
     
     // 1st float
     cudaSafe(AT,cudaMemcpy(mf_device , chi_packed , 
                                    size*sizeof(float), cudaMemcpyHostToDevice), "cudaMemcpy");
     cudaSafe(AT,cudaMemcpy(mf_device +   size , chi_packed +   size , 
                                    size*sizeof(float), cudaMemcpyHostToDevice), "cudaMemcpy");
     cudaSafe(AT,cudaMemcpy(mf_device + 2*size , chi_packed + 2*size , 
                                    size*sizeof(float), cudaMemcpyHostToDevice), "cudaMemcpy");

     // 2nd float
     cudaSafe(AT,cudaMemcpy(mf_device +          3*size, chi_packed +          3*size,  size*sizeof(float), cudaMemcpyHostToDevice), "cudaMemcpy");
     cudaSafe(AT,cudaMemcpy(mf_device +   size + 3*size, chi_packed +   size + 3*size,  size*sizeof(float), cudaMemcpyHostToDevice), "cudaMemcpy");
     cudaSafe(AT,cudaMemcpy(mf_device + 2*size + 3*size, chi_packed + 2*size + 3*size,  size*sizeof(float), cudaMemcpyHostToDevice), "cudaMemcpy");

     

  // allocate & initialize to zero smf_device (even & odd)
  // again 2 since 1double~2float
  cudaSafe(AT,cudaMalloc((void**)&smf_device, 2*vector_size_f), "cudaMalloc"); 
  cudaSafe(AT,cudaMemset(smf_device, 0, 2*vector_size_f), "cudaMemset"); 

  #ifdef DEBUG_MODE
  printf("\tterminated cuda_init1\n");
  #endif
  }


extern "C" void cuda_end(void)
  {
  #ifdef DEBUG_MODE
  printf("DEBUG: inside cuda_end ...\n");
  #endif

  cudaSafe(AT,cudaMemcpy(gauge_field_packed, gauge_field_device, 2*3*no_links*sizeof(float4), cudaMemcpyDeviceToHost), 
                                                                                                            "cudaMemcpy");

  cudaSafe(AT,cudaFree(gauge_field_device), "cudaFree");
  cudaSafe(AT,cudaFree(device_table), "cudaFree");
  cudaSafe(AT,cudaFree(device_phases), "cudaFree");

  cudaSafe(AT,cudaFree(mf_device), "cudaFree");
  cudaSafe(AT,cudaFree(smf_device), "cudaFree");

  #ifdef DEBUG_MODE
  printf("\tterminated cuda_end\n");
  #endif
  }


extern "C" void cuda_get_conf(void)
 {
 #ifdef DEBUG_MODE
 printf("DEBUG: inside cuda_get_conf ...\n");
 #endif

 cudaSafe(AT,cudaMemcpy(gauge_field_packed, gauge_field_device, 2*12*no_links*sizeof(float), cudaMemcpyDeviceToHost), 
       "cudaMemcpy");

 #ifdef DEBUG_MODE
 printf("\tterminated cuda_get_conf ...\n");
 #endif
 }


