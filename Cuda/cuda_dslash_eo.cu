// SINGLE PRECISION KERNEL for even/odd fermions

__global__ void switchParityKernel(float2 *out,
                                     float2 *in)
{

    int idx = blockIdx.x * blockDim.x + threadIdx.x ; 

    float2 temp;
    //1st float
    for(int color = 0 ; color < 3 ; color++){
        temp = in[idx+(color*2+1)*size_dev_h];
    	out[idx+(color*2+1)*size_dev_h] = in[idx+color*2*size_dev_h];
        out[idx+color*2*size_dev_h] = temp ;
    }
    //2nd float
    for(int color = 3 ; color < 6 ; color++){
        temp = in[idx+(color*2+1)*size_dev_h];
    	out[idx+(color*2+1)*size_dev_h] = in[idx+color*2*size_dev_h];
        out[idx+color*2*size_dev_h] = temp ;
    }



}

__global__ void copyKernel(float2 *out,
                                     float2 *in)
{

    int idx = blockIdx.x * blockDim.x + threadIdx.x ; 

    
    //1st float
    for(int color = 0 ; color < 3 ; color++){
        out[idx+color*2*size_dev_h] = in[idx+color*2*size_dev_h];
        out[idx+(color*2+1)*size_dev_h] = in[idx+(color*2+1)*size_dev_h];
    }
    //2nd float
    for(int color = 3 ; color < 6 ; color++){
        out[idx+color*2*size_dev_h] = in[idx+color*2*size_dev_h];
        out[idx+(color*2+1)*size_dev_h] = in[idx+(color*2+1)*size_dev_h];
    }

}



__global__ void DslashKernelEO(float2  *out,
                                 float2  *in,
                                 int *tables, 
                                 int *phases, 
                                 size_t gauge_offset)
{
  int idx = blockIdx.x * blockDim.x + threadIdx.x + size_dev_h;  // idx>sizeh, ODD

  float stag_phase = 1.0;
  //if(idx < 192) return;

  //Store result in sharedMem
  __shared__ float ferm_out[3][2][NUM_THREADS];

  //New tables indexing (index fastest)
  __shared__ int site_table[NUM_THREADS];

  //Load link matrix U_mu(ix) in registers
  float link0x, link0y, link0z, link0w, 
         link1x, link1y, link1z, link1w, 
         link2x, link2y, link2z, link2w;   
  float4 auxlink;

  float2 ferm_in_0, ferm_in_1, ferm_in_2;


  // DIRECTION 0
  site_table[threadIdx.x]  = tables[idx+4*size_dev];

  ferm_in_0 = in[              site_table[threadIdx.x]];
  ferm_in_1 = in[   size_dev + site_table[threadIdx.x]];
  ferm_in_2 = in[ 2*size_dev + site_table[threadIdx.x]];

  // 1st float 
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(0+3*0));
  link0x=(float) auxlink.x;
  link0y=(float) auxlink.y;
  link0z=(float) auxlink.z;
  link0w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(1+3*0));
  link1x=(float) auxlink.x;
  link1y=(float) auxlink.y;
  link1z=(float) auxlink.z;
  link1w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(2+3*0));
  link2x=(float) auxlink.x;
  link2y=(float) auxlink.y;
  link2z=(float) auxlink.z;
  link2w=(float) auxlink.w;
  

  ferm_out[0][0][threadIdx.x] = link0x*ferm_in_0.x-link0y*ferm_in_0.y+  
                                link0z*ferm_in_1.x-link0w*ferm_in_1.y+ 
                                link1x*ferm_in_2.x-link1y*ferm_in_2.y; 
  ferm_out[0][1][threadIdx.x] = link0x*ferm_in_0.y+link0y*ferm_in_0.x+ 
                                link0z*ferm_in_1.y+link0w*ferm_in_1.x+ 
                                link1x*ferm_in_2.y+link1y*ferm_in_2.x; 

  ferm_out[1][0][threadIdx.x] = link1z*ferm_in_0.x-link1w*ferm_in_0.y+  
                                link2x*ferm_in_1.x-link2y*ferm_in_1.y+ 
                                link2z*ferm_in_2.x-link2w*ferm_in_2.y; 
  ferm_out[1][1][threadIdx.x] = link1z*ferm_in_0.y+link1w*ferm_in_0.x+ 
                                link2x*ferm_in_1.y+link2y*ferm_in_1.x+ 
                                link2z*ferm_in_2.y+link2w*ferm_in_2.x; 

  ferm_out[2][0][threadIdx.x] = C1RED*ferm_in_0.x-C1IMD*ferm_in_0.y+  
                                C2RED*ferm_in_1.x-C2IMD*ferm_in_1.y+ 
                                C3RED*ferm_in_2.x-C3IMD*ferm_in_2.y; 
  ferm_out[2][1][threadIdx.x] = C1RED*ferm_in_0.y+C1IMD*ferm_in_0.x+ 
                                C2RED*ferm_in_1.y+C2IMD*ferm_in_1.x+ 
                                C3RED*ferm_in_2.y+C3IMD*ferm_in_2.x; 

  //DIRECTION 1
  site_table[threadIdx.x] = tables[idx+5*size_dev];
  stag_phase              = (float) phases[idx+size_dev];

  ferm_in_0 = in[              site_table[threadIdx.x]];
  ferm_in_1 = in[   size_dev + site_table[threadIdx.x]];
  ferm_in_2 = in[ 2*size_dev + site_table[threadIdx.x]];


  // 1st float 
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(0+3*1));
  link0x=(float) auxlink.x;
  link0y=(float) auxlink.y;
  link0z=(float) auxlink.z;
  link0w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(1+3*1));
  link1x=(float) auxlink.x;
  link1y=(float) auxlink.y;
  link1z=(float) auxlink.z;
  link1w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(2+3*1));
  link2x=(float) auxlink.x;
  link2y=(float) auxlink.y;
  link2z=(float) auxlink.z;
  link2w=(float) auxlink.w;
  

ferm_out[0][0][threadIdx.x] += link0x*ferm_in_0.x-link0y*ferm_in_0.y+  
                                 link0z*ferm_in_1.x-link0w*ferm_in_1.y+ 
                                 link1x*ferm_in_2.x-link1y*ferm_in_2.y; 
  ferm_out[0][1][threadIdx.x] += link0x*ferm_in_0.y+link0y*ferm_in_0.x+ 
                                 link0z*ferm_in_1.y+link0w*ferm_in_1.x+ 
                                 link1x*ferm_in_2.y+link1y*ferm_in_2.x; 

  ferm_out[1][0][threadIdx.x] += link1z*ferm_in_0.x-link1w*ferm_in_0.y+  
                                 link2x*ferm_in_1.x-link2y*ferm_in_1.y+ 
                                 link2z*ferm_in_2.x-link2w*ferm_in_2.y; 
  ferm_out[1][1][threadIdx.x] += link1z*ferm_in_0.y+link1w*ferm_in_0.x+ 
                                 link2x*ferm_in_1.y+link2y*ferm_in_1.x+ 
                                 link2z*ferm_in_2.y+link2w*ferm_in_2.x; 

  ferm_out[2][0][threadIdx.x] += stag_phase*(C1RED*ferm_in_0.x-C1IMD*ferm_in_0.y+  
					     C2RED*ferm_in_1.x-C2IMD*ferm_in_1.y+ 
					     C3RED*ferm_in_2.x-C3IMD*ferm_in_2.y); 
  ferm_out[2][1][threadIdx.x] += stag_phase*(C1RED*ferm_in_0.y+C1IMD*ferm_in_0.x+ 
					     C2RED*ferm_in_1.y+C2IMD*ferm_in_1.x+ 
					     C3RED*ferm_in_2.y+C3IMD*ferm_in_2.x); 


  //DIRECTION 2
  site_table[threadIdx.x] = tables[idx+6*size_dev];
  stag_phase              = (float) phases[idx+2*size_dev];

  ferm_in_0 = in[              site_table[threadIdx.x]];
  ferm_in_1 = in[   size_dev + site_table[threadIdx.x]];
  ferm_in_2 = in[ 2*size_dev + site_table[threadIdx.x]];


  // 1st float 
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(0+3*2));
  link0x=(float) auxlink.x;
  link0y=(float) auxlink.y;
  link0z=(float) auxlink.z;
  link0w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(1+3*2));
  link1x=(float) auxlink.x;
  link1y=(float) auxlink.y;
  link1z=(float) auxlink.z;
  link1w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(2+3*2));
  link2x=(float) auxlink.x;
  link2y=(float) auxlink.y;
  link2z=(float) auxlink.z;
  link2w=(float) auxlink.w;

  ferm_out[0][0][threadIdx.x] += link0x*ferm_in_0.x-link0y*ferm_in_0.y+  
                                 link0z*ferm_in_1.x-link0w*ferm_in_1.y+ 
                                 link1x*ferm_in_2.x-link1y*ferm_in_2.y; 
  ferm_out[0][1][threadIdx.x] += link0x*ferm_in_0.y+link0y*ferm_in_0.x+ 
                                 link0z*ferm_in_1.y+link0w*ferm_in_1.x+ 
                                 link1x*ferm_in_2.y+link1y*ferm_in_2.x; 

  ferm_out[1][0][threadIdx.x] += link1z*ferm_in_0.x-link1w*ferm_in_0.y+  
                                 link2x*ferm_in_1.x-link2y*ferm_in_1.y+ 
                                 link2z*ferm_in_2.x-link2w*ferm_in_2.y; 
  ferm_out[1][1][threadIdx.x] += link1z*ferm_in_0.y+link1w*ferm_in_0.x+ 
                                 link2x*ferm_in_1.y+link2y*ferm_in_1.x+ 
                                 link2z*ferm_in_2.y+link2w*ferm_in_2.x; 

  ferm_out[2][0][threadIdx.x] += stag_phase*(C1RED*ferm_in_0.x-C1IMD*ferm_in_0.y+  
					     C2RED*ferm_in_1.x-C2IMD*ferm_in_1.y+ 
					     C3RED*ferm_in_2.x-C3IMD*ferm_in_2.y); 
  ferm_out[2][1][threadIdx.x] += stag_phase*(C1RED*ferm_in_0.y+C1IMD*ferm_in_0.x+ 
					     C2RED*ferm_in_1.y+C2IMD*ferm_in_1.x+ 
					     C3RED*ferm_in_2.y+C3IMD*ferm_in_2.x); 



  //DIRECTION 3
  site_table[threadIdx.x]  = tables[idx+7*size_dev];
  stag_phase               = (float) phases[idx+3*size_dev];

  ferm_in_0 = in[              site_table[threadIdx.x]];
  ferm_in_1 = in[   size_dev + site_table[threadIdx.x]];
  ferm_in_2 = in[ 2*size_dev + site_table[threadIdx.x]];


  // 1st float 
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(0+3*3));
  link0x=(float) auxlink.x;
  link0y=(float) auxlink.y;
  link0z=(float) auxlink.z;
  link0w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(1+3*3));
  link1x=(float) auxlink.x;
  link1y=(float) auxlink.y;
  link1z=(float) auxlink.z;
  link1w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(2+3*3));
  link2x=(float) auxlink.x;
  link2y=(float) auxlink.y;
  link2z=(float) auxlink.z;
  link2w=(float) auxlink.w;

  ferm_out[0][0][threadIdx.x] += link0x*ferm_in_0.x-link0y*ferm_in_0.y+  
                                 link0z*ferm_in_1.x-link0w*ferm_in_1.y+ 
                                 link1x*ferm_in_2.x-link1y*ferm_in_2.y; 
  ferm_out[0][1][threadIdx.x] += link0x*ferm_in_0.y+link0y*ferm_in_0.x+ 
                                 link0z*ferm_in_1.y+link0w*ferm_in_1.x+ 
                                 link1x*ferm_in_2.y+link1y*ferm_in_2.x; 

  ferm_out[1][0][threadIdx.x] += link1z*ferm_in_0.x-link1w*ferm_in_0.y+  
                                 link2x*ferm_in_1.x-link2y*ferm_in_1.y+ 
                                 link2z*ferm_in_2.x-link2w*ferm_in_2.y; 
  ferm_out[1][1][threadIdx.x] += link1z*ferm_in_0.y+link1w*ferm_in_0.x+ 
                                 link2x*ferm_in_1.y+link2y*ferm_in_1.x+ 
                                 link2z*ferm_in_2.y+link2w*ferm_in_2.x; 

  ferm_out[2][0][threadIdx.x] += stag_phase*(C1RED*ferm_in_0.x-C1IMD*ferm_in_0.y+  
					     C2RED*ferm_in_1.x-C2IMD*ferm_in_1.y+ 
					     C3RED*ferm_in_2.x-C3IMD*ferm_in_2.y); 
  ferm_out[2][1][threadIdx.x] += stag_phase*(C1RED*ferm_in_0.y+C1IMD*ferm_in_0.x+ 
					     C2RED*ferm_in_1.y+C2IMD*ferm_in_1.x+ 
					     C3RED*ferm_in_2.y+C3IMD*ferm_in_2.x); 
  

  
  //---------------------------------------------------end of first block

  //DIRECTION 0
  site_table[threadIdx.x] = tables[idx];

  ferm_in_0 = in[              site_table[threadIdx.x]];
  ferm_in_1 = in[   size_dev + site_table[threadIdx.x]];
  ferm_in_2 = in[ 2*size_dev + site_table[threadIdx.x]];


  // 1st float 
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(0+3*0));
  link0x=(float) auxlink.x;
  link0y=(float) auxlink.y;
  link0z=(float) auxlink.z;
  link0w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(1+3*0));
  link1x=(float) auxlink.x;
  link1y=(float) auxlink.y;
  link1z=(float) auxlink.z;
  link1w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(2+3*0));
  link2x=(float) auxlink.x;
  link2y=(float) auxlink.y;
  link2z=(float) auxlink.z;
  link2w=(float) auxlink.w;

  ferm_out[0][0][threadIdx.x] -= link0x*ferm_in_0.x+link0y*ferm_in_0.y +
              			 link1z*ferm_in_1.x+link1w*ferm_in_1.y +
				 C1RED*ferm_in_2.x   +C1IMD*ferm_in_2.y; 
  
  ferm_out[0][1][threadIdx.x] -= link0x*ferm_in_0.y-link0y*ferm_in_0.x +
                                 link1z*ferm_in_1.y-link1w*ferm_in_1.x +
                                 C1RED*ferm_in_2.y   -C1IMD*ferm_in_2.x; 

  ferm_out[1][0][threadIdx.x] -= link0z*ferm_in_0.x+link0w*ferm_in_0.y +
                                 link2x*ferm_in_1.x+link2y*ferm_in_1.y +
                                 C2RED*ferm_in_2.x   +C2IMD*ferm_in_2.y; 

  ferm_out[1][1][threadIdx.x] -= link0z*ferm_in_0.y-link0w*ferm_in_0.x +
                                 link2x*ferm_in_1.y-link2y*ferm_in_1.x +
                                 C2RED*ferm_in_2.y   -C2IMD*ferm_in_2.x; 

  ferm_out[2][0][threadIdx.x] -= link1x*ferm_in_0.x+link1y*ferm_in_0.y +
                                 link2z*ferm_in_1.x+link2w*ferm_in_1.y +
                                 C3RED*ferm_in_2.x   +C3IMD*ferm_in_2.y; 

  ferm_out[2][1][threadIdx.x] -= link1x*ferm_in_0.y-link1y*ferm_in_0.x +
                                 link2z*ferm_in_1.y-link2w*ferm_in_1.x +
                                 C3RED*ferm_in_2.y   -C3IMD*ferm_in_2.x; 



  
  //DIRECTION 1
  site_table[threadIdx.x] = tables[idx+size_dev];
  stag_phase              = (float) phases[site_table[threadIdx.x]+size_dev];

  ferm_in_0 = in[              site_table[threadIdx.x]];
  ferm_in_1 = in[   size_dev + site_table[threadIdx.x]];
  ferm_in_2 = in[ 2*size_dev + site_table[threadIdx.x]];


  // 1st float 
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(0+3*1));
  link0x=(float) auxlink.x;
  link0y=(float) auxlink.y;
  link0z=(float) auxlink.z;
  link0w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(1+3*1));
  link1x=(float) auxlink.x;
  link1y=(float) auxlink.y;
  link1z=(float) auxlink.z;
  link1w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(2+3*1));
  link2x=(float) auxlink.x;
  link2y=(float) auxlink.y;
  link2z=(float) auxlink.z;
  link2w=(float) auxlink.w;


  ferm_out[0][0][threadIdx.x] -= link0x*ferm_in_0.x+link0y*ferm_in_0.y +
                                 link1z*ferm_in_1.x+link1w*ferm_in_1.y +
                                 stag_phase*(C1RED*ferm_in_2.x+C1IMD*ferm_in_2.y); 

  ferm_out[0][1][threadIdx.x] -= link0x*ferm_in_0.y-link0y*ferm_in_0.x +
                                 link1z*ferm_in_1.y-link1w*ferm_in_1.x +
                                 stag_phase*(C1RED*ferm_in_2.y-C1IMD*ferm_in_2.x); 

  ferm_out[1][0][threadIdx.x] -= link0z*ferm_in_0.x+link0w*ferm_in_0.y +
                                 link2x*ferm_in_1.x+link2y*ferm_in_1.y +
                                 stag_phase*(C2RED*ferm_in_2.x+C2IMD*ferm_in_2.y); 

  ferm_out[1][1][threadIdx.x] -= link0z*ferm_in_0.y-link0w*ferm_in_0.x +
                                 link2x*ferm_in_1.y-link2y*ferm_in_1.x +
                                 stag_phase*(C2RED*ferm_in_2.y-C2IMD*ferm_in_2.x); 

  ferm_out[2][0][threadIdx.x] -= link1x*ferm_in_0.x+link1y*ferm_in_0.y +
                                 link2z*ferm_in_1.x+link2w*ferm_in_1.y +
                                 stag_phase*(C3RED*ferm_in_2.x+C3IMD*ferm_in_2.y); 

  ferm_out[2][1][threadIdx.x] -= link1x*ferm_in_0.y-link1y*ferm_in_0.x +
                                 link2z*ferm_in_1.y-link2w*ferm_in_1.x +
                                 stag_phase*(C3RED*ferm_in_2.y- C3IMD*ferm_in_2.x); 



  //DIRECTION 2
  site_table[threadIdx.x] = tables[idx+2*size_dev];
  stag_phase              = (float) phases[site_table[threadIdx.x]+2*size_dev];

  ferm_in_0 = in[              site_table[threadIdx.x]];
  ferm_in_1 = in[   size_dev + site_table[threadIdx.x]];
  ferm_in_2 = in[ 2*size_dev + site_table[threadIdx.x]];

  // 1st float 
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(0+3*2));
  link0x=(float) auxlink.x;
  link0y=(float) auxlink.y;
  link0z=(float) auxlink.z;
  link0w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(1+3*2));
  link1x=(float) auxlink.x;
  link1y=(float) auxlink.y;
  link1z=(float) auxlink.z;
  link1w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(2+3*2));
  link2x=(float) auxlink.x;
  link2y=(float) auxlink.y;
  link2z=(float) auxlink.z;
  link2w=(float) auxlink.w;


  ferm_out[0][0][threadIdx.x] -= link0x*ferm_in_0.x+link0y*ferm_in_0.y +
                                 link1z*ferm_in_1.x+link1w*ferm_in_1.y +
                                 stag_phase*(C1RED*ferm_in_2.x+ C1IMD*ferm_in_2.y); 

  ferm_out[0][1][threadIdx.x] -= link0x*ferm_in_0.y-link0y*ferm_in_0.x +
                                 link1z*ferm_in_1.y-link1w*ferm_in_1.x +
                                 stag_phase*(C1RED*ferm_in_2.y- C1IMD*ferm_in_2.x); 

  ferm_out[1][0][threadIdx.x] -= link0z*ferm_in_0.x+link0w*ferm_in_0.y +
                                 link2x*ferm_in_1.x+link2y*ferm_in_1.y +
                                 stag_phase*(C2RED*ferm_in_2.x+ C2IMD*ferm_in_2.y); 

  ferm_out[1][1][threadIdx.x] -= link0z*ferm_in_0.y-link0w*ferm_in_0.x +
                                 link2x*ferm_in_1.y-link2y*ferm_in_1.x +
                                 stag_phase*(C2RED*ferm_in_2.y- C2IMD*ferm_in_2.x); 

  ferm_out[2][0][threadIdx.x] -= link1x*ferm_in_0.x+link1y*ferm_in_0.y +
                                 link2z*ferm_in_1.x+link2w*ferm_in_1.y +
                                 stag_phase*(C3RED*ferm_in_2.x+ C3IMD*ferm_in_2.y); 

  ferm_out[2][1][threadIdx.x] -= link1x*ferm_in_0.y-link1y*ferm_in_0.x +
                                 link2z*ferm_in_1.y-link2w*ferm_in_1.x +
                                 stag_phase*(C3RED*ferm_in_2.y- C3IMD*ferm_in_2.x); 



  //DIRECTION 3
  site_table[threadIdx.x] = tables[idx+3*size_dev];
  stag_phase              = (float) phases[site_table[threadIdx.x]+3*size_dev];

  ferm_in_0 = in[              site_table[threadIdx.x]];
  ferm_in_1 = in[   size_dev + site_table[threadIdx.x]];
  ferm_in_2 = in[ 2*size_dev + site_table[threadIdx.x]];


  // 1st float 
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(0+3*3));
  link0x=(float) auxlink.x;
  link0y=(float) auxlink.y;
  link0z=(float) auxlink.z;
  link0w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(1+3*3));
  link1x=(float) auxlink.x;
  link1y=(float) auxlink.y;
  link1z=(float) auxlink.z;
  link1w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(2+3*3));
  link2x=(float) auxlink.x;
  link2y=(float) auxlink.y;
  link2z=(float) auxlink.z;
  link2w=(float) auxlink.w;


  ferm_out[0][0][threadIdx.x] -= link0x*ferm_in_0.x+link0y*ferm_in_0.y +
                                 link1z*ferm_in_1.x+link1w*ferm_in_1.y +
                                 stag_phase*(C1RED*ferm_in_2.x+  C1IMD*ferm_in_2.y); 

  ferm_out[0][1][threadIdx.x] -= link0x*ferm_in_0.y-link0y*ferm_in_0.x +
                                 link1z*ferm_in_1.y-link1w*ferm_in_1.x +
                                 stag_phase*(C1RED*ferm_in_2.y- C1IMD*ferm_in_2.x); 

  ferm_out[1][0][threadIdx.x] -= link0z*ferm_in_0.x+link0w*ferm_in_0.y +
                                 link2x*ferm_in_1.x+link2y*ferm_in_1.y +
                                 stag_phase*(C2RED*ferm_in_2.x+ C2IMD*ferm_in_2.y); 

  ferm_out[1][1][threadIdx.x] -= link0z*ferm_in_0.y-link0w*ferm_in_0.x +
                                 link2x*ferm_in_1.y-link2y*ferm_in_1.x +
                                 stag_phase*(C2RED*ferm_in_2.y- C2IMD*ferm_in_2.x); 

  ferm_out[2][0][threadIdx.x] -= link1x*ferm_in_0.x+link1y*ferm_in_0.y +
                                 link2z*ferm_in_1.x+link2w*ferm_in_1.y +
                                 stag_phase*(C3RED*ferm_in_2.x+ C3IMD*ferm_in_2.y); 

  ferm_out[2][1][threadIdx.x] -= link1x*ferm_in_0.y-link1y*ferm_in_0.x +
                                 link2z*ferm_in_1.y-link2w*ferm_in_1.x +
                                 stag_phase*(C3RED*ferm_in_2.y- C3IMD*ferm_in_2.x); 
  
  //-------------------------------------------------end of second block

  // even
  ferm_in_0 = in[              idx - size_dev_h];
  ferm_in_1 = in[   size_dev + idx - size_dev_h];
  ferm_in_2 = in[ 2*size_dev + idx - size_dev_h];

//  out[idx              - size_dev_h ].x = 0 ; // mass_d_dev*ferm_in_0.x;
//  out[idx              - size_dev_h ].y = 0 ; // mass_d_dev*ferm_in_0.y;
//  out[idx +   size_dev - size_dev_h ].x = 0 ; // mass_d_dev*ferm_in_1.x;
//  out[idx +   size_dev - size_dev_h ].y = 0 ; // mass_d_dev*ferm_in_1.y;
//  out[idx + 2*size_dev - size_dev_h ].x = 0 ; // mass_d_dev*ferm_in_2.x;
//  out[idx + 2*size_dev - size_dev_h ].y = 0 ; // mass_d_dev*ferm_in_2.y;

  //odd
  out[idx               ].x =ferm_out[0][0][threadIdx.x]*(float)0.5;
  out[idx               ].y =ferm_out[0][1][threadIdx.x]*(float)0.5;
  out[idx +   size_dev  ].x =ferm_out[1][0][threadIdx.x]*(float)0.5;
  out[idx +   size_dev  ].y =ferm_out[1][1][threadIdx.x]*(float)0.5;
  out[idx + 2*size_dev  ].x =ferm_out[2][0][threadIdx.x]*(float)0.5;
  out[idx + 2*size_dev  ].y =ferm_out[2][1][threadIdx.x]*(float)0.5;

  //-------------------------------------------------end of Dslash
  }







__global__ void DslashDaggerKernelEO(float2 *out,
                                       float2 *in,
                                       int *tables, 
                                       int *phases,
                                       size_t gauge_offset) 
  { 
  int idx = blockIdx.x*blockDim.x + threadIdx.x;     // idx< sizeh, EVEN!!
  float stag_phase = 1.0;

  //Store result in sharedMem
  __shared__ float ferm_out[3][2][NUM_THREADS];
 
  //New tables indexing (index fastest)
  __shared__ int site_table[NUM_THREADS];

  //Load link matrix U_mu(ix) in registers
  float link0x, link0y, link0z, link0w, 
         link1x, link1y, link1z, link1w, 
         link2x, link2y, link2z, link2w;   
  float4 auxlink;

  float2 ferm_in_0, ferm_in_1, ferm_in_2;
  
  // DIRECTION 0
  site_table[threadIdx.x] = tables[idx+4*size_dev];

  ferm_in_0 = in[              site_table[threadIdx.x]];
  ferm_in_1 = in[   size_dev + site_table[threadIdx.x]];
  ferm_in_2 = in[ 2*size_dev + site_table[threadIdx.x]];



 
  // 1st float 
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(0+3*0));
  link0x=(float) auxlink.x;
  link0y=(float) auxlink.y;
  link0z=(float) auxlink.z;
  link0w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(1+3*0));
  link1x=(float) auxlink.x;
  link1y=(float) auxlink.y;
  link1z=(float) auxlink.z;
  link1w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(2+3*0));
  link2x=(float) auxlink.x;
  link2y=(float) auxlink.y;
  link2z=(float) auxlink.z;
  link2w=(float) auxlink.w;


  ferm_out[0][0][threadIdx.x] = link0x*ferm_in_0.x-link0y*ferm_in_0.y+  
                                link0z*ferm_in_1.x-link0w*ferm_in_1.y+ 
                                link1x*ferm_in_2.x-link1y*ferm_in_2.y; 
  ferm_out[0][1][threadIdx.x] = link0x*ferm_in_0.y+link0y*ferm_in_0.x+ 
                                link0z*ferm_in_1.y+link0w*ferm_in_1.x+ 
                                link1x*ferm_in_2.y+link1y*ferm_in_2.x; 

  ferm_out[1][0][threadIdx.x] = link1z*ferm_in_0.x-link1w*ferm_in_0.y+  
                                link2x*ferm_in_1.x-link2y*ferm_in_1.y+ 
                                link2z*ferm_in_2.x-link2w*ferm_in_2.y; 
  ferm_out[1][1][threadIdx.x] = link1z*ferm_in_0.y+link1w*ferm_in_0.x+ 
                                link2x*ferm_in_1.y+link2y*ferm_in_1.x+ 
                                link2z*ferm_in_2.y+link2w*ferm_in_2.x; 

  ferm_out[2][0][threadIdx.x] = C1RED*ferm_in_0.x-C1IMD*ferm_in_0.y+  
                                C2RED*ferm_in_1.x-C2IMD*ferm_in_1.y+ 
                                C3RED*ferm_in_2.x-C3IMD*ferm_in_2.y; 
  ferm_out[2][1][threadIdx.x] = C1RED*ferm_in_0.y+C1IMD*ferm_in_0.x+ 
                                C2RED*ferm_in_1.y+C2IMD*ferm_in_1.x+ 
                                C3RED*ferm_in_2.y+C3IMD*ferm_in_2.x; 



  //DIRECTION 1
  site_table[threadIdx.x] = tables[idx+5*size_dev];
  stag_phase              = (float) phases[idx+size_dev];

  ferm_in_0 = in[              site_table[threadIdx.x]];
  ferm_in_1 = in[   size_dev + site_table[threadIdx.x]];
  ferm_in_2 = in[ 2*size_dev + site_table[threadIdx.x]];


  // 1st float 
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(0+3*1));
  link0x=(float) auxlink.x;
  link0y=(float) auxlink.y;
  link0z=(float) auxlink.z;
  link0w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(1+3*1));
  link1x=(float) auxlink.x;
  link1y=(float) auxlink.y;
  link1z=(float) auxlink.z;
  link1w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(2+3*1));
  link2x=(float) auxlink.x;
  link2y=(float) auxlink.y;
  link2z=(float) auxlink.z;
  link2w=(float) auxlink.w;


  ferm_out[0][0][threadIdx.x] += link0x*ferm_in_0.x-link0y*ferm_in_0.y+  
                                 link0z*ferm_in_1.x-link0w*ferm_in_1.y+ 
                                 link1x*ferm_in_2.x-link1y*ferm_in_2.y; 
  ferm_out[0][1][threadIdx.x] += link0x*ferm_in_0.y+link0y*ferm_in_0.x+ 
                                 link0z*ferm_in_1.y+link0w*ferm_in_1.x+ 
                                 link1x*ferm_in_2.y+link1y*ferm_in_2.x; 

  ferm_out[1][0][threadIdx.x] += link1z*ferm_in_0.x-link1w*ferm_in_0.y+  
                                 link2x*ferm_in_1.x-link2y*ferm_in_1.y+ 
                                 link2z*ferm_in_2.x-link2w*ferm_in_2.y; 
  ferm_out[1][1][threadIdx.x] += link1z*ferm_in_0.y+link1w*ferm_in_0.x+ 
                                 link2x*ferm_in_1.y+link2y*ferm_in_1.x+ 
                                 link2z*ferm_in_2.y+link2w*ferm_in_2.x; 

  ferm_out[2][0][threadIdx.x] += stag_phase*(C1RED*ferm_in_0.x-C1IMD*ferm_in_0.y+  
					     C2RED*ferm_in_1.x-C2IMD*ferm_in_1.y+ 
					     C3RED*ferm_in_2.x-C3IMD*ferm_in_2.y); 
  ferm_out[2][1][threadIdx.x] += stag_phase*(C1RED*ferm_in_0.y+C1IMD*ferm_in_0.x+ 
					     C2RED*ferm_in_1.y+C2IMD*ferm_in_1.x+ 
					     C3RED*ferm_in_2.y+C3IMD*ferm_in_2.x); 
   


  //DIRECTION 2
  site_table[threadIdx.x] = tables[idx+6*size_dev];
  stag_phase              = (float) phases[idx+2*size_dev];

  ferm_in_0 = in[              site_table[threadIdx.x]];
  ferm_in_1 = in[   size_dev + site_table[threadIdx.x]];
  ferm_in_2 = in[ 2*size_dev + site_table[threadIdx.x]];

  // 1st float 
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(0+3*2));
  link0x=(float) auxlink.x;
  link0y=(float) auxlink.y;
  link0z=(float) auxlink.z;
  link0w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(1+3*2));
  link1x=(float) auxlink.x;
  link1y=(float) auxlink.y;
  link1z=(float) auxlink.z;
  link1w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(2+3*2));
  link2x=(float) auxlink.x;
  link2y=(float) auxlink.y;
  link2z=(float) auxlink.z;
  link2w=(float) auxlink.w;


  ferm_out[0][0][threadIdx.x] += link0x*ferm_in_0.x-link0y*ferm_in_0.y+  
                                 link0z*ferm_in_1.x-link0w*ferm_in_1.y+ 
                                 link1x*ferm_in_2.x-link1y*ferm_in_2.y; 
  ferm_out[0][1][threadIdx.x] += link0x*ferm_in_0.y+link0y*ferm_in_0.x+ 
                                 link0z*ferm_in_1.y+link0w*ferm_in_1.x+ 
                                 link1x*ferm_in_2.y+link1y*ferm_in_2.x; 

  ferm_out[1][0][threadIdx.x] += link1z*ferm_in_0.x-link1w*ferm_in_0.y+  
                                 link2x*ferm_in_1.x-link2y*ferm_in_1.y+ 
                                 link2z*ferm_in_2.x-link2w*ferm_in_2.y; 
  ferm_out[1][1][threadIdx.x] += link1z*ferm_in_0.y+link1w*ferm_in_0.x+ 
                                 link2x*ferm_in_1.y+link2y*ferm_in_1.x+ 
                                 link2z*ferm_in_2.y+link2w*ferm_in_2.x; 

  ferm_out[2][0][threadIdx.x] += stag_phase*(C1RED*ferm_in_0.x-C1IMD*ferm_in_0.y+  
					     C2RED*ferm_in_1.x-C2IMD*ferm_in_1.y+ 
					     C3RED*ferm_in_2.x-C3IMD*ferm_in_2.y); 
  ferm_out[2][1][threadIdx.x] += stag_phase*(C1RED*ferm_in_0.y+C1IMD*ferm_in_0.x+ 
					     C2RED*ferm_in_1.y+C2IMD*ferm_in_1.x+ 
					     C3RED*ferm_in_2.y+C3IMD*ferm_in_2.x); 

  
  //DIRECTION 3
  site_table[threadIdx.x] = tables[idx+7*size_dev];
   stag_phase              = (float) phases[idx+3*size_dev];

  ferm_in_0 = in[              site_table[threadIdx.x]];
  ferm_in_1 = in[   size_dev + site_table[threadIdx.x]];
  ferm_in_2 = in[ 2*size_dev + site_table[threadIdx.x]];


  // 1st float 
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(0+3*3));
  link0x=(float) auxlink.x;
  link0y=(float) auxlink.y;
  link0z=(float) auxlink.z;
  link0w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(1+3*3));
  link1x=(float) auxlink.x;
  link1y=(float) auxlink.y;
  link1z=(float) auxlink.z;
  link1w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, idx + gauge_offset + size_dev*(2+3*3));
  link2x=(float) auxlink.x;
  link2y=(float) auxlink.y;
  link2z=(float) auxlink.z;
  link2w=(float) auxlink.w;

  ferm_out[0][0][threadIdx.x] += link0x*ferm_in_0.x-link0y*ferm_in_0.y+  
                                 link0z*ferm_in_1.x-link0w*ferm_in_1.y+ 
                                 link1x*ferm_in_2.x-link1y*ferm_in_2.y; 
  ferm_out[0][1][threadIdx.x] += link0x*ferm_in_0.y+link0y*ferm_in_0.x+ 
                                 link0z*ferm_in_1.y+link0w*ferm_in_1.x+ 
                                 link1x*ferm_in_2.y+link1y*ferm_in_2.x; 

  ferm_out[1][0][threadIdx.x] += link1z*ferm_in_0.x-link1w*ferm_in_0.y+  
                                 link2x*ferm_in_1.x-link2y*ferm_in_1.y+ 
                                 link2z*ferm_in_2.x-link2w*ferm_in_2.y; 
  ferm_out[1][1][threadIdx.x] += link1z*ferm_in_0.y+link1w*ferm_in_0.x+ 
                                 link2x*ferm_in_1.y+link2y*ferm_in_1.x+ 
                                 link2z*ferm_in_2.y+link2w*ferm_in_2.x; 

  ferm_out[2][0][threadIdx.x] += stag_phase*(C1RED*ferm_in_0.x-C1IMD*ferm_in_0.y+  
					     C2RED*ferm_in_1.x-C2IMD*ferm_in_1.y+ 
					     C3RED*ferm_in_2.x-C3IMD*ferm_in_2.y); 
  ferm_out[2][1][threadIdx.x] += stag_phase*(C1RED*ferm_in_0.y+C1IMD*ferm_in_0.x+ 
					     C2RED*ferm_in_1.y+C2IMD*ferm_in_1.x+ 
					     C3RED*ferm_in_2.y+C3IMD*ferm_in_2.x); 


  //---------------------------------------------------end of first block
 
  //DIRECTION 0
  site_table[threadIdx.x] = tables[idx];
 
  ferm_in_0 = in[              site_table[threadIdx.x]];
  ferm_in_1 = in[   size_dev + site_table[threadIdx.x]];
  ferm_in_2 = in[ 2*size_dev + site_table[threadIdx.x]];


  // 1st float 
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(0+3*0));
  link0x=(float) auxlink.x;
  link0y=(float) auxlink.y;
  link0z=(float) auxlink.z;
  link0w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(1+3*0));
  link1x=(float) auxlink.x;
  link1y=(float) auxlink.y;
  link1z=(float) auxlink.z;
  link1w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(2+3*0));
  link2x=(float) auxlink.x;
  link2y=(float) auxlink.y;
  link2z=(float) auxlink.z;
  link2w=(float) auxlink.w;


  ferm_out[0][0][threadIdx.x] -= link0x*ferm_in_0.x+link0y*ferm_in_0.y +
              			 link1z*ferm_in_1.x+link1w*ferm_in_1.y +
				 C1RED*ferm_in_2.x   +C1IMD*ferm_in_2.y; 
  
  ferm_out[0][1][threadIdx.x] -= link0x*ferm_in_0.y-link0y*ferm_in_0.x +
                                 link1z*ferm_in_1.y-link1w*ferm_in_1.x +
                                 C1RED*ferm_in_2.y   -C1IMD*ferm_in_2.x; 

  ferm_out[1][0][threadIdx.x] -= link0z*ferm_in_0.x+link0w*ferm_in_0.y +
                                 link2x*ferm_in_1.x+link2y*ferm_in_1.y +
                                 C2RED*ferm_in_2.x   +C2IMD*ferm_in_2.y; 

  ferm_out[1][1][threadIdx.x] -= link0z*ferm_in_0.y-link0w*ferm_in_0.x +
                                 link2x*ferm_in_1.y-link2y*ferm_in_1.x +
                                 C2RED*ferm_in_2.y   -C2IMD*ferm_in_2.x; 

  ferm_out[2][0][threadIdx.x] -= link1x*ferm_in_0.x+link1y*ferm_in_0.y +
                                 link2z*ferm_in_1.x+link2w*ferm_in_1.y +
                                 C3RED*ferm_in_2.x   +C3IMD*ferm_in_2.y; 

  ferm_out[2][1][threadIdx.x] -= link1x*ferm_in_0.y-link1y*ferm_in_0.x +
                                 link2z*ferm_in_1.y-link2w*ferm_in_1.x +
                                 C3RED*ferm_in_2.y   -C3IMD*ferm_in_2.x; 
  

  //DIRECTION 1
  site_table[threadIdx.x] = tables[idx+size_dev];
  stag_phase              = (float) phases[site_table[threadIdx.x]+size_dev];

  ferm_in_0 = in[              site_table[threadIdx.x]];
  ferm_in_1 = in[   size_dev + site_table[threadIdx.x]];
  ferm_in_2 = in[ 2*size_dev + site_table[threadIdx.x]];

  // 1st float 
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(0+3*1));
  link0x=(float) auxlink.x;
  link0y=(float) auxlink.y;
  link0z=(float) auxlink.z;
  link0w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(1+3*1));
  link1x=(float) auxlink.x;
  link1y=(float) auxlink.y;
  link1z=(float) auxlink.z;
  link1w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(2+3*1));
  link2x=(float) auxlink.x;
  link2y=(float) auxlink.y;
  link2z=(float) auxlink.z;
  link2w=(float) auxlink.w;


  ferm_out[0][0][threadIdx.x] -= link0x*ferm_in_0.x+link0y*ferm_in_0.y +
                                 link1z*ferm_in_1.x+link1w*ferm_in_1.y +
                                 stag_phase*(C1RED*ferm_in_2.x+C1IMD*ferm_in_2.y); 

  ferm_out[0][1][threadIdx.x] -= link0x*ferm_in_0.y-link0y*ferm_in_0.x +
                                 link1z*ferm_in_1.y-link1w*ferm_in_1.x +
                                 stag_phase*(C1RED*ferm_in_2.y-C1IMD*ferm_in_2.x); 

  ferm_out[1][0][threadIdx.x] -= link0z*ferm_in_0.x+link0w*ferm_in_0.y +
                                 link2x*ferm_in_1.x+link2y*ferm_in_1.y +
                                 stag_phase*(C2RED*ferm_in_2.x+C2IMD*ferm_in_2.y); 

  ferm_out[1][1][threadIdx.x] -= link0z*ferm_in_0.y-link0w*ferm_in_0.x +
                                 link2x*ferm_in_1.y-link2y*ferm_in_1.x +
                                 stag_phase*(C2RED*ferm_in_2.y-C2IMD*ferm_in_2.x); 

  ferm_out[2][0][threadIdx.x] -= link1x*ferm_in_0.x+link1y*ferm_in_0.y +
                                 link2z*ferm_in_1.x+link2w*ferm_in_1.y +
                                 stag_phase*(C3RED*ferm_in_2.x+C3IMD*ferm_in_2.y); 

  ferm_out[2][1][threadIdx.x] -= link1x*ferm_in_0.y-link1y*ferm_in_0.x +
                                 link2z*ferm_in_1.y-link2w*ferm_in_1.x +
                                 stag_phase*(C3RED*ferm_in_2.y- C3IMD*ferm_in_2.x); 

  //DIRECTION 2
  site_table[threadIdx.x] = tables[idx+2*size_dev];
  stag_phase              = (float) phases[site_table[threadIdx.x]+2*size_dev];

  ferm_in_0 = in[              site_table[threadIdx.x]];
  ferm_in_1 = in[   size_dev + site_table[threadIdx.x]];
  ferm_in_2 = in[ 2*size_dev + site_table[threadIdx.x]];

  // 1st float
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(0+3*2));
  link0x=(float) auxlink.x;
  link0y=(float) auxlink.y;
  link0z=(float) auxlink.z;
  link0w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(1+3*2));
  link1x=(float) auxlink.x;
  link1y=(float) auxlink.y;
  link1z=(float) auxlink.z;
  link1w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(2+3*2));
  link2x=(float) auxlink.x;
  link2y=(float) auxlink.y;
  link2z=(float) auxlink.z;
  link2w=(float) auxlink.w;


  ferm_out[0][0][threadIdx.x] -= link0x*ferm_in_0.x+link0y*ferm_in_0.y +
                                 link1z*ferm_in_1.x+link1w*ferm_in_1.y +
                                 stag_phase*(C1RED*ferm_in_2.x+ C1IMD*ferm_in_2.y); 

  ferm_out[0][1][threadIdx.x] -= link0x*ferm_in_0.y-link0y*ferm_in_0.x +
                                 link1z*ferm_in_1.y-link1w*ferm_in_1.x +
                                 stag_phase*(C1RED*ferm_in_2.y- C1IMD*ferm_in_2.x); 

  ferm_out[1][0][threadIdx.x] -= link0z*ferm_in_0.x+link0w*ferm_in_0.y +
                                 link2x*ferm_in_1.x+link2y*ferm_in_1.y +
                                 stag_phase*(C2RED*ferm_in_2.x+ C2IMD*ferm_in_2.y); 

  ferm_out[1][1][threadIdx.x] -= link0z*ferm_in_0.y-link0w*ferm_in_0.x +
                                 link2x*ferm_in_1.y-link2y*ferm_in_1.x +
                                 stag_phase*(C2RED*ferm_in_2.y- C2IMD*ferm_in_2.x); 

  ferm_out[2][0][threadIdx.x] -= link1x*ferm_in_0.x+link1y*ferm_in_0.y +
                                 link2z*ferm_in_1.x+link2w*ferm_in_1.y +
                                 stag_phase*(C3RED*ferm_in_2.x+ C3IMD*ferm_in_2.y); 

  ferm_out[2][1][threadIdx.x] -= link1x*ferm_in_0.y-link1y*ferm_in_0.x +
                                 link2z*ferm_in_1.y-link2w*ferm_in_1.x +
                                 stag_phase*(C3RED*ferm_in_2.y- C3IMD*ferm_in_2.x); 



  //DIRECTION 3

 site_table[threadIdx.x] = tables[idx+3*size_dev];
  stag_phase              = (float) phases[site_table[threadIdx.x]+3*size_dev];

  ferm_in_0 = in[              site_table[threadIdx.x]];
  ferm_in_1 = in[   size_dev + site_table[threadIdx.x]];
  ferm_in_2 = in[ 2*size_dev + site_table[threadIdx.x]];

  // 1st float
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(0+3*3));
  link0x=(float) auxlink.x;
  link0y=(float) auxlink.y;
  link0z=(float) auxlink.z;
  link0w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(1+3*3));
  link1x=(float) auxlink.x;
  link1y=(float) auxlink.y;
  link1z=(float) auxlink.z;
  link1w=(float) auxlink.w;
  auxlink = tex1Dfetch(gauge_texRef, site_table[threadIdx.x] + gauge_offset + size_dev*(2+3*3));
  link2x=(float) auxlink.x;
  link2y=(float) auxlink.y;
  link2z=(float) auxlink.z;
  link2w=(float) auxlink.w;

  ferm_out[0][0][threadIdx.x] -= link0x*ferm_in_0.x+link0y*ferm_in_0.y +
                                 link1z*ferm_in_1.x+link1w*ferm_in_1.y +
                                 stag_phase*(C1RED*ferm_in_2.x+  C1IMD*ferm_in_2.y); 

  ferm_out[0][1][threadIdx.x] -= link0x*ferm_in_0.y-link0y*ferm_in_0.x +
                                 link1z*ferm_in_1.y-link1w*ferm_in_1.x +
                                 stag_phase*(C1RED*ferm_in_2.y- C1IMD*ferm_in_2.x); 

  ferm_out[1][0][threadIdx.x] -= link0z*ferm_in_0.x+link0w*ferm_in_0.y +
                                 link2x*ferm_in_1.x+link2y*ferm_in_1.y +
                                 stag_phase*(C2RED*ferm_in_2.x+ C2IMD*ferm_in_2.y); 

  ferm_out[1][1][threadIdx.x] -= link0z*ferm_in_0.y-link0w*ferm_in_0.x +
                                 link2x*ferm_in_1.y-link2y*ferm_in_1.x +
                                 stag_phase*(C2RED*ferm_in_2.y- C2IMD*ferm_in_2.x); 

  ferm_out[2][0][threadIdx.x] -= link1x*ferm_in_0.x+link1y*ferm_in_0.y +
                                 link2z*ferm_in_1.x+link2w*ferm_in_1.y +
                                 stag_phase*(C3RED*ferm_in_2.x+ C3IMD*ferm_in_2.y); 

  ferm_out[2][1][threadIdx.x] -= link1x*ferm_in_0.y-link1y*ferm_in_0.x +
                                 link2z*ferm_in_1.y-link2w*ferm_in_1.x +
                                 stag_phase*(C3RED*ferm_in_2.y- C3IMD*ferm_in_2.x); 

  //-------------------------------------------------end of second block

  // even   
  ferm_in_0 = in[              idx];
  ferm_in_1 = in[   size_dev + idx];
  ferm_in_2 = in[ 2*size_dev + idx];
/*
  out[idx               ].x = mass_d_dev*ferm_in_0.x - ferm_out[0][0][threadIdx.x]*(float)0.5;
  out[idx               ].y = mass_d_dev*ferm_in_0.y - ferm_out[0][1][threadIdx.x]*(float)0.5;
  out[idx +   size_dev  ].x = mass_d_dev*ferm_in_1.x - ferm_out[1][0][threadIdx.x]*(float)0.5;
  out[idx +   size_dev  ].y = mass_d_dev*ferm_in_1.y - ferm_out[1][1][threadIdx.x]*(float)0.5;
  out[idx + 2*size_dev  ].x = mass_d_dev*ferm_in_2.x - ferm_out[2][0][threadIdx.x]*(float)0.5;
  out[idx + 2*size_dev  ].y = mass_d_dev*ferm_in_2.y - ferm_out[2][1][threadIdx.x]*(float)0.5;
*/
  out[idx               ].x =ferm_out[0][0][threadIdx.x]*(float)0.5;
  out[idx               ].y =ferm_out[0][1][threadIdx.x]*(float)0.5;
  out[idx +   size_dev  ].x =ferm_out[1][0][threadIdx.x]*(float)0.5;
  out[idx +   size_dev  ].y =ferm_out[1][1][threadIdx.x]*(float)0.5;
  out[idx + 2*size_dev  ].x =ferm_out[2][0][threadIdx.x]*(float)0.5;
  out[idx + 2*size_dev  ].y =ferm_out[2][1][threadIdx.x]*(float)0.5;


  // odd
  out[idx              + size_dev_h ].x = (float)0.0;
  out[idx              + size_dev_h ].y = (float)0.0;
  out[idx +   size_dev + size_dev_h ].x = (float)0.0;
  out[idx +   size_dev + size_dev_h ].y = (float)0.0;
  out[idx + 2*size_dev + size_dev_h ].x = (float)0.0;
  out[idx + 2*size_dev + size_dev_h ].y = (float)0.0;

  //-------------------------------------------------end of DslashDagger
  }






/*
================================================================= EXTERNAL C FUNCTION
*/

void switchParity(float2 *out, float2* in)
{

  dim3 BlockDimension(NUM_THREADS);
  dim3 GridDimension(sizeh/BlockDimension.x);  //Half sites

  switchParityKernel<<<GridDimension,BlockDimension>>>(out,in);


}
void copy(float2 *out, float2* in)
{

  dim3 BlockDimension(NUM_THREADS);
  dim3 GridDimension(sizeh/BlockDimension.x);  //Half sites

  copyKernel<<<GridDimension,BlockDimension>>>(out,in);


}

void DslashOperatorEO(float2 *out, 
 		        float2 *in, 
 		        const int isign)
  {
  #ifdef DEBUG_MODE_2
  printf("\033[32mDEBUG: inside DslashOperatorDDEO ...\033[0m\n");
  #endif

  dim3 BlockDimension(NUM_THREADS);
  dim3 GridDimension(sizeh/BlockDimension.x);  //Half sites

  size_t gauge_field_size = sizeof(float4)*size*12;

  size_t offset_g;
  cudaSafe(AT,cudaBindTexture(&offset_g, gauge_texRef, gauge_field_device, 2*gauge_field_size), "cudaBindTexture");  
  offset_g/=sizeof(float4);

  if(isign == PLUS) 
    {
    DslashKernelEO<<<GridDimension,BlockDimension>>>(out, in, device_table, device_phases, offset_g); 
    cudaCheckError(AT,"DslashDDKernelEO"); 
    }
  
  if(isign == MINUS) 
    {
    DslashDaggerKernelEO<<<GridDimension,BlockDimension>>>(out, in, device_table, device_phases, offset_g); 
    cudaCheckError(AT,"DslashDaggerDDKernelEO"); 
    }

  cudaSafe(AT,cudaUnbindTexture(gauge_texRef), "cudaUnbindTexture");

  #ifdef DEBUG_MODE_2
  printf("\033[32m\tterminated DslashOperatorDDEO \033[0m\n");
  #endif
  }

