
#ifndef PACKER_CC_
#define PACKER_CC_

#include <complex>

// configuration packer
// 12 = 2(row)*3(complex)*2(real)
void smartpack_gauge(float out[2*12*no_links] , const Conf *in)
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside smartpack_gauge ..."<<endl;
 #endif

 long int i, mu;
 long int offs=12*no_links;
 double aux;

 for(mu=0; mu<4; mu++)
    {
    for(i=0; i<size; i++)
       {
       // 1st float
       out[4*i              + 12*mu*size] = (float) real(in->u_work[i + mu*size].comp[0][0]);
       out[4*i + 1          + 12*mu*size] = (float) imag(in->u_work[i + mu*size].comp[0][0]);
       out[4*i + 2          + 12*mu*size] = (float) real(in->u_work[i + mu*size].comp[0][1]);
       out[4*i + 3          + 12*mu*size] = (float) imag(in->u_work[i + mu*size].comp[0][1]);
	 
       out[4*i     + 4*size + 12*mu*size] = (float) real(in->u_work[i + mu*size].comp[0][2]);
       out[4*i + 1 + 4*size + 12*mu*size] = (float) imag(in->u_work[i + mu*size].comp[0][2]);
       out[4*i + 2 + 4*size + 12*mu*size] = (float) real(in->u_work[i + mu*size].comp[1][0]);
       out[4*i + 3 + 4*size + 12*mu*size] = (float) imag(in->u_work[i + mu*size].comp[1][0]);
	 
       out[4*i     + 8*size + 12*mu*size] = (float) real(in->u_work[i + mu*size].comp[1][1]);
       out[4*i + 1 + 8*size + 12*mu*size] = (float) imag(in->u_work[i + mu*size].comp[1][1]);
       out[4*i + 2 + 8*size + 12*mu*size] = (float) real(in->u_work[i + mu*size].comp[1][2]);
       out[4*i + 3 + 8*size + 12*mu*size] = (float) imag(in->u_work[i + mu*size].comp[1][2]);

       // 2nd float
       aux = real(in->u_work[i + mu*size].comp[0][0]) - (double) out[4*i              + 12*mu*size];
       out[offs + 4*i              + 12*mu*size] = (float) aux;
       aux = imag(in->u_work[i + mu*size].comp[0][0]) - (double) out[4*i + 1          + 12*mu*size]; 
       out[offs + 4*i + 1          + 12*mu*size] = (float) aux;
       aux = real(in->u_work[i + mu*size].comp[0][1]) - (double) out[4*i + 2          + 12*mu*size];
       out[offs + 4*i + 2          + 12*mu*size] = (float) aux;
       aux = imag(in->u_work[i + mu*size].comp[0][1]) - (double) out[4*i + 3          + 12*mu*size];
       out[offs + 4*i + 3          + 12*mu*size] = (float) aux;

       aux = real(in->u_work[i + mu*size].comp[0][2]) - (double) out[4*i     + 4*size + 12*mu*size];
       out[offs + 4*i     + 4*size + 12*mu*size] = (float) aux;
       aux = imag(in->u_work[i + mu*size].comp[0][2]) - (double) out[4*i + 1 + 4*size + 12*mu*size];
       out[offs + 4*i + 1 + 4*size + 12*mu*size] = (float) aux;
       aux = real(in->u_work[i + mu*size].comp[1][0]) - (double) out[4*i + 2 + 4*size + 12*mu*size];
       out[offs + 4*i + 2 + 4*size + 12*mu*size] = (float) aux;
       aux = imag(in->u_work[i + mu*size].comp[1][0]) - (double) out[4*i + 3 + 4*size + 12*mu*size];
       out[offs + 4*i + 3 + 4*size + 12*mu*size] = (float) aux;

       aux = real(in->u_work[i + mu*size].comp[1][1]) - (double) out[4*i     + 8*size + 12*mu*size];
       out[offs + 4*i     + 8*size + 12*mu*size] = (float) aux;
       aux = imag(in->u_work[i + mu*size].comp[1][1]) - (double) out[4*i + 1 + 8*size + 12*mu*size];
       out[offs + 4*i + 1 + 8*size + 12*mu*size] = (float) aux;
       aux = real(in->u_work[i + mu*size].comp[1][2]) - (double) out[4*i + 2 + 8*size + 12*mu*size];
       out[offs + 4*i + 2 + 8*size + 12*mu*size] = (float) aux;
       aux = imag(in->u_work[i + mu*size].comp[1][2]) - (double) out[4*i + 3 + 8*size + 12*mu*size];
       out[offs + 4*i + 3 + 8*size + 12*mu*size] = (float) aux;
       }
    }

 #ifdef DEBUG_MODE
 cout << "\tterminated smartpack_gauge"<<endl;
 #endif
 }


// Unpack gauge configuration
void smartunpack_gauge(Conf *out, const float in[2*12*no_links])
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside smartunpack_gauge ..."<<endl;
 #endif

 long int i;
 int mu;
 long int offs=12*no_links;
 REAL tempR, tempI;
 complex<REAL> aux[3][3]; 

 for(mu=0; mu<4; mu++)
    {
    for(i=0; i<size;i++)
       {
       tempR  = (double) in[       4*i              + 12*mu*size];
       tempR += (double) in[offs + 4*i              + 12*mu*size];
       tempI  = (double) in[       4*i + 1          + 12*mu*size];
       tempI += (double) in[offs + 4*i + 1          + 12*mu*size];
       aux[0][0]=complex<REAL>(tempR, tempI);

       tempR  = (double) in[       4*i + 2          + 12*mu*size];
       tempR += (double) in[offs + 4*i + 2          + 12*mu*size];
       tempI  = (double) in[       4*i + 3          + 12*mu*size];
       tempI += (double) in[offs + 4*i + 3          + 12*mu*size];
       aux[0][1]=complex<REAL>(tempR, tempI);

       tempR  = (double) in[       4*i     + 4*size + 12*mu*size];
       tempR += (double) in[offs + 4*i     + 4*size + 12*mu*size];
       tempI  = (double) in[       4*i + 1 + 4*size + 12*mu*size];
       tempI += (double) in[offs + 4*i + 1 + 4*size + 12*mu*size];
       aux[0][2]=complex<REAL>(tempR, tempI);

       tempR  = (double) in[       4*i + 2 + 4*size + 12*mu*size];
       tempR += (double) in[offs + 4*i + 2 + 4*size + 12*mu*size];
       tempI  = (double) in[       4*i + 3 + 4*size + 12*mu*size];
       tempI += (double) in[offs + 4*i + 3 + 4*size + 12*mu*size];
       aux[1][0]=complex<REAL>(tempR, tempI);

       tempR  = (double) in[       4*i     + 8*size + 12*mu*size];
       tempR += (double) in[offs + 4*i     + 8*size + 12*mu*size];
       tempI  = (double) in[       4*i + 1 + 8*size + 12*mu*size];
       tempI += (double) in[offs + 4*i + 1 + 8*size + 12*mu*size];
       aux[1][1]=complex<REAL>(tempR, tempI);

       tempR  = (double) in[       4*i + 2 + 8*size + 12*mu*size];
       tempR += (double) in[offs + 4*i + 2 + 8*size + 12*mu*size];
       tempI  = (double) in[       4*i + 3 + 8*size + 12*mu*size];
       tempI += (double) in[offs + 4*i + 3 + 8*size + 12*mu*size];
       aux[1][2]=complex<REAL>(tempR, tempI);

       aux[2][0]=complex<REAL>(0.0,0.0);
       aux[2][1]=complex<REAL>(0.0,0.0);
       aux[2][2]=complex<REAL>(0.0,0.0);

       out->u_work[i+mu*size]=Su3(aux);      
       }
    }
 out->unitarize_with_eta();

 #ifdef DEBUG_MODE
 cout << "\tterminated smartunpack_gauge"<<endl;
 #endif
 }




// fermion packer
//                                                            size
//                                                           |----|        v2
// Fermion are only on even sites on CPU                               |^^^^^^^^^|
// are on even and odd sites on GPU                            e    o    e    o    e    o
//                  |v1|                      fermion_device |----|----|----|----|----|----|
//  Fermion(site_i)=|v2|                                     |_________|         |_________|
//                  |v3|                                          v1                 v3
//
//  |--v1--|=|v1(site0).re,v1(site0).im,v1(site1).re,v1(site1).im-----v1(sizeh).re, v1(sizeh).im|
//  and similarly for the other components 
 
void smartpack_fermion(float out[6*sizeh], const Fermion *in)
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside smartpack_fermion ..."<<endl;
 #endif

 long int i;

 for(i=0; i<sizeh; i++)
    {
    out[2*i            ] = (float) real(in->fermion[i].comp[0]);
    out[2*i          +1] = (float) imag(in->fermion[i].comp[0]);
    out[2*i +   size   ] = (float) real(in->fermion[i].comp[1]);
    out[2*i +   size +1] = (float) imag(in->fermion[i].comp[1]);
    out[2*i + 2*size   ] = (float) real(in->fermion[i].comp[2]);
    out[2*i + 2*size +1] = (float) imag(in->fermion[i].comp[2]);
    }

 #ifdef DEBUG_MODE
 cout << "\tterminated smartpack_fermion"<<endl;
 #endif
 }


void smartpack_fermion_d(float out[6*sizeh*2], const Fermion *in)
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside smartpack_fermion_d ..."<<endl;
 #endif

 long int i;
 long int offs=6*sizeh;
 double aux;

 for(i=0; i<sizeh; i++)
    {
    // 1st float
    out[2*i            ] = (float) real(in->fermion[i].comp[0]);
    out[2*i          +1] = (float) imag(in->fermion[i].comp[0]);
    out[2*i +   size   ] = (float) real(in->fermion[i].comp[1]);
    out[2*i +   size +1] = (float) imag(in->fermion[i].comp[1]);
    out[2*i + 2*size   ] = (float) real(in->fermion[i].comp[2]);
    out[2*i + 2*size +1] = (float) imag(in->fermion[i].comp[2]);

    //2nd float
    aux = real(in->fermion[i].comp[0]) - (double)out[2*i            ];
    out[offs + 2*i            ] = (float) aux;
    aux = imag(in->fermion[i].comp[0]) - (double)out[2*i          +1];
    out[offs + 2*i          +1] = (float) aux;
    aux = real(in->fermion[i].comp[1]) - (double)out[2*i +   size   ];
    out[offs + 2*i +   size   ] = (float) aux;
    aux = imag(in->fermion[i].comp[1]) - (double)out[2*i +   size +1];
    out[offs + 2*i +   size +1] = (float) aux;
    aux = real(in->fermion[i].comp[2]) - (double)out[2*i + 2*size   ];
    out[offs + 2*i + 2*size   ] = (float) aux;
    aux = imag(in->fermion[i].comp[2]) - (double)out[2*i + 2*size +1];
    out[offs + 2*i + 2*size +1] = (float) aux;
    }

 #ifdef DEBUG_MODE
 cout << "\tterminated smartpack_fermion_d"<<endl;
 #endif
 }

/*
void smartpack_fermion(float *out, const Fermion *in)
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside smartpack_fermion_d ..."<<endl;
 #endif

 long int i;
 double aux;

 for(i=0; i<sizeh; i++)
    {
    // 1st float
    out[2*i            ] = (float) real(in->fermion[i].comp[0]);
    out[2*i          +1] = (float) imag(in->fermion[i].comp[0]);
    out[2*i +   size   ] = (float) real(in->fermion[i].comp[1]);
    out[2*i +   size +1] = (float) imag(in->fermion[i].comp[1]);
    out[2*i + 2*size   ] = (float) real(in->fermion[i].comp[2]);
    out[2*i + 2*size +1] = (float) imag(in->fermion[i].comp[2]);

 }

 #ifdef DEBUG_MODE
 cout << "\tterminated smartpack_fermion"<<endl;
 #endif
 }
*/


void smartunpack_fermion_d(Fermion *out, const float in[6*sizeh*2])
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside smartunpack_fermion_d ..."<<endl;
 #endif

 long int i;
 long int offs=6*sizeh;

 double d_re, d_im;
 complex<REAL> aux[3];

 for(i=0; i<sizeh; i++)
    {
    d_re = (double) in[2*i            ];
    d_re+= (double) in[2*i             +offs];
    d_im = (double) in[2*i          +1];
    d_im+= (double) in[2*i          +1 +offs];
    aux[0]=complex<REAL>(d_re, d_im);

    d_re = (double) in[2*i +   size   ];
    d_re+= (double) in[2*i +   size    +offs];
    d_im = (double) in[2*i +   size +1];
    d_im+= (double) in[2*i +   size +1 +offs];
    aux[1]=complex<REAL>(d_re, d_im);

    d_re = (double) in[2*i + 2*size   ];
    d_re+= (double) in[2*i + 2*size    +offs];
    d_im = (double) in[2*i + 2*size +1];
    d_im+= (double) in[2*i + 2*size +1 +offs];
    aux[2]=complex<REAL>(d_re, d_im);

    out->fermion[i]=Vec3(aux);
    }

 #ifdef DEBUG_MODE
 cout << "\tterminated smartunpack_fermion_d"<<endl;
 #endif
 }

void smartunpack_fermion(Fermion *out, const float *in)
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside smartunpack_fermion_d ..."<<endl;
 #endif

 long int i;

 double d_re, d_im;
 complex<REAL> aux[3];

 for(i=0; i<sizeh; i++)
    {
    d_re = (double) in[2*i            ];
    d_im = (double) in[2*i          +1];
    aux[0]=complex<REAL>(d_re, d_im);

    d_re = (double) in[2*i +   size   ];
    d_im = (double) in[2*i +   size +1];
    aux[1]=complex<REAL>(d_re, d_im);

    d_re = (double) in[2*i + 2*size   ];
    d_im = (double) in[2*i + 2*size +1];
    aux[2]=complex<REAL>(d_re, d_im);

    out->fermion[i]=Vec3(aux);
    }

 #ifdef DEBUG_MODE
 cout << "\tterminated smartunpack_fermion_d"<<endl;
 #endif
 }


// create shift table
void make_shift_table(int table[8*size])
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside make_shift_table ..."<<endl;
 #endif

 long int i;

 for(i=0; i<size; i++)
    {
    table[i         ]= nnm[i][0];
    table[i +   size]= nnm[i][1];
    table[i + 2*size]= nnm[i][2];
    table[i + 3*size]= nnm[i][3];

    table[i + 4*size]= nnp[i][0];
    table[i + 5*size]= nnp[i][1];
    table[i + 6*size]= nnp[i][2];
    table[i + 7*size]= nnp[i][3];
    }

 #ifdef DEBUG_MODE
 cout << "\tterminated make_shift_table"<<endl;
 #endif
 }



#ifdef __CUDACC__ 


void smartPackFermionOnDevice(float2* outDevice, Fermion* inHost){
    

    float* TempPacked = new float[2*sizeh];


    // 1st float
    for(int color = 0; color <3 ; color++){
        for(int i =0; i < sizeh ; i++){
            //1st float
            TempPacked[2*i] = (float) real(inHost->fermion[i].comp[color]);
            TempPacked[2*i+1] = (float) imag(inHost->fermion[i].comp[color]);
        }
        cudaSafe(AT, cudaMemcpy(outDevice+color*2*sizeh,TempPacked,2*sizeh*sizeof(float),cudaMemcpyHostToDevice ),"cudaMemcpy");
    }
    
    delete[] TempPacked;
}

void smartPackFermionOnDeviceD(float2* outDevice, Fermion* inHost){
    
    float* TempPacked = new float[3*2*sizeh];

    for(int color = 0; color <3 ; color++){
        for(int i =0; i < sizeh ; i++){
            //1st float
            TempPacked[color*2*sizeh + 2*i  ] = (float) real(inHost->fermion[i].comp[color]);
            TempPacked[color*2*sizeh + 2*i+1] = (float) imag(inHost->fermion[i].comp[color]);
        }
        cudaSafe(AT, cudaMemcpy(outDevice+color*2*sizeh,TempPacked+color*2*sizeh,2*sizeh*sizeof(float), cudaMemcpyHostToDevice ),"cudaMemcpy");
//        cudaSafe(AT, cudaMemset(outDevice+color*2*sizeh,0,2*sizeh*sizeof(float)),"cudaMemcpy");
    }

    for(int color = 0; color <3 ; color++){
        for(int i =0; i < sizeh ; i++){
            //2nd float
            TempPacked[color*2*sizeh + 2*i  ] = (float) ( real(inHost->fermion[i].comp[color]) - (double) TempPacked[color*2*sizeh + 2*i  ] );
            TempPacked[color*2*sizeh + 2*i+1] = (float) ( imag(inHost->fermion[i].comp[color]) - (double) TempPacked[color*2*sizeh + 2*i+1] );
        }
        cudaSafe(AT, cudaMemcpy(outDevice+(color+3)*2*sizeh,TempPacked + color*2*sizeh,2*sizeh*sizeof(float), cudaMemcpyHostToDevice ),"cudaMemcpy");
  //      cudaSafe(AT, cudaMemset(outDevice+(color+3)*2*sizeh,0,2*sizeh*sizeof(float)),"cudaMemcpy");
        
    }
   
    delete[] TempPacked;

}
void smartPackFermionOnDeviceD(double2* outDevice, Fermion* inHost){
    
    double* TempPacked = new double[2*sizeh];

    for(int color = 0; color <3 ; color++){
        for(int i =0; i < sizeh ; i++){
            //1st and only double
            TempPacked[2*i  ] = real(inHost->fermion[i].comp[color]);
            TempPacked[2*i+1] = imag(inHost->fermion[i].comp[color]);
        }
        cudaSafe(AT, cudaMemcpy(outDevice+color*2*sizeh,TempPacked,2*sizeh*sizeof(double), cudaMemcpyHostToDevice ),"cudaMemcpy");
    }

    delete[] TempPacked;

}

void smartUnpackFermionFromDevice(Fermion *outHost, float2 *inDevice){

    float* TempPacked = new float[2*sizeh];

    for(int color =0; color < 3 ; color++){
            //1st float
            cudaSafe(AT, cudaMemcpy(TempPacked, inDevice + color*2*sizeh, 2*sizeh*sizeof(float), cudaMemcpyDeviceToHost), "cudaMemcpy");
            for(int i =0; i < sizeh ; i++)
                outHost->fermion[i].comp[color] = complex<REAL> ((double)TempPacked[2*i], (double)TempPacked[2*i+1]);
    }
    
    for(int color =0; color < 3 ; color++){
            //2nd float
            cudaSafe(AT, cudaMemcpy(TempPacked, inDevice + (color+3)*2*sizeh, 2*sizeh*sizeof(float), cudaMemcpyDeviceToHost), "cudaMemcpy");
            for(int i =0; i < sizeh ; i++)
                outHost->fermion[i].comp[color]+= complex<REAL> ((double)TempPacked[2*i], (double)TempPacked[2*i+1]);
    }


    delete[] TempPacked;
}

void smartUnpackFermionFromDeviceOdd(Fermion *outHost, float2 *inDevice){

    float* TempPacked = new float[2*sizeh];

    for(int color =0; color < 3 ; color++){
            //1st float
            cudaSafe(AT, cudaMemcpy(TempPacked, inDevice + sizeh + color*2*sizeh, 2*sizeh*sizeof(float), cudaMemcpyDeviceToHost), "cudaMemcpy");
            for(int i =0; i < sizeh ; i++)
                outHost->fermion[i].comp[color] = complex<REAL> ((double)TempPacked[2*i], (double)TempPacked[2*i+1]);
    }
    
    for(int color =0; color < 3 ; color++){
            //2nd float
            cudaSafe(AT, cudaMemcpy(TempPacked, inDevice + sizeh + (color+3)*2*sizeh, 2*sizeh*sizeof(float), cudaMemcpyDeviceToHost), "cudaMemcpy");
            for(int i =0; i < sizeh ; i++)
                outHost->fermion[i].comp[color] += complex<REAL> ((double)TempPacked[2*i], (double)TempPacked[2*i+1]);

    }

    delete[] TempPacked;
}
void smartUnpackFermionFromDevice(Fermion* outHost, double2 *inDevice){

     double* TempPacked = new double[2*sizeh];

     for(int color =0; color < 3 ; color++){
            //1st and only double   
            cudaSafe(AT, cudaMemcpy(TempPacked, inDevice + color*2*sizeh, 2*sizeh*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
            for(int i = 0; i < sizeh ; i++)
                outHost->fermion[i].comp[color] = complex<REAL> (TempPacked[2*i],TempPacked[2*i+1]);
    }

}

void smartUnpackFermionFromDeviceOdd(Fermion* outHost, double2 *inDevice){

     double* TempPacked = new double[2*sizeh];

     for(int color =0; color < 3 ; color++){
            //1st and only double   
            cudaSafe(AT, cudaMemcpy(TempPacked, inDevice + (color*2+1)*sizeh, 2*sizeh*sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy");
            for(int i = 0; i < sizeh ; i++)
                outHost->fermion[i].comp[color] = complex<REAL> (TempPacked[2*i],TempPacked[2*i+1]);
    }

}


#endif 




#endif
