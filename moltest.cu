#define USEGPU //necessary for all the inclusions to be performed correctly
//#define DOUBLE_PREC

#define NITERMAX 1000
//#include <sys/time.h>
#include <ctime>

#include <iostream>

//constants, global variables, geometric parameters et cetera 
#include "./Include/global_const.cc"
#include "./Include/global_macro.cc"
#include "./Include/global_var.cc"
#include "./Init/init.cc"
#include "./Fermions/fermions.cc"




#include "./DeviceInfo/device_info.cu"
#include "./Cuda/cuda_globalstuff.cu"
#include "./Packer/packer.cc"
#ifndef GPU_DOUBLE_PREC
 #include "./Cuda/cuda_dslash_eo.cu"
#else 
 #include "./Cuda/cuda_dslash_dd_eo.cu"
#endif

using namespace std;

int main(){

/**************************************
Matrix vector multiplication
***************************************/

  //chooseGPU(0);


//All SU(3) links set to identity
//A Random SU(3) field can be produced with init(1)
   init(1);
   cout << "Initialized Random Geuge Matrix.\n";
//If needed the gauge field can be loaded from file
//it can also be loaded using Init(2), in that case a 
//file named 'config' will be looked for.
//The gauge field can also be saved :
//gauge_conf->saveToFile("TestConf.cnf");
   Fermion* tempFermion = new Fermion();//initialized to 0
   //tempFermion->gauss(0); //produces a gaussian noise on the first component 
                            //the other components are set to constant
   cout << "Loading fermion from StartFermion.fer...\n";
   tempFermion->loadFromFile("StartFermion.fer");//choose a name

     
// GPU gauge and geometry setup
    smartpack_gauge( gauge_field_packed , gauge_conf );
    make_shift_table(shift_table);//writes the geometry info in shift_table
    chooseGPU(gpu_device_to_use);
    cuda_init_all();//moves the packed conf and shift_table onto the gpu
// END GPU gauge and geometry setup

// creating GPU arrays for the fermions
// actually for the current implementation of DslashOperatorEO()
// one would suffice

#ifndef DOUBLE_PREC
    size_t vector_size = sizeof(float2)*6*size ;//in both cases,same dimension
    float2 *Mp, *MMp;
#else
    //in this case the double precision fermion vector must be reconstructed.
    size_t vector_size = sizeof(double2)*3*size ;//in both cases,same dimension
    double2 *Mp, *MMp;
    float2 *startFermGPU;//fermion vector will initially be loaded here
#endif
   
 
    cudaSafe(AT,cudaMalloc((void**)&Mp ,vector_size), "cudaMalloc");
    cudaSafe(AT,cudaMalloc((void**)&MMp,vector_size), "cudaMalloc");
    cudaSafe(AT,cudaMemset(Mp,0,vector_size),"cudaMemset");
    cudaSafe(AT,cudaMemset(MMp,0,vector_size),"cudaMemset");
#ifdef DOUBLE_PREC    
    cudaSafe(AT,cudaMalloc((void**)&startFermGPU,vector_size), "cudaMalloc");
    cudaSafe(AT,cudaMemset(startFermGPU,0,vector_size),"cudaMemset");
#endif
// this function takes a fermion on the host 
// and puts it into an array directly on the gpu


// DslashOperatorEO called with 1 as the third argument reads the EVEN part 
// of the 2nd argument and writes the result in the ODD part of the 1st 
// DslashOperatorEO called with -1 as the third argument reads the ODD part 
// of the 2nd argument and writes the result in the EVEN part of the 1st 



#ifndef DOUBLE_PREC
   smartPackFermionOnDeviceD(MMp,tempFermion);//packing the fermion vector into
                                              // MMp
#else
   smartPackFermionOnDeviceD(startFermGPU,tempFermion);//loading fermion into
                                                       // startFermGPU
   InitR(MMp, startFermGPU ,3*size);// Initializing the double precision 
                                    // fermion vector, into MMp
#endif


//  cudaEvent_t start, stop;
//  float time;
//  cudaEventCreate(&start);
//  cudaEventCreate(&stop);
     
//  cudaEventRecord(start, 0);
    
//  struct timeval tim;  
//  gettimeofday(&tim, NULL);  
//  double t1=tim.tv_sec+(tim.tv_usec/1000000.0);  

    clock_t start, end;
    start = clock();

    for(int iter = 0; iter < NITERMAX ; iter++){
#ifndef DOUBLE_PREC
//      cout << "Multiplying by Doe, on device.\n";
        DslashOperatorEO(Mp,MMp, 1 );
//      cout << "Multiplying by Deo, on device.\n";
        DslashOperatorEO(MMp, Mp, -1 ); 
#else
//      cout << "Multiplying by DDDoe, on device.\n";
        DslashOperatorDDEO(Mp, MMp, 1 );
//      cout << "Multiplying by DDDeo, on device.\n";
        DslashOperatorDDEO(MMp, Mp, -1 );
#endif
   }

   end = clock();

//  gettimeofday(&tim, NULL);  
//  double t2=tim.tv_sec+(tim.tv_usec/1000000.0); 

    
//  cudaEventRecord(stop, 0);
//  cudaEventSynchronize(stop);
//  cudaEventElapsedTime(&time, start, stop);

//  cout << "TEST ran in "<< time / NITERMAX << "ms. \n";
//  cout << "Test run (in "<< (t2-t1)/ NITERMAX<< " sec.)\n";
    cout << "Test run (in "<< (double)(end-start)/( NITERMAX *CLOCKS_PER_SEC) << " sec.)\n";





// this function takes an array on the gpu 
// and puts it into a fermion on the host
   smartUnpackFermionFromDevice(tempFermion,MMp);

#ifndef DOUBLE_PREC
    const char* endFermionFilename = "EndFermionGPU_SinglePrec.fer"; 
#else
    const char* endFermionFilename = "EndFermionGPU_DoublePrec.fer";
#endif
   
    cout << "Saving fermion in " << endFermionFilename << endl; 
    tempFermion->saveToFile(endFermionFilename);
 
    cudaFree(MMp);
    cudaFree(Mp);
#ifdef DOUBLE_PREC
    cudaFree(startFermGPU);
#endif
    cuda_end();
 
    delete tempFermion;
    return 0;
 
}
