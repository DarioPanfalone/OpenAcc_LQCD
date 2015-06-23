#source setenv-pgi-14.9
source setenv-pgi-14.7
#source setenv-pgi-15.1
#source setenv-cuda.6.5
source setenv-cuda.5.5

#bsub -q gpuday ./compile_multishift_and_gauge.sh
#bsub -q gpuday ./compile_multishift_and_gauge_CUDAVER.sh
bsub -q gpuday ./compile_update_solopartecpu.sh
