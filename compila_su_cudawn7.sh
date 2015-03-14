source setenv-pgi-14.9
#source setenv-pgi-15.1
source setenv-cuda.6.5

bsub -q gpuday ./compile_multishift_and_gauge.sh
