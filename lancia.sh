module load pgi
module load cuda/6.0
export PGI_ACC_BUFFERSIZE=3221225472
env > ambiente

nvprof ./prog > output 2> errput 

