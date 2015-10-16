module load pgi/15.9
module load cuda/7.0
export PGI_ACC_BUFFERSIZE=3221225472
env > ambiente

nvprof ./prog_solo > output 2> errput 
#./prog > output 2> errput 

