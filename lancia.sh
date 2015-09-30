module load pgi
module load cuda
export PGI_ACC_BUFFERSIZE=3221225472
env > ambiente

nvprof ./prog_solo > output 2> errput 
#./prog > output 2> errput 

