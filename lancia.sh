module load pgi
module load cuda/6.5
export PGI_ACC_BUFFERSIZE=3221225472
env > ambiente

./prog > output 2> errput 

