module load pgi
module load cuda

export PGI_ACC_BUFFERSIZE=3221225472

rm prog
rm *.o

pgcc  -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35 -c  Rand/random.c  2> msg_err_0 #-Mcuda=maxregcount:128 
pgcc  -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35 -c  OpenAcc/include_all.c  2> msg_err_1 #-Mcuda=maxregcount:128 
pgcc  *.o -o prog -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35  2> msg_err_3 #-Mcuda=maxregcount:128 

