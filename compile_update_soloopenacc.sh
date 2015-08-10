module load pgi
module load cuda

export PGI_ACC_BUFFERSIZE=3221225472

rm prog_solo
rm *.o

#pgcc  -acc=noautopar -Mlarge_arrays -Minfo=accel -O3  -v -ta=tesla:cc35 -c  Rand/random.c                   2> msg_err_0 #-Mcuda=maxregcount:128 
#pgcc  -acc=noautopar -Mlarge_arrays -Minfo=accel -O3  -v -ta=tesla:cc35 -c  OpenAcc/include_all_main.c      2> msg_err_1 #-Mcuda=maxregcount:128 
#pgcc  *.o -o prog_solo -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35                     2> msg_err_3 #-Mcuda=maxregcount:128  

pgcc  -O3  -c  Rand/random.c                   2> msg_err_0 #-Mcuda=maxregcount:128 
pgcc  -O3  -c  OpenAcc/include_all_main.c      2> msg_err_1 #-Mcuda=maxregcount:128 
pgcc  -O3  *.o -o prog_solo                    2> msg_err_3 #-Mcuda=maxregcount:128  
