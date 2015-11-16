module load pgi/15.9
module load cuda/7.0


export PGI_ACC_BUFFERSIZE=3221225472

rm prog_solo_GPU
rm *.o

pgcc  -acc  -Minfo=accel -O3  -v -ta=tesla:cc35,cuda7.0 -c  Rand/random.c                   2> msg_err_0 #-Mcuda=maxregcount:128 
pgcc  -acc  -Minfo=accel -O3  -v -ta=tesla:cc35,cuda7.0 -c  OpenAcc/include_all_main.c      2> msg_err_1 #-Mcuda=maxregcount:128 
pgcc  *.o -o prog_solo_GPU -acc  -Minfo=accel -O3 -v -ta=tesla:cc35,cuda7.0                     2> msg_err_3 #-Mcuda=maxregcount:128  

rm *.o

#pgcc  -acc  -Minfo=accel -O1  -v -ta=tesla:nollvm,cuda6.0,time,keep,ptxinfo -c  Rand/random.c                   2> msg_err_0 #-Mcuda=maxregcount:128 
#pgcc  -acc  -Minfo=accel -O1  -v -ta=tesla:nollvm,cuda6.0,time,keep,ptxinfo -c  OpenAcc/include_all_main.c      2> msg_err_1 #-Mcuda=maxregcount:128 
#pgcc  *.o -o prog_solo -acc  -Minfo=accel -O1 -v -ta=tesla:nollvm,cuda6.0,time,keep,ptxinfo                     2> msg_err_3 #-Mcuda=maxregcount:128  

#-ta=tesla:cc35

#pgcc  -O3  -c  Rand/random.c                   2> msg_err_0 #-Mcuda=maxregcount:128 
#pgcc  -O3  -c  OpenAcc/include_all_main.c      2> msg_err_1 #-Mcuda=maxregcount:128 
#pgcc  -O3  *.o -o prog_solo                    2> msg_err_3 #-Mcuda=maxregcount:128  
