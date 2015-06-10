#module load pgi
#module load cuda/6.5

#module load pgi/15.1
#module load cuda/6.0

#module load pgi/14.7
#module load cuda/6.0

module load pgi
module load cuda

export PGI_ACC_BUFFERSIZE=3221225472

rm prog
rm *.o

#pgcc -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda6.5,time,keep,ptxinfo -c  Rand/random.c  2> msg_err_0 #-Mcuda=maxregcount:128 
#pgcc -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda6.5,time,keep,ptxinfo -c  OpenAcc/include_all.c  2> msg_err_1 #-Mcuda=maxregcount:128 
#pgcpp -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda6.5,time,keep,ptxinfo -c update_test_openacc_singlestep.cpp 2> msg_err_2 #-Mcuda=maxregcount:128 
#pgcpp *.o -o prog -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda6.5,time,keep,ptxinfo  2> msg_err_3 #-Mcuda=maxregcount:128 

#pgcc -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda6.0,time,keep,ptxinfo -c  Rand/random.c  2> msg_err_0 #-Mcuda=maxregcount:128 
#pgcc -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda6.0,time,keep,ptxinfo -c  OpenAcc/include_all.c  2> msg_err_1 #-Mcuda=maxregcount:128 
#pgcpp -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda6.0,time,keep,ptxinfo -c update_test_openacc_singlestep.cpp 2> msg_err_2 #-Mcuda=maxregcount:128 
#pgcpp *.o -o prog -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda6.0,time,keep,ptxinfo  2> msg_err_3 #-Mcuda=maxregcount:128 


#pgcc -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda5.5,time,keep,ptxinfo -c  Rand/random.c  2> msg_err_0 #-Mcuda=maxregcount:128 
#pgcc -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda5.5,time,keep,ptxinfo -c  OpenAcc/include_all.c  2> msg_err_1 #-Mcuda=maxregcount:128 
#pgcpp -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda5.5,time,keep,ptxinfo -c update_test_openacc_singlestep.cpp 2> msg_err_2 #-Mcuda=maxregcount:128 
#pgcpp *.o -o prog -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda5.5,time,keep,ptxinfo  2> msg_err_3 #-Mcuda=maxregcount:128 


#pgcc -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda7.0,time,keep,ptxinfo -c  Rand/random.c  2> msg_err_0 #-Mcuda=maxregcount:128 
#pgcc -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda7.0,time,keep,ptxinfo -c  OpenAcc/include_all.c  2> msg_err_1 #-Mcuda=maxregcount:128 
#pgcpp -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda7.0,time,keep,ptxinfo -c update_test_openacc_singlestep.cpp 2> msg_err_2 #-Mcuda=maxregcount:128
#pgcpp *.o -o prog -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda7.0,time,keep,ptxinfo  2> msg_err_3 #-Mcuda=maxregcount:128 

pgcc -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35 -c  Rand/random.c  2> msg_err_0 #-Mcuda=maxregcount:128 
pgcc -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35 -c  OpenAcc/include_all.c  2> msg_err_1 #-Mcuda=maxregcount:128 
pgcpp -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35 -c update_test_openacc_singlestep.cpp 2> msg_err_2 #-Mcuda=maxregcount:128
pgcpp *.o -o prog -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35  2> msg_err_3 #-Mcuda=maxregcount:128 



#rm *.o

#if [ -f prog ]
#then
#touch successful_compilation
#else
#touch unsuccessful_compilation
#fi
