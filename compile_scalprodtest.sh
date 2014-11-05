
rm prog
pgcc -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda5.5,time,keep,ptxinfo -c OpenAcc/fermion_matrix.c  #-Mcuda=maxregcount:128
pgcpp -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda5.5,time,keep,ptxinfo -c scal_prod_check.cpp  #-Mcuda=maxregcount:128 
pgcpp *.o -o prog -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda5.5,time,keep,ptxinfo   #-Mcuda=maxregcount:128 
rm *.o


