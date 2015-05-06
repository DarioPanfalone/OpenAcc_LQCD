
rm prog
rm *.o
pgcc -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda5.5,time,keep,ptxinfo -c OpenAcc/include_all.c  2> msg_err_1 #-Mcuda=maxregcount:128 
pgcpp -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda5.5,time,keep,ptxinfo -c update_test_openacc.cpp 2> msg_err_2 #-Mcuda=maxregcount:128 
pgcpp *.o -o prog -acc=noautopar -Mlarge_arrays -Minfo=accel -O3 -v -ta=tesla:cc35,cuda5.5,time,keep,ptxinfo  2> msg_err_3 #-Mcuda=maxregcount:128 
#rm *.o

if [ -f prog ]
then
touch successful_compilation
else
touch unsuccessful_compilation
fi
