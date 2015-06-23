
rm prog
rm *.o
pgcc -O3 -c Rand/random.c  2> msg_err_1 #-Mcuda=maxregcount:128 
pgcc -O3 -c OpenAcc/include_all.c  2> msg_err_1 #-Mcuda=maxregcount:128 
pgcpp -O3 -c update_test_openacc.cpp 2> msg_err_2 #-Mcuda=maxregcount:128 
pgcpp *.o -o prog -O3   2> msg_err_3 #-Mcuda=maxregcount:128 
#rm *.o

if [ -f prog ]
then
touch successful_compilation
else
touch unsuccessful_compilation
fi
