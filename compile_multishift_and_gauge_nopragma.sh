
rm prog
rm *.o
pgcc -O3 -c OpenAcc/inverter_multishift_full.c  #-Mcuda=maxregcount:128
pgcc -O3 -c OpenAcc/su3_utilities.c  #-Mcuda=maxregcount:128
pgcpp -O3 -c invtest_openacc.cpp  #-Mcuda=maxregcount:128 
pgcpp *.o -o prog -O3
rm *.o

