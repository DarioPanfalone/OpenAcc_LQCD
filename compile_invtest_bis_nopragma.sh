
rm prog
rm *.o
pgcc  -O3   -c OpenAcc/inverter_full.c
pgcpp -O3   -c invtest_openacc.cpp
pgcpp -O3  *.o -o prog 
rm *.o

