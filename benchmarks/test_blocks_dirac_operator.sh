#!/bin/bash

# LINES CONTAINING '#' ARE IGNORED
BLOCKS_FILE=dirac_blocks.txt

while read DDT0 DDT1 DDT2 DDG3
do
   echo 
   echo Compiling dirac test for 
   echo DEODOETILE0=$DDT0 DEODOETILE1=$DDT1  DEODOETILE2=$DDT2 DEODOEGANG3=$DDG3
   echo 
   rm fermion_matrix.o sp_fermion_matrix.o
   make DEODOETILE0=$DDT0 DEODOETILE1=$DDT1  DEODOETILE2=$DDT2 DEODOEGANG3=$DDG3 deo_doe_test
   echo mv deo_doe_test deo_doe_test_$DDT0\_$DDT1\_$DDT2\_$DDG3
   mv deo_doe_test deo_doe_test_$DDT0\_$DDT1\_$DDT2\_$DDG3
   echo Compiled deo_doe_test_$DDT0\_$DDT1\_$DDT2\_$DDG3

done < <(grep -v '#' $BLOCKS_FILE)
