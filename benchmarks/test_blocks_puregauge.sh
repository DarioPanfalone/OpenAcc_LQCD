#!/bin/bash

# LINES CONTAINING '#' ARE IGNORED
BLOCKS_FILE=puregauge_blocks.txt

while read IST0 IST1 IST2 ISG3 ST0 ST1 ST2 SG3
do
   echo 
   echo Compiling pure gauge test for 
   echo IMPSTAPTILE0=$IST0  IMPSTAPTILE1=$IST1  IMPSTAPTILE2=$IST2  IMPSTAPGANG3=$ISG3
   echo STAPTILE0=$ST0  STAPTILE1=$ST1  STAPTILE2=$ST2  STAPGANG3=$SG3

   echo 
   rm rettangoli.o ipdot_gauge.o plaquettes.o su3_utilities.o
   rm sp_rettangoli.o sp_ipdot_gauge.o sp_plaquettes.o sp_su3_utilities.o
   make IMPSTAPTILE0=$IST0  IMPSTAPTILE1=$IST1  IMPSTAPTILE2=$IST2  IMPSTAPGANG3=$ISG3\
       STAPTILE0=$ST0  STAPTILE1=$ST1  STAPTILE2=$ST2  STAPGANG3=$SG3 main
   echo mv main pg_$IST0\_$IST1\_$IST2\_$ISG3\_$ST0\_$ST1\_$ST2\_$SG3 
   mv main pg_$IST0\_$IST1\_$IST2\_$ISG3\_$ST0\_$ST1\_$ST2\_$SG3 
   echo mv run/main run/pg_$IST0\_$IST1\_$IST2\_$ISG3\_$ST0\_$ST1\_$ST2\_$SG3 
   mv run/main run/pg_$IST0\_$IST1\_$IST2\_$ISG3\_$ST0\_$ST1\_$ST2\_$SG3 
   echo Compiled pg_$IST0\_$IST1\_$IST2\_$ISG3\_$ST0\_$ST1\_$ST2\_$SG3 
done < <(grep -v '#' $BLOCKS_FILE)
