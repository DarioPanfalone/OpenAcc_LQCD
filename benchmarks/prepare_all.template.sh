#/usr/bin/bash

module purge
# too change
module load cuda/8.0 openmpi/2.1.1-cuda8.0 pgi/17.3

module list

grep -v '#' lattices.txt | while read L T N
do
     echo L $L "," T $T "," N $N 
     LOCT=$((T/N)) 
     DIRNAME=$L\_$T\_$N

     SIZE=$((576*L**3*(LOCT+4)))
     
     # directory creation    
     if [ ! -d $DIRNAME ]
     then 
     mkdir $DIRNAME
     fi 
     # preparing make settings
     sed 's/SED_LOCT/'$LOCT'/;s/SED_L/'$L'/;s/SED_N/'$N'/' \
            make_settings.set.proto > $DIRNAME/make_settings.set

     # compiling   
     cd $DIRNAME
     #make clean
     #touch makefile # to trigger ricompilation by higher-level makefile
     #touch make_settings.set
     #make -f ../makefile ROOTDIR=~/OPENACCREPO_TEST BENCHDIR=benchmarks_3 COMMTI_HASH=$(git rev-parse HEAD)
     cd -

     # preparing slurm scripts ans submitting them
     for PROTOSLURMSCRIPT in test.deodoe.slurm.proto test.inverter.slurm.proto test.pg.slurm.proto  #test.main.slurm.proto 
     do
         SLURMSCRIPT=${PROTOSLURMSCRIPT%.proto}
         sed 's/SEDN/'$N'/;s/SEDDIR/'$DIRNAME'/;s/SED_SIZE/'$SIZE'/' \
                         $PROTOSLURMSCRIPT > $DIRNAME/run/$SLURMSCRIPT
         # submitting
         cd $DIRNAME/run
         #sbatch ./$SLURMSCRIPT
         cd -
     done 

     for PROTOSLURMSCRIPT in test.deodoeprof.slurm.proto test.inverterprof.slurm.proto test.pgprof.slurm.proto
     do
         SLURMSCRIPT=${PROTOSLURMSCRIPT%.proto}
         sed 's/SEDN/'$N'/;s/SEDDIR/'$DIRNAME'/;s/SED_SIZE/'$SIZE'/' \
                         $PROTOSLURMSCRIPT > $DIRNAME/run/$SLURMSCRIPT
         # submitting
         cd $DIRNAME/run
         #sbatch ./$SLURMSCRIPT
         cd -
     done 

     # copying rational approximations
     cp approx* $DIRNAME/run/

done 

