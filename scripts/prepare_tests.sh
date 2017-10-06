COMPILE_INFO_FILENAME=$1
SCRIPTSDIR=$2
if test ! -f $COMPILE_INFO_FILENAME
then
    echo File $COMPILE_INFO_FILENAME does not exist! 
    exit
else
    echo Reading compilation/geometry parameters from $COMPILE_INFO_FILENAME ...
fi

if test $# -eq 1
then 
    SCRIPTSDIR=$(dirname $BASH_SOURCE)
    echo SCRIPTSDIR not set, setting to $SCRIPTSDIR
fi

FILETEMPLATEDIR=$SCRIPTSDIR


if test ! -f "$FILETEMPLATEDIR/test.pg.DPSP.set.proto" -a \
    -f "$FILETEMPLATEDIR/fermion_parameters.set"
then
    echo Directory $SCRIPTSDIR does not exist or not contains benchmarks templates.
    exit
else
    echo Copying template benchmark input files from directory $SCRIPTSDIR
fi



L0=$( grep -E "^\s*LOC_N0\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
L1=$( grep -E "^\s*LOC_N1\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
L2=$( grep -E "^\s*LOC_N2\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')
L3=$( grep -E "^\s*LOC_N3\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')

TASKS=$( grep -E "^\s*NRANKS_D3\s+" $COMPILE_INFO_FILENAME | awk '{print $2}')


# size of the local gauge configuration in bytes.
# the '+4' is needed only on the direction which is cut (L3)
SIZE=$((576*L0*L1*L2*(L3+4)))

LABEL=$L0$L1$L2$L3\_$TASKS


echo Reading modules to load in slurm scripts...
MODULES_TO_LOAD=$(module list -t 2>&1 | tail -n+2)
# creating a string with all the modules in one line
for MODULE_TO_LOAD in $MODULES_TO_LOAD
do
    MODULES_TO_LOAD_TEMP="$MODULES_TO_LOAD_TEMP $MODULE_TO_LOAD"
done
MODULES_TO_LOAD=$MODULES_TO_LOAD_TEMP

if test "$MODULES_TO_LOAD" == "" 
then
    echo "WARNING: No module loaded."
    echo "Preparing only shell scripts and input files, no slurm scripts"
    PREPARESLURM=no
else 
    PREPARESLURM=yes
fi

echo Modules to load:$MODULES_TO_LOAD....

# preparing base benchmark file
FILEBASE=test.pg.DPSP.set.proto
for MODE in dp sp
do
    if test $MODE == sp
    then 
        SINGLEPRECISIONMD=1
    elif test $MODE == dp
    then 
        SINGLEPRECISIONMD=0
    fi
    OUTPROTOSETFILE=$(echo ${FILEBASE%.proto} | sed 's/DPSP/'$MODE'/')
    echo $FILEBASE | sed 's/DPSP/'$MODE'/'
    echo Writing $OUTPROTOSETFILE
    sed 's/SEDNRANKS/'$TASKS'/' $FILETEMPLATEDIR/$FILEBASE |\
        sed 's/SEDNX/'$L0'/' | sed 's/SEDNY/'$L1'/' | sed 's/SEDNZ/'$L2'/' |\
        sed 's/SEDNODEDIM/16/' |  sed 's/SEDNT/'$((L3*TASKS))'/' |\
        sed 's/DPSP/'$MODE'/' |  sed 's/SINGLEPRECISIONMD/'$SINGLEPRECISIONMD'/'  >\
        $OUTPROTOSETFILE


    echo Writing test.$MODE.set
    cat $FILETEMPLATEDIR/fermion_parameters.set $OUTPROTOSETFILE > test.main.$MODE.set
done

echo Writing test.deo_doe_test.set
cp test.main.dp.set test.deo_doe_test.set
echo Writing test.inverter_multishift_test.set
cp test.main.dp.set test.inverter_multishift_test.set



EXECUTABLES="main.dp main.sp deo_doe_test inverter_multishift_test pg.dp pg.sp"

# executable 'pg' is actually 'main'
rm ./bin/pg
ln -s ./main ./bin/pg

for EXECUTABLE in $EXECUTABLES
do
    INPUTFILENAME="test.$EXECUTABLE.set"
    PROFILINGSTR=""
    OUTPUTFILENAME="out.$EXECUTABLE.txt"
    SLURMFILENAME="test.$EXECUTABLE.slurm"
    BASHFILENAME="test.$EXECUTABLE.sh"

    if test $PREPARESLURM == yes
    then
        cat > $SLURMFILENAME << EOF
#!/bin/bash
#SBATCH --job-name=test.${EXECUTABLE}_$LABEL
#SBATCH --ntasks=$TASKS
#SBATCH --cpus-per-task=1
#SBATCH --error=test.${EXECUTABLE}.%J.err 
#SBATCH --output=test.${EXECUTABLE}.%J.out
#SBATCH --gres=gpu:16
#SBATCH --partition=shortrun
#SBATCH --mem-per-cpu=12000

module purge
module load $MODULES_TO_LOAD
export PGI_ACC_BUFFERSIZE=$SIZE

rm stop

srun --cpu_bind=v,sockets ./bin/$EXECUTABLE ./$INPUTFILENAME > $OUTPUTFILENAME

EOF
    fi

    cat > $BASHFILENAME << EOF
#!/bin/bash

mpirun -n $TASKS ./bin/$EXECUTABLE ./$INPUTFILENAME > $OUTPUTFILENAME

EOF




done
