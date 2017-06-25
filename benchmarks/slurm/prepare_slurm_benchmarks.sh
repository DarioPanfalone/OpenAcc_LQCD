

COMPILE_INFO_FILENAME=$1
BENCHDIR=$2
if test ! -f $COMPILE_INFO_FILENAME
then
    echo File $COMPILE_INFO_FILENAME does not exist! 
    exit
else
    echo Reading compilation/geometry parameters from $COMPILE_INFO_FILENAME ...
fi

if test $# -eq 1
then 
    BENCHDIR=$(dirname $BASH_SOURCE)
    BENCHDIR=$(dirname $BENCHDIR)
    echo BENCHDIR not set, setting to $BENCHDIR
fi

FILETEMPLATEDIR=$BENCHDIR/setting_file_templates


if test ! -f "$FILETEMPLATEDIR/benchmark.pg.set.proto" -a \
    -f "$FILETEMPLATEDIR/benchmark_for_profiling.pg.set.proto"  -a \
    -f "$FILETEMPLATEDIR/fermion_parameters.set"
then
    echo Directory $BENCHDIR does not exist or not contains benchmarks templates.
    exit
else
    echo Copying template benchmark input files from directory $BENCHDIR
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
echo Modules to load: $MODULES_TO_LOAD

# preparing base benchmark file
for FILEBASE in benchmark.pg.set.proto benchmark_profiling.pg.set.proto
do
echo $FILEBASE
sed 's/SEDNRANKS/'$TASKS'/' $FILETEMPLATEDIR/$FILEBASE |\
    sed 's/SEDNX/'$L0'/' | sed 's/SEDNY/'$L1'/' | sed 's/SEDNZ/'$L2'/' |\
    sed 's/SEDNODEDIM/16/' |  sed 's/SEDNT/'$((L3*TASKS))'/' > ${FILEBASE%.proto}
done

cat $FILETEMPLATEDIR/fermion_parameters.set benchmark.pg.set >> benchmark.set
cat $FILETEMPLATEDIR/fermion_parameters.set benchmark_profiling.pg.set >>\
    benchmark_profiling.set

cp benchmark.set benchmark.deodoe.set
cp benchmark.set benchmark.inverter.set
cp benchmark_profiling.set benchmark_profiling.deodoe.set
cp benchmark_profiling.set benchmark_profiling.inverter.set


EXECUTABLES="deodoe inverter pg"
for EXECUTABLE in $EXECUTABLES
do
    for PROFILING in YES NO
    do
        
        if test $PROFILING == "YES"
        then
            INPUTFILENAME="benchmark_profiling.$EXECUTABLE.set"
            PROFILINGSTR="nvprof --log-file \"$EXECUTABLE.%p"\"
            OUTPUTFILENAME="out_profiling.$EXECUTABLE.txt"
            SLURMFILENAME="test_profiling.$EXECUTABLE.slurm"
        else
            INPUTFILENAME="benchmark.$EXECUTABLE.set"
            PROFILINGSTR=""
            OUTPUTFILENAME="out.$EXECUTABLE.txt"
            SLURMFILENAME="test.$EXECUTABLE.slurm"
        fi


        cat > $SLURMFILENAME << EOF
#!/bin/bash
#SBATCH --job-name=$EXECUTABLE_$LABEL
#SBATCH --ntasks=$TASKS
#SBATCH --cpus-per-task=1
#SBATCH --error=$EXECUTABLE_.%J.err 
#SBATCH --output=$EXECUTABLE_.%J.out
#SBATCH --gres=gpu:16
#SBATCH --partition=shortrun
#SBATCH --mem-per-cpu=12000

module purge
module load $MODULES_TO_LOAD
export PGI_ACC_BUFFERSIZE=$SIZE

rm stop

srun --cpu_bind=v,sockets $PROFILINGSTR $EXECUTABLE ./$INPUTFILENAME > $OUTPUTFILENAME
 
EOF
    done

done
