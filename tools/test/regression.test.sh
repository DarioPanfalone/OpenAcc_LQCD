#!/usr/bin/env bash
# This script runs the tests comparing the results between different commits.
# Options:
#   -s, --scheduler: The job scheduler that will be used to launch the tests.
#                    If 'bash' is chosen, everything will be run in serial on 
#                    the machine this scripts runs on.
#   -o, --oldcommit: The old commit to test against.
#   -n, --newcommit: The new commit to test.
#   -g, --geomfile : The file containing the geometry information for the 
#                    compilation.
#   -m, --modules  : Comma-separated-value list of modules to load (will be 
#                    propagated to the rest of the code.
#   -c, --config-options: Comma-separated-value list of options to pass to the
#                    configure-wrapper script for compilation.

# Default values
# Scheduler can be either 'bash' or a proper scheduling system. 
SCHEDULER=bash
OLDCOMMIT=8d4dcd6cbf91be1354a471665ec71a1db4756628
NEWCOMMIT=$(git rev-parse HEAD)
GEOMFILE=
GEOMFILE1=
GEOMFILE2=
MODULES_TO_LOAD_CSV=""
WORKDIR=$PWD
SCRIPTSDIR=$PWD/$(dirname $BASH_SOURCE)
REPODIR=$(dirname $(dirname $SCRIPTSDIR))
BUILDJOBS=4
SLURMPARTITION=
CONFIGOPTIONS_CSV=""
CONFIGOPTIONS_CSV1=""
CONFIGOPTIONS_CSV2=""


# Types of test
DOMAINTEST=no
DOMAINREVTEST=no
DOPGTEST=no
DOPGREVTEST=no
DODIRACTEST=no
DOINVTEST=no

CURRENTCOMMIT=$(git rev-parse --abbrev-ref HEAD)


CLEAREVERYTHING(){
    git checkout -f $CURRENTCOMMIT
    exit
}

trap CLEAREVERYTHING EXIT


# function for compilation and setting up of all directories
PREPARE_TEST(){
    # variables that are expected to change form call to call are passed as arguments.
    COMMIT=$1
    TEMPGEOMFILE=$2
    TEMPCONFIGOPTIONS_CSV=$3 # compiler versions may change
    # the other are taken as global.

    cd $REPODIR
    autoreconf --install
    git checkout $COMMIT
    if [ $? -ne 0 ]
    then 
        echo "ERROR: could not checkout to commit $COMMIT"
        echo "Exiting."
        exit 1
    fi
    mkdir $WORKDIR/$COMMIT
    # configure_wrapper and configure need the geom_defines.txt file to be named
    # exactly geom_defines.txt
    cp $WORKDIR/$TEMPGEOMFILE $WORKDIR/$COMMIT/geom_defines.txt
    cd $WORKDIR/$COMMIT
    yes | $REPODIR/configure_wrapper $( echo $TEMPCONFIGOPTIONS_CSV | sed 's/,/ /g')
    make -j$BUILDJOBS &&  make install 
    if [ $? -ne 0 ]
    then 
        echo "Error: make failed, in $PWD"
        CLEAREVERYTHING
    fi

    $WORKDIR/test/prepare_tbps.sh -c geom_defines.txt -p test $SLURMFLAGS $MODULESFLAGS
    if [ $? -ne 0 ]
    then 
        echo "Error: prepare_tbps.sh failed, in $PWD"
        CLEAREVERYTHING
    fi

    cd -
}



while (( "$#" ))
do
    case $1 in
    -s|--scheduler)
        SCHEDULER=$2
        shift
        ;;

    -o|--oldcommit) 
        OLDCOMMIT=$2
        shift
        ;;
    -n|--newcommit)
        NEWCOMMIT=$2
        shift
        ;;
    -g|--geomfile)
        GEOMFILE=$2
        shift
        ;;
    -g1|--geomfile1)
        GEOMFILE1=$2
        shift
        ;;
    -g2|--geomfile2)
        GEOMFILE2=$2
        shift
        ;;
    -m|--modules)
        MODULES_TO_LOAD_CSV=$2
        shift
        ;;
    -c|--config-options)
        CONFIGOPTIONS_CSV=$2
        shift
        ;;
    -c1|--config-options1)
        CONFIGOPTIONS_CSV1=$2
        shift
        ;;
    -c2|--config-options2)
        CONFIGOPTIONS_CSV2=$2
        shift
        ;;
    -bj|--build-jobs)
        BUILDJOBS=$2
        shift
        ;;
    -t|--test)
        TESTSCSV=$2
        shift
        ;;
    -p|--partition)
        SLURMPARTITION=$2
        shift
        ;;
    esac
    shift 
done

if [ ! -z $GEOMFILE1 ] || [ ! -z $GEOMFILE2 ] &&  [ ! -z $GEOMFILE ]  
then
    echo "$0 Error: Incompatible options: -g and -g1 or -g2"
    echo $GEOMFILE1
    echo $GEOMFILE2
    echo $GEOMFILE
    exit
fi

if [ -z $GEOMFILE ] && [ -z $GEOMFILE1 ] && [ -z $GEOMFILE2 ]
then 
    GEOMFILE=geom_defines.txt
    echo "No geom_defines file names provided. Using default $GEOMFILE"
fi

if [ -z $GEOMFILE1 ] && [ ! -z $GEOMFILE ] 
then 
    GEOMFILE1=$GEOMFILE
fi
if [ -z $GEOMFILE2 ] && [ ! -z $GEOMFILE ] 
then 
    GEOMFILE2=$GEOMFILE
fi

echo "Selecting  scheduler $SCHEDULER (choose with -s)"
echo "Selecting old commit $OLDCOMMIT (choose with -o)"
echo "Selecting new commit $NEWCOMMIT (choose with -n)" 
echo "Selecting geom_defines,1 file $GEOMFILE1 (choose with -g1)" 
echo "Selecting geom_defines,2 file $GEOMFILE2 (choose with -g2)" 
echo "Selecting modules $MODULES_TO_LOAD_CSV (choose with -m)" 
echo "Selecting config wrapper options $CONFIGOPTIONS_CSV (choose with -c)" 



if [ ! -z $MODULES_TO_LOAD_CSV ]
then
    echo 'Loading modules...'
    for MODULE in $(echo $MODULES_TO_LOAD_CSV | sed 's/,/ /g')
    do
        echo "module load $MODULE"
        module load $MODULE
    done
    MODULESFLAGS="-m $MODULES_TO_LOAD_CSV"
else 
    MODULESFLAGS=""
fi

if [ $SCHEDULER == "slurm" ]
then
    if [ -z $SLURMPARTITION ] 
    then
        echo "$0 Error: Slurm Partition must be selected with -a or --partition."
        exit
    fi
    SLURMFLAGS="-s -a $SLURMPARTITION"
else 
    SLURMFLAGS=""
fi




TESTSLIST=(${TESTCSV//,/ }) # bash arrays
for TEST in ${TESTSLIST[@]}
do
    case $TEST in
        main)
            DOMAINTEST=yes
            ;;
        mainrev)
            DOMAINREVTEST=yes
            ;;
        pg)
            DOPGTEST=yes
            ;;
        pgrev)
            DOPGREVTEST=yes
            ;;
        dirac)
            DODIRACTEST=yes
            ;;
        inv|inverter)
            DOINVTEST=yes
            ;;
    esac
done

cp -r $SCRIPTSDIR test

# compilation and setting up of all directories
PREPARE_TEST $OLDCOMMIT $GEOMFILE1 $CONFIGOPTIONS_CSV1
PREPARE_TEST $NEWCOMMIT $GEOMFILE2 $CONFIGOPTIONS_CSV2

# setting up all the tests, submitting all the jobs, linking the files
# that must be the same

if [ $DOMAINTEST == "yes" ] 
then
    # both single and double precision tests will be run on both commits.
    for DPSP in dp sp
    do 
        if [  $SCHEDULER == "slurm" ]
        then 
            diff $WORKDIR/$OLDCOMMIT/test.main.$DPSP\.set $WORKDIR/$NEWCOMMIT/test.main.$DPSP\.set
            # 1. job on first commit
            cd $WORKDIR/$OLDCOMMIT/test.main.$DPSP
            CAJOB=$(sbatch test.main.$DPSP\.slurm | cut -d' ' -f4)
            # 2. "connection" job
            cd $WORKDIR
            # WHAT EXACTLY NEEDS TO BE COPIED-LINKED? TO DEFINE
            cat > test.main.$DPSP\.connection.slurm << EOF
#!/bin/bash
#SBATCH --job-name=main.$DPSP\.connection
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --error=main.$DPSP\.connection.%J.err 
#SBATCH --output=main.$DPSP\.connection.%J.out
#SBATCH --partition=$SLURMPARTITION

for file in $WORKDIR/$OLDCOMMIT/test.main.$DPSP/*norndtest* 
do 
    ln $file $WORKDIR/$OLDCOMMIT/test.main.$DPSP
done

EOF
            CONNJOB=$(\
                sbatch test.main.$DPSP\.connection.slurm --dependency=afterok:$CAJOB\
                | cut -d' ' -f4)
            # 3. job on second commit
            cd $WORKDIR/$NEWCOMMIT/test.main.$DPSP
            CBJOB=$(sbatch test.main.$DPSP\.slurm --dependency=afterok:$CONNJOB\
                | cut -d' ' -f4)

            # 4. final-check job to compare the results
            cd $WORKDIR
            # WHAT EXACTLY NEED TO BE CHECKED? TO DEFINE
            cat > test.main.$DPSP\.finalcheck.slurm << EOF
#!/bin/bash
#SBATCH --job-name=main.$DPSP\.finalcheck
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --error=main.$DPSP\.finalcheck.%J.err 
#SBATCH --output=main.$DPSP\.finalcheck.%J.out
#SBATCH --partition=$SLURMPARTITION

for file in $WORKDIR/$OLDCOMMIT/test.main.$DPSP/* 
do 
    ./test/diff_ascii_files.py $file $WORKDIR/$OLDCOMMIT/test.main.$DPSP/$(basename $file) 
done
EOF           
            sbatch test.main.$DPSP\.finalcheck.slurm --dependency=afterok:$CBJOB

        else
           # BASH VERSION TO WRITE 
        fi

    done
fi

if [ $DOMAINREVTEST == "yes" ] 
then
    # both single and double precision tests will be run on both commits.
    for $DPSP in dp sp
    do
    done
fi

if [ $DOPGTEST == "yes" ] 
then
    # both single and double precision tests will be run on both commits.
fi

if [ $DOPGREVTEST == "yes" ] 
then
    # both single and double precision tests will be run on both commits.
fi

if [ $DODIRACTEST == "yes" ] 
then
fi

if [ $DOINVTEST == "yes" ] 
then
fi



git checkout $CURRENTCOMMIT
