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
###################################################
####   READING COMMAND LINE OPTIONS   #############
###################################################
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
    -a|--partition) # -a is choosen for uniformity with prepare_tbps.sh
        SLURMPARTITION=$2
        shift
        ;;
    esac
    shift 
done

#########################################################
#######   LOGIC OPERATIONS ON OPTIONS   #################
#########################################################

# this function checks if an option has been set equal for both runs, e.g. GEOMFILE, or
# CONFIGOPTIONS_CSV, or MODULES_TO_LOAD_CSV. In that case sets the run-specific options 
# equal to the common value.
# If instead it has been set separately for the two runs (e.g. using GEOMFILE1 and GEOMFILE2)
# basically nothing is done, except some checking (you cannot specify GEOMFILE and GEOMFILE1 
# at the same time, for example).

OPTION12CHECK(){

    VARNAME=$1
    VARDEFAULT=$2
    CLOPTION=$3
    VARNAME1=${VARNAME}1
    VARNAME2=${VARNAME}2


    if [ ! -z ${!VARNAME1} ] || [ ! -z ${!VARNAME2} ] &&  [ ! -z ${!VARNAME} ]  
    then
        echo "$0 Error: Incompatible options: -$CLOPTION and -$CLOPTION\1 or -$CLOPTION\2"
        echo ${!VARNAME1}
        echo ${!VARNAME2}
        echo ${!VARNAME}
        exit
    fi

    if [ -z ${!VARNAME} ] && [ -z ${!VARNAME1} ] && [ -z ${!VARNAME2} ]
    then 
            eval $VARNAME=$VARDEFAULT
            echo "No geom_defines file names provided. Using default ${!VARNAME}"
        fi

        if [ -z ${!VARNAME1} ] && [ ! -z ${!VARNAME} ] 
        then 
            eval $VARNAME1=${!VARNAME}
        fi
        if [ -z ${!VARNAME2} ] && [ ! -z ${!VARNAME} ] 
        then 
            eval $VARNAME2=${!VARNAME}
        fi

}

OPTION12CHECK GEOMFILE geom_defines.txt -g
OPTION12CHECK MODULES_TO_LOAD_CSV " " -m
OPTION12CHECK CONFIGOPTIONS_CSV "gcc" -c

# SLURM-related checks
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

# Showing which options we got.
echo "Selecting  scheduler $SCHEDULER (choose with -s)"
echo "Selecting commit $OLDCOMMIT for (1) (choose with -o)"
echo "Selecting commit $NEWCOMMIT for (2) (choose with -n)" 
echo "Selecting geom_defines file $GEOMFILE1 for (1) (choose with -g1 or -g)" 
echo "Selecting geom_defines file $GEOMFILE2 for (2) (choose with -g2 or -g)" 
echo "Selecting config wrapper options $CONFIGOPTIONS_CSV1 for (1) (choose with -c1 or -c)" 
echo "Selecting config wrapper options $CONFIGOPTIONS_CSV2 for (2) (choose with -c2 or -c)" 
echo "Selecting modules \"$MODULES_TO_LOAD_CSV1\" for (1) (choose with -m1 or -m)" 
echo "Selecting modules \"$MODULES_TO_LOAD_CSV2\" for (1) (choose with -m1 or -m)" 

# creating directory with all the test scripts, which will not be changed by git
# when we checkout to another commit.
echo "Using regression test scripts at commit $CURRENTCOMMIT"
cp -r $SCRIPTSDIR test

########################################################
########  COMPILING AND PREPARING RUNS   ###############
########################################################
# function for compilation and setting up of all directories
PREPARE_RUN(){
    # variables that are expected to change form call to call are passed as arguments.
    COMMIT=$1
    TEMPGEOMFILE=$2
    TEMPCONFIGOPTIONS_CSV=$3 # compiler versions may change
    TEMPMODULES_TO_LOAD_CSV=$4
    # the other are taken as global.

    if [ ! -z $TEMPMODULES_TO_LOAD_CSV ]
    then
        echo 'Loading modules...'
        for MODULE in $(echo $TEMPMODULES_TO_LOAD_CSV | sed 's/,/ /g')
        do
            echo "module load $MODULE"
            module load $MODULE
        done
        TEMPMODULESFLAGS="-m $TEMPMODULES_TO_LOAD_CSV"
    else 
        TEMPMODULESFLAGS=""
    fi

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

    $WORKDIR/test/prepare_tbps.sh -c geom_defines.txt -p test $SLURMFLAGS $TEMPMODULESFLAGS
    if [ $? -ne 0 ]
    then 
        echo "Error: prepare_tbps.sh failed, in $PWD"
        CLEAREVERYTHING
    fi

    cd -
}

# compilation and setting up of all directories
PREPARE_RUN $OLDCOMMIT $GEOMFILE1 $CONFIGOPTIONS_CSV1 $MODULES_TO_LOAD_CSV1
PREPARE_RUN $NEWCOMMIT $GEOMFILE2 $CONFIGOPTIONS_CSV2 $MODULES_TO_LOAD_CSV2


#####################################################
####### SETTING UP ALL THE TESTS ####################
#####################################################
# setting up all the tests, submitting all the jobs, linking the files
# that must be the same

TESTSLIST=(${TESTCSV//,/ }) # using a bash array
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



if [ $DOMAINTEST == "yes" ] 
then
    # both single and double precision tests will be run on both commits.
    for DPSP in dp sp
    do 
            # 1. job on first commit
            cd $WORKDIR/$OLDCOMMIT/test.main.$DPSP
            if [  $SCHEDULER == "slurm" ]
            then 
            CAJOB=$(sbatch test.main.$DPSP\.slurm | cut -d' ' -f4)
            else
               bash ./test.main.$DPSP\.sh
            fi
            # 2. "connection" job
            cd $WORKDIR
            # WHAT EXACTLY NEEDS TO BE COPIED-LINKED? TO DEFINE
            cat > test.main.$DPSP\.connection.sh-slurm << EOF
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
            if [  $SCHEDULER == "slurm" ]
            then
            CONNJOB=$(\
                sbatch test.main.$DPSP\.connection.sh-slurm --dependency=afterok:$CAJOB\
                | cut -d' ' -f4)
            else
                bash ./test.main.$DPSP\.connection.sh-slurm
            fi
            # 3. job on second commit
            cd $WORKDIR/$NEWCOMMIT/test.main.$DPSP
            if [  $SCHEDULER == "slurm" ]
            then
            CBJOB=$(sbatch test.main.$DPSP\.slurm --dependency=afterok:$CONNJOB\
                | cut -d' ' -f4)
            else
                bash ./test.main.$DPSP\.sh
            fi
            # 4. final-check job to compare the results
            cd $WORKDIR
            # WHAT EXACTLY NEED TO BE CHECKED? TO DEFINE
            cat > test.main.$DPSP\.finalcheck.sh-slurm << EOF
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

            if [  $SCHEDULER == "slurm" ]
            then
                sbatch test.main.$DPSP\.finalcheck.sh-slurm --dependency=afterok:$CBJOB
            else
                bash test.main.$DPSP\.finalcheck.sh
            fi
    done
fi

if [ $DOMAINREVTEST == "yes" ] 
then
    # both single and double precision tests will be run on both commits.
    for $DPSP in dp sp
    do
        echo "To implement"
    done
fi

if [ $DOPGTEST == "yes" ] 
then
    # both single and double precision tests will be run on both commits.
    for $DPSP in dp sp
    do
        echo "To implement"
    done
fi

if [ $DOPGREVTEST == "yes" ] 
then
    # both single and double precision tests will be run on both commits.
    for $DPSP in dp sp
    do
        echo "To implement"
    done
fi

if [ $DODIRACTEST == "yes" ] 
then
fi

if [ $DOINVTEST == "yes" ] 
then
fi



git checkout $CURRENTCOMMIT
