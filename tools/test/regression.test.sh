#!/bin/bash
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
GEOMFILE1=geom_defines.txt
GEOMFILE2=geom_defines.txt
MODULES_TO_LOAD_CSV=""

CURRENTCOMMIT=$(git rev-parse --abbrev-ref HEAD)


CLEAREVERYTHING(){
    git checkout -f $CURRENTCOMMIT
    exit
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
    esac
    shift 
done

if [ -z $GEOMFILE ] && [ -z $GEOMFILE1  -o -z $GEOMFILE2 ] 
then
    echo "$0 Error: Incompatible options: -g and -g1 or -g2"
    exit
fi

if [ ! -z $GEOMFILE ] && [ ! -z $GEOMFILE1 ] && [ ! -z $GEOMFILE2 ]
then 
    GEOMFILE=geom_defines.txt
    echo "No geom_defines file names provided. Using default $GEOMFILE"
fi

if [ ! -z $GEOMFILE1 ] && [ -z $GEOMFILE ] 
then 
    GEOMFILE1=$GEOMFILE
fi
if [ ! -z $GEOMFILE2 ] && [ -z $GEOMFILE ] 
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

if test -z $MODULES_TO_LOAD_CSV
then
   for MODULE in $(echo $MODULES_TO_LOAD_CSV | sed 's/,/ /g')
   do
       module load $MODULE
   done
   MODULESFLAGS="-m $MODULES_TO_LOAD_CSV"
else 
   MODULESFLAGS=""
fi

if test $SCHEDULER == "slurm"
then
   SLURMFLAGS="-s"
else 
   SLURMFLAGS=""
fi

WORKDIR=$PWD
SCRIPTSDIR=$PWD/$(dirname $BASH_SOURCE)
REPODIR=$(dirname $(dirname $SCRIPTSDIR))

cp -r $SCRIPTSDIR test

for COMMIT in $OLDCOMMIT $NEWCOMMIT
do 
   cd $REPODIR
   autoreconf --install
   git checkout $COMMIT
   if test $? -ne 0
   then 
       echo "ERROR: could not checkout to commit $COMMIT"
       echo "Exiting."
       exit 1
   fi
   mkdir $WORKDIR/$COMMIT
   cp $WORKDIR/$GEOMFILE $WORKDIR/$COMMIT/
   cd $WORKDIR/$COMMIT
   yes | $REPODIR/configure_wrapper $( echo $CONFIGOPTIONS_CSV | sed 's/,/ /g')
   make -j4 &&  make install 
   if test $? -ne 0 ; then 
       echo "Error: make failed, in $PWD"
       CLEAREVERYTHING
   fi
   $WORKDIR/test/prepare_tbps.sh -c $GEOMFILE -p test $SLURMFLAGS $MODULESFLAGS
   if test $? -ne 0 ; then 
       echo "Error: prepare_tbps.sh failed, in $PWD"
       CLEAREVERYTHING
   fi
   cd -

done

git checkout $CURRENTCOMMIT
