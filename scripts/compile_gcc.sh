
rm prog_GCC
rm *.o

SOURCES="\
OpenAcc/alloc_vars.c                \
OpenAcc/backfield.c                \
OpenAcc/cayley_hamilton.c          \
OpenAcc/cooling.c                  \
OpenAcc/deviceinit.c               \
OpenAcc/fermion_force.c            \
OpenAcc/fermion_force_utilities.c  \
OpenAcc/fermion_matrix.c           \
OpenAcc/fermionic_utilities.c      \
OpenAcc/find_min_max.c             \
OpenAcc/geometry.c                 \
OpenAcc/inverter_multishift_full.c \
OpenAcc/inverter_full.c \
OpenAcc/ipdot_gauge.c              \
OpenAcc/md_integrator.c            \
OpenAcc/random_assignement.c       \
OpenAcc/rettangoli.c               \
OpenAcc/single_types.c             \
OpenAcc/stouting.c                 \
Meas/gauge_meas.c                 \
Meas/ferm_meas.c                 \
OpenAcc/su3_utilities.c            \
OpenAcc/update_versatile.c         \
RationalApprox/rationalapprox.c         \
Include/fermion_parameters.c \
"


gcc -O3  -c  Rand/random.c                 2> gccmsg_err_0  
echo Compiling OpenAcc/main.c ...
gcc -std=c99 -O3  -c  OpenAcc/main.c  2> gccmsg_err_1
for sourcefile in $SOURCES
do 
echo Compiling $sourcefile ...
gcc -std=c99 -O3  -c  $sourcefile           2>> gccmsg_err_1
done 
gcc -O3  *.o -o prog_GCC  -lm              2> gccmsg_err_3






















