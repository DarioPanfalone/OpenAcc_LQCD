
rm prog_PGCC
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
OpenAcc/su3_utilities.c            \
OpenAcc/update_versatile.c         \
RationalApprox/rationalapprox.c         \
Include/fermion_parameters.c \
"


pgcc -O0  -c  Rand/random.c                 2> pgccmsg_err_0  
echo Compiling OpenAcc/include_all_main.c ...
pgcc -O0  -c  OpenAcc/include_all_main.c  2> pgccmsg_err_1
for sourcefile in $SOURCES
do 
echo Compiling $sourcefile ...
echo Compiling $sourcefile ... >> pgccmsg_err_1
pgcc -O0  -c  $sourcefile           2>> pgccmsg_err_1
done 
pgcc -O0  *.o -o prog_PGCC  -lm              2> pgccmsg_err_3






















