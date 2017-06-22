AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = main 


# single precision code generation
#
#__sp_sources = OpenAcc/sp_alloc_vars.c  \
#	OpenAcc/sp_alloc_vars.h  \
#	OpenAcc/sp_backfield.c  \
#	OpenAcc/sp_backfield.h  \
#	OpenAcc/sp_cayley_hamilton.h  \
#	OpenAcc/sp_fermion_force.c  \
#	OpenAcc/sp_fermion_force.h  \
#	OpenAcc/sp_fermion_force_utilities.c  \
#	OpenAcc/sp_fermion_force_utilities.h  \
#	OpenAcc/sp_fermionic_utilities.c  \
#	OpenAcc/sp_fermionic_utilities.h  \
#	OpenAcc/sp_fermion_matrix.c  \
#	OpenAcc/sp_fermion_matrix.h  \
#	OpenAcc/sp_inverter_full.c  \
#	OpenAcc/sp_inverter_full.h  \
#	OpenAcc/sp_inverter_multishift_full.c  \
#	OpenAcc/sp_inverter_multishift_full.h  \
#	OpenAcc/sp_ipdot_gauge.c  \
#	OpenAcc/sp_ipdot_gauge.h  \
#	OpenAcc/sp_matvecmul.h  \
#	OpenAcc/sp_md_integrator.c  \
#	OpenAcc/sp_md_integrator.h  \
#	OpenAcc/sp_plaquettes.c  \
#	OpenAcc/sp_plaquettes.h  \
#	OpenAcc/sp_rettangoli.c  \
#	OpenAcc/sp_rettangoli.h  \
#	OpenAcc/sp_single_types.h  \
#	OpenAcc/sp_stouting.c  \
#	OpenAcc/sp_stouting.h  \
#	OpenAcc/sp_struct_c_def.h  \
#	OpenAcc/sp_su3_measurements.c  \
#	OpenAcc/sp_su3_measurements.h  \
#	OpenAcc/sp_su3_utilities.c  \
#	OpenAcc/sp_su3_utilities.h  \
#	DbgTools/sp_dbgtools.h  \
#	DbgTools/sp_dbgtools.c  \
#	Mpi/sp_communications.h  \
#	Mpi/sp_communications.c 
#
#BUILT_SOURCES = $(__sp_sources)
#
#__dp_sources = OpenAcc/alloc_vars.c  \
#	OpenAcc/alloc_vars.h  \
#	OpenAcc/backfield.c  \
#	OpenAcc/backfield.h  \
#	OpenAcc/cayley_hamilton.h  \
#	OpenAcc/fermion_force.c  \
#	OpenAcc/fermion_force.h  \
#	OpenAcc/fermion_force_utilities.c  \
#	OpenAcc/fermion_force_utilities.h  \
#	OpenAcc/fermionic_utilities.c  \
#	OpenAcc/fermionic_utilities.h  \
#	OpenAcc/fermion_matrix.c  \
#	OpenAcc/fermion_matrix.h  \
#	OpenAcc/inverter_full.c  \
#	OpenAcc/inverter_full.h  \
#	OpenAcc/inverter_multishift_full.c  \
#	OpenAcc/inverter_multishift_full.h  \
#	OpenAcc/ipdot_gauge.c  \
#	OpenAcc/ipdot_gauge.h  \
#	OpenAcc/matvecmul.h  \
#	OpenAcc/md_integrator.c  \
#	OpenAcc/md_integrator.h  \
#	OpenAcc/plaquettes.c  \
#	OpenAcc/plaquettes.h  \
#	OpenAcc/rettangoli.c  \
#	OpenAcc/rettangoli.h  \
#	OpenAcc/single_types.h  \
#	OpenAcc/stouting.c  \
#	OpenAcc/stouting.h  \
#	OpenAcc/struct_c_def.h  \
#	OpenAcc/su3_measurements.c  \
#	OpenAcc/su3_measurements.h  \
#	OpenAcc/su3_utilities.c  \
#	OpenAcc/su3_utilities.h  \
#	DbgTools/dbgtools.h  \
#	DbgTools/dbgtools.c  \
#	Mpi/communications.h  \
#	Mpi/communications.c 
#
#$(__sp_sources) : $(__dp_sources)
#	./double_to_single_transformer.py autoMode
#

main_SOURCES =	DbgTools/dbgtools.c \
	DbgTools/sp_dbgtools.c \
	Include/debug.c \
	Include/fermion_parameters.c \
	Include/hash.c \
	Include/inverter_tricks.c \
	Include/montecarlo_parameters.c \
	Include/setting_file_parser.c \
	Include/tell_geom_defines.c \
	Meas/baryon_number_utilities.c \
	Meas/ferm_meas.c \
	Meas/gauge_meas.c \
	Meas/magnetic_susceptibility_utilities.c \
	Meas/polyakov.c \
	Mpi/communications.c \
	Mpi/multidev.c \
	Mpi/sp_communications.c \
	OpenAcc/action.c \
	OpenAcc/alloc_settings.c \
	OpenAcc/alloc_vars.c \
	OpenAcc/backfield.c \
	OpenAcc/backfield_parameters.c \
	OpenAcc/binary.c \
	OpenAcc/cooling.c \
	OpenAcc/deviceinit.c \
	OpenAcc/fermion_force.c \
	OpenAcc/fermion_force_utilities.c \
	OpenAcc/fermionic_utilities.c \
	OpenAcc/fermion_matrix.c \
	OpenAcc/fermion_matrix_indirect.c \
	OpenAcc/field_times_fermion_matrix.c \
	OpenAcc/find_min_max.c \
	OpenAcc/float_double_conv.c \
	OpenAcc/geometry.c \
	OpenAcc/inverter_full.c \
	OpenAcc/inverter_mixedp.c \
	OpenAcc/inverter_multishift_full.c \
	OpenAcc/inverter_package.c \
	OpenAcc/inverter_wrappers.c \
	OpenAcc/io.c \
	OpenAcc/ipdot_gauge.c \
	OpenAcc/md_integrator.c \
	OpenAcc/md_parameters.c \
	OpenAcc/plaquettes.c \
	OpenAcc/random_assignement.c \
	OpenAcc/rettangoli.c \
	OpenAcc/sp_alloc_vars.c \
	OpenAcc/sp_backfield.c \
	OpenAcc/sp_fermion_force.c \
	OpenAcc/sp_fermion_force_utilities.c \
	OpenAcc/sp_fermionic_utilities.c \
	OpenAcc/sp_fermion_matrix.c \
	OpenAcc/sp_inverter_full.c \
	OpenAcc/sp_inverter_multishift_full.c \
	OpenAcc/sp_ipdot_gauge.c \
	OpenAcc/sp_md_integrator.c \
	OpenAcc/sp_plaquettes.c \
	OpenAcc/sp_rettangoli.c \
	OpenAcc/sp_stouting.c \
	OpenAcc/sp_su3_measurements.c \
	OpenAcc/sp_su3_utilities.c \
	OpenAcc/stouting.c \
	OpenAcc/su3_measurements.c \
	OpenAcc/su3_utilities.c \
	OpenAcc/update_versatile.c \
	Rand/random.c \
	RationalApprox/rationalapprox_andreoli.c \
	RationalApprox/rationalapprox.c \
	tests_and_benchmarks/test_and_benchmarks.c \
	OpenAcc/main.c 
