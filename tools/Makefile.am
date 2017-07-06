
AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = rgen rgenfo eval_ratapprox ratnormalizer

__common_sources = \
	../src/RationalApprox/rationalapprox.c\
	../src/RationalApprox/rationalapprox.h\
	../src/RationalApprox/Remez/alg_remez.h\
	../src/RationalApprox/Remez/bigfloat.h

rgen_SOURCES = ../src/RationalApprox/Remez/alg_remez.cc\
	../src/RationalApprox/rationalapprox_calc.cc\
	$(__common_sources)

rgenfo_SOURCES = ../src/RationalApprox/Remez/alg_remez.cc\
	../src/RationalApprox/rationalapprox_calc_fixed_order.cc\
	$(__common_sources)

eval_ratapprox_SOURCES = eval_ratapprox.c\
	$(__common_sources)

ratnormalizer_SOURCES = normalize_ratapprox.c\
	$(__common_sources)

