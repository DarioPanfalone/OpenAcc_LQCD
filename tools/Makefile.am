
AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = rgen rgenfo ratnormalizer

__common_sources = \
	./rationalapprox.h\
	./rationalapprox.cpp\
	../src/RationalApprox/Remez/alg_remez.h\
	../src/RationalApprox/Remez/bigfloat.h

rgen_SOURCES = ./rationalapprox_calc.cc\
	../src/RationalApprox/Remez/alg_remez.cc\
	$(__common_sources)

rgenfo_SOURCES = ./rationalapprox_calc_fixed_order.cc\
	../src/RationalApprox/Remez/alg_remez.cc\
	$(__common_sources)

ratnormalizer_SOURCES = normalize_ratapprox.cc\
	$(__common_sources)

