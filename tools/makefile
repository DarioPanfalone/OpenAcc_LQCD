
rationalapprox.o: ../src/RationalApprox/rationalapprox.c ../src/RationalApprox/rationalapprox.h
	g++ -c ../src/RationalApprox/rationalapprox.c

normalize_ratapprox.o: normalize_ratapprox.c ../src/RationalApprox/rationalapprox.h
	g++ -c normalize_ratapprox.c

ratnormalizer: normalize_ratapprox.o rationalapprox.o
	g++ -o ratnormalizer normalize_ratapprox.o rationalapprox.o -lm


rationalapprox_calc.o: ../src/RationalApprox/rationalapprox.h  ../src/RationalApprox/rationalapprox_calc.cc ../src/RationalApprox/Remez/alg_remez.h ../src/RationalApprox/Remez/bigfloat.h
	g++ -c ../src/RationalApprox/rationalapprox_calc.cc

alg_remez.o: ../src/RationalApprox/Remez/alg_remez.cc ../src/RationalApprox/Remez/alg_remez.h ../src/RationalApprox/Remez/bigfloat.h
	g++ -c ../src/RationalApprox/Remez/alg_remez.cc

rgen: rationalapprox_calc.o rationalapprox.o alg_remez.o
	g++ -o rgen rationalapprox_calc.o rationalapprox.o alg_remez.o -lgmp -lmpfr


eval_ratapprox.o: ../src/RationalApprox/rationalapprox.h eval_ratapprox.c
	g++ -c eval_ratapprox.c

eval_ratapprox: eval_ratapprox.o rationalapprox.o
	g++ -o eval_ratapprox eval_ratapprox.o rationalapprox.o -lm


clean: 
	rm -rf normalize_ratapprox.o rationalapprox.o ratnormalizer rgen alg_remez.o rationalapprox_calc.o


