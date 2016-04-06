
rationalapprox.o: ../src/RationalApprox/rationalapprox.c ../src/RationalApprox/rationalapprox.h
	g++ -O3 -c ../src/RationalApprox/rationalapprox.c -DNOMULTI

normalize_ratapprox.o: normalize_ratapprox.c ../src/RationalApprox/rationalapprox.h
	g++ -O3 -c normalize_ratapprox.c

ratnormalizer: normalize_ratapprox.o rationalapprox.o
	g++ -O3 -o ratnormalizer normalize_ratapprox.o rationalapprox.o -lm


rationalapprox_calc.o: ../src/RationalApprox/rationalapprox.h  ../src/RationalApprox/rationalapprox_calc.cc ../src/RationalApprox/Remez/alg_remez.h ../src/RationalApprox/Remez/bigfloat.h
	g++ -O3  -c ../src/RationalApprox/rationalapprox_calc.cc

rationalapprox_calc_fixed_order.o: ../src/RationalApprox/rationalapprox.h  ../src/RationalApprox/rationalapprox_calc_fixed_order.cc ../src/RationalApprox/Remez/alg_remez.h ../src/RationalApprox/Remez/bigfloat.h
	g++ -O3 -c ../src/RationalApprox/rationalapprox_calc_fixed_order.cc

alg_remez.o: ../src/RationalApprox/Remez/alg_remez.cc ../src/RationalApprox/Remez/alg_remez.h ../src/RationalApprox/Remez/bigfloat.h
	g++ -O3 -c ../src/RationalApprox/Remez/alg_remez.cc

rgen: rationalapprox_calc.o rationalapprox.o alg_remez.o
	g++ -O3 -o rgen rationalapprox_calc.o rationalapprox.o alg_remez.o -lgmp -lmpfr

rgenfo: rationalapprox_calc_fixed_order.o rationalapprox.o alg_remez.o
	g++ -O3 -o rgenfo rationalapprox_calc_fixed_order.o  rationalapprox.o alg_remez.o -lgmp -lmpfr



eval_ratapprox.o: ../src/RationalApprox/rationalapprox.h eval_ratapprox.c
	g++ -O3 -c eval_ratapprox.c

eval_ratapprox: eval_ratapprox.o rationalapprox.o
	g++ -O3 -o eval_ratapprox eval_ratapprox.o rationalapprox.o -lm


clean: 
	rm -rf *.o rgen rgenfo


