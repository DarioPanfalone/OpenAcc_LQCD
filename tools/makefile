
rationalapprox.o: ../src/RationalApprox/rationalapprox.c ../src/RationalApprox/rationalapprox.h
	g++ -O3 -c ../src/RationalApprox/rationalapprox.c -DNOMULTI  -DLOC_N0=4  -DLOC_N1=4 -DLOC_N2=4 -DLOC_N3=4 -DNRANKS_D3=1

normalize_ratapprox.o: normalize_ratapprox.c ../src/RationalApprox/rationalapprox.h
	g++ -O3 -c normalize_ratapprox.c -DNOMULTI  -DLOC_N0=4  -DLOC_N1=4 -DLOC_N2=4 -DLOC_N3=4 -DNRANKS_D3=1

ratnormalizer: normalize_ratapprox.o rationalapprox.o
	g++ -O3 -o ratnormalizer normalize_ratapprox.o rationalapprox.o -lm 


rationalapprox_calc.o: ../src/RationalApprox/rationalapprox.h  ../src/RationalApprox/rationalapprox_calc.cc ../src/RationalApprox/Remez/alg_remez.h ../src/RationalApprox/Remez/bigfloat.h
	g++ -O3  -c ../src/RationalApprox/rationalapprox_calc.cc -DNOMULTI  -DLOC_N0=4  -DLOC_N1=4 -DLOC_N2=4 -DLOC_N3=4 -DNRANKS_D3=1

rationalapprox_calc_fixed_order.o: ../src/RationalApprox/rationalapprox.h  ../src/RationalApprox/rationalapprox_calc_fixed_order.cc ../src/RationalApprox/Remez/alg_remez.h ../src/RationalApprox/Remez/bigfloat.h
	g++ -O3 -c ../src/RationalApprox/rationalapprox_calc_fixed_order.cc -DNOMULTI  -DLOC_N0=4  -DLOC_N1=4 -DLOC_N2=4 -DLOC_N3=4 -DNRANKS_D3=1

alg_remez.o: ../src/RationalApprox/Remez/alg_remez.cc ../src/RationalApprox/Remez/alg_remez.h ../src/RationalApprox/Remez/bigfloat.h
	g++ -O3 -c ../src/RationalApprox/Remez/alg_remez.cc  -DNOMULTI  -DLOC_N0=4  -DLOC_N1=4 -DLOC_N2=4 -DLOC_N3=4 -DNRANKS_D3=1

rgen: rationalapprox_calc.o rationalapprox.o alg_remez.o
	g++ -O3 -o rgen rationalapprox_calc.o rationalapprox.o alg_remez.o  -lgmp -lmpfr

rgenfo: rationalapprox_calc_fixed_order.o rationalapprox.o alg_remez.o
	g++ -O3 -o rgenfo rationalapprox_calc_fixed_order.o  rationalapprox.o alg_remez.o -lgmp -lmpfr



eval_ratapprox.o: ../src/RationalApprox/rationalapprox.h eval_ratapprox.c
	g++ -O3 -c eval_ratapprox.c -DNOMULTI  -DLOC_N0=4  -DLOC_N1=4 -DLOC_N2=4 -DLOC_N3=4 -DNRANKS_D3=1

eval_ratapprox: eval_ratapprox.o rationalapprox.o
	g++ -O3 -o eval_ratapprox eval_ratapprox.o rationalapprox.o -lm


clean: 
	rm -rf *.o rgen rgenfo


