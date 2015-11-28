
rationalapprox.o: ../src/RationalApprox/rationalapprox.c
	gcc -c -std=c99 ../src/RationalApprox/rationalapprox.c

normalize_ratapprox.o: normalize_ratapprox.c
	gcc -c normalize_ratapprox.c

ratnormalizer: normalize_ratapprox.o rationalapprox.o
	gcc -o ratnormalizer normalize_ratapprox.o rationalapprox.o -lm

clean: 
	rm -rf normalize_ratapprox.o rationalapprox.o ratnormalizer

