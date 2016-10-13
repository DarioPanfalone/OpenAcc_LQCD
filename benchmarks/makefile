

CHOICE=GNU # or PGI, or MPI

# DIMENSIONS OF THE LOCAL LATTICE
L=32
T=8

# NUMBER OF RANKS
NR3=1


all: run/main run/deo_doe_test run/inverter_multishift_test run/benchmark.set run/benchmark.puregauge.set

../build/makefile: ../build/lattice_dimensions.txt
	cd ../src && ./double_to_single_transformer.py autoMode && cd ../build && echo 'LATTICE DIMENSIONS:'  &&  echo PORCODIO && cat ../build/lattice_dimensions.txt && \
	./generate_makefile.py ../src/OpenAcc/main.c ../src/tests_and_benchmarks/deo_doe_test.c ../src/tests_and_benchmarks/inverter_multishift_test.c ../src/tools/checkinput.c  ../src/tools/confrestart.c $(CHOICE) > ../build/makefile 
	
run:
	if [ ! -d run ] ; then mkdir run ; fi

../build/lattice_dimensions.txt:
	echo N0 $(L) "\n"N1 $(L) "\n"N2 $(L) "\n"N3 $(T) "\n"NR3 $(NR3) > ../build/lattice_dimensions.txt

run/rgen run/main: ../build/makefile run
	cd ../build && make COMMIT_HASH=$(COMMIT_HASH) main &&cp run/main run/rgen ../benchmarks/run && cd ../benchmarks

run/deo_doe_test: ../build/makefile run
	cd ../build && make COMMIT_HASH=$(COMMIT_HASH) deo_doe_test && cp deo_doe_test ../benchmarks/run && cd ../benchmarks

run/inverter_multishift_test: ../build/makefile run
	cd ../build && make COMMIT_HASH=$(COMMIT_HASH) inverter_multishift_test && cp inverter_multishift_test ../benchmarks/run && cd ../benchmarks


run/benchmark.set: benchmark.set.proto fermion_parameters.set run
	cat fermion_parameters.set > run/benchmark.set && sed 's/SEDNRANKS/'$(NR3)'/' benchmark.set.proto | sed 's/SEDL/'$(L)'/' | sed 's/SEDT/'$(T)'/' >> run/benchmark.set 

run/benchmark.puregauge.set: benchmark.set.proto run

	sed 's/SEDNRANKS/'$(NR3)'/' benchmark.set.proto | sed 's/SEDL/'$(L)'/' | sed 's/SEDT/'$(T)'/' > run/benchmark.puregauge.set 


test: all
	cd run && ./deo_doe_test ./benchmark.set | tee  deo_doe_test_output.txt
	cd run && ./inverter_multishift_test ./benchmark.set | tee  inverter_multishift_test_output.txt
	cd run && ./main ./benchmark.puregauge.set | tee  puregauge_output.txt
	cd run && ./main ./benchmark.set | tee  full_md_output.txt

clean: 
	cd run && rm -f main rgen deo_doe_test inverter_multishift_test benchmark.set benchmark.puregauge.set && cd ..
