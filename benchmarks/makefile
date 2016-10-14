



# DIMENSIONS OF THE LOCAL LATTICE
N0:=$(shell grep N0 make_settings.set | awk '{print $$2}')
N1:=$(shell grep N1 make_settings.set | awk '{print $$2}')
N2:=$(shell grep N2 make_settings.set | awk '{print $$2}')
N3:=$(shell grep N3 make_settings.set | awk '{print $$2}')
# NUMBER OF RANKS
NR3:=$(shell grep NR3 make_settings.set | awk '{print $$2}')

# COMPILER CHOICE: PGI, GNU or MPI
COMPILER:=$(shell grep COMPILER make_settings.set | awk '{print $$2}')


ROOTDIR=~/OPENACCREPO




all: run/main run/deo_doe_test run/inverter_multishift_test run/benchmark.set run/benchmark.puregauge.set run/rgen make_settings.set

make_settings.set:
	$(error The file make_settings.set must be present.)


lattice_dimensions.txt:
	grep "N0\|N1\|N2\|N3\|NR3" make_settings.set > lattice_dimensions.txt


makefile: make_settings.set lattice_dimensions.txt
	cd $(ROOTDIR)/src && ./double_to_single_transformer.py autoMode && cd - && echo 'MAKE SETTINGS:'  &&  cat make_settings.set && \
	$(ROOTDIR)/build/generate_makefile.py $(ROOTDIR)/src/OpenAcc/main.c $(ROOTDIR)/src/tests_and_benchmarks/deo_doe_test.c $(ROOTDIR)/src/tests_and_benchmarks/inverter_multishift_test.c $(ROOTDIR)/src/tools/checkinput.c  $(ROOTDIR)/src/tools/confrestart.c $(COMPILER) > ./makefile 
	
run:
	if [ ! -d run ] ; then mkdir run ; fi

rgen run/rgen : makefile run
	make COMMIT_HASH=$(COMMIT_HASH) rgen

main run/main: makefile run
	make COMMIT_HASH=$(COMMIT_HASH) main -j32

deo_doe_test run/deo_doe_test: makefile run
	make COMMIT_HASH=$(COMMIT_HASH) deo_doe_test -j32

inverter_multishift_test run/inverter_multishift_test: makefile run
	make COMMIT_HASH=$(COMMIT_HASH) inverter_multishift_test -j32


run/benchmark.set: run/benchmark.puregauge.set fermion_parameters.set run
	cat fermion_parameters.set > run/benchmark.set && cat run/benchmark.puregauge.set >> run/benchmark.set 

run/benchmark.puregauge.set: benchmark.set.proto run

	sed 's/SEDNRANKS/'$(NR3)'/' benchmark.set.proto | sed 's/SEDNX/'$(N0)'/' |\
	   	sed 's/SEDNY/'$(N1)'/' | sed 's/SEDNZ/'$(N2)'/' | sed 's/SEDNT/'$$(($(N3)*$(NR3)))'/' > run/benchmark.puregauge.set 


test: testdeodoe testmultishift testmain testmainpg

testdeodoe: run/deo_doe_test run/benchmark.set
	cd run && mpirun -n $(NR3) ./deo_doe_test ./benchmark.set | tee -a  deo_doe_test_output.txt

testmultishift: run/inverter_multishift_test run/benchmark.set
	cd run && mpirun -n $(NR3) ./inverter_multishift_test ./benchmark.set | tee -a  inverter_multishift_test_output.txt

testmainpg: run/main run/benchmark.puregauge.set
	cd run && mpirun -n $(NR3) ./main ./benchmark.puregauge.set | tee -a  puregauge_output.txt

testmain: run/main run/benchmark.set run/rgen
	cd run && mpirun -n $(NR3) ./main ./benchmark.set | tee -a  full_md_output.txt

clean: 
	cd run && rm -f main rgen deo_doe_test inverter_multishift_test benchmark.set benchmark.puregauge.set 	&& cd ..
