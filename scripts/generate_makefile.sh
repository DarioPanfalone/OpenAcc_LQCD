#!/bin/bash
# launch this from the build directory
./generate_makefile.py $@ ../src/{OpenAcc,RationalApprox,Include,Rand,Meas}/*.[ch] > makefile

