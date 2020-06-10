#!/bin/bash

/afs/pi.infn.it/project/hpc/software/CUDA/CUDA6/bin/nvcc -o test_device ./DeviceInfo/test_device.cu

./test_device

