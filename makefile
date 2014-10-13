
GLOBAL_STUFF=Include/global_const.cc Include/global_macro.cc Include/global_var.cc Include/parameters.cc
COMMON_STUFF=Fermions/fermions.cc Init/init.cc Exception/exception.cc Su3/su3.cc Vec3/vec3.cc Conf/conf.cc Geometry/geometry.cc Gauss/gauss.cc 
CPU_ONLY_STUFF=Inverter/inverter.cc FermionMatrix/fermionmatrix.cc
GPU_ONLY_STUFF=Cuda/cuda_init_all.cu Packer/packer.cc DeviceInfo/device_info.cu
GPU_ONLY_STUFF_ONLY_FLOAT=Cuda/cuda_dslash_eo.cu 
GPU_ONLY_STUFF_ONLY_DOUBLE=Cuda/cuda_dslash_dd_eo.cu 

CC=g++
NVCC=nvcc
ARCH=sm_20

all: moltestCPU moltestGPU moltestGPU_D

moltestGPU: $(GLOBAL_STUFF) $(COMMON_STUFF) $(GPU_ONLY_STUFF) $(GPU_ONLY_STUFF_ONLY_FLOAT) moltest.cu
	$(NVCC) -O3 -arch $(ARCH) -o moltest_CU moltest.cu

moltestCPU: $(GLOBAL_STUFF) $(COMMON_STUFF) $(CPU_ONLY_STUFF) moltest.cpp
	$(CC) -O3 -o moltestCPU moltest.cpp

moltestGPU_D:$(GLOBAL_STUFF) $(COMMON_STUFF) $(GPU_ONLY_STUFF) $(GPU_ONLY_STUFF_ONLY_DOUBLE) moltest.cu
	$(NVCC) -O3 -arch $(ARCH) -DDOUBLE_PREC -o moltest_CU_D moltest.cu


