#ifndef CUDA_INVERTER_H
#define CUDA_INVERTER_H

extern "C" void cuda_inverter(const double residual, 
			              int *ncount);

extern "C" void cuda_inverter_d(const double residual, 
			              int *ncount);

#endif
