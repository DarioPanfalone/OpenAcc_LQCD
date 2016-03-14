#ifndef DEVICE_INIT_H_
#define DEVICE_INIT_H_

typedef struct device_param_t{

  int device_choice;


} device_param;

extern device_param dev_settings;



#ifndef __GNUC__  // assuming then PGI is used for compilation on accelerator
       // at the moment only nvidia gpus are supported.

#include "openacc.h"
void select_init_acc_device(acc_device_t my_device_type, int dev_index);
void shutdown_acc_device(acc_device_t my_device_type);


#endif

#endif
