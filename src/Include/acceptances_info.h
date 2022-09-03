#ifndef ACC_INFO_H
#define ACC_INFO_H

typedef struct accept_info_t {
	char hmc_file_name[500];
	char swap_file_name[500];
	char file_label_name[500];
} accept_info;

extern accept_info *acc_info;

#endif
