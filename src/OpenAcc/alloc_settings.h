#ifndef ALLOC_SETTINGS_H_
#define ALLOC_SETTINGS_H_

typedef struct alloc_settings_t{

    int NDiffFlavs;
    int NPS_tot;
    int conf_acc_size; 
    int maxNeededShifts;
    int maxApproxOrder;

    int revTestAllocations;
    int diagnosticsAllocations;

    int singlePrecCoreAllocations;
    int singlePrecExtendedAllocations;
    int stoutAllocations;
    
    int num_replicas;

}alloc_settings;

extern alloc_settings alloc_info;

#endif
