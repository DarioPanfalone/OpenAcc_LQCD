#ifndef ALLOC_SETTINGS_H_
#define ALLOC_SETTINGS_H_



typedef struct alloc_settings_t{

    int NDiffFlavs; // set in Include/setting_file_parser.c, from input file 
    int NPS_tot; // set in Include/fermion_parameters.c 
    int conf_acc_size; 
    int maxNeededShifts;
    int maxApproxOrder;

    int revTestAllocations;
    int diagnosticsAllocations;

    int singlePrecCoreAllocations;
    int singlePrecExtendedAllocations;
    int stoutAllocations;
    
    int num_replicas; //numero di repliche dell'hasenbush

}alloc_settings;

extern alloc_settings alloc_info;




#endif
