//
//  defect_info.h
//  
//
//  Created by Luca Parente on 01/10/21.
//

#ifndef defect_info_h
#define defect_info_h




typedef struct defect_info_t{
    
    int ** defect_swap_min;
    int ** defect_swap_max;
    #ifdef GAUGE_ACT_TLSM
    int ** defect_swap_min_TLSM;
    int ** defect_swap_max_TLSM;
    #endif
    int def_axis_mapped;
    int nu_vector[3];

    
    
}defect_info;
extern defect_info *def;

#endif /* defect_info_h */

