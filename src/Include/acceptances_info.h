//
//  acceptances_info.h
//  
//
//  Created by Luca Parente on 27/08/21.
//

#ifndef acceptances_info_h
#define acceptances_info_h


typedef struct accept_info_t{
    
    char hmc_file_name[500];
    char swap_file_name[500];
    char file_label_name[500];
   
    
    
}accept_info;
extern accept_info *acc_info;



#endif /* acceptances_info_h */
