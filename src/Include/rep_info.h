#ifndef rep_info_h
#define rep_info_h



typedef struct rep_info_t{
    
    int replicas_total_number;
    int defect_boundary;
    int defect_coordinates[3];
    double *cr_vet;
		int *label;   
    
}rep_info;
extern rep_info *rep;


#endif /* rep_info_h */
