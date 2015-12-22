#ifndef UPDATE_VERSATILE_H_
#define UPDATE_VERSATILE_H_

#include "./struct_c_def.h"

// se metro==1 allora fa il test di metropolis
// se metro==0 allora non fa il test di metropolis --> termalizzazione


int UPDATE_SOLOACC_UNOSTEP_VERSATILE(su3_soa *tconf_acc,
        double res_metro, double res_md, int id_iter,int acc,int metro);


#endif
