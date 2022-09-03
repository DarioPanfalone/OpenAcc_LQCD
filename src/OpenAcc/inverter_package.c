#ifndef INVERTER_PACKAGE_C_
#define INVERTER_PACKAGE_C_

#include "./sp_struct_c_def.h"
#include "./struct_c_def.h"
#include "./inverter_package.h"
#include <stdio.h>
#include <stdlib.h>

// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
#define __restrict
#endif



void setup_inverter_package_dp(inverter_package * ip, 
															 su3_soa * u, 
															 vec3_soa * ferm_shift_temp,
															 int nshifts,
															 vec3_soa * loc_r,
															 vec3_soa * loc_h,
															 vec3_soa * loc_s,
															 vec3_soa * loc_p)
{
	ip->u =  u;
	ip->ferm_shift_temp = ferm_shift_temp;
	ip->nshifts = nshifts;
	ip->loc_r =  loc_r;
	ip->loc_h =  loc_h;
	ip->loc_s =  loc_s;
	ip->loc_p =  loc_p;

	vec3_soa * allptrs[] = {loc_r,loc_h,loc_s, loc_p,ferm_shift_temp};
	int nptrs = sizeof(allptrs)/sizeof(vec3_soa*);
	int i,j;
	// check 
	for(i=0;i<nptrs;i++) for(j=i+1;j<nptrs;j++)
												 if(allptrs[i] == allptrs[j] && 0 != allptrs[i]){
													 printf("BAD SETUP OF INVERTER PACKAGE! (%s:%d)\n", __FILE__, __LINE__);
													 printf("Pointer %p used twice (%d == %d).\n",allptrs[i],i,j);
													 exit(1);
												 }

}

void setup_inverter_package_sp(inverter_package * ip, 
															 su3_soa_f * u_f, 
															 vec3_soa_f * ferm_shift_temp_f,
															 int nshifts,
															 vec3_soa_f * loc_r_f,
															 vec3_soa_f * loc_h_f,
															 vec3_soa_f * loc_s_f,
															 vec3_soa_f * loc_p_f,
															 vec3_soa_f * out_f)
{
	ip->u_f =  u_f;
	ip->ferm_shift_temp_f = ferm_shift_temp_f;
	ip->nshifts = nshifts;
	ip->loc_r_f =  loc_r_f;
	ip->loc_h_f =  loc_h_f;
	ip->loc_s_f =  loc_s_f;
	ip->loc_p_f =  loc_p_f;
	ip->out_f   =    out_f;
	vec3_soa_f * allptrs[] = {loc_r_f,loc_h_f,loc_s_f, loc_p_f,ferm_shift_temp_f, out_f};
    
	int nptrs = sizeof(allptrs)/sizeof(vec3_soa_f*);
	int i,j;
	// check 
	for(i=0;i<nptrs;i++) for(j=i+1;j<nptrs;j++)
												 if(allptrs[i] == allptrs[j] && 0 != allptrs[i]){
													 printf("BAD SETUP OF INVERTER PACKAGE! (%s:%d)\n", __FILE__, __LINE__);
													 printf("Pointer %p used twice (%d == %d).\n",allptrs[i],i,j);
													 exit(1);
												 }

}


#endif
