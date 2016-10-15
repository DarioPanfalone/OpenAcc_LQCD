#ifndef _HASH_C_
#define _HASH_C_

#include "./hash.h"

#include "./montecarlo_parameters.h"
#include "./fermion_parameters.h"
#include "../OpenAcc/action.h"
#include "../OpenAcc/md_integrator.h"
#include "../OpenAcc/backfield.h"
#include "../OpenAcc/geometry.h"
#include "../OpenAcc/alloc_settings.h"
#include "../Meas/ferm_meas.h"

#include <stdint.h>

// modified version of
// https://en.wikipedia.org/wiki/Jenkins_hash_function
void hashfun(uint32_t *hash, char *key, size_t len)
{
    uint32_t i;
    for(i = 0; i < len; ++i)
    {
        *hash += key[i];
        *hash += (*hash << 10);
        *hash ^= (*hash >> 6);
    }
    *hash += (*hash << 3);
    *hash ^= (*hash >> 11);
    *hash += (*hash << 15);
}


uint32_t hash_settings_explicit(
        ferm_param *flpar,
        int nflavs,
        action_param *act_par,
        bf_param *bfpar,
        md_param *mdpar,
        mc_params_t *mcpar,
        geom_parameters *gpar,
        ferm_meas_params *fmpar)
{


    uint32_t hash=0;

    // fermion parameters
    int iflav;
    for(iflav=0; iflav < nflavs; iflav++){

        hashfun(&hash,(char*)&(flpar[iflav].ferm_mass       ),sizeof(double));
        hashfun(&hash,(char*)&(flpar[iflav].degeneracy      ),sizeof(int));
        hashfun(&hash,(char*)&(flpar[iflav].number_of_ps    ),sizeof(int));
        hashfun(&hash,(char*)&(flpar[iflav].ferm_charge     ),sizeof(double));
        hashfun(&hash,(char*)&(flpar[iflav].ferm_im_chem_pot),sizeof(double));

    }
    // action parameters
    hashfun(&hash,(char*)&(act_par->beta),sizeof(double));
    hashfun(&hash,(char*)&(act_par->stout_steps),sizeof(int));
    hashfun(&hash,(char*)&(act_par->stout_rho),sizeof(double));
    // backfield info
    hashfun(&hash,(char*)&(bfpar->ex),sizeof(double));
    hashfun(&hash,(char*)&(bfpar->ey),sizeof(double));
    hashfun(&hash,(char*)&(bfpar->ez),sizeof(double));
    hashfun(&hash,(char*)&(bfpar->bx),sizeof(double));
    hashfun(&hash,(char*)&(bfpar->by),sizeof(double));
    hashfun(&hash,(char*)&(bfpar->bz),sizeof(double));
    // molecular dynamics info
    hashfun(&hash,(char*)&(mdpar->no_md      ),sizeof(int));
    hashfun(&hash,(char*)&(mdpar->gauge_scale),sizeof(int));
    hashfun(&hash,(char*)&(mdpar->t          ),sizeof(double));
    hashfun(&hash,(char*)&(mdpar->residue_md ),sizeof(double));
    hashfun(&hash,(char*)&(mdpar->residue_metro),sizeof(double));
    hashfun(&hash,(char*)&(mdpar->singlePrecMD),sizeof(int));
    hashfun(&hash,(char*)&(mdpar->expected_max_eigenvalue),
            sizeof(double));
    // ONLY RELEVANT geometry info
    hashfun(&hash,(char*)&(gpar->gnx),sizeof(int));
    hashfun(&hash,(char*)&(gpar->gny),sizeof(int));
    hashfun(&hash,(char*)&(gpar->gnz),sizeof(int));
    hashfun(&hash,(char*)&(gpar->gnt),sizeof(int));
    // ONLY RELEVANT mcpar info
    hashfun(&hash,(char*)&(mcpar->JarzynskiMode),sizeof(int));
    // ONLY RELEVANT fermion measurements info
    hashfun(&hash,(char*)&(fmpar->measEvery),sizeof(int));
    hashfun(&hash,(char*)&(fmpar->SingleInvNVectors),sizeof(int));
    hashfun(&hash,(char*)&(fmpar->DoubleInvNVectorsChiral),sizeof(int));
    hashfun(&hash,(char*)&(fmpar->DoubleInvNVectorsQuarkNumber),sizeof(int));

    return hash;

}


uint32_t hash_settings(){

    return hash_settings_explicit(fermions_parameters,alloc_info.NDiffFlavs,
                &act_params,&backfield_parameters,&md_parameters,
                &mc_params,&geom_par,&fm_par);


}






#endif
