#ifndef _SETTING_FILE_PARSER_C_
#define _SETTING_FILE_PARSER_C_

#include "../Meas/ferm_meas.h"
#include "../Meas/gauge_meas.h"
#include "../Mpi/multidev.h"
#include "../OpenAcc/action.h"
#include "../OpenAcc/backfield.h"
#include "../OpenAcc/geometry.h"
#include "../OpenAcc/geometry.h" // for MULTIDEVICE to be defined or not
#include "../OpenAcc/md_integrator.h"
#include "../OpenAcc/su3_measurements.h"
#include "../RationalApprox/rationalapprox.h"
#include "./common_defines.h"
#include "./debug.h"
#include "./fermion_parameters.h"
#include "./hash.h"
#include "./montecarlo_parameters.h"
#include "./setting_file_parser.h"

#include <stdio.h>
#include <strings.h>
#include <stdlib.h>

#include <inttypes.h>
#include <stdint.h>

#define MAXLINES 300
#define MAXLINELENGTH 500 // pretty long to accomodate all the comments

char input_file_str[MAXLINES*MAXLINELENGTH];

//types that can be found in an input file
#define TYPE_INT 0
#define TYPE_DOUBLE 1
#define TYPE_STR 2
const char * type_strings[]={"(int)", "(double)", "(string)" };
typedef struct par_info_t{

    void* par;
    int type;
    char* name;
    int is_optional;
    const void* default_value;

}par_info;

#define NPMGTYPES 10
#define MAXPMG 20
const char * par_macro_groups_names[] ={ // NPMTYPES strings, NO SPACES!
    "ActionParameters",         // 0 
    "FlavourParameters",        // 1   
    "BackgroundFieldParameters",// 2            
    "MDParameters",             // 3
    "MontecarloParameters",     // 4       
    "GaugeMeasuresSettings",    // 5
    "FermionMeasuresSettings",  // 6
    "DeviceSettings"         ,  // 7
    "Geometry"               ,  // 8
    "DebugSettings"             // 9
};
#define PMG_ACTION        0
#define PMG_FERMION       1
#define PMG_BACKGROUND    2
#define PMG_MD            3
#define PMG_MC            4
#define PMG_GMEAS         5
#define PMG_FMEAS         6
#define PMG_DEVICE        7
#define PMG_GEOMETRY      8
#define PMG_DEBUG         9
// last number should be NPMGTYPES - 1 !!


FILE * helpfile;
char IGNORE_IT[50];

// just to save it in the conf
int prepare_string_from_stringarray(char file_lines[MAXLINES][MAXLINELENGTH], 
        int maxlines,
        char* input_filename_to_write){

    int iline;int totlen = 0;
    for(iline = 0; iline < maxlines ; iline++){
        strcat(input_filename_to_write,file_lines[iline]);
        totlen += strlen(file_lines[iline]);
    }
    return totlen;

}

void erase_comments(char file_lines[MAXLINES][MAXLINELENGTH], int maxlines)
{
    for(int i = 0; i<maxlines; i++){
        char* start_comment_ptr = index(file_lines[i],'#' );
        if(start_comment_ptr){
            int comment_start = start_comment_ptr - file_lines[i];
            file_lines[i][comment_start] = '\n';
            file_lines[i][comment_start+1] = '\0';
            for(int j = comment_start+2; j < MAXLINELENGTH; j++)
                file_lines[i][j] = ' ';
        }
    }
}

void reorder_par_infos(int npar, par_info * par_infos ){
    // this is necessary if a parameter name is contained in another
    // parameter name, for example "Ntraj" is contained in "ThermNtraj".
    // If this happens, if the parser looks first for "Ntraj", we'll have 
    // 2 matches for "Ntraj" and 0 matches for "ThermNtraj" (see scan_group_NV()).
    // We are safe instead if "ThermNtraj" comes first.
    // So, if i<j and name[i] is contained in name[j], the two parameters
    // are exchanged.

    int allright = 0;
    while(! allright && npar > 1){
        for(int i = 0; i < npar; i++){
            for(int j = i+1; j < npar; j++){
                allright = 1;
                char* match = strstr(par_infos[j].name,par_infos[i].name);
                // if(0==devinfo.myrank) printf("DEBUG (%d,%d)params %s and %s, match %p \n",i,j,par_infos[j].name,par_infos[i].name, match );
                if (match){
                    if(0==devinfo.myrank)
                        printf("  <<Reordering %s and %s>>\n",par_infos[j].name,par_infos[i].name);
                    par_info tmp = par_infos[i];
                    par_infos[i] = par_infos[j];
                    par_infos[j] = tmp;
                    allright = 0;
                    break;
                }           
            }
        }
        if(allright) break;
    }
}


int scan_group_V(int ntagstofind, const char **strtofind, 
        int *tagcount,
        int *taglines, int *tagtypes, int maxnres,
        char filelines[MAXLINES][MAXLINELENGTH], 
        int startline, int endline)
{
    // returns number of groups found, writes arrays
    // scans for lines in the format
    // (stuff that will be ignored )NAME (stuff that will be ignored)
    // the arrays 'taglines' and 'types' and must be maxnres long

    if(0==devinfo.myrank)
        printf("Scanning for parameter macro groups..\n");
    for(int i =0; i<ntagstofind; i++) tagcount[i] = 0;
    for(int i =0; i<maxnres; i++){// initializing output arrays
        taglines[i] = -1;
        tagtypes[i] = -1;
    }
    int nres = 0;
    int iline = startline;
    while(iline<endline && nres < maxnres ){
        char * found_something;
        for(int itype =0; itype <ntagstofind; itype++){
            found_something = strstr(filelines[iline],strtofind[itype]);
            if(found_something){
                //printf("Found group %s on line %d\n",strtofind[itype], iline);
                taglines[nres] = iline;
                tagtypes[nres] = itype;
                nres++;
                tagcount[itype]++;
            }
        } 
        iline++;
    }
    for(int i =0; i<ntagstofind; i++)
        if(0==devinfo.myrank)
            printf("Found %s %d times.\n", strtofind[i], tagcount[i]); 

    return nres;
}

int scan_group_NV(int npars,par_info* par_infos,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{   
    if(startline >= endline){ // goes into 'help mode'
        for(int ipar = 0;ipar< npars ; ipar++){
            if(par_infos[ipar].is_optional) switch(par_infos[ipar].type){
                case TYPE_INT:
                    if(0==devinfo.myrank) fprintf(
                            helpfile,"%-30s%-20d#%s\n",
                            par_infos[ipar].name,
                            *((const int *) par_infos[ipar].default_value),
                            type_strings[par_infos[ipar].type]);
                    break;
                case TYPE_DOUBLE:
                    if(0==devinfo.myrank) fprintf(
                            helpfile,"%-30s%-20e#%s\n",
                            par_infos[ipar].name,
                            *((const double *) par_infos[ipar].default_value),
                            type_strings[par_infos[ipar].type]);
                    break;
                case TYPE_STR:
                    if(0==devinfo.myrank) fprintf(
                            helpfile,"%-30s%-20s#%s\n",
                            par_infos[ipar].name,
                            ((const char *) par_infos[ipar].default_value),
                            type_strings[par_infos[ipar].type]);


                    break;
            }
            else 
                if(0==devinfo.myrank) fprintf(
                        helpfile,"%-50s#%s\n",
                        par_infos[ipar].name,
                        type_strings[par_infos[ipar].type]);

        }
        return 0;
    }
    else // 'normal mode'
    {   // scans lines in the form 
        // NAME VALUE  +stuff which will be ignored
        // rc = readcheck 

        // necessary step if a parameter name is contained 
        // in another parameter name.

        reorder_par_infos(npars, par_infos);


        int * rc = (int *) malloc(npars*sizeof(int));
        for(int i =0; i<npars; i++) rc[i] = 0;

        int res = 0 ;

        int iline = startline+1;

        while(! res && iline < endline)
        {
            char * found_something;
            for(int i =0; i<npars; i++){
                found_something = strstr(filelines[iline],par_infos[i].name);
                if(found_something){ // looks at the beginning of the line.
                    // found parameter
                    if(par_infos[i].par == IGNORE_IT){
                        if(0==devinfo.myrank)
                            printf("WARNING, LINE %-3d: IGNORING %s\n",
                                    iline+1,par_infos[i].name);
                        break;
                    }
                    else {
                        if(0==devinfo.myrank)
                            printf("%-3d  %s\r\t\t\t\t ",iline+1,par_infos[i].name);
                        int reads = 0;
                        char parname[50];

                        switch(par_infos[i].type){
                            case TYPE_INT: 
                                reads = sscanf(filelines[iline],
                                        "%s %d",parname,(int*)par_infos[i].par);
                                if(reads == 2) 
                                    if(0==devinfo.myrank)
                                        printf("%d\n", *((int*)par_infos[i].par));
                                break;
                            case TYPE_DOUBLE:

                                reads = sscanf(filelines[iline],
                                        "%s %lf",parname,(double*) par_infos[i].par);
                                if(reads == 2)
                                    if(0==devinfo.myrank)
                                        printf("%e\n", *((double*)par_infos[i].par));
                                break;
                            case TYPE_STR: 
                                reads = sscanf(filelines[iline],
                                        "%s %s",parname,(char*) par_infos[i].par);
                                if(reads == 2) 
                                    if(0==devinfo.myrank)
                                        printf("\"%s\"\n", ((char*)par_infos[i].par));
                                break;
                            default: 
                                if(0==devinfo.myrank)
                                    printf("WARNING, variable type not set in sourcecode.\n");
                                break;

                        }
                        if(reads == 2)rc[i]++;
                        else if(0==devinfo.myrank) printf("WARNING, NO VALUE READ!");
                        break;
                    }
                };
            }
            if(!found_something){
                char word[50];
                int reads = sscanf(filelines[iline],"%s", word);
                if(reads==1){
                    if(0==devinfo.myrank)
                        printf("line: %d, ERROR, parameter %s not recognized\n",iline+1,word);
                    printf("%s\n", filelines[iline]);
                    return 1;
                }
            }


            iline++; 
            res = 1;
            for(int i =0; i<npars; i++)
                res = res && rc[i];
        }



        if(! res){
            res = 1;
            for(int i =0; i<npars; i++){
                if (rc[i]==0){
                    if(0==devinfo.myrank)
                        printf("  (NOTICE) Parameter %s not set in input file.",par_infos[i].name);
                    if(par_infos[i].is_optional==1){
                        if(0==devinfo.myrank)
                            printf(" Parameter is optional. Default value: ");
                        switch(par_infos[i].type){
                            case TYPE_INT:
                                *((int*)(par_infos[i].par)) = 
                                    *((const int*) (par_infos[i].default_value));
                                if(0==devinfo.myrank)
                                    printf(" %d", *((int*)(par_infos[i].par)));
                                break;
                            case TYPE_DOUBLE:
                                *((double*)(par_infos[i].par)) = 
                                    *((const double*) (par_infos[i].default_value));
                                if(0==devinfo.myrank)
                                    printf(" %f", *((double*)(par_infos[i].par)));
                                break;
                            case TYPE_STR:
                                strcpy((char*)par_infos[i].par,(const char*) par_infos[i].default_value);
                                if(0==devinfo.myrank)
                                    printf(" %s", (char*)(par_infos[i].par));
                                break;
                        }
                        if(0==devinfo.myrank) printf("\n");
                    }
                    else{
                        if(0==devinfo.myrank)
                            printf("\n\nERROR: Parameter is NOT optional!\n\n");
                        res = 0;
                    }
                }
            }
        } // if(!res)
        free(rc);
        if(!res){

            if(0==devinfo.myrank)
                printf("NON OPTIONAL PARAMETERS NOT FOUND IN INPUT FILE. PLEASE CHECK.\n");
            return 1;
        }
        else return 0;
    }// else normal mode
}

int read_flavour_info(ferm_param *flpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    // see /Include/fermion_parameters.h

    par_info fp[]={
        (par_info){(void*) &(flpar->ferm_mass       ),TYPE_DOUBLE, "Mass"          , 0 , NULL},
        (par_info){(void*) &(flpar->degeneracy      ),TYPE_INT   , "Degeneracy"    , 0 , NULL},
        (par_info){(void*) &(flpar->number_of_ps    ),TYPE_INT   , "PseudoFermions", 0 , NULL},
        (par_info){(void*) &(flpar->name            ),TYPE_STR   , "Name"          , 0 , NULL},
        (par_info){(void*) &(flpar->ferm_charge     ),TYPE_DOUBLE, "Charge"        , 0 , NULL},
        (par_info){(void*) &(flpar->ferm_im_chem_pot),TYPE_DOUBLE, "MuOverPiT"     , 0 , NULL}};

    // from here on, you should not have to modify anything.
    return scan_group_NV(sizeof(fp)/sizeof(par_info),fp, filelines, startline, endline);

}

int read_action_info(action_param *act_par,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{
    // see OpenAcc/su3_measurements.h
    const double stout_rho_def = RHO;

    par_info ap[]={
        (par_info){(void*) &(act_par->beta)       ,TYPE_DOUBLE,"Beta"      , 0 , NULL},
        (par_info){(void*) &(act_par->stout_steps),TYPE_INT   ,"StoutSteps", 0 , NULL},
        (par_info){(void*) &(act_par->stout_rho)  ,TYPE_DOUBLE,"StoutRho"  , 1 , (const void*) &stout_rho_def}};


    // from here on, you should not have to modify anything.
    int res = scan_group_NV(sizeof(ap)/sizeof(par_info),ap, filelines, startline, endline);

    if(startline<endline)
        if(act_par->stout_rho != RHO ){ 

            if(0==devinfo.myrank){
                printf("Error, input file stout_rho != RHO \n");
                printf("  Either modify the input file, or recompile changing RHO\n");
                printf(" (input) stout_rho = %f, (code) RHO = %f\n", act_par->stout_rho,RHO);
            }
            exit(1);

        }

    return res;



}
int read_backfield_info(bf_param *bfpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    // see /OpenAcc/backfield.h
    const double default_value = 0.0;
    par_info bfp[]={
        (par_info){(void*) &(bfpar->ex ),TYPE_DOUBLE, "ex", 1 ,(const void*) &default_value  },
        (par_info){(void*) &(bfpar->ey ),TYPE_DOUBLE, "ey", 1 ,(const void*) &default_value  },
        (par_info){(void*) &(bfpar->ez ),TYPE_DOUBLE, "ez", 1 ,(const void*) &default_value  },
        (par_info){(void*) &(bfpar->bx ),TYPE_DOUBLE, "bx", 1 ,(const void*) &default_value  },
        (par_info){(void*) &(bfpar->by ),TYPE_DOUBLE, "by", 1 ,(const void*) &default_value  },
        (par_info){(void*) &(bfpar->bz ),TYPE_DOUBLE, "bz", 1 ,(const void*) &default_value  }};

    // from here on, you should not have to modify anything.
    return scan_group_NV(sizeof(bfp)/sizeof(par_info),bfp, filelines, startline, endline);

}
int read_md_info(md_param *mdpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    // see /OpenAcc/md_integrator.h
    const double tlendef = 1.0;
    const double expmaxeigenv_def = 5.5 ; 
    const int max_cg_iterations_def = 10000;

    par_info mdp[]={
        (par_info){(void*) &(mdpar->no_md ),       TYPE_INT, "NmdSteps"     , 0 , NULL},
        (par_info){(void*) &(mdpar->gauge_scale ), TYPE_INT, "GaugeSubSteps", 0 , NULL},
        (par_info){(void*) &(mdpar->t ),        TYPE_DOUBLE, "TrajLength"   , 1 , (const void*) &tlendef},
        (par_info){(void*) &(mdpar->residue_metro),       TYPE_DOUBLE,   "residue_metro"          , 0, NULL},
        (par_info){(void*) &(mdpar->expected_max_eigenvalue),TYPE_DOUBLE,"ExpMaxEigenvalue"       , 1,(const void*) &expmaxeigenv_def},
        (par_info){(void*) &(mdpar->residue_md),TYPE_DOUBLE, "residue_md"   , 0 , NULL},
        (par_info){(void*) &(mdpar->max_cg_iterations),TYPE_INT, "MaxCGIterations"   , 1 , (const void*) &max_cg_iterations_def}};

    // from here on, you should not have to modify anything.
    return scan_group_NV(sizeof(mdp)/sizeof(par_info),mdp, filelines, startline, endline);

}
int read_mc_info(mc_params_t *mcpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    // see /OpenAcc/md_integrator.h
    const int seed_def = 0;  // which means time()
    const double epsgen_def = 0.1 ; 
    const char RandGenStatusFilename_def[] = "rgstatus.bin"; 
    const double MaxRunTimeS_def = 1.0e9; // 30 years should be enough
    const int MaxConfIdIter_def = 1000000; 

    par_info mcp[]={
        (par_info){(void*) &(mcpar->ntraj                  ),TYPE_INT,   "Ntraj"                  , 0, NULL},
        (par_info){(void*) &(mcpar->therm_ntraj            ),TYPE_INT,   "ThermNtraj"             , 0, NULL},
        (par_info){(void*) &(mcpar->storeconfinterval      ),TYPE_INT,   "StoreConfInterval"      , 0, NULL},
        (par_info){(void*) &(mcpar->saveconfinterval),       TYPE_INT,   "SaveConfInterval"       , 0, NULL},
        (par_info){(void*) &(mcpar->store_conf_name),        TYPE_STR,   "StoreConfName"          , 0, NULL},
        (par_info){(void*) &(mcpar->save_conf_name),         TYPE_STR,   "SaveConfName"           , 0, NULL},
        (par_info){(void*) &(mcpar->MaxConfIdIter),          TYPE_INT,   "MaxConfIdIter"          , 1,(const void*) &MaxConfIdIter_def},
        (par_info){(void*) &(mcpar->RandGenStatusFilename),  TYPE_STR,   "RandGenStatusFilename"  , 1,(const void*) &RandGenStatusFilename_def},
        (par_info){(void*) &(mcpar->MaxRunTimeS),         TYPE_DOUBLE,   "MaxRunTimeS"            , 1,(const void*) &MaxRunTimeS_def},
        (par_info){(void*) &(mcpar->seed),                   TYPE_INT,   "Seed"                   , 1,(const void*) &seed_def},
        (par_info){(void*) &(mcpar->eps_gen),             TYPE_DOUBLE,   "EpsGen"                 , 1,(const void*) &epsgen_def}};

    // from here on, you should not have to modify anything.

    return scan_group_NV(sizeof(mcp)/sizeof(par_info),mcp, filelines, startline, endline);

}

int read_debug_info(debug_settings_t * dbg_settings,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    const int useildg_def = 1;
    const int input_vbl_def = 2;
    const int SaveAllAtEnd_def = 1;
    const int save_diagnostics_def = 0;
    const char diagnostics_filename_def[] = "md_diagnostics.dat"; 
    const int do_reversibility_test_def = 0; 
    const int do_norandom_test_def = 0; 
    const int rng_fakeness_level_def = 0 ;
    const int md_dbg_print_max_count_def = 0;
    // see /Meas
    par_info gmp[]= {
        (par_info){(void*) &(dbg_settings->use_ildg),              TYPE_INT,"UseILDG"                , 1,(const void*) &useildg_def},
        (par_info){(void*) &(dbg_settings->input_vbl),             TYPE_INT,"VerbosityLv"            , 1,(const void*) &input_vbl_def},
        (par_info){(void*) &(dbg_settings->SaveAllAtEnd),          TYPE_INT,"SaveAllAtEnd"           , 1,(const void*) &SaveAllAtEnd_def},
        (par_info){(void*) &(dbg_settings->save_diagnostics),      TYPE_INT,"SaveDiagnostics"        , 1,(const void*) &save_diagnostics_def},
        (par_info){(void*) &(dbg_settings->do_reversibility_test), TYPE_INT,"DoRevTest"              , 1,(const void*) &do_reversibility_test_def},
        (par_info){(void*) &(dbg_settings->do_norandom_test),      TYPE_INT,"DoNoRandomTest"         , 1,(const void*) &do_norandom_test_def},
        (par_info){(void*) &(dbg_settings->rng_fakeness_level),    TYPE_INT,"RngFakenessLevel"       , 1,(const void*) &rng_fakeness_level_def},
        (par_info){(void*) &(dbg_settings->md_dbg_print_max_count),TYPE_INT,"MDDbgPrintMaxCount"     , 1,(const void*) &md_dbg_print_max_count_def},
        (par_info){(void*) &(dbg_settings->diagnostics_filename),  TYPE_STR,"SaveDiagnosticsFilename", 1,(const void*) &diagnostics_filename_def}};




    // from here on, you should not have to modify anything.
    return scan_group_NV(sizeof(gmp)/sizeof(par_info),gmp, filelines, startline, endline);

}



int read_gaugemeas_info(char *outfilename,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    // see /Meas
    par_info gmp[]= {
        (par_info){(void*) outfilename ,TYPE_STR, "GaugeOutfilename", 0 , NULL}};


    // from here on, you should not have to modify anything.
    return scan_group_NV(sizeof(gmp)/sizeof(par_info),gmp, filelines, startline, endline);

}
int read_fermmeas_info(ferm_meas_params * fmpars,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    const int doubleinv_def = 0;
    par_info fmp[]={
        (par_info){(void*) &(fmpars->fermionic_outfilename),TYPE_STR,"FermionicOutfilename", 0 , NULL},
        (par_info){(void*) &(fmpars->SingleInvNVectors),TYPE_INT,    "SingleInvNVectors"   , 0, NULL },
        (par_info){(void*) &(fmpars->DoubleInvNVectorsChiral),TYPE_INT,    "DoubleInvNVectorsChiral", 1, (const void*) &doubleinv_def  },
        (par_info){(void*) &(fmpars->DoubleInvNVectorsQuarkNumber),TYPE_INT,    "DoubleInvNVectorsQuarkNumber", 1, (const void*) &doubleinv_def   }};


    // from here on, you should not have to modify anything.
    return scan_group_NV(sizeof(fmp)/sizeof(par_info),fmp, filelines, startline, endline);

}
int read_device_setting(dev_info * di,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    // notice that pre_init_multidev1D or any relevant function
    // must have been called before, in order to get nranks
    // from MPI_Init()
    const int single_dev_choice_def = 0;
    const int async_comm_fermion_def = 0;
    const int async_comm_gauge_def   = 0;

    int helpmode = (int) (startline == endline);
#ifndef MULTIDEVICE        
    const int ignored_def = 1; // ignored, actually, but necessary
#endif

    par_info tp[]= {
#ifndef MULTIDEVICE        
        (par_info){(void*) &(di->single_dev_choice),TYPE_INT,"device_choice", 0 , NULL  },
        (par_info){(void*) IGNORE_IT,TYPE_INT,"AsyncFermionComms",1,(const void*) &async_comm_fermion_def},
        (par_info){(void*) IGNORE_IT,TYPE_INT,"AsyncGaugeComms"  ,1,(const void*) &async_comm_gauge_def  },
        (par_info){(void*) IGNORE_IT,TYPE_INT,"NProcPerNode",   1,(const void*) &ignored_def},
#else
        (par_info){(void*) &(di->single_dev_choice), TYPE_INT,"device_choice"    ,1,(const void*) &single_dev_choice_def},
        (par_info){(void*) &(di->async_comm_fermion),TYPE_INT,"AsyncFermionComms",1,(const void*) &async_comm_fermion_def},
        (par_info){(void*) &(di->async_comm_gauge),  TYPE_INT,"AsyncGaugeComms"  ,1,(const void*) &async_comm_gauge_def  },
        (par_info){(void*) &(di->proc_per_node),TYPE_INT,"NProcPerNode", 0 , NULL},
#endif
        (par_info){(void*) &(di->nranks_read),TYPE_INT,"NRanks", 0 , NULL}
    };

    // from here on, you should not have to modify anything.
    int res = scan_group_NV(sizeof(tp)/sizeof(par_info),tp, filelines, startline, endline);

    if(!helpmode){
#ifdef MULTIDEVICE
        if(di->nranks_read != di->nranks){
            printf("MPI%02d: ERROR: nranks from settings file ", di->myrank);
            printf("and from MPI_Init() DIFFER!\n");
            printf("settings: %d , MPI_Init(): %d\n",
                    di->nranks_read, di->nranks);
            res = 1;
        }
#else 
        if(di->nranks_read != 1){
            printf("ERROR: \'Nranks\' from setting file is %d,", di->nranks_read);
            printf(" but code is not compiled for muiltidevice\n");
            res = 1;
        }
#endif
    }
    return res;

}

int read_geometry(geom_parameters *gpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    // see /OpenAcc/backfield.h
    const int nx_def   = nd0; 
    const int ny_def   = nd1; 
    const int nz_def   = nd2; 
    const int nt_def   = nd3; 
    const int xmap_def = 0;
    const int ymap_def = 1;
    const int zmap_def = 2;
    const int tmap_def = 3;

    par_info gp[]={
        (par_info){(void*) &(gpar->gnx ),TYPE_INT,  "nx"  ,1,(const void*) &nx_def  },
        (par_info){(void*) &(gpar->gny ),TYPE_INT,  "ny"  ,1,(const void*) &ny_def  },
        (par_info){(void*) &(gpar->gnz ),TYPE_INT,  "nz"  ,1,(const void*) &nz_def  },
        (par_info){(void*) &(gpar->gnt ),TYPE_INT,  "nt"  ,1,(const void*) &nt_def  },
        (par_info){(void*) &(gpar->xmap ),TYPE_INT, "xmap",1,(const void*) &xmap_def},
        (par_info){(void*) &(gpar->ymap ),TYPE_INT, "ymap",1,(const void*) &ymap_def},
        (par_info){(void*) &(gpar->zmap ),TYPE_INT, "zmap",1,(const void*) &zmap_def},
        (par_info){(void*) &(gpar->tmap ),TYPE_INT, "tmap",1,(const void*) &tmap_def}};

    int res = scan_group_NV(sizeof(gp)/sizeof(par_info),gp, filelines, startline, endline);

    set_geom_glv(gpar);

    if(startline<endline){

        int expnx =(gpar->nd[gpar->xmap]-2*gpar->halos[gpar->xmap])*
            gpar->nranks[gpar->xmap]; 
        int expny =(gpar->nd[gpar->ymap]-2*gpar->halos[gpar->ymap])*
            gpar->nranks[gpar->ymap]; 
        int expnz =(gpar->nd[gpar->zmap]-2*gpar->halos[gpar->zmap])*
            gpar->nranks[gpar->zmap]; 
        int expnt =(gpar->nd[gpar->tmap]-2*gpar->halos[gpar->tmap])*
            gpar->nranks[gpar->tmap]; 

        if(gpar->gnx != expnx || gpar->gny != expny ||
                gpar->gnz != expnz  || gpar->gnt != expnt ){ 

            if(0==devinfo.myrank){
                printf("Error, input file lattice dimensions are not compatible\n");
                printf("       with the lattice dimensions written in geometry.h.\n");
                printf("       Either modify the input file, or recompile,\n");
                printf("(input) nx=%d\tny=%d\tnz=%d\tnt=%d\n",
                        gpar->gnx,gpar->gny,gpar->gnz,gpar->gnt);
                printf("(code)  nx=%dx%d\tny=%dx%d\tnz=%dx%d\tnt=%dx%d\n",
                        gpar->nd[gpar->xmap], gpar->nranks[gpar->xmap],
                        gpar->nd[gpar->ymap], gpar->nranks[gpar->ymap],
                        gpar->nd[gpar->zmap], gpar->nranks[gpar->zmap],
                        gpar->nd[gpar->tmap], gpar->nranks[gpar->tmap]);
            }
            res = 1;
        }
        int maps[4] = {gpar->xmap,gpar->ymap,gpar->zmap,gpar->tmap};
        int stop = 0;
        int imap,jmap;
        for(imap = 0 ; imap<3; imap++) for(jmap = imap+1 ; jmap<4; jmap++)
            stop = stop || (maps[imap] == maps[jmap]);

        if(stop){

            if(0==devinfo.myrank)
                printf("ERROR: found two equal direction mappings (%s:%d)\n",
                        __FILE__,__LINE__);
            res = 1;

        }
    }
    return res;

}


void set_global_vars_and_fermions_from_input_file(const char* input_filename)
{

    // Opening filenames and reading it
    int helpmode = 0;
    char filelines[MAXLINES][MAXLINELENGTH];
    FILE *input = fopen(input_filename,"r");
    if (input == NULL){

        if(0==devinfo.myrank){
            printf("Could not open file %s \n",input_filename );
            printf("writing an template_input file for your convenience.\n" );
        }
        helpmode = 1;
    }

    int lines_read = 0;

    int tagpositions[MAXPMG], tagtypes[MAXPMG],tagcounts[NPMGTYPES];
    int found_tags = 0;
    fermions_parameters = NULL;
    if (! helpmode){
        char *readcheck = filelines[0];
        while(readcheck != NULL){
            readcheck = fgets(filelines[lines_read],MAXLINELENGTH,input);
            if(readcheck != NULL) lines_read++;
        }
        fclose(input);

        if(0==devinfo.myrank)
            printf("lines read: %d\n", lines_read);

        int totlen = prepare_string_from_stringarray(filelines,lines_read,input_file_str);

        if(0==devinfo.myrank)
            printf("Written %d characters into input_file_str\n",totlen);


        // erasing comments 
        erase_comments(filelines,lines_read);
        //for(int dbgline = 0 ; dbgline < lines_read; dbgline++) // DBG
        //  printf("%d %s",dbgline, filelines[dbgline]);

        //scanning for macro parameter families

        found_tags = scan_group_V(NPMGTYPES,par_macro_groups_names,
                tagcounts,
                tagpositions,tagtypes,MAXPMG,
                filelines,0,lines_read );

        // see global var in /Include/fermion_parameters.
        // setting NDiffFlavs first
        NDiffFlavs = tagcounts[PMG_FERMION];
        if(NDiffFlavs==0)
        {
            fermions_parameters = NULL;
            if(0==devinfo.myrank){
                printf("NO FERMIONS FOUND, ");
                printf("SIMULATING PURE GAUGE THEORY...\n");
            }
        }
        else fermions_parameters = (ferm_param*) malloc(NDiffFlavs*sizeof(ferm_param));
    }
    else
    {   // goes into help mode
        found_tags = NPMGTYPES;
        for(int ifake_tag = 0; ifake_tag < found_tags; ifake_tag ++){
            tagpositions[ifake_tag] = 0; // so that all scan_group_NV() 
            // will go into 'help mode'
            tagcounts[ifake_tag] = 1;
            tagtypes[ifake_tag] = ifake_tag ; // so we have a tag for each 
            // type anyway

        }
        if(0==devinfo.myrank)
            helpfile = fopen("template_input", "w");
    }


    // check if all parameter groups were found
    int check = 1;
    if(!helpmode){
        for(int igrouptype  = 0 ; igrouptype < NPMGTYPES; igrouptype++)
            if(igrouptype != PMG_FERMION  && igrouptype != PMG_DEBUG )  check *= tagcounts[igrouptype];
        if(!check){
            for(int igrouptype  = 0 ; igrouptype < NPMGTYPES; igrouptype++)
                if (!tagcounts[igrouptype]) if(0==devinfo.myrank)
                    printf("\"%s\"  parameter group not found!\n",
                            par_macro_groups_names[igrouptype]);
            exit(1);
        }
    }


    int fermion_count = 0;
    // note 'check' is reused here
    int totcheck = 0;
    for(int igroup  = 0 ; igroup < found_tags; igroup++){
        int startline = tagpositions[igroup];
        int endline = (igroup<found_tags-1)?tagpositions[igroup+1]:lines_read;

        if(helpmode){
            if(0==devinfo.myrank){
                fprintf(helpfile,"\n\n%s\n",  par_macro_groups_names[tagtypes[igroup]]);
                printf("Writing %s...\n",  par_macro_groups_names[tagtypes[igroup]]);
            }
        }
        else if(0==devinfo.myrank){
            printf("\nReading %s, lines %d - %d ...\n", 
                    par_macro_groups_names[tagtypes[igroup]],
                    startline, endline);
        }

        switch(tagtypes[igroup]){
            case PMG_ACTION     :
                check = read_action_info(&act_params,filelines,startline,endline);
                break; 
            case PMG_FERMION   : 
                check = read_flavour_info(&(fermions_parameters[fermion_count]),
                        filelines,startline,endline);
                fermion_count++;
                break; 
            case PMG_BACKGROUND: 
                check = read_backfield_info(&backfield_parameters, 
                        filelines,startline,endline);
                break; 
            case PMG_MD        : 
                check = read_md_info(&md_parameters,filelines,startline,endline);
                break; 
            case PMG_MC        : 
                check = read_mc_info(&mc_params,filelines,startline,endline);
                break; 
            case PMG_GMEAS     : 
                check = read_gaugemeas_info(gauge_outfilename,
                        filelines,startline,endline);
                break; 
            case PMG_FMEAS     : 
                check = read_fermmeas_info(&fm_par,
                        filelines,startline,endline);
                break; 
            case PMG_DEVICE     : 
                check = read_device_setting(&devinfo,
                        filelines,startline,endline);
                break; 
            case PMG_GEOMETRY   : 
                check = read_geometry(&geom_par, 
                        filelines,startline,endline);
                break; 
            case PMG_DEBUG   : 
                check = read_debug_info(&debug_settings, 
                        filelines,startline,endline);
                break; 

        }
        if(check)
            if(0==devinfo.myrank)
                printf("Problem in group %s\n", par_macro_groups_names[tagtypes[igroup]]);
        totcheck += check;
    }
    if(tagcounts[PMG_DEBUG]==0)
        read_debug_info(&debug_settings,filelines,0,1);// Just to set default values

    // check == 1 means at least a parameter was not found.
    if(helpmode) exit(1);
    if(!helpmode){

        uint32_t hash = hash_settings();
        char hash_string[32];
        sprintf(hash_string,"%" PRIu32, hash);
        if(0==devinfo.myrank)
            printf("Hash of all relevant settings: %s\n", hash_string );
        strcat(gauge_outfilename, hash_string);
        strcat(fm_par.fermionic_outfilename, hash_string);
        strcat(debug_settings.diagnostics_filename,hash_string);

    }
    if(totcheck!=0){

        if(0==devinfo.myrank)
            printf("There are errors in some groups, exiting.\n")   ;
        exit(1);

    }


}




#endif


