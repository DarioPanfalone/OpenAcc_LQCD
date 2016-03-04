#ifndef _INIT_C_
#define _INIT_C_

#include "./init.h"
#include "./common_defines.h"
#include "./markowchain.h"
#include "./fermion_parameters.h"
#include "../OpenAcc/action.h"
#include "../OpenAcc/md_integrator.h"
#include "../OpenAcc/backfield.h"
#include "../OpenAcc/su3_measurements.h"
#include "../OpenAcc/deviceinit.h"
#include "../OpenAcc/geometry.h"
#include "../RationalApprox/rationalapprox.h"
#include "../Meas/ferm_meas.h"
#include "../Meas/gauge_meas.h"

#include <stdio.h>
#include <strings.h>
#include <stdlib.h>

#define MAXLINES 300
#define MAXLINELENGTH 500 // pretty long to accomodate all the comments

char input_file_str[MAXLINES*MAXLINELENGTH];

//types that can be found in an input file
#define TYPE_INT 0
#define TYPE_DOUBLE 1
#define TYPE_STR 2
const char * type_strings[]={"(int)", "(double)", "(string)"};

typedef struct par_info_t{

    void* par;
    int type;
    char* name;

}par_info;

#define NPMGTYPES 9
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
    "Geometry"                  // 8
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
// last number should be NPMGTYPES - 1 !!


FILE * helpfile;

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
                //printf("DEBUG (%d,%d)params %s and %s, match %p \n",i,j,par_infos[j].name,par_infos[i].name, match );
                if (match){
                    printf("Reordering %s and %s\n",par_infos[j].name,par_infos[i].name);
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
        printf("Found %s %d times.\n", strtofind[i], tagcount[i]); 

    return nres;
}

int scan_group_NV(int npars,par_info* par_infos,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{   
    if(startline >= endline){ // goes into 'help mode'
        for(int ipar = 0;ipar< npars ; ipar++)
            fprintf(
                    helpfile,"%s\t%s\n",
                    par_infos[ipar].name,
                    type_strings[par_infos[ipar].type]);
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

        int iline = startline;

        while(! res && iline < endline)
        {
            char * found_something;
            for(int i =0; i<npars; i++){
                found_something = strstr(filelines[iline],par_infos[i].name);
                if(found_something){ // looks at the beginning of the line.
                    // found parameter
                    printf("  %s\r\t\t\t\t ",par_infos[i].name);
                    int reads = 0;
                    char parname[50];
                    switch(par_infos[i].type){
                        case TYPE_INT: 
                            reads = sscanf(filelines[iline],
                                    "%s %d",parname,(int*)par_infos[i].par);
                            if(reads == 2) 
                                printf("%d\n", *((int*)par_infos[i].par));
                            break;
                        case TYPE_DOUBLE:

                            reads = sscanf(filelines[iline],
                                    "%s %lf",parname,(double*) par_infos[i].par);
                            if(reads == 2)
                                printf("%e\n", *((double*)par_infos[i].par));
                            break;
                        case TYPE_STR: 
                            reads = sscanf(filelines[iline],
                                    "%s %s",parname,(char*) par_infos[i].par);
                            if(reads == 2) 
                                printf("\"%s\"\n", ((char*)par_infos[i].par));
                            break;
                        default: 
                            printf("WARNING, variable type not set in sourcecode.\n");
                            break;

                    }
                    if(reads == 2)rc[i]++;
                    else printf("WARNING, NO VALUE READ!");
                    break;
                };
            }
            iline++; 
            res = 1;
            for(int i =0; i<npars; i++)
                res = res && rc[i];
        }



        if(! res){
            printf("ERROR: not all parameters needed read!");
            for(int i =0; i<npars; i++) 
                if (rc[i]==0) printf("Parameter %s not set!\n",par_infos[i].name);
            free(rc);
            return 1;
        }
        else{
            free(rc);
            return 0;
        }
    }
}

int read_flavour_info(ferm_param *flpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    // see /Include/fermion_parameters.h

    par_info fp[]={
    (par_info){(void*) &(flpar->ferm_mass       ),TYPE_DOUBLE, "Mass"          },
    (par_info){(void*) &(flpar->degeneracy      ),TYPE_INT   , "Degeneracy"    },
    (par_info){(void*) &(flpar->number_of_ps    ),TYPE_INT   , "PseudoFermions"},
    (par_info){(void*) &(flpar->name            ),TYPE_STR   , "Name"          },
    (par_info){(void*) &(flpar->ferm_charge     ),TYPE_DOUBLE, "Charge"        },
    (par_info){(void*) &(flpar->ferm_im_chem_pot),TYPE_DOUBLE, "MuOverPiT"     }};

    // from here on, you should not have to modify anything.
    return scan_group_NV(sizeof(fp)/sizeof(par_info),fp, filelines, startline, endline);
    
}

int read_action_info(action_param *act_par,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{
    // see OpenAcc/su3_measurements.h

    par_info ap[]={
    (par_info){(void*) &(act_par->beta)       ,TYPE_DOUBLE,"Beta"      },
    (par_info){(void*) &(act_par->stout_steps),TYPE_INT   ,"StoutSteps"},
    (par_info){(void*) &(act_par->stout_rho)  ,TYPE_DOUBLE,"StoutRho"  }};


    // from here on, you should not have to modify anything.
    int res = scan_group_NV(sizeof(ap)/sizeof(par_info),ap, filelines, startline, endline);

    if(startline<endline)
    if(act_par->stout_rho != RHO ){ 
        printf("Error, input file stout_rho != RHO \n");
        printf("  Either modify the input file, or recompile changing RHO\n");
        printf(" (input) stout_rho = %f, (code) RHO = %f\n", act_par->stout_rho,RHO);
        exit(1);

    }

    return res;



}
int read_backfield_info(bf_param *bfpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    // see /OpenAcc/backfield.h
    par_info bfp[]={
    (par_info){(void*) &(bfpar->ex ),TYPE_DOUBLE, "ex" },
    (par_info){(void*) &(bfpar->ey ),TYPE_DOUBLE, "ey" },
    (par_info){(void*) &(bfpar->ez ),TYPE_DOUBLE, "ez" },
    (par_info){(void*) &(bfpar->bx ),TYPE_DOUBLE, "bx" },
    (par_info){(void*) &(bfpar->by ),TYPE_DOUBLE, "by" },
    (par_info){(void*) &(bfpar->bz ),TYPE_DOUBLE, "bz" }};

    // from here on, you should not have to modify anything.
    return scan_group_NV(sizeof(bfp)/sizeof(par_info),bfp, filelines, startline, endline);

}
int read_md_info(md_param *mdpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    // see /OpenAcc/md_integrator.h

    par_info mdp[]={
    (par_info){(void*) &(mdpar->no_md ),       TYPE_INT, "NmdSteps"     },
    (par_info){(void*) &(mdpar->gauge_scale ), TYPE_INT, "GaugeSubSteps"},
    (par_info){(void*) &(mdpar->t ),        TYPE_DOUBLE, "TrajLength"   },
    (par_info){(void*) &(mdpar->residue_md),TYPE_DOUBLE, "residue_md"   }};


    // from here on, you should not have to modify anything.
    return scan_group_NV(sizeof(mdp)/sizeof(par_info),mdp, filelines, startline, endline);

}
int read_mc_info(mc_param *mcpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    // see /OpenAcc/md_integrator.h
    par_info mcp[]={
    (par_info){(void*) &(mcpar->ntraj                  ),TYPE_INT,   "Ntraj"                  },
    (par_info){(void*) &(mcpar->therm_ntraj            ),TYPE_INT,   "ThermNtraj"             },
    (par_info){(void*) &(mcpar->storeconfinterval      ),TYPE_INT,   "StoreConfInterval"      },
    (par_info){(void*) &(mcpar->saveconfinterval),TYPE_INT,          "SaveConfInterval"       },
    (par_info){(void*) &(mcpar->residue_metro    ),TYPE_DOUBLE,      "residue_metro"          },
    (par_info){(void*) &(mcpar->store_conf_name  ),TYPE_STR,         "StoreConfName"          },
    (par_info){(void*) &(mcpar->save_conf_name   ),TYPE_STR,         "SaveConfName"           },
    (par_info){(void*) &(mcpar->use_ildg),TYPE_INT,                  "UseILDG"                },
    (par_info){(void*) &(mcpar->seed   ),TYPE_INT,                   "Seed"                   },
    (par_info){(void*) &(mcpar->eps_gen  ),TYPE_DOUBLE,              "EpsGen"                 },
    (par_info){(void*) &(mcpar->input_vbl  ),TYPE_INT,               "VerbosityLv"            },
    (par_info){(void*) &(mcpar->expected_max_eigenvalue),TYPE_DOUBLE,"ExpMaxEigenvalue"       },
    (par_info){(void*) &(mcpar->save_diagnostics),TYPE_INT,          "SaveDiagnostics"        },
    (par_info){(void*) &(mcpar->diagnostics_filename),TYPE_STR,      "SaveDiagnosticsFilename"}};

    // from here on, you should not have to modify anything.

    return scan_group_NV(sizeof(mcp)/sizeof(par_info),mcp, filelines, startline, endline);

}
int read_gaugemeas_info(char *outfilename,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    // see /Meas
    par_info gmp[]= {
    (par_info){(void*) outfilename ,TYPE_STR, "GaugeOutfilename"}};


    // from here on, you should not have to modify anything.
   return scan_group_NV(sizeof(gmp)/sizeof(par_info),gmp, filelines, startline, endline);

}
int read_fermmeas_info(ferm_meas_params * fmpars,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    par_info fmp[]={
    (par_info){(void*) &(fmpars->fermionic_outfilename),TYPE_STR,"FermionicOutfilename"},
    (par_info){(void*) &(fmpars->SingleInvNVectors),TYPE_INT,    "SingleInvNVectors"   },
    (par_info){(void*) &(fmpars->DoubleInvNVectors),TYPE_INT,    "DoubleInvNVectors"   }};

    // from here on, you should not have to modify anything.
    return scan_group_NV(sizeof(fmp)/sizeof(par_info),fmp, filelines, startline, endline);

}
int read_device_setting(device_param *device_settings,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    par_info tp[]= {
    (par_info){(void*) &(device_settings->device_choice),TYPE_INT,"device_choice" }};


    // from here on, you should not have to modify anything.
    return scan_group_NV(sizeof(tp)/sizeof(par_info),tp, filelines, startline, endline);

}

int read_geometry(geom_parameters *gpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    // see /OpenAcc/backfield.h
    par_info gp[]={
    (par_info){(void*) &(gpar->gnx ),TYPE_INT,  "nx" },
    (par_info){(void*) &(gpar->gny ),TYPE_INT,  "ny" },
    (par_info){(void*) &(gpar->gnz ),TYPE_INT,  "nz" },
    (par_info){(void*) &(gpar->gnt ),TYPE_INT,  "nt" },
    (par_info){(void*) &(gpar->xmap ),TYPE_INT, "xmap" },
    (par_info){(void*) &(gpar->ymap ),TYPE_INT, "ymap" },
    (par_info){(void*) &(gpar->zmap ),TYPE_INT, "zmap" },
    (par_info){(void*) &(gpar->tmap ),TYPE_INT, "tmap" }};

    int res = scan_group_NV(sizeof(gp)/sizeof(par_info),gp, filelines, startline, endline);

    set_geom_glv(gpar);

    if(startline<endline){

        int expnx = gpar->nd[gpar->xmap] * gpar->nranks[gpar->xmap]; 
        int expny = gpar->nd[gpar->ymap] * gpar->nranks[gpar->ymap]; 
        int expnz = gpar->nd[gpar->zmap] * gpar->nranks[gpar->zmap]; 
        int expnt = gpar->nd[gpar->tmap] * gpar->nranks[gpar->tmap]; 

    if(gpar->gnx != expnx || gpar->gny != expny ||
            gpar->gnz != expnz  || gpar->gnt != expnt ){ 
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
        exit(1);
    }
    int maps[4] = {gpar->xmap,gpar->ymap,gpar->zmap,gpar->tmap};
    int stop = 0;
    int imap,jmap;
    for(imap = 0 ; imap<3; imap++) for(jmap = imap+1 ; jmap<4; jmap++)
        stop = stop || (maps[imap] == maps[jmap]);

    if(stop){
        printf("ERROR: found two equal direction mappings (%s:%d)\n",
                __FILE__,__LINE__);
        exit(1);
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
    if (input == NULL) {
        printf("Could not open file %s \n",input_filename );
        printf("writing an input_template for your convenience.\n" );
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

        printf("lines read: %d\n", lines_read);

        int totlen = prepare_string_from_stringarray(filelines,lines_read,input_file_str);
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
        fermions_parameters = (ferm_param*) malloc(NDiffFlavs*sizeof(ferm_param));
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
        helpfile = fopen("input_template", "w");
    }


    // check if all parameter groups were found
    int check = 1;
    for(int igrouptype  = 0 ; igrouptype < NPMGTYPES; igrouptype++)
        check *= tagcounts[igrouptype];
    if(!check){
        for(int igrouptype  = 0 ; igrouptype < NPMGTYPES; igrouptype++)
            if (!tagcounts[igrouptype]) printf("\"%s\"  parameter group not found!\n",
                    par_macro_groups_names[igrouptype]);
        exit(1);
    }


    int fermion_count = 0;
    // note 'check' is reused here
    int totcheck = 0;
    for(int igroup  = 0 ; igroup < found_tags; igroup++){
        int startline = tagpositions[igroup];
        int endline = (igroup<found_tags-1)?tagpositions[igroup+1]:lines_read;

        if(helpmode) fprintf(helpfile,"\n\n%s\n",  par_macro_groups_names[tagtypes[igroup]]);
        else printf("Reading %s...\n", par_macro_groups_names[tagtypes[igroup]]);
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
                check = read_mc_info(&mkwch_pars,filelines,startline,endline);
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
                check = read_device_setting(&dev_settings,
                        filelines,startline,endline);
                break; 
            case PMG_GEOMETRY   : 
                check = read_geometry(&geom_par, 
                        filelines,startline,endline);
                break; 

        }
        if(check)
            printf("Problem in group %s\n", par_macro_groups_names[tagtypes[igroup]]);
        totcheck += check;
    }

    // check == 1 means at least a parameter was not found.
    if(helpmode) exit(1);
    if(totcheck!=0){
    
        printf("There are errors in some groups, exiting.\n")   ;
        exit(1);

    }

    
}




#endif


