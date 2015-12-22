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
#include "../RationalApprox/rationalapprox.h"
#include "../Meas/ferm_meas.h"
#include "../Meas/gauge_meas.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAXLINES 300
#define MAXLINELENGTH 500 // pretty long to accomodate all the comments

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

#define NPMGTYPES 8
#define MAXPMG 20
const char * par_macro_groups_name[] ={ // NPMTYPES strings, NO SPACES!
    "ActionParameters",         // 0 
    "FlavourParameters",        // 1   
    "BackgroundFieldParameters",// 2            
    "MDParameters",             // 3
    "MontecarloParameters",     // 4       
    "GaugeMeasuresSettings",    // 5
    "FermionMeasuresSettings",  // 6
    "DeviceSettings"            // 7
};
#define PMG_ACTION        0
#define PMG_FERMION       1
#define PMG_BACKGROUND    2
#define PMG_MD            3
#define PMG_MC            4
#define PMG_GMEAS         5
#define PMG_FMEAS         6
#define PMG_DEVICE        7
// last number should be NPMGTYPES - 1 !!


FILE * helpfile;

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

int check_parinfo_sequence(int npar, par_info * par_infos){

    for(int i = 0; i < npar; i++)
        for(int j = i+1; j < npar; j++)
            if (strstr(par_infos[j].name,par_infos[i].name)) return 0;
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
                if (strstr(par_infos[j].name,par_infos[i].name)){
                    printf("Reordering %s and %s \n",par_infos[j].name,par_infos[i].name );
                    par_info tmp = par_infos[i];
                    par_infos[i] = par_infos[j];
                    par_infos[j] = tmp;
                    allright = 0;
                    break;
                }           
            }
            if(allright) break;
        }
    }
}

int scan_group_V(int ntagstofind, const char **strtofind, 
        int *tagcount,
        int *taglines, int *tagtypes, int maxnres,
        char filelines[MAXLINES][MAXLINELENGTH], 
        int startline, int endline)
    // returns number of groups found, writes arrays
{   // scans for lines in the format
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
                printf("Found group %s on line %d\n",strtofind[itype], iline);
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

    void scan_group_NV(int npars,par_info* par_infos,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{   if(startline >= endline) // goes into 'help mode'
    for(int ipar = 0;ipar< npars ; ipar++)
        fprintf(
                helpfile,"%s\t%s\n",
                par_infos[ipar].name,
                type_strings[par_infos[ipar].type]);
    else // 'normal mode'
    {   // scans lines in the form 
        // NAME VALUE  +stuff which will be ignored
        // rc = readcheck 

        // necessary step if a parameter name is contained 
        // in another parameter name.
        reorder_par_infos(npars, par_infos);


        int * rc = (int *) malloc(npars*sizeof(int));
        for(int i =0; i<npars; i++) rc[i] = 0;

        int rc_all() {
            int res = 1 ; for(int i =0; i<npars; i++) res = res && rc[i];
            return res;
        }

        int iline = startline;

        while(! rc_all() && iline < endline)
        {
            char * found_something;
            for(int i =0; i<npars; i++){
                found_something = strstr(filelines[iline],par_infos[i].name);
                //                if(found_something == filelines[iline]){ // looks at the beginning of the line.
                if(found_something){ // looks at the beginning of the line.
                    // found parameter
                    printf("  %s\r\t\t\t ",par_infos[i].name);
                    int reads = 0;
                    char parname[20];
                    double tempdouble; int tempint ; char tempstr[20];
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
                                printf("%f\n", *((double*)par_infos[i].par));
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
            }


            if(! rc_all()){
                printf("ERROR: not all parameters needed read!");
                for(int i =0; i<npars; i++) 
                    if (rc[i]==0) printf("Parameter %s not set!\n",par_infos[i].name);
                exit(1);
            }

        }
    }

    void read_flavour_info(ferm_param *flpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
    {

        // see /Include/fermion_parameters.h
        const unsigned int  npar_fermions = 6;
        par_info fp[npar_fermions];

        // all names
        char sferm_mass[]        = "Mass"        ;
        char sdegeneracy[]       = "Degeneracy"       ;      
        char snumber_of_ps[]     = "PseudoFermions"     ;
        char sname[]             = "Name"             ;
        char sferm_charge[]      = "Charge"      ;
        char sferm_im_chem_pot[] = "Mu" ;

        fp[0]=(par_info){(void*) &(flpar->ferm_mass       ),TYPE_DOUBLE, sferm_mass       };
        fp[1]=(par_info){(void*) &(flpar->degeneracy      ),TYPE_INT   , sdegeneracy      };
        fp[2]=(par_info){(void*) &(flpar->number_of_ps    ),TYPE_INT   , snumber_of_ps    };
        fp[3]=(par_info){(void*) &(flpar->name            ),TYPE_STR   , sname            };
        fp[4]=(par_info){(void*) &(flpar->ferm_charge     ),TYPE_DOUBLE, sferm_charge     };
        fp[5]=(par_info){(void*) &(flpar->ferm_im_chem_pot),TYPE_DOUBLE, sferm_im_chem_pot};


        // from then on, you should not have to modify anything.
        scan_group_NV(npar_fermions,fp, filelines, startline, endline);

#ifndef IMCHEMPOT
        if(flpar)
            if(flpar->ferm_im_chem_pot != 0) {
                printf("ERROR! Found a nonzero chemical potential for quark %s!\n", flpar->name);
                printf("ERROR! This program was not compiled with \'IMCHEMPOT\'.\n");
                printf("ERROR! This program cannot perform simulations at nonzero \n");
                printf("       chemical potential. \n");
                exit(1);
            }
#endif
    }

    void read_action_info(action_param *act_par,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
    {
        // see OpenAcc/su3_measurements.h
        const unsigned int npar_action = 2;
        par_info ap[npar_action];

        char sbeta[]        = "Beta" ;
        char sstoutsteps[]  = "StoutSteps"  ;
        ap[0]=(par_info){(void*) &(act_par->beta),TYPE_DOUBLE, sbeta };
        ap[1]=(par_info){(void*) &(act_par->stout_steps),TYPE_INT, sstoutsteps };

        // from then on, you should not have to modify anything.
        scan_group_NV(npar_action,ap, filelines, startline, endline);

    }
    void read_backfield_info(bf_param *bfpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
    {

        // see /OpenAcc/backfield.h
        const unsigned int npar_background = 6;
        par_info bfp[npar_background];

        char sex[] = "ex" ;
        char sey[] = "ey" ;
        char sez[] = "ez" ;
        char sbx[] = "bx" ;
        char sby[] = "by" ;
        char sbz[] = "bz" ;
        bfp[0]=(par_info){(void*) &(bfpar->ex ),TYPE_DOUBLE, sex };
        bfp[1]=(par_info){(void*) &(bfpar->ey ),TYPE_DOUBLE, sey };
        bfp[2]=(par_info){(void*) &(bfpar->ez ),TYPE_DOUBLE, sez };
        bfp[3]=(par_info){(void*) &(bfpar->bx ),TYPE_DOUBLE, sbx };
        bfp[4]=(par_info){(void*) &(bfpar->by ),TYPE_DOUBLE, sby };
        bfp[5]=(par_info){(void*) &(bfpar->bz ),TYPE_DOUBLE, sbz };


        // from then on, you should not have to modify anything.
        scan_group_NV(npar_background,bfp, filelines, startline, endline);

#ifndef BACKFIELD
        if(bfpar->ex || bfpar->ey || bfpar->ez || bfpar->bx || bfpar->by || bfpar->bz ){
            printf("ERROR! Found nonzero field value, please check.\n");
            printf("ERROR! \'BACKFIELD\' not defined at compile time.\n");
            printf("ERROR! This program DOES NOT run simulation with a\n");
            printf("       background field.\n");
            exit(1);
        }
#endif
    }
    void read_md_info(md_param *mdpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
    {

        // see /OpenAcc/md_integrator.h
        const unsigned int  npar_md =  3;
        par_info mdp[npar_md];

        char snomd[] = "NmdSteps" ;
        char sgs[] = "GaugeSubSteps" ;
        char st[] = "TrajLength" ;

        mdp[0]=(par_info){(void*) &(mdpar->no_md ),TYPE_INT, snomd };
        mdp[1]=(par_info){(void*) &(mdpar->gauge_scale ),TYPE_INT, sgs };
        mdp[2]=(par_info){(void*) &(mdpar->t ),TYPE_DOUBLE, st};


        // from then on, you should not have to modify anything.
        scan_group_NV(npar_md,mdp, filelines, startline, endline);

    }
    void read_mc_info(mc_param *mcpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
    {

        // see /OpenAcc/md_integrator.h
        const unsigned int  npar_mc =  10;
        par_info mcp[npar_mc];

        char sntraj[] = "Ntraj" ;
        char stherm_ntraj[] = "ThermNtraj" ;
        char ssaveconfinterval[] = "SaveConfInterval" ;
        char ssaverunningconfinterval[] = "SaveRunningConfInterval";
        char sresidue_metro[] = "residue_metro";
        char sresidue_md[] = "residue_md";
        char sstore_conf_name[] = "StoreConfName";
        char ssave_conf_name[] = "SaveConfName";
        char sseed[] = "Seed";
        char sinput_vbl[] = "VerbosityLv";

        mcp[0]=(par_info){(void*) &(mcpar->ntraj                  ),TYPE_INT,sntraj          };
        mcp[1]=(par_info){(void*) &(mcpar->therm_ntraj            ),TYPE_INT,stherm_ntraj    };
        mcp[2]=(par_info){(void*) &(mcpar->saveconfinterval       ),TYPE_INT, ssaveconfinterval};
        mcp[3]=(par_info){(void*) &(mcpar->saverunningconfinterval),TYPE_INT, ssaverunningconfinterval};
        mcp[4]=(par_info){(void*) &(mcpar->residue_metro    ),TYPE_DOUBLE,sresidue_metro};
        mcp[5]=(par_info){(void*) &(mcpar->residue_md       ),TYPE_DOUBLE,sresidue_md};
        mcp[6]=(par_info){(void*) &(mcpar->store_conf_name  ),TYPE_STR,sstore_conf_name};
        mcp[7]=(par_info){(void*) &(mcpar->save_conf_name   ),TYPE_STR,ssave_conf_name};
        mcp[8]=(par_info){(void*) &(mcpar->seed   ),TYPE_INT,sseed};

        mcp[9]=(par_info){(void*) &(mcpar->input_vbl   ),TYPE_INT,sinput_vbl};

        // from then on, you should not have to modify anything.

        scan_group_NV(npar_mc,mcp, filelines, startline, endline);

    }
    void read_gaugemeas_info(char *outfilename,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
    {

        // see /Meas
        const unsigned int npar_gaugemeas = 1 ; 
        par_info gmp[npar_gaugemeas];

        char soutfilename[] = "GaugeOutfilename" ;

        gmp[0]=(par_info){(void*) outfilename ,TYPE_STR,soutfilename };


        // from then on, you should not have to modify anything.
        scan_group_NV(npar_gaugemeas,gmp, filelines, startline, endline);

    }
    void read_fermmeas_info(char * outfilename,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
    {

        const unsigned int npar_fermmeas = 1 ; 
        par_info fmp[npar_fermmeas];

        char soutfilename[] = "FermionicOutfilename" ;

        fmp[0]=(par_info){(void*) outfilename ,TYPE_STR,soutfilename };


        // from then on, you should not have to modify anything.
        scan_group_NV(npar_fermmeas,fmp, filelines, startline, endline);

    }
    void read_device_setting(device_param *device_settings,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
    {

        const unsigned int npar_device_settings = 1 ; 
        par_info tp[npar_device_settings];

        char sdevice_choice[] = "device_choice" ;

        tp[0]=(par_info){(void*) &(device_settings->device_choice),TYPE_INT,sdevice_choice };


        // from then on, you should not have to modify anything.
        scan_group_NV(npar_device_settings,tp, filelines, startline, endline);

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


            // erasing comments 
            erase_comments(filelines,lines_read);
            //for(int dbgline = 0 ; dbgline < lines_read; dbgline++) // DBG
            //  printf("%d %s",dbgline, filelines[dbgline]);




            //scanning for macro parameter families

            found_tags = scan_group_V(NPMGTYPES,par_macro_groups_name,
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
                if (!tagcounts[igrouptype]) printf("\"%s\" not found!\n",par_macro_groups_name[igrouptype]);
            exit(1);
        }


        int fermion_count = 0;
        for(int igroup  = 0 ; igroup < found_tags; igroup++){
            int startline = tagpositions[igroup];
            int endline = (igroup<found_tags-1)?tagpositions[igroup+1]:lines_read;

            if(helpmode) fprintf(helpfile,"\n\n%s\n",  par_macro_groups_name[tagtypes[igroup]]);
            else printf("Reading %s...\n", par_macro_groups_name[tagtypes[igroup]]);
            switch(tagtypes[igroup]){
                case PMG_ACTION     :
                    read_action_info(&act_params,filelines,startline,endline);
                    break; 
                case PMG_FERMION   : 
                    read_flavour_info(&(fermions_parameters[fermion_count]),
                            filelines,startline,endline);
                    fermion_count++;
                    break; 
                case PMG_BACKGROUND: 
                    read_backfield_info(&backfield_parameters, 
                            filelines,startline,endline);
                    break; 
                case PMG_MD        : 
                    read_md_info(&md_parameters,filelines,startline,endline);
                    break; 
                case PMG_MC        : 
                    read_mc_info(&mkwch_pars,filelines,startline,endline);
                    break; 
                case PMG_GMEAS     : 
                    read_gaugemeas_info(gauge_outfilename,
                            filelines,startline,endline);
                    break; 
                case PMG_FMEAS     : 
                    read_fermmeas_info(fermionic_outfilename,
                            filelines,startline,endline);
                    break; 
                case PMG_DEVICE     : 
                    read_device_setting(&dev_settings,
                            filelines,startline,endline);
                    break; 

            }
        }

        if(helpmode) exit(1);
    }




#endif


