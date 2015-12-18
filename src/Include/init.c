#ifndef _INIT_C_
#define _INIT_C_

#include "./init.h"
#include "./fermion_parameters.h"
#include "../OpenAcc/md_integrator.h"
#include "../RationalApprox/rationalapprox.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAXLINES 300
#define MAXLINELENGTH 50

//types that can be found in an input file
#define TYPE_INT 0
#define TYPE_DOUBLE 1
#define TYPE_STR 2
typedef struct par_info_t{

    void* par;
    int type;
    char* name;

}par_info;

#define NPMGTYPES 7
#define MAXPMG 20
char par_macro_groups_name[][] ={
    "Gauge Parameters",           // 0 
    "Flavour Parameters",         // 1   
    "Background Field Parameters",// 2            
    "MD parameters",              // 3
    "Montecarlo Parameters",      // 4       
    "Gauge Measures Settings",    // 5
    "Fermion Measures Settings"   // 6
}
#define PMG_GAUGE      0
#define PMG_FERMION    1
#define PMG_BACKGROUND 2
#define PMG_MD         3
#define PMG_MC         4
#define PMG_GMEAS      5
#define PMG_FMEAS      6



void erase_comments(char **file_lines, int maxlines)
{
    for(int i = 0; i<maxlines; i++){
       char* start_comment_ptr = index(file_lines[i],'#' );
       if(start_comment_ptr){
           int comment_start = start_comment_ptr - file_lines[i];
           file_lines[i][comment_start] = '\0';
           for(int j = comment_start+1; j < MAXLINELENGTH; j++)
               file_lines[i][j] = ' ';
       }
    }
}

int find_parameter_macro_groups(int *parameter_lines, int *types, char **filelines, int maxlines)// returns number of groups found, writes arrays
{
    // the arrays 'parameter_lines' and 'types' and must be MAXPMG long
    printf("Scanning for parameter macro groups..\n");
    
    int * rc = (int *) malloc(NPMG*sizeof(int));//read check
    for(int i =0; i<NPMGTYPES; i++) rc[i] = 0;
    
    for(int i =0; i<MAXPMG; i++){// initializing output arrays
        parameter_lines[i] = -1;
        types[i] = -1;
    }
    int impg = 0;
    for(int iline = 0; iline<maxlines; iline++){
        char *found_something;
        for(int itype =0; itype <NPMGTYPES; itype++){
            found_something = 
                strcasestr(filelines[iline],par_macro_groups_name[itype]);
            if(found_something){
                printf("Found group %s \n",par_macro_groups_name[itype]);
                parameter_lines[impg] = iline;
                types[impg] = itype;
                impg++;
                if(rc[itype]!=0 && itype != PMG_FERMION ){
                    printf("Error, found %s more than once in input file!\n",
                            par_macro_groups_name[itype]);
                    printf("See line %d\n", iline);
               //     exit(1);
                }
                else rc[itype]++;
            } 
        }
    }
    for(int i =0; i<NPMGTYPES; i++)
        printf("Found %s %d times.\n", par_macro_groups_name[i], rc[i]); 

    return impg;
}

void scan_group_NV(int npars,par_info* par_infos, char **filelines, int startline, int endline)
{
// scans lines in the form 
// NAME VALUE  +stuff which will be ignored
    // rc = readcheck 
    int * rc = (int *) malloc(npars*sizeof(int));
    for(int i =0; i<npars; i++) rc[i] = 0;

    int rc_all() {
        int res = 1 ; for(int i =0; i<npars; i++) res = res && rc[i];
        return res;
    }

    int iline = startline;

    while(! rc_all() && iline < endline){
        char *found_something;
        for(int i =0; i<npars; i++){
            found_something = strcasestr(filelines[iline],par_infos[i].name);
            if(found_something){
               // found parameter
               printf("Found %s = ",par_infos[i].name);
               int reads = 0;
               char parname[20];
               if(partypes[i] == TYPE_INT){ 
                   reads = sscanf(filelines[iline],
                           "%s %d",parname,(int*)par_infos[i].par);
                   if(reads == 1) printf("%d\n", *((int*)par_infos[i].par));
               }
               else if(partypes[i] == TYPE_DOUBLE){
                   reads = sscanf(filelines[iline],
                           "%s %f",parname,(double*) par_infos[i].par);
                   if(reads == 1) printf("%f\n", *((double*)par_infos[i].par));
               }
               else if(partypes[i] == TYPE_STR){
                   reads = sscanf(filelines[iline],
                           "%s %s",parname,(char*) par_infos[i].par);
                   if(reads == 1) printf("%s\n", *((char*)par_infos[i].par));
               }
               else printf("WARNING, variable type not set in sourcecode.\n");

               if(reads == 1) rc[i] = 1;
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



void read_flavour_info(ferm_param *flpar, char** filelines, int startline, int endline)
{

  // see src/Include/fermion_parameters.h
  par_info fp[6];

  char sferm_mass[]        = "ferm_mass"        ;
  char sdegeneracy[]       = "degeneracy"       ;      
  char snumber_of_ps[]     = "number_of_ps"     ;
  char sname[]             = "name"             ;
  char sferm_charge[]      = "ferm_charge"      ;
  char sferm_im_chem_pot[] = "ferm_im_chem_pot" ;

  fp[0]={(void*) &(flpar->ferm_mass       ),TYPE_DOUBLE, sferm_mass       };
  fp[1]={(void*) &(flpar->degeneracy      ),TYPE_INT   , sdegeneracy      };
  fp[2]={(void*) &(flpar->number_of_ps    ),TYPE_INT   , snumber_of_ps    };
  fp[3]={(void*) &(flpar->name            ),TYPE_STR   , sname            };
  fp[4]={(void*) &(flpar->ferm_charge     ),TYPE_DOUBLE, sferm_charge     };
  fp[5]={(void*) &(flpar->ferm_im_chem_pot),TYPE_DOUBLE, sferm_im_chem_pot};

  scan_group_NV(6,fp, filelines, startline, endline);


}
void read_gauge_info(gauge_param *gpar, char** filelines, int startline, int endline)
{

  par_info gp[1];

  char sbeta[]        = "beta"        ;
  gp[0]={(void*) &(gpar->beta       ),TYPE_DOUBLE, sbeta };
  scan_group_NV(1,gp, filelines, startline, endline);

}
void read_backfield_info(bf_param *bfpar, char** filelines, int startline, int endline)
{

  // see src/OpenAcc/backfield.h
  par_info bfp[6];

  char sex[] = "ex" ;
  char sey[] = "ey" ;
  char sez[] = "ez" ;
  char sbx[] = "bx" ;
  char sby[] = "by" ;
  char sbz[] = "bz" ;
  bfp[0]={(void*) &(bf_param->ex ),TYPE_DOUBLE, sex };
  bfp[1]={(void*) &(bf_param->ey ),TYPE_DOUBLE, sey };
  bfp[2]={(void*) &(bf_param->ez ),TYPE_DOUBLE, sez };
  bfp[3]={(void*) &(bf_param->bx ),TYPE_DOUBLE, sbx };
  bfp[4]={(void*) &(bf_param->by ),TYPE_DOUBLE, sby };
  bfp[5]={(void*) &(bf_param->bz ),TYPE_DOUBLE, sbz };

  scan_group_NV(6,bfp, filelines, startline, endline);

}
void read_md_info(md_param *mdpar, char** filelines, int startline, int endline)
{

  // see src/OpenAcc/md_integrator.h
  par_info mdp[3];

  char snomd[] = "NmdSteps" ;
  char sgs[] = "GaugeSteps" ;
  char st[] = "t" ;

  mdp[0]={(void*) &(mdpar->nomd ),TYPE_INT, snomd };
  mdp[1]={(void*) &(mdpar->gs ),TYPE_INT, sgs };
  mdp[2]={(void*) &(mdpar->t ),TYPE_DOUBLE, st};

  scan_group_NV(3,mdp, filelines, startline, endline);

}
void read_mc_info(mc_param *mcpar, char** filelines, int startline, int endline)
{

  // see src/OpenAcc/md_integrator.h
  par_info mcp[2];

  char sntraj[] = "Ntraj" ;
  char ssaveconfinterval[] = "SaveConfInterval" ;

  mcp[0]={(void*) &(mcpar->ntraj           ),TYPE_INT,sntraj           };
  mcp[1]={(void*) &(mcpar->saveconfinterval),TYPE_INT,ssaveconfinterval};

  scan_group_NV(2,mcp, filelines, startline, endline);

}
void read_gaugemeas_info(gm_param *gmpar, char** filelines, int startline, int endline)
{

  // see src/OpenAcc/md_integrator.h
  par_info gmp[1];

  char soutfilename[] = "Gauge Outfilename" ;

  mcp[0]={(void*) &(gmpar->outfilename ),TYPE_STR,soutfilename };

  scan_group_NV(1,gmp, filelines, startline, endline);

}






void set_global_vars_and_fermions_from_input_file(const char* input_filename){

    FILE *input = fopen(input_filename);
    if (input == NULL) {
        printf("Could not open file %s \n",nomefile );
        exit(1);
    }

    char filelines[MAXLINES][MAXLINELENGTH];
    char *readcheck = filelines[0]; int i = 0
    while(readcheck != NULL){
        readcheck = fgets(&filelines[i],MAXLINELENGTH,input);
        if(readcheck != NULL) i++;
    }

    fclose(input);







}




#endif


