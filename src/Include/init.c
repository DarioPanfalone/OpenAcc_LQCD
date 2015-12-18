#ifndef _INIT_C_
#define _INIT_C_


#include "./init.h"
#include "./common_defines.h"
#include "./markowchain.h"
#include "./fermion_parameters.h"
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

#define NPMGTYPES 8
#define MAXPMG 20
const char * par_macro_groups_name[] ={
    "Gauge Parameters",           // 0 
    "Flavour Parameters",         // 1   
    "Background Field Parameters",// 2            
    "MD parameters",              // 3
    "Montecarlo Parameters",      // 4       
    "Gauge Measures Settings",    // 5
    "Fermion Measures Settings"   // 6
};
#define PMG_GAUGE      0
#define PMG_FERMION    1
#define PMG_BACKGROUND 2
#define PMG_MD         3
#define PMG_MC         4
#define PMG_GMEAS      5
#define PMG_FMEAS      6
#define PMG_DEVICE     7



void erase_comments(char file_lines[MAXLINES][MAXLINELENGTH], int maxlines)
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

int scan_group_V(int ntagstofind, const char **strtofind, 
        int *tagcount,
        int *taglines, int *tagtypes, int maxnres,
        char filelines[MAXLINES][MAXLINELENGTH], 
        int startline, int endline)
    // returns number of groups found, writes arrays
{
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
        int found_something;
        for(int itype =0; itype <ntagstofind; itype++){
            found_something = 
                strcasestr(filelines[iline],strtofind[itype]);
            if(found_something){
                printf("Found group %s \n",strtofind[itype]);
                taglines[nres] = iline;
                tagtypes[nres] = itype;
                nres++;
                }
                else tagcount[itype]++;
            } 
        iline++;
        }
    for(int i =0; i<ntagstofind; i++)
        printf("Found %s %d times.\n", strtofind[i], tagcount[i]); 

    return nres;
}

void scan_group_NV(int npars,par_info* par_infos,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
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
        int found_something;
        for(int i =0; i<npars; i++){
            found_something = strcasestr(filelines[iline],par_infos[i].name);
            if(found_something){
               // found parameter
               printf("Found %s = ",par_infos[i].name);
               int reads = 0;
               char parname[20];
               switch(par_infos[i].type){
                   case TYPE_INT: 
                       reads = sscanf(filelines[iline],
                               "%s %d",parname,(int*)par_infos[i].par);
                       if(reads == 2) 
                           printf("%d\n", *((int*)par_infos[i].par));
                       break;
                   case TYPE_DOUBLE:
                       reads = sscanf(filelines[iline],
                               "%s %f",parname,(double*) par_infos[i].par);
                       if(reads == 2)
                           printf("%f\n", *((double*)par_infos[i].par));
                       break;
                   case TYPE_STR: 
                       reads = sscanf(filelines[iline],
                               "%s %s",parname,(char*) par_infos[i].par);
                       if(reads == 2) 
                           printf("%s\n", *((char*)par_infos[i].par));
                       break;
                   default: 
                       printf("WARNING, variable type not set in sourcecode.\n");
                       break;

                       if(reads == 1) rc[i] = 1;
                       else printf("WARNING, NO VALUE READ!");

                       break;
               }

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


void read_flavour_info(ferm_param *flpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

  // see /Include/fermion_parameters.h
  par_info fp[6];

  char sferm_mass[]        = "ferm_mass"        ;
  char sdegeneracy[]       = "degeneracy"       ;      
  char snumber_of_ps[]     = "number_of_ps"     ;
  char sname[]             = "name"             ;
  char sferm_charge[]      = "ferm_charge"      ;
  char sferm_im_chem_pot[] = "ferm_im_chem_pot" ;

  fp[0]=(par_info){(void*) &(flpar->ferm_mass       ),TYPE_DOUBLE, sferm_mass       };
  fp[1]=(par_info){(void*) &(flpar->degeneracy      ),TYPE_INT   , sdegeneracy      };
  fp[2]=(par_info){(void*) &(flpar->number_of_ps    ),TYPE_INT   , snumber_of_ps    };
  fp[3]=(par_info){(void*) &(flpar->name            ),TYPE_STR   , sname            };
  fp[4]=(par_info){(void*) &(flpar->ferm_charge     ),TYPE_DOUBLE, sferm_charge     };
  fp[5]=(par_info){(void*) &(flpar->ferm_im_chem_pot),TYPE_DOUBLE, sferm_im_chem_pot};

  scan_group_NV(6,fp, filelines, startline, endline);


}
void read_gauge_info(double *beta,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{
  // see OpenAcc/su3_measurements.h
  par_info gp[1];

  char sbeta[]        = "beta"        ;
  gp[0]=(par_info){(void*) beta,TYPE_DOUBLE, sbeta };
  scan_group_NV(1,gp, filelines, startline, endline);

}
void read_backfield_info(bf_param *bfpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

  // see /OpenAcc/backfield.h
  par_info bfp[6];

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

  scan_group_NV(6,bfp, filelines, startline, endline);

}
void read_md_info(md_param *mdpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

  // see /OpenAcc/md_integrator.h
  par_info mdp[3];

  char snomd[] = "NmdSteps" ;
  char sgs[] = "GaugeSubSteps" ;
  char st[] = "TrajLength" ;

  mdp[0]=(par_info){(void*) &(mdpar->no_md ),TYPE_INT, snomd };
  mdp[1]=(par_info){(void*) &(mdpar->gauge_scale ),TYPE_INT, sgs };
  mdp[2]=(par_info){(void*) &(mdpar->t ),TYPE_DOUBLE, st};

  scan_group_NV(3,mdp, filelines, startline, endline);

}
void read_mc_info(mc_param *mcpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

  // see /OpenAcc/md_integrator.h
  par_info mcp[9];

  char sntraj[] = "Ntraj" ;
  char stherm_ntraj[] = "ThermNtraj" ;
  char ssaveconfinterval[] = "SaveConfInterval" ;
  char ssaverunningconfinterval[] = "SaveRunningConfInterval";
  char sresidue_metro[] = "residue_metro";
  char sresidue_md[] = "residue_md";
  char sstore_conf_name[] = "StoreConfName";
  char ssave_conf_name[] = "SaveConfName";
  char sseed[] = "Seed";

  mcp[0]=(par_info){(void*) &(mcpar->ntraj                  ),TYPE_INT,sntraj          };
  mcp[1]=(par_info){(void*) &(mcpar->therm_ntraj            ),TYPE_INT,stherm_ntraj    };
  mcp[2]=(par_info){(void*) &(mcpar->saveconfinterval       ),TYPE_INT, ssaveconfinterval};
  mcp[3]=(par_info){(void*) &(mcpar->saverunningconfinterval),TYPE_INT, ssaverunningconfinterval};
  mcp[4]=(par_info){(void*) &(mcpar->residue_metro    ),TYPE_DOUBLE,sresidue_metro};
  mcp[5]=(par_info){(void*) &(mcpar->residue_md       ),TYPE_DOUBLE,sresidue_md};
  mcp[6]=(par_info){(void*) &(mcpar->store_conf_name  ),TYPE_STR,sstore_conf_name};
  mcp[7]=(par_info){(void*) &(mcpar->save_conf_name   ),TYPE_STR,ssave_conf_name};
  mcp[8]=(par_info){(void*) &(mcpar->seed   ),TYPE_INT,sseed};

  scan_group_NV(9,mcp, filelines, startline, endline);

}
void read_gaugemeas_info(char *outfilename,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

  // see /Meas
  par_info gmp[1];

  char soutfilename[] = "Gauge Outfilename" ;

  gmp[0]=(par_info){(void*) outfilename ,TYPE_STR,soutfilename };

  scan_group_NV(1,gmp, filelines, startline, endline);

}
void read_fermmeas_info(char * outfilename,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

  par_info fmp[1];

  char soutfilename[] = "Fermionic Outfilename" ;

  fmp[0]=(par_info){(void*) outfilename ,TYPE_STR,soutfilename };

  scan_group_NV(1,fmp, filelines, startline, endline);

}

void read_tech_setting(tech_param *tech_settings,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

  par_info tp[1];

  char sdevice_choice[] = "device_choice" ;

  tp[0]=(par_info){(void*) &(tech_settings->device_choice),TYPE_INT,sdevice_choice };

  scan_group_NV(1,tp, filelines, startline, endline);

}


void set_global_vars_and_fermions_from_input_file(const char* input_filename)
{

    // Opening filenames and reading it
    FILE *input = fopen(input_filename,"r");
    if (input == NULL) {
        printf("Could not open file %s \n",input_filename );
        exit(1);
    }

    char filelines[MAXLINES][MAXLINELENGTH];
    char *readcheck = filelines[0]; int lines_read = 0;
    while(readcheck != NULL){
        readcheck = fgets(filelines[lines_read],MAXLINELENGTH,input);
        if(readcheck != NULL) lines_read++;
    }
    
    // erasing comments 
    erase_comments(filelines,lines_read);

    //scanning for macro parameter families

    int tagpositions[MAXPMG], tagtypes[MAXPMG],tagcounts[MAXPMG];
    int found_tags = scan_group_V(NPMGTYPES,par_macro_groups_name,
           tagcounts,
           tagpositions,tagtypes,MAXPMG,
           filelines,0,lines_read );

    // see global var in /Include/fermion_parameters.
    // setting NDiffFlavs first
    NDiffFlavs = tagcounts[PMG_FERMION];
    fermions_parameters = (ferm_param*) malloc(NDiffFlavs*sizeof(ferm_param));
    int fermion_count = 0;
    for(int igroup  = 0 ; igroup < found_tags; igroup++){
        int startline = tagpositions[igroup];
        int endline = (igroup<found_tags-1)?tagpositions[igroup+1]:lines_read;

        printf("Reading %s...\n", par_macro_groups_name[tagtypes[igroup]]);
        switch(tagtypes[igroup]){
            case PMG_GAUGE     :
                read_gauge_info(&beta,filelines,startline,endline);
                break; 
            case PMG_FERMION   : 
                read_flavour_info(&(fermions_parameters[fermion_count]),
                            filelines,startline,endline);
                fermion_count++;
                break; 
#ifdef BACKFIELD
            case PMG_BACKGROUND: 
                read_backfield_info(&backfield_parameters, 
                        filelines,startline,endline);
                break; 
#endif
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
                read_tech_setting(&dev_settings,
                        filelines,startline,endline);
                break; 

        }
    }

    fclose(input);

}




#endif


