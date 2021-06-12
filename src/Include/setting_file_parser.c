#ifndef _SETTING_FILE_PARSER_C_
#define _SETTING_FILE_PARSER_C_

#include "../Meas/ferm_meas.h"
#include "../Meas/gauge_meas.h"
#include "../Meas/measure_topo.h"
#include "../Mpi/multidev.h"
#include "../OpenAcc/action.h"
#include "../OpenAcc/md_parameters.h"
#include "../OpenAcc/backfield.h"
#include "../OpenAcc/geometry.h" // for MULTIDEVICE to be defined or not
#include "../OpenAcc/su3_measurements.h"
#include "../RationalApprox/rationalapprox.h"
#include "./common_defines.h"
#include "./debug.h"
#include "./fermion_parameters.h"
#include "./hash.h"
#include "./rep_info.h"
#include "./montecarlo_parameters.h"
#include "./setting_file_parser.h"
#include "./inverter_tricks.h"
#include "../tests_and_benchmarks/test_and_benchmarks.h"
#include "../OpenAcc/alloc_settings.h"

#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <math.h>

#include <inttypes.h>
#include <stdint.h>

#define MAXLINES 300
#define MAXLINELENGTH 500 // pretty long to accomodate all the comments

#define MAXCRLENGTH 300

char input_file_str[MAXLINES*MAXLINELENGTH];

//types that can be found in an input file
enum dtype {TYPE_INT, TYPE_DOUBLE, TYPE_STR, NUM_TYPES}; //i tipi  che si trovano nel file. é un enum (0,3)
const char * type_strings[]={"(int)", "(double)", "(string)" }; //vettore di stringhe non modificabile.


//qui è definita questa struct, che verrà usata per leggere i valori dal file.
typedef struct par_info_t{

    void* par;  //The void pointer in C is a pointer which is not associated with any data types. It points to some data location in the storage means points to the address of variables.
    enum dtype type; //da il tipo
    const char* name; // nome della variabile
    const void* default_value;//il valore di default
    const char* comment; //eventuale commento

}par_info;

#define MAXPMG 20

//definisce il nome dei gruppi di variabili.(MOD)
const char * par_macro_groups_names[] ={ // NPMTYPES strings, NO SPACES!
    "ActionParameters",         // 0 
    "FlavourParameters",        // 1   
    "BackgroundFieldParameters",// 2            
    "MDParameters",             // 3
    "MontecarloParameters",     // 4       
    "GaugeMeasuresSettings",    // 5
    "FermionMeasuresSettings",  // 6
    "TopoMeasuresSettings"   ,  // 7
    "DeviceSettings"         ,  // 8
    "Geometry"               ,  // 9
    "DebugSettings"          ,  // 10
    "InverterTricks"         ,  // 11
    "TestSettings"           ,   // 12
    "Replicas numbers"          //13
};
enum pmg_types {  //ricorda enum gli da valore da 0 a n-1. Perciò l'ultimo è davvero il numero di tutti i gruppi. qui va aggiunto replicas. (MOD)
    PMG_ACTION         ,
    PMG_FERMION        ,
    PMG_BACKGROUND     ,
    PMG_MD             ,
    PMG_MC             ,
    PMG_GMEAS          ,
    PMG_FMEAS          ,
    PMG_TMEAS	       ,
    PMG_DEVICE         ,
    PMG_GEOMETRY       ,
    PMG_DEBUG          ,
    PMG_INVERTER_TRICKS,
    PMG_TESTS          ,
    PMG_REPLICAS       ,
    NPMGTYPES};
    

FILE * helpfile; //definisce l'help file.
char IGNORE_IT[50];

//FUNZIONE 1: appiccica tutto e ridà la launghezza totale del file in caratteri. (NON MOD)

// just to save it in the conf
//i parametri sono la matrice con dentro il contenuto del file input, un int  del numero massimo di linee che è di interesse scrivere, il nome del file da mandare in input


int prepare_string_from_stringarray(char file_lines[MAXLINES][MAXLINELENGTH], 
        int maxlines,
        char* input_filename_to_write){

    int iline;int totlen = 0; //contatori
    for(iline = 0; iline < maxlines ; iline++){
        strcat(input_filename_to_write,file_lines[iline]); //aggiunge le righe una a una.
        totlen += strlen(file_lines[iline]);
    }
    return totlen;

}

//FUNZIONE 2: Funzione che cancella i commenti.(trova # ed elimina cioò che viene dopo.) (NON MOD)

void erase_comments(char file_lines[MAXLINES][MAXLINELENGTH], int maxlines)
{
    //cerca le righe fino a maxlines.
    for(int i = 0; i<maxlines; i++){
        char* start_comment_ptr = index(file_lines[i],'#' ); //questa funzione index appartiene a strings.h e restituisce la stringa che inizia con #
        //immagino che l'if cosi voglia dire, se la funzione è definita.
        if(start_comment_ptr){
            int comment_start = start_comment_ptr - file_lines[i]; //se il carattere # si trova al 15 carattere della riga, il risulatato sarà 15
            file_lines[i][comment_start] = '\n';
            file_lines[i][comment_start+1] = '\0'; //it stands for null.
            for(int j = comment_start+2; j < MAXLINELENGTH; j++)
                file_lines[i][j] = ' ';
        }
    }
}

//FUNZIONE 3: funzione che serve per l'omonimia tra parametri. Non mi serve. (NON MOD)

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

//FUNZIONE 4: CERCA I GRUPPI DI PARAMETRI E RESTITUISCE IL NUMERO DI GRUPPI TROVATI.(NO MOD.)

//ovviamente a parametro gli passa par_macro_groups_names[], quindi mi basta modificare questo parametro. e anche ntagstofind.
int scan_group_V(int ntagstofind, const char **strtofind, 
        int *tagcount, //Vettore definito per contare i gruppi.
        int *taglines, int *tagtypes, int maxnres, //maxres massimo numero dei risultati
        char filelines[MAXLINES][MAXLINELENGTH], 
        int startline, int endline)
{
    // returns number of groups found, writes arrays
    // scans for lines in the format
    // (stuff that will be ignored )NAME (stuff that will be ignored)
    // the arrays 'taglines' and 'types' and must be maxnres long

    if(0==devinfo.myrank)
        printf("Scanning for parameter macro groups..\n");
    for(int i =0; i<ntagstofind; i++) tagcount[i] = 0; //manda a zero gli elementi di quel vettore.
    for(int i =0; i<maxnres; i++){// initializing output arrays
        taglines[i] = -1; //inizializza i vettori in output a -1.
        tagtypes[i] = -1;
    }
    int nres = 0;
    int iline = startline; //linea da cui parte la ricerca. //ricorda filelines è la matrice in cui viene copiato il file.
    while(iline<endline && nres < maxnres ){
        char * found_something;
        for(int itype =0; itype <ntagstofind; itype++){
            found_something = strstr(filelines[iline],strtofind[itype]); //strtofind è un'altra matrice che viene data da input.
            if(found_something){
                printf("Found group %s on line %d\n",strtofind[itype], iline);
                taglines[nres] = iline; //alla fine tagline è un vettore che segna la posizione dei gruppi di macro.
                tagtypes[nres] = itype; //tagtypes è una roba che segna il type segnato.
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

//FUNZIONE 5

int scan_group_NV(int npars,par_info* par_infos,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{   
    if(startline >= endline){ // goes into 'help mode'
        if(0==devinfo.myrank)for(int ipar = 0;ipar< npars ; ipar++){
            // printing name
            fprintf(helpfile,"%-30s",par_infos[ipar].name);
           
            // printing default value (if present)
            if(par_infos[ipar].default_value) switch(par_infos[ipar].type){
                case TYPE_INT:
                    fprintf(helpfile,"%-20d",*((const int *) par_infos[ipar].default_value ));
                    break;
                case TYPE_DOUBLE:
                    fprintf(helpfile,"%-20e",*((const double *) par_infos[ipar].default_value));
                    break;
                case TYPE_STR:
                    fprintf(helpfile,"%-20s",((const char *) par_infos[ipar].default_value));
                    break;
            }
            else fprintf(helpfile,"%-20s","(set yours!)");

            // printing type
            fprintf(helpfile,"#%s\n",type_strings[par_infos[ipar].type]);
            if(par_infos[ipar].comment)
                fprintf(helpfile,"%s\n\n",par_infos[ipar].comment);

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

        while(iline < endline )
        {
            char * found_something = NULL;
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
                            case TYPE_INT: //il case che gli da scannerizza la striga con il nome e il codice.
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
                    if(par_infos[i].default_value){
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


//FUNZIONE 6: Legge le informazioni relative ai flavour. //questo ferm_param non si capisce dove lo definisce va beh meno male che per ora i fermioni non servono.

int read_flavour_info(ferm_param *flpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    // see /Include/fermion_parameters.h
    
    //vediamo che definisce par_info utilizzando la struttura flpar
    //basta passargli dentro filelines, inizio e fine righe.

    par_info fp[]={
        (par_info){(void*) &(flpar->ferm_mass       ),TYPE_DOUBLE, "Mass"          , NULL,NULL},
        (par_info){(void*) &(flpar->degeneracy      ),TYPE_INT   , "Degeneracy"    , NULL,NULL},
        (par_info){(void*) &(flpar->number_of_ps    ),TYPE_INT   , "PseudoFermions", NULL,NULL},
        (par_info){(void*) &(flpar->name            ),TYPE_STR   , "Name"          , NULL,NULL},
        (par_info){(void*) &(flpar->ferm_charge     ),TYPE_DOUBLE, "Charge"        , NULL,NULL},
        (par_info){(void*) &(flpar->ferm_im_chem_pot),TYPE_DOUBLE, "MuOverPiT"     , NULL,NULL}};

    // from here on, you should not have to modify anything.
    return scan_group_NV(sizeof(fp)/sizeof(par_info),fp, filelines, startline, endline);

}












//anche qui questo action_param non ho capito dove lo definisce.
int read_action_info(action_param *act_par,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{
    // see OpenAcc/su3_measurements.h
    const double stout_rho_def = RHO;
	const int topo_action_def = 0;
	const double barrier_def = 2.0;
	const double width_def = 0.1;
	const char topo_file_path_def[] = "ToPotential";
	const int topo_stout_steps_def = 0;
	const char topo_act_comment[] = "#Set to 1 to enable multicanonical topological potential.";
    par_info ap[]={
        (par_info){(void*) &(act_par->beta)            , TYPE_DOUBLE, "Beta"          , NULL, NULL},
        (par_info){(void*) &(act_par->stout_steps)     , TYPE_INT   , "StoutSteps"    , NULL, NULL},
        (par_info){(void*) &(act_par->stout_rho)       , TYPE_DOUBLE, "StoutRho"      , (const void*) &stout_rho_def, NULL},
		(par_info){(void*) &(act_par->topo_action)     , TYPE_INT   , "TopoAct"       , (const void*) &topo_action_def, topo_act_comment},
		(par_info){(void*) &(act_par->barrier)         , TYPE_DOUBLE, "Barrier"       , (const void*) &barrier_def, NULL},
		(par_info){(void*) &(act_par->width)           , TYPE_DOUBLE, "Width"         , (const void*) &width_def, NULL},
		(par_info){(void*) &(act_par->topo_file_path)  , TYPE_STR   , "TopoPath"      , (const void*) &topo_file_path_def, NULL},
		(par_info){(void*) &(act_par->topo_stout_steps), TYPE_INT   , "TopoStoutSteps", (const void*) &topo_stout_steps_def, NULL},
		(par_info){(void*) &(act_par->topo_rho)        , TYPE_DOUBLE, "TopoRho"       , NULL, NULL}
	};			   
				   
    // from here on, you should not have to modify anything.
    int res = scan_group_NV(sizeof(ap)/sizeof(par_info),ap, filelines, startline, endline);

    if(startline<endline){
        if(act_par->stout_rho != RHO ){ 

            if(0==devinfo.myrank){
                printf("Error, input file stout_rho != RHO \n");
                printf("  Either modify the input file, or recompile changing RHO\n");
                printf(" (input) stout_rho = %f, (code) RHO = %f\n", act_par->stout_rho,RHO);
            }
            res = 1;

        }

        if(act_par->topo_rho != TOPO_RHO ){ 

            if(0==devinfo.myrank){
                printf("Error, input file topo_rho != TOPO_RHO \n");
                printf("  Either modify the input file, or recompile changing TOPO_RHO\n");
                printf(" (input) topo_rho = %f, (code) TOPO_RHO = %f\n", act_par->topo_rho,TOPO_RHO);
            }
            res = 1;

        }
    }
    return res;



}
int read_backfield_info(bf_param *bfpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    // see /OpenAcc/backfield.h
    const double default_value = 0.0;
    par_info bfp[]={
        (par_info){(void*) &(bfpar->ex ),TYPE_DOUBLE, "ex", (const void*) &default_value  ,NULL},
        (par_info){(void*) &(bfpar->ey ),TYPE_DOUBLE, "ey", (const void*) &default_value  ,NULL},
        (par_info){(void*) &(bfpar->ez ),TYPE_DOUBLE, "ez", (const void*) &default_value  ,NULL},
        (par_info){(void*) &(bfpar->bx ),TYPE_DOUBLE, "bx", (const void*) &default_value  ,NULL},
        (par_info){(void*) &(bfpar->by ),TYPE_DOUBLE, "by", (const void*) &default_value  ,NULL},
        (par_info){(void*) &(bfpar->bz ),TYPE_DOUBLE, "bz", (const void*) &default_value  ,NULL}};

    // from here on, you should not have to modify anything.
    return scan_group_NV(sizeof(bfp)/sizeof(par_info),bfp, filelines, startline, endline);

}
int read_md_info(md_param *mdpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    // see /OpenAcc/md_integrator.h
    const double tlendef = 1.0;
    const double expmaxeigenv_def = 5.5 ; 
    const int singlePrecMDdef = 0;
    const int max_cg_iterations_def = 10000;
    const int recycleInvsForceDef = 0;
    const int extrapolateInvsForceDef = 0;
    
    const char singlePrecMD_comment[] = "#For volumes large as or larger than about 24^4 double precision MD is needed. ";
    const char residue_md_comment[] = "#If using single precision MD, the residue can't reliably be set to less than 1e-4" ;
    const char recycleInvsForce_comment[] = "# For chronological inverters - not implemented - keep 0 or remove";
    const char extrapolateInvsForce_comment[] = "# For chronological inverters - not implemented - keep 0 or remove";
    const char expmaxeigenv_comment[] = "# 5.5 is ok for 2 stout, with less stouting this number must be higher. " ; 

    par_info mdp[]={
        (par_info){(void*) &(mdpar->no_md ),       TYPE_INT, "NmdSteps"     , NULL,NULL},
        (par_info){(void*) &(mdpar->gauge_scale ), TYPE_INT, "GaugeSubSteps", NULL,NULL},
        (par_info){(void*) &(mdpar->t ),        TYPE_DOUBLE, "TrajLength"   , (const void*) &tlendef,NULL},
        (par_info){(void*) &(mdpar->residue_metro),       TYPE_DOUBLE,   "residue_metro"          , NULL,NULL},
        (par_info){(void*) &(mdpar->expected_max_eigenvalue),TYPE_DOUBLE,"ExpMaxEigenvalue"       , (const void*) &expmaxeigenv_def,expmaxeigenv_comment},
        (par_info){(void*) &(mdpar->singlePrecMD),TYPE_INT , "SinglePrecMD", (const void*) &singlePrecMDdef,singlePrecMD_comment},
        (par_info){(void*) &(mdpar->residue_md),TYPE_DOUBLE, "residue_md"   , NULL,residue_md_comment},
        (par_info){(void*) &(mdpar->max_cg_iterations),TYPE_INT, "MaxCGIterations" , (const void*) &max_cg_iterations_def,NULL},
        (par_info){(void*) &(mdpar->recycleInvsForce),TYPE_INT, "recycleInvsForce" , (const void*) &recycleInvsForceDef,recycleInvsForce_comment},
        (par_info){(void*) &(mdpar->extrapolateInvsForce),TYPE_INT, "extrapolateInvsForce", (const void*) &extrapolateInvsForceDef,extrapolateInvsForce_comment}};

    // from here on, you should not have to modify anything.
    return scan_group_NV(sizeof(mdp)/sizeof(par_info),mdp, filelines, startline, endline);

}
int read_inv_tricks_info(inv_tricks *invinfo,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline){

    const int singlePInvAccelMultiInvDef = 0;
    const int useMixedPrecisionDef = 0;
    const int restartingEveryDef = 10000;
    const double mixedPrecisionDeltaDef = 0.1; 

    const char singlePInvAccelMultiInv_comment[] = "# At present, multishift inverter is faster. Keep 0 or delete.";

    par_info iip[]={
        (par_info){(void*) &(invinfo->singlePInvAccelMultiInv),TYPE_INT,"singlePInvAccelMultiInv", (const void*) &singlePInvAccelMultiInvDef,singlePInvAccelMultiInv_comment},
        (par_info){(void*) &(invinfo->useMixedPrecision),TYPE_INT,"useMixedPrecision", (const void*) &useMixedPrecisionDef,NULL},
        (par_info){(void*) &(invinfo->restartingEvery),TYPE_INT,"restartingEvery", (const void*) &restartingEveryDef,NULL},
        (par_info){(void*) &(invinfo->mixedPrecisionDelta), TYPE_DOUBLE, "mixedPrecisionDelta", (const void*)&mixedPrecisionDeltaDef,NULL}};
    return scan_group_NV(sizeof(iip)/sizeof(par_info),iip, filelines, startline, endline);

}




int read_mc_info(mc_params_t *mcpar,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    // see /OpenAcc/md_integrator.h
    const int seed_def = 0;  // which means time()
    const double epsgen_def = 0.1 ; 
    const char RandGenStatusFilename_def[] = "rgstatus.bin"; 
    const double MaxRunTimeS_def = 1.0e9; // 30 years should be enough
    const int MaxConfIdIter_def = 1000000;
    const int JarzynskiMode_def = 0;

    const char JarzynskiMode_comment[]="#0 - normal operation,\n#1 - from bz(as set in input file) to bz+1\n#-1 - from bz to bz-1";
    const char MaxConfIdIter_comment[]="#In Jarzynski mode, this number will be taken as the number of steps\n#to go from the starting value of bz to the last.";
    const char therm_ntraj_comment[] = 
        "# When starting from the identity links, metropolis test will likely fail for a few trajectories.\n#Set this number to a few to allow the run to thermalise.";
    const char eps_gen_comment[] = "#The level of randomness in creating a new random gauge conf";
    const char store_conf_name_comment[] = "#This is the prefix of the filenames of the gauge conf files which will be saved (with an index).";
    const char save_conf_name_comment[] = "#This is the name of the starting gauge conf file. THIS FILE WILL BE OVERWRITTEN.";
    const char seed_comment[] = 
        "# set this to 42 \n\
#(https://en.wikipedia.org/wiki/Phrases_from_The_Hitchhiker\%27s_Guide_to_the_Galaxy#Answer_to_the_Ultimate_Question_of_Life.2C_the_Universe.2C_and_Everything_.2842.29)";
    const char RandGenStatusFilename_comment[] = "# The status of the random number generator will be saved in a file named like this (plus suffixes).\n\
# In correspondence of every stored gauge conf, the random number generator statuses will be saved for each MPI rank.\n\
# With these files, the gauge conf and the right setting file (like this one) reproducibility should be assured.";

    const char statusFileName_def[] = "program_status.txt"; 
    const char statusFileName_comment[] = "# The status of the program at the end will be saved here.\n\
# This is useful to run the program on \"short\" queues."; 

    par_info mcp[]={
        (par_info){(void*) &(mcpar->ntraj                  ),TYPE_INT,   "Ntraj"                  , NULL,NULL},
        (par_info){(void*) &(mcpar->therm_ntraj            ),TYPE_INT,   "ThermNtraj"             , NULL,therm_ntraj_comment},
        (par_info){(void*) &(mcpar->storeconfinterval      ),TYPE_INT,   "StoreConfInterval"      , NULL,NULL},
        (par_info){(void*) &(mcpar->saveconfinterval),       TYPE_INT,   "SaveConfInterval"       , NULL,NULL},
        (par_info){(void*) &(mcpar->store_conf_name),        TYPE_STR,   "StoreConfName"          , NULL,store_conf_name_comment},
        (par_info){(void*) &(mcpar->save_conf_name),         TYPE_STR,   "SaveConfName"           , NULL,save_conf_name_comment},
        (par_info){(void*) &(mcpar->MaxConfIdIter),          TYPE_INT,   "MaxConfIdIter"          , (const void*) &MaxConfIdIter_def,MaxConfIdIter_comment},
        (par_info){(void*) &(mcpar->RandGenStatusFilename),  TYPE_STR,   "RandGenStatusFilename"  , (const void*) &RandGenStatusFilename_def,RandGenStatusFilename_comment},
        (par_info){(void*) &(mcpar->MaxRunTimeS),         TYPE_DOUBLE,   "MaxRunTimeS"            , (const void*) &MaxRunTimeS_def,NULL},
        (par_info){(void*) &(mcpar->seed),                   TYPE_INT,   "Seed"                   , (const void*) &seed_def,seed_comment},
        (par_info){(void*) &(mcpar->eps_gen),             TYPE_DOUBLE,   "EpsGen"                 , (const void*) &epsgen_def,eps_gen_comment},
        (par_info){(void*) &(mcpar->JarzynskiMode),             TYPE_INT,   "JarzynskiMode"       , (const void*) &JarzynskiMode_def,JarzynskiMode_comment},
        (par_info){(void*) &(mcpar->statusFileName),            TYPE_STR,   "StatusFileName"      , (const void*) &statusFileName_def,statusFileName_comment}
    };

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
    const int print_bfield_dbginfo_def = 0;
    const int md_diag_print_every_def = 10000;

    const char SaveAllAtEnd_comment[] = "# Set this to 0 if you want the program not to save its state at the end\n\
# (or it may overwrite some files and you won't be able to reproduce bugs/do other things)";
    // see /Meas
    
    //è un vettore di par_info.
    par_info gmp[]= {
        (par_info){(void*) &(dbg_settings->use_ildg),              TYPE_INT,"UseILDG"                , (const void*) &useildg_def,NULL},
        (par_info){(void*) &(dbg_settings->input_vbl),             TYPE_INT,"VerbosityLv"            , (const void*) &input_vbl_def,NULL},
        (par_info){(void*) &(dbg_settings->SaveAllAtEnd),          TYPE_INT,"SaveAllAtEnd"           , (const void*) &SaveAllAtEnd_def,SaveAllAtEnd_comment},
        (par_info){(void*) &(dbg_settings->print_bfield_dbginfo),  TYPE_INT,"PrintBackFieldDbgInfo"  , (const void*) &print_bfield_dbginfo_def,NULL},
        (par_info){(void*) &(dbg_settings->save_diagnostics),      TYPE_INT,"SaveDiagnostics"        , (const void*) &save_diagnostics_def,NULL},
        (par_info){(void*) &(dbg_settings->do_reversibility_test), TYPE_INT,"DoRevTest"              , (const void*) &do_reversibility_test_def,NULL},
        (par_info){(void*) &(dbg_settings->do_norandom_test),      TYPE_INT,"DoNoRandomTest"         , (const void*) &do_norandom_test_def,NULL},
        (par_info){(void*) &(dbg_settings->rng_fakeness_level),    TYPE_INT,"RngFakenessLevel"       , (const void*) &rng_fakeness_level_def,NULL},
        (par_info){(void*) &(dbg_settings->md_dbg_print_max_count),TYPE_INT,"MDDbgPrintMaxCount"     , (const void*) &md_dbg_print_max_count_def,NULL},
        (par_info){(void*) &(dbg_settings->diagnostics_filename),  TYPE_STR,"SaveDiagnosticsFilename", (const void*) &diagnostics_filename_def,NULL},
        (par_info){(void*) &(dbg_settings->md_diag_print_every) ,  TYPE_INT,"PrintDiagInfoEvery",      (const void*) &md_diag_print_every_def,NULL},
            
    };




    // from here on, you should not have to modify anything.
    return scan_group_NV(sizeof(gmp)/sizeof(par_info),gmp, filelines, startline, endline);

}





int read_gaugemeas_info(char *outfilename,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    // see /Meas
    par_info gmp[]= {
        (par_info){(void*) outfilename ,TYPE_STR, "GaugeOutfilename", NULL, NULL}};


    // from here on, you should not have to modify anything.
    return scan_group_NV(sizeof(gmp)/sizeof(par_info),gmp, filelines, startline, endline);

}
int read_fermmeas_info(ferm_meas_params * fmpars,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    const int doubleinv_def = 0;
    const int measEvery_def = 1;
    const char measEvery_comment[] = "#Fermionic measurements will be performed once very MeasEvery times.";
    const int printPlaqAndRect_def = 0;
    const char printPlaqAndRect_comment[] = "# Each measurement line will contain the value of the plaquette and the rectangle - reweighting will then be easier.";



    par_info fmp[]={
        (par_info){(void*) &(fmpars->fermionic_outfilename),       TYPE_STR,"FermionicOutfilename",        NULL ,                       NULL},
        (par_info){(void*) &(fmpars->measEvery),                   TYPE_INT,"MeasEvery"           ,(const void*) &measEvery_def,measEvery_comment},
        (par_info){(void*) &(fmpars->SingleInvNVectors),           TYPE_INT,"SingleInvNVectors"   ,        NULL ,                       NULL},
        (par_info){(void*) &(fmpars->DoubleInvNVectorsChiral),     TYPE_INT,"DoubleInvNVectorsChiral",     (const void*) &doubleinv_def,NULL},
        (par_info){(void*) &(fmpars->DoubleInvNVectorsQuarkNumber),TYPE_INT,"DoubleInvNVectorsQuarkNumber",(const void*) &doubleinv_def,NULL},
        (par_info){(void*) &(fmpars->printPlaqAndRect),TYPE_INT,"PrintPlaqAndRect",(const void*) &printPlaqAndRect_def,printPlaqAndRect_comment},
    
    };


    // from here on, you should not have to modify anything.
    return scan_group_NV(sizeof(fmp)/sizeof(par_info),fmp, filelines, startline, endline);

}


int read_topomeas_info(meastopo_param * meastopars,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    const int meascool_def = 0;
    const char pathcool_def[] = "TopoCool";
    const int coolmeasstep_def = 1;
    const int cool_measinterval_def = 1;
    const char cool_measinterval_comment[] = "# Cooled Topological charge will be measured every CoolMeasInterval cooling steps.";
    const int cooleach_def = 1;
    const char cooleach_comment[] = "# Cooled Topological charge will be measured every CoolMeasEach MC step.";
    
    const int measstout_def = 0;
    const char pathstout_def[] = "TopoStout";
    const double measrhostout_def = TOPO_RHO;
    const char measrhostout_comment[]="# MeasStoutRho can't be changed, it have to be equal to StoutRho, which is setted in src/Include/common_defines.h";
    const int stoutmeasstep_def = 1;
    const int stout_measinterval_def = 1;
    const char stout_measinterval_comment[] = "# Stouted Topological charge will be measured every StoutMeasInterval cooling steps.";
    const int stouteach_def = 1;
    const char stouteach_comment[] = "# Stouted Topological charge will be measured every StoutMeasEach MC step.";






    par_info tomp[]={
        (par_info){(void*) &(meastopars->meascool),       	       TYPE_INT,"MeasCool"            ,(const void*) &meascool_def,         NULL},
        (par_info){(void*) &(meastopars->pathcool),                    TYPE_STR,"PathCool"            ,(const void*) &pathcool_def,         NULL},
        (par_info){(void*) &(meastopars->coolmeasstep),                    TYPE_INT,"CoolMeasSteps"       ,(const void*) &coolmeasstep_def,         NULL},
        (par_info){(void*) &(meastopars->cool_measinterval),           TYPE_INT,"CoolMeasInterval"    ,(const void*) &cool_measinterval_def,cool_measinterval_comment},
        (par_info){(void*) &(meastopars->cooleach),                TYPE_INT,"CoolMeasEach"        ,(const void*) &cooleach_def,cooleach_comment},

        (par_info){(void*) &(meastopars->measstout),       	       TYPE_INT,"MeasStout"           ,(const void*) &measstout_def,         NULL},
        (par_info){(void*) &(meastopars->pathstout),                   TYPE_STR,"PathStout"           ,(const void*) &pathstout_def,         NULL},
        (par_info){(void*) &(meastopars->measrhostout),       	    TYPE_DOUBLE,"MeasStoutRho"        ,(const void*) &measrhostout_def,measrhostout_comment},
        (par_info){(void*) &(meastopars->stoutmeasstep),                   TYPE_INT,"StoutMeasSteps"      ,(const void*) &stoutmeasstep_def,         NULL},
        (par_info){(void*) &(meastopars->stout_measinterval),          TYPE_INT,"StoutMeasInterval"   ,(const void*) &stout_measinterval_def,stout_measinterval_comment},
        (par_info){(void*) &(meastopars->stouteach),               TYPE_INT,"StoutMeasEach"       ,(const void*) &stouteach_def,stouteach_comment},
    };


    // from here on, you should not have to modify anything.
    int res = scan_group_NV(sizeof(tomp)/sizeof(par_info),tomp, filelines, startline, endline);
    if(startline<endline){
        if(meastopars->measrhostout != TOPO_RHO ){ 

            if(0==devinfo.myrank){
                printf("Error, input file MeasStoutRho != TOPO_RHO \n");
                printf("  Either modify the input file, or recompile changing TOPO_RHO\n");
                printf(" (input) topo_rho = %f, (code) TOPO_RHO = %f\n", meastopars->measrhostout,TOPO_RHO);
            }
            res = 1;
        }
    }
    return res;
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
    
    const char nranks_read_comment[] = "# NRanks has been set at make time with the NR3 variable.";


    par_info tp[]= {
#ifndef MULTIDEVICE        
        (par_info){(void*) &(di->single_dev_choice),TYPE_INT,"device_choice",    NULL,                                 NULL},
        (par_info){(void*) IGNORE_IT,               TYPE_INT,"AsyncFermionComms",(const void*) &async_comm_fermion_def,NULL},
        (par_info){(void*) IGNORE_IT,               TYPE_INT,"AsyncGaugeComms",  (const void*) &async_comm_gauge_def,  NULL},
        (par_info){(void*) IGNORE_IT,               TYPE_INT,"NProcPerNode",     (const void*) &ignored_def,           NULL},
#else
        (par_info){(void*) &(di->single_dev_choice), TYPE_INT,"device_choice",    (const void*) &single_dev_choice_def, NULL},
        (par_info){(void*) &(di->async_comm_fermion),TYPE_INT,"AsyncFermionComms",(const void*) &async_comm_fermion_def,NULL},
        (par_info){(void*) &(di->async_comm_gauge),  TYPE_INT,"AsyncGaugeComms",  (const void*) &async_comm_gauge_def  ,NULL},
        (par_info){(void*) &(di->proc_per_node),     TYPE_INT,"NProcPerNode",     NULL, NULL},
#endif
        (par_info){(void*) &(di->nranks_read),       TYPE_INT,"NRanks",           NULL, nranks_read_comment}};

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

    const char dimensions_comment[] = "# Lattice dimensions have been fixed at make time through N0,N1,N2,N3 and NR3.";
    const char dimension_mapping_content[] = "#Notice that only the dimension mapped as 3 will be parallelized";

    par_info gp[]={
        (par_info){(void*) &(gpar->gnx ),TYPE_INT,  "nx"  ,(const void*) &nx_def  ,NULL},
        (par_info){(void*) &(gpar->gny ),TYPE_INT,  "ny"  ,(const void*) &ny_def  ,NULL},
        (par_info){(void*) &(gpar->gnz ),TYPE_INT,  "nz"  ,(const void*) &nz_def  ,NULL},
        (par_info){(void*) &(gpar->gnt ),TYPE_INT,  "nt"  ,(const void*) &nt_def  ,dimensions_comment},
        (par_info){(void*) &(gpar->xmap ),TYPE_INT, "xmap",(const void*) &xmap_def,NULL},
        (par_info){(void*) &(gpar->ymap ),TYPE_INT, "ymap",(const void*) &ymap_def,NULL},
        (par_info){(void*) &(gpar->zmap ),TYPE_INT, "zmap",(const void*) &zmap_def,NULL},
        (par_info){(void*) &(gpar->tmap ),TYPE_INT, "tmap",(const void*) &tmap_def,dimension_mapping_content}};

    int res = scan_group_NV(sizeof(gp)/sizeof(par_info),gp, filelines, startline, endline);

    if(startline < endline)
        if (1 == set_geom_glv(gpar) ) res = 1 ;

    return res;

}
//devo solo capire dove test_info e questi type cosi sono definiti.



int read_test_setting(test_info * ti,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline)
{

    int helpmode = (int) (startline == endline);

    const int saveResults_def = 0;
    const char saveResults_comment[] = "\
# set to 1 if you want to save all results (the results of Deo Doe or the result of \n\
# individual shifts)\n";

    const char fakeShift_comment[]="\
# tHe shift (equal for all terms) in the fake rational approximations in benchmarkMode for CG-M.";
    const char benchMarkMode_comment[]="\
# In benchmark mode, a fake rationan approximation will be used in CG-M.\n";

    par_info tp[]= {
        //ti è un puntatore a test_info.
        
        //qui si definisce questo vettore tp, che  ha come elementi, struct par_info e come elementi, 1 indirizzo in cui ti punta a deodoeiteration. 2 Type int, nome del coso, NULL(dovrebbe essere il default value), NULL(un eventuale commento).
        (par_info){(void*) &(ti->deoDoeIterations),       TYPE_INT,"DeoDoeIterations",NULL, NULL},
        (par_info){(void*) &(ti->multiShiftInverterRepetitions),TYPE_INT,"MultiShiftInverterRepetitions",NULL, NULL},
        (par_info){(void*) &(ti->fakeShift),TYPE_DOUBLE,"FakeShift",NULL,fakeShift_comment},
        (par_info){(void*) &(ti->benchmarkMode),TYPE_INT,"BenchmarkMode",NULL,benchMarkMode_comment},
        (par_info){(void*) &(ti->saveResults),TYPE_INT,"SaveResults",(const void*) &saveResults_def,saveResults_comment},
    
    };

    // from here on, you should not have to modify anything.
 
    int res = scan_group_NV(sizeof(tp)/sizeof(par_info),tp, filelines, startline, endline);

    if(!res) ti->parametersAreSet = 1; 

    return res;

}
//********************************************************************************************************************************************//
//********************************************************************************************************************************************//




//qui sto solo definendo una funziona. Here i modified the function by defining a new function "read_replicas_numbers", that will read the replicas_info from the input file.



int read_replicas_numbers(rep_info * re,char filelines[MAXLINES][MAXLINELENGTH], int startline, int endline){
    
    
    int i2=0;

    par_info rp[]={
        
        (par_info){(void*) &(re->replicas_total_number),TYPE_INT,"totalnumber",NULL, NULL},
        (par_info){(void*) &(re->defect_boundary),TYPE_INT,"defect_boundary",NULL, NULL},
        

        
        
    };
    printf("startline: %d  endline %d\n",startline,startline+3);
    
   int res = scan_group_NV(sizeof(rp)/sizeof(par_info),rp, filelines, startline, startline+3);
  
    
    
  
    alloc_info.num_replicas=re->replicas_total_number;

    
    par_info rp1[6];
    /*rp2=malloc(alloc_info.num_replicas*sizeof(par_info));*/
    
    
    char str0[13];//uno in più perchè deve avere il carattere di fine stringa.
    
    for(i2=0;i2<6;i2++){
    
        
        rp1[i2].par= &(re->defect_coordinates[i2]);
        
        
        rp1[i2].type=TYPE_INT;
        
       
        snprintf(str0,13,"defect_size%d",i2);
        
        
        rp1[i2].name=malloc(sizeof(char)*20);
        //non puoi eguagliare due stringhe perchè sono due puntatori e tu vuoi modificare invece il contenuto della stringa-
        strcpy(rp1[i2].name,str0);
        
        rp1[i2].default_value=NULL;
        rp1[i2].comment=NULL;
        
    }
     startline=startline+2;
    printf("startline: %d  endline %d\n",startline,startline+7);
    
    
    int res1 = scan_group_NV(6,rp1, filelines, startline, startline+7);
    if(res1!=0){ res=res1;}
    
    re->cr_vet=malloc(alloc_info.num_replicas*sizeof(double));
    
  
    par_info *rp2 ;
    rp2=malloc(alloc_info.num_replicas*sizeof(par_info));
   
    char str[10];//uno in più perchè deve avere il carattere di fine stringa.
    
       for(i2=0;i2<(alloc_info.num_replicas);i2++){
           
           rp2[i2].par= &(re->cr_vet[i2]);
          
           rp2[i2].type=TYPE_DOUBLE;
           
    
           snprintf(str,10,"cr_value%d",i2);
           
           
            rp2[i2].name=malloc(sizeof(char)*20);
         //non puoi eguagliare due stringhe perchè sono due puntatori e tu vuoi modificare invece il contenuto della stringa-
           strcpy(rp2[i2].name,str);
          
           rp2[i2].default_value=NULL;
           rp2[i2].comment=NULL;
           
       }
    startline=startline+6;
    
    
 
    
  int res2 = scan_group_NV(alloc_info.num_replicas,rp2, filelines, startline, endline);
  if(res2!=0){ res=res2;}
   

        
  
   


    free(rp2);
    //checking right choice of axis.
    if(re->defect_boundary>3 || re->defect_boundary<0){
        printf("wrong defect axis choice!\n");
        res=1;
    }
    int target_def_dim[3];
    
    switch (re->defect_boundary) {
        case 0:
            target_def_dim[0]=nd1;
            target_def_dim[1]=nd2;
            target_def_dim[2]=nd3;
            
            break;
    
        case 1:
            target_def_dim[0]=nd0;
            target_def_dim[1]=nd2;
            target_def_dim[2]=nd3;
    
            break;
            
        case 2:
            target_def_dim[0]=nd0;
            target_def_dim[1]=nd1;
            target_def_dim[2]=nd3;
            
            break;

        case 3:
            target_def_dim[0]=nd0;
            target_def_dim[1]=nd1;
            target_def_dim[2]=nd2;
            
            break;

        default:
            printf("ERROR WRONG AXIS CHOICE!\n");
            break;
            
    
    
    }
    
    for(i2=0;i2<3;i2++){
        
        
        
        
    if( re->defect_coordinates[2*i2]<0 || re->defect_coordinates[2*i2+1]>target_def_dim[i2] ){
        printf("wrong defect coordinates choice!\n");
        res=1;
    }
    }
    
    
    return res;
}

//********************************************************************************************************************************************//
//********************************************************************************************************************************************//




//LE FUNZIONI PRECEDENTI INVECE LEGGONO E STAMPANO A VIDEO.
//QUESTA FUNZIONE LEGGE DA FILE E DA IL VALORE AD ALLOC_INFO.

int set_global_vars_and_fermions_from_input_file(const char* input_filename)
{

    // Opening filenames and reading it
    int helpmode = 0;
    int i2=0;
    char filelines[MAXLINES][MAXLINELENGTH]; //questi sono delle stringhe.
    FILE *input = fopen(input_filename,"r"); //questo ovviamente apre il file
    if (input == NULL){

        printf("MPI%02d: Could not open file %s \n",devinfo.myrank,input_filename );
        if(0==devinfo.myrank){
            printf("writing an template_input file for your convenience.\n" );
        }
        helpmode = 1;
    }

    
    //in questa parte legge il file  e lo copia in taglines.
    int lines_read = 0;

    int tagpositions[MAXPMG], tagtypes[MAXPMG],tagcounts[NPMGTYPES]; //qui definisce
    int found_tags = 0;
    fermions_parameters = NULL;
    if (! helpmode){//reading input file 
        char *readcheck = filelines[0];
        while(readcheck != NULL){ // saving input file in array of lines
            readcheck = fgets(filelines[lines_read],MAXLINELENGTH,input);
            if(readcheck != NULL){ lines_read++;
                printf("linea %d",lines_read);
        
            }
        }
        
        fclose(input); //qui chiude il file

        if(0==devinfo.myrank)
            printf("lines read: %d\n", lines_read); //printa il numero di linee lette.
       

        int totlen = prepare_string_from_stringarray(filelines,lines_read,input_file_str); // mi da la lunghezza totale di quanto letto.

        if(0==devinfo.myrank)
            printf("Written %d characters into input_file_str\n",totlen);


        // erasing comments 
        erase_comments(filelines,lines_read);  //cancella i commmenti da tutte le linee.

        //scanning for macro parameter families //cerca tutti itag
        found_tags = scan_group_V(NPMGTYPES,par_macro_groups_names,  //trova i tag. //il numero di type da trovare viene definito dall'enum di sopra.
                tagcounts,
                tagpositions,tagtypes,MAXPMG, //MAXPMG é definito 20
                filelines,0,lines_read );

        
        printf("last\n");
        rep =(rep_info*) malloc(sizeof(rep_info)); //funziona?
        
       
        
        
        // see global var in /Include/fermion_parameters.
        // setting alloc_info.NDiffFlavs first
        alloc_info.NDiffFlavs = tagcounts[PMG_FERMION]; //tagcount PMG_Fermions è il numero di flavour
        
        if(alloc_info.NDiffFlavs==0)
        {
            fermions_parameters = NULL; //se non trova i parametri fermionici
            if(0==devinfo.myrank){
                printf("NO FERMIONS FOUND, ");
                printf("SIMULATING PURE GAUGE THEORY...\n");
            }
        }
        else fermions_parameters = (ferm_param*) malloc(alloc_info.NDiffFlavs*sizeof(ferm_param)); //altrimenti nullo.
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


    // check if all parameter groups were found (neglecting optional goups) //qui controlla che tutti i parametri siano stati trovati.
    int check = 1;
    if(!helpmode){
        for(int igrouptype  = 0 ; igrouptype < NPMGTYPES; igrouptype++)
            if(igrouptype != PMG_FERMION  && igrouptype != PMG_DEBUG 
             && igrouptype != PMG_TESTS )  check *= tagcounts[igrouptype];
        if(!check){
            for(int igrouptype  = 0 ; igrouptype < NPMGTYPES; igrouptype++)
                if (!tagcounts[igrouptype]) if(0==devinfo.myrank)
                    printf("\"%s\"  parameter group not found!\n",
                            par_macro_groups_names[igrouptype]);
            return 1;
        }
    }


    int fermion_count = 0;
    // note 'check' is reused here
    int totcheck = 0;
    
    //qui inizia un grande for dove si analizzano tutti i gruppi di variabili.
    for(int igroup  = 0 ; igroup < found_tags; igroup++){
        int startline = tagpositions[igroup]; //la startline ovviamente è data dalle posizioni dei gruppi. //ovvero dalla riga dove quel gruppo si trova. A ogni iterazione la startline cambia.
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
            case PMG_TMEAS     : 
                check = read_topomeas_info(&meastopo_params,
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
            case PMG_INVERTER_TRICKS   : 
                check = read_inv_tricks_info(&inverter_tricks, 
                        filelines,startline,endline);
                break;
            case PMG_TESTS           :
                check = read_test_setting(&test_settings,
                        filelines,startline,endline);
                break;
            case PMG_REPLICAS       :
                check = read_replicas_numbers(rep,filelines,startline,endline) ;
              
                break;
            default:
                printf("TAG TYPE NOT RECOGNIZED\n");
                return 1;
                break;
        }
        if(check) //check diverso da 0.
            if(0==devinfo.myrank)
                printf("Problem in group %s\n", par_macro_groups_names[tagtypes[igroup]]);
        totcheck += check;
    }
    

    if(tagcounts[PMG_DEBUG]==0 ){
        if(!helpmode) read_debug_info(&debug_settings,filelines,0,1);// Just to set default values
        else read_debug_info(&debug_settings,filelines,0,0);// Just to set default values
    }


    alloc_info.stoutAllocations = act_params.stout_steps > 0 || act_params.topo_stout_steps > 0 || meastopo_params.stoutmeasstep > 0;
    if(devinfo.myrank == 0 )
        printf("Set alloc_info.stoutAllocations to %d\n",  alloc_info.stoutAllocations );


    // check == 1 means at least a parameter was not found.
    if(totcheck!=0 || helpmode ){
        
        if(helpmode && 0 == devinfo.myrank) fclose(helpfile);

        if(0==devinfo.myrank && totcheck )
            printf("There are errors in some groups, exiting.\n");
        printf("MPI%02d: Returning error 1...\n",devinfo.myrank);
        return 1;

    }
    else{

        uint32_t hash = hash_settings();
        char hash_string[32];
        char input_to_rename_filename[200];
        sprintf(hash_string,"%" PRIu32, hash);

        strcat(gauge_outfilename, hash_string);
        strcat(fm_par.fermionic_outfilename, hash_string);
        strcat(debug_settings.diagnostics_filename,hash_string);

        //if(0==devinfo.myrank){
        //    printf("Hash of all relevant settings: %s\n", hash_string );

        //    strcpy(input_to_rename_filename, input_filename);
        //    strcat(input_to_rename_filename, hash_string );
        //    
        //    FILE * input2 = fopen(input_filename,"r");
        //    FILE * input_renamed = fopen(input_to_rename_filename,"w"); 
        //    printf("Creating copy of run configuration file (%s) \nwith hash appended in the name (%s)...\n", input_filename, input_to_rename_filename);
        //    char ch;
        //    while(( ch = fgetc(input2) ) != EOF)
        //        fputc(ch, input_renamed);

        //    fclose(input2);
        //    fclose(input_renamed);

        //}


    }
   
  
    
  

    return 0;

}




#endif


