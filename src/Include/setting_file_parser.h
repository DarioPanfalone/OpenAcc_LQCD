#ifndef _SETTING_FILE_PARSER_H_
#define _SETTING_FILE_PARSER_H_

#define MAXLINES 300
#define MAXLINELENGTH 500 // pretty long to accomodate all the comments

extern char input_file_str[MAXLINES*MAXLINELENGTH];

int set_global_vars_and_fermions_from_input_file(const char* input_filename);

#endif



