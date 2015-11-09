
#ifndef DEBUG_MACROS_GLVARCHECK
#define DEBUG_MACROS_GLVARCHECK
// machinery to check for statuses/ use conflicts

#define FREE 0 
#define REQUESTED 1
#define IN_USE 2
#define PRINT(x) (#x)

#define SETINUSE(var)\
  if(var->status & IN_USE )\
  { printf("%s:%s, var %s is in use!(%p)\n",__FILE__,__LINE__,PRINT(var)),var; }\
  else {printf("%s:%s, var %s is not in use, using it!\n",\
          __FILE__,__LINE__, PRINT(var));\
      var->status = var->status | IN_USE;}

#define SETREQUESTED(var)\
  if(var->status == REQUESTED )\
 printf("%s:%s, var %s is already requested (%p).\n",\
        __FILE__,__LINE__,PRINT(var),var);\
else if(var->status | IN_USE)\
 printf("%s:%s, var %s is in use! (%p)\n",\
        __FILE__,__LINE__,PRINT(var),var);\
else printf("%s:%s, var %s is not in use or requested, requiring it!(%p)\n"\
        ,__FILE__,__LINE__, PRINT(var),var);\
    var->status = REQUESTED;}

#define SETFREE(var) {printf("%s:%s, freeing %s (%p)\n",__FILE__,__LINE__,PRINT(var),var); var->status = FREE;}

#endif

