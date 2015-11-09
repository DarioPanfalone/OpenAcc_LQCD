
#ifndef DEBUG_MACROS_GLVARCHECK
#define DEBUG_MACROS_GLVARCHECK
// machinery to check for statuses/ use conflicts

#define FREE 0 
#define REQUESTED 1
#define IN_USE 2
#define PRINT(x) (#x)
#define SETINUSE(var)  if(var->status & IN_USE ) { printf("%s:%s, var %s is in use!\n",__FILE__,__LINE__,PRINT(var)); } else {printf("%s:%s, var %s is not in use, using it!\n",__FILE__,__LINE__, PRINT(var)); var->status = var->status | IN_USE;}
#define SETREQUESTED(var)  if(var->status != FREE ) { printf("%s:%s, var %s is in use or requested!(status %d)\n",__FILE__,__LINE__,PRINT(var),var->status); }else {printf("%s:%s, var %s is not in use or requested, requiring it!\n",__FILE__,__LINE__, PRINT(var)); var->status = var->status | REQUESTED;}
#define SETFREE(var) {printf("%s:%s, freeing %s\n",__FILE__,__LINE__,PRINT(var)); var->status = FREE;}

#endif

