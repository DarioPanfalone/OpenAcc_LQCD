
#ifndef DEBUG_MACROS_GLVARCHECK
#define DEBUG_MACROS_GLVARCHECK
// machinery to check for statuses/ use conflicts

#define FREE 0 
#define REQUESTED 1
#define IN_USE 2
#define PRINT(x) (#x)

// (i1 & i2 ) e' un operazione logica & bit a bit sui due interi i1 ed i2
// esempio: i1 = 0000111001
//          i2 = 0010010011
//   risultato = 0000010001

#define SETINUSE(var)\
if(var->status & IN_USE )\
  { printf("SETINUSE %s:%d, var %s is in use! (%p)\n",__FILE__,__LINE__,PRINT(var)),var; }\
  else {printf("SETINUSE %s:%d, var %s is not in use, using it! (%p)\n",\
	       __FILE__,__LINE__, PRINT(var),var);				\
      var->status = var->status | IN_USE;}

#define SETREQUESTED(var)\
{if(var->status == REQUESTED )	   /* richiesta ma non in uso */\
 printf("SETREQUESTED %s:%d, var %s is already requested (%p).\n",\
        __FILE__,__LINE__,PRINT(var),var);\
else if(var->status & IN_USE)    /* in uso (sia gia' richiesta che no)*/ \
 printf("SETREQUESTED %s:%d, var %s is in use! (%p)\n",\
        __FILE__,__LINE__,PRINT(var),var); /* ne' richiesta ne' in uso*/ \
else printf("SETREQUESTED %s:%d, var %s is not in use or requested, requiring it!(%p)\n"\
        ,__FILE__,__LINE__, PRINT(var),var);\
 var->status = REQUESTED;}

#define SETFREE(var)\
  {printf("SETFREE %s:%d, freeing %s (%p)  [previous status = %d] \n",__FILE__,__LINE__,PRINT(var),var,var->status); var->status = FREE;}


#define CHECKSTATUS(var)\
  {printf("CHECKSTATUS %s:%d, var %s (%p) has status %d.\n",__FILE__,__LINE__,PRINT(var),var,var->status);}

#endif

