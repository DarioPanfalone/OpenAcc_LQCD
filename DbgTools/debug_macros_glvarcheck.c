
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


#define SETINUSE(var)

#define SETREQUESTED(var)

#define SETFREE(var)

#define CHECKSTATUS(var)




#endif

