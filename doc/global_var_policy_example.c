/*************************************************
 * This file contains a description              *
 * of the policies about the use of              *
 * global parking variables in various functions *
 ************************************************/

// EXAMPLE

#define FREE 0 
#define REQUESTED 1
#define IN_USE 2
#define PRINT(x) (#x)
#define SETINUSE(var)  if(var->status & IN_USE ) { printf("%s is in use!\n",PRINT(var)); } else {printf("%s is not in use, using it!\n", PRINT(var)); var->status = var->status | IN_USE;}
#define SETREQUESTED(var)  if(var->status != FREE ) { printf("%s is in use or requested!(%d)\n",PRINT(var),var->status); }else {printf("%s is not in use or requested, requiring it!\n", PRINT(var)); var->status = var->status | REQUESTED;}
#define SETFREE(var) {printf("freeing %s\n",PRINT(var)); var->status = FREE;}


typedef struct global_parking_t{

  int status;

  // stuff

} global_parking;


global_parking * global_thingy1, global_thingy2;

//allocation
//...


void b(global_parking * tglobal_thing){// will be probably global_thingy2


//    if(tglobal_thing->in_use) printf("ARGH111!!!\n");
//    tglobal_thing-> status = 1 ; 
     SETINUSE(tglobal_thing);


     //do stuff
     //..

     // so long and thanks for all the fish
     tglobal_thing->in_use = 0 ; 
     
     
}


void a(tglobal_thing * tglobal_thing){ // will be probably global_thingy1


//   if(tglobal_thing->in_use) printf("ARGH111!!!\n");
//   tglobal_thing->in_use = 1 ; 
     SETINUSE(tglobal_thing);

     //do stuff
     //..
     // here it is: the call of b()!!

     // perhaps...
     SETREQUESTED(global_thingy2);
     b(global_thingy2);

     
    

     // so long and thanks for all the fish

    SETFREE(tglobal_thing);
}



