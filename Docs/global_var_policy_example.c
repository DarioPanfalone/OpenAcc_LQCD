/*************************************************
 * This file contains a description              *
 * of the policies about the use of              *
 * global parking variables in various functions *
 ************************************************/

// EXAMPLE


typedef struct global_parking_t{

  int in_use;

  // stuff

} global_parking;


global_parking * global_thingy1, global_thingy2;

//allocation
//...


void b(global_parking * tglobal_thing){// will be probably global_thingy2


     if(tglobal_thing->in_use) printf("ARGH111!!!\n");
     
     global_thing->in_use = 1 ; 

     //do stuff
     //..

     // so long and thanks for all the fish
     global_thing->in_use = 0 ; 
     

}


void a(tglobal_thing * tglobal_thing){ // will be probably global_thingy1


     if(tglobal_thing->in_use) printf("ARGH111!!!\n");
     
     global_thing->in_use = 1 ; 

     //do stuff
     //..
     // here it is: the call of b()!!

     b(global_thingy2);

     
    

     // so long and thanks for all the fish
     global_thing->in_use = 0 ; 
 


     
}



