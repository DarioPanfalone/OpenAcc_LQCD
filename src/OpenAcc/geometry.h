#ifndef GEOMETRY_H_
#define GEOMETRY_H_


// if using GCC, there are some problems with __restrict.
#ifdef __GNUC__
 #define __restrict
#endif

// lattice dimensions
#define nd0 48  
#define nd1 48  
#define nd2 48  
#define nd3 16

#define nd0h (nd0 >> 1) // nx/2
#define nd1h (nd1 >> 1)
#define nd2h (nd2 >> 1)
#define nd3h (nd3 >> 1)

// for legacy - to remove after refactoring is complete
#define nx nd0    
#define ny nd1
#define nz nd2
#define nt nd3

#define nxh (nx >> 1) // nx/2
#define nyh (ny >> 1)
#define nzh (nz >> 1)
#define nth (nt >> 1)



#define sizehh nd0*nd1*nd2*nd3/2 

#define ANTIPERIODIC_T_BC  // else periodic time bc are taken

#define vol1 nd0
#define vol2 (nd1 * vol1)
#define vol3 (nd2 * vol2)
#define vol4 (nd3 * vol3)
#define size vol4
#define size2 (2*size)
#define size3 (3*size)
#define sizeh (size / 2)
#define no_links (4 * vol4)


typedef struct geom_parameters_t{

    int gnx,gny,gnz,gnt;
    // map of physical directions onto logical directions
    int xmap,ymap,zmap,tmap;// tmap is the "antiperiodic direction" for 
                            // fermions
    

} geom_parameters;

extern geom_parameters geom_par;


#pragma acc routine seq
static inline int snum_acc(int d0, int d1, int d2, int d3) {
  int ris;
  ris = d0 + (d1*vol1) + (d2*vol2) + (d3*vol3);
  return ris/2;   // <---  /2 Pay attention to even/odd 
}

// TABLES FOR THE NEAREST NEIGHBOURS       
// nnp[site_half][dir][par] = nearest neighbour in the positive direction "dir"            
//                            starting from the site "site_half" (from 0 to sizeh) of parity "par"         
// nnm[site_half][dir][par] = nearest neighbour in the negative direction "dir"            
//                            starting from the site "site_half" (from 0 to sizeh) of parity "par"         
// temporarily defined and computed here (should be moved elsewhere!)      
int nnp_openacc[sizeh][4][2];
int nnm_openacc[sizeh][4][2];



void compute_nnp_and_nnm_openacc(void);

#endif
