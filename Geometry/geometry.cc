#ifndef GEOMETRY_CC_
#define GEOMETRY_CC_

#include "../Include/global_const.cc"
#include "../Include/global_var.cc"
#include "../Include/parameters.cc"


// superindex

long int snum(int x, int y, int z, int t)
  {
  long int ris;

  ris=x+y*vol1+z*vol2+t*vol3;
  return ris/2;   // <---  /2 Pay attention to even/odd  (see init_geo) 
  }

int srotola(int x,int y,int z, int t) { //evolution of snum, gives directly the index (TO TEST)

	int xi[4];
	xi[0] = x;	xi[1] = y;	xi[2] = z;	xi[3] = t;

	int xiRestr[4];

	xiRestr[0] = xi[0]%nx;	if(xiRestr[0] < 0) xiRestr[0] += nx;
	xiRestr[1] = xi[1]%ny;	if(xiRestr[1] < 0) xiRestr[1] += ny;
	xiRestr[2] = xi[2]%nz;	if(xiRestr[2] < 0) xiRestr[2] += nz;
	xiRestr[3] = xi[3]%nt;	if(xiRestr[3] < 0) xiRestr[3] += nt;


	return (int) (xiRestr[0] +
			xiRestr[1] * vol1 +
			xiRestr[2] * vol2 +
			xiRestr[3] * vol3 )/2  +
			((xiRestr[0]+xiRestr[1]+xiRestr[2]+xiRestr[3]) % 2 )* sizeh ;

}



//   links used according to this scheme
//
//   0            size         size2         size3         no_links
//   |------|------|------|------|------|------|------|------|
//      e      o      e       o     e      o      e      o
//        x-dir         y-dir         z-dir         t-dir

// initialize geometry
// periodic spatial bc are always assumed



void init_geo(void)
  {
  #ifdef DEBUG_MODE
  cout << "DEBUG: inside init_geo ..."<<endl;
  #endif

  int even;
  int x, y, z, t, xm, ym, zm, tm, xp, yp, zp, tp;
  long int num;
  int sum, rest;
  
  // staggered phases
  eta = new int[no_links];

  // allocate nnp[size][4]
  nnp = new long int * [size];
  for(num=0; num<size; num++)
     {
     nnp[num] = new long int [4];
     }

  // allocate nnm[size][4]
  nnm = new long int * [size];
  for(num=0; num<size; num++)
     {
     nnm[num] = new long int [4];
     }

  for(t=0; t<nt; t++)
     {
     tp=t+1;
     tm=t-1;
     if(t==nt-1) tp=0;
     if(t==0) tm=nt-1;

     for(z=0; z<nz; z++)
        {
        zp=z+1;
        zm=z-1;
        if(z==nz-1) zp=0;
        if(z==0) zm=nz-1;

        for(y=0; y<ny; y++)
           {
           yp=y+1;
           ym=y-1;
           if(y==ny-1) yp=0;
           if(y==0) ym=ny-1;

           for(x=0; x<nx; x++)
              {
              xp=x+1;
              xm=x-1;
              if(x==nx-1) xp=0;
              if(x==0) xm=nx-1;

              // the even sites get the lower half (     0 --> sizeh-1 )
              // the odd  sites get the upper half ( sizeh --> size -1 )
              
              sum = x+y+z+t;
              even = sum%2;        // even=0 for even sites
                                   // even=1 for odd sites

              num = even*sizeh + snum(x,y,z,t);

              // NEXT-NEIGHBOURS DEFINITION

              // x-dir             
              nnp[num][0]=(1-even)*sizeh + snum(xp,y,z,t);
              nnm[num][0]=(1-even)*sizeh + snum(xm,y,z,t);
             
              //y-dir
              nnp[num][1]=(1-even)*sizeh + snum(x,yp,z,t);
              nnm[num][1]=(1-even)*sizeh + snum(x,ym,z,t);

              //z-dir
              nnp[num][2]=(1-even)*sizeh + snum(x,y,zp,t);
              nnm[num][2]=(1-even)*sizeh + snum(x,y,zm,t);
 
              //t-dir
              nnp[num][3]=(1-even)*sizeh + snum(x,y,z,tp);
              nnm[num][3]=(1-even)*sizeh + snum(x,y,z,tm);


              // ETA definitions
              
              // x-dir
              eta[num]=1;

              // y-dir
              sum=x;
              rest=sum%2;
              if(rest==0)  
                {
                eta[num+size]=1;
                }
              else
                {
                eta[num+size]=-1;
                }

              // z-dir
              sum=x+y;
              rest=sum%2;
              if(rest==0)  
                {
                eta[num+size2]=1;
                }
              else
                {
                eta[num+size2]=-1;
                }

              // t-dir
              sum=x+y+z;
              rest=sum%2;
              if(rest==0)  
                {
                eta[num+size3]=1;
                }
              else
                {
                eta[num+size3]=-1;
                }
              if(ferm_temp_bc==0)
                {
                if(t==nt-1)        // antiperiodic temporal b.c. for fermions
                  {
                  eta[num+size3]=-eta[num+size3];
                  }
                }
              }
           } 
        }
     }
  #ifdef DEBUG_MODE
  cout << "\tterminated init_geo"<<endl;
  #endif
  }



// clear geometry
void end_geo(void)
  {
  #ifdef DEBUG_MODE
  cout << "DEBUG: inside end_geo ..."<<endl;
  #endif

  long int num;

  for(num=0; num<size; num++)
     {
     delete [] nnm[num];
     }
  delete [] nnm;

  for(num=0; num<size; num++)
     {
     delete [] nnp[num];
     }
  delete [] nnp;

   delete [] eta;

  #ifdef DEBUG_MODE
  cout << "\tterminated end_geo"<<endl;
  #endif
  }

void coord(long int index,int &x,int &y,int &z,int &t){    //RESTITUISCE IL VALORE DELLE COORDINATE SPAZIO-TEMPORALI ASSOCIATE AL MEGA INDICE index                
                                                           // index VA DA 0 FINO A size-1                                                                          
  int even;
  long int temp_index;

  temp_index = index;
  even = 1;
  if(temp_index < sizeh) even = 0;
  temp_index = temp_index - even * sizeh;
  temp_index = temp_index * 2;
  x = (((temp_index % vol3) % vol2) % vol1);
  y = (((temp_index % vol3) % vol2)-x) / vol1;
  z = ((temp_index % vol3) - x - y * vol1) / vol2;
  t = (temp_index -x - y * vol1 - z * vol2) / vol3;
  if(((x+y+z+t)%2)!=even) x = x + 1;
}


#endif
