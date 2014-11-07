#ifndef INVERTER_CC_ 
#define INVERTER_CC_ 

#include <time.h> 
#include <fstream>
#include "../FermionMatrix/fermionmatrix.cc"
#include "../Include/global_const.cc"
#include "../Include/parameters.cc"
#include "../Include/global_var.cc"


#include "../Global_sum/global_sum.cc"


using namespace std;

// Calculate "out" from  (M^dag M) * out = in
// using the CG-algorithm
// (M^dag M) = mass2 - Deo*Doe
// Using: loc_r, loc_h, loc_s, loc_p   

// Doe e Deo sono definite in '../FermionMatrix/fermionmatrix.cc'
// Doe 

void invert (Fermion *out, const Fermion *in, REAL res, const Fermion *trialSolution = NULL ){
  #if ((defined DEBUG_MODE) || (defined DEBUG_INVERTER))
  cerr <<  "DEBUG: inside invert ..."<< endl;
  cerr.flush();
  #endif

  clock_t time_start, time_finish;
  time_start=clock();

  static Fermion vloc_r;
  static Fermion vloc_h;
  static Fermion vloc_s;
  static Fermion vloc_p;  
 
  /****************************************************************
  I puntatori di seguito sono stati definiti per pura pigrizia, 
  altrimenti avrei dovuto cambiare tutti i '->' in '.' nel seguito
  ****************************************************************/
  Fermion* loc_r= &vloc_r;
  Fermion* loc_h= &vloc_h;
  Fermion* loc_s= &vloc_s;
  Fermion* loc_p= &vloc_p;  


  int cg;
  long int i;
  double delta, alpha, lambda, omega, gammag;
  Vec3 vr_1, vr_2, vr_3, vr_4, vr_5;

  // initializations
  if (trialSolution == NULL) 
  out->gauss();      // trial solution for out
  else (*out) = Fermion(*trialSolution);

  // s=(M^dagM)out
  Doe(loc_h, out);
  Deo(loc_s, loc_h);


  for(i=0; i<sizeh; i++)
     {
     vr_1=(out->fermion[i]);
     vr_2=(loc_s->fermion[i]);
     loc_s->fermion[i]=mass2*vr_1-vr_2;
     }

  // r=in-s, p=r, delta=(r,r)
  for(i=0; i<sizeh; i++)
     {
     vr_3=(in->fermion[i])-(loc_s->fermion[i]);

     loc_r->fermion[i]=vr_3;
     loc_p->fermion[i]=vr_3;
     d_vector1[i]=vr_3.l2norm2();
     }
  global_sum(d_vector1, sizeh);
  delta=d_vector1[0];

  clock_t start, ttm;
  start = clock();

  // loop over cg iterations
  cg=0;
  do {
     cg++;
    
     // s=(M^dag M)p    alpha=(p,s)
     Doe(loc_h, loc_p);
     Deo(loc_s, loc_h);

     for(i=0; i<sizeh; i++)
        {
        vr_1=(loc_p->fermion[i]);
        vr_2=(loc_s->fermion[i]);
        vr_3=mass2*vr_1 - vr_2;
        (loc_s->fermion[i])=vr_3;
        d_vector1[i]=r_scalprod(vr_1,vr_3);
        }
     global_sum(d_vector1, sizeh);
     alpha=d_vector1[0];

     omega=delta/alpha;
     
     // out+=omega*p  r-=omega*s
     // lambda=(r,r);
     for(i=0; i<sizeh; i++)
        {
        vr_1=(out->fermion[i]);
        vr_2=(loc_p->fermion[i]);
        vr_3=vr_1 + omega*vr_2;
        (out->fermion[i])=vr_3;

        vr_3=(loc_r->fermion[i]);
        vr_4=(loc_s->fermion[i]);
        vr_5=vr_3 - omega*vr_4;
        (loc_r->fermion[i])=vr_5;

        d_vector1[i]=vr_5.l2norm2();
        }
     global_sum(d_vector1, sizeh);
     lambda=d_vector1[0];

     gammag=lambda/delta;
     delta=lambda;

     // p=r+gammag*p
     for(i=0; i<sizeh; i++)
        {
        vr_1=(loc_r->fermion[i]);
        vr_2=(loc_p->fermion[i]);
        vr_3=vr_1+gammag*vr_2;
        (loc_p->fermion[i])=vr_3;
        }
    cerr << "Inversion iter: " << cg << "  residuo " << sqrt(lambda) << " (target=" << res << ")  " << endl;    

     } while( (sqrt(lambda)>res) && cg<max_cg);
  
  //  #if ((defined DEBUG_MODE) || (defined DEBUG_INVERTER))
  cout << "\tterminated invert terminated in "<<cg<<" iterations [";
  // test
  Doe(loc_h, out);
  Deo(loc_s, loc_h);


  for(i=0; i<sizeh; i++)
     {
     vr_1=(out->fermion[i]);
     vr_2=(loc_s->fermion[i]);
     vr_3=(in->fermion[i]);
     vr_4=mass2*vr_1 -vr_2 -vr_3;
     d_vector1[i]=vr_4.l2norm2();
     }
  global_sum(d_vector1, sizeh);
  cout << " res/stop_res="<< sqrt(d_vector1[0])/res << " ,  stop_res="<< res << " ]"<<endl;
  //  #endif

  if(cg==max_cg)
    {
    ofstream err_file;
    err_file.open(QUOTEME(ERROR_FILE), ios::app);   
    err_file << "WARNING: maximum number of iterations reached in invert\n";
    err_file.close();
    }

  time_finish=clock();
  cout << "CPU INVERSION times:        Tot time: " << ((REAL)(time_finish)-(REAL)(time_start))/CLOCKS_PER_SEC << " sec     AvgTime/cg_iter: " << ((REAL)(time_finish)-(REAL)(time_start))/CLOCKS_PER_SEC/cg << endl << endl ;



  }


#endif
