class Staples {
  //private:
public:

  Su3 staples[no_links];                      // ho tante staples quanti link ---> ogni staple è la somma di 6 matrici di SU(3)
                                              // ognuna delle quali è il prodotto di tre matrici di SU(3) (quelle che formano [ )
  Staples(void);
 
  friend void calc_staples(int aggiungo);

  // defined in Ipdot/ipdot.cc
  /*
  friend void calc_ipdot_gauge(void);
  friend void calc_ipdot(void);
  */
};

// base constructor
Staples::Staples(void)                          // il costruttore base con argomento void (-> niente)
 {
 for(long int i=0; i<no_links; i++)
    {
    staples[i].zero();                                 // tutte le staples sono delle matrici piene di zeri
    }
 }


// calculate staples 
void calc_staples(int aggiungo)
 {
 #ifdef DEBUG_MODE
 cout <<"DEBUG: inside calc_staples ..."<<endl;
 #endif
 int mu, nu,rho,sigma,k;
 long int pos, index_mu, index_nu, helper;
 Su3 aux,aux1,aux2,aux3;
  
 for(pos=0; pos < size; pos++)
   {
     for(mu=0; mu<4; mu++)
       {
	 index_mu=mu*size;
	 
	 (gauge_staples->staples[index_mu+pos]).zero();
	 
	 for(nu=0; nu<4; nu++)
	   {
	     if(nu!=mu)
	       {
		 index_nu=nu*size;
		 
		 //             (1)
		 //           +-->--+
		 //           |     |
		 // ^mu       |     V  (2)   to calculate (1)*(2)*(3)
		 // |         |     |
		 // |         +--<--+
		 // -->nu   pos  (3)
		 
		 aux1 = (gauge_conf->u_work[index_nu + nnp[pos][mu]]);  // 1
		 aux2=~(gauge_conf->u_work[index_mu + nnp[pos][nu]]);  // 2 
		 aux3=~(gauge_conf->u_work[index_nu + pos]);           // 3
		 
		 (gauge_staples->staples[index_mu+pos])+=aux1*aux2*aux3;

		 //             (1)
		 //           +--<--+
		 //           |     |
		 // ^mu   (2) V     |   to calculate (1)*(2)*(3)
		 // |         |     |
		 // |         +-->--+
		 // -->nu  help  (3)  pos
		 
		 helper=nnm[pos][nu];
		 aux1=~(gauge_conf->u_work[index_nu + nnp[helper][mu]]); // 1
		 aux2=~(gauge_conf->u_work[index_mu + helper]);          // 2
		 aux3= (gauge_conf->u_work[index_nu + helper]);          // 3
		 
		 (gauge_staples->staples[index_mu+pos])+=aux1*aux2*aux3;
	       }
	   }
       }
   }
 

#ifdef DEBUG_MODE
 cout <<"\tterminated calc_staples"<<endl;
#endif
 }


