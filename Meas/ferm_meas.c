// use z2 noise instead of gaussian noise (see hep-lat/9308015)
// use the global defined fermions loc_chi, loc_phi, rnd_o, rnd_e, chi_o and loc_h

#ifndef FERM_MEAS_C 
#define FERM_MEAS_C

// vedi tesi F.Negro per ragguagli
void eo_inversion(su3_soa *tconf_acc,
		  double_soa * tbackfield,
		  ferm_param * tfermions_parameters,
                  double res,
		  vec3_soa * in_e,     // z2 noise
		  vec3_soa * in_o,     // z2 noise
		  vec3_soa * out_e,
		  vec3_soa * out_o,
		  vec3_soa * phi_e,    // parking variable
		  vec3_soa * phi_o,    // parking variable
                  vec3_soa * trialSolution,       // initial vector for the inversion 
		  vec3_soa * tloc_r,    // parking variable for the inverter      
		  vec3_soa * tloc_h,    // parking variable for the inverter
		  vec3_soa * tloc_s,    // parking variable for the inverter
                  vec3_soa * tloc_p){   // parking variable for the inverter




	    acc_Deo(tconf_acc, phi_e, in_o,tfermions_parameters,tbackfield);
	    combine_in1_x_fact1_minus_in2_back_into_in2(in_e, tfermions_parameters->ferm_mass , phi_e);
	    ker_invert_openacc(tconf_acc,tbackfield,tfermions_parameters,
				      out_e,phi_e,res,trialSolution,
				      tloc_r,tloc_h,tloc_s,tloc_p);
	    acc_Doe(tconf_acc, phi_o, out_e,tfermions_parameters,tbackfield);
	    combine_in1_minus_in2_allxfact(in_o,phi_o,(double)1/tfermions_parameters->ferm_mass,out_o);

    
}// end eo_inversion


d_complex chiral_condensate(vec3_soa * rnd_e,
	       vec3_soa * rnd_o,
	       vec3_soa * chi_e,
	       vec3_soa * chi_o){

  return scal_prod_global(rnd_o,chi_o) + scal_prod_global(rnd_e,chi_e);
     
}

void perform_chiral_measures(){


}




/*


void chiral_meas( su3_soa *tconf_acc,
		  double_soa * backfield,
		  ferm_param * tfermions_parameters,
		  vec3_soa * tkloc_r,///  forse ne servono anche altri vedere un po ... cmq parking variables
                  d_complex *chiral,
                  d_complex *e_dens,
                  d_complex *b_dens,
                  d_complex *p_dens) {
    

#ifdef IMCHEMPOT
  d_complex  eim=cos(tfermions_parameters->ferm_im_chem_pot)/((double)(nt))+I*sin(tfermions_parameters->ferm_im_chem_pot)/((double)(nt));
  d_complex  emim=conj(eim);
#endif
    
  // complex auxiliary vectors to be used in energy and baryon 
  // density when imaginary chemical potential is used
  dcomplex_soa *c_vector1,*c_vector2; // devono essere lunghi size = 2*sizeh!!!
  
    if(rand_vect>0)
      {
	int k, iter;
	long int i, j, index1, index2;
	complex<double> loc_chiral, loc_e_dens, loc_b_dens, loc_p_dens;
	complex<double> p1, p2;
	Vec3 v_e, v_o, w_e, w_o;
	
	loc_chiral=complex<double>(0.0,0.0);
	loc_e_dens=complex<double>(0.0,0.0);
	loc_b_dens=complex<double>(0.0,0.0);
	loc_p_dens=complex<double>(0.0,0.0);
	
	for(iter=0; iter<rand_vect; iter++)
	  {
	    rnd_e->z2noise();
	    rnd_o->z2noise();
	    
#ifdef MAGN
	    Deo(phi_e, rnd_o,no_ps);  // gli passo la carica; l'ultima componente del vettore di cariche che è nulla---> da qualche parte i link sono moltiplicati per 1
#else
	    Deo(phi_e, rnd_o);
#endif   //  ----> in realtà quello che viene fatto nei due casi è del tutto equivalente , controllare!!!!
	    for(i=0; i<sizeh; i++)
	      {
		(phi_e->fermion[i])=mass*(rnd_e->fermion[i])-(phi_e->fermion[i]);
	      }
	    
#ifndef USE_GPU
	    invert(chi_e, phi_e, residue_metro,no_ps);  // inverto la matrice mettendo come carica quella nulla --> passo come indice del campo l'ultimo (no_ps) che ha tutti 1
#else
	    smartpack_fermion_d(simple_fermion_packed, phi_e);
	    cuda_meas_init();
	    
	    int cg, order;
	    double *shifts = new double[max_approx_order];
	    int psferm=0;
	    get_order(order, *meas_inv_coeff);  // order=1
	    get_shifts(shifts, *meas_inv_coeff); // just one shift equal to 0.0
	    
	    cuda_shifted_inverter_d(residue_metro, shifts, order, psferm, &cg, no_ps); //anche qui inverto mettendo come carica Q=0 (quindi passo no_ps)
	    if(cg==max_cg)
	      {
		ofstream err_file;
		err_file.open(QUOTEME(ERROR_FILE), ios::app);
		err_file  << "WARNING: maximum number of iterations reached in cuda_multips_shift_invert in chiralmeas\n";
		err_file.close();
	      }
	    delete [] shifts;
	    
	    cuda_meas_end();
	    smartunpack_fermion_d(chi_e, simple_fermion_packed);
#endif
	    
#ifdef MAGN
	    Doe(phi_o, chi_e, no_ps); 
#else
	    Doe(phi_o, chi_e);
#endif
	    for(i=0; i<sizeh; i++)
	      { 
		chi_o->fermion[i]=one_by_mass*(rnd_o->fermion[i] - phi_o->fermion[i]);
	      } 
	    
	    // chiral condensate calculation
	    for(i=0; i<sizeh; i++)
	      { 
		c_vector1[i]=c_scalprod(rnd_o->fermion[i], chi_o->fermion[i]); 
		c_vector1[i]+=c_scalprod(rnd_e->fermion[i], chi_e->fermion[i]);
		//		if(i%200==0) cout << c_vector1[i] << endl; ----------------uncomment to check---------------------------
	      }
	    global_sum(c_vector1,sizeh);
	    loc_chiral+=c_vector1[0]*complex<double>(no_flavours*0.25*inv_size,0.0); 

	    
#ifndef IM_CHEM_POT
	    // energy and baryon density calculation
	    for(i=0; i<sizeh; i++)   // i=even index
	      {
		j=i+sizeh;            // j=odd
		index1=i+size3;
		index2=nnp[i][3]-sizeh;    // nnp[even,3]=odd   sizeh<=odd<size
		
		v_e=(gauge_conf->u_work[index1])*(chi_o->fermion[index2]);
		w_e=(gauge_conf->u_work[index1])*(rnd_o->fermion[index2]);
		
		p1=c_scalprod(rnd_e->fermion[i], v_e);
		p2=c_scalprod(w_e, chi_e->fermion[i]);
		
		c_vector1[i]=p1-p2;
		c_vector2[i]=p1+p2;
		
		index1=j+size3;
		index2=nnp[j][3];
		
		v_o=(gauge_conf->u_work[index1])*(chi_e->fermion[index2]);
		w_o=(gauge_conf->u_work[index1])*(rnd_e->fermion[index2]);
		
		p1=c_scalprod(rnd_o->fermion[i], v_o);
		p2=c_scalprod(w_o, chi_o->fermion[i]);
		
		c_vector1[i]+=p1-p2;
		c_vector2[i]+=p1+p2;
	      }
	    global_sum(c_vector1,sizeh);
	    loc_e_dens+=c_vector1[0]*complex<double>(0.5*no_flavours*0.25*inv_size,0.0);  // energy density
	    global_sum(c_vector2,sizeh);
	    loc_b_dens+=c_vector2[0]*complex<double>(0.5*no_flavours*0.25*inv_size,0.0);  // barion density
#else
	    // energy and baryon density calculation
	    for(i=0; i<sizeh; i++)   // i=even index
	      {
		j=i+sizeh;            // j=odd
		index1=i+size3;
		index2=nnp[i][3]-sizeh;    // nnp[even,3]=odd   sizeh<=odd<size
		
		v_e=(gauge_conf->u_work[index1])*(chi_o->fermion[index2]);
		w_e=(gauge_conf->u_work[index1])*(rnd_o->fermion[index2]);
		
		c_vector1[i]=c_scalprod(rnd_e->fermion[i], v_e);
		c_vector2[i]=c_scalprod(w_e, chi_e->fermion[i]);
		
		index1=j+size3;
		index2=nnp[j][3];
		
		v_o=(gauge_conf->u_work[index1])*(chi_e->fermion[index2]);
		w_o=(gauge_conf->u_work[index1])*(rnd_e->fermion[index2]);
		
		c_vector1[i]+=c_scalprod(rnd_o->fermion[i], v_o);
		c_vector2[i]+=c_scalprod(w_o, chi_o->fermion[i]);
	      }
	    global_sum(c_vector1,sizeh);
	    global_sum(c_vector2,sizeh);
	    
	    loc_e_dens+=(eim*c_vector1[0]-emim*c_vector2[0])*complex<double>(0.5*no_flavours*0.25*inv_size,0.0);  // energy density
	    loc_b_dens+=(eim*c_vector1[0]+emim*c_vector2[0])*complex<double>(0.5*no_flavours*0.25*inv_size,0.0);  // barion density
#endif
	    
	    // pressure density calculation
	    for(i=0; i<sizeh; i++)    // i=even
	      {
		j=i+sizeh;             // j=odd
		c_vector1[i]=0.0;
		for(k=0; k<3; k++)
		  {
		    index1=i+k*size;
		    index2=nnp[i][k]-sizeh;
		    
		    v_e=(gauge_conf->u_work[index1])*(chi_o->fermion[index2]);
		    w_e=(gauge_conf->u_work[index1])*(rnd_o->fermion[index2]);
		    
		    c_vector1[i]+=c_scalprod(rnd_e->fermion[i], v_e);
		    c_vector1[i]-=c_scalprod(w_e, chi_e->fermion[i]);
		    
		    index1=j+k*size;
		    index2=nnp[j][k];
		    
		    v_o=(gauge_conf->u_work[index1])*(chi_e->fermion[index2]);
		    w_o=(gauge_conf->u_work[index1])*(rnd_e->fermion[index2]);
		    
		    c_vector1[i]+=c_scalprod(rnd_o->fermion[i], v_o);
		    c_vector1[i]-=c_scalprod(w_o, chi_o->fermion[i]);
		  }
	      }
	    global_sum(c_vector1,sizeh);
	    loc_p_dens+=c_vector1[0]*complex<double>(0.5*no_flavours*0.25*inv_size,0.0);     // pressure density
	  }
	
	p1=complex<double>(1.0/(double) rand_vect, 0.0);
	
	chiral=loc_chiral*p1;
	e_dens=loc_e_dens*p1;
	b_dens=loc_b_dens*p1;
	p_dens=loc_p_dens*p1;
      }
    
    delete [] c_vector1;
    delete [] c_vector2;
    
#ifdef DEBUG_MODE
    cout << "\tterminated chiral_meas"<<endl;
#endif
  }

////////////////////////////////////////////////////////////////////////////////
/////////// Ridefinizione di chiral meas per il caso B =/= 0 ///////////////////
////////////////////////////////////////////////////////////////////////////////

#ifdef MAGN
void chiral_meas_M(complex<REAL> *chiral, complex<REAL> *e_dens, complex<REAL> *b_dens, complex<REAL> *p_dens)
  {
  #ifdef DEBUG
  cout << "DEBUG: inside chiral_meas ..."<<endl;
  #endif

  #ifdef IM_CHEM_POT
  const complex<REAL> eim=complex<REAL>(eim_cos, eim_sin);
  const complex<REAL> emim=complex<REAL>(eim_cos, -eim_sin);
  #endif

  complex<double> *c_vector1, *c_vector2;
  c_vector1=new complex<double>[size];    // complex auxiliary vectors to be used in energy and baryon 
  c_vector2=new complex<double>[size];    // density when imaginary chemical potential is used

  if(rand_vect>0)
    {
    int k, iter;
    long int i, j, index1, index2;
    complex<double> loc_chiral[no_ps], loc_e_dens[no_ps], loc_b_dens[no_ps], loc_p_dens[no_ps];
    complex<double> p1, p2;
    Vec3 v_e, v_o, w_e, w_o;
    int pseudofermion;

    for(int ips=0;ips<no_ps;ips++){
    loc_chiral[ips]=complex<double>(0.0,0.0);
    loc_e_dens[ips]=complex<double>(0.0,0.0);
    loc_b_dens[ips]=complex<double>(0.0,0.0);
    loc_p_dens[ips]=complex<double>(0.0,0.0);
    chiral[ips]=complex<double>(0.0,0.0);
    e_dens[ips]=complex<double>(0.0,0.0);
    b_dens[ips]=complex<double>(0.0,0.0);
    p_dens[ips]=complex<double>(0.0,0.0);
    }
    
    for(iter=0; iter<rand_vect; iter++)
      {
	rnd_e->z2noise();
	rnd_o->z2noise();
	
	for(pseudofermion=0;pseudofermion<no_ps;pseudofermion++)
	  {
	    
	    
	    Deo(phi_e, rnd_o, pseudofermion);
	    for(i=0; i<sizeh; i++)
	      {
		(phi_e->fermion[i])=mass*(rnd_e->fermion[i])-(phi_e->fermion[i]);
	      }

#ifndef USE_GPU
	    invert(chi_e, phi_e, residue_metro, pseudofermion);
#else
	    smartpack_fermion_d(simple_fermion_packed, phi_e);
	    cuda_meas_init();
	    
	    int cg, order;
	    double *shifts = new double[max_approx_order];
	    int psferm=0;
	    get_order(order, *meas_inv_coeff);  // order=1
	    get_shifts(shifts, *meas_inv_coeff); // just one shift equal to 0.0
	    
	    cuda_shifted_inverter_d(residue_metro, shifts, order, psferm, &cg, pseudofermion); // psferm =0 --> serve x l'indicizzazione e l'offset, ps.._index --> serve x la carica!!!
	    if(cg==max_cg)
	      {
		ofstream err_file;
		err_file.open(QUOTEME(ERROR_FILE), ios::app);   
		err_file  << "WARNING: maximum number of iterations reached in cuda_multips_shift_invert in chiralmeas\n";
		err_file.close();
	      }
	    delete [] shifts;
	    
	    cuda_meas_end();
	    smartunpack_fermion_d(chi_e, simple_fermion_packed);
#endif
	    
	    Doe(phi_o, chi_e, pseudofermion);
	    for(i=0; i<sizeh; i++)
	      { 
		chi_o->fermion[i]=one_by_mass*(rnd_o->fermion[i] - phi_o->fermion[i]);
	      } 
	    
	    // chiral condensate calculation
	    for(i=0; i<sizeh; i++)
	      { 
		c_vector1[i]=c_scalprod(rnd_o->fermion[i], chi_o->fermion[i]); 
		c_vector1[i]+=c_scalprod(rnd_e->fermion[i], chi_e->fermion[i]);
		//		if(i%200==0) cout << c_vector1[i] << endl; ---------------uncomment to check------------------------------------------
	      }
	    global_sum(c_vector1,sizeh);
	    loc_chiral[pseudofermion] +=c_vector1[0]*complex<double>(no_flavours*0.25*inv_size,0.0); 
#ifndef IM_CHEM_POT
	    // energy and baryon density calculation
	    for(i=0; i<sizeh; i++)   // i=even index
	      {
		j=i+sizeh;            // j=odd
		index1=i+size3;
		index2=nnp[i][3]-sizeh;    // nnp[even,3]=odd   sizeh<=odd<size
		
		v_e=(gauge_conf->u_work[index1])*(chi_o->fermion[index2]);
		v_e *= b_field[pseudofermion][index1];
		w_e=(gauge_conf->u_work[index1])*(rnd_o->fermion[index2]);
		w_e *= b_field[pseudofermion][index1];
		
		p1=c_scalprod(rnd_e->fermion[i], v_e);
		p2=c_scalprod(w_e, chi_e->fermion[i]);
		
		c_vector1[i]=p1-p2;
		c_vector2[i]=p1+p2;
		
		index1=j+size3;
		index2=nnp[j][3];
		
		v_o=(gauge_conf->u_work[index1])*(chi_e->fermion[index2]);
		v_o *= b_field[pseudofermion][index1];
		w_o=(gauge_conf->u_work[index1])*(rnd_e->fermion[index2]);
		w_o *= b_field[pseudofermion][index1];
		
		p1=c_scalprod(rnd_o->fermion[i], v_o);
		p2=c_scalprod(w_o, chi_o->fermion[i]);
		
		c_vector1[i]+=p1-p2;
		c_vector2[i]+=p1+p2;
	      }
	    global_sum(c_vector1,sizeh);
	    loc_e_dens[pseudofermion]+=c_vector1[0]*complex<double>(0.5*no_flavours*0.25*inv_size,0.0);  // energy density
	    global_sum(c_vector2,sizeh);
	    loc_b_dens[pseudofermion]+=c_vector2[0]*complex<double>(0.5*no_flavours*0.25*inv_size,0.0);  // barion density
#else
	    // energy and baryon density calculation
	    for(i=0; i<sizeh; i++)   // i=even index
	      {
		j=i+sizeh;            // j=odd
		index1=i+size3;
		index2=nnp[i][3]-sizeh;    // nnp[even,3]=odd   sizeh<=odd<size
		
		v_e=(gauge_conf->u_work[index1])*(chi_o->fermion[index2]);
		v_e *= b_field[pseudofermion][index1];
		w_e=(gauge_conf->u_work[index1])*(rnd_o->fermion[index2]);
		w_e *= b_field[pseudofermion][index1];
		
		c_vector1[i]=c_scalprod(rnd_e->fermion[i], v_e);
		c_vector2[i]=c_scalprod(w_e, chi_e->fermion[i]);
		
		index1=j+size3;
		index2=nnp[j][3];
		
		v_o=(gauge_conf->u_work[index1])*(chi_e->fermion[index2]);
		v_o *= b_field[pseudofermion][index1];
		w_o=(gauge_conf->u_work[index1])*(rnd_e->fermion[index2]);
		w_o *= b_field[pseudofermion][index1];
		
		c_vector1[i]+=c_scalprod(rnd_o->fermion[i], v_o);
		c_vector2[i]+=c_scalprod(w_o, chi_o->fermion[i]);
	      }
	    global_sum(c_vector1,sizeh);
	    global_sum(c_vector2,sizeh);
	    
	    loc_e_dens[pseudofermion]+=(eim*c_vector1[0]-emim*c_vector2[0])*complex<double>(0.5*no_flavours*0.25*inv_size,0.0);  // energy density
	    loc_b_dens[pseudofermion]+=(eim*c_vector1[0]+emim*c_vector2[0])*complex<double>(0.5*no_flavours*0.25*inv_size,0.0);  // barion density  
#endif
	    
	    // pressure density calculation
	    for(i=0; i<sizeh; i++)    // i=even
	      {
		j=i+sizeh;             // j=odd
		c_vector1[i]=0.0;
		for(k=0; k<3; k++)
		  {
		    index1=i+k*size;
		    index2=nnp[i][k]-sizeh;
		    
		    v_e=(gauge_conf->u_work[index1])*(chi_o->fermion[index2]);
		    v_e *= b_field[pseudofermion][index1];
		    w_e=(gauge_conf->u_work[index1])*(rnd_o->fermion[index2]);
		    w_e *= b_field[pseudofermion][index1];
		    
		    c_vector1[i]+=c_scalprod(rnd_e->fermion[i], v_e);
		    c_vector1[i]-=c_scalprod(w_e, chi_e->fermion[i]);
		    
		    index1=j+k*size;
		    index2=nnp[j][k];
		    
		    v_o=(gauge_conf->u_work[index1])*(chi_e->fermion[index2]);
		    v_o *= b_field[pseudofermion][index1];
		    w_o=(gauge_conf->u_work[index1])*(rnd_e->fermion[index2]);
		    w_o *= b_field[pseudofermion][index1];
		    
		    c_vector1[i]+=c_scalprod(rnd_o->fermion[i], v_o);
		    c_vector1[i]-=c_scalprod(w_o, chi_o->fermion[i]);
		  }
	      }
	    global_sum(c_vector1,sizeh);
	    loc_p_dens[pseudofermion]+=c_vector1[0]*complex<double>(0.5*no_flavours*0.25*inv_size,0.0);     // pressure density
	  }
      }
    p1=complex<double>(1.0/(double) rand_vect, 0.0);
    
    for(int ips=0;ips<no_ps;ips++){
      chiral[ips]=loc_chiral[ips]*p1;
      e_dens[ips]=loc_e_dens[ips]*p1;
      b_dens[ips]=loc_b_dens[ips]*p1;
      p_dens[ips]=loc_p_dens[ips]*p1;
    }
    
    }
  
  delete [] c_vector1;
  delete [] c_vector2;
#ifdef DEBUG_MODE
  cout << "\tterminated chiral_meas"<<endl;
#endif
  }
#endif



void calc_Q(REAL &tch0,REAL &tch1,REAL &tch2,REAL &tch3,REAL &tch4){
  tch0=0.0;
  tch1=0.0;
  tch2=0.0;
  tch3=0.0;
  tch4=0.0;
#ifdef TOPOLOGIC
  gauge_conf->save();                //salva la configurazione per futuro utilizzo
  long double* carica_top=new long double[size]; //sulla configurazione senza cooling
  gauge_conf->Q_topologica(carica_top);
  tch0=(double)carica_top[0];
  for(int k=0; k<30; k++)  gauge_conf->cooling();
  gauge_conf->Q_topologica(carica_top);
  tch1=(double)carica_top[0];
  for(int k=0; k<30; k++)  gauge_conf->cooling();
  gauge_conf->Q_topologica(carica_top);
  tch2=(double)carica_top[0];
  for(int k=0; k<30; k++)  gauge_conf->cooling();
  gauge_conf->Q_topologica(carica_top);
  tch3=(double)carica_top[0];
  for(int k=0; k<30; k++)  gauge_conf->cooling();
  gauge_conf->Q_topologica(carica_top);
  tch4=(double)carica_top[0];
  //tch4=tch3;
  delete[]carica_top;                                                                                                                                                  
  gauge_conf->copy_saved();                //riprende la config prima del coolig
#endif
}

void meas(ofstream &out,int what)
  {
  int i;
  REAL ps, pt, pr, pi;
  REAL pxyr,pxyi,pxzr,pxzi,pxtr,pxti,pyzr,pyzi,pytr,pyti,pztr,pzti;
  REAL tch0,tch1,tch2,tch3,tch4;

#ifdef MAGN
  complex<REAL> chiral_cond[no_ps], energy_density[no_ps], barion_density[no_ps], pressure_density[no_ps];
  complex<REAL> chiral_cond_dyn, energy_density_dyn, barion_density_dyn, pressure_density_dyn;
#endif

#ifndef MAGN
  complex<REAL> chiral_cond, energy_density, barion_density, pressure_density;
#endif


  if(what==0){
    gauge_conf->calc_plaq(ps,pt);
  }
  if(what==1){
    gauge_conf->calc_all_plaq(pxyr,pxyi,pxzr,pxzi,pxtr,pxti,pyzr,pyzi,pytr,pyti,pztr,pzti);
#ifdef TOPOLOGIC
  calc_Q(tch0,tch1,tch2,tch3,tch4);
#endif
  }
  gauge_conf->calc_poly(pr,pi);


  if(what==1){
    out << update_iteration << " "<< pxyr << "  " << pxyi << "  " << pxzr << "  " << pxzi << "  " << pxtr << "  " << pxti << "  " << pyzr << "  " << pyzi << "  " << pytr << "  " << pyti << "  " << pztr << "  " << pzti << "  " << pr << "  " << pi;
#ifndef TOPOLOGIC
    out << endl;
#endif

#ifdef TOPOLOGIC
    out << "   " << tch0 << "  " << tch1 << "  " << tch2 << "  " << tch3 << "  " << tch4  << endl;
#endif
  }

  if(what==0){
  for(i=0; i<4; i++)
    {
      out << update_iteration << "  ";
      //  plaqs        plaqt
      out << ps << "  " << pt << "  ";
      
      //  poly_re       poly_im 
      out << pr << "  " << pi << "  ";
      //     out << endl;
#ifndef MAGN
	chiral_meas(chiral_cond, energy_density, barion_density, pressure_density);
	out << real(chiral_cond) << "  " << imag(chiral_cond) << "  ";
	out << real(energy_density) << "  " << imag(energy_density) << "  ";
	out << real(barion_density) << "  " << imag(barion_density) << "  ";
	out << real(pressure_density) << "  " << imag(pressure_density) << "   ";
#endif
	
	
#ifdef MAGN
	// Misura delle osservabili totali per tutti i quarks
	//     cout << endl << "Misura del condensato con campo magnetico" << endl << endl;
	
	chiral_meas_M(chiral_cond, energy_density, barion_density, pressure_density);
	for(int ips=0;ips<no_ps;ips++){
	  out << real(chiral_cond[ips]) << "  " << imag(chiral_cond[ips]) << "  ";
	  out << real(energy_density[ips]) << "  " << imag(energy_density[ips]) << "  ";
	  out << real(barion_density[ips]) << "  " << imag(barion_density[ips]) << "  ";
	  out << real(pressure_density[ips]) << "  " << imag(pressure_density[ips]) << "  ";
	}
	//Misura del contributo dinamico delle osservabili (q=0)
	//     cout << endl << "Misura del condensato senza campo magnetico" << endl << endl;
	
	chiral_meas(chiral_cond_dyn, energy_density_dyn, barion_density_dyn, pressure_density_dyn);
	out << real(chiral_cond_dyn) << "  " << imag(chiral_cond_dyn) << "  ";
	out << real(energy_density_dyn) << "  " << imag(energy_density_dyn) << "  ";
	out << real(barion_density_dyn) << "  " << imag(barion_density_dyn) << "  ";
	out << real(pressure_density_dyn) << "  " << imag(pressure_density_dyn) << "   ";
#endif

    out << endl;
    }
  }
  }
*/

#endif
