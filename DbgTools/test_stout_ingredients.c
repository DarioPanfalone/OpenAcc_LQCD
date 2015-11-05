

#include <complex.h>
#include "../OpenAcc/struct_c_def.c"
#include "../OpenAcc/single_types.c"
#include "../OpenAcc/cayley_hamilton.c"
#include "../OpenAcc/stouting_deottimizzato.c"

struct anti_hermitian_exp_ingredients
{
    int sign;
    double c0,c1;
    double cu,su;
    double c2u,s2u;
    single_su3 Q,Q2;
    double u,w,theta;
    double xi0w;
    double cw;
    d_complex f[3];
};



void anti_hermitian_exact_i_exponentiate_ingredients(anti_hermitian_exp_ingredients &out,single_tamat * QA)
  {
    //copy Q

    tamat_to_su3(&(out.Q),QA);

    double c0 = out.c0 = det_i_times_QA(QA); //(14) // cosi calcolo direttamente il det(Q)
    
    //takes the square of Q
    single_su3xsu3(&out.Q2,&out.Q,&out.Q);
    
    //takes 1/2 of the real part of the trace of Q2 (eq. 15)
    double c1=out.c1=0.5 * Tr_i_times_QA_sq(QA); // (15)
    //compute c0_max (eq. 17)
    double c0_max=2*pow(c1/3,1.5);
    
    //consider the case in which c1<4*10^-3 apart, as done in MILC
    if(c1<4e-3)
      {
	out.f[0]=     1-c0*c0/720;
	out.f[0]+=  -c0*(1-c1*(1-c1/42)/20)/6 * I;
	out.f[1]=c0*(1-c1*(1-3*c1/112)/15)/24;
	out.f[1]+=1-c1*(1-c1*(1-c1/42)/20)/6-c0*c0/5040*I;
	out.f[2]=0.5*(-1+c1*(1-c1*(1-c1/56)/30)/12+c0*c0/20160);
	out.f[2]+=0.5*(c0*(1-c1*(1-c1/48)/21)/60)*I;
      }
    else
      {
	//take c0 module and write separately its sign (see note after eq. 34)
	out.sign=0;
	if(c0<0)
	  {
	    out.sign=1;
	    c0=-c0;
	  }
	
	//check rounding error
	double eps=(c0_max-c0)/c0_max;
	
	//(eqs. 23-24)
	double theta;
	if(eps<0) theta=out.theta=0; //only possible as an effect of rounding error when c0/c0_max=1
	else
	  if(eps<1e-3) theta=out.theta=sqrt(2*eps)*(1+(1.0/12+(3.0/160+(5.0/896+(35.0/18432+63.0/90112*eps)*eps)*eps)*eps)*eps);
	  else theta=out.theta=acos(c0/c0_max);
	double u=out.u=sqrt(c1/3)*cos(theta/3);
	double w=out.w=sqrt(c1)*sin(theta/3);
	
	//auxiliary variables for the computation of h0, h1, h2
	double u2=u*u,w2=w*w,u2mw2=u2-w2,w2p3u2=w2+3*u2,w2m3u2=w2-3*u2;
	double cu=out.cu=cos(u),c2u=out.c2u=cos(2*u);
	double su=out.su=sin(u),s2u=out.s2u=sin(2*u);
	double cw=out.cw=cos(w);
	
	//compute xi function defined after (eq. 33)
	double xi0w;
	if(fabs(w)<0.05)
	  {
	    double temp0=w*w,temp1=1-temp0/42,temp2=1.0-temp0/20*temp1;
	    xi0w=1-temp0/6*temp2;
	  }
	else xi0w=sin(w)/w;
	out.xi0w=xi0w;
	
	//computation of h0, h1, h2 (eqs. 30-32)
	d_complex h0=
	  u2mw2*c2u+ //(u2-w2)*cos(2u)
	  cu*8*u2*cw+ //cos(u)*8*u2*cos(w)
	  2*su*u*w2p3u2*xi0w+ //sin(u)*2*mu*(3*u2+w2)*xi0(w)
	  I*(u2mw2*s2u+ //(u2-w2)*sin(2u)
	  -su*8*u2*cw+ //-sin(u)*8*u2*cos(w)
	  cu*2*u*w2p3u2*xi0w); //cos(u)*2*u*(3*u2+w2)*xi0(w)
	d_complex h1=
	  2*u*c2u+ //2*u*cos(2u)
	  -cu*2*u*cw+ //cos(u)*2*u*cos(w)
	  -su*w2m3u2*xi0w + //sin(u)*(u2-3*w2)*xi0(w)
	  I * (2*u*s2u+ //2*u*sin(2u)
	  su*2*u*cos(w)+ //sin(u)*2*u*cos(w)
	  -cu*w2m3u2*xi0w);//cos(u)*(3*u2-w2)*xi0(w)
	d_complex h2=
	  c2u+ //cos(2u)
	  -cu*cw+ //-cos(u)*cos(w)
	  -3*su*u*xi0w +  //-3*sin(u)*u*xi0(w)
	  I * (s2u+ //sin(2u)
	  su*cw+ //sin(w)*cos(w)
	  -cu*3*u*xi0w)};//-cos(u)*3*u*xi0(w)
	
	//build f (eq. 29)
	double fact=1/(9*u*u-w*w);
	out.f[0] = h0*fact;
	out.f[1] = h1*fact;
	out.f[2] = h2*fact;
	
	//change sign to f according to (eq. 34)
	if(out.sign!=0)
	  {
	    //out.f[0][IM]*=-1;
	    out.f[0]=conj(out.f[0]);
	    //out.f[1][RE]*=-1;
	    out.f[1]=-conj(out.f[1]);
	    //out.f[2][IM]*=-1;
	    out.f[2]=-conj(out.f[2]);
	  }
      }

void stouted_force_compute_coeffs_NISSA(anti_hermitian_exp_ingredients *ing, single_tamat *QA)
  {
    d_complex b[2][3];
    
    if(ing->c1<4e-3)
      {
	double c0=ing->c0,c1=ing->c1;
	b[1][0]=-c0/360;
	b[1][0]+=(-1.0/6*(1-c1/20*(1-c1/42))) * I;

	b[0][0]=0;
	b[0][0]+=(c0/120*(1-c1/21))*I;

	b[1][1]=1.0/24*(1-c1/15*(1-3*c1/112));
	b[1][1]+=(-c0/2520)*I;

	b[0][1]=-c0/360*(1-3*c1/56);
	b[0][1]+=(-1.0/6*(1-c1/10*(1.0-c1/28)))*I;

	b[1][2]=0.5*c0/10080;
	b[1][2]+=(0.5*(1.0/60*(1-c1/21*(1-c1/48))))*I;

	b[0][2]=0.5*(1.0/12*(1-2*c1/30*(1-3*c1/112))); 
	b[0][2]+=(0.5*(-c0/1260*(1-c1/24)))*I;
      }
    else
      {
	//copy back the stored variables
	double u=ing->u;
	double w=ing->w;
	double xi0w=ing->xi0w;
	double cu=ing->cu,c2u=ing->c2u;
	double su=ing->su,s2u=ing->s2u;
	double cw=ing->cw;
	
	//compute additional variables
	double u2=u*u,w2=w*w;
	double xi1w; //eq. (67)
	if(fabs(w)<0.05) xi1w=-(1-w2*(1-w2*(1-w2/54)/28)/10)/3;
	else xi1w=cw/w2-sin(w)/(w2*w);
	
	//eqs. (60-65)
	d_complex r[2][3]=
	  {{2*c2u*u+s2u*(-2*u2+2*w2)+2*cu*u*(8*cw+3*u2*xi0w+w2*xi0w)+su*(-8*cw*u2+18*u2*xi0w+2*w2*xi0w) +
	    I* (-8*cw*(2*su*u+cu*u2)+2*(s2u*u+c2u*u2-c2u*w2)+2*(9*cu*u2-3*su*u*u2+cu*w2 -su*u*w2)*xi0w),
	    2*c2u-4*s2u*u+su*(2*cw*u+6*u*xi0w)+cu*(-2*cw+3*u2*xi0w-w2*xi0w) +
	     I* (2*s2u+4*c2u*u+2*cw*(su+cu*u)+(6*cu*u-3*su*u2+su*w2)*xi0w),
	    -2*s2u+cw*su-3*(su+cu*u)*xi0w +
	     I* (2*c2u+cu*cw+(-3*cu+3*su*u)*xi0w)},
	   {-2*c2u+2*cw*su*u+2*su*u*xi0w-8*cu*u2*xi0w+6*su*u*u2*xi1w +
	     I*(2*(-s2u+4*su*u2*xi0w+cu*u*(cw+xi0w+3*u2*xi1w))),
	    2*cu*u*xi0w+su*(-cw-xi0w+3*u2*xi1w) +
	     I*(-2*su*u*xi0w-cu*(cw+xi0w-3*u2*xi1w)),
	    cu*xi0w-3*su*u*xi1w +
	     I*(-(su*xi0w)-3*cu*u*xi1w)}};
	
	//change sign to the f if needed, eq. (34) - changed again below
	if(ing->sign!=0)
	  {
	    //ing->f[0][IM]*=-1;
	    ing->f[0]=conj(ing->f[0]);
        //ing->f[1][RE]*=-1;
        ing->f[1]= -conj(ing->f[1]);
	    //ing->f[2][IM]*=-1;
	    ing->f[2]=conj(ing->f[2]);
	  }
	
	//compute b
	double t1=9*u2-w2,t2=1/(2*t1*t1);
	for(int j=0;j<3;j++)
	  {
	    //eq. (57-58)
		b[0][j]=(2*u*r[0][j]+(3*u2-w2)*r[1][j]-2*(15*u2+w2)*ing->f[j])*t2;
		b[1][j]=(r[0][j]-3*u*r[1][j]-24*u*ing->f[j])*t2;
	    
	    //take into account the sign of c0, eq. (70)
	    if(ing->sign!=0)
	      for(int i=0;i<2;i++)
		{
		  //change the sign to real or imag part
		  int ri=(i+j+1)%2;
          // b[i][j][ri]=-b[i][j][ri]; 
          if(ri==0)
              b[i][j]=-conj(b[i][j]);
          if(ri==1)
              b[i][j]=conj(b[i][j]);
		}
	  }
	
	//change back sign to the f if needed, eq. (34)
	if(ing->sign!=0)
	  {
	    //ing->f[0][IM]*=-1;
	    ing->f[0]=conj(ing->f[0]);
        //ing->f[1][RE]*=-1;
        ing->f[1]= -conj(ing->f[1]);
	    //ing->f[2][IM]*=-1;
	    ing->f[2]=conj(ing->f[2]);

	  }
      }
    
    //compute B eq. (69)
    single_su3 B[2];
    for(int i=0;i<2;i++)
      {
	
         Itamat_2ndDeg_poly(b[i][0],b[i][1],b[i][2],QA,&B[i]) ;
//    su3_put_to_diag(B[i],b[i][0]);
//  su3_summ_the_prod_complex(B[i],ing->Q, b[i][1]);
//	su3_summ_the_prod_complex(B[i],ing->Q2,b[i][2]);
      }



    printf("c0 = %.18lf\n",ing->c0);
    printf("c1 = %.18lf\n",ing->c1);
    printf("f0 = %.18lf\n",ing->f[0]);
    printf("f1 = %.18lf\n",ing->f[1]);
    printf("f2 = %.18lf\n",ing->f[2]);
    //       printf("tr1 = (%.18lf) + (%.18lf)*I\n",creal(tr1),cimag(tr1));
    //       printf("tr2 = (%.18lf) + (%.18lf)*I\n",creal(tr2),cimag(tr2));
    printf("b10 = (%.18lf) + (%.18lf)*I\n",creal(b[0][0]),cimag(b[0][0]));
    printf("b20 = (%.18lf) + (%.18lf)*I\n",creal(b[1][0]),cimag(b[1][0]));
    printf("b11 = (%.18lf) + (%.18lf)*I\n",creal(b[0][1]),cimag(b[0][1]));
    printf("b21 = (%.18lf) + (%.18lf)*I\n",creal(b[1][1]),cimag(b[1][1]));
    printf("b12 = (%.18lf) + (%.18lf)*I\n",creal(b[0][2]),cimag(b[0][2]));
    printf("b22 = (%.18lf) + (%.18lf)*I\n",creal(b[1][2]),cimag(b[1][2]));





// segate (MI)
/*    
    //compute Gamma, eq. (74)
    su3 Gamma;
    //compute U*Sigma', to be used many times
    su3 aux;
    unsafe_su3_prod_su3(aux,U,F);
    //compute the trace of U*Sigma'*B[j]
    complex we[2];
    for(int j=0;j<2;j++) trace_su3_prod_su3(we[j],aux,B[j]);
    //first term
    unsafe_su3_prod_complex(Gamma,ing->Q,we[0]);
    //first term
    su3_summ_the_prod_complex(Gamma,ing->Q2,we[1]);
    //third term
    su3_summ_the_prod_complex(Gamma,aux,ing->f[1]);
    //fourth and fifth term
    su3 temp;
    unsafe_su3_prod_su3  (temp,ing->Q,aux);
    su3_summ_the_prod_su3(temp,aux,ing->Q);
    su3_summ_the_prod_complex(Gamma,temp,ing->f[2]);
    
    //compute Lambda eq. (73)
    unsafe_su3_traceless_hermitian_part(Lambda,Gamma);

    */
  }


int main(){

   // Lambda 
    single_su3 gl3_temp0;
    single_su3 sSP;
    single_tamat sQA;
    su3_soa * U =   (su3_soa * ) malloc(sizeof(su3_soa)); 
    su3_soa * SP =  (su3_soa * ) malloc(sizeof(su3_soa));
    su3_soa * tmp = (su3_soa * ) malloc(sizeof(su3_soa));
    thmat_soa * L =  (thmat_soa * ) malloc(sizeof(thmat_soa));
    tamat_soa * QA = (tamat_soa * ) malloc(sizeof(tamat_soa));
    int idx;


    // FOR Lambda testing
    idx = 0;

    sQA.rc00 = -0.5;
    sQA.rc11 = 0.5+1.2*I;
    sQA.c01 = -1.0*I;
    sQA.c02 = 1.2+2*I;
    sQA.c12 = 2.3*I;

    CH_exponential_antihermitian_nissalike(&gl3_temp0,&sQA);
    // sSP for lambda
    sSP.comp[0][0] =  0.4-4.4*I;
    sSP.comp[0][1] =  1.5+5.1*I;
    sSP.comp[0][2] = -0.2+0.5*I;
    sSP.comp[1][0] =  0.7-1.7*I;
    sSP.comp[1][1] = -1.8+0.8*I;
    sSP.comp[1][2] =  1.9+1.5*I;
    sSP.comp[2][0] =  0.6-0.4*I;
    sSP.comp[2][1] = -1.5+1.1*I;
    sSP.comp[2][2] =  2.4+5.2*I;
    single_gl3_into_su3_soa(SP,idx,&sSP);
    // U for lambda
    single_su3_into_su3_soa(U,idx,&gl3_temp0);
    // QA for lambda
    single_tamat_into_tamat_soa(QA,idx,&sQA);
    // temp should not necessitate any initialization, but...
    // ad cazzum initialization
    single_gl3_into_su3_soa(tmp,idx,&sSP);


    DEOTT_compute_loc_Lambda(
            L, 
            SP,
            U,
            QA,
            tmp,
            idx );


 anti_hermitian_exp_ingredients ing;
 anti_hermitian_exact_i_exponentiate_ingredients(&ing, &sQA);
 stouted_force_compute_coeffs_NISSA(&ing,&sQA);

}
