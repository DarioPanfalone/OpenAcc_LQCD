#ifndef HOME_ACOS_C
#define HOME_ACOS_C

#define FC00 1.0000000000000000000
#define FC01 0.1666666666666666667
#define FC02 0.4500000000000000000
#define FC03 0.5952380952380951324
#define FC04 0.6805555555555559528
#define FC05 0.7363636363636368121
#define FC06 0.7756410256410246628
#define FC07 0.8047619047619061329
#define FC08 0.8272058823529423090
#define FC09 0.8450292397660798621
#define FC10 0.8595238095238085963
#define FC11 0.8715415019762845588
#define FC12 0.8816666666666654268
#define FC13 0.8903133903133937477
#define FC14 0.8977832512315269865
#define FC15 0.9043010752688167351
#define FC16 0.9100378787878786412
#define FC17 0.9151260504201654769
#define FC18 0.9196696696696787060
#define FC19 0.9237516869095804500
#define FC20 0.9274390243902376965
#define FC21 0.9307862679955640940
#define FC22 0.9338383838383934450
#define FC23 0.9366327474560528359
#define FC24 0.9392006802721121926
#define FC25 0.9415686274509820285
#define FC26 0.9437590711175601763
#define FC27 0.9457912457912391206
#define FC28 0.9476817042606531270
#define FC29 0.9494447691408663556

#define FSD00 1.0000000000000000000
#define FSD01 0.0833333333333333333
#define FSD02 0.2250000000000000000
#define FSD03 0.2976190476190475662
#define FSD04 0.3402777777777779764
#define FSD05 0.3681818181818184061
#define FSD06 0.3878205128205123314
#define FSD07 0.4023809523809530664
#define FSD08 0.4136029411764711545
#define FSD09 0.4225146198830399311
#define FSD10 0.4297619047619042981
#define FSD11 0.4357707509881422794
#define FSD12 0.4408333333333327134
#define FSD13 0.4451566951566968738
#define FSD14 0.4488916256157634932
#define FSD15 0.4521505376344083676
#define FSD16 0.4550189393939393206
#define FSD17 0.4575630252100827384
#define FSD18 0.4598348348348393530
#define FSD19 0.4618758434547902250
#define FSD20 0.4637195121951188483
#define FSD21 0.4653931339977820470
#define FSD22 0.4669191919191967225
#define FSD23 0.4683163737280264179
#define FSD24 0.4696003401360560963
#define FSD25 0.4707843137254910143
#define FSD26 0.4718795355587800881
#define FSD27 0.4728956228956195603
#define FSD28 0.4738408521303265635
#define FSD29 0.4747223845704331778

#pragma acc routine seq
double homebrew_acos(double x){
  double x2=x*x;
  double y=(double)1.0-fabs(x);
  double w=sqrt(2.0*y);
  double r;
  if(x2<=0.25){
  return M_PI*0.5-
         FC00*x*((double)1.0 +x2*FC01*(
		  (double)1.0 +x2*FC02*(
		  (double)1.0 +x2*FC03*(
		  (double)1.0 +x2*FC04*(
		  (double)1.0 +x2*FC05*(
		  (double)1.0 +x2*FC06*(
		  (double)1.0 +x2*FC07*(
		  (double)1.0 +x2*FC08*(
		  (double)1.0 +x2*FC09*(
		  (double)1.0 +x2*FC10*(
		  (double)1.0 +x2*FC11*(
		  (double)1.0 +x2*FC12*(
		  (double)1.0 +x2*FC13*(
		  (double)1.0 +x2*FC14*(
		  (double)1.0 +x2*FC15*(
		  (double)1.0 +x2*FC16*(
		  (double)1.0 +x2*FC17*(
		  (double)1.0 +x2*FC18*(
		  (double)1.0 +x2*FC19*(
		  (double)1.0 +x2*FC20*(
		  (double)1.0 +x2*FC21*(
		  (double)1.0 +x2*FC22*(
		  (double)1.0 +x2*FC23)))))))))))))))))))))));
  }else{
    r =  FSD00*w*((double)1.0 +y*FSD01*(
		  (double)1.0 +y*FSD02*(
	 	  (double)1.0 +y*FSD03*(
		  (double)1.0 +y*FSD04*(
		  (double)1.0 +y*FSD05*(
		  (double)1.0 +y*FSD06*(
		  (double)1.0 +y*FSD07*(
		  (double)1.0 +y*FSD08*(
		  (double)1.0 +y*FSD09*(
		  (double)1.0 +y*FSD10*(
		  (double)1.0 +y*FSD11*(
		  (double)1.0 +y*FSD12*(
		  (double)1.0 +y*FSD13*(
		  (double)1.0 +y*FSD14*(
		  (double)1.0 +y*FSD15*(
		  (double)1.0 +y*FSD16*(
		  (double)1.0 +y*FSD17*(
		  (double)1.0 +y*FSD18*(
		  (double)1.0 +y*FSD19*(
		  (double)1.0 +y*FSD20*(
		  (double)1.0 +y*FSD21*(
		  (double)1.0 +y*FSD22*(
		  (double)1.0 +y*FSD23)))))))))))))))))))))));
  }

  if(x>0.5){
    return r;
  }else{
    return M_PI-r;
  }
}

#endif
