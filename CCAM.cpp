// Censored Catch Assessment model
//  a Statistical catch-at-age model with (un)censored catch
// Authors:
// 		*Elisabeth Van Beveren <Elisabeth.vanBeveren@dfo-mpo.gc.ca>
// 		*Daniel Duplisea
// 		*Martin Castonguay
// 		*Thomas Doniol-Valcroze
// 		*Stephane Plourde
// 		*Noel Cadigan
// Journal:  (submitted)
// Title: 	How catch underreporting can bias stock assessment and advice in northwest Atlantic mackerel 
// 		and a possible resolution using censored catch
//
// © 2017 E. Van Beveren, D. Duplisea, M. Castonguay, T. Doniol-Valcroze, S. Plourde, N. Cadigan 	  
// All rights reserved.

#include <TMB.hpp>
#include <algorithm>

template<class Type>
Type objective_function<Type>::operator() () {

//---------------------------------------
//------------ data ---------------------
//---------------------------------------
DATA_VECTOR(y);  		// vector with years
DATA_VECTOR(a); 		// vector with ages
DATA_VECTOR(logClower);		// lower boundary of catches in kg
DATA_VECTOR(logCmean);		// average catches in kg (lower+upper)/2
DATA_VECTOR(logCupper);  	// upper boundary of catches in kg
DATA_ARRAY(crl);		// continous-ratio logit catch composition
DATA_ARRAY(M);			// natural mortality by age and year
DATA_ARRAY(Weight);		// weight by age and year
DATA_ARRAY(propMature);	// proportion mature by age and year

DATA_VECTOR(ys);		// years of the survey (continuous)
DATA_VECTOR(as);		// ages of survey (ages summed to obtain the annual index)
DATA_SCALAR(Stime);		// proportion of year
DATA_VECTOR(S); 		// survey SSB values (kg)  (negative values if NA)
DATA_ARRAY(WeightS);		// weight by age and year, at index time (Stime)

DATA_VECTOR(env1);		// environmental factor 1 (zeros if not used)
DATA_VECTOR(env2);		// environmental factor 2 (zeros if not used)

DATA_INTEGER(state);		// 1: uncensored: fit to logCmean, 2: censored model
DATA_INTEGER(rec);		// recruitment options
DATA_INTEGER(recDist);		// recruitment error distribution: 0: (log)normal, 1:t , 2: robuust (mixture), 3: robuust + (mixture)

DATA_INTEGER(logsd_logC);	// standard deviation catches

DATA_VECTOR(f); 		// vector of mortalities

//---------------------------------------
//------------ parameters ---------------
//---------------------------------------
PARAMETER_ARRAY(logN); 		// abundance
PARAMETER_VECTOR(logFa);	// Fishing mortality per age
PARAMETER_VECTOR(logFy);	// Fishing mortality per year
PARAMETER(logsd_logS);		// standard deviation survey
PARAMETER(logsd_logFy);		// standard deviation Fy
PARAMETER(logsd_logrec);	// standard deviation recruitment 
PARAMETER_VECTOR(logsd_crl);	// standard deviation age composition
PARAMETER(logrecruitMean);      // mean recruitment (if rec = 1)

PARAMETER(logQ);		// survey catchability

PARAMETER(R1);			// parameter stock-recruitment
PARAMETER(R2);			// parameter stock-recruitment
PARAMETER(E1);			// parameter stock-recruitment environment
PARAMETER(E2);			// parameter stock-recruitment environment

PARAMETER(log_std_pe);   	// standard deviation process error         
PARAMETER(logit_ar_pe_age);     // autocorrelation PE age
PARAMETER(logit_ar_pe_year);    // autocorrelation PE year

//------- simple objects and transformations --------
//- lengths
int ny = y.size(); 
int na = a.size();
int nys = ys.size(); 
int nas = as.size();

//- mortality
vector<Type> Fy=exp(logFy);  
vector<Type> Fa=exp(logFa);
array<Type> F(na,ny);		
array<Type> ZZ(na,ny);	

//- Numbers and masses
array<Type> N(na,ny); 
vector<Type> ssb(ny);
vector<Type> tsb(ny);

//- Catches
array<Type> logCpred(na,ny);
array<Type> Cpred(na,ny);
array<Type> CpredW(na,ny);
vector<Type> logCpredtot(ny);
vector<Type> logCpredtotW(ny);
vector<Type> Cpredtot(ny);
vector<Type> CpredtotW(ny);
array<Type> propCpred(na,ny);
array<Type> propCpred_cond(na-1,ny);
array<Type> pred_crl(na-1,ny);

//- Residuals
vector<Type> RESstd_rec(ny);
vector<Type> RESraw_rec(ny);
vector<Type> RESraw_Spred(nys);
vector<Type> RESstd_Spred(nys);
array<Type> RESstd_crl(na-1,ny);
vector<Type> RESstd_Cobs(ny);
array<Type> RESraw_crl(na-1,ny);
vector<Type> RESraw_Cobs(ny);

//- other
Type one = 1.0; 
Type zero = 0.0;

Type nll = 0.0;

//---------------------------------------
//------- START calculations ------------
//---------------------------------------

//-- mortality -----------
for(int i=0;i<na;i++){ 
 for(int j=0;j<ny;j++){ 
   F(i,j)=Fy(j)*Fa(i);		// separable mortality
   ZZ(i,j)=F(i,j)+M(i,j);
 }
}

for(int j = 1 ; j < ny; j++){
 nll -= dnorm(logFy(j), logFy(j-1), exp(logsd_logFy), true); 	// random walk for Fy
}

//-- SSB -----------
for(int j=0;j<ny;j++){ 
 ssb(j)=0; 
 tsb(j)=0;       
  for(int i=0; i<na; ++i){
   ssb(j)+=exp(logN(i,j))*Weight(i,j)*propMature(i,j);
   tsb(j)+=exp(logN(i,j))*Weight(i,j);
  }
}

//----- PE-----------
array<Type> pe(na,ny);  
Type phi_age = exp(logit_ar_pe_age)/(one + exp(logit_ar_pe_age));
Type phi_year = exp(logit_ar_pe_year)/(one + exp(logit_ar_pe_year));                
Type sigmape = exp(log_std_pe);

Type mZ;
 
pe(0,0)=zero;

//----- N: Recruitment process -----------
vector<Type> predR(ny);

if(rec==1){  // Centered around mean
for(int j=0;j<ny;j++){ 
  RESstd_rec(j)=(logN(0,j)-logrecruitMean)/exp(logsd_logrec); 
  RESraw_rec(j)=(logN(0,j)-logrecruitMean);
  }
}

if(rec==2){  // Centered around linear model (fails)
for(int j=0;j<ny;j++){ 
  RESstd_rec(j)=(logN(0,j)-log(R1+R2*j))/exp(logsd_logrec); 
  RESraw_rec(j)=(logN(0,j)-log(R1+R2*j));
  }
}

if(rec==3){  // Random walk
for(int j=1;j<ny;j++){ 
  predR(j)=logN(0,j-1);
  RESstd_rec(j)=(logN(0,j)-predR(j))/exp(logsd_logrec); 
  RESraw_rec(j)=(logN(0,j)-predR(j));
 }
}

if(rec>3){ //SR relationship
for(int j=1;j<ny;j++){ 
   if(rec==4){predR(j) = R1+log(ssb(j-1))-exp(R2)*ssb(j-1);}; //ricker
   if(rec==5){predR(j) = R1+log(ssb(j-1))-log(1.0+exp(R2)*ssb(j-1));}; // basic BH
   if(env1.sum()!=0){predR(j) += E1*env1(j-1);};         // add environmental factor one
   if(env2.sum()!=0){predR(j) += E2*env2(j-1);};         // add environmental factor two
   if(recDist==0){predR(j) += Type(0.5)*exp(logsd_logrec)*exp(logsd_logrec);};   // lognormal distribution
   RESstd_rec(j)=(logN(0,j)-predR(j))/exp(logsd_logrec); 
   RESraw_rec(j)=(logN(0,j)-predR(j));
 }
}

if(rec>2){
for(int j=1;j<ny;j++){
if(recDist==0){nll -= dnorm(RESraw_rec(j),zero,exp(logsd_logrec),true);};
if(recDist==1){nll -= dt(RESstd_rec(j),Type(3),true);};
if(recDist==2){nll -=(Type(0.95)*dnorm(RESraw_rec(j),zero,exp(logsd_logrec),true)+Type(0.05)*dt(RESstd_rec(j),Type(3),true));};
if(recDist==3){nll -=(Type(0.95)*dnorm(RESraw_rec(j),zero,exp(logsd_logrec),true)+Type(0.05)*dt(RESstd_rec(j),Type(1),true));};
  }
}else{
for(int j=0;j<ny;j++){
if(recDist==0){nll -= dnorm(RESraw_rec(j),zero,exp(logsd_logrec),true);};
if(recDist==1){nll -= dt(RESstd_rec(j),Type(3),true);};
if(recDist==2){nll -=(Type(0.95)*dnorm(RESraw_rec(j),zero,exp(logsd_logrec),true)+Type(0.05)*dt(RESstd_rec(j),Type(3),true));};
if(recDist==3){nll -=(Type(0.95)*dnorm(RESraw_rec(j),zero,exp(logsd_logrec),true)+Type(0.05)*dt(RESstd_rec(j),Type(1),true));};
  }
}
//----- N:  non-recruits -----------
  vector<Type> predN(na-1);
  for(int i=1; i<na; ++i){
   pe(i,0)=zero;
      for(int j=1;j<ny;j++){  
      predN(i-1)=logN(i-1,j-1)-ZZ(i-1,j-1);	
        if(i==(na-1)){
           predN(i-1)=log(exp(logN(i,j))+exp(logN(i,j-1)-ZZ(i,j-1)));	// plus group
        }
      pe(i,j)=logN(i,j)-predN(i-1);
      mZ = phi_year*pe(i,j-1) + phi_age*(pe(i-1,j) - phi_year*pe(i-1,j-1));
      nll -= dnorm(pe(i,j),mZ,sigmape,true);    
     }
  }

//------- Catch -----------
//--- Estimation
for(int j = 0 ; j < ny; j++){
 for(int i = 0 ; i < na; i++){
  N(i,j)=exp(logN(i,j));
  logCpred(i,j)=log(F(i,j))-log(ZZ(i,j))+log(one-exp(-ZZ(i,j)))+logN(i,j); // Baranov catch equation
  Cpred(i,j)=exp(logCpred(i,j));
  CpredW(i,j)=Cpred(i,j)*Weight(i,j);
 }
 Cpredtot(j) = Cpred.col(j).sum();
 CpredtotW(j) = CpredW.col(j).sum();  // for fit to total catch 
 logCpredtotW(j) = log(CpredtotW(j));
 for(int i = 0 ; i < na; i++){
  propCpred(i,j) = Cpred(i,j)/Cpredtot(j);  // for at age proportions
 }
}

//--- Match total landings
if(state==1){ // match to the mean
 for(int j = 0 ; j < ny; j++){
  RESraw_Cobs(j) = (logCpredtotW(j) - logCmean(j));
  RESstd_Cobs(j) = (logCpredtotW(j) - logCmean(j))/exp(logsd_logC);
  nll -= (dnorm(RESstd_Cobs(j),zero,one,true)-logsd_logC);
 }
}
if(state==2){ // censored
vector<Type> ZU(ny);
vector<Type> ZL(ny);
 for(int j = 0 ; j < ny; j++){
  ZU(j) = (logCupper(j) - logCpredtotW(j))/exp(logsd_logC);
  ZL(j) = (logClower(j) - logCpredtotW(j))/exp(logsd_logC);
  nll -= log(pnorm(ZU(j)) - pnorm(ZL(j)));
 }
}

//--- Match proportions
Type total;  
for(int j = 0 ;j <ny; j++){
  propCpred_cond(0,j) = propCpred(0,j); 
  for(int i = 1 ;i < na-1; i++){
    total=0.0;
    for(int k = 0 ;k < i; k++){total += propCpred(k,j);}  
    propCpred_cond(i,j) = propCpred(i,j)/(one - total);
  }
}

for(int j = 0 ;j <ny; j++){
 for(int i = 0 ; i < na-1; i++){
   pred_crl(i,j) = log(propCpred_cond(i,j)/(one - propCpred_cond(i,j)));
   RESraw_crl(i,j) = (crl(i,j) - pred_crl(i,j));
   RESstd_crl(i,j) = (crl(i,j) - pred_crl(i,j))/exp(logsd_crl(i));
   nll -= (dnorm(RESstd_crl(i,j),zero,one,true)-logsd_crl(i));
 }
}   

//------- Estimate Survey SSB---------------
array<Type> logSpredN(nas,nys);
vector<Type> SpredW(nys);
array<Type> propMatureS=propMature.matrix().block(na-nas,ny-nys,nas,nys); 

for(int j = 0 ; j < nys; j++){
 SpredW(j)=0.0; 
  for(int i = 0 ; i < nas; i++){
   logSpredN(i,j)=exp(logQ)-ZZ(i,ny-nys+j)*Stime+logN(i,ny-nys+j);      
    SpredW(j)+=exp(logSpredN(i,j))*propMatureS(i,j)*WeightS(i,j);
  }
}

//------- Match Survey---------------
for(int j = 0; j < nys; j++){ 
  if(S(j)>0){
   nll -= dnorm(log(S(j))-log(SpredW(j)),zero,exp(logsd_logS),true); 
   RESraw_Spred(j)=log(S(j))-log(SpredW(j));
   RESstd_Spred(j)=(log(S(j))-log(SpredW(j)))/exp(logsd_logS);
   }else{
    RESraw_Spred(j)=-999999;
    RESstd_Spred(j)=-999999;
 }
}

//---------------------------------------
//------- Reference Points ------------
//---------------------------------------
vector<Type> RP(19);  

//- Catch selectivity
vector<Type> Sel(na); 
for(int i = 0 ; i < na; i++){Sel(i)=F(i,ny-1)/F.col(ny-1).maxCoeff();};

//- spr and ypr
int nf = f.size(); 
vector<Type> lx(na);
vector<Type> spr(nf);
vector<Type> spr_inv(nf);
vector<Type> ypr(nf);
for(int k = 0 ;k <nf; k++){
  lx(0)=one;
  spr(k)=lx(0)*Weight(0,ny-1)*propMature(0,ny-1);
  ypr(k)=zero;
  for(int i = 1 ; i < na; i++){
    lx(i)=lx(i-1)*exp(-(M(i-1,ny-1)+Sel(i-1)*f(k))+pe(i,ny-1)); 
    if(i==(na-1)){lx(i)=lx(i)*(1+exp((M(i,ny-1)+Sel(i)*f(k))+pe(i,ny-1)));};
    spr(k)+=lx(i)*Weight(i,ny-1)*propMature(i,ny-1);
  }
  spr_inv(k)=one/spr(k);
  for(int i = 0 ; i < na; i++){
    ypr(k)+=lx(i)*Weight(i,ny-1)*(1-exp(-Sel(i)*f(k)));
  }
}

//- sprF0, spr20, spr30, spr40,Fmax		could be looped
vector<Type> sprPercent(nf);
vector<Type> sprPercent20(nf);
vector<Type> sprPercent30(nf);
vector<Type> sprPercent40(nf);
for(int k = 0 ;k <nf; k++){
  sprPercent(k)=spr(k)/spr.maxCoeff();
  sprPercent20(k)=abs(sprPercent(k)-Type(0.20));
  sprPercent30(k)=abs(sprPercent(k)-Type(0.30));
  sprPercent40(k)=abs(sprPercent(k)-Type(0.40));
}

for(int k = 0 ;k <nf; k++){
  if(sprPercent20(k)==sprPercent20.minCoeff()){RP(0)=f(k);};		//spr20
  if(sprPercent30(k)==sprPercent30.minCoeff()){RP(1)=f(k);};		//spr30
  if(sprPercent40(k)==sprPercent40.minCoeff()){RP(2)=f(k);};		//spr40
  if(ypr(k)==ypr.maxCoeff()){RP(3)=f(k);};				//Fmax
}

RP(4)=spr(0);								//sprF0

//- F01
vector<Type> slope(nf-1);
vector<Type> slopePercent(nf-1);
vector<Type> diff10(nf-1);
for(int k = 1 ;k <nf; k++){
  slope(k-1)=(ypr(k)-ypr(k-1))/(f(k)-f(k-1));
  slopePercent(k-1)=slope(k-1)/slope(0);
  diff10(k-1)=abs(slopePercent(k-1)-Type(0.1));
}
for(int k = 0 ;k < (nf-1); k++){
  if(diff10(k)==diff10.minCoeff()){RP(5)=f(k);};
}

//- Fhigh, Flow, Fmed       		      could be looped
vector<Type> recssb(ny);
for(int j = 0 ; j < ny; j++){recssb(j)=N(0,j)/ssb(j);};

std::sort(recssb.data(),recssb.data()+recssb.size());

vector<Type> Q01_diff(nf);
vector<Type> Q05_diff(nf);
vector<Type> Q09_diff(nf);
for(int k = 0 ;k < nf; k++){
  Q01_diff(k)=abs(spr_inv(k)-recssb(round(0.1*ny+1)));
  Q05_diff(k)=abs(spr_inv(k)-recssb(round(0.5*ny+1)));
  Q09_diff(k)=abs(spr_inv(k)-recssb(round(0.9*ny+1)));
}

for(int k = 0 ;k < nf; k++){
  if(Q01_diff(k)==Q01_diff.minCoeff()){RP(6)=f(k);};		//Flow
  if(Q05_diff(k)==Q05_diff.minCoeff()){RP(7)=f(k);};		//Fmed
  if(Q09_diff(k)==Q09_diff.minCoeff()){RP(8)=f(k);};		//Fhigh
}

//- Fmsy, SSBmsy, RECmsy : SO FAR ONLY FOR BH and 1 environmental effect!!!  
vector<Type> equiS(nf);
vector<Type> equiR(nf); 
vector<Type> yield(nf); 

vector<Type> equiS_noENV(nf);
vector<Type> equiR_noENV(nf); 
vector<Type> yield_noENV(nf); 

if(rec==5){
  Type alpha=exp(R1);
  Type K=one/exp(R2);
  
  for(int k = 0 ;k < nf; k++){
    equiS_noENV(k)=(alpha*spr(k)-one)*K; 
      equiS(k)=(alpha*spr(k)*exp(E1*env1(ny-1)-(Type(0.5)*exp(logsd_logrec)*exp(logsd_logrec)))-one)*K;
      equiR_noENV(k)=(alpha*equiS_noENV(k))/(one+(equiS_noENV(k)/K));
      equiR(k)=(alpha*equiS(k)*exp(E1*env1(ny-1)-(Type(0.5)*exp(logsd_logrec)*exp(logsd_logrec))))/(one+(equiS(k)/K));
      yield_noENV(k)=equiR_noENV(k)*ypr(k);
      yield(k)=equiR(k)*ypr(k);
  }
  
  for(int k = 0 ;k < nf; k++){
    if(yield(k)==yield.maxCoeff()){RP(9)=f(k);        //Fmsy
    RP(10)=equiS(k);   //SSBmsy
    RP(11)=equiR(k);};   //RECmsy  
    if(yield_noENV(k)==yield_noENV.maxCoeff()){RP(14)=f(k);        //Fmsy without environmental effect 
    RP(15)=equiS_noENV(k);   //SSBmsy without environmental effect 
    RP(16)=equiR_noENV(k);};   //RECmsy without environmental effect                                		
  }
  
  RP(12)=RP(10)*0.4;						//SSBmsy40
  
  vector<Type> alpha_diff(nf);
  for(int k = 0 ;k <nf; k++){
    alpha_diff(k)=abs(spr_inv(k)-alpha);
  }
  
  for(int k = 0 ;k < nf; k++){
    if(alpha_diff(k)==alpha_diff.minCoeff()){RP(13)=f(k);};		//Fcol=Fcrash
  }
}

// Brecovery
vector<Type> smaller20(ny-1);
vector<Type> smaller40(ny-1);
for(int j = 0 ; j < ny-1; j++){
  if(ssb.tail(ny-j).maxCoeff()>ssb.maxCoeff()*Type(0.2)){smaller20(j)=ssb(j);}else{smaller20(j)=99999999999;};
  if(ssb.tail(ny-j).maxCoeff()>ssb.maxCoeff()*Type(0.4)){smaller40(j)=ssb(j);}else{smaller40(j)=99999999999;};
}

RP(17)=smaller20.minCoeff();			// Brec 20
RP(18)=smaller40.minCoeff();			// Brec 40
    
//------- Report---------------

// random effects
REPORT(logN);
REPORT(logFy);

// of interest
REPORT(ssb);
REPORT(RP);
REPORT(CpredtotW);
REPORT(F);
REPORT(pe);
REPORT(logFa);
REPORT(SpredW);
REPORT(ypr);
REPORT(spr);

// sd
ADREPORT(ssb);

return nll;
}


