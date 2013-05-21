//  Basic RBC Model with volatility shocks
// 
//  Antti Ripatti, 8.2.2011

//  modified by Markus Haavio 13.4.2011
//  - model version where the bank's  incentive constraint does not 
//   depend on ra, the return to bank capital
//  - endogenous labour supply, disutilty function v(L)=(ksi/(1+fi)*L^(1+fi))
//  - steady state distortion (too little capital + investments) corrected by
//  an investment subsidy
//  steady identical to the steady state of the corresponding standard RBC
//  model 
//  ...
//  many modifications by Markus 
//  ...
// 
//  modified by Antti Ripatti, 15.3.2012 
//  - Gertler-Karadi setup
//  modified by Antti Ripatti 3.4.2012
//  - analytical steady state matches that of the dynamic model
//  modified by Markus Haavio, May
//  - analytical calibration of parameters. 
//  modified by Antti Ripatti 7.6.2012
//  - corrections found by Markus
//  - goverment capital injection
//  - plotting stuff added, Antti 13.6.2012
// 
//close all;
// //  Define variables, parameters, shocks, etc.




var 
  mU //  marginal utility of consumption
  Y  //  output
  I //  investments
  C //  consumption
  K //  capital stock
  L //  hours/employment, here =1 
  Z //  technology shock
  Zc //  preference shock
  rK //  capital return
  w //  real wage
  e_I
  sigma_I
  q
;

//  exogenous shock
varexo e_, e_c, ee_I, e_sigma;
//  list of parameters
parameters 
  sigma2 //  shock variance
  sigma //  risk aversion
  delta //  depreciation rate
  beta //  discount rate
  alpha //  capital share
  rho //  shock persistence
  rhoc //preference shock persistence
  ksi //  marginal disutility of labour 
  fi //  inverse of frich elasticity of labour supply
  rho_eI
  b //  habit persistence
  cstar // monitoring cost
;


// //  Parameter values, King-Rebelo (2000)
beta = 1/1.02^0.25; 
alpha = 1/3;
delta = 0.025;
ksi = 2; //  2
fi = 0.1; //  3
rho = 0.979; 
rhoc = 0.6;
sigma = 1;
sigma2 = ((0.007/(1-alpha)))^2;
cstar = 0.000825;
//cstar=0.008;

//  Hansen (1985) parameterization
//beta = 0.99;
//alpha = 0.36;
//delta = 0.0025;
//rho = 0.95;
//sigma = 1;


rr=1/beta-1;
rho_eI=.9;
b=0;

// 
//  Model code: start equation at line x1 in order to be able to map the
//  equation numbers that dynare gives with the code below
model; 
// (use_dll);
q=1+cstar;
  mU = Zc*((C-b*C(-1))/(1-b))^(-sigma); 
  mU = beta*mU(+1)*(q*(1-delta) + rK(+1))/q; 
//   log(Z) = rho*log(Z(-1)) + e_; 
  Z = rho*Z(-1) + e_; 
  log(Zc)=rhoc*log(Zc(-1)) + e_c; 
  Y = K(-1)^alpha*(exp(Z)*L)^(1-alpha); 
  rK/w = (alpha/(1-alpha))*L/K(-1); 
  rK = alpha*Y/K(-1);
  K = (1-delta)*K(-1) + I*(1+e_I); 
  L = (w*mU/ksi)^(1/fi); 
  (C + I*(1+cstar)) = Y;  //  resource constraint corrected for monitoring costs
  sigma_I = (1-rho_eI) + rho_eI*sigma_I(-1) + e_sigma; 
  e_I=sigma_I*ee_I; 
end;




// 
//  Compute stuff
initval; //  this is the analytical steady state
  Z = 0; 
  Zc = 1;
  e_ = 0;
  e_I = 0;
  sigma_I = 1; 
 q=1+cstar;
  //  return to capital
 
  rK = q*(1/beta-1+delta);
  //  real wage
  w=(1-alpha)*(rK/alpha)^(-alpha/(1-alpha));
  //  capital stock
  //K =((1-alpha)/ksi)^(1/(sigma+fi))*(rK/alpha)^(-(alpha+fi)/((1-alpha)*(sigma+fi)))*(rK/alpha-delta*q)^(-sigma/(sigma+fi));
  K =((1-alpha)/ksi)^(1/(sigma+fi))*(rK/alpha)^(-(alpha+fi)/((1-alpha)*(sigma+fi)))*(rK/alpha-delta*(1+cstar))^(-sigma/(sigma+fi));
  //  labour 
  L=K*(rK/alpha)^(1/(1-alpha));
  //  output
  Y = K*rK/alpha;
  //  investments
  I = delta*K;
  //  insider wealth, J = A + N
  C = Y - I*(1+cstar); // -I*chi/pH;
  //  consumption of entrepreneurs
  mU =C^(-sigma);
end;

//steady;
//resid;
//e_, e_A, e_N, e_mu, e_c, e_m, ee_I, e_sigma;

//  shocks;
//    var e_ = 0.01^2;
//    var e_c = 0;
//    var ee_I = 0.01^2;
//    var e_sigma = 0.01^2;
//  end;
vcov = [0.00005184 0 0 0 ,
          0 0 0 0 ,
          0 0 0 0,
          0 0 0 0 ];

vcov = [0.0000168 0 0 0 ;
          0 0 0 0 ;
          0 0 0.0001 0;
          0 0 0 0.0001 ];
           
order = 3;
           

 //resid(1);
// steady(solve_algo = 2);
// check;
// stoch_simul(order=3, irf=40);
