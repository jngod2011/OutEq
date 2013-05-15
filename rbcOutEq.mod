// rbcOutEq.mod
// Financial contract with macro investment shocks and outside equity
// 
//  (c) Antti Ripatti and Markus Haavio, 2011-2013
//
//  builds on ./code/GK/buffer.mod and RBC.mod by Haavio and Ripatti
// 
// 


//  Define variables, parameters, shocks, etc.
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
  e_I // investment shock 
  sigma_I // uncertainty shock (stochastic volatility)
  PI // Pledgeable income 
  q // Tobin's q, price of capital
  phie // outside equity investment ratio, firms
  phib // outside equity investment ratio, banks
  // marginal value of capital 
  ve // entrepreneurial
  vb // banks 
  vw // worker-owned
  Qbhat // risk-adjusted price of banks' outside equity
  ratilde // return on banks capital
  retilde // return on entrepreneurial capital
  rhobhat // risk-adjusted NPB of the investment project
  rd // deposit rate
  etae // risk adjustment 
  etab
  etaw 
  G // inverse of the leverage
  A // bankers' capital
  N // firms' capital
  M // informed capital
  omega // firms capital relative to informed capital
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
  pH
  dp 
  bstar //  private benefit
  cstar //  banks monitoring cost
  kappa0e
  kappa0b
  kappa1e
  kappa1b // parameters of the efficiency costs, ie kappa0*phi^kappa1
  R // size of the cake
  lambdab // bankers' survival probability
  lambdae // firms' survival probability
  s // investment subsidy
  Rstar
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

//  Hansen (1985) parameterization
//beta = 0.99;
//alpha = 0.36;
//delta = 0.0025;
//rho = 0.95;
//sigma = 1;


rr=1/beta-1;
rho_eI=.9;
b=0;

//  Financial sector calibration
//  Based on the following data
//  riskless rate: 1+rbar = 
rbar = 1.02^0.25; 
//  excess return on banks capital: 1+rabar = 1.20^(1/4)
rabar = 1.13^0.25;
rabar = 1.20^0.25;
//rabar = 1.10^0.25;
//  excess return on 1+rebar = (1.065-1.02)^0.25
rebar = (1.065-rbar^4+1)^0.25;
//rebar = (1.08-rbar^4+1)^0.25;
//  Entrepreneurs capital ratio N/I = 0.3
CRF = 0.3;
CRF = 0.45;
//  Banks capital ratio CRB= A/(I-N) = 0.08
CRB = 0.08;
//  Banks  monitoring costs 0.03 relative to assets (per annum / 4)
CORB = 0.01/4;
CORB = 0.006/4;
//  calculated parameter values
lambdab = beta/rabar;
lambdae = beta/rebar;
lambdab = 1-lambdab;
lambdae = 1-lambdae;
pH = 0.95;
rdp = CORB/(CRB*rabar);
dp = pH*rdp;                 
// dp = .5;
pL = pH - dp;
R = 1/pH;
bstar = 0.00815075;
cstar = 0.000825;
// cstar = 0.825;
kappa0e = 1;
kappa0b = 1;
kappa1e = 0.1;
kappa1b = 0.1;
// kappa0e = 0.24;
// kappa0b = 0.002;
// kappa1e = 7;
// kappa1b = 11;
s = ((pH/dp)*(1 - ((1-lambdae)/beta))*bstar + (pH/dp)*(1 - ((1-lambdab)/beta))*cstar + cstar - cstar)/(1+cstar);
rhobhatbar =   (pH/dp)*(1 - ((1-lambdae)/beta))*bstar + (pH/dp)*(1 - ((1-lambdab)/beta))*cstar + cstar;
Rstar = R*(1+rhobhatbar)/(1+cstar);





// 
//  Model code: start equation at line x1 in order to be able to map the
//  equation numbers that dynare gives with the code below
model; 
  mU = Zc*((C-b*C(-1))/(1-b))^(-sigma); 
  mU = beta*mU(+1)*((q(+1)*(1-delta) + rK(+1))/q); 
  Z = rho*Z(-1) + e_; 
  log(Zc)=rhoc*log(Zc(-1)) + e_c; 
  Y = K(-1)^alpha*(exp(Z)*L)^(1-alpha); 
  rK/w = (alpha/(1-alpha))*L/K(-1); 
  rK = alpha*Y/K(-1);
  K = (1-delta)*K(-1) + I*(1+e_I); 
  L = (w*mU/ksi)^(1/fi);
  (C + I) = Y; // 10 
  sigma_I = (1-rho_eI) + rho_eI*sigma_I(-1) + e_sigma; 
  e_I=sigma_I*ee_I;
  rd = 0;
  PI = pH*q*R/(1+rd) - pH*etae*bstar/(dp*(1+rd)); // (15)
  (1/etaw - 1/etab)*(PI - etaw*phie) = (pH*cstar/dp + phib/Qbhat)*(Qbhat-1);
     kappa0e*(exp(kappa1e*phie)-1) + kappa0e*kappa1e*exp(kappa1e*phie)*phie = 1 - etaw/etab;
     kappa0b*(exp(kappa1b*phib)-1) + kappa0b*kappa1b*exp(kappa1b*phib)*phib = (pH/dp)*cstar*(Qbhat - 1)/((pH/dp)*cstar*Qbhat + phib/Qbhat);
  G = (pH/dp)*(etae*bstar/(etab*(1+rd))) + (1+pH/dp)*cstar - rhobhat 
    + (kappa0e*(exp(kappa1e*phie)-1) - 1 + etaw/etab)*phie 
    + (kappa0b*(exp(kappa1b*phib)-1) - 1 + 1/Qbhat)*phib;
//    + (kappa0e*kappa1e*exp(kappa1e*phie) - 1 + etaw/etab)*phie 
//    + (kappa0b*kappa1b*exp(kappa1b*phib) - 1 + 1/Qbhat)*phib;
  I = M/G;
  rhobhat = pH*q*R*(1+s)/(etab*(1+rd))-1;  // note that there is no Rstar
  ve = beta*(((rK(+1)+(1-delta)*q(+1))/q)*(mU(+1)/mU)*(lambdae+(1-lambdae)*(1+retilde(+1))*ve(+1)));
  vb = beta*(((rK(+1)+(1-delta)*q(+1))/q)*(mU(+1)/mU)*(lambdab+(1-lambdab)*(1+ratilde(+1))*vb(+1)));
  vw = beta*((rK(+1)+(1-delta)*q(+1))/q)*(mU(+1)/mU);
  (1+e_I(+1))*ve(+1)*etae = ve(+1);
  (1+e_I(+1))*vb(+1)*etab = vb(+1);
  (1+e_I(+1))*vw(+1)*etaw = vw(+1);
  1 + retilde = etae*(pH/dp)*(bstar/G)*(A + N)*(1+e_I(+1))/N;
  1 + ratilde = ((1+rd)*(pH/dp)*cstar*(A + N)/(A*G))*(1 + (PI-etaw*phie)*(1+e_I(+1)-1/etab)/((pH/dp)*cstar + phib/Qbhat));
  M = A + N;
  omega = N/M;
  M = (1+rd)*(pH/dp)*(M(-1)/G)*((rK + (1-delta)*q(+1))/q)*((1-lambdae)*etae*(bstar/(1+rd))*(1+e_I) + (1-lambdab)*cstar*(1 + (PI - etaw*phie)*(1+e_I-1/etab)/((pH/dp)*cstar + phib/Qbhat)));
  omega = (1-lambdae)*etae*(bstar/(1+rd))*(1+e_I)/((1-lambdae)*etae*(bstar/(1+rd))*(1+e_I)+(1-lambdab)*cstar*(1 + (PI - etaw*phie)*(1+e_I-1/etab)/((pH/dp)*cstar + phib/Qbhat)));
end;

initval; //  this is the analytical steady state
  Z = 0; 
  Zc = 1;
  e_ = 0;
  e_I = 0;
  ee_I = 0;
  e_sigma = 0;
  sigma_I = 1;
  rhobhat = (pH/dp)*(1 - ((1-lambdae)/beta))*bstar + (pH/dp)*(1 - ((1-lambdab)/beta))*cstar + cstar;
  q = (1+cstar)/(pH*R);
//  q = 1;
  //  return to capital
  rK = q*(1/beta-1+delta);
  //  real wage
  w=(1-alpha)*(rK/alpha)^(-alpha/(1-alpha));
  //  capital stock
  K =((1-alpha)/ksi)^(1/(sigma+fi))*(rK/alpha)^(-(alpha+fi)/((1-alpha)*(sigma+fi)))*(rK/alpha-delta)^(-sigma/(sigma+fi));
  //  labour 
  L=K*(rK/alpha)^(1/(1-alpha));
  //  output
  Y = K*rK/alpha;
  //  investments
  I = delta*K;
  //  insider wealth, J = A + N
//  C = Y - (1+cstar)*I; // -I*chi/pH;
  C = Y - I;
  //  consumption of entrepreneurs
  mU =C^(-sigma);
  G = (1/beta)*(pH/dp)*((1-lambdae)*bstar+(1-lambdab)*cstar);
  M = G*I;
  A = (1-lambdab)*cstar*M/((1-lambdab)*cstar+(1-lambdae)*bstar);
  N = (1-lambdae)*bstar*M/((1-lambdab)*cstar+(1-lambdae)*bstar);
  etab = 1;
  etae = 1;
  etaw = 1;
  ratilde = (pH/dp)*((1-lambdab)*cstar+(1-lambdae)*bstar)/((1-lambdab)*G)-1;
  retilde = etae*(pH/dp)*(bstar/G)*(A + N)/N-1;
  ve = lambdae/(1-(1-lambdae)*(1+retilde));
  vb = lambdab/(1-(1-lambdab)*(1+ratilde));
  vw = 1;
  Qbhat = 1;
  rd = 0;
  phie = 0;
  phib = 0;
  PI = pH*q*R/(1+rd) - pH*etae*bstar/(dp*(1+rd));
  omega = N/M;
end;

//   shocks;
//     var e_ = 0.01^2;
//    var e_c = 0;
//     var ee_I = 0.01^2;
//    var e_sigma = 0.01^2;
//   end;
vcov = [0.00005184 0 0 0 ,
          0 0 0 0 ,
          0 0 0 0,
          0 0 0 0 ];

vcov = [0.0000168 0 0 0 ;
          0 0 0 0 ;
          0 0 0.0001 0;
          0 0 0 0.0001 ];
order = 3;
// resid(1);
// steady(solve_algo = 2);
// check;
// for i=1:length(oo_.steady_state);
//   fprintf(1,'%-s %10.6f\n',M_.endo_names(i,:),oo_.steady_state(i));
// end;
// stoch_simul(order=3, irf=400) ;
