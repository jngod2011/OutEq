//% rbcOutEq.mod
//% Financial contract with macro investment shocks and outside equity
//% 
//%  (c) Antti Ripatti and Markus Haavio, 2011-2013
//%
//%  builds on ./code/GK/buffer.mod and RBC.mod by Haavio and Ripatti
//%
//%  CHANGE LOG:
//%  Markus 20.5.2013
//%
//% PI = pH*q*Rstar/(1+rd) - pH*etae*bstar/(dp*(1+rd))
//% aiemmin
//% PI = pH*q*R/(1+rd) - pH*etae*bstar/(dp*(1+rd))
//% Siis R =Rstar (missä Rstar = R*(1+rhobhatbar)/(1+cstar);
//% ii)    rhobhat = pH*q*Rstar/(etab*(1+rd))-1;
//% aiemmin
//% rhobhat = pH*q*R*(1+s)/(etab*(1+rd))-1;
//% Siis  R*(1+s)=> Rstar
//% Alkoarvoihin korjaus iii)
//% C = Y - I*(1+cstar)/pH*R;
//% Koodeista löytyi vielä bugeja.  Resurssirajoitteen 
//% C+I*(1+cstar+osakerahoitukseen liittyvät kustannukset )  =Y
//% kaikkia implikaatioita ei ollut viety systemaattisesti dynaamisiin 
//% mallikoodeihin (HT + rbc) ja steady statea kuvaaviin yhtälöihin 
//% (jälleen sekä HT että rbc): alkuarvot eivät ihan täsmänneet kaikkien 
//% yhtälöiden kohdalla.   Nyt sain mielestäni bugit korjattua: 
//% residuaalit menevät nollille. 
//% 


//%  Define variables, parameters, shocks, etc.
var 
  mU //%  marginal utility of consumption
  Y  //%  output
  I //%  investments
  C //%  consumption
  K //%  capital stock
  L //%  hours/employment, here =1 
  Z //%  technology shock
  Zc //%  preference shock
  rK //%  capital return
  w //%  real wage
  e_I //% investment shock 
  sigma_I //% uncertainty shock (stochastic volatility)
  PI //% Pledgeable income 
  q //% Tobin's q, price of capital
  phie //% outside equity investment ratio, firms
  phib //% outside equity investment ratio, banks
  //% marginal value of capital 
  ve //% entrepreneurial
  vb //% banks 
  vw //% worker-owned
  Qbhat //% risk-adjusted price of banks' outside equity
  ratilde //% return on banks capital
  retilde //% return on entrepreneurial capital
  rhobhat //% risk-adjusted NPB of the investment project
  rd //% deposit rate
  etae //% risk adjustment 
  etab
  etaw 
  G //% inverse of the leverage
  A //% bankers' capital
  N //% firms' capital
  M //% informed capital
  omega //% firms capital relative to informed capital
;

//%  exogenous shock
varexo e_, e_c, ee_I, e_sigma;
//%  list of parameters
parameters 
  sigma2 //%  shock variance
  sigma //%  risk aversion
  delta //%  depreciation rate
  beta //%  discount rate
  alpha //%  capital share
  rho //%  shock persistence
  rhoc //%preference shock persistence
  ksi //%  marginal disutility of labour 
  fi //%  inverse of frich elasticity of labour supply
  rho_eI
  b //%  habit persistence
  pH
  dp 
  bstar //%  private benefit
  cstar //%  banks monitoring cost
  kappa0e
  kappa0b
  kappa1e
  kappa1b //% parameters of the efficiency costs, ie kappa0*phi^kappa1
  R //% size of the cake
  lambdab //% bankers' survival probability
  lambdae //% firms' survival probability
  s //% investment subsidy
  Rstar
  rhobhatbar 
;


//% //%  Parameter values, King-Rebelo (2000)
beta = 1/1.02^0.25; 
alpha = 1/3;
delta = 0.025;
ksi = 2; //%  2
fi = 0.1; //%  3
rho = 0.979; 
rhoc = 0.6;
sigma = 1;
sigma2 = ((0.007/(1-alpha)))^2;

//%  Hansen (1985) parameterization
//%beta = 0.99;
//%alpha = 0.36;
//%delta = 0.0025;
//%rho = 0.95;
//%sigma = 1;


rr=1/beta-1;
rho_eI=.9;
b=0;

//%  Financial sector calibration
//%  Based on the following data
//%  riskless rate: 1+rbar = 
rbar = 1.02^0.25; 
//%  excess return on banks capital: 1+rabar = 1.20^(1/4)
rabar = 1.13^0.25;
rabar = 1.20^0.25;
//%rabar = 1.10^0.25;
//%  excess return on 1+rebar = (1.065-1.02)^0.25
rebar = (1.065-rbar^4+1)^0.25;
//%rebar = (1.08-rbar^4+1)^0.25;
//%  Entrepreneurs capital ratio N/I = 0.3
CRF = 0.3;
CRF = 0.45;
//%  Banks capital ratio CRB= A/(I-N) = 0.08
CRB = 0.08;
//%  Banks  monitoring costs 0.03 relative to assets (per annum / 4)
CORB = 0.01/4;
CORB = 0.006/4;
//%  calculated parameter values
lambdab = beta/rabar;
lambdae = beta/rebar;
lambdab = 1-lambdab;
lambdae = 1-lambdae;
pH = 0.95;
rdp = CORB/(CRB*rabar);
dp = pH*rdp;                 
//% dp = .5;
pL = pH - dp;
R = 1/pH;
bstar = 0.00815075;
//%cstar=0;
cstar = 0.000825;
//%cstar=.04;
//% cstar = 0.825;
kappa0e = 1;
kappa0b = 1;
kappa1e = 100;
kappa1b = 5;


//% kappa1e = 10;
//% kappa1b = 5;
//% kappa0e = 0.24;
//% kappa0b = 0.002;
//% kappa1e = 7;
//% kappa1b = 11;
s = ((pH/dp)*(1 - ((1-lambdae)/beta))*bstar + (pH/dp)*(1 - ((1-lambdab)/beta))*cstar + cstar - cstar)/(1+cstar);
rhobhatbar =   (pH/dp)*(1 - ((1-lambdae)/beta))*bstar + (pH/dp)*(1 - ((1-lambdab)/beta))*cstar + cstar;
Rstar = R*(1+rhobhatbar)/(1+cstar);
//%Rstar=R;

//% 
//%  Model code: start equation at line x1 in order to be able to map the
//%  equation numbers that dynare gives with the code below
model; 
  mU = Zc*((C-b*C(-1))/(1-b))^(-sigma); 
  mU = beta*mU(+1)*((q(+1)*(1-delta) + rK(+1))/q); 
  Z = rho*Z(-1) + e_; 
  log(Zc)=rhoc*log(Zc(-1)) + e_c; 
  Y = K^alpha*(exp(Z)*L)^(1-alpha); 
  rK/w = (alpha/(1-alpha))*L/K; 
  rK = alpha*Y/K;
  K = (1-delta)*K(-1) + pH*R*I(-1)*(1+e_I); 
  L = (w*mU/ksi)^(1/fi);
  C + (1+cstar+kappa0b*(exp(kappa1b*phib)-1)*phib+kappa0e*(exp(kappa1e*phie)-1)*phie)*I = Y; //% 10 
  sigma_I = (1-rho_eI) + rho_eI*sigma_I(-1) + e_sigma; 
  e_I=  sigma_I*ee_I;
  rd = 0;
  PI = pH*q*Rstar/(1+rd) - pH*etae*bstar/(dp*(1+rd)); //% (15)
  (1/etaw - 1/etab)*(PI - etaw*phie) = (pH*cstar/dp + phib/Qbhat)*(Qbhat-1);
  phie = 0; //% only one instrument to the planner; it is phib
 //%    kappa0e*(exp(kappa1e*phie)-1) + kappa0e*kappa1e*exp(kappa1e*phie)*phie = 1 - etaw/etab;
 //%    kappa0b*(exp(kappa1b*phib)-1) + kappa0b*kappa1b*exp(kappa1b*phib)*phib = (pH/dp)*cstar*(Qbhat - 1)/((pH/dp)*cstar*Qbhat + phib/Qbhat);
  G = (pH/dp)*(etae*bstar/(etab*(1+rd))) + (1+pH/dp)*cstar - rhobhat 
    + (kappa0e*(exp(kappa1e*phie)-1) - 1 + etaw/etab)*phie 
    + (kappa0b*(exp(kappa1b*phib)-1) - 1 + 1/Qbhat)*phib;
  I = M/G;
  rhobhat = pH*q*Rstar/(etab*(1+rd))-1;  //% note that there is no Rstar
  ve = beta*(((rK(+1)+(1-delta)*q(+1))/q)*(mU(+1)/mU)*(lambdae+(1-lambdae)*(1+retilde(+1))*ve(+1)));
  vb = beta*(((rK(+1)+(1-delta)*q(+1))/q)*(mU(+1)/mU)*(lambdab+(1-lambdab)*(1+ratilde(+1))*vb(+1)));
  vw = beta*((rK(+1)+(1-delta)*q(+1))/q)*(mU(+1)/mU);
  etae*(1+e_I(+1))*beta*(((rK(+1)+(1-delta)*q(+1))/q)*(mU(+1)/mU)*(lambdae+(1-lambdae)*(1+retilde(+1))*ve(+1))) = ve;
  etab*(1+e_I(+1))*beta*(((rK(+1)+(1-delta)*q(+1))/q)*(mU(+1)/mU)*(lambdab+(1-lambdab)*(1+ratilde(+1))*vb(+1))) = vb;
  etaw*(1+e_I(+1))*beta*((rK(+1)+(1-delta)*q(+1))/q)*(mU(+1)/mU) = vw;
  1 + retilde = etae*(pH/dp)*(bstar/G)*M*(1+e_I(+1))/N;
  1 + ratilde = ((1+rd)*(pH/dp)*cstar*M/(A*G))*(1 + (PI-etaw*phie)*(1+e_I(+1)-1/etab)/((pH/dp)*cstar + phib/Qbhat));
  N=omega*M;
  A=(1-omega)*M;
  //%M = A + N;
  //%omega = N/M;
  M = (1+rd(-1))*(pH/dp)*(M(-1)/G(-1))*((rK + (1-delta)*q)/q(-1))*((1-lambdae)*etae(-1)*(bstar/(1+rd(-1)))*(1+e_I) + (1-lambdab)*cstar*(1 + (PI(-1) - etaw(-1)*phie(-1))*(1+e_I-1/etab(-1))/((pH/dp)*cstar + phib(-1)/Qbhat(-1))));
  omega = (1-lambdae)*etae(-1)*(bstar/(1+rd(-1)))*(1+e_I)/((1-lambdae)*etae(-1)*(bstar/(1+rd(-1)))*(1+e_I)+(1-lambdab)*cstar*(1 + (PI(-1) - etaw(-1)*phie(-1))*(1+e_I-1/etab(-1))/((pH/dp)*cstar + phib(-1)/Qbhat(-1))));
end;

initval; //%  this is the analytical steady state
//%steady_state_model
  Z = 0; 
  Zc = 1;
  e_ = 0;
  e_I = 0;
  ee_I = 0;
  e_sigma = 0;
  sigma_I = 1;
  rhobhat = (pH/dp)*(1 - ((1-lambdae)/beta))*bstar + (pH/dp)*(1 - ((1-lambdab)/beta))*cstar + cstar;
  q = (1+cstar)/(pH*R);
//%  q = 1;
 
  rK = q*(1/beta-1+delta);
  //%  real wage
  w=(1-alpha)*(rK/alpha)^(-alpha/(1-alpha));
  //%  capital stock
  //%K =((1-alpha)/ksi)^(1/(sigma+fi))*(rK/alpha)^(-(alpha+fi)/((1-alpha)*(sigma+fi)))*(rK/alpha-delta)^(-sigma/(sigma+fi));
  K =((1-alpha)/ksi)^(1/(sigma+fi))*(rK/alpha)^(-(alpha+fi)/((1-alpha)*(sigma+fi)))*(rK/alpha-delta*(1+cstar)/(pH*R))^(-sigma/(sigma+fi));
  //%  labour 
  L=K*(rK/alpha)^(1/(1-alpha));
  //%  output
  Y = K*rK/alpha;
  //%  investments
  I = delta*K;
  //%  insider wealth, J = A + N
 C = Y - (1+cstar)*I; 

 //% C = Y - I;
  
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
  PI = pH*q*Rstar/(1+rd) - pH*etae*bstar/(dp*(1+rd));
  omega = N/M;
end;

//%steady; 

//%resid;

//%   shocks;
//%     var e_ = 0.01^2;
//%    var e_c = 0;
//%     var ee_I = 0.01^2;
//%var e_sigma = 0.01
//%  %  var e_sigma = 0.01^2;
//%   end;


vcov = [0.0000168 0 0 0 ;
          0 0 0 0 ;
          0 0 .00008 0;
          0 0 0 0.0 ];
order = 3;

planner_objective (((C)/(1))^(1-sigma))/(1-sigma)-(ksi/(1+fi))*L^(1+fi);
planner_discount beta; 

//% resid(1);
 //%steady(solve_algo = 1);
//% check;
//% for i=1:length(oo_.steady_state);
//%   fprintf(1,'%-s %10.6f\n',M_.endo_names(i,:),oo_.steady_state(i));
//% end;
//% stoch_simul(order=3, irf=400) ;
