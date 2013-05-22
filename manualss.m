% manualss.m
% Calculates the steady state "manually", ie using fsolve separately
%
% (c) Antti Ripatti, 2013-
%
%% Parameter values (copied from the mod file)
beta = 1/1.02^0.25; 
alpha = 1/3;
delta = 0.025;
ksi = 2; %  2
fi = 0.1; %  3
rho = 0.979; 
rhoc = 0.6;
sigma = 1;
sigma2 = ((0.007/(1-alpha)))^2;

rr=1/beta-1;
rho_eI=.9;
b=0;

%  Financial sector calibration
%  Based on the following data
%  riskless rate: 1+rbar = 
rbar = 1.02^0.25; 
%  excess return on banks capital: 1+rabar = 1.20^(1/4)
rabar = 1.13^0.25;
rabar = 1.20^0.25;
%rabar = 1.10^0.25;
%  excess return on 1+rebar = (1.065-1.02)^0.25
rebar = (1.065-rbar^4+1)^0.25;
%rebar = (1.08-rbar^4+1)^0.25;
%  Entrepreneurs capital ratio N/I = 0.3
CRF = 0.3;
CRF = 0.45;
%  Banks capital ratio CRB= A/(I-N) = 0.08
CRB = 0.08;
%  Banks  monitoring costs 0.03 relative to assets (per annum / 4)
CORB = 0.01/4;
CORB = 0.006/4;
%  calculated parameter values
lambdab = beta/rabar;
lambdae = beta/rebar;
lambdab = 1-lambdab;
lambdae = 1-lambdae;
pH = 0.95;
rdp = CORB/(CRB*rabar);
dp = pH*rdp;                 
% dp = .5;
pL = pH - dp;
R = 1/pH;
bstar = 0.00815075;
%cstar=0;
cstar = 0.000825;
%cstar=.04;
% cstar = 0.825;
kappa0e = 1;
kappa0b = 1;
kappa1e = 100;
kappa1b = 100;


% kappa1e = 10;
% kappa1b = 5;
% kappa0e = 0.24;
% kappa0b = 0.002;
% kappa1e = 7;
% kappa1b = 11;
s = ((pH/dp)*(1 - ((1-lambdae)/beta))*bstar + (pH/dp)*(1 - ((1-lambdab)/beta))*cstar + cstar - cstar)/(1+cstar);
rhobhatbar =   (pH/dp)*(1 - ((1-lambdae)/beta))*bstar + (pH/dp)*(1 - ((1-lambdab)/beta))*cstar + cstar;
Rstar = R*(1+rhobhatbar)/(1+cstar);
%Rstar=R;
 
%% Parameter vector
p = [sigma2, sigma, delta, beta, alpha, rho, rhoc, ksi, fi, rho_eI, b, pH, dp, bstar, cstar, kappa0e, kappa0b, kappa1e, kappa1b, R, lambdab, lambdae, s, Rstar, rhobhatbar]';
%% Analytical steady-state
  Z = 0; 
  Zc = 1;
  e_ = 0;
  e_I = 0;
  ee_I = 0;
  e_sigma = 0;
  sigma_I = 1;
  rhobhat = (pH/dp)*(1 - ((1-lambdae)/beta))*bstar + (pH/dp)*(1 - ((1-lambdab)/beta))*cstar + cstar;
  q = (1+cstar)/(pH*R);
%  q = 1;
 
  rK = q*(1/beta-1+delta);
  %  real wage
  w=(1-alpha)*(rK/alpha)^(-alpha/(1-alpha));
  %  capital stock
  %K =((1-alpha)/ksi)^(1/(sigma+fi))*(rK/alpha)^(-(alpha+fi)/((1-alpha)*(sigma+fi)))*(rK/alpha-delta)^(-sigma/(sigma+fi));
  K =((1-alpha)/ksi)^(1/(sigma+fi))*(rK/alpha)^(-(alpha+fi)/((1-alpha)*(sigma+fi)))*(rK/alpha-delta*(1+cstar)/(pH*R))^(-sigma/(sigma+fi));
  %  labour 
  L=K*(rK/alpha)^(1/(1-alpha));
  %  output
  Y = K*rK/alpha;
  %  investments
  I = delta*K;
  %  insider wealth, J = A + N
 C = Y - (1+cstar)*I; 

 % C = Y - I;
  
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

init_y = [mU, Y, I, C, K, L, Z, Zc, rK, w, e_I, sigma_I, PI, q, phie, phib, ve, vb, vw, Qbhat, ratilde, retilde, rhobhat, rd, etae, etab, etaw, G, A, N, M, omega]';
%% Computations
disp('Residuals of the steady-state equations:');
rbcOutEq_f(p,init_y)
opt = optimset('Jacobian','on','Display','iter','TolFun',1e-14);
y = fsolve(@(y) rbcOutEq_fsolve(p,y),init_y, opt);
sslabels = {'mU','Y','I','C','K','L','Z','Zc','rK','w','e_I','sigma_I','PI','q','phie','phib','ve','vb','vw','Qbhat','ratilde','retilde','rhobhat','rd','etae','etab','etaw','G','A','N','M','omega'}';
disp('Calculated steady-state');
disp('                     iterated      initial');
horzcat(sslabels,num2cell([y init_y]))
