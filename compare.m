% compare.m 
% plots the basic simulation results of the outside equity model
%
% (c) Antti Ripatti, 2013 -
% 
%% addpath and clear workspace
%addpath('c:\hy-data\ajripatt\mallit\matlab\iris');
%addpath('c:\hy-data\ajripatt\\mallit\matlab\bofiris');
%irisstartup
%addpath('c:\hy-data\AJRIPATT\mallit\matlab\dynare\dynare++');
clear('all');close('all');clear('struct');
%% User defined control data
rundynarepp = true(); % run mod files using dynare++ 
myMod1 = 'rbcOutEq.mat';
myMod2 = 'rbc.mat';
korder = 3; % order of the Taylor approximation
nIrf = 40; % length of the irfs
plotName = 'OutEq_Vola_Shock'; 
legendtxt = {'RBC with H&T and outside equity'; 'RBC'};
% list of variables to be plotted
% list of variables where percent from the steady-state is computed
myList = {'Y','I','C','A','N','q','ve','vb','vw','w','L','etab','etae','etaw','phib','phie','Z'};
myListLabels = {'"Output" Y','"Investments" I','"Consumption" C','"Bank capital" A',...
  '"Entrepreneurial capital" N','"$q$" q','"$v^e$" ve',...
  '"$v^b$" vb','"$v^w$" vw','"Real wages" w','"Hours" L',...
  '"$\eta^b$" etab','"$\eta^e$" etae','"$\eta^w$" etaw',...
  '"Technology shock, $Z_t$" Z'};
% list of variables, where the difference is computed
myDiffList = {'retilde','ratilde','e_I','etab','etae','etaw','phib','phie'};
myDiffListLabels = {'"$\tilde r^e$" retilde','"$\tilde r^a$" ratilde',...
  '"Investment shock" e_I','"$\eta^b$" etab','"$\eta^e$" etae','"$\eta^w$" etaw','"$\phi_t^b$" phib','"$\phi_t^e$" phie',};
verysmallnumber = 1e-12;
%% run dynare++ if needed
if rundynarepp
  disp('Running dynare++ ...');
  [status,result] = system('dynare++ --threads 6 --no-irfs rbc.mod','-echo');
  if ~isempty(result)
    fprintf(1,'dynare++ gives the following error message:_\n');
    error('fake','%s\n',result);
  end;
  [status,result] = system('dynare++ --threads 6 --no-irfs rbcOutEq.mod','-echo');
  if ~isempty(result)
    fprintf(1,'dynare++ gives the following error message while running rbcOutEq.mod\n');
%    error('%s\n',result);
  end;
  disp('...passed!');
end;
%% Load first model
load(myMod1);
nShocks = length(dyn_vcov_exo);
%% Then the simulations using dynare++'s own dynare_simul function
%  investment shock
ex_=zeros(nShocks,nIrf);
%ex_(dyn_i_ee_I,1) = -0.1;
%ex_(dyn_i_e_,1) = 0.01;
ex_(dyn_i_e_sigma,1) = 1;
irf1=dynare_simul(myMod1,ex_);
irf1(irf1 < verysmallnumber) = 0; dyn_ss(dyn_ss < verysmallnumber) = 0;
if any(any(isnan(irf1)))
  error('Explosive system');
end;
db1 = irfpp2db([],irf1,dyn_ss,'eI',cellstr(dyn_vars));
db1diff = irfpp2db([],irf1,dyn_ss,'eI',cellstr(dyn_vars),'relative',false);
ss1 = dyn_ss;
disp('Holmstrom -- Tirole in Macro with Outside Equity-----------------');
disp('          Deterministic Stoch. fix-point Uncoditional mean');
horzcat(cellstr(dyn_vars),num2cell([dyn_steady_states dyn_ss dyn_mean]))
%% Next scenario
clear dyn* ex_;
load(myMod2);
nShocks = length(dyn_vcov_exo);
disp('Basic RBC---------------');
disp('          Deterministic Stoch. fix-point Uncoditional mean');
horzcat(cellstr(dyn_vars),num2cell([dyn_steady_states dyn_ss dyn_mean]))
%%  ex post capitalization of investment shock
ex_=zeros(nShocks,nIrf);
%ex_(dyn_i_ee_I,1) = -0.1;
%ex_(dyn_i_e_,1) = 0.01;
ex_(dyn_i_e_sigma,1) = 1;
irf3=dynare_simul(myMod2,ex_);
irf3(irf3 < verysmallnumber) = 0; dyn_ss(dyn_ss < verysmallnumber) = 0;
db2 = irfpp2db([],irf3,dyn_ss,'eI',cellstr(dyn_vars));
db2diff = irfpp2db([],irf3,dyn_ss,'eI',cellstr(dyn_vars),'relative',false);
% concatenate the dbases of the above simulations
db = dbconcat(db1.eI,db2.eI); 
dbdiff = dbconcat(db1diff.eI,db2diff.eI); 
% plot stuff
fig = dbplot(db,myListLabels,range(db.(myList{1})),'tight',true,'zeroline',true,'interpreter=','latex'); % ,'subplot',[3 4]
fig = updateStyleLegend(fig,plotName,legendtxt);
fig = dbplot(dbdiff,myDiffListLabels,range(db.(myList{1})),'tight',true,'zeroline',true,'interpreter=','latex'); % ','subplot',[3 4]
fig = updateStyleLegend(fig,[plotName 'Diff'],legendtxt);
