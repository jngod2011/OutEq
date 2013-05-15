function db = irfpp2db(mydb,irf,ss,sname,vnames,varargin)
% store dynare++ irfs to Iris dbase structure
% Input
%   mydb dbase object to which the new shocks are added
%   irf nxT matrix of impulse responses 
%   sname string name of shock
%   vnames n cell vector of variable names
% Variable input
%  'relative' boolean (default true) to determine whether we compute
%  percent from the steady-state
%
% Output
%   db Iris toolbox dbase structure organized such that 
%       db.sname.vnames{i}
% Created
%   Antti Ripatti, 7.6.2012
%
% (c) Antti Ripatti, 2012-
%
if isempty(mydb)
  db = struct();
else
  db = mydb;
end;
if nargin < 4
    error('Too few parameters. See: help irfpp2db');
end;
default = {...
    'relative',true,@islogical,...
};
options = passvalopt(default,varargin{:});
nvnames = length(vnames);
irflen = size(irf,2);
if options.relative
  for j = 1:size(irf,1);
    if ss(j) ~= 0;
      relss(j,:) = 100*(irf(j,:)./ss(j)-1);
    else
      warning('no:no','Steady-state of %s is zero. The original path will be returned.',vnames{j});
      relss(j,:) = irf(j,:);
    end;
  end;
else
  for j = 1:size(irf,1);
    relss(j,:) = irf(j,:)-ss(j);
  end;
end;
if ~iscellstr(vnames)
  disp('the last argument is not a cell vector of strings! I am trying to transform it.');
  vnames = cellstr(vnames);
end;
if sum(isspace(sname))
  warning('antti:stuff','shock name contains spaces, I am getting rid of them!');
  sname = strrep(sname,' ','');
end;
for i=1:nvnames;
  shockname = sname;
  varname = vnames{i};
  scomment = sprintf('%s to %s',varname, shockname);
  db.(shockname).(varname) = tseries(1:irflen,relss(i,:),scomment);
end;