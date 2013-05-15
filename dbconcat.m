function db = dbconcat(d1,d2)
% concatenate tseries of d1 and d2 into db
% Input
%  d1 db structure
%  d2 db structure
%
% Output
%  db db structure
% Example:
% db = dbconcat(d1,d2);
%
% Created
%   Antti Ripatti, 7.6.2012
% 
% (c) Antti Ripatti, 2012-
%
if nargin == 0
  error('You must provide at least one dbase object file');
end;
fnv = fieldnames(d1);
nv = length(fnv);
db = struct();
for i = 1:nv;
  if isfield(d2,fnv{i}) % check if the field exists in the second dbase
    if length(d1.(fnv{i})) == length(d2.(fnv{i}))
      data = [get(d1.(fnv{i}),'data') get(d2.(fnv{i}),'data')];
      db.(fnv{i}) = tseries(range(d1.(fnv{i})),data);
    else
      warning('no:no','tseries %s are of different range in dbases',fnv{i});
    end;
  else
    warning('no:no','tseries %s is not present in the second database.',fnv{i});
    db.(fnv{i}) = d1.(fnv{i});
  end;
end;
