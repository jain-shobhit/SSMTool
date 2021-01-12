function varargout = coco_read_complementary(fid, runid, lab, varargin)
% COCO_READ_COMPLEMENTARY   Extract complementary solution information from file
%
% VARARGIN  = { [NAME...] }
% NAME      = 'data' | 'chart' | 'uidx' | 'lidx' | 'vidx'
%
% If no names are given, {NAME...} = {'data' 'chart'}.

% Copyright (C) Harry Dankowicz
% $Id: coco_read_solution.m 2972 2017-01-10 19:44:34Z hdankowicz $

persistent cache

if nargin==0
  % clear cache and exit
  cache = [];
  return
end

s = coco_stream(varargin{:});

[func, chart] = coco_read_solution(fid, runid, lab, 'data', 'chart');

sfname = sprintf('sol%d.mat', lab);
sfname = coco_fname(runid, sfname);

lst = dir(sfname);
if isempty(cache) || ~strcmp(sfname,cache.sfname) || cache.datenum<lst.datenum
  cache.sfname  = sfname;
  cache.datenum = lst.datenum;
  cache.vars    = who('-file', sfname);
  cache.cont    = load(sfname);
end

nout = max(1,nargout);
vars = cache.vars;
if any(strcmp('version', vars))
  version = cache.cont.version;
  switch version{1}
    case 1
      [varargout{1:nout}] = read_sol(func, chart, version{2}, s);
    otherwise
      error('%s: unknown format of solution file', mfilename);
  end
else
  error('%s: unknown format of solution file', mfilename);
end

end


function varargout = read_sol(func, chart, format, s)

list = {};
while ~isempty(s)
  oname = s.peek;
  switch lower(oname)
    case {'data' 'chart' 'uidx' 'lidx' 'vidx'}
      list = [ list s.get ]; %#ok<AGROW>
    otherwise
      if ischar(oname)
        error('%s: unrecognised property ''%s''', mfilename, oname);
      else
        error('%s: options must be strings', mfilename);
      end
  end
end
if isempty(list)
  list = {'data' 'chart'};
end

assert(nargout<=numel(list), '%s: too few input arguments', mfilename);

switch format
  case 'full'
    soldata  = func.data;
    chart.v  = chart.x(func.vidx);
    chart.tv = chart.t(func.vidx);
  case 'reduced'
    error('not supported');
end

for i=1:nargout
  switch list{i}
    case 'data'
      varargout{i} = soldata; %#ok<AGROW>
    case 'chart'
      varargout{i} = chart; %#ok<AGROW>
    case 'uidx'
      varargout{i} = func.u_idx; %#ok<AGROW>
    case 'lidx'
      varargout{i} = func.l_idx; %#ok<AGROW>
    case 'vidx'
      varargout{i} = func.v_idx; %#ok<AGROW>
  end
end

end
