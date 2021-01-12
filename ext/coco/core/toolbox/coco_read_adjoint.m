function varargout = coco_read_adjoint(fid, runid, lab, varargin)
% COCO_READ_ADJOINT   Extract adjoint information from file
%
% VARARGIN  = { [NAME...] }
% NAME      = 'data' | 'chart' | 'lidx'
%
% If no names are given, {NAME...} = {'data' 'chart'}.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_read_adjoint.m 2839 2015-03-05 17:09:01Z fschild $

persistent cache

if nargin==0
  % clear cache and exit
  cache = [];
  return
end

s = coco_stream(varargin{:});

[data, chart] = coco_read_solution('coco_adjoint', runid, lab, 'data', ...
  'chart');

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
      [varargout{1:nout}] = read_sol(data, chart, cache.cont.adata, fid, version{2}, s);
    otherwise
      error('%s: unknown format of solution file', mfilename);
  end
else
  error('%s: unknown format of solution file', mfilename);
end

end

function varargout = read_sol(data, chart, adata, fid, format, s)

idx = find(strcmpi(fid, adata(:,1)));
assert(~isempty(idx), '%s: no such adjoint', mfilename);

list = {};
while ~isempty(s)
  oname = s.peek;
  switch lower(oname)
    case {'data' 'chart' 'lidx'}
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
    adjoint = data.adjoint;
    soldata = adjoint.funcs(idx).data;
    lidx    = data.l_idx(coco_slot_data(fid, adata));
    chart.x = chart.x(lidx);
    chart.t = chart.t(lidx);
  case 'reduced'
    error('not supported');
end

for i=1:nargout
  switch list{i}
    case 'data'
      varargout{i} = soldata; %#ok<AGROW>
    case 'chart'
      varargout{i} = chart; %#ok<AGROW>
    case 'lidx'
      varargout{i} = lidx; %#ok<AGROW>
  end
end

end
