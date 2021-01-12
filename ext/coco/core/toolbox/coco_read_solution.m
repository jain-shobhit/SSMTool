function varargout = coco_read_solution(varargin)
% COCO_READ_SOLUTION   Extract solution information from file
%
% VARARGIN   = { [FID] RUNID LAB [OPTIONS] [NAME...] }
% OPTIONS   = '-no-extract-xp'
% NAME      = 'data' | 'run' | 'chart' | 'chart1' | 'uidx' | 'xidx'
%
% If no names are given, {NAME...} = {'data' 'chart' 'run' 'chart1'}.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_read_solution.m 2972 2017-01-10 19:44:34Z hdankowicz $

persistent cache

if nargin==0
  % clear cache and exit
  cache = [];
  return
end

s = coco_stream(varargin{:});
[fid runid] = s.get;
if ischar(runid) || iscellstr(runid)
  lab = s.get;
else
  [fid runid lab] = coco_deal('', fid, runid);
end

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
      % varargout = {soldata chart run chart1}
      %           | defined by varargin
      [varargout{1:nout}] = read_sol_v1(cache, fid, version{2}, s);
    otherwise
      error('%s: unknown format of solution file', mfilename);
  end
else
  % varargout = {soldata chart run chart1}
  [varargout{1:nout}] = read_sol_v0(cache, fid);
end

end

function varargout = read_sol_v1(cache, fid, format, s)

extract_xp = true;

list = {};
while ~isempty(s)
  oname = s.peek;
  switch lower(oname)
    case 'no-extract-xp'
      fprintf(2, '%s: option ''no-extract-xp''deprecated, use ''-no-extract-xp'' instead', mfilename);
      extract_xp = false;
      s.skip;
    case '-no-extract-xp'
      extract_xp = false;
      s.skip;
    case {'data' 'run' 'chart' 'chart1' 'uidx' 'xidx'}
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
  list = {'data' 'chart' 'run' 'chart1'};
end

assert(nargout<=numel(list), '%s: too few input arguments', mfilename);

switch format
  case 'full'
    % save(fname, 'version', 'run', 'data', ['chart', 'chart1', 'fdata'])
    run    = cache.cont.run;
    data   = cache.cont.data;
    chart  = cache.cont.chart;
    chart1 = cache.cont.chart1;
    fdata  = cache.cont.fdata;
    if isempty(fid)
      soldata = data;
    else
      if isempty(data)
        soldata = data;
      else
        soldata = coco_slot_data(fid, data);
      end
      if extract_xp
        xidx = coco_slot_data(fid, fdata);
        if ~isempty(xidx)
          chart.x  = chart.x(xidx);
          if isfield(chart, 't')
            chart.t = chart.t(xidx);
          end
          if isfield(chart, 'TS')
            chart.TS = chart.TS(xidx,:);
          end
          
          if isfield(chart1, 'x')
            chart1.x = chart1.x(xidx);
          end
          if isfield(chart1, 't')
            chart1.t = chart1.t(xidx);
          end
          if isfield(chart1, 'TS')
            chart1.TS = chart1.TS(xidx,:);
          end
        end
      end
    end
    
  case 'reduced'
    % save(sfname, 'run', 'data');
    run  = cache.cont.run;
    data = cache.cont.data;
    if isempty(fid) || isempty(data)
      soldata = data;
    else
      soldata = coco_slot_data(fid, data);
    end
    chart  = struct();
    chart1 = struct();
end

for i=1:nargout
  switch list{i}
    case 'data'
      varargout{i} = soldata;
    case 'run'
      varargout{i} = run;
    case 'chart'
      varargout{i} = chart;
    case 'chart1'
      varargout{i} = chart1;
    case {'uidx' 'xidx'}
      varargout{i} = xidx;
  end
end

end

function [soldata chart run chart1] = read_sol_v0(cache, fid)
% save(fname, 'data', 'sol', 'run', ['pt0'])
vars = cache.vars;
if any(strcmp('pt0', vars))
  % load(sfname, 'data', 'sol', 'run', 'pt0');
  run    = cache.cont.run;
  data   = cache.cont.data;
  chart  = cache.cont.sol;
  chart1 = cache.cont.pt0;
else
  % load(sfname, 'data', 'sol', 'run');
  run    = cache.cont.run;
  data   = cache.cont.data;
  chart  = cache.cont.sol;
  chart1 = struct();
end
if isempty(fid) || isempty(data)
  soldata = data;
else
  soldata = coco_slot_data(fid, data);
end
end
