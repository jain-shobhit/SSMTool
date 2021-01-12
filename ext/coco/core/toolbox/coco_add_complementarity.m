function prob = coco_add_complementarity(prob, fid, varargin)
%COCO_ADD_COMPLEMENTARITY   Add nonlinear complimentarity condition as complementary monitor function
%
% PROB = COCO_ADD_COMPLEMENTARITY(PROB, FID, VARARGIN)
%
% VARARGIN = { [ FHAN [ DFHAN ] V0 [ TV0 ] ] PAR_NAMES }

% Copyright (C) Harry Dankowicz, Mingwu Li
% $Id: coco_add_functionals.m 2969 2017-01-10 19:13:51Z hdankowicz $

%% parse input arguments

assert(isstruct(prob) && isfield(prob, 'opts') && isa(prob.opts, 'coco_opts_tree'), ...
  '%s: the first argument must be a continuation problem structure', mfilename);
assert(isfield(prob, 'efunc'), ...
  '%s: the continuation problem structure cannot be empty', mfilename);
efunc = prob.efunc;
assert(isfield(prob, 'adjoint'), ...
  '%s: the continuation problem structure must contain adjoint functions', mfilename);
adjoint = prob.adjoint;

eidx = find(strcmpi(fid, efunc.identifyers), 1);
aidx = find(strcmpi(fid, adjoint.identifyers), 1);
assert(~isempty(eidx) && ~isempty(aidx), ...
  '%s: the function ''%s'' has not been defined', mfilename, fid);

data.fid  = fid;

argidx = 1;

if isa(varargin{argidx}, 'function_handle')
  data.fhan = varargin{argidx};
  argidx = argidx + 1;
  if isa(varargin{argidx}, 'function_handle')
    data.dfhan = varargin{argidx};
    argidx = argidx + 1;
  else
    data.dfhan = [];
  end
  v0 = varargin{argidx};
  tv0 = [];
  argidx = argidx + 1;
  if isnumeric(varargin{argidx})
    tv0 = varargin{argidx};
    argidx = argidx + 1;
  end
else
  data.fhan  = @(a,b,~) sqrt(a.^2+b.^2)-a-b;
  data.dfhan = @(a,b,~) [ a./sqrt(a.^2+b.^2)-1 b./sqrt(a.^2+b.^2)-1 ];
  v0  = [];
  tv0 = [];
end

pnames = varargin{argidx};
if ~iscell(pnames)
  pnames = { pnames };
end
assert(iscellstr(pnames), ...
  '%s: argument %d must be a cell array of strings', mfilename, 2+argidx);
pnum   = numel(pnames);
pnames = reshape(pnames, 1, pnum);

[data, uidx, lidx] = init_data_ncp(prob, data, fid);
assert(data.l_dim==pnum, '%s: %s %s', mfilename, ...
  'number of parameter labels must equal', ...
  'the number of continuation multipliers');

cfid = coco_get_id('ncp', fid);
prob = coco_add_comp(prob, cfid, @func_ncp, data, 'inactive', pnames, ...
  'uidx', uidx, 'lidx', lidx, 'v0', v0, 'tv0', tv0, 'f+df', ...
  'remesh', @remesh_ncp);

end

function [data, uidx, lidx] = init_data_ncp(prob, data, fid)

efunc   = prob.efunc;
adjoint = prob.adjoint;

idx  = find(strcmpi(fid, efunc.identifyers), 1);
ncp  = efunc.funcs(idx);
uidx = ncp.x_idx;

idx = find(strcmpi(fid, adjoint.identifyers), 1);
lidx = adjoint.funcs(idx).af_idx;

data.u_dim = numel(uidx);
data.l_dim = numel(lidx);

data.ncp = ncp;

end

function [data, y, J] = func_ncp(prob, data, u, l, v)

ncp = data.ncp;
switch ncp.call_mode
  case 1 % [d f]=F(o,d,x); [d J]=DFDX(o,d,x)
    [ncp.data, g] = ncp.F(prob, ncp.data, u);
    if nargout==3
      [ncp.data, J] = ncp.DFDX(prob, ncp.data, u);
    end
  case 2 % [d f J]=FDF(o,d,x)
    [ncp.data, g, J] = ncp.F(prob, ncp.data, u);
end

v = repmat(v', [data.l_dim 1]);
y = data.fhan(l, -g, v);

if nargout==3
  jac = data.dfhan(l, -g, v);
  Ju  = repmat(-jac(:,2), [1 data.u_dim]).*J;
  Jl  = jac(:,1);
  Jv  = jac(:,3:end);
  J   = {Ju, diag(Jl), Jv};
end

data.ncp = ncp;

end

function [prob, stat, xtr] = remesh_ncp(prob, data, chart, ub, Vb) %#ok<INUSD,INUSL>

data = data.data;
[data, uidx, lidx] = init_data_ncp(prob, data, data.fid);
prob = coco_change_comp(prob, data, 'uidx', uidx, 'lidx' , lidx);
xtr  = 1:numel(ub);
stat = 'success';

end
