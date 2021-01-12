function [opts, chart, J] = efunc_DFDU(opts, chart, x)
%EFUNC_DFDU  Evaluate Jacobian of extended system
%
%   [OPTS CHART J] = EFUNC_DFDU(OPTS, CHART, X) 
%
%   Evaluates the Jacobian of the extended system with respect to
%   the continuation variables.
%
% See also: EFUNC_CALL_DF, EFUNC_DFDX
%

% create cache for variables that grow
persistent urows ucols uvals

if nargin==0
  % clear cache
  urows = [];
  ucols = [];
  uvals = [];
  return;
end

%% extract u from xp=[u;p]

efunc = opts.efunc;
nrows = efunc.f_dim;
ncols = efunc.xp_dim;
idx_MX = 0;

%% compute Jacobian of extended continuation problem

for i = [ efunc.zero efunc.embedded ]
  func  = efunc.funcs(i);
  f_idx = func.f_idx;
  x_idx = func.x_idx;
  data  = func.data;
  
  % evaluate Jacobian of function
  [opts, data, chart, DFDX] = efunc_call_DF(opts, data, chart, func, x(x_idx));
  
  % merge Jacobians
  [r, c, v]   = find(DFDX);
  idx       = idx_MX+(1:numel(r));
  idx_MX    = idx_MX+numel(r);
  urows(idx) = f_idx(r)';
  ucols(idx) = x_idx(c)';
  uvals(idx) = v;
  
  % allow modification of associated toolbox data
  % for example, counting evaluations of F and DFDX
  opts.efunc.funcs(i).data = data;
  
end

%% apply permutations

idx = 1:idx_MX;
J = sparse(efunc.f_idx(urows(idx)), efunc.xp_idx(ucols(idx)), uvals(idx), nrows, ncols);
