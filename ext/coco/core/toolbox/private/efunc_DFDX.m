function [opts, chart, J] = efunc_DFDX(opts, chart, xp)
%EFUNC_DFDX  Evaluate Jacobian of extended system at UP=[U;P].
%
%   [OPTS J] = EFUNC_DFDX(OPTS, CHART, X) 
%
%   Evaluates the Jacobian of the extended continuation problem with
%   respect to the continuation variables and all active continuation
%   parameters.
%
% See also: EFUNC_DFDU
%

% create cache for variables that grow
persistent rows cols vals

if nargin==0
  % clear cache
  rows = [];
  cols = [];
  vals = [];
  return;
end

%% extract x from xp=[x;p] and
%  initialise derivatives with respect to x(efunc.p_idx,1)

efunc          = opts.efunc;
x              = xp(efunc.xp_idx,1);
J              = efunc.J_init;
[r, c, v]      = find(J);
[nrows, ncols] = size(J);
idx_MX         = numel(r);
idx            = 1:idx_MX;
rows(idx)      = r;
cols(idx)      = c;
vals(idx)      = v;

%% compute Jacobian of extended continuation problem

[opts, chart, J] = efunc_DFDU(opts, chart, x);
[r, c, v] = find(J);
idx       = idx_MX+(1:numel(r));
idx_MX    = idx_MX+numel(r);
rows(idx) = r;
cols(idx) = c;
vals(idx) = v;

%% apply permutations
% [opts chart J1] = coco_ezDFDX('f(o,c,x)', opts, @efunc_F, chart, xp);

idx = 1:idx_MX;
J = sparse(efunc.f_idx(rows(idx)), efunc.xp_idx(cols(idx)), vals(idx), nrows, ncols);
