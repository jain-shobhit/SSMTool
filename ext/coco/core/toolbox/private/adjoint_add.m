function prob = adjoint_add(prob)
% ADJOINT_ADD   Add zero and default monitor functions corresponding to the generalized adjoint problem
%
% This function appends the generalized adjoint zero function 
%   (u,lambda,eta,sigma)->(lambda^T, eta^T, sigma^T)*C(u)
% and, when string labels are included in the call to coco_add_adjt, the
% corresponding monitor functions
%   (u,lambda,eta,sigma)->eta
% or
%   (u,lambda,eta,sigma)->sqrt(sigma.^2+G(u).^2)-sigma+G(u)
% with associated continuation parameters.

% Copyright (C) Harry Dankowicz, Mingwu Li
% $Id: adjoint_add.m 3085 2018-12-10 21:21:13Z hdankowicz $

data    = coco_func_data;
data    = init_data_adjoint(prob, data);
adjoint = data.adjoint;
prob    = coco_add_func(prob, 'coco_adjoint', @adeqs, data, 'zero', ...
  'uidx', data.x_idx, 'x0', adjoint.l0, 't0', adjoint.tl0, ...
  'fdim', adjoint.a_dim(2), 'f+df', 'remesh', @remesh_adjoint);
prob    = coco_add_slot(prob, 'coco_adjoint', @coco_save_data, data, ...
  'save_full');

if ~isempty(adjoint.lactive)
  prob = coco_add_comp_pars(prob, 'coco_adjoint_active_pars', ...
    adjoint.l_idx(adjoint.lactive), adjoint.lnames(adjoint.lactive), ...
    'active');
end
if ~isempty(adjoint.linactive)
  prob = coco_add_comp_pars(prob, 'coco_adjoint_inactive_pars', ...
    adjoint.l_idx(adjoint.linactive), adjoint.lnames(adjoint.linactive), ...
    'inactive');
end

for i = 1:numel(adjoint.identifyers)
  func = adjoint.funcs(i);
  if strcmpi(func.type, 'inequality') && ~isempty(func.snames)
    prob = coco_add_complementarity(prob, func.identifyer, func.snames);
  end
end

end

function data = init_data_adjoint(prob, data)

efunc      = prob.efunc;
adjoint    = prob.adjoint;
data.x_dim = efunc.x_dim;
data.x_idx = 1:data.x_dim;
data.l_dim = adjoint.a_dim(1);
data.l_idx = data.x_dim + (1:data.l_dim);

data.xtr     = zeros(data.l_dim,1);
data.xtrstay = data.l_idx([adjoint.l_idx; adjoint.s_idx]) - data.x_dim;
data.xtr(data.xtrstay) = data.xtrstay;

ncols    = adjoint.a_dim(2);
data.shp = [ncols, data.x_dim];
adjoint.dax_idx = 1:(ncols*data.x_dim);
for i=1:numel(adjoint.identifyers)
  func   = adjoint.funcs(i);
  ax_idx = func.ax_idx;
  x_idx  = func.x_idx;
  
  dax_idx = kron(ones(1,numel(x_idx)), ax_idx);
  dx_idx  = kron(x_idx, ones(1,numel(ax_idx)));
  dax_idx = dax_idx+ncols*(dx_idx-1);
  
  adjoint.funcs(i).dax_idx = dax_idx;
end

data.adjoint = adjoint;

end

function [data, f, J] = adeqs(prob, data, u)

pr = data.pr;

x = u(pr.x_idx);
l = u(pr.l_idx);

[pr, C] = adfunc(prob, pr, x);
f = l'*C;
f = f(:);

if nargout>=3
  [pr, dC] = adfunc_DU(prob, pr, x);
  J = [reshape(l'*dC, pr.shp),  C'];
end

data.pr = pr;

end

function [data, C] = adfunc(prob, data, x)

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

adjoint = data.adjoint;
nrows   = adjoint.a_dim(1);
ncols   = adjoint.a_dim(2);
idx_MX  = 0;

%% compute adjoint coefficient matrix

for i = 1:numel(adjoint.identifyers)
  func   = adjoint.funcs(i);
  af_idx = func.af_idx;
  ax_idx = func.ax_idx;
  x_idx  = func.x_idx;
  fdata  = func.data;
  
  % evaluate adjoint of function
  [fdata, A] = func.F(prob, fdata, x(x_idx));
  
  % merge adjoints
  [r, c, v]  = find(A);
  idx        = idx_MX+(1:numel(r));
  idx_MX     = idx_MX+numel(r);
  urows(idx) = af_idx(r)';
  ucols(idx) = ax_idx(c)';
  uvals(idx) = v;
  
  adjoint.funcs(i).data = fdata;
end

data.adjoint = adjoint;
%% apply permutations

idx = 1:idx_MX;
C = sparse(adjoint.af_idx(urows(idx)), adjoint.ax_idx(ucols(idx)), uvals(idx), nrows, ncols);

end

function [data, dC] = adfunc_DU(prob, data, x)

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

adjoint = data.adjoint;
nrows   = adjoint.a_dim(1);
ncols   = adjoint.a_dim(2);
idx_MX  = 0;

%% compute Jacobian of adjoint matrix

for i = 1:numel(adjoint.identifyers)
  func    = adjoint.funcs(i);
  af_idx  = func.af_idx;
  dax_idx = func.dax_idx;
  x_idx   = func.x_idx;
  fdata   = func.data;
  
  % evaluate jacobian of adjoint of function
  if ~isempty(func.DFDX)
    [fdata, dA] = func.DFDX(prob, fdata, x(x_idx));
  else
    [fdata, dA] = coco_ezDFDX('f(o,d,x)', prob, fdata, @func.F, x(x_idx));
  end
  
  % reshape jacobian
  dA = reshape(dA, [numel(af_idx), numel(dA)/numel(af_idx)]);
    
  % merge adjoints
  [r, c, v]  = find(dA);
  idx        = idx_MX+(1:numel(r));
  idx_MX     = idx_MX+numel(r);
  urows(idx) = af_idx(r)';
  ucols(idx) = dax_idx(c)';
  uvals(idx) = v;
  
  adjoint.funcs(i).data = fdata;
end

data.adjoint = adjoint;

%% apply permutations

idx = 1:idx_MX;
dC = sparse(adjoint.af_idx(urows(idx)), adjoint.dax_idx(ucols(idx)), uvals(idx), nrows, ncols*numel(data.x_idx));

end

function [prob, stat, xtr] = remesh_adjoint(prob, data, chart, ub, Vb) %#ok<INUSD>
% Remesh function for adjoint problem.

xtr     = data.xtr;
xtrstay = data.xtrstay;
data    = init_data_adjoint(prob, data);
adjoint = data.adjoint;
xtr(xtrstay) = data.xtr(data.xtrstay);
prob    = coco_change_func(prob, data, 'uidx', data.x_idx, ...
  'u0', adjoint.l0, 'fdim', adjoint.a_dim(2), 'vecs', adjoint.Vl0);
stat = 'success';

end
