function [prob, data] = ep_var_construct_eqn(prob, data, vecs)
%EP_VAR_CONSTRUCT_EQN   Add EP variational equation.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ep_var_construct_eqn.m 2897 2015-10-07 17:43:39Z hdankowicz $

[data, sol] = init_data(data, vecs);
fid  = data.ep_var.fid;
tbid = data.ep_eqn.fid;
uidx = coco_get_func_data(prob, tbid, 'uidx');
prob = coco_add_func(prob, fid, @FDF, data, 'zero', 'uidx', uidx, ...
  'u0', sol.var.u0, 't0', sol.var.t0, 'F+DF', 'requires', tbid);

end

function [data, sol] = init_data(data, vecs)
%INIT_DATA   Initialize data for variational equation.

eqn  = data.ep_eqn;
tbid = eqn.fid;
fid  = coco_get_id(tbid, 'var');

assert(isfield(data, 'ep_Jx'), ...
  '%s: %s: cannot add variational equation, ''%s'' not constructed with option ''-cache-jac''', ...
  mfilename, fid, tbid);

xdim = data.xdim;
nvec = size(vecs,2);
var.nvec  = nvec;
var.fid   = fid;
idx  = reshape(1:xdim*nvec, [xdim nvec]);
var.x_idx = eqn.x_idx;
var.p_idx = eqn.p_idx;
var.v_idx = var.p_idx(end)+idx;
var.w_idx = var.v_idx(end)+idx;

var.rows  = repmat(idx, [xdim 1]);
var.cols  = repmat(1:nvec*xdim, [xdim 1]);
var.J     = sparse(xdim*nvec, var.p_idx(end));
var.Jw    = -speye(xdim*nvec);

v = [vecs data.ep_Jx*vecs];
sol.var.u0 = v(:);
sol.var.t0 = [];

data.ep_var  = var;
data.no_save = union(data.no_save, ...
  {'ep_var.Jw' 'ep_var.J' 'ep_var.rows', 'ep_var.cols'});

end

function [data, y, J] = FDF(~, data, u)
%FDF   Variational equation and linearization.

pr  = data.pr;
var = pr.ep_var;

x = u(var.x_idx);
p = u(var.p_idx);
v = u(var.v_idx);
w = u(var.w_idx);

Jx = data.ep_Jx;
y  = Jx*v - w;
y  = y(:);

if nargout>=3
  xdim = pr.xdim;
  J = var.J;
  h = 1.0e-4*(1+norm(x));
  dfdx = pr.ode_DFDX;
  dfdp = pr.ode_DFDP;
  for i=1:var.nvec
    J0 = dfdx(pr, [], x-h*v(:,i), p);
    J1 = dfdx(pr, [], x+h*v(:,i), p);
    Jxx = (0.5/h)*(J1-J0);
    
    J0 = dfdp(pr, [], x-h*v(:,i), p);
    J1 = dfdp(pr, [], x+h*v(:,i), p);
    Jpx = (0.5/h)*(J1-J0);
    J((i-1)*xdim+(1:xdim),:) = [Jxx Jpx];
  end
  Jx = sparse(var.rows, var.cols, repmat(Jx, [1 var.nvec]));
  J = [J Jx var.Jw];
end

end
