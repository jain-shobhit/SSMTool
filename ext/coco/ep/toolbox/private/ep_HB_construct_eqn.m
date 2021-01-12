function [prob, data] = ep_HB_construct_eqn(prob, data, sol)
%EP_HB_CONSTRUCT_EQN   Add Hopf system.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ep_HB_construct_eqn.m 2897 2015-10-07 17:43:39Z hdankowicz $

vfid = data.ep_var.fid;

data = init_data(data, sol);
fid  = data.ep_hb.fid;
uidx = coco_get_func_data(prob, vfid, 'uidx');
uidx = uidx([data.ep_var.v_idx(:,1) data.ep_var.w_idx(:,2)]);

prob = coco_add_func(prob, fid, @FDF, data, 'zero', ...
  'uidx', uidx, 'u0', sol.hb.u0, 't0', sol.hb.t0, 'F+dF', ...
  'requires', vfid);
prob = coco_add_slot(prob, fid, @update, data, 'update');

end

function data = init_data(data, sol)
%INIT_DATA   Initialize data for EP Hopf system.

xdim     = data.xdim;
hb.nv    = sol.hb.nv';
hb.v_idx = 1:xdim;
hb.w_idx = xdim+hb.v_idx;
hb.k_idx = 2*xdim+1;
hb.m     = xdim+2;
hb.n     = 2*xdim+1;
o        = ones(1,xdim);
hb.rows  = [1:xdim 1:xdim 1:xdim (xdim+1)*o (xdim+2)*o];
hb.cols  = [1:2*xdim (2*xdim+1)*o 1:xdim 1:xdim];
hb.vals  = o;
hb.fid   = coco_get_id(data.ep_eqn.fid, 'HB');

data.ep_hb   = hb;
data.no_save = [ data.no_save ...
  { 'ep_hb.rows' 'ep_hb.cols' 'ep_hb.vals' 'ep_hb.m' 'ep_hb.n' } ];

end

function [data, y, J] = FDF(~, data, u)
%FDF   Hopf zero problem.

hb  = data.ep_hb;

v = u(hb.v_idx);
w = u(hb.w_idx);
k = u(hb.k_idx);

y = [k*v+w ; v'*v-1 ; hb.nv*v ];

if nargout>=3
  vals = [k*hb.vals hb.vals v' 2*v' hb.nv];
  J = sparse(hb.rows, hb.cols, vals, hb.m, hb.n);
end

end

function data = update(prob, data, cseg, varargin)
%UPDATE   Update orthogonality vector data.hb.nv.

pr     = data.pr;
hb     = pr.ep_hb;
ep_var = pr.ep_var;

chart = cseg.src_chart;
uidx  = coco_get_func_data(prob, ep_var.fid, 'uidx');
u     = chart.x(uidx);
v     = u(ep_var.v_idx(:,1));
w     = u(ep_var.w_idx(:,1));
al    = [-w'*v v'*v];
al    = al/norm(al);
nv    = al(1)*v + al(2)*w;
nv    = nv/norm(nv);
hb.nv = nv';

pr.ep_hb = hb;
data.pr = pr;

end
