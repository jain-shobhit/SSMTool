function [prob, data] = ep_SN_construct_eqn(prob, data, sol)
%EP_SN_CONSTRUCT_EQN   Add Moore-Spence system for SN points.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ep_SN_construct_eqn.m 2872 2015-08-06 20:12:06Z hdankowicz $

vfid = data.ep_var.fid;

data = init_data(data);
fid  = data.ep_sn.fid;
uidx = coco_get_func_data(prob, vfid, 'uidx');
uidx = uidx([data.ep_var.v_idx ; data.ep_var.w_idx]);
prob = coco_add_func(prob, fid, @FDF, data, 'zero', ...
  'uidx', uidx, 'u0', sol.sn.u0, 't0', sol.sn.t0, 'F+dF', ...
  'requires', vfid);

end

function data = init_data(data)
%INIT_DATA   Initialize data for Moore-Spence system.

xdim     = data.xdim;
sn.v_idx = (1:xdim)';
sn.w_idx = xdim+sn.v_idx;
sn.m     = xdim+1;
sn.n     = 2*xdim;
sn.rows  = [sn.m*ones(1,xdim) 1:xdim];
sn.cols  = 1:sn.n;
sn.vals  = ones(1,xdim);
sn.fid   = coco_get_id(data.ep_eqn.fid, 'SN');

data.ep_sn   = sn;
data.no_save = [ data.no_save { 'ep_sn.rows' 'ep_sn.cols' 'ep_sn.vals' } ];

end

function [data, y, J] = FDF(~, data, u)
%FDF   Moore-Spence system and linearization

sn  = data.ep_sn;

v = u(sn.v_idx);
w = u(sn.w_idx);

y  = [ w ; v'*v-1 ];

if nargout>=3
  J = sparse(sn.rows, sn.cols, [2*v' sn.vals], sn.m, sn.n);
end

end
