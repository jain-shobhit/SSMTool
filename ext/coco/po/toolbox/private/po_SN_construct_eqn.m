function [prob, data] = po_SN_construct_eqn(prob, data, sol)
%PO_SN_CONSTRUCT_EQN   Add saddle-node zero problem.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_SN_construct_eqn.m 2904 2015-10-10 02:25:50Z hdankowicz $

vfid = coco_get_id(data.cid, 'var');
[fdata, uidx] = coco_get_func_data(prob, vfid, 'data', 'uidx');
var  = fdata.coll_var;
if data.ode.autonomous  
  maps = fdata.coll_seg.maps;
  uidx = uidx([maps.x0_idx; maps.p_idx; var.v0_idx; var.v1_idx]);
  data = init_data1(prob, data);
  fid  = data.po_sn.fid;
  prob = coco_add_func(prob, fid, @FDF1, data, 'zero', ...
    'uidx', uidx, 'u0', sol.sn.u0, 't0', sol.sn.t0, ...
    'requires', vfid, 'F+DF');
else
  uidx = uidx([var.v0_idx; var.v1_idx]);
  data = init_data2(prob, data);
  fid  = data.po_sn.fid;
  prob = coco_add_func(prob, fid, @FDF2, data, 'zero', ...
    'uidx', uidx, 'u0', sol.sn.u0, 't0', sol.sn.t0, ...
    'requires', vfid, 'F+DF');
end

end

function data = init_data1(prob, data)
%INIT_DATA   Initialize data for autonomous saddle-node bifurcation zero problem.

fdata = coco_get_func_data(prob, data.cid, 'data', 'u0');
seg  = fdata.coll_seg;
dim  = seg.int.dim;
pdim = seg.maps.pdim;

sn.ode_F    = fdata.ode_F;
sn.ode_DFDX = fdata.ode_DFDX;
sn.ode_DFDP = fdata.ode_DFDP;
sn.fhan     = fdata.fhan;
sn.dfdxhan  = fdata.dfdxhan;
sn.dfdphan  = fdata.dfdphan;
sn.ode      = fdata.ode;

sn.x_idx = 1:dim;
sn.p_idx = dim+(1:pdim);
sn.v_idx = dim+pdim+(1:dim);
sn.w_idx = 2*dim+pdim+(1:dim);
sn.b_idx = 3*dim+pdim+1;

sn.m     = dim+2;
sn.n     = 3*dim+pdim+1;
sn.rows  = [repmat(1:dim, [1 dim+pdim+3]) (dim+1)*ones(1,2*dim+pdim) ...
  (dim+2)*ones(1,dim)];
sn.cols  = [reshape(repmat(sn.x_idx, [dim 1]), [1 dim^2])...
  reshape(repmat(sn.p_idx, [dim 1]), [1 dim*pdim]) ...
  sn.v_idx sn.w_idx sn.b_idx*ones(1,dim) ...
  sn.x_idx sn.p_idx sn.v_idx sn.v_idx];
sn.vals  = [ones(dim,1); -ones(dim,1)];

sn.fid   = coco_get_id(data.po_orb.fid, 'SN');

data.po_sn   = sn;
data.no_save = [ data.no_save ...
  { 'po_sn.ode' 'po_sn.ode_F' 'po_sn.ode_DFDX' 'po_sn.ode_DFDP' ...
  'po_sn.fhan' 'po_sn.dfdxhan' 'po_sn.dfdphan'...
  'po_sn.m' 'po_sn.n' 'po_sn.rows' 'po_sn.cols' 'po_sn.vals' }];

end

function data = init_data2(prob, data)
%INIT_DATA   Initialize data for non-autonomous saddle-node bifurcation zero problem.

fdata = coco_get_func_data(prob, data.cid, 'data');
seg  = fdata.coll_seg;
dim  = seg.int.dim;

sn.v_idx = 1:dim;
sn.w_idx = dim+(1:dim);

sn.m     = dim+1;
sn.n     = 2*dim;
sn.rows  = [repmat(1:dim, [1 2]) (dim+1)*ones(1,dim)];
sn.cols  = [sn.v_idx sn.w_idx sn.v_idx];
sn.vals  = [ones(dim,1); -ones(dim,1)];

sn.fid   = coco_get_id(data.po_orb.fid, 'SN');

data.po_sn   = sn;
data.no_save = [ data.no_save ...
  { 'po_sn.m' 'po_sn.n' 'po_sn.rows' 'po_sn.cols' 'po_sn.vals' } ];

end

function [data, y, J] = FDF1(~, data, u)
%FDF   Saddle-node zero problem and linearization for autonomous vector field

sn  = data.po_sn;

x = u(sn.x_idx);
p = u(sn.p_idx);
v = u(sn.v_idx);
w = u(sn.w_idx);
b = u(sn.b_idx);

f = sn.ode_F(sn, 0, x, p);
y = [ v+b*f-w ; f'*v; v'*v-1 ];
if nargout>=3
  dfdx = sn.ode_DFDX(sn, 0, x, p);
  dfdp = sn.ode_DFDP(sn, 0, x, p);
  J = sparse(sn.rows, sn.cols, [b*dfdx(:); b*dfdp(:); sn.vals; f; ...
    dfdx'*v; dfdp'*v; f; 2*v], sn.m, sn.n);
end

end

function [data, y, J] = FDF2(~, data, u)
%FDF   Saddle-node zero problem and linearization for non-autonomous vector field

sn  = data.po_sn;

v = u(sn.v_idx);
w = u(sn.w_idx);

y  = [ v-w ; v'*v-1 ];
if nargout>=3
  J    = sparse(sn.rows, sn.cols, [sn.vals; 2*v], sn.m, sn.n);
end

end
