function prob = po_construct_adjt(prob, data, sol)
%PO_CONSTRUCT_ADJT   Add PO adjoint problem.

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: po_construct_adjt.m 2872 2015-08-06 20:12:06Z hdankowicz $

data = init_data(data);
opt  = data.po_opt;
fid  = opt.fid;

[fdata, axidx] = coco_get_adjt_data(prob, data.cid, 'data', 'axidx');
fopt = fdata.coll_opt;

if data.ode.autonomous
  aidx   = axidx([fopt.xcn_idx; fopt.x0_idx; fopt.x1_idx]);
  prob   = coco_add_adjt(prob, fid, @adj, @adj_DU, data, 'aidx', aidx, ....
    'l0', sol.l0, 'tl0', sol.tl0, 'remesh', @adj_remesh);
  up_fid = coco_get_id(fid, 'up');
  prob   = coco_add_slot(prob, up_fid, @update, data, 'update');
else
  aidx = axidx([fopt.x0_idx; fopt.x1_idx]);
  prob = coco_add_adjt(prob, fid, @adj, @adj_DU, data, 'aidx', aidx, ....
  'l0', sol.l0, 'tl0', sol.tl0);
  sfid = coco_get_id(fid, 'tinit');
  dpar = coco_get_id('d', sfid);
  prob = coco_add_adjt(prob, sfid, dpar, 'aidx',  axidx(fopt.T0_idx), ...
    'l0', sol.tinit_l0, 'tl0', sol.tinit_tl0);
end
pfid   = coco_get_id(fid, 'period');
dpar   = coco_get_id('d', pfid);
prob   = coco_add_adjt(prob, pfid, dpar, 'aidx', axidx(fopt.T_idx), ...
  'l0', sol.T_l0, 'tl0', sol.T_tl0);

end

function data = init_data(data)
%INIT_DATA   Initialize PO data with phase condition

opt.fid = coco_get_id(data.oid, 'po');
data.po_opt = opt;

end

function [data, J] = adj(prob, data, u) %#ok<INUSD,INUSL>
%FDF   PO zero problem with Jacobian

J = data.po_orb.J_adj;

end

function [data, dJ] = adj_DU(prob, data, u) %#ok<INUSD,INUSL>

[m,n] = size(data.po_orb.J_adj);
dJ = zeros(m,n,n);

end

function data = update(prob, data, cseg, varargin) %#ok<INUSD>

fdata = coco_get_func_data(prob, data.po_opt.fid, 'data');
data.po_orb = fdata.po_orb;

end

function [prob, status, xtr, ftr] = adj_remesh(prob, data, chart, lb, Vlb)  %#ok<INUSL>

fdata = coco_get_func_data(prob, data.po_opt.fid, 'data');
data.po_orb = fdata.po_orb;
[fdata, axidx] = coco_get_adjt_data(prob, fdata.cid, 'data', 'axidx');
fopt = fdata.coll_opt;

xtr    = [];
aidx = axidx([fopt.xcn_idx; fopt.x0_idx; fopt.x1_idx]);
adim = [numel(fopt.x0_idx)+1, numel(aidx)];
ftr  = (1:numel(fopt.x0_idx)+1)';

prob   = coco_change_adjt(prob, data, 'l0', lb, 'aidx', aidx, ...
  'adim', adim, 'vecs', Vlb);
status = 'success';

end
