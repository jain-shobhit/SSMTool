function [prob, data] = po_TR_construct_eqn(prob, data, sol)
%PO_TR_CONSTRUCT_EQN   Add torus bifurcation zero problem.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_TR_construct_eqn.m 2904 2015-10-10 02:25:50Z hdankowicz $

vfid = coco_get_id(data.cid, 'var');
[fdata, uidx] = coco_get_func_data(prob, vfid, 'data', 'uidx');
var  = fdata.coll_var;
uidx = uidx([var.v0_idx; var.v1_idx]);

data = init_data(prob, data);
fid  = data.po_tr.fid;
prob = coco_add_func(prob, fid, @FDF, data, 'zero', 'uidx', uidx, ...
  'u0', sol.tr.u0, 't0', sol.tr.t0, 'requires', vfid, 'F+DF');

end

function data = init_data(prob, data)
%INIT_DATA   Initialize data for torus bifurcation zero problem.

fdata = coco_get_func_data(prob, data.cid, 'data');
seg  = fdata.coll_seg;
dim  = seg.int.dim;

tr.v1_idx = 1:dim;
tr.w1_idx = dim+(1:dim);
tr.v2_idx = 2*dim+(1:dim);
tr.w2_idx = 3*dim+(1:dim);
tr.a_idx  = 4*dim+1;
tr.b_idx  = 4*dim+2;

tr.m      = 2*dim+3;
tr.n      = 4*dim+2;
tr.rows   = [repmat(1:dim, [1 5]), repmat(dim+(1:dim), [1 5]), ...
  (2*dim+1)*ones(1,2*dim), (2*dim+2)*ones(1,2*dim), 2*dim+3, 2*dim+3];
tr.cols   = [tr.v1_idx, tr.v2_idx, tr.w1_idx, ...
  tr.a_idx*ones(1,dim), tr.b_idx*ones(1,dim), ...
  tr.v1_idx, tr.v2_idx, tr.w2_idx, ...
  tr.a_idx*ones(1,dim), tr.b_idx*ones(1,dim), ...
  tr.v1_idx, tr.v2_idx, tr.v1_idx, tr.v2_idx, tr.a_idx, tr.b_idx];
tr.vals   = ones(dim,1);

tr.fid    = coco_get_id(data.po_orb.fid, 'TR');

data.po_tr   = tr;
data.no_save = [ data.no_save ...
  { 'po_tr.m' 'po_tr.n' 'po_tr.rows' 'po_tr.cols' 'po_tr.vals' } ];

end

function [data, y, J] = FDF(~, data, u)
%FDF   Torus bifurcation system and linearization

tr  = data.po_tr;
dim = numel(tr.v1_idx);

v1 = u(tr.v1_idx);
v2 = u(tr.v2_idx);
w1 = u(tr.w1_idx);
w2 = u(tr.w2_idx);
a  = u(tr.a_idx);
b  = u(tr.b_idx);

y  = [ w1-a*v1+b*v2; w2-a*v2-b*v1 ; v1'*v1+v2'*v2-1; v1'*v2 ; a^2+b^2-1];
if nargout>=3
  vals = [-a*ones(dim,1); b*ones(dim,1); tr.vals; -v1; v2; ...
    -b*ones(dim,1); -a*ones(dim,1); tr.vals; -v2; -v1; ...
    2*v1; 2*v2; v2; v1; 2*a; 2*b];
  J    = sparse(tr.rows, tr.cols, vals, tr.m, tr.n);
end

end
