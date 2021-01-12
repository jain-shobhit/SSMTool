function [prob, data] = po_close_orb(prob, data, PhaseCondition)

if PhaseCondition % autonomous vector field
  data = init_data1(prob, data);
  orb  = data.po_orb;
  fid  = orb.fid;
  
  [fdata, uidx] = coco_get_func_data(prob, data.cid, 'data', 'uidx');
  seg  = fdata.coll_seg;
  maps = seg.maps;
  
  prob = coco_add_func(prob, fid, @FDF, data, 'zero', ...
    'uidx', uidx(maps.xbp_idx), 'fdim', seg.int.dim + 1, ...
    'remesh', @remesh, 'F+DF');
  prob = coco_add_slot(prob, fid, @update, data, 'update');  
  pfid = coco_get_id(fid, 'period');
  prob = coco_add_pars(prob, pfid, uidx(maps.T_idx), pfid, 'active');
else % non-autonomous vector field
  data = init_data2(prob, data);
  orb  = data.po_orb;
  fid  = orb.fid;
  
  [fdata, uidx] = coco_get_func_data(prob, data.cid, 'data', 'uidx');
  seg  = fdata.coll_seg;
  maps = seg.maps;
  
  prob = coco_add_func(prob, fid, @FDF, data, 'zero', ...
    'uidx', uidx([maps.x0_idx; maps.x1_idx]), 'fdim', seg.int.dim, ...
    'F+DF');
  sfid = coco_get_id(fid, 'tinit');
  prob = coco_add_pars(prob, sfid, uidx(maps.T0_idx), sfid);
  pfid = coco_get_id(fid, 'period');
  prob = coco_add_pars(prob, pfid, uidx(maps.T_idx), pfid);
end

end

function data = init_data1(prob, data)
%INIT_DATA1   Initialize PO data with phase condition

[fdata, u0] = coco_get_func_data(prob, data.cid, 'data', 'u0');

seg  = fdata.coll_seg;
dim  = seg.int.dim;
maps = seg.maps;
mesh = seg.mesh;

po.dim    = dim-1;
rows      = [1:dim 1:dim];
cols      = [maps.x0_idx maps.x1_idx];
vals      = [ones(1,dim) -ones(1,dim)];
J         = sparse(rows, cols, vals, dim, numel(maps.xbp_idx)); % periodic bound.
x0        = u0(maps.xbp_idx);
po.intfac = maps.Wp'*mesh.wts2*maps.W; % integral phase cond.
po.J      = [J; x0'*po.intfac];

po.J_pbc  = [eye(dim) -eye(dim)];
po.J_phfc = maps.Wp';
po.J_pha  = x0'*po.J_phfc;
po.J_adj  = [[zeros(dim, size(po.J_pha,2)) po.J_pbc]; [po.J_pha zeros(1,2*dim)]];

po.fid    = coco_get_id(data.oid, 'po');

data.po_orb = po;
data.no_save = [ data.no_save { 'po_orb.J' 'po_orb.intfac' } ];

end

function data = init_data2(prob, data)
%INIT_DATA2   Initialize PO data with no phase condition

fdata = coco_get_func_data(prob, data.cid, 'data');

seg  = fdata.coll_seg;
dim  = seg.int.dim;

rows = [1:dim 1:dim];
cols = 1:2*dim;
vals = [ones(1,dim) -ones(1,dim)];
po.J = sparse(rows, cols, vals, dim, 2*dim);
po.J_adj = po.J;

po.dim = dim;
po.fid = coco_get_id(data.oid, 'po');

data.po_orb = po;
data.no_save = [ data.no_save { 'po_orb.J' } ];

end

function [data, y, J] = FDF(prob, data, u) %#ok<INUSL>
%FDF   PO zero problem with Jacobian

J = data.po_orb.J;
y = J*u;

end

function data = update(prob, data, cseg, varargin)

[fdata, uidx] = coco_get_func_data(prob, data.cid, 'data', 'uidx');
seg  = fdata.coll_seg;
dim  = seg.int.dim;
maps = seg.maps;
po   = data.po_orb;

u           = cseg.src_chart.x(uidx);
x0          = u(maps.xbp_idx);
po.J(end,:) = x0'*po.intfac;
po.J_pha    = x0'*po.J_phfc;
po.J_adj    = [[zeros(dim, size(po.J_pha,2)) po.J_pbc]; [po.J_pha zeros(1,2*dim)]];

data.po_orb = po;

end

function [prob, status, xtr] = remesh(prob, data, chart, old_u, old_V) %#ok<INUSD>

[fdata, uidx] = coco_get_func_data(prob, data.cid, 'data', 'uidx');
seg  = fdata.coll_seg;
dim  = seg.int.dim;
maps = seg.maps;
mesh = seg.mesh;
po   = data.po_orb;

rows      = [1:dim 1:dim];
cols      = [maps.x0_idx maps.x1_idx];
vals      = [ones(1,dim) -ones(1,dim)];
J         = sparse(rows, cols, vals, dim, numel(maps.xbp_idx));
po.J      = [J; zeros(1, numel(maps.xbp_idx))];
po.intfac = maps.Wp'*mesh.wts2*maps.W;
po.J_phfc = maps.Wp';

data.po_orb = po;

xtr       = [];
prob      = coco_change_func(prob, data, 'uidx', uidx(maps.xbp_idx), ...
  'fdim', dim + 1);
status    = 'success';

end
