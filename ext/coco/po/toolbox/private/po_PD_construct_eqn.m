function [prob, data] = po_PD_construct_eqn(prob, data, sol)
%PO_PD_CONSTRUCT_EQN   Add period-doubling zero problem.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_PD_construct_eqn.m 2904 2015-10-10 02:25:50Z hdankowicz $

data = init_data(prob, data);

vfid = coco_get_id(data.cid, 'var');
[fdata, uidx_var] = coco_get_func_data(prob, vfid, 'data', 'uidx');
var  = fdata.coll_var;

fid  = data.po_pd.fid;
uidx = uidx_var([var.v0_idx; var.v1_idx]);
prob = coco_add_func(prob, fid, @FDF, data, 'zero', 'uidx', uidx, ...
  'u0', sol.pd.u0, 't0', sol.pd.t0, 'requires', vfid, 'F+DF');

prob = coco_add_chart_data(prob, fid, [], []);

cfid = coco_get_id(fid, 'dbl');
prob = coco_add_func(prob, cfid, @double, data, 'regular', {}, ...
  'uidx', uidx_var, 'remesh', @remesh, 'passChart');

end

function data = init_data(prob, data)
%INIT_DATA   Initialize data for period-doubling zero problem.

fdata = coco_get_func_data(prob, data.cid, 'data');
seg  = fdata.coll_seg;
maps = seg.maps;
mesh = seg.mesh;
dim  = seg.int.dim;

pd.v_idx = 1:dim;
pd.w_idx = dim+(1:dim);

pd.xbp_idx = maps.xbp_idx;
pd.T_idx   = maps.T_idx;
pd.p_idx   = maps.p_idx;
pd.var_idx = numel(maps.xbp_idx)+1+numel(maps.p_idx)+maps.xbp_idx;

pd.maps  = maps;
pd.mesh  = mesh;

pd.m     = dim+1;
pd.n     = 2*dim;
pd.rows  = [repmat(1:dim, [1 2]) (dim+1)*ones(1,dim)];
pd.cols  = [pd.v_idx pd.w_idx pd.v_idx];
pd.vals  = [ones(dim,1); ones(dim,1)];

pd.fid   = coco_get_id(data.po_orb.fid, 'PD');

data.po_pd   = pd;
data.no_save = [ data.no_save ...
  { 'po_pd.maps' 'po_pd.mesh' 'po_pd.xbp_idx' ...
  'po_pd.T_idx' 'po_pd.p_idx' 'po_pd.var_idx' ...
  'po_pd.m' 'po_pd.n' 'po_pd.rows' 'po_pd.cols' 'po_pd.vals' } ];

end

function [data, y, J] = FDF(~, data, u)
%FDF   Period-doubling system and linearization

pd = data.po_pd;

v = u(pd.v_idx);
w = u(pd.w_idx);

y = [ v+w ; v'*v-1 ];
if nargout>=3
  J = sparse(pd.rows, pd.cols, [pd.vals; 2*v], pd.m, pd.n);
end

end

function [data, chart, y] = double(prob, data, chart, u) %#ok<INUSL>

pd = data.po_pd;
maps = pd.maps;
mesh = pd.mesh;

x = u(pd.xbp_idx);
T = u(pd.T_idx);
p = u(pd.p_idx);
v = u(pd.var_idx);

% Construct initial solution guess for period-doubled orbit
xp1   = reshape(x+0.01*v, maps.xbp_shp)'; 
xp1   = xp1(maps.tbp_idx,:);
t1    = mesh.tbp(maps.tbp_idx)*T;
xp2   = reshape(x-0.01*v, maps.xbp_shp)';
xp2   = xp2(maps.tbp_idx,:);
t2    = mesh.tbp(maps.tbp_idx)*T;
x0    = [xp1; xp2(2:end,:)];
t0    = [t1; T+t2(2:end)];
cdata = coco_get_chart_data(chart, pd.fid);
cdata.pd.x0 = x0;
cdata.pd.t0 = t0;
cdata.pd.p0 = p;
chart = coco_set_chart_data(chart, pd.fid, cdata);

y = [];

end

function [prob, status, xtr] = remesh(prob, data, chart, old_u, old_V) %#ok<INUSD>

pd = data.po_pd;

[fdata, uidx] = coco_get_func_data(prob, data.cid, 'data', 'uidx');
seg  = fdata.coll_seg;
maps = seg.maps;
mesh = seg.mesh;

pd.xbp_idx = maps.xbp_idx;
pd.T_idx   = maps.T_idx;
pd.p_idx   = maps.p_idx;
pd.var_idx = numel(maps.xbp_idx)+maps.xbp_idx;

pd.maps = maps;
pd.mesh = mesh;

vfid = coco_get_id(data.cid, 'var');
[fdata, uidx_var] = coco_get_func_data(prob, vfid, 'data', 'uidx');
var  = fdata.coll_var;

uidx   = [uidx; uidx_var(var.v_idx);];
xtr    = [];
prob   = coco_change_func(prob, data, 'uidx', uidx);
status = 'success';

end
