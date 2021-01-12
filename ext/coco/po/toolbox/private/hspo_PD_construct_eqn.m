function [prob, data] = hspo_PD_construct_eqn(prob, data, sol)
%HSPO_PD_CONSTRUCT_EQN   Add period-doubling zero problem.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_PD_construct_eqn.m 2849 2015-05-17 20:32:46Z hdankowicz $

data = init_data(prob, data);
pd   = data.hspo_pd;
fid  = pd.fid;

didx = [];
xidx = [];
vidx = [];
widx = [];
for i=1:pd.nsegs
  [fdata, uidx_var] = coco_get_func_data(prob, pd.vids{i}, 'data', 'uidx');
  didx = [didx; uidx_var]; %#ok<AGROW>
  maps = fdata.coll_seg.maps;
  var  = fdata.coll_var;
  xidx = [xidx; uidx_var(maps.x1_idx)]; %#ok<AGROW>
  vidx = [vidx; uidx_var(var.v0_idx)]; %#ok<AGROW>
  widx = [widx; uidx_var(var.v1_idx)]; %#ok<AGROW>
end
uidx = [uidx_var(maps.p_idx); xidx; vidx; widx];

prob = coco_add_func(prob, fid, @FDF, data, 'zero', 'uidx', uidx, ...
  'u0', sol.pd.u0, 't0', sol.pd.t0, 'requires', pd.vids, 'F+dF');

prob = coco_add_chart_data(prob, fid, [], []);

cfid = coco_get_id(fid, 'dbl');
prob = coco_add_func(prob, cfid, @double, data, 'regular', {}, ...
  'uidx', didx, 'remesh', @remesh, 'passChart');

end

function data = init_data(prob, data)
%INIT_DATA   Initialize data for period-doubling zero problem.

po = data.hspo_orb;
cdim = po.cdim;
pdim = po.pdim;
pd.p_idx = 1:pdim;
pd.x_idx = pdim + (1:cdim);
pd.v_idx = pdim + cdim + (1:cdim);
pd.w_idx = pdim + 2*cdim + (1:cdim);

fdata = coco_get_func_data(prob, data.bvid, 'data');
nsegs = fdata.nsegs;
vids = cell(1,nsegs);
off  = 0;
sidx = cell(1,nsegs);
for i=1:nsegs
  cdata   = coco_get_func_data(prob, fdata.cids{i}, 'data');
  vids{i} = cdata.coll_var.fid;
  sidx{i} = off + (1:po.dim(i));
  off     = off + po.dim(i);
end
pd.nsegs  = nsegs;
pd.cids   = fdata.cids;
pd.vids   = vids;
pd.sidx   = sidx;

contmat = eye(cdim);
frstdim = po.dim(1);
contmat(1:frstdim,1:frstdim)=-eye(frstdim);
pd.cont = circshift(contmat, -frstdim);
pd.zrs1 = zeros(1,pdim+cdim);
pd.zrs2 = zeros(1,2*cdim-frstdim);

segs = cell(1,nsegs);
off  = 0;
for i=1:nsegs
  cdata = coco_get_func_data(prob, vids{i}, 'data');
  mesh  = cdata.coll_seg.mesh;
  var   = cdata.coll_var;
  segs{i}.maps.tbp_idx = var.tbp_idx;
  segs{i}.maps.xbp_idx = off + var.xbp_idx;
  segs{i}.maps.xbp_shp = var.xbp_shp;
  segs{i}.maps.T_idx = off + var.T_idx;
  segs{i}.maps.p_idx = off + var.p_idx;
  segs{i}.maps.var_idx = off + var.v_idx;
  segs{i}.mesh.tbp = mesh.tbp;
  off   = off + var.v_idx(end);
end
pd.segs = segs;

pd.fid   = coco_get_id(data.oid, 'hspo.PD');

data.hspo_pd   = pd;
data.no_save = [ data.no_save ...
  { 'hspo_pd.zrs1' 'hspo_pd.zrs2' 'hspo_pd.segs' } ];

end

function [data, y, J] = FDF(prob, data, u)
%FDF   Period-doubling system and linearization

function [data, f] = temp(prob, data, u, w)
  f = hspo_proj(prob, data, u)*w;
end

pr = data.pr;
pd = pr.hspo_pd;

p = u(pd.p_idx);
x = u(pd.x_idx);
v = u(pd.v_idx);
w = u(pd.w_idx);

Px = hspo_proj(prob, pr, [p; x]);

cont = Px*w - pd.cont*v;
unit = v(pd.sidx{1})'*v(pd.sidx{1}) - 1;

y = [cont(:); unit];

if nargout>=3
  [data, Jpx] = coco_ezDFDX('f(o,d,x)', prob, data, ...
    @(o,d,u) temp(o,d,u,w), [p; x]);
  J = [Jpx, -pd.cont, Px; pd.zrs1, 2*v(pd.sidx{1})', pd.zrs2];
end

end

function Px = hspo_proj(prob, data, u)
%HSPO_PROJ   Compute Jacobian of event projections
%
% For each segment, extract Jacobian of time-T map and premultiply by
% saltation matrix (correcting for difference in time-of-flight to event
% surface, and the imposition of the reset) to obtain transfer matrix.
%
% P = HSPO_PROJ(PROB, DATA, U)
%
% P    - Cell array of transfer matrices.
% PROB - Continuation problem structure.
% DATA - hspo_mult_eigs_bddat function data.
% U    - Array of end points of segments at t=1 and problem parameters.

pd = data.hspo_pd;

Px = cell(1,pd.nsegs);
p = u(pd.p_idx); % Problem parameters
for i=1:pd.nsegs
  fdata = coco_get_func_data(prob, pd.cids{i}, 'data');
  x     = u(pd.p_idx(end) + pd.sidx{i}); % Segment end point at t=1
  fs    = fdata.ode_F(fdata, 0, x, p);
  dhdx  = data.event_DFDX(data, x, p, data.events{i});
  dgdx  = data.reset_DFDX(data, x, p, data.resets{i});
  Px{i} = dgdx - dgdx*fs*dhdx/(dhdx*fs);
end
Px = blkdiag(Px{:});

end

function [data, chart, y] = double(prob, data, chart, u) %#ok<INUSL>

pr = data.pr;
pd = pr.hspo_pd;
segs = pd.segs;

nsegs = pd.nsegs;
t0   = cell(1,nsegs);
x10  = cell(1,nsegs);
x20  = cell(1,nsegs);
for i=1:nsegs
  maps = segs{i}.maps;
  mesh = segs{i}.mesh;
  x = u(maps.xbp_idx);
  T = u(maps.T_idx);
  v = u(maps.var_idx);
  
  t0{i}   = mesh.tbp(maps.tbp_idx)*T;
  x1      = reshape(x+0.01*v, maps.xbp_shp)'; % First loop
  x10{i}  = x1(maps.tbp_idx,:);
  x2      = reshape(x-0.01*v, maps.xbp_shp)'; % Second loop
  x20{i}  = x2(maps.tbp_idx,:);
end
p = u(maps.p_idx);
cdata = coco_get_chart_data(chart, pd.fid);
cdata.pd.x0     = [x10 x20];
cdata.pd.t0     = [t0 t0];
cdata.pd.p0     = p;
cdata.pd.modes  = [pr.modes  pr.modes];
cdata.pd.events = [pr.events pr.events];
cdata.pd.resets = [pr.resets pr.resets];
chart = coco_set_chart_data(chart, pd.fid, cdata);

y = [];

end

function [prob, status, xtr] = remesh(prob, data, chart, old_u, old_V) %#ok<INUSD>

pd = data.hspo_pd;

fdata = coco_get_func_data(prob, data.bvid, 'data');
nsegs = fdata.nsegs;

segs = cell(1,nsegs);
off  = 0;
for i=1:nsegs
  cdata = coco_get_func_data(prob, pd.vids{i}, 'data');
  mesh  = cdata.coll_seg.mesh;
  var   = cdata.coll_var;
  segs{i}.maps.tbp_idx = var.tbp_idx;
  segs{i}.maps.xbp_idx = off + var.xbp_idx;
  segs{i}.maps.xbp_shp = var.xbp_shp;
  segs{i}.maps.T_idx = off + var.T_idx;
  segs{i}.maps.p_idx = off + var.p_idx;
  segs{i}.maps.var_idx = off + var.v_idx;
  segs{i}.mesh.tbp = mesh.tbp;
  off   = off + var.v_idx(end);
end
pd.segs = segs;
data.hspo_pd = pd;

uidx = [];
for i=1:nsegs
  uidx_var = coco_get_func_data(prob, pd.vids{i}, 'uidx');
  uidx = [uidx; uidx_var]; %#ok<AGROW>
end

xtr    = [];
prob   = coco_change_func(prob, data, 'uidx', uidx);
status = 'success';

end
