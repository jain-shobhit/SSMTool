function [prob, data] = hspo_TR_construct_eqn(prob, data, sol)
%HSPO_TR_CONSTRUCT_EQN   Add torus bifurcation zero problem.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_PD_construct_eqn.m 2849 2015-05-17 20:32:46Z hdankowicz $

data = init_data(prob, data);
tr   = data.hspo_tr;
fid  = tr.fid;

xidx = [];
vidx = [];
widx = [];
for i=1:tr.nsegs
  [fdata, uidx_var] = coco_get_func_data(prob, tr.vids{i}, 'data', 'uidx');
  maps = fdata.coll_seg.maps;
  var  = fdata.coll_var;
  xidx = [xidx; uidx_var(maps.x1_idx)]; %#ok<AGROW>
  vidx = [vidx; uidx_var(var.v0_idx)]; %#ok<AGROW>
  widx = [widx; uidx_var(var.v1_idx)]; %#ok<AGROW>
end
uidx = [uidx_var(maps.p_idx); xidx; vidx(:); widx(:)];

prob = coco_add_func(prob, fid, @FDF, data, 'zero', 'uidx', uidx, ...
  'u0', sol.tr.u0, 't0', sol.tr.t0, 'requires', tr.vids, 'F+dF');

end

function data = init_data(prob, data)
%INIT_DATA   Initialize data for torus bifurcation zero problem.

po = data.hspo_orb;
cdim = po.cdim;
pdim = po.pdim;
tr.p_idx = 1:pdim;
tr.x_idx = pdim + (1:cdim);
tr.v_idx = pdim + cdim + (1:2*cdim);
tr.w_idx = pdim + 3*cdim + (1:2*cdim);
tr.a_idx = pdim + 5*cdim + 1;
tr.b_idx = pdim + 5*cdim + 2;

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
tr.nsegs  = nsegs;
tr.cids   = fdata.cids;
tr.vids   = vids;
tr.sidx   = sidx;

frstdim = po.dim(1);
lastdim = po.dim(end);
contmat = circshift(eye(cdim), -frstdim);
tr.cont = contmat(1:cdim-lastdim,:);
tr.zeros = zeros(cdim-lastdim,cdim);
tr.zero_row = zeros(frstdim, cdim-frstdim);
tr.eye  = eye(frstdim);

tr.perm_a = [0*tr.cont tr.zeros; ...
  tr.eye tr.zero_row -0*tr.eye tr.zero_row;
  tr.zeros 0*tr.cont; ...
  0*tr.eye tr.zero_row tr.eye tr.zero_row];
tr.perm_b = [0*tr.cont tr.zeros; ...
  0*tr.eye tr.zero_row -tr.eye tr.zero_row;
  tr.zeros 0*tr.cont; ...
  tr.eye tr.zero_row 0*tr.eye tr.zero_row];

tr.zrs1 = zeros(1,pdim+cdim);
tr.zrs2 = zeros(1,cdim-frstdim);
tr.zrs3 = zeros(1,cdim-frstdim+2*cdim+2);
tr.zrs4 = zeros(1,pdim+5*cdim);

tr.fid   = coco_get_id(data.oid, 'hspo.TR');

data.hspo_tr   = tr;
data.no_save = [ data.no_save { 'hspo_pd.rows' 'hspo_pd.cols'} ];

end

function [data, y, J] = FDF(prob, data, u)
%FDF   Torus bifurcation system and linearization

function [data, f] = temp(prob, data, u, w)
  P = hspo_proj(prob, data, u);
  f = blkdiag(P,P)*w;
end

pr = data.pr;
tr = pr.hspo_tr;
po = pr.hspo_orb;

p = u(tr.p_idx);
x = u(tr.x_idx);
v = u(tr.v_idx);
w = u(tr.w_idx);
a = u(tr.a_idx);
b = u(tr.b_idx);

Px = hspo_proj(prob, pr, [p; x]);
Px = blkdiag(Px,Px);

perm = [tr.cont tr.zeros; ...
  a*tr.eye tr.zero_row -b*tr.eye tr.zero_row;
  tr.zeros tr.cont; ...
  b*tr.eye tr.zero_row a*tr.eye tr.zero_row];
cont = Px*w - perm*v;
unit = [v(tr.sidx{1})' * v(tr.sidx{1}) + ...
  v(po.cdim+tr.sidx{1})' * v(po.cdim+tr.sidx{1}) - 1; ...
  v(tr.sidx{1})' * v(po.cdim+tr.sidx{1}); a^2 + b^2 - 1];

y = [cont(:); unit];

if nargout>=3
  [data, Jpx] = coco_ezDFDX('f(o,d,x)', prob, data, ...
    @(o,d,u) temp(o,d,u,w), [p; x]);
  J = [Jpx, -perm, Px, -tr.perm_a*v -tr.perm_b*v; ...
    tr.zrs1, 2*v(tr.sidx{1})', tr.zrs2, 2*v(po.cdim+tr.sidx{1})', tr.zrs3; ...
    tr.zrs1, v(po.cdim+tr.sidx{1})', tr.zrs2, v(tr.sidx{1})', tr.zrs3; ...
    tr.zrs4, 2*a, 2*b];
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

tr = data.hspo_tr;

Px = cell(1,tr.nsegs);
p = u(tr.p_idx); % Problem parameters
for i=1:tr.nsegs
  fdata = coco_get_func_data(prob, tr.cids{i}, 'data');
  x     = u(tr.p_idx(end)+tr.sidx{i}); % Segment end point at t=1
  fs    = fdata.ode_F(fdata, 0, x, p);
  dhdx  = data.event_DFDX(data, x, p, data.events{i});
  dgdx  = data.reset_DFDX(data, x, p, data.resets{i});
  Px{i} = dgdx - dgdx*fs*dhdx/(dhdx*fs);
end
Px = blkdiag(Px{:});

end
