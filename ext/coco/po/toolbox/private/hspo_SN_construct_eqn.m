function [prob, data] = hspo_SN_construct_eqn(prob, data, sol)
%HSPO_SN_CONSTRUCT_EQN   Add saddle-node zero problem.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_PD_construct_eqn.m 2849 2015-05-17 20:32:46Z hdankowicz $

data = init_data(prob, data);
sn   = data.hspo_sn;
fid  = sn.fid;

didx = [];
xidx = [];
vidx = [];
widx = [];
for i=1:sn.nsegs
  [fdata, uidx_var] = coco_get_func_data(prob, sn.vids{i}, 'data', 'uidx');
  didx = [didx; uidx_var]; %#ok<AGROW>
  maps = fdata.coll_seg.maps;
  var  = fdata.coll_var;
  xidx = [xidx; uidx_var(maps.x1_idx)]; %#ok<AGROW>
  vidx = [vidx; uidx_var(var.v0_idx)]; %#ok<AGROW>
  widx = [widx; uidx_var(var.v1_idx)]; %#ok<AGROW>
end
uidx = [uidx_var(maps.p_idx); xidx; vidx; widx];

prob = coco_add_func(prob, fid, @FDF, data, 'zero', 'uidx', uidx, ...
  'u0', sol.sn.u0, 't0', sol.sn.t0, 'requires', sn.vids, 'F+dF');

end

function data = init_data(prob, data)
%INIT_DATA   Initialize data for saddle-node zero problem.

po = data.hspo_orb;
cdim = po.cdim;
pdim = po.pdim;
sn.p_idx = 1:pdim;
sn.x_idx = pdim + (1:cdim);
sn.v_idx = pdim + cdim + (1:cdim);
sn.w_idx = pdim + 2*cdim + (1:cdim);

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
sn.nsegs  = nsegs;
sn.cids   = fdata.cids;
sn.vids   = vids;
sn.sidx   = sidx;

contmat = eye(cdim);
frstdim = po.dim(1);
sn.cont = circshift(contmat, -frstdim);
sn.zrs1 = zeros(1,pdim+cdim);
sn.zrs2 = zeros(1,2*cdim-frstdim);

sn.fid   = coco_get_id(data.oid, 'hspo.SN');

data.hspo_sn = sn;
data.no_save = [ data.no_save { 'hspo_sn.zrs1' 'hspo_sn.zrs2' } ];

end

function [data, y, J] = FDF(prob, data, u)
%FDF   Saddle-node system and linearization

function [data, f] = temp(prob, data, u, w)
  f = hspo_proj(prob, data, u)*w;
end

pr = data.pr;
sn = pr.hspo_sn;

p = u(sn.p_idx);
x = u(sn.x_idx);
v = u(sn.v_idx);
w = u(sn.w_idx);

Px = hspo_proj(prob, pr, [p; x]);

cont = Px*w - sn.cont*v;
unit = v(sn.sidx{1})'*v(sn.sidx{1}) - 1;

y = [cont(:); unit];

if nargout>=3
  [data, Jpx] = coco_ezDFDX('f(o,d,x)', prob, data, ...
    @(o,d,u) temp(o,d,u,w), [p; x]);
  J = [Jpx, -sn.cont, Px; sn.zrs1, 2*v(sn.sidx{1})', sn.zrs2];
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

sn = data.hspo_sn;

Px = cell(1,sn.nsegs);
p = u(sn.p_idx); % Problem parameters
for i=1:sn.nsegs
  fdata = coco_get_func_data(prob, sn.cids{i}, 'data');
  x     = u(sn.p_idx(end) + sn.sidx{i}); % Segment end point at t=1
  fs    = fdata.ode_F(fdata, 0, x, p);
  dhdx  = data.event_DFDX(data, x, p, data.events{i});
  dgdx  = data.reset_DFDX(data, x, p, data.resets{i});
  Px{i} = dgdx - dgdx*fs*dhdx/(dhdx*fs);
end
Px = blkdiag(Px{:});

end
