function prob = bvp_close_segs(prob, data)
%BVP_CLOSE_SEGS   Append an instance of 'bvp' to problem.
%
% PROB = BVP_CLOSE_SEGS(PROB, DATA)
%
% PROB - Continuation problem structure.
% DATA - Toolbox data structure.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: bvp_close_segs.m 2839 2015-03-05 17:09:01Z fschild $

data = init_data(prob, data);
bc   = data.bvp_bc;
fid  = bc.fid;

if ~isempty(data.bc_update) % Optional inclusion of update slot
  prob = coco_add_slot(prob, fid, @update, data, 'update');
end
nsegs  = data.nsegs;
T0_idx = zeros(nsegs,1);
T_idx  = zeros(nsegs,1);
x0_idx = zeros(bc.cdim,1);
x1_idx = zeros(bc.cdim,1);
s_idx  = cell(1, nsegs);
off = 0;
for i=1:nsegs
  [fdata, uidx] = coco_get_func_data(prob, data.cids{i}, 'data', 'uidx');
  dim  = fdata.xdim;
  maps = fdata.coll_seg.maps;
  T0_idx(i) = uidx(maps.T0_idx);
  T_idx(i)  = uidx(maps.T_idx);
  x0_idx(off + (1:dim)') = uidx(maps.x0_idx);
  x1_idx(off + (1:dim)') = uidx(maps.x1_idx);
  s_idx{i} = uidx(maps.p_idx);
  off  = off + dim;
end
uidx = [T0_idx; T_idx; x0_idx; x1_idx; s_idx{1}];
switch data.basemode % Add boundary conditions zero functions
  case 1
    prob = coco_add_func(prob, fid, @FDF1, data, 'zero', 'uidx', uidx, ...
      'F+DF');
  case 2
    prob = coco_add_func(prob, fid, @FDF2, data, 'zero', 'uidx', uidx, ...
      'F+DF');
  otherwise
    prob = coco_add_func(prob, fid, @FDF2, data, 'zero', 'uidx', uidx);
end
for i=2:nsegs % Glue redundant copies of problem parameters
  sfid  = coco_get_id(fid, sprintf('shared%d', i-1));
  prob = coco_add_glue(prob, sfid, s_idx{1}, s_idx{i});
end
if ~isempty(data.pnames) % Optional monitor functions
  pfid  = coco_get_id(fid, 'pars');
  prob = coco_add_pars(prob, pfid, s_idx{1}, data.pnames);
end
prob = coco_add_slot(prob, fid, @coco_save_data, data, 'save_full');

end

function data = init_data(prob, data)
%INIT_DATA   Initialize toolbox data for an instance of 'bvp'.

nsegs = data.nsegs;
cdim  = 0;
for i=1:nsegs
  fdata = coco_get_func_data(prob, data.cids{i}, 'data');
  cdim  = cdim + fdata.xdim; % Track total dimension
end
bc.cdim = cdim;

bc.T0_idx = (1:nsegs)';
bc.T_idx  = nsegs + (1:nsegs)';
bc.x0_idx = 2*nsegs + (1:cdim)';
bc.x1_idx = 2*nsegs + cdim + (1:cdim)';
bc.p_idx  = 2*nsegs + 2*cdim + (1:fdata.pdim)';
bc.nargs  = (nargin(data.fhan)==5);

bc.fid    = coco_get_id(data.oid, 'bvp');

data.bvp_bc = bc;
data.no_save = [ data.no_save ...
  { 'bc.T0_idx' 'bc.T_idx' 'bc.x0_idx' 'bc.x1_idx' ...
  'bc.p_idx' 'bc.cdim' } ];

end

function [data, y, J] = FDF1(prob, data, u) %#ok<INUSL>
%FDF1   COCO-compatible wrapper to boundary conditions, including Jacobian

pr = data.pr;
bc = pr.bvp_bc;
  
T0 = u(bc.T0_idx);
T  = u(bc.T_idx);
x0 = u(bc.x0_idx);
x1 = u(bc.x1_idx);
p  = u(bc.p_idx);

args = {T0, T, x0, x1, p};

if nargout==2
  y = pr.fhan(pr.bc_data, args{bc.nargs+1:end});
else
  [y, J] = pr.fhan(pr.bc_data, args{bc.nargs+1:end});
  J = [zeros(numel(y), bc.nargs*numel(T0)) J];
end

end

function [data, y, J] = FDF2(prob, data, u) %#ok<INUSL>
%FDF2   COCO-compatible wrapper to boundary conditions and optional Jacobian

pr = data.pr;
bc = pr.bvp_bc;
  
T0 = u(bc.T0_idx);
T  = u(bc.T_idx);
x0 = u(bc.x0_idx);
x1 = u(bc.x1_idx);
p  = u(bc.p_idx);

args = {T0, T, x0, x1, p};

y  = pr.fhan(pr.bc_data, args{bc.nargs+1:end});

if nargout>=3
  J  = [zeros(numel(y), bc.nargs*numel(T0)) ...
    pr.dfdxhan(pr.bc_data, args{bc.nargs+1:end})];
end

end

function data = update(prob, data, cseg, varargin)
%UPDATE   COCO-compatible wrapper to boundary conditions update function.

pr = data.pr;
bc = pr.bvp_bc;

uidx = coco_get_func_data(prob, bc.fid, 'uidx');
u    = cseg.src_chart.x(uidx);
T0   = u(bc.T0_idx);
T    = u(bc.T_idx);
x0   = u(bc.x0_idx);
x1   = u(bc.x1_idx);
p    = u(bc.p_idx);
args = {T0, T, x0, x1, p};

data.bc_data = pr.bc_update(pr.bc_data, args{bc.nargs+1:end});

end
