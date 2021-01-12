function prob = bvp_construct_adjt(prob, data, sol)
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
  up_fid = coco_get_id(fid, 'up');
  prob = coco_add_slot(prob, up_fid, @update, data, 'update');
end
nsegs  = data.nsegs;
T0_idx = zeros(nsegs,1);
T_idx  = zeros(nsegs,1);
x0_idx = zeros(bc.cdim,1);
x1_idx = zeros(bc.cdim,1);
s_idx  = cell(1, nsegs);
off = 0;
for i=1:nsegs
  [fdata, aidx] = coco_get_adjt_data(prob, data.cids{i}, 'data', 'axidx');
  dim = fdata.xdim;
  opt = fdata.coll_opt;
  T0_idx(i) = aidx(opt.T0_idx);
  T_idx(i)  = aidx(opt.T_idx);
  x0_idx(off + (1:dim)') = aidx(opt.x0_idx);
  x1_idx(off + (1:dim)') = aidx(opt.x1_idx);
  s_idx{i} = aidx(opt.p_idx);
  off  = off + dim;
end
aidx = [T0_idx; T_idx; x0_idx; x1_idx; s_idx{1}];
prob = coco_add_adjt(prob, fid, @F, @DF, data, 'aidx', aidx, ...
  'l0', sol.l0, 'tl0', sol.tl0);
for i=2:nsegs
  shid = sprintf('shared%d', i-1);
  sfid = coco_get_id(fid, shid);
  prob = coco_add_adjt(prob, sfid, 'aidx', [s_idx{1}; s_idx{i}], ...
    'l0', sol.([shid, '_l0']), 'tl0', sol.([shid, '_tl0']));
end
if ~isempty(data.pnames)
  pfid   = coco_get_id(fid, 'pars');
  dnames = coco_get_id('d', data.pnames);
  prob   = coco_add_adjt(prob, pfid, dnames, 'aidx', s_idx{1}, ...
    'l0', sol.pars_l0, 'tl0', sol.pars_tl0);
end

end

function data = init_data(prob, data)
%INIT_DATA   Initialize toolbox data for an instance of 'bvp'.

bc = data.bvp_bc;

fid  = coco_get_id(data.oid, 'bvp');
fidx = coco_get_func_data(prob, fid, 'fidx');

rownum = numel(fidx);
colnum = (2-bc.nargs)*data.nsegs + 2*bc.cdim + numel(bc.p_idx);
opt.dfdxrows = repmat(1:rownum, [1 colnum]);
opt.dfdxcols = repmat(bc.nargs*data.nsegs+1:bc.p_idx(end), [rownum 1]);

opt.dfdxdxrows = repmat(1:rownum, [1 colnum^2]);
cols = repmat(1:colnum, [rownum 1]);
opt.dfdxdxcols = repmat(cols(:), [1 colnum]) + ...
  repmat(bc.nargs*data.nsegs:colnum:(colnum-1)*colnum+bc.nargs*data.nsegs, ...
  [rownum*colnum 1]);

opt.fid = fid;
data.bvp_opt = opt;

end

function [data, J] = F(prob, data, u) %#ok<INUSL>
%F   COCO-compatible wrapper to adjoint boundary conditions

pr  = data.pr;
bc  = pr.bvp_bc;
opt = pr.bvp_opt;
  
T0 = u(bc.T0_idx);
T  = u(bc.T_idx);
x0 = u(bc.x0_idx);
x1 = u(bc.x1_idx);
p  = u(bc.p_idx);

args = {T0, T, x0, x1, p};
dfdx = pr.dfdxhan(pr.bc_data, args{bc.nargs+1:end});

J = sparse(opt.dfdxrows, opt.dfdxcols, dfdx(:));

end

function [data, dJ] = DF(prob, data, u) %#ok<INUSL>
%DF   COCO-compatible wrapper to Jacobian of adjoint to boundary conditions

pr  = data.pr;
bc  = pr.bvp_bc;
opt = pr.bvp_opt;
  
T0 = u(bc.T0_idx);
T  = u(bc.T_idx);
x0 = u(bc.x0_idx);
x1 = u(bc.x1_idx);
p  = u(bc.p_idx);

args = {T0, T, x0, x1, p};

dfdxdx = pr.dfdxdxhan(pr.bc_data, args{bc.nargs+1:end});

dJ = sparse(opt.dfdxdxrows, opt.dfdxdxcols, dfdxdx(:));

end

function data = update(prob, data, cseg, varargin) %#ok<INUSD>
%UPDATE   COCO-compatible wrapper to boundary conditions update function.

fdata = coco_get_func_data(prob, data.bvp_opt.fid, 'data');
data.bc_data = fdata.bc_data;

end
