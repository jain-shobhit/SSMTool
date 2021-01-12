function prob = riess_close_het_1(prob, data, vec, lam, eps)
%RIESS_CLOSE_HET_1   Append boundary conditions to heteroclinic problem.
%
% Glue redundant copies of parameters, impose eigenspace and eigenvector
% conditions, and append hyperplane projection conditions.
%
% PROB = RIESS_CLOSE_HET_1(PROB, DATA, VEC, LAM, EPS)
%
% PROB - Continuation problem structure.
% DATA - Problem data structure.
% VEC  - Eigenvector.
% LAM  - Floquet multiplier.
% EPS  - Array of separations at end points from orbit and equilibrium.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: riess_close_het_1.m 2839 2015-03-05 17:09:01Z fschild $

[data1 uidx1] = coco_get_func_data(prob, 'col1.coll', 'data', 'uidx'); % Data and index array for segment in unstable manifold of origin
[data2 uidx2] = coco_get_func_data(prob, 'col2.coll', 'data', 'uidx'); % Data and index array for segment in stable manifold of orbit
[data3 uidx3] = coco_get_func_data(prob, 'po.seg.var', 'data', 'uidx'); % Data and index array for fundamental solution associated with orbit
[data4 uidx4] = coco_get_func_data(prob, data3.coll_id, 'data', 'uidx'); % Data and index array for periodic orbit segment

prob = coco_add_glue(prob, 'shared', ...
  [uidx1(data1.p_idx); uidx1(data1.p_idx)], ...
  [uidx2(data2.p_idx); uidx4(data4.p_idx)]); % Glue redundant copies of problem parameters

evsdata = struct('tbid', 'po.seg.var');
prob = coco_add_func(prob, 'evs', @var_evs, evsdata, 'zero', ...
  'uidx', uidx3(data3.ubp_idx), 'u0', [vec; lam]); % Append eigenvector conditions
uidx = coco_get_func_data(prob, 'evs', 'uidx');
data.vec_idx = uidx(end-3:end-1); % Store index arrays
data.lam_idx = uidx(end);         % Store index arrays

prob = coco_add_func(prob, 'bcs1', @eig_bcs, [], 'zero', 'uidx', ...
    [uidx1(data1.x0_idx); uidx2(data2.x1_idx); ...
    uidx4(data4.xbp_idx(end-data.dim+1:end)); ...
    uidx4(data4.p_idx(1:2)); data.vec_idx], 'x0', eps); % Append eigenspace conditions
uidx = coco_get_func_data(prob, 'bcs1', 'uidx');
data.eps_idx = uidx(end-1:end);   % Store index arrays

prob = coco_add_func(prob, 'bcs2', @proj_bcs, data, 'inactive', ...
  {'sg1' 'sg2'}, 'uidx', [uidx1(data1.x1_idx); uidx2(data2.x0_idx)]); % Append hyperplane conditions

prob = coco_add_pars(prob, 'pars', ...
  [data.eps_idx; uidx1(data1.T_idx); uidx2(data2.T_idx)], ...
  {'eps1' 'eps2' 'T1' 'T2'});

prob = coco_add_slot(prob, 'riess_save_1', @coco_save_data, data, ...
  'save_full');

end
