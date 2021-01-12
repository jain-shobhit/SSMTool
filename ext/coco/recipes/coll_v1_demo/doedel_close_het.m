function prob = doedel_close_het(prob, eps, th, vec, lam)
%DOEDEL_CLOSE_HET   Glue together two instances of the 'coll' toolbox and two instances of the 'alg' toolbox.
%
% PROB = DOEDEL_CLOSE_HET(PROB, EPS, th, vec, lam)
%
% PROB - continuation problem structure.
% EPS  - array of deviations of heteroclinic trajectory end points from
%        corresponding equilibria.
% TH   - angle relative to horizontal
% VEC  - eigenvector
% LAM  - eigenvalue

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: doedel_close_het.m 2839 2015-03-05 17:09:01Z fschild $

[data1 uidx1] = coco_get_func_data(prob, 'doedel1.coll', ...
  'data', 'uidx'); % Extract toolbox data and context-dependent index set
[data2 uidx2] = coco_get_func_data(prob, 'doedel2.coll', ...
  'data', 'uidx'); % Extract toolbox data and context-dependent index set
[data3 uidx3] = coco_get_func_data(prob, 'doedel3.alg', ...
  'data', 'uidx'); % Extract toolbox data and context-dependent index set
[data4 uidx4] = coco_get_func_data(prob, 'doedel4.alg', ...
  'data', 'uidx'); % Extract toolbox data and context-dependent index set

prob = coco_add_glue(prob, 'shared', ...
  [uidx1(data1.p_idx); uidx1(data1.p_idx); uidx1(data1.p_idx)], ...
  [uidx2(data2.p_idx); uidx3(data3.p_idx); uidx4(data4.p_idx)]); % Glue parameters

prob = coco_add_func(prob, 'evs', @doedel_evs, [], 'zero', ...
  'uidx', [uidx4(data4.p_idx); uidx4(data4.x_idx)], 'u0', [vec; lam]); % Apply eigenspace conditions
uidx_evs = coco_get_func_data(prob, 'evs', 'uidx');
data_evs.vec_idx = [numel(uidx_evs)-2; numel(uidx_evs)-1]; % Store index array
data_evs.lam_idx  = numel(uidx_evs); % Store index array
prob = coco_add_slot(prob, 'evs', @coco_save_data, data_evs, ...
  'save_full');

prob = coco_add_func(prob, 'bcs', @doedel_bcs, [], 'zero', 'uidx', ...
    [uidx1(data1.x0_idx); uidx2(data2.x1_idx); ...
    uidx3(data3.x_idx); uidx4(data4.x_idx); ...
    uidx_evs(data_evs.vec_idx)], 'u0', [eps; th]); % Apply boundary conditions
uidx_bcs = coco_get_func_data(prob, 'bcs', 'uidx');
data_bcs.eps_idx = [numel(uidx_bcs)-2; numel(uidx_bcs)-1]; % Store index array
data_bcs.th_idx  = numel(uidx_bcs); % Store index array
prob = coco_add_slot(prob, 'bcs', @coco_save_data, data_bcs, ...
  'save_full');

prob = coco_add_glue(prob, 'lin', uidx1(data1.x1_idx(1)), ...
  uidx2(data2.x0_idx(1)), 'gap', 'inactive'); % Monitor horizontal gap

prob = coco_add_pars(prob, 'pars', ...
  [uidx1(data1.p_idx); uidx_bcs(data_bcs.eps_idx); ...
  uidx1(data1.x1_idx(2)); uidx2(data2.x0_idx(2))], ...
    {'p1' 'p2' 'eps1' 'eps2' 'y12e' 'y22e'});

end
