function prob = doedel_close_het(prob, dev, th, lam)
%DOEDEL_CLOSE_HET   Glue together two instances of the 'coll' toolbox and two instances of the ‘ep’ toolbox.
%
% PROB = DOEDEL_CLOSE_HET(PROB, DEV, TH, LAM)
%
% PROB - continuation problem structure.
% DEV  - array of deviations of heteroclinic trajectory end points from
%        corresponding equilibria.
% TH   - angle relative to horizontal
% LAM  - eigenvalue

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: doedel_close_het.m 2896 2015-10-05 20:08:38Z hdankowicz $

% Extract toolbox data and context-dependent index sets
[data1, uidx1] = coco_get_func_data(prob, 'doedel1.coll', ...
  'data', 'uidx'); 
[data2, uidx2] = coco_get_func_data(prob, 'doedel2.coll', ...
  'data', 'uidx');
[data3, uidx3] = coco_get_func_data(prob, 'doedel3.ep', ...
  'data', 'uidx');
[data4, uidx4] = coco_get_func_data(prob, 'doedel4.ep.var', ...
  'data', 'uidx');

maps1 = data1.coll_seg.maps;
maps2 = data2.coll_seg.maps;
ep1   = data3.ep_eqn;
ep2   = data4.ep_var;

% Glue parameters
prob = coco_add_glue(prob, 'shared', ...
  [uidx1(maps1.p_idx); uidx1(maps1.p_idx); uidx1(maps1.p_idx)], ...
  [uidx2(maps2.p_idx); uidx3(ep1.p_idx); uidx4(ep2.p_idx)]);

% Apply eigenvalue condition
prob = coco_add_func(prob, 'evs', @doedel_evs, [], 'zero', ...
  'uidx', uidx4([ep2.v_idx; ep2.w_idx]), 'u0', lam);
uidx_evs = coco_get_func_data(prob, 'evs', 'uidx');

% Store context-independent index array for eigenvalue
data_evs.lam_idx  = numel(uidx_evs); 
prob = coco_add_slot(prob, 'evs', @coco_save_data, data_evs, ...
  'save_full');

% Apply boundary conditions
prob = coco_add_func(prob, 'bcs', @doedel_bcs, [], 'zero', 'uidx', ...
    [uidx1(maps1.x0_idx); uidx2(maps2.x1_idx); ...
    uidx3(ep1.x_idx); uidx4([ep2.x_idx; ep2.v_idx])], 'u0', [dev; th]);
uidx_bcs = coco_get_func_data(prob, 'bcs', 'uidx');

% Store context-independent index array for deviations and
% angle relative to horizontal
data_bcs.dev_idx = [numel(uidx_bcs)-2; numel(uidx_bcs)-1]; 
data_bcs.th_idx  = numel(uidx_bcs);
prob = coco_add_slot(prob, 'bcs', @coco_save_data, data_bcs, ...
  'save_full');

% Monitor horizontal gap
prob = coco_add_glue(prob, 'lin', uidx1(maps1.x1_idx(1)), ...
  uidx2(maps2.x0_idx(1)), 'gap', 'inactive');

% Introduce continuation parameters
prob = coco_add_pars(prob, 'pars', ...
  [uidx1(maps1.p_idx); uidx_bcs(data_bcs.dev_idx); ...
  uidx1(maps1.x1_idx(2)); uidx2(maps2.x0_idx(2))], ...
    {'p1' 'p2' 'dev1' 'dev2' 'y12e' 'y22e'});

end
