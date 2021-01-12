function prob = huxley_close_het(prob, devs)
%HUXLEY_CLOSE_HET   Glue together two instances of the 'coll' toolbox.
%
% PROB = HUXLEY_CLOSE_HET(PROB, DEVS)
%
% PROB - continuation problem structure.
% DEVS - array of deviations of heteroclinic trajectory end points from
%        corresponding equilibria.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: huxley_close_het.m 2862 2015-07-26 22:11:32Z hdankowicz $

[data1, uidx1] = coco_get_func_data(prob, 'huxley1.coll', ...
  'data', 'uidx'); % Extract toolbox data and context-dependent index set
[data2, uidx2] = coco_get_func_data(prob, 'huxley2.coll', ...
  'data', 'uidx'); % Extract toolbox data and context-dependent index set

maps1 = data1.coll_seg.maps;
maps2 = data2.coll_seg.maps;

prob = coco_add_glue(prob, 'shared', uidx1(maps1.p_idx), ...
  uidx2(maps2.p_idx)); % Glue parameters

prob = coco_add_func(prob, 'bcs', @huxley_bcs, [], 'zero', 'uidx', ...
  [uidx1(maps1.x0_idx); uidx2(maps2.x1_idx); uidx1(maps1.p_idx)], ...
  'u0', devs); % Apply boundary conditions
uidx = coco_get_func_data(prob, 'bcs', 'uidx');
data.dev_idx = [numel(uidx)-1; numel(uidx)]; % store index array
prob = coco_add_slot(prob, 'bcs', @coco_save_data, data, 'save_full');

prob = coco_add_glue(prob, 'gap', uidx1(maps1.x1_idx(2)), ...
  uidx2(maps2.x0_idx(2)), 'gap', 'inactive'); % Monitor vertical gap

prob = coco_add_pars(prob, 'pars', ...
  [uidx1(maps1.p_idx); uidx(data.dev_idx); ...
  uidx1(maps1.x1_idx(1)); uidx2(maps2.x0_idx(1))], ...
  {'p1' 'p2' 'dev1' 'dev2' 'y11e' 'y21e'});

end
