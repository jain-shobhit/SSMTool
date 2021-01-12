function prob = riess_close_het_2(prob, data)
%RIESS_CLOSE_HET_2   Append additional boundary conditions to heteroclinic problem.
%
% Append lin gap and phase conditions.
%
% PROB = RIESS_CLOSE_HET_2(PROB, DATA)
%
% PROB - Continuation problem structure.
% DATA - Problem data structure.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: riess_close_het_2.m 2839 2015-03-05 17:09:01Z fschild $

[data1 uidx1] = coco_get_func_data(prob, 'col1.coll', 'data', 'uidx'); % Data and index array for segment in unstable manifold of origin
[data2 uidx2] = coco_get_func_data(prob, 'col2.coll', 'data', 'uidx'); % Data and index array for segment in stable manifold of orbit

prob = coco_add_func(prob, 'gap', @lingap, data, 'inactive', ...
  'lingap', 'uidx', [uidx1(data1.x1_idx); uidx2(data2.x0_idx)]); % Append lin gap condition
prob = coco_add_func(prob, 'phase', @linphase, data, 'zero', ...
  'uidx', [uidx1(data1.x1_idx); uidx2(data2.x0_idx)]); % Append lin phase condition
prob = coco_add_slot(prob, 'riess_save_2', @coco_save_data, data, ...
  'save_full');

end
