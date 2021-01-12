coco_use_recipes_toolbox atlas1d_v6 % Add atlas1d_v6 atlas algorithm to search path

%% Example 16.1

% The continuation problem structure encoded below corresponds to a single
% zero function in two continuation variables, a family of two monitor
% functions that evaluate to the continuation variables, and the
% corresponding inactive continuation parameters 'x' and 'y'. Its
% dimensional deficit is -1. A one-dimensional family of solutions is
% obtained by releasing 'x' and 'y' and allowing these to vary during
% continuation. Changes in the sign of the atlas.test.FP and atlas.test.BP
% continuation parameters corresponds to special points with event type
% 'FP' and 'BP'.

prob = coco_add_func(coco_prob(), 'lemniscate', @lemniscate, ...
  [], 'zero', 'u0', [0.6; 0.4]);
prob = coco_add_pars(prob, '', [1 2], {'x' 'y'});
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
prob = coco_set(prob, 'cont', 'almax', 15);
cont_pars = {'x' 'y' 'atlas.test.FP' 'atlas.test.BP'};
coco(prob, 'lemniscate', [], 1, cont_pars);

coco_use_recipes_toolbox % Remove atlas1d_v6 atlas algorithm from search path
