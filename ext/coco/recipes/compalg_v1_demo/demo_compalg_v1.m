coco_use_recipes_toolbox alg_v5 compalg_v1 % add the alg_v5 and compalg_v1 toolboxes to the search path

%% Example 5.1

% The continuation problem encoded below corresponds to a family of three
% zero function (two instances of 'alg' and one gluing condition) in terms
% of four continuation variables (two problem variables and two copies of
% the problem parameter), a single monitor function that evaluates to the
% problem parameter, and a corresponding inactive continuation parameter
% 'y'. Its dimensional deficit equals 0. The call to the coco entry-point
% function indicates a desired manifold dimension of 1. To this end, the
% continuation parameter 'y' is released and allowed to vary during
% continuation.

funs = {@(x,p) x^2+(p-1)^2-1, @(x,p) x^2+(p-1)^2-1};
prob = coco_prob();
prob = compalg_isol2sys(prob, 'sys1', funs, {1 -1}, 'y', 1);
coco(prob, 'run1', [], 1, 'y', [0.5 1.5]);

%% Example 5.2

% The compalg_read_solution solution extractor returns information
% compatible with reconstructing the continuation problem from saved data.

[sol data] = compalg_read_solution('sys1', 'run1', 4)

%% Restarted run [Not in Recipes for Continuation]

prob = coco_prob();
prob = compalg_sol2sys(prob, 'sys2', 'run1', 'sys1', 4);
coco(prob, 'run2', [], 1, 'y', [0.5 1.5]);

%% Test for single zero functions [Not in Recipes for Continuation]

prob = coco_prob();
prob = compalg_isol2sys(prob, '', funs{2}, -1, 'y', 1);
coco(prob, 'run3', [], 1, 'y', [0.5 1.5]);

coco_use_recipes_toolbox % remove the alg_v5 and compalg_v1 toolboxes from the search path
