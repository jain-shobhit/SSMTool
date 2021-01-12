coco_use_recipes_toolbox alg_v5 compalg_v2 % add the alg_v5 and compalg_v2 toolboxes to the search path

%% Example 5.3

% The continuation problem encoded below corresponds to a single zero
% function in terms of two continuation variables, a single monitor function that
% evaluates to the problem parameter, and a corresponding inactive
% continuation parameter 'y'. Its dimensional deficit equals 0. Each call to
% the coco entry-point function indicates a desired manifold dimension of
% 1. To this end, the continuation parameter 'y' is released and allowed to
% vary during continuation.

fun = @(x,p) x^2+(p-1)^2-1;
alg_args = {fun, 0.9, 'y', 1};
coco('run1', @alg_isol2eqn, alg_args{:}, 1, 'y', [0.5 1.5]);
coco('run2', @alg_sol2eqn, 'run1', 4, 1, 'y', [0.1 1.9]);

%% Example 5.4

% The continuation problem encoded below corresponds to a family of three
% zero function (two instances of 'alg' and one gluing condition) in terms
% of four continuation variables (two problem variables and two copies of
% the problem parameter), a single monitor function that evaluates to the
% problem parameter, and a corresponding inactive continuation parameter
% 'y'. Its dimensional deficit equals 0. Each call to the coco entry-point
% function indicates a desired manifold dimension of 1. To this end, the
% continuation parameter 'y' is released and allowed to vary during
% continuation.

fun = @(x,p) x^2+(p-1)^2-1;
alg1_args = {fun, 1.1, 1};
alg2_args = {fun, 0.9, 1};
compalg_args = [alg1_args, alg2_args, 'y'];
coco('run3', @compalg_isol2sys, compalg_args{:}, 1, 'y', [0.5 1.5]);
coco('run4', @compalg_sol2sys, 'run3', 2, 1, 'y', [0.5 1.5]);

sol = compalg_read_solution('', 'run4', 2)

coco_use_recipes_toolbox % remove the alg_v5 and compalg_v2 toolboxes from the search path
