coco_use_recipes_toolbox alg_v6 compalg_v2 % add the alg_v6 and compalg_v2 toolboxes to the search path

%% Example 5.5

% The continuation problem encoded below corresponds to a single zero
% function in terms of two continuation variables, a single monitor
% function that evaluates to the problem parameter, and a
% corresponding inactive continuation parameter 'mu'. Its dimensional
% deficit equals 0. The call to the coco entry-point function indicates a
% desired manifold dimension of 1. To this end, the continuation parameter
% 'mu' is released and allowed to vary during continuation. The absolute
% value of the problem variable is appended to the bifurcation data stored
% to disk during continuation.

fun = @(x,p) x^2+(p-1)^2-1;
alg_args = {fun, 1.1, 'mu', 1};
prob = coco_prob();
prob = coco_set(prob, 'alg', 'norm', true);
bd1 = coco(prob, 'run1', @alg_isol2eqn, alg_args{:}, 1, ...
  'mu', [0.5 1.5]);

%% Example 5.6

% The continuation problem encoded below corresponds to a family of three
% zero function in terms of four continuation variables (two problem
% variables and a redundant duplication of a single problem parameter), a
% single monitor function that evaluates to the problem parameter, and a
% corresponding inactive continuation parameter 'mu'. Its dimensional
% deficit equals 0. The call to the coco entry-point function indicates a
% desired manifold dimension of 1. To this end, the continuation parameter
% 'mu' is released and allowed to vary during continuation. The absolute
% value of the first (but not the second) problem variable is appended to
% the bifurcation data stored to disk during continuation.

fun = @(x,p) x^2+(p-1)^2-1;
alg1_args = {fun, 1.1, 1};
alg2_args = {fun, 0.9, 1};
compalg_args = [alg1_args, alg2_args, 'mu'];
prob = coco_prob();
prob = coco_set(prob, 'compalg.alg', 'norm', true);
prob = coco_set(prob, 'compalg.eqn2.alg', 'norm', false);
bd2 = coco(prob, 'run2', @compalg_isol2sys, compalg_args{:}, ...
  1, 'mu', [0.5 1.5]);

%% Restarted run [not in Recipes for Continuation]

% The continuation problem encoded below is identical to that above, but
% constructed from a stored solution from the previous run.

bd3 = coco(prob, 'run3', @compalg_sol2sys, 'run2', 2, ...
  1, 'mu', [0.5 1.5]);

coco_use_recipes_toolbox % remove the alg_v6 and compalg_v2 toolboxes from the search path
