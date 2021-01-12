coco_use_recipes_toolbox alg_v3 % add the alg_v3 toolbox to the search path

%% Example 4.3

% The continuation problem encoded below corresponds to a single zero
% function in terms of two continuation variables, a single monitor
% function that evaluates to the second continuation variables, and a
% corresponding inactive continuation parameter 'y'. Its dimensional
% deficit equals 0. The call to the coco entry-point function indicates a
% desired manifold dimension of 1. To this end, the continuation parameter
% 'y' is released and allowed to vary during continuation.

funs = {@(x,p) x^2+(p-1)^2-1, @(x,p) 2*x, @(x,p) 2*(p-1)};
prob = alg_construct_eqn(funs{:}, 1, 'y', 1.1);
coco(prob, 'run', [], 1, 'y', [0.1 5]);

data  = coco_read_solution('', 'run', 3)
data  = coco_read_solution('alg', 'run', 3)
[x p] = alg_read_solution('run', 4)

coco_use_recipes_toolbox % remove the alg_v3 toolbox from the search path
