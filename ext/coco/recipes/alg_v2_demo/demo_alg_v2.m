coco_use_recipes_toolbox alg_v2 % add the alg_v2 toolbox to the search path

%% Example 4.2

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

%% Example 4.2 without explicit Jacobian with respect to problem parameters

% The continuation problem encoded below is identical to that above, but
% the Jacobian with respect to problem parameters is approximated using
% numerical differentiation.

funs = {@(x,p) x^2+(p-1)^2-1, @(x,p) 2*x, []};
prob = alg_construct_eqn(funs{:}, 1, 'y', 1.1);
coco(prob, 'run', [], 1, 'y', [0.1 5]);

coco_use_recipes_toolbox % remove the alg_v2 toolbox from the search path
