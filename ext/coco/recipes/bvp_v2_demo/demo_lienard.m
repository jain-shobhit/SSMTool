coco_use_recipes_toolbox coll_v1 bvp_v2 % add the coll_v1 and bvp_v2 toolboxes to the search path

%% Example 8.5

% The continuation problem structure encoded below corresponds to a family
% of 101 zero functions (80 collocation conditions, 18 continuity
% conditions, two periodicity conditions, and one Poincare section
% condition) in terms of 102 continuation variables (100 basepoint values,
% one interval length, and one problem parameter), a single monitor
% function that evaluates to the problem parameter, and a corresponding
% inactive continuation parameter 'p'. Its dimensional deficit is 0. The
% calls to the coco entry-point function indicate a desired manifold
% dimensionality of 1. To this end, the continuation parameter 'p' is
% released and allowed to vary during continuation.

p0 = 1;
x0 = [0.4; -1.2];
f  = @(t,x) lienard(x, p0);
[t0 x0] = ode45(f, [0 6.7], x0);
coll_args = {@lienard, t0, x0, 'p', p0};

data = struct();
data.fhan = @lienard; % Used by per_bc_update for computing normal vector to section
data = per_bc_update(data, [], x0(1,:)', [], p0); % Initialize Poincare section

bvp_args = {@per_bc, @per_bc_DFDX, data, @per_bc_update};
prob = bvp_isol2seg(coco_prob(), '', coll_args{:}, bvp_args{:});
prob = coco_set(prob, 'cont', 'NPR', 1); % More frequent outputs than in Recipes for Continuation
coco(prob, 'run_moving', [], 1, 'p', [-1 1]); % Run with updates to section

bvp_args = {@per_bc, @per_bc_DFDX, data};
prob = bvp_isol2seg(coco_prob(), '', coll_args{:}, bvp_args{:});
prob = coco_set(prob, 'cont', 'NPR', 1); % More frequent outputs than in Recipes for Continuations
coco(prob, 'run_fixed', [], 1, 'p', [-1 1]); % Run without updates to section

coco_use_recipes_toolbox % remove the coll_v1 and bvp_v2 toolboxes from the search path
