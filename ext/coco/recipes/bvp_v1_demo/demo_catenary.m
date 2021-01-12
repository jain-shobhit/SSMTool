coco_use_recipes_toolbox coll_v1 bvp_v1 % add the coll_v1 and bvp_v1 toolboxes to the search path

%% Example 8.1

% The continuation problem structure encoded below consists of a family of 101
% zero functions (80 collocation conditions, 18 continuity conditions, 3
% boundary conditions) in 102 continuation variables (100 basepoint values,
% 1 interval length, 1 problem parameter), a single monitor function, and
% the corresponding inactive continuation parameter 'Y'. Its dimensional
% deficit equals 0. The call to the coco entry-point function indicates a
% desired manifold dimension of 1. To this end, the continuation parameter
% 'Y' is released and allowed to vary during continuation.

t0 = [0; 1];
x0 = [1 0; 1 0];
Y0 = 1;
coll_args = {@catn, t0, x0, 'Y', Y0};
bvp_args  = [coll_args, {@catn_bc, @catn_bc_DFDX}];
prob      = bvp_isol2seg(coco_prob(), '', bvp_args{:});
coco(prob, 'catn', [], 1, 'Y', [0 3]);

coco_use_recipes_toolbox % remove the coll_v1 and bvp_v1 toolboxes from the search path
