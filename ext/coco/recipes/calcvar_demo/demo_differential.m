coco_use_recipes_toolbox coll_v1 bvp_v1 % add the coll_v1 and bvp_v1 toolboxes to the search path

% The continuation problem encoded below corresponds to a family of 101
% zero functions (80 collocation conditions, 18 continuity conditions, 3
% boundary conditions) in terms of 102 continuation variables (100
% basepoint values, one interval length, and one problem parameter), a
% single monitor function that evaluates to the problem parameter, and the
% corresponding inactive continuation parameter 'Y'. Its dimensional
% deficit equals 0. The call to the coco entry-point function indicates a
% desired manifold dimension of 1. To this end, the continuation parameter
% 'Y' is released and allowed to vary during continuation. Special points
% associated with specific values of 'Y' and denoted by 'UZ' are detected
% and located using a Newton method. The family of solutions correspond to
% a parameterization of solutions to the Euler-Lagrange equations of the
% catenary problem, in terms of continuous, piecewise polynomial functions
% (f(t), f'(t)) that satisfy the boundary conditions f(0)=1 and f(1)=Y.

t0 = (0:0.01:1)';
x0 = [cosh(t0) sinh(t0)];
prob = coco_prob();
prob = bvp_isol2seg(prob, '', @caty_ode, t0, x0, 'Y', cosh(1), ...
    @caty_bc, @caty_bc_DFDX);
prob = coco_add_event(prob, 'UZ', 'Y', [0.6 0.7 1]);
coco(prob, 'run_diff', [], 1, 'Y', [0 2]);

coco_use_recipes_toolbox % remove the coll_v1 and bvp_v1 toolboxes from the search path
