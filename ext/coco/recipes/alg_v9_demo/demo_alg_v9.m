coco_use_recipes_toolbox alg_v9 atlas2d_v6 % add the alg_v9 toolbox and atlas2d_v6 atlas algorithm to the search path

%% Example 17.1

% The continuation problem encoded below consists of two zero functions in
% terms of four continuation variables (two problem variables and two
% problem parameters), four monitor functions that evaluate to the problem
% variables and parameters, respectively, and the corresponding inactive
% continuation parameters 'p1', 'p2', 'x', and 'y'. Its dimensional deficit
% is -2. The call to the coco entry-point function indicates a desired
% manifold dimensionality of 2. To this end, all four continuation
% parameters are released and allowed to vary during continuation. By
% default Hopf bifurcation and fold point detection are enabled and the
% associated special points are denoted by the HB and FO point types. The
% algorithm does not distinguish between Hopf bifurcations and neutral
% saddle points.

alg_args = {@popul, [1.76; 1.52], {'p1' 'p2'}, [0.3; 0.1]};
prob = alg_isol2eqn(coco_prob(), '', alg_args{:});
prob = coco_add_pars(prob, 'pars', [1 2], {'x', 'y'});
prob = coco_set(prob, 'cont', 'atlas', @atlas_2d_min.create);
prob = coco_set(prob, 'cont', 'h', 0.15, 'PtMX', 1000);
coco(prob, 'run', [], 2, {'p1' 'p2' 'x' 'y' }, ...
  {[0 0.5], [0 0.25], [0 10], [0 10]});

coco_use_recipes_toolbox % remove the alg_v9 toolbox and atlas2d_v6 atlas algorithm from the search path
