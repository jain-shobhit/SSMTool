coco_use_recipes_toolbox atlas2d_v6 bvp_v1 coll_v1 % Add atlas2d_v6 atlas algorithm and bvp_v1 and coll_v1 toolboxes to search path

% The continuation problem structure encoded below corresponds to 151 zero
% functions (120 collocation conditions, 27 continuity conditions, and four
% boundary conditions) in terms of 154 continuation variables (150
% basepoint values, one interval length, and three problem parameters), a
% family of three monitor functions that evaluate to the problem
% parameters, and the corresponding inactive continuation parameters 'om',
% 'ro', and 'eps'. Its dimensional deficit equals 0. For 'eps'=0, there
% exists a periodic orbit for the reduced system, for which the rotation
% number is approximately three. A two-dimensional family of such 1:3
% phase-locked periodic orbits is obtained by releasing 'ro' and 'eps' and
% allowing these to vary during continuation.

p0 = [3.5; 0.35; 0];
[t x0] = ode45(@(t,x) lang(x, p0), [0 5.3], [0.3; 0; 0.4]); % Approximate periodic orbit

prob = coco_prob();
prob = bvp_isol2seg(prob, '', @lang, @lang_DFDX, @lang_DFDP,  ...
  t, x0, {'om' 'ro' 'eps'}, p0, @po_bc, @po_bc_DFDX);
prob = coco_set(prob, 'cont', 'atlas', @atlas_2d_min.create);
prob = coco_set(prob, 'cont', 'h', .25, 'PtMX', 5000, 'NPR', 500);
coco(prob, 'run1', [], 2, {'ro' 'eps'}, {[] [-0.5 0.5]});

coco_use_recipes_toolbox % Remove atlas2d_v6 atlas algorithm and bvp_v1 and coll_v1 toolboxes from search path
